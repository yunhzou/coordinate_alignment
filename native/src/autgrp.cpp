// Native replacement for ``pynauty.autgrp(graph)[0]`` on a fixed graph that is
// recoloured per call (rxn_core.matcher.canonical.atom_generators).
//
// The generators nauty emits -- and their order -- are a function of the
// graph, the initial ``lab``/``ptn`` (cell order and vertex order inside each
// cell), the option block and the workspace size.  Everything pynauty 2.8.8.1
// (nautywrap.c, dense nauty 2.8.8, the same sources as native/nauty) sets is
// reproduced here:
//
//   options   DEFAULTOPTIONS_GRAPH (getcanon 0, digraph FALSE, writeautoms
//             FALSE, writemarkers FALSE, defaultptn TRUE, cartesian FALSE,
//             linelength CONSOLWIDTH, no user procs, invarproc NULL,
//             tc_level 100, mininvarlevel 0, maxinvarlevel 1, invararg 0,
//             dispatch_graph, schreier FALSE), then create_nygraph /
//             graph_autgrp override: digraph FALSE, getcanon FALSE,
//             writeautoms FALSE, cartesian TRUE, linelength 0,
//             userautomproc = store_generator, defaultptn TRUE when the
//             vertex colouring is empty and FALSE otherwise.
//   m         (n + WORDSIZE - 1) / WORDSIZE   (create_nygraph: no_setwords)
//   worksize  WORKSPACE_FACTOR (66) * m setwords
//   adjacency every (x, y) of the adjacency dict sets bit y of row x and,
//             for an undirected graph, bit x of row y (ADDELEMENT both ways;
//             duplicates are idempotent, (x, x) sets the diagonal bit).
//   lab/ptn   set_partition: cells concatenated in list order, each cell in
//             the order the C loop iterates the Python set object (the caller
//             passes ``list(cell)``, which is that order); ptn is 1 inside a
//             cell and 0 at the last vertex of each non-empty cell.
//   output    a copy of perm[0..n-1] for every userautomproc call, in call
//             order (nauty.c: gca cases 1/2 and extra_autom).
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <stdexcept>
#include <algorithm>
#include <map>
#include <set>
#include <tuple>
#include <unordered_map>
#include <string>
#include <utility>
#include <vector>

extern "C" {
#include "nauty.h"
}

namespace py = pybind11;

namespace autgrp_native {

// nauty's userautomproc carries no user-data pointer; nautywrap.c uses a
// static global (GRAPH_PTR) for the same reason.
static thread_local std::vector<std::vector<int>>* g_sink = nullptr;

static void store_generator(int /*count*/, int* perm, int* /*orbits*/,
                            int /*numorbits*/, int /*stabvertex*/, int n) {
    if (g_sink != nullptr) g_sink->emplace_back(perm, perm + n);
}

struct AutGraph {
    int n = 0;
    int m = 0;
    std::vector<setword> adjacency;
    std::vector<int> lab, ptn, orbits;
    std::vector<setword> workspace;
    optionblk options;
    statsblk stats;

    AutGraph(int vertices, const std::vector<std::pair<int, int>>& edges) : n(vertices) {
        if (n < 1) throw std::invalid_argument("AutGraph needs at least one vertex");
        m = (n + WORDSIZE - 1) / WORDSIZE;
        nauty_check(WORDSIZE, m, n, NAUTYVERSIONID);
        adjacency.assign((size_t)n * m, 0);
        for (const auto& e : edges) {
            check_vertex(e.first);
            check_vertex(e.second);
            ADDELEMENT(GRAPHROW(adjacency.data(), e.first, m), e.second);
            ADDELEMENT(GRAPHROW(adjacency.data(), e.second, m), e.first);
        }
        lab.assign(n, 0);
        ptn.assign(n, 0);
        orbits.assign(n, 0);
        workspace.assign((size_t)66 * m, 0);
        static DEFAULTOPTIONS_GRAPH(defaults);
        options = defaults;
        options.digraph = FALSE;
        options.getcanon = FALSE;
        options.defaultptn = TRUE;
        options.writeautoms = FALSE;
        options.cartesian = TRUE;
        options.linelength = 0;
        options.userautomproc = store_generator;
    }

    void check_vertex(int v) const {
        if (v < 0 || v >= n)
            throw std::invalid_argument("vertex " + std::to_string(v) +
                                        " conflicts with n_vertices=" + std::to_string(n));
    }

    // ``cells``: the final pynauty ``vertex_coloring`` (after its completion
    // and single-cell dropping rules), each set already converted with
    // ``list(cell)``.  Empty -> nauty's default single-cell partition.
    std::vector<std::vector<int>> generators(const std::vector<std::vector<int>>& cells) {
        if (cells.empty()) {
            options.defaultptn = TRUE;
        } else {
            options.defaultptn = FALSE;
            std::vector<char> seen((size_t)n, 0);
            int k = 0;
            for (const auto& cell : cells) {
                for (int v : cell) {
                    check_vertex(v);
                    if (seen[(size_t)v])
                        throw std::invalid_argument("vertex " + std::to_string(v) +
                                                    " appears in two cells");
                    seen[(size_t)v] = 1;
                    lab[(size_t)k] = v;
                    ptn[(size_t)k] = 1;
                    ++k;
                }
                if (k > 0) ptn[(size_t)k - 1] = 0;
            }
            if (k != n) throw std::invalid_argument("colouring does not cover every vertex");
        }
        std::vector<std::vector<int>> out;
        g_sink = &out;
        nauty(adjacency.data(), lab.data(), ptn.data(), nullptr, orbits.data(), &options,
              &stats, workspace.data(), (int)workspace.size(), m, n, nullptr);
        g_sink = nullptr;
        return out;
    }
};

// Fixed topology, reusable atom-role colors and exact conditioned partitions.
// The partition and generator order match the Python implementation; only
// construction/caching move here. No symmetry group is expanded into elements.
class ConditionedAutGraph {
    AutGraph graph;
    std::vector<int> atoms, vertices;
    std::vector<std::string> base_colors;
    std::map<int,int> atom_index;
    std::vector<std::pair<std::string,std::vector<int>>> edge_cells;
    using ColorKey=std::tuple<std::string,std::string,std::string>;
    std::map<ColorKey,std::string> color_cache;
    std::unordered_map<std::string,std::vector<std::vector<int>>> group_cache;
    size_t budget, group_bytes=0, color_bytes=0;
    long calls=0,hits=0,solves=0,color_hits=0,clears=0;
public:
    ConditionedAutGraph(int n,const std::vector<std::pair<int,int>>& edges,
        std::vector<int> atom_labels,std::vector<int> vertex_indices,
        std::vector<std::string> colors,
        std::vector<std::pair<std::string,std::vector<int>>> fixed_cells,size_t cache_bytes)
        :graph(n,edges),atoms(std::move(atom_labels)),vertices(std::move(vertex_indices)),
         base_colors(std::move(colors)),edge_cells(std::move(fixed_cells)),budget(cache_bytes) {
        for(size_t i=0;i<atoms.size();++i) atom_index.emplace(atoms[i],i);
    }
    std::vector<std::vector<int>> generators(const std::vector<std::pair<int,int>>& locked,
                                             const std::map<int,std::string>& roles) {
        ++calls;
        std::vector<std::vector<int>> lock_roles(atoms.size());
        for(auto [r,p]:locked) {
            auto found=atom_index.find(p);
            if(found!=atom_index.end()) lock_roles[found->second].push_back(r);
        }
        std::map<std::string,std::vector<int>> colored;
        for(size_t i=0;i<atoms.size();++i) {
            auto& locks=lock_roles[i];std::sort(locks.begin(),locks.end());
            std::string locked_repr="(";
            for(size_t j=0;j<locks.size();++j) {
                if(j) locked_repr+=", ";
                locked_repr+="('locked', "+std::to_string(locks[j])+")";
            }
            if(locks.size()==1) locked_repr+=",";
            locked_repr+=")";
            auto role=roles.find(atoms[i]);
            ColorKey key{base_colors[i],locked_repr,role==roles.end()?"()":role->second};
            auto found=color_cache.find(key);
            if(found==color_cache.end()) {
                std::string value="('atom', "+std::get<0>(key)+", ("+std::get<1>(key)+", "+std::get<2>(key)+"))";
                size_t bytes=sizeof(ColorKey)+sizeof(std::string)+128+value.size()
                    +std::get<0>(key).size()+std::get<1>(key).size()+std::get<2>(key).size();
                if(color_bytes+bytes>budget/4) {color_cache.clear();color_bytes=0;++clears;}
                color_bytes+=bytes;
                found=color_cache.emplace(std::move(key),std::move(value)).first;
            } else ++color_hits;
            colored[found->second].push_back(vertices[i]);
        }
        for(const auto& [name,vs]:edge_cells) colored[name]=vs;
        std::vector<std::vector<int>> cells;
        std::string key;
        for(const auto& [name,vs]:colored) {
            // Reproduce Python set iteration, including non-contiguous vertex IDs.
            py::set cell;
            for(int v:vs) cell.add(py::int_(v));
            std::vector<int> ordered;
            for(auto v:cell) ordered.push_back(py::cast<int>(v));
            for(int v:ordered) key.append(reinterpret_cast<const char*>(&v),sizeof(v));
            int sep=-1;key.append(reinterpret_cast<const char*>(&sep),sizeof(sep));
            cells.push_back(std::move(ordered));
        }
        auto cached=group_cache.find(key);
        if(cached!=group_cache.end()) {++hits;return cached->second;}
        ++solves;
        if(cells.size()==1) cells.clear(); // pynauty's full-cell coloring rule
        auto raw=graph.generators(cells);
        std::map<int,int> vertex_to_atom;
        for(size_t i=0;i<atoms.size();++i) vertex_to_atom.emplace(vertices[i],atoms[i]);
        int degree=atoms.empty()?0:*std::max_element(atoms.begin(),atoms.end())+1;
        std::vector<int> identity(degree);
        for(int i=0;i<degree;++i) identity[i]=i;
        std::set<std::vector<int>> seen;
        std::vector<std::vector<int>> out;
        for(const auto& g:raw) {
            auto projected=identity;
            for(size_t i=0;i<atoms.size();++i) projected[atoms[i]]=vertex_to_atom.at(g[vertices[i]]);
            if(projected!=identity && seen.insert(projected).second) out.push_back(std::move(projected));
        }
        size_t bytes=key.size()+128+out.size()*(sizeof(std::vector<int>)+degree*sizeof(int));
        if(group_bytes+bytes>budget*3/4) {group_cache.clear();group_bytes=0;++clears;}
        if(bytes<=budget*3/4) {group_cache.emplace(std::move(key),out);group_bytes+=bytes;}
        return out;
    }
    py::dict stats() const {
        py::dict d;d["calls"]=calls;d["partition_hits"]=hits;d["group_solves"]=solves;
        d["color_hits"]=color_hits;d["cache_clears"]=clears;
        d["resident_bytes"]=color_bytes+group_bytes;d["cache_budget_bytes"]=budget;
        return d;
    }
};

}  // namespace autgrp_native

void register_autgrp(py::module_& mod) {
    using autgrp_native::AutGraph;
    py::class_<AutGraph>(mod, "AutGraph",
                         "Fixed undirected graph whose automorphism generators are computed\n"
                         "per vertex colouring exactly as pynauty.autgrp(graph)[0] would\n"
                         "(same nauty options, workspace and lab/ptn construction).")
        .def(py::init<int, const std::vector<std::pair<int, int>>&>(), py::arg("n_vertices"),
             py::arg("edges"))
        .def("generators", &AutGraph::generators, py::arg("cells"),
             "Generators (list of permutations) in nauty's emission order for the\n"
             "ordered partition ``cells`` (list of vertex lists; [] for no colouring).")
        .def_readonly("n_vertices", &AutGraph::n);
    using autgrp_native::ConditionedAutGraph;
    py::class_<ConditionedAutGraph>(mod,"ConditionedAutGraph")
        .def(py::init<int,const std::vector<std::pair<int,int>>&,std::vector<int>,std::vector<int>,
             std::vector<std::string>,std::vector<std::pair<std::string,std::vector<int>>>,size_t>())
        .def("generators",&ConditionedAutGraph::generators)
        .def("stats",&ConditionedAutGraph::stats);
}
