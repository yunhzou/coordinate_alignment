// Native AAM growth engine.
//
// A C++ port of rxn_core.growth.island.grow_island together with the matcher
// step it drives (rxn_core.matcher: state, extend, support, dedupe,
// canonical).  Every rule below mirrors one Python function; the Python
// engine remains the reference and the default, and bench/compare_grow_calls.py
// replays recorded Python calls through this module to prove identical
// outputs (mappings, deferred edges, fragments, symmetry states, cap
// behaviour) call by call.
//
// Only the element node policy and orbit maps carrying a structural zero
// bucket are supported; the Python wrapper falls back otherwise.
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <functional>
#include <map>
#include <memory>
#include <queue>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

extern "C" {
#include "nauty.h"
}

namespace py = pybind11;

namespace {

using Pair = std::pair<int, int>;

// ---------------------------------------------------------------------------
// nauty
// ---------------------------------------------------------------------------

struct DenseGraph {
    int n = 0;
    int m = 0;
    std::vector<setword> adjacency;
    std::vector<int> lab, ptn, orbits;
    std::vector<setword> canon, workspace;
    optionblk options;
    statsblk stats;

    DenseGraph(int vertices, const std::vector<Pair>& edges)
        : n(vertices), m(SETWORDSNEEDED(vertices)) {
        nauty_check(WORDSIZE, m, n, NAUTYVERSIONID);
        adjacency.assign((size_t)n * m, 0);
        for (const auto& e : edges) ADDONEEDGE(adjacency.data(), e.first, e.second, m);
        lab.assign(n, 0);
        ptn.assign(n, 0);
        orbits.assign(n, 0);
        canon.assign((size_t)n * m, 0);
        workspace.assign((size_t)66 * m, 0);
        static DEFAULTOPTIONS_GRAPH(defaults);
        options = defaults;
        options.getcanon = TRUE;
        options.defaultptn = FALSE;
        options.digraph = FALSE;
        options.writeautoms = FALSE;
    }

    // Exact memo: the graph is fixed, so the canonical form is a pure
    // function of the ordered cell partition.  Keyed by the partition itself
    // (no hashing shortcut); bounded so a pathological run cannot grow it
    // without limit.
    std::unordered_map<std::string, std::string> memo;
    size_t memo_limit = 500000;
    long memo_hits = 0, memo_misses = 0;

    // cells: vertex lists in cell order (must cover every vertex).
    std::string certificate(const std::vector<std::vector<int>>& cells) {
        std::string key;
        key.reserve((size_t)n * 4 + cells.size() * 4);
        int k = 0;
        for (const auto& cell : cells) {
            for (int v : cell) {
                lab[k] = v;
                ptn[k] = 1;
                ++k;
                key.append(reinterpret_cast<const char*>(&v), sizeof(int));
            }
            if (!cell.empty()) ptn[k - 1] = 0;
            int sep = -1;
            key.append(reinterpret_cast<const char*>(&sep), sizeof(int));
        }
        if (k != n) throw std::runtime_error("colouring does not cover every vertex");
        auto found = memo.find(key);
        if (found != memo.end()) {
            ++memo_hits;
            return found->second;
        }
        ++memo_misses;
        nauty(adjacency.data(), lab.data(), ptn.data(), nullptr, orbits.data(), &options,
              &stats, workspace.data(), (int)workspace.size(), m, n, canon.data());
        std::string out(reinterpret_cast<const char*>(canon.data()),
                        canon.size() * sizeof(setword));
        if (memo.size() < memo_limit) memo.emplace(std::move(key), out);
        return out;
    }
};

// ---------------------------------------------------------------------------
// graphs
// ---------------------------------------------------------------------------

struct Graph {
    int n = 0;
    std::vector<std::string> elem;
    std::vector<int> ecode;             // element code (shared table)
    std::vector<std::vector<int>> adj;  // ascending neighbour lists (active edges)
    std::vector<uint8_t> edge;          // n*n adjacency flags
    std::vector<double> w;              // n*n full WBO matrix
    double bond_cut = 0.2;

    bool has_edge(int a, int b) const { return edge[(size_t)a * n + b] != 0; }
    // primitives._edge_wbo: 0.0 on the diagonal, else the matrix value
    double wbo(int a, int b) const { return a == b ? 0.0 : w[(size_t)a * n + b]; }
};

// primitives._wbo_bucket: int(round(w * 5)) with Python's half-to-even round
inline long wbo_bucket(double w) { return (long)std::nearbyint(w * 5.0); }

// primitives._growth_edge_supported
inline bool growth_edge_supported(double w_R, double w_P, double iso_tol, double floor) {
    return w_P >= floor && std::fabs(w_R - w_P) <= iso_tol;
}

struct VecHash {
    size_t operator()(const std::vector<int64_t>& v) const noexcept {
        uint64_t h = 1469598103934665603ull;
        for (int64_t x : v) {
            h ^= (uint64_t)x;
            h *= 1099511628211ull;
            h ^= h >> 29;
        }
        return (size_t)h;
    }
};

// Injective integer ids for role / colour keys, assigned in first-seen
// order.  Certificates and role keys are only ever compared for equality
// inside one dedupe call, and an id is a pure function of its key for the
// lifetime of the table, so replacing the key vectors by ids preserves every
// equality decision while removing the per-atom vector allocations.
struct Interner {
    std::unordered_map<std::vector<int64_t>, int, VecHash> ids;
    int intern(const std::vector<int64_t>& key) {
        auto it = ids.find(key);
        if (it != ids.end()) return it->second;
        int id = (int)ids.size();
        ids.emplace(key, id);
        return id;
    }
};

struct Target : Graph {
    std::vector<int> orbit;        // p_orbits[v]
    std::vector<int> orbit_size;   // by orbit id
    std::vector<int> bucket;       // n*n pair bucket table (symmetric)
    int zero_bucket = 0;
    // nauty base graph: atoms 0..n-1, one vertex per pair with nonzero bucket
    int n_vertices = 0;
    std::vector<int> edge_vertex_bucket;                 // per edge vertex
    std::map<int, std::vector<int>> edge_cells;          // bucket -> edge vertices
    std::vector<int> edge_cell_colour;                   // colour id per edge cell (edge_cells order)
    std::unique_ptr<DenseGraph> nauty;
    mutable Interner roles;    // role keys -> ids
    mutable Interner colours;  // atom / edge colour keys -> ids
    std::vector<std::vector<int>> same_element;          // ecode -> atoms ascending

    int pair_bucket(int a, int b) const { return bucket[(size_t)a * n + b]; }
};

struct Source : Graph {};

// ---------------------------------------------------------------------------
// candidate state (rxn_core.matcher.state._SymCand / _SymBlock)
// ---------------------------------------------------------------------------

struct Block {
    std::vector<int> r;  // sorted unique
    std::vector<int> p;  // sorted unique
    bool extendable = true;
    bool open() const { return r.size() < p.size(); }
    bool complete() const { return r.size() >= p.size(); }
};

inline void sort_unique(std::vector<int>& v) {
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
}

inline Block make_block(std::vector<int> r, std::vector<int> p, bool extendable) {
    sort_unique(r);
    sort_unique(p);
    return Block{std::move(r), std::move(p), extendable};
}

inline bool contains(const std::vector<int>& sorted, int x) {
    return std::binary_search(sorted.begin(), sorted.end(), x);
}

struct Cand {
    std::vector<int> img;     // NR entries, -1 = unmapped
    std::vector<int> mapped;  // ascending mapped R atoms
    std::vector<Block> blocks;
    std::vector<int> exact_fixed;  // sorted unique
    long long mult = 1;
    std::vector<Block> automorph;

    bool has(int r) const { return img[r] >= 0; }
    bool has_open_choice() const {
        for (const auto& b : blocks)
            if (b.open()) return true;
        return false;
    }
};

struct BlockIndex {
    std::vector<int> r_to_block;  // NR, -1
    std::vector<int> p_to_block;  // NP, -1
};

BlockIndex block_indexes(const Cand& c, int NR, int NP) {
    BlockIndex bi;
    bi.r_to_block.assign(NR, -1);
    bi.p_to_block.assign(NP, -1);
    for (size_t i = 0; i < c.blocks.size(); ++i) {
        for (int r : c.blocks[i].r) bi.r_to_block[r] = (int)i;  // later blocks win
        for (int p : c.blocks[i].p) bi.p_to_block[p] = (int)i;
    }
    return bi;
}

// _SymCand.__init__: validate + back-fill the witness inside blocks; nullopt
// stands for the Python ValueError (callers return None).
bool make_cand(const std::vector<int>& raw, const std::vector<Block>& blocks,
               std::vector<int> exact_fixed, long long mult,
               const std::vector<Block>& automorph, int NP, Cand& out) {
    const int NR = (int)raw.size();
    std::vector<uint8_t> block_r(NR, 0);
    for (const auto& b : blocks)
        for (int r : b.r) block_r[r] = 1;
    std::vector<int> img(NR, -1);
    std::vector<uint8_t> used(NP, 0);
    for (int r = 0; r < NR; ++r) {
        if (raw[r] >= 0 && !block_r[r]) {
            img[r] = raw[r];
            used[raw[r]] = 1;
        }
    }
    for (const auto& b : blocks) {
        if (b.r.size() > b.p.size()) return false;
        for (int r : b.r) {
            if (raw[r] < 0) continue;
            int p = raw[r];
            if (!contains(b.p, p) || used[p]) return false;
            img[r] = p;
            used[p] = 1;
        }
        std::vector<int> available;
        for (int p : b.p)
            if (!used[p]) available.push_back(p);
        std::vector<int> missing;
        for (int r : b.r)
            if (img[r] < 0) missing.push_back(r);
        if (available.size() < missing.size()) return false;
        for (size_t i = 0; i < missing.size(); ++i) {
            img[missing[i]] = available[i];
            used[available[i]] = 1;
        }
    }
    out.img = std::move(img);
    out.mapped.clear();
    for (int r = 0; r < NR; ++r)
        if (out.img[r] >= 0) out.mapped.push_back(r);
    out.blocks = blocks;
    sort_unique(exact_fixed);
    out.exact_fixed = std::move(exact_fixed);
    out.mult = mult;
    out.automorph.clear();
    for (const auto& block : automorph) {
        std::vector<int> r_atoms = block.r;
        std::vector<uint8_t> in_pool(NP, 0);
        for (int p : block.p) in_pool[p] = 1;
        for (int r : out.mapped)
            if (in_pool[out.img[r]]) r_atoms.push_back(r);
        out.automorph.push_back(make_block(std::move(r_atoms), block.p, false));
    }
    return true;
}

bool rebuild(const Cand& c, const std::vector<int>& raw, const std::vector<Block>& blocks,
             int NP, Cand& out) {
    return make_cand(raw, blocks, c.exact_fixed, c.mult, c.automorph, NP, out);
}

bool with_fixed(const Cand& c, int r, int p, int NP, Cand& out) {
    if (c.has(r)) {
        if (c.img[r] != p) return false;
        out = c;
        return true;
    }
    for (int q : c.mapped)
        if (c.img[q] == p) return false;
    std::vector<int> raw = c.img;
    raw[r] = p;
    return rebuild(c, raw, c.blocks, NP, out);
}

bool with_new_block(const Cand& c, int r, std::vector<int> p_atoms, bool extendable, int NP,
                    Cand& out) {
    sort_unique(p_atoms);
    if (p_atoms.empty()) return false;
    if (p_atoms.size() == 1) return with_fixed(c, r, p_atoms[0], NP, out);
    if (c.has(r)) return false;
    std::vector<Block> blocks = c.blocks;
    blocks.push_back(make_block({r}, std::move(p_atoms), extendable));
    return rebuild(c, c.img, blocks, NP, out);
}

bool with_extended_block(const Cand& c, int idx, int r, int NP, Cand& out) {
    if (c.has(r)) return false;
    const Block& b = c.blocks[idx];
    if (!b.extendable || b.complete()) return false;
    std::vector<Block> blocks = c.blocks;
    std::vector<int> r_atoms = b.r;
    r_atoms.push_back(r);
    blocks[idx] = make_block(std::move(r_atoms), b.p, b.extendable);
    return rebuild(c, c.img, blocks, NP, out);
}

// assignments: (r, p) pairs, unique r
bool with_witness(const Cand& c, const std::vector<Pair>& assignments, int NR, int NP,
                  Cand& out) {
    std::vector<int> r_to_block(NR, -1);
    for (size_t i = 0; i < c.blocks.size(); ++i)
        for (int r : c.blocks[i].r) r_to_block[r] = (int)i;
    std::vector<uint8_t> touched(c.blocks.size() + 1, 0);
    for (const auto& a : assignments)
        if (r_to_block[a.first] >= 0) touched[r_to_block[a.first]] = 1;
    std::vector<int> raw(NR, -1);
    for (int r : c.mapped) {
        int bidx = r_to_block[r];
        if (bidx < 0 || !touched[bidx]) raw[r] = c.img[r];
    }
    for (const auto& a : assignments) raw[a.first] = a.second;
    return rebuild(c, raw, c.blocks, NP, out);
}

Cand with_multiplicity(const Cand& c, long long mult) {
    Cand out = c;
    out.mult = mult;
    return out;
}

// state._SymCand.with_automorph_equivalent
bool with_automorph_equivalent(const Cand& self, const Cand& other, int NP, Cand& out) {
    std::vector<int> varying_r;
    for (int r : self.mapped)
        if (other.has(r) && other.img[r] != self.img[r]) varying_r.push_back(r);
    std::vector<Block> blocks = self.automorph;
    blocks.insert(blocks.end(), other.automorph.begin(), other.automorph.end());
    if (!varying_r.empty()) {
        std::vector<int> p_atoms;
        for (int r : varying_r) {
            p_atoms.push_back(self.img[r]);
            p_atoms.push_back(other.img[r]);
        }
        blocks.push_back(make_block(varying_r, std::move(p_atoms), false));
    }
    std::vector<Block> merged;
    for (const auto& block : blocks) {
        std::vector<int> r_set = block.r, p_set = block.p;
        bool changed = true;
        while (changed) {
            changed = false;
            std::vector<Block> keep;
            for (const auto& prior : merged) {
                bool r_hit = false, p_hit = false;
                for (int r : r_set)
                    if (contains(prior.r, r)) { r_hit = true; break; }
                if (!r_hit)
                    for (int p : p_set)
                        if (contains(prior.p, p)) { p_hit = true; break; }
                if (r_hit || p_hit) {
                    r_set.insert(r_set.end(), prior.r.begin(), prior.r.end());
                    p_set.insert(p_set.end(), prior.p.begin(), prior.p.end());
                    sort_unique(r_set);
                    sort_unique(p_set);
                    changed = true;
                } else {
                    keep.push_back(prior);
                }
            }
            merged = std::move(keep);
        }
        merged.push_back(make_block(r_set, p_set, false));
    }
    std::sort(merged.begin(), merged.end(), [](const Block& a, const Block& b) {
        if (a.r != b.r) return a.r < b.r;
        return a.p < b.p;
    });
    return make_cand(self.img, self.blocks, self.exact_fixed, self.mult + other.mult, merged,
                     NP, out);
}

// support._refine_sym_assignments (the _SymCand branch)
bool refine_sym_assignments(const Cand& c, const std::vector<Pair>& assignments_in, int NR,
                            int NP, Cand& out) {
    // dict(assignments): later duplicates of r overwrite earlier ones
    std::vector<int> assign(NR, -1);
    std::vector<int> assigned_r;
    for (const auto& a : assignments_in) {
        if (assign[a.first] < 0) assigned_r.push_back(a.first);
        assign[a.first] = a.second;
    }
    std::sort(assigned_r.begin(), assigned_r.end());
    {
        std::vector<int> values;
        for (int r : assigned_r) values.push_back(assign[r]);
        std::sort(values.begin(), values.end());
        if (std::adjacent_find(values.begin(), values.end()) != values.end()) return false;
    }
    std::vector<uint8_t> assigned_value(NP, 0);
    for (int r : assigned_r) assigned_value[assign[r]] = 1;
    std::vector<int> r_to_block(NR, -1);
    for (size_t i = 0; i < c.blocks.size(); ++i)
        for (int r : c.blocks[i].r) r_to_block[r] = (int)i;
    std::vector<int> m(NR, -1);
    std::vector<uint8_t> used(NP, 0);
    for (int r : c.mapped)
        if (r_to_block[r] < 0) {
            m[r] = c.img[r];
            used[c.img[r]] = 1;
        }
    for (int r : assigned_r) {
        int p = assign[r];
        if (used[p] && m[r] != p) return false;
        m[r] = p;
        used[p] = 1;
    }
    std::vector<Block> new_blocks;
    for (const auto& block : c.blocks) {
        std::vector<int> remaining_r, remaining_p;
        for (int r : block.r)
            if (assign[r] < 0) remaining_r.push_back(r);
        for (int p : block.p)
            if (!assigned_value[p]) remaining_p.push_back(p);
        if (remaining_r.empty()) continue;
        if (remaining_r.size() > remaining_p.size()) return false;
        for (int r : remaining_r) {
            int p = c.img[r];
            if (p >= 0 && contains(remaining_p, p) && !used[p]) {
                m[r] = p;
                used[p] = 1;
            }
        }
        if (remaining_p.size() == 1) {
            int p = remaining_p[0];
            for (int r : remaining_r) {
                if (m[r] >= 0) continue;
                if (used[p]) return false;
                m[r] = p;
                used[p] = 1;
            }
            continue;
        }
        new_blocks.push_back(make_block(remaining_r, remaining_p, block.extendable));
    }
    return make_cand(m, new_blocks, c.exact_fixed, c.mult, c.automorph, NP, out);
}

// support._force_sym_value (the _SymCand branch)
bool force_sym_value(const Cand& c, int r, int p, int NR, int NP, Cand& out) {
    BlockIndex bi = block_indexes(c, NR, NP);
    if (c.has(r)) {
        int current = c.img[r];
        if (current == p) {
            out = c;
            return true;
        }
        int bidx = bi.r_to_block[r];
        if (bidx < 0 || !contains(c.blocks[bidx].p, p)) return false;
        std::vector<Pair> assignment{{r, p}};
        for (int other : c.blocks[bidx].r) {
            if (other != r && c.img[other] == p) {
                assignment.push_back({other, current});
                break;
            }
        }
        return refine_sym_assignments(c, assignment, NR, NP, out);
    }
    int bidx = bi.p_to_block[p];
    if (bidx >= 0) {
        Cand nc;
        if (!with_extended_block(c, bidx, r, NP, nc)) return false;
        return refine_sym_assignments(nc, {{r, p}}, NR, NP, out);
    }
    for (int q : c.mapped)
        if (bi.r_to_block[q] < 0 && c.img[q] == p) return false;
    return with_fixed(c, r, p, NP, out);
}

// ---------------------------------------------------------------------------
// extension context (matcher.extend._ExtensionContext)
// ---------------------------------------------------------------------------

struct Context {
    const Source* R;
    const Target* P;
    std::vector<uint8_t> fragment_old;   // NR flags
    int n;
    const std::vector<int>* mapping;     // locked mapping image, NR, -1
    std::vector<uint8_t> locked_p;       // NP flags
    double iso_tol;
    const std::vector<Pair>* islands;    // islands_R items in dict order, or null
    const std::vector<int>* island_of;   // NR, -1
    std::vector<Pair> deferred_edges;
    int anchor_u;
    double anchor_wbo;
    bool has_anchor_wbo;
    std::vector<Pair> dedupe_edges;      // boundary evidence for dedupe
    std::vector<int> bonded_in_frag;     // ascending
    std::vector<double> r_wbos;          // parallel to bonded_in_frag
    int strict_r;                        // anchor_u when it lies in fragment_old, else -1
    double strict_w;
    std::vector<int> island_atoms;
    std::vector<uint8_t> sig_fragment;   // fragment_old | {n} | island_atoms

    bool is_merge() const { return (*mapping)[n] >= 0; }
};

// support._support_witness_for_value; returns false for None, else fills
// `support` with (u, p) block assignments.
bool support_witness_for_value(const Cand& c, const Context& ctx, int v_n, int join_block_idx,
                               const BlockIndex& bi, std::vector<Pair>& support,
                               int max_states = 4096) {
    const Target& P = *ctx.P;
    const double floor = P.bond_cut;
    auto pair_ok = [&](int r_atom, double w_R, int p_atom) {
        double w_P = P.wbo(p_atom, v_n);
        if (r_atom == ctx.strict_r) return growth_edge_supported(ctx.strict_w, w_P, ctx.iso_tol, floor);
        return growth_edge_supported(w_R, w_P, ctx.iso_tol, floor);
    };
    // fixed_used: images of atoms outside every block
    for (int r : c.mapped)
        if (bi.r_to_block[r] < 0 && c.img[r] == v_n && c.img[ctx.n] != v_n) return false;
    int value_block = bi.p_to_block[v_n];
    if (value_block >= 0 && value_block != join_block_idx) return false;
    if (join_block_idx >= 0) {
        const Block& block = c.blocks[join_block_idx];
        if (!contains(block.p, v_n) || !block.open()) return false;
    }
    support.clear();
    std::vector<int> block_order;
    std::vector<std::vector<std::pair<int, double>>> by_block(c.blocks.size());
    for (size_t i = 0; i < ctx.bonded_in_frag.size(); ++i) {
        int u = ctx.bonded_in_frag[i];
        double w = ctx.r_wbos[i];
        if (!c.has(u)) return false;
        int bidx = bi.r_to_block[u];
        if (bidx < 0) {
            int p_u = c.img[u];
            if (p_u == v_n) return false;
            if (!pair_ok(u, w, p_u)) return false;
            continue;
        }
        if (by_block[bidx].empty()) block_order.push_back(bidx);
        by_block[bidx].push_back({u, w});
    }
    for (int bidx : block_order) {
        const Block& block = c.blocks[bidx];
        std::vector<int> available_pool;
        for (int p : block.p)
            if (!(bidx == join_block_idx && p == v_n)) available_pool.push_back(p);
        std::vector<std::pair<int, std::vector<int>>> domains;
        for (const auto& item : by_block[bidx]) {
            std::vector<int> vals;
            for (int p : available_pool)
                if (pair_ok(item.first, item.second, p)) vals.push_back(p);
            if (vals.empty()) return false;
            domains.push_back({item.first, std::move(vals)});
        }
        bool all_full = true;
        for (const auto& d : domains)
            if (d.second.size() != available_pool.size()) { all_full = false; break; }
        if (all_full) continue;
        std::stable_sort(domains.begin(), domains.end(),
                         [](const auto& a, const auto& b) {
                             if (a.second.size() != b.second.size())
                                 return a.second.size() < b.second.size();
                             return a.first < b.first;
                         });
        std::vector<Pair> chosen;
        std::vector<uint8_t> used(P.n, 0);
        if (bidx == join_block_idx) used[v_n] = 1;
        int states = 0;
        std::function<bool(size_t)> backtrack = [&](size_t pos) -> bool {
            ++states;
            if (states > max_states) return false;
            if (pos == domains.size()) return true;
            int u = domains[pos].first;
            for (int p : domains[pos].second) {
                if (used[p]) continue;
                used[p] = 1;
                chosen.push_back({u, p});
                if (backtrack(pos + 1)) return true;
                chosen.pop_back();
                used[p] = 0;
            }
            return false;
        };
        if (!backtrack(0)) return false;
        support.insert(support.end(), chosen.begin(), chosen.end());
    }
    return true;
}

// support._r_compatible_with_block
bool r_compatible_with_block(const Cand& c, int bidx, const Context& ctx) {
    const Block& b = c.blocks[bidx];
    if (!b.extendable || b.complete()) return false;
    const Source& R = *ctx.R;
    std::vector<int> outside;
    for (int r = 0; r < R.n; ++r)
        if (ctx.fragment_old[r] && !contains(b.r, r) && c.has(r)) outside.push_back(r);
    auto rel_sig = [&](int r) {
        std::vector<std::pair<int, long>> out;
        for (int x : outside)
            if (R.has_edge(r, x)) out.push_back({x, wbo_bucket(R.wbo(r, x))});
        return out;
    };
    auto n_sig = rel_sig(ctx.n);
    for (int r : b.r)
        if (rel_sig(r) != n_sig) return false;
    return true;
}

// ---------------------------------------------------------------------------
// relation signatures (matcher.dedupe._p_relation_signature_from_parts)
// ---------------------------------------------------------------------------

// Python repr of the dense signature tuple, used to order target groups.
std::string python_repr_int(long v) { return std::to_string(v); }

struct DenseSignature {
    std::string elem;
    int orbit;
    std::vector<std::pair<int, int>> rel;                 // (r, bucket)
    std::vector<std::tuple<int, bool, std::vector<int>>> block_rel;

    std::string repr() const {
        std::string s = "('" + elem + "', " + std::to_string(orbit) + ", ";
        // rel tuple
        if (rel.empty()) s += "()";
        else {
            s += "(";
            for (size_t i = 0; i < rel.size(); ++i) {
                if (i) s += ", ";
                s += "(" + std::to_string(rel[i].first) + ", " + std::to_string(rel[i].second) + ")";
            }
            if (rel.size() == 1) s += ",";
            s += ")";
        }
        s += ", ";
        if (block_rel.empty()) s += "()";
        else {
            s += "(";
            for (size_t i = 0; i < block_rel.size(); ++i) {
                if (i) s += ", ";
                const auto& br = block_rel[i];
                s += "(" + std::to_string(std::get<0>(br)) + ", " +
                     (std::get<1>(br) ? "True" : "False") + ", ";
                const auto& ews = std::get<2>(br);
                if (ews.empty()) s += "()";
                else {
                    s += "(";
                    for (size_t j = 0; j < ews.size(); ++j) {
                        if (j) s += ", ";
                        s += std::to_string(ews[j]);
                    }
                    if (ews.size() == 1) s += ",";
                    s += ")";
                }
                s += ")";
            }
            if (block_rel.size() == 1) s += ",";
            s += ")";
        }
        s += ")";
        return s;
    }
    bool operator==(const DenseSignature& o) const {
        return elem == o.elem && orbit == o.orbit && rel == o.rel && block_rel == o.block_rel;
    }
};

DenseSignature dense_signature(const Cand& c, int v, const Target& P) {
    DenseSignature sig;
    sig.elem = P.elem[v];
    sig.orbit = P.orbit[v];
    for (int r : c.mapped) {
        int p = c.img[r];
        if (p == v) continue;
        sig.rel.push_back({r, P.pair_bucket(p, v)});
    }
    for (size_t i = 0; i < c.blocks.size(); ++i) {
        const Block& b = c.blocks[i];
        std::vector<int> ews;
        for (int p : b.p) {
            if (p == v) continue;
            ews.push_back(P.pair_bucket(p, v));
        }
        std::sort(ews.begin(), ews.end());
        sig.block_rel.emplace_back((int)i, contains(b.p, v), std::move(ews));
    }
    return sig;
}

// Compact signature encoded as an integer blob (equality only).  The compact
// and dense signatures induce the same partition (structural-zero pairs are
// exactly the non-edges of the bucket table).
void compact_signature_blob(const Cand& c, int v, const Target& P, std::vector<int64_t>& blob) {
    blob.clear();
    blob.push_back(P.ecode[v]);
    blob.push_back(P.orbit[v]);
    // mapped_r sequence excluding r with image v
    blob.push_back(-1);
    for (int r : c.mapped)
        if (c.img[r] != v) blob.push_back(r);
    // nonzero relations (r, bucket) sorted
    std::vector<std::pair<int, int>> rel;
    for (int p : P.adj[v]) {
        int bucket = P.pair_bucket(p, v);
        if (bucket == P.zero_bucket) continue;
        // r with image p
        for (int r : c.mapped) {
            if (c.img[r] == p && p != v) { rel.push_back({r, bucket}); break; }
        }
    }
    std::sort(rel.begin(), rel.end());
    blob.push_back(-2);
    for (const auto& x : rel) { blob.push_back(x.first); blob.push_back(x.second); }
    blob.push_back(-3);
    for (size_t i = 0; i < c.blocks.size(); ++i) {
        const Block& b = c.blocks[i];
        bool member = contains(b.p, v);
        std::vector<int> ews;
        for (int p : b.p) {
            if (p == v) continue;
            if (!P.has_edge(p, v)) continue;
            int bucket = P.pair_bucket(p, v);
            if (bucket != P.zero_bucket) ews.push_back(bucket);
        }
        std::sort(ews.begin(), ews.end());
        blob.push_back(-4);
        blob.push_back((int64_t)i);
        blob.push_back(member ? 1 : 0);
        blob.push_back((int64_t)b.p.size() - (member ? 1 : 0));
        for (int e : ews) blob.push_back(e);
    }
}

// ---------------------------------------------------------------------------
// certificates (matcher.canonical._CandidateAutomorphismCanonicalizer)
// ---------------------------------------------------------------------------

using Color = std::vector<int64_t>;

struct Certificate {
    std::string canon;
    std::string profile;  // ordered colour keys with cell sizes
    bool operator==(const Certificate& o) const { return canon == o.canon && profile == o.profile; }
};
struct CertificateHash {
    size_t operator()(const Certificate& c) const {
        return std::hash<std::string>()(c.canon) ^ (std::hash<std::string>()(c.profile) << 1);
    }
};

void append_blob(std::string& out, const int64_t* data, size_t count) {
    out.append(reinterpret_cast<const char*>(data), count * sizeof(int64_t));
}

// roles per target atom (matcher.state._cand_roles_from_scratch) as a CSR
// table of interned role ids, ascending within each atom.  Python sorts an
// atom's roles by repr so its role tuple is a function of the role multiset;
// the ascending id list is the same function of the same multiset.
struct CandRoles {
    std::vector<int> start;  // NP + 1
    std::vector<int> ids;
    bool any(int p) const { return start[p + 1] > start[p]; }
};

struct Canonicalizer {
    const Target* P;
    std::vector<int> locked_r_of_p;  // NP, -1: locked R atom whose image is p
    mutable std::vector<int> mapped_role_id;  // per r, -1 until interned

    Canonicalizer(const Target* target, const std::vector<int>& mapping) : P(target) {
        locked_r_of_p.assign(target->n, -1);
        for (int r = 0; r < (int)mapping.size(); ++r)
            if (mapping[r] >= 0) locked_r_of_p[mapping[r]] = r;
        mapped_role_id.assign(mapping.size(), -1);
    }

    int mapped_role(int r) const {
        if (r >= (int)mapped_role_id.size()) mapped_role_id.resize(r + 1, -1);
        int& id = mapped_role_id[r];
        if (id < 0) id = P->roles.intern({0, (int64_t)r});
        return id;
    }

    void candidate_roles(const Cand& c, CandRoles& out) const {
        static thread_local std::vector<std::pair<int, int>> pairs;  // (p, role id)
        static thread_local std::vector<uint8_t> block_r;
        static thread_local std::vector<int64_t> key;
        pairs.clear();
        block_r.assign(c.img.size(), 0);
        for (const auto& b : c.blocks)
            for (int r : b.r) block_r[r] = 1;
        for (int r : c.mapped)
            if (!block_r[r]) pairs.push_back({c.img[r], mapped_role(r)});
        for (const auto& b : c.blocks) {
            key.clear();
            key.push_back(1);
            key.push_back(b.extendable ? 1 : 0);
            for (int r : b.r) key.push_back(r);
            int id = P->roles.intern(key);
            for (int p : b.p) pairs.push_back({p, id});
        }
        for (const auto& b : c.automorph) {
            key.clear();
            key.push_back(2);
            for (int r : b.r) key.push_back(r);
            int id = P->roles.intern(key);
            for (int p : b.p) pairs.push_back({p, id});
        }
        std::sort(pairs.begin(), pairs.end());
        const int NP = P->n;
        out.start.assign(NP + 1, 0);
        for (const auto& pr : pairs) out.start[pr.first + 1]++;
        for (int p = 0; p < NP; ++p) out.start[p + 1] += out.start[p];
        out.ids.resize(pairs.size());
        for (size_t i = 0; i < pairs.size(); ++i) out.ids[i] = pairs[i].second;
    }

    // Certificate of the candidate colouring: atom colour (0, ecode,
    // locked_r, role ids...) and edge colour (1, bucket) keys are interned,
    // cells are ordered by colour id (a function of the key alone), and the
    // profile records (colour id, cell size) per cell.  Two candidates get
    // equal certificates iff an automorphism of the target maps one colouring
    // onto the other, exactly as with the Python colour-key ordering.
    Certificate certificate(const CandRoles& roles) const {
        static thread_local std::vector<int64_t> key;
        static thread_local std::vector<std::pair<int, int>> coloured;  // (colour id, atom)
        static thread_local std::vector<std::vector<int>> ordered;
        static thread_local std::vector<std::pair<int, int>> cell_order;  // (colour id, cell index)
        const int NP = P->n;
        coloured.clear();
        for (int p = 0; p < NP; ++p) {
            key.clear();
            key.push_back(0);
            key.push_back(P->ecode[p]);
            key.push_back(locked_r_of_p[p]);
            for (int k = roles.start[p]; k < roles.start[p + 1]; ++k) key.push_back(roles.ids[k]);
            coloured.push_back({P->colours.intern(key), p});
        }
        std::sort(coloured.begin(), coloured.end());
        // cells: atom cells (grouped by colour id) plus the edge cells
        size_t n_cells = 0;
        cell_order.clear();
        auto cell_at = [&](size_t index) -> std::vector<int>& {
            if (ordered.size() <= index) ordered.emplace_back();
            ordered[index].clear();
            return ordered[index];
        };
        for (size_t i = 0; i < coloured.size();) {
            size_t j = i;
            std::vector<int>& cell = cell_at(n_cells);
            while (j < coloured.size() && coloured[j].first == coloured[i].first) {
                cell.push_back(coloured[j].second);
                ++j;
            }
            cell_order.push_back({coloured[i].first, (int)n_cells});
            ++n_cells;
            i = j;
        }
        size_t edge_index = 0;
        for (const auto& item : P->edge_cells) {
            std::vector<int>& cell = cell_at(n_cells);
            cell.assign(item.second.begin(), item.second.end());
            cell_order.push_back({P->edge_cell_colour[edge_index++], (int)n_cells});
            ++n_cells;
        }
        std::sort(cell_order.begin(), cell_order.end());
        static thread_local std::vector<std::vector<int>> cells;
        cells.resize(n_cells);
        Certificate cert;
        for (size_t k = 0; k < n_cells; ++k) {
            cells[k].swap(ordered[cell_order[k].second]);
            int64_t colour = cell_order[k].first;
            int64_t size = (int64_t)cells[k].size();
            append_blob(cert.profile, &colour, 1);
            append_blob(cert.profile, &size, 1);
        }
        cert.canon = P->nauty->certificate(cells);
        // give the buffers back so `ordered` keeps its capacity
        for (size_t k = 0; k < n_cells; ++k) ordered[cell_order[k].second].swap(cells[k]);
        return cert;
    }

    // orbit-role key: multiset of (orbit, locked_r, role ids...) over role
    // atoms; `singleton` is true when every role atom sits in a singleton orbit.
    std::string role_key(const CandRoles& roles, bool& singleton) const {
        static thread_local std::vector<std::vector<int64_t>> items;
        singleton = true;
        size_t count = 0;
        for (int p = 0; p < P->n; ++p) {
            if (!roles.any(p)) continue;
            if (P->orbit_size[P->orbit[p]] > 1) singleton = false;
            if (items.size() <= count) items.emplace_back();
            auto& item = items[count++];
            item.clear();
            item.push_back((int64_t)P->orbit[p]);
            item.push_back((int64_t)locked_r_of_p[p]);
            for (int k = roles.start[p]; k < roles.start[p + 1]; ++k) item.push_back(roles.ids[k]);
        }
        std::sort(items.begin(), items.begin() + count);
        std::string blob;
        for (size_t i = 0; i < count; ++i) {
            append_blob(blob, items[i].data(), items[i].size());
            int64_t sep = INT64_MIN;
            append_blob(blob, &sep, 1);
        }
        return blob;
    }
};

// ---------------------------------------------------------------------------
// boundary signature (matcher.dedupe._boundary_signature, candidate part)
// ---------------------------------------------------------------------------

struct BoundaryContext {
    // pools in first-seen order of boundary atoms (sorted): one per element
    std::vector<int> pool_ecodes;                // distinct element codes in order
    std::vector<std::vector<int>> pools;         // compatible unlocked targets
    bool active = false;
};

BoundaryContext boundary_context(const Context& ctx, const std::vector<uint8_t>& fragment,
                                 const std::vector<Pair>& deferred_edges) {
    BoundaryContext bc;
    if (deferred_edges.empty()) return bc;
    bool any_fragment = false;
    for (uint8_t f : fragment)
        if (f) { any_fragment = true; break; }
    if (!any_fragment) return bc;
    std::vector<int> boundary;
    for (const auto& e : deferred_edges) {
        bool a_in = fragment[e.first], b_in = fragment[e.second];
        if (a_in == b_in) continue;
        boundary.push_back(a_in ? e.second : e.first);
    }
    sort_unique(boundary);
    if (boundary.empty()) return bc;
    bc.active = true;
    for (int x : boundary) {
        int code = ctx.R->ecode[x];
        if (std::find(bc.pool_ecodes.begin(), bc.pool_ecodes.end(), code) != bc.pool_ecodes.end())
            continue;
        bc.pool_ecodes.push_back(code);
        std::vector<int> pool;
        for (int v = 0; v < ctx.P->n; ++v)
            if (!ctx.locked_p[v] && ctx.P->ecode[v] == code) pool.push_back(v);
        bc.pools.push_back(std::move(pool));
    }
    return bc;
}

std::string boundary_signature(const Cand& c, const Context& ctx, const BoundaryContext& bc) {
    std::string out;
    if (!bc.active) return out;
    std::vector<uint8_t> used_possible(ctx.P->n, 0);
    for (int r : c.mapped) used_possible[c.img[r]] = 1;
    for (const auto& b : c.blocks)
        for (int p : b.p) used_possible[p] = 1;
    std::vector<int64_t> blob;
    for (const auto& pool : bc.pools) {
        std::vector<std::vector<int64_t>> sigs;
        for (int v : pool) {
            if (used_possible[v]) continue;
            compact_signature_blob(c, v, *ctx.P, blob);
            sigs.push_back(blob);
        }
        std::sort(sigs.begin(), sigs.end());
        for (const auto& s : sigs) {
            append_blob(out, s.data(), s.size());
            int64_t sep = INT64_MIN;
            append_blob(out, &sep, 1);
        }
        int64_t sep2 = INT64_MIN + 1;
        append_blob(out, &sep2, 1);
    }
    return out;
}

// ---------------------------------------------------------------------------
// dedupe (matcher.dedupe._dedup_sym_cands)
// ---------------------------------------------------------------------------

struct Engine;

std::vector<Cand> dedup_sym_cands(std::vector<Cand>& cands, const Context& ctx,
                                  const Canonicalizer& canon, const std::vector<uint8_t>& fragment,
                                  const std::vector<Pair>& deferred_edges, long& certificate_calls) {
    if (cands.size() <= 1) return cands;
    const int NP = ctx.P->n;
    // roles and orbit-role keys
    std::vector<CandRoles> roles(cands.size());
    std::vector<std::string> keys(cands.size());
    std::vector<bool> singleton(cands.size());
    std::unordered_map<std::string, std::vector<int>> classes;
    std::vector<std::string> class_order;
    for (size_t i = 0; i < cands.size(); ++i) {
        canon.candidate_roles(cands[i], roles[i]);
        bool single;
        keys[i] = canon.role_key(roles[i], single);
        singleton[i] = single;
        auto it = classes.find(keys[i]);
        if (it == classes.end()) {
            classes.emplace(keys[i], std::vector<int>{(int)i});
            class_order.push_back(keys[i]);
        } else {
            it->second.push_back((int)i);
        }
    }
    // certificates are interned to ids for this call: stand-ins for
    // singleton classes and singleton-orbit classes, nauty otherwise
    // (memoised per identical roles)
    std::vector<Certificate> cert_pool;
    std::unordered_map<Certificate, int, CertificateHash> cert_ids;
    auto cert_id = [&](Certificate&& cert) {
        auto it = cert_ids.find(cert);
        if (it != cert_ids.end()) return it->second;
        int id = (int)cert_pool.size();
        cert_pool.push_back(cert);
        cert_ids.emplace(std::move(cert), id);
        return id;
    };
    std::vector<int> certs(cands.size(), -1);
    std::unordered_map<std::string, int> by_roles;
    for (const auto& key : class_order) {
        const auto& members = classes[key];
        if (members.size() == 1 || singleton[members[0]]) {
            Certificate stand_in;
            stand_in.profile = key;  // unique per class; canon empty
            int id = cert_id(std::move(stand_in));
            for (int i : members) certs[i] = id;
        } else {
            for (int i : members) {
                const CandRoles& cr = roles[i];
                std::string roles_blob(reinterpret_cast<const char*>(cr.start.data()),
                                       cr.start.size() * sizeof(int));
                roles_blob.append(reinterpret_cast<const char*>(cr.ids.data()),
                                  cr.ids.size() * sizeof(int));
                auto it = by_roles.find(roles_blob);
                if (it == by_roles.end()) {
                    ++certificate_calls;
                    int id = cert_id(canon.certificate(cr));
                    it = by_roles.emplace(std::move(roles_blob), id).first;
                }
                certs[i] = it->second;
            }
        }
    }
    // counts per certificate
    std::vector<int> counts(cert_pool.size(), 0);
    for (int id : certs) counts[id] += 1;
    BoundaryContext bc;
    bool bc_built = false;
    struct Key {
        int cert;
        std::string boundary;
        bool operator==(const Key& o) const { return cert == o.cert && boundary == o.boundary; }
    };
    struct KeyHash {
        size_t operator()(const Key& k) const {
            return std::hash<int>()(k.cert) ^ (std::hash<std::string>()(k.boundary) << 2);
        }
    };
    std::unordered_map<Key, int, KeyHash> seen;  // key -> index into out
    std::vector<Cand> out;
    std::unordered_map<std::string, std::string> boundary_memo;
    for (size_t i = 0; i < cands.size(); ++i) {
        std::string boundary;
        if (counts[certs[i]] > 1 && !deferred_edges.empty()) {
            if (!bc_built) {
                bc = boundary_context(ctx, fragment, deferred_edges);
                bc_built = true;
            }
            boundary = boundary_signature(cands[i], ctx, bc);
        }
        Key key{certs[i], boundary};
        auto it = seen.find(key);
        if (it == seen.end()) {
            seen.emplace(std::move(key), (int)out.size());
            out.push_back(std::move(cands[i]));
        } else {
            Cand merged;
            if (!with_automorph_equivalent(out[it->second], cands[i], NP, merged))
                throw std::runtime_error("automorph merge rejected by candidate validation");
            out[it->second] = std::move(merged);
        }
    }
    return out;
}

// ---------------------------------------------------------------------------
// extension (matcher.extend._extend_sym_cands)
// ---------------------------------------------------------------------------

struct TargetEntry {
    int v;
    std::vector<Pair> support;
    bool can_extend;
};

bool island_merge_wbo_consistent(const Cand& child, const Context& ctx) {
    const Source& R = *ctx.R;
    const Target& P = *ctx.P;
    const double floor = P.bond_cut;
    std::vector<uint8_t> is_island(R.n, 0);
    for (int r : ctx.island_atoms) is_island[r] = 1;
    for (int r : ctx.island_atoms) {
        for (int r2 : child.mapped) {  // sorted(check_set)
            if (r2 == r) continue;
            if (r >= r2 && is_island[r2]) continue;
            if (!R.has_edge(r, r2)) continue;
            double w_r = R.wbo(r, r2);
            int p = child.img[r], p2 = child.img[r2];
            if (!growth_edge_supported(w_r, P.wbo(p, p2), ctx.iso_tol, floor)) return false;
        }
    }
    return true;
}

bool extend_locked_merge(const Cand& c, const Context& ctx, Cand& out) {
    const int NR = ctx.R->n, NP = ctx.P->n;
    int v_n = (*ctx.mapping)[ctx.n];
    BlockIndex bi = block_indexes(c, NR, NP);
    int join_idx = bi.p_to_block[v_n];
    if (join_idx >= 0 && !r_compatible_with_block(c, join_idx, ctx)) return false;
    std::vector<Pair> support;
    if (!support_witness_for_value(c, ctx, v_n, join_idx, bi, support)) return false;
    Cand child;
    if (!force_sym_value(c, ctx.n, v_n, NR, NP, child)) return false;
    support.push_back({ctx.n, v_n});
    Cand witnessed;
    if (!with_witness(child, support, NR, NP, witnessed)) return false;
    child = std::move(witnessed);
    for (int r : ctx.island_atoms) {
        if (r == ctx.n) continue;
        Cand forced;
        if (!force_sym_value(child, r, (*ctx.mapping)[r], NR, NP, forced)) return false;
        child = std::move(forced);
    }
    if (!island_merge_wbo_consistent(child, ctx)) return false;
    out = std::move(child);
    return true;
}

void extend_free_atom(const Cand& c, const Context& ctx, std::vector<Cand>& children) {
    const Source& R = *ctx.R;
    const Target& P = *ctx.P;
    const int NR = R.n, NP = P.n;
    BlockIndex bi = block_indexes(c, NR, NP);
    // admissible targets: neighbours of every fixed image / of some pool atom
    std::vector<uint8_t> admissible;
    bool first = true;
    for (int u : ctx.bonded_in_frag) {
        std::vector<uint8_t> reach(NP, 0);
        int bidx = bi.r_to_block[u];
        if (bidx < 0) {
            for (int q : P.adj[c.img[u]]) reach[q] = 1;
        } else {
            for (int p : c.blocks[bidx].p)
                for (int q : P.adj[p]) reach[q] = 1;
        }
        if (first) { admissible = reach; first = false; }
        else for (int q = 0; q < NP; ++q) admissible[q] &= reach[q];
    }
    if (first) admissible.assign(NP, 1);
    std::map<int, std::vector<TargetEntry>> block_join;
    struct Group { DenseSignature sig; std::vector<TargetEntry> entries; std::string repr; };
    std::vector<Group> groups;
    std::vector<int> compat_cache(c.blocks.size(), -1);
    const int n_code = R.ecode[ctx.n];
    if (n_code < 0 || n_code >= (int)P.same_element.size()) return;
    for (int v : P.same_element[n_code]) {
        if (!admissible[v] || ctx.locked_p[v]) continue;
        int join_idx = bi.p_to_block[v];
        bool can_extend = false;
        if (join_idx >= 0) {
            if (!c.blocks[join_idx].open()) join_idx = -1;
            else {
                if (compat_cache[join_idx] < 0)
                    compat_cache[join_idx] = r_compatible_with_block(c, join_idx, ctx) ? 1 : 0;
                can_extend = compat_cache[join_idx] == 1;
            }
        }
        std::vector<Pair> support;
        if (!support_witness_for_value(c, ctx, v, join_idx, bi, support)) continue;
        if (join_idx >= 0) {
            block_join[join_idx].push_back({v, std::move(support), can_extend});
        } else {
            DenseSignature sig = dense_signature(c, v, P);
            bool placed = false;
            for (auto& g : groups)
                if (g.sig == sig) { g.entries.push_back({v, std::move(support), true}); placed = true; break; }
            if (!placed) {
                Group g;
                g.sig = sig;
                g.repr = g.sig.repr();
                g.entries.push_back({v, std::move(support), true});
                groups.push_back(std::move(g));
            }
        }
    }
    // block joins in ascending block index
    for (auto& item : block_join) {
        int join_idx = item.first;
        auto& entries = item.second;
        std::vector<int> free_vs;
        for (const auto& e : entries)
            if (e.support.empty() && e.can_extend) free_vs.push_back(e.v);
        if (!free_vs.empty()) {
            int witness_v = *std::min_element(free_vs.begin(), free_vs.end());
            Cand child;
            if (with_extended_block(c, join_idx, ctx.n, NP, child)) {
                Cand witnessed;
                if (with_witness(child, {{ctx.n, witness_v}}, NR, NP, witnessed))
                    children.push_back(std::move(witnessed));
            }
        }
        for (const auto& e : entries) {
            if (e.support.empty() && e.can_extend) continue;
            std::vector<Pair> fixed = e.support;
            fixed.push_back({ctx.n, e.v});
            Cand child;
            if (refine_sym_assignments(c, fixed, NR, NP, child)) children.push_back(std::move(child));
        }
    }
    // context groups ordered by Python str() of the dense signature
    std::vector<size_t> order(groups.size());
    for (size_t i = 0; i < order.size(); ++i) order[i] = i;
    std::stable_sort(order.begin(), order.end(),
                     [&](size_t a, size_t b) { return groups[a].repr < groups[b].repr; });
    for (size_t gi : order) {
        const auto& entries = groups[gi].entries;
        for (const auto& e : entries) {
            if (e.support.empty()) continue;
            std::vector<Pair> fixed = e.support;
            fixed.push_back({ctx.n, e.v});
            Cand child;
            if (refine_sym_assignments(c, fixed, NR, NP, child)) children.push_back(std::move(child));
        }
        std::vector<int> group;
        for (const auto& e : entries)
            if (e.support.empty()) group.push_back(e.v);
        if (group.empty()) continue;
        std::sort(group.begin(), group.end());
        int witness_v = group[0];
        Cand child;
        bool ok;
        if (group.size() > 1) ok = with_new_block(c, ctx.n, group, true, NP, child);
        else ok = with_fixed(c, ctx.n, witness_v, NP, child);
        if (!ok) continue;
        Cand witnessed;
        if (with_witness(child, {{ctx.n, witness_v}}, NR, NP, witnessed))
            children.push_back(std::move(witnessed));
    }
}

std::vector<Cand> extend_sym_cands(const std::vector<Cand>& cands, const Context& ctx,
                                   const Canonicalizer& canon, long& certificate_calls) {
    std::vector<Cand> children;
    for (const auto& c : cands) {
        bool covers = true;
        for (int u : ctx.bonded_in_frag)
            if (!c.has(u)) { covers = false; break; }
        if (!covers) continue;
        if (ctx.is_merge()) {
            Cand child;
            if (extend_locked_merge(c, ctx, child)) children.push_back(std::move(child));
        } else {
            extend_free_atom(c, ctx, children);
        }
    }
    return dedup_sym_cands(children, ctx, canon, ctx.sig_fragment, ctx.dedupe_edges,
                           certificate_calls);
}

// ---------------------------------------------------------------------------
// growth loop (rxn_core.growth.island.grow_island)
// ---------------------------------------------------------------------------

struct HeapItem {
    double neg_w;
    int u, n;
    bool operator>(const HeapItem& o) const {
        if (neg_w != o.neg_w) return neg_w > o.neg_w;
        if (u != o.u) return u > o.u;
        return n > o.n;
    }
};

struct GrowProfile {
    long seed_targets = 0, seed_groups = 0, init_cands = 0, initial_heap = 0, heap_pops = 0,
         stale_pops = 0, fragment_skip_pops = 0, extend_calls = 0, max_cands_before = 0,
         max_cands_after = 0, max_fragment_size = 0, max_heap_len = 0, commits = 0,
         deferred = 0, merge_calls = 0, free_extend_calls = 0, certificate_calls = 0,
         final_cands = 0, final_fragment_size = 0, branches = 0;
    std::string result;
};

struct BranchCap {
    long count, limit;
};

struct IsoOut {
    std::vector<Pair> mapping;             // sorted (r, p)
    std::vector<Pair> deferred_edges;      // sorted
    std::vector<int> fragment;             // sorted
    Cand cand;
};

struct GrowResult {
    std::vector<IsoOut> isos;
    GrowProfile profile;
    bool capped = false;
    BranchCap cap{0, 0};
};

GrowResult grow_island(const Source& R, const Target& P, int seed, const std::vector<int>& mapping,
                       double graph_floor, double iso_tol, int min_lock_size, long max_branches,
                       const std::vector<Pair>* islands, const std::vector<Pair>& prior_deferred,
                       bool allow_mapped_seed) {
    GrowResult result;
    GrowProfile& prof = result.profile;
    const int NR = R.n, NP = P.n;
    std::vector<uint8_t> locked_p(NP, 0);
    for (int r = 0; r < NR; ++r)
        if (mapping[r] >= 0) locked_p[mapping[r]] = 1;
    std::vector<int> island_of(NR, -1);
    if (islands)
        for (const auto& item : *islands) island_of[item.first] = item.second;

    std::vector<Cand> cands;
    if (mapping[seed] >= 0) {
        if (!allow_mapped_seed) {
            prof.result = "already_mapped";
            return result;
        }
        prof.seed_targets = 1;
        prof.seed_groups = 1;
        std::vector<int> raw(NR, -1);
        raw[seed] = mapping[seed];
        Cand c;
        make_cand(raw, {}, {seed}, 1, {}, NP, c);
        cands.push_back(std::move(c));
    } else {
        std::vector<int> targets;
        for (int v = 0; v < NP; ++v)
            if (!locked_p[v] && P.ecode[v] == R.ecode[seed]) targets.push_back(v);
        prof.seed_targets = (long)targets.size();
        // groups keyed by (element, orbit), ordered by Python str(key)
        std::map<int, std::vector<int>> by_orbit;
        for (int v : targets) by_orbit[P.orbit[v]].push_back(v);
        std::vector<std::pair<std::string, std::vector<int>>> groups;
        for (auto& item : by_orbit) {
            std::string repr = "('" + R.elem[seed] + "', " + std::to_string(item.first) + ")";
            groups.push_back({repr, item.second});
        }
        std::stable_sort(groups.begin(), groups.end(),
                         [](const auto& a, const auto& b) { return a.first < b.first; });
        prof.seed_groups = (long)groups.size();
        for (auto& g : groups) {
            std::vector<int> raw(NR, -1);
            raw[seed] = g.second[0];
            Cand c;
            if (g.second.size() > 1)
                make_cand(raw, {make_block({seed}, g.second, false)}, {}, 1, {}, NP, c);
            else
                make_cand(raw, {}, {}, 1, {}, NP, c);
            cands.push_back(std::move(c));
        }
    }
    if (cands.empty()) {
        prof.result = "no_initial_cands";
        return result;
    }
    std::vector<uint8_t> fragment(NR, 0);
    fragment[seed] = 1;
    int fragment_size = 1;
    std::vector<uint8_t> used_edge(NR * (size_t)NR, 0);
    auto edge_used = [&](int a, int b) -> uint8_t& { return used_edge[(size_t)std::min(a, b) * NR + std::max(a, b)]; };
    std::vector<Pair> deferred = prior_deferred;
    for (auto& e : deferred)
        if (e.first > e.second) std::swap(e.first, e.second);
    std::sort(deferred.begin(), deferred.end());
    deferred.erase(std::unique(deferred.begin(), deferred.end()), deferred.end());

    std::priority_queue<HeapItem, std::vector<HeapItem>, std::greater<HeapItem>> heap;
    auto push_edges_from = [&](int atom) {
        for (int nb : R.adj[atom]) {
            if (fragment[nb]) continue;
            if (edge_used(atom, nb)) continue;
            double w = R.wbo(atom, nb);
            if (w >= graph_floor) heap.push({-w, atom, nb});
        }
    };
    push_edges_from(seed);
    prof.init_cands = (long)cands.size();
    prof.initial_heap = (long)heap.size();
    prof.max_cands_before = prof.max_cands_after = (long)cands.size();
    prof.max_fragment_size = 1;
    prof.max_heap_len = (long)heap.size();
    if ((long)cands.size() > max_branches) {
        prof.result = "live_branch_cap";
        prof.final_cands = (long)cands.size();
        prof.final_fragment_size = fragment_size;
        prof.branches = (long)cands.size();
        result.capped = true;
        result.cap = {(long)cands.size(), max_branches};
        return result;
    }
    Canonicalizer canon(&P, mapping);
    long certificate_calls = 0;

    while (!heap.empty()) {
        HeapItem item = heap.top();
        heap.pop();
        prof.heap_pops++;
        int u = item.u, n = item.n;
        double wbo = -item.neg_w;
        if (edge_used(u, n)) { prof.stale_pops++; continue; }
        edge_used(u, n) = 1;
        if (fragment[n]) { prof.fragment_skip_pops++; continue; }
        bool n_in_mapping = mapping[n] >= 0;
        prof.extend_calls++;
        if (n_in_mapping) prof.merge_calls++; else prof.free_extend_calls++;
        long old_count = (long)cands.size();
        prof.max_cands_before = std::max(prof.max_cands_before, old_count);

        // dedupe fragment: fragment | {n} | n's whole island (if any)
        std::vector<uint8_t> dedupe_fragment = fragment;
        dedupe_fragment[n] = 1;
        if (n_in_mapping && islands && island_of[n] >= 0) {
            int iid = island_of[n];
            for (const auto& it : *islands)
                if (it.second == iid) dedupe_fragment[it.first] = 1;
        }
        // dedupe edges: deferred | one-hop frontier of dedupe_fragment
        std::vector<Pair> dedupe_edges = deferred;
        for (int a = 0; a < NR; ++a) {
            if (!dedupe_fragment[a]) continue;
            for (int nb : R.adj[a]) {
                if (dedupe_fragment[nb]) continue;
                if (R.wbo(a, nb) >= graph_floor)
                    dedupe_edges.push_back({std::min(a, nb), std::max(a, nb)});
            }
        }
        std::sort(dedupe_edges.begin(), dedupe_edges.end());
        dedupe_edges.erase(std::unique(dedupe_edges.begin(), dedupe_edges.end()), dedupe_edges.end());

        // context
        Context ctx;
        ctx.R = &R;
        ctx.P = &P;
        ctx.fragment_old = fragment;
        ctx.n = n;
        ctx.mapping = &mapping;
        ctx.locked_p = locked_p;
        ctx.iso_tol = iso_tol;
        ctx.islands = islands;
        ctx.island_of = &island_of;
        ctx.deferred_edges = deferred;
        ctx.anchor_u = u;
        ctx.anchor_wbo = wbo;
        ctx.has_anchor_wbo = true;
        ctx.dedupe_edges = std::move(dedupe_edges);
        for (int a = 0; a < NR; ++a)
            if (fragment[a] && R.has_edge(a, n)) ctx.bonded_in_frag.push_back(a);
        std::vector<Cand> new_cands;
        if (!ctx.bonded_in_frag.empty()) {
            for (int a : ctx.bonded_in_frag) ctx.r_wbos.push_back(R.wbo(a, n));
            ctx.strict_r = fragment[u] ? u : -1;
            ctx.strict_w = wbo;
            if (!n_in_mapping || !islands || island_of[n] < 0) {
                ctx.island_atoms = {n};
            } else {
                int iid = island_of[n];
                for (const auto& it : *islands)
                    if (it.second == iid && !fragment[it.first]) ctx.island_atoms.push_back(it.first);
            }
            ctx.sig_fragment = fragment;
            ctx.sig_fragment[n] = 1;
            for (int a : ctx.island_atoms) ctx.sig_fragment[a] = 1;
            new_cands = extend_sym_cands(cands, ctx, canon, certificate_calls);
        }
        if ((long)new_cands.size() > max_branches) {
            prof.result = "live_branch_cap";
            prof.final_cands = (long)new_cands.size();
            prof.final_fragment_size = fragment_size + 1;
            prof.branches = (long)new_cands.size();
            prof.certificate_calls = certificate_calls;
            result.capped = true;
            result.cap = {(long)new_cands.size(), max_branches};
            return result;
        }
        if (!new_cands.empty()) {
            prof.commits++;
            prof.max_cands_after = std::max(prof.max_cands_after, (long)new_cands.size());
            cands = std::move(new_cands);
            std::vector<int> added{n};
            if (n_in_mapping && islands && island_of[n] >= 0) {
                int iid = island_of[n];
                for (const auto& it : *islands)
                    if (it.second == iid && !fragment[it.first] && it.first != n) added.push_back(it.first);
            }
            for (int a : added) {
                if (!fragment[a]) {
                    fragment[a] = 1;
                    fragment_size++;
                }
            }
            for (int a : added) push_edges_from(a);
            prof.max_fragment_size = std::max(prof.max_fragment_size, (long)fragment_size);
            prof.max_heap_len = std::max(prof.max_heap_len, (long)heap.size());
        } else {
            prof.deferred++;
            Pair e{std::min(u, n), std::max(u, n)};
            auto pos = std::lower_bound(deferred.begin(), deferred.end(), e);
            if (pos == deferred.end() || *pos != e) deferred.insert(pos, e);
        }
    }
    // saturation quotient
    {
        Context ctx;
        ctx.R = &R;
        ctx.P = &P;
        ctx.mapping = &mapping;
        ctx.locked_p = locked_p;
        ctx.iso_tol = iso_tol;
        cands = dedup_sym_cands(cands, ctx, canon, fragment, deferred, certificate_calls);
    }
    prof.certificate_calls = certificate_calls;
    prof.final_fragment_size = fragment_size;
    if (cands.empty() || fragment_size < min_lock_size) {
        prof.result = cands.empty() ? "no_cands" : "too_small";
        prof.final_cands = (long)cands.size();
        return result;
    }
    auto emit = [&](const Cand& c) {
        IsoOut iso;
        for (int r : c.mapped) iso.mapping.push_back({r, c.img[r]});
        iso.deferred_edges = deferred;
        for (int r = 0; r < NR; ++r)
            if (fragment[r]) iso.fragment.push_back(r);
        iso.cand = c;
        result.isos.push_back(std::move(iso));
    };
    // _set_unique
    bool set_unique = true;
    for (const auto& c : cands)
        if (c.has_open_choice()) { set_unique = false; break; }
    if (set_unique && cands.size() > 1) {
        for (size_t i = 1; i < cands.size(); ++i)
            if (cands[i].img != cands[0].img) { set_unique = false; break; }
    }
    prof.final_cands = (long)cands.size();
    if (set_unique) {
        prof.result = "success";
        prof.branches = 1;
        emit(cands[0]);
        return result;
    }
    if ((long)cands.size() > max_branches) {
        prof.result = "subtree_branch_cap";
        prof.branches = (long)cands.size();
        result.capped = true;
        result.cap = {(long)cands.size(), max_branches};
        return result;
    }
    prof.result = "branched";
    prof.branches = (long)cands.size();
    for (const auto& c : cands) emit(c);
    return result;
}

// ---------------------------------------------------------------------------
// Python conversions
// ---------------------------------------------------------------------------

std::unordered_map<std::string, int>& element_table() {
    static std::unordered_map<std::string, int> table;
    return table;
}

int element_code(const std::string& e) {
    auto& table = element_table();
    auto it = table.find(e);
    if (it == table.end()) it = table.emplace(e, (int)table.size()).first;
    return it->second;
}

void fill_graph(Graph& g, const std::vector<std::string>& elements,
                const std::vector<std::vector<double>>& wbo, double bond_cut,
                const std::vector<Pair>& edges) {
    g.n = (int)elements.size();
    g.elem = elements;
    g.ecode.resize(g.n);
    for (int i = 0; i < g.n; ++i) g.ecode[i] = element_code(elements[i]);
    g.adj.assign(g.n, {});
    g.edge.assign((size_t)g.n * g.n, 0);
    g.w.assign((size_t)g.n * g.n, 0.0);
    for (int i = 0; i < g.n; ++i)
        for (int j = 0; j < g.n; ++j) g.w[(size_t)i * g.n + j] = wbo[i][j];
    g.bond_cut = bond_cut;
    for (const auto& e : edges) {
        g.edge[(size_t)e.first * g.n + e.second] = 1;
        g.edge[(size_t)e.second * g.n + e.first] = 1;
        g.adj[e.first].push_back(e.second);
        g.adj[e.second].push_back(e.first);
    }
    for (auto& list : g.adj) std::sort(list.begin(), list.end());
}

struct PySource {
    Source g;
    PySource(const std::vector<std::string>& elements, const std::vector<std::vector<double>>& wbo,
             double bond_cut, const std::vector<Pair>& edges) {
        fill_graph(g, elements, wbo, bond_cut, edges);
    }
};

struct PyTarget {
    Target g;
    PyTarget(const std::vector<std::string>& elements, const std::vector<std::vector<double>>& wbo,
             double bond_cut, const std::vector<Pair>& edges, const std::vector<int>& orbits,
             const std::vector<std::tuple<int, int, int>>& pair_buckets, int zero_bucket) {
        fill_graph(g, elements, wbo, bond_cut, edges);
        g.orbit = orbits;
        int max_orbit = 0;
        for (int o : orbits) max_orbit = std::max(max_orbit, o);
        g.orbit_size.assign(max_orbit + 1, 0);
        for (int o : orbits) g.orbit_size[o]++;
        g.zero_bucket = zero_bucket;
        g.bucket.assign((size_t)g.n * g.n, zero_bucket);
        std::vector<std::tuple<int, int, int>> sorted_pairs = pair_buckets;
        for (auto& t : sorted_pairs)
            if (std::get<0>(t) > std::get<1>(t)) std::swap(std::get<0>(t), std::get<1>(t));
        std::sort(sorted_pairs.begin(), sorted_pairs.end());
        for (const auto& t : sorted_pairs) {
            int a = std::get<0>(t), b = std::get<1>(t), bucket = std::get<2>(t);
            g.bucket[(size_t)a * g.n + b] = bucket;
            g.bucket[(size_t)b * g.n + a] = bucket;
        }
        // nauty base graph: atoms then one vertex per nonzero-bucket pair in
        // sorted pair order (canonical._CandidateAutomorphismCanonicalizer)
        std::vector<Pair> nauty_edges;
        int next = g.n;
        for (const auto& t : sorted_pairs) {
            int a = std::get<0>(t), b = std::get<1>(t), bucket = std::get<2>(t);
            if (bucket == zero_bucket) continue;
            nauty_edges.push_back({a, next});
            nauty_edges.push_back({b, next});
            g.edge_vertex_bucket.push_back(bucket);
            g.edge_cells[bucket].push_back(next);
            ++next;
        }
        g.n_vertices = next;
        g.nauty = std::make_unique<DenseGraph>(next, nauty_edges);
        for (const auto& item : g.edge_cells)
            g.edge_cell_colour.push_back(g.colours.intern({1, (int64_t)item.first}));
        int max_code = 0;
        for (int c : g.ecode) max_code = std::max(max_code, c);
        g.same_element.assign((size_t)max_code + 1, {});
        for (int v = 0; v < g.n; ++v) g.same_element[g.ecode[v]].push_back(v);
    }
};

py::dict block_dict(const Block& b, bool automorph) {
    py::dict d;
    py::list r, p;
    for (int x : b.r) r.append(x);
    for (int x : b.p) p.append(x);
    d["r_atoms"] = r;
    d["p_atoms"] = p;
    if (automorph) {
        d["extendable"] = false;
        d["open"] = false;
        d["assignments"] = "exact_group";
        d["source"] = "exact_automorph_group";
    } else {
        d["extendable"] = b.extendable;
        d["open"] = b.open();
        size_t n = b.p.size(), k = b.r.size();
        std::string expr;
        if (k <= 0 || n <= 1) expr = "1";
        else if (k == 1) expr = std::to_string(n);
        else if (k == n) expr = std::to_string(n) + "!";
        else expr = "P(" + std::to_string(n) + "," + std::to_string(k) + ")";
        d["assignments"] = expr;
    }
    return d;
}

py::dict symmetry_state(const Cand& c) {
    py::dict item;
    py::dict witness;
    for (int r : c.mapped) witness[py::int_(r)] = c.img[r];
    item["witness"] = witness;
    py::list blocks;
    item["blocks"] = blocks;
    py::list exact;
    for (int x : c.exact_fixed) exact.append(x);
    item["exact_fixed"] = exact;
    item["multiplicity"] = c.mult;
    py::list automorph;
    for (const auto& b : c.automorph) automorph.append(block_dict(b, true));
    item["automorph_blocks"] = automorph;
    for (const auto& b : c.blocks) blocks.append(block_dict(b, false));
    for (auto entry : automorph) blocks.append(entry);
    return item;
}

py::object py_grow_island(const PySource& source, const PyTarget& target, int seed,
                          const std::vector<int>& mapping, double graph_floor, double iso_tol,
                          int min_lock_size, long max_branches, py::object islands_obj,
                          const std::vector<Pair>& prior_deferred, bool allow_mapped_seed) {
    std::vector<Pair> islands;
    const std::vector<Pair>* islands_ptr = nullptr;
    if (!islands_obj.is_none()) {
        islands = islands_obj.cast<std::vector<Pair>>();
        islands_ptr = &islands;
    }
    GrowResult res;
    {
        py::gil_scoped_release release;
        res = grow_island(source.g, target.g, seed, mapping, graph_floor, iso_tol, min_lock_size,
                          max_branches, islands_ptr, prior_deferred, allow_mapped_seed);
    }
    py::dict out;
    py::list isos;
    for (const auto& iso : res.isos) {
        py::dict d;
        py::dict m;
        for (const auto& pr : iso.mapping) m[py::int_(pr.first)] = pr.second;
        d["mapping"] = m;
        d["deferred_edges"] = iso.deferred_edges;
        d["fragment"] = iso.fragment;
        d["symmetry"] = symmetry_state(iso.cand);
        isos.append(d);
    }
    out["isos"] = isos;
    out["capped"] = res.capped;
    out["cap_count"] = res.cap.count;
    out["cap_limit"] = res.cap.limit;
    py::dict prof;
    const GrowProfile& p = res.profile;
    prof["seed_targets"] = p.seed_targets;
    prof["seed_groups"] = p.seed_groups;
    prof["init_cands"] = p.init_cands;
    prof["initial_heap"] = p.initial_heap;
    prof["heap_pops"] = p.heap_pops;
    prof["stale_pops"] = p.stale_pops;
    prof["fragment_skip_pops"] = p.fragment_skip_pops;
    prof["extend_calls"] = p.extend_calls;
    prof["max_cands_before"] = p.max_cands_before;
    prof["max_cands_after"] = p.max_cands_after;
    prof["max_fragment_size"] = p.max_fragment_size;
    prof["max_heap_len"] = p.max_heap_len;
    prof["commits"] = p.commits;
    prof["deferred"] = p.deferred;
    prof["merge_calls"] = p.merge_calls;
    prof["free_extend_calls"] = p.free_extend_calls;
    prof["certificate_calls"] = p.certificate_calls;
    prof["final_cands"] = p.final_cands;
    prof["final_fragment_size"] = p.final_fragment_size;
    prof["branches"] = p.branches;
    prof["result"] = p.result;
    out["profile"] = prof;
    return out;
}

}  // namespace

// native/src/freeze.cpp
void register_freeze(pybind11::module_&);
// native/src/autgrp.cpp
void register_autgrp(pybind11::module_&);
// native/src/repair.cpp
void register_repair(pybind11::module_&);

PYBIND11_MODULE(_engine, mod) {
    mod.doc() = "Native AAM growth engine";
    register_freeze(mod);
    register_autgrp(mod);
    register_repair(mod);
    py::class_<DenseGraph>(mod, "DenseGraph")
        .def(py::init<int, const std::vector<Pair>&>())
        .def("certificate", [](DenseGraph& g, const std::vector<std::vector<int>>& cells) {
            std::string s = g.certificate(cells);
            return py::bytes(s);
        })
        .def_readonly("n", &DenseGraph::n);
    py::class_<PySource>(mod, "SourceGraph")
        .def(py::init<const std::vector<std::string>&, const std::vector<std::vector<double>>&,
                      double, const std::vector<Pair>&>());
    py::class_<PyTarget>(mod, "TargetGraph")
        .def(py::init<const std::vector<std::string>&, const std::vector<std::vector<double>>&,
                      double, const std::vector<Pair>&, const std::vector<int>&,
                      const std::vector<std::tuple<int, int, int>>&, int>())
        .def("memo_stats", [](const PyTarget& t) {
            return py::make_tuple(t.g.nauty->memo_hits, t.g.nauty->memo_misses,
                                  (long)t.g.nauty->memo.size());
        });
    mod.def("grow_island", &py_grow_island, py::arg("source"), py::arg("target"), py::arg("seed"),
            py::arg("mapping"), py::arg("graph_floor"), py::arg("iso_tol"),
            py::arg("min_lock_size"), py::arg("max_branches"), py::arg("islands"),
            py::arg("prior_deferred_edges"), py::arg("allow_mapped_seed"));
}
