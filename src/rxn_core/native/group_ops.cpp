// Exact permutation bookkeeping. Search choices and chemical semantics live in Python.
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <algorithm>
#include <map>
#include <numeric>
#include <optional>
#include <set>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace py = pybind11;
using Images = std::vector<int>;
using Fragment = std::pair<int, Images>;
using Bond = std::pair<int, int>;
using Key = std::tuple<Images, Images, std::vector<Fragment>, std::vector<Bond>>;
using State = std::pair<Images, Images>;

// Only a closed, permanently unobservable set can be collapsed. Include every
// point that any recorded generator can transport to an observed point. This
// over-approximates chronological reachability and is safe for every stage.
static Images observation_labels(int degree,
        const std::vector<std::vector<Images>>& stages,
        const std::optional<Images>& observed) {
    Images labels(degree);
    std::iota(labels.begin(), labels.end(), 0);
    if (!observed) return labels;
    Images parent = labels;
    auto root = [&](int a) {
        while (parent[a] != a) {
            parent[a] = parent[parent[a]];
            a = parent[a];
        }
        return a;
    };
    for (const auto& stage : stages) for (const auto& g : stage)
        for (int a = 0; a < degree; ++a) parent[root(a)] = root(g[a]);
    std::vector<bool> visible(degree, false);
    for (int a : *observed) {
        if (a < 0 || a >= degree) throw std::invalid_argument("observed atom outside frame");
        visible[root(a)] = true;
    }
    for (int a = 0; a < degree; ++a)
        if (!visible[root(a)]) labels[a] = -1;
    return labels;
}

static void check_permutation(const Images& p) {
    Images ordered = p;
    std::sort(ordered.begin(), ordered.end());
    for (size_t i = 0; i < ordered.size(); ++i)
        if (ordered[i] != static_cast<int>(i))
            throw std::invalid_argument("expected a permutation of its complete frame");
}

static Images select_sorted(const Images& images, const Images& positions) {
    Images selected;
    selected.reserve(positions.size());
    for (int position : positions) selected.push_back(images[position]);
    std::sort(selected.begin(), selected.end());
    return selected;
}

static Key occupation_key(const Images& images, const Images& attachments,
                          const std::vector<Fragment>& fragments,
                          const std::vector<Bond>& bonds) {
    Images covered = images;
    std::sort(covered.begin(), covered.end());
    std::vector<Fragment> parts;
    parts.reserve(fragments.size());
    for (const auto& [label, positions] : fragments)
        parts.emplace_back(label, select_sorted(images, positions));
    std::sort(parts.begin(), parts.end());
    std::vector<Bond> edges;
    edges.reserve(bonds.size());
    for (const auto& [a, b] : bonds)
        edges.push_back(std::minmax(images[a], images[b]));
    std::sort(edges.begin(), edges.end());
    return {std::move(covered), select_sorted(images, attachments), std::move(parts), std::move(edges)};
}

// Same chronological stage closure and first witness as the Python algorithm.
// This enumerates distinct fragment relations, never the complete bijection group.
static std::vector<State> occupation_orbit(
        const Images& witness, int degree, const std::vector<std::vector<Images>>& stages,
        const Images& attachments, const std::vector<Fragment>& fragments,
        const std::vector<Bond>& bonds, int limit,
        const std::optional<Images>& observed = std::nullopt) {
    if (degree < 0) throw std::invalid_argument("negative permutation degree");
    for (int image : witness)
        if (image < 0 || image >= degree) throw std::invalid_argument("image outside frame");
    auto check_position = [&](int i) {
        if (i < 0 || static_cast<size_t>(i) >= witness.size())
            throw std::invalid_argument("source position outside witness");
    };
    for (int i : attachments) check_position(i);
    for (const auto& fragment : fragments) for (int i : fragment.second) check_position(i);
    for (const auto& bond : bonds) { check_position(bond.first); check_position(bond.second); }
    for (const auto& stage : stages) for (const auto& g : stage) {
        if (g.size() != static_cast<size_t>(degree)) throw std::invalid_argument("generator frame mismatch");
        check_permutation(g);
    }
    Images identity(degree);
    std::iota(identity.begin(), identity.end(), 0);
    const Images labels = observation_labels(degree, stages, observed);
    auto projected = [&](const Images& images) {
        Images value;
        value.reserve(images.size());
        for (int image : images) value.push_back(labels[image]);
        return value;
    };
    std::vector<State> states{{witness, identity}};
    std::map<Key, size_t> seen;
    seen.emplace(occupation_key(projected(witness), attachments, fragments, bonds), 0);
    std::map<Images, size_t> evaluated;
    for (const auto& stage : stages) {
        std::vector<size_t> already;
        for (const auto& g : stage) already.push_back(evaluated[g]);
        for (size_t current = 0; current < states.size(); ++current) {
            // Appending states may reallocate; do not hold references into states.
            const State state = states[current];
            const Images state_projection = projected(state.first);
            for (size_t gi = 0; gi < stage.size(); ++gi) {
                if (current < already[gi]) continue;
                const auto& g = stage[gi];
                Images moved;
                moved.reserve(witness.size());
                for (int image : state.first) moved.push_back(g[image]);
                Images moved_projection = projected(moved);
                if (moved_projection == state_projection) continue;
                Key key = occupation_key(moved_projection, attachments, fragments, bonds);
                if (seen.find(key) != seen.end()) continue;
                if (limit >= 0 && states.size() >= static_cast<size_t>(limit))
                    throw std::length_error("fragment occupation limit exceeded");
                Images action;
                action.reserve(degree);
                for (int image : state.second) action.push_back(g[image]);
                seen.emplace(std::move(key), states.size());
                states.emplace_back(std::move(moved), std::move(action));
            }
        }
        for (const auto& g : stage) evaluated[g] = states.size();
    }
    // Preserve discovery order; the Python caller merges families and sorts final keys.
    return states;
}

static py::tuple project_generators(const std::vector<Images>& raw,
                                    const Images& atom_labels, const Images& vertex_indices) {
    if (atom_labels.size() != vertex_indices.size()) throw std::invalid_argument("atom index size mismatch");
    const int degree = atom_labels.empty() ? 0 : *std::max_element(atom_labels.begin(), atom_labels.end()) + 1;
    std::map<int, int> atoms;
    for (size_t i = 0; i < atom_labels.size(); ++i) atoms.emplace(vertex_indices[i], atom_labels[i]);
    Images identity(degree);
    std::iota(identity.begin(), identity.end(), 0);
    std::set<Images> seen;
    std::vector<Images> output;
    {
        py::gil_scoped_release release;
        for (const auto& g : raw) {
            Images p = identity;
            for (size_t i = 0; i < atom_labels.size(); ++i) {
                const int vertex = vertex_indices[i];
                if (vertex < 0 || static_cast<size_t>(vertex) >= g.size())
                    throw std::invalid_argument("atom vertex outside generator frame");
                auto image = atoms.find(g[vertex]);
                if (image == atoms.end()) throw std::invalid_argument("candidate automorphism mixed atom/edge vertices");
                p.at(atom_labels[i]) = image->second;
            }
            if (p != identity && seen.insert(p).second) output.push_back(std::move(p));
        }
    }
    py::tuple result(output.size());
    for (size_t i = 0; i < output.size(); ++i) result[i] = py::tuple(py::cast(output[i]));
    return result;
}

static py::tuple conjugate_generators(const std::vector<Images>& generators, const Images& action) {
    check_permutation(action);
    std::vector<Images> output;
    {
        py::gil_scoped_release release;
        for (const auto& g : generators) {
            check_permutation(g);
            size_t degree = std::max(g.size(), action.size());
            Images moved(degree);
            auto image = [&](size_t i) { return i < action.size() ? action[i] : static_cast<int>(i); };
            for (size_t i = 0; i < degree; ++i)
                moved[image(i)] = image(i < g.size() ? g[i] : i);
            output.push_back(std::move(moved));
        }
    }
    py::tuple result(output.size());
    for (size_t i = 0; i < output.size(); ++i) result[i] = py::tuple(py::cast(output[i]));
    return result;
}

PYBIND11_MODULE(_group_ops, m) {
    m.def("occupation_orbit", &occupation_orbit,
          py::arg("witness"), py::arg("degree"), py::arg("stages"),
          py::arg("attachments"), py::arg("fragments"), py::arg("bonds"),
          py::arg("limit"), py::arg("observed_atoms") = py::none(),
          py::call_guard<py::gil_scoped_release>());
    m.def("project_generators", &project_generators);
    m.def("conjugate_generators", &conjugate_generators);
    py::register_exception<std::length_error>(m, "OccupationLimitExceeded");
}
