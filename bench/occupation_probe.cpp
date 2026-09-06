// Exact permutation bookkeeping. Search choices and chemical semantics live in Python.
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <algorithm>
#include <map>
#include <numeric>
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

#include <chrono>

// Diagnostic replica of the production orbit walk. The optional cache only
// skips state/generator pairs evaluated in an earlier stage; never search states.
static py::dict probe(const Images& witness, int degree,
        const std::vector<std::vector<Images>>& stages, const Images& attachments,
        const std::vector<Fragment>& fragments, const std::vector<Bond>& bonds,
        bool cache_edges) {
    Images identity(degree);
    std::iota(identity.begin(), identity.end(), 0);
    std::vector<State> states{{witness, identity}};
    std::map<Key, size_t> seen;
    seen.emplace(occupation_key(witness, attachments, fragments, bonds), 0);
    std::map<Images, size_t> evaluated;
    py::list measurements;
    for (size_t si = 0; si < stages.size(); ++si) {
        const auto& stage = stages[si];
        auto started = std::chrono::steady_clock::now();
        size_t before = states.size(), attempts = 0, moved_equal = 0, repeated_key = 0, skipped = 0;
        std::vector<size_t> already;
        for (const auto& g : stage) already.push_back(evaluated[g]);
        for (size_t current = 0; current < states.size(); ++current) {
            const State state = states[current];
            for (size_t gi = 0; gi < stage.size(); ++gi) {
                if (cache_edges && current < already[gi]) { ++skipped; continue; }
                const auto& g = stage[gi];
                ++attempts;
                Images moved;
                moved.reserve(witness.size());
                for (int image : state.first) moved.push_back(g[image]);
                if (moved == state.first) { ++moved_equal; continue; }
                Key key = occupation_key(moved, attachments, fragments, bonds);
                if (seen.find(key) != seen.end()) { ++repeated_key; continue; }
                Images action;
                action.reserve(degree);
                for (int image : state.second) action.push_back(g[image]);
                seen.emplace(std::move(key), states.size());
                states.emplace_back(std::move(moved), std::move(action));
            }
        }
        for (const auto& g : stage) evaluated[g] = states.size();
        py::dict row;
        row["stage"] = si;
        row["states_before"] = before;
        row["states_after"] = states.size();
        row["generators"] = stage.size();
        row["attempts"] = attempts;
        row["unchanged_witness"] = moved_equal;
        row["repeated_occupation"] = repeated_key;
        row["cached_edges"] = skipped;
        row["seconds"] = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - started).count();
        measurements.append(row);
    }
    py::dict result;
    result["stages"] = measurements;
    result["states"] = py::cast(states);
    return result;
}

PYBIND11_MODULE(_occupation_probe, m) { m.def("probe", &probe); }

