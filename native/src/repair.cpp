// Exact kernel for the per-group search inside
// rxn_core.alignment.branch.symmetry_repair_mapping.
//
// One call performs, for one target group, exactly what the Python loop
// body does: the breadth-first enumeration of the group's target tuples under
// the supplied generators (with the same evaluation cap bookkeeping), the
// lexicographic ordering of the states, the scoring of the first `budget`
// states with the same elementwise IEEE operations and the same strict
// left-to-right double summation, and the `score < best_score` scan.  The
// 12-place rounding is CPython's own float rounding (dtoa mode 3 followed by
// a correctly rounded strtod), so every returned score equals the Python
// `round(total, 12)`.
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <deque>
#include <string>
#include <unordered_set>
#include <vector>

namespace py = pybind11;

namespace {

struct StateHash {
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

// Python's round(x, 12) for finite doubles: CPython formats with
// _Py_dg_dtoa(x, 3, 12) and re-reads the decimal with _Py_dg_strtod;
// PyOS_double_to_string('f', 12) / PyOS_string_to_double perform the same two
// correctly rounded conversions.
double round12(double x) {
    if (!std::isfinite(x)) return x;
    char* text = PyOS_double_to_string(x, 'f', 12, 0, nullptr);
    if (text == nullptr) throw py::error_already_set();
    double out = PyOS_string_to_double(text, nullptr, nullptr);
    PyMem_Free(text);
    if (out == -1.0 && PyErr_Occurred()) throw py::error_already_set();
    return out;
}

template <typename T>
using Arr = py::array_t<T, py::array::c_style | py::array::forcecast>;

py::tuple repair_group(Arr<int64_t> current_images, Arr<int64_t> rs, Arr<int64_t> generators,
                       Arr<int64_t> pair_i, Arr<int64_t> pair_j, Arr<double> pair_wbo_R,
                       Arr<double> pair_threshold, Arr<bool> pair_r_active, Arr<double> wbo_P,
                       double bond_floor, long evaluated, long max_evals, long best_count,
                       double best_total) {
    auto img0 = current_images.unchecked<1>();
    auto rs_v = rs.unchecked<1>();
    auto gens = generators.unchecked<2>();
    auto pi = pair_i.unchecked<1>();
    auto pj = pair_j.unchecked<1>();
    auto wR = pair_wbo_R.unchecked<1>();
    auto thr = pair_threshold.unchecked<1>();
    auto active = pair_r_active.unchecked<1>();
    auto wP = wbo_P.unchecked<2>();
    const size_t k = (size_t)rs_v.shape(0);
    const size_t n_gen = (size_t)gens.shape(0);
    const size_t n_pairs = (size_t)pi.shape(0);

    // breadth-first orbit enumeration (states is a set; order fixed below)
    std::vector<int64_t> state(k);
    for (size_t t = 0; t < k; ++t) state[t] = img0(rs_v(t));
    std::unordered_set<std::vector<int64_t>, StateHash> states;
    states.insert(state);
    std::deque<std::vector<int64_t>> queue;
    queue.push_back(state);
    bool capped = false;
    std::vector<int64_t> image(k);
    while (!queue.empty() && evaluated < max_evals) {
        state = std::move(queue.front());
        queue.pop_front();
        for (size_t g = 0; g < n_gen; ++g) {
            for (size_t t = 0; t < k; ++t) image[t] = gens(g, state[t]);
            if (states.count(image)) continue;
            states.insert(image);
            queue.push_back(image);
            if ((long)states.size() + evaluated >= max_evals) {
                capped = !queue.empty();
                break;
            }
        }
    }
    std::vector<std::vector<int64_t>> ordered(states.begin(), states.end());
    std::sort(ordered.begin(), ordered.end());
    long budget = max_evals - evaluated;
    if ((long)ordered.size() > budget) capped = true;
    size_t to_score = (size_t)std::max<long>(budget, 0);
    if (to_score > ordered.size()) to_score = ordered.size();

    // scoring: the same per-pair operations as local_scores, summed strictly
    // left to right (np.add.accumulate's last column)
    std::vector<int64_t> images(img0.data(0), img0.data(0) + img0.shape(0));
    bool have_best = false;
    std::vector<int64_t> best_state;
    for (size_t s = 0; s < to_score; ++s) {
        const auto& st = ordered[s];
        for (size_t t = 0; t < k; ++t) images[rs_v(t)] = st[t];
        long count = 0;
        double total = 0.0;
        for (size_t q = 0; q < n_pairs; ++q) {
            double w_p = wP(images[pi(q)], images[pj(q)]);
            double difference = wR(q) - w_p;
            double magnitude = std::fabs(difference);
            bool changed = magnitude >= thr(q);
            double contribution;
            if (changed) {
                contribution = magnitude;
            } else if (active(q) || w_p >= bond_floor) {
                double scaled = magnitude * 0.01;
                contribution = scaled;
            } else {
                contribution = 0.0;
            }
            count += changed ? 1 : 0;
            if (q == 0) {
                total = contribution;
            } else {
                double next = total + contribution;
                total = next;
            }
        }
        if (n_pairs == 0) total = 0.0;
        double rounded = round12(total);
        // Python tuple comparison: count first, then the rounded total
        if (count < best_count || (count == best_count && rounded < best_total)) {
            best_count = count;
            best_total = rounded;
            best_state = st;
            have_best = true;
        }
    }
    py::object best = py::none();
    if (have_best) {
        py::tuple t(best_state.size());
        for (size_t i = 0; i < best_state.size(); ++i) t[i] = py::int_(best_state[i]);
        best = std::move(t);
    }
    return py::make_tuple(best, py::make_tuple(py::int_(best_count), py::float_(best_total)),
                          py::int_((long)to_score), py::bool_(capped));
}

// Scores of explicit states (differential testing against local_scores).
py::list repair_scores(Arr<int64_t> current_images, Arr<int64_t> rs, Arr<int64_t> states,
                       Arr<int64_t> pair_i, Arr<int64_t> pair_j, Arr<double> pair_wbo_R,
                       Arr<double> pair_threshold, Arr<bool> pair_r_active, Arr<double> wbo_P,
                       double bond_floor) {
    auto img0 = current_images.unchecked<1>();
    auto rs_v = rs.unchecked<1>();
    auto st = states.unchecked<2>();
    auto pi = pair_i.unchecked<1>();
    auto pj = pair_j.unchecked<1>();
    auto wR = pair_wbo_R.unchecked<1>();
    auto thr = pair_threshold.unchecked<1>();
    auto active = pair_r_active.unchecked<1>();
    auto wP = wbo_P.unchecked<2>();
    const size_t k = (size_t)rs_v.shape(0);
    const size_t n_pairs = (size_t)pi.shape(0);
    std::vector<int64_t> images(img0.data(0), img0.data(0) + img0.shape(0));
    py::list out;
    for (py::ssize_t s = 0; s < st.shape(0); ++s) {
        for (size_t t = 0; t < k; ++t) images[rs_v(t)] = st(s, t);
        long count = 0;
        double total = 0.0;
        for (size_t q = 0; q < n_pairs; ++q) {
            double w_p = wP(images[pi(q)], images[pj(q)]);
            double difference = wR(q) - w_p;
            double magnitude = std::fabs(difference);
            bool changed = magnitude >= thr(q);
            double contribution;
            if (changed) {
                contribution = magnitude;
            } else if (active(q) || w_p >= bond_floor) {
                double scaled = magnitude * 0.01;
                contribution = scaled;
            } else {
                contribution = 0.0;
            }
            count += changed ? 1 : 0;
            if (q == 0) {
                total = contribution;
            } else {
                double next = total + contribution;
                total = next;
            }
        }
        if (n_pairs == 0) total = 0.0;
        out.append(py::make_tuple(py::int_(count), py::float_(round12(total))));
    }
    return out;
}

}  // namespace

void register_repair(py::module_& m) {
    m.def("repair_group", &repair_group, py::arg("current_images"), py::arg("rs"),
          py::arg("generators"), py::arg("pair_i"), py::arg("pair_j"), py::arg("pair_wbo_R"),
          py::arg("pair_threshold"), py::arg("pair_r_active"), py::arg("wbo_P"),
          py::arg("bond_floor"), py::arg("evaluated"), py::arg("max_evals"),
          py::arg("best_count"), py::arg("best_total"));
    m.def("repair_scores", &repair_scores);
    m.def("round12", &round12);
}
