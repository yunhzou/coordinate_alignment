#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <algorithm>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace {

using Integer = long long;
using Relation = std::pair<Integer, Integer>;
using RawColor = std::pair<Integer, Integer>;
using NeighborhoodKey = std::tuple<Relation, int>;
using RelationCountKey = std::tuple<int, int, Relation>;

bool read_integer(PyObject* value, Integer& output) {
    output = PyLong_AsLongLong(value);
    return !(output == -1 && PyErr_Occurred());
}

bool read_pair(PyObject* value, std::pair<Integer, Integer>& output) {
    PyObject* sequence = PySequence_Fast(value, "expected an integer pair");
    if (sequence == nullptr) {
        return false;
    }
    if (PySequence_Fast_GET_SIZE(sequence) != 2) {
        Py_DECREF(sequence);
        PyErr_SetString(PyExc_ValueError, "expected an integer pair");
        return false;
    }
    const bool valid =
        read_integer(PySequence_Fast_GET_ITEM(sequence, 0), output.first) &&
        read_integer(PySequence_Fast_GET_ITEM(sequence, 1), output.second);
    Py_DECREF(sequence);
    return valid;
}

std::string pair_repr(const std::pair<Integer, Integer>& value) {
    return "(" + std::to_string(value.first) + ", " +
        std::to_string(value.second) + ")";
}

std::string neighborhood_entry_repr(
        const NeighborhoodKey& key, int count) {
    const auto& [relation, color] = key;
    return "((" + std::to_string(relation.first) + ", " +
        std::to_string(relation.second) + "), " +
        std::to_string(color) + "), " + std::to_string(count) + ")";
}

std::string tuple_repr(const std::vector<std::string>& entries) {
    if (entries.empty()) {
        return "()";
    }
    std::string result = "(";
    for (std::size_t index = 0; index < entries.size(); ++index) {
        if (index != 0) {
            result += ", ";
        }
        result += entries[index];
    }
    if (entries.size() == 1) {
        result += ",";
    }
    result += ")";
    return result;
}

std::vector<int> compact(const std::vector<std::string>& values) {
    std::set<std::string> ordered(values.begin(), values.end());
    std::map<std::string, int> classes;
    int index = 0;
    for (const auto& value : ordered) {
        classes.emplace(value, index++);
    }
    std::vector<int> colors;
    colors.reserve(values.size());
    for (const auto& value : values) {
        colors.push_back(classes.at(value));
    }
    return colors;
}

PyObject* integer_pair(const Integer left, const Integer right) {
    return Py_BuildValue("(LL)", left, right);
}

PyObject* color_counts_object(const std::map<int, long long>& counts) {
    PyObject* result = PyTuple_New(static_cast<Py_ssize_t>(counts.size()));
    if (result == nullptr) {
        return nullptr;
    }
    Py_ssize_t index = 0;
    for (const auto& [color, count] : counts) {
        PyObject* item = Py_BuildValue("(iL)", color, count);
        if (item == nullptr) {
            Py_DECREF(result);
            return nullptr;
        }
        PyTuple_SET_ITEM(result, index++, item);
    }
    return result;
}

PyObject* paired_mapping_invariant(PyObject*, PyObject* args) {
    PyObject* raw_colors_object = nullptr;
    PyObject* zero_relation_object = nullptr;
    PyObject* active_relations_object = nullptr;
    if (!PyArg_ParseTuple(
            args, "OOO", &raw_colors_object, &zero_relation_object,
            &active_relations_object)) {
        return nullptr;
    }

    PyObject* raw_sequence = PySequence_Fast(
        raw_colors_object, "raw colors must be a sequence");
    if (raw_sequence == nullptr) {
        return nullptr;
    }
    const Py_ssize_t size = PySequence_Fast_GET_SIZE(raw_sequence);
    std::vector<RawColor> raw_colors(static_cast<std::size_t>(size));
    for (Py_ssize_t index = 0; index < size; ++index) {
        if (!read_pair(
                PySequence_Fast_GET_ITEM(raw_sequence, index),
                raw_colors[static_cast<std::size_t>(index)])) {
            Py_DECREF(raw_sequence);
            return nullptr;
        }
    }
    Py_DECREF(raw_sequence);

    Relation zero_relation;
    if (!read_pair(zero_relation_object, zero_relation)) {
        return nullptr;
    }

    std::vector<std::vector<std::pair<int, Relation>>> active_by_position(
        static_cast<std::size_t>(size));
    std::map<std::pair<int, int>, Relation> relations;
    PyObject* active_sequence = PySequence_Fast(
        active_relations_object, "active relations must be a sequence");
    if (active_sequence == nullptr) {
        return nullptr;
    }
    for (Py_ssize_t index = 0;
         index < PySequence_Fast_GET_SIZE(active_sequence); ++index) {
        PyObject* item = PySequence_Fast(
            PySequence_Fast_GET_ITEM(active_sequence, index),
            "active relation must contain four integers");
        if (item == nullptr) {
            Py_DECREF(active_sequence);
            return nullptr;
        }
        if (PySequence_Fast_GET_SIZE(item) != 4) {
            Py_DECREF(item);
            Py_DECREF(active_sequence);
            PyErr_SetString(
                PyExc_ValueError, "active relation must contain four integers");
            return nullptr;
        }
        Integer raw_left, raw_right, source_relation, target_relation;
        const bool valid =
            read_integer(PySequence_Fast_GET_ITEM(item, 0), raw_left) &&
            read_integer(PySequence_Fast_GET_ITEM(item, 1), raw_right) &&
            read_integer(PySequence_Fast_GET_ITEM(item, 2), source_relation) &&
            read_integer(PySequence_Fast_GET_ITEM(item, 3), target_relation);
        Py_DECREF(item);
        if (!valid) {
            Py_DECREF(active_sequence);
            return nullptr;
        }
        if (raw_left < 0 || raw_right < 0 || raw_left >= size ||
                raw_right >= size || raw_left == raw_right) {
            Py_DECREF(active_sequence);
            PyErr_SetString(PyExc_ValueError, "invalid active relation endpoints");
            return nullptr;
        }
        const int left = static_cast<int>(std::min(raw_left, raw_right));
        const int right = static_cast<int>(std::max(raw_left, raw_right));
        const Relation relation{source_relation, target_relation};
        relations[{left, right}] = relation;
    }
    Py_DECREF(active_sequence);
    for (const auto& [endpoints, relation] : relations) {
        const auto [left, right] = endpoints;
        active_by_position[static_cast<std::size_t>(left)].push_back(
            {right, relation});
        active_by_position[static_cast<std::size_t>(right)].push_back(
            {left, relation});
    }

    std::map<RawColor, long long> initial_counts;
    std::vector<std::string> raw_representations;
    raw_representations.reserve(raw_colors.size());
    for (const auto& raw_color : raw_colors) {
        ++initial_counts[raw_color];
        raw_representations.push_back(pair_repr(raw_color));
    }
    std::vector<int> colors = compact(raw_representations);

    for (Py_ssize_t iteration = 0; iteration < size; ++iteration) {
        std::map<int, int> class_counts;
        for (const int color : colors) {
            ++class_counts[color];
        }
        std::vector<std::string> signatures;
        signatures.reserve(static_cast<std::size_t>(size));
        for (Py_ssize_t left = 0; left < size; ++left) {
            std::map<NeighborhoodKey, int> neighborhood;
            for (const auto& [color, count] : class_counts) {
                neighborhood[{zero_relation, color}] = count;
            }
            --neighborhood[{zero_relation, colors[static_cast<std::size_t>(left)]}];
            for (const auto& [right, relation] :
                    active_by_position[static_cast<std::size_t>(left)]) {
                const int neighbor_color = colors[static_cast<std::size_t>(right)];
                --neighborhood[{zero_relation, neighbor_color}];
                ++neighborhood[{relation, neighbor_color}];
            }
            std::vector<std::string> entries;
            entries.reserve(neighborhood.size());
            for (const auto& [key, count] : neighborhood) {
                if (count != 0) {
                    entries.push_back(neighborhood_entry_repr(key, count));
                }
            }
            signatures.push_back(
                "(" + std::to_string(colors[static_cast<std::size_t>(left)]) +
                ", " + tuple_repr(entries) + ")");
        }
        std::vector<int> refined = compact(signatures);
        if (refined == colors) {
            break;
        }
        colors = std::move(refined);
    }

    std::map<int, long long> final_color_counts;
    for (const int color : colors) {
        ++final_color_counts[color];
    }
    std::map<RelationCountKey, long long> relation_counts;
    for (const auto& [left_color, left_count] : final_color_counts) {
        relation_counts[{left_color, left_color, zero_relation}] +=
            left_count * (left_count - 1) / 2;
        for (const auto& [right_color, right_count] : final_color_counts) {
            if (left_color < right_color) {
                relation_counts[{left_color, right_color, zero_relation}] +=
                    left_count * right_count;
            }
        }
    }
    for (const auto& [endpoints, relation] : relations) {
        const auto [left, right] = endpoints;
        const int left_color = std::min(colors[left], colors[right]);
        const int right_color = std::max(colors[left], colors[right]);
        --relation_counts[{left_color, right_color, zero_relation}];
        ++relation_counts[{left_color, right_color, relation}];
    }

    PyObject* initial_result = PyTuple_New(
        static_cast<Py_ssize_t>(initial_counts.size()));
    if (initial_result == nullptr) {
        return nullptr;
    }
    Py_ssize_t output_index = 0;
    for (const auto& [raw_color, count] : initial_counts) {
        PyObject* color = integer_pair(raw_color.first, raw_color.second);
        PyObject* item = color == nullptr
            ? nullptr : Py_BuildValue("(NL)", color, count);
        if (item == nullptr) {
            Py_DECREF(initial_result);
            return nullptr;
        }
        PyTuple_SET_ITEM(initial_result, output_index++, item);
    }
    PyObject* color_result = color_counts_object(final_color_counts);
    if (color_result == nullptr) {
        Py_DECREF(initial_result);
        return nullptr;
    }
    std::size_t nonzero_relation_count = 0;
    for (const auto& [key, count] : relation_counts) {
        if (count != 0) {
            ++nonzero_relation_count;
        }
    }
    PyObject* relation_result = PyTuple_New(
        static_cast<Py_ssize_t>(nonzero_relation_count));
    if (relation_result == nullptr) {
        Py_DECREF(initial_result);
        Py_DECREF(color_result);
        return nullptr;
    }
    output_index = 0;
    for (const auto& [key, count] : relation_counts) {
        if (count == 0) {
            continue;
        }
        const auto& [left_color, right_color, relation] = key;
        PyObject* relation_object = integer_pair(relation.first, relation.second);
        PyObject* relation_key = relation_object == nullptr
            ? nullptr
            : Py_BuildValue("(iiN)", left_color, right_color, relation_object);
        PyObject* item = relation_key == nullptr
            ? nullptr : Py_BuildValue("(NL)", relation_key, count);
        if (item == nullptr) {
            Py_DECREF(initial_result);
            Py_DECREF(color_result);
            Py_DECREF(relation_result);
            return nullptr;
        }
        PyTuple_SET_ITEM(relation_result, output_index++, item);
    }
    return Py_BuildValue("(NNN)", initial_result, color_result, relation_result);
}

PyMethodDef methods[] = {
    {"paired_mapping_invariant", paired_mapping_invariant, METH_VARARGS,
     "Compute the exact sparse paired-mapping refinement."},
    {nullptr, nullptr, 0, nullptr},
};

PyModuleDef module = {
    PyModuleDef_HEAD_INIT,
    "_native",
    "Native exact kernels for rxn_core.",
    -1,
    methods,
};

}  // namespace

PyMODINIT_FUNC PyInit__native() {
    return PyModule_Create(&module);
}
