// Native port of rxn_core.alignment.sweep._freeze_analytical_py.
//
// The semantics must match the pure-Python original exactly:
//
//   dict            -> tuple(sorted((str(key), freeze(item))
//                                   for key, item in value.items()))
//   list / tuple    -> tuple(freeze(item) for item in value)
//   set / frozenset -> tuple(sorted((freeze(item) for item in value), key=repr))
//   anything else   -> the very same object (no copy, no conversion)
//
// isinstance() semantics are kept: PyDict_Check / PyList_Check /
// PyTuple_Check / PyAnySet_Check all accept subclasses.  Exact dicts, lists
// and tuples take direct C-API paths; subclasses go through the same protocol
// Python uses (dict.items() / __iter__), because e.g. an OrderedDict that was
// reordered with move_to_end() is *not* in PyDict_Next order, and the stable
// sort keeps the input order between equal pairs.  Sorting itself is done by
// Python's own list.sort / builtins.sorted so comparison semantics (and the
// TypeError raised for unorderable pairs) are identical.
#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace {

py::object freeze_analytical(py::handle value);

inline PyObject* check(PyObject* p) {
    if (p == nullptr) throw py::error_already_set();
    return p;
}

inline py::object steal(PyObject* p) {
    return py::reinterpret_steal<py::object>(check(p));
}

// Mirror Python's recursion limit so a pathologically deep input raises
// RecursionError (as the Python version does) instead of overflowing the
// C stack.  Py_EnterRecursiveCall undoes its increment when it fails, so a
// throwing constructor (whose destructor never runs) is correct.
struct RecursionGuard {
    RecursionGuard() {
        if (Py_EnterRecursiveCall(" in freeze_analytical")) throw py::error_already_set();
    }
    ~RecursionGuard() { Py_LeaveRecursiveCall(); }
    RecursionGuard(const RecursionGuard&) = delete;
    RecursionGuard& operator=(const RecursionGuard&) = delete;
};

inline PyObject* make_pair(py::object first, py::object second) {
    PyObject* pair = check(PyTuple_New(2));
    PyTuple_SET_ITEM(pair, 0, first.release().ptr());
    PyTuple_SET_ITEM(pair, 1, second.release().ptr());
    return pair;
}

inline void append_steal(PyObject* list, PyObject* item) {
    int rc = PyList_Append(list, item);
    Py_DECREF(item);
    if (rc < 0) throw py::error_already_set();
}

// list.sort() then tuple(): exactly what tuple(sorted(<list>)) does.
py::object sort_to_tuple(py::object list) {
    if (PyList_Sort(list.ptr()) < 0) throw py::error_already_set();
    return steal(PyList_AsTuple(list.ptr()));
}

// [freeze(item) for item in iterable] through the iterator protocol; this is
// the path Python's own `for item in value` takes (including e.g. the
// "Set changed size during iteration" check of the set iterator).
py::object freeze_iterable_to_list(py::handle obj) {
    py::object out = steal(PyList_New(0));
    py::object iter = steal(PyObject_GetIter(obj.ptr()));
    while (PyObject* raw = PyIter_Next(iter.ptr())) {
        py::object item = py::reinterpret_steal<py::object>(raw);
        py::object frozen = freeze_analytical(item);
        if (PyList_Append(out.ptr(), frozen.ptr()) < 0) throw py::error_already_set();
    }
    if (PyErr_Occurred()) throw py::error_already_set();
    return out;
}

// --- dict -------------------------------------------------------------------

// Exact dict: PyDict_Next is items() order.  Every callback into Python
// (str(key), freeze(item)) may run arbitrary __str__/__repr__ code, so both
// the key and the value are owned before calling out, and a size change is
// reported the way Python's dict iterator reports it.
py::object freeze_dict_exact(PyObject* obj) {
    const Py_ssize_t n = PyDict_Size(obj);
    py::object out = steal(PyList_New(n));
    Py_ssize_t pos = 0, i = 0;
    PyObject *key, *item;
    while (PyDict_Next(obj, &pos, &key, &item)) {
        py::object key_ref = py::reinterpret_borrow<py::object>(key);
        py::object item_ref = py::reinterpret_borrow<py::object>(item);
        py::object skey = steal(PyObject_Str(key_ref.ptr()));  // str(key) first, as in the original
        py::object frozen = freeze_analytical(item_ref);
        if (PyDict_Size(obj) != n) {
            PyErr_SetString(PyExc_RuntimeError, "dictionary changed size during iteration");
            throw py::error_already_set();
        }
        PyObject* pair = make_pair(std::move(skey), std::move(frozen));
        if (i < n) {
            PyList_SET_ITEM(out.ptr(), i, pair);
        } else {
            // A same-size mutation (del + insert) can make PyDict_Next yield
            // more entries than the dict has; keep everything it yields.
            append_steal(out.ptr(), pair);
        }
        ++i;
    }
    if (i < n) {
        // ... or fewer: drop the unused (NULL) slots before sorting.
        if (PyList_SetSlice(out.ptr(), i, n, nullptr) < 0) throw py::error_already_set();
    }
    return sort_to_tuple(std::move(out));
}

// dict subclass: `for key, item in value.items()` verbatim.
py::object freeze_dict_items(py::handle obj) {
    py::object out = steal(PyList_New(0));
    py::object items = obj.attr("items")();
    py::object iter = steal(PyObject_GetIter(items.ptr()));
    while (PyObject* raw = PyIter_Next(iter.ptr())) {
        py::object entry = py::reinterpret_steal<py::object>(raw);
        // `key, item = entry`: any iterable of exactly two elements.
        py::object seq = steal(PySequence_Tuple(entry.ptr()));
        const Py_ssize_t size = PyTuple_GET_SIZE(seq.ptr());
        if (size < 2) {
            PyErr_Format(PyExc_ValueError,
                         "not enough values to unpack (expected 2, got %zd)", size);
            throw py::error_already_set();
        }
        if (size > 2) {
            PyErr_SetString(PyExc_ValueError, "too many values to unpack (expected 2)");
            throw py::error_already_set();
        }
        py::object skey = steal(PyObject_Str(PyTuple_GET_ITEM(seq.ptr(), 0)));
        py::object frozen = freeze_analytical(PyTuple_GET_ITEM(seq.ptr(), 1));
        append_steal(out.ptr(), make_pair(std::move(skey), std::move(frozen)));
    }
    if (PyErr_Occurred()) throw py::error_already_set();
    return sort_to_tuple(std::move(out));
}

// --- list / tuple -----------------------------------------------------------

// Exact list.  Python's list iterator re-reads the length on every step and
// keeps going if the list grew (or stops early if it shrank) under it, so the
// result tuple is resized in those cases - the same thing PySequence_Tuple
// does for an iterator of unknown length.
py::object freeze_list_exact(PyObject* obj) {
    py::object out = steal(PyTuple_New(PyList_GET_SIZE(obj)));
    Py_ssize_t i = 0;
    for (; i < PyList_GET_SIZE(obj); ++i) {
        py::object item = py::reinterpret_borrow<py::object>(PyList_GET_ITEM(obj, i));
        py::object frozen = freeze_analytical(item);
        if (i >= PyTuple_GET_SIZE(out.ptr())) {
            PyObject* raw = out.release().ptr();
            if (_PyTuple_Resize(&raw, i + 1) < 0) throw py::error_already_set();  // raw freed, NULL
            out = py::reinterpret_steal<py::object>(raw);
        }
        PyTuple_SET_ITEM(out.ptr(), i, frozen.release().ptr());
    }
    if (i < PyTuple_GET_SIZE(out.ptr())) {
        PyObject* raw = out.release().ptr();
        if (_PyTuple_Resize(&raw, i) < 0) throw py::error_already_set();
        out = py::reinterpret_steal<py::object>(raw);
    }
    return out;
}

// Exact tuple: immutable, so the size is fixed and items cannot go away.
py::object freeze_tuple_exact(PyObject* obj) {
    const Py_ssize_t n = PyTuple_GET_SIZE(obj);
    py::object out = steal(PyTuple_New(n));
    for (Py_ssize_t i = 0; i < n; ++i) {
        py::object frozen = freeze_analytical(PyTuple_GET_ITEM(obj, i));
        PyTuple_SET_ITEM(out.ptr(), i, frozen.release().ptr());
    }
    return out;
}

// --- set / frozenset --------------------------------------------------------

py::object freeze_set(py::handle obj) {
    py::object frozen = freeze_iterable_to_list(obj);
    // sorted(<list>, key=repr): Python computes every key first, then does a
    // stable sort with default str comparison.  Calling the builtin keeps
    // that (and its error behaviour) identical by construction.
    py::module_ builtins = py::module_::import("builtins");
    py::object sorted_list =
        builtins.attr("sorted")(frozen, py::arg("key") = builtins.attr("repr"));
    return steal(PyList_AsTuple(sorted_list.ptr()));
}

// --- entry point ------------------------------------------------------------

py::object freeze_analytical(py::handle value) {
    PyObject* obj = value.ptr();
    if (PyDict_Check(obj)) {
        RecursionGuard guard;
        return PyDict_CheckExact(obj) ? freeze_dict_exact(obj) : freeze_dict_items(value);
    }
    if (PyList_Check(obj)) {
        RecursionGuard guard;
        if (PyList_CheckExact(obj)) return freeze_list_exact(obj);
        return steal(PyList_AsTuple(freeze_iterable_to_list(value).ptr()));
    }
    if (PyTuple_Check(obj)) {
        RecursionGuard guard;
        if (PyTuple_CheckExact(obj)) return freeze_tuple_exact(obj);
        return steal(PyList_AsTuple(freeze_iterable_to_list(value).ptr()));
    }
    if (PyAnySet_Check(obj)) {
        RecursionGuard guard;
        return freeze_set(value);
    }
    return py::reinterpret_borrow<py::object>(value);  // leaf: the same object
}

}  // namespace

void register_freeze(py::module_& m) {
    m.def("freeze_analytical", &freeze_analytical, py::arg("value"),
          "Exact native port of rxn_core.alignment.sweep._freeze_analytical_py.\n\n"
          "dict -> tuple(sorted((str(k), freeze(v)) for k, v in value.items()));\n"
          "list/tuple -> tuple of frozen items; set/frozenset -> tuple(sorted(..., key=repr));\n"
          "anything else is returned unchanged (same object).");
}
