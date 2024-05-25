#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <ndarraytypes.h>

/*
Provides an implementation of count_collisions().

Expects arguments:
  xs (np.ndarray, dtype=float64, shape=(N,3)): positions of points
  vs (np.ndarray, dtype=float64, shape=(N,3)): velocities of points
  rs (np.ndarray, dtype=float64, shape=(N)): radii of points

Returns:
  count: int64

May raise an exception if internal memory allocation fails.
*/
static PyObject* pcoll_count(PyObject *self, PyObject *args) {
    // raw arg
    PyObject *arg1=NULL, *arg2=NULL, *arg3=NULL;
    // numpy array objects
    PyObject *arr_xs=NULL, *arr_vs=NULL, *arr_rs=NULL;
    // dimension counts
    int ndims_xs, ndims_vs, ndims_rs;
    // dimensions
    npy_intp *dims_xs, *dims_vs, *dims_rs;
    // number of points
    size_t n_pts;
    // pointers to actual data
    float *xs, *vs, *rs;
    // collision count
    int count;

    // get raw python objects
    if (!PyArg_ParseTuple(args, "OOO", &arg1, &arg2, &arg3)) {
        return NULL;
    }

    // conversion
    arr_xs = PyArray_FROM_OTF(arg1, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (arr_xs == NULL) return NULL;
    arr_vs = PyArray_FROM_OTF(arg2, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (arr_vs == NULL) goto fail;
    arr_rs = PyArray_FROM_OTF(arg3, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (arr_rs == NULL) goto fail;

    // get count of dimensions for all arrays
    ndims_xs = PyArray_NDIM(arr_xs);
    ndims_vs = PyArray_NDIM(arr_vs);
    ndims_rs = PyArray_NDIM(arr_rs);

    // check array sizes
    if (ndims_xs != 2) goto fail;
    if (ndims_vs != 2) goto fail;
    if (ndims_rs != 1) goto fail;

    // get dimensions of all arrays
    dims_xs = PyArray_DIMS(arr_xs);
    dims_vs = PyArray_DIMS(arr_vs);
    dims_rs = PyArray_DIMS(arr_rs);

    // assert shapes are correct
    if ((dims_rs[0] != dims_vs[0])) || (dims_rs[0] != dims_xs[0])) {
        // array shape mismatch
        // TODO: raise an exception
        goto fail;
    }
    if ((dims_vs[1] != 3) || (dims_xs[1] != 3)) {
        // invlaid array shape for vs or xs
        // TODO: raise an exception
        goto fail;
    }

    // get number of points
    n_pts = dims_rs[0];

    // get data
    xs = (double*) PyArray_DATA(arr_xs);
    vs = (double*) PyArray_DATA(arr_vs);
    rs = (double*) PyArray_DATA(arr_rs);

    // perform the actual operation
    count = count_collisions_fast(n_pts, xs, vs, rs);

    return PyLong_FromLong(count);

fail:
    Py_XDECREF(arr_xs);
    Py_XDECREF(arr_vs);
    Py_XDECREF(arr_rs);
    return NULL;
}

static PyMethodDef PointCollisionMethods[] = {
    ...
    {"count_collisions",  pcoll_count, METH_VARARGS,
     "Execute a shell command."},
    ...
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

static struct PyModuleDef pointcollision = {
    PyModuleDef_HEAD_INIT,
    "pointcollision",   /* name of module */
    NULL,     /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    PointCollisionMethods
};

PyMODINIT_FUNC PyInit_pointcollision(void)
{
    return PyModule_Create(&pointcollision);
}
