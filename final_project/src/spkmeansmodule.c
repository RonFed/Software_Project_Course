#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "spkmeans.h"

static PyObject *convert_mat_to_PyObject(matrix *mat)
{
    int i, j;
    unsigned int rows, cols;
    rows = mat->rows;
    cols = mat->cols;
    PyObject *obj = PyList_New(rows);
    for (i = 0; i < rows; i++)
    {
        PyObject *currnet = PyList_New(cols);
        for (j = 0; j < cols; j++)
        {
            PyObject *val = Py_BuildValue("d", (mat->data)[i][j]);
            PyList_SetItem(currnet, j, val);
        }
        PyList_SetItem(obj, i, currnet);
    }
    return obj;
}

static PyObject *convert_sym_mat_to_PyObject(sym_matrix *mat)
{
    int i, j;
    unsigned int dim;
    dim = mat->dim;
    PyObject *obj = PyList_New(dim);
    for (i = 0; i < dim; i++)
    {
        PyObject *currnet = PyList_New(dim);
        for (j = 0; j < dim; j++)
        {
            PyObject *val = Py_BuildValue("d", get_val_sym(mat, i, j));
            PyList_SetItem(currnet, j, val);
        }
        PyList_SetItem(obj, i, currnet);
    }
    return obj;
}

static matrix *convert_PyObject_to_mat(PyObject *mat_py_obj)
{
    int i, j;
    unsigned int rows, cols;
    matrix *result;
    rows = PyObject_Length(mat_py_obj);
    cols = PyList_Size(PyList_GET_ITEM(mat_py_obj, 0));
    result = init_mat(rows, cols);
    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < cols; j++)
        {
            (result->data)[i][j] = PyFloat_AS_DOUBLE(PyList_GET_ITEM(PyList_GET_ITEM(mat_py_obj, i), j));
        }
    }
    result->rows = rows;
    result->cols = cols;
    return result;
}

static PyObject *weights_mat_c_api(PyObject *self, PyObject *args)
{
    PyObject *data_arr_list;
    matrix *data_mat;
    sym_matrix *weights_mat_sym;

    if (!PyArg_ParseTuple(args, "O", &data_arr_list))
    {
        return NULL;
    }

    data_mat = convert_PyObject_to_mat(data_arr_list);
    weights_mat_sym = weights_mat(data_mat);

    return convert_sym_mat_to_PyObject(weights_mat_sym);
}

static PyMethodDef capiMethods[] = {
    {"weights_mat",
     (PyCFunction)weights_mat_c_api,
     METH_VARARGS,
     PyDoc_STR("Calculate the weights matrix of the given matrix")},
    {NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "spkmeans",
    NULL,
    -1,
    capiMethods};

PyMODINIT_FUNC
PyInit_spkmeans(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m)
    {
        return NULL;
    }
    return m;
}