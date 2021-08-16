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

static PyObject *convert_diag_mat_to_PyObject(diag_matrix *mat)
{
    int i, j;
    unsigned int dim;
    PyObject *val;
    dim = mat->dim;
    PyObject *obj = PyList_New(dim);
    for (i = 0; i < dim; i++)
    {
        PyObject *currnet = PyList_New(dim);
        for (j = 0; j < dim; j++)
        {
            if (i == j)
            {
                val = Py_BuildValue("d", (mat->data)[i]);
            } else {
                val = Py_BuildValue("d", 0.0);
            }
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
    PyObject *weights_mat_py;
    matrix *data_mat;
    sym_matrix *weights_mat_sym;

    if (!PyArg_ParseTuple(args, "O", &data_arr_list))
    {
        return NULL;
    }

    data_mat = convert_PyObject_to_mat(data_arr_list);
    weights_mat_sym = weights_mat(data_mat);
    weights_mat_py = convert_sym_mat_to_PyObject(weights_mat_sym);
    free_sym_mat(weights_mat_sym);
    return weights_mat_py;
}

static PyObject *degree_mat_c_api(PyObject *self, PyObject *args)
{
    PyObject *data_arr_list;
    PyObject *diagonal_mat_py;
    matrix *data_mat;
    diag_matrix *diagonal_mat;

    if (!PyArg_ParseTuple(args, "O", &data_arr_list))
    {
        return NULL;
    }

    data_mat = convert_PyObject_to_mat(data_arr_list);
    diagonal_mat = degree_mat(data_mat);
    diagonal_mat_py = convert_diag_mat_to_PyObject(diagonal_mat);
    free_diag_mat(diagonal_mat);
    return diagonal_mat_py;
}

static PyObject *l_norm_mat_c_api(PyObject *self, PyObject *args)
{
    PyObject *data_arr_list;
    PyObject *lnorm_mat_py;
    matrix *data_mat;
    sym_matrix *lnorm_m;

    if (!PyArg_ParseTuple(args, "O", &data_arr_list))
    {
        return NULL;
    }

    data_mat = convert_PyObject_to_mat(data_arr_list);
    lnorm_m = l_norm_mat(data_mat);
    lnorm_mat_py = convert_sym_mat_to_PyObject(lnorm_m);
    free_sym_mat(lnorm_m);
    return lnorm_mat_py;
}

/*
Functions exported to the extrenal API
Can be called by Python
*/
static PyMethodDef capiMethods[] = {
    {"weights_mat",
     (PyCFunction)weights_mat_c_api,
     METH_VARARGS,
     PyDoc_STR("Calculate the weights matrix of a given matrix")
     },
     {
    "degree_mat",
     (PyCFunction)degree_mat_c_api,
     METH_VARARGS,
     PyDoc_STR("Calculate the degree matrix of a given matrix")
     },
     {
    "l_norm_mat",
     (PyCFunction)l_norm_mat_c_api,
     METH_VARARGS,
     PyDoc_STR("Calculate the normalized laplacian matrix of a given matrix")
     },
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