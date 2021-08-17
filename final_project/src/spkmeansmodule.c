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

/*The returned value is a 1D list containing the diagonal elements */
static PyObject *convert_diag_mat_to_PyObject(diag_matrix *mat)
{
    int i, j;
    unsigned int dim;
    dim = mat->dim;
    PyObject *vector = PyList_New(dim);
    for (i = 0; i < dim; i++)
    {
        PyObject *val = Py_BuildValue("d", (mat->data)[i]);
        PyList_SetItem(vector, i, val);
    }
    return vector;
}

/* The returned representation of the Jacobi matrix is a Python list :
    first element is a list of eigan values
    second element is a matrix (list of lists) containing the eiganvectors */
static PyObject *convert_jacobi_mat_to_PyObject(jacobi_matrix *j_mat)
{
    int i, j;
    unsigned int dim;
    dim = j_mat->mat->dim;
    PyObject *eigan_values_arr_py;
    PyObject* eigan_mat_py;
    PyObject *result = PyList_New(2);
    
    /* building the 1D array with eiganvalues*/
    eigan_values_arr_py = PyList_New(dim);
    for (i = 0; i < dim; i++)
    {
        PyObject *val = Py_BuildValue("d", (j_mat->e_mat)[i].e_val);
        PyList_SetItem(eigan_values_arr_py, i, val);
    }
    PyList_SetItem(result, 0, eigan_values_arr_py);

     /* building the 2D array with eiganvectors*/
    eigan_mat_py = PyList_New(dim);
    for (i = 0; i < dim; i++)
    {
        PyObject *currnet_row = PyList_New(dim);
        for (j = 0; j < dim; j++)
        {
            PyObject *val = Py_BuildValue("d", ((j_mat->e_mat)[i].vec)[j]);
            PyList_SetItem(currnet_row, j, val);
        }
        PyList_SetItem(eigan_mat_py, i, currnet_row);
    }
    PyList_SetItem(result, 1, eigan_mat_py);

    return result;
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

static PyObject *jacobi_mat_c_api(PyObject *self, PyObject *args)
{
    PyObject *data_arr_list;
    PyObject *jacobi_mat_data_py;
    matrix *data_mat;
    sym_matrix *data_sym_mat;
    jacobi_matrix *jacobi_m;

    if (!PyArg_ParseTuple(args, "O", &data_arr_list))
    {
        return NULL;
    }

    data_mat = convert_PyObject_to_mat(data_arr_list);
    data_sym_mat = matrix_to_sym_matrix(data_mat);
    jacobi_m = init_jac_mat(data_sym_mat);
    jacobi(jacobi_m);
    jacobi_mat_data_py = convert_jacobi_mat_to_PyObject(jacobi_m);
    return jacobi_mat_data_py;
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
     PyDoc_STR("Calculate the degree matrix of a given matrix (return a 1D array of diagonal elements values")
     },
     {
    "l_norm_mat",
     (PyCFunction)l_norm_mat_c_api,
     METH_VARARGS,
     PyDoc_STR("Calculate the normalized laplacian matrix of a given matrix")
     },
    {
    "jacobi_mat",
     (PyCFunction)jacobi_mat_c_api,
     METH_VARARGS,
     PyDoc_STR("Calculate the eiganvalues and eiganvectoes using the jacobi algorithm")
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