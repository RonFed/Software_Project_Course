#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#define BUFFERSIZE 1000
#define BASE_ARR_SIZE 100
#define ARR_SIZE_MULTIPLY 2
#define COMMA ','
#define MALLOC_FOR_STRING(str) (char *)(malloc((strlen(str) + 1) * sizeof(char)))
#define MALLOC_ARR(type, size) (type *)malloc((size) * (sizeof(type)))
#define ASSERT_MALLOC(ptr) assert((ptr) != NULL)
#define MY_MALLOC(ptr, type, size)      \
    do                                  \
    {                                   \
        (ptr) = MALLOC_ARR(type, size); \
        ASSERT_MALLOC(ptr);             \
    } while (0)

int d = 0;
int k = 0;
double minus_one = -1;
int max_iter = 200;
int lines_count;
double distacne(double **a, double **b);
int closest_cluster_index(double *vector, double **centers);
void multiply_scalar(double *vector, double *scalar, double *result);
void add_vectors(double *a, double *b, double *result);
void print_final(double **mat);
void malloc_for_mat_and_set_zeros(double ***mat, int lines, int cols);
void k_means(double ***centers, double ***vectors);
void convert_PyObject_to_mat(PyObject *obj, double ***mat, int lines, int cols);
static PyObject *convert_mat_to_PyObject(double ***mat, int lines, int cols);

int main(int argc, char const *argv[])
{
    return 0;
}

// ********************* Python stuff ***********************

static PyObject *fit(PyObject *self, PyObject *args)
{
    PyObject *data_list;
    PyObject *initial_centers;
    double **vectors;
    double **centeroids;

    if (!PyArg_ParseTuple(args, "OOi", &initial_centers, &data_list, &max_iter))
    {
        return NULL;
    }
    k = PyObject_Length(initial_centers);
    lines_count = PyObject_Length(data_list);
    d = PyList_Size(PyList_GET_ITEM(initial_centers, 0));

    malloc_for_mat_and_set_zeros(&vectors, lines_count, d);
    malloc_for_mat_and_set_zeros(&centeroids, k, d);

    convert_PyObject_to_mat(data_list, &vectors, lines_count, d);
    convert_PyObject_to_mat(initial_centers, &centeroids, k, d);

    k_means(&centeroids, &vectors);

    return convert_mat_to_PyObject(&centeroids, k, d);
}

void convert_PyObject_to_mat(PyObject *obj, double ***mat, int lines, int cols)
{
    int i, j;
    for (i = 0; i < lines; i++)
    {
        for (j = 0; j < cols; j++)
        {
            (*mat)[i][j] = PyFloat_AS_DOUBLE(PyList_GET_ITEM(PyList_GET_ITEM(obj, i), j));
        }
    }
}

static PyObject *convert_mat_to_PyObject(double ***mat, int lines, int cols)
{
    int i, j;
    PyObject *obj = PyList_New(lines);
    for (i = 0; i < lines; i++)
    {
        PyObject *currnet = PyList_New(cols);
        for (j = 0; j < cols; j++)
        {
            PyObject *val = Py_BuildValue("d", (*mat)[i][j]);
            PyList_SetItem(currnet, j, val);
        }
        PyList_SetItem(obj, i, currnet);
    }
    return obj;
}

static PyMethodDef capiMethods[] = {
    {"fit",
     (PyCFunction)fit,
     METH_VARARGS,
     PyDoc_STR("K means algorithm with initial values provided")},
    {NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "mykmeanssp",
    NULL,
    -1,
    capiMethods};

PyMODINIT_FUNC
PyInit_mykmeanssp(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m)
    {
        return NULL;
    }
    return m;
}

// ********************* end of Python stuff ***********************


double distacne(double **a, double **b)
{
    double sum = 0;
    int i;
    for (i = 0; i < d; i++)
    {
        sum += ((*a)[i] - (*b)[i]) * ((*a)[i] - (*b)[i]);
    }
    return sum;
}

int closest_cluster_index(double *vector, double **centers)
{
    int min_dist_index = 0, i;
    double min_dist;
    min_dist = distacne(&vector, &centers[0]);
    for (i = 0; i < k; i++)
    {
        double current_dist;
        current_dist = distacne(&vector, &centers[i]);
        if (current_dist < min_dist)
        {
            min_dist = current_dist;
            min_dist_index = i;
        }
        if (min_dist == 0)
        {
            break;
        }
    }
    return min_dist_index;
}

void multiply_scalar(double *vector, double *scalar, double *result)
{
    int i;
    for (i = 0; i < d; i++)
    {
        result[i] = vector[i] * (*scalar);
    }
}

void add_vectors(double *a, double *b, double *result)
{
    int i;
    for (i = 0; i < d; i++)
    {
        result[i] = a[i] + b[i];
    }
}

void print_final(double **mat)
{
    int i, j;
    for (i = 0; i < k; i++)
    {
        for (j = 0; j < d - 1; j++)
        {
            printf("%.4f,", mat[i][j]);
        }
        printf("%.4f", mat[i][d - 1]);
        printf("\n");
    }
}

void malloc_for_mat_and_set_zeros(double ***mat, int lines, int cols)
{
    int i;
    MY_MALLOC(*mat, double *, lines);
    for (i = 0; i < lines; i++)
    {
        (*mat)[i] = (double *)calloc(cols, sizeof(double));
        ASSERT_MALLOC((*mat)[i]);
    }
}

void k_means(double ***centers, double ***vectors)
{
    int i, j, l;
    double **cluster_sums;
    double *negative_vec;
    double inverse_size;
    int *cluster_size;
    int *which_cluster;
    int new_closest_cluster;
    MY_MALLOC(negative_vec, double, d);
    MY_MALLOC(cluster_size, int, k);
    MY_MALLOC(which_cluster, int, lines_count);
    malloc_for_mat_and_set_zeros(&cluster_sums, k, d);
    for (i = 0; i < k; i++)
    {
        cluster_size[i] = 0;
    }
    for (i = 0; i < lines_count; i++)
    {
        which_cluster[i] = -1;
    }
    for (i = 0; i < max_iter; i++)
    {
        int differ = 0;
        for (j = 0; j < lines_count; j++)
        {
            new_closest_cluster = closest_cluster_index((*vectors)[j], *centers);
            if (new_closest_cluster != which_cluster[j])
            {
                if (which_cluster[j] != -1)
                {
                    cluster_size[which_cluster[j]] -= 1;
                    multiply_scalar((*vectors)[j], &minus_one, negative_vec);
                    add_vectors(cluster_sums[which_cluster[j]], negative_vec, cluster_sums[which_cluster[j]]);
                }
                cluster_size[new_closest_cluster] += 1;
                add_vectors(cluster_sums[new_closest_cluster], (*vectors)[j], cluster_sums[new_closest_cluster]);
                which_cluster[j] = new_closest_cluster;
                differ = 1;
            }
        }
        if (differ == 0)
        {
            break;
        }
        for (l = 0; l < k; l++)
        {
            if (cluster_size[l] != 0)
            {
                inverse_size = (double)1 / cluster_size[l];
                multiply_scalar(cluster_sums[l], &inverse_size, (*centers)[l]);
            }
        }
    }
}