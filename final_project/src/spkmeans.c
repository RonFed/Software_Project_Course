#include "spkmeans.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

/*
GENERAL MATRIX FUNCTIONS
*/
matrix *init_mat(unsigned int rows, unsigned int cols)
{
    int i;
    matrix *mat = (matrix *)malloc(sizeof(matrix));
    mat->rows = rows;
    mat->cols = cols;
    MALLOC_ARR_ASSERT(mat->data, double *, rows);
    for (i = 0; i < rows; i++)
    {
        CALLOC_ARR_ASSERT((mat->data)[i], double, cols);
    }
    return mat;
}

void print_mat(matrix *mat)
{
    int i, j;
    for (i = 0; i < mat->rows; i++)
    {
        for (j = 0; j < mat->cols; j++)
        {
            printf("%.4f,", (mat->data)[i][j]);
        }
        printf("\n");
    }
}

sym_matrix *init_sym_mat(unsigned int dim)
{
    int i;
    sym_matrix *mat = (sym_matrix *)malloc(sizeof(sym_matrix));
    mat->dim = dim;
    MALLOC_ARR_ASSERT(mat->data, double *, dim);
    for (i = 0; i < dim; i++)
    {
        CALLOC_ARR_ASSERT((mat->data)[i], double, i + 1);
    }
    return mat;
}

void print_sym_mat(sym_matrix *mat)
{
    int i, j;
    unsigned int dim = mat->dim;
    for (i = 0; i < dim; i++)
    {
        for (j = 0; j < i + 1; j++)
        {
            printf("%.4f,", (mat->data)[i][j]);
        }
        for (j = i + 1; j < dim; j++)
        {
            printf("%.4f,", (mat->data)[j][i]);
        }
        printf("\n");
    }
}

diag_matrix *init_diag_mat(unsigned int dim)
{
    diag_matrix *mat = (diag_matrix *)malloc(sizeof(diag_matrix));
    mat->dim = dim;
    MALLOC_ARR_ASSERT(mat->data, double, dim);
    return mat;
}

void print_diag_mat(diag_matrix *mat)
{
    int i, j;
    unsigned int dim = mat->dim;
    for (i = 0; i < dim; i++)
    {
        for (j = 0; j < dim; j++)
        {
            if (j == i)
            {
                printf("%.4f", (mat->data)[i]);
            }
            else
            {
                printf(ZERO_4F);
            }
            if (j != dim - 1)
            {
                printf(",");
            }
        }
        printf("\n");
    }
}

/*
WEIGHTS MATRIX FUNCTIONS
*/

double l2_norm_vectors(matrix *mat, unsigned int row1, unsigned int row2)
{
    assert((row1 < mat->rows) && (row2 < mat->rows));
    double result = 0, dif;
    int i;
    for (i = 0; i < mat->cols; i++)
    {
        dif = ((mat->data)[row1][i]) - ((mat->data)[row2][i]);
        result += (pow(dif, 2));
    }
    result = pow(result, 0.5);
    return result;
}

double exp_norm_vectors(matrix *mat, unsigned int row1, unsigned int row2)
{
    double l2_norm = l2_norm_vectors(mat, row1, row2);
    double result;
    result = exp(-l2_norm / 2);
    return result;
}

sym_matrix *weights_mat(matrix *mat)
{
    sym_matrix *sm = init_sym_mat(mat->rows);
    int i, j;
    for (i = 0; i < mat->rows; i++)
    {
        for (j = 0; j < i; j++)
        {
            /* each cell in the weights matrix is the result
            of exp norm between two vectrose */
            (sm->data)[i][j] = exp_norm_vectors(mat, i, j);
        }
    }
    return sm;
}

/*
DEGREE MATRIX FUNCTIONS
*/

double row_sum_sym_mat(sym_matrix *mat, unsigned int row)
{
    assert(row < mat->dim);
    double sum = 0;
    int j;
    for (j = 0; j < row + 1; j++)
    {
        sum += (mat->data)[row][j];
    }
    for (j = row + 1; j < mat->dim; j++)
    {
        sum += (mat->data)[j][row];
    }
    return sum;
}

diag_matrix *degree_mat(matrix *mat)
{
    int i;
    // allocating memory for degree mat as diagonal matrix
    diag_matrix *degree_m = init_diag_mat(mat->rows);
    sym_matrix *weight_m = weights_mat(mat);
    for (i = 0; i < mat->rows; i++)
    {
        // sum of a row in the weights matrix
        (degree_m->data)[i] = row_sum_sym_mat(weight_m, i);
    }
    return degree_m;
}

sym_matrix *l_norm_mat(matrix *mat)
{
    sym_matrix *l_norm_m = init_sym_mat(mat->rows);
    sym_matrix *weight_m = weights_mat(mat);
    diag_matrix *degree_m = degree_mat(mat);
    double curr_weight, curr_d_i, curr_d_j, d_sqrt;
    int i, j;
    for (i = 0; i < l_norm_m->dim; i++)
    {
        curr_d_i = degree_m->data[i];
        for (j = 0; j <= i; j++)
        {
            curr_d_j = degree_m->data[j];
            curr_weight = weight_m->data[i][j];
            if (j == i)
            {
                // on the diagonal w(i,i) = 0 => l_norm(i,i) = 1
                (l_norm_m->data)[i][j] = 1;
            }
            else
            {
                d_sqrt = pow(curr_d_i * curr_d_j, -0.5);
                (l_norm_m->data)[i][j] = -d_sqrt * curr_weight;
            }
        }
    }

    return l_norm_m;
}