#include "spkmeans.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

/*
GENERAL MATRIX FUNCTIONS
*/
matrix* init_mat(unsigned int rows, unsigned int cols) {
    int i;
    matrix * mat = (matrix *) malloc(sizeof(matrix));
    mat->rows = rows;
    mat->cols = cols;
    MALLOC_ARR_ASSERT(mat->data, double*, rows);
    for (i = 0; i < rows; i++)
    {
       CALLOC_ARR_ASSERT((mat->data)[i], double, cols);
    }
    return mat;
}

void print_mat(matrix *mat) {
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

sym_matrix * init_sym_mat(unsigned int dim) {
    int i;
    sym_matrix * mat = (sym_matrix *) malloc(sizeof(sym_matrix));
    mat->dim = dim;
    MALLOC_ARR_ASSERT(mat->data, double*, dim);
    for (i = 0; i < dim; i++)
    {
       CALLOC_ARR_ASSERT((mat->data)[i], double, i+1);
    }
    return mat;
}

void print_sym_mat(sym_matrix *mat) {
    int i, j;
    unsigned int dim = mat->dim;
    for (i = 0; i < dim; i++)
    {
        for (j = 0; j < i+1; j++)
        {
            printf("%.4f,", (mat->data)[i][j]);
        }
        for (j = i+1; j < dim; j++)
        {
            printf("%.4f,", (mat->data)[j][i]);
        }
        printf("\n");
    }
}

/*
WEIGHTS MATRIX FUNCTIONS
*/

double l2_norm_vectors(matrix * mat, unsigned int row1, unsigned int row2) {
    assert((row1 < mat->rows) && (row2 < mat->rows));
    double result = 0 ,dif;
    int i;
    for (i = 0; i < mat->cols; i++)
    {
        dif = ((mat->data)[row1][i]) - ((mat->data)[row2][i]);
        result += (pow(dif, 2));
    }
    result = pow(result, 0.5);
    return result;
}

double exp_norm_vectors(matrix * mat, unsigned int row1, unsigned int row2) {
    double l2_norm = l2_norm_vectors(mat, row1, row2);
    double result;
    result = exp(-l2_norm/2);
    return result;
}

sym_matrix * weights_mat(matrix * mat) {
    sym_matrix * sm = init_sym_mat(mat->rows);
    int i,j;
    for (i = 0; i < mat->rows; i++)
    {
        for (j = 0; j < i; j++)
        {
            /* each cell in the weights matrix is the resul
            t of exp norm between two vectrose */
        (sm->data)[i][j] = exp_norm_vectors(mat,i,j);
        }
    }
    return sm;
}

/*
DEGREE MATRIX FUNCTIONS
*/

double row_sum_sym_mat(sym_matrix * mat, unsigned int row) {
    assert(row < mat->dim);
    double sum = 0;
    int j;
    for (j = 0; j < row+1; j++)
    {
        sum += (mat->data)[row][j];
    }
    for (j = row+1; j < mat->dim; j++)
    {
        sum += (mat->data)[j][row];
    }
    return sum;
}