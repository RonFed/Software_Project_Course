#include <stdlib.h>
#include <assert.h>

/*
GENERAL USE MACROS
*/
#define MALLOC_ARR(type, size) (type *)malloc((size) * (sizeof(type)))
#define CALLOC_ARR(type, size) (type *)calloc(size, sizeof(type))
#define ASSERT_MALLOC(ptr) assert((ptr) != NULL)
#define MALLOC_ARR_ASSERT(ptr, type, size) \
    do                                     \
    {                                      \
        (ptr) = MALLOC_ARR(type, size);    \
        ASSERT_MALLOC(ptr);                \
    } while (0)
#define CALLOC_ARR_ASSERT(ptr, type, size) \
    do                                     \
    {                                      \
        (ptr) = CALLOC_ARR(type, size);    \
        ASSERT_MALLOC(ptr);                \
    } while (0)

/*
GENERAL MATRIX
*/
typedef struct
{
    unsigned int rows;
    unsigned int cols;
    double **data;
} matrix;

matrix * init_mat(unsigned int rows, unsigned int cols);
void print_mat(matrix *mat); 

/*
SYMMETRIC MATRIX
*/
typedef struct
{
    // takes (dim^2 / 2) memory
    // stored as a "triangle"
    unsigned int dim;
    double **data;
} sym_matrix;

sym_matrix * init_sym_mat(unsigned int dim); 
void print_sym_mat(sym_matrix *mat); 

/*
DIAGONAL MATRIX
*/
typedef struct
{
    // takes (dim) memory
    unsigned int dim;
    double *data;
} diag_matrix;

/*
WEIGHT MATRIX RELATED FUNCTIONS 
*/

// L2 Distance between two rows in matrix
double l2_norm_vectors(matrix * mat, unsigned int row1, unsigned int row2);
// exp norm as defined in the weights matrix
double exp_norm_vectors(matrix * mat, unsigned int row1, unsigned int row2);
// Calculate weight matrix of a given matrix
sym_matrix * weights_mat(matrix * mat);

/*
DEGREE MATRIX RELATED FUNCTION  
*/

// Sum of specific row in a symmetric matrix
double row_sum_sym_mat(sym_matrix * mat, unsigned int row);
// Calculate the degree matrix
diag_matrix * degree_mat(matrix * mat); // TODO !!