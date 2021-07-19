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
#define ZERO_4F "0.0000"
#define SIGN(X) (((X) >= 0) ? 1 : -1)

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
double get_val_sym(sym_matrix *mat, unsigned int row, unsigned col);
void set_val_sym(sym_matrix *mat, unsigned int row, unsigned col, double val);

/*
DIAGONAL MATRIX
*/
typedef struct
{
    // takes (dim) memory
    unsigned int dim;
    double *data;
} diag_matrix;

diag_matrix * init_diag_mat(unsigned int dim); 
void print_diag_mat(diag_matrix * mat);

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
diag_matrix * degree_mat(matrix * mat); 

/*
NORMALIZED LAPLACIAN MATRIX FUNCTION
*/

// Comput the normalized laplacian matrix - symmetric matrix
sym_matrix * l_norm_mat (matrix * mat);


/*
JACOBI ALGORITHM FUNCTIONS AND STRUCTURES
*/

typedef struct
{
    // the matrix being diagonalized (symmetric)
   sym_matrix * mat;
    // rotation params, holds the index for the element
    // with max abs value
   unsigned int row_ind;
   unsigned int col_ind;
   // rotation params
   double c;
   double s;
    // array containing the column index of the element with the
    // highest abs value (off diagonal) i.e max_indes[i] = j means mat[i][j] has 
    // max abs value in row i (right to the diagonal)
   unsigned int * max_inds;
} jacobi_matrix;

// Initialize Jacobi matrix from a given symmetric matrix
// (not allocating new space, ovverides s_mat)
// (find inital max element loc, c and s params)
jacobi_matrix * init_jac_mat(sym_matrix * s_mat);

// Initail abs value lookup
void max_abs_val_initial(jacobi_matrix * j_mat);
// getter for jacobi matrix
double get_val(jacobi_matrix * j_mat, unsigned int row, unsigned int col);

void set_val(jacobi_matrix * j_mat, unsigned int row, unsigned int col, double val);
// update the c and s params after new max is found
void update_c_s_params(jacobi_matrix * j_mat);
// rotate the jacobi matrix using c and s parmas
void rotate_jac(jacobi_matrix * j_mat);
// find new max efiicianly after rotation
void update_max_jac(jacobi_matrix * j_mat);
// update j_mat->max_inds for row i
void update_max_in_row(jacobi_matrix * j_mat, unsigned int row);
// upatate the toal max and row_ind, col_ind fields
void update_total_max_inds(jacobi_matrix * j_mat);

void jacobi(jacobi_matrix * j_mat);
