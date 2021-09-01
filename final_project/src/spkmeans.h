#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/*
GENERAL USE MACROS
*/
#define INVALID_INPUT_MSG "Invalid Input!\n"
#define ERROR_MSG "An Error Has Occured\n"
#define MALLOC_ARR(type, size) (type *)malloc((size) * (sizeof(type)))
#define CALLOC_ARR(type, size) (type *)calloc(size, sizeof(type))
#define ASSERT_MALLOC(ptr) assert((ptr) != NULL)
#define ASSERT_WITH_MSG(cond, msg) \
    do                             \
    {                              \
        if (!(cond))               \
        {                          \
            printf(msg);           \
            exit(0);               \
        }                          \
    } while (0)
#define MALLOC_ARR_ASSERT(ptr, type, size)       \
    do                                           \
    {                                            \
        (ptr) = MALLOC_ARR(type, size);          \
        ASSERT_WITH_MSG(ptr != NULL, ERROR_MSG); \
    } while (0)
#define CALLOC_ARR_ASSERT(ptr, type, size)       \
    do                                           \
    {                                            \
        (ptr) = CALLOC_ARR(type, size);          \
        ASSERT_WITH_MSG(ptr != NULL, ERROR_MSG); \
    } while (0)
#define ZERO_4F "0.0000"
#define SIGN(X) (((X) >= 0) ? 1 : -1)
#define SQUARE(X) ((X) * (X))
#define EPSILON 1.0e-15
#define BASE_ARR_SIZE 50
#define ARR_SIZE_MULTIPLY 2
#define MAX_KMEANS_ITER 300
#define MAX_JACOBI_ITERS 100

/*
GENERAL REAL VALUES MATRIX
*/
typedef struct
{
    unsigned int rows;
    unsigned int cols;
    double **data;
} matrix;

/* allocating memory for new matrix and settig dimensions*/
matrix *init_mat(unsigned int rows, unsigned int cols);
void print_mat(matrix *mat);
void free_mat(matrix *mat);

/*
SYMMETRIC MATRIX
*/
typedef struct
{
    /* takes (dim^2 / 2) memory
    / stored as a "triangle"*/
    unsigned int dim;
    double **data;
} sym_matrix;

/* allocating memory for symmetric matrix*/
sym_matrix *init_sym_mat(unsigned int dim);
/* free symmetric matrix */
void free_sym_mat(sym_matrix *mat);
void print_sym_mat(sym_matrix *mat);
double get_val_sym(sym_matrix *mat, unsigned int row, unsigned col);
void set_val_sym(sym_matrix *mat, unsigned int row, unsigned col, double val);
sym_matrix *matrix_to_sym_matrix(matrix *mat);

/*
DIAGONAL MATRIX
*/
typedef struct
{
    /* takes (dim) memory */
    unsigned int dim;
    double *data;
} diag_matrix;

/* allocating memory (1d array) for diagonal matrix */
diag_matrix *init_diag_mat(unsigned int dim);
/* free diagonal matrix */
void free_diag_mat(diag_matrix *mat);
void print_diag_mat(diag_matrix *mat);

/*
WEIGHT MATRIX RELATED FUNCTIONS 
*/

/* Calculate weight matrix of a given matrix (allocating memory for weights matrix) */
sym_matrix *weights_mat(matrix *mat);

/*
DEGREE MATRIX RELATED FUNCTION  
*/

/* Calculate the degree matrix */
diag_matrix *degree_mat_from_data(matrix *mat);

/*
NORMALIZED LAPLACIAN MATRIX FUNCTION
*/

/* Comput the normalized laplacian matrix - symmetric matrix */
sym_matrix *l_norm_mat(matrix *mat);

/*
JACOBI ALGORITHM FUNCTIONS AND STRUCTURES
*/

/* Eiganvector */
typedef struct
{
    /*eiganvalue*/
    double e_val;
    /*eiganvector*/
    double *vec;
    /* index of vector in matrix*/
    unsigned int ind;
} e_vector;

typedef struct
{
    /* the matrix being diagonalized (symmetric) */
    sym_matrix *mat;
    /* rotation params, holds the index for the element
    with max abs value*/
    unsigned int row_ind;
    unsigned int col_ind;
    /* rotation params */
    double c;
    double s;
    /* array containing the column index of the element with the
    highest abs value (off diagonal) i.e max_indes[i] = j means mat[i][j] has
    max abs value in row i (right to the diagonal) */
    unsigned int *max_inds;
    /* the off diagonal diffrence between rotations, used to determine convergence */
    double off_diff;
    /* matrix containing the eiganvectors as columns */
    e_vector *e_mat;
} jacobi_matrix;

/* Initialize Jacobi matrix from a given symmetric matrix
 (not allocating new space for main matrix, ovverides s_mat)
 alocating memory for max_inds and e_mat */
jacobi_matrix *init_jac_mat(sym_matrix *s_mat);
/* free the allocated memory in jacobi matrix */
void free_jacobi(jacobi_matrix *j_mat);
/* getter for jacobi matrix */
double get_val(jacobi_matrix *j_mat, unsigned int row, unsigned int col);
/* set a specifi cell in jacobi matrix to given value */
void set_val(jacobi_matrix *j_mat, unsigned int row, unsigned int col, double val);
/* initialize eigan-vectors matrix (memeory allocation and set to identity matrix) */
void init_eigan_mat(jacobi_matrix *j_mat);
/* print the eigan-vectors matrix (each eigan-vector as a row) */
void print_e_mat(jacobi_matrix *j_mat);
/* print eigan-values as first row and then eigan-vectors as rows*/
void print_eigan_vectors_values(jacobi_matrix *j_mat);
/* compare two vectors to be used by qsort */
int cmp_vecs(const void *vec1, const void *vec2);
/* main jacobi algorithm */
void jacobi(jacobi_matrix *j_mat);

/*
THE EIGANGAP HEURISTIC
*/

unsigned int find_k(jacobi_matrix *jm);

/*
TEXT PARSING
*/

matrix *read_file_to_mat(FILE *file_pointer);

/*
NORMLIZING THE EIGAN MAT BEFORE K-MEAN
*/

/* step 4 in the Normalized Spectral Clustering Algorithm */
matrix *create_u_matrix(jacobi_matrix *j_mat, unsigned int k);
/* step 5 in the Normalized Spectral Clustering Algorithm - normlize each row */
void normlize_rows(matrix *mat);

typedef struct
{
    /* The data matrix - each row is a data point */
    matrix *vectors;
    /* Centroids matrix - each row is a centroid */
    matrix *centroids;
    /* Each row is a vector sum of all the data points assoiciated to the releveant
    cluster i.e if cluster j contains vec1, vec2, vec3 then the j-th row in clusters_sums is
    a vector sum of vec1 + vec2 + vec3 */
    matrix *clusters_sums;
    /* 1D array sized as number of clusters containing their size */
    unsigned int *clusters_size;
    /* 1D array sized as number of data points containing the cluster index for each data-point
    for an un-associated data point j whick_cluster[j] = -1 (Only initial value) */
    int *which_cluster;
} kmeans_data;

/* Main K-Means Algorithm :
    centroids matrix is the inital centroids (every row a centroid) K rows dim columns
        the centroids must be initialized by an external source
    vectors matrix is the data matrix, every row is a data-point
    centroids and vectors matrixes are asuumed to be initialized and allocated.
    vectors matrix doesn't contain duplicate rows (from the course forum)
    The centroids matrix will contain the final centroids */
void k_means(matrix *centroids, matrix *vectors, unsigned int k);

/* SPECTRAL CLUSTERING FUNCTCION
INPUT: k (may be 0 to trigger the The Eigengap Heuristic) 
       filename to csv formatted .txt or .csv files
OUTPUT: Centroids matrix containing the final centroids
*/
matrix *spectral_clustering(unsigned int k, matrix *data_mat);

/*
Handlers functions used in main C program to wrap each goal
*/

void handle_spectral_clustering(unsigned int k, matrix *data_mat);

void handle_weight_matrix(matrix *data_mat);

void handle_degree_matrix(matrix *data_mat);

void handle_lnorm_matrix(matrix *data_mat);

void handle_jacobi_matrix(matrix *data_mat);
