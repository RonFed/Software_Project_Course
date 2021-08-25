#include "spkmeans.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

static int is_legal_goal(const char *str)
{
    return (!strcmp(str, "spk") ||
            !strcmp(str, "wam") ||
            !strcmp(str, "ddg") ||
            !strcmp(str, "lnorm") ||
            !strcmp(str, "jacobi"));
}

static int is_non_negative_int(const char *str)
{
    char digit = str[0];
    size_t i;
    for (i = 0; i < strlen(str); i++)
    {
        digit = str[i];
        if (digit < '0' || digit > '9')
        {
            return 0;
        }
    }
    return 1;
}

static void check_command_line_args(int argc, char const *argv[])
{
    /* 4 arguments including program name */
    ASSERT_WITH_MSG(argc == 4, INVALID_INPUT_MSG);
    /* k is a non- negative integer */
    ASSERT_WITH_MSG(is_non_negative_int(argv[1]), INVALID_INPUT_MSG);
    /* valid goal */
    ASSERT_WITH_MSG(is_legal_goal(argv[2]), INVALID_INPUT_MSG);
}

/* Main C interface */
int main(int argc, char const *argv[])
{
    FILE *file_pointer;
    char const *goal;
    char const *filename;
    unsigned int k;
    matrix *data_mat;

    /* Validate the command line arguments */
    check_command_line_args(argc, argv);

    k = atoi(argv[1]);
    goal = argv[2];
    filename = argv[3];

    /* Read the file to matrix */
    file_pointer = fopen(filename, "r");
    ASSERT_WITH_MSG(file_pointer != NULL, ERROR_MSG);

    data_mat = read_file_to_mat(file_pointer);
    fclose(file_pointer);

    /* Handle the goal */
    if (strcmp(goal, "spk") == 0)
    {
        handle_spectral_clustering(k, data_mat);
    }
    else if (strcmp(goal, "wam") == 0)
    {
        handle_weight_matrix(data_mat);
    }
    else if (strcmp(goal, "ddg") == 0)
    {
        handle_degree_matrix(data_mat);
    }
    else if (strcmp(goal, "lnorm") == 0)
    {
        handle_lnorm_matrix(data_mat);
    }
    else if (strcmp(goal, "jacobi") == 0)
    {
        handle_jacobi_matrix(data_mat);
    }
    else
    {
        ASSERT_WITH_MSG(0, INVALID_INPUT_MSG);
    }
    return 0;
}

void handle_spectral_clustering(unsigned int k, matrix *data_mat)
{
    matrix *centroids = spectral_clustering(k, data_mat);
    print_mat(centroids);
    free_mat(centroids);
}

void handle_weight_matrix(matrix *data_mat)
{
    sym_matrix *weights_m = weights_mat(data_mat);
    print_sym_mat(weights_m);
    free_sym_mat(weights_m);
    free_mat(data_mat);
}

void handle_degree_matrix(matrix *data_mat)
{
    diag_matrix *degree_m = degree_mat_from_data(data_mat);
    print_diag_mat(degree_m);
    free_diag_mat(degree_m);
    free_mat(data_mat);
}

void handle_lnorm_matrix(matrix *data_mat)
{
    sym_matrix *l_norm_m = l_norm_mat(data_mat);
    print_sym_mat(l_norm_m);
    free_sym_mat(l_norm_m);
    free_mat(data_mat);
}

void handle_jacobi_matrix(matrix *data_mat)
{
    sym_matrix *sym_m = matrix_to_sym_matrix(data_mat);
    jacobi_matrix *jac_m = init_jac_mat(sym_m);
    jacobi(jac_m);
    print_eigan_vectors_values(jac_m);
    free_jacobi(jac_m);
    free_mat(data_mat);
}

/*
GENERAL MATRIX FUNCTIONS
*/

/* Allocate memory for a new matrix object sized rows * cols of double */
matrix *init_mat(unsigned int rows, unsigned int cols)
{
    unsigned int i;
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

/* Print matrix object in CSV format */
void print_mat(matrix *mat)
{
    unsigned int i, j;
    for (i = 0; i < mat->rows; i++)
    {
        for (j = 0; j < mat->cols - 1; j++)
        {
            printf("%.4f,", (mat->data)[i][j]);
        }
        printf("%.4f", (mat->data)[i][mat->cols - 1]);
        printf("\n");
    }
}

void free_mat(matrix *mat)
{
    unsigned int i;
    for (i = 0; i < mat->rows; i++)
    {
        free((mat->data)[i]);
    }
    free(mat->data);
    free(mat);
}

sym_matrix *init_sym_mat(unsigned int dim)
{
    unsigned int i;
    sym_matrix *mat = (sym_matrix *)malloc(sizeof(sym_matrix));
    mat->dim = dim;
    MALLOC_ARR_ASSERT(mat->data, double *, dim);
    for (i = 0; i < dim; i++)
    {
        CALLOC_ARR_ASSERT((mat->data)[i], double, dim - i);
    }
    return mat;
}

void free_sym_mat(sym_matrix *mat)
{
    unsigned int i;
    for (i = 0; i < mat->dim; i++)
    {
        free((mat->data)[i]);
    }
    free(mat->data);
    free(mat);
}

/* Converts a matrix object to symmetric matrix object by allocating memory for new
    symmetric matrix and taking the upper triangle elements of the given matrix
    For correctness, Assumes the original matrix is indeed symmetric */
sym_matrix *matrix_to_sym_matrix(matrix *mat)
{
    sym_matrix *result;
    unsigned int i, j;
    result = init_sym_mat(mat->rows);
    for (i = 0; i < mat->rows; i++)
    {
        for (j = i; j < mat->rows; j++)
        {
            set_val_sym(result, i, j, (mat->data)[i][j]);
        }
    }
    return result;
}

void print_sym_mat(sym_matrix *mat)
{
    unsigned int i, j;
    unsigned int dim = mat->dim;
    for (i = 0; i < dim; i++)
    {
        for (j = 0; j < i; j++)
        {
            printf("%.4f,", (mat->data)[j][i - j]);
        }
        for (j = i; j < dim - 1; j++)
        {
            printf("%.4f,", (mat->data)[i][j - i]);
        }
        printf("%.4f", (mat->data)[i][dim - 1 - i]);
        printf("\n");
    }
}
/* Getter for symmetric matrix 
logic is requiered becouse of the memory allocatin for symmetric matrx 
*/
double get_val_sym(sym_matrix *mat, unsigned int row, unsigned col)
{
    if (col >= row)
    {
        return (mat->data)[row][col - row];
    }
    return (mat->data)[col][row - col];
}
void set_val_sym(sym_matrix *mat, unsigned int row, unsigned col, double val)
{
    if (col >= row)
    {
        (mat->data)[row][col - row] = val;
    }
    else
    {
        (mat->data)[col][row - col] = val;
    }
}

diag_matrix *init_diag_mat(unsigned int dim)
{
    diag_matrix *mat = (diag_matrix *)malloc(sizeof(diag_matrix));
    mat->dim = dim;
    MALLOC_ARR_ASSERT(mat->data, double, dim);
    return mat;
}

void free_diag_mat(diag_matrix *mat)
{
    free(mat->data);
    free(mat);
}

/* Print diagonal matrix, all elments except the diagonal elements are 0 */
void print_diag_mat(diag_matrix *mat)
{
    unsigned int i, j;
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

/* L2 Distance between two rows in matrix */
static double l2_norm_vectors(matrix *mat, unsigned int row1, unsigned int row2)
{
    double result, dif;
    unsigned int i;
    result = 0;
    assert((row1 < mat->rows) && (row2 < mat->rows));
    for (i = 0; i < mat->cols; i++)
    {
        dif = ((mat->data)[row1][i]) - ((mat->data)[row2][i]);
        result += (dif * dif);
    }
    result = sqrt(result);
    return result;
}

/* exp norm as defined in the weights matrix */
static double exp_norm_vectors(matrix *mat, unsigned int row1, unsigned int row2)
{
    double l2_norm = l2_norm_vectors(mat, row1, row2);
    double result;
    result = exp(-l2_norm / 2);
    return result;
}

sym_matrix *weights_mat(matrix *mat)
{
    sym_matrix *sm = init_sym_mat(mat->rows);
    unsigned int i, j;
    for (i = 0; i < mat->rows; i++)
    {
        for (j = 0; j < i; j++)
        {
            /* each cell in the weights matrix is the result
            of exp norm between two vectrose */
            set_val_sym(sm, i, j, exp_norm_vectors(mat, i, j));
        }
    }
    return sm;
}

/*
DEGREE MATRIX FUNCTIONS
*/

double row_sum_sym_mat(sym_matrix *mat, unsigned int row)
{
    double sum;
    unsigned int j;
    sum = 0;
    assert(row < mat->dim);
    for (j = 0; j < mat->dim; j++)
    {
        sum += get_val_sym(mat, row, j);
    }
    return sum;
}

/* Compute the degree matrix from data matrix (including weight matrix inside this function) */
diag_matrix *degree_mat_from_data(matrix *mat)
{
    unsigned int i;
    /* allocating memory for degree mat as diagonal matrix */
    diag_matrix *degree_m = init_diag_mat(mat->rows);
    sym_matrix *weight_m = weights_mat(mat);
    for (i = 0; i < mat->rows; i++)
    {
        /* sum of a row in the weights matrix */
        (degree_m->data)[i] = row_sum_sym_mat(weight_m, i);
    }
    free_sym_mat(weight_m);
    return degree_m;
}

/* Compute the degree matrix from a given weights matrix */
static diag_matrix *degree_mat_from_wieght_mat(sym_matrix *weights_m)
{
    unsigned int i;
    /* allocating memory for degree mat as diagonal matrix */
    diag_matrix *degree_m = init_diag_mat(weights_m->dim);
    for (i = 0; i < weights_m->dim; i++)
    {
        /* sum of a row in the weights matrix */
        (degree_m->data)[i] = row_sum_sym_mat(weights_m, i);
    }
    return degree_m;
}

/* Compute normalized Laplacian matrix from data matrix */
sym_matrix *l_norm_mat(matrix *mat)
{
    sym_matrix *l_norm_m = init_sym_mat(mat->rows);
    sym_matrix *weight_m = weights_mat(mat);
    diag_matrix *degree_m = degree_mat_from_wieght_mat(weight_m);
    double curr_weight, curr_d_i, curr_d_j, d_sqrt;
    unsigned int i, j;
    for (i = 0; i < l_norm_m->dim; i++)
    {
        curr_d_i = degree_m->data[i];
        for (j = 0; j <= i; j++)
        {
            curr_d_j = degree_m->data[j];
            curr_weight = get_val_sym(weight_m, i, j);
            if (j == i)
            {
                /* on the diagonal w(i,i) = 0 => l_norm(i,i) = 1 */
                set_val_sym(l_norm_m, i, j, 1);
            }
            else
            {
                d_sqrt = 1 / (sqrt(curr_d_i * curr_d_j));
                set_val_sym(l_norm_m, i, j, -d_sqrt * curr_weight);
            }
        }
    }
    free_diag_mat(degree_m);
    free_sym_mat(weight_m);
    return l_norm_m;
}

/* ----------------------------------------------------------------------------------//
///////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  Jacobi Algorithm   //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------------------*/

void init_eigan_mat(jacobi_matrix *j_mat)
{
    unsigned int i, n;
    n = (j_mat->mat)->dim;
    CALLOC_ARR_ASSERT(j_mat->e_mat, e_vector, n);
    for (i = 0; i < n; i++)
    {
        /* Allocate memory for each eiganvector */
        CALLOC_ARR_ASSERT((j_mat->e_mat)[i].vec, double, n);
        /* setting the index of eigan-vector */
        (j_mat->e_mat)[i].ind = i;
        /* initialize the matrix to the identity */
        ((j_mat->e_mat)[i].vec)[i] = 1;
    }
}

jacobi_matrix *init_jac_mat(sym_matrix *s_mat)
{
    jacobi_matrix *jac_mat = (jacobi_matrix *)malloc(sizeof(jacobi_matrix));
    jac_mat->mat = s_mat;
    init_eigan_mat(jac_mat);
    /* calloc for max inds array */
    CALLOC_ARR_ASSERT(jac_mat->max_inds, unsigned int, s_mat->dim);
    jac_mat->off_diff = 0;
    return jac_mat;
}

void free_jacobi(jacobi_matrix *j_mat)
{
    unsigned int i, n;
    n = (j_mat->mat)->dim;
    /* free eigan-mat */
    for (i = 0; i < n; i++)
    {
        free((j_mat->e_mat)[i].vec);
    }
    free(j_mat->e_mat);
    /* free max-inds array */
    free(j_mat->max_inds);
    free_sym_mat(j_mat->mat);
    free(j_mat);
}

/* Initail abs value lookup */
static void max_abs_val_initial(jacobi_matrix *j_mat)
{
    double max_total = -1, max_in_row, curr_abs;
    unsigned int i, j;
    sym_matrix *sm = j_mat->mat;
    unsigned int max_row = 0;
    for (i = 0; i < sm->dim; i++)
    {
        /* for each row find the max element */
        max_in_row = -1;
        /* max abs value off - diagonal */
        for (j = i + 1; j < (sm->dim); j++)
        {
            curr_abs = fabs(get_val_sym(sm, i, j));
            if (curr_abs > max_in_row)
            {
                max_in_row = curr_abs;
                (j_mat->max_inds)[i] = j;
            }
        }
        if (max_in_row > max_total)
        {
            max_total = max_in_row;
            max_row = i;
        }
    }
    j_mat->row_ind = max_row;
    j_mat->col_ind = (j_mat->max_inds)[max_row];
}

double get_val(jacobi_matrix *j_mat, unsigned int row, unsigned int col)
{
    return get_val_sym(j_mat->mat, row, col);
}

void set_val(jacobi_matrix *j_mat, unsigned int row, unsigned int col, double val)
{
    set_val_sym(j_mat->mat, row, col, val);
}

/* this fuction is called after the max value row and col indexes are updated */
void update_c_s_params(jacobi_matrix *j_mat)
{
    unsigned int i, j;
    double theta, c, s, t, a_ii, a_jj, a_ij;
    i = j_mat->row_ind;
    j = j_mat->col_ind;
    a_ii = get_val(j_mat, i, i);
    a_jj = get_val(j_mat, j, j);
    a_ij = get_val(j_mat, i, j);
    theta = (a_jj - a_ii) / (2 * a_ij);
    t = SIGN(theta) / (fabs(theta) + sqrt((theta * theta) + 1));
    c = 1 / (sqrt((t * t) + 1));
    s = t * c;
    j_mat->c = c;
    j_mat->s = s;
}

void update_eigan_mat(jacobi_matrix *j_mat)
{
    e_vector *e_mat_data = j_mat->e_mat;
    unsigned int i, j, r;
    double c, s, curr_i, curr_j;
    e_vector vec_i, vec_j;
    i = j_mat->row_ind;
    j = j_mat->col_ind;
    c = j_mat->c;
    s = j_mat->s;
    /* update only the i-th and j-th eigan-vectors */
    vec_i = e_mat_data[i];
    vec_j = e_mat_data[j];
    for (r = 0; r < (j_mat->mat->dim); r++)
    {
        /* the update according to the rotation matrix
        update row i - the i-th and j-th eigan-vectors */
        curr_i = (vec_i.vec)[r];
        curr_j = (vec_j.vec)[r];
        (vec_i.vec)[r] = (c * curr_i) - (s * curr_j);
        (vec_j.vec)[r] = (s * curr_i) + (c * curr_j);
    }
}

void print_e_mat(jacobi_matrix *j_mat)
{
    e_vector *e_mat_data = j_mat->e_mat;
    unsigned int n = (j_mat->mat)->dim;
    unsigned int i, j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n - 1; j++)
        {
            printf("%.4f,", (e_mat_data[i].vec)[j]);
        }
        printf("%.4f", (e_mat_data[i].vec)[n - 1]);
        printf("\n");
    }
}

void print_eigan_vectors_values(jacobi_matrix *j_mat)
{
    e_vector *e_mat_data = j_mat->e_mat;
    unsigned int n = (j_mat->mat)->dim;
    unsigned int i;
    for (i = 0; i < n - 1; i++)
    {
        printf("%.4f,", (e_mat_data[i].e_val));
    }
    printf("%.4f\n", (e_mat_data[n - 1].e_val));
    print_e_mat(j_mat);
}

void rotate_jacobi(jacobi_matrix *j_mat)
{
    sym_matrix *sm_data = j_mat->mat;
    unsigned int i, j, r;
    double curr_i, curr_j, new_i, new_j, new_ii, new_jj, c, s;
    double positive_diff = 0.0, negative_diff = 0.0;
    i = j_mat->row_ind;
    j = j_mat->col_ind;
    c = j_mat->c;
    s = j_mat->s;
    /* update the i-th and j-th columns */
    for (r = 0; r < sm_data->dim; r++)
    {
        if (r != j && r != i)
        {
            curr_i = get_val(j_mat, r, i);
            curr_j = get_val(j_mat, r, j);
            new_i = (c * curr_i) - (s * curr_j);
            new_j = (c * curr_j) + (s * curr_i);
            set_val(j_mat, r, i, new_i);
            set_val(j_mat, r, j, new_j);

            /* updating the diffrence in off diagonal elements field */
            positive_diff += (SQUARE(curr_i) + SQUARE(curr_j));
            negative_diff += (SQUARE(new_i) + SQUARE(new_j));
        }
    }
    /* update relevant diagonal elements */
    curr_i = get_val(j_mat, i, i);
    curr_j = get_val(j_mat, j, j);
    new_ii = (c * c * curr_i) + (s * s * curr_j) - (2 * s * c * (get_val(j_mat, i, j)));
    new_jj = (s * s * curr_i) + (c * c * curr_j) + (2 * s * c * (get_val(j_mat, i, j)));
    set_val(j_mat, i, i, new_ii);
    set_val(j_mat, j, j, new_jj);
    positive_diff += SQUARE(get_val(j_mat, i, j));
    set_val(j_mat, i, j, 0);
    /* multiply by 2 since the matrix is symmetric */
    j_mat->off_diff = 2 * (positive_diff - negative_diff);
}

/* update j_mat->max_inds for row i */
static void update_max_in_row(jacobi_matrix *j_mat, unsigned int row)
{
    unsigned int j, old_max_col;
    double curr, max;
    old_max_col = (j_mat->max_inds)[row];
    max = fabs(get_val(j_mat, row, old_max_col));
    /* only checking elements right to the diagonal */
    for (j = row + 1; j < (j_mat->mat)->dim; j++)
    {
        curr = fabs(get_val(j_mat, row, j));
        if (curr > max)
        {
            max = curr;
            (j_mat->max_inds)[row] = j;
        }
    }
}

/* upatate the toal max and row_ind, col_ind fields */
static void update_total_max_inds(jacobi_matrix *j_mat)
{
    double max = -1, curr;
    unsigned int r, max_row = 0, max_col = 0;
    /* last row doesn't contain relevant elements (just one element on diagonal) */
    for (r = 0; r < (j_mat->mat)->dim - 1; r++)
    {
        /* max value in line r */
        curr = fabs(get_val(j_mat, r, (j_mat->max_inds)[r]));
        if (curr > max)
        {
            max = curr;
            max_row = r;
            max_col = (j_mat->max_inds)[r];
        }
    }
    j_mat->row_ind = max_row;
    j_mat->col_ind = max_col;
}

/* find new max efiicianly after rotation */
static void update_max_jac(jacobi_matrix *j_mat, int is_not_first)
{
    if (is_not_first)
    {
        /* finding the new max off-diagonal abs value */
        unsigned int r, curr_row_max_ind, curr_col_max_ind;
        unsigned int *max_arr;
        curr_row_max_ind = j_mat->row_ind;
        curr_col_max_ind = j_mat->col_ind;
        /* first check the rows being rotated */
        update_max_in_row(j_mat, curr_row_max_ind);
        update_max_in_row(j_mat, curr_col_max_ind);
        max_arr = j_mat->max_inds;
        for (r = 0; r < (j_mat->mat)->dim; r++)
        {
            /* all rows excpet the two that were updated */
            if ((r != curr_row_max_ind) && (r != curr_col_max_ind))
            {
                /* the max in the row r is in a cloumn that wasn't changed in the last rotation */
                if ((max_arr[r] != curr_row_max_ind) && (max_arr[r] != curr_col_max_ind))
                {
                    /* only two columns have changed in this row to check */
                    if (fabs(get_val(j_mat, r, curr_row_max_ind)) > fabs(get_val(j_mat, r, max_arr[r])))
                    {
                        max_arr[r] = curr_row_max_ind;
                    }
                    if (fabs(get_val(j_mat, r, curr_col_max_ind)) > fabs(get_val(j_mat, r, max_arr[r])))
                    {
                        max_arr[r] = curr_col_max_ind;
                    }
                }
                else
                {
                    /* the column in row r containing the max has been changed in the rotation
                    find the max in row r */
                    update_max_in_row(j_mat, r);
                }
            }
        }
        update_total_max_inds(j_mat);
    }
    else
    {
        /* only occurs once in first iteration of jacobi */
        max_abs_val_initial(j_mat);
    }
}

/* get the value of max abs element off-diagonal */
static double get_val_max_off_diagonal(jacobi_matrix *j_mat)
{
    return get_val(j_mat, j_mat->row_ind, j_mat->col_ind);
}

/* comparison function for qsort */
int cmp_vecs(const void *vec1, const void *vec2)
{
    /* if the two vectors have the same eigan-value, sort them by index (the smaller index first)
    else sort them by the eigan-value - smaller one comes first */
    double val1, val2;
    val1 = ((e_vector *)vec1)->e_val;
    val2 = ((e_vector *)vec2)->e_val;
    if (val1 == val2)
    {
        return ((e_vector *)vec1)->ind - ((e_vector *)vec2)->ind;
    }
    else
    {
        if (val1 > val2)
        {
            return 1;
        }
    }
    return -1;
}

void set_eigan_values(jacobi_matrix *j_mat)
{
    unsigned int n, i;
    n = (j_mat->mat)->dim;
    for (i = 0; i < n; i++)
    {
        (j_mat->e_mat)[i].e_val = get_val(j_mat, i, i);
    }
}

/* main Jacobi algorithm */
void jacobi(jacobi_matrix *j_mat)
{
    int iters = 0;
    do
    {
        update_max_jac(j_mat, iters);
        /* if diagonal then stop, avoiding divide by zero */
        if (get_val_max_off_diagonal(j_mat) == 0)
        {
            break;
        }
        update_c_s_params(j_mat);
        rotate_jacobi(j_mat);
        update_eigan_mat(j_mat);
        iters++;
    } while (iters < MAX_JACOBI_ITERS && j_mat->off_diff > EPSILON);
    /* setting each eigan-value to the coresponding eigan-vector */
    set_eigan_values(j_mat);
}

/* ----------------------------------------------------------------------------------//
///////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  Text Parsing   //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------------------*/

/* File is open, read the first line into first_line - looking for CSV format
first line is a pointer to a pointer that may change using realloc (since dimension is not known)
return the dimension of the data i.e how many elements found in the first row (assuming each row will have
the same dimension) */
static unsigned int find_dimension_from_first_line(FILE *file_pointer, double **first_line)
{
    double current_element;
    char current_delimeter;
    unsigned int dimension = 0;
    unsigned int current_arr_len = BASE_ARR_SIZE;

    while (fscanf(file_pointer, "%lf%c", &current_element, &current_delimeter) == 2)
    {
        /* rarely happens - first line is over BASE_ARR_SIZE elements */
        if (dimension == current_arr_len)
        {
            current_arr_len *= ARR_SIZE_MULTIPLY;
            (*first_line) = realloc(*first_line, current_arr_len * sizeof(double));
            assert(*first_line);
        }
        (*first_line)[dimension] = current_element;
        dimension++;
        if (current_delimeter == '\n')
        {
            break;
        }
    }
    (*first_line) = realloc(*first_line, dimension * sizeof(double));
    return dimension;
}

/* read data from string buffer to line (assuming dimension elements in line) */
static void read_line_to_row(char *buffer, double *line, unsigned int dimension)
{
    double current_element;
    unsigned int i;
    int chars_read;

    for (i = 0; i < dimension; i++)
    {
        sscanf(buffer, "%lf,%n", &current_element, &chars_read);
        line[i] = current_element;
        /* update the pointer to the buffer */
        buffer += chars_read;
    }
}

/* main parsing - return matrix object from csv formatted file */
matrix *read_file_to_mat(FILE *file_pointer)
{
    matrix *data_mat;
    double **data;
    char *buffer;
    unsigned int lines_count = 0;
    unsigned int current_arr_len = BASE_ARR_SIZE;
    unsigned int buffer_size = 0;
    unsigned int dimension = 0;
    data_mat = (matrix *)malloc(sizeof(matrix));

    MALLOC_ARR_ASSERT(data, double *, BASE_ARR_SIZE);

    /* first line analysis to determine the dimension of the csv format (columns number) */
    CALLOC_ARR_ASSERT(data[0], double, BASE_ARR_SIZE);
    /* passing pointer to first line since it might change using realloc (dimension is not known yet) */
    dimension = find_dimension_from_first_line(file_pointer, &data[0]);
    lines_count++;

    /* assuming each numeber is less than 100 digits */
    buffer_size = dimension * 100;
    MALLOC_ARR_ASSERT(buffer, char, buffer_size);

    while (fgets(buffer, buffer_size, file_pointer))
    {
        /* multiply the size of the array when it's full */
        if (lines_count == current_arr_len)
        {
            current_arr_len *= ARR_SIZE_MULTIPLY;
            data = (double **)realloc(data, current_arr_len * sizeof(double *));
            ASSERT_WITH_MSG(data != NULL, ERROR_MSG);
        }
        /* Allocate memory for current row*/
        CALLOC_ARR_ASSERT(data[lines_count], double, dimension);
        /* Read data from buffer to current row and increase row counter*/
        read_line_to_row(buffer, data[lines_count++], dimension);
    }
    free(buffer);
    /* triming the final array */
    data = (double **)realloc(data, lines_count * sizeof(double *));
    ASSERT_WITH_MSG(data != NULL, ERROR_MSG);
    data_mat->cols = dimension;
    data_mat->rows = lines_count;
    data_mat->data = data;
    return data_mat;
}

/* ----------------------------------------------------------------------------------//
///////////////////////////////////////////////////////////////////////////////////////
//////////////////////////  The Eigengap Heuristic   //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------------------*/

unsigned int find_k(jacobi_matrix *j_mat)
{
    /* the eigan vectors matrix */
    e_vector *eigan_mat = j_mat->e_mat;
    unsigned int top_limit = ((j_mat->mat)->dim) / 2;
    unsigned int i, max_gap_index = 0;
    double current_gap, max_gap = 0;

    for (i = 0; i < top_limit; i++)
    {
        /* gap between two eigan-values */
        current_gap = fabs(eigan_mat[i + 1].e_val - eigan_mat[i].e_val);
        /* in case of equality between gaps take the lowest index (done by weak inequality) */
        if (current_gap > max_gap)
        {
            max_gap = current_gap;
            max_gap_index = i;
        }
    }
    return max_gap_index + 1;
}

/* ----------------------------------------------------------------------------------//
///////////////////////////////////////////////////////////////////////////////////////
///////// Steps 4-5 in Normalized Spectral Clustering Algorithm  //////////////////////
///////////////////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------------------*/

/* form the matrix containing the first k eignan-vectors as colmuns 
(allocates memory for new matrix) */
matrix *create_u_matrix(jacobi_matrix *j_mat, unsigned int k)
{
    unsigned int i, j;
    matrix *u_mat;
    unsigned int n;
    n = (j_mat->mat)->dim;
    u_mat = init_mat(n, k);
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < k; j++)
        {
            (u_mat->data)[i][j] = (j_mat->e_mat)[j].vec[i];
        }
    }
    return u_mat;
}

/* calculate the l2 norm of a given row in matrix */
static double l2_norm_row_in_mat(matrix *mat, unsigned int row)
{
    double result = 0;
    unsigned int i;
    for (i = 0; i < mat->cols; i++)
    {
        result += SQUARE((mat->data)[row][i]);
    }
    result = sqrt(result);
    return result;
}

/* normlaize each row */
void normlize_rows(matrix *mat)
{
    unsigned int i, j;
    double norm;
    for (i = 0; i < mat->rows; i++)
    {
        norm = l2_norm_row_in_mat(mat, i);
        if (norm != 0)
        {
            for (j = 0; j < mat->cols; j++)
            {
                (mat->data)[i][j] /= norm;
            }
        }
    }
}

/* ----------------------------------------------------------------------------------//
///////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  K-Means  //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------------------*/

/* Compute the euclidian distance between two vectors */
static double distacne(double *a, double *b, unsigned int dim)
{
    double sum = 0;
    unsigned int i;
    for (i = 0; i < dim; i++)
    {
        sum += SQUARE(a[i] - b[i]);
    }
    return sum;
}

/* Find the closest cluster to a given vector */
static int closest_cluster_index(kmeans_data *kmeans_data, unsigned int vector_ind)
{
    unsigned int i, dim, k;
    int min_dist_index = 0;
    double min_dist;
    double *vector;
    double *centroied;
    dim = (kmeans_data->vectors)->cols;
    k = (kmeans_data->centroids)->rows;
    vector = ((kmeans_data->vectors)->data)[vector_ind];
    centroied = ((kmeans_data->centroids)->data)[0];
    min_dist = distacne(vector, centroied, dim);
    for (i = 0; i < k; i++)
    {
        double current_dist;
        centroied = ((kmeans_data->centroids)->data)[i];
        current_dist = distacne(vector, centroied, k);
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

/* Initialize kmeans_data :
    vectors and centroids matrixes are given
    allocate memory for clusters_sums matrix, cluster_size array and which_cluster array */
static kmeans_data *init_kmeans(matrix *vectors, matrix *centroids, unsigned int k)
{
    unsigned int i;
    kmeans_data *result;
    matrix *clusters_sum;
    unsigned int *clusters_size;
    int *which_cluster;
    clusters_sum = init_mat(k, vectors->cols);
    CALLOC_ARR_ASSERT(clusters_size, unsigned int, k);
    MALLOC_ARR_ASSERT(which_cluster, int, vectors->rows);
    for (i = 0; i < vectors->rows; i++)
    {
        which_cluster[i] = -1;
    }

    result = (kmeans_data *)malloc(sizeof(kmeans_data));
    ASSERT_WITH_MSG(result != NULL, ERROR_MSG);

    result->vectors = vectors;
    result->centroids = centroids;
    result->clusters_size = clusters_size;
    result->which_cluster = which_cluster;
    result->clusters_sums = clusters_sum;
    return result;
}

/* Remove a given vector to a target cluster, 
update sum and cluster size accordingly */
static void remove_vector_from_cluster(kmeans_data *kmeans_data, unsigned int vector_ind)
{
    unsigned int current_cluster, i;
    matrix *clusters_sum;
    double *vector = ((kmeans_data->vectors)->data)[vector_ind];
    current_cluster = (kmeans_data->which_cluster)[vector_ind];

    (kmeans_data->clusters_size)[current_cluster] -= 1;

    clusters_sum = kmeans_data->clusters_sums;
    for (i = 0; i < clusters_sum->cols; i++)
    {
        (clusters_sum->data)[current_cluster][i] -= vector[i];
    }
}

/* Add a given vector to a target cluster, 
update sum and cluster size accordingly */
static void add_vector_to_cluster(kmeans_data *kmeans_data, unsigned int vector_ind, unsigned int target_cluster)
{
    unsigned int i;
    matrix *clusters_sum;
    double *vector = ((kmeans_data->vectors)->data)[vector_ind];

    (kmeans_data->clusters_size)[target_cluster] += 1;

    clusters_sum = kmeans_data->clusters_sums;
    for (i = 0; i < clusters_sum->cols; i++)
    {
        (clusters_sum->data)[target_cluster][i] += vector[i];
    }

    (kmeans_data->which_cluster)[vector_ind] = target_cluster;
}

/* Update the centroids :
    given the vector sum of each cluster and number of vectors in each cluster
    computer the updated centroid of each cluster by taking the mean */
static void update_centroids(kmeans_data *kmeans_data)
{
    unsigned int i, j;
    matrix *centroids;
    matrix *centroids_sums;
    unsigned int *clusters_size;
    centroids = kmeans_data->centroids;
    centroids_sums = kmeans_data->clusters_sums;
    clusters_size = kmeans_data->clusters_size;
    for (i = 0; i < centroids->rows; i++)
    {
        if (clusters_size[i] != 0)
        {
            for (j = 0; j < centroids->cols; j++)
            {
                (centroids->data)[i][j] = ((centroids_sums->data)[i][j]) / (clusters_size[i]);
            }
        }
    }
}

static void free_kmeans(kmeans_data *kmeans_data)
{
    /* Free the un-relevant data and keep the centroids*/
    free_mat(kmeans_data->clusters_sums);
    free(kmeans_data->which_cluster);
    free(kmeans_data->clusters_size);
    free(kmeans_data);
}

/* Main K-Means Algorithm :
    centroids matrix is the inital centroids (every row a centroid) K rows dim columns
        the centroids must be initialized by an external source
    vectors matrix is the data matrix, every row is a data-point
    centroids and vectors matrixes are asuumed to be initialized and allocates
    The centroids matrix will contain the final centroids */
void k_means(matrix *centroids, matrix *vectors, unsigned int k)
{
    unsigned int i, current_vector, vectors_num;
    int new_closest_cluster;
    kmeans_data *kmeans_data;
    kmeans_data = init_kmeans(vectors, centroids, k);
    vectors_num = vectors->rows;

    for (i = 0; i < MAX_ITER; i++)
    {
        int not_converged = 0;
        /* Iterate over all data points*/
        for (current_vector = 0; current_vector < vectors_num; current_vector++)
        {
            /* find the closest cluster*/
            new_closest_cluster = closest_cluster_index(kmeans_data, current_vector);
            if (new_closest_cluster != (kmeans_data->which_cluster)[current_vector])
            {
                /* Found a closer cluster */
                if ((kmeans_data->which_cluster)[current_vector] != -1)
                {
                    /* Data point was already associated with a cluster
                    Remove it from cluster */
                    remove_vector_from_cluster(kmeans_data, current_vector);
                }
                /* Add the data point to new cluster*/
                add_vector_to_cluster(kmeans_data, current_vector, new_closest_cluster);
                not_converged = 1;
            }
        }
        /* If last main iteration didn't move any data point to a new cluster --> converged */
        if (not_converged == 0)
        {
            break;
        }
        update_centroids(kmeans_data);
    }
    /* Free the un-relevant data and keep the centroids and vectors data */
    free_kmeans(kmeans_data);
}

/* Simple init for centroids - the first k vectors in the vectors data
   Allocates memory for the returned matrix
   This should be called within the C implementation only */
static matrix *centroids_init_simple(matrix *vectors, unsigned int k)
{
    unsigned int i, j;
    matrix *centeroids = init_mat(k, vectors->cols);
    for (i = 0; i < k; i++)
    {
        for (j = 0; j < vectors->cols; j++)
        {
            (centeroids->data)[i][j] = (vectors->data)[i][j];
        }
    }
    return centeroids;
}

matrix *spectral_clustering(unsigned int k, matrix *data_mat)
{
    sym_matrix *l_norm;
    jacobi_matrix *jacobi_mat;
    matrix *u_mat;
    matrix *centeroids;

    /* the input must be in the form of k < N 
    (number of clusters is less then number of data-points */
    ASSERT_WITH_MSG(k < data_mat->rows, INVALID_INPUT_MSG);

    /* compute normalized laplacian */
    l_norm = l_norm_mat(data_mat);
    free_mat(data_mat);

    /* find eiganvectors and eiganvalues using the Jacobi algorithm */
    jacobi_mat = init_jac_mat(l_norm);
    jacobi(jacobi_mat);

    /* sorting the eigan-vectors */
    qsort(jacobi_mat->e_mat, (jacobi_mat->mat)->dim, sizeof(e_vector), cmp_vecs);

    if (k == 0)
    {
        k = find_k(jacobi_mat);
    }
    /* Step 4 - U matrix containing the first k eiganvectors as columns */
    u_mat = create_u_matrix(jacobi_mat, k);

    free_jacobi(jacobi_mat);
    /* Step 5 - Create T matrix by normalize each row in U
    (not allocating new memory) - done in-place */
    normlize_rows(u_mat);
    /* Simple init for centroids - Only for C interface */
    centeroids = centroids_init_simple(u_mat, k);
    /* Step 6 K-Means Algorithm */
    k_means(centeroids, u_mat, k);

    free_mat(u_mat);
    return centeroids;
}
