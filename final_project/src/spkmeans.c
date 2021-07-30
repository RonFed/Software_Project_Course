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
        for (j = 0; j < mat->cols - 1; j++)
        {
            printf("%.4f,", (mat->data)[i][j]);
        }
        printf("%.4f", (mat->data)[i][mat->cols - 1]);
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
        CALLOC_ARR_ASSERT((mat->data)[i], double, dim - i);
    }
    return mat;
}

void print_sym_mat(sym_matrix *mat)
{
    int i, j;
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
        result += (dif * dif);
    }
    result = sqrt(result);
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
    assert(row < mat->dim);
    double sum = 0;
    int j;
    for (j = 0; j < mat->dim; j++)
    {
        sum += get_val_sym(mat, row, j);
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
            curr_weight = get_val_sym(weight_m, i, j);
            if (j == i)
            {
                // on the diagonal w(i,i) = 0 => l_norm(i,i) = 1
                set_val_sym(l_norm_m, i, j, 1);
            }
            else
            {
                d_sqrt = 1 / (sqrt(curr_d_i * curr_d_j));
                set_val_sym(l_norm_m, i, j, -d_sqrt * curr_weight);
            }
        }
    }

    return l_norm_m;
}

// ----------------------------------------------------------------------------------//
///////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  Jacobi Algorithm   //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------------------//

void init_e_mat(jacobi_matrix *j_mat)
{
    int i, n;
    n = (j_mat->mat)->dim;
    CALLOC_ARR_ASSERT(j_mat->e_mat, e_vector, n);
    for (i = 0; i < n; i++)
    {
        CALLOC_ARR_ASSERT((j_mat->e_mat)[i].vec, double, n);
        // setting the index of eigan-vector
        (j_mat->e_mat)[i].ind = i;
        // initialize the matrix to the identity
        ((j_mat->e_mat)[i].vec)[i] = 1;
    }
}

jacobi_matrix *init_jac_mat(sym_matrix *s_mat)
{
    jacobi_matrix *jac_mat = (jacobi_matrix *)malloc(sizeof(jacobi_matrix));
    jac_mat->mat = s_mat;
    init_e_mat(jac_mat);
    // calloc for max inds array
    CALLOC_ARR_ASSERT(jac_mat->max_inds, unsigned int, s_mat->dim);
    jac_mat->off_diff = 0;
    return jac_mat;
}

void max_abs_val_initial(jacobi_matrix *j_mat)
{
    double max_total = -1, max_in_row, curr_abs;
    unsigned int i, j;
    sym_matrix *sm = j_mat->mat;
    unsigned int max_row;
    for (i = 0; i < sm->dim; i++)
    {
        // for each row find the max element
        max_in_row = -1;
        // max abs value off - diagonal
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
    // symmetric matrix so we can change the indexes here
    // since the implementation of symmetric matrix used top triangle representation
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

// this fuction is called after the max value row and col indexes are updated
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

void update_e_mat(jacobi_matrix *j_mat)
{
    e_vector *e_mat_data = j_mat->e_mat;
    unsigned int i, j, r;
    double c, s, curr_i, curr_j;
    e_vector vec_i, vec_j;
    i = j_mat->row_ind;
    j = j_mat->col_ind;
    c = j_mat->c;
    s = j_mat->s;
    // update only the i-th and j-th eigan-vectors
    vec_i = e_mat_data[i];
    vec_j = e_mat_data[j];
    for (r = 0; r < (j_mat->mat->dim); r++)
    {
        // the update according to the rotation matrix
        // update row i - the i-th and j-th eigan-vectors
        curr_i = (vec_i.vec)[r];
        curr_j = (vec_j.vec)[r];
        (vec_i.vec)[r] = (c * curr_i) - (s * curr_j);
        (vec_j.vec)[r] = (s * curr_i) + (c * curr_j);
    }
}

void print_e_mat(jacobi_matrix *j_mat)
{
    e_vector *e_mat_data = j_mat->e_mat;
    int n = (j_mat->mat)->dim;
    int i, j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n - 1; j++)
        {
            printf("%.4f,", (e_mat_data[j].vec)[i]);
        }
        printf("%.4f", (e_mat_data[n - 1].vec)[i]);
        printf("\n");
    }
}

void rotate_jac(jacobi_matrix *j_mat)
{
    sym_matrix *sm_data = j_mat->mat;
    unsigned int i, j, r;
    double curr_i, curr_j, new_i, new_j, new_ii, new_jj, c, s;
    double pos_diff = 0, neg_diff = 0;
    i = j_mat->row_ind;
    j = j_mat->col_ind;
    c = j_mat->c;
    s = j_mat->s;
    // update the i-th and j-th columns
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

            // updating the diffrence in off diagonal elements field
            pos_diff += (SQUARE(curr_i) + SQUARE(curr_j));
            neg_diff += (SQUARE(new_i) + SQUARE(new_j));
        }
    }
    // update relevant diagonal elements
    curr_i = get_val(j_mat, i, i);
    curr_j = get_val(j_mat, j, j);
    new_ii = (c * c * curr_i) + (s * s * curr_j) - (2 * s * c * (get_val(j_mat, i, j)));
    new_jj = (s * s * curr_i) + (c * c * curr_j) + (2 * s * c * (get_val(j_mat, i, j)));
    set_val(j_mat, i, i, new_ii);
    set_val(j_mat, j, j, new_jj);
    pos_diff += SQUARE(get_val(j_mat, i, j));
    set_val(j_mat, i, j, 0);
    // multiply by 2 since the matrix is symmetric
    j_mat->off_diff = 2 * (pos_diff - neg_diff);
}

void update_max_jac(jacobi_matrix *j_mat, int is_not_first)
{
    if (is_not_first)
    {
        // finding the new max off-diagonal abs value
        unsigned int r, curr_row_max_ind, curr_col_max_ind;
        curr_row_max_ind = j_mat->row_ind;
        curr_col_max_ind = j_mat->col_ind;
        // first check the rows being rotated
        update_max_in_row(j_mat, curr_row_max_ind);
        update_max_in_row(j_mat, curr_col_max_ind);
        unsigned int *max_arr = j_mat->max_inds;
        for (r = 1; r < (j_mat->mat)->dim; r++)
        {
            // all rows excpet the two that were updated
            if ((r != curr_row_max_ind) && (r != curr_col_max_ind))
            {
                // the max in the row r is in a cloumn that wasn't changed in the last rotation
                if ((max_arr[r] != curr_row_max_ind) && (max_arr[r] != curr_col_max_ind))
                {
                    // only two columns have changed in this row to check
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
                    // the column in row r containing the max has been changed in the rotation
                    // find the max in row r
                    update_max_in_row(j_mat, r);
                }
            }
        }
        update_total_max_inds(j_mat);
    }
    else
    {
        // only occurs once in first iteration of jacobi
        max_abs_val_initial(j_mat);
    }
}

void update_max_in_row(jacobi_matrix *j_mat, unsigned int row)
{
    unsigned int j, max_ind, old_max_col;
    double curr, max;
    old_max_col = (j_mat->max_inds)[row];
    max = fabs(get_val(j_mat, row, old_max_col));
    // only checking elements right to the diagonal
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

void update_total_max_inds(jacobi_matrix *j_mat)
{
    double max = -1, curr;
    unsigned int r, max_row, max_col;
    // last row doesn't contain relevant elements (just one element on diagonal)
    for (r = 0; r < (j_mat->mat)->dim - 1; r++)
    {
        // max value in line r
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

int cmp_vecs(const void *vec1, const void *vec2)
{
    // if the two vectors have the same eigan-value, sort them by index (the smaller index first)
    // else sort them by the eigan-value - smaller one comes firdt
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
    int n, i;
    n = (j_mat->mat)->dim;
    for (i = 0; i < n; i++)
    {
        (j_mat->e_mat)[i].e_val = get_val(j_mat, i, i);
    }
}

// main Jacobi algorithm
void jacobi(jacobi_matrix *j_mat)
{
    int iters = 0;
    do
    {
        update_max_jac(j_mat, iters);
        update_c_s_params(j_mat);
        rotate_jac(j_mat);
        update_e_mat(j_mat);
        iters++;
    } while (iters < 100 && j_mat->off_diff > EPSILON);
    // setting each eigan-value to the coresponding eigan-vector
    set_eigan_values(j_mat);
    // sorting the eigan-vectors
    qsort(j_mat->e_mat, (j_mat->mat)->dim, sizeof(e_vector), cmp_vecs);
}

// ----------------------------------------------------------------------------------//
///////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  Text Parsing   //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------------------//

matrix *read_file_to_mat(FILE *file_pointer)
{
    matrix *data_mat;
    double **data;
    double current_element;
    char current_delimeter;
    unsigned int lines_count = 0;
    unsigned int current_arr_len = BASE_ARR_SIZE;
    unsigned int dimension = 0;

    MALLOC_ARR_ASSERT(data, double *, BASE_ARR_SIZE);
    // first line analysis to determine the dimension of the csv format (columns number)
    dimension = find_dimension_from_first_line(file_pointer, data[0]);
    lines_count++;

    // TODO !!

    data_mat->cols = dimension;
    data_mat->rows = lines_count;
    data_mat->data = data;
    return data_mat;
}

unsigned int find_dimension_from_first_line(FILE * file_pointer, double * first_line) {
    double current_element;
    char current_delimeter;
    unsigned int dimension = 0;
    unsigned int current_arr_len = BASE_ARR_SIZE;

    CALLOC_ARR_ASSERT(first_line, double, BASE_ARR_SIZE);

    while (fscanf(file_pointer, "%lf%c", &current_element, &current_delimeter) == 2)
    {
        // rarely happens - first line is over BASE_ARR_SIZE elements
        if (dimension == current_arr_len)
        {
            current_arr_len *= ARR_SIZE_MULTIPLY;
            first_line = realloc(first_line, current_arr_len * sizeof(double));
            assert(first_line);
        }
        first_line[dimension] = current_element;
        dimension++;
        if (current_delimeter == '\n'){
            break;
        }
    }
    first_line = realloc(first_line, dimension * sizeof(double));
    return dimension;
}
