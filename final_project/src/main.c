#include <stdio.h>
#include <math.h>
#include "spkmeans.h"

int main(int argc, char const *argv[])
{
    int i,j;
    unsigned int rows = 5;
    unsigned int cols = 5;
    matrix * m = init_mat(rows, cols);
     for (i = 0; i < m->rows; i++)
    {
        for (j = 0; j < m->cols; j++) 
        {
           (m->data)[i][j] = i + j;
        }
    }
    print_mat(m);

    // sym_matrix * sm;
    // (sm->data)[1][0] = 3;
    // print_sym_mat(sm);

    // double x = l2_norm_vectors(m, 1, 2);
    // printf("%f\n", x);
    // x = exp_norm_vectors(m, 1, 2);
    // printf("%.10f\n", x);

    // sm = weights_mat(m);
    // print_sym_mat(sm);
    // x = row_sum_sym_mat(sm,1);
    // printf("%.10f\n", x);

    // diag_matrix * dm;
    // dm = degree_mat(m);
    // // print_diag_mat(dm);

    sym_matrix * l_norm;
    l_norm = l_norm_mat(m);
    print_sym_mat(l_norm);
    return 0;
}
