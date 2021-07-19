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
    // print_mat(m);

    // sym_matrix * sm = init_sym_mat(4);
    // set_val_sym(sm, 1,2,1.21);
    // set_val_sym(sm, 3,0,3.14);
    // set_val_sym(sm, 3,3,22);
    // print_sym_mat(sm);
    // printf("%.4f\n", get_val_sym(sm,2,1));
    // printf("%.4f\n", get_val_sym(sm,0,3));
    // double x = l2_norm_vectors(m, 1, 2);
    // printf("%f\n", x);
    // x = exp_norm_vectors(m, 1, 2);
    // printf("%.10f\n", x);

    // sym_matrix * sm = weights_mat(m);
    // print_sym_mat(sm);
    // double x = row_sum_sym_mat(sm,1);
    // printf("%.10f\n", x);

    // diag_matrix * dm;
    // dm = degree_mat(m);
    // // print_diag_mat(dm);

    // sym_matrix * l_norm;
    // l_norm = l_norm_mat(m);
    // print_sym_mat(l_norm);
    sym_matrix * rami;
    rami = init_sym_mat(3);
    set_val_sym(rami,0,0,3);
    set_val_sym(rami,1,0,2);
    set_val_sym(rami,2,0,4);
    set_val_sym(rami,2,1,2);
    set_val_sym(rami,2,2,3);
    // (rami->data)[0][0] = 3;
    // (rami->data)[1][0] = 2;
    // (rami->data)[2][0] = 4;
    // (rami->data)[2][1] = 2;
    // (rami->data)[2][2] = 3;
    // print_sym_mat(rami);
    jacobi_matrix * jm = init_jac_mat(rami);
    // printf("%.15f\n", l_norm[2,4]);
    jacobi(jm);
    return 0;
}
