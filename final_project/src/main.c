#include <stdio.h>
#include <math.h>
#include "spkmeans.h"

int main(int argc, char const *argv[])
{
    FILE * fp;
    char * filename = "final_project\\src\\python_scripts\\data1_data.csv";
    fp = fopen(filename, "r");
    matrix * mm = read_file_to_mat(fp);
    print_mat(mm);

    // int i, j;
    // unsigned int rows = 10;
    // unsigned int cols = 5;
    // matrix *m = init_mat(rows, cols);
    // for (i = 0; i < m->rows; i++)
    // {
    //     for (j = 0; j < m->cols; j++)
    //     {
    //         (m->data)[i][j] = i + j*j -i*i + j*i;
    //     }
    // }
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
    // // print_sym_mat(l_norm);
    // jacobi_matrix * jm = init_jac_mat(l_norm);
    // jacobi(jm);

    // sym_matrix * rami;
    // rami = init_sym_mat(3);
    // set_val_sym(rami,0,0,3);
    // set_val_sym(rami,1,0,2);
    // set_val_sym(rami,2,0,4);
    // set_val_sym(rami,2,1,2);
    // set_val_sym(rami,2,2,3);
    // jacobi_matrix * jm = init_jac_mat(rami);
    // jacobi(jm);

    // python example
    // sym_matrix *rami;
    // rami = init_sym_mat(4);
    // set_val_sym(rami, 0, 0, 0.91);
    // set_val_sym(rami, 1, 0, 0.18);
    // set_val_sym(rami, 2, 0, 0.48);
    // set_val_sym(rami, 3, 0, 0.99);
    // set_val_sym(rami, 1, 1, 0.63);
    // set_val_sym(rami, 2, 1, 0.18);
    // set_val_sym(rami, 3, 1, 0.1);
    // set_val_sym(rami, 2, 2, 0.22);
    // set_val_sym(rami, 3, 2, 0.13);
    // set_val_sym(rami, 3, 3, 0.48);
    // print_sym_mat(rami);
    // jacobi_matrix *jm = init_jac_mat(rami);
    // jacobi(jm);
    // print_sym_mat(jm->mat);
    // printf("-------------e_mat------\n");
    // print_e_mat(jm);

    // sym_matrix * sm2 = init_sym_mat(rows);
    //  for (i = 0; i < m->rows; i++)
    // {
    //     for (j = 0; j <=i; j++)
    //     {
    //        set_val_sym(sm2, i, j, i+j);
    //     }
    // }
    // print_sym_mat(sm2);
    // jacobi_matrix * jm2 = init_jac_mat(sm2);
    // jacobi(jm2);
    return 0;
}
