import numpy as np
import sys
import spkmeans as sp

def print_list(lst):
    for i in range(len(lst) - 1):
        print("{:.4f}".format(lst[i]), end=",")
    print("{:.4f}".format(lst[len(lst) - 1]))

def print_mat(mat):
    for line in mat:
        print_list(line)

def print_jacobi(j_mat):
    # print the first line of eiganvalues
    print_list(j_mat[0])
    # print the eiganvectors matrix
    print_mat(j_mat[1])


def main():
    file_path = sys.argv[1]
    data = np.genfromtxt(file_path, delimiter=',')
    weights = sp.weights_mat(data.tolist())
    print_mat(weights)
    degree_1d = sp.degree_mat(data.tolist())
    degree_2d = [[degree_1d[j] if i==j else 0 for j in range(len(degree_1d))] for i in range(len(degree_1d))]
    print_mat(degree_2d)
    lnorm = sp.l_norm_mat(data.tolist())
    print_mat(lnorm)
    jacobi_m = sp.jacobi_mat(data.tolist())
    print_jacobi(jacobi_m)

if __name__ == '__main__':
    main()
