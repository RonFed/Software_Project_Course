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


def main():
    file_path = sys.argv[1]
    data = np.genfromtxt(file_path, delimiter=',')
    weights = sp.weights_mat(data.tolist())
    print_mat(weights)


if __name__ == '__main__':
    main()
