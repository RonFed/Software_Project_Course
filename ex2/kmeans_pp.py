import numpy as np
import sys
import pandas as pd


def k_means_pp_alg(data, k):
    n, d = data.shape
    centers = np.full((n, d), 0.0, dtype = float)
    rand_index = np.random.choice(n)
    centers[0] = data[rand_index]
    z = 1
    while z < k:
        Distances = np.full(n, 0.0, dtype=float)
        for i in range(n):
            Distances[i] = min([sum((data[i] - centers[j]) ** 2) for j in range(z)])
        z += 1
        sum_dists = sum(Distances)
        prob = Distances / sum_dists
        prob_index = np.random.choice(n, p=prob)
        centers[z - 1] = data[prob_index]
    return centers[:k]


def read_args(argv):
    num_args = len(argv)
    try:
        k = int(argv[1])
    except ValueError:
        print("Invalid value for k")
        sys.exit(1)
    if num_args == 4:
        max_iter = 300
        path1 = argv[2]
        path2 = argv[3]
    elif num_args == 5:
        try:
            max_iter = int(argv[2])
        except ValueError:
            print("Invalid value for max_iter")
            sys.exit(1)
        path1 = argv[3]
        path2 = argv[4]
    else:
        print("Number of variables is not valid")
        sys.exit(1)
    return k, max_iter, path1, path2


def main():
    k, max_iter, path1, path2 = read_args(sys.argv)
    data1 = pd.read_csv(path1, index_col=0, header=None)
    data2 = pd.read_csv(path2, index_col=0, header=None)
    d = pd.merge(data1, data2, left_index=True, right_index=True).values
    centers = k_means_pp_alg(d, k)
    print(centers)


if __name__ == '__main__':
    main()

