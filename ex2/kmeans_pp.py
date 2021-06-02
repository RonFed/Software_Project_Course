import numpy as np
import sys
import pandas as pd
import mykmeanssp as km


def k_means_pp_alg(data, k):
    n, d = data.shape
    np.random.seed(0)
    centers = np.full((n, d), 0.0, dtype=float)
    centers_ids = np.full(k, -1)
    centers_ids[0] = np.random.choice(n)
    centers[0] = data[centers_ids[0]]
    z = 1
    while z < k:
        Distances = np.full(n, 0.0, dtype=float)
        for i in range(n):
            Distances[i] = min([sum((data[i] - centers[j]) ** 2) for j in range(z)])
        z += 1
        sum_dists = sum(Distances)
        prob = Distances / sum_dists
        centers_ids[z - 1] = np.random.choice(n, p=prob)
        centers[z - 1] = data[centers_ids[z - 1]]
    return centers[:k], centers_ids


def read_args(argv):
    num_args = len(argv)
    if num_args == 1:
        print("No arguments provided, exiting...")
        sys.exit(1)
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

    # input arguments validation
    k, max_iter, path1, path2 = read_args(sys.argv)

    # getting the data from files and merging
    data1 = pd.read_csv(path1, index_col=0, header=None).sort_index()
    data2 = pd.read_csv(path2, index_col=0, header=None).sort_index()
    data_total = pd.merge(data1, data2, left_index=True, right_index=True).values

    # calculation of initial centroids
    centers, centers_inds = k_means_pp_alg(data_total, k)
    print(*centers_inds, sep=",")

    # calling the C API to calculate k means using the initial vectors calculated
    final_result = km.fit(centers.tolist(), data_total.tolist(), max_iter)
    for i in range(len(final_result)-1):
        print(*np.round(final_result[i],4), sep=",")
    print(*np.round(final_result[i+1],4), sep=",", end='')

if __name__ == '__main__':
    main()
