from os import error
import numpy as np
import sys
import spkmeans as sp
from enum import Enum

invalid_input_msg = "Invalid Input!"
class Goal(Enum):
    SPK = "spk"
    WAM = "wam"
    DDG = "ddg"
    LNORM = "lnorm"
    JACOBI = "jacobi"

def myAssert(cond, print_msg):
    if not cond:
        print(print_msg)
        exit()

def isLegalGoal(x):
    try:
        Goal(x)
        return True
    except Exception:
        return False

def isInt(x):
    try:
        int(x)
        return True
    except ValueError:
        return False

# check command line arguments and return them if they are valid
# k - number of clusters
# goal - enum string
# file_path path to csv or txt file containing the data
def read_args(argv):
    num_args = len(argv)
    myAssert(num_args == 4, invalid_input_msg)
    k = argv[1]
    myAssert(isInt(k), invalid_input_msg)
    k = int(k)
    myAssert(k >= 0, invalid_input_msg)
    goal = argv[2]
    myAssert(isLegalGoal(goal), invalid_input_msg)
    file_path = argv[3]
   
    return k, goal, file_path

# Avoiding "-0.0000" in the output printouts 
# according to the project forum 
def negative_zero_fix(num):
    if num >= 0 or num <= -5e-5:
        return num
    return -num

# print list in csv format with 4 decimal places
def print_list(lst, last_line = False):
    for i in range(len(lst) - 1):
        print("{:.4f}".format(negative_zero_fix(lst[i])), end=",")
    if last_line:
         end_char = "" 
    else: 
        end_char = "\n"
    print("{:.4f}".format(negative_zero_fix(lst[len(lst) - 1])), end=end_char)

# print matrix (list of lists) in csv format
def print_mat(mat):
    for i in range(len(mat) - 1):
        print_list(mat[i])
    print_list(mat[len(mat) - 1], last_line=True)

def print_jacobi(j_mat):
    # print the first line of eiganvalues
    print_list(j_mat[0])
    # print the eiganvectors matrix
    print_mat(j_mat[1])

# k-means++ algorithm
#   INPUT: list of lists representing tha data
#          k - number of clusters
#   OUTOUT: tuple of 2 elements : first is a matrix containing the inital centroids
#          second is a list containing the initial centroids indices
def k_means_pp_alg(data, k):
    data = np.array(data)
    n, d = data.shape
    # setting seed to be consistent 
    np.random.seed(0)
    centers = np.full((n, d), 0.0, dtype=float)
    centers_ids = np.full(k, -1)
    # choosing the first centroid randomly
    centers_ids[0] = np.random.choice(n)
    centers[0] = data[centers_ids[0]]
    z = 1
    # distances array to hold distances between data vectors and centroids
    Distances = np.full(n, 'Inf', dtype=float)
    while z < k:
        for i in range(n):
            # dynamic programming, reducing complexity
            Distances[i] = min([sum((data[i] - centers[z-1]) ** 2), Distances[i]])
        z += 1
        sum_dists = sum(Distances)
        prob = Distances / sum_dists
        centers_ids[z - 1] = np.random.choice(n, p=prob)
        centers[z - 1] = data[centers_ids[z - 1]]
    return centers[:k], centers_ids

# spectral clustrering goal has been requested
def handle_spk(data, k):
    # get k and T matrix from C API 
    k_T_lst = sp.k_and_T_mat(data, k)
    k = k_T_lst[0]
    T_mat = k_T_lst[1]
    # use k-means++ algorithm to determine initial centroids
    centers, centers_inds = k_means_pp_alg(T_mat, k)
    # call C API to cluster the T matrix using k-means (with initial centroids)
    final_centers = sp.kmeans_from_centroids(T_mat, centers.tolist())
    # print as first line the centroids indices that were taken
    print(*centers_inds, sep=",")
    # print the final centroids matrix
    print_mat(final_centers)

def handle_goal(goal, data, k):
    if (goal == Goal.SPK.value):
        handle_spk(data, k)
    elif (goal == Goal.WAM.value):
        weights = sp.weights_mat(data)
        print_mat(weights)
    elif (goal == Goal.DDG.value):
        degree_1d = sp.degree_mat(data)
        # converting the 1d representation of 
        # diagonal matrix to 2d representation for printing
        degree_2d = [[degree_1d[j] if i==j else 0 for j in range(len(degree_1d))] for i in range(len(degree_1d))]
        print_mat(degree_2d)
    elif (goal == Goal.LNORM.value):
        lnorm = sp.l_norm_mat(data)
        print_mat(lnorm)
    elif (goal == Goal.JACOBI.value):
        jacobi_m = sp.jacobi_mat(data)
        print_jacobi(jacobi_m)
    else :
        myAssert(False, invalid_input_msg)

def main():
    # check validity of command line arguments
    k, goal, file_path = read_args(sys.argv)

    # read text to python's list of lists
    try:
        data = np.genfromtxt(file_path, delimiter=',').tolist()
    except OSError:
        myAssert(False, invalid_input_msg)

    # validate k if it is given (positive)
    if k > 0 and goal == Goal.SPK.value:
        myAssert(k < len(data), invalid_input_msg)

    handle_goal(goal, data, k)

if __name__ == '__main__':
    main()
