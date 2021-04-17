import sys


##########################################################
#####################  Helper functions ##################
##########################################################


def split_to_nums(str):
    return tuple(float(x) for x in str.split(','))


def dist(x, y):
    sumo = 0
    for i in range(len(x)):
        sumo += (x[i] - y[i]) ** 2
    return sumo


def add_tuples(x, y):
    return tuple(x[t] + y[t] for t in range(len(x)))


def multiply_scalar(x, scalar):
    return tuple(scalar * x[t] for t in range(len(x)))


def print_tuple(t):
    for i in range(len(t) - 1):
        print("{:.4f}".format(t[i]), end=",")
    print("{:.4f}".format(t[len(t) - 1]))


def closest_cluster_index(x):
    min_dist = float('Inf')
    min_dist_index = -1
    for i in range(len(centers)):
        d = dist(x, centers[i])
        if d < min_dist:
            min_dist = d
            min_dist_index = i
    return min_dist_index


##########################################################
#####################  Text Parsing ######################
##########################################################


n = len(sys.argv)
# testing for correct command line arguments
# k must be provided as second argument
# max_iter is optional as third argument
if n < 2:
    print("K not provided")
    exit()
elif n == 2:
    max_iter = 200
elif n > 3:
    print("Too many arguments")
    exit()
else:
    try:
        max_iter = int(sys.argv[2])
    except ValueError:
        print("Invalid value for max_iter")
        exit()
try:
    k = int(sys.argv[1])
except ValueError:
    print("Invalid value for k")
    exit()

data = [split_to_nums(input())]
d = len(data[0])
while True:
    try:
        data.append(split_to_nums(input()))
    except EOFError:
        break

##########################################################
#####################  K-mean algorithm ##################
##########################################################


centers = data[0:k]
which_cluster_ls = [-1] * len(data)
clusters_sums = [tuple([0] * d)] * k
clusters_size = [0] * k
for i in range(max_iter):
    differ = False
    for j in range(len(data)):
        new_closest_cluster = closest_cluster_index(data[j])
        if new_closest_cluster != which_cluster_ls[j]:
            # updating the closest cluster for data[j]
            if which_cluster_ls[j] != -1:
                # decrease 1 from the old cluster size
                clusters_size[which_cluster_ls[j]] -= 1
                # decrease the amount of data[j] from old cluster sum
                clusters_sums[which_cluster_ls[j]] = add_tuples(clusters_sums[which_cluster_ls[j]],
                                                                multiply_scalar(data[j], -1))
            # add 1 to the new cluster size
            clusters_size[new_closest_cluster] += 1
            # update the sum of the new cluster
            clusters_sums[new_closest_cluster] = add_tuples(clusters_sums[new_closest_cluster], data[j])
            # update the closest cluster field
            which_cluster_ls[j] = new_closest_cluster
            differ = True

    if not differ:
        # converged
        break
    for l in range(k):
        if clusters_size[l] != 0:
            centers[l] = multiply_scalar(clusters_sums[l], 1 / clusters_size[l])

for x in centers:
    print_tuple(x)
