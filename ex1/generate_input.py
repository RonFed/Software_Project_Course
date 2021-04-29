import random
from numpy import random
import sys

if len(sys.argv) >=2 :
    k = int(sys.argv[1])
else:
    k = 3
d = 5
n = 10000
mean = [10 * i for i in range(k)]
std_div = 5

for i in range(n):
    x = random.normal(loc = random.choice(mean), scale = std_div, size=(d))
    for j in range(d-1):
        print("{:.4f}".format(x[j]), end=",")
    print("{:.4f}".format(x[d-1]))
