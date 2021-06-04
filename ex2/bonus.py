import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse
from sklearn import datasets, cluster


def differ_arr(x):
    return [(x[i + 1] - x[i]) / (x[i] - x[i - 1]) for i in range(1, len(x) - 1)]

print('a')
iris = datasets.load_iris()

print('b')
inertia = np.full(10, 0)
k = np.arange(1, 11)
for j in k:
    print('c')
    inertia[j - 1] = cluster.KMeans(n_clusters=j, random_state=0, init='k-means++').fit(iris.data).inertia_
print('eeee')
ax = plt.subplot(1,1,1)
print('d')
plt.plot(k, inertia)
print('d1')
plt.title('Elbow method for optimal K in k-means algorithm', fontsize=12)
plt.xlabel('k - number of clusters')
plt.ylabel('Inertia')

differs = differ_arr(inertia)
print('d2')
elbow_k = differs.index(min(differs)) + 1
print('e')

ax.add_artist(Ellipse((elbow_k + 1, inertia[elbow_k]), width=0.5, height=50, facecolor='none', edgecolor='b', lw=1))
print('f')
plt.savefig('elbow.png')

