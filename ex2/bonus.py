import matplotlib.pyplot as plt
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib.patches import Ellipse
from sklearn import datasets, cluster


def differ_arr(x):
    return [(x[i + 1] - x[i]) / (x[i] - x[i - 1]) for i in range(1, len(x) - 1)]

iris = datasets.load_iris()

inertia = np.full(10, 0)
k = np.arange(1, 11)
for j in k:
    inertia[j - 1] = cluster.KMeans(n_clusters=j, random_state=0, init='k-means++').fit(iris.data).inertia_
ax = plt.subplot(1,1,1)
plt.plot(k, inertia)
plt.title('Elbow method for optimal K in k-means algorithm', fontsize=12)
plt.xlabel('k - number of clusters')
plt.ylabel('Inertia')

differs = differ_arr(inertia)

elbow_k = differs.index(min(differs)) + 1


ax.add_artist(Ellipse((elbow_k + 1, inertia[elbow_k]), width=0.5, height=50, facecolor='none', edgecolor='b', lw=1))

plt.savefig('elbow.png')

