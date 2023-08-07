# %%
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
# from sklearn.cluster import KMeans
from matplotlib.colors import ListedColormap


CMAP = ListedColormap(['#2342cd', '#df4951', '#f4c82a'])
EDGECOLOR = 'k'
N_FRAMES = 10

# Generate some random data
np.random.seed(43)
X = np.random.randn(300, 2)
# add random x and y offset to each 100 points
# offsets = np.random.uniform(-4, 4, size=(3, 2))
# X[:100] += offsets[0]
# X[100:200] += offsets[1]
# X[200:] += offsets[2]


# Initialize the k-means model
class KMeans:
    def __init__(self, data, n_clusters):
        self.data = data
        self.n_clusters = n_clusters
        self.iteration = 0
        # pick random points as initial centroids
        self.centers = [data[np.random.choice(range(len(data)), n_clusters, replace=False)]]
        self.labels = np.zeros(len(data))
        # self.update_labels()
        self.fresh = True

    def update_labels(self):
        for i, x in enumerate(self.data):
            self.labels[i] = np.argmin(np.linalg.norm(self.centers[-1] - x, axis=1))
    
    def update_centers(self):
        tmp_centers = np.zeros((self.n_clusters, 2))
        for i in range(self.n_clusters):
            tmp_centers[i] = np.mean(self.data[self.labels == i], axis=0)
        self.centers.append(tmp_centers)

    def fit(self):
        if self.fresh:
            self.update_labels()
            self.fresh = False
        else:
            self.update_labels()
            self.update_centers()

# %%
kmeans = KMeans(X, 3)

# plt.scatter(X[:, 0], X[:, 1], c=kmeans.labels)
# plt.scatter(kmeans.centers[:, 0], kmeans.centers[:, 1], marker='x', s=200, linewidths=3, color='r')
# plt.title('Iteration 0')


# Define the update function for the animation
def update(frame):
    # Fit the k-means model to the data
    kmeans.fit()

    # Plot the data points and cluster centers
    plt.clf()
    plt.scatter(X[:, 0], X[:, 1], c=kmeans.labels, cmap=CMAP)
    plt.scatter(kmeans.centers[-1][:, 0], kmeans.centers[-1][:, 1], marker='s', s=200, linewidths=3, c=range(kmeans.n_clusters), edgecolors=EDGECOLOR, cmap=CMAP)
    for i in range(kmeans.n_clusters):
        center_x = [x[i, 0] for x in kmeans.centers]
        center_y = [x[i, 1] for x in kmeans.centers]
        center_xy = list(zip(center_x, center_y))
        # plt.plot(center_x, center_y, c=EDGECOLOR)
        for i, (start, stop) in enumerate(zip(center_xy[:-1], center_xy[1:])):
            plt.plot([start[0], stop[0]], [start[1], stop[1]], c=EDGECOLOR, alpha=((i+1)/len(center_xy)))
    plt.title(f'Iteration {frame+1}')

# Create the animation
anim = FuncAnimation(plt.gcf(), update, frames=N_FRAMES, repeat=True)

# save the animation as a good quality gif
anim.save('kmeans.gif', writer='pillow', fps=1, dpi=600)

# Show the animation
plt.show()

# %%
