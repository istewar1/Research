import numpy as np
from numpy.random import uniform,randn
import matplotlib.pyplot as plt

def create_uniform_particles(x_range, y_range, N):
    particles = np.empty((N, 3))
    particles[:, 0] = uniform(x_range[0], x_range[1], size=N)
    particles[:, 1] = uniform(y_range[0], y_range[1], size=N)
    return particles

def create_gaussian_particles(mean, std, N):
    particles = np.empty((N, 3))
    particles[:, 0] = mean[0] + (randn(N) * std[0])
    particles[:, 1] = mean[1] + (randn(N) * std[1])
    return particles

p = create_uniform_particles([0,1],[0,1],10000)
plt.scatter(p[:,0],p[:,1],c='k');plt.show()