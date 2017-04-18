import numpy as np
import tpcf
from matplotlib import pyplot as plt

def uniform_sphere(N, dim):
    data = []
    center = np.array([.5,.5,.5])

    while N > 0:
        point = np.random.rand(3)
        norm2 = np.sum(np.multiply(point-center,point-center))
        if norm2 < .25:
            data.append(point)
            N -= 1

    return np.array(data).T
        
dsize = 1000
rsize = dsize*20

np.random.seed(3)
X_data = uniform_sphere(dsize,3)
X_random = uniform_sphere(rsize,3)

N = 10
gen = tpcf._jackknife(X_data, X_random, N=N)

fig = plt.figure()
plt.title("Jackknife error data subset selection")
for i, (data, rand) in enumerate(gen):
    ax = plt.subplot(2,(N/2),i+1)
    ax.plot(rand[0],rand[1], '.', markersize=1, label='Random')
    ax.plot(data[0],data[1], 'r.', markersize=3, label='Data')
    ax.set_title("i = %d/%d" % (i,N))
    ax.set_xlim([0,1])
    ax.set_ylim([0,1])

plt.show()
