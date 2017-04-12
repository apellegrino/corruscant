import numpy as np
import tpcf
from matplotlib import pyplot as plt

dsize = 2000
rsize = dsize*20

#X_data = np.vstack([np.random.rand(dsize/2,3),np.random.rand(dsize/2,3)*.25])
np.random.seed(3)
X_data = np.random.rand(dsize,3).T
X_random = np.random.rand(rsize,3).T

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
