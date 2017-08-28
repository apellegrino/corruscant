import timeit
import numpy as np

setup = """
from corruscant import twopoint
from corruscant import clustering

import numpy as np

# data set sizes
dsize = 20000
rsize = dsize*20

# generate random points in a cube
np.random.seed(3)
X_data = np.random.rand(dsize, 3)
X_rand = np.random.rand(rsize, 3)

# choose radius bin edges
radii = np.logspace(-2.5,-1.,6)

# give each data point and random point a field ID
# we will split up the data into fourths in the x dimension
data_fields = (X_data[:,0] * {N_fields}).astype('int32')
rand_fields = (X_rand[:,0] * {N_fields}).astype('int32')

# generate K-d trees
dtree = clustering.tree(X_data, data_fields)
rtree = clustering.tree(X_rand, rand_fields)
"""

stmt = "twopoint.threedim.autocorr(dtree, rtree, radii, num_threads=4)"

tot = []
for N_fields in range(10,210,10):
    times=timeit.repeat(stmt, setup.format(N_fields=N_fields), repeat=5, number=1)
    print times
    tot.append(times)

with open("results.csv", 'w') as f:
    f.write("n_fields,avg(sec),stddev(sec)\n")
    for n, times in zip(range(10,210,10), tot):
        f.write("{:d},{:f},{:f}\n".format(n, np.average(times), np.std(times)))
