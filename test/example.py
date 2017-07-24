import tpcf
import numpy as np

# data set sizes
dsize = 4000
rsize = dsize*20

# generate random points in a cube
np.random.seed(3)
X_data = np.random.rand(dsize, 3)
X_random = np.random.rand(rsize, 3)

# choose radius bin edges
radii = np.logspace(-2.5,-1.,6)

# give each data point and random point a field ID
# for a simple example, call the x < 0.5 region `field 1` and the x > 0.5
# region `field 2`
data_fields = np.where(X_data[:,0] < 0.5, 1, 2)
rand_fields = np.where(X_random[:,0] < 0.5, 1, 2)

# generate K-d trees
dtree = tpcf.tree(X_data, data_fields)
rtree = tpcf.tree(X_random, rand_fields)

# get the correlation function results
results = tpcf.twopoint(dtree, rtree, radii,
                           est_type="landy-szalay",
                           err_type='jackknife', num_threads=4)

print(results)
