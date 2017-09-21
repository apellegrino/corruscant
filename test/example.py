from corruscant import twopoint
from corruscant import clustering

import numpy as np

# data set sizes
dsize = 4000
rsize = dsize*20

# generate random points in a cube
np.random.seed(3)
X_data = np.random.rand(dsize, 3)
X_rand = np.random.rand(rsize, 3)

# choose radius bin edges
radii = np.logspace(-2.5,-1.,6)

# give each data point and random point a field ID
# we will split up the data into fourths in the x dimension
data_fields = (X_data[:,0] * 4).astype('int32')
rand_fields = (X_rand[:,0] * 4).astype('int32')

# generate K-d trees
dtree = clustering.tree(X_data, data_fields)
rtree = clustering.tree(X_rand, rand_fields)

# get the correlation function results
results = twopoint.threedim.autocorr(dtree, rtree, radii,
                                    est_type="landy-szalay", num_threads=4)

print(results.ftf_error())
print(results)
