import twopoint
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
dtree = twopoint.clustering.tree(X_data, data_fields)
rtree = twopoint.clustering.tree(X_rand, rand_fields)

# get the correlation function results

from twopoint.clustering import est_standard, est_hamilton, est_landy_szalay
err_types = [None, "poisson", "ftf", "jackknife", "bootstrap"]
est_types = [est_standard, est_hamilton, est_landy_szalay]

results = twopoint.threedim.autocorr(dtree, rtree, radii,
                           num_threads=4)

for est in est_types:
    for err in err_types:
        results.estimator = est
        results.error_type = err
        print("\n")
        print("Estimator type = {:s}, Error type = {:s}".format(est, err))
        print(results)
