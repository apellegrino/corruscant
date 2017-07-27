import twopoint
import numpy as np

# data set sizes
dsize = 10000
rsize = dsize*20

# generate uniform sphere points in RA, DEC
np.random.seed(3)

data_ra = np.random.rand(dsize) * 360.
data_dec = np.arcsin(np.random.rand(dsize) * 2. - 1.) * 180. / np.pi

rand_ra = np.random.rand(rsize) * 360.
rand_dec = np.arcsin(np.random.rand(rsize) * 2. - 1.) * 180. / np.pi

# choose radius bin edges
radii_deg = np.logspace(-2.,0.,6)

# give each data point and random point a field ID
# for a simple example, call the x < 0.5 region `field 1` and the x > 0.5
# region `field 2`
data_fields = np.where(data_ra < 180., 1, 2)
rand_fields = np.where(rand_ra < 180., 1, 2)

# transform to cartesian coords
X_data = twopoint.angular.cartesian(data_ra, data_dec)
X_rand = twopoint.angular.cartesian(rand_ra, rand_dec)

# generate K-d trees
dtree = twopoint.clustering.tree(X_data, data_fields)
rtree = twopoint.clustering.tree(X_rand, rand_fields)

# get the correlation function results
results = twopoint.angular.autocorr(dtree, rtree, radii_deg,
                           est_type="landy-szalay",
                           err_type='jackknife', num_threads=4)

print(results)
