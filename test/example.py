import tpcf
import numpy as np

dsize = 4000
rsize = dsize*20

np.random.seed(3)
X_data = np.random.rand(dsize, 3)
X_random = np.random.rand(rsize, 3)

radii = np.logspace(-2.5,-1.,6)

data_fields = np.where(X_data.T[0] < 0.5, 1, 2)
rand_fields = np.where(X_random.T[0] < 0.5, 1, 2)

dtree = tpcf.tree(X_data, data_fields, 2)
rtree = tpcf.tree(X_random, rand_fields, 2)

results = tpcf.twopoint(dtree, rtree, radii,
                           est_type="landy-szalay",
                           err_type='jackknife', num_threads=2)

print "------- Using Landy-Szalay estimator, jackknife error -------"

print "r = "
print radii

print "DD, DR, RR differential counts"
dd, dr, rr = results.total_pair_counts()
print np.vstack([dd,dr,rr]).T

print "Estimated Xi(r)"
print results.estimate()
print "Error of Xi(r)"
print results.error()
print "Covariance matrix"
print results.covariance()
