import tpcf
import numpy as np

dsize = 4000
rsize = dsize*20

np.random.seed(3)
X_data = np.random.rand(dsize, 3)
X_random = np.random.rand(rsize, 3)

radii = np.logspace(-2.5,-1.,6)

results = tpcf.pair_counts(X_data, X_random, radii,
                           xi_estimator_type="landy-szalay",
                           xi_error_type="field-to-field", num_threads=2)

print "------- Using Landy-Szalay estimator, field-to-field error -------"

print "r = "
print radii

print "DD, DR, RR differential counts"
dd = results["DD"]
dr = results["DR"]
rr = results["RR"]
print np.vstack([dd,dr,rr]).T

print "Estimated Xi(r)"
print results["estimator"]
print "sigma^2 of Xi(r)"
print results["error"]
