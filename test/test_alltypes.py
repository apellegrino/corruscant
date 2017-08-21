from corruscant import twopoint
from corruscant import clustering
import numpy as np
import testing

dsize = 4000
rsize = dsize*20
radii = np.logspace(-2.5,-1.,6)
radii_ang = np.logspace(-1.5,0.5,6)
N_fields = 4
num_threads = 4

np.random.seed(3)

from corruscant.clustering import est_standard, est_hamilton, est_landy_szalay

def print_alltypes(results):
    err_types = [None, "poisson", "ftf", "jackknife", "bootstrap"]
    est_types = [est_standard, est_hamilton, est_landy_szalay]
    est_names = ["Standard", "Hamilton", "Landy-Szalay"]

    for est,est_name in zip(est_types,est_names):
        for err in err_types:
            results.estimator = est
            results.error_type = err
            print("\n")
            print("Estimator type = {:s}, Error type = {:s}".format(est_name, str(err)))
            print(results)

results = twopoint.threedim.autocorr(
                clustering.tree(*testing.random_cube(dsize, N_fields)),
                clustering.tree(*testing.random_cube(rsize, N_fields)),
                radii, num_threads=num_threads
                                     )
print_alltypes(results)

results = twopoint.threedim.crosscorr(
                clustering.tree(*testing.random_cube(dsize, N_fields)),
                clustering.tree(*testing.random_cube(dsize, N_fields)),
                clustering.tree(*testing.random_cube(rsize, N_fields)),
                clustering.tree(*testing.random_cube(rsize, N_fields)),
                radii, num_threads=num_threads
                                     )
print_alltypes(results)

results = twopoint.angular.autocorr(
                clustering.tree(*testing.random_sphere(dsize, N_fields)),
                clustering.tree(*testing.random_sphere(rsize, N_fields)),
                radii_ang, num_threads=num_threads
                                     )
print_alltypes(results)

results = twopoint.angular.crosscorr(
                clustering.tree(*testing.random_sphere(dsize, N_fields)),
                clustering.tree(*testing.random_sphere(dsize, N_fields)),
                clustering.tree(*testing.random_sphere(rsize, N_fields)),
                clustering.tree(*testing.random_sphere(rsize, N_fields)),
                radii_ang, num_threads=num_threads
                                     )
print_alltypes(results)
