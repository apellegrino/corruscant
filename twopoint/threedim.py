from ctypes import CDLL, POINTER, c_double, c_int
import os
import numpy as np

from . import clustering

path_here = os.path.abspath(__file__)
path_dir = os.path.dirname(path_here)
coordlib = CDLL("{:s}/../bin/libcoords.so".format(path_dir))

coordlib.radecdist2cart64.restype = None
coordlib.radecdist2cart64.argtypes = [
                                POINTER(c_double), # ra
                                POINTER(c_double), # dec
                                POINTER(c_double), # dist
                                POINTER(c_double), # x,y,z
                                c_int, # N
                                ]

def autocorr(data_tree, rand_tree, radii, est_type="landy-szalay",
             err_type="jackknife", num_threads=4):

    results = clustering._autocorr(data_tree, rand_tree, radii,
                            est_type=est_type, err_type=err_type,
                            num_threads=num_threads)

    results.radii_nominal = results.radii_euclidean

    return results

def cartesian(ra, dec, distance=None, z=None, cosmology=None):
    """Returns array of cartesian points given ra, dec, distance arrays or 
    ra, dec, redshift arrays with an Astropy cosmology. Provided distances
    take precedence over provided redshift and cosmology.
    """

    if distance is None:
        # calculate comoving distance
        distance = cosmology.comoving_distance(z)

    #sph = np.require(
    #                np.vstack([ra, dec, distance]).T,
    #                requirements = 'COA'
    #                 )

    #cartesian = np.empty(sph.shape)

    #for in_pt, out_pt in zip(sph, cartesian):
    #    coordlib.radecdist2cart64(
    #                        in_pt.ctypes.data_as(POINTER(c_double)),
    #                        out_pt.ctypes.data_as(POINTER(c_double)),
    #                        )

    N = ra.size

    ra = np.require(ra, requirements='AC')
    dec = np.require(dec, requirements='AC')
    distance = np.require(distance, requirements='AC')

    cartesian = np.empty((N, 3))

    coordlib.radec2cart64(
                        ra.ctypes.data_as(POINTER(c_double)),
                        dec.ctypes.data_as(POINTER(c_double)),
                        distance.ctypes.data_as(POINTER(c_double)),
                        cartesian.ctypes.data_as(POINTER(c_double)),
                        c_int(N)
                        )

    return cartesian
