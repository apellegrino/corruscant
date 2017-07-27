from ctypes import CDLL, POINTER, c_double, c_int
import os
import numpy as np

from . import clustering

path_here = os.path.abspath(__file__)
path_dir = os.path.dirname(path_here)
coordlib = CDLL("{:s}/../bin/libcoords.so".format(path_dir))

coordlib.radec2cart64.restype = None
coordlib.radec2cart64.argtypes = [
                                POINTER(c_double), # ra, dec
                                POINTER(c_double), # x,y,z
                                POINTER(c_double),
                                c_int,
                                ]

#def autocorrelation(data_tree, rand_tree, radii, est_type="landy-szalay",
#             err_type="jackknife", num_threads=4):
#    # TODO
#    # call clustering.auto with *args, **kwargs, maybe
#    # indicate in results that the 3D function is being calculated
#
#    pass

def autocorr(data_tree, rand_tree, radii_deg, est_type="landy-szalay",
             err_type="jackknife", num_threads=4):

    radii_euclidean = 2. * np.sin( np.array(radii_deg) * (np.pi / 180.) / 2. )

    results = clustering._autocorr(data_tree, rand_tree, radii_euclidean,
                                   est_type=est_type, err_type=err_type,
                                   num_threads=4)

    results.radii_nominal = np.array(radii_deg)
    results.radii_units = "degrees"

    return results

def cartesian(ra, dec):
    """Returns array of cartesian points projected onto the unit sphere given
    ra, dec arrays.
    """

    #sph = np.require(
    #                np.vstack([ra, dec]).T,
    #                requirements = 'AC'
    #                 )
    #N = sph.shape[0]

    #cartesian = np.empty((N, 3))

    #for in_pt, out_pt in zip(sph, cartesian):
    #    coordlib.radec2cart64(
    #                        in_pt.ctypes.data_as(POINTER(c_double)),
    #                        out_pt.ctypes.data_as(POINTER(c_double)),
    #                        )

    N = ra.size

    ra = np.require(ra, requirements='AC')
    dec = np.require(dec, requirements='AC')

    cartesian = np.empty((N, 3))

    coordlib.radec2cart64(
                        ra.ctypes.data_as(POINTER(c_double)),
                        dec.ctypes.data_as(POINTER(c_double)),
                        cartesian.ctypes.data_as(POINTER(c_double)),
                        c_int(N)
                        )

    return cartesian
