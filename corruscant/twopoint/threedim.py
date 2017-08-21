"""
Copyright (C) 2016-2017 Andrew Pellegrino

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
from __future__ import division, print_function

from ctypes import CDLL, POINTER, c_double, c_int
import ctypes.util

import numpy as np
from os.path import dirname, abspath

import corruscant
from corruscant import clustering

coordlib = corruscant.coordlib

coordlib.radecdist2cart64.restype = None
coordlib.radecdist2cart64.argtypes = [
                                POINTER(c_double), # ra
                                POINTER(c_double), # dec
                                POINTER(c_double), # dist
                                POINTER(c_double), # x,y,z
                                c_int, # N
                                ]

def autocorr(data_tree, rand_tree, radii, **kwargs):

    results = clustering._autocorr(data_tree, rand_tree, radii,
                                   **kwargs)
    results.radii_nominal = results.radii_euclidean
    results.radii_units = "Data units"

    return results

def crosscorr(data_tree_1, data_tree_2, rand_tree_1, rand_tree_2, radii,
              **kwargs):
    
    results = clustering._crosscorr(data_tree_1, data_tree_2, rand_tree_1,
                                    rand_tree_2, radii, **kwargs)
    results.radii_nominal = results.radii_euclidean
    results.radii_units = "Data units"

    return results

def cartesian(ra, dec, distance=None, z=None, cosmology=None):
    """Creates array of cartesian points given ra, dec, distance arrays or ra,
    dec, redshift arrays with an Astropy cosmology. Provided distances take
    precedence over provided redshift and cosmology. Either `distance` or `z`
    and `cosmology` must be provided.

    Parameters
    ----------
    ra (array-like): RA input values
    dec (array-like): DEC input values
    distance (array-like, optional): distance input values
    z (array-like, optional): redshift input values
    cosmology (astropy.cosmology object, optional): a cosmology to calculate
    comoving distances from redshifts

    Returns
    -------
    Numpy array of shape (N,3) of x,y,z points where N is the length of the
    input arrays
    """

    if distance is None:
        # calculate comoving distance
        z = np.array(z)
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

    # require float64 before passing to double C function
    ra = np.require(ra, requirements='AC', dtype="float64")
    dec = np.require(dec, requirements='AC', dtype="float64")
    distance = np.require(distance, requirements='AC', dtype="float64")

    N = ra.size
    cartesian = np.empty((N, 3))

    coordlib.radecdist2cart64(
                        ra.ctypes.data_as(POINTER(c_double)),
                        dec.ctypes.data_as(POINTER(c_double)),
                        distance.ctypes.data_as(POINTER(c_double)),
                        cartesian.ctypes.data_as(POINTER(c_double)),
                        c_int(N)
                        )

    return cartesian
