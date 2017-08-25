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

import corruscant
from corruscant import clustering

coordlib = corruscant.coordlib

coordlib.radec2cart64.restype = None
coordlib.radec2cart64.argtypes = [
                                POINTER(c_double), # ra, dec
                                POINTER(c_double), # x,y,z
                                POINTER(c_double),
                                c_int,
                                ]

def autocorr(data_tree, rand_tree, radii_deg, **kwargs):

    radii_deg = np.array(radii_deg)
    radii_euclidean = 2. * np.sin( radii_deg * (np.pi / 180.) / 2. )

    results = clustering._autocorr(data_tree, rand_tree, radii_euclidean,
                                   **kwargs)

    results.radii_nominal = radii_deg
    results.radii_units = "degrees"

    return results

def crosscorr(data_tree_1, data_tree_2, rand_tree_1, rand_tree_2, radii_deg,
              **kwargs):

    radii_deg = np.array(radii_deg)
    radii_euclidean = 2. * np.sin( radii_deg * (np.pi / 180.) / 2. )

    results = clustering._crosscorr(data_tree_1, data_tree_2, rand_tree_1,
                                    rand_tree_2, radii_euclidean, **kwargs)

    results.radii_nominal = radii_deg
    results.radii_units = "degrees"

    return results

def cartesian(ra, dec):
    """Convert points described by arrays of RA and DEC into cartesian points
    projected onto the unit sphere.

    Parameters
    ----------
    ra (array-like): RA values in degrees
    dec (array-like): DEC values in degrees

    Returns
    -------
    Numpy array of shape (N,3) of x,y,z points where N is the length of the
    input arrays
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

    ra = np.require(ra, requirements='AC', dtype='float64')
    dec = np.require(dec, requirements='AC', dtype='float64')

    cartesian = np.empty((N, 3), dtype='float64')

    coordlib.radec2cart64(
                        ra.ctypes.data_as(POINTER(c_double)),
                        dec.ctypes.data_as(POINTER(c_double)),
                        cartesian.ctypes.data_as(POINTER(c_double)),
                        c_int(N)
                        )

    return cartesian
