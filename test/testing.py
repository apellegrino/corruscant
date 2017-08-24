import numpy as np
from corruscant import twopoint

def random_sphere(N_data, N_fields):
    ra = np.random.rand(N_data) * 360.
    dec = np.arcsin(np.random.rand(N_data) * 2. - 1.) * 180. / np.pi
    xyz = twopoint.angular.cartesian(ra, dec)
    fields = (ra / 360. * N_fields).astype('int32')
    return xyz, fields

def random_cube(N_data, N_fields):
    xyz = np.random.rand(N_data, 3)
    fields = (xyz[:,0] * N_fields).astype('int32')
    return xyz, fields
