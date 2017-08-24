import numpy as np
from corruscant import twopoint

def random_sphere(N_data, N_fields=None):
    ra = np.random.rand(N_data) * 360.
    dec = np.arcsin(np.random.rand(N_data) * 2. - 1.) * 180. / np.pi
    xyz = twopoint.angular.cartesian(ra, dec)

    if N_fields is not None:
        fields = (ra / 360. * N_fields).astype('int32')
        return xyz, fields
    else:
        return xyz

def random_cube(N_data, N_fields=None):
    xyz = np.random.rand(N_data, 3)

    if N_fields is not None:
        fields = (xyz[:,0] * N_fields).astype('int32')
        return xyz, fields
    else:
        return xyz
