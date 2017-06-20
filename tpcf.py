from ctypes import CDLL, Structure, POINTER, c_double, c_int, c_longlong
from os.path import abspath, dirname
import numpy as np

class node(Structure):
    pass

class array3d(Structure):
    _fields_ = [
        ("x",POINTER(c_double)),
        ("y",POINTER(c_double)),
        ("z",POINTER(c_double)),
        ("fields",POINTER(c_int)),
        ("num_fields",c_int),
        ("size",c_int),
        ]

class argarray3d(Structure):
    _fields_ = [
        ("x",POINTER(c_int)),
        ("y",POINTER(c_int)),
        ("z",POINTER(c_int)),
        ("size",c_int),
        ]

class kdtree(Structure):
    _fields_ = [
        ("node_data",POINTER(node)),
        ("size",c_int),
        ("memsize",c_int),
        ("data",array3d),
        ("arg_data",argarray3d),
        ]

path_here = abspath(__file__)
path_dir = dirname(path_here)
kdlib = CDLL("%s/bin/libkdtree.so" % path_dir)

kdlib.tree_construct.restype = kdtree
kdlib.tree_construct.argtypes = [array3d]

kdlib.destroy.restype = None
kdlib.destroy.argtypes = [kdtree]

# returning array w/ numbers of pair counts
kdlib.pair_count_jackknife.restype = np.ctypeslib.ndpointer(dtype=c_longlong, shape=(255,))

kdlib.pair_count_jackknife.argtypes = [
                            kdtree, # tree to query
                            array3d, # data to query with
                            c_double, # radius
                            c_int, # num_threads
                            ]

kdlib.pair_count_ftf.restype = np.ctypeslib.ndpointer(dtype=c_longlong, shape=(255,))

kdlib.pair_count_ftf.argtypes = [
                            kdtree, # tree to query
                            array3d, # data to query with
                            c_double, # radius
                            c_int, # num_threads
                            ]

kdlib.pair_count_noerr.restype = np.ctypeslib.ndpointer(dtype=c_longlong, shape=(255,))

kdlib.pair_count_noerr.argtypes = [
                            kdtree, # tree to query
                            array3d, # data to query with
                            c_double, # radius
                            c_int, # num_threads
                            ]
def _unpack(data):
    return tuple([np.copy(row) for row in data])

# whatever arrays that ctypes.data_as() is called on must remain outside the
# scope of this function because .ctypes uses a reference to it, otherwise the
# references will be dangling
def make_clike_array(x, y, z, fields, N_fields):
    if fields is None:
        fields = np.zeros_like(x, dtype=np.int32)

    if not fields.dtype == np.int32:
        raise ValueError("Fields must be of dtype int32")

    array = array3d()
    array.x = x.ctypes.data_as(POINTER(c_double))
    array.y = y.ctypes.data_as(POINTER(c_double))
    array.z = z.ctypes.data_as(POINTER(c_double))
    array.fields = fields.ctypes.data_as(POINTER(c_int))
    array.num_fields = c_int(N_fields)
    array.size = c_int(x.shape[0])
    return array

def _make_tree(points, fields, N_fields=1):
    points_x, points_y, points_z = _unpack(points)

    # numpy will use 64-bit ints unless otherwise specified
    # type change must happen on ITS OWN LINE!
    fields = fields.astype('int32')

    array = make_clike_array(points_x,points_y,points_z,fields,N_fields)
    tree = kdlib.tree_construct(array)

    return tree

def _query_tree(tree, points, radius, num_threads, errtype, fields=None, N_fields=0):
    points_x, points_y, points_z = _unpack(points)

    N_fields = tree.N_fields
    
    # numpy will use 64-bit ints unless otherwise specified
    # type change must happen on ITS OWN LINE!
    fields = fields.astype('int32')

    array = make_clike_array(points_x,points_y,points_z,fields,N_fields)
    counts = None
    if errtype == 'jackknife':
        counts = kdlib.pair_count_jackknife(tree.ctree,array,radius,num_threads)
    elif errtype == 'field-to-field':
        counts = kdlib.pair_count_ftf(tree.ctree,array,radius,num_threads)
    else:
        counts = kdlib.pair_count_noerr(tree.ctree,array,radius,num_threads)

    return counts[:N_fields+1]

def validate_points(arr):
    if not type(arr) is np.ndarray:
        arr = np.array(arr)

    if not (arr.ndim == 2):
        raise ValueError("Point array must be two-dimensional (i.e. a \
                        list of points)")

    # try to make array of shape (3, N) if shape is (N, 3)
    if (arr.shape[1] == 3):
        arr = arr.T

    if not(arr.shape[0] == 3):
        raise ValueError("Array must be of shape (3, N) or (N, 3). The \
                        provided array has shape %s" % str(arr.shape))

    return arr

def validate_fields(fields, points):
    fields = np.array(fields)
    N_points = max(points.shape)

    if fields.shape != (N_points,):
        raise ValueError("Field IDs must be a 1-d array of length equal to \
                        the number of points")

    return fields

def est_landy_szalay(dd,dr,rr,dsize,rsize):
    f = float(rsize)/dsize
    return ( f*f*np.array(dd) - 2*f*np.array(dr) + np.array(rr) ) / np.array(rr)

def est_hamilton(dd,dr,rr,dsize,rsize):
    return np.divide( np.multiply(dd,rr).astype("float64"),
                      np.multiply(dr,dr).astype("float64") ) - 1.

def est_standard(dd,dr,rr,dsize,rsize):
    return float(rsize)/dsize * np.divide( dd.astype("float64"), dr ) - 1.

def twopoint(data_tree, rand_tree, radii, est_type="landy-szalay",
             err_type='jackknife', num_threads=4):
    """Given a set of 3D cartesian data points and random points, calculate the
    estimated two-point correlation function with error estimation.

    Arguments:
    data_tree (tree) -- tree of data points
    rand_tree (tree) -- tree of random points
    radii (array-like) -- floats which define the radius bin sizes

    Keyword arguments:
    est_type -- the type of estimator to use for the correlation function.
    (default "landy-szalay")
        Possible values: "standard", "landy-szalay", "hamilton"

    err_type (str) -- what type of error to calculate. (default None)
        Possible values: "jackknife", "field-to-field", "poisson", None

    num_threads -- the maximum number of threads that the C code will create.
    For max performance, set this to the number of logical cores (threads) on
    your machine. (default 4)

    Return:
    A twopoint_data object containing pair count arrays, estimator, error and
    covariance matrix functions and the error type
    """

    if data_tree.N_fields != rand_tree.N_fields:
        raise ValueError("data and random trees must have same number of fields")

    if est_type == "landy-szalay":
        estimator = est_landy_szalay
    elif est_type == "hamilton":
        estimator = est_hamilton
    elif est_type == "standard":
        estimator = est_standard
    else:
        raise ValueError("Estimator type for Xi %s not valid" % est_type)

    valid_err_types = ["jackknife", "field-to-field", "poisson", None]
    if err_type not in valid_err_types:
        raise ValueError("Estimator error type %s not valid" % err_type)
    if data_tree.fields is None and err_type is not None:
        raise ValueError("Error cannot be calculated when data tree has no fields")
    if rand_tree.fields is None and err_type is not None:
        raise ValueError("Error cannot be calculated when random tree has no fields")

    dd = lambda r: _query_tree(data_tree, data_tree.points, r, num_threads,
                                err_type, data_tree.fields)
    dd_array = np.diff([dd(r) for r in radii], axis=0)

    dr = lambda r: _query_tree(rand_tree, data_tree.points, r, num_threads,
                                err_type, data_tree.fields)
    dr_array = np.diff([dr(r) for r in radii], axis=0)

    rr_array = None
    if not est_type == "standard":
        rr = lambda r: _query_tree(rand_tree, rand_tree.points, r,
                                    num_threads, err_type, rand_tree.fields)
        rr_array = np.diff([rr(r) for r in radii], axis=0)
    

    data = twopoint_data(dd_array, dr_array, rr_array, data_tree, rand_tree, estimator, radii)
    data.error_type = err_type

    return data

class tree:
    def __init__(self, points, fields=None):
        points = validate_points(points)

        if fields is None:
            self.fields = None
            self.N_fields = 0
        else: 
            fields = validate_fields(fields, points)

            self.N_fields = np.unique(fields).size

            if self.N_fields < 2:
                raise ValueError("At least two unique field IDs must be provided")

            min_field, max_field = np.min(fields), np.max(fields)

            if min_field != 1:
                raise ValueError("Minimum field ID must be 1")
            if max_field != self.N_fields:
                raise ValueError("Maximum field ID must be the same as the number \
                                        of unique field IDs (%d)" % self.N_fields)
            self.fields = fields
            self.field_sizes = self._calc_field_sizes()

        self.points = points
        self.size = points.shape[1]

        self.ctree = _make_tree(points, fields, self.N_fields)

    def _calc_field_sizes(self):
        fields = self.fields
        return np.array([fields[fields == id].size for id in range(1,self.N_fields+1)])

    def __del__(self):
        kdlib.destroy(self.ctree)

class twopoint_data:
    def __init__(self, dd, dr, rr, dtree, rtree, estimator, radii):
        self.dd = dd
        self.dr = dr
        self.rr = rr
        self.dtree = dtree
        self.rtree = rtree
        self.estimator = estimator
        self.radii = radii
        self.error_type = None

    def total_pair_counts(self):
        if self.rr is None:
            return self.dd[:,0], self.dr[:,0], None
        else:
            return self.dd[:,0], self.dr[:,0], self.rr[:,0]

    def field_pair_counts(self, fid):
        if self.rr is None:
            return self.dd[:,fid], self.dr[:,fid], None
        else:
            return self.dd[:,fid], self.dr[:,fid], self.rr[:,fid]

    def estimate(self):
        estfunc = self.estimator
        dd, dr, rr = self.total_pair_counts()
        self.estimation = estfunc(dd, dr, rr, self.dtree.size, self.rtree.size)
        return self.estimation

    def covariance(self):

        dd_tot, dr_tot, rr_tot = self.total_pair_counts()

        if self.error_type is None:
            return None, None
    
        elif self.error_type == "poisson":
            with np.errstate(divide='ignore', invalid='ignore'):
                error = np.divide(1 + estimation,np.sqrt(dd_tot))
                error[np.isneginf(error)] = np.nan
                error[np.isinf(error)] = np.nan
            return error, None

        elif self.error_type == "jackknife":
            # sizes with one field out
            dfield_sizes = self.dtree.size - self.dtree.field_sizes
            rfield_sizes = self.rtree.size - self.rtree.field_sizes

        elif self.error_type == "field-to-field":
            # sizes with one field in
            dfield_sizes = self.dtree.field_sizes
            rfield_sizes = self.rtree.field_sizes

        est_func = self.estimator

        cov = np.zeros((self.radii.size-1,self.radii.size-1))

        for fid in range(1,self.dtree.N_fields+1):
            dd, dr, rr = self.field_pair_counts(fid)
            est_per_field = est_func(dd,dr,rr,dfield_sizes[fid-1],rfield_sizes[fid-1])

            diff = est_per_field - self.estimation

            # covariance matrix
            rr_quotient = np.sqrt(rr.astype('float64')/rr_tot)
            rr_subi, rr_subj = np.meshgrid(rr_quotient,rr_quotient)
            
            xi_subi, xi_subj = np.meshgrid(diff,diff)
            cov += rr_subi*xi_subi*rr_subj*xi_subj

        if self.error_type == "field-to-field":
            cov /= float(self.dtree.N_fields - 1)

        return cov

    def error(self):
        cov = self.covariance()

        return np.sqrt(np.diagonal(cov))

    def normalized_covariance(self):
        # normalize covariance matrix to get the regression matrix
        cov = self.covariance()

        sigmas = np.sqrt(np.diagonal(cov))
        i_divisor, j_divisor = np.meshgrid(sigmas,sigmas)
        total_divisor = np.multiply(i_divisor,j_divisor)
        return np.divide(cov,total_divisor)
