from ctypes import CDLL, Structure, POINTER, c_double, c_int, c_longlong
from os.path import abspath, dirname
import numpy as np

class node(Structure):
    pass

class vpoint(Structure):
    pass

class vptree(Structure):
    _fields_ = [
        ("node_data",POINTER(node)),
        ("size",c_int),
        ("memsize",c_int),
        ]

path_here = abspath(__file__)
path_dir = dirname(path_here)
vplib = CDLL("%s/bin/libvptree.so" % path_dir)

vplib.make_vp_tree.restype = vptree
vplib.make_vp_tree.argtypes = [POINTER(c_double), POINTER(c_double), POINTER(c_int), c_int]

vplib.destroy.restype = None
vplib.destroy.argtypes = [vptree]

# returning array w/ numbers of pair counts
vplib.pair_count_noerr.restype = np.ctypeslib.ndpointer(dtype=c_longlong, shape=(255,))

vplib.pair_count_noerr.argtypes = [
                            vptree, # tree to query
                            POINTER(c_double), # data to query with
                            POINTER(c_double), # data to query with
                            c_int, # array size
                            c_double, # radius
                            c_int, # num. threads
                            ]

vplib.pair_count_jackknife.restype = np.ctypeslib.ndpointer(dtype=c_longlong, shape=(255,))

vplib.pair_count_jackknife.argtypes = [
                            vptree, # tree to query
                            POINTER(c_double), # data to query with
                            POINTER(c_double), # data to query with
                            POINTER(c_int), # array of fields
                            c_int, # array size
                            c_int, # num. of fields
                            c_double, # radius
                            c_int, # num. threads
                            ]

vplib.pair_count_ftf.restype = np.ctypeslib.ndpointer(dtype=c_longlong, shape=(255,))

vplib.pair_count_ftf.argtypes = [
                            vptree, # tree to query
                            POINTER(c_double), # data to query with
                            POINTER(c_double), # data to query with
                            POINTER(c_int), # array of fields
                            c_int, # array size
                            c_int, # num. of fields
                            c_double, # radius
                            c_int, # num. threads
                            ]

def validate_points(points, ndims):
    """Check that the point array has the proper dimensions, and make the view
    correspond to parallel arrays. `ndims` represents the dimensionality of the
    space which the points lie in, not the dimensionality of the array
    """

    points = np.asarray(points)
 
    if not type(points) is np.ndarray:
        points = np.array(points)

    if not points.ndim == 2:
        raise ValueError("Point array must be two-dimensional (i.e. a list of points)")

    # first dimension should be ndims, second can be anything. If first is not
    # ndims, transpose to handle (N, ndims) case before raising an exception
    oldshape = points.shape
    if points.shape[1] == ndims:
        points = points.T

    if not points.shape[0] == ndims:
        raise ValueError("Array must be of shape (%d, N) or (N, %d). The provided array has shape %s" % (ndims, ndims, str(oldshape)))

    return points

def validate_fields(fields, points):
    fields = np.asarray(fields)
    N_points = max(points.shape)

    if fields.shape != (N_points,):
        raise ValueError("Field IDs must be a 1-d array of length equal to the number of points (shape %s, actual shape %s)" % (str((N_points,)), str(fields.shape)))

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
             err_type="jackknife", num_threads=4):
    """Given a tree of data points and a tree of random points, estimate the
    two-point correlation function with error estimation.

    Arguments:
    data_tree (tree) -- tree of data points
    rand_tree (tree) -- tree of random points
    radii (array-like) -- floats which define the radius bin edges

    Keyword arguments:
    est_type -- the type of estimator to use for the correlation function.
    (default "landy-szalay")
        Possible values: "standard", "landy-szalay", "hamilton"

    err_type (str) -- what type of error to calculate. (default jackknife)
        Possible values: "jackknife", "field-to-field", "poisson", None

    num_threads -- the maximum number of threads that the C code will create.
    For max performance, set this to the number of logical cores (threads) on
    your machine. (default 4)

    Return:
    A twopoint_data object containing pair count arrays, estimator, error and
    covariance matrix functions and the error type
    """

    if data_tree.N_fields != rand_tree.N_fields:
        raise ValueError("data, random trees must have same number of fields")

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
    if data_tree.fields is None and err_type == "jackknife" or err_type == "field-to-field":
        raise ValueError("%s error cannot be calculated when data tree has no fields" % err_type)
    if rand_tree.fields is None and err_type == "jackknife" or err_type == "field-to-field":
        raise ValueError("%s error cannot be calculated when random tree has no fields" % err_type)

    dd = lambda r: data_tree._query(data_tree, r, num_threads, err_type)
    dd_array = np.diff([dd(r) for r in radii], axis=0)

    dr = lambda r: rand_tree._query(data_tree, r, num_threads, err_type)
    dr_array = np.diff([dr(r) for r in radii], axis=0)

    rr_array = None
    if not est_type == "standard":
        rr = lambda r: rand_tree._query(rand_tree, r, num_threads, err_type)
        rr_array = np.diff([rr(r) for r in radii], axis=0)

    data = twopoint_data(dd_array, dr_array, rr_array, data_tree, rand_tree,
                            estimator, radii)
    data.error_type = err_type

    return data

class tree:
    def __init__(self, points, fields=None):
        points = validate_points(points, 2)

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
                raise ValueError("Maximum field ID must be the same as the \
                            number of unique field IDs (%d)" % self.N_fields)

            self.fields = np.require(fields, dtype='int32', requirements='CAWO')
            self.field_sizes = self._calc_field_sizes()

        # precondition: parallel arrays
        self.points = np.require(points, dtype='float64', requirements='CAWO')
        self.size = points.shape[1]

        self._make_ctree()

    def _calc_field_sizes(self):
        """Calculate the number of points in each field"""
        fields = self.fields
        return np.array([fields[fields == id].size for id in range(1,self.N_fields+1)])

    def _make_ctree(self):
        """Call the C function to build the tree"""
        lat_ptr, lon_ptr = [dim.ctypes.data_as(POINTER(c_double)) for dim in self.points]
        fields_ptr = self.fields.ctypes.data_as(POINTER(c_int))
        self.ctree = vplib.make_vp_tree(lat_ptr, lon_ptr, fields_ptr, c_int(self.size))

    def _query(self, qtree, radius, num_threads, err_type):
        """Get the number of points on this tree inside a distance radius for
        all points on qtree"""
        lat_ptr, lon_ptr = [dim.ctypes.data_as(POINTER(c_double)) for dim in qtree.points]

        if qtree.fields is not None:
            fields = qtree.fields.ctypes.data_as(POINTER(c_int))

        if err_type == "jackknife":
            counts = vplib.pair_count_jackknife(
                                        self.ctree,lat_ptr,lon_ptr,fields,
                                        c_int(qtree.size),c_int(self.N_fields),
                                        c_double(radius),c_int(num_threads)
                                        )
        elif err_type == "field-to-field":
            counts = vplib.pair_count_ftf(
                                        self.ctree,lat_ptr,lon_ptr,fields,
                                        c_int(qtree.size),c_int(self.N_fields),
                                        c_double(radius),c_int(num_threads)
                                        )
        else:
            counts = vplib.pair_count_noerr(
                                        self.ctree,lat_ptr,lon_ptr,
                                        c_int(qtree.size),
                                        c_double(radius),c_int(num_threads)
                                        )

        return counts[:self.N_fields+1]

    #def __del__(self):
    #    vplib.destroy(self.ctree)

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

        if self.error_type is None:
            return None
        elif self.error_type == "poisson":
            return None

        dd_tot, dr_tot, rr_tot = self.total_pair_counts()

        if self.error_type == "jackknife":
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
        if self.error_type is None:
            return None
        elif self.error_type == "poisson":
            with np.errstate(divide='ignore', invalid='ignore'):
                error = np.divide(1 + self.estimation,np.sqrt(dd_tot))
                error[np.isneginf(error)] = np.nan
                error[np.isinf(error)] = np.nan
            return error
        else:
            cov = self.covariance()
            return np.sqrt(np.diagonal(cov))

    def normalized_covariance(self):
        # normalize covariance matrix to get the regression matrix
        cov = self.covariance()

        if cov is None:
            return None

        sigmas = np.sqrt(np.diagonal(cov))
        i_divisor, j_divisor = np.meshgrid(sigmas,sigmas)
        total_divisor = np.multiply(i_divisor,j_divisor)
        return np.divide(cov,total_divisor)
