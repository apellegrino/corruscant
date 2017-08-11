from ctypes import CDLL, Structure, POINTER, c_double, c_int, c_longlong
from os.path import abspath, dirname
import numpy as np
import warnings

class node(Structure):
    pass

class datum(Structure):
    pass

class kdtree(Structure):
    _fields_ = [
        ("node_data",POINTER(node)),
        ("size",c_int),
        ("memsize",c_int),
        ("num_fields",c_int),
        ("data",POINTER(datum)),
        ("fields",POINTER(c_int)),
        ("args",POINTER(c_int)),
                ]

PROJECT_PATH = dirname(dirname(abspath(__file__)))

try:
    kdlib = CDLL("{:s}/bin/libkdtree.so".format(PROJECT_PATH))
except OSError:
    raise OSError("Shared library not found. Did you remember to \"make\"?")

kdlib.tree_construct.restype = kdtree
kdlib.tree_construct.argtypes = [
                                POINTER(c_double), # array of double arrays
                                POINTER(c_int), # field array
                                c_int, # length
                                c_int, # num_fields
                                ]

kdlib.destroy.restype = None
kdlib.destroy.argtypes = [kdtree]

# returning array w/ numbers of pair counts
kdlib.pair_count.restype = np.ctypeslib.ndpointer(dtype=c_longlong,
                                                            shape=(255,) )

kdlib.pair_count.argtypes = [
                            kdtree, # tree to query
                            POINTER(c_double), # data to query with
                            POINTER(c_int), # fields to query with
                            c_int, # data array size
                            c_int, # num_fields
                            c_double, # radius
                            c_int, # num_threads
                            ]

def _query_tree(tree, points, radius, num_threads, fields=None,
                N_fields=1):

    fields_arg = None
    if fields is not None:
        fields_arg = fields.ctypes.data_as(POINTER(c_int))

    counts = kdlib.pair_count(tree.ctree,
                                    points.ctypes.data_as(POINTER(c_double)),
                                    fields_arg,
                                    c_int(points.shape[0]),
                                    c_int(N_fields),
                                    c_double(radius),
                                    c_int(num_threads)
                                        )

    return counts[:N_fields*tree.N_fields].reshape((N_fields,tree.N_fields))

def validate_points(points):
    points = np.array(points)
    N_points = max(points.shape)

    # try to make array of shape (N, 3) if shape is (3, N)
    if not points.shape == (N_points, 3):
        points = points.T

    if not points.shape == (N_points, 3):
        raise ValueError("Array must be of shape (N, 3) "
                         "or (3, N). The provided array has "
                         "shape {:s}".format(str(points.shape)))

    newpoints = np.require(points, requirements='AC')

    if newpoints is not points:
        warnings.warn(
                "Data array is being copied "
                "to have shape {:s}".format(str(newpoints.shape)),
                RuntimeWarning
                      )

    return newpoints

def validate_fields(fields, points):
    fields = np.array(fields)
    N_points = max(points.shape)

    if not fields.shape == (N_points,):
        raise ValueError("Field IDs must be a 1-d array of length equal to "
                         "the number of points")

    newfields = np.require(fields, dtype='int32', requirements='AC')

    if newfields is not fields:
        warnings.warn("Field array is being copied to have type int32",
                       RuntimeWarning)

    return newfields

def est_landy_szalay(dd,dr,rr,dsize,rsize):
    try:
        f = float(rsize) / dsize
    except TypeError:
        f = rsize.astype("float64") / dsize

    return (f*f*np.array(dd) - 2*f*np.array(dr) + np.array(rr)) / np.array(rr)

def est_hamilton(dd,dr,rr,dsize,rsize):
    return np.divide( np.multiply(dd,rr).astype("float64"),
                      np.multiply(dr,dr).astype("float64") ) - 1.

def est_standard(dd,dr,rr,dsize,rsize):
    try:
        f = float(rsize) / dsize
    except TypeError:
        f = rsize.astype("float64") / dsize

    return f * np.divide( dd.astype("float64"), dr ) - 1.

def _autocorr(data_tree, rand_tree, radii, est_type="landy-szalay",
             err_type="jackknife", num_threads=4):
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

    err_type (str) -- what type of error to calculate. (default "jackknife")
        Possible values: "jackknife", "ftf", "poisson", "bootstrap", None

    num_threads -- the maximum number of threads that the C code will create.
    For max performance, set this to the number of logical cores (threads) on
    your machine. (default 4)

    Return:
    A twopoint_data object containing pair count arrays, estimator, error and
    covariance matrix functions and the error type
    """

    if data_tree.N_fields != rand_tree.N_fields:
        raise ValueError("data and random trees must have same number of "
                         "fields")

    if est_type == "landy-szalay":
        estimator = est_landy_szalay
    elif est_type == "hamilton":
        estimator = est_hamilton
    elif est_type == "standard":
        estimator = est_standard
    else:
        raise ValueError("Estimator type for Xi \"{:s}\" "
                         "not valid. Try \"standard\", \"hamilton\", or "
                         "\"landy-szalay\"".format(est_type))

    #valid_err_types = ["jackknife", "ftf", "poisson", "bootstrap", None]

    #if err_type not in valid_err_types:
    #    raise ValueError("Estimator error type \"{:s}\" not "
    #                     "valid".format(err_type))

    #if data_tree.fields is None and err_type is not None:
    #    raise ValueError("Error cannot be calculated when data tree has no "
    #                     "fields")

    #if rand_tree.fields is None and err_type is not None:
    #    raise ValueError("Error cannot be calculated when random tree has no "
    #                     "fields")

    dd = lambda r: _query_tree(data_tree, data_tree.points, r, num_threads,
                                data_tree.fields, data_tree.N_fields)
    dd_array = np.diff([dd(r) for r in radii], axis=0)

    dr = lambda r: _query_tree(rand_tree, data_tree.points, r, num_threads,
                                data_tree.fields, data_tree.N_fields)
    dr_array = np.diff([dr(r) for r in radii], axis=0)

    rr = lambda r: _query_tree(rand_tree, rand_tree.points, r,
                                num_threads, rand_tree.fields, rand_tree.N_fields)
    rr_array = np.diff([rr(r) for r in radii], axis=0)

    data = twopoint_data(dd_array, dr_array, rr_array, data_tree, rand_tree,
                         estimator, radii)

    data.error_type = err_type
    return data

class tree:
    def __init__(self, points, fields=None):

        points = validate_points(points)

        if fields is None:
            self.fields = None
            self.N_fields = 1
        else: 
            fields = validate_fields(fields, points)

            unique = np.unique(fields)
            if not np.all(unique == range(np.max(fields)+1)):
                raise ValueError("fields must be a 1-D array of integers from "
                                 "1 to the number of unique fields minus one")

            self.N_fields = unique.size

            if self.N_fields < 2:
                raise ValueError("At least two unique field IDs must be "
                                 "provided if using error fields")

            #min_field, max_field = np.min(fields), np.max(fields)

            #if min_field != 1:
            #    raise ValueError("Minimum field ID must be 1")
            #if max_field != self.N_fields:
            #    raise ValueError("Maximum field ID must be the same as the "
            #                     "number of unique field "
            #                     "IDs ({:d})".format(self.N_fields))
            self.fields = fields
            self.field_sizes = self._calc_field_sizes()

        self.points = points

        self._make_tree()

        self.size = self.ctree.size

    def _calc_field_sizes(self):
        f = self.fields
        return np.array([f[f == id].size for id in range(0,self.N_fields)])

    def _make_tree(self):
        data = self.points.ctypes.data_as(POINTER(c_double))

        f = None
        if self.fields is not None:
            f = self.fields.ctypes.data_as(POINTER(c_int))
        self.ctree = kdlib.tree_construct(data, f, 
                                          c_int(self.points.shape[0]),
                                          c_int(self.N_fields) )

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
        self.radii_euclidean = radii
        self.radii_nominal = None # to be set by wrappers of `_autocorr()`
        self.radii_units = None # to be set by wrappers of `_autocorr()`
        self.error_type = None

    def total_pair_counts(self):
        # sum each matrix across field indices 1 and 2, leaving bin index 0
        return (
                np.sum(self.dd, axis=(1,2)),
                np.sum(self.dr, axis=(1,2)),
                np.sum(self.rr, axis=(1,2)),
                )

    def ftf_pair_counts(self, fid):
        # get the diagonal element with index `fid`
        return self.dd[:,fid,fid], self.dr[:,fid,fid], self.rr[:,fid,fid]

    def jackknife_pair_counts(self, fid):
        # ignore all pair counts between the field `fld` and any other, and
        # sum the whole matrix
        dd_copy = np.copy(self.dd)
        dr_copy = np.copy(self.dr)
        rr_copy = np.copy(self.rr)

        dd_copy[:,fid,:] = 0
        dd_copy[:,:,fid] = 0

        dr_copy[:,fid,:] = 0
        dr_copy[:,:,fid] = 0

        rr_copy[:,fid,:] = 0
        rr_copy[:,:,fid] = 0

        return (
                np.sum(dd_copy, axis=(1,2)),
                np.sum(dr_copy, axis=(1,2)),
                np.sum(rr_copy, axis=(1,2)),
                )

    def bootstrap_pair_counts(self, N_trials):
        """Perform bootstrap resampling on the results, and return the
        estimator value.

        Parameters
        ----------
        N_trials (int): The number of resamples to perform.

        Returns
        -------
        A 2D numpy array of shape (R, N_trials) where R is the number of radii
        bins, i.e. the number of radii values minus one.
        """
        Nfld = self.dtree.N_fields

        # random field data assumed to be the same
        resample = np.random.randint(
                                0, Nfld, N_trials*Nfld
                                     ).reshape(N_trials, Nfld)

        hists = np.apply_along_axis(np.bincount, 1, resample, minlength=Nfld)

        # Create a matrix M where M(i,j) is product of frequencies of ith, jth
        # fields
        func = lambda x: np.multiply(*np.meshgrid(x, x))
        freqs = np.apply_along_axis(func, 1, hists)

        # compute new field sizes
        dsizes = np.sum(hists * self.dtree.field_sizes, axis=1)
        rsizes = np.sum(hists * self.rtree.field_sizes, axis=1)

        # Dot product over dimensions representing different fields only,
        # leaving num. of trials and radii intact
        dd = np.tensordot(self.dd, freqs, axes=((1,2),(1,2)))
        dr = np.tensordot(self.dr, freqs, axes=((1,2),(1,2)))
        rr = np.tensordot(self.rr, freqs, axes=((1,2),(1,2)))

        return dd, dr, rr, dsizes, rsizes

    def bootstrap_error(self, N_trials=1000):
        dd, dr, rr, dsizes, rsizes = self.bootstrap_pair_counts(N_trials)
        est = self.estimator(dd, dr, rr, dsizes, rsizes)
        return np.std(est, axis=1)

    def estimate(self):
        estfunc = self.estimator
        dd, dr, rr = self.total_pair_counts()
        return estfunc(dd, dr, rr, self.dtree.size, self.rtree.size)

    def covariance(self):

        dd_tot, dr_tot, rr_tot = self.total_pair_counts()

        if self.error_type is None:
            return None
    
        elif self.error_type == "poisson":
            return None

        est_func = self.estimator

        cov = np.zeros((self.radii_euclidean.size-1,
                        self.radii_euclidean.size-1)
                       )

        estimation = self.estimate()

        if self.error_type == "jackknife":
            # sizes with one field out
            dfield_sizes = self.dtree.size - self.dtree.field_sizes
            rfield_sizes = self.rtree.size - self.rtree.field_sizes

            for fid in range(self.dtree.N_fields):
                dd, dr, rr = self.jackknife_pair_counts(fid)

                est_per_field = est_func(dd, dr, rr, 
                                         dfield_sizes[fid],rfield_sizes[fid])

                diff = est_per_field - estimation

                # covariance matrix
                rr_quotient = np.sqrt(rr.astype('float64') / rr_tot)
                rr_subi, rr_subj = np.meshgrid(rr_quotient,rr_quotient)
                
                xi_subi, xi_subj = np.meshgrid(diff,diff)
                cov += rr_subi*xi_subi*rr_subj*xi_subj

        elif self.error_type == "ftf":
            # sizes with one field in
            dfield_sizes = self.dtree.field_sizes
            rfield_sizes = self.rtree.field_sizes

            for fid in range(self.dtree.N_fields):
                dd, dr, rr = self.ftf_pair_counts(fid)

                est_per_field = est_func(dd, dr, rr, 
                                         dfield_sizes[fid],rfield_sizes[fid])

                diff = est_per_field - estimation

                # covariance matrix
                rr_quotient = np.sqrt(rr.astype('float64') / rr_tot)
                rr_subi, rr_subj = np.meshgrid(rr_quotient,rr_quotient)
                
                xi_subi, xi_subj = np.meshgrid(diff,diff)
                cov += rr_subi*xi_subi*rr_subj*xi_subj

            cov /= float(self.dtree.N_fields - 1)

        return cov

    def error(self):
        if self.error_type is None:
            return None
        elif self.error_type == "poisson":
            dd_tot, dr_tot, rr_tot = self.total_pair_counts()

            estimation = self.estimate()
            with np.errstate(divide='ignore', invalid='ignore'):
                error = np.divide(1 + estimation,np.sqrt(dd_tot))
                error[np.isneginf(error)] = np.nan
                error[np.isinf(error)] = np.nan

            return error
        elif self.error_type == "bootstrap":
            return self.bootstrap_error()
        else:
            cov = self.covariance()

            return np.sqrt(np.diagonal(cov))

    def normalized_covariance(self):
        # normalize covariance matrix to get the regression matrix
        cov = self.covariance()

        sigmas = np.sqrt(np.diagonal(cov))
        i_divisor, j_divisor = np.meshgrid(sigmas,sigmas)
        total_divisor = np.multiply(i_divisor,j_divisor)
        return np.divide(cov,total_divisor)

    def __str__(self):
        r_lower = self.radii_nominal[:-1]
        r_upper = self.radii_nominal[1:]
        dd_tot, dr_tot, rr_tot = self.total_pair_counts()
        est = self.estimate()
        err = self.error()

        if err is None:
            err = [None] * len(est)

        # calc column widths for pair counts
        dd_width = max([len(str(count)) for count in dd_tot]) + 2
        dr_width = max([len(str(count)) for count in dr_tot]) + 2
        rr_width = max([len(str(count)) for count in rr_tot]) + 2

        lines = [ "\n" ]
        labels = ["Bin L", "Bin R", "DD",     "DR",     "RR",     "Estimator", "Error"]
        sizes =  [9,       9,       dd_width, dr_width, rr_width, 11,          10]

        header = [label.ljust(size) for label, size in zip(labels, sizes)]

        lines.append(''.join(header))

        # add units to radii column, e.g. degrees
        lines.append("{:^14}".format("({:s})".format(self.radii_units)))


        lines.append("-" * 79)

        ########## end of table header ##########

        for rl, ru, dd, dr, rr, estv, errv in \
                zip(r_lower, r_upper, dd_tot, dr_tot, rr_tot, est, err):

            rl_s =  "{:.2E}".format(rl)
            ru_s =  "{:.2E}".format(ru)
            estv_s =  "{:+.2E}".format(estv)

            if errv is None:
                errv_s = str(None)
            else:
                errv_s =  "{:.2E}".format(errv)

            s = [
                    "{:<9}".format(rl_s),
                    "{:<9}".format(ru_s),
                    "{:<{:d}}".format(dd, dd_width),
                    "{:<{:d}}".format(dr, dr_width),
                    "{:<{:d}}".format(rr, rr_width),
                    "{:<11}".format(estv_s),
                    "{:<10}".format(errv_s),
                 ]

            lines.append(''.join(s))

        return '\n'.join(lines)

