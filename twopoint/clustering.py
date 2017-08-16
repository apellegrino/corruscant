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
kdlib.pair_count.restype = POINTER(c_longlong)

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

    return np.ctypeslib.as_array(counts, shape=(N_fields, tree.N_fields))

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

def est_landy_szalay(*args):
    if len(args) == 3:
        d1d2, d1r2, r1r2 = args
        d2r1 = d1r2
    elif len(args) == 4:
        d1d2, d1r2, d2r1, r1r2 = args

    return (d1d2 - d1r2 - d2r1 + r1r2) / r1r2

def est_hamilton(*args):
    if len(args) == 3:
        d1d2, d1r2, r1r2 = args
        d2r1 = d1r2

    elif len(args) == 4:
        d1d2, d1r2, d2r1, r1r2 = args

    return ( d1d2 * r1r2 / (d1r2 * d2r1) ) - 1. 

def est_standard(*args):
    if len(args) == 3:
        d1d2, d1r2, r1r2 = args
        d2r1 = d1r2

    elif len(args) == 4:
        d1d2, d1r2, d2r1, r1r2 = args

    return d1d2 / d1r2 - 1.

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

    if not (data_tree.N_fields == rand_tree.N_fields):
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

    # TODO simplify call with .fields, .N_fields
    dd = lambda r: _query_tree(data_tree, data_tree.points, r, num_threads,
                                data_tree.fields, data_tree.N_fields)
    dd_array = np.diff([dd(r) for r in radii], axis=0)

    dr = lambda r: _query_tree(rand_tree, data_tree.points, r, num_threads,
                                data_tree.fields, data_tree.N_fields)
    dr_array = np.diff([dr(r) for r in radii], axis=0)

    rr = lambda r: _query_tree(rand_tree, rand_tree.points, r,
                               num_threads, rand_tree.fields,
                               rand_tree.N_fields)

    rr_array = np.diff([rr(r) for r in radii], axis=0)

    results = twopoint_autocorr_results(dd_array, dr_array, rr_array, data_tree,
                                     rand_tree, estimator, radii)

    results.error_type = err_type
    return results

def _crosscorr(data_tree_1, data_tree_2, rand_tree_1, rand_tree_2, radii,
               est_type="landy-szalay", err_type="jackknife", num_threads=4):
    
    if not (data_tree_1.N_fields == data_tree_2.N_fields == \
            rand_tree_1.N_fields == rand_tree_2.N_fields):
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

    if data_tree_1.size > data_tree_2.size:
        d1d2 = lambda r: _query_tree(data_tree_1, data_tree_2.points, r,
                                     num_threads, data_tree_2.fields,
                                     data_tree_2.N_fields)
    else:
        d1d2 = lambda r: _query_tree(data_tree_2, data_tree_1.points, r,
                                     num_threads, data_tree_1.fields,
                                     data_tree_1.N_fields)

    d1d2_array = np.diff([d1d2(r) for r in radii], axis=0)

    d1r2 = lambda r: _query_tree(rand_tree_2, data_tree_1.points, r, num_threads,
                                data_tree_1.fields, data_tree_1.N_fields)
    d1r2_array = np.diff([d1r2(r) for r in radii], axis=0)

    d2r1 = lambda r: _query_tree(rand_tree_1, data_tree_2.points, r, num_threads,
                                data_tree_2.fields, data_tree_2.N_fields)
    d2r1_array = np.diff([d2r1(r) for r in radii], axis=0)

    if rand_tree_1.size > rand_tree_2.size:
        r1r2 = lambda r: _query_tree(rand_tree_1, rand_tree_2.points, r,
                                     num_threads, rand_tree_2.fields,
                                     rand_tree_2.N_fields)
    else:
        r1r2 = lambda r: _query_tree(rand_tree_2, rand_tree_1.points, r,
                                     num_threads, rand_tree_1.fields,
                                     rand_tree_1.N_fields)

    r1r2_array = np.diff([r1r2(r) for r in radii], axis=0)

    results = twopoint_crosscorr_results(d1d2_array, d1r2_array, d2r1_array,
                                     r1r2_array, data_tree_1, data_tree_2,
                                     rand_tree_1, rand_tree_2, estimator,
                                     radii)

    results.error_type = err_type
    return results

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

class twopoint_results(object):
    def __init__(self, estimator, radii):
        self.estimator = estimator
        self.nbins = len(radii)-1
        self.radii_euclidean = radii
        self.radii_nominal = None # to be set by wrappers of `_autocorr()`
        self.radii_units = None # to be set by wrappers of `_autocorr()`
        self.error_type = None

    def total_pair_counts(self, normalized):
        # sum each matrix across field indices 1 and 2, leaving bin index 0

        cts = [np.sum(arr, axis=(1,2)) for arr in self.count_arrays]

        if normalized:

            cts = [arr.astype('float64') for arr in cts]
            matrix = lambda x, y: np.multiply(*np.meshgrid(x, y))
            for arr,(asizes,bsizes) in zip(cts,self.count_sizes):
                arr /= np.sum(matrix(asizes, bsizes))

        return cts

    def ftf_pair_counts(self, fid, normalized):
        # get the diagonal element with index `fid`
        cts = [np.copy(arr[:,fid,fid]) for arr in self.count_arrays]

        if normalized:
            cts = [arr.astype('float64') for arr in cts]
            for arr,(asizes,bsizes) in zip(cts,self.count_sizes):
                arr /= float(asizes[fid] * bsizes[fid])

        return cts

    def jackknife_pair_counts(self, fid, normalized):
        # ignore all pair counts between the field `fld` and any other, and
        # sum the whole matrix

        count_copies = [np.copy(arr) for arr in self.count_arrays]
            
        for arr in count_copies:
            arr[:,fid,:] = 0
            arr[:,:,fid] = 0

        jk_cts = [np.sum(arr, axis=(1,2)) for arr in count_copies]

        if normalized:
            jk_cts = [arr.astype('float64') for arr in jk_cts]

            sizes = [(np.copy(a), np.copy(b)) for a,b in self.count_sizes]

            for a,b in sizes:
                a[fid] = 0
                b[fid] = 0

            matrix = lambda x, y: np.multiply(*np.meshgrid(x, y))

            for cts,(asize,bsize) in zip(jk_cts,sizes):
                cts /= np.sum(matrix(asize,bsize)).astype('float64')

        return jk_cts

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

        # random field data assumed to be the same
        Nfld = self.N_fields

        resample = np.random.randint(
                                0, Nfld, N_trials*Nfld
                                     ).reshape(N_trials, Nfld)

        hists = np.apply_along_axis(np.bincount, 1, resample, minlength=Nfld)

        # Create a matrix M where M(i,j) is product of frequencies of ith, jth
        # fields
        matrix = lambda x: np.multiply(*np.meshgrid(x, x))
        freqs = np.apply_along_axis(matrix, 1, hists)

        matrix = lambda x, y: np.multiply(*np.meshgrid(x, y))
        cts = [arr / np.sum(matrix(asize, bsize)).astype('float64') for arr,(asize,bsize) in zip(self.count_arrays,self.count_sizes)]

        # Dot product over dimensions representing different fields only,
        # leaving num. of trials and radii intact

        cts = [np.tensordot(arr, freqs, axes=((1,2),(1,2))) for arr in cts]

        return cts

    def bootstrap_error(self, N_trials=1000):
        cts = self.bootstrap_pair_counts(N_trials)
        est = self.estimator(*cts)
        return np.std(est, axis=1)

    def estimate(self):
        cts = self.total_pair_counts(normalized=True)
        return self.estimator(*cts)

    def covariance(self):
        if self.error_type is None:
            return None
    
        elif self.error_type == "poisson":
            return None

        cov = np.zeros((self.nbins,self.nbins))

        estimation = self.estimate()
        cts_tot = self.total_pair_counts(normalized=False)

        if self.error_type == "jackknife":

            for fid in range(self.N_fields):
                cts_norm = self.jackknife_pair_counts(fid, normalized=True)
                cts = self.jackknife_pair_counts(fid, normalized=False)

                est_per_field = self.estimator(*cts_norm)
                diff = est_per_field - estimation

                # covariance matrix
                rr_quotient = np.sqrt(cts[-1] / cts_tot[-1].astype('float64'))
                rr_subi, rr_subj = np.meshgrid(rr_quotient,rr_quotient)
                
                xi_subi, xi_subj = np.meshgrid(diff,diff)
                cov += rr_subi*xi_subi*rr_subj*xi_subj

        elif self.error_type == "ftf":
            for fid in range(self.N_fields):
                cts_norm = self.ftf_pair_counts(fid, normalized=True)
                cts = self.ftf_pair_counts(fid, normalized=False)

                est_per_field = self.estimator(*cts_norm)
                diff = est_per_field - estimation

                # covariance matrix
                rr_quotient = np.sqrt(cts[-1] / cts_tot[-1].astype('float64'))
                rr_subi, rr_subj = np.meshgrid(rr_quotient,rr_quotient)
                
                xi_subi, xi_subj = np.meshgrid(diff,diff)
                cov += rr_subi*xi_subi*rr_subj*xi_subj

            cov /= float(self.N_fields - 1)

        return cov

    def error(self):
        if self.error_type is None:
            return None
        elif self.error_type == "poisson":
            cts_tot = self.total_pair_counts(normalized=False)

            estimation = self.estimate()
            with np.errstate(divide='ignore', invalid='ignore'):
                error = np.divide(1 + estimation,np.sqrt(cts_tot[0]))
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
        cts = self.total_pair_counts(normalized=False)
        est = self.estimate()
        err = self.error()

        if err is None:
            err = [None] * len(est)

        # calc column widths for pair counts
        ct_widths = [max([len(str(count)) for count in ct] + [4]) + 2 for ct in cts]

        lines = [ "\n" ]

        if len(cts) == 3:
            labels = ["Bin L", "Bin R", "DD", "DR", "RR", "Estimator", "Error"]
        elif len(cts) == 4:
            labels = ["Bin L", "Bin R", "D1D2", "D1R2", "D2R1", "R1R2", "Estimator", "Error"]

        sizes =  [9,       9,   ] + ct_widths + [11,          10]

        header = [label.ljust(size) for label, size in zip(labels, sizes)]

        lines.append(''.join(header))

        # add units to radii column, e.g. degrees
        lines.append("{:^14}".format("({:s})".format(self.radii_units)))


        lines.append("-" * 79)

        ########## end of table header ##########

        for rl, ru, estv, errv, ct in \
                zip(r_lower, r_upper, est, err, zip(*cts)):

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
                 ] + ["{:<{:d}}".format(c, w) for c,w in zip(ct,ct_widths)] + \
                [
                    "{:<11}".format(estv_s),
                    "{:<10}".format(errv_s),
                 ]

            lines.append(''.join(s))

        return '\n'.join(lines)

class twopoint_autocorr_results(twopoint_results):
    def __init__(self, dd, dr, rr, dtree, rtree, estimator, radii):
        self.dd = dd
        self.dr = dr
        self.rr = rr

        # TODO more consistent way to determine N_fields
        self.N_fields = dtree.N_fields

        self.count_arrays = [self.dd, self.dr, self.rr]
        self.count_sizes = [
                (np.copy(dtree.field_sizes), np.copy(dtree.field_sizes)),
                (np.copy(dtree.field_sizes), np.copy(rtree.field_sizes)),
                (np.copy(rtree.field_sizes), np.copy(rtree.field_sizes)),
                            ]
        super(twopoint_autocorr_results, self).__init__(estimator, radii)

class twopoint_crosscorr_results(twopoint_results):
    def __init__(self, d1d2, d1r2, d2r1, r1r2, dtree_1, dtree_2, rtree_1,
                 rtree_2, estimator, radii):
        self.d1d2 = d1d2
        self.d1r2 = d1r2
        self.d2r1 = d2r1
        self.r1r2 = r1r2

        # TODO more consistent way to determine N_fields
        self.N_fields = dtree_1.N_fields

        self.count_arrays = [self.d1d2, self.d1r2, self.d2r1, self.r1r2]
        self.count_sizes = [
                (np.copy(dtree_1.field_sizes), np.copy(dtree_2.field_sizes)),
                (np.copy(dtree_1.field_sizes), np.copy(rtree_2.field_sizes)),
                (np.copy(dtree_2.field_sizes), np.copy(rtree_1.field_sizes)),
                (np.copy(rtree_1.field_sizes), np.copy(rtree_2.field_sizes)),
                            ]

        super(twopoint_crosscorr_results, self).__init__(estimator, radii)

