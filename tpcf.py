from ctypes import *
from os.path import abspath, dirname
import numpy as np

class node(Structure):
    pass

class kdtree(Structure):
    _fields_ = [
        ("root",POINTER(node)),
        ("size",c_int),
        ("memsize",c_int),
        ("x",POINTER(c_double)),
        ("y",POINTER(c_double)),
        ("z",POINTER(c_double)),
        ]

path_here = abspath(__file__)
path_dir = dirname(path_here)
kdlib = CDLL("%s/bin/libkdtree.so" % path_dir)

kdlib.tree_construct.restype = kdtree
kdlib.tree_construct.argtypes = [
                            c_int, # array size
                            POINTER(c_double), # x array
                            POINTER(c_double), # y array
                            POINTER(c_double), # z array
                            ]

kdlib.destroy.restype = None
kdlib.destroy.argtypes = [ POINTER(kdtree) ] # tree to destroy

kdlib.pair_count.restype = c_longlong # num. of pair counts
kdlib.pair_count.argtypes = [
                            kdtree, # tree to query
                            POINTER(c_double), # x in array of querying pts.
                            POINTER(c_double), # y in array of querying pts.
                            POINTER(c_double), # z in array of querying pts.
                            c_int, # array size
                            c_double, # radius
                            ]

def _unpack(data):
    #return tuple([np.require(row, requirements=reqs) for row in data])
    return tuple([np.copy(row) for row in data])

def _make_tree(data):

    # tree_construct is inefficient for sorted data due to quicksort.
    # Maybe use mergesort instead
    np.random.shuffle(data.T)
    data_x, data_y, data_z = _unpack(data)

    tree = kdlib.tree_construct(
                    c_int(data.shape[1]),
                    data_x.ctypes.data_as(POINTER(c_double)),
                    data_y.ctypes.data_as(POINTER(c_double)),
                    data_z.ctypes.data_as(POINTER(c_double))
                    )

    return tree

# exclude 1/N th of the data and random sets in a dimension
def _jackknife(data, random, N, dim=0):
    nd = data.shape[1]

    # sort data and random in the (dim) dimension
    # X=0, Y=1, Z=2
    data = data[:,np.argsort(data[dim,:])]
    random = random[:,np.argsort(random[dim,:])]


    for i in range(N):
        start = i * nd/N + min(i,nd % N)
        stop = start + nd/N + (nd % N > i)

        # exclude the range from start to stop
        data_sub = np.concatenate([data[:,:start], data[:,stop:]], axis=1)

        # counter bias due to splitting based on data position
        if start <= 0:
            data_min = -np.inf
        else:
            data_min = (data[dim,start]+data[dim,start-1])/2

        if stop >= nd:
            data_max = np.inf
        else:
            data_max = (data[dim,stop-1]+data[dim,stop])/2

        random_sub = np.concatenate([random[:,random[dim] < data_min],
                                    random[:,random[dim] > data_max]], axis=1)

        yield data_sub, random_sub

def _ftf(data, random, N=10, dim=0):
    nd = data.shape[1]

    # sort data and random in dim
    data = data[:,np.argsort(data[dim,:])]
    random = random[:,np.argsort(random[dim,:])]


    for i in range(N):
        start = i * nd/N + min(i,nd % N)
        stop = start + nd/N + (nd % N > i)

        # include from start to stop
        data_sub = data[:,start:stop]

        # counter bias due to splitting based on data position
        if start <= 0:
            data_min = -np.inf
        else:
            data_min = (data[dim,start]+data[dim,start+1])/2

        if stop >= nd:
            data_max = np.inf
        else:
            data_max = (data[dim,stop-1]+data[dim,stop-2])/2

        random_sub = random[:,random[dim] > data_min]
        random_sub = random_sub[:,random_sub[dim] < data_max]
        yield data_sub, random_sub

def _query_tree(tree, data, radius, num_threads):

    data_x, data_y, data_z = _unpack(data)
    nd = data.shape[1]

    count = kdlib.pair_count(tree,
                             data_x.ctypes.data_as(POINTER(c_double)),
                             data_y.ctypes.data_as(POINTER(c_double)),
                             data_z.ctypes.data_as(POINTER(c_double)),
                             nd,radius,num_threads)

    return count

class InputError:
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

def validate_array(arr):
    if not type(arr) is np.ndarray:
        raise InputError("Array must be a Numpy array")

    if not (arr.shape[0] == 3 and arr.ndim == 2):
        raise InputError("Array must be of shape (3, N). The provided array has shape %s" % str(arr.shape))

def est_landy_szalay(dd,dr,rr,dsize,rsize):
    f = float(rsize)/dsize
    return ( f*f*np.array(dd) - 2*f*np.array(dr) + np.array(rr) ) / np.array(rr)

def est_hamilton(dd,dr,rr):
    return np.divide( np.multiply(dd,rr).astype("float64"),
                      np.multiply(dr,dr).astype("float64") ) - 1.

def est_standard(dd,dr,dsize,rsize):
    return float(rsize)/dsize * np.divide( dd.astype("float64"), dr ) - 1.

def pair_counts(data, rand, radii, xi_error_type=None,
                xi_estimator_type="landy-szalay", N_error=10, num_threads=4):

    est_types = [None, "landy-szalay", "hamilton", "standard"]
    if not xi_estimator_type in est_types:
        raise InputError("Estimator type for Xi %s not valid" % xi_estimator_type)

    validate_array(data)
    validate_array(rand)

    data_tree = _make_tree(data)
    rand_tree = _make_tree(rand)

    dd_array = np.diff([_query_tree(data_tree, data, r, num_threads) for r in radii])
    dr_array = np.diff([_query_tree(rand_tree, data, r, num_threads) for r in radii])

    if not xi_estimator_type == "standard":
        rr_array = np.diff([_query_tree(rand_tree, rand, r, num_threads) for r in radii])
    else:
        rr_array = None

    if xi_estimator_type == "landy-szalay":
        est = est_landy_szalay(dd_array,dr_array,rr_array,data.size,rand.size)
    elif xi_estimator_type == "hamilton":
        est = est_hamilton(dd_array,dr_array,rr_array)
    elif xi_estimator_type == "standard":
        est = est_standard(dd_array,dr_array,data.size,rand.size)

    kdlib.destroy(data_tree)
    kdlib.destroy(rand_tree)

    #error = np.array([np.inf] * len(radii-1))
    error = None

    if xi_error_type == "jackknife" or xi_error_type == "field-to-field":
        if xi_error_type == "jackknife":
            gen = _jackknife(data, rand, N_error)
        elif xi_error_type == "field-to-field":
            gen = _ftf(data, rand, N_error)

        error = np.zeros(len(radii)-1)

        for data_subset, rand_subset in gen:
            data_tree = _make_tree(data_subset)
            rand_tree = _make_tree(rand_subset)

            dd_err = np.diff([_query_tree(data_tree, data, r, num_threads) for r in radii])
            dr_err = np.diff([_query_tree(rand_tree, data, r, num_threads) for r in radii])
            rr_err = np.diff([_query_tree(rand_tree, rand, r, num_threads) for r in radii])
            
            if xi_estimator_type == "landy-szalay":
                est_sub = est_landy_szalay(dd_err,dr_err,rr_err,data_subset.size,rand_subset.size)
            elif xi_estimator_type == "hamilton":
                est_sub = est_hamilton(dd_err,dr_err,rr_err)
            elif xi_estimator_type == "standard":
                est_sub = est_standard(dd_err,dr_err,data_subset.size,rand_subset.size)

            diff = est - est_sub
            error += np.divide(dr_err.astype("float64"),dr_array)*diff*diff

            kdlib.destroy(data_tree)
            kdlib.destroy(rand_tree)

        if xi_error_type == "field-to-field":
            error /= float(N_error - 1)

    elif xi_error_type == "poisson":
        with np.errstate(divide='ignore', invalid='ignore'):
            error = np.divide(1 + est,np.sqrt(dd_array))
            error[np.isneginf(error)] = np.nan
            error[np.isinf(error)] = np.nan

    output = {  
                "xi_error_type":xi_error_type,
                "radii":radii,
                "DD":dd_array,
                "DR":dr_array,
                "RR":rr_array,
                "estimator":est,
                "error":error,
            }
    return output
