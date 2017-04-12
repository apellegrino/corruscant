import os
from ctypes import *
import numpy as np

'''
# Cosmology library for computing comoving distances quickly, to be used with
# scipy.integrate

coslib = CDLL(os.path.abspath("cosmology.so"))
coslib.f.argtypes = (c_int,c_double)
coslib.f.restype = c_double
'''

# To define a ctypes struct that contains pointers to the same struct type, we
# must declare the class first and then add references
class node(Structure):
    pass

'''
# python never manipulates nodes directly so we may not need
# fields for the ctpyes structs
node._fields_ = [

    ("x", c_double),
    ("y", c_double),
    ("z", c_double),
    ("left_child", POINTER(node) ),
    ("left_child", POINTER(node) ),
    ]
'''

# accessing ctypes data:
# print cast(self.data_tree.y, POINTER(c_double)).contents


class kdtree(Structure):
    _fields_ = [
        ("root",POINTER(node)),
        ("size",c_int),
        ("memsize",c_int),
        ("x",POINTER(c_double)),
        ("y",POINTER(c_double)),
        ("z",POINTER(c_double)),
        ]

kdlib = CDLL(os.path.abspath("libkdtree.so"))

kdlib.tree_construct.restype = kdtree
kdlib.tree_construct.argtypes = [
                            c_int,
                            POINTER(c_double),
                            POINTER(c_double),
                            POINTER(c_double),
                            ]

kdlib.destroy.restype = None
kdlib.destroy.argtypes = [ POINTER(kdtree) ]

kdlib.pair_count.restype = c_longlong
kdlib.pair_count.argtypes = [
                            kdtree,
                            POINTER(c_double),
                            POINTER(c_double),
                            POINTER(c_double),
                            c_int,
                            c_double,
                            ]

def _unpack(data):
    #return tuple([np.require(row, requirements=reqs) for row in data])
    return tuple([np.copy(row) for row in data])

def _make_tree(data):

    # tree_construct is inefficient for sorted data
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

    # sort data and random in dim
    data = data[:,np.argsort(data[dim,:])]
    random = random[:,np.argsort(random[dim,:])]


    for i in range(N):
        start = i * nd/N + min(i,nd % N)
        stop = start + nd/N + (nd % N > i)

        # exclude from start to stop
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
        yield np.copy(data_sub), np.copy(random_sub)

def landy_szalay(dd,dr,rr,f):
    return ( f*f*np.array(dd) - 2*f*np.array(dr) + np.array(rr) ) / np.array(rr)

def _query_tree(data, tree, radius, num_threads):

    data_x, data_y, data_z = _unpack(data)
    nd = data.shape[1]

    count = kdlib.pair_count(tree,
                             data_x.ctypes.data_as(POINTER(c_double)),
                             data_y.ctypes.data_as(POINTER(c_double)),
                             data_z.ctypes.data_as(POINTER(c_double)),
                             nd,radius,num_threads)

    return count

def _pc(data_tree, X_data, rand_tree, X_rand, radii, num_threads):
    # Num. random points
    nr = rand_tree.size
    # Num. data points
    nd = data_tree.size

    # Data-data pair counts as function of radius
    dd_array = []
    # Data-random pair counts
    dr_array = []
    # Random-random pair counts
    rr_array = []

    for r in radii:
        dd = _query_tree(X_data, data_tree, r, num_threads)
        dr = _query_tree(X_data, rand_tree, r, num_threads)
        rr = _query_tree(X_rand, rand_tree, r, num_threads)

        dd_array.append(dd)
        dr_array.append(dr)
        rr_array.append(rr)

    return np.diff(dd_array), np.diff(dr_array), np.diff(rr_array)

def pair_counts(data, rand, radii, xi_error_type=None,
                xi_estimator_type='landy-szalay', N_error=10, num_threads=32):

    data_tree = _make_tree(data)
    rand_tree = _make_tree(rand)

    dd_array, dr_array, rr_array = _pc(data_tree, data, rand_tree, rand, radii, num_threads)
    f = float(rand_tree.size)/data_tree.size
    ls = landy_szalay(dd_array,dr_array,rr_array,f)

    kdlib.destroy(data_tree)
    kdlib.destroy(rand_tree)

    # Error as specified by error type:
    error = []

    if xi_error_type == 'jackknife':
        gen = _jackknife(data, rand, N_error)
        error = np.zeros(len(radii)-1)

        for data_subset, rand_subset in gen:
            f_sub = float(rand_subset.shape[1])/float(data_subset.shape[1])
            data_tree = _make_tree(data_subset)
            rand_tree = _make_tree(rand_subset)

            dd_err, dr_err, rr_err = _pc(data_tree, data_subset, rand_tree, rand_subset, radii, num_threads)
            ls_sub = landy_szalay(dd_err,dr_err,rr_err,f_sub)

            diff = ls - ls_sub
            error += np.divide(dr_err.astype('float64'),dr_array)*diff*diff

            kdlib.destroy(data_tree)
            kdlib.destroy(rand_tree)

    elif xi_error_type == 'field-to-field':
        gen = _ftf(data, rand, N_error)
        error = np.zeros(len(radii)-1)

        for data_subset, rand_subset in gen:
            f_sub = float(rand_subset.shape[1])/float(data_subset.shape[1])
            data_tree = _make_tree(data_subset)
            rand_tree = _make_tree(rand_subset)

            dd_err, dr_err, rr_err = _pc(data_tree, data_subset, rand_tree, rand_subset, radii, num_threads)
            ls_sub = landy_szalay(dd_err,dr_err,rr_err,f_sub)

            diff = ls - ls_sub
            error += np.divide(dr_err.astype('float64'),dr_array)*diff*diff

            kdlib.destroy(data_tree)
            kdlib.destroy(rand_tree)

        error /= float(N_error - 1)

    elif xi_error_type == 'poisson':
        error = np.divide(1 + ls,np.sqrt(dd_array))

    return dd_array, dr_array, rr_array, ls, error
