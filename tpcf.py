import ctypes
from os.path import abspath, dirname
import numpy as np

class node(ctypes.Structure):
    pass

class array3d(ctypes.Structure):
    _fields_ = [
        ("x",ctypes.POINTER(ctypes.c_double)),
        ("y",ctypes.POINTER(ctypes.c_double)),
        ("z",ctypes.POINTER(ctypes.c_double)),
        ("fields",ctypes.POINTER(ctypes.c_int)),
        ("num_fields",ctypes.c_int),
        ("size",ctypes.c_int),
        ]

class argarray3d(ctypes.Structure):
    _fields_ = [
        ("x",ctypes.POINTER(ctypes.c_int)),
        ("y",ctypes.POINTER(ctypes.c_int)),
        ("z",ctypes.POINTER(ctypes.c_int)),
        ("size",ctypes.c_int),
        ]

class kdtree(ctypes.Structure):
    _fields_ = [
        ("node_data",ctypes.POINTER(node)),
        ("size",ctypes.c_int),
        ("memsize",ctypes.c_int),
        ("data",array3d),
        ("arg_data",argarray3d),
        ]

path_here = abspath(__file__)
path_dir = dirname(path_here)
kdlib = ctypes.CDLL("%s/bin/libkdtree.so" % path_dir)

kdlib.tree_construct.restype = kdtree
kdlib.tree_construct.argtypes = [array3d]

kdlib.form_array.restype = array3d
kdlib.form_array.argtypes = [
                            ctypes.POINTER(ctypes.c_double), # x
                            ctypes.POINTER(ctypes.c_double), # y
                            ctypes.POINTER(ctypes.c_double), # z
                            ctypes.POINTER(ctypes.c_int), # field ids
                            ctypes.c_int, #size
                            ]

kdlib.destroy.restype = None
kdlib.destroy.argtypes = [kdtree]

#                                    array w/ numbers of pair counts
kdlib.pair_count_jackknife.restype = np.ctypeslib.ndpointer(dtype=ctypes.c_longlong, shape=(255,))

kdlib.pair_count_jackknife.argtypes = [
                            kdtree, # tree to query
                            array3d, # data to query with
                            ctypes.c_double, # radius
                            ctypes.c_int, # num_threads
                            ]

kdlib.pair_count_ftf.restype = np.ctypeslib.ndpointer(dtype=ctypes.c_longlong, shape=(255,))

kdlib.pair_count_ftf.argtypes = [
                            kdtree, # tree to query
                            array3d, # data to query with
                            ctypes.c_double, # radius
                            ctypes.c_int, # num_threads
                            ]

kdlib.pair_count_noerr.restype = np.ctypeslib.ndpointer(dtype=ctypes.c_longlong, shape=(255,))

kdlib.pair_count_noerr.argtypes = [
                            kdtree, # tree to query
                            array3d, # data to query with
                            ctypes.c_double, # radius
                            ctypes.c_int, # num_threads
                            ]
# class Dataset:
#     class Data:
#         def __init__(self, type, dim1, dim2, dim3=None):
# 
#             if type == 'cartesian':
#                 self.x = dim1
#                 self.y = dim2
#                 self.z = dim3
#             else if type == '2d_spherical_polar':
#                 self.az = 
# 
#     def __init__(self):
# 
#         self.datatype = 'none'
#         self.data = None
#         self.rand = None
#         self.cos = None
#         self.error_type = None
#         self.error_bounds
# 
#     def add_cartesian(self, data, random):
#         self.data = data
#         self.rand = random
#         self.datatype = 'cartesian'
# 
#     def add_2d_spherical(self, data, random):
#         self.data = data
#         self.rand = random
#         self.datatype = '2d_spherical'
# 
#     def add_3d_spherical(self, data, random):
#         self.data = data
#         self.rand = random
#         self.datatype = '3d_spherical'
# 
#     def set_cosmology(self, cos):
#         self.cos = cos
# 
#     def coord_transform(self):
#         if self.datatype = '3d_spherical':
#             
# 
#     def sph_to_cart(self):
#         if self.cos is None:
#             raise InputError("Provide a cosmology with set_cosmology() before converting to cartesian coordinates")
# 
#         r = self.cos.comoving_distance(self.coord_redshift)
# 
#         if self.spherical_type == 'latitude':
#             theta = np.pi / 2 - self.coord_el
#         else: # polar
#             theta = self.coord_el
# 
#         self.coord_x = r * np.sin(theta) * np.sin(phi)
#         self.coord_y = r * np.sin(theta) * np.cos(phi)
#         self.coord_z = r * np.cos(theta)

# input: ra, dec in degrees, redshift
# output: x, y, z
def sph_to_cart(data):
    ra = data[0]
    dec = data[1]
    redshift = data[2]

    from astropy import cosmology

    H_0 = 100.0
    Om0 = .237
    Ode0 = .763

    cos = cosmology.LambdaCDM(H_0, Om0, Ode0)
    r = cos.comoving_distance(redshift)

    x = r * np.cos(np.pi / 180. * dec) * np.cos(np.pi / 180. * ra)
    y = r * np.cos(np.pi / 180. * dec) * np.sin(np.pi / 180. * ra)
    z = r * np.sin(np.pi / 180. * dec)

    return np.array([x,y,z])

def _unpack(data):
    return tuple([np.copy(row) for row in data])

# array breaks when returned from this function ???
def make_clike_array(data, fields=None):
    # tree_construct is inefficient for sorted data due to quicksort.
    # Maybe use mergesort instead
    #np.random.shuffle(data.T)
    data_x, data_y, data_z = _unpack(data)

    if fields is None:
        fields = np.zeros_like(data_z)
    # numpy will use 64-bit ints unless otherwise specified
    # type change must happen on ITS OWN LINE!
    fields = fields.astype('int32')

    array = kdlib.form_array(
                        data_x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                        data_y.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                        data_z.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                        fields.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
                        ctypes.c_int(data.shape[1])
                        )
    return array

def _make_tree(data, fields, N_fields=1):
    # tree_construct is inefficient for sorted data due to quicksort.
    # Maybe use mergesort instead
    #np.random.shuffle(data.T)
    data_x, data_y, data_z = _unpack(data)

    if fields is None:
        fields = np.zeros_like(data_x)
    # numpy will use 64-bit ints unless otherwise specified
    # type change must happen on ITS OWN LINE!
    fields = fields.astype('int32')

    array = kdlib.form_array(
                        data_x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                        data_y.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                        data_z.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                        fields.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
                        ctypes.c_int(N_fields),
                        ctypes.c_int(data.shape[1]),
                        )
    #array = make_clike_array(data,fields)

    tree = kdlib.tree_construct(array)

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

def cart_to_sph(data):
    x = data[0]
    y = data[1]
    z = data[2]

    offset = -np.sign(x)*(np.sign(y)+1)/2.+1

    ra = np.arctan2(y,x)/np.pi * 180.
    if ra < 0:
        ra = ra + 360.

    dec = np.sign(z) * np.arctan(np.sqrt(z*z/(x*x+y*y))) / np.pi * 180.

    return ra, dec

# user-defined boundaries for jackknife error
# bounds is an array of boundaries, where each boundary defines a min and max
# value for each coordinate e.g. [min_ra, max_ra, min_dec, max_dec]
def _jackknife_bounds(data, random, bounds):
    for b in bounds:
        print b
        
        data_sub = []
        random_sub = []
    
        for point in data.T:
            ra, dec = cart_to_sph(point)
            if not b[0] < ra < b[1]:
                data_sub.append(point)
                continue
            if not b[2] < dec < b[3]:
                data_sub.append(point)
                continue

        for point in random.T:
            ra, dec = cart_to_sph(point)
            if not b[0] < ra < b[1]:
                random_sub.append(point)
                continue
            if not b[2] < dec < b[3]:
                random_sub.append(point)
                continue

        data_sub = np.copy(np.array(data_sub).T)
        random_sub = np.copy(np.array(random_sub).T)

        print "%d out of %d" % (data_sub.shape[1], data.shape[1])
        yield data_sub, random_sub

# user-defined boundaries for field-to-field error
def _ftf_bounds(data, random, bounds):
    for b in bounds:
        
        data_sub = []
        random_sub = []
    
        for point in data.T:
            ra, dec = cart_to_sph(point)
            if not bounds[0][0] < ra < bounds[0][1]: continue
            if not bounds[1][0] < dec < bounds[1][1]: continue
            data_sub.append(point)

        for point in random.T:
            ra, dec = cart_to_sph(point)
            if not bounds[0][0] < ra < bounds[0][1]: continue
            if not bounds[1][0] < dec < bounds[1][1]: continue
            random_sub.append(point)

        data_sub = np.copy(np.array(data_sub).T)
        random_sub = np.copy(np.array(random_sub).T)

        print data_sub.shape, random_sub.shape
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

def _query_tree(tree, data, radius, num_threads, errtype, fields=None, N_fields=0):
    data_x, data_y, data_z = _unpack(data)
    if fields is None:
        fields = np.zeros_like(data_x)
    # numpy will use 64-bit ints unless otherwise specified
    # type change must happen on ITS OWN LINE!
    fields = fields.astype('int32')

    array = kdlib.form_array(
                        data_x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                        data_y.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                        data_z.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                        fields.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
                        ctypes.c_int(N_fields),
                        ctypes.c_int(data.shape[1]),
                        )

    #array = make_clike_array(data, fields)
    counts = None
    if errtype == 'jackknife':
        counts = kdlib.pair_count_jackknife(tree,array,radius,num_threads)
    elif errtype == 'field-to-field':
        counts = kdlib.pair_count_ftf(tree,array,radius,num_threads)
    else:
        counts = kdlib.pair_count_noerr(tree,array,radius,num_threads)

    return counts[:N_fields+1]

class InputError:
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

def validate_array(arr):
    if not type(arr) is np.ndarray:
        arr = np.array(arr)

    if not (arr.ndim == 2):
        raise InputError("Array must be two-dimensional (i.e. a list of points)")

    # try to make array of shape (3, N) if shape is (N, 3)
    if (arr.shape[1] == 3):
        arr = arr.T

    if not(arr.shape[0] == 3):
        raise InputError("Array must be of shape (3, N) or (N, 3). The provided array has shape %s" % str(arr.shape))

    return arr

def est_landy_szalay(dd,dr,rr,dsize,rsize):
    f = float(rsize)/dsize
    return ( f*f*np.array(dd) - 2*f*np.array(dr) + np.array(rr) ) / np.array(rr)

def est_hamilton(dd,dr,rr):
    return np.divide( np.multiply(dd,rr).astype("float64"),
                      np.multiply(dr,dr).astype("float64") ) - 1.

def est_standard(dd,dr,dsize,rsize):
    return float(rsize)/dsize * np.divide( dd.astype("float64"), dr ) - 1.

def pair_counts(data, rand, radii, 
                data_fields=None, rand_fields=None,
                xi_estimator_type="landy-szalay",
                xi_error_type=None, N_error=10, num_threads=4,
                bounds=None, N_fields=0, datatype='cartesian'):
    """Given a set of 3D cartesian data points and random points, calculate the
    estimated two-point correlation function with error estimation.

    Arguments:
    data (np.ndarray) -- a numpy array of data points with shape (3, N_data).
    rand (np.ndarray) -- a numpy array of random points with shape (3, N_random).
    radii (array-like) -- an array-like of floats which define the radius bin sizes.

    Keyword arguments:
    xi_estimator_type -- the type of estimator to use for the correlation
    function. (default "landy-szalay") Possible values in order of decreasing
    speed: "standard" > "landy-szalay" = "hamilton"

    xi_error_type -- a string defining what type of error to calculate. (default None)
        Possible values in order of decreasing speed:
            None > "poisson" > "field-to-field" > "jackknife"

    N_error -- an integer describing the number of bins to use when calculating
    jackknife or field-to-field errors. Lower is faster, higher is more
    accurate. (default 10)

    num_threads -- an integer describing the number of threads that the C code
    will create. For max performance, set this to the number of logical cores
    on your machine. (default 4)

    Return:
    A dictionary of pair counts, the requested estimator values for input
    radii, the errors of estimator values if specified, and the input radii
    list and estimator and error types.
    """
    est_types = ["landy-szalay", "hamilton", "standard"]
    if not xi_estimator_type in est_types:
        raise InputError("Estimator type for Xi %s not valid" % xi_estimator_type)

    err_types = [None, "poisson", "field-to-field", "jackknife"]
    if not xi_error_type in err_types:
        raise InputError("Estimator error type %s not valid" % xi_error_type)

    data = validate_array(data)
    rand = validate_array(rand)

    if datatype == 'spherical':
        data = sph_to_cart(data)
        rand = sph_to_cart(rand)

    data_tree = _make_tree(data, data_fields, N_fields)
    rand_tree = _make_tree(rand, rand_fields, N_fields)

    dd_array = np.diff([_query_tree(data_tree, data, r, num_threads, xi_error_type, data_fields, N_fields) for r in radii], axis=0)
    dr_array = np.diff([_query_tree(rand_tree, data, r, num_threads, xi_error_type, data_fields, N_fields) for r in radii], axis=0)

    rr_array = None
    if not xi_estimator_type == "standard":
        rr_array = np.diff([_query_tree(rand_tree, rand, r, num_threads, xi_error_type, rand_fields, N_fields) for r in radii], axis=0)

    dd_total = dd_array[:,0]
    dr_total = dr_array[:,0]
    rr_total = rr_array[:,0]
    if xi_estimator_type == "landy-szalay":
        est = est_landy_szalay(dd_total,dr_total,rr_total,data.size,rand.size)
    elif xi_estimator_type == "hamilton":
        est = est_hamilton(dd_total,dr_total,rr_total)
    elif xi_estimator_type == "standard":
        est = est_standard(dd_total,dr_total,data.size,rand.size)
    error = None

    if xi_error_type == "jackknife" or xi_error_type == "field-to-field":
        error = np.zeros(len(radii)-1)
        
        dd_err = dd_array[:,1:N_fields+1]
        dr_err = dr_array[:,1:N_fields+1]
        rr_err = rr_array[:,1:N_fields+1]

        for sub_i in range(N_fields):
            this_dd_err = dd_err[:,sub_i]
            print this_dd_err
            this_dr_err = dr_err[:,sub_i]
            this_rr_err = rr_err[:,sub_i]

            id = sub_i + 1
            if xi_error_type == "jackknife":
                nd = data_fields[data_fields != id].size
                nr = rand_fields[rand_fields != id].size
            elif xi_error_type == "field-to-field":
                nd = data_fields[data_fields == id].size
                nr = rand_fields[rand_fields == id].size
       
            print nd, nr
            if xi_estimator_type == "landy-szalay":
                est_sub = est_landy_szalay(this_dd_err,this_dr_err,this_rr_err,nd,nr)
            elif xi_estimator_type == "hamilton":
                est_sub = est_hamilton(this_dd_err,this_dr_err,this_rr_err)
            elif xi_estimator_type == "standard":
                est_sub = est_standard(this_dd_err,this_dr_err,nd,nr)

            diff = est - est_sub
            error += np.divide(this_dr_err.astype("float64"),dr_total)*diff*diff

        if xi_error_type == "field-to-field":
            error /= float(N_error - 1)

    elif xi_error_type == "poisson":
        with np.errstate(divide='ignore', invalid='ignore'):
            error = np.divide(1 + est,np.sqrt(dd_array))
            error[np.isneginf(error)] = np.nan
            error[np.isinf(error)] = np.nan

    kdlib.destroy(data_tree)
    kdlib.destroy(rand_tree)

    if error is not None: error = np.sqrt(error)
    output = {  
                "radii":radii,
                "DD":dd_total,
                "DR":dr_total,
                "RR":rr_total,
                "estimator":est,
                "error": error,
                "xi_error_type":xi_error_type,
            }
    return output
