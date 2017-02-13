import os
from ctypes import *
import numpy as np
from mpi4py import MPI

'''
coslib = CDLL(os.path.abspath("cosmology.so"))
coslib.f.argtypes = (c_int,c_double)
coslib.f.restype = c_double
'''

# To define a ctypes struct that contains pointers to the same struct type, we
# must declare the class first and then add references
class node(Structure):
    pass

node._fields_ = [
    # python never manipulates nodes directly so we may not need
    # fields for the ctpyes structs

    #("x", c_double),
    #("y", c_double),
    #("z", c_double),
    #("left_child", POINTER(node) ),
    #("left_child", POINTER(node) ),
    ]

class kdtree(Structure):
    _fields_ = [
    	("root",POINTER(node)),
    	("size",c_int),
    	("x",POINTER(c_double)),
    	("y",POINTER(c_double)),
    	("z",POINTER(c_double)),
        ]

kdlib = CDLL(os.path.abspath("driver.so"))

if MPI._sizeof(MPI.Comm) == sizeof(c_int):
    MPI_Comm = c_int
else:
    MPI_Comm = c_void_p

kdlib.tree_construct.restype = kdtree
kdlib.tree_construct.argtypes = [
    							c_int,
    							POINTER(c_double),
    							POINTER(c_double),
    							POINTER(c_double),
                                ]
'''
kdlib.landy_szalay.restype = c_double
kdlib.landy_szalay.argtypes = [
    						kdtree,
    						kdtree,
    						c_double ]
'''
kdlib.two_point_correlation.restype = c_longlong
kdlib.two_point_correlation.argtypes = [
    						kdtree,
    						POINTER(c_double),
    						POINTER(c_double),
    						POINTER(c_double),
    						c_int,
    						c_double,
                            MPI_Comm,
                            ]

class twopoint_constructor:
    def __init__(self, X_data, X_random):
    	try:
    		self.X_data = np.array(X_data)
    		self.X_random = np.array(X_random)
    	except:
    		raise

    def _unpack(self,x,y,z):
    	return (x,y,z)

    def construct(self):
    	data_x, data_y, data_z = self._unpack(*self.X_data)

    	self.data_tree = kdlib.tree_construct(c_int(self.X_data.shape[1]),
    					data_x.ctypes.data_as(POINTER(c_double)),
    					data_y.ctypes.data_as(POINTER(c_double)),
    					data_z.ctypes.data_as(POINTER(c_double)) )
    	print cast(self.data_tree.y, POINTER(c_double)).contents
    	rand_x, rand_y, rand_z = self._unpack(*self.X_random)

    	self.rand_tree = kdlib.tree_construct(c_int(self.X_random.shape[1]),
    					rand_x.ctypes.data_as(POINTER(c_double)),
    					rand_y.ctypes.data_as(POINTER(c_double)),
    					rand_z.ctypes.data_as(POINTER(c_double)) )

    def landy_szalay(self, radii):
    	data_x, data_y, data_z = self._unpack(*self.X_data)
    	rand_x, rand_y, rand_z = self._unpack(*self.X_random)
    	nr = self.X_random.shape[1]
    	nd = self.X_data.shape[1]
    	f = float(nr)/nd

    	out = []

        comm_ptr = MPI._addressof(MPI.COMM_WORLD)
        comm_val = MPI_Comm.from_address(comm_ptr)

    	for r in radii:
    		dd = kdlib.two_point_correlation(self.data_tree,
    								data_x.ctypes.data_as(POINTER(c_double)),
    								data_y.ctypes.data_as(POINTER(c_double)),
    								data_z.ctypes.data_as(POINTER(c_double)),
    								nd,r,comm_val)
    		dr = kdlib.two_point_correlation(self.data_tree,
    								rand_x.ctypes.data_as(POINTER(c_double)),
    								rand_y.ctypes.data_as(POINTER(c_double)),
    								rand_z.ctypes.data_as(POINTER(c_double)),
    								nr,r,comm_val)
    		rr = kdlib.two_point_correlation(self.rand_tree,
    								rand_x.ctypes.data_as(POINTER(c_double)),
    								rand_y.ctypes.data_as(POINTER(c_double)),
    								rand_z.ctypes.data_as(POINTER(c_double)),
    								nr,r,comm_val)

    		out.append((f*f*dd-2*f*dr+rr)/float(rr))

    	return out
