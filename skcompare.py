from mpi4py import MPI
import twopoint
import math
from scipy import integrate
import numpy as np
from sklearn import neighbors
import time

# constants
c = 299792.4580 # km/s
H_0 = 100.0 #h^-1 km/s/Mpc


def integrand(a, omega_lambda, omega_k, omega_matter, omega_rel):
	return 1.0/math.sqrt(omega_rel + omega_matter*a + omega_k*a*a + 
										omega_lambda*a*a*a*a)

# return floats from string of float literals separated by spaces
def parse_line(string):
	return [float(i) for i in filter(None, string.strip('\r\n').split(' '))]

def readfile(path):
	x = []
	y = []
	z = []
	with open(path, "r") as f:
		for line in f:
			ra, dec, redshift = parse_line(line)

			# convert from degrees to radians
			ra *= math.pi/180.0
			dec = math.pi*(.5 - dec/180.0)

			#integral, error = integrate.quad(coslib.f,1/(1+redshift),1.0)
			integral, error = integrate.quad(integrand,1/(1+redshift),1.0,
												args = (0.7, 0., 0.3, 0.))
			d_comoving = c / H_0 * integral #Mpc

			# transform to cartesian coordinates
			xv = d_comoving * math.sin(dec) * math.sin(ra)
			yv = d_comoving * math.sin(dec) * math.cos(ra)
			zv = d_comoving * math.cos(dec)

			x.append(xv)
			y.append(yv)
			z.append(zv)

	return x, y, z

#r = np.logspace(-.7,1.6,25)
r = np.array([200.])
print r,r.shape

#data_file = "data/2SLAQ_clustering/testset_data_4.dat"
#rand_file = "data/2SLAQ_clustering/testset_rand_4.dat"

data_file = "data/2SLAQ_clustering/3yr_obj_v4_Sample8_nota01ao2s01_radecz.cat"
rand_file = "data/2SLAQ_clustering/rand_LRG_Sm8_spe.dat"

data_x, data_y, data_z = readfile(data_file)
rand_x, rand_y, rand_z = readfile(rand_file)
X_data = np.array([data_x,data_y,data_z])
X_random = np.array([rand_x,rand_y,rand_z])

# --------- C library ---------

#initialize object
cons = twopoint.twopoint_constructor(X_data,X_random)

#construct trees
start = time.clock()
cons.construct()
end = time.clock()
test_buildtime = end - start

#calculate correlation function
start = time.clock()
ls_test = cons.landy_szalay(r)
end = time.clock()
test_runtime = end - start

# --------- SKLearn library ---------

#construct trees
start = time.clock()
dtree = neighbors.KDTree(X_data.T)
rtree = neighbors.KDTree(X_random.T, leaf_size = 40)
end = time.clock()
sk_buildtime = end - start

#calclate correlation function
start = time.clock()
dd = dtree.two_point_correlation(X_data.T,r)
dr = rtree.two_point_correlation(X_data.T,r)
rr = rtree.two_point_correlation(X_random.T,r)

f = float(X_random.shape[1])/X_data.shape[1]

ls_scikit = (f*f*dd-2*f*dr+rr)/[float(i) for i in rr]
end = time.clock()
sk_runtime = end - start

'''
plt.loglog(r,ls_test,r,ls_scikit)
plt.xlabel("$h^-1 Mpc$")
plt.ylabel("$\\psi(s)$")
plt.figsave("compare_function.png", format='png')
plt.clf()
plt.xticks([1,2],('C library','scikit-learn'))
plt.ylabel("Calculation time (sec)")
plt.bar([1,2],[test_runtime,sk_runtime], align='center')
plt.figsave("compare_runtime.png", format='png')
'''

print "C Library took %f sec to build, %f sec to run" % (test_buildtime, test_runtime)
print "scikit-learn took %f sec to build, %f sec to run" % (sk_buildtime, sk_runtime)
print "C Library output:"
print ls_test
print "scikit-learn output:"
print ls_scikit
