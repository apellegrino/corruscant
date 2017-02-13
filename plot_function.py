import math
from scipy import integrate
import numpy as np
from sklearn import neighbors
import matplotlib
#matplotlib.use('Agg')
from matplotlib import pyplot as plt
from astroML import correlation

# constants
c = 299792.4580 # km/s
H_0 = 100.0 #h^-1 km/s/Mpc


def integrand(a, omega_lambda, omega_k, omega_matter, omega_rel):
	return 1.0/math.sqrt(omega_rel + omega_matter*a + omega_k*a*a + 
										omega_lambda*a*a*a*a)

# return floats from string of float literals separated by spaces
def parse_line(string):
	return [float(i) for i in filter(None, string.strip('\r\n').split(' '))]

def get_comoving(redshift):
    integral, error = integrate.quad(integrand, # function
                                    1/(1.+redshift), 1.0, # limits
                                    args = (0.7, 0., 0.3, 0.))
    return c / H_0 * integral

def readfile_for_astroML(path):
    x = []
    y = []
    z = []
    with open(path, "r") as f:
        for line in f:
            ra, dec, redshift = parse_line(line)

            integral, error = integrate.quad(integrand,1/(1+redshift),1.0,
												args = (0.7, 0., 0.3, 0.))
            d_comoving = c / H_0 * integral # units of Mpc
            xv, yv, zv = tuple([d_comoving*i for i in  correlation.ra_dec_to_xyz(ra,dec)])
            x.append(xv)
            y.append(yv)
            z.append(zv)

    return x, y, z

def readfile(path):
	x = []
	y = []
	z = []
	with open(path, "r") as f:
		for line in f:
			ra, dec, redshift = parse_line(line)

			# convert from degrees to radians
			ra *= math.pi/180.0
			dec *= math.pi/180.0

			integral, error = integrate.quad(integrand,1/(1+redshift),1.0,
												args = (0.7, 0., 0.3, 0.))
			d_comoving = c / H_0 * integral # units of Mpc

			# transform to cartesian coordinates
			xv = d_comoving * math.sin(math.pi/2 - dec) * math.sin(ra)
			yv = d_comoving * math.sin(math.pi/2 - dec) * math.cos(ra)
			zv = d_comoving * math.cos(math.pi/2 - dec)

			x.append(xv)
			y.append(yv)
			z.append(zv)

	return x, y, z


#data_file = "data/2SLAQ_clustering/testset_data_4.dat"
#rand_file = "data/2SLAQ_clustering/testset_rand_4.dat"

data_file = "data/2SLAQ_clustering/3yr_obj_v4_Sample8_nota01ao2s01_radecz.cat"
rand_file = "data/2SLAQ_clustering/rand_LRG_Sm8_spe.dat"
old_calc_file = "data/2SLAQ_clustering/k_output_jack_perl_full_newcor_v2_ed.dat"

# Fields:
# [0] log10(s (h^-1 Mpc))
# [1] s (h^-1 Mpc)
# [2]
# [3]
# [4] DD(s)
# [5] DR(s)
# [6] Xi(s)
# [7]
# [8] RR(s)
# [9]
# [10]
Slaq=open('data/2SLAQ_clustering/k_output_jack_perl_full_newcor_v2_ed.dat','r')
slaqdat=Slaq.readlines()
logR = []
R=[]
xislaq=[]
r2xislaq=[]
DR=[]
DD=[]
RR=[]
for i in range(1,len(slaqdat)-1):
	#print slaqdat[i].split()[1], slaqdat[i].split()[6]
	logR.append(float(slaqdat[i].split()[0]))
	R.append(float(slaqdat[i].split()[1]))
	xislaq.append(float(slaqdat[i].split()[6]))
	r2xislaq.append(float(slaqdat[i].split()[1])**2*float(slaqdat[i].split()[6]))
	DD.append(float(slaqdat[i].split()[4]))
	DR.append(float(slaqdat[i].split()[5]))
	RR.append(float(slaqdat[i].split()[8]))

r = 10**np.arange(-.75,2.35,.1)

#astroML
data_x_a, data_y_a, data_z_a = readfile_for_astroML(data_file)
rand_x_a, rand_y_a, rand_z_a = readfile_for_astroML(rand_file)
X_data_a = np.array([data_x_a,data_y_a,data_z_a])
X_random_a = np.array([rand_x_a,rand_y_a,rand_z_a])

ls_astro = correlation.two_point(X_data_a.T, r, method='landy-szalay', data_R=X_random_a.T)

#sklearn
data_x, data_y, data_z = readfile(data_file)
rand_x, rand_y, rand_z = readfile(rand_file)
X_data = np.array([data_x,data_y,data_z])
X_random = np.array([rand_x,rand_y,rand_z])

#ls_astro = correlation.two_point(X_data.T, r, method='landy-szalay', data_R=X_random.T)

#construct trees
dtree = neighbors.KDTree(X_data.T)
rtree = neighbors.KDTree(X_random.T)

# calclate correlation function
# The arrays dd,dr,rr will have length len(r) - 1
dd = np.diff(dtree.two_point_correlation(X_data.T,r))
dr = np.diff(dtree.two_point_correlation(X_random.T,r))
rr = np.diff(rtree.two_point_correlation(X_random.T,r))

f = float(X_random.shape[1])/X_data.shape[1]

ls_estimator = (f*f*dd-2*f*dr+rr)/[float(i) for i in rr]
#ls_estimator = np.delete(ls_estimator,1)

# Match length of dd,dr,rr
r_diff = r[:-1]
#r_diff = r[1:]
#r_diff = np.delete(r_diff,1)

print R
print r_diff
print ls_estimator

#plt.loglog(r,ls_estimator,'o')
plt.loglog(r_diff,ls_estimator, 'o', label='Andrew')
plt.loglog(R,xislaq, label='Old')
plt.loglog(r_diff,ls_astro, '+', label='AstroML')
plt.legend()
plt.xlim([1e-1,1e3])
plt.ylim([1e-6,1e3])
plt.xlabel("$s (h^{-1} Mpc)$")
plt.ylabel("$\\xi(s)$")
plt.savefig("compare_function_diff.png", format='png')
