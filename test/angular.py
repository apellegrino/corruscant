import tpcf
import tpcf_angular
import tpcf_helper as thelp
import numpy as np
import timer
t = timer.timer()

nd = 1000
nr = 20 * nd

def sph_to_cart(lat, lon):
    x = np.cos(lat)*np.cos(lon)
    y = np.cos(lat)*np.sin(lon)
    z = np.sin(lat)

    return x,y,z

data = np.random.rand(2,nd)
data[0] = np.arcsin(data[0] * 2. - 1.)
data[1] = data[1] * 2 * np.pi

rand = np.random.rand(2,nr)
rand[0] = np.arcsin(rand[0] * 2. - 1.)
rand[1] = rand[1] * 2 * np.pi

dfields = np.where(data[0] < 0., 1, 2)
rfields = np.where(rand[0] < 0., 1, 2)

dtree = tpcf_angular.tree(data, dfields)
t.start()
rtree = tpcf_angular.tree(rand, rfields)
print "VP build time:"
t.printout()

#radii = np.logspace(-3,-1.5,5)
radii=[0,1.0]
eu = [2*np.sin(r/2) for r in radii]
t.start()
results = tpcf_angular.twopoint(dtree, rtree, radii=radii, num_threads=8)
print "VP query time:"
t.printout()
print "------- Using Landy-Szalay estimator, jackknife error -------"

print "theta = "
print radii

print "DD, DR, RR differential counts"
dd, dr, rr = results.total_pair_counts()
print np.vstack([dd,dr,rr]).T

datacart = np.array(sph_to_cart(data[0], data[1]))
randcart = np.array(sph_to_cart(rand[0], rand[1]))

dtree = tpcf.tree(datacart, dfields)
t.start()
rtree = tpcf.tree(randcart, rfields)
print "KD build time:"
t.printout()

t.start()
results = tpcf.twopoint_angular(dtree, rtree, radii=radii, num_threads=8)
print "KD query time:"
t.printout()
dd, dr, rr = results.total_pair_counts()
print np.vstack([dd,dr,rr]).T

#print "Estimated w(theta)"
#print results.estimate()
#print "Error of w(theta)"
#print results.error()
#print "Covariance matrix"
#print results.covariance()
