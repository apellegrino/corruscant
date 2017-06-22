import tpcf_angular
import numpy as np

nd = 10000
nr = 20 * nd

data = np.random.rand(2,nd)
data[0] = np.arcsin(data[0] * 2. - 1.)
data[1] = data[1] * 2 * np.pi

rand = np.random.rand(2,nr)
rand[0] = np.arcsin(rand[0] * 2. - 1.)
rand[1] = rand[1] * 2 * np.pi

dtree = tpcf_angular.tree(data)
rtree = tpcf_angular.tree(rand)
results = tpcf_angular.twopoint(dtree, rtree, radii=[0,0.05], err_type=None)
