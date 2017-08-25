import numpy as np

from corruscant import twopoint
from corruscant import clustering

x = range(2,201)

jk, bs = [], []
for N_fields in x:
    # data set sizes
    dsize = 20000
    rsize = dsize*20

    # generate random points in a cube
    np.random.seed(3)
    X_data = np.random.rand(dsize, 3)
    X_rand = np.random.rand(rsize, 3)

    # choose radius bin edges
    radii = np.logspace(-2.5,-1.,6)

    # give each data point and random point a field ID
    # we will split up the data into fourths in the x dimension

    data_fields = (X_data[:,0] * N_fields).astype('int32')
    rand_fields = (X_rand[:,0] * N_fields).astype('int32')

    # generate K-d trees
    dtree = clustering.tree(X_data, data_fields)
    rtree = clustering.tree(X_rand, rand_fields)

    results = twopoint.threedim.autocorr(dtree, rtree, radii, num_threads=26)
    jk.append(results.error())
    bs.append(results.bootstrap_error())
    print results
    #times=timeit.repeat(stmt, setup.format(N_fields=N_fields), repeat=5, number=1)
    #print times
    #tot.append(times)

jk = np.array(jk)
bs = np.array(bs)

import matplotlib.pyplot as plt
import richardsplot

for i,(lo,hi) in enumerate(zip(radii[:-1],radii[1:])):
    
    plt.figure()
    ax = plt.gca()

    jk_slice = jk[:,i]
    bs_slice = bs[:,i]
    plt.plot(x,jk_slice,'o',label='Jackknife Error',ls='none',markersize=2)
    plt.plot(x,bs_slice,'o',label='Bootstrap Error',ls='none',markersize=2)
    plt.xlabel("Number of fields")
    plt.ylabel("Estimated Error")
    plt.title("Error for bin r = [{:.3f}, {:.3f}]".format(lo,hi))
    plt.ylim((0.0,ax.get_ylim()[1]))
    ax.set_yticklabels([format(label, '.04f') for label in ax.get_yticks()])
    plt.tight_layout()
    plt.legend(loc='lower right')
    plt.savefig("linear/err_{:d}.png".format(i+1), dpi=200)
    #plt.show()
exit()

with open("results.csv", 'w') as f:
    f.write("n_fields,avg(sec),stddev(sec)\n")
    for n, times in zip(range(10,210,10), tot):
        f.write("{:d},{:f},{:f}\n".format(n, np.average(times), np.std(times)))
