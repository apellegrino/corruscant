from corruscant import clustering
import numpy as np

def ndupes(arr):
    arrnew = np.sort(arr, kind='mergesort')

    N = len(arr)
    c = 0
    for i in range(N-1):
        j = i+1
        if np.all(arrnew[i] == arrnew[j]):
            c += 1
    return c

N = 30000
X_data = np.repeat(np.random.rand(N, 3).astype('float64'), 3, axis=0)

print("Building tree with points of shape {:s}".format(str(X_data.shape)))
dfields = (X_data[:,0] * 4).astype('int32')
t = clustering.tree(X_data, dfields)
print("Tree successfully built with {:d} unique points".format(t.size))
