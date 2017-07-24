# README #

This is a C and Python library for calculation of the two-point correlation function on 3D cartesian data.

### Installation instructions for Unix-based systems ###

```
#!bash

cd ~
git clone https://[your Bitbucket username]@bitbucket.org/apellegrino/clustering-tree.git
cd clustering-tree
make python
```
Then, add the following line to the end of your ~/.bash_profile or ~/.profile file:
```
#!bash
export PYTHONPATH="$PYTHONPATH:/home/[your username here]/clustering-tree"
```
Then finally:
```
#!bash
source ~/.bash_profile
```
or
```
#!bash
source ~/.profile
```
You can now use the library in Python with
```
#!python
import tpcf
```

Example usage can be found in test/example.py:

```
#!python
import tpcf
import numpy as np

# data set sizes
dsize = 4000
rsize = dsize*20

# generate random points in a cube
np.random.seed(3)
X_data = np.random.rand(dsize, 3)
X_random = np.random.rand(rsize, 3)

# choose radius bin edges
radii = np.logspace(-2.5,-1.,6)

# give each data point and random point a field ID
# for a simple example, call the x < 0.5 region `field 1` and the x > 0.5
# region `field 2`
data_fields = np.where(X_data[:,0] < 0.5, 1, 2)
rand_fields = np.where(X_random[:,0] < 0.5, 1, 2)

# generate K-d trees
dtree = tpcf.tree(X_data, data_fields)
rtree = tpcf.tree(X_random, rand_fields)

# get the correlation function results
results = tpcf.twopoint(dtree, rtree, radii,
                           est_type="landy-szalay",
                           err_type='jackknife', num_threads=4)

print(results)
```

### Memory Efficiency ###
`clustering-tree` requires that input arrays be C-contiguous and that all values in each dimension come before any values in the next dimensions. Additionally, the array of integers describing field IDs must be of Numpy datatype `int32`. If the user's input does not match these requirements, the input will be copied to new arrays as necessary. If conserving memory is critical, the user can update their data to be properly formatted using the `validate_points` and `validate_fields` functions:

```
#!python
data = tpcf.validate_points(data)
fields = tpcf.validate_fields(fields)
```

This assures that the input arrays will not be copied when constructing a tree.
