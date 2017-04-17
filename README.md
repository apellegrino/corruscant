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

dsize = 4000
rsize = dsize*20

np.random.seed(3)
X_data = np.random.rand(dsize,3)
X_random = np.random.rand(rsize,3)

radii = np.logspace(-2.5,-1.,6)

results = tpcf.pair_counts(X_data, X_random, radii,
                           xi_estimator_type="landy-szalay",
                           xi_error_type="field-to-field", num_threads=2)

print "------- Using Landy-Szalay estimator, field-to-field error -------"

print "r = "
print radii

print "DD, DR, RR differential counts"
dd = results["DD"]
dr = results["DR"]
rr = results["RR"]
print np.vstack([dd,dr,rr]).T

print "Estimated Xi(r)"
print results["estimator"]
print "sigma^2 of Xi(r)"
print results["error"]
```
