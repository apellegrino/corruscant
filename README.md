# README #

This is a C and Python library for calculation of the two-point correlation function on 3D cartesian data.

### Installation instructions ###

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
source ~/.bash_profile *(or ~/.profile)*
```
You can now use the library in Python with
```
#!python
import tpcf
```