======
README
======

Corruscant is a C and Python library for calculation of the two-point
correlation function on astronomical data. It currently includes the
following features:

* Auto- or cross- correlation
* 3-dimensional (\xi(r)) or angular (w(\theta)) functions
* Jackknife, bootstrap, field-to-field, and Poisson errors
* Automatic division of (\theta, \phi) or (RA, DEC) data into fields

Installation
------------

Corruscant is currently packaged with pip, although it is not yet on the
PyPI public repository. This means installation is somewhat more complicated
for the moment. To install, distutils must be used to create a wheel file,
which can then be installed with pip like a regular package.

With Python 2:

.. code-block:: console

    python setup.py bdist_wheel
    pip install dist/*cp2*.whl

With Python 3:

.. code-block:: console

    python3 setup.py bdist_wheel
    pip3 install dist/*cp3*.whl

When updating, uninstall the current version before reinstalling:

With Python 2:

.. code-block:: console

    sudo -H pip uninstall corruscant <<< y

With Python 3:

.. code-block:: console

    sudo -H pip3 uninstall corruscant <<< y


Usage
-----
You can use the library in Python with

.. code-block:: python

    import tpcf

Example usage can be found in test/example.py. Because of how Python looks
for packages, be sure to change into the test/ directory before trying to
run it. `python test/example.py` will not work from the root directory.

example.py:

.. code-block:: python

    from corruscant import twopoint
    from corruscant import clustering

    import numpy as np

    # data set sizes
    dsize = 4000
    rsize = dsize*20

    # generate random points in a cube
    np.random.seed(3)
    X_data = np.random.rand(dsize, 3)
    X_rand = np.random.rand(rsize, 3)

    # choose radius bin edges
    radii = np.logspace(-2.5,-1.,6)

    # give each data point and random point a field ID
    # we will split up the data into fourths in the x dimension
    data_fields = (X_data[:,0] * 4).astype('int32')
    rand_fields = (X_rand[:,0] * 4).astype('int32')

    # generate K-d trees
    dtree = clustering.tree(X_data, data_fields)
    rtree = clustering.tree(X_rand, rand_fields)

    # get the correlation function results
    results = twopoint.threedim.autocorr(dtree, rtree, radii,
                                        est_type="landy-szalay", num_threads=4)

    results.error_type = "jackknife"
    print(results)

Memory Efficiency
-----------------

Corruscant requires that input arrays be C-contiguous and that all values in
each dimension come before any values in the next dimensions. Additionally,
the array of integers describing field IDs must be of Numpy datatype ``int32``.
If the user's input does not match these requirements, the input will be copied
to new arrays as necessary. If conserving memory is critical, the user can
update their data to be properly formatted using the ``validate_points()`` and
``validate_fields()`` functions:

.. code-block:: python

    data = tpcf.validate_points(data)
    fields = tpcf.validate_fields(fields)

This assures that the input arrays will not be copied when constructing a tree.
