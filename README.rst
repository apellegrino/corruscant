======
README
======

Corruscant is a C and Python library for calculation of the two-point
correlation function on astronomical data.

Installation instructions for Unix-based systems
------------------------------------------------

.. code-block:: console
    $ cd ~
    $ git clone https://[your Bitbucket username]@bitbucket.org/apellegrino/clustering-tree.git
    $ cd clustering-tree
    $ make python


Then, add the following line to the end of your ~/.bash_profile or ~/.profile file:

.. code-block:: console
    export PYTHONPATH="$PYTHONPATH:/home/[your username here]/clustering-tree"

Then finally:

.. code-block:: console
    source ~/.bash_profile

or

.. code-block:: console
    source ~/.profile

You can now use the library in Python with

.. code-block:: console
    import tpcf

Example usage can be found in test/example.py:

.. code-block:: console
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
                               err_type="jackknife", est_type="landy-szalay",
                               num_threads=4)

    print(results)

Memory Efficiency
-----------------

Corruscant requires that input arrays be C-contiguous and that all values in
each dimension come before any values in the next dimensions. Additionally,
the array of integers describing field IDs must be of Numpy datatype `int32`.
If the user's input does not match these requirements, the input will be copied
to new arrays as necessary. If conserving memory is critical, the user can
update their data to be properly formatted using the `validate_points` and
`validate_fields` functions:

.. code-block:: console
    data = tpcf.validate_points(data)
    fields = tpcf.validate_fields(fields)

This assures that the input arrays will not be copied when constructing a tree.
