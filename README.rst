.. image:: https://raw.githubusercontent.com/solvebio/veppy/master/images/veppy.png
    :target: https://github.com/solvebio/veppy
    :alt: Veppy: Variant Effect Prediction for Python

.. image:: https://travis-ci.org/solvebio/veppy.svg?branch=master
       :target: https://travis-ci.org/solvebio/veppy

About Veppy
-----------

Veppy is a genetic variant effect predictor for Python. Inspired by SnpEff and VEP.

**WARNING:** This code is an alpha release and not production-ready. APIs may change at any time.


.. image:: https://img.shields.io/pypi/v/veppy.svg
    :target: https://pypi.python.org/pypi/veppy


Installation
------------

.. code-block:: bash

    $ pip install veppy


Installation from source
------------------------

.. code-block:: bash

    $ git clone git@github.com:solvebio/veppy.git
    $ cd veppy
    $ python setup.py install


Setup
-----

**Step 1 (OPTIONAL):** Prepare a directory for veppy data

The default data path is: ``./data``

You can override this by setting ``$VEPPY_DATA_DIR``.

.. code-block:: bash

    export VEPPY_DATA_DIR=/opt/veppy


**Step 2:** Download source data and build indexes

**NOTE:** This step downloads about 1gb of data. After indexing, the data directory will consume about 8gb of disk space.

.. code-block:: bash

    ./scripts/download_data_GRCh37.sh


**Step 3:** Index the source data

.. code-block:: bash

    python ./run_index.py


Example Usage
-------------

.. code-block:: python

    >>> from veppy.veppy import calculate_consequences
    >>> variant = ('1', 8025384, 'A', 'T')
    >>> result = calculate_consequences('GRCh37', *variant)
    >>> print result.results


Testing
-------

Tests are currently based on chr1 versions of input data.
Full genome tests are coming soon!

.. code-block:: bash

    # pip install -r requirements.txt
    $ nosetests


Coverage:


.. code-block:: bash

    $ nosetests --with-coverage --cover-package=veppy


About SolveBio
--------------

SolveBio is a genomics company based in New York City.

.. image:: https://raw.githubusercontent.com/solvebio/veppy/master/images/solvebio.png
    :target: https://www.solvebio.com/
    :alt: SolveBio
