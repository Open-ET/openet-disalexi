===============================
Contributing to OpenET DisALEXI
===============================

Thank you for your interesting in support the OpenET DisALEXI project.

Versioning
==========

The OpenET DisALEXI project is currently in Beta and the version numbers should be "0.0.X" until a non-Beta release is made.

Coding Conventions
==================

OpenET DisALEXI primarily supports Python 3.6 at this time.

The code will likely work on other version of Python 3 but there are no plans to support Python 2.7.

All code should follow the `PEP8
<https://www.python.org/dev/peps/pep-0008/>`__ style guide.

Docstrings should be written for all functions that follow the `NumPy docstring format <https://numpydoc.readthedocs.io/en/latest/format.html>`__.

Testing
=======

PyTest
------

Testing is done using `pytest <https://docs.pytest.org/en/latest/>`__.

.. code-block:: console

    python -m pytest

Detailed testing results can be obtained using the "-v" and/or "-s" tags.

.. code-block:: console

    python -m pytest -v -s

Approach
--------

The current test suite is comparing precomputed values at three sites for 1 day (the AmeriFlux sites in Nebraska for the Landsat image LC08_028031_20140708).  The approach for testing is to generate separate input images for each site that are ee.Image.constants() and then extract the output value at the sites using a simple reduceRegion call (see the constant_image_value function in disalexi/tests/test_utils.py constant_image_value).  These calculations are relatively fast but they do not provide full coverage of the possible inputs and outputs.
