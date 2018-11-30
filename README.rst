=================
OpenET - DisALEXI
=================

|version| |build|

This repository provides an Earth Engine Python API based implementation of the DisALEXI ET model.

Input Collections
=================

DisALEXI ET can currently only be computed for Landsat Collection 1 TOA image from the following Earth Engine image collections:

 * LANDSAT/LC08/C01/T1_RT_TOA or LANDSAT/LC08/C01/T1_TOA or LANDSAT/LC08/C01/T1_SR
 * LANDSAT/LE07/C01/T1_RT_TOA or LANDSAT/LE07/C01/T1_TOA or LANDSAT/LE07/C01/T1_SR
 * LANDSAT/LT05/C01/T1_TOA

Model Design
============

The primary component of the SSEBop model is the Image() class.  The Image class should generally be instantiated from an Earth Engine Landsat image using the collection specific methods listed below.  ET image collections can be built by computing ET in a function that is mapped over a collection of input images.  Please see the `Example Notebooks`_ for more details.

Landsat Collection 1 TOA Input Image
------------------------------------

To instantiate the class for a Landsat Collection 1 TOA image, use the Image().from_landsat_c1_toa() method.

The input Landsat image must have the following bands and properties:

=================  ======================================
SPACECRAFT_ID      Band Names
=================  ======================================
LANDSAT_5          B1, B2, B3, B4, B5, B7, B6, BQA
LANDSAT_7          B1, B2, B3, B4, B5, B7, B6_VCID_1, BQA
LANDSAT_8          B2, B3, B4, B5, B6, B7, B10, BQA
=================  ======================================

=================  =============================================
Property           Description
=================  =============================================
system:index       - Landsat Scene ID
                   - Must be in the Earth Engine format (e.g. LC08_044033_20170716)
                   - Used to lookup the scene specific c-factor
system:time_start  Image datetime in milliseconds since 1970
SPACECRAFT_ID      - Used to determine which Landsat type
                   - Must be: LANDSAT_5, LANDSAT_7, or LANDSAT_8
=================  =============================================

Model Output
------------

Example Notebooks
=================

Jupyter notebooks are provided in the "examples" folder that show various approaches for calling the OpenET DisALEXI model.

+ `Computing daily ET for a single Landsat image <examples/single_image.ipynb>`__
+ `Computing a daily ET image collection from Landsat image collection <examples/collection.ipynb>`__
+ `Computing annual ET from a collection <examples/interpolate.ipynb>`__

Ancillary Datasets
==================

ALEXI ET
------------------------------------
The daily maximum air temperature (Tmax) is essential for establishing the maximum ET limit (cold boundary) as explained in Anderson2007JGR_.

Default Asset ID: xxxx

Land Surface Temperature
------------------------
Land Surface Temperature (LST) is currently calculated in the SSEBop approach from Landsat Top-of-Atmosphere images by including commonly used calibration steps and atmospheric correction techniques. These include calculations for: (1) spectral radiance conversion to the at-sensor brightness temperature; (2) atmospheric absorption and re-emission value; (3) surface emissivity; and (4) land surface temperature. For additional information, users can refer to section 3.2 of the Methodology in Senay2016_.

dT
--
The SSEBop ET model uses dT as a predefined temperature difference between Thot and Tcold for each pixel.
In SSEBop formulation, hot and cold limits are defined on the same pixel; therefore, dT actually represents the vertical temperature difference between the surface temperature of a theoretical bare/dry condition of a given pixel and the air temperature at the canopy level of the same pixel as explained in Senay2013_. The input dT is calculated under average-sky conditions and assumed not to change from year to year, but is unique for each day and location.

Default Asset ID: projects/usgs-ssebop/dt/daymet_median_v1

Elevation
---------
The default elevation dataset is the USGS SRTM global image asset.

Default Asset ID: `USGS/SRTMGL1_003 <https://developers.google.com/earth-engine/datasets/catalog/USGS_SRTMGL1_003>`__

The elevation parameter will accept any Earth Engine image.

Installation
============

The OpenET DisALEXI python module can be installed via pip:

.. code-block:: console

    pip install openet-disalexi

Dependencies
============

Modules needed to run the model:

 * `earthengine-api <https://github.com/google/earthengine-api>`__
 * `openet <https://github.com/Open-ET/openet-core-beta>`__

OpenET Namespace Package
========================

Each OpenET model should be stored in the "openet" folder (namespace).  The benefit of the namespace package is that each ET model can be tracked in separate repositories but called as a "dot" submodule of the main openet module,

.. code-block:: console

    import openet.disalexi as disalexi


References
==========

.. _references:

.. [Anderson2007JGR] Anderson, M. C., J. M. Norman, J. R. Mecikalski, J. A. Otkin, and W. P. Kustas (2007), A climatological study of evapotranspiration and moisture stress across the continental United States based on thermal remote sensing: 1. Model formulation, J. Geophys. Res., 112, D10117, doi:10.1029/2006JD007506.

.. |build| image:: https://travis-ci.org/Open-ET/openet-disalexi-beta.svg?branch=master
   :alt: Build status
   :target: https://travis-ci.org/Open-ET/openet-disalexi-beta
.. |version| image:: https://badge.fury.io/py/openet-disalexi.svg
   :alt: Latest version on PyPI
   :target: https://badge.fury.io/py/openet-disalexi