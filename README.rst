
=================
OpenET - DisALEXI
=================

|version| |build|

This repository provides Google Earth Engine Python API based implementation of the DisALEXI ET model.

DisALEXI is the disaggregation component of a multi-scale system for modeling actual evapotranspiration (ETa) at field to global scales.  DisALEXI spatially downscales regional gridded ET output from the Atmosphere-Landsat Exchange Inverse (ALEXI) model to finer scales using moderate to high resolution remotely sensed land-surface temperature data.  Both ALEXI and DisALEXI are based on the Two Source Energy Balance (TSEB) land-surface representation originally developed by Norman et al., (1995).

Input Collections
=================

DisALEXI ET can currently only be computed for Landsat Collection 1 TOA and Surface Reflection (SR) image from the following Earth Engine image collections:
    * LANDSAT/LC08/C01/T1_RT_TOA or LANDSAT/LC08/C01/T1_TOA or LANDSAT/LC08/C01/T1_SR
    * LANDSAT/LE07/C01/T1_RT_TOA or LANDSAT/LE07/C01/T1_TOA or LANDSAT/LE07/C01/T1_SR
    * LANDSAT/LT05/C01/T1_TOA

Note that this version of DisALEXI can only be run over CONUS.

Model Design
============

The primary component of the DisALEXI model is the Image() class. The Image class should generally be instantiated from an Earth Engine Landsat image using the collection specific methods listed below. ET image collections can be built by computing ET in a function that is mapped over a collection of input images. Please see the `Example Notebooks` for more details.

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

system:time_start  Image datetime in milliseconds since 1970
SPACECRAFT_ID      - Used to determine which Landsat type
                   - Must be: LANDSAT_5, LANDSAT_7, or LANDSAT_8
=================  =============================================

Landsat Collection 1 SR Input Image
-----------------------------------

To instantiate the class for a Landsat Collection 1 TOA image, use the Image().from_landsat_c1_sr() method.

The input Landsat image must have the following bands and properties:

=================  ======================================
SPACECRAFT_ID      Band Names
=================  ======================================
LANDSAT_5          B1, B2, B3, B4, B5, B7, B6, pixel_qa
LANDSAT_7          B1, B2, B3, B4, B5, B7, B6_VCID_1, pixel_qa
LANDSAT_8          B2, B3, B4, B5, B6, B7, B10, pixel_qa
=================  ======================================

=================  =============================================
Property           Description
=================  =============================================
system:index       - Landsat Scene ID
                   - Must be in the Earth Engine format (e.g. LC08_044033_20170716)
system:time_start  Image datetime in milliseconds since 1970
SPACECRAFT_ID      - Used to determine which Landsat type
                   - Must be: LANDSAT_5, LANDSAT_7, or LANDSAT_8
=================  =============================================

Model Output
------------

The primary output of the DisALEXI model is daily ETa.  Internally this is partitioned to contributions from soil evaporation (E) and canopy transpiration (T).

Examples
--------

Jupyter notebooks are provided in the "examples" folder that show various approaches for calling the OpenET DisALEXI model.

Example Notebooks
=================

Jupyter notebooks are provided in the "examples" folder that show various approaches for calling the OpenET DisALEXI model.


* `computing daily ET for a single Landsat TOA image`
* `computing daily ET for a single Landsat SR image`

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
.. [Anderson2012a] Anderson, M. C., R. G. Allen, A. Morse, W. P. Kustas (2012a), Use of Landsat thermal imagery in monitoring evapotranspiration and managing water resources, Remote Sens. Environ. 122, 50-65.
.. [Anderson2018] Anderson, M. C., F. Gao, K. Knipper, C. Hain, W. Dulaney, D. D. Baldocchi, E. Eichelmann, K. S. Hemes, Y. Yang, J. Medellin-Azuara, W. P. Kustas (2018), Field-scale assessment of land and water use change over the California Delta using remote sensing. Remote Sens. 10:889.
.. [Norman1995] Norman, J. M., W. P. Kustas, K. S. Humes (1995), A two-source approach for estimating soil and vegetation energy fluxes from observations of directional radiometric surface temperature. Agric. For. Meteorol. 77:263-293.
.. [Anderson2007] Anderson, M. C., J. M. Norman, J. R. Mecikalski, J. A. Otkin, and W. P. Kustas (2007), A climatological study of evapotranspiration and moisture stress across the continental United States based on thermal remote sensing: 1. Model formulation, J. Geophys. Res., 112, D10117, doi:10.1029/2006JD007506.
.. [Anderson1997] Anderson, M. C., J. M. Norman, G. R. Diak, W. P. Kustas, J. R. Mecikalski (1997), A two-source time integrated model for estimating surface fluxes using thermal infrared remote sensing, Remote Sens. Environ. 60, 195-216.
.. [Anderson2004] Anderson, M. C., J. M. Norman, J. R. Mecikalski, R. D. Torn, W. P. Kustas, J. B. Basara (2004), A multiscale remote sensing model for disaggregating regional fluxes to micrometeorological scales, J. Hydrometeorol. 5, 343-363.
.. [Anderson2012b] Anderson, M. C., W.P. Kustas, J. G. Alfieri, F. Gao, C. Hain, J. H. Prueger, S. Evett, P. Colaizzi, T. Howell, J. L. Chavez (2012b), Mapping daily evapotranspiration at Landsat spatial scales during the BEAREX'08 field campaign (2012b), Adv. Water Resour, 50, 162-177.

.. |build| image:: https://travis-ci.org/Open-ET/openet-disalexi-beta.svg?branch=master
   :alt: Build status
   :target: https://travis-ci.org/Open-ET/openet-disalexi-beta
.. |version| image:: https://badge.fury.io/py/openet-disalexi.svg
   :alt: Latest version on PyPI
