=================
OpenET - DisALEXI
=================

|version| |build|

This repository provides Google Earth Engine Python API based implementation of the DisALEXI ET model.

DisALEXI is the disaggregation component of a multi-scale system for modeling actual evapotranspiration (ETa) at field to global scales.  DisALEXI spatially downscales regional gridded ET output from the Atmosphere-Landsat Exchange Inverse (ALEXI) model to finer scales using moderate to high resolution remotely sensed land-surface temperature data.  Both ALEXI and DisALEXI are based on the Two Source Energy Balance (TSEB) land-surface representation originally developed by Norman et al., (1995).

Input Collections
=================

DisALEXI ET can currently be computed for Landsat Collection 2 Level 2 (surface reflectance) images from the following Earth Engine image collections:

 * LANDSAT/LC09/C02/T1_L2
 * LANDSAT/LC08/C02/T1_L2
 * LANDSAT/LE07/C02/T1_L2
 * LANDSAT/LT05/C02/T1_L2

**Note that this version of DisALEXI can only be run over the conterminous United States (CONUS)**

Model Design
============

The primary component of the DisALEXI model is the Image() class. The Image class should generally be instantiated from an Earth Engine Landsat image using the collection specific methods listed below. ET image collections can be built by computing ET in a function that is mapped over a collection of input images. Please see the `Example Notebooks` for more details.

Landsat Collection 2 SR Input Image
-----------------------------------

To instantiate the class for a Landsat Collection 2 level 2 image, use the Image.from_landsat_c2_sr() method.

The input Landsat image must have the following bands and properties:

=================  ===========================================
SPACECRAFT_ID      Band Names
=================  ===========================================
LANDSAT_5          SR_B1, SR_B2, SR_B3, SR_B4, SR_B5, SR_B7, ST_B6, QA_PIXEL
LANDSAT_7          SR_B1, SR_B2, SR_B3, SR_B4, SR_B5, SR_B7, ST_B6, QA_PIXEL
LANDSAT_8          SR_B1, SR_B2, SR_B3, SR_B4, SR_B5, SR_B6, SR_B7, ST_B10, QA_PIXEL
LANDSAT_9          SR_B1, SR_B2, SR_B3, SR_B4, SR_B5, SR_B6, SR_B7, ST_B10, QA_PIXEL
=================  ===========================================

=================  =============================================
Property           Description
=================  =============================================
system:id          Earth Engine Asset ID (e.g. LANDSAT/LC08/C02/T1_L2/LC08_044033_20170716)
system:index       - Landsat Scene ID
                   - Must be in the Earth Engine format (e.g. LC08_044033_20170716)
system:time_start  Image datetime in milliseconds since 1970
SPACECRAFT_ID      - Used to determine which Landsat type
                   - Must be: LANDSAT_5, LANDSAT_7, LANDSAT_8, LANDSAT_9
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
 * `openet <https://github.com/Open-ET/openet-core>`__

OpenET Namespace Package
========================

Each OpenET model should be stored in the "openet" folder (namespace).  The benefit of the namespace package is that each ET model can be tracked in separate repositories but called as a "dot" submodule of the main openet module.

.. code-block:: console

    import openet.disalexi as disalexi

References
==========

.. _references:
.. [Anderson2012a] Anderson, M. C., R. G. Allen, A. Morse, W. P. Kustas (2012a), Use of Landsat thermal imagery in monitoring evapotranspiration and managing water resources, Remote Sens. Environ. 122, 50-65.    `https://doi.org/10.1016/j.rse.2011.08.025 <https://doi.org/10.1016/j.rse.2011.08.025>`__

.. [Anderson2018a] Anderson, M. C., F. Gao, K. Knipper, C. Hain, W. Dulaney, D. D. Baldocchi, E. Eichelmann, K. S. Hemes, Y. Yang, J. Medellin-Azuara, W. P. Kustas (2018), Field-scale assessment of land and water use change over the California Delta using remote sensing. Remote Sens. 10:889. `https://doi.org/10.3390/rs10060889 <https://doi.org/10.3390/rs10060889>`__

.. [Norman1995] Norman, J. M., W. P. Kustas, K. S. Humes (1995), A two-source approach for estimating soil and vegetation energy fluxes from observations of directional radiometric surface temperature. Agric. For. Meteorol. 77:263-293.  `https://doi.org/10.1016/0168-1923(95)02265-Y <https://doi.org/10.1016/0168-1923(95)02265-Y>`__

.. [Anderson2007] Anderson, M. C., J. M. Norman, J. R. Mecikalski, J. A. Otkin, and W. P. Kustas (2007), A climatological study of evapotranspiration and moisture stress across the continental United States based on thermal remote sensing: 1. Model formulation, J. Geophys. Res., 112, D10117. `https://doi.org/10.1029/2006JD007506 <https://doi.org/10.1029/2006JD007506>`__

.. [Anderson1997] Anderson, M. C., J. M. Norman, G. R. Diak, W. P. Kustas, J. R. Mecikalski (1997), A two-source time integrated model for estimating surface fluxes using thermal infrared remote sensing, Remote Sens. Environ. 60, 195-216. `https://doi.org/10.1029/2006JD007507 <https://doi.org/10.1029/2006JD007507>`__

.. [Anderson2004] Anderson, M. C., J. M. Norman, J. R. Mecikalski, R. D. Torn, W. P. Kustas, J. B. Basara (2004), A multiscale remote sensing model for disaggregating regional fluxes to micrometeorological scales, J. Hydrometeorol. 5, 343-363. `https://doi.org/10.1175/1525-7541(2004)005<0343:AMRSMF>2.0.CO;2 <https://doi.org/10.1175/1525-7541(2004)005%3C0343:AMRSMF%3E2.0.CO;2>`__

.. [Anderson2012b] Anderson, M. C., W.P. Kustas, J. G. Alfieri, F. Gao, C. Hain, J. H. Prueger, S. Evett, P. Colaizzi, T. Howell, J. L. Chavez (2012b), Mapping daily evapotranspiration at Landsat spatial scales during the BEAREX'08 field campaign (2012b), Adv. Water Resour, 50, 162-177. `https://doi.org/10.1016/j.advwatres.2012.06.005 <https://doi.org/10.1016/j.advwatres.2012.06.005>`__

.. [Cammalleri2014] Cammalleri, C., M.C. Anderson, F. Gao, C.R. Hain, W.P. Kustas, Mapping daily evapotranspiration at field scales over rainfed and irrigated agricultural areas using remote sensing data fusion, Agricultural and Forest Meteorology, Volume 186, 2014, Pages 1-11, ISSN 0168-1923, https://doi.org/10.1016/j.agrformet.2013.11.001.

.. [Semmens2016] Semmens, K. A., M. C. Anderson, W. P. Kustas, F. Gao, J. G. Alfieri, L. McKee, J. H. Prueger, C. R. Hain, C. Cammalleri, Y. Yang and T. Xia (2016), Monitoring daily evapotranspiration over two California vineyards using Landsat 8 in a multi-sensor data fusion approach, Remote Sens. Environ., 185, 155–170.

.. [Yang2017] Yang, Y., M. C. Anderson, F. Gao, C. R. Hain, K. A. Semmens, W. P. Kustas, A. Noormets, R. H. Wynne, V. A. Thomas, and G. Sun (2017), Daily Landsat-scale evapotranspiration estimation over a forested landscape in North Carolina, USA using multi-satellite data fusion, Hydrol. Earth Syst. Sci., 21, 1017–1037.  `doi:doi:10.5194/hess-21-1017-2017 <doi:doi:10.5194/hess-21-1017-2017>`__

.. [Anderson2018b] Anderson, M., G. Diak, F. Gao, K. Knipper, C. Hain, E. Eichelmann, K.S. Hemes, D. Baldocchi, W. Kustas, Y. Yang (2018), Impact of Insolation Data Source on Remote Sensing Retrievals of Evapotranspiration over the California Delta, Remote Sensing, 11(3): 216, doi: 10.3390/rs11030216

.. [Yang2020] Yang, Y., M. Anderson, F. Gao, C. Hain, A. Noormets, G. Sun, R. Wynne and V. Thomas (2020), Investigating impacts of drought and disturbance on evapotranspiration over a forested landscape in North Carolina, USA using high spatiotemporal resolution remotely sensed data, Remote Sensing of Environment, 238, p. 111018

.. |build| image:: https://github.com/Open-ET/openet-disalexi/actions/workflows/build.yml/badge.svg
   :alt: Build status
   :target: https://travis-ci.org/Open-ET/openet-disalexi
.. |version| image:: https://badge.fury.io/py/openet-disalexi.svg
   :alt: Latest version on PyPI
   :target: https://badge.fury.io/py/openet-disalexi
