# from builtins import input
# import datetime
import pprint
import re

import ee
# import numpy as np

# Why can't these be imported directly
from .lc_properties import remaps
from . import landsat
from . import tseb
from . import tseb_utils
from . import utils


def lazy_property(fn):
    """Decorator that makes a property lazy-evaluated

    https://stevenloria.com/lazy-properties/
    """
    attr_name = '_lazy_' + fn.__name__

    @property
    def _lazy_property(self):
        if not hasattr(self, attr_name):
            setattr(self, attr_name, fn(self))
        return getattr(self, attr_name)

    return _lazy_property


class Image(object):
    """Earth Engine DisALEXI"""

    def __init__(
            self,
            image,
            ta_source='CONUS_V006',
            alexi_source='CONUS_V006',
            lai_source='projects/earthengine-legacy/assets/projects/openet/lai/landsat/c02_unsat',
            tir_source='projects/earthengine-legacy/assets/projects/openet/tir/landsat/c02',
            elevation_source='USGS/SRTMGL1_003',
            landcover_source='USGS/NLCD_RELEASES/2019_REL/NLCD',
            ta0_source='CFSR',
            rs_daily_source='CFSR',
            rs_hourly_source='CFSR',
            windspeed_source='CFSR',
            vp_source='CFSR',
            airpressure_source='CFSR',
            stability_iterations=None,
            albedo_iterations=10,
            rs_interp_flag=True,
            ta_interp_flag=True,
            ta_smooth_flag=True,
            lat=None,
            lon=None,
            et_min=0.01,
            **kwargs
    ):
        """Initialize an image for computing DisALEXI

        Parameters
        ----------
        image : ee.Image
            Prepped image
        ta_source : {'CONUS_V003', 'CONUS_V004', 'CONUS_V005', 'CONUS_V006'}
            ALEXI scale air temperature image collection ID.
            The default is 'CONUS_V005'.
        alexi_source : {'CONUS_V003', 'CONUS_V004', 'CONUS_V005', 'CONUS_V006'}
            ALEXI ET image collection ID (the default is 'CONUS_V005').
        lai_source : string
            LAI image collection ID.
        tir_source : string
            Sharpened thermal infrared image collection ID.
        elevation_source: str, ee.Image
            Elevation source keyword or asset (the default is USGS/SRTMGL1_003).
            Units must be in meters.
        landcover_source : {'USGS/NLCD_RELEASES/2019_REL/NLCD',
                            'USGS/NLCD_RELEASES/2019_REL/NLCD/2016',
                            'USGS/NLCD_RELEASES/2016_REL',
                            'USGS/NLCD_RELEASES/2016_REL/2011',
                            'GLOBELAND30'}
            Land cover source keyword
            (the default is 'USGS/NLCD_RELEASES/2019_REL/NLCD').
        ta0_source : {'CFSR'}
            Air temperature source keyword (the default is 'CFSR').
        rs_daily_source : {'MERRA2', 'CFSR'}
            Daily solar insolation source keyword (the default is 'CFSR').
        rs_hourly_source : {'MERRA2', 'CFSR'}
            Hourly solar insolation source keyword (the default is 'CFSR').
        windspeed_source : {'CFSV2', 'CFSR'}
            Windspeed source keyword (the default is 'CFSR').
        vp_source : {'CFSR'}
            Vapour pressure source keyword (the default is 'CFSR').
        airpressure_source: {'CFSR'}
            Air pressure source keyword (the default is 'CFSR').
        stability_iterations : int, optional
            Number of istability calculation iterations.  If not set, the
            number will be computed dynamically.
        albedo_iterations : int, optional
            Number albedo separation iterations (the default is 10).
        rs_interp_flag : bool, optional
            If True, interpolate incoming solar radiation.
            If False, select image with same date and hour.
            The default is True.
        ta_interp_flag : bool, optional
            If True, interpolate between Ta step images (the default is True).
        ta_smooth_flag : bool, optional
            If True, smooth and resample Ta imagee (the default is True).
        lat : ee.Image, optional
            Latitude [deg].  If not set will default to ee.Image.pixelLonLat().
        lon : ee.Image, optional
            Longitude [deg].  If not set will default to ee.Image.pixelLonLat().
        et_min: float, optional
            Minimum output ET value. (the default is 0.01).
        kwargs: dict, optional
            et_reference_source : str, float
                Reference ET source (the default is None).
                Parameter is required if computing 'et_fraction' or 'et_reference'.
            et_reference_band : str
                Reference ET band name (the default is None).
                Parameter is required if computing 'et_fraction' or 'et_reference'.
            et_reference_factor : float, None
                Reference ET scaling factor.  The default is None which is
                equivalent to 1.0 (or no scaling).
            et_reference_resample : {'nearest', 'bilinear', 'bicubic', None}
                Reference ET resampling.  The default is None which is
                equivalent to nearest neighbor resampling.

        Notes
        -----
        For now defaulting all inputs to the CONUS based inputs.

        FilterDate looks at the time_starts, so if the ALEXI image
            has a start time of 0 UTC, to get the ALEXI image for the
            image date, you may need to move the image date back a day.

        """
        self.image = ee.Image(image)

        # Copy system properties
        self.id = self.image.get('system:id')
        self.index = self.image.get('system:index')
        self.time_start = self.image.get('system:time_start')
        self.properties = {
            'system:index': self.index,
            'system:time_start': self.time_start,
            'image_id': self.id,
        }

        # # Build SCENE_ID from the (possibly merged) system:index
        # sceneid = ee.List(ee.String(self.index).split('_')).slice(-3)
        # self.scene_id = ee.String(scene_id.get(0)).cat('_') \
        #     .cat(ee.String(scene_id.get(1))).cat('_') \
        #     .cat(ee.String(scene_id.get(2)))

        # # Build WRS2_TILE from the scene_id
        # self.wrs2_tile = ee.String('p').cat(self.scene_id.slice(5, 8)) \
        #     .cat('r').cat(self.scene_id.slice(8, 11))

        # Set server side date/time properties using the 'system:time_start'
        self.datetime = ee.Date(self.time_start)
        self.year = ee.Number(self.datetime.get('year'))
        self.month = ee.Number(self.datetime.get('month'))
        self.start_date = ee.Date(utils.date_to_time_0utc(self.datetime))
        self.end_date = self.start_date.advance(1, 'day')
        self.doy = ee.Number(self.datetime.getRelative('day', 'year')).add(1).int()
        # self.cycle_day = self.start_date.difference(
        #     ee.Date.fromYMD(1970, 1, 3), 'day').mod(8).add(1).int()

        # # Set server side date/time properties using the 'system:time_start'
        # self.datetime = ee.Date(image.get('system:time_start'))
        # self.date = ee.Date(self.datetime.format('yyyy-MM-dd'))
        # self.doy = ee.Number(self.datetime.getRelative('day', 'year')).add(1).double()
        self.hour = ee.Number(self.datetime.getFraction('day')).multiply(24)
        self.hour_int = self.hour.floor()
        # Time used in IDL is hours and fractional minutes (no seconds)
        # self.time = ee.Date(self.datetime).get('hour').add(
        #     ee.Date(self.datetime).get('minute').divide(60))

        # CGM - Applying cloud mask directly to input image
        #   In the IDL code the mask was applied to a_pt in main TSEB function
        # CGM - Buffering clouds a small amount
        # Buffer distances are currently hardcoded at -30 and +90
        #   but these could be changed to input parameters
        self.cloud_mask = image.select('cfmask').gte(1) \
            .reduceNeighborhood(ee.Reducer.min(), ee.Kernel.circle(30, 'meters')) \
            .reduceNeighborhood(ee.Reducer.max(), ee.Kernel.circle(90, 'meters')) \
            .eq(0) \
            .rename(['cloud_mask'])
        # self.cfmask = image.select('cfmask')
        # self.cloud_mask = self.cfmask.eq(0)
        input_image = image.updateMask(self.cloud_mask)

        # Get input bands from the image
        self.albedo = input_image.select('albedo')
        self.ndvi = input_image.select('ndvi')
        # LAT and LST are being read from source image collections
        # self.lai = input_image.select('lai')
        # self.lai = self.lai.where(lai.mask(), 0.01)
        # self.lst = input_image.select('lst')

        # Set input parameters
        self.ta_source = ta_source
        self.alexi_source = alexi_source
        self.lai_source = lai_source
        self.tir_source = tir_source
        self.elevation_source = elevation_source
        self.landcover_source = landcover_source
        self.ta0_source = ta0_source
        self.rs_daily_source = rs_daily_source
        self.rs_hourly_source = rs_hourly_source
        self.windspeed_source = windspeed_source
        self.vp_source = vp_source
        self.airpressure_source = airpressure_source
        if stability_iterations:
            self.stabil_iter = int(stability_iterations + 0.5)
        else:
            self.stabil_iter = None
        self.albedo_iter = int(albedo_iterations + 0.5)
        self.rs_interp_flag = utils.boolean(rs_interp_flag)
        self.ta_interp_flag = utils.boolean(ta_interp_flag)
        self.ta_smooth_flag = utils.boolean(ta_smooth_flag)
        self.et_min = et_min

        # Reference ET parameters
        try:
            self.et_reference_source = kwargs['et_reference_source']
        except:
            self.et_reference_source = None
        try:
            self.et_reference_band = kwargs['et_reference_band']
        except:
            self.et_reference_band = None
        try:
            self.et_reference_factor = kwargs['et_reference_factor']
        except:
            self.et_reference_factor = None
        try:
            self.et_reference_resample = kwargs['et_reference_resample']
        except:
            self.et_reference_resample = None

        # Check reference ET parameters
        if (self.et_reference_factor and
                not utils.is_number(self.et_reference_factor)):
            raise ValueError('et_reference_factor must be a number')
        if self.et_reference_factor and self.et_reference_factor < 0:
            raise ValueError('et_reference_factor must be greater than zero')
        resample_methods = ['nearest', 'bilinear', 'bicubic']
        if (self.et_reference_resample and
                self.et_reference_resample.lower() not in resample_methods):
            raise ValueError('Unsupported et_reference_resample method\n')

        if lat is None:
            self.lat = self.ndvi.multiply(0).add(
                ee.Image.pixelLonLat().select(['latitude']))
        elif utils.is_number(lat):
            self.lat = ee.Image.constant(lat)
        elif isinstance(lat, ee.computedobject.ComputedObject):
            self.lat = lat
        else:
            raise ValueError('invalid lat parameter')

        if lon is None:
            self.lon = self.ndvi.multiply(0).add(
                ee.Image.pixelLonLat().select(['longitude']))
        elif utils.is_number(lon):
            self.lon = ee.Image.constant(lon)
        elif isinstance(lon, ee.computedobject.ComputedObject):
            self.lon = lon
        else:
            raise ValueError('invalid lon parameter')

        # TODO: Move this into a landcover source function
        # Set default land cover image and type
        # For now default to CONUS and use default if image and type were not set
        if utils.is_number(self.landcover_source):
            self.lc_source = ee.Image.constant(int(self.landcover_source)) \
                .rename(['landcover'])
            self.lc_type = 'NLCD'
        elif isinstance(self.landcover_source, ee.computedobject.ComputedObject):
            # If the source is an ee.Image assume it is an NLCD image
            self.lc_source = self.landcover_source.rename(['landcover'])
            self.lc_type = 'NLCD'
        elif self.landcover_source == 'USGS/NLCD_RELEASES/2019_REL/NLCD':
            # If the image collection is passed in, assume the target is CONUS
            #   and select a close year
            # Making this a server side call to avoid a getInfo call on self.year
            # Clamp the year to 1999-2019 for now
            year_remap = ee.Dictionary({
                '1999': '2001', '2000': '2001', '2001': '2001', '2002': '2001',
                '2003': '2004', '2004': '2004', '2005': '2004',
                '2006': '2006', '2007': '2006',
                '2008': '2008', '2009': '2008',
                '2010': '2011', '2011': '2011', '2012': '2011',
                '2013': '2013', '2014': '2013',
                '2015': '2016', '2016': '2016', '2017': '2016',
                '2018': '2019', '2019': '2019',
            })
            nlcd_year = year_remap.get(self.year.min(2019).max(1999).format('%d'))
            # The timestart didn't seem quite right on the NLCD assets,
            #   so filtering using the system:index instead
            self.lc_source = ee.ImageCollection(self.landcover_source) \
                .filter(ee.Filter.equals('system:index', nlcd_year))\
                .first().select(['landcover'])
            self.lc_type = 'NLCD'
        elif self.landcover_source.upper().startswith('USGS/NLCD_RELEASES/2019_REL/NLCD/'):
            # Assume an image was passed in and use it directly
            self.lc_source = ee.Image(self.landcover_source.upper()) \
                .select(['landcover'])
            self.lc_type = 'NLCD'
        elif self.landcover_source in ['USGS/NLCD_RELEASES/2016_REL']:
            # If the image collection is passed in, assume the target is CONUS
            #   and select a close year
            year_remap = ee.Dictionary({
                '1999': 2001, '2000': 2001, '2001': 2001, '2002': 2001,
                '2003': 2004, '2004': 2004, '2005': 2004,
                '2006': 2006, '2007': 2006,
                '2008': 2008, '2009': 2008,
                '2010': 2011, '2011': 2011, '2012': 2011,
                '2013': 2013, '2014': 2013,
                '2015': 2016, '2016': 2016,
            })
            nlcd_year = year_remap.get(self.year.min(2016).max(1999).format('%d'))
            nlcd_date = ee.Date.fromYMD(nlcd_year, 1, 1)
            self.lc_source = ee.ImageCollection(self.landcover_source)\
                .filterDate(nlcd_date, nlcd_date.advance(1, 'year'))\
                .first().select(['landcover'])
            self.lc_type = 'NLCD'
        elif self.landcover_source.upper().startswith('USGS/NLCD_RELEASES/2016_REL/'):
            # Assume an image was passed in and use it directly
            self.lc_source = ee.Image(self.landcover_source.upper())\
                .select(['landcover'])
            self.lc_type = 'NLCD'
        # TODO: Test out the ESA 10m Landcover asset
        #   Collection ID: ESA/WorldCover/v100
        #   Only 2020 is currently available
        # elif self.landcover_source == 'ESA/WorldCover/v100/2020':
        #     self.lc_source = ee.Image(self.landcover_source).select(['Map'])
        #     self.lc_type = 'ESA?'
        # TODO: Test out the Copernicus 100m landcover asset
        #   Collection ID: COPERNICUS/Landcover/100m/Proba-V-C3/Global
        #   Years: 2015-2019
        # elif self.landcover_source == 'COPERNICUS/Landcover/100m/Proba-V-C3/Global':
        # elif self.landcover_source.startswith('COPERNICUS/Landcover/100m/Proba-V-C3/Global/2'):
        #     self.lc_source = ee.ImageCollection(self.landcover_source)
        #         .filterDate(ee.Date(self.year, 1, 1),
        #                     ee.Date(self.year, 12, 31))\
        #         .select(['discrete_classification])
        #     self.lc_type = ''
        # TODO: Eventually remove this option since the image collection is deprecated
        elif self.landcover_source.upper().startswith('USGS/NLCD/'):
            self.lc_source = ee.Image(self.landcover_source.upper())\
                .select(['landcover'])
            self.lc_type = 'NLCD'
        # TODO: Eventually remove the NLCD  keyword sources
        elif self.landcover_source.upper() == 'NLCD2016':
            self.lc_source = ee.Image('USGS/NLCD/NLCD2016').select(['landcover'])
            self.lc_type = 'NLCD'
        elif self.landcover_source.upper() == 'NLCD2011':
            self.lc_source = ee.Image('USGS/NLCD/NLCD2011').select(['landcover'])
            self.lc_type = 'NLCD'
        elif self.landcover_source.upper() == 'NLCD2006':
            self.lc_source = ee.Image('USGS/NLCD/NLCD2006').select(['landcover'])
            self.lc_type = 'NLCD'
        # TODO: Eventually remove this option
        elif self.landcover_source.upper() == 'GLOBELAND30':
            # GlobeLand30 values need to be set to the lowest even multiple of 10,
            #   since that is currently what is in the landcover.xlsx file.
            # http://www.globallandcover.com/GLC30Download/index.aspx
            lc_coll = ee.ImageCollection('users/cgmorton/GlobeLand30')\
                .filterBounds(image.geometry().bounds(1))
            self.lc_source = ee.Image(lc_coll.mosaic())\
                .divide(10).floor().multiply(10)\
                .rename(['landcover'])
            self.lc_type = 'GLOBELAND30'
        else:
            raise ValueError(f'Unsupported landcover_source: '
                             f'{self.landcover_source}\n')
        self.set_landcover_vars()

        # Image projection and geotransform
        self.crs = image.projection().crs()
        self.transform = ee.List(ee.Dictionary(
            ee.Algorithms.Describe(image.projection())).get('transform'))
        # self.crs = image.select([0]).projection().getInfo()['crs']
        # self.transform = image.select([0]).projection().getInfo()['transform']

        # CGM - This should probably be set in et_alexi() but that wasn't working
        if type(self.alexi_source) is str:
            if self.alexi_source.upper() == 'CONUS_V006':
                self.alexi_geo = [0.04, 0, -125.02, 0, -0.04, 49.78]
                self.alexi_crs = 'EPSG:4326'
            elif self.alexi_source.upper() == 'CONUS_V005':
                self.alexi_geo = [0.04, 0, -125.02, 0, -0.04, 49.78]
                self.alexi_crs = 'EPSG:4326'
            elif self.alexi_source.upper() == 'CONUS_V004':
                self.alexi_geo = [0.04, 0, -125.02, 0, -0.04, 49.78]
                self.alexi_crs = 'EPSG:4326'
            elif self.alexi_source.upper() == 'CONUS_V003':
                self.alexi_geo = [0.04, 0, -125.04, 0, -0.04, 49.8]
                self.alexi_crs = 'EPSG:4326'
            else:
                # Assume ALEXI source is an image collection ID if it is a string
                #   but doesn't match on any of the keywords.
                alexi_img = ee.Image(ee.ImageCollection(self.alexi_source).first())
                self.alexi_geo = ee.List(ee.Dictionary(
                    ee.Algorithms.Describe(alexi_img.projection())).get('transform'))
                self.alexi_crs = alexi_img.projection().crs()
        else:
            self.alexi_geo = [0.04, 0, -125.04, 0, -0.04, 49.8]
            self.alexi_crs = 'EPSG:4326'

    @classmethod
    def from_image_id(cls, image_id, **kwargs):
        """Constructs a DisALEXI Image instance from an image ID

        Parameters
        ----------
        image_id : str
            An earth engine image ID.
            (i.e. 'LANDSAT/LC08/C01/T1_SR/LC08_044033_20170716')
        kwargs
            Keyword arguments to pass through to model init.

        Returns
        -------
        new instance of Image class

        """
        # CGM - Should the supported image collection IDs and helper
        # function mappings be set in a property or method of the Image class?
        collection_methods = {
            'LANDSAT/LC09/C02/T1_L2': 'from_landsat_c2_sr',
            'LANDSAT/LC08/C02/T1_L2': 'from_landsat_c2_sr',
            'LANDSAT/LE07/C02/T1_L2': 'from_landsat_c2_sr',
            'LANDSAT/LT05/C02/T1_L2': 'from_landsat_c2_sr',
            'LANDSAT/LT04/C02/T1_L2': 'from_landsat_c2_sr',
            'LANDSAT/LC08/C01/T1_SR': 'from_landsat_c1_sr',
            'LANDSAT/LE07/C01/T1_SR': 'from_landsat_c1_sr',
            'LANDSAT/LT05/C01/T1_SR': 'from_landsat_c1_sr',
            'LANDSAT/LT04/C01/T1_SR': 'from_landsat_c1_sr',
        }

        try:
            method_name = collection_methods[image_id.rsplit('/', 1)[0]]
        except KeyError:
            raise ValueError(f'Unsupported collection ID: {image_id}\n')
        except Exception as e:
            raise Exception(f'unhandled exception: {e}')

        method = getattr(Image, method_name)

        return method(ee.Image(image_id), **kwargs)

    @classmethod
    def from_landsat_c1_sr(cls, sr_image, **kwargs):
        """Returns a DisALEXI Image instance from a Landsat Collection 1 SR image

        Parameters
        ----------
        sr_image : ee.Image
            A raw Landsat Collection 1 SR image.
        cloudmask_args : dict
            Keyword arguments to pass through to cloud mask function.
        kwargs : dict
            Keyword arguments to pass through to Image init function.

        Returns
        -------
        Instance of Image class

        """
        return cls(landsat.Landsat_C01_SR(sr_image).prep(), **kwargs)

    @classmethod
    def from_landsat_c2_sr(cls, sr_image, **kwargs):
        """Returns a DisALEXI Image instance from a Landsat Collection 2 SR image

        Parameters
        ----------
        sr_image : ee.Image
            A raw Landsat Collection 2 SR image.
        cloudmask_args : dict
            Keyword arguments to pass through to cloud mask function.
        kwargs : dict
            Keyword arguments to pass through to Image init function.

        Returns
        -------
        Instance of Image class

        """
        return cls(landsat.Landsat_C02_SR(sr_image).prep(), **kwargs)

    def calculate(self, variables=['et']):
        """Return a multiband image of calculated variables

        Parameters
        ----------
        variables : list

        Returns
        -------
        ee.Image

        """
        output_images = []
        for v in variables:
            if v.lower() == 'et':
                output_images.append(self.et.float())
            # elif v.lower() == 'et_fraction':
            #     output_images.append(self.et_fraction.float())
            elif v.lower() == 'et_reference':
                output_images.append(self.et_reference.float())
            elif v.lower() == 'lst':
                output_images.append(self.tir.float())
            elif v.lower() == 'mask':
                output_images.append(self.mask)
            elif v.lower() == 'ndvi':
                output_images.append(self.ndvi.float())
            # elif v.lower() == 'qa':
            #     output_images.append(self.qa)
            # elif v.lower() == 'quality':
            #     output_images.append(self.quality)
            elif v.lower() == 'time':
                output_images.append(self.time)
            else:
                raise ValueError(f'Unsupported variable: {v}\n')

        return ee.Image(output_images).set(self.properties)

    @lazy_property
    def et(self):
        """Compute Landsat scale DisALEXI ET

        Returns
        -------
        image : ee.Image
            DisALEXI ET image

        """
        et = tseb.tseb_pt(
            t_air=self.ta, t_rad=self.tir, t_air0=self.t_air0, e_air=self.vp,
            u=self.windspeed, p=self.pressure, z=self.elevation,
            rs_1=self.rs1, rs24=self.rs24, vza=0,
            aleafv=self.aleafv, aleafn=self.aleafn, aleafl=self.aleafl,
            adeadv=self.adeadv, adeadn=self.adeadn, adeadl=self.adeadl,
            albedo=self.albedo, ndvi=self.ndvi, lai=self.lai,
            clump=self.clump, leaf_width=self.leaf_width,
            hc_min=self.hc_min, hc_max=self.hc_max,
            datetime=self.datetime, lat=self.lat, lon=self.lon,
            stabil_iter=self.stabil_iter, albedo_iter=self.albedo_iter,
            et_min=self.et_min,
        )
        return et.rename(['et']).set(self.properties)

    @lazy_property
    def et_reference(self):
        """Compute reference ET for the image date"""
        if utils.is_number(self.et_reference_source):
            # Interpret numbers as constant images
            # CGM - Should we use the ee_types here instead?
            #   i.e. ee.ee_types.isNumber(self.etr_source)
            et_reference_img = ee.Image.constant(float(self.et_reference_source))
        elif type(self.et_reference_source) is str:
            # Assume a string source is an image collection ID (not an image ID)
            et_reference_img = ee.Image(
                ee.ImageCollection(self.et_reference_source)\
                    .filterDate(self.start_date, self.end_date)\
                    .select([self.et_reference_band])\
                    .first())
        else:
            raise ValueError(f'Unsupported etr_source: {self.etr_source}\n')

        # Map ETr values directly to the input (i.e. Landsat) image pixels
        # The benefit of this is the ETr image is now in the same crs as the
        #   input image.  Not all models may want this though.
        # CGM - Should the output band name match the input ETr band name?
        return self.ndvi.multiply(0).add(et_reference_img)\
            .multiply(self.et_reference_factor)\
            .rename(['et_reference']).set(self.properties)

    # @lazy_property
    # def et_fraction(self):
    #     """Compute ET fraction as actual ET divided by the reference ET"""
    #     return self.et.divide(self.et_reference) \
    #         .rename(['et_fraction']).set(self.properties)

    @lazy_property
    def et_alexi(self):
        """Extract ALEXI ET image for the target image time

        ALEXI ET is converted from MJ m-2 d-1 to mm d-1

        """
        alexi_keyword_sources = {
            'CONUS_V003': 'projects/disalexi/alexi/CONUS_V003',
            'CONUS_V004': 'projects/disalexi/alexi/CONUS_V004',
            'CONUS_V005': 'projects/ee-tulipyangyun-2/assets/alexi/ALEXI_V005',
            'CONUS_V006': 'projects/ee-tulipyangyun-2/assets/alexi/ALEXI_V006',
        }
        alexi_re = re.compile('(projects/earthengine-legacy/assets/)?'
                              'projects/disalexi/alexi/CONUS_V\\w+')

        if utils.is_number(self.alexi_source):
            # Interpret numbers as constant images
            alexi_img = ee.Image.constant(float(self.alexi_source))
        elif isinstance(self.alexi_source, ee.computedobject.ComputedObject):
            alexi_img = self.alexi_source
        elif self.alexi_source.upper() in alexi_keyword_sources.keys():
            alexi_coll_id = alexi_keyword_sources[self.alexi_source.upper()]
            alexi_coll = ee.ImageCollection(alexi_coll_id)\
                .filterDate(self.start_date, self.end_date)
            # TODO: Check if collection size is 0
            alexi_img = ee.Image(alexi_coll.first()).multiply(0.408)
        elif alexi_re.match(self.alexi_source):
            alexi_coll = ee.ImageCollection(self.alexi_source)\
                .filterDate(self.start_date, self.end_date)
            alexi_img = ee.Image(alexi_coll.first()).multiply(0.408)
        elif self.alexi_source in alexi_keyword_sources.values():
            # CGM - Quick fix for catching if the alexi_source was to as the
            #   collection ID, specifically for V005 since it is currently in a
            #   different project and won't get matched by the regex.
            alexi_coll = ee.ImageCollection(self.alexi_source)\
                .filterDate(self.start_date, self.end_date)
            alexi_img = ee.Image(alexi_coll.first()).multiply(0.408)
        else:
            raise ValueError(f'Unsupported alexi_source: {self.alexi_source}\n')

        return alexi_img.max(0.1).rename(['et_alexi'])

    @lazy_property
    def elevation(self):
        """Elevation"""
        if utils.is_number(self.elevation_source):
            elev_img = ee.Image.constant(float(self.elevation_source))
        elif isinstance(self.elevation_source, ee.computedobject.ComputedObject):
            elev_img = self.elevation_source
        elif type(self.elevation_source) is str:
            elev_img = ee.Image(self.elevation_source)
        else:
            raise ValueError(f'Unsupported elev_source: {self.elevation_source}\n')

        return elev_img.select([0], ['elevation'])

    @lazy_property
    def lai(self):
        """Leaf Area Index (LAI)"""
        if utils.is_number(self.lai_source):
             lai_img = ee.Image.constant(float(self.lai_source))
        # elif isinstance(self.lai_source, ee.computedobject.ComputedObject):
        #     lai_img = self.lai_source
        elif type(self.lai_source) is str:
             # Assumptions (for now)
             #   String lai_source is an image collection ID
             #   Images are single band and don't need a select()
             #   LAI images always need to be scaled
             # CGM - This will raise a .get() error if the image doesn't exist
             lai_coll = ee.ImageCollection(self.lai_source)\
                .filterMetadata('scene_id', 'equals', self.index)  # yun modified
             lai_img = ee.Image(lai_coll.first()).select(['LAI'])
             lai_img = lai_img.multiply(ee.Number(lai_img.get('scale_factor')))\
                 .set({'landsat_lai_version': lai_img.get('landsat_lai_version')})
             #lai_img = lai_img.multiply(ee.Number(0.01))\
             #    .set({'landsat_lai_version': lai_img.get('landsat_version')})  #yun to test nosat lai
             self.landsat_lai_version = lai_img.get('landsat_lai_version')
        else:
             raise ValueError(f'Unsupported lai_source: {self.lai_source}\n')

        return lai_img.select([0], ['lai'])

    # TODO: Rename "tir" to "lst" once Collection 2 is fully loaded
    @lazy_property
    def tir(self):
        """Sharpened thermal infrared (TIR)"""
        if utils.is_number(self.tir_source):
            tir_img = ee.Image.constant(float(self.tir_source))
        # elif isinstance(self.tir_source, ee.computedobject.ComputedObject):
        #     tir_img = self.tir_source
        elif type(self.tir_source) is str:
            # Assumptions (for now)
            #   String tir_source is an image collection ID
            #   Images are single band and don't need a select()
            #   TIR images always need to be scaled
            # CGM - This will raise a .get() error if the image doesn't exist
            tir_coll = ee.ImageCollection(self.tir_source)\
                .filterMetadata('scene_id', 'equals', self.index)
            tir_img = ee.Image(tir_coll.first())
            tir_img = tir_img.multiply(ee.Number(tir_img.get('scale_factor')))\
                .set({'sharpen_version': tir_img.get('sharpen_version')})
            self.sharpen_version = tir_img.get('sharpen_version')
        else:
            raise ValueError(f'Unsupported tir_source: {self.tir_source}\n')

        # tir_img = ee.Image('users/tulipyangyun/LST_20162091720').add(273.16)
        return tir_img.select([0], ['tir'])

    # CGM - This is not being called but maybe should be applied to just the
    #   Collection 1 TIR image (since it is not LST)
    # @lazy_property
    # def lst(self):
    #     """Return land surface temperature (LST) image"""
    #     emissivity = self.lai.divide(300).add(0.97)\
    #         .where(self.ndvi.lte(0), 0.99)\
    #         .where(self.ndvi.gt(0).And(self.lai.gt(3)), 0.98)
    #
    #     # Get properties from image
    #     k1 = ee.Number(ee.Image(self.image).get('k1_constant'))
    #     k2 = ee.Number(ee.Image(self.image).get('k2_constant'))
    #
    #     # First back out radiance from brightness temperature
    #     # Then recalculate emissivity corrected Ts
    #     thermal_rad_toa = self.tir.expression(
    #         'k1 / (exp(k2 / ts_brightness) - 1)',
    #         {'ts_brightness': self.tir, 'k1': k1, 'k2': k2},
    #     )
    #
    #     tnb = 0.866  # narrow band transmissivity of air
    #     rp = 0.91  # path radiance
    #     rsky = 1.32  # narrow band clear sky downward thermal radiation
    #     rc = thermal_rad_toa.expression(
    #         '((thermal_rad_toa - rp) / tnb) - ((1. - emiss) * rsky)',
    #         {
    #             'thermal_rad_toa': thermal_rad_toa,
    #             'emiss': emissivity,
    #             'rp': rp, 'tnb': tnb, 'rsky': rsky
    #         },
    #     )
    #     lst = rc.expression(
    #         'k2 / log(emiss * k1 / rc + 1)',
    #         {'emiss': emissivity, 'rc': rc, 'k1': k1, 'k2': k2},
    #     )
    #
    #     return lst.set(self.properties).rename(['lst'])

    @lazy_property
    def mask(self):
        """Mask of all active pixels (based on the final et)"""
        return self.et.multiply(0).add(1).updateMask(1)\
            .rename(['mask']).set(self.properties).uint8()

    @lazy_property
    def time(self):
        """Return an image of the 0 UTC time (in milliseconds)"""
        return self.mask\
            .double().multiply(0).add(utils.date_to_time_0utc(self.datetime))\
            .rename(['time']).set(self.properties)

    # @lazy_property
    # def albedo(self):
    #     """Return albedo image"""
    #     return self.image.select(['albedo']).set(self.properties)

    # @lazy_property
    # def ndvi(self):
    #     """Return NDVI image"""
    #     return self.image.select(['ndvi']).set(self.properties)

    # @lazy_property
    # def quality(self):
    #     """Set quality to 1 for all active pixels (for now)"""
    #     return self.mask\
    #         .rename(['quality']).set(self.properties)

    @lazy_property
    def ta(self):
        """Return the precomputed Ta asset for the target image

        Returns
        -------
        image : ee.Image
            Ta image

        """
        ta_keyword_sources = {
            'CONUS_V003': 'projects/openet/disalexi/tair/conus_v003_1k',
            'CONUS_V004': 'projects/openet/disalexi/tair/conus_v004_1k',
            'CONUS_V005': 'projects/openet/disalexi/tair/conus_v005_1k',
            'CONUS_V006': 'projects/openet/disalexi/tair/conus_v006_1k',
            'CONUS_V006B': 'projects/openet/disalexi/tair/conus_v006b_1k',
        }
        ta_source_re = re.compile(
            '(projects/earthengine-legacy/assets/)?'
            'projects/(\w+/)?disalexi/ta(ir)?/(conus|global)_v\d{3}\w?(_\w+)',
            re.IGNORECASE
        )

        # if self.ta_source is None:
        #     raise ValueError('ta_source must be set to compute et')
        if utils.is_number(self.ta_source):
            ta_img = ee.Image.constant(float(self.ta_source))
            #     .set({'ta_iteration': 'constant'})
        elif isinstance(self.ta_source, ee.computedobject.ComputedObject):
            ta_img = ee.Image(self.ta_source)
            #     .set({'ta_iteration': 'image'})
        elif self.ta_source.upper() in ta_keyword_sources.keys():
            ta_coll_id = ta_keyword_sources[self.ta_source.upper()]
            ta_coll = ee.ImageCollection(ta_coll_id)\
                .filterMetadata('image_id', 'equals', self.id)\
                .limit(1, 'step_size', False)
            input_img = ee.Image(ta_coll.first())

            # Select the Ta image with the minimum bias
            ta_array = input_img.select('step_\\d+_ta').toArray()
            bias_array = input_img.select('step_\\d+_bias').toArray()
            diff = bias_array.arraySlice(0, 1).subtract(bias_array.arraySlice(0, 0, -1))
            bias_array_mask = diff.abs().lt(0.001)
            # Make the array the same length by repeating the first value
            bias_array_mask = bias_array_mask.arraySlice(0, 0, 1).arrayCat(bias_array_mask, 0)
            bias_array_new = bias_array.add(bias_array_mask.multiply(99))
            index = bias_array_new.abs().multiply(-1).arrayArgmax()\
                .arraySlice(0, 0, 1).arrayFlatten([['array']])
            ta_img = ta_array.arrayGet(index)

            if self.ta_smooth_flag:
                # CGM - The v003 calculations used a radius of 2 instead of 1
                ta_img = ta_img.focal_mean(1, 'circle', 'pixels')\
                    .reproject(crs=self.alexi_crs, crsTransform=self.alexi_geo)\
                    .resample('bilinear')\
                    .reproject(crs=self.crs,crsTransform=self.transform)

        elif ta_source_re.match(self.ta_source):
            # CGM - How can I ensure it is an image collection ID?
            ta_coll = ee.ImageCollection(self.ta_source)\
                .filterMetadata('image_id', 'equals', self.id)\
                .limit(1, 'step_size', False)
            input_img = ee.Image(ta_coll.first())

            # Select the Ta image with the minimum bias
            ta_array = input_img.select('step_\\d+_ta').toArray()
            bias_array = input_img.select('step_\\d+_bias').toArray()
            diff = bias_array.arraySlice(0, 1).subtract(bias_array.arraySlice(0, 0, -1))
            bias_array_mask = diff.abs().lt(0.001)
            # Make the array the same length by repeating the first value
            bias_array_mask = bias_array_mask.arraySlice(0, 0, 1).arrayCat(bias_array_mask, 0)
            bias_array_new = bias_array.add(bias_array_mask.multiply(99))
            index = bias_array_new.abs().multiply(-1).arrayArgmax()\
                .arraySlice(0, 0, 1).arrayFlatten([['array']])
            ta_img = ta_array.arrayGet(index)

            if self.ta_smooth_flag:
                # CGM - The v003 calculations used a radius of 2 instead of 1
                ta_img = ta_img.focal_mean(1, 'circle', 'pixels')\
                    .reproject(crs=self.alexi_crs, crsTransform=self.alexi_geo)\
                    .resample('bilinear')\
                    .reproject(crs=self.crs,crsTransform=self.transform)

        else:
            raise ValueError(f'Unsupported ta_source: {self.ta_source}\n')

        # elif re.match('projects/(\w+/)?disalexi/ta/conus_v', self.ta_source, re.I):

        return ta_img.rename(['ta']).set(self.properties)

    @lazy_property
    def pressure(self):
        """Air pressure [kPa]"""
        if utils.is_number(self.airpressure_source):
            ap_img = ee.Image.constant(float(self.airpressure_source))
        elif isinstance(self.airpressure_source, ee.computedobject.ComputedObject):
            ap_img = self.airpressure_source
        elif self.airpressure_source.upper() == 'ESTIMATE':
            ap_img = self.elevation.expression(
                '101.3 * (((293.0 - 0.0065 * z) / 293.0) ** 5.26)',
                {'z': self.elevation})
        elif self.airpressure_source.upper() == 'CFSR':
            ap_coll_id = 'projects/disalexi/meteo_data/airpressure/global_v001_3hour'
            ap_coll = ee.ImageCollection(ap_coll_id).select(['airpressure'])
            ap_img = utils.interpolate(ap_coll, self.datetime, timestep=3)
            ap_img = ee.Image(ap_img)\
                .multiply(self.elevation.multiply(-0.0065).add(293)
                          .divide(293.0).pow(5.26))\
                .resample('bicubic')\
                .reproject(crs=self.alexi_crs, crsTransform=self.alexi_geo)\
                .resample('bilinear')\
                .reproject(crs=self.crs, crsTransform=self.transform)
        else:
            raise ValueError(f'Invalid airpressure_source: '
                             f'{self.airpressure_source}\n')

        return ee.Image(ap_img).rename(['pressure'])

    @lazy_property
    def rs1(self):
        """Hourly solar insolation [W m-2]"""
        if utils.is_number(self.rs_hourly_source):
            rs1_img = ee.Image.constant(float(self.rs_hourly_source))
        elif isinstance(self.rs_hourly_source, ee.computedobject.ComputedObject):
            rs1_img = self.rs_hourly_source
        elif self.rs_hourly_source.upper() == 'CFSR':
            rs1_coll_id = 'projects/disalexi/insol_data/global_v001_hourly'
            rs1_coll = ee.ImageCollection(rs1_coll_id).select(['insolation'])
            if self.rs_interp_flag:
                # TODO: Check if the CFSR instolation are instantaneous or accumulations
                rs1_img = utils.interpolate(rs1_coll, self.datetime, timestep=1)
            else:
                # Select the source image before the image time
                # This would be off if image time was "exactly" on the hour
                # Could change the offsets to -0.5 hours to select the closest image
                rs1_img = rs1_coll\
                    .filterDate(self.datetime.advance(-1, 'hour'),
                                self.datetime.advance(0, 'hour'))\
                    .first()
            rs1_img = ee.Image(rs1_img)\
                .resample('bicubic')\
                .reproject(crs=self.alexi_crs, crsTransform=self.alexi_geo)\
                .resample('bilinear')\
                .reproject(crs=self.crs, crsTransform=self.transform)
        elif self.rs_hourly_source.upper() == 'MERRA2':
            # Extract MERRA2 solar images for the target image time
            # Interpolate rs hourly image at image time
            # Hourly Rs is time average so time starts are 30 minutes early
            # Move image time 30 minutes earlier to simplify filtering/interpolation
            # This interpolation scheme will only work for hourly data
            rs1_coll = ee.ImageCollection('projects/disalexi/merra2/hourly')\
                .select(['SWGDNCLR'])
            if self.rs_interp_flag:
                # rs1_img = utils.interpolate(
                #     rs1_coll, self.datetime.advance(-0.5, 'hour'), timestep=1)
                interp_dt = self.datetime.advance(-0.5, 'hour')
                rs_a_img = ee.Image(rs1_coll.filterDate(
                    interp_dt.advance(-1, 'hour'), interp_dt).first())
                rs_b_img = ee.Image(rs1_coll.filterDate(
                    interp_dt, interp_dt.advance(1, 'hour')).first())
                t_a = ee.Number(rs_a_img.get('system:time_start'))
                t_b = ee.Number(rs_b_img.get('system:time_start'))
                rs1_img = rs_b_img.subtract(rs_a_img)\
                    .multiply(interp_dt.millis().subtract(t_a)
                              .divide(t_b.subtract(t_a)))\
                    .add(rs_a_img)
            else:
                rs1_img = ee.Image(
                    ee.ImageCollection(rs1_coll)\
                        .filterDate(self.start_date, self.end_date)\
                        .filter(ee.Filter.calendarRange(
                            self.hour_int, self.hour_int, 'hour'))
                        .first())
        else:
            raise ValueError(f'Unsupported rs_hourly_source: '
                             f'{self.rs_hourly_source}\n')

        return ee.Image(rs1_img).rename(['rs'])

    @lazy_property
    def rs24(self):
        """Daily (24 hour) solar insolation [W m-2]"""
        if utils.is_number(self.rs_daily_source):
            rs24_img = ee.Image.constant(float(self.rs_daily_source))
        elif isinstance(self.rs_daily_source, ee.computedobject.ComputedObject):
            rs24_img = self.rs_daily_source
        elif self.rs_daily_source.upper() == 'MERRA2':
            rs24_coll = ee.ImageCollection('projects/climate-engine/merra2/daily')\
                .select(['SWGDNCLR'])\
                .filterDate(self.start_date, self.end_date)
            rs24_img = ee.Image(rs24_coll.first())
        elif self.rs_daily_source.upper() == 'CFSR':
            rs24_coll_id = 'projects/disalexi/insol_data/global_v001_hourly'
            rs24_coll = ee.ImageCollection(rs24_coll_id)\
                .filterDate(self.datetime.advance(-8, 'hours'),
                            self.datetime.advance(12, 'hours'))
            rs24_img = ee.Image(rs24_coll.sum())\
                .resample('bicubic')\
                .reproject(crs=self.alexi_crs, crsTransform=self.alexi_geo)\
                .resample('bilinear')\
                .reproject(crs=self.crs, crsTransform=self.transform)
        else:
            raise ValueError(f'Unsupported rs_daily_source: {self.rs_daily_source}\n')

        return ee.Image(rs24_img).rename(['rs'])

    # TODO: Change the naming of this parameter, it is a little inconsistent.
    #   The band is called "tair0", but the property is "t_air0".
    # It might make more sense to change "ta" to something else and rename
    #   this to "ta".
    @lazy_property
    def t_air0(self):
        """Air temperature [K]"""
        if utils.is_number(self.ta0_source):
            tair0_img = ee.Image.constant(float(self.ta0_source))
        # elif isinstance(self.ta0_source, ee.computedobject.ComputedObject):
        #     tair0_img = self.ta0_source
        elif self.ta0_source.upper() == 'CFSR':
            tair0_coll_id = 'projects/disalexi/meteo_data/airtemperature/global_v001_3hour'
            tair0_coll = ee.ImageCollection(tair0_coll_id).select(['temperature'])
            tair0_img = utils.interpolate(tair0_coll, self.datetime, timestep=3)
            tair0_img = ee.Image(tair0_img)\
                .resample('bicubic')\
                .reproject(crs=self.alexi_crs, crsTransform=self.alexi_geo)\
                .resample('bilinear')\
                .reproject(crs=self.crs, crsTransform=self.transform)
        else:
            raise ValueError(f'Unsupported ta0_source: {self.ta0_source}\n')

        return ee.Image(tair0_img).rename(['tair0'])

    @lazy_property
    def windspeed(self):
        """Windspeed [m/s]"""
        if utils.is_number(self.windspeed_source):
            windspeed_img = ee.Image.constant(float(self.windspeed_source))
        elif isinstance(self.windspeed_source, ee.computedobject.ComputedObject):
            windspeed_img = self.windspeed_source
        elif self.windspeed_source.upper() == 'CFSV2':
            # It would be more correct to compute the magnitude for each image,
            #   then compute the average.
            # Do we need daily, 6hr, or interpolated instantaneous data?
            windspeed_coll = ee.ImageCollection('NOAA/CFSV2/FOR6H')\
                .select(['u-component_of_wind_height_above_ground',
                         'v-component_of_wind_height_above_ground'])\
                .filterDate(self.start_date, self.end_date)
            windspeed_img = windspeed_coll.mean()\
                .expression('sqrt(b(0) ** 2 + b(1) ** 2)')
        elif self.windspeed_source.upper() == 'CFSR':
            windspeed_coll_id = 'projects/disalexi/meteo_data/windspeed/global_v001_3hour'
            windspeed_coll = ee.ImageCollection(windspeed_coll_id)\
                .select(['windspeed'])
            windspeed_img = utils.interpolate(windspeed_coll, self.datetime, timestep=3)
            windspeed_img = ee.Image(windspeed_img)\
                .resample('bicubic')\
                .reproject(crs=self.alexi_crs, crsTransform=self.alexi_geo)\
                .resample('bilinear')\
                .reproject(crs=self.crs, crsTransform=self.transform)
        else:
            raise ValueError(f'Unsupported windspeed_source: '
                             f'{self.windspeed_source}\n')

        return ee.Image(windspeed_img).max(2).min(20).rename(['windspeed'])

    @lazy_property
    def vp(self):
        """Vapor pressure [kPa]"""
        if utils.is_number(self.vp_source):
            vp_img = ee.Image.constant(float(self.vp_source))
        elif isinstance(self.vp_source, ee.computedobject.ComputedObject):
            vp_img = self.vp_source
        elif self.vp_source.upper() == 'CFSR':
            vp_coll_id = 'projects/disalexi/meteo_data/vp/global_v001_3hour'
            vp_coll = ee.ImageCollection(vp_coll_id).select(['vp'])
            vp_img = utils.interpolate(vp_coll, self.datetime, timestep=3)
            vp_img = ee.Image(vp_img)\
                .resample('bicubic')\
                .reproject(crs=self.alexi_crs, crsTransform=self.alexi_geo)\
                .resample('bilinear')\
                .reproject(crs=self.crs, crsTransform=self.transform)
        else:
            raise ValueError(f'Unsupported vp_source: {self.vp_source}\n')

        return ee.Image(vp_img).rename(['vp'])

    def set_landcover_vars(self):
        """Compute Land Cover / LAI derived variables

        Eventually add code to fall back on default values
            aleafv: 0.9, aleafn: 0.9, aleafl: 0.9
            adeadv: 0.2, adeadn: 0.2, adeadl: 0.2
            hc_min: 0.1, hc_max: 0.5, xl: 0.5, clump: 0.99

        Parameters
        ----------
        lai : ee.Image
            Leaf area index
        landcover : ee.Image
            Landcover
        lc_type : string
            Landcover type (choices are "NLCD" or "GlobeLand30")

        """
        if self.lc_type.upper() not in remaps.keys():
            raise KeyError(f'Invalid lc_type (choices are {", ".join(remaps.keys())})')

        def lc_remap(landcover, lc_type, lc_var):
            """Remap land cover values to target values

            Parameters
            ----------
            landcover : ee.Image
                Land cover
            lc_type : string
                Land cover type (choices are "NLCD" or "GlobeLand30")
            lc_var: string
                Target variable

            Returns
            -------
            remap_img : ee.Image
            """
            lc_items = sorted(remaps[lc_type.upper()][lc_var].items())
            input_values = [k for k, v in lc_items]
            # Scale output values by 100 since remap values must be integer
            output_values = [v * 100 for k, v in lc_items]

            # Get the remap values from the dataframe
            #   and apply to the land cover image
            return ee.Image(landcover)\
                .remap(input_values, output_values)\
                .divide(100)\
                .rename([lc_var])

        # Get LC based variables
        self.aleafv = lc_remap(self.lc_source, self.lc_type, 'aleafv')
        self.aleafn = lc_remap(self.lc_source, self.lc_type, 'aleafn')
        self.aleafl = lc_remap(self.lc_source, self.lc_type, 'aleafl')
        self.adeadv = lc_remap(self.lc_source, self.lc_type, 'adeadv')
        self.adeadn = lc_remap(self.lc_source, self.lc_type, 'adeadn')
        self.adeadl = lc_remap(self.lc_source, self.lc_type, 'adeadl')
        # Yun modified to use hmax for hmin too!!!
        self.hc_min = lc_remap(self.lc_source, self.lc_type, 'hmax')
        self.hc_max = lc_remap(self.lc_source, self.lc_type, 'hmax')
        self.leaf_width = lc_remap(self.lc_source, self.lc_type, 'xl')
        self.clump = lc_remap(self.lc_source, self.lc_type, 'omega')

    # Helper functions for computing Ta (and ET)
    def ta_coarse(self, ta_img, threshold=0.5):
        """Compute the air temperature for each ALEXI ET cell that minimizes
        the bias between Landsat scale ET and ALEXI ET for a fixed range of
        air temperature values

        Parameters
        ----------
        ta_img : ee.Image
        threshold : float, optional

        Returns
        -------
        image : ee.Image
            ALEXI scale air temperature image

        """
        et_fine = tseb.tseb_pt(
            t_air=ta_img, t_rad=self.tir, t_air0=self.t_air0,
            e_air=self.vp, u=self.windspeed, p=self.pressure, z=self.elevation,
            rs_1=self.rs1, rs24=self.rs24, vza=0,
            aleafv=self.aleafv, aleafn=self.aleafn, aleafl=self.aleafl,
            adeadv=self.adeadv, adeadn=self.adeadn, adeadl=self.adeadl,
            albedo=self.albedo, ndvi=self.ndvi, lai=self.lai,
            clump=self.clump, leaf_width=self.leaf_width,
            hc_min=self.hc_min, hc_max=self.hc_max,
            datetime=self.datetime, lat=self.lat, lon=self.lon,
            stabil_iter=self.stabil_iter, albedo_iter=self.albedo_iter,
        )

        # Aggregate the Landsat scale ET up to the ALEXI scale
        et_coarse = ee.Image(et_fine)\
            .reproject(crs=self.crs, crsTransform=self.transform)\
            .reduceResolution(reducer=ee.Reducer.mean().unweighted(),
                              maxPixels=30000)\
            .reproject(crs=self.alexi_crs, crsTransform=self.alexi_geo)
        et_mask = et_coarse.mask().gte(threshold)
        bias = et_coarse.subtract(self.et_alexi)
        return ee.Image([ta_img, bias, et_coarse])\
            .updateMask(et_mask)\
            .rename(['ta', 'bias', 'et'])

    def ta_qm(self, ta_img, step_size, step_count, threshold=0.5):
        """Compute the air temperature for each ALEXI ET cell that minimizes
        the bias between Landsat scale ET and ALEXI ET for a fixed range of
        air temperature values

        Parameters
        ----------
        ta_img : ee.Image
        step_size : float
            The size of each Ta step.
        step_count : int
            The number of Ta steps.
        threshold : float, optional

        Returns
        -------
        image : ee.Image
            ALEXI scale air temperature image

        """

        # CGM - This is basically identical to ta_coarse()
        # Can't use ta_coarse since a mapped function can only have one input
        # Redefining it here allows threshold to be passed in
        def ta_func(ta):
            """Compute TSEB ET for the target t_air value"""
            et_fine = tseb.tseb_pt(
                t_air=ta, t_rad=self.tir, t_air0=self.t_air0,
                e_air=self.vp, u=self.windspeed, p=self.pressure,
                z=self.elevation, rs_1=self.rs1, rs24=self.rs24, vza=0,
                aleafv=self.aleafv, aleafn=self.aleafn, aleafl=self.aleafl,
                adeadv=self.adeadv, adeadn=self.adeadn, adeadl=self.adeadl,
                albedo=self.albedo, ndvi=self.ndvi, lai=self.lai,
                clump=self.clump, leaf_width=self.leaf_width,
                hc_min=self.hc_min, hc_max=self.hc_max,
                datetime=self.datetime, lat=self.lat, lon=self.lon,
                stabil_iter=self.stabil_iter, albedo_iter=self.albedo_iter,
            )

            # Aggregate the Landsat scale ET up to the ALEXI scale
            et_coarse = ee.Image(et_fine)\
                .reproject(crs=self.crs, crsTransform=self.transform)\
                .reduceResolution(reducer=ee.Reducer.mean().unweighted(),
                                  maxPixels=30000)\
                .reproject(crs=self.alexi_crs, crsTransform=self.alexi_geo)
            et_mask = et_coarse.mask().gte(threshold)
            bias = et_coarse.subtract(self.alexi_et)
            # Invert the abs(bias) since quality mosaic sorts descending
            qm_bias = bias.abs().multiply(-1)

            return ee.Image([qm_bias, ta, bias, et_coarse])\
                .updateMask(et_mask)\
                .rename(['qm_bias', 'ta', 'bias', 'et'])

        # Get output for a range of Ta values
        # Mapping over the list seemed a little slower than the FC
        ta_coll = ee.ImageCollection([
            ta_img.add(i * step_size)
            for i in range(step_count + 1)])

        # Return the air temperature associated with the lowest difference
        #   in ET from the ALEXI ET value
        return ee.ImageCollection(ta_coll.map(ta_func))\
            .qualityMosaic('qm_bias')\
            .select(['ta', 'bias', 'et'])

    def ta_mosaic(self, ta_img, step_size, step_count, threshold=0.5):
        """

        Parameters
        ----------
        ta_img : ee.Image
        step_size : float
            The size of each Ta step.
        step_count : int
            The number of Ta steps.
        threshold : float, optional

        Returns
        -------
        image : ee.Image
            ALEXI scale air temperature image

        """

        # CGM - This is basically identical to ta_coarse()
        # Can't use ta_coarse since a mapped function can only have one input
        # Redefining it here allows threshold to be passed in
        def ta_func(ta):
            """Compute TSEB ET for the target t_air value"""
            et_fine = tseb.tseb_pt(
                t_air=ta, t_rad=self.tir, t_air0=self.t_air0,
                e_air=self.vp, u=self.windspeed, p=self.pressure,
                z=self.elevation, rs_1=self.rs1, rs24=self.rs24, vza=0,
                aleafv=self.aleafv, aleafn=self.aleafn, aleafl=self.aleafl,
                adeadv=self.adeadv, adeadn=self.adeadn, adeadl=self.adeadl,
                albedo=self.albedo, ndvi=self.ndvi, lai=self.lai,
                clump=self.clump, leaf_width=self.leaf_width,
                hc_min=self.hc_min, hc_max=self.hc_max,
                datetime=self.datetime, lat=self.lat, lon=self.lon,
                stabil_iter=self.stabil_iter, albedo_iter=self.albedo_iter,
            )

            # Aggregate the Landsat scale ET up to the ALEXI scale
            et_coarse = ee.Image(et_fine)\
                .reproject(crs=self.crs, crsTransform=self.transform)\
                .reduceResolution(reducer=ee.Reducer.mean().unweighted(),
                                  maxPixels=30000)\
                .reproject(crs=self.alexi_crs, crsTransform=self.alexi_geo)
            et_agg = ee.Image(et_fine)\
                .reproject(crs=self.crs, crsTransform=self.transform)\
                .reduceResolution(reducer=ee.Reducer.count()
                                      .combine(reducer2=ee.Reducer.countEvery(),
                                               outputPrefix='all',
                                               sharedInputs=False),
                                  maxPixels=30000)\
                .reproject(crs=self.alexi_crs, crsTransform=self.alexi_geo)
            et_perc = et_agg.select(['et_count']).divide(et_agg.select(['et_allcount']))
            et_mask = et_perc.gt(threshold)
            #et_mask = et_coarse.mask().gte(threshold)
            bias = et_coarse.subtract(self.et_alexi)
            return ee.Image([ta, bias]).updateMask(et_mask)\
                .rename(['ta', 'bias'])\
                .set({'system:index': ee.Image(ta).get('system:index')})

        # CGM - Adding one extra Ta step and beginning and end to handle
        #   rounding or reduceResolution/projection error
        ta_coll = ee.ImageCollection([
            ta_img.add(j * step_size).set({'system:index': 'step_{}'.format(i)})
            for i, j in enumerate(range(-step_count // 2, step_count // 2 + 1))])
        # ta_coll = ee.ImageCollection([
        #     ta_img.add(i * step_size).set({'system:index': 'step_{}'.format(i)})
        #     for i in range(step_count + 1)])

        return ee.ImageCollection(ta_coll.map(ta_func)).toBands()

    def et_coarse(self, ta_img, threshold=0.5):
        """Compute the Landsat ET summed to the ALEXI grid for the
        specified air temperature image

        Parameters
        ----------
        ta_img : ee.Image
        threshold : float, optional

        Returns
        -------
        image : ee.Image
            ALEXI scale ET image

        """

        # ta_fine = self.lst.multiply(0).add(ta).rename(['ta'])
        # ta_coarse = self.alexi_et.multiply(0).add(ta).rename(['ta'])

        # Aggregate the Landsat scale ET up to the ALEXI scale
        et_fine = tseb.tseb_pt(
            t_air=ta_img, t_rad=self.tir, t_air0=self.t_air0,
            e_air=self.vp, u=self.windspeed, p=self.pressure, z=self.elevation,
            rs_1=self.rs1, rs24=self.rs24, vza=0,
            aleafv=self.aleafv, aleafn=self.aleafn, aleafl=self.aleafl,
            adeadv=self.adeadv, adeadn=self.adeadn, adeadl=self.adeadl,
            albedo=self.albedo, ndvi=self.ndvi, lai=self.lai,
            clump=self.clump, leaf_width=self.leaf_width,
            hc_min=self.hc_min, hc_max=self.hc_max,
            datetime=self.datetime, lat=self.lat, lon=self.lon,
            stabil_iter=self.stabil_iter, albedo_iter=self.albedo_iter,
        )
        et_coarse = et_fine\
            .reproject(crs=self.crs, crsTransform=self.transform)\
            .reduceResolution(reducer=ee.Reducer.mean().unweighted(),
                              maxPixels=30000)\
            .reproject(crs=self.alexi_crs, crsTransform=self.alexi_geo)\
            .rename(['et'])
        et_mask = et_coarse.mask().gte(threshold)
        return et_coarse.updateMask(et_mask)

    def et_bias(self, et_coarse_img):
        """Compute the bias between the computed ET and the ALEXI ET"""
        return et_coarse_img.subtract(self.et_alexi).rename(['bias'])

    def et_fine(self, ta_img):
        """Compute Landsat scale DisALEXI ET

        Parameters
        ----------
        ta_img : ee.Images

        Returns
        -------
        image : ee.Image
            DisALEXI ET image

        """
        et = tseb.tseb_pt(
            t_air=ta_img, t_rad=self.tir, t_air0=self.t_air0,
            e_air=self.vp, u=self.windspeed, p=self.pressure, z=self.elevation,
            rs_1=self.rs1, rs24=self.rs24, vza=0,
            aleafv=self.aleafv, aleafn=self.aleafn, aleafl=self.aleafl,
            adeadv=self.adeadv, adeadn=self.adeadn, adeadl=self.adeadl,
            albedo=self.albedo, ndvi=self.ndvi, lai=self.lai,
            clump=self.clump, leaf_width=self.leaf_width,
            hc_min=self.hc_min, hc_max=self.hc_max,
            datetime=self.datetime, lat=self.lat, lon=self.lon,
            stabil_iter=self.stabil_iter, albedo_iter=self.albedo_iter,
        )
        return et.rename(['et']).set(self.properties)
        #     .reproject(crs=self.crs, crsTransform=self.transform)
