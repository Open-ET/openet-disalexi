# from builtins import input
# import datetime
import pprint

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
            ta_source='CONUS_V001',
            alexi_source='CONUS_V001',
            elevation_source='USGS/SRTMGL1_003',
            landcover_source='NLCD2011',
            rs_daily_source='MERRA2',
            rs_hourly_source='MERRA2',
            windspeed_source='CFSV2',
            stabil_iterations=6,
            albedo_iterations=3,
            ta_interp_flag=True,
            rs_interp_flag=True,
        ):
        """Initialize an image for computing DisALEXI

        Parameters
        ----------
        image : ee.Image
            Prepped image
        ta_source : {'CONUS_V001'}
            ALEXI scale air temperature image collection ID (the default is
            'CONUS_V001').
        alexi_source : {'CONUS_V001'}
            ALEXI ET image collection ID (the default is 'CONUS_V001').
        elevation_source: str, ee.Image
            Elevation source keyword or asset (the default is USGS/SRTMGL1_003).
            Units must be in meters.
        landcover_source : {'NLCD2011', 'NLCD2006', 'GLOBELAND30'}
            Land cover source keyword (the default is 'NLCD2011').
        rs_daily_source : {'MERRA2'}
            Daily solar insolation source keyword (the default is 'MERRA2').
        rs_hourly_source : {'MERRA2'}
            Hourly solar insolation source keyword (the default is 'MERRA2').
        windspeed_source : {'CFSV2}
            Windspeed source keyword (the default is 'CFSV2').
        stabil_iterations : int
            Number of istability calculation iterations (the default is 10).
        albedo_iterations : int
            Number albedo separation iterations (the default is 10).
        ta_interp_flag : bool
        rs_interp_flag : bool

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
        self.date = ee.Date(self.time_start)
        # self.year = ee.Number(self.date.get('year'))
        # self.month = ee.Number(self.date.get('month'))
        # self.start_date = ee.Date(utils.date_to_time_0utc(self.date))
        # self.end_date = self.start_date.advance(1, 'day')
        # self.doy = ee.Number(self.date.getRelative('day', 'year')).add(1).int()
        # self.cycle_day = self.start_date.difference(
        #     ee.Date.fromYMD(1970, 1, 3), 'day').mod(8).add(1).int()

        # # Set server side date/time properties using the 'system:time_start'
        # self.datetime = ee.Date(image.get('system:time_start'))
        # self.date = ee.Date(self.datetime.format('yyyy-MM-dd'))
        # self.doy = ee.Number(self.datetime.getRelative('day', 'year')).add(1).double()
        self.hour = ee.Number(self.date.getFraction('day')).multiply(24)
        self.hour_int = self.hour.floor()
        # Time used in IDL is hours and fractional minutes (no seconds)
        # self.time = ee.Date(self.datetime).get('hour').add(
        #     ee.Date(self.datetime).get('minute').divide(60))

        # CGM - Applying cloud mask directly to input image
        #   In the IDL code the mask was applied to a_pt in main TSEB function
        # CGM - Buffering clouds a small amount
        # Buffer distances are currently hardcoded at -30 and +90
        #   but these could be changed to input parameters
        self.mask = image.select('cfmask').gte(1)\
            .reduceNeighborhood(ee.Reducer.min(), ee.Kernel.circle(30, 'meters'))\
            .reduceNeighborhood(ee.Reducer.max(), ee.Kernel.circle(90, 'meters'))\
            .eq(0)\
            .rename(['mask'])
        # self.cfmask = image.select('cfmask')
        # self.mask = self.cfmask.eq(0)
        input_image = image.updateMask(self.mask)

        # Get input bands from the image
        self.albedo = input_image.select('albedo')
        self.lai = input_image.select('lai')
        # self.lai = self.lai.where(lai.mask(), 0.01)
        self.lst = input_image.select('lst')
        self.ndvi = input_image.select('ndvi')

        # Set input parameters
        self.ta_source = ta_source
        self.alexi_source = alexi_source
        self.elevation_source = elevation_source
        self.landcover_source = landcover_source
        self.rs_daily_source = rs_daily_source
        self.rs_hourly_source = rs_hourly_source
        self.windspeed_source = windspeed_source
        self.stabil_iter = int(stabil_iterations)
        self.albedo_iter = int(albedo_iterations)
        self.ta_interp_flag = ta_interp_flag
        self.rs_interp_flag = rs_interp_flag

        # Convert flags from string to bool if necessary
        # This should be moved to a utils function
        if type(self.ta_interp_flag) is str:
            if self.ta_interp_flag.upper() in ['TRUE', 'T']:
                self.ta_interp_flag = True
            elif self.ta_interp_flag.upper() in ['FALSE', 'F']:
                self.ta_interp_flag = False
            else:
                raise ValueError(
                    'ta_interp_flag "{}" could not be interpreted as '
                    'bool'.format(self.ta_interp_flag))
        if type(self.rs_interp_flag) is str:
            if self.rs_interp_flag.upper() in ['TRUE', 'T']:
                self.rs_interp_flag = True
            elif self.rs_interp_flag.upper() in ['FALSE', 'F']:
                self.rs_interp_flag = False
            else:
                raise ValueError(
                    'rs_interp_flag "{}" could not be interpreted as '
                    'bool'.format(self.rs_interp_flag))

        # Set default land cover image and type
        # For now default to CONUS and use default if image and type were not set
        # GlobeLand30 values need to be set to the lowest even multiple of 10,
        #   since that is currently what is in the landcover.xlsx file.
        # http://www.globallandcover.com/GLC30Download/index.aspx
        if utils.is_number(self.landcover_source):
            self.lc_source = ee.Image.constant(int(self.landcover_source))\
                .rename(['landcover'])
            self.lc_type = 'NLCD'
        elif isinstance(self.landcover_source, ee.computedobject.ComputedObject):
            self.lc_source = self.landcover_source.rename(['landcover'])
            self.lc_type = 'NLCD'
        elif self.landcover_source.upper() == 'NLCD2011':
            self.lc_source = ee.Image('USGS/NLCD/NLCD2011').select(['landcover'])
            self.lc_type = 'NLCD'
        elif self.landcover_source.upper() == 'NLCD2006':
            self.lc_source = ee.Image('USGS/NLCD/NLCD2006').select(['landcover'])
            self.lc_type = 'NLCD'
        elif self.landcover_source.upper() == 'GLOBELAND30':
            lc_coll =  ee.ImageCollection('users/cgmorton/GlobeLand30')\
                .filterBounds(image.geometry().bounds(1))
            self.lc_source = ee.Image(lc_coll.mosaic())\
                .divide(10).floor().multiply(10)\
                .rename(['landcover'])
            self.lc_type = 'GLOBELAND30'
        else:
            raise ValueError('unsupported landcover_source: {}'.format(
                self.landcover_source))
        self.set_landcover_vars()

        # Image projection and geotransform
        self.crs = image.projection().crs()
        self.transform = ee.List(ee.Dictionary(
            ee.Algorithms.Describe(image.projection())).get('transform'))
        # self.crs = image.select([0]).projection().getInfo()['crs']
        # self.transform = image.select([0]).projection().getInfo()['transform']

        # CGM - This should probably be set in et_alexi() but that wasn't working
        if self.alexi_source.upper() == 'CONUS_V001':
            self.alexi_geo = [0.04, 0, -125.04, 0, -0.04, 49.8]
            self.alexi_crs = 'EPSG:4326'
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
        # DEADBEEF - Should the supported image collection IDs and helper
        # function mappings be set in a property or method of the Image class?
        collection_methods = {
            'LANDSAT/LC08/C01/T1_SR': 'from_landsat_c1_sr',
            'LANDSAT/LE07/C01/T1_SR': 'from_landsat_c1_sr',
            'LANDSAT/LT05/C01/T1_SR': 'from_landsat_c1_sr',
            # 'LANDSAT/LT04/C01/T1_SR': 'from_landsat_c1_sr',
            # 'LANDSAT/LC08/C01/T1_RT_TOA': 'from_landsat_c1_toa',
            # 'LANDSAT/LE07/C01/T1_RT_TOA': 'from_landsat_c1_toa',
            # 'LANDSAT/LC08/C01/T1_TOA': 'from_landsat_c1_toa',
            # 'LANDSAT/LE07/C01/T1_TOA': 'from_landsat_c1_toa',
            # 'LANDSAT/LT05/C01/T1_TOA': 'from_landsat_c1_toa',
            # # 'LANDSAT/LT04/C01/T1_TOA': 'from_landsat_c1_toa',
        }

        try:
            method_name = collection_methods[image_id.rsplit('/', 1)[0]]
        except KeyError:
            raise ValueError('unsupported collection ID: {}'.format(image_id))
        except Exception as e:
            raise Exception('unhandled exception: {}'.format(e))

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
        return cls(landsat.LandsatSR(sr_image).prep(), **kwargs)

    # @classmethod
    # def from_landsat_c1_toa(cls, toa_image, **kwargs):
    #     """Returns a DisALEXI Image instance from a Landsat Collection 1 TOA image
    #
    #     Parameters
    #     ----------
    #     toa_image : ee.Image
    #         A raw Landsat Collection 1 TOA image.
    #     kwargs : dict
    #         Keyword arguments to pass through to Image init function.
    #
    #     Returns
    #     -------
    #     Instance of Image class
    #
    #     """
    #     return cls(landsat.LandsatTOA(toa_image).prep(), **kwargs)

    # def calculate(self, variables=['et', 'etr', 'etf']):
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
                output_images.append(self.et)
            # elif v.lower() == 'etf':
            #     output_images.append(self.etf)
            # elif v.lower() == 'etr':
            #     output_images.append(self.etr)
            elif v.lower() == 'lst':
                output_images.append(self.lst)
            elif v.lower() == 'mask':
                output_images.append(self.mask)
            elif v.lower() == 'ndvi':
                output_images.append(self.ndvi)
            # elif v.lower() == 'qa':
            #     output_images.append(self.qa)
            # elif v.lower() == 'quality':
            #     output_images.append(self.quality)
            elif v.lower() == 'time':
                output_images.append(self.time)
            else:
                raise ValueError('unsupported variable: {}'.format(v))

        return ee.Image(output_images).set(self.properties)

    # @lazy_property
    # def albedo(self):
    #     """Return albedo image"""
    #     return self.image.select(['albedo']).set(self.properties).double()

    @lazy_property
    def et(self):
        """Compute Landsat scale DisALEXI ET

        Returns
        -------
        image : ee.Image
            DisALEXI ET image

        """
        et = tseb.tseb_pt(
            T_air=self.ta, T_rad=self.lst,
            u=self.windspeed, p=self.pressure, z=self.elevation,
            Rs_1=self.rs1, Rs24=self.rs24, vza=0,
            # CGM - Need to add GEE gaussian_filter call to Rs24
            # Rs24=ndimage.gaussian_filter(self.Rs24, sigma=5),
            aleafv=self.aleafv, aleafn=self.aleafn, aleafl=self.aleafl,
            adeadv=self.adeadv, adeadn=self.adeadn, adeadl=self.adeadl,
            albedo=self.albedo, ndvi=self.ndvi, lai=self.lai,
            clump=self.clump, leaf_width=self.leaf_width,
            hc_min=self.hc_min, hc_max=self.hc_max,
            datetime=self.date, stabil_iter=self.stabil_iter,
            albedo_iter=self.albedo_iter,
        )
        return et.rename(['et']).set(self.properties).double()

    # @lazy_property
    # def etr(self):
    #     """Compute reference ET for the image date"""
    #     return True

    # @lazy_property
    # def etf(self):
    #     """Compute ET fraction as actual ET divided by the reference ET"""
    #     return self.et.divide(self.etr)\
    #         .rename(['etf']).set(self.properties).double()

    @lazy_property
    def et_alexi(self):
        """Extract ALEXI ET image for the target image time"""
        if utils.is_number(self.alexi_source):
            # Interpret numbers as constant images
            alexi_img = ee.Image.constant(float(self.alexi_source))
        elif isinstance(self.alexi_source, ee.computedobject.ComputedObject):
            alexi_img = self.alexi_source
        elif self.alexi_source.upper() == 'CONUS_V001':
            alexi_coll_id = 'projects/disalexi/alexi/CONUS_V001'
            alexi_coll = ee.ImageCollection(alexi_coll_id) \
                .filterDate(self.date, self.date.advance(1, 'day'))
            alexi_img = ee.Image(alexi_coll.first())
            # self.alexi_geo = [0.04, 0, -125.04, 0, -0.04, 49.8]
            # self.alexi_crs = 'EPSG:4326'
        else:
            raise ValueError('unsupported alexi_source: {}'.format(
                self.alexi_source))
        return alexi_img.rename(['et_alexi'])

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
            raise ValueError('Unsupported elev_source: {}\n'.format(
                self.elevation_source))
        return elev_img.select([0], ['elevation'])

    # @lazy_property
    # def lai(self):
    #     """Return LAI image"""
    #     return self.image.select(['lai']).set(self.properties).double()

    # @lazy_property
    # def lst(self):
    #     """Return land surface temperature (LST) image"""
    #     return self.image.select(['lst']).set(self.properties).double()

    # CGM - Using CFMask for now
    # @lazy_property
    # def mask(self):
    #     """Mask of all active pixels (based on the final et)"""
    #     return self.et.multiply(0).add(1).updateMask(1)\
    #         .rename(['mask']).set(self.properties).uint8()

    # @lazy_property
    # def ndvi(self):
    #     """Return NDVI image"""
    #     return self.image.select(['ndvi']).set(self.properties).double()

    # @lazy_property
    # def quality(self):
    #     """Set quality to 1 for all active pixels (for now)"""
    #     return self.mask\
    #         .rename(['quality']).set(self.properties)

    @lazy_property
    def pressure(self):
        """Air pressure [kPa]"""
        return self.elevation \
            .expression(
                '101.3 * (((293.0 - 0.0065 * z) / 293.0) ** 5.26)',
                {'z': self.elevation}) \
            .rename(['pressure'])

    @lazy_property
    def rs1(self):
        """Extract MERRA2 solar images for the target image time"""
        # Interpolate rs hourly image at image time
        # Hourly Rs is time average so time starts are 30 minutes early
        # Move image time 30 minutes earlier to simplify filtering/interpolation
        # This interpolation scheme will only work for hourly data
        if utils.is_number(self.rs_hourly_source):
            rs1_img = ee.Image.constant(float(self.rs_hourly_source))
        elif isinstance(self.rs_hourly_source, ee.computedobject.ComputedObject):
            rs1_img = ee.Image(self.rs_hourly_source)
        elif self.rs_hourly_source.upper() == 'MERRA2':
            rs_coll = ee.ImageCollection('projects/disalexi/merra2/hourly')\
                .select(['SWGDNCLR'])

            if self.rs_interp_flag:
                interp_dt = self.date.advance(-0.5, 'hour')
                rs_a_img = ee.Image(rs_coll \
                    .filterDate(interp_dt.advance(-1, 'hour'), interp_dt).first())
                rs_b_img = ee.Image(rs_coll \
                    .filterDate(interp_dt, interp_dt.advance(1, 'hour')).first())
                t_a = ee.Number(rs_a_img.get('system:time_start'))
                t_b = ee.Number(rs_b_img.get('system:time_start'))
                rs1_img = rs_b_img.subtract(rs_a_img) \
                    .multiply(interp_dt.millis().subtract(t_a).divide(t_b.subtract(t_a))) \
                    .add(rs_a_img) \
                    .rename(['rs'])
            else:
                rs1_img = ee.Image(
                    ee.ImageCollection(rs_coll) \
                        .filterDate(self.date, self.date.advance(1, 'day'))\
                        .filter(ee.Filter.calendarRange(
                            self.hour_int, self.hour_int, 'hour'))
                        .first()) \
                    .rename(['rs'])
        else:
            raise ValueError('Unsupported rs_hourly_source: {}\n'.format(
                self.rs_hourly_source))

        return rs1_img.rename(['rs'])

    @lazy_property
    def rs24(self):
        if utils.is_number(self.rs_daily_source):
            rs24_img = ee.Image.constant(float(self.rs_daily_source))
        elif isinstance(self.rs_daily_source, ee.computedobject.ComputedObject):
            rs24_img = self.rs_daily_source
        elif self.rs_daily_source.upper() == 'MERRA2':
            rs_coll = ee.ImageCollection('projects/climate-engine/merra2/daily')\
                .select(['SWGDNCLR'])\
                .filterDate(self.date, self.date.advance(1, 'day'))
            rs24_img = ee.Image(rs_coll.first())
        else:
            raise ValueError('Unsupported rs_daily_source: {}\n'.format(
                self.rs_daily_source))
        return rs24_img.rename(['rs'])

    @lazy_property
    def ta(self):
        """Return the precomputed Ta asset for the target image

        Returns
        -------
        image : ee.Image
            Ta image

        """
        # if self.ta_source is None:
        #     raise ValueError('ta_source must be set to compute et')
        if utils.is_number(self.ta_source):
            ta_img = ee.Image.constant(float(self.ta_source))
        elif isinstance(self.ta_source, ee.computedobject.ComputedObject):
            ta_img = self.ta_source
        elif self.ta_source.upper() == 'CONUS_V001':
            ta_coll_id = 'projects/disalexi/ta/CONUS_V001'
            ta_coll = ee.ImageCollection(ta_coll_id)\
                .filterMetadata('id', 'equals', self.id)\
                .limit(1, 'iteration', False)
            input_img = ee.Image(ta_coll.first())

            # Compute new Ta as mean of bracketing values
            # ta_img = input_img.select(['ta_a', 'ta_b']).reduce(ee.Reducer.mean())
            # Interpolate new Ta based on bias values
            bias_a = input_img.select(['bias_a'])
            bias_b = input_img.select(['bias_b'])
            ta_a = input_img.select(['ta_a'])
            ta_b = input_img.select(['ta_b'])
            ta_img = ta_a.multiply(bias_b).subtract(ta_b.multiply(bias_a))\
                .divide(bias_b.subtract(bias_a))
            ta_img = ee.Image(ta_img.copyProperties(input_img))
        else:
            raise ValueError('Unsupported ta_source: {}\n'.format(
                self.ta_source))
        return ta_img.rename(['ta'])

    @lazy_property
    def time(self):
        """Return an image of the 0 UTC time (in milliseconds)"""
        return self.mask\
            .double().multiply(0).add(utils.date_to_time_0utc(self.date))\
            .rename(['time']).set(self.properties)
        # return ee.Image.constant(utils.date_to_time_0utc(self.date)) \
        #     .double().rename(['time']).set(self.properties)

    @lazy_property
    def windspeed(self):
        """Windspeed"""
        if utils.is_number(self.windspeed_source):
            windspeed_img = ee.Image.constant(float(self.windspeed_source))
        elif isinstance(self.windspeed_source, ee.computedobject.ComputedObject):
            windspeed_img = self.windspeed_source
        elif self.windspeed_source.upper() == 'CFSV2':
            # It would be more correct to compute the magnitude for each image,
            #   then compute the average.
            # Do we need daily, 6hr, or interpolated instantaneous data?
            windspeed_coll = ee.ImageCollection('NOAA/CFSV2/FOR6H') \
                .select([
                    'u-component_of_wind_height_above_ground',
                    'v-component_of_wind_height_above_ground']) \
                .filterDate(self.date, self.date.advance(1, 'day'))
            windspeed_img = windspeed_coll.mean() \
                .expression('sqrt(b(0) ** 2 + b(1) ** 2)')
        else:
            raise ValueError('Invalid windspeed_source: {}\n'.format(
                self.windspeed_source))
        return windspeed_img.rename(['windspeed'])

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
            raise KeyError('Invalid lc_type (choices are {})'.format(
                ', '.join(remaps.keys())))

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
            return ee.Image(landcover) \
                .remap(input_values, output_values) \
                .divide(100) \
                .rename([lc_var])

        # Get LC based variables
        self.aleafv = lc_remap(self.lc_source, self.lc_type, 'aleafv')
        self.aleafn = lc_remap(self.lc_source, self.lc_type, 'aleafn')
        self.aleafl = lc_remap(self.lc_source, self.lc_type, 'aleafl')
        self.adeadv = lc_remap(self.lc_source, self.lc_type, 'adeadv')
        self.adeadn = lc_remap(self.lc_source, self.lc_type, 'adeadn')
        self.adeadl = lc_remap(self.lc_source, self.lc_type, 'adeadl')
        self.hc_min = lc_remap(self.lc_source, self.lc_type, 'hmin')
        self.hc_max = lc_remap(self.lc_source, self.lc_type, 'hmax')
        self.leaf_width = lc_remap(self.lc_source, self.lc_type, 'xl')
        self.clump = lc_remap(self.lc_source, self.lc_type, 'omega')

    # Helper functions for compute Ta
    def ta_coarse(self, ta_img, cell_size=30, threshold=0.5):
        """Compute the air temperature for each ALEXI ET cell that minimizes
        the bias between Landsat scale ET and ALEXI ET for a fixed range of
        air temperature values

        Parameters
        ----------
        ta_img : ee.Image
        cell_size : float, optional
            Cellsize for "fine" scale calculation and aggregation.
        threshold : float, optional

        Returns
        -------
        image : ee.Image
            ALEXI scale air temperature image

        """
        et_fine = tseb.tseb_pt(
            T_air=ta_img, T_rad=self.lst,
            u=self.windspeed, p=self.pressure, z=self.elevation,
            Rs_1=self.rs1, Rs24=self.rs24, vza=0,
            aleafv=self.aleafv, aleafn=self.aleafn, aleafl=self.aleafl,
            adeadv=self.adeadv, adeadn=self.adeadn, adeadl=self.adeadl,
            albedo=self.albedo, ndvi=self.ndvi, lai=self.lai,
            clump=self.clump, leaf_width=self.leaf_width,
            hc_min=self.hc_min, hc_max=self.hc_max,
            datetime=self.date, a_PT_in=1.32,
            stabil_iter=self.stabil_iter, albedo_iter=self.albedo_iter,
        )

        # Aggregate the Landsat scale ET up to the ALEXI scale
        et_coarse = ee.Image(et_fine)\
            .reproject(crs=self.crs,
                       crsTransform=[cell_size, 0, 15, 0, -cell_size, 15])\
            .reduceResolution(reducer=ee.Reducer.mean(), maxPixels=20000)\
            .reproject(crs=self.alexi_crs, crsTransform=self.alexi_geo)
        et_mask = et_coarse.mask().gt(threshold)
        bias = et_coarse.subtract(self.et_alexi)
        # Invert the abs(bias) since quality mosaic sorts descending
        # qm_bias = bias.abs().multiply(-1)

        return ee.Image([ta_img, bias, et_coarse])\
            .updateMask(et_mask)\
            .rename(['ta', 'bias', 'et'])

    def ta_qm(self, ta_img, step_size, step_count, cell_size=30, threshold=0.5):
        """Compute the air temperature for each ALEXI ET cell that minimizes
        the bias between Landsat scale ET and ALEXI ET for a fixed range of
        air temperature values

        Parameters
        ----------
        ta_img : ee.Image
        step_size : float
            The size of each Ta step.
        step_count : int
            The number of + & - steps, so the total number of air temperatures
            images/values will be ta_steps * 2 + 1.
        cell_size : float, optional
            Cellsize for "fine" scale calculation and aggregation.
        threshold : float, optional

        Returns
        -------
        image : ee.Image
            ALEXI scale air temperature image

        """
        # CGM - This is basically identical to ta_coarse()
        # Can't use ta_coarse since a mapped function can only have one input
        # Redefining it here allows cell_size and threshold to be passed in
        def ta_func(ta):
            """Compute TSEB ET for the target T_air value"""
            et_fine = tseb.tseb_pt(
                T_air=ta, T_rad=self.lst,
                u=self.windspeed, p=self.pressure, z=self.elevation,
                Rs_1=self.rs1, Rs24=self.rs24, vza=0,
                aleafv=self.aleafv, aleafn=self.aleafn, aleafl=self.aleafl,
                adeadv=self.adeadv, adeadn=self.adeadn, adeadl=self.adeadl,
                albedo=self.albedo, ndvi=self.ndvi, lai=self.lai,
                clump=self.clump, leaf_width=self.leaf_width,
                hc_min=self.hc_min, hc_max=self.hc_max,
                datetime=self.date, a_PT_in=1.32,
                stabil_iter=self.stabil_iter, albedo_iter=self.albedo_iter,
            )

            # Aggregate the Landsat scale ET up to the ALEXI scale
            et_coarse = ee.Image(et_fine)\
                .reproject(crs=self.crs,
                           crsTransform=[cell_size, 0, 15, 0, -cell_size, 15])\
                .reduceResolution(reducer=ee.Reducer.mean(), maxPixels=20000)\
                .reproject(crs=self.alexi_crs, crsTransform=self.alexi_geo)
            et_mask = et_coarse.mask().gt(threshold)
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
            for i in range(-step_count, step_count + 1)])

        # Return the air temperature associated with the lowest difference
        #   in ET from the ALEXI ET value
        return ee.ImageCollection(ta_coll.map(ta_func))\
            .qualityMosaic('qm_bias')\
            .select(['ta', 'bias', 'et'])

    def et_coarse(self, ta_img, cell_size=30, threshold=0.5):
        """Compute the Landsat ET summed to the ALEXI grid for the
        specified air temperature image

        Parameters
        ----------
        ta_img : ee.Image
        cell_size : float, optional
            Cellsize for "fine" scale calculation and aggregation.
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
            T_air=ta_img, T_rad=self.lst,
            u=self.windspeed, p=self.pressure, z=self.elevation,
            Rs_1=self.rs1, Rs24=self.rs24, vza=0,
            # CGM - Need to add GEE gaussian_filter call to Rs24
            # Rs24=ndimage.gaussian_filter(self.Rs24, sigma=5),
            aleafv=self.aleafv, aleafn=self.aleafn, aleafl=self.aleafl,
            adeadv=self.adeadv, adeadn=self.adeadn, adeadl=self.adeadl,
            albedo=self.albedo, ndvi=self.ndvi, lai=self.lai,
            clump=self.clump, leaf_width=self.leaf_width,
            hc_min=self.hc_min, hc_max=self.hc_max,
            datetime=self.date, stabil_iter=self.stabil_iter,
            albedo_iter=self.albedo_iter,
        )
        et_coarse = et_fine\
            .reproject(crs=self.crs,
                       crsTransform=[cell_size, 0, 15, 0, -cell_size, 15])\
            .reduceResolution(reducer=ee.Reducer.mean(), maxPixels=20000)\
            .reproject(crs=self.alexi_crs, crsTransform=self.alexi_geo)\
            .rename(['et'])
        et_mask = et_coarse.mask().gt(threshold)
        return et_coarse.updateMask(et_mask)

    def et_bias(self, et_coarse_img):
        """Compute the bias between the computed ET and the ALEXI ET"""
        return et_coarse_img.subtract(self.et_alexi).rename(['bias'])

    # def smooth(self, T_air):
    #     """Resample image
    #
    #     Parameters
    #     ----------
    #     T_air
    #
    #     Returns
    #     -------
    #     image : ee.Image
    #
    #     """
    #     T_air = ee.Image(T_air) \
    #         .resample('bilinear') \
    #         .reproject(crs=self.alexi_crs, crsTransform=self.alexi_geo)
    #     return T_air
