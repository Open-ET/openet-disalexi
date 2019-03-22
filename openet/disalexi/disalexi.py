# from builtins import input
# import datetime
import pprint

import ee
# import numpy as np

# Why can't these be imported directly
from .lc_properties import remaps
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
            # ta_values=None,
            # ta_source=None,
            alexi_et_source='CONUS_V001',
            elevation_source='USGS/SRTMGL1_003',
            landcover_source='NLCD2011',
            rs_daily_source='MERRA2',
            rs_hourly_source='MERRA2',
            windspeed_source='CFSV2',
            stabil_iterations=10,
            albedo_iterations=4,
        ):
        """Initialize an image for computing DisALEXI

        Parameters
        ----------
        image : ee.Image
            Prepped image
        ta_values : list
            Air temperature values to test(the default is None).
            Must be set to compute "ta".
        ta_source : str, float
            ALEXI scale air temperature image collection ID (the default is None).
            Must be set to compute "et".
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

        Notes
        -----
        For now defaulting all inputs to the CONUS based inputs.

        FilterDate looks at the time_starts, so if the ALEXI image
            has a start time of 0 UTC, to get the ALEXI image for the
            image date, you may need to move the image date back a day.

        """
        # self.image = ee.Image(image)

        # Set server side date/time properties using the 'system:time_start'
        self.datetime = ee.Date(ee.Image(image).get('system:time_start'))
        self.date = ee.Date(self.datetime.format('yyyy-MM-dd'))
        self.doy = ee.Number(self.datetime.getRelative('day', 'year')).add(1).double()
        self.hour = ee.Number(self.datetime.getFraction('day')).multiply(24)
        self.hour_int = self.hour.floor()
        # Time used in IDL is hours and fractional minutes (no seconds)
        # self.time = ee.Date(self.datetime).get('hour').add(
        #     ee.Date(self.datetime).get('minute').divide(60))

        # Client side date/time properties can't be used in mapped functions
        # self.dt_client = datetime.datetime.utcfromtimestamp(
        #     self.datetime.millis().getInfo() / 1000)
        # self.date_client = self.datetime.strftime('%Y-%m-%d')
        # self.doy_client = int(self.datetime.strftime('%j'))
        # self.hour_client = self.datetime.hour

        # CGM - Applying cloud mask directly to input image
        #   instead of to a_pt in main TSEB function
        self.cfmask = image.select('cfmask')
        self.mask = self.cfmask.eq(0)
        input_image = ee.Image(image).updateMask(self.mask)

        # input_image = ee.ImageCollection(image).map()

        # Get input bands from the image
        self.albedo = input_image.select('albedo')
        self.lai = input_image.select('lai')
        # self.lai = self.lai.where(lai.mask(), 0.01)
        self.lst = input_image.select('lst')
        self.ndvi = input_image.select('ndvi')

        # Set input parameters
        # self.ta_values = ta_values
        # self.ta_source = ta_source
        self.alexi_et_source = alexi_et_source
        self.landcover_source = landcover_source
        self.elevation_source = elevation_source
        self.windspeed_source = windspeed_source
        self.rs_daily_source = rs_daily_source
        self.rs_hourly_source = rs_hourly_source
        self.stabil_iter = int(stabil_iterations)
        self.albedo_iter = int(albedo_iterations)

        if self.alexi_et_source.upper() == 'CONUS_V001':
            alexi_id = 'projects/disalexi/alexi/CONUS_V001'
            alexi_coll = ee.ImageCollection(alexi_id) \
                .filterDate(self.date, self.date.advance(1, 'day'))
            self.alexi_et = ee.Image(alexi_coll.first()).rename(['alexi_et'])
            self.alexi_mask = self.alexi_et.multiply(0).rename(['alexi_mask'])
            self.alexi_geo = [0.04, 0, -125.04, 0, -0.04, 49.8]
            self.alexi_crs = 'EPSG:4326'
        else:
            raise ValueError('Unsupported alexi_et_source: {}\n'.format(
                self.alexi_et_source))

        # Set default land cover image and type
        # For now default to CONUS and use default if image and type were not set
        # GlobeLand30 values need to be set to the lowest even multiple of 10,
        #   since that is currently what is in the landcover.xlsx file.
        # http://www.globallandcover.com/GLC30Download/index.aspx
        if self.landcover_source.upper() == 'NLCD2011':
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

        if self.rs_daily_source.upper() == 'MERRA2':
            self.rs_daily_coll = ee.ImageCollection(
                    'projects/climate-engine/merra2/daily')\
                .select(['SWGDNCLR'])
        #         .filterDate(self.date, self.date.advance(1, 'day')) \
        # #       'projects/climate-engine/merra2/daily').select(['SWGDN'])
        # elif self.rs_daily_source.upper() == 'MERRA2':
        #     self.rs_daily_coll = ee.ImageCollection(
        #         'projects/climate-engine/merra2/rs_daily')
        else:
            raise ValueError('unsupported rs_daily_source: {}'.format(
                self.rs_daily_source))

        if self.rs_hourly_source.upper() == 'MERRA2':
            self.rs_hourly_coll = ee.ImageCollection(
                    'projects/disalexi/merra2/hourly')\
                .select(['SWGDNCLR'])
        #         'projects/climate-engine/merra2/hourly').select(['SWGDN'])
        #         .filterDate(self.date, self.date.advance(1, 'day')) \
        #     self.rs_hourly_coll = ee.ImageCollection(
        #         'projects/climate-engine/merra2/rs_hourly')
        # elif self.rs_hourly_source.upper() == 'CFSV2':
        #     self.rs_hourly_coll = ee.ImageCollection('NOAA/CFSV2/FOR6H')
        else:
            raise ValueError('unsupported rs_hourly_source: {}'.format(
                self.rs_hourly_source))

        # Image projection and geotransform
        self.crs = image.projection().crs()
        self.transform = ee.List(ee.Dictionary(
            ee.Algorithms.Describe(image.projection())).get('transform'))
        # self.crs = image.select([0]).projection().getInfo()['crs']
        # self.transform = image.select([0]).projection().getInfo()['transform']

        self._set_landcover_vars()
        self._set_solar_vars()

    def ta_qm(self, cellsize=120, ta_values=list(range(250, 351, 1)),
              threshold=0.5):
        """Compute the air temperature for each ALEXI ET cell that minimizes
        the bias between Landsat scale ET and ALEXI ET for a fixed range of
        air temperature values

        Parameters
        ----------
        cellsize : float, optional
            Cellsize for "fine" scale calculation and aggregation.
        ta_values : list, optional
            List of air temperature values that will be used in calculation.
        threshold : float, optional

        Returns
        -------
        image : ee.Image
            ALEXI scale air temperature image

        """
        def ta_func(ta_ftr):
            """Compute TSEB ET for the target T_air value"""
            ta = ee.Number(ta_ftr.get('ta'))
            ta_fine = self.lst.multiply(0).add(ta).rename(['ta'])
            ta_coarse = self.alexi_et.multiply(0).add(ta).rename(['ta'])

            et_fine = tseb.tseb_pt(
                T_air=ta_fine, T_rad=self.lst,
                u=self.windspeed, p=self.pressure, z=self.elevation,
                Rs_1=self.rs1, Rs24=self.rs24, vza=0,
                aleafv=self.aleafv, aleafn=self.aleafn, aleafl=self.aleafl,
                adeadv=self.adeadv, adeadn=self.adeadn, adeadl=self.adeadl,
                albedo=self.albedo, ndvi=self.ndvi, lai=self.lai,
                clump=self.clump, leaf_width=self.leaf_width,
                hc_min=self.hc_min, hc_max=self.hc_max,
                datetime=self.datetime, a_PT_in=1.32,
                stabil_iter=self.stabil_iter, albedo_iter=self.albedo_iter,
            )

            # Aggregate the Landsat scale ET up to the ALEXI scale
            et_coarse = ee.Image(et_fine)\
                .reproject(crs=self.crs,
                           crsTransform=[cellsize, 0, 15, 0, -cellsize, 15])\
                .reduceResolution(reducer=ee.Reducer.mean(), maxPixels=20000)\
                .reproject(crs=self.alexi_crs, crsTransform=self.alexi_geo)
            et_mask = et_coarse.mask().gt(threshold)
            bias = et_coarse.subtract(self.alexi_et)
            # Invert the abs(bias) since quality mosaic sorts descending
            qm_bias = bias.abs().multiply(-1)

            return ee.Image([qm_bias, ta_coarse, bias, et_coarse])\
                .updateMask(et_mask)\
                .rename(['qm_bias', 'ta', 'bias', 'et'])

        # Get output for a range of Ta values
        # Mapping over the list seemed a little slower than the FC
        ta_coll = ee.FeatureCollection([
            ee.Feature(None, {'ta': ta}) for ta in ta_values])
        # ta_coll = ee.List(ta_values)
        # ta_coll = ee.List.sequence(250, 351, 1)

        # Return the air temperature associated with the lowest difference
        #   in ET from the ALEXI ET value
        return ee.ImageCollection(ta_coll.map(ta_func))\
            .qualityMosaic('qm_bias') \
            .select(['ta', 'bias', 'et'])

    def ta_iter(self, cellsize=1920, iterations=2, method='golden',
                ta_init=[250, 350]):
        """Compute the air temperature for each ALEXI ET cell that minimizes
        the bias between Landsat scale ET and ALEXI ET over a range of air
        temperature values.

        Parameters
        ----------
        cellsize : float, optional
            Cellsize for "fine" scale calculation and aggregation.
        iterations : int, optional
            Number of iterations (the default is 2).
        method : str, optional
            Iteration/optimization method (the default is 'golden').
        ta_init : list, optional
            The initial Ta air lower and upper bound values.

        Returns
        -------
        image : ee.Image
            ALEXI scale air temperature image

        """

        def et_coarse(ta_img):
            """Compute the Landsat ET summed to the ALEXI grid for the
            specified air temperature image"""

            # ta_fine = self.lst.multiply(0).add(ta).rename(['ta'])
            # ta_coarse = self.alexi_et.multiply(0).add(ta).rename(['ta'])

            et_fine = tseb.tseb_pt(
                T_air=ta_img, T_rad=self.lst,
                u=self.windspeed, p=self.pressure, z=self.elevation,
                Rs_1=self.rs1, Rs24=self.rs24, vza=0,
                aleafv=self.aleafv, aleafn=self.aleafn, aleafl=self.aleafl,
                adeadv=self.adeadv, adeadn=self.adeadn, adeadl=self.adeadl,
                albedo=self.albedo, ndvi=self.ndvi, lai=self.lai,
                clump=self.clump, leaf_width=self.leaf_width,
                hc_min=self.hc_min, hc_max=self.hc_max,
                datetime=self.datetime, stabil_iter=self.stabil_iter,
                albedo_iter=self.albedo_iter,
            )

            # Aggregate the Landsat scale ET up to the ALEXI scale
            return ee.Image(et_fine)\
                .reproject(crs=self.crs,
                           crsTransform=[cellsize, 0, 15, 0, -cellsize, 15])\
                .reduceResolution(reducer=ee.Reducer.mean(), maxPixels=20000)\
                .reproject(crs=self.alexi_crs, crsTransform=self.alexi_geo)\
                .rename(['et'])

        def et_bias(et_coarse_img):
            return et_coarse_img.subtract(self.alexi_et).rename(['bias'])

        def bisect_func(n, prev):
            # A is the lower bound, B is the upper bound
            # The minimum should lie within this bound
            ta_a = ee.Image(ee.Dictionary(prev).get('ta_a'))
            ta_b = ee.Image(ee.Dictionary(prev).get('ta_b'))
            bias_a = ee.Image(ee.Dictionary(prev).get('bias_a'))
            bias_b = ee.Image(ee.Dictionary(prev).get('bias_b'))

            # Generate new guess
            ta_c = ta_a.add(ta_b).multiply(0.5)
            # et_c = et_coarse(ta_c)
            bias_c = et_bias(et_coarse(ta_c))

            # Only use the new value if it causes [A, B] to bracket 0
            mask1 = bias_c.lte(bias_b).And(bias_c.gte(0))
            mask2 = bias_c.gte(bias_a).And(bias_c.lt(0))

            ta_b = ta_b.where(mask1, ta_c)
            bias_b = bias_b.where(mask1, bias_c)
            ta_a = ta_a.where(mask2, ta_c)
            bias_a = bias_a.where(mask2, bias_c)

            return ee.Dictionary({
                'ta_a': ta_a, 'bias_a': bias_a, 'ta_b': ta_b, 'bias_b': bias_b
            })

        # def false_pos_func(n, prev):
        #     # A is the lower bound, B is the upper bound
        #     # The minimum should lie within this bound
        #     ta_a = ee.Image(ee.Dictionary(prev).get('ta_a'))
        #     ta_b = ee.Image(ee.Dictionary(prev).get('ta_b'))
        #     bias_a = ee.Image(ee.Dictionary(prev).get('bias_a'))
        #     bias_b = ee.Image(ee.Dictionary(prev).get('bias_b'))
        #
        #     # Generate new guess
        #     ta_c = ta_a.multiply(bias_b).subtract(ta_b.multiply(bias_a)) \
        #         .divide(bias_b.subtract(bias_a))
        #     # et_c = et_coarse(ta_c)
        #     bias_c = et_bias(et_coarse(ta_c))
        #
        #     # Only use the new value if it causes [A, B] to bracket 0
        #     mask1 = bias_c.lte(bias_b).And(bias_c.gte(0))
        #     mask2 = bias_c.gte(bias_a).And(bias_c.lt(0))
        #
        #     ta_b = ta_b.where(mask1, ta_c)
        #     bias_b = bias_b.where(mask1, bias_c)
        #     ta_a = ta_a.where(mask2, ta_c)
        #     bias_a = bias_a.where(mask2, bias_c)
        #
        #     return ee.Dictionary({
        #         'ta_a': ta_a, 'bias_a': bias_a, 'ta_b': ta_b, 'bias_b': bias_b
        #     })

        def golden_func(n, prev):
            """Golden section search function"""
            # A is the lower bound, B is the upper bound
            # The minimum should lie within this bound
            ta_a = ee.Image(ee.Dictionary(prev).get('ta_a'))
            ta_b = ee.Image(ee.Dictionary(prev).get('ta_b'))
            ta_c = ee.Image(ee.Dictionary(prev).get('ta_c'))
            ta_d = ee.Image(ee.Dictionary(prev).get('ta_d'))
            bias_a = ee.Image(ee.Dictionary(prev).get('bias_a'))
            bias_b = ee.Image(ee.Dictionary(prev).get('bias_b'))
            bias_c = ee.Image(ee.Dictionary(prev).get('bias_c'))
            bias_d = ee.Image(ee.Dictionary(prev).get('bias_d'))

            # If f(c) < f(d): move the data from d to b and c to d
            mask1 = bias_c.lt(bias_d)
            # If f(c) > f(d): move the data from c to a and d to c
            mask2 = bias_c.gte(bias_d)

            ta_b = ta_b.where(mask1, ta_d)
            bias_b = bias_b.where(mask1, bias_d)
            # CGM - This isn't working but (re)computing below is?
            #ta_d = ta_d.where(mask1, ta_c);
            #bias_d = bias_d.where(mask1, bias_c);

            ta_a = ta_a.where(mask2, ta_c)
            bias_a = bias_a.where(mask2, bias_c)
            # ta_c = ta_d.where(mask2, ta_d);
            # bias_c = bias_d.where(mask2, bias_d)

            # Compute new test points
            ta_c = ta_b.subtract(ta_b.subtract(ta_a).multiply(0.618034))
            ta_d = ta_a.add(ta_b.subtract(ta_a).multiply(0.618034))

            # Compute the bias at the new test points
            bias_c = et_bias(et_coarse(ta_c))
            bias_d = et_bias(et_coarse(ta_d))

            return ee.Dictionary({
                'ta_a': ta_a, 'bias_a': bias_a,
                'ta_b': ta_b, 'bias_b': bias_b,
                'ta_c': ta_c, 'bias_c': bias_c,
                'ta_d': ta_d, 'bias_d': bias_d,
            })

        if method == 'bisect':
            input_images = ee.Dictionary({
                'ta_a': self.alexi_mask.add(ta_init[0]),
                'ta_b': self.alexi_mask.add(ta_init[1]),
                'bias_a': et_bias(et_coarse(self.alexi_mask.add(ta_init[0]))),
                'bias_b': et_bias(et_coarse(self.alexi_mask.add(ta_init[1]))),
            })
            output = ee.Dictionary(ee.List.sequence(1, iterations).iterate(
                bisect_func, input_images))

            # Use the average of the final bracketing values
            ta = ee.Image(output.get('ta_a')).add(
                ee.Image(output.get('ta_b'))).multiply(0.5)

        elif method == 'golden':
            gr = 0.618034 * (ta_init[1] - ta_init[0])
            input_images = ee.Dictionary({
                'ta_a': self.alexi_mask.add(ta_init[0]),
                'ta_b': self.alexi_mask.add(ta_init[1]),
                'ta_c': self.alexi_mask.add(ta_init[1] - gr),
                'ta_d': self.alexi_mask.add(ta_init[0] + gr),
                'bias_a': et_bias(et_coarse(self.alexi_mask.add(ta_init[0]))),
                'bias_b': et_bias(et_coarse(self.alexi_mask.add(ta_init[1]))),
                'bias_c': et_bias(et_coarse(self.alexi_mask.add(ta_init[1] - gr))),
                'bias_d': et_bias(et_coarse(self.alexi_mask.add(ta_init[0] + gr))),
            })

            output = ee.Dictionary(ee.List.sequence(1, iterations).iterate(
                golden_func, input_images))

            # Use the average of the final bracketing values
            ta = ee.Image(output.get('ta_c')).add(
                ee.Image(output.get('ta_d'))).multiply(0.5)

        else:
            raise ValueError('unsupported Ta iteration method')

        et_img = et_coarse(ta)
        return ee.Image([ta, et_bias(et_img), et_img])\
            .rename(['ta', 'bias', 'et'])

    def ta_mask(self, threshold=None, cellsize=30, ta_value=290):
        """Identify ALEXI scale pixels with sufficient cloud free non-water
        Landsat pixels to compute ALEXI scale air temperature.

        Parameters
        ----------
        threshold : float, optional
            Fraction of Landsat pixels that must be clear (the default is None).
            If a threshold is set, cells below the threshold will be masked out.
            If not set (or set to None) the mask percentage will be returned.
        cellsize : float, optional
            Cellsize for "fine" scale calculation and aggregation
            (the default is 120).
        ta_value : list, optional
            Air temperature value to use for the ET calculation
            (the default is 290).

        Returns
        -------
        image : ee.Image

        """
        et_fine = tseb.tseb_pt(
            T_air=ee.Image.constant(ta_value), T_rad=self.lst,
            u=self.windspeed, p=self.pressure, z=self.elevation,
            Rs_1=self.rs1, Rs24=self.rs24, vza=0,
            aleafv=self.aleafv, aleafn=self.aleafn, aleafl=self.aleafl,
            adeadv=self.adeadv, adeadn=self.adeadn, adeadl=self.adeadl,
            albedo=self.albedo, ndvi=self.ndvi, lai=self.lai,
            clump=self.clump, leaf_width=self.leaf_width,
            hc_min=self.hc_min, hc_max=self.hc_max,
            datetime=self.datetime, stabil_iter=self.stabil_iter,
            albedo_iter=self.albedo_iter,
        )
        mask_img = ee.Image(et_fine) \
            .reproject(crs=self.crs,
                       crsTransform=[cellsize, 0, 15, 0, -cellsize, 15]) \
            .reduceResolution(reducer=ee.Reducer.count(), maxPixels=20000) \
            .reproject(crs=self.alexi_crs, crsTransform=self.alexi_geo) \
            .mask()

        if threshold is not None:
            mask_img = mask_img.gt(threshold)
            mask_img = mask_img.updateMask(mask_img).updateMask(1)

        return mask_img.rename(['mask'])

        # DEADBEEF
        # return mask_img.updateMask(mask_img).updateMask(1).rename(['mask'])

    @lazy_property
    def et(self, ta_img):
        """Compute Landsat scale DisALEXI ET

        Returns
        -------
        image : ee.Image
            DisALEXI ET image

        """
        # if self.ta_source is None:
        #     raise ValueError('ta_source must be set to compute et')
        # elif utils.is_number(self.ta_source):
        #     ta_img = ee.Image.constant(float(self.ta_source))
        # elif type(self.ta_source) is str:
        #     ta_img = ee.Image(ee.ImageCollection(self.ta_source) \
        #         .filterDate(self.date, self.date.advance(1, 'day')).first())

        return tseb.tseb_pt(
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
            datetime=self.datetime, stabil_iter=self.stabil_iter,
            albedo_iter=self.albedo_iter,
        )

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

    # @lazy_property
    # def alexi_et(self):
    #     """Extract ALEXI ET image for the target image time"""
    #     if self.alexi_et_source.upper() == 'CONUS_V001':
    #         alexi_et_id = 'projects/disalexi/alexi/CONUS_V001'
    #         alexi_et_coll = ee.ImageCollection(alexi_et_id) \
    #             .filterDate(self.date, self.date.advance(1, 'day'))
    #         alexi_et_img = ee.Image(alexi_et_coll.first())
    #         self.alexi_geo = [0.04, 0, -125.04, 0, -0.04, 49.8]
    #         self.alexi_crs = 'EPSG:4326'
    #     else:
    #         raise ValueError('Unsupported alexi_et_source: {}\n'.format(
    #             self.alexi_et_source))
    #
    #     return alexi_et_img.rename(['alexi_et'])

    @lazy_property
    def elevation(self):
        """Elevation"""
        if utils.is_number(self.elevation_source):
            elev_img = ee.Image.constant(float(self.elevation_source))
        # elif self.elevation_source.upper() == 'GTOPO':
        #     elev_image = ee.Image('USGS/GTOPO30')
        # elif self.elevation_source.upper() == 'NED':
        #     elev_image = ee.Image('USGS/NED')
        # elif self.elevation_source.upper() == 'SRTM':
        #     elev_image = ee.Image('USGS/SRTMGL1_003')
        # elif (self.elevation_source.lower().startswith('projects/') or
        #       self.elevation_source.lower().startswith('users/')):
        #     elev_image = ee.Image(self.elevation_source)
        elif type(self.elevation_source) is str:
            elev_img = ee.Image(self.elevation_source)
        else:
            raise ValueError('Unsupported elev_source: {}\n'.format(
                self.elevation_source))

        return elev_img.select([0], ['elevation'])

    @lazy_property
    def pressure(self):
        """Air pressure"""
        return self.elevation \
            .expression(
                '101.3 * (((293.0 - 0.0065 * z) / 293.0) ** 5.26)',
                {'z': self.elevation}) \
            .rename(['pressure'])

    @lazy_property
    def windspeed(self):
        """Windspeed"""

        # if utils.is_number(self.windspeed_source):
        #     dt_img = ee.Image.constant(float(self._dt_source))
        if self.windspeed_source.upper() == 'CFSV2':
            # It would be more correct to compute the magnitude for each image,
            #   then compute the average.
            # Do we need daily, 6hr, or interpolated instantaneous data?
            windspeed_coll = ee.ImageCollection('NOAA/CFSV2/FOR6H') \
                .select([
                    'u-component_of_wind_height_above_ground',
                    'v-component_of_wind_height_above_ground']) \
                .filterDate(self.date, self.date.advance(1, 'day'))
            windspeed_img = windspeed_coll.mean() \
                .expression('sqrt(b(0) ** 2 + b(1) ** 2)') \
                .rename(['windspeed'])

        else:
            raise ValueError('Invalid windspeed_source: {}\n'.format(
                self.windspeed_source))

        return windspeed_img

    def _set_landcover_vars(self):
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

    def _set_solar_vars(self, interpolate_flag=True):
        """Extract MERRA2 solar images for the target image time"""
        # Interpolate rs hourly image at image time
        # Hourly Rs is time average so time starts are 30 minutes early
        # Move image time 30 minutes earlier to simplify filtering/interpolation
        # This interpolation scheme will only work for hourly data
        if interpolate_flag:
            interp_dt = self.datetime.advance(-0.5, 'hour')
            # time_a = interp_time
            # time_b = interp_time
            rs_a_img = ee.Image(self.rs_hourly_coll \
                .filterDate(interp_dt.advance(-1, 'hour'), interp_dt).first())
            rs_b_img = ee.Image(self.rs_hourly_coll \
                .filterDate(interp_dt, interp_dt.advance(1, 'hour')).first())
            t_a = ee.Number(rs_a_img.get('system:time_start'))
            t_b = ee.Number(rs_b_img.get('system:time_start'))
            self.rs1 = rs_b_img.subtract(rs_a_img) \
                .multiply(interp_dt.millis().subtract(t_a).divide(t_b.subtract(t_a))) \
                .add(rs_a_img) \
                .rename(['rs'])
        else:
            self.rs1 = ee.Image(
                ee.ImageCollection(self.rs_hourly_coll) \
                    .filterDate(self.date, self.date.advance(1, 'day')) \
                    .filter(ee.Filter.calendarRange(
                        self.hour_int, self.hour_int, 'hour'))
                    .first()) \
                .rename(['rs'])
        self.rs24 = ee.Image(
            ee.ImageCollection(self.rs_daily_coll) \
                .filterDate(self.date, self.date.advance(1, 'day')) \
                .first()) \
            .rename(['rs'])
