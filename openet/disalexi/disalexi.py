# from builtins import input
# import datetime
import math

import ee
# import numpy as np

# Why can't these be imported directly
from .lc_properties import remaps
from . import tseb
from . import tseb_utils
from . import utils


class Image(object):
    """Earth Engine DisALEXI"""

    def __init__(
            self,
            image,
            elevation=None,
            landcover=None,
            lc_type=None,
            iterations=10,
            tair_values=list(range(273, 321, 1)),
        ):
        """Initialize an image for computing DisALEXI

        FilterDate looks at the time_starts, so if the Alexi image
            has a start time of 0 UTC, to get the Alexi image for the
            image date, you may need to move the image date back a day.

        Parameters
        ----------
        image : ee.Image
            Prepped image
        elevation: ee.Image
            Elevation asset.  Units must be in meters.
            Will default to ee.Image('USGS/NED') if not set or None.
        lc_type : str
            Land cover type.  Choices are 'NLCD' or 'GlobeLand30'.  Will
            default to 'NLCD' if not set and 'landcover' parameter not set.
        landcover : ee.Image
            Land cover asset.  Will default to ee.Image('USGS/NLCD/NLCD2011')
            if not set and 'lc_type' parameter not set.
        iterations : int
            Number of iterations of main calculation (the default is 10).
        """
        # self.image = ee.Image(image)
        self.iterations = iterations
        self.tair_values = tair_values

        # Set server side date/time properties using the 'system:time_start'
        self.datetime = ee.Date(ee.Image(image).get('system:time_start'))
        self.date = ee.Date(self.datetime.format('yyyy-MM-dd'))
        self.doy = ee.Number(self.datetime.getRelative('day', 'year')).add(1).double()
        self.hour = ee.Number(self.datetime.getFraction('day')).multiply(24)
        self.hour_int = self.hour.floor()
        # Time used in IDL is hours and fractional minutes (no seconds)
        self.time = ee.Date(self.datetime).get('hour').add(
            ee.Date(self.datetime).get('minute').divide(60))

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

        # Get input bands from the image
        self.albedo = input_image.select('albedo')
        self.lai = input_image.select('lai')
        # self.lai = self.lai.where(lai.mask(), 0.01)
        self.lst = input_image.select('lst')
        self.ndvi = input_image.select('ndvi')

        # Elevation [m]
        if elevation is None:
            self.elevation = ee.Image('USGS/NED').rename(['elevation'])
            # self.elevation = ee.Image('USGS/SRTMGL1_003').rename(['elevation'])
            # self.elevation = ee.Image('USGS/GTOPO30').rename(['elevation'])
            # self.elevation = ee.Image('USGS/GMTED2010').rename(['elevation'])
            # self.elevation = ee.Image('CGIAR/SRTM90_V4').rename(['elevation'])
            # self.elevation = ee.Image.constant(350.0).rename(['elevation'])
        else:
            self.elevation = elevation.rename(['elevation'])

        # Set default land cover image and type
        # For now default to CONUS and use default if image and type were not set
        # GlobeLand30 values need to be set to the lowest even multiple of 10,
        #   since that is currently what is in the landcover.xlsx file.
        # http://www.globallandcover.com/GLC30Download/index.aspx
        if landcover is None and lc_type is None:
            # Using NLCD as default land cover and type
            self.landcover = ee.Image('USGS/NLCD/NLCD2011').select(['landcover'])
            self.lc_type = 'NLCD'
            # Using GlobeLand30 land cover and type
            # lc_coll =  ee.ImageCollection('users/cgmorton/GlobeLand30')\
            #     .filterBounds(image.geometry().bounds(1))
            # self.landcover = ee.Image(lc_coll.mosaic()) \
            #     .divide(10).floor().multiply(10) \
            #     .rename(['landcover'])
            # self.lc_type = 'GLOBELAND30'
        elif landcover is None:
            # What should happen if on only the land cover image is set?
            raise ValueError('landcover must be set if lc_type is set')
        elif lc_type is None:
            # What should happen if on only the land cover type is set?
            # The images could be looked up by type from a default LC image dict
            raise ValueError('lc_type must be set if landcover is set')
        else:
            self.landcover = landcover
            self.lc_type = lc_type

        # ALEXI ET - CONUS
        self.et_coll = ee.ImageCollection('projects/disalexi/alexi/CONUS')
        self.et_transform = [0.04, 0, -125.04, 0, -0.04, 49.82]
        self.et_crs = 'EPSG:4326'
        # ALEXI ET - Global
        # self.et_coll = ee.ImageCollection(
        #     'projects/climate-engine/alexi/global/daily/et')
        # self.et_transform = [0.05, 0, -180.0, 0, -0.05, 90]
        # self.et_crs = 'EPSG:4326'

        # Hard coding using CFSR for wind speed
        self.windspeed_coll = ee.ImageCollection('NOAA/CFSV2/FOR6H') \
            .select([
                'u-component_of_wind_height_above_ground',
                'v-component_of_wind_height_above_ground'])

        # Hard coding using MERRA2 solar insolation [W m-2]
        self.rs_hourly_coll = ee.ImageCollection(
            'projects/climate-engine/merra2/rs_hourly')
        self.rs_daily_coll = ee.ImageCollection(
            'projects/climate-engine/merra2/rs_daily')

    def compute_ta_test(self):
        """"""

        self._set_alexi_et_vars()
        self._set_elevation_vars()
        self._set_landcover_vars()
        self._set_solar_vars()
        self._set_time_vars()
        self._set_weather_vars()

        def t_air_func(t_air):
            """Compute TSEB ET for each air temperature value"""

            # Start with ALEXI scale air temperature
            t_air_coarse_img = self.alexi_et.multiply(0) \
                .add(ee.Number(t_air.get('t_air')))

            # Does the Ta image need to be projected to the Landsat image?
            # .reproject(crs=self.image_crs, crsTransform=self.image_transform)

            et_img = tseb.tseb_pt(
                T_air=t_air_coarse_img,
                T_rad=self.lst,
                u=self.windspeed,
                p=self.pressure,
                z=self.elevation,
                Rs_1=self.rs1,
                Rs24=self.rs24,
                vza=0,
                zs=self.sol_zenith,
                aleafv=self.aleafv,
                aleafn=self.aleafn,
                aleafl=self.aleafl,
                adeadv=self.adeadv,
                adeadn=self.adeadn,
                adeadl=self.adeadl,
                albedo=self.albedo,
                ndvi=self.ndvi,
                lai=self.lai,
                clump=self.clump,
                hc=self.hc,
                time=self.time,
                t_rise=self.t_rise,
                t_end=self.t_end,
                leaf_width=self.leaf_width,
                a_PT_in=1.32,
                iterations=self.iterations)

            et_coarse_img = ee.Image(et_img) \
                .reduceResolution(reducer=ee.Reducer.mean(), maxPixels=30000) \
                .reproject(crs=self.et_crs, crsTransform=self.et_transform) \
                .updateMask(1)

            # Invert the abs(bias) since quality mosaic sorts descending
            bias = et_coarse_img.subtract(self.alexi_et).abs().multiply(-1)
            # bias = self.et_img.subtract(ET).abs().multiply(-1)
            return ee.Image([bias, et_coarse_img, t_air_coarse_img]) \
                .rename(['bias', 'et', 't_air'])

        t_air_values = ee.FeatureCollection([
            ee.Feature(None, {'t_air': t_air})
            for t_air in [x * 0.1 for x in range(2950, 3050, 1)]])
        # DEADBEEF - NumPy is only being used for this line
        #     for t_air in np.arange(295, 305, 0.1)])
        output_coll = ee.ImageCollection(t_air_values.map(t_air_func))

        # Return the air temperature associated with the lowest difference
        #   in ET from the ALEXI ET value
        t_air = ee.Image(output_coll.qualityMosaic('bias'))
        #     .select(['t_air'])

        # # Project the output back to the ALEXI pixels
        # #     .reproject(crs=self.et_crs, crsTransform=self.et_transform) \
        # #     .rename(['t_air'])

        return t_air

    def compute_ta(self):
        """Compute Landsat scale air temperature that minimizes bias between
        Landsat scale ET and ALEXI ET

        Parameters
        ----------

        Returns
        -------
        image : ee.Image
            Landsat scale air temperature image

        """
        self._set_alexi_et_vars()
        self._set_elevation_vars()
        self._set_landcover_vars()
        self._set_solar_vars()
        self._set_time_vars()
        self._set_weather_vars()

        def t_air_func(t_air):
            """Compute TSEB ET for each T_air value

            Assume the function is being mapped over a FC with a T_air property
            Return bias (ALEXI ET - TSEB ET) as first band for quality mosaic
            """

            # Apply T_air values to Landsat pixels
            t_air_img = self.lst.multiply(0).add(ee.Number(t_air.get('t_air')))
            # T_air_img = ee.Image.constant(ee.Number(t_air.get('t_air'))).double()
            # T_air_img = ee.Image.constant(ee.Number(t_air)).double()
            # return t_air_img

            et = tseb.tseb_pt(
                T_air=t_air_img,
                T_rad=self.lst,
                u=self.windspeed,
                p=self.pressure,
                z=self.elevation,
                Rs_1=self.rs1,
                Rs24=self.rs24,
                vza=0,
                zs=self.sol_zenith,
                aleafv=self.aleafv,
                aleafn=self.aleafn,
                aleafl=self.aleafl,
                adeadv=self.adeadv,
                adeadn=self.adeadn,
                adeadl=self.adeadl,
                albedo=self.albedo,
                ndvi=self.ndvi,
                lai=self.lai,
                clump=self.clump,
                hc=self.hc,
                time=self.time,
                t_rise=self.t_rise,
                t_end=self.t_end,
                leaf_width=self.leaf_width,
                a_PT_in=1.32,
                iterations=self.iterations)

            # Invert the abs(bias) since quality mosaic sorts descending
            bias = ee.Image(et).subtract(self.alexi_et).abs().multiply(-1)
            # bias = self.et_img.subtract(ET).abs().multiply(-1)
            return ee.Image([bias, ee.Image(et), t_air_img]) \
                .rename(['bias', 'et', 't_air'])

        # Get output for a single Tair value
        # output_img = t_air_func(ee.Feature(None, {'t_air': 285}))
        # # output_img = t_air_func(ee.Number(285))
        # print('Output: {}'.format(utils.image_value(output_img).values()))

        # Get output for a range of Tair values
        # Mapping over the list seemed a little slower than the FC
        t_air_values = ee.FeatureCollection([
            ee.Feature(None, {'t_air': t_air}) for t_air in self.tair_values])
        # T_air_values = ee.List.sequence(275, 335, 5)
        output_coll = ee.ImageCollection(t_air_values.map(t_air_func))

        # Return the air temperature associated with the lowest difference
        #   in ET from the ALEXI ET value
        t_air = ee.Image(output_coll.qualityMosaic('bias')) \
            .select(['t_air'])

        # Project the output back to the landsat scene
        #     .reproject(crs=, crsTransform=) \
        #     .rename(['t_air'])

        return t_air

    def compute_et(self, T_air):
        """Compute Landsat scale DisALEXI ET

        Parameters
        ----------
        T_air: ee.Image
            ALEXI ET scale air temperature (K)

        Returns
        -------
        image : ee.Image
            DisALEXI ET image

        """
        self._set_elevation_vars()
        self._set_alexi_et_vars()
        self._set_landcover_vars()
        self._set_solar_vars()
        self._set_time_vars()
        self._set_weather_vars()

        et = tseb.tseb_pt(
            T_air=T_air,
            T_rad=self.lst,
            u=self.windspeed,
            p=self.pressure,
            z=self.elevation,
            Rs_1=self.rs1,
            Rs24=self.rs24,
            # CGM - Need to add GEE gaussian_filter call to Rs24
            # Rs24=ndimage.gaussian_filter(self.Rs24, sigma=5),
            vza=0,
            zs=self.sol_zenith,
            aleafv=self.aleafv,
            aleafn=self.aleafn,
            aleafl=self.aleafl,
            adeadv=self.adeadv,
            adeadn=self.adeadn,
            adeadl=self.adeadl,
            albedo=self.albedo,
            ndvi=self.ndvi,
            lai=self.lai,
            clump=self.clump,
            hc=self.hc,
            time=self.time,
            t_rise=self.t_rise,
            t_end=self.t_end,
            leaf_width=self.leaf_width,
            a_PT_in=1.32,
            iterations=self.iterations)
        return ee.Image(et).rename(['et'])

    def aggregate(self, T_air):
        """EE approach for aggregating Landsat Ta to the ALEXI grid

        30000 pixels should works for aggregating Landsat to global ET (5 km)
        The ".updateMask(1)" call is used to remove the transparency mask
        Should we use the ".mean().unweighted()" reducer?

        Parameters
        ----------
        T_air: ee.Image
            Landsat scale air temperature (K)

        Returns
        -------
        image : ee.Image

        """
        T_air = ee.Image(T_air) \
            .reduceResolution(reducer=ee.Reducer.mean(), maxPixels=30000) \
            .reproject(crs=self.et_crs, crsTransform=self.et_transform) \
            .updateMask(1)
        return T_air

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
    #         .reproject(crs=self.et_crs, crsTransform=self.et_transform)
    #     return T_air

    def _set_alexi_et_vars(self):
        """Extract ALEXI ET image for the target image time"""
        self.alexi_et = ee.Image(ee.ImageCollection(self.et_coll) \
            .filterDate(self.date, self.date.advance(1, 'day')).first())
        self.alexi_et = self.alexi_et.rename(['alexi_et'])

    def _set_elevation_vars(self):
        """Compute elevation derived variables"""
        self.pressure = self.elevation \
            .expression(
                '101.3 * (((293.0 - 0.0065 * z) / 293.0) ** 5.26)',
                {'z': self.elevation}) \
            .rename(['pressure'])

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
        self.aleafv = lc_remap(self.landcover, self.lc_type, 'aleafv')
        self.aleafn = lc_remap(self.landcover, self.lc_type, 'aleafn')
        self.aleafl = lc_remap(self.landcover, self.lc_type, 'aleafl')
        self.adeadv = lc_remap(self.landcover, self.lc_type, 'adeadv')
        self.adeadn = lc_remap(self.landcover, self.lc_type, 'adeadn')
        self.adeadl = lc_remap(self.landcover, self.lc_type, 'adeadl')
        hc_min = lc_remap(self.landcover, self.lc_type, 'hmin')
        hc_max = lc_remap(self.landcover, self.lc_type, 'hmax')
        self.leaf_width = lc_remap(self.landcover, self.lc_type, 'xl')
        self.clump = lc_remap(self.landcover, self.lc_type, 'omega')

        # LAI for leafs spherical distribution
        F = self.lai.multiply(self.clump).rename(['F'])

        # Fraction cover at nadir (view=0)
        f_c = self.lai.expression('1.0 - exp(-0.5 * F)', {'F': F}) \
            .clamp(0.01, 0.9) \
            .rename(['f_c'])

        # ======================================================================
        # Compute canopy height and roughness parameters
        self.hc = self.lai \
            .expression(
                'hc_min + ((hc_max - hc_min) * f_c)',
                {'hc_min': hc_min, 'hc_max': hc_max, 'f_c': f_c}) \
            .rename(['hc'])

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
                    .filter(ee.Filter.calendarRange(self.hour_int, self.hour_int, 'hour'))
                    .first()) \
                .rename(['rs'])
        self.rs24 = ee.Image(
            ee.ImageCollection(self.rs_daily_coll) \
                .filterDate(self.date, self.date.advance(1, 'day')) \
                .first()) \
            .rename(['rs'])

    def _set_time_vars(self):
        """Compute time and position related variables

        CGM - The zs returned by this function is not used in the original
            Python code.
        The hour in the original call was the integer hour, even though it
            seems like it should be the float hour.

        """
        self.t_rise, self.t_end = tseb_utils.sunrise_sunset(
            date=self.datetime,
            lon=ee.Image.pixelLonLat().select(['longitude']).multiply(math.pi / 180),
            lat=ee.Image.pixelLonLat().select(['latitude']).multiply(math.pi / 180))
        self.sol_zenith = tseb_utils.solar_zenith(
            date=self.datetime,
            lon=ee.Image.pixelLonLat().select(['longitude']).multiply(math.pi / 180),
            lat=ee.Image.pixelLonLat().select(['latitude']).multiply(math.pi / 180))

    def _set_weather_vars(self):
        """Compute weather derived variables (only wind from CFSv2 for now)

        Assume input image is a CFSv2 collection
        Assume wind speed image has two bands and compute magnitude
        It would probably make more sense to compute a single windspeed image
            in the input collection.
        For simplicity, computing the mean of all images in the UTC day

        Do we need daily, 6hr, or interpolated instantaneous data?

        CFSv2 units q: kg/kg, p: Pa, ta: K, wind: m/s
        """
        windspeed_img = ee.Image(self.windspeed_coll \
            .filterDate(self.date, self.date.advance(1, 'day')).mean())
        self.windspeed = windspeed_img \
            .expression('sqrt(b(0) ** 2 + b(1) ** 2)') \
            .rename(['windspeed'])

        # CGM - Per Mitch, q2, Ta, pressure, and ea are not used in the current
        #   implementation of DisALEXI.
        # self.q2 = cfsr_img.select(['Specific_humidity_height_above_ground']) \
        #     .rename(['q2'])
        # self.q2 = cfsr_img.expression(
        #     '0.5 * (a + b)',
        #     {'a': cfsr_img.select(['Minimum_specific_humidity_at_2m_height_above_ground_6_Hour_Interval'])
        #      'b': cfsr_img.select(['Maximum_specific_humidity_at_2m_height_above_ground_6_Hour_Interval'])})

        # self.ta = cfsr_img.select(['Temperature_height_above_ground']) \
        #     .rename(['ta'])
        # self.ta = cfsr_img.expression(
        #     '0.5 * (a + b)',
        #     {'a': cfsr_img.select(['Minimum_temperature_height_above_ground_6_Hour_Interval'])
        #      'b': cfsr_img.select(['Maximum_temperature_height_above_ground_6_Hour_Interval'])})

        # Get modeled air pressure (convert Pa to kPa)
        # self.pressure = cfsr_img.select(['Pressure_surface']) \
        #     .multiply(0.001).rename(['pressure'])

