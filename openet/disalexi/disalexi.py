from importlib import metadata
import pprint
import re
import warnings

import ee
import openet.lai

# Why can't these be imported directly
from .lc_properties import lc_remaps
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
            lai_source='openet-landsat-lai',
            # lai_source='projects/openet/assets/lai/landsat/c02',
            lst_source='projects/openet/assets/lst/landsat/c02',
            elevation_source='USGS/SRTMGL1_003',
            # elevation_source='NASA/NASADEM_HGT/001',
            landcover_source='USGS/NLCD_RELEASES/2021_REL/NLCD',
            air_pres_source='CFSR',
            air_temp_source='CFSR',
            rs_daily_source='CFSR',
            rs_hourly_source='CFSR',
            vapor_pres_source='CFSR',
            wind_speed_source='CFSR',
            stability_iterations=None,
            albedo_iterations=10,
            rs_interp_flag=True,
            ta_smooth_flag=True,
            latitude=None,
            longitude=None,
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
            The default is 'CONUS_V006'.
        alexi_source : {'CONUS_V004', 'CONUS_V005', 'CONUS_V006'}
            ALEXI ET image collection ID (the default is 'CONUS_V006').
        lai_source : string
            LAI image collection ID or "openet-landsat-lai" to compute dynamically
            using the openet-landsat-lai module.
        lst_source : string
            LST image collection ID.
        elevation_source : str, ee.Image
            Elevation source keyword or asset (the default is USGS/SRTMGL1_003).
            Units must be in meters.
        landcover_source : {'USGS/NLCD_RELEASES/2021_REL/NLCD',
                            'USGS/NLCD_RELEASES/2019_REL/NLCD',
                            'USGS/NLCD_RELEASES/2019_REL/NLCD/2016',
                            'USGS/NLCD_RELEASES/2021_REL/NLCD/2021'}
            Land cover source collection or image ID
            (the default is 'USGS/NLCD_RELEASES/2021_REL/NLCD').
        air_pres_source: {'CFSR'}
            Air pressure source keyword (the default is 'CFSR').
        air_temp_source : {'CFSR'}
            Air temperature source keyword (the default is 'CFSR').
        rs_daily_source : {'CFSR'}
            Daily solar insolation source keyword (the default is 'CFSR').
        rs_hourly_source : {'CFSR'}
            Hourly solar insolation source keyword (the default is 'CFSR').
        vapor_pres_source : {'CFSR'}
            Vapour pressure source keyword (the default is 'CFSR').
        wind_speed_source : {'CFSR'}
            Wind speed source keyword (the default is 'CFSR').
        stability_iterations : int, optional
            Number of stability calculation iterations.  If not set, the
            number will be computed dynamically.
        albedo_iterations : int, optional
            Number albedo separation iterations (the default is 10).
        rs_interp_flag : bool, optional
            If True, interpolate incoming solar radiation.
            If False, select image with same date and hour.
            The default is True.
        ta_smooth_flag : bool, optional
            If True, smooth and resample Ta image (the default is True).
        latitude : ee.Image, optional
            Latitude [deg].  If not set will default to ee.Image.pixelLonLat().
        longitude : ee.Image, optional
            Longitude [deg].  If not set will default to ee.Image.pixelLonLat().
        et_min : float, optional
            Minimum output ET value. (the default is 0.01).
        kwargs : dict, optional
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
        # scene_id = ee.List(ee.String(self.index).split('_')).slice(-3)
        # self.scene_id = (
        #     ee.String(scene_id.get(0)).cat('_')
        #     .cat(ee.String(scene_id.get(1))).cat('_')
        #     .cat(ee.String(scene_id.get(2)))
        # )

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
        self.cloud_mask = (
            image.select('cfmask').gte(1)
            .reduceNeighborhood(ee.Reducer.min(), ee.Kernel.circle(30, 'meters'))
            .reduceNeighborhood(ee.Reducer.max(), ee.Kernel.circle(90, 'meters'))
            .eq(0)
            .rename(['cloud_mask'])
        )
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
        self.lst_source = lst_source
        self.elevation_source = elevation_source
        self.landcover_source = landcover_source
        self.air_pres_source = air_pres_source
        self.air_temp_source = air_temp_source
        self.rs_daily_source = rs_daily_source
        self.rs_hourly_source = rs_hourly_source
        self.vapor_pres_source = vapor_pres_source
        self.wind_speed_source = wind_speed_source
        if stability_iterations:
            self.stabil_iter = int(stability_iterations + 0.5)
        else:
            self.stabil_iter = None
        self.albedo_iter = int(albedo_iterations + 0.5)
        self.rs_interp_flag = utils.boolean(rs_interp_flag)
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
        if self.et_reference_factor and not utils.is_number(self.et_reference_factor):
            raise ValueError('et_reference_factor must be a number')
        if self.et_reference_factor and self.et_reference_factor < 0:
            raise ValueError('et_reference_factor must be greater than zero')
        resample_methods = ['nearest', 'bilinear', 'bicubic']
        if (self.et_reference_resample and
                self.et_reference_resample.lower() not in resample_methods):
            raise ValueError('Unsupported et_reference_resample method\n')

        if latitude is None:
            self.latitude = self.ndvi.multiply(0).add(ee.Image.pixelLonLat().select(['latitude']))
        elif utils.is_number(latitude):
            self.latitude = ee.Image.constant(latitude)
        elif isinstance(latitude, ee.computedobject.ComputedObject):
            self.latitude = latitude
        else:
            raise ValueError('invalid lat parameter')

        if longitude is None:
            self.longitude = self.ndvi.multiply(0).add(ee.Image.pixelLonLat().select(['longitude']))
        elif utils.is_number(longitude):
            self.longitude = ee.Image.constant(longitude)
        elif isinstance(longitude, ee.computedobject.ComputedObject):
            self.longitude = longitude
        else:
            raise ValueError('invalid lon parameter')

        # TODO: This is now getting set in landcover function
        #   Only NLCD types are currently supported
        self.lc_type = 'NLCD'
        self.set_landcover_vars()

        # Image projection and geotransform
        self.crs = image.projection().crs()
        self.transform = ee.List(ee.Dictionary(
            ee.Algorithms.Describe(image.projection())).get('transform'))
        # self.crs = image.select([0]).projection().getInfo()['crs']
        # self.transform = image.select([0]).projection().getInfo()['transform']

        # CGM - Hardcoding to the CONUS transform and ancillary assets since they
        #   are identical for all versions and global is not yet supported
        self.alexi_geo = [0.04, 0, -125.04, 0, -0.04, 49.8]
        self.alexi_crs = 'EPSG:4326'
        self.alexi_elev = ee.Image('projects/openet/assets/alexi/ancillary/conus/v006/elevation')
        self.alexi_lat = ee.Image('projects/openet/assets/alexi/ancillary/conus/v006/latitude')
        self.alexi_lon = ee.Image('projects/openet/assets/alexi/ancillary/conus/v006/longitude')

        # # CGM - This should probably be set in et_alexi() but that wasn't working
        # if type(self.alexi_source) is str:
        #     if self.alexi_source.upper() == 'CONUS_V006':
        #         self.alexi_geo = [0.04, 0, -125.02, 0, -0.04, 49.78]
        #         self.alexi_crs = 'EPSG:4326'
        #     elif self.alexi_source.upper() == 'CONUS_V005':
        #         self.alexi_geo = [0.04, 0, -125.02, 0, -0.04, 49.78]
        #         self.alexi_crs = 'EPSG:4326'
        #     elif self.alexi_source.upper() == 'CONUS_V004':
        #         self.alexi_geo = [0.04, 0, -125.02, 0, -0.04, 49.78]
        #         self.alexi_crs = 'EPSG:4326'
        #     else:
        #         # Assume ALEXI source is an image collection ID if it is a string
        #         #   but doesn't match on any of the keywords.
        #         alexi_img = ee.Image(ee.ImageCollection(self.alexi_source).first())
        #         self.alexi_geo = ee.List(ee.Dictionary(
        #             ee.Algorithms.Describe(alexi_img.projection())).get('transform'))
        #         self.alexi_crs = alexi_img.projection().crs()
        # else:
        #     self.alexi_geo = [0.04, 0, -125.04, 0, -0.04, 49.8]
        #     self.alexi_crs = 'EPSG:4326'

    @classmethod
    def from_image_id(cls, image_id, **kwargs):
        """Constructs a DisALEXI Image instance from an image ID

        Parameters
        ----------
        image_id : str
            An earth engine image ID.
            (i.e. 'LANDSAT/LC08/C02/T1_L2/LC08_044033_20170716')
        kwargs
            Keyword arguments to pass through to model init.

        Returns
        -------
        new instance of Image class

        """
        # CGM - Should the supported image collection IDs and helper
        # function mappings be set in a property or method of the Image class?
        collection_methods = {
            'LANDSAT/LC09/C02/T1_L2': 'from_landsat_c02_l2',
            'LANDSAT/LC08/C02/T1_L2': 'from_landsat_c02_l2',
            'LANDSAT/LE07/C02/T1_L2': 'from_landsat_c02_l2',
            'LANDSAT/LT05/C02/T1_L2': 'from_landsat_c02_l2',
            'LANDSAT/LT04/C02/T1_L2': 'from_landsat_c02_l2',
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
    def from_landsat_c02_l2(cls, sr_image, **kwargs):
        """Returns a DisALEXI Image instance from a Landsat C02 Level 2 (SR) image

        Parameters
        ----------
        sr_image : ee.Image
            A raw Landsat Collection 2 SR image.
        kwargs : dict
            Keyword arguments to pass through to Image init function.

        Returns
        -------
        Instance of Image class

        """
        return cls(landsat.Landsat_C02_L2(sr_image).prep(), **kwargs)

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
                output_images.append(self.lst.float())
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
            t_air=self.ta, t_rad=self.lst, t_air0=self.air_temperature,
            e_air=self.vapor_pressure, u=self.wind_speed, p=self.air_pressure,
            z=self.elevation, rs_1=self.rs1, rs24=self.rs24, vza=0,
            aleafv=self.aleafv, aleafn=self.aleafn, aleafl=self.aleafl,
            adeadv=self.adeadv, adeadn=self.adeadn, adeadl=self.adeadl,
            albedo=self.albedo, ndvi=self.ndvi, lai=self.lai,
            clump=self.clump, leaf_width=self.leaf_width,
            hc_min=self.hc_min, hc_max=self.hc_max,
            datetime=self.datetime, lat=self.latitude, lon=self.longitude,
            stabil_iter=self.stabil_iter, albedo_iter=self.albedo_iter,
            et_min=self.et_min,
        )

        # Filter out pixels that are very different than ALEXI values
        # Calculate the relative difference between aggregated ET and ALEXI as a percentage value
        # Apply different threshold of the relative difference for different ALEXI range.
        #   ALEXI <= 0.1, we use 200% as threshold
        #   0.1 < ALEXI <= 1.0, we use 150% as threshold
        #   1.0 < ALEXI <= 8.0, we use 50% as threshold
        #   ALEXI >= 8.0, we use 30% as threshold
        et_coarse_new = (
            et.reproject(crs=self.crs, crsTransform=self.transform)
            .reduceResolution(reducer=ee.Reducer.mean().unweighted(), maxPixels=30000)
            .reproject(crs=self.alexi_crs, crsTransform=self.alexi_geo)
        )
        alexi_diff = et_coarse_new.subtract(self.et_alexi).abs().divide(self.et_alexi)
        combined_et_mask = (
            alexi_diff.lt(2.0).And(self.et_alexi.lte(0.1))
            .Or(alexi_diff.lt(1.5).And(self.et_alexi.gt(0.1)).And(self.et_alexi.lte(1.0)))
            .Or(alexi_diff.lt(0.5).And(self.et_alexi.gt(1.0)).And(self.et_alexi.lte(8.0)))
            .Or(alexi_diff.lt(0.3).And(self.et_alexi.gt(8.0)))
        )

        # Remove retile call if et masking above is not applied!
        return (
            et.rename(['et'])
            .updateMask(combined_et_mask)
            .retile(128)
            .set(self.properties)
        )

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
                ee.ImageCollection(self.et_reference_source)
                .filterDate(self.start_date, self.end_date)
                .select([self.et_reference_band])
                .first()
            )
        else:
            raise ValueError(f'Unsupported etr_source: {self.etr_source}\n')

        # Map ETr values directly to the input (i.e. Landsat) image pixels
        # The benefit of this is the ETr image is now in the same crs as the
        #   input image.  Not all models may want this though.
        # CGM - Should the output band name match the input ETr band name?
        return (
            self.ndvi.multiply(0).add(et_reference_img)
            .multiply(self.et_reference_factor)
            .rename(['et_reference']).set(self.properties)
        )

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
    def landcover(self):
        """Landcover"""
        if utils.is_number(self.landcover_source):
            lc_img = ee.Image.constant(int(self.landcover_source)).rename(['landcover'])
            self.lc_type = 'NLCD'
        elif isinstance(self.landcover_source, ee.computedobject.ComputedObject):
            # If the source is an ee.Image assume it is an NLCD image
            lc_img = self.landcover_source.rename(['landcover'])
            self.lc_type = 'NLCD'
        elif (re.match('USGS/NLCD_RELEASES/2021_REL/NLCD/\\d{4}', self.landcover_source, re.I) or
              re.match('USGS/NLCD_RELEASES/2019_REL/NLCD/\\d{4}', self.landcover_source, re.I)):
            # Assume an NLCD image ID was passed in and use it directly
            lc_img = ee.Image(self.landcover_source.upper()).select(['landcover'])
            self.lc_type = 'NLCD'
        elif re.match('USGS/NLCD_RELEASES/2016_REL/\\d{4}', self.landcover_source.upper()):
            warnings.warn(
                'The NLCD 2016 release is deprecated and support will be removed in a future version',
                FutureWarning
            )
            lc_img = ee.Image(self.landcover_source.upper()).select(['landcover'])
            self.lc_type = 'NLCD'
        elif self.landcover_source.upper() == 'USGS/NLCD_RELEASES/2021_REL/NLCD':
            # Automatically switch to the 2019 release for all years before 2020
            #   since the 2021 Release does not currently contain earlier images
            # Use the 2021 release and image for all years on or after 2020
            #   but intentionally limit end year to 2021 for now
            year_remap = ee.Dictionary({
                '1999': 2001, '2000': 2001, '2001': 2001, '2002': 2001,
                '2003': 2004, '2004': 2004, '2005': 2004,
                '2006': 2006, '2007': 2006,
                '2008': 2008, '2009': 2008,
                '2010': 2011, '2011': 2011, '2012': 2011,
                '2013': 2013, '2014': 2013,
                '2015': 2016, '2016': 2016, '2017': 2016,
                '2018': 2019, '2019': 2019,
                '2020': 2021, '2021': 2021,
            })
            nlcd_year = year_remap.get(self.year.min(2021).max(1999).format('%d'))
            # The filterDate calls are probably not needed,
            # but adding just in case additional years get added to either collection
            lc_img = (
                ee.ImageCollection(self.landcover_source)
                .filterDate('2020-01-01', '2022-01-01')
                .select(['landcover'])
                .merge(
                    ee.ImageCollection('USGS/NLCD_RELEASES/2019_REL/NLCD')
                    .filterDate('2001-01-01', '2020-01-01')
                    .select(['landcover'])
                )
                .filter(ee.Filter.calendarRange(nlcd_year, nlcd_year, 'year'))
                .first()
            )
            self.lc_type = 'NLCD'
            # TODO: Save the actual land cover source image as a property?
            # self.landcover_source = lc_img.get('system:index')
        elif self.landcover_source.upper() == 'USGS/NLCD_RELEASES/2019_REL/NLCD':
            year_remap = ee.Dictionary({
                '1999': 2001, '2000': 2001, '2001': 2001, '2002': 2001,
                '2003': 2004, '2004': 2004, '2005': 2004,
                '2006': 2006, '2007': 2006,
                '2008': 2008, '2009': 2008,
                '2010': 2011, '2011': 2011, '2012': 2011,
                '2013': 2013, '2014': 2013,
                '2015': 2016, '2016': 2016, '2017': 2016,
                '2018': 2019, '2019': 2019,
            })
            nlcd_year = year_remap.get(self.year.min(2019).max(1999).format('%d'))
            lc_img = (
                ee.ImageCollection(self.landcover_source)
                .filter(ee.Filter.calendarRange(nlcd_year, nlcd_year, 'year'))
                .first().select(['landcover'])
            )
            self.lc_type = 'NLCD'
        elif self.landcover_source.upper() == 'USGS/NLCD_RELEASES/2016_REL':
            warnings.warn(
                'The NLCD 2016 release is deprecated and support will be removed in a future version',
                FutureWarning
            )
            year_remap = ee.Dictionary({
                '1999': '2001', '2000': '2001', '2001': '2001', '2002': '2001',
                '2003': '2004', '2004': '2004', '2005': '2004',
                '2006': '2006', '2007': '2006',
                '2008': '2008', '2009': '2008',
                '2010': '2011', '2011': '2011', '2012': '2011',
                '2013': '2013', '2014': '2013',
                '2015': '2016', '2016': '2016',
            })
            nlcd_year = year_remap.get(self.year.min(2016).max(1999).format('%d'))
            lc_img = (
                ee.ImageCollection(self.landcover_source)
                .filter(ee.Filter.equals('system:index', nlcd_year))
                .first().select(['landcover'])
            )
            self.lc_type = 'NLCD'
        else:
            raise ValueError(f'Unsupported landcover_source: {self.landcover_source}\n')
        # TODO: Test out the ESA 10m Landcover asset
        #   Collection ID: ESA/WorldCover/v100
        #   Only 2020 is currently available
        # elif self.landcover_source == 'ESA/WorldCover/v100/2020':
        #     lc_img = ee.Image(self.landcover_source).select(['Map'], ['landcover'])
        #     self.lc_type = 'ESA?'
        # TODO: Test out the Copernicus 100m landcover asset
        #   Collection ID: COPERNICUS/Landcover/100m/Proba-V-C3/Global
        #   Years: 2015-2019
        # elif self.landcover_source == 'COPERNICUS/Landcover/100m/Proba-V-C3/Global':
        # elif self.landcover_source.startswith('COPERNICUS/Landcover/100m/Proba-V-C3/Global/2'):
        #     lc_img = (
        #         ee.ImageCollection(self.landcover_source)
        #         .filterDate(ee.Date(self.year, 1, 1), ee.Date(self.year, 12, 31))
        #         .select(['discrete_classification])
        #         .rename(['landcover'])
        #     self.lc_type = ''

        return lc_img

    @lazy_property
    def lai(self):
        """Leaf Area Index (LAI)"""
        if utils.is_number(self.lai_source):
            lai_img = ee.Image.constant(float(self.lai_source))
        elif (type(self.lai_source) is str and
              self.lai_source.lower() in ['openet-landsat-lai', 'openet-lai']):
            # The LAI module only accepts python strings as input for the image ID
            # Calling .getInfo() here to avoid needing to modify the LAI module (for now)
            # This will make it NOT possible to map this function
            #   and in the long run the LAI module should be modified to support either
            #   ee.Strings or ee.Images
            lai_img = openet.lai.Landsat(utils.getinfo(self.id)).lai(nonveg=True)
            # lai_img = openet.lai.Landsat(utils.getinfo(self.index)).lai(nonveg=True)

            # This approach for getting the version will only work with the new LAI v0.2.0
            #   that was built with a pyproject.toml and pins DisALEXI to Python 3.8+
            self.landsat_lai_version = metadata.metadata('openet-landsat-lai')['Version']
        elif type(self.lai_source) is str:
            # Assumptions (for now)
            #   String lai_source is an image collection ID
            #   Images are single band and don't need a select()
            #   LAI images always need to be scaled
            # CGM - This will raise a .get() error if the image doesn't exist
            lai_coll = ee.ImageCollection(self.lai_source)\
               .filterMetadata('scene_id', 'equals', self.index)
            lai_img = ee.Image(lai_coll.first()).select(['LAI'])
            lai_img = lai_img.multiply(ee.Number(lai_img.get('scale_factor')))\
                .set({'landsat_lai_version': lai_img.get('landsat_lai_version')})
            self.landsat_lai_version = lai_img.get('landsat_lai_version')
        else:
            raise ValueError(f'Unsupported lai_source: {self.lai_source}\n')
        # elif isinstance(self.lai_source, ee.computedobject.ComputedObject):
        #     lai_img = self.lai_source

        # Apply the cloud mask since LAI was not read from the masked input image
        return lai_img.select([0], ['lai']).updateMask(self.cloud_mask)

    @lazy_property
    def lst(self):
        """Sharpened land surface temperature (LST)"""
        if utils.is_number(self.lst_source):
            lst_img = ee.Image.constant(float(self.lst_source))
        elif type(self.lst_source) is str:
            # Assumptions (for now)
            #   String lst_source is an image collection ID
            #   Images are single band and don't need a select()
            #   LST images always need to be scaled
            # CGM - This will raise a .get() error if the image doesn't exist
            lst_coll = ee.ImageCollection(self.lst_source)\
                .filterMetadata('scene_id', 'equals', self.index)
            lst_img = ee.Image(lst_coll.first())
            lst_img = lst_img.multiply(ee.Number(lst_img.get('scale_factor')))\
                .set({'landsat_lst_version': lst_img.get('landsat_lst_version')})
            self.landsat_lst_version = lst_img.get('landsat_lst_version')
        else:
            raise ValueError(f'Unsupported lst_source: {self.lst_source}\n')
        # elif isinstance(self.lst_source, ee.computedobject.ComputedObject):
        #     lst_img = self.lst_source

        # Apply the cloud mask since LST was not read from the masked input image
        return lst_img.select([0], ['lst']).updateMask(self.cloud_mask)

    @lazy_property
    def mask(self):
        """Mask of all active pixels (based on the final et)"""
        return (
            self.et.multiply(0).add(1).uint8().updateMask(1)
            .rename(['mask']).set(self.properties)
        )

    @lazy_property
    def time(self):
        """Return an image of the 0 UTC time (in milliseconds)"""
        return (
            self.mask
            .double().multiply(0).add(utils.date_to_time_0utc(self.datetime))
            .rename(['time']).set(self.properties)
        )

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
            'CONUS_V004': 'projects/openet/disalexi/tair/conus_v004_1k',
            'CONUS_V005': 'projects/openet/disalexi/tair/conus_v005_1k',
            'CONUS_V006': 'projects/openet/disalexi/tair/conus_v006_1k',
            # 'CONUS_V006': 'projects/openet/disalexi/tair/conus_v006',
            # 'CONUS_V007': 'projects/openet/disalexi/tair/conus_v007',
        }
        ta_source_re = re.compile(
            '(projects/earthengine-legacy/assets/)?'
            'projects/(\\w+/)?(assets/)?'
            'disalexi/ta(ir)?/(conus|global)_v\\d{3}\\w?(_\\w+)?',
            re.IGNORECASE
        )

        if utils.is_number(self.ta_source):
            ta_img = ee.Image.constant(float(self.ta_source))
            #     .set({'ta_iteration': 'constant'})
        elif isinstance(self.ta_source, ee.computedobject.ComputedObject):
            ta_img = ee.Image(self.ta_source)
        elif ((self.ta_source.upper() in ta_keyword_sources.keys()) and
              (('1k' in ta_keyword_sources[self.ta_source.upper()]) or
               ('10k' in ta_keyword_sources[self.ta_source.upper()]))):
            ta_coll_id = ta_keyword_sources[self.ta_source.upper()]
            ta_coll = (
                ee.ImageCollection(ta_coll_id)
                .filterMetadata('image_id', 'equals', self.id)
                .limit(1, 'step_size', False)
            )
            ta_img = ta_mosaic_min_bias(ee.Image(ta_coll.first()))
            if self.ta_smooth_flag:
                ta_img = (
                    ta_img.focal_mean(1, 'circle', 'pixels')
                    .reproject(crs=self.alexi_crs, crsTransform=self.alexi_geo)
                    .resample('bilinear')
                    .reproject(crs=self.crs, crsTransform=self.transform)
                )
        elif (ta_source_re.match(self.ta_source) and
              (('1k' in self.ta_source) or ('10k' in self.ta_source))):
            ta_coll = (
                ee.ImageCollection(self.ta_source)
                .filterMetadata('image_id', 'equals', self.id)
                .limit(1, 'step_size', False)
            )
            ta_img = ta_mosaic_min_bias(ee.Image(ta_coll.first()))
            if self.ta_smooth_flag:
                ta_img = (
                    ta_img.focal_mean(1, 'circle', 'pixels')
                    .reproject(crs=self.alexi_crs, crsTransform=self.alexi_geo)
                    .resample('bilinear')
                    .reproject(crs=self.crs, crsTransform=self.transform)
                )
        elif ta_source_re.match(self.ta_source):
            # For now assuming Ta source has the correct band (ta_smooth or ta_interp)
            #   and an image_id property
            ta_img = ee.Image(
                ee.ImageCollection(self.ta_source)
                .filterMetadata('image_id', 'equals', self.id)
                .first()
                .select('ta_smooth' if self.ta_smooth_flag else 'ta_interp')
                .resample('bilinear')
                .reproject(crs=self.crs, crsTransform=self.transform)
            )
        else:
            raise ValueError(f'Unsupported ta_source: {self.ta_source}\n')
        # if self.ta_source is None:
        #     raise ValueError('ta_source must be set to compute et')

        return ta_img.rename(['ta']).set(self.properties)

    @lazy_property
    def air_pressure_coarse(self):
        """Air pressure [kPa] at coarse (ALEXI) scale"""
        if utils.is_number(self.air_pres_source):
            ap_img = ee.Image.constant(float(self.air_pres_source))
        elif isinstance(self.air_pres_source, ee.computedobject.ComputedObject):
            ap_img = self.air_pres_source
        elif self.air_pres_source.upper() == 'ESTIMATE':
            # pressure = 101.3 * (((293.0 - 0.0065 * elevation) / 293.0) ** 5.26)
            ap_img = self.alexi_elev.multiply(-0.0065).add(293).divide(293.0).pow(5.26).multiply(101.3)
        elif self.air_pres_source.upper() == 'CFSR':
            ap_coll_id = 'projects/disalexi/meteo_data/airpressure/global_v001_3hour'
            ap_coll = ee.ImageCollection(ap_coll_id).select(['airpressure'])
            ap_img = utils.interpolate(ap_coll, self.datetime, timestep=3)
            ap_img = (
                ee.Image(ap_img)
                .resample('bicubic')
                .reproject(crs=self.alexi_crs, crsTransform=self.alexi_geo)
                .multiply(self.alexi_elev.multiply(-0.0065).add(293).divide(293.0).pow(5.26))
                # # Resample/reproject to the Landsat scale will happen in non-coarse function
                # .resample('bilinear')
                # .reproject(crs=self.crs, crsTransform=self.transform)
            )
        else:
            raise ValueError(f'Invalid air_pres_source: {self.air_pres_source}\n')

        return ee.Image(ap_img).rename(['air_pressure'])

    @lazy_property
    def air_pressure(self):
        """Air pressure [kPa] at fine (Landsat) scale"""
        return (
            self.air_pressure_coarse
            .resample('bilinear')
            .reproject(crs=self.crs, crsTransform=self.transform)
        )

    @lazy_property
    def air_temperature_coarse(self):
        """Air temperature [K] at coarse (ALEXI) scale"""
        if utils.is_number(self.air_temp_source):
            at_img = ee.Image.constant(float(self.air_temp_source))
        # elif isinstance(self.air_temp_source, ee.computedobject.ComputedObject):
        #     air_temp_img = self.air_temp_source
        elif self.air_temp_source.upper() == 'CFSR':
            at_coll_id = 'projects/disalexi/meteo_data/airtemperature/global_v001_3hour'
            at_coll = ee.ImageCollection(at_coll_id).select(['temperature'])
            at_img = utils.interpolate(at_coll, self.datetime, timestep=3)
            at_img = (
                ee.Image(at_img)
                .resample('bicubic')
                .reproject(crs=self.alexi_crs, crsTransform=self.alexi_geo)
                # # Resample/reproject to the Landsat scale will happen in non-coarse function
                #.resample('bilinear')
                # .reproject(crs=self.crs, crsTransform=self.transform)
            )
        else:
            raise ValueError(f'Unsupported air_temp_source: {self.air_temp_source}\n')

        return ee.Image(at_img).rename(['air_temperature'])

    @lazy_property
    def air_temperature(self):
        """Air temperature [K] at fine (Landsat) scale"""
        return (
            self.air_temperature_coarse
            .resample('bilinear')
            .reproject(crs=self.crs, crsTransform=self.transform)
        )
    
    @lazy_property
    def rs1_coarse(self):
        """Hourly solar insolation [W m-2] at coarse (ALEXI) scale"""
        if utils.is_number(self.rs_hourly_source):
            rs1_img = ee.Image.constant(float(self.rs_hourly_source))
        elif isinstance(self.rs_hourly_source, ee.computedobject.ComputedObject):
            rs1_img = self.rs_hourly_source
        elif self.rs_hourly_source.upper() == 'CFSR':
            rs1_coll_id = 'projects/disalexi/insol_data/global_v001_hourly'
            rs1_coll = ee.ImageCollection(rs1_coll_id).select(['insolation'])
            if self.rs_interp_flag:
                # TODO: Check if the CFSR insolation are instantaneous or accumulations
                rs1_img = utils.interpolate(rs1_coll, self.datetime, timestep=1)
            else:
                # Select the source image before the image time
                # This would be off if image time was "exactly" on the hour
                # Could change the offsets to -0.5 hours to select the closest image
                rs1_img = (
                    rs1_coll
                    .filterDate(self.datetime.advance(-1, 'hour'),
                                self.datetime.advance(0, 'hour'))
                    .first()
                )
            rs1_img = (
                ee.Image(rs1_img)
                .resample('bicubic')
                .reproject(crs=self.alexi_crs, crsTransform=self.alexi_geo)
                # # Resample/reproject to the Landsat scale will happen in non-coarse function
                #.resample('bilinear')
                # .reproject(crs=self.crs, crsTransform=self.transform)
            )
        else:
            raise ValueError(f'Unsupported rs_hourly_source: {self.rs_hourly_source}\n')

        return ee.Image(rs1_img).rename(['rs'])

    @lazy_property
    def rs1(self):
        """Hourly solar insolation [W m-2] at fine (Landsat) scale"""
        return (
            self.rs1_coarse
            .resample('bilinear')
            .reproject(crs=self.crs, crsTransform=self.transform)
        )

    @lazy_property
    def rs24_coarse(self):
        """Daily (24 hour) solar insolation [W m-2] at coarse (ALEXI) scale"""
        if utils.is_number(self.rs_daily_source):
            rs24_img = ee.Image.constant(float(self.rs_daily_source))
        elif isinstance(self.rs_daily_source, ee.computedobject.ComputedObject):
            rs24_img = self.rs_daily_source
        elif self.rs_daily_source.upper() == 'CFSR':
            rs24_coll_id = 'projects/disalexi/insol_data/global_v001_hourly'
            rs24_coll = (
                ee.ImageCollection(rs24_coll_id)
                .filterDate(self.datetime.advance(-8, 'hours'),
                            self.datetime.advance(12, 'hours'))
            )
            rs24_img = (
                ee.Image(rs24_coll.sum())
                .resample('bicubic')
                .reproject(crs=self.alexi_crs, crsTransform=self.alexi_geo)
                # # Resample/reproject to the Landsat scale will happen in non-coarse function
                #.resample('bilinear')
                # .reproject(crs=self.crs, crsTransform=self.transform)
            )
        else:
            raise ValueError(f'Unsupported rs_daily_source: {self.rs_daily_source}\n')

        return ee.Image(rs24_img).rename(['rs'])

    @lazy_property
    def rs24(self):
        """Daily (24 hour) solar insolation [W m-2] at fine (Landsat) scale"""
        return (
            self.rs24_coarse
            .resample('bilinear')
            .reproject(crs=self.crs, crsTransform=self.transform)
        )

    @lazy_property
    def wind_speed_coarse(self):
        """Wind speed [m/s] at coarse (ALEXI) scale"""
        if utils.is_number(self.wind_speed_source):
            ws_img = ee.Image.constant(float(self.wind_speed_source))
        elif isinstance(self.wind_speed_source, ee.computedobject.ComputedObject):
            ws_img = self.wind_speed_source
        # elif self.wind_speed_source.upper() == 'CFSV2':
        #     # It would be more correct to compute the magnitude for each image,
        #     #   then compute the average.
        #     # Do we need daily, 6hr, or interpolated instantaneous data?
        #     ws_coll = (
        #         ee.ImageCollection('NOAA/CFSV2/FOR6H')
        #         .select(['u-component_of_wind_height_above_ground',
        #                  'v-component_of_wind_height_above_ground'])
        #         .filterDate(self.start_date, self.end_date)
        #     )
        #     ws_img = ws_coll.mean().expression('sqrt(b(0) ** 2 + b(1) ** 2)')
        elif self.wind_speed_source.upper() == 'CFSR':
            ws_coll_id = 'projects/disalexi/meteo_data/windspeed/global_v001_3hour'
            ws_coll = ee.ImageCollection(ws_coll_id).select(['windspeed'])
            ws_img = utils.interpolate(ws_coll, self.datetime, timestep=3)
            ws_img = (
                ee.Image(ws_img)
                .resample('bicubic')
                .reproject(crs=self.alexi_crs, crsTransform=self.alexi_geo)
                # # Resample/reproject to the Landsat scale will happen in non-coarse function
                # .resample('bilinear')
                # .reproject(crs=self.crs, crsTransform=self.transform)
            )
        else:
            raise ValueError(f'Unsupported wind_speed_source: {self.wind_speed_source}\n')

        # TODO: Clamping here will cause slight differences with the existing approach
        #   It might be more consistent to only clamp the fine scale image
        return ee.Image(ws_img).max(2).min(20).rename(['wind_speed'])

    @lazy_property
    def wind_speed(self):
        """Wind speed [m/s] at fine (Landsat) scale"""
        return (
            self.wind_speed_coarse
            .resample('bilinear')
            .reproject(crs=self.crs, crsTransform=self.transform)
        )

    @lazy_property
    def vapor_pressure_coarse(self):
        """Vapor pressure [kPa] at coarse (ALEXI) scale"""
        if utils.is_number(self.vapor_pres_source):
            vp_img = ee.Image.constant(float(self.vapor_pres_source))
        elif isinstance(self.vapor_pres_source, ee.computedobject.ComputedObject):
            vp_img = self.vapor_pres_source
        elif self.vapor_pres_source.upper() == 'CFSR':
            vp_coll_id = 'projects/disalexi/meteo_data/vp/global_v001_3hour'
            vp_coll = ee.ImageCollection(vp_coll_id).select(['vp'])
            vp_img = utils.interpolate(vp_coll, self.datetime, timestep=3)
            vp_img = (
                ee.Image(vp_img)
                .resample('bicubic')
                .reproject(crs=self.alexi_crs, crsTransform=self.alexi_geo)
                # # Resample/reproject to the Landsat scale will happen in non-coarse function
                # .resample('bilinear')
                # .reproject(crs=self.crs, crsTransform=self.transform)
            )
        else:
            raise ValueError(f'Unsupported vapor_pres_source: {self.vapor_pres_source}\n')

        return ee.Image(vp_img).rename(['vapor_pressure'])

    @lazy_property
    def vapor_pressure(self):
        """Vapor pressure [kPa] at fine (Landsat) scale"""
        return (
            self.vapor_pressure_coarse
            .resample('bilinear')
            .reproject(crs=self.crs, crsTransform=self.transform)
        )

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

        # Call landcover function to trigger building lc_type
        lc_img = self.landcover

        if self.lc_type.upper() not in lc_remaps.keys():
            raise KeyError(f'Invalid lc_type: {self.lc_type.upper()}\n'
                           f'Choices are {", ".join(lc_remaps.keys())})')

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
            lc_items = sorted(lc_remaps[lc_type.upper()][lc_var].items())
            input_values = [k for k, v in lc_items]
            # Scale output values by 100 since remap values must be integer
            output_values = [v * 100 for k, v in lc_items]

            # Get the remap values from the dataframe and apply to the land cover image
            return (
                ee.Image(landcover)
                .remap(input_values, output_values)
                .divide(100)
                .rename([lc_var])
            )

        # Get LC based variables
        self.aleafv = lc_remap(lc_img, self.lc_type, 'aleafv')
        self.aleafn = lc_remap(lc_img, self.lc_type, 'aleafn')
        self.aleafl = lc_remap(lc_img, self.lc_type, 'aleafl')
        self.adeadv = lc_remap(lc_img, self.lc_type, 'adeadv')
        self.adeadn = lc_remap(lc_img, self.lc_type, 'adeadn')
        self.adeadl = lc_remap(lc_img, self.lc_type, 'adeadl')
        # Yun modified to use hmax for hmin too
        self.hc_min = lc_remap(lc_img, self.lc_type, 'hmax')
        self.hc_max = lc_remap(lc_img, self.lc_type, 'hmax')
        self.leaf_width = lc_remap(lc_img, self.lc_type, 'xl')
        self.clump = lc_remap(lc_img, self.lc_type, 'omega')

    def ta_coarse(
            self,
            offsets=[-12, -7, -4, -2, -1, 0, 1, 2, 4, 7, 12],
            #offsets=[-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5],
    ):
        """Compute coarse scale air temperature estimate

        Parameters
        ----------
        offsets : list
            Air temperature offset values to use when generating mosaic

        Returns
        -------
        image : ee.Image
            ALEXI scale air temperature image

        """
        # Compute the initial Ta from the meteorology
        ta_initial_img = self.ta_coarse_initial().rename(['ta_initial'])

        # Compute the Ta mosaic from the initial Ta image
        ta_mosaic_img = self.ta_mosaic(ta_img=ta_initial_img, offsets=offsets)

        # Interpolate the minimum bias Ta from the 1k steps
        ta_interp_img = ta_mosaic_interpolate(ta_mosaic_img)

        # TODO: Should the retile parameter be hardcoded?
        #   If so, where should it go in this function?
        # ta_interp_img = ta_interp_img.retile(4)

        # Apply simple smoothing to the interpolated Ta band and save as "ta_smooth"
        # This will fill small 1 pixel holes
        if self.ta_smooth_flag:
            ta_interp_img = ta_interp_img.addBands(
                ta_interp_img.select('ta_interp')
                .focal_mean(1, 'circle', 'pixels')
                # CGM - Testing without reproject call, but it may be needed
                # .reproject(crs=self.alexi_crs, crsTransform=self.alexi_geo)
                .multiply(10).round().divide(10)
                .rename('ta_smooth')
            )

        return ta_initial_img.addBands(ta_interp_img).set(self.properties)

    def ta_coarse_initial(self):
        """Compute initial coarse scale air temperature estimate from meteorology

        Parameters
        ----------

        Returns
        -------
        image : ee.Image
            ALEXI scale air temperature image

        """

        # TODO: Test making a single multi-band image of all of these and making
        #   a single reduceResolution call
        # TODO: Test out using a median reducer instead the mean
        rr_params = {'reducer': ee.Reducer.mean().unweighted(), 'maxPixels': 30000}
        # rr_params = {'reducer': ee.Reducer.median().unweighted(), 'maxPixels': 30000}
        proj_params = {'crs': self.alexi_crs, 'crsTransform': self.alexi_geo}

        # Intentionally passing measured air temperature in for both t_air and t_air0
        # The "coarse" images are already at the ALEXI scale and shouldn't need to be reduced
        ta_invert = tseb.tseb_invert(
            et_alexi=self.et_alexi,
            t_air=self.air_temperature_coarse,
            t_air0=self.air_temperature_coarse,
            e_air=self.vapor_pressure_coarse,
            u=self.wind_speed_coarse,
            p=self.air_pressure_coarse,
            rs_1=self.rs1_coarse,
            rs24=self.rs24_coarse,
            t_rad=self.lst.reduceResolution(**rr_params).reproject(**proj_params),
            z=self.alexi_elev,
            vza=0,
            aleafv=self.aleafv.reduceResolution(**rr_params).reproject(**proj_params),
            aleafn=self.aleafn.reduceResolution(**rr_params).reproject(**proj_params),
            aleafl=self.aleafl.reduceResolution(**rr_params).reproject(**proj_params),
            adeadv=self.adeadv.reduceResolution(**rr_params).reproject(**proj_params),
            adeadn=self.adeadn.reduceResolution(**rr_params).reproject(**proj_params),
            adeadl=self.adeadl.reduceResolution(**rr_params).reproject(**proj_params),
            albedo=self.albedo.reduceResolution(**rr_params).reproject(**proj_params),
            ndvi=self.ndvi.reduceResolution(**rr_params).reproject(**proj_params),
            lai=self.lai.reduceResolution(**rr_params).reproject(**proj_params),
            clump=self.clump.reduceResolution(**rr_params).reproject(**proj_params),
            leaf_width=self.leaf_width.reduceResolution(**rr_params).reproject(**proj_params),
            hc_min=self.hc_min.reduceResolution(**rr_params).reproject(**proj_params),
            hc_max=self.hc_max.reduceResolution(**rr_params).reproject(**proj_params),
            datetime=self.datetime,
            lat=self.alexi_lat,
            lon=self.alexi_lon,
            stabil_iter=self.stabil_iter,
            albedo_iter=self.albedo_iter,
        )

        # CGM - Does the output also need the reduceResolution().reproject() call
        #   if all the inputs are already at the coarse/ALEXI scale?
        return (
            ta_invert
            .round()
            # # Round to nearest tenth
            # .multiply(10).round().divide(10)
            .rename(['ta_initial'])
            .set(self.properties)
        )

    # TODO: Maybe rename as ta_coarse_mosaic()
    def ta_mosaic(
            self,
            ta_img,
            offsets=[-12, -7, -4, -2, -1, 0, 1, 2, 4, 7, 12],
            #offsets=[-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5],
            threshold=0.5
    ):
        """Compute the air temperature for each ALEXI ET cell that minimizes
        the bias between Landsat scale ET and ALEXI ET for a fixed range of
        air temperature values

        Parameters
        ----------
        ta_img : ee.Image
        offsets : list
            Air temperature offset values to use when generating mosaic
        threshold : float, optional

        Returns
        -------
        image : ee.Image
            ALEXI scale multi-band air temperature image

        """

        # Redefining ta_func() here allows threshold to be passed into function
        def ta_func(ta):
            """Compute TSEB ET for the target t_air value"""
            et_fine = tseb.tseb_pt(
                t_air=ta, t_rad=self.lst, t_air0=self.air_temperature,
                e_air=self.vapor_pressure, u=self.wind_speed, p=self.air_pressure,
                z=self.elevation, rs_1=self.rs1, rs24=self.rs24, vza=0,
                aleafv=self.aleafv, aleafn=self.aleafn, aleafl=self.aleafl,
                adeadv=self.adeadv, adeadn=self.adeadn, adeadl=self.adeadl,
                albedo=self.albedo, ndvi=self.ndvi, lai=self.lai,
                clump=self.clump, leaf_width=self.leaf_width,
                hc_min=self.hc_min, hc_max=self.hc_max,
                datetime=self.datetime, lat=self.latitude, lon=self.longitude,
                stabil_iter=self.stabil_iter, albedo_iter=self.albedo_iter,
                et_min=self.et_min,
            )

            # TODO: Test if these two calls can be combined into one
            # Aggregate the Landsat scale ET up to the ALEXI scale
            et_coarse = (
                ee.Image(et_fine)
                .reproject(crs=self.crs, crsTransform=self.transform)
                .reduceResolution(reducer=ee.Reducer.mean().unweighted(), maxPixels=30000)
                .reproject(crs=self.alexi_crs, crsTransform=self.alexi_geo)
            )

            # Count the number of Landsat pixels in the aggregation
            agg_count_reducer = (
                ee.Reducer.count()
                .combine(reducer2=ee.Reducer.countEvery(), outputPrefix='all', sharedInputs=False)
            )
            et_agg = (
                ee.Image(et_fine)
                .reproject(crs=self.crs, crsTransform=self.transform)
                .reduceResolution(reducer=agg_count_reducer, maxPixels=30000)
                .reproject(crs=self.alexi_crs, crsTransform=self.alexi_geo)
            )
            et_perc = et_agg.select(['et_count']).divide(et_agg.select(['et_allcount']))

            bias = et_coarse.subtract(self.et_alexi)
            return (
                ee.Image([ta, bias])
                .updateMask(et_perc.gt(threshold))
                .rename(['ta', 'bias'])
                .set({'system:index': ee.Image(ta).get('system:index')})
            )

        ta_coll = ee.ImageCollection([
            ta_img.add(j).set({'system:index': 'step_{:02d}'.format(i)})
            for i, j in enumerate(offsets)
        ])

        return ee.ImageCollection(ta_coll.map(ta_func)).toBands()


# TODO: Maybe rename also, ta_mosaic_zero_bias()?
def ta_mosaic_interpolate(ta_mosaic_img):
    """Interpolate the air temperature for a bias of 0 from a mosaic stack

    Parameters
    ----------
    ta_mosaic_img : ee.Image

    Returns
    -------
    image : ee.Image
        ALEXI scale single band air temperature image

    """
    # Reverse the band order so that we can find the last transition
    #   from decreasing to increasing with a positive bias
    ta_bands = ta_mosaic_img.select('step_\\d+_ta').bandNames().reverse()
    bias_bands = ta_mosaic_img.select('step_\\d+_bias').bandNames().reverse()
    ta_array = ta_mosaic_img.select(ta_bands).toArray()
    bias_array = ta_mosaic_img.select(bias_bands).toArray()

    # Assign the bias that are very similar a very large value so that they will not be selected
    diff_array = bias_array.arraySlice(0, 1).subtract(bias_array.arraySlice(0, 0, -1))
    adj_bias_mask = diff_array.abs().lt(0.001)
    # repeat the last value to make the array the same length. array is reversed order.
    adj_bias_mask = adj_bias_mask.arrayCat(adj_bias_mask.arraySlice(0, -1), 0)
    adj_bias_array = bias_array.add(adj_bias_mask.multiply(99))

    # Identify the "last" transition from a negative to positive bias
    # CGM - Having problems with .ceil() limiting result to the image data range
    #   Multiplying by a big number seemed to fix the issue but this could still
    #     be a problem with the bias ranges get really small
    sign_array = bias_array.multiply(1000).ceil().max(0).min(1).int()
    transition_array = sign_array.arraySlice(0, 0, -1).subtract(sign_array.arraySlice(0, 1))
    # Insert an extra value at the beginning (of reverse, so actually at end)
    #   of the transition array so the indexing lines up for all steps
    transition_array = bias_array.arraySlice(0, 0, 1).multiply(0).arrayCat(transition_array, 0)
    transition_index = transition_array.arrayArgmax().arrayFlatten([['index']])
    # Get the max transition value in order to know if there was a transition
    transition_max = transition_array.arrayReduce(ee.Reducer.max(), [0]).arrayFlatten([['max']])

    # Identify the position of minimum absolute bias
    min_bias_index = adj_bias_array.abs().multiply(-1).arrayArgmax().arrayFlatten([['index']])

    # Identify the "bracketing" Ta and bias values
    # If there is a transition, use the "last" transition
    # If there is not a transition, use the minimum absolute bias for both
    # Note, the index is for the reversed arrays
    # B is the "high" value, A is the "low value"
    index_b = transition_index.subtract(1).max(0).where(transition_max.eq(0), min_bias_index)
    index_a = (
        transition_index.min(ta_bands.size().subtract(1))
        .where(transition_max.eq(0), min_bias_index)
    )
    ta_b = ta_array.arrayGet(index_b)
    ta_a = ta_array.arrayGet(index_a)
    bias_b = bias_array.arrayGet(index_b)
    bias_a = bias_array.arrayGet(index_a)

    # Linearly interpolate Ta
    # Limit the interpolated value to the bracketing values (don't extrapolate)
    ta_img = (
        ta_b.subtract(ta_a).divide(bias_b.subtract(bias_a))
        .multiply(bias_a.multiply(-1)).add(ta_a)
        .max(ta_a).min(ta_b)
    )
    # # Compute the target Ta as the average of the bracketing Ta values
    # #   instead of interpolating
    # ta_img = ta_a.add(ta_b).multiply(0.5)

    # Mask out Ta cells outside the interpolation range
    # if extrapolate_mask:
    ta_img = ta_img.updateMask(ta_a.lt(ta_b))

    # # CGM - This is mostly needed at the 10k step size
    # #   Commenting out for now
    # # Mask out Ta cells with all negative biases
    # ta_img = ta_img.updateMask(bias_b.lt(0).And(bias_a.lt(0)).Not())

    # Round to the nearest tenth (should it be hundredth?)
    return (
        ta_img.multiply(10).round().divide(10)
        .addBands([ta_a, bias_a, ta_b, bias_b])
        .rename(['ta_interp', 'ta_a', 'bias_a', 'ta_b', 'bias_b'])
    )


def ta_mosaic_min_bias(ta_mosaic_img):
    """Return the air temperature with the minimum bias from a mosaic stack

    Parameters
    ----------
    ta_mosaic_img : ee.Image

    Returns
    -------
    image : ee.Image
        ALEXI scale single band air temperature image

    """
    # Select the Ta image with the minimum bias
    ta_array = ta_mosaic_img.select('step_\\d+_ta').toArray()
    bias_array = ta_mosaic_img.select('step_\\d+_bias').toArray()
    diff = bias_array.arraySlice(0, 1).subtract(bias_array.arraySlice(0, 0, -1))
    bias_array_mask = diff.abs().lt(0.001)
    # Make the array the same length by repeating the first value
    bias_array_mask = bias_array_mask.arraySlice(0, 0, 1).arrayCat(bias_array_mask, 0)
    bias_array_new = bias_array.add(bias_array_mask.multiply(99))
    index = (
        bias_array_new.abs().multiply(-1).arrayArgmax()
        .arraySlice(0, 0, 1).arrayFlatten([['array']])
    )
    ta_img = ta_array.arrayGet(index)

    return ta_img.rename(['ta'])
