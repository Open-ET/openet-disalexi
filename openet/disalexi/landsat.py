import ee

# import openet.core.common as common


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


class Landsat(object):
    """"""
    def __init__(self):
        """"""
        pass

    @lazy_property
    def _albedo(self):
        """Compute total shortwave broadband albedo following [Liang2001]

        The Python code had the following line and comment:
            "bands = [1, 3, 4, 5, 7]  # dont use blue"
        IDL code and [Liang2001] indicate that the green band is not used.
        Coefficients were derived for Landsat 7 ETM+, but were found to be
            "suitable" to Landsat 4/5 TM also.

        Parameters
        ----------
        self.input_image : ee.Image
            "Prepped" Landsat image with standardized band names.

        Returns
        -------
        albedo : ee.Image

        References
        ----------
        .. [Liang2001] Shunlin Liang (2001),
            Narrowband to broadband conversions of land surface albedo -
            I Algorithms, Remote Sensing of Environment,
            Volume 76, Issue2, Pages 213-238,
            http://doi.org/10.1016/S0034-4257(00)00205-4
        """
        return ee.Image(self.input_image) \
            .select(['blue', 'red', 'nir', 'swir1', 'swir2']) \
            .multiply([0.356, 0.130, 0.373, 0.085, 0.072]) \
            .reduce(ee.Reducer.sum()) \
            .subtract(0.0018) \
            .rename(['albedo'])

    @lazy_property
    def _lai(self):
        """Compute LAI using METRIC NDVI / LAI empirical equation

        Parameters
        ----------
        self.input_image : ee.Image
            "Prepped" Landsat image with standardized band names.

        Returns
        -------
        lai : ee.Image

        """
        return self._ndvi.pow(3).multiply(7.0).clamp(0, 6).rename(['lai'])

    @lazy_property
    def _emissivity(self):
        """METRIC narrowband emissivity"""
        lai = self._lai
        ndvi = self._ndvi
        # Initial values are for NDVI > 0 and LAI <= 3
        return lai.divide(300).add(0.97) \
            .where(ndvi.lte(0), 0.99) \
            .where(ndvi.gt(0).And(lai.gt(3)), 0.98)

    @lazy_property
    def _ndvi(self):
        """Compute NDVI

        Parameters
        ----------
        self.input_image : ee.Image
            "Prepped" Landsat image with standardized band names.

        Returns
        -------
        ndvi : ee.Image

        """
        return ee.Image(self.input_image).normalizedDifference(['nir', 'red']) \
            .rename(['ndvi'])


class LandsatTOA(Landsat):
    def __init__(self, raw_image):
        """Initialize a Landsat Collection 1 image

        Parameters
        ----------
        raw_image : ee.Image
            Landsat 5/7/8 Collection 1 TOA/SR image
            (i.e. from the "LANDSAT/X/C01/T1_XX" collection)

        """
        self.raw_image = ee.Image(raw_image)
        self._index = self.raw_image.get('system:index')
        self._time_start = self.raw_image.get('system:time_start')

        # Use the SPACECRAFT_ID property identify each Landsat type
        self._spacecraft_id = ee.String(self.raw_image.get('SPACECRAFT_ID'))

        input_bands = ee.Dictionary({
            'LANDSAT_5': ['B1', 'B2', 'B3', 'B4', 'B5', 'B7', 'B6', 'BQA'],
            'LANDSAT_7': ['B1', 'B2', 'B3', 'B4', 'B5', 'B7', 'B6_VCID_1', 'BQA'],
            'LANDSAT_8': ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B10', 'BQA']})
        # Rename bands to generic names
        output_bands = ['blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'lst',
                        'bqa']

        # Rename thermal band "k" coefficients to generic names
        input_k1 = ee.Dictionary({
            'LANDSAT_5': 'K1_CONSTANT_BAND_6',
            'LANDSAT_7': 'K1_CONSTANT_BAND_6_VCID_1',
            'LANDSAT_8': 'K1_CONSTANT_BAND_10'})
        input_k2 = ee.Dictionary({
            'LANDSAT_5': 'K2_CONSTANT_BAND_6',
            'LANDSAT_7': 'K2_CONSTANT_BAND_6_VCID_1',
            'LANDSAT_8': 'K2_CONSTANT_BAND_10'})
        output_k1 = self.raw_image.get(input_k1.get(self._spacecraft_id))
        output_k2 = self.raw_image.get(input_k2.get(self._spacecraft_id))

        self.input_image = ee.Image(self.raw_image) \
            .select(input_bands.get(self._spacecraft_id), output_bands) \
            .setMulti({
                'system:time_start': self._time_start,
                'system:index': self._index,
                'k1_constant': ee.Number(output_k1),
                'k2_constant': ee.Number(output_k2)})
        super()

    def prep(self):
        """Return an image with the bands/products needed to run EE DisALEXI

        Parameters
        ----------
        self.raw_image : ee.Image

        Returns
        -------
        prep_image : ee.Image

        """
        self.prep_image = ee.Image([
            self._albedo,
            self._cfmask,
            self._lai,
            self._lst,
            self._ndvi])
        self.prep_image = self.prep_image.setMulti({
            'system:time_start': self._time_start,
            'system:index': self._index,
        })
        return self.prep_image

    @lazy_property
    def _lst(self):
        """Compute emissivity corrected land surface temperature (LST)
        from brightness temperature.

        Note, the coefficients were derived from a small number of scenes in
        southern Idaho [Allen2007] and may not be appropriate for other areas.

        Parameters
        ----------
        self.input_image : ee.Image

        Returns
        -------
        lst : ee.Image

        References
        ----------
        .. [Allen2007a] R. Allen, M. Tasumi, R. Trezza (2007),
            Satellite-Based Energy Balance for Mapping Evapotranspiration with
            Internalized Calibration (METRIC) Model,
            Journal of Irrigation and Drainage Engineering, Vol 133(4),
            http://dx.doi.org/10.1061/(ASCE)0733-9437(2007)133:4(380)

        """
        # Get properties from image
        k1 = ee.Number(ee.Image(self.input_image).get('k1_constant'))
        k2 = ee.Number(ee.Image(self.input_image).get('k2_constant'))

        ts_brightness = ee.Image(self.input_image).select(['lst'])
        emissivity = self._emissivity

        # First back out radiance from brightness temperature
        # Then recalculate emissivity corrected Ts
        thermal_rad_toa = ts_brightness.expression(
            'k1 / (exp(k2 / ts_brightness) - 1)',
            {'ts_brightness': ts_brightness, 'k1': k1, 'k2': k2})

        # tnb = 0.866   # narrow band transmissivity of air
        # rp = 0.91     # path radiance
        # rsky = 1.32   # narrow band clear sky downward thermal radiation
        rc = thermal_rad_toa.expression(
            '((thermal_rad_toa - rp) / tnb) - ((1. - emiss) * rsky)',
            {
                'thermal_rad_toa': thermal_rad_toa,
                'emiss': emissivity,
                'rp': 0.91, 'tnb': 0.866, 'rsky': 1.32})
        lst = rc.expression(
            'k2 / log(emiss * k1 / rc + 1)',
            {'emiss': emissivity, 'rc': rc, 'k1': k1, 'k2': k2})

        return lst.rename(['lst'])

    @lazy_property
    def _cfmask(self):
        """Extract CFmask from Landsat Collection 1 BQA band

        https://landsat.usgs.gov/collectionqualityband

        Confidence values
        00 = "Not Determined" = Algorithm did not determine the status of this condition
        01 = "No" = Algorithm has low to no confidence that this condition exists
            (0-33 percent confidence)
        10 = "Maybe" = Algorithm has medium confidence that this condition exists
            (34-66 percent confidence)
        11 = "Yes" = Algorithm has high confidence that this condition exists
            (67-100 percent confidence


        Parameters
        ----------
        self.input_image : ee.Image

        Returns
        -------
        cfmask : ee.Image

        """
        bqa_image = ee.Image(self.input_image).select(['bqa'])

        def getQABits(bqa_image, start, end, newName):
            """
                From Tyler's function
                https://ee-api.appspot.com/#97ab9a8f694b28128a5a5ca2e2df7841
                """
            pattern = 0
            for i in range(start, end + 1):
                pattern += int(2 ** i)
            return bqa_image.select([0], [newName]) \
                .bitwise_and(pattern).right_shift(start)

        # Extract the various masks from the QA band
        fill_mask = getQABits(bqa_image, 0, 0, 'designated_fill')
        cloud_mask = getQABits(bqa_image, 5, 6, 'cloud_confidence').gte(2)
        shadow_mask = getQABits(bqa_image, 7, 8, 'shadow_confidence').gte(3)
        snow_mask = getQABits(bqa_image, 9, 10, 'snow_confidence').gte(3)
        # Landsat 8 only
        # cirrus_mask = getQABits(bqa_image, 11, 12, 'cirrus_confidence').gte(3)

        # Convert masks to old style Fmask values
        # 0 - Clear land
        # 1 - Clear water
        # 2 - Cloud shadow
        # 3 - Snow
        # 4 - Cloud
        return fill_mask \
            .add(shadow_mask.multiply(2)) \
            .add(snow_mask.multiply(3)) \
            .add(cloud_mask.multiply(4)) \
            .rename(['cfmask'])


class LandsatSR(Landsat):
    def __init__(self, raw_image):
        """Initialize a Landsat Collection 1 image

        Parameters
        ----------
        raw_image : ee.Image
            Landsat 5/7/8 Collection 1 SR image
            (i.e. from the "LANDSAT/X/C01/T1_XX" collection)

        """
        scalars = [0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.1, 1]

        self.raw_image = ee.Image(raw_image)
        self._index = self.raw_image.get('system:index')
        self._time_start = self.raw_image.get('system:time_start')

        # Use the SATELLITE property to identify each Landsat type
        self._spacecraft_id = ee.String(self.raw_image.get('SATELLITE'))

        input_bands = ee.Dictionary({
            'LANDSAT_5': ['B1', 'B2', 'B3', 'B4', 'B5', 'B7', 'B6', 'pixel_qa'],
            'LANDSAT_7': ['B1', 'B2', 'B3', 'B4', 'B5', 'B7', 'B6',
                          'pixel_qa'],
            'LANDSAT_8': ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B10', 'pixel_qa']})
        # Rename bands to generic names
        output_bands = [
            'blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'lst', 'pixel_qa']

        # TODO: Follow up with Simon about adding K1/K2 to SR collection
        k1 = ee.Dictionary({
            'LANDSAT_5': 607.76, 'LANDSAT_7': 666.09, 'LANDSAT_8': 774.8853})
        k2 = ee.Dictionary({
            'LANDSAT_5': 1260.56, 'LANDSAT_7': 1282.71, 'LANDSAT_8': 1321.0789})

        self.input_image = ee.Image(self.raw_image) \
            .select(input_bands.get(self._spacecraft_id), output_bands) \
            .multiply(scalars) \
            .setMulti({
                'system:time_start': self._time_start,
                'system:index': self._index,
                'k1_constant': ee.Number(k1.get(self._spacecraft_id)),
                'k2_constant': ee.Number(k2.get(self._spacecraft_id))})
        super()

    def prep(self):
        """Return an image with the bands/products needed to run EE DisALEXI

        Parameters
        ----------
        self.raw_image : ee.Image

        Returns
        -------
        prep_image : ee.Image

        """

        self.prep_image = ee.Image([
            self._albedo,
            self._cfmask,
            self._lai,
            self._lst,
            self._ndvi])

        self.prep_image = ee.Image(self.prep_image.setMulti({
            'system:time_start': self._time_start,
            'system:index': self._index,
        }))
        return self.prep_image

        # The cloud mask could be applied here
        # mask_img = self.prep_image.select('Mask').eq(0)
        # return self.prep_image.updateMask(mask_img)

    @lazy_property
    def _cfmask(self):
        """Extract CFmask from Landsat Collection 1 pixel_qa band

        https://landsat.usgs.gov/collectionqualityband

        Confidence values
        00 = "Not Determined" = Algorithm did not determine the status of this condition
        01 = "No" = Algorithm has low to no confidence that this condition exists
            (0-33 percent confidence)
        10 = "Maybe" = Algorithm has medium confidence that this condition exists
            (34-66 percent confidence)
        11 = "Yes" = Algorithm has high confidence that this condition exists
            (67-100 percent confidence


        Parameters
        ----------
        self.input_image : ee.Image

        Returns
        -------
        cfmask : ee.Image

        """
        bqa_image = ee.Image(self.input_image).select(['pixel_qa'])

        def getQABits(bqa_image, start, end, newName):
            """
                From Tyler's function
                https://ee-api.appspot.com/#97ab9a8f694b28128a5a5ca2e2df7841
                """
            pattern = 0
            for i in range(start, end + 1):
                pattern += int(2 ** i)
            return bqa_image.select([0], [newName]) \
                .bitwise_and(pattern).right_shift(start)

        # Extract the various masks from the QA band
        fill_mask = getQABits(bqa_image, 0, 0, 'fill')
        cloud_mask = getQABits(bqa_image, 6, 7, 'cloud_confidence').gte(2)
        shadow_mask = getQABits(bqa_image, 3, 3, 'shadow')
        snow_mask = getQABits(bqa_image, 4, 4, 'snow')
        # Landsat 8 only
        # cirrus_mask = getQABits(bqa_image, 11, 12, 'cirrus_confidence').gte(3)

        # Convert masks to old style Fmask values
        # 0 - Clear land
        # 1 - Clear water
        # 2 - Cloud shadow
        # 3 - Snow
        # 4 - Cloud
        return fill_mask \
            .add(shadow_mask.multiply(2)) \
            .add(snow_mask.multiply(3)) \
            .add(cloud_mask.multiply(4)) \
            .rename(['cfmask'])

    @lazy_property
    def _lst(self):
        """Compute emissivity corrected land surface temperature (LST)
        from brightness temperature.

        Note, the coefficients were derived from a small number of scenes in
        southern Idaho [Allen2007] and may not be appropriate for other areas.

        Parameters
        ----------
        self.input_image : ee.Image

        Returns
        -------
        lst : ee.Image

        References
        ----------
        .. [Allen2007a] R. Allen, M. Tasumi, R. Trezza (2007),
            Satellite-Based Energy Balance for Mapping Evapotranspiration with
            Internalized Calibration (METRIC) Model,
            Journal of Irrigation and Drainage Engineering, Vol 133(4),
            http://dx.doi.org/10.1061/(ASCE)0733-9437(2007)133:4(380)

        """
        # Get properties from image
        k1 = ee.Number(ee.Image(self.input_image).get('k1_constant'))
        k2 = ee.Number(ee.Image(self.input_image).get('k2_constant'))

        ts_brightness = ee.Image(self.input_image).select(['lst'])
        emissivity = self._emissivity

        # First back out radiance from brightness temperature
        # Then recalculate emissivity corrected Ts
        thermal_rad_toa = ts_brightness.expression(
            'k1 / (exp(k2 / ts_brightness) - 1)',
            {'ts_brightness': ts_brightness, 'k1': k1, 'k2': k2})

        # tnb = 0.866   # narrow band transmissivity of air
        # rp = 0.91     # path radiance
        # rsky = 1.32   # narrow band clear sky downward thermal radiation
        rc = thermal_rad_toa.expression(
            '((thermal_rad_toa - rp) / tnb) - ((1. - emiss) * rsky)',
            {
                'thermal_rad_toa': thermal_rad_toa,
                'emiss': emissivity,
                'rp': 0.91, 'tnb': 0.866, 'rsky': 1.32})
        lst = rc.expression(
            'k2 / log(emiss * k1 / rc + 1)',
            {'emiss': emissivity, 'rc': rc, 'k1': k1, 'k2': k2})

        return lst.rename(['lst'])

    # @lazy_property
    # def _lst(self):
    #     """just return the Brightness Temperature (BT) as land surface temperature (LST) for now
    #
    #     Parameters
    #     ----------
    #     self.input_image : ee.Image
    #
    #     Returns
    #     -------
    #     lst : ee.Image
    #
    #     """
    #     ts_brightness = ee.Image(self.input_image).select(['lst'])
    #     return ts_brightness.rename(['lst'])
