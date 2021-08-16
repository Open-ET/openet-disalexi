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
        """Total shortwave broadband albedo following [Liang2001]

        Parameters
        ----------
        self.input_image : ee.Image
            "Prepped" Landsat image with standardized band names.

        Returns
        -------
        albedo : ee.Image

        Notes
        -----
        The Python DisALEXI code had the following line and comment:
            "bands = [1, 3, 4, 5, 7]  # dont use blue"
        IDL code and [Liang2001] indicate that the green band is not used.
        Coefficients were derived for Landsat 7 ETM+, but were found to be
            "suitable" to Landsat 4/5 TM also.

        References
        ----------
        .. [Liang2001] Shunlin Liang (2001),
            Narrowband to broadband conversions of land surface albedo -
            I Algorithms, Remote Sensing of Environment,
            Volume 76, Issue2, Pages 213-238,
            http://doi.org/10.1016/S0034-4257(00)00205-4

        """
        albedo = self.input_image\
            .select(['blue', 'red', 'nir', 'swir1', 'swir2'])\
            .multiply([0.356, 0.130, 0.373, 0.085, 0.072])
        return albedo.select([0])\
            .add(albedo.select([1])).add(albedo.select([2]))\
            .add(albedo.select([3])).add(albedo.select([4]))\
            .subtract(0.0018)\
            .rename(['albedo'])

        # # Using a sum reducer was returning an unbounded image
        # return ee.Image(self.input_image)\
        #     .select(['blue', 'red', 'nir', 'swir1', 'swir2'])\
        #     .multiply([0.356, 0.130, 0.373, 0.085, 0.072])\
        #     .reduce(ee.Reducer.sum())\
        #     .subtract(0.0018)\
        #     .rename(['albedo'])

    # DEADBEEF - LAI is being read from a source image collection
    #   Leaving this method since it is currently used in emissivity calculation
    @lazy_property
    def _lai(self):
        """Leaf Area Index (LAI) computed from METRIC NDVI / LAI equation

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
        """Normalized difference vegetation index

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


class Landsat_C01_SR(Landsat):
    def __init__(self, raw_image):
        """Initialize a Landsat Collection 1 SR image

        Parameters
        ----------
        raw_image : ee.Image
            Landsat 5/7/8 Collection 1 SR image
            (i.e. from the "LANDSAT/X/C01/T1_SR" collection)

        """
        scalars = [0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.1, 1]

        self.raw_image = ee.Image(raw_image)
        self._id = self.raw_image.get('system:id')
        self._index = self.raw_image.get('system:index')
        self._time_start = self.raw_image.get('system:time_start')

        # Use the SATELLITE property to identify each Landsat type
        self._spacecraft_id = ee.String(self.raw_image.get('SATELLITE'))

        # Rename bands to generic names
        input_bands = ee.Dictionary({
            'LANDSAT_4': ['B1', 'B2', 'B3', 'B4', 'B5', 'B7', 'B6', 'pixel_qa'],
            'LANDSAT_5': ['B1', 'B2', 'B3', 'B4', 'B5', 'B7', 'B6', 'pixel_qa'],
            'LANDSAT_7': ['B1', 'B2', 'B3', 'B4', 'B5', 'B7', 'B6', 'pixel_qa'],
            'LANDSAT_8': ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B10', 'pixel_qa']})
        output_bands = [
            'blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'tir', 'pixel_qa']

        self.input_image = ee.Image(self.raw_image) \
            .select(input_bands.get(self._spacecraft_id), output_bands) \
            .multiply(scalars) \
            .set({
                'system:time_start': self._time_start,
                'system:index': self._index,
                'system:id': self._id,
                'SATELLITE': self._spacecraft_id,
                'SPACECRAFT_ID': self._spacecraft_id,
            })
        super()

    def prep(self):
        """Return an image with the bands/products needed to run EE DisALEXI

        Parameters
        ----------

        Returns
        -------
        prep_image : ee.Image

        """

        # CGM - TIR/LST is being read from a source image collection
        self.prep_image = ee.Image([
            self._albedo,
            self._cfmask,
            # self._lai,
            # self._lst,
            self._ndvi,
        ])
        self.prep_image = self.prep_image.set({
            'system:time_start': self._time_start,
            'system:index': self._index,
            'system:id': self._id,
        })
        self.prep_image = ee.Image(self.prep_image.copyProperties(self.input_image))

        return self.prep_image

    @lazy_property
    def _cfmask(self):
        """Extract CFmask like image from Landsat Collection 1 SR pixel_qa band
        Parameters
        ----------
        self.input_image : ee.Image
        Returns
        -------
        cfmask : ee.Image
        Notes
        -----
        https://landsat.usgs.gov/collectionqualityband
        Confidence values
        00 = "Not Determined" = Algorithm did not determine the status of this condition
        01 = "No" = Algorithm has low to no confidence that this condition exists
            (0-33 percent confidence)
        10 = "Maybe" = Algorithm has medium confidence that this condition exists
            (34-66 percent confidence)
        11 = "Yes" = Algorithm has high confidence that this condition exists
            (67-100 percent confidence
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
        #cloud_mask = getQABits(bqa_image, 6, 7, 'cloud_confidence').gte(2)
        #Yun modified the cloud mask. Instead of using the confidence level, directly using cloud bit.
        #This is a relatively stricter cloud mask, but matching with our inhouse standard.
        cloud_mask = getQABits(bqa_image, 5, 5,'cloud_confidence')
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


class Landsat_C02_SR(Landsat):
    def __init__(self, raw_image):
        """Initialize a Landsat Collection 2 SR image

        Parameters
        ----------
        raw_image : ee.Image
            Landsat 5/7/8 Collection 2 SR image
            (i.e. from the "LANDSAT/X/C02/T1_L2" collection)

        """
        scalars_multi = [
            0.0000275, 0.0000275, 0.0000275, 0.0000275, 0.0000275, 0.0000275,
            0.00341802, 1]
        scalars_add = [-0.2, -0.2, -0.2, -0.2, -0.2, -0.2, 149.0, 1]

        self.raw_image = ee.Image(raw_image)
        self._id = self.raw_image.get('system:id')
        self._index = self.raw_image.get('system:index')
        self._time_start = self.raw_image.get('system:time_start')

        # Use the SATELLITE property to identify each Landsat type
        self._spacecraft_id = ee.String(self.raw_image.get('SPACECRAFT_ID'))

        input_bands = ee.Dictionary({
            'LANDSAT_4': ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7',
                          'ST_B6', 'QA_PIXEL'],
            'LANDSAT_5': ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7',
                          'ST_B6', 'QA_PIXEL'],
            'LANDSAT_7': ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7',
                          'ST_B6', 'QA_PIXEL'],
            'LANDSAT_8': ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7',
                          'ST_B10', 'QA_PIXEL']})
        # Rename bands to generic names
        output_bands = [
            'blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'lst', 'QA_PIXEL']

        self.input_image = ee.Image(self.raw_image) \
            .select(input_bands.get(self._spacecraft_id), output_bands) \
            .multiply(scalars_multi) \
            .add(scalars_add)\
            .set({
                'system:time_start': self._time_start,
                'system:index': self._index,
                'system:id': self._id,
                'SATELLITE': self._spacecraft_id,
                'SPACECRAFT_ID': self._spacecraft_id,
            })
        super()

    def prep(self):
        """Return an image with the bands/products needed to run EE DisALEXI

        Parameters
        ----------

        Returns
        -------
        prep_image : ee.Image

        """

        # CGM - TIR/LST is being read from a source image collection
        self.prep_image = ee.Image([
            self._albedo,
            self._cfmask,
            # self._lai,
            # self._lst,
            self._ndvi,
        ])
        self.prep_image = self.prep_image.set({
            'system:time_start': self._time_start,
            'system:index': self._index,
            'system:id': self._id,
        })
        self.prep_image = ee.Image(self.prep_image.copyProperties(self.input_image))

        return self.prep_image

    @lazy_property
    def _cfmask(self, cirrus_flag=False, dilate_flag=False,
                shadow_flag=True, snow_flag=True,
                ):
        """Extract cloud mask from the Landsat Collection 2 SR QA_PIXEL band

        Parameters
        ----------
        img : ee.Image
            Image from a Landsat Collection 2 SR image collection with a QA_PIXEL
            band (e.g. LANDSAT/LC08/C02/T1_L2).
        cirrus_flag : bool
            If true, mask cirrus pixels (the default is False).
            Note, cirrus bits are only set for Landsat 8 (OLI) images.
        dilate_flag : bool
            If true, mask dilated cloud pixels (the default is False).
        shadow_flag : bool
            If true, mask shadow pixels (the default is True).
        snow_flag : bool
            If true, mask snow pixels (the default is True).

        Returns
        -------
        ee.Image

        Notes
        -----
        Output image is structured to be applied directly with updateMask()
            i.e. 0 is cloud/masked, 1 is clear/unmasked
        Assuming Cloud must be set to check Cloud Confidence
        Bits
            0: Fill
                0 for image data
                1 for fill data
            1: Dilated Cloud
                0 for cloud is not dilated or no cloud
                1 for cloud dilation
            2: Cirrus
                0 for no confidence level set or low confidence
                1 for high confidence cirrus
            3: Cloud
                0 for cloud confidence is not high
                1 for high confidence cloud
            4: Cloud Shadow
                0 for Cloud Shadow Confidence is not high
                1 for high confidence cloud shadow
            5: Snow
                0 for Snow/Ice Confidence is not high
                1 for high confidence snow cover
            6: Clear
                0 if Cloud or Dilated Cloud bits are set
                1 if Cloud and Dilated Cloud bits are not set
            7: Water
                0 for land or cloud
                1 for water
            8-9: Cloud Confidence
            10-11: Cloud Shadow Confidence
            12-13: Snow/Ice Confidence
            14-15: Cirrus Confidence
        Confidence values
            00: "No confidence level set"
            01: "Low confidence"
            10: "Medium confidence" (for Cloud Confidence only, otherwise "Reserved")
            11: "High confidence"

        References
        ----------
        https://prd-wret.s3.us-west-2.amazonaws.com/assets/palladium/production/atoms/files/LSDS-1328_Landsat8-9-OLI-TIRS-C2-L2-DFCB-v6.pdf

        """
        qa_img = ee.Image(self.input_image).select(['QA_PIXEL'])
        cloud_mask = qa_img.rightShift(3).bitwiseAnd(1).neq(0).multiply(4)
        #     .And(qa_img.rightShift(6).bitwiseAnd(3).gte(cloud_confidence))
        if cirrus_flag:
            cirrus_mask = qa_img.rightShift(2).bitwiseAnd(1).neq(0)
            cloud_mask = cloud_mask.add(cirrus_mask.multiply(5))
        if dilate_flag:
            dilate_mask = qa_img.rightShift(1).bitwiseAnd(1).neq(0)
            cloud_mask = cloud_mask.add(dilate_mask.multiply(6))
        if shadow_flag:
            shadow_mask = qa_img.rightShift(4).bitwiseAnd(1).neq(0)
            cloud_mask = cloud_mask.add(shadow_mask.multiply(2))
        if snow_flag:
            snow_mask = qa_img.rightShift(5).bitwiseAnd(1).neq(0)
            cloud_mask = cloud_mask.add(snow_mask.multiply(3))

        return cloud_mask.rename(['cfmask'])
