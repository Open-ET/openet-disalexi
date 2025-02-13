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
        Using a sum reducer was returning an unbounded image,
            but this could probably be fixed with a setDefaultProjection() call
        Clamping the reflectance bands instead of the output
            since the saturated bands have a significant jump from 0.6-1.0 up to 1.6
            but also adding a max(0) on the output to avoid any possibility
            of a negative output

        References
        ----------
        .. [Liang2001] Shunlin Liang (2001),
            Narrowband to broadband conversions of land surface albedo -
            I Algorithms, Remote Sensing of Environment,
            Volume 76, Issue2, Pages 213-238,
            http://doi.org/10.1016/S0034-4257(00)00205-4

        """
        albedo = (
            self.input_image
            .select(['blue', 'red', 'nir', 'swir1', 'swir2'])
            .clamp(0, 1)
            .multiply([0.356, 0.130, 0.373, 0.085, 0.072])
        )
        return (
            albedo.select([0]).add(albedo.select([1])).add(albedo.select([2]))
            .add(albedo.select([3])).add(albedo.select([4])).subtract(0.0018)
            .max(0)
            .rename(['albedo'])
        )

    # LAI is being read from a source image collection
    #   but keeping this method since it is currently used in emissivity calculation
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
        # Initial values are for NDVI >= 0 and LAI <= 3
        # CGM - This function originally set water as .lte(0)
        return (
            self._lai.divide(300).add(0.97)
            .where(self._ndvi.lt(0), 0.99)
            .where(self._ndvi.gt(0).And(self._lai.gt(3)), 0.98)
            .rename(['emissivity'])
        )

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
        # DEADBEEF - Original NDVI calculation that masks negative values
        # return self.input_image.normalizedDifference(['nir', 'red']).rename(['ndvi'])

        # Force the input values to be at greater than or equal to zero
        #   since C02 surface reflectance values can be negative
        #   but the normalizedDifference function will return nodata
        ndvi_img = self.input_image.select(['nir', 'red']).max(0).normalizedDifference(['nir', 'red'])

        # Assume that very high reflectance values are unreliable for computing the index
        #   and set the output value to 0
        # Threshold value could be set lower, but for now only trying to catch saturated pixels
        b1 = self.input_image.select(['nir'])
        b2 = self.input_image.select(['red'])
        ndvi_img = ndvi_img.where(b1.gte(1).Or(b2.gte(1)), 0)

        # Assume that low reflectance values are unreliable for computing the index
        # If both reflectance values are below the threshold,
        #   and if the pixel is flagged as water, set the output to -0.1 (should this be -1?)
        #   otherwise set the output to 0
        ndvi_img = ndvi_img.where(b1.lt(0.01).And(b2.lt(0.01)), 0)

        # Including the global surface water maximum extent to help remove shadows that
        #   are misclassified as water
        # The flag is needed so that the image can be bypassed during testing with constant images
        if self.gsw_extent_flag:
            qa_water_mask = self.input_image.select(['QA_PIXEL']).rightShift(7).bitwiseAnd(1).neq(0)
            gsw_mask = ee.Image('JRC/GSW1_4/GlobalSurfaceWater').select(['max_extent']).gte(1)
            qa_water_mask = qa_water_mask.And(gsw_mask)
            ndvi_img = ndvi_img.where(b1.lt(0.01).And(b2.lt(0.01)).And(qa_water_mask), -0.1)

        # Should there be an additional check for if either value was negative?
        # ndvi_img = ndvi_img.where(b1.lt(0).Or(b2.lt(0)), 0)

        return ndvi_img.rename(['ndvi'])


class Landsat_C02_L2(Landsat):
    def __init__(self, raw_image, gsw_extent_flag=True):
        """Initialize a Landsat Collection 2 SR image

        Parameters
        ----------
        raw_image : ee.Image, str
            Landsat 5/7/8/9 Collection 2 SR image or image ID
            (i.e. from the "LANDSAT/X/C02/T1_L2" collection)

        """
        scalars_multi = [0.0000275, 0.0000275, 0.0000275, 0.0000275, 0.0000275, 0.0000275, 0.00341802, 1]
        scalars_add = [-0.2, -0.2, -0.2, -0.2, -0.2, -0.2, 149.0, 0]

        self.raw_image = ee.Image(raw_image)
        self._id = self.raw_image.get('system:id')
        self._index = self.raw_image.get('system:index')
        self._time_start = self.raw_image.get('system:time_start')

        # Use the SATELLITE property to identify each Landsat type
        self._spacecraft_id = ee.String(self.raw_image.get('SPACECRAFT_ID'))

        input_bands = ee.Dictionary({
            'LANDSAT_4': ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7', 'ST_B6', 'QA_PIXEL'],
            'LANDSAT_5': ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7', 'ST_B6', 'QA_PIXEL'],
            'LANDSAT_7': ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7', 'ST_B6', 'QA_PIXEL'],
            'LANDSAT_8': ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7', 'ST_B10', 'QA_PIXEL'],
            'LANDSAT_9': ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7', 'ST_B10', 'QA_PIXEL'],
        })
        # Rename bands to generic names
        output_bands = ['blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'lst', 'QA_PIXEL']

        self.input_image = (
            ee.Image(self.raw_image)
            .select(input_bands.get(self._spacecraft_id), output_bands)
            .multiply(scalars_multi)
            .add(scalars_add)
            .set({
                'system:time_start': self._time_start,
                'system:index': self._index,
                'system:id': self._id,
                'SATELLITE': self._spacecraft_id,
                'SPACECRAFT_ID': self._spacecraft_id,
            })
        )
        self.gsw_extent_flag = gsw_extent_flag
        super()

    def prep(self):
        """Return an image with the bands/products needed to run EE DisALEXI

        Parameters
        ----------

        Returns
        -------
        prep_image : ee.Image

        """

        # TIR/LST is being read from a source image collection and does not need to be passed in
        self.prep_image = ee.Image([
            self._albedo,
            self._cfmask,
            self._ndvi,
            # self._lai,
            # self._lst,
        ])
        self.prep_image = self.prep_image.set({
            'system:time_start': self._time_start,
            'system:index': self._index,
            'system:id': self._id,
        })
        self.prep_image = ee.Image(self.prep_image.copyProperties(self.input_image))

        return self.prep_image

    @lazy_property
    def _cfmask(
            self,
            cirrus_flag=False,
            dilate_flag=False,
            shadow_flag=True,
            snow_flag=True,
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
        # The following line could be added to the mask_img call
        #     to include the cloud confidence bits
        # .Or(qa_img.rightShift(8).bitwiseAnd(3).gte(cloud_confidence))

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
