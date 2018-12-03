import logging

import ee
import pytest

import openet.disalexi.disalexi as disalexi
import openet.disalexi.utils as utils

# AmeriFlux sites adjusted to nearest Landsat cell centroid
ne1_xy = [-96.47672812080845, 41.16506126041818]
ne2_xy = [-96.46994024736414, 41.16491226772292]
ne3_xy = [-96.43968912903934, 41.17964494123755]

# Initialize a "prepped" Landsat image
asset_ws = 'users/cgmorton/disalexi/LC08_028031_20140708/'

img_date_str = '2014-07-08'
img_hour = 17
img_doy = 189
img_time_start = 1404839150550
# img_time_start = ee.Date('2014-05-21T17:00', 'GMT').millis()
img_hour_start = ee.Date('2014-07-08T17:00', 'GMT').millis()
img_date_start = ee.Date('2014-07-08', 'GMT').millis()
# img_hour = ee.Date('2014-05-21T17:05', 'GMT').millis()
img_date = ee.Date('2014-07-08', 'GMT').millis()
# Time used in IDL is chours and fractional minutes (no seconds)
img_time = 17 + 5.0/60

test_img = ee.Image([
    ee.Image(asset_ws + 'albedo'),
    ee.Image(asset_ws + 'cfmask'),
    ee.Image(asset_ws + 'lai'),
    ee.Image(asset_ws + 'lst'),
    ee.Image(asset_ws + 'ndvi')
])
test_img = ee.Image(test_img \
    .rename(['albedo', 'cfmask', 'lai', 'lst', 'ndvi']) \
    .setMulti({'system:time_start': img_time_start}))


def test_Image_init_dates():
    d_obj = disalexi.Image(test_img)
    assert d_obj.date.format('yyyy-MM-dd').getInfo() == img_date_str
    assert d_obj.doy.getInfo() == img_doy
    assert int(d_obj.hour.getInfo()) == img_hour
    assert float(d_obj.time.getInfo()) == img_time


def test_Image_init_missing_landcover():
    """Test that only setting a land cover type raises a ValueError exception"""
    with pytest.raises(ValueError) as e_info:
        d_obj = disalexi.Image(test_img, lc_type='NLCD')


def test_Image_init_missing_lc_type():
    """Test that only setting a land cover image raises a ValueError exception"""
    with pytest.raises(ValueError) as e_info:
        d_obj = disalexi.Image(
            test_img, landcover=ee.Image('USGS/NLCD/NLCD2011'))


@pytest.mark.parametrize(
    'xy,et',
    [
        # CONUS ALEXI ET values
        [ne1_xy, 7.287301],
        [ne2_xy, 7.287301],
        [ne3_xy, 7.477266],
    ]
)
def test_Image_set_alexi_et_vars_defaults(xy, et, tol=1E-6):
    d_obj = disalexi.Image(test_img)
    d_obj._set_alexi_et_vars()
    assert abs(utils.image_value(
        ee.Image(d_obj.alexi_et), xy=xy, scale=0.1)['alexi_et'] - et) <= tol


@pytest.mark.parametrize(
    'xy,et',
    [
        [ne1_xy, 5.995176],
        [ne2_xy, 5.995176],
        [ne3_xy, 6.320339],
    ]
)
def test_Image_set_alexi_et_vars_assets(xy, et, tol=1E-6):
    """

    Don't use scale parameter in image_value since ALEXI ET assets are already
        resampled to the Landsat grid.
    Add separate check that band name is set correctly?
    """
    d_obj = disalexi.Image(test_img)
    d_obj.et_coll = ee.ImageCollection([
        ee.Image(asset_ws + 'alexiET') \
            .setMulti({'system:time_start': img_date_start})])
    d_obj.et_transform = [0.04, 0, -96.442, 0, -0.04, 41.297]
    d_obj._set_alexi_et_vars()
    assert abs(utils.image_value(
        ee.Image(d_obj.alexi_et), xy=xy)['alexi_et'] - et) <= tol


@pytest.mark.parametrize(
    'xy,elevation,pressure',
    [
        [ne1_xy, 350, 97.2306250000000],
        [ne2_xy, 350, 97.2306250000000],
        [ne3_xy, 350, 97.2306250000000],
    ]
)
def test_Image_set_elevation_vars(xy, elevation, pressure, tol=1E-6):
    """"""
    d_obj = disalexi.Image(test_img)
    d_obj.elevation = ee.Image.constant(elevation)
    d_obj._set_elevation_vars()
    assert abs(utils.image_value(ee.Image(d_obj.pressure), xy=xy)['pressure'] -
               pressure) <= tol


def test_Image_set_landcover_vars_invalid_lc_type():
    """Test that setting an invalid lc_type value raises a KeyError exception"""
    d_obj = disalexi.Image(
        test_img, landcover=ee.Image('USGS/NLCD/NLCD2011'), lc_type='DEADBEEF')
    # print(d_obj._set_landcover_vars())
    with pytest.raises(KeyError) as e_info:
        d_obj._set_landcover_vars()


def test_Image_set_landcover_vars_default(tol=1E-6):
    """Test default land cover image and type

    It might make more sense to just test that the value at the test pixel
    is 82 (for NLCD) for 10 (for GLC30)"""
    d_obj = disalexi.Image(test_img)
    d_obj._set_landcover_vars()
    assert utils.image_value(ee.Image(d_obj.aleafv))['aleafv'] == 0.83
    assert utils.image_value(ee.Image(d_obj.aleafn))['aleafn'] == 0.35
    assert utils.image_value(ee.Image(d_obj.aleafl))['aleafl'] == 0.95
    assert utils.image_value(ee.Image(d_obj.adeadv))['adeadv'] == 0.49
    assert utils.image_value(ee.Image(d_obj.adeadn))['adeadn'] == 0.13
    assert utils.image_value(ee.Image(d_obj.adeadl))['adeadl'] == 0.95
    assert utils.image_value(ee.Image(d_obj.leaf_width))['xl'] == 0.05
    assert utils.image_value(ee.Image(d_obj.clump))['omega'] == 0.83
    # assert abs(utils.image_value(ee.Image(d_obj.hc))['hc'] -
    #            0.44531310527793) <= tol


def test_Image_set_landcover_vars_init_asset(tol=1E-6):
    """Test setting the land cover image and type as the object is initialized"""
    d_obj = disalexi.Image(
        test_img, lc_type='NLCD',
        landcover=ee.Image(asset_ws + 'landcover'))
    d_obj._set_landcover_vars()
    assert utils.image_value(ee.Image(d_obj.aleafv))['aleafv'] == 0.83
    # assert utils.image_value(ee.Image(d_obj.aleafn))['aleafn'] == 0.35
    # assert utils.image_value(ee.Image(d_obj.aleafl))['aleafl'] == 0.95
    # assert utils.image_value(ee.Image(d_obj.adeadv))['adeadv'] == 0.49
    # assert utils.image_value(ee.Image(d_obj.adeadn))['adeadn'] == 0.13
    # assert utils.image_value(ee.Image(d_obj.adeadl))['adeadl'] == 0.95
    # assert utils.image_value(ee.Image(d_obj.leaf_width))['xl'] == 0.05
    # assert utils.image_value(ee.Image(d_obj.clump))['omega'] == 0.83
    # assert abs(utils.image_value(ee.Image(d_obj.hc))['hc'] -
    #            0.410955099827) <= tol


def test_Image_set_landcover_vars_set_asset(tol=1E-6):
    """Test setting the land cover image and type directly on the object"""
    d_obj = disalexi.Image(test_img)
    d_obj.landcover = ee.Image(asset_ws + 'landcover'),
    d_obj.lc_type = 'NLCD'
    d_obj._set_landcover_vars()
    assert utils.image_value(ee.Image(d_obj.aleafv))['aleafv'] == 0.83
    # assert utils.image_value(ee.Image(d_obj.aleafn))['aleafn'] == 0.35
    # assert utils.image_value(ee.Image(d_obj.aleafl))['aleafl'] == 0.95
    # assert utils.image_value(ee.Image(d_obj.adeadv))['adeadv'] == 0.49
    # assert utils.image_value(ee.Image(d_obj.adeadn))['adeadn'] == 0.13
    # assert utils.image_value(ee.Image(d_obj.adeadl))['adeadl'] == 0.95
    # assert utils.image_value(ee.Image(d_obj.leaf_width))['xl'] == 0.05
    # assert utils.image_value(ee.Image(d_obj.clump))['omega'] == 0.83
    # assert abs(utils.image_value(ee.Image(d_obj.hc))['hc'] -
    #            0.410955099827) <= tol


@pytest.mark.parametrize(
    'xy,interp,rs1,rs24',
    [
        [ne1_xy, True, 880.75 + (35.8425 / 60) * (956.25 - 880.75), 8506.97168],
        [ne2_xy, True, 880.75 + (35.8425 / 60) * (956.25 - 880.75), 8506.97168],
        [ne3_xy, True, 880.75 + (35.8425 / 60) * (956.25 - 880.75), 8506.97168],
        [ne1_xy, False, 956.25, 8506.97168],
        [ne2_xy, False, 956.25, 8506.97168],
        [ne3_xy, False, 956.25, 8506.97168],
    ]
)
def test_Image_set_solar_vars_defaults(xy, interp, rs1, rs24, tol=1E-4):
    """Test that the default MERRA2 Rs values are returned"""
    d_obj = disalexi.Image(test_img)
    d_obj._set_solar_vars(interpolate_flag=interp)
    assert abs(utils.image_value(
        ee.Image(d_obj.rs1), xy, scale=0.1)['rs'] - rs1) <= tol
    assert abs(utils.image_value(
        ee.Image(d_obj.rs24), xy, scale=0.1)['rs'] - rs24) <= tol


@pytest.mark.parametrize(
    'xy,rs1,rs24',
    [
        [ne1_xy, 917.87845865885413, 8559.066406],
        [ne2_xy, 917.87921142578125, 8558.934570],
        [ne3_xy, 917.72406005859375, 8557.359375],
    ]
)
def test_Image_set_solar_vars_assets_no_interp(xy, rs1, rs24, tol=1E-4):
    """Test that the default MERRA2 Rs values are returned"""
    d_obj = disalexi.Image(test_img)
    d_obj.rs_hourly_coll = ee.ImageCollection([
        ee.Image(asset_ws + 'Insol1')
            .setMulti({'system:time_start': img_hour_start})])
    d_obj.rs_daily_coll = ee.ImageCollection([
        ee.Image(asset_ws + 'Insol24') \
            .setMulti({'system:time_start': img_date_start})])
    d_obj._set_solar_vars(interpolate_flag=False)
    assert abs(utils.image_value(
        ee.Image(d_obj.rs1), xy)['rs'] - rs1) <= tol
    assert abs(utils.image_value(
        ee.Image(d_obj.rs24), xy)['rs'] - rs24) <= tol


@pytest.mark.parametrize(
    'xy,rs1,rs24',
    [
        [ne1_xy, 917.87845865885413, 8559.066406],
        [ne2_xy, 917.87921142578125, 8558.934570],
        [ne3_xy, 917.72406005859375, 8557.359375],
    ]
)
def test_Image_set_solar_vars_assets_interp(xy, rs1, rs24, tol=1E-4):
    """Test that the default MERRA2 Rs values are returned"""
    d_obj = disalexi.Image(test_img)
    d_obj.rs_hourly_coll = ee.ImageCollection([
        ee.Image(asset_ws + 'Insol1')
            .setMulti({'system:time_start': img_hour_start.subtract(3600000)}),
        ee.Image(asset_ws + 'Insol1')
            .setMulti({'system:time_start': img_hour_start}),
        ee.Image(asset_ws + 'Insol1')
            .setMulti({'system:time_start': img_hour_start.add(3600000)})
    ])
    d_obj.rs_daily_coll = ee.ImageCollection([
        ee.Image(asset_ws + 'Insol24') \
            .setMulti({'system:time_start': img_date_start})])
    d_obj._set_solar_vars(interpolate_flag=True)
    assert abs(utils.image_value(
        ee.Image(d_obj.rs1), xy)['rs'] - rs1) <= tol
    assert abs(utils.image_value(
        ee.Image(d_obj.rs24), xy)['rs'] - rs24) <= tol


@pytest.mark.parametrize(
    'xy,t_rise,t_end',
    [
        [ne1_xy, 11.02448526443880, 26.01087501850882],
        [ne2_xy, 11.02404074400912, 26.01041448914592],
        [ne3_xy, 11.02123229380177, 26.00918945690997],
    ]
)
def test_Image_set_time_vars_defaults(xy, t_rise, t_end, tol=1E-8):
    """Test setting the land cover image and type directly on the object

    High NDVI test point values
    CGM - Should probably switch this to a constant image test
    """
    d_obj = disalexi.Image(test_img)
    d_obj._set_time_vars()
    assert abs(utils.image_value(
        ee.Image(d_obj.t_rise), xy)['t_rise'] - t_rise) <= tol
    assert abs(utils.image_value(
        ee.Image(d_obj.t_end), xy)['t_end'] - t_end) <= tol


def test_Image_set_weather_vars_defaults(tol=0.01):
    d_obj = disalexi.Image(test_img)
    d_obj._set_weather_vars()
    assert abs(utils.image_value(
        ee.Image(d_obj.windspeed))['windspeed'] - 4.12) <= tol


def test_Image_set_weather_var_assets(tol=0.01):
    d_obj = disalexi.Image(test_img)
    d_obj.windspeed_coll = ee.ImageCollection([
        ee.Image([ee.Image(asset_ws + 'u'), ee.Image(asset_ws + 'u').multiply(0)]) \
            .setMulti({'system:time_start': img_date_start})])
    d_obj._set_weather_vars()
    assert abs(utils.image_value(
        ee.Image(d_obj.windspeed))['windspeed'] - 7.02662301063538) <= tol


@pytest.mark.parametrize(
    'xy,iterations,expected',
    [
        [ne1_xy, 10, 297.00],
        # [ne2_xy, 10, 301.00],
        # [ne3_xy, 10, 302.00],
    ]
)
def test_Image_compute_ta_asset(xy, iterations, expected, tol=0.01):
    """Test fine scale air temperature at a single point using the test assets"""
    d_obj = disalexi.Image(
        test_img,
        elevation=ee.Image.constant(350.0),
        iterations=iterations,
        lc_type='NLCD',
        landcover=ee.Image(asset_ws + 'landcover')
    )

    # Overwrite the default ancillary images with the test assets
    d_obj.windspeed_coll = ee.ImageCollection([
        ee.Image([
                ee.Image(asset_ws + 'u'),
                ee.Image(asset_ws + 'u').multiply(0)]) \
            .setMulti({'system:time_start': img_date_start})])
    d_obj.rs_hourly_coll = ee.ImageCollection([
        ee.Image(asset_ws + 'Insol1')
            .setMulti({'system:time_start': img_hour_start.subtract(3600000)}),
        ee.Image(asset_ws + 'Insol1')
            .setMulti({'system:time_start': img_hour_start}),
        ee.Image(asset_ws + 'Insol1')
            .setMulti({'system:time_start': img_hour_start.add(3600000)})
    ])
    d_obj.rs_daily_coll = ee.ImageCollection([
        ee.Image(asset_ws + 'Insol24')  \
            .setMulti({'system:time_start': img_date_start})])
    d_obj.et_coll = ee.ImageCollection([
        ee.Image(asset_ws + 'alexiET') \
            .setMulti({'system:time_start': img_date_start})])
    d_obj.et_transform = [0.04, 0, -96.442, 0, -0.04, 41.297]

    # Get the spatial reference and geoTransform of the assets
    asset_crs = ee.Image(asset_ws + 'albedo').projection().crs().getInfo()
    asset_transform = ee.Image(asset_ws + 'albedo') \
        .projection().getInfo()['transform']

    # Compute Tair
    ta_img = d_obj.compute_ta() \
        .reproject(crs=asset_crs, crsTransform=asset_transform)

    # Extract image values at a point using reduceRegion (with point geom)
    output = list(utils.image_value(ta_img, xy=xy).values())[0]
    # output = utils.image_value(ta_img)['t_air']

    logging.debug('  Target values: {}'.format(expected))
    logging.debug('  Output values: {}'.format(output))
    assert abs(output - expected) <= tol



@pytest.mark.parametrize(
    'xy,iterations,expected',
    [
        [ne1_xy, 10, {'t_air': 298.0, 'et': 5.995176}],
    ]
)
def test_Image_compute_ta_test_asset(xy, iterations, expected, tol=0.01):
    """Test coarse scale air temperature at a single point using the test assets"""
    d_obj = disalexi.Image(
        test_img,
        elevation=ee.Image.constant(350.0),
        iterations=iterations,
        lc_type='NLCD',
        landcover=ee.Image(asset_ws + 'landcover')
    )

    # Overwrite the default ancillary images with the test assets
    d_obj.windspeed_coll = ee.ImageCollection([
        ee.Image([
            ee.Image(asset_ws + 'u'),
            ee.Image(asset_ws + 'u').multiply(0)]) \
            .setMulti({'system:time_start': img_date_start})])
    d_obj.rs_hourly_coll = ee.ImageCollection([
        ee.Image(asset_ws + 'Insol1')
            .setMulti({'system:time_start': img_hour_start.subtract(3600000)}),
        ee.Image(asset_ws + 'Insol1')
            .setMulti({'system:time_start': img_hour_start}),
        ee.Image(asset_ws + 'Insol1')
            .setMulti({'system:time_start': img_hour_start.add(3600000)})
    ])
    d_obj.rs_daily_coll = ee.ImageCollection([
        ee.Image(asset_ws + 'Insol24')  \
            .setMulti({'system:time_start': img_date_start})])
    d_obj.et_coll = ee.ImageCollection([
        ee.Image(asset_ws + 'alexiET') \
            .setMulti({'system:time_start': img_date_start})])
    d_obj.et_transform = [0.04, 0, -96.442, 0, -0.04, 41.297]

    d_obj._set_solar_vars(interpolate_flag=False)
    d_obj._set_weather_vars()

    # Get the spatial reference and geoTransform of the assets
    # asset_crs = ee.Image(asset_ws + 'albedo').projection().crs().getInfo()
    # asset_transform = ee.Image(asset_ws + 'albedo') \
    #     .projection().getInfo()['transform']

    # Compute ALEXI scale air temperature
    ta_coarse_img = d_obj.compute_ta_test()
    #     .reproject(crs=asset_crs, crsTransform=asset_transform)
    #     .reproject(crs='EPSG:4326', crsTransform=d_obj.et_transform)

    # Extract image values at a point using reduceRegion (with point geom)
    output = utils.image_value(ta_coarse_img.select(['t_air']), tile_scale=4)
    logging.debug('  Target values: {}'.format(expected['t_air']))
    logging.debug('  Output values: {}'.format(output['t_air']))
    assert abs(output['t_air'] - expected['t_air']) <= tol

    # # This test consistently times out, commenting out for now
    # output = utils.image_value(ta_coarse_img.select(['et']))
    # logging.debug('  Target values: {}'.format(expected['et']))
    # logging.debug('  Output values: {}'.format(output['et']))
    # assert abs(output['et'] - expected['et']) <= tol


# @pytest.mark.skip(reason="Skipping until TSEB is working")
# @pytest.mark.parametrize(
#     'xy,iterations,expected',
#     [
#         [ne1_xy, 10, 298.01],
#         [ne1_xy, 20, 298.01],
#         [ne2_xy, 10, 298.01],
#         [ne3_xy, 10, 298.01],
#     ]
# )
# def test_Image_compute_ta_coarse_asset(xy, iterations, expected, tol=0.01):
#     """Test coarse scale air temperature at a single point using the test assets"""
#     d_obj = disalexi.Image(
#         test_img,
#         elevation=ee.Image.constant(350.0),
#         iterations=iterations,
#         lc_type='NLCD',
#         landcover=ee.Image(asset_ws + 'landcover')
#     )
#
#     # Overwrite the default ancillary images with the test assets
#     d_obj.windspeed_coll = ee.ImageCollection([
#         ee.Image([
#             ee.Image(asset_ws + 'u'),
#             ee.Image(asset_ws + 'u').multiply(0)]) \
#             .setMulti({'system:time_start': img_date_start})])
#     d_obj.rs_hourly_coll = ee.ImageCollection([
#         ee.Image(asset_ws + 'Insol1') \
#             .setMulti({'system:time_start': img_hour_start})])
#     d_obj.rs_daily_coll = ee.ImageCollection([
#         ee.Image(asset_ws + 'Insol24')  \
#             .setMulti({'system:time_start': img_date_start})])
#     d_obj.et_coll = ee.ImageCollection([
#         ee.Image(asset_ws + 'alexiET') \
#             .setMulti({'system:time_start': img_date_start})])
#     d_obj.et_transform = [0.04, 0, -96.442, 0, -0.04, 41.297]  # ETd test asset
#
#     # Get the spatial reference and geoTransform of the assets
#     asset_crs = ee.Image(asset_ws + 'albedo').projection().crs().getInfo()
#     asset_transform = ee.Image(asset_ws + 'albedo') \
#         .projection().getInfo()['transform']
#
#     # Compute coarsened Tair
#     # CGM - Would there be some way to chain the methods together?
#     ta_img = d_obj.compute_ta() \
#         .reproject(crs=asset_crs, crsTransform=asset_transform)
#     ta_coarse_img = d_obj.aggregate(ta_img) \
#         .reproject(crs='EPSG:4326', crsTransform=d_obj.et_transform)
#
#     # Extract image values at a point using reduceRegion (with point geom)
#     output = list(utils.image_value(ta_coarse_img).values())[0]
#     # output = utils.image_value(ta_img)['t_air']
#
#     logging.debug('  Target values: {}'.format(expected))
#     logging.debug('  Output values: {}'.format(output))
#     assert abs(output - expected) <= tol




# @pytest.mark.skip(reason="Skipping until TSEB is working")
# def test_disalexi_constant(tol=0.001):
#     """Check output values at a single point for constant image inputs
#
#
#     """
#     logging.debug(
#         '\nGeneric EE DisALEXI testing function using constant image')
#
#     expected = 0.002
#
#     # We could rename the output image to have a common band name so that
#     #   we wouldn't need to know the output band name.
#     # Maybe the band name should be explicitly tested for.
#     output_band = 'et'
#
#     # Generate a totally fake Landsat image using constant images
#     # What is a reasonable BQA value?
#     landsat_img = ee.Image.constant([0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 300, 0]) \
#         .rename(['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B10', 'BQA']) \
#         .setMulti({'system:time_start': ee.Date('2015-08-05', 'GMT').millis()})
#
#     # Generate a totally fake set of inputs using constant images
#     alexi_et_coll = ee.ImageCollection([
#         ee.Image.constant(10).setMulti({
#             'system:time_start': ee.Date('2015-08-05', 'GMT').millis()}),
#         ee.Image.constant(10).setMulti({
#             'system:time_start': ee.Date('2015-08-06', 'GMT').millis()})])
#     landcover_img = ee.Image.constant(11)
#     landcover_type = 'NLCD'
#
#     # Prepare the Landsat 8 image for DisALEXI
#     input_img = ee.Image(landsat.Landsat(landsat_img).prep_landsat8())
#     # logging.info(input_img.getInfo())
#
#     # Initialize DisALEXI object and compute ET
#     def run_model(img):
#         return disalexi.Image(
#                 img, et_coll=alexi_et_coll,
#                 landcover=landcover_img, lc_type=landcover_type) \
#             .run_disalexi()
#     et_img = run_model(input_img).rename(['et'])
#     # logging.info(et_img.getInfo())
#
#     # Extract image values at a point using reduceRegion (with point geom)
#     output = utils.constant_image_value(et_img)['et']
#
#     logging.debug('  Target values: {}'.format(expected))
#     logging.debug('  Output values: {}'.format(output))
#     assert abs(output - expected) <= tol
#
#
# @pytest.mark.skip(reason="Skipping until TSEB is working")
# def test_disalexi_full(tol=0.01):
#     """Check output values at a single point for a real set of inputs"""
#     logging.debug(
#         '\nGeneric EE DisALEXI testing function using constant image')
#
#     expected = 0.089
#
#     # Initialize DisALEXI object
#     # Use fake Alexi ET collection for now
#     alexi_et_coll = ee.ImageCollection([
#         ee.Image.constant(10).setMulti({
#             'system:time_start': ee.Date('2015-08-05', 'GMT').millis()}),
#         ee.Image.constant(10).setMulti({
#             'system:time_start': ee.Date('2015-08-06', 'GMT').millis()})])
#     landcover_img = ee.Image('USGS/NLCD/NLCD2011')
#     landcover_type = 'NLCD'
#
#     # Initialize the Landsat collection (with a single image)
#     landsat_coll = ee.ImageCollection(
#         ee.Image('LANDSAT/LC08/C01/T1_TOA/LC08_043033_20150805'))
#
#     # Prep each image in the Landsat collection
#     def run_prep(img):
#         return ee.Image(landsat.Landsat(img).prep_landsat8())
#     input_coll = ee.ImageCollection(landsat_coll.map(run_prep))
#
#     # CGM - Example of applying prep with lambda call
#     # input_coll = ee.ImageCollection(landsat_coll.map(
#     #     lambda x: ee.Image(landsat.Landsat(ee.Image(x)).prep_landsat8())))
#
#     # Compute ET for each image in the Landsat collection
#     def run_model(img):
#         return disalexi.Image(
#                 img, et_coll=alexi_et_coll,
#                 landcover=landcover_img, lc_type=landcover_type) \
#             .run_disalexi()
#     et_coll = ee.ImageCollection(input_coll.map(run_model))
#     et_img = ee.Image(et_coll.first()).rename(['et'])
#
#     # Extract image values at a point using reduceRegion (with point geom)
#     output = utils.image_value(et_img, xy=(-96.509, 41.2835))['et']
#
#     logging.debug('  Target values: {}'.format(expected))
#     logging.debug('  Output values: {}'.format(output))
#     assert abs(output - expected) <= tol


if __name__ == "__main__":
    pass
    # test_disalexi_asset()
    # test_disalexi_constant()
    # test_disalexi_full()
