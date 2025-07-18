import logging
import pprint

import ee
import pytest

import openet.disalexi as disalexi
import openet.disalexi.utils as utils

COLL_ID = 'LANDSAT/LC08/C02/T1_L2'
SCENE_ID = 'LC08_044033_20170716'
SCENE_TIME = 1500230731090
SCENE_POINT = (-119.5, 36.0)
TEST_POINT = (-121.5265, 38.7399)

# AmeriFlux sites adjusted to nearest Landsat cell centroid
ne1_xy = [-96.47672812080845, 41.16506126041818]
ne2_xy = [-96.46994024736414, 41.16491226772292]
ne3_xy = [-96.43968912903934, 41.17964494123755]


# DisALEXI needs a non-constant image since it computes lat/lon
# LAI and LST are read from sources and net set in the "prep" functions anymore
def default_image(
        albedo=0.125,
        cfmask=0,
        ndvi=0.875,
        # lai=4.2,
        # lst=306.5,
):
    return (
        ee.Image(f'{COLL_ID}/{SCENE_ID}')
        .select(['SR_B2', 'SR_B3', 'SR_B4'])
        .multiply([0, 0, 0]).add([albedo, cfmask, ndvi])
        .rename(['albedo', 'cfmask', 'ndvi'])
        .set({
            'system:index': SCENE_ID,
            'system:time_start': SCENE_TIME,
            'system:id': f'{COLL_ID}/{SCENE_ID}',
        })
)

def default_image_args(
        albedo=0.125,
        cfmask=0,
        ndvi=0.875,
        ta_source='projects/openet/assets/disalexi/tair/conus_v006_1k',
        alexi_source='CONUS_V006',
        # lai_source='openet-landsat-lai',
        lai_source='projects/openet/assets/lai/landsat/c02',
        lst_source='projects/openet/assets/lst/landsat/c02',
        # lai_source=4.2,
        # lst_source=306.5,
        et_reference_source=10,
        et_reference_band=None,
        et_reference_factor=1.0,
        et_reference_resample='nearest',
        stability_iterations=10,
        albedo_iterations=10,
):
    return {
        'image': default_image(albedo=albedo, cfmask=cfmask, ndvi=ndvi),
        'ta_source': ta_source,
        'alexi_source': alexi_source,
        'lai_source': lai_source,
        'lst_source': lst_source,
        'et_reference_source': et_reference_source,
        'et_reference_band': et_reference_band,
        'et_reference_factor': et_reference_factor,
        'et_reference_resample': et_reference_resample,
        'stability_iterations': stability_iterations,
        'albedo_iterations': albedo_iterations,
    }

def default_image_obj(
        albedo=0.125,
        cfmask=0,
        ndvi=0.875,
        ta_source='projects/openet/assets/disalexi/tair/conus_v006_1k',
        alexi_source='CONUS_V006',
        #lai_source='openet-landsat-lai',
        lai_source='projects/openet/assets/lai/landsat/c02',
        lst_source='projects/openet/assets/lst/landsat/c02',
        # lai_source=4.2,
        # lst_source=306.5,
        et_reference_source=10,
        et_reference_band=None,
        et_reference_factor=1.0,
        et_reference_resample='nearest',
        stability_iterations=10,
        albedo_iterations=10,
):
    return disalexi.Image(**default_image_args(
        albedo=albedo,
        cfmask=cfmask,
        ndvi=ndvi,
        ta_source=ta_source,
        alexi_source=alexi_source,
        lai_source=lai_source,
        lst_source=lst_source,
        et_reference_source=et_reference_source,
        et_reference_band=et_reference_band,
        et_reference_factor=et_reference_factor,
        et_reference_resample=et_reference_resample,
        stability_iterations=stability_iterations,
        albedo_iterations=albedo_iterations,
    ))


def test_Image_ta_coarse_initial_properties():
    d_obj = disalexi.Image(**default_image_args())
    output = utils.getinfo(d_obj.ta_coarse_initial())
    assert output['bands'][0]['id'] == 'ta_initial'
    assert output['bands'][0]['crs'] == 'EPSG:4326'
    assert output['bands'][0]['crs_transform'][0] == 0.04
    assert output['bands'][0]['crs_transform'][4] == -0.04


def test_Image_ta_coarse_initial_values():
    d_obj = disalexi.Image(**default_image_args())
    ta_img = d_obj.ta_coarse_initial()
    output = utils.point_image_value(ta_img, TEST_POINT)
    assert output['ta_initial'] > 0


def test_Image_ta_coarse_mosaic_step_count(step_size=5, step_count=4):
    # Check that the output step counts match the input (+1)
    # Use the ALEXI image as the air temperature mask
    d_obj = disalexi.Image(**default_image_args())
    ta_img = d_obj.et_alexi.multiply(0).add(290)
    ta_mosaic_img = d_obj.ta_coarse_mosaic(ta_img, offsets=[-10, -5, 0, 5, 10])
    output = utils.getinfo(ta_mosaic_img)
    ta_bands = [b for b in output['bands'] if '_ta' in b['id']]
    bias_bands = [b for b in output['bands'] if '_bias' in b['id']]
    assert len(ta_bands) == (step_count + 1)
    assert len(bias_bands) == (step_count + 1)


def test_Image_ta_coarse_mosaic_step_size(ta_init=290, step_size=5, step_count=4):
    # Check that the output step sizes match the input
    # Use the ALEXI image as the air temperature mask
    d_obj = disalexi.Image(**default_image_args())
    ta_img = d_obj.et_alexi.multiply(0).add(ta_init)
    ta_mosaic_img = d_obj.ta_coarse_mosaic(ta_img, offsets=[-10, -5, 0, 5, 10])
    output = utils.point_image_value(ta_mosaic_img, TEST_POINT)
    assert (abs(output['step_01_ta'] - output['step_00_ta']) - 5) < 0.000001


def test_Image_ta_coarse_mosaic_min_bias(ta_init=290):
    # Check that a reasonable air temperature value is returned
    d_obj = disalexi.Image(**default_image_args())
    ta_mosaic_img = d_obj.ta_coarse_mosaic(
        d_obj.et_alexi.multiply(0).add(ta_init), offsets=[-10, -5, 0, 5, 10]
    )
    ta_min_bias = disalexi.ta_mosaic_min_bias(ta_mosaic_img)
    output = utils.point_image_value(ta_min_bias, TEST_POINT)
    assert abs(output['ta'] - 290) < 10


def test_Image_ta_mosaic_interpolate(ta_init=290):
    # Check that a reasonable air temperature value is returned
    d_obj = disalexi.Image(**default_image_args())
    ta_mosaic_img = d_obj.ta_coarse_mosaic(
        d_obj.et_alexi.multiply(0).add(ta_init), offsets=[-10, -5, 0, 5, 10]
    )
    ta_interp = disalexi.ta_mosaic_interpolate(ta_mosaic_img)
    output = utils.point_image_value(ta_interp, TEST_POINT)
    assert abs(output['ta_interp'] - 290) < 10


@pytest.mark.parametrize(
    'ta_values, bias_values, expected',
    [
        [[280, 285, 290, 295, 300], [-4, -3, -2, 1, 3], 293.33],
        [[280, 285, 290, 295, 300], [-2, -2, -1, 1, 3], 292.5],
        # Constant negative bias values should pick the last one before going positive
        [[280, 285, 290, 295, 300, 305, 310], [-2, -2, -2, 1, 3, 4, 5], 293.33],
        [[270, 275, 280, 285, 290, 295, 300, 305], [-2, -2, -2, -2, -2, 1, 3, 4], 293.33],
        # Test multiple negative to positive bias transitions
        [[270, 275, 280, 285, 290, 295, 300, 305], [-2, -0.1, 0.1, -2, -1, 1, 3, 4], 293.33],
        # Test that pixels with all positive or negative biases are masked out
        # This is different than the original implementation
        # The original unmasked test values are commented out below
        [[280, 285, 290, 295, 300], [1, 2, 3, 4, 5], None],
        [[270, 275, 280, 285, 290, 295, 300], [2, 2, 2, 2, 3, 4, 5], None],
        [[280, 285, 290, 295, 300], [-5, -4, -3, -2, -1], None],
        # # DEADBEEF
        # #   The following tests can eventually be removed
        # #   if pixels that have all positive or negative bias values are masked
        # # All positive biases
        # [[280, 285, 290, 295, 300], [1, 2, 3, 4, 5], 280],
        # # Test that the Ta for the first "increasing" bias is used
        # [[270, 275, 280, 285, 290, 295, 300], [2, 2, 2, 2, 3, 4, 5], 290],
        # # In the original 10k processing, pixels with all negative biases were masked
        # # This is commented out for now but may be added back in and should be tested for
        # [[280, 285, 290, 295, 300], [-5, -4, -3, -2, -1], 300],
    ]
)
def test_Image_ta_mosaic_interpolate_values(ta_values, bias_values, expected):
    # Check that a reasonable air temperature value is returned
    d_obj = disalexi.Image(**default_image_args())
    mask_img = d_obj.et_alexi.multiply(0)
    ta_images = [mask_img.add(x).rename(f'step_{i:02d}_ta') for i, x in enumerate(ta_values)]
    bias_images = [mask_img.add(x).rename(f'step_{i:02d}_bias') for i, x in enumerate(bias_values)]
    ta_mosaic_img = ee.Image(ta_images + bias_images)
    ta_interp = disalexi.ta_mosaic_interpolate(ta_mosaic_img)
    output = utils.point_image_value(ta_interp, TEST_POINT)
    if expected is None:
        assert output['ta_interp'] is None
    else:
        assert abs(output['ta_interp'] - expected) < 1


@pytest.mark.parametrize(
    'ta_values, bias_values, expected',
    [
        [[280, 285, 290, 295, 300], [-4, -3, -2, 1, 3], 295],    # Normal bias profile
        [[280, 285, 290, 295, 300], [-2, -2, -2, 1, 3], 295],    # Constant negative biases
        [[280, 285, 290, 295, 300], [1, 2, 3, 4, 5], 280],       # All positive biases
        [[280, 285, 290, 295, 300], [-5, -4, -3, -2, -1], 300],  # All negative biases
    ]
)
def test_Image_ta_mosaic_min_bias_values(ta_values, bias_values, expected):
    d_obj = disalexi.Image(**default_image_args())
    mask_img = d_obj.et_alexi.multiply(0)
    ta_images = [mask_img.add(x).rename(f'step_{i:02d}_ta') for i, x in enumerate(ta_values)]
    bias_images = [mask_img.add(x).rename(f'step_{i:02d}_bias') for i, x in enumerate(bias_values)]
    ta_mosaic_img = ee.Image(ta_images + bias_images)
    ta_interp = disalexi.ta_mosaic_min_bias(ta_mosaic_img)
    output = utils.point_image_value(ta_interp, TEST_POINT)
    assert abs(output['ta'] - expected) < 1


def test_Image_set_landcover_vars_default(tol=1E-6):
    """Test default land cover image and type

    It might make more sense to just test that the value at the test pixel
    is 82 (for NLCD) for 10 (for GLC30)"""
    d_obj = disalexi.Image(**default_image_args())
    d_obj.set_landcover_vars()
    assert utils.point_image_value(d_obj.aleafv, ne1_xy, scale=30)['aleafv'] == 0.83
    assert utils.point_image_value(d_obj.aleafn, ne1_xy, scale=30)['aleafn'] == 0.35
    assert utils.point_image_value(d_obj.aleafl, ne1_xy, scale=30)['aleafl'] == 0.95
    assert utils.point_image_value(d_obj.adeadv, ne1_xy, scale=30)['adeadv'] == 0.49
    assert utils.point_image_value(d_obj.adeadn, ne1_xy, scale=30)['adeadn'] == 0.13
    assert utils.point_image_value(d_obj.adeadl, ne1_xy, scale=30)['adeadl'] == 0.95
    assert utils.point_image_value(d_obj.leaf_width, ne1_xy, scale=30)['xl'] == 0.05
    assert utils.point_image_value(d_obj.clump, ne1_xy, scale=30)['omega'] == 0.83


def test_Image_set_landcover_vars_init_asset(tol=1E-6):
    """Test setting the land cover image and type as the object is initialized"""
    d_obj = disalexi.Image(landcover_source='USGS/NLCD_RELEASES/2019_REL/NLCD/2011', **default_image_args())
    d_obj.set_landcover_vars()
    # Only need to check the first test value
    assert utils.point_image_value(d_obj.aleafv, ne1_xy, scale=30)['aleafv'] == 0.83


def test_Image_set_landcover_vars_set_asset(tol=1E-6):
    """Test setting the land cover image and type directly on the object"""
    d_obj = disalexi.Image(**default_image_args())
    d_obj.landcover_source = 'USGS/NLCD_RELEASES/2019_REL/NLCD/2011'
    d_obj.set_landcover_vars()
    # Only need to check the first test value
    assert utils.point_image_value(d_obj.aleafv, ne1_xy, scale=30)['aleafv'] == 0.83
