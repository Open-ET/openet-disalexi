import pprint

import ee
import pytest

import openet.disalexi.landsat as landsat
import openet.disalexi.utils as utils


# Default property values
l8_properties = {
    'system:time_start': 1404839150550,
    'system:index': 'LC08_028031_20140708',
    'SPACECRAFT_ID': 'LANDSAT_8',
    'SATELLITE': 'LANDSAT_8',
}
l7_properties = {
    'system:time_start': 992192137355,
    'system:index': 'LE07_028031_20010610',
    'SPACECRAFT_ID': 'LANDSAT_7',
    'SATELLITE': 'LANDSAT_7',
}
l5_properties = {
    'system:time_start': 988735571649,
    'system:index': 'LT05_028031_20010501',
    'SPACECRAFT_ID': 'LANDSAT_5',
    'SATELLITE': 'LANDSAT_5',
}

Landsat_C02_L2_scalars_multi = [
    0.0000275, 0.0000275, 0.0000275, 0.0000275, 0.0000275, 0.0000275, 0.00341802, 1
]
Landsat_C02_L2_scalars_add = [-0.2, -0.2, -0.2, -0.2, -0.2, -0.2, 149.0, 1]


# Landsat C02 SR tests
@pytest.mark.parametrize(
    'img_id',
    [
        'LANDSAT/LC09/C02/T1_L2/LC09_044033_20220127',
        'LANDSAT/LC08/C02/T1_L2/LC08_028031_20140708',
        'LANDSAT/LE07/C02/T1_L2/LE07_028031_20010610',
        'LANDSAT/LT05/C02/T1_L2/LT05_028031_20010501',
    ]
)
def test_Landsat_C02_L2_init(img_id):
    """Test that the Landsat bands are renamed and properties are copied"""
    l_info = ee.Image(img_id).getInfo()['properties']
    input_img = landsat.Landsat_C02_L2(ee.Image(img_id)).input_image
    input_info = ee.Image(input_img).getInfo()
    bands = {'blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'lst', 'QA_PIXEL'}
    assert {b['id'] for b in input_info['bands']} == bands
    assert input_info['properties']['system:time_start'] == l_info['system:time_start']
    assert input_info['properties']['system:index'] == img_id.split('/')[-1]


def test_Landsat_C02_L2_prep(img_id='LANDSAT/LC08/C02/T1_L2/LC08_028031_20140708'):
    """Test that the prepped image has the target bands and properties"""
    l_info = ee.Image(img_id).getInfo()['properties']
    prepped_img = landsat.Landsat_C02_L2(ee.Image(img_id)).prep()
    prepped_info = ee.Image(prepped_img).getInfo()
    assert {b['id'] for b in prepped_info['bands']} == {'albedo', 'cfmask', 'ndvi'}
    assert prepped_info['properties']['system:time_start'] == l_info['system:time_start']
    assert prepped_info['properties']['system:index'] == img_id.split('/')[-1]


@pytest.mark.parametrize(
    'blue, green, red, nir, swir1, swir2',
    [
        [0.2, 0.1, 0.2, 0.2, 0.2, 0.2],
        [0.2, 0.9, 0.2, 0.2, 0.2, 0.2],
        # Check albedo calculation for saturated pixels
        [1.6, 1.6, 1.6, 1.0, 1.0, 1.0],
        # Check albedo calculation for negative pixels
        [-0.1, -0.1, -0.1, 0.0, 0.0, 0.0],
    ]
)
def test_Landsat_C02_L2_albedo(blue, green, red, nir, swir1, swir2, tol=0.000001):
    """Test the albedo calculation

    Note that the Green band is not being used to compute albedo
    """
    expected = sum([max(min(a, 1), 0) * b for a, b in zip(
        [blue, red, nir, swir1, swir2, 1],
        [0.356, 0.130, 0.373, 0.085, 0.072, -0.0018]
    )])
    expected = max(0, expected)

    # Undo the scalars that are applied in the init
    raw_img = (
        ee.Image.constant([blue, green, red, nir, swir1, swir2, 300, 0])
        .rename(['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7', 'ST_B10', 'QA_PIXEL'])
        .subtract(Landsat_C02_L2_scalars_add)
        .divide(Landsat_C02_L2_scalars_multi)
        .set(l8_properties)
    )
    albedo = ee.Image(landsat.Landsat_C02_L2(ee.Image(raw_img))._albedo)
    assert abs(utils.constant_image_value(albedo)['albedo'] - expected) <= tol


@pytest.mark.parametrize(
    'pixel_qa, expected',
    [
        ['0000000000000000', 0],  # Designated Fill
        ['0000000000000001', 0],
        ['0000000000000010', 0],  # Dilated Cloud
        ['0000000000000100', 0],  # Cirrus
        ['0000000000001000', 4],  # Cloud
        ['0000000000010000', 2],  # Cloud Shadow
        ['0000000000100000', 3],  # Snow
        ['0000000001000000', 0],  # Clear
        ['0000000010000000', 0],  # Water
    ]
)
def test_Landsat_C02_L2_cfmask(pixel_qa, expected):
    input_img = (
        ee.Image.constant([0.2, 0, 0, 0, 0, 0, 300, int(pixel_qa, 2)])
        .rename(['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7', 'ST_B10', 'QA_PIXEL'])
        .set(l8_properties)
    )
    cfmask = ee.Image(landsat.Landsat_C02_L2(input_img)._cfmask)
    assert utils.constant_image_value(cfmask)['cfmask'] == expected


def test_Landsat_C02_L2_lai(red=0.2, nir=0.7, expected=1.200, tol=0.001):
    # Undo the scalars that are applied in the init
    input_img = (
        ee.Image.constant([0.2, 0.2, red, nir, 0.2, 0.2, 300, 0])
        .rename(['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7', 'ST_B10', 'QA_PIXEL'])
        .subtract(Landsat_C02_L2_scalars_add)
        .divide(Landsat_C02_L2_scalars_multi)
        .set(l8_properties)
    )
    lai = ee.Image(landsat.Landsat_C02_L2(input_img, gsw_extent_flag=False)._lai)
    assert abs(utils.constant_image_value(lai)['lai'] - expected) <= tol


def test_Landsat_C02_L2_ndvi(red=0.2, nir=0.7, expected=0.5556, tol=0.001):
    input_img = (
        ee.Image.constant([0.2, 0.2, red, nir, 0.2, 0.2, 300, 0])
        .rename(['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7', 'ST_B10', 'QA_PIXEL'])
        .subtract(Landsat_C02_L2_scalars_add)
        .divide(Landsat_C02_L2_scalars_multi)
        .set(l8_properties)
    )
    ndvi = ee.Image(landsat.Landsat_C02_L2(input_img, gsw_extent_flag=False)._ndvi)
    assert abs(utils.constant_image_value(ndvi)['ndvi'] - expected) <= tol


@pytest.mark.parametrize(
    'red, nir, expected',
    [
        [0.2, 9.0 / 55, -0.1],
        [0.2, 0.2, 0.0],
        [0.1, 11.0 / 90,  0.1],
        [0.2, 0.3, 0.2],
        [0.1, 13.0 / 70, 0.3],
        [0.3, 0.7, 0.4],
        [0.2, 0.6, 0.5],
        [0.2, 0.8, 0.6],
        [0.1, 17.0 / 30, 0.7],
        [0.2, 0.7, 0.55555555],
        # Check that negative values are not masked
        [-0.01, 0.1, 1.0],
        [0.1, -0.01, -1.0],
        # Check that low values are set to 0
        [-0.1, -0.1, 0.0],
        [0.0, 0.0, 0.0],
        [0.009, 0.009, 0.0],
        [0.009, -0.01, 0.0],
        [-0.01, 0.009, 0.0],
        # Don't adjust NDVI if only one reflectance value is low
        [0.005, 0.1, 0.9047619104385376],
    ]
)
def test_Landsat_C02_L2_ndvi_calculation(red, nir, expected, tol=0.000001):
    input_img = (
        ee.Image.constant([0.2, 0.2, red, nir, 0.2, 0.2, 300, 0])
        .rename(['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7', 'ST_B10', 'QA_PIXEL'])
        .subtract(Landsat_C02_L2_scalars_add)
        .divide(Landsat_C02_L2_scalars_multi)
        .set(l8_properties)
    )
    ndvi = ee.Image(landsat.Landsat_C02_L2(input_img, gsw_extent_flag=False)._ndvi)
    assert abs(utils.constant_image_value(ndvi)['ndvi'] - expected) <= tol
