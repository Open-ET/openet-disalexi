import pprint

import ee
import pytest

import openet.disalexi.landsat as landsat
import openet.disalexi.utils as utils


# Default property values
l8_properties = {
    'system:time_start': 1404839150550, 'system:index': 'LC08_028031_20140708',
    'SPACECRAFT_ID': 'LANDSAT_8', 'SATELLITE': 'LANDSAT_8',
}
l7_properties = {
    'system:time_start': 992192137355, 'system:index': 'LE07_028031_20010610',
    'SPACECRAFT_ID': 'LANDSAT_7', 'SATELLITE': 'LANDSAT_7',
}
l5_properties = {
    'system:time_start': 988735571649, 'system:index': 'LT05_028031_20010501',
    'SPACECRAFT_ID': 'LANDSAT_5', 'SATELLITE': 'LANDSAT_5',
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
    assert set([b['id'] for b in input_info['bands']]) == set([
        'blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'lst', 'QA_PIXEL'])
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
        [0.2, 0.9, 0.2, 0.2, 0.2, 0.2]
    ]
)
def test_Landsat_C02_L2_albedo(blue, green, red, nir, swir1, swir2, tol=0.000001):
    """Test the albedo calculation

    Ensure that the Green band is not being used to compute albedo
    """
    expected = sum([a * b for a, b in zip(
        [blue, red, nir, swir1, swir2, 1],
        [0.356, 0.130, 0.373, 0.085, 0.072, -0.0018]
    )])

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


# DEADBEEF - LAI is being read from a source image collection
#   Leaving this test since existing lai method is used in emissivity calculation
def test_Landsat_C02_L2_lai(red=0.2, nir=0.7, expected=1.200, tol=0.001):
    # Undo the scalars that are applied in the init
    input_img = (
        ee.Image.constant([0.2, 0.2, red, nir, 0.2, 0.2, 300, 0])
        .rename(['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7', 'ST_B10', 'QA_PIXEL'])
        .subtract(Landsat_C02_L2_scalars_add)
        .divide(Landsat_C02_L2_scalars_multi)
        .set(l8_properties)
    )
    lai = ee.Image(landsat.Landsat_C02_L2(input_img)._lai)
    assert abs(utils.constant_image_value(lai)['lai'] - expected) <= tol


# # DEADBEEF - LST is being read from a source image collection
# @pytest.mark.parametrize(
#     'red, nir, bt, expected',
#     [
#         [0.2, 0.7, 300, 304.5866],
#         [0.2, 0.3, 300, 304.8238],
#         [0.2, 0.1, 300, 303.6067],
#     ]
# )
# def test_Landsat_C02_L2_lst(red, nir, bt, expected, tol=0.001):
#     """Test that different emissivity values (from NDVI & LAI) change LST"""
#     input_img = (
#         ee.Image.constant([0.2, 0.2, red, nir, 0.2, 0.2, bt, 0])
#         .rename(['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7', 'ST_B10', 'QA_PIXEL'])
#         .subtract(Landsat_C02_L2_scalars_add)
#         .divide(Landsat_C02_L2_scalars_multi)
#         .set(l8_properties)
#     )
#     lst = ee.Image(landsat.Landsat_C02_L2(input_img)._lst)
#     assert abs(utils.constant_image_value(lst)['lst'] - expected) <= tol


def test_Landsat_C02_L2_ndvi(red=0.2, nir=0.7, expected=0.5556, tol=0.001):
    input_img = (
        ee.Image.constant([0.2, 0.2, red, nir, 0.2, 0.2, 300, 0])
        .rename(['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7', 'ST_B10', 'QA_PIXEL'])
        .subtract(Landsat_C02_L2_scalars_add)
        .divide(Landsat_C02_L2_scalars_multi)
        .set(l8_properties)
    )
    ndvi = ee.Image(landsat.Landsat_C02_L2(input_img)._ndvi)
    assert abs(utils.constant_image_value(ndvi)['ndvi'] - expected) <= tol