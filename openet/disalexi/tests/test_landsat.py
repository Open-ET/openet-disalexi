import pprint

import ee
import pytest

import openet.disalexi.landsat as landsat
import openet.disalexi.utils as utils


# Default property values
l8_properties = {
    'system:time_start': 1404839150550, 'system:index': 'LC08_028031_20140708',
    'SPACECRAFT_ID': 'LANDSAT_8', 'SATELLITE': 'LANDSAT_8',
    'K1_CONSTANT_BAND_10': 774.8853, 'K2_CONSTANT_BAND_10': 1321.0789}
l7_properties = {
    'system:time_start': 992192137355, 'system:index': 'LE07_028031_20010610',
    'SPACECRAFT_ID': 'LANDSAT_7', 'SATELLITE': 'LANDSAT_7',
    'K1_CONSTANT_BAND_6_VCID_1': 666.09, 'K2_CONSTANT_BAND_6_VCID_1': 1282.71}
l5_properties = {
    'system:time_start': 988735571649, 'system:index': 'LT05_028031_20010501',
    'SPACECRAFT_ID': 'LANDSAT_5', 'SATELLITE': 'LANDSAT_5',
    'K1_CONSTANT_BAND_6': 607.76, 'K2_CONSTANT_BAND_6': 1260.56}
sr_scalars = [0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.1, 1]


# LandsatTOA tests
@pytest.mark.parametrize(
    'img_id',
    [
        'LANDSAT/LC08/C01/T1_RT_TOA/LC08_028031_20140708',
        'LANDSAT/LC08/C01/T1_TOA/LC08_028031_20140708',
        'LANDSAT/LE07/C01/T1_RT_TOA/LE07_028031_20010610',
        'LANDSAT/LE07/C01/T1_TOA/LE07_028031_20010610',
        'LANDSAT/LT05/C01/T1_TOA/LT05_028031_20010501',
    ]
)
def test_LandsatTOA_init(img_id):
    """Test that the Landsat bands are renamed and properties are copied"""
    l_info = ee.Image(img_id).getInfo()['properties']
    input_img = landsat.LandsatTOA(ee.Image(img_id)).input_image
    input_info = ee.Image(input_img).getInfo()
    assert [b['id'] for b in input_info['bands']] == [
        'blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'lst', 'bqa']
    assert input_info['properties']['system:time_start'] == l_info['system:time_start']
    assert input_info['properties']['system:index'] == img_id.split('/')[-1]
    assert input_info['properties']['k1_constant']
    assert input_info['properties']['k2_constant']


def test_LandsatTOA_prep(img_id='LANDSAT/LC08/C01/T1_RT_TOA/LC08_028031_20140708'):
    """Test that the prepped image has the target bands and properties"""
    l_info = ee.Image(img_id).getInfo()['properties']
    prepped_img = landsat.LandsatTOA(ee.Image(img_id)).prep()
    prepped_info = ee.Image(prepped_img).getInfo()
    assert sorted([b['id'] for b in prepped_info['bands']]) == [
        'albedo', 'cfmask', 'lai', 'lst', 'ndvi']
    assert prepped_info['properties']['system:time_start'] == l_info['system:time_start']
    assert prepped_info['properties']['system:index'] == img_id.split('/')[-1]


@pytest.mark.parametrize(
    'blue, green, red, nir, swir1, swir2',
    [
        [0.2, 0.1, 0.2, 0.2, 0.2, 0.2],
        [0.2, 0.9, 0.2, 0.2, 0.2, 0.2]
    ]
)
def test_LandsatTOA_albedo(blue, green, red, nir, swir1, swir2, tol=0.000001):
    """Test the albedo calculation

    Ensure that the Green band is not being used to compute albedo
    """
    expected = sum([a * b for a, b in zip(
        [blue, red, nir, swir1, swir2, 1],
        [0.356, 0.130, 0.373, 0.085, 0.072, -0.0018])])

    raw_img = ee.Image.constant([blue, green, red, nir, swir1, swir2, 300, 0]) \
        .rename(['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B10', 'BQA']) \
        .setMulti(l8_properties)
    albedo = ee.Image(landsat.LandsatTOA(ee.Image(raw_img))._albedo)
    assert abs(utils.image_value(albedo)['albedo'] - expected) <= tol


@pytest.mark.parametrize(
    'bqa, expected',
    [
        ['0000000000000000', 0],
        ['0000000000000001', 1],
        ['0000000110000000', 2],
        ['0000011000000000', 3],
        ['0000000000000000', 0],
        ['0000000000100000', 0],
        ['0000000001000000', 4],
        ['0000000001100000', 4]
    ]
)
def test_LandsatTOA_cfmask(bqa, expected):
    input_img = ee.Image.constant([0.2, 0, 0, 0, 0, 0, 300, int(bqa, 2)]) \
        .rename(['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B10', 'BQA']) \
        .setMulti(l8_properties)
    cfmask = ee.Image(landsat.LandsatTOA(input_img)._cfmask)
    assert utils.image_value(cfmask)['cfmask'] == expected


def test_LandsatTOA_lai(red=0.2, nir=0.7, expected=1.200, tol=0.001):
    input_img = ee.Image.constant([0.2, 0.2, red, nir, 0.2, 0.2, 300, 0]) \
        .rename(['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B10', 'BQA']) \
        .setMulti(l8_properties)
    lai = ee.Image(landsat.LandsatTOA(input_img)._lai)
    assert abs(utils.image_value(lai)['lai'] - expected) <= tol


@pytest.mark.parametrize(
    'red, nir, bt, expected',
    [
        [0.2, 0.7, 300, 304.5866],
        [0.2, 0.3, 300, 304.8238],
        [0.2, 0.1, 300, 303.6067],
    ]
)
def test_LandsatTOA_lst(red, nir, bt, expected, tol=0.001):
    """Test that different emissivity values (from NDVI & LAI) change LST"""
    input_img = ee.Image.constant([0.2, 0.2, red, nir, 0.2, 0.2, bt, 0]) \
        .rename(['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B10', 'BQA']) \
        .setMulti(l8_properties)
    lst = ee.Image(landsat.LandsatTOA(input_img)._lst)
    assert abs(utils.image_value(lst)['lst'] - expected) <= tol


def test_LandsatTOA_ndvi(red=0.2, nir=0.7, expected=0.5556, tol=0.001):
    input_img = ee.Image.constant([0.2, 0.2, red, nir, 0.2, 0.2, 300, 0]) \
        .rename(['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B10', 'BQA']) \
        .setMulti(l8_properties)
    ndvi = ee.Image(landsat.LandsatTOA(input_img)._ndvi)
    assert abs(utils.constant_image_value(ndvi)['ndvi'] - expected) <= tol


# LandsatSR tests
@pytest.mark.parametrize(
    'img_id',
    [
        'LANDSAT/LC08/C01/T1_SR/LC08_028031_20140708',
        'LANDSAT/LE07/C01/T1_SR/LE07_028031_20010610',
        'LANDSAT/LT05/C01/T1_SR/LT05_028031_20010501',
    ]
)
def test_LandsatSR_init(img_id):
    """Test that the Landsat bands are renamed and properties are copied"""
    l_info = ee.Image(img_id).getInfo()['properties']
    input_img = landsat.LandsatSR(ee.Image(img_id)).input_image
    input_info = ee.Image(input_img).getInfo()
    assert [b['id'] for b in input_info['bands']] == [
        'blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'lst', 'pixel_qa']
    assert input_info['properties']['system:time_start'] == l_info['system:time_start']
    assert input_info['properties']['system:index'] == img_id.split('/')[-1]
    assert input_info['properties']['k1_constant']
    assert input_info['properties']['k2_constant']


def test_LandsatSR_prep(img_id='LANDSAT/LC08/C01/T1_SR/LC08_028031_20140708'):
    """Test that the prepped image has the target bands and properties"""
    l_info = ee.Image(img_id).getInfo()['properties']

    prepped_img = landsat.LandsatSR(ee.Image(img_id)).prep()
    prepped_info = ee.Image(prepped_img).getInfo()
    assert sorted([b['id'] for b in prepped_info['bands']]) == [
        'albedo', 'cfmask', 'lai', 'lst', 'ndvi']
    assert prepped_info['properties']['system:time_start'] == l_info['system:time_start']
    assert prepped_info['properties']['system:index'] == img_id.split('/')[-1]


@pytest.mark.parametrize(
    'blue, green, red, nir, swir1, swir2',
    [
        [0.2, 0.1, 0.2, 0.2, 0.2, 0.2],
        [0.2, 0.9, 0.2, 0.2, 0.2, 0.2]
    ]
)
def test_LandsatSR_albedo(blue, green, red, nir, swir1, swir2, tol=0.000001):
    """Test the albedo calculation

    Ensure that the Green band is not being used to compute albedo
    """
    expected = sum([a * b for a, b in zip(
        [blue, red, nir, swir1, swir2, 1],
        [0.356, 0.130, 0.373, 0.085, 0.072, -0.0018])])

    # Undo the scalars that are applied in the init
    raw_img = ee.Image.constant([blue, green, red, nir, swir1, swir2, 300, 0]) \
        .rename(['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B10', 'pixel_qa']) \
        .divide(sr_scalars) \
        .setMulti(l8_properties)
    albedo = ee.Image(landsat.LandsatSR(ee.Image(raw_img))._albedo)
    assert abs(utils.image_value(albedo)['albedo'] - expected) <= tol


@pytest.mark.parametrize(
    'pixel_qa, expected',
    [
        # Note - The cfmask function is only check cloud confidence bits (6-7),
        #   not the cloud bit (5) itself
        ['0000000000000000', 0],
        ['0000000000000001', 1],
        ['0000000000001000', 2],
        ['0000000000010000', 3],
        ['0000000001000000', 0],
        ['0000000010000000', 4],
        ['0000000011000000', 4],
        ['0000000000100000', 0],
        ['0000000001100000', 0],
        ['0000000010100000', 4],
        ['0000000011100000', 4],
    ]
)
def test_LandsatSR_cfmask(pixel_qa, expected):
    input_img = ee.Image.constant([0.2, 0, 0, 0, 0, 0, 300, int(pixel_qa, 2)]) \
        .rename(['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B10', 'pixel_qa']) \
        .setMulti(l8_properties)
    cfmask = ee.Image(landsat.LandsatSR(input_img)._cfmask)
    assert utils.image_value(cfmask)['cfmask'] == expected


def test_LandsatSR_lai(red=0.2, nir=0.7, expected=1.200, tol=0.001):
    # Undo the scalars that are applied in the init
    input_img = ee.Image.constant([0.2, 0.2, red, nir, 0.2, 0.2, 300, 0]) \
        .rename(['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B10', 'pixel_qa']) \
        .divide(sr_scalars) \
        .setMulti(l8_properties)
    lai = ee.Image(landsat.LandsatSR(input_img)._lai)
    assert abs(utils.image_value(lai)['lai'] - expected) <= tol


@pytest.mark.parametrize(
    'red, nir, bt, expected',
    [
        [0.2, 0.7, 300, 304.5866],
        [0.2, 0.3, 300, 304.8238],
        [0.2, 0.1, 300, 303.6067],
        # # Current LandsatSR LST function just returns brightness Ts
        # [0.2, 0.7, 300, 300],
    ]
)
def test_LandsatSR_lst(red, nir, bt, expected, tol=0.001):
    """Test that different emissivity values (from NDVI & LAI) change LST"""
    input_img = ee.Image.constant([0.2, 0.2, red, nir, 0.2, 0.2, bt, 0]) \
        .rename(['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B10', 'pixel_qa']) \
        .divide(sr_scalars) \
        .setMulti(l8_properties)
    lst = ee.Image(landsat.LandsatSR(input_img)._lst)
    assert abs(utils.image_value(lst)['lst'] - expected) <= tol


def test_LandsatSR_ndvi(red=0.2, nir=0.7, expected=0.5556, tol=0.001):
    input_img = ee.Image.constant([0.2, 0.2, red, nir, 0.2, 0.2, 300, 0]) \
        .rename(['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B10', 'pixel_qa']) \
        .divide(sr_scalars) \
        .setMulti(l8_properties)
    ndvi = ee.Image(landsat.LandsatSR(input_img)._ndvi)
    assert abs(utils.constant_image_value(ndvi)['ndvi'] - expected) <= tol
