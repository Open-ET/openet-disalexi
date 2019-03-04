import ee
import pytest

import openet.disalexi.utils as utils


@pytest.fixture
def constant_image_value(image):
    return image.reduceRegion(
        reducer=ee.Reducer.mean(),
        geometry=ee.Geometry.Rectangle([0, 0, 10, 10], 'EPSG:32610', False),
        scale=1).getInfo()


def test_constant_image_value(tol=0.000001):
    expected = 10.123456789
    input = ee.Image.constant(expected)
    output = utils.constant_image_value(input)['constant']
    assert abs(output - expected) <= tol


@pytest.mark.parametrize(
    'lon, lat, expected',
    [
        # US-NE1
        [-96.47672812080845, 41.16506126041818, 356.7449951171875],
        # US-NE2
        # [-96.46994024736414, 41.16491226772292, 356.7449951171875],
        # US-NE3
        [-96.43968912903934, 41.17964494123755, 357.64923095703125],
    ]
)
def test_image_value_default(lon, lat, expected, tol=0.000001):
    input = ee.Image('USGS/NED')
    output = utils.image_value(input, xy=(lon, lat))['elevation']
    assert abs(output - expected) <= tol


@pytest.mark.parametrize(
    'lon, lat, expected',
    [
        # US-NE1
        [-96.47672812080845, 41.16506126041818, 356.7449951171875],
        # US-NE2
        # [-96.46994024736414, 41.16491226772292, 356.7449951171875],
        # US-NE3
        [-96.43968912903934, 41.17964494123755, 357.6092224121094],
    ]
)
def test_image_value_scale(lon, lat, expected, tol=0.000001):
    input = ee.Image('USGS/NED')
    output = utils.image_value(input, xy=(lon, lat), scale=0.1)['elevation']
    assert abs(output - expected) <= tol


@pytest.mark.parametrize(
    'lon, lat, expected',
    [
        # US-NE1
        [-96.47672812080845, 41.16506126041818, 356.7449951171875],
        # US-NE2
        # [-96.46994024736414, 41.16491226772292, 356.7449951171875],
        # US-NE3
        [-96.43968912903934, 41.17964494123755, 357.6092224121094],
    ]
)
def test_coll_value(lon, lat, expected, tol=0.000001):
    input = ee.ImageCollection([ee.Image('USGS/NED'), ee.Image('USGS/NED')])
    output = utils.coll_value(input, xy=(lon, lat))
    assert abs(output[1][4] - expected) <= tol
    assert abs(output[2][4] - expected) <= tol
