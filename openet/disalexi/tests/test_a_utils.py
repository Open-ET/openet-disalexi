import datetime

import ee
import pytest

import openet.disalexi.utils as utils


def test_getinfo():
    assert utils.getinfo(ee.Number(1)) == 1


def test_constant_image_value(tol=0.000001):
    expected = 10.123456789
    input_img = ee.Image.constant(expected)
    output = utils.constant_image_value(input_img)
    assert abs(output['constant'] - expected) <= tol


def test_point_image_value(tol=0.001):
    expected = 2362
    output = utils.point_image_value(ee.Image('USGS/SRTMGL1_003'), [-106.03249, 37.17777])
    assert abs(output['elevation'] - expected) <= tol


def point_coll_value(tol=0.001):
    expected = 2362
    output = utils.point_coll_value(
        ee.ImageCollection([ee.Image('USGS/SRTMGL1_003')]), [-106.03249, 37.17777]
    )
    assert abs(output['elevation']['2012-04-04'] - expected) <= tol


@pytest.mark.parametrize(
    # Note: These are made up values
    'input, expected',
    [
        ['True', True],
        ['true', True],
        ['t', True],
        ['False', False],
        ['false', False],
        ['f', False],
        [True, True],
        [False, False],
    ]
)
def test_boolean(input, expected):
    assert utils.boolean(input) == expected


def test_boolean_exception():
    with pytest.raises(ValueError):
        utils.boolean('DEADBEEF')


@pytest.mark.parametrize(
    'input, expected',
    [
        ['2015-07-13T18:33:39', 1436745600000],
        ['2015-07-13T00:00:00', 1436745600000],

    ]
)
def test_date_to_time_0utc(input, expected):
    assert utils.getinfo(utils.date_to_time_0utc(ee.Date(input))) == expected


# TODO: Test a few more parameterization and timesteps
def test_interpolate():
    a_dt = ee.Date('2020-07-01').advance(10, 'hour')
    b_dt = ee.Date('2020-07-01').advance(11, 'hour')
    coll = ee.ImageCollection([
        ee.Image.constant(10).set('system:time_start', a_dt.millis()),
        ee.Image.constant(20).set('system:time_start', b_dt.millis()),
    ])
    output = utils.constant_image_value(utils.interpolate(coll, a_dt.advance(0.5, 'hour'), timestep=1))
    assert output['constant'] == 15


@pytest.mark.parametrize(
    # Note: These are made up values
    'value, expected',
    [
        [300, True],
        ['300', True],
        [300.25, True],
        ['300.25', True],
        ['a', False],
    ]
)
def test_is_number(value, expected):
    assert utils.is_number(value) == expected


def test_millis():
    assert utils.millis(datetime.datetime(2015, 7, 13)) == 1436745600000


def test_valid_date():
    assert utils.valid_date('2015-07-13') == True
    assert utils.valid_date('2015-02-30') == False
    assert utils.valid_date('20150713') == False
    assert utils.valid_date('07/13/2015') == False
    assert utils.valid_date('07-13-2015', '%m-%d-%Y') == True
