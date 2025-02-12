from datetime import datetime, timedelta, UTC

import ee
import pytest

import openet.disalexi as disalexi
import openet.disalexi.utils as utils
# TODO: import utils from openet.core
# import openet.core.utils as utils


COLL_ID = 'LANDSAT/LC08/C02/T1_L2'
SCENE_ID = 'LC08_044033_20170716'
IMAGE_ID = f'{COLL_ID}/{SCENE_ID}'
SCENE_TIME = 1500230731090
SCENE_DT = datetime.fromtimestamp(SCENE_TIME / 1000.0, UTC)
SCENE_DATE = SCENE_DT.strftime('%Y-%m-%d')
SCENE_DOY = int(SCENE_DT.strftime('%j'))
SCENE_HOUR = (
    SCENE_DT.hour + SCENE_DT.minute / 60.0 +
    (SCENE_DT.second + SCENE_DT.microsecond / 1000000.0) / 3600.0
)
SCENE_0UTC_DT = datetime.strptime(SCENE_DATE, '%Y-%m-%d')
SCENE_POINT = (-119.5, 36.0)
TEST_POINT = (-121.5265, 38.7399)


# DisALEXI needs a non-constant image since it computes lat/lon
# LAI and LST are read from sources and net set in the "prep" functions anymore
def default_image(
        albedo=0.125,
        cfmask=0,
        ndvi=0.875,
        # lai=4.2,
        # lst=306.5,
        coll_id=None,
        scene_id=None,
        scene_time=None,
):
    if coll_id is None:
        coll_id = COLL_ID
    if scene_id is None:
        scene_id = SCENE_ID
    if scene_time is None:
        scene_time = SCENE_TIME
    mask_img = ee.Image(f'{coll_id}/{scene_id}').select(['SR_B1']).multiply(0)
    return (
        ee.Image([mask_img.add(albedo), mask_img.add(cfmask), mask_img.add(ndvi)])
        .rename(['albedo', 'cfmask', 'ndvi'])
        .set({
            'system:index': scene_id,
            'system:time_start': scene_time,
            # 'system:time_start': ee.Date(SCENE_DATE).millis(),
            'system:id': f'{coll_id}/{scene_id}',
        })
    )

def default_image_args(
        albedo=0.125,
        cfmask=0, ndvi=0.875,
        ta_source='CONUS_V006',
        alexi_source='CONUS_V006',
        lai_source=4.2,
        lst_source=306.5,
        et_reference_source=10,
        et_reference_band=None,
        et_reference_factor=1.0,
        et_reference_resample='nearest',
        stability_iterations=25,
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
        ta_source='CONUS_V006',
        alexi_source='CONUS_V006',
        lai_source=4.2,
        lst_source=306.5,
        et_reference_source=10,
        et_reference_band=None,
        et_reference_factor=1.0,
        et_reference_resample='nearest',
        stability_iterations=25,
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


def test_Image_init_default_parameters():
    m = disalexi.Image(default_image())
    assert m.ta_source == 'CONUS_V006'
    assert m.alexi_source == 'CONUS_V006'
    assert m.lai_source == 'openet-landsat-lai'
    #assert m.lai_source == 'projects/openet/assets/lai/landsat/c02'
    assert m.lst_source == 'projects/openet/assets/lst/landsat/c02'
    assert m.elevation_source == 'USGS/SRTMGL1_003'
    assert m.landcover_source == 'projects/sat-io/open-datasets/USGS/ANNUAL_NLCD/LANDCOVER'
    assert m.air_pres_source == 'CFSR'
    assert m.air_temp_source == 'CFSR'
    assert m.rs_daily_source == 'CFSR'
    assert m.rs_hourly_source == 'CFSR'
    assert m.vapor_pres_source == 'CFSR'
    assert m.wind_speed_source == 'CFSR'
    assert m.stabil_iter is None
    assert m.albedo_iter == 10
    assert m.rs_interp_flag is True
    assert m.ta_smooth_flag is True
    assert m.et_reference_source is None
    assert m.et_reference_band is None
    assert m.et_reference_factor is None
    assert m.et_reference_resample is None


# Todo: Break these up into separate functions?
def test_Image_init_calculated_properties():
    d = default_image_obj()
    assert utils.getinfo(d.time_start) == SCENE_TIME


def test_Image_init_date_properties():
    d = default_image_obj()
    assert utils.getinfo(d.datetime)['value'] == SCENE_TIME
    assert utils.getinfo(d.start_date)['value'] == utils.millis(SCENE_0UTC_DT)
    assert utils.getinfo(d.end_date)['value'] == utils.millis(SCENE_0UTC_DT + timedelta(days=1))
    assert utils.getinfo(d.hour) == SCENE_HOUR
    assert utils.getinfo(d.hour_int) == int(SCENE_HOUR)


@pytest.mark.parametrize(
    'source, xy, expected',
    [
        # Direct & variable step assets
        ['projects/openet/assets/disalexi/tair/conus_v006', TEST_POINT, 299.68001353628875],
        ['projects/openet/assets/disalexi/tair/conus_v006', [-121.50822, 38.71776], 298.61358548716754],
        # 1k step assets
        ['projects/openet/assets/disalexi/tair/conus_v006_1k', TEST_POINT, 298.047307556939],
        ['projects/openet/assets/disalexi/tair/conus_v006_1k', [-121.50822, 38.71776], 297.32998269291346],
        # The openet legacy assets will eventually be removed
        ['projects/openet/disalexi/tair/conus_v006_1k', TEST_POINT, 298.047307556939],
        ['projects/openet/disalexi/tair/conus_v006_1k', [-121.50822, 38.71776], 297.32998269291346],
        # The CONUS_V006 keyword is currently pointed at the 1k steps assets
        ['CONUS_V006', TEST_POINT, 298.047307556939],
        ['CONUS_V006', [-121.50822, 38.71776], 297.32998269291346],
        [ee.Image('USGS/SRTMGL1_003').multiply(0).add(10), TEST_POINT, 10],
        ['294.8', TEST_POINT, 294.8],  # Check constant values
        [294.8, TEST_POINT, 294.8],    # Check constant values
    ]
)
def test_Image_ta_source(source, xy, expected, tol=0.01):
    m = disalexi.Image(default_image(), ta_source=source, ta_smooth_flag=True)
    output = utils.point_image_value(ee.Image(m.ta), xy)
    assert abs(output['ta'] - expected) <= tol


# CGM - This test should not have an identical value to the one above but the
#   smoothing radius is currently too small to change the data
def test_Image_ta_smooth_flag(tol=0.01):
    """The smooth flag defaults to True, so set False to test"""
    m = disalexi.Image(
        default_image(),
        ta_source='projects/openet/assets/disalexi/tair/conus_v006_1k',
        ta_smooth_flag=False
    )
    output = utils.point_image_value(ee.Image(m.ta), TEST_POINT)
    assert abs(output['ta'] - 290) <= tol


def test_Image_ta_interp_flag(tol=0.01):
    """The interp flag defaults to True, so set False to test"""
    m = disalexi.Image(
        default_image(),
        ta_source='projects/openet/assets/disalexi/tair/conus_v006_1k',
        ta_interp_flag=False
    )
    output = utils.point_image_value(ee.Image(m.ta), TEST_POINT)
    assert abs(output['ta'] - 298.047307556939) <= tol


def test_Image_ta_source_exception():
    with pytest.raises(ValueError):
        utils.getinfo(disalexi.Image(default_image(), ta_source='').ta)


def test_Image_ta_properties():
    """Test if properties are set on the ET image"""
    output = utils.getinfo(default_image_obj().ta)
    assert output['bands'][0]['id'] == 'ta'
    assert output['properties']['system:index'] == SCENE_ID
    assert output['properties']['system:time_start'] == SCENE_TIME
    assert output['properties']['image_id'] == IMAGE_ID


# CGM - Added the scene_id for this test to make it easy to switch the datetime
# I also modified the default_image to take the scene_id & scene_time parameters
@pytest.mark.parametrize(
    'scene_id, source, xy, expected',
    [
        # ALEXI ET is currently in MJ m-2 d-1
        ['LC08_044033_20200708', 'CONUS_V006', TEST_POINT, 17.999324798583984 * 0.408],
        ['LC08_044033_20200724', 'CONUS_V006', TEST_POINT, 17.819684982299805 * 0.408],
        [None, 'projects/ee-tulipyangyun-2/assets/alexi/ALEXI_V006', TEST_POINT, 12.765579223632812 * 0.408],
        [None, ee.Image('USGS/SRTMGL1_003').multiply(0).add(10), TEST_POINT, 10],
        [None, '10.382039', TEST_POINT, 10.382039],
        [None, 10.382039, TEST_POINT, 10.382039],
    ]
)
def test_Image_alexi_source(scene_id, source, xy, expected, tol=0.0001):
    if scene_id:
        scene_time = ee.Image(f'{COLL_ID}/{scene_id}').get('system:time_start').getInfo()
    else:
        scene_time = None
    m = disalexi.Image(default_image(scene_id=scene_id, scene_time=scene_time), alexi_source=source)
    output = utils.point_image_value(m.et_alexi, xy)
    assert abs(output['et_alexi'] - expected) <= tol


def test_Image_alexi_source_exception():
    with pytest.raises(ValueError):
        utils.getinfo(disalexi.Image(default_image(), alexi_source='').et_alexi)


def test_Image_alexi_band_name():
    output = utils.getinfo(disalexi.Image(default_image()).et_alexi)['bands'][0]['id']
    assert output == 'et_alexi'


@pytest.mark.parametrize(
    'source, xy, expected',
    [
        ['USGS/SRTMGL1_003', TEST_POINT, 3],
        [ee.Image('USGS/SRTMGL1_003'), TEST_POINT, 3],
        ['2364.351', TEST_POINT, 2364.351],
        [2364.351, TEST_POINT, 2364.351],
    ]
)
def test_Image_elevation_source(source, xy, expected, tol=0.001):
    output = utils.point_image_value(
        disalexi.Image(default_image(), elevation_source=source).elevation, xy
    )
    assert abs(output['elevation'] - expected) <= tol


# def test_Image_elevation_source_exception():
#     with pytest.raises(ValueError):
#         utils.getinfo(disalexi.Image(default_image(), elevation_source='').elevation)


def test_Image_elevation_band_name():
    output = utils.getinfo(disalexi.Image(default_image()).elevation)['bands'][0]['id']
    assert output == 'elevation'


@pytest.mark.parametrize(
    'source, xy, expected',
    [
        [
            'projects/sat-io/open-datasets/USGS/ANNUAL_NLCD/LANDCOVER/Annual_NLCD_LndCov_2021_CU_C1V0',
            TEST_POINT, 82
        ],
        ['projects/sat-io/open-datasets/USGS/ANNUAL_NLCD/LANDCOVER', TEST_POINT, 82],
        ['USGS/NLCD_RELEASES/2021_REL/NLCD/2021', TEST_POINT, 82],
        ['USGS/NLCD_RELEASES/2021_REL/NLCD', TEST_POINT, 82],
        ['USGS/NLCD_RELEASES/2019_REL/NLCD', TEST_POINT, 82],
        ['USGS/NLCD_RELEASES/2019_REL/NLCD/2016', TEST_POINT, 82],
        [ee.Image('USGS/SRTMGL1_003').multiply(0).add(82), TEST_POINT, 82],
        ['82', TEST_POINT, 82],
        [82, TEST_POINT, 82],
    ]
)
def test_Image_landcover_source(source, xy, expected, tol=0.001):
    output = utils.point_image_value(
        disalexi.Image(default_image(), landcover_source=source).landcover, xy
    )
    assert abs(output['landcover'] - expected) <= tol


@pytest.mark.parametrize(
    'source',
    [
        # No source
        '',
        # Fail on collection ID with trailing slash
        'USGS/NLCD_RELEASES/2019_REL/NLCD/',
    ]
)
def test_Image_landcover_source_invalid(source):
    with pytest.raises(ValueError):
        utils.getinfo(disalexi.Image(default_image(), landcover_source=source).landcover)


@pytest.mark.parametrize(
    'source',
    [
        # Fail on missing year images
        # Not sure why these tests aren't working
        'USGS/NLCD_RELEASES/2019_REL/NLCD/2021',
        'USGS/NLCD_RELEASES/2021_REL/NLCD/2019',
    ]
)
def test_Image_landcover_source_invalid(source):
    with pytest.raises(Exception):
        utils.getinfo(disalexi.Image(default_image(), landcover_source=source).landcover)


def test_Image_landcover_band_name():
    output = utils.getinfo(disalexi.Image(default_image()).landcover)['bands'][0]['id']
    assert output == 'landcover'


@pytest.mark.parametrize(
    'source, xy, expected',
    [
        ['CFSR', TEST_POINT, 100.48347318265635],
        ['ESTIMATE', TEST_POINT, 101.2395327760238],
        # ['ESTIMATE', TEST_POINT, 101.26454311195941],
        ['100.41653321557092', TEST_POINT, 100.41653321557092],
        [100.41653321557092, TEST_POINT, 100.41653321557092],
    ]
)
def test_Image_air_pressure_source(source, xy, expected, tol=0.0001):
    output = utils.point_image_value(
        disalexi.Image(default_image(), air_pres_source=source).air_pressure, xy
    )
    assert abs(output['air_pressure'] - expected) <= tol


def test_Image_air_pressure_source_exception():
    with pytest.raises(ValueError):
        utils.getinfo(disalexi.Image(default_image(), air_pres_source='').air_pressure)


def test_Image_air_pressure_band_name():
    output = utils.getinfo(disalexi.Image(default_image()).air_pressure)['bands'][0]['id']
    assert output == 'air_pressure'


@pytest.mark.parametrize(
    'source, xy, expected',
    [
        ['CFSR', TEST_POINT, 306.3861390180871],
        ['306', [-121.50822, 38.71776], 306],
        [306, [-121.50822, 38.71776], 306],
    ]
)
def test_Image_air_temperature_source(source, xy, expected, tol=0.0001):
    output = utils.point_image_value(
        disalexi.Image(default_image(), air_temp_source=source).air_temperature, xy
    )
    assert abs(output['air_temperature'] - expected) <= tol


def test_Image_air_temperature_source_exception():
    with pytest.raises(ValueError):
        utils.getinfo(disalexi.Image(default_image(), air_temp_source='').air_temperature)


def test_Image_air_temperature_band_name():
    output = utils.getinfo(disalexi.Image(default_image()).air_temperature)['bands'][0]['id']
    assert output == 'air_temperature'


@pytest.mark.parametrize(
    'source, xy, expected',
    [
        # CGM - I'm not sure why these two values are different
        #   The intermediate values are identical but end up different after smoothing
        ['CFSR', TEST_POINT, 8594.84901460842],
        [ee.Image('USGS/SRTMGL1_003').multiply(0).add(10), TEST_POINT, 10],
        ['8587.4091796875', TEST_POINT, 8587.4091796875],
        [8587.4091796875, TEST_POINT, 8587.4091796875],
    ]
)
def test_Image_rs_daily_source(source, xy, expected, tol=0.0001):
    output = utils.point_image_value(disalexi.Image(default_image(), rs_daily_source=source).rs24, xy)
    assert abs(output['rs'] - expected) <= tol


def test_Image_rs_daily_source_exception():
    with pytest.raises(ValueError):
        utils.getinfo(disalexi.Image(default_image(), rs_daily_source='').rs24)


def test_Image_rs_daily_band_name():
    output = utils.getinfo(disalexi.Image(default_image()).rs24)['bands'][0]['id']
    assert output == 'rs'


@pytest.mark.parametrize(
    'source, xy, expected',
    [
        ['CFSR', TEST_POINT, 936.7493916867356],
        [ee.Image('USGS/SRTMGL1_003').multiply(0).add(10), TEST_POINT, 10],
        ['946.6906', TEST_POINT, 946.6906],
        [946.6906, TEST_POINT, 946.6906],
    ]
)
def test_Image_rs_hourly_source(source, xy, expected, tol=0.0001):
    output = utils.point_image_value(disalexi.Image(default_image(), rs_hourly_source=source).rs1, xy)
    assert abs(output['rs'] - expected) <= tol


def test_Image_rs_hourly_interp(tol=0.001):
    output = utils.point_image_value(
        disalexi.Image(default_image(), rs_hourly_source='CFSR', rs_interp_flag=True).rs1, TEST_POINT
    )
    assert abs(output['rs'] - 936.7493916867356) <= tol


def test_Image_rs_hourly_no_interp(tol=0.001):
    output = utils.point_image_value(
        disalexi.Image( default_image(), rs_hourly_source='CFSR', rs_interp_flag=False).rs1, TEST_POINT
    )
    assert abs(output['rs'] - 872.6140747070312) <= tol


def test_Image_rs_hourly_source_exception():
    with pytest.raises(ValueError):
        utils.getinfo(disalexi.Image(default_image(), rs_hourly_source='').rs1)


def test_Image_rs_hourly_band_name():
    output = utils.getinfo(disalexi.Image(default_image()).rs1)['bands'][0]['id']
    assert output == 'rs'


@pytest.mark.parametrize(
    'source, xy, expected',
    [
        ['CFSR', TEST_POINT, 0.8393702479538856],
        ['2.001476', TEST_POINT, 2.001476],
        [2.001476, TEST_POINT, 2.001476],
    ]
)
def test_Image_vapor_pressure_source(source, xy, expected, tol=0.0001):
    output = utils.point_image_value(
        disalexi.Image(default_image(), vapor_pres_source=source).vapor_pressure, xy
    )
    assert abs(output['vapor_pressure'] - expected) <= tol


def test_Image_vapor_pressure_source_exception():
    with pytest.raises(ValueError):
        utils.getinfo(disalexi.Image(default_image(), vapor_pres_source='').vapor_pressure)


def test_Image_vapor_pressure_band_name():
    output = utils.getinfo(disalexi.Image(default_image()).vapor_pressure)['bands'][0]['id']
    assert output == 'vapor_pressure'


@pytest.mark.parametrize(
    'source, xy, expected',
    [
        ['CFSR', TEST_POINT, 2.814758509561662],
        [ee.Image('USGS/SRTMGL1_003').multiply(0).add(10), TEST_POINT, 10],
        ['2.001476', TEST_POINT, 2.001476],
        [2.001476, TEST_POINT, 2.001476],
    ]
)
def test_Image_wind_speed_source(source, xy, expected, tol=0.0001):
    output = utils.point_image_value(
        disalexi.Image(default_image(), wind_speed_source=source).wind_speed, xy
    )
    assert abs(output['wind_speed'] - expected) <= tol


@pytest.mark.parametrize(
    'source, xy, expected',
    [
        [1.9999, TEST_POINT, 2.0000],
        [20.0001, TEST_POINT, 20.0000],
    ]
)
def test_Image_wind_speed_clamping(source, xy, expected, tol=0.0001):
    output = utils.point_image_value(
        disalexi.Image(default_image(), wind_speed_source=source).wind_speed, xy
    )
    assert abs(output['wind_speed'] - expected) <= tol


def test_Image_wind_speed_source_exception():
    with pytest.raises(ValueError):
        utils.getinfo(disalexi.Image(default_image(), wind_speed_source='').wind_speed)


def test_Image_wind_speed_band_name():
    output = utils.getinfo(disalexi.Image(default_image()).wind_speed)['bands'][0]['id']
    assert output == 'wind_speed'


def test_Image_lai_band_name():
    output = utils.getinfo(disalexi.Image(default_image()).lai)['bands'][0]['id']
    assert output == 'lai'


@pytest.mark.parametrize(
    'source, xy, expected',
    [
        ['openet-landsat-lai', TEST_POINT, 3.6301],
        ['openet-landsat-lai', [-121.50822, 38.71776], 0.5750],
        # # DEADBEEF - These collections will be removed in the future
        # ['projects/openet/assets/lai/landsat/c02', TEST_POINT, 3.6301],
        # ['projects/openet/assets/lai/landsat/c02', [-121.50822, 38.71776], 0.5750],
        ['OPENET-LAI', TEST_POINT, 3.6301],
        ['4.123', TEST_POINT, 4.123],
        [4.123, TEST_POINT, 4.123],
    ]
)
def test_Image_lai_source(source, xy, expected, tol=0.0001):
    output = utils.point_image_value(disalexi.Image(default_image(), lai_source=source).lai, xy)
    assert abs(output['lai'] - expected) <= tol


@pytest.mark.parametrize(
    'source',
    [
        'openet-landsat-lai',
        # # DEADBEEF
        # 'projects/openet/assets/lai/landsat/c02',
        # 4.123,
    ]
)
def test_Image_lai_version_property(source):
    model_obj = disalexi.Image(default_image(), lai_source=source)
    # Trigger building the LAI in order to set the landsat_lai_version
    lai_img = model_obj.lai
    assert model_obj.landsat_lai_version


def test_Image_lst_band_name():
    output = utils.getinfo(disalexi.Image(default_image()).lst)['bands'][0]['id']
    assert output == 'lst'


@pytest.mark.parametrize(
    ' source, xy, expected',
    [
        ['projects/openet/assets/lst/landsat/c02', TEST_POINT, 304.5],
        ['projects/openet/assets/lst/landsat/c02', [-121.50822, 38.71776], 323.8],
        ['300', TEST_POINT, 300],
        [300, TEST_POINT, 300],
    ]
)
def test_Image_lst_source(source, xy, expected, tol=0.0001):
    output = utils.point_image_value(disalexi.Image(default_image(), lst_source=source).lst, xy)
    assert abs(output['lst'] - expected) <= tol


def test_Image_et_default_values(expected=5.6248832186825375, tol=0.0001):
    output = utils.point_image_value(default_image_obj().et, TEST_POINT)
    assert abs(output['et'] - expected) <= tol


def test_Image_et_properties():
    """Test if properties are set on the ET image"""
    output = utils.getinfo(default_image_obj().et)
    assert output['bands'][0]['id'] == 'et'
    assert output['properties']['system:index'] == SCENE_ID
    assert output['properties']['system:time_start'] == SCENE_TIME
    assert output['properties']['image_id'] == IMAGE_ID


def test_Image_et_reference_properties():
    """Test if properties are set on the ETr image"""
    output = utils.getinfo(default_image_obj().et_reference)
    assert output['bands'][0]['id'] == 'et_reference'
    assert output['properties']['system:index'] == SCENE_ID
    assert output['properties']['system:time_start'] == SCENE_TIME
    assert output['properties']['image_id'] == IMAGE_ID


def test_Image_et_reference_default_values(expected=10.0, tol=0.0001):
    output = utils.point_image_value(
        default_image_obj(et_reference_source=expected).et_reference, TEST_POINT
    )
    assert abs(output['et_reference'] - expected) <= tol


@pytest.mark.parametrize(
    'source, band, xy, expected',
    [
        ['IDAHO_EPSCOR/GRIDMET', 'etr', TEST_POINT, 11.2],
        ['projects/openet/assets/reference_et/conus/gridmet/daily/v1', 'etr', TEST_POINT, 9.3565],
        ['projects/openet/assets/reference_et/california/cimis/daily/v1', 'etr', TEST_POINT, 10.174],
        ['10.8985', None, TEST_POINT, 10.8985],
        [10.8985, None, TEST_POINT, 10.8985],
    ]
)
def test_Image_et_reference_source(source, band, xy, expected, tol=0.001):
    output = utils.point_image_value(
        default_image_obj(et_reference_source=source, et_reference_band=band).et_reference, xy
    )
    assert abs(output['et_reference'] - expected) <= tol


def test_Image_mask_values():
    # Mask is 1 for active pixels and nodata for cloudy pixels (cfmask >= 1)
    output = utils.point_image_value(default_image_obj(cfmask=1).mask, TEST_POINT)
    assert output['mask'] == None
    output = utils.point_image_value(default_image_obj(cfmask=0).mask, TEST_POINT)
    assert output['mask'] == 1


def test_Image_mask_properties():
    """Test if properties are set on the mask image"""
    output = utils.getinfo(default_image_obj().mask)
    assert output['bands'][0]['id'] == 'mask'
    assert output['properties']['system:index'] == SCENE_ID
    assert output['properties']['system:time_start'] == SCENE_TIME
    assert output['properties']['image_id'] == IMAGE_ID


# TODO: Debug "Computation is too complex." error
# def test_Image_time_point_values():
#     output = utils.point_image_value(default_image_obj().time, TEST_POINT)
#     assert output['time'] == utils.millis(SCENE_0UTC_DT)


def test_Image_time_properties():
    """Test if properties are set on the time image"""
    output = utils.getinfo(default_image_obj().time)
    assert output['bands'][0]['id'] == 'time'
    assert output['properties']['system:index'] == SCENE_ID
    assert output['properties']['system:time_start'] == SCENE_TIME
    assert output['properties']['image_id'] == IMAGE_ID


def test_Image_calculate_variables_default():
    output = utils.getinfo(default_image_obj().calculate())
    assert {x['id'] for x in output['bands']} == {'et'}


@pytest.mark.parametrize(
    'variables, expected',
    [
        [['ndvi'], ['ndvi']],
        [['lst'], ['lst']],
        [['ndvi', 'lst'], ['ndvi', 'lst']],
        [['mask'], ['mask']],
        [['et_reference'], ['et_reference']],
        # [['et_fraction'], ['et_fraction']],
    ]
)
def test_Image_calculate_variables_custom(variables, expected):
    output = utils.getinfo(default_image_obj().calculate(variables))
    assert {x['id'] for x in output['bands']} == set(expected)


def test_Image_calculate_variables_all():
    variables = {'et', 'mask', 'ndvi', 'time'}
    output = utils.getinfo(default_image_obj().calculate(variables=list(variables)))
    assert {x['id'] for x in output['bands']} == variables


def test_Image_calculate_properties():
    """Test if properties are set on the output image"""
    output = utils.getinfo(default_image_obj().calculate(['ndvi']))
    assert output['properties']['system:index'] == SCENE_ID
    assert output['properties']['system:time_start'] == SCENE_TIME
    assert output['properties']['image_id'] == IMAGE_ID


def test_Image_calculate_values(tol=0.0001):
    """Test if the calculate method returns values"""
    output_img = disalexi.Image(
            default_image(albedo=0.125, cfmask=0, ndvi=0.875),
            ta_source=300, alexi_source=10, lai_source=4.7, lst_source=306,
            elevation_source=10, landcover_source=82,
            rs_daily_source=8600, rs_hourly_source=950, wind_speed_source=3.25,
            stability_iterations=6, albedo_iterations=3)\
        .calculate(['et'])
    #     .calculate(['et', 'et_reference', 'et_fraction'])
    output = utils.point_image_value(output_img, TEST_POINT)
    assert abs(output['et'] - 7.392495155334473) <= tol
    # assert abs(output['et_reference'] - 10) <= tol
    # assert abs(output['et_fraction'] - 0.58) <= tol


def test_Image_calculate_variables_valueerror():
    """Test if calculate method raises a valueerror for invalid variables"""
    with pytest.raises(ValueError):
        utils.getinfo(default_image_obj().calculate(['FOO']))


# Landsat Collection 2 Level 2
def test_Image_from_landsat_c02_l2_default_image():
    """Test that the classmethod is returning a class object"""
    output = disalexi.Image.from_landsat_c02_l2(ee.Image(IMAGE_ID))
    assert type(output) == type(disalexi.Image(default_image()))


@pytest.mark.parametrize(
    'image_id',
    [
        'LANDSAT/LC09/C02/T1_L2/LC09_044033_20220127',
        'LANDSAT/LC08/C02/T1_L2/LC08_044033_20170716',
        'LANDSAT/LE07/C02/T1_L2/LE07_044033_20170708',
        'LANDSAT/LT05/C02/T1_L2/LT05_044033_20110716',
    ]
)
def test_Image_from_landsat_c02_l2_image_id(image_id):
    """Test instantiating the class from a Landsat image ID"""
    output = utils.getinfo(disalexi.Image.from_landsat_c02_l2(image_id).ndvi)
    assert output['properties']['system:index'] == image_id.split('/')[-1]


def test_Image_from_landsat_c02_l2_image():
    """Test instantiating the class from a Landsat ee.Image"""
    image_id = 'LANDSAT/LC08/C02/T1_L2/LC08_044033_20170716'
    output = utils.getinfo(disalexi.Image.from_landsat_c02_l2(ee.Image(image_id)).ndvi)
    assert output['properties']['system:index'] == image_id.split('/')[-1]


def test_Image_from_landsat_c02_l2_et():
    """Test if ET can be built for a Landsat images"""
    image_id = 'LANDSAT/LC08/C02/T1_L2/LC08_044033_20170716'
    output = utils.getinfo(disalexi.Image.from_landsat_c02_l2(image_id).et)
    assert output['properties']['system:index'] == image_id.split('/')[-1]


def test_Image_from_landsat_c02_l2_exception():
    """Test instantiating the class for an invalid image ID"""
    with pytest.raises(Exception):
        utils.getinfo(disalexi.Image.from_landsat_c02_l2(ee.Image('FOO')).ndvi)


def test_Image_from_landsat_c02_l2_scaling():
    """Test if Landsat SR images images are being scaled"""
    sr_img = ee.Image('LANDSAT/LC08/C02/T1_L2/LC08_044033_20170716')
    input_img = (
        ee.Image.constant([10909, 10909, 10909, 14545, 10909, 10909, 44177.6, 21824, 0])
        .rename(['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7', 'ST_B10', 'QA_PIXEL', 'QA_RADSAT'])
        .set({'SPACECRAFT_ID': ee.String(sr_img.get('SPACECRAFT_ID')),
              'system:id': ee.String(sr_img.get('system:id')),
              'system:index': ee.String(sr_img.get('system:index')),
              'system:time_start': ee.Number(sr_img.get('system:time_start'))})
    )

    output = utils.constant_image_value(disalexi.Image.from_landsat_c02_l2(input_img).ndvi)
    assert abs(output['ndvi'] - 0.333) <= 0.001

    output = utils.constant_image_value(disalexi.Image.from_landsat_c02_l2(input_img).albedo)
    assert abs(output['albedo'] - 0.137) <= 0.001


@pytest.mark.parametrize(
    'image_id',
    [
        'LANDSAT/LC08/C02/T1_L2/LC08_044033_20170716',
    ]
)
def test_Image_from_image_id(image_id):
    """Test instantiating the class using the from_image_id method"""
    output = utils.getinfo(disalexi.Image.from_image_id(image_id).ndvi)
    assert output['properties']['system:index'] == image_id.split('/')[-1]


def test_Image_from_method_kwargs():
    """Test that the init parameters can be passed through the helper methods"""
    assert disalexi.Image.from_landsat_c02_l2(
        'LANDSAT/LC08/C02/T1_L2/LC08_042035_20150713',
        elevation_source='DEADBEEF').elevation_source == 'DEADBEEF'
