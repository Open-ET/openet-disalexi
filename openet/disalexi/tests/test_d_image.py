import datetime

import ee
import pytest

import openet.disalexi as disalexi
import openet.disalexi.utils as utils
# TODO: import utils from openet.core
# import openet.core.utils as utils


COLL_ID = 'LANDSAT/LC08/C02/T1_L2/'
SCENE_ID = 'LC08_044033_20170716'
SCENE_TIME = 1500230731090
# SCENE_ID = 'LC08_042035_20150713'
# SCENE_TIME = 1436812419150
SCENE_DT = datetime.datetime.utcfromtimestamp(SCENE_TIME / 1000.0)
# SCENE_DT = datetime.datetime.strptime(SCENE_ID[-8:], '%Y%m%d')
SCENE_DATE = SCENE_DT.strftime('%Y-%m-%d')
SCENE_DOY = int(SCENE_DT.strftime('%j'))
SCENE_HOUR = (SCENE_DT.hour + SCENE_DT.minute / 60.0 +
              (SCENE_DT.second + SCENE_DT.microsecond / 1000000.0) / 3600.0)
# SCENE_HOUR = 18 + 33.0/60 + 39.15/3600
# SCENE_TIME = utils.millis(SCENE_DT)
SCENE_0UTC_DT = datetime.datetime.strptime(SCENE_DATE, '%Y-%m-%d')
SCENE_POINT = (-119.5, 36.0)
TEST_POINT = (-121.5265, 38.7399)
# TEST_POINT = (-19.44252382373145, 36.04047742246546)

# SCENE_TIME = utils.getinfo(ee.Date(SCENE_DATE).millis())
# SCENE_POINT = (-19.44252382373145, 36.04047742246546)
# SCENE_POINT = utils.getinfo(
#     ee.Image(COLL_ID + SCENE_ID).geometry().centroid())['coordinates']


# # Should these be test fixtures instead?
# # I'm not sure how to make them fixtures and allow input parameters
# def toa_image(red=0.1, nir=0.9, bt=305):
#     """Construct a fake Landsat 8 TOA image with renamed bands"""
#     return ee.Image.constant([red, nir, bt])\
#         .rename(['red', 'nir', 'lst']) \
#         .set({
#             'system:time_start': ee.Date(SCENE_DATE).millis(),
#             'k1_constant': ee.Number(607.76),
#             'k2_constant': ee.Number(1260.56),
#         })


# DisALEXI needs a non-constant image since it compute lat/lon
def default_image(albedo=0.125, cfmask=0, lai=4.2, lst=306.5, ndvi=0.875,
                  coll_id=None, scene_id=None, scene_time=None):
    if coll_id is None:
        coll_id = COLL_ID
    if scene_id is None:
        scene_id = SCENE_ID
    if scene_time is None:
        scene_time = SCENE_TIME
    mask_img = ee.Image(coll_id + scene_id).select(['SR_B1']).multiply(0)
    # mask_img = ee.Image('projects/disalexi/ta/CONUS_V001/{}_0'.format(SCENE_ID))\
    #     .select(0).multiply(0)
    return ee.Image([mask_img.add(albedo), mask_img.add(cfmask),
                     mask_img.add(lai), mask_img.add(lst),
                     mask_img.add(ndvi)])\
        .rename(['albedo', 'cfmask', 'lai', 'lst', 'ndvi'])\
        .set({
            'system:index': scene_id,
            'system:time_start': scene_time,
            # 'system:time_start': ee.Date(SCENE_DATE).millis(),
            'system:id': coll_id + scene_id,
        })

def default_image_args(albedo=0.125, cfmask=0, lai=4.2, lst=306.5, ndvi=0.875,
                       ta_source='CONUS_V005', alexi_source='CONUS_V005',
                       ta_smooth_flag=False,
                       et_reference_source=10, et_reference_band=None,
                       et_reference_factor=1.0, et_reference_resample='nearest',
                       stability_iterations=25, albedo_iterations=10,
                       ):
    return {
        'image': default_image(albedo=albedo, cfmask=cfmask, lai=lai, lst=lst,
                               ndvi=ndvi),
        'ta_source': ta_source,
        'alexi_source': alexi_source,
        'et_reference_source': et_reference_source,
        'et_reference_band': et_reference_band,
        'et_reference_factor': et_reference_factor,
        'et_reference_resample': et_reference_resample,
        'stability_iterations': stability_iterations,
        'albedo_iterations': albedo_iterations,
    }

def default_image_obj(albedo=0.125, cfmask=0, lai=4.2, lst=306.5, ndvi=0.875,
                      ta_source='CONUS_V005', alexi_source='CONUS_V005',
                      et_reference_source=10, et_reference_band=None,
                      et_reference_factor=1.0, et_reference_resample='nearest',
                      stability_iterations=25, albedo_iterations=10,
                      ):
    return disalexi.Image(**default_image_args(
        albedo=albedo, cfmask=cfmask, lai=lai, lst=lst, ndvi=ndvi,
        ta_source=ta_source,
        alexi_source=alexi_source,
        et_reference_source=et_reference_source,
        et_reference_band=et_reference_band,
        et_reference_factor=et_reference_factor,
        et_reference_resample=et_reference_resample,
        stability_iterations=stability_iterations,
        albedo_iterations=albedo_iterations,
    ))


def test_Image_init_default_parameters():
    m = disalexi.Image(default_image())
    assert m.ta_source == 'CONUS_V005'
    assert m.alexi_source == 'CONUS_V005'
    assert m.lai_source == 'projects/earthengine-legacy/assets/projects/openet/lai/landsat/c02'
    assert m.tir_source == 'projects/earthengine-legacy/assets/projects/openet/tir/landsat/c02'
    assert m.elevation_source == 'USGS/SRTMGL1_003'
    assert m.landcover_source == 'USGS/NLCD_RELEASES/2019_REL/NLCD'
    assert m.airpressure_source == 'CFSR'
    assert m.rs_daily_source == 'CFSR'
    assert m.rs_hourly_source == 'CFSR'
    assert m.ta0_source == 'CFSR'
    assert m.vp_source == 'CFSR'
    assert m.windspeed_source == 'CFSR'
    # assert m.stabil_iter == 36
    assert m.albedo_iter == 10
    assert m.rs_interp_flag == True
    # assert m.ta_interp_flag == True
    assert m.ta_smooth_flag == True
    assert m.et_reference_source == None
    assert m.et_reference_band == None
    assert m.et_reference_factor == None
    assert m.et_reference_resample == None


# Todo: Break these up into separate functions?
def test_Image_init_calculated_properties():
    d = default_image_obj()
    assert utils.getinfo(d.time_start) == SCENE_TIME
    # assert utils.getinfo(d.scene_id) == SCENE_ID
    # assert utils.getinfo(d.wrs2_tile) == 'p{}r{}'.format(
    #     SCENE_ID.split('_')[1][:3], SCENE_ID.split('_')[1][3:])


def test_Image_init_date_properties():
    d = default_image_obj()
    assert utils.getinfo(d.datetime)['value'] == SCENE_TIME
    # assert utils.getinfo(d.year) == int(SCENE_DATE.split('-')[0])
    # assert utils.getinfo(d.month) == int(SCENE_DATE.split('-')[1])
    assert utils.getinfo(d.start_date)['value'] == utils.millis(SCENE_0UTC_DT)
    assert utils.getinfo(d.end_date)['value'] == utils.millis(
        SCENE_0UTC_DT + datetime.timedelta(days=1))
    # assert utils.getinfo(d.doy) == SCENE_DOY
    # assert utils.getinfo(d.cycle_day) == int(
    #     (SCENE_DT - datetime.datetime(1970, 1, 3)).days % 8 + 1)
    assert utils.getinfo(d.hour) == SCENE_HOUR
    assert utils.getinfo(d.hour_int) == int(SCENE_HOUR)


# def test_Image_init_scene_id_property():
#     """Test that the system:index from a merged collection is parsed"""
#     input_img = default_image()
#     m = disalexi.Image(input_img.set('system:index', '1_2_' + SCENE_ID))
#     assert utils.getinfo(m.scene_id) == SCENE_ID


# CGM - These aren't lazy properties and don't have properties being set on them
# def test_Image_ndvi_properties():
#     """Test if properties are set on the NDVI image"""
#     output = utils.getinfo(disalexi.Image(default_image()).ndvi)
#     assert output['bands'][0]['id'] == 'ndvi'
#     assert output['properties']['system:index'] == SCENE_ID
#     assert output['properties']['system:time_start'] == SCENE_TIME
#     assert output['properties']['image_id'] == COLL_ID + SCENE_ID
#
#
# def test_Image_lst_properties():
#     """Test if properties are set on the LST image"""
#     output = utils.getinfo(disalexi.Image(default_image()).lst)
#     assert output['bands'][0]['id'] == 'lst'
#     assert output['properties']['system:index'] == SCENE_ID
#     assert output['properties']['system:time_start'] == SCENE_TIME
#     assert output['properties']['image_id'] == COLL_ID + SCENE_ID


@pytest.mark.parametrize(
    'source, xy, expected',
    [
        ['CONUS_V005', TEST_POINT, 297.4977364906148],
        ['CONUS_V005', [-121.50822, 38.71776], 296.99708009212264],
        # ['CONUS_V004', TEST_POINT, 298.1013811163322],
        # ['CONUS_V004', [-121.50822, 38.71776], 300.04066476402267],
        # ['CONUS_V003', TEST_POINT, 300.99823651442586],
        # ['CONUS_V003', [-121.50822, 38.71776], 300.4560056192083],
        # ['projects/openet/disalexi/tair/conus_v003_1k', TEST_POINT, 300.99823651442586],
        # ['projects/earthengine-legacy/assets/projects/openet/disalexi/tair/conus_v003_1k',
        # CGM - I'm not sure why the test values below don't work anymore when
        #   Ta shouldn't have changed.
        # ['CONUS_V003', TEST_POINT, 298.1013811163322],
        # ['CONUS_V003', [-121.50822, 38.71776], 300.04066476402267],
        # TEST_POINT, 300.99823651442586],
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
    m = disalexi.Image(default_image(), ta_source='CONUS_V005',
                       ta_smooth_flag=False)
    output = utils.point_image_value(ee.Image(m.ta), TEST_POINT)
    assert abs(output['ta'] - 290) <= tol
    # assert abs(output['ta'] - 298.1013811163322) <= tol


def test_Image_ta_interp_flag(tol=0.01):
    """The interp flag defaults to True, so set False to test"""
    m = disalexi.Image(default_image(), ta_source='CONUS_V005',
                       ta_interp_flag=False)
    output = utils.point_image_value(ee.Image(m.ta), TEST_POINT)
    assert abs(output['ta'] - 297.5) <= tol


def test_Image_ta_source_exception():
    with pytest.raises(ValueError):
        utils.getinfo(disalexi.Image(default_image(), ta_source='').ta)


def test_Image_ta_properties():
    """Test if properties are set on the ET image"""
    output =  utils.getinfo(default_image_obj().ta)
    assert output['bands'][0]['id'] == 'ta'
    assert output['properties']['system:index'] == SCENE_ID
    assert output['properties']['system:time_start'] == SCENE_TIME
    assert output['properties']['image_id'] == COLL_ID + SCENE_ID
   # assert output['properties']['ta_iteration'] > 0


# CGM - Added the scene_id for this test since the ALEXI_V004 is only currently
#   loaded for 2020 and 2021 and I needed a quick way to switch the datetime.
# I also modified the default_image to take the scene_id & scene_time parameters
# This could all be removed once ALEXI V004 is loaded for more years.
@pytest.mark.parametrize(
    'scene_id, source, xy, expected',
    [
        # ALEXI ET is currently in MJ m-2 d-1
        ['LC08_044033_20200708', 'CONUS_V005', TEST_POINT, 17.999324798583984 * 0.408],
        ['LC08_044033_20200724', 'CONUS_V005', TEST_POINT, 17.819684982299805 * 0.408],
        [None, 'CONUS_V005', TEST_POINT, 12.765579223632812 * 0.408],
        [None, 'projects/ee-tulipyangyun-2/assets/alexi/ALEXI_V005',
         TEST_POINT, 12.765579223632812 * 0.408],
        ['LC08_044033_20200708', 'CONUS_V004', TEST_POINT, 16.58306312561035 * 0.408],
        ['LC08_044033_20200724', 'CONUS_V004', TEST_POINT, 16.664167404174805 * 0.408],
        ['LC08_044033_20200708', 'CONUS_V003', TEST_POINT, 15.465087890625 * 0.408],
        ['LC08_044033_20200724', 'CONUS_V003', TEST_POINT, 14.467782974243164 * 0.408],
        ['LC08_044033_20170716', 'CONUS_V003', TEST_POINT, 12.03880786895752 * 0.408],
        [None, 'projects/disalexi/alexi/CONUS_V003', TEST_POINT, 12.03880786895752 * 0.408],
        [None, 'projects/earthengine-legacy/assets/projects/disalexi/alexi/CONUS_V003',
         TEST_POINT, 12.03880786895752 * 0.408],
        [None, ee.Image('USGS/SRTMGL1_003').multiply(0).add(10), TEST_POINT, 10],
        [None, '10.382039', TEST_POINT, 10.382039],
        [None, 10.382039, TEST_POINT, 10.382039],
    ]
)
def test_Image_alexi_source(scene_id, source, xy, expected, tol=0.0001):
    if scene_id:
        scene_time = ee.Image(COLL_ID + scene_id).get('system:time_start').getInfo()
    else:
        scene_time = None
    output = utils.point_image_value(disalexi.Image(
        default_image(scene_id=scene_id, scene_time=scene_time),
        alexi_source=source).et_alexi, xy)
    print(output['et_alexi']/0.408)
    assert abs(output['et_alexi'] - expected) <= tol


def test_Image_alexi_source_exception():
    with pytest.raises(ValueError):
        utils.getinfo(disalexi.Image(default_image(), alexi_source='').et_alexi)


def test_Image_alexi_band_name():
    output = utils.getinfo(
        disalexi.Image(default_image()).et_alexi)['bands'][0]['id']
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
    output = utils.point_image_value(disalexi.Image(
        default_image(), elevation_source=source).elevation, xy)
    assert abs(output['elevation'] - expected) <= tol


# def test_Image_elevation_source_exception():
#     with pytest.raises(ValueError):
#         utils.getinfo(disalexi.Image(
#             default_image(), elevation_source='').elevation)


def test_Image_elevation_band_name():
    output = utils.getinfo(
        disalexi.Image(default_image()).elevation)['bands'][0]['id']
    assert output == 'elevation'


@pytest.mark.parametrize(
    'source, xy, expected',
    [
        ['USGS/NLCD_RELEASES/2019_REL/NLCD', TEST_POINT, 82],
        ['USGS/NLCD_RELEASES/2019_REL/NLCD/2016', TEST_POINT, 82],
        ['USGS/NLCD_RELEASES/2016_REL', TEST_POINT, 82],
        ['USGS/NLCD_RELEASES/2016_REL/2016', TEST_POINT, 82],
        ['USGS/NLCD/NLCD2016', TEST_POINT, 82],
        ['GLOBELAND30', TEST_POINT, 10],
        [ee.Image('USGS/SRTMGL1_003').multiply(0).add(82), TEST_POINT, 82],
        ['82', TEST_POINT, 82],
        [82, TEST_POINT, 82],
    ]
)
def test_Image_landcover_source(source, xy, expected, tol=0.001):
    output = utils.point_image_value(disalexi.Image(
        default_image(), landcover_source=source).lc_source, xy)
    assert abs(output['landcover'] - expected) <= tol


def test_Image_landcover_source_exception():
    with pytest.raises(ValueError):
        utils.getinfo(disalexi.Image(
            default_image(), landcover_source='').landcover)


def test_Image_landcover_band_name():
    output = utils.getinfo(
        disalexi.Image(default_image()).lc_source)['bands'][0]['id']
    assert output == 'landcover'


@pytest.mark.parametrize(
    'source, xy, expected',
    [
        ['CFSR_NEW', TEST_POINT, 306.3861390180871],
        ['CFSR', TEST_POINT, 306.3861390180871],
        ['CFSR', [-121.5265, 38.7399], 306.3861390180871],
        ['CFSR', [-121.50822, 38.71776], 306.3840757488451],
        ['306', [-121.50822, 38.71776], 306],
        [306, [-121.50822, 38.71776], 306],
    ]
)
def test_Image_ta0_source(source, xy, expected, tol=0.0001):
    output = utils.point_image_value(disalexi.Image(
        default_image(), ta0_source=source).t_air0, xy)
    assert abs(output['tair0'] - expected) <= tol


def test_Image_ta0_source_exception():
    with pytest.raises(ValueError):
        utils.getinfo(disalexi.Image(
            default_image(), ta0_source='').t_air0)


def test_Image_t_air0_band_name():
    output = utils.getinfo(disalexi.Image(default_image()).t_air0)['bands'][0]['id']
    assert output == 'tair0'


@pytest.mark.parametrize(
    'source, xy, expected',
    [
        ['CFSR_NEW', TEST_POINT, 100.41703705955021],
        ['CFSR', TEST_POINT, 100.41703705955021],
        ['ESTIMATE', TEST_POINT, 101.26454311195941],
        ['100.41653321557092', TEST_POINT, 100.41653321557092],
        [100.41653321557092, TEST_POINT, 100.41653321557092],
    ]
)
def test_Image_pressure_source(source, xy, expected, tol=0.0001):
    output = utils.point_image_value(disalexi.Image(
        default_image(), airpressure_source=source).pressure, xy)
    assert abs(output['pressure'] - expected) <= tol


def test_Image_pressure_source_exception():
    with pytest.raises(ValueError):
        utils.getinfo(disalexi.Image(
            default_image(), airpressure_source='').pressure)


def test_Image_pressure_band_name():
    output = utils.getinfo(disalexi.Image(default_image()).pressure)['bands'][0]['id']
    assert output == 'pressure'


@pytest.mark.parametrize(
    'source, xy, expected',
    [
        ['CFSR_NEW', TEST_POINT, 8594.84901460842],
        ['CFSR', TEST_POINT, 8589.374265816517],
        ['MERRA2', TEST_POINT, 8587.4091796875],
        [ee.Image('USGS/SRTMGL1_003').multiply(0).add(10), TEST_POINT, 10],
        ['8587.4091796875', TEST_POINT, 8587.4091796875],
        [8587.4091796875, TEST_POINT, 8587.4091796875],
    ]
)
def test_Image_rs_daily_source(source, xy, expected, tol=0.0001):
    output = utils.point_image_value(disalexi.Image(
        default_image(), rs_daily_source=source).rs24, xy)
    assert abs(output['rs'] - expected) <= tol


def test_Image_rs_daily_source_exception():
    with pytest.raises(ValueError):
        utils.getinfo(disalexi.Image(
            default_image(), rs_daily_source='').rs24)


def test_Image_rs_daily_band_name():
    output = utils.getinfo(disalexi.Image(default_image()).rs24)['bands'][0]['id']
    assert output == 'rs'


@pytest.mark.parametrize(
    'source, xy, expected',
    [
        ['CFSR_NEW', TEST_POINT, 936.7493916867359],
        ['CFSR', TEST_POINT, 936.7493916867359],
        # ['MERRA2', TEST_POINT, 946.6906],
        [ee.Image('USGS/SRTMGL1_003').multiply(0).add(10), TEST_POINT, 10],
        ['946.6906', TEST_POINT, 946.6906],
        [946.6906, TEST_POINT, 946.6906],
    ]
)
def test_Image_rs_hourly_source(source, xy, expected, tol=0.0001):
    output = utils.point_image_value(disalexi.Image(
        default_image(), rs_hourly_source=source).rs1, xy)
    assert abs(output['rs'] - expected) <= tol


def test_Image_rs_hourly_interp(tol=0.001):
    output = utils.point_image_value(disalexi.Image(
        default_image(), rs_hourly_source='CFSR',
        rs_interp_flag=True).rs1, TEST_POINT)
    assert abs(output['rs'] - 936.7493916867359) <= tol


def test_Image_rs_hourly_no_interp(tol=0.001):
    output = utils.point_image_value(disalexi.Image(
        default_image(), rs_hourly_source='CFSR',
        rs_interp_flag=False).rs1, TEST_POINT)
    assert abs(output['rs'] - 936.7493916867359) <= tol


def test_Image_rs_hourly_source_exception():
    with pytest.raises(ValueError):
        utils.getinfo(disalexi.Image(default_image(), rs_hourly_source='').rs1)


def test_Image_rs_hourly_band_name():
    output = utils.getinfo(disalexi.Image(default_image()).rs1)['bands'][0]['id']
    assert output == 'rs'


@pytest.mark.parametrize(
    'source, xy, expected',
    [
        ['CFSR_NEW', TEST_POINT, 0.8393702479538856],
        ['CFSR', TEST_POINT, 0.8393702479538856],
        ['2.001476', TEST_POINT, 2.001476],
        [2.001476, TEST_POINT, 2.001476],
    ]
)
def test_Image_vp_source(source, xy, expected, tol=0.0001):
    output = utils.point_image_value(disalexi.Image(
        default_image(), vp_source=source).vp, xy)
    assert abs(output['vp'] - expected) <= tol


def test_Image_vp_source_exception():
    with pytest.raises(ValueError):
        utils.getinfo(disalexi.Image(default_image(), vp_source='').vp)


def test_Image_vp_band_name():
    output = utils.getinfo(disalexi.Image(default_image()).vp)['bands'][0]['id']
    assert output == 'vp'


@pytest.mark.parametrize(
    'source, xy, expected',
    [
        ['CFSR_NEW', TEST_POINT, 2.433912],
        ['CFSR', TEST_POINT, 2.433912],
        ['CFSV2', TEST_POINT, 2.169500185267768],
        [ee.Image('USGS/SRTMGL1_003').multiply(0).add(10), TEST_POINT, 10],
        ['2.001476', TEST_POINT, 2.001476],
        [2.001476, TEST_POINT, 2.001476],
    ]
)
def test_Image_windspeed_source(source, xy, expected, tol=0.0001):
    output = utils.point_image_value(disalexi.Image(
        default_image(), windspeed_source=source).windspeed, xy)
    assert abs(output['windspeed'] - expected) <= tol


def test_Image_windspeed_source_exception():
    with pytest.raises(ValueError):
        utils.getinfo(disalexi.Image(default_image(), windspeed_source='').windspeed)


def test_Image_windspeed_band_name():
    output = utils.getinfo(disalexi.Image(default_image()).windspeed)['bands'][0]['id']
    assert output == 'windspeed'


def test_Image_et_default_values(expected=5.497461682063846, tol=0.0001):
    output = utils.point_image_value(default_image_obj().et, TEST_POINT)
    assert abs(output['et'] - expected) <= tol


# CGM - What is this test doing?
# It seems be very sensitive to the number of iterations
#   but I can't set it higher than ~25
@pytest.mark.parametrize(
    'ta, iterations, expected',
    [
        # [298.427, 10, 7.602516963336294],
        # [298.427, 20, 6.750980231303747],
        [298.427, 25, 5.900454963186251],
        [300, 25, 6.7716727187269585],
    ]
)
def test_Image_et_fixed_source(ta, iterations, expected, tol=0.001):
    m = disalexi.Image(
        default_image(), ta_source=ta, stability_iterations=iterations)
    output = utils.point_image_value(m.et, TEST_POINT)
    assert abs(output['et'] - expected) <= tol


@pytest.mark.parametrize(
    'albedo, cfmask, lai, lst, ndvi, '
    'ta, alexi, elevation, landcover, '
    'ta0, rs_daily, rs_hourly, windspeed, '
    'stabil_iter, albedo_iter, lat, lon, expected',
    [
        # Rounded values based on high NDVI site (below)
        [0.125, 0, 4.2, 306.5, 0.875,
         298.5, 3.25, 10.0, 82, 306, 8600, 950, 3.25,
         25, 10, 38.7399, -121.5265, 5.7136059507579615],
        # # Same as above but with fewer iterations
        # [0.125, 0, 4.2, 306.5, 0.875,
        #  298.5, 3.25, 10.0, 82, 306, 8600, 950, 3.25,
        #  20, 10, 38.7399, -121.5265, 6.575051230124045],
        # # Same as above but with fewer iterations
        # [0.125, 0, 4.2, 306.5, 0.875,
        #  298.5, 3.25, 10.0, 82, 306, 8600, 950, 3.25,
        #  10, 10, 38.7399, -121.5265, 7.664097170794699],
        # # Same as above but with fewer iterations
        # [0.125, 0, 4.2, 306.5, 0.875,
        #  298.5, 3.25, 10.0, 82, 306, 8600, 950, 3.25,
        #  6, 3, 38.7399, -121.5265, 6.9778049553556825],
        # High NDVI site in LC08_044033_20170716
        [0.1259961, 0, 4.234, 306.5, 0.8743930074457752,
         298.4269785619036, 3.35975885, 3.0, 82,
         306.3861390180871, 8603.212890625, 946.69066527778, 3.2665367230039,
         25, 10, 38.7399, -121.5265, 5.723137358281673],
        # Low NDVI site in LC08_044033_20170716
        [0.17163020, 0, 0.5791, 323.6, 0.16195230171935662,
         300.36076433814867, 3.35975885, 4.0, 82,
         306.3840757488451, 8603.212890625, 946.69066527778, 3.2665367230039,
         25, 10, 38.71776, -121.50822, 1.3412559416959626],
    ]
)
def test_Image_et_values(albedo, cfmask, lai, lst, ndvi,
                         ta, alexi, elevation, landcover,
                         ta0, rs_daily, rs_hourly, windspeed,
                         stabil_iter, albedo_iter, lat, lon,
                         expected, tol=0.001):
    output_img = disalexi.Image(
        default_image(albedo=albedo, cfmask=cfmask, lai=lai, lst=lst, ndvi=ndvi),
        ta_source=ta, alexi_source=alexi,
        elevation_source=elevation, landcover_source=landcover,
        ta0_source=ta0, rs_daily_source=rs_daily, rs_hourly_source=rs_hourly,
        windspeed_source=windspeed, lat=lat, lon=lon,
        stability_iterations=stabil_iter, albedo_iterations=albedo_iter).et
    output = utils.point_image_value(output_img, [lon, lat])
    assert abs(output['et'] - expected) <= tol


def test_Image_et_properties():
    """Test if properties are set on the ET image"""
    output =  utils.getinfo(default_image_obj().et)
    assert output['bands'][0]['id'] == 'et'
    assert output['properties']['system:index'] == SCENE_ID
    assert output['properties']['system:time_start'] == SCENE_TIME
    assert output['properties']['image_id'] == COLL_ID + SCENE_ID
    # assert output['properties']['ta_iteration'] > 0


def test_Image_et_reference_properties():
    """Test if properties are set on the ETr image"""
    output =  utils.getinfo(default_image_obj().et_reference)
    assert output['bands'][0]['id'] == 'et_reference'
    assert output['properties']['system:index'] == SCENE_ID
    assert output['properties']['system:time_start'] == SCENE_TIME
    assert output['properties']['image_id'] == COLL_ID + SCENE_ID


def test_Image_et_reference_default_values(expected=10.0, tol=0.0001):
    output = utils.point_image_value(default_image_obj(
        et_reference_source=expected).et_reference, TEST_POINT)
    assert abs(output['et_reference'] - expected) <= tol


@pytest.mark.parametrize(
    'source, band, xy, expected',
    [
        ['IDAHO_EPSCOR/GRIDMET', 'etr', TEST_POINT, 11.2],
        ['projects/climate-engine/cimis/daily', 'ETr_ASCE', TEST_POINT, 10.124],
        ['10.8985', None, TEST_POINT, 10.8985],
        [10.8985, None, TEST_POINT, 10.8985],
    ]
)
def test_Image_et_reference_source(source, band, xy, expected, tol=0.001):
    output = utils.point_image_value(default_image_obj(
        et_reference_source=source, et_reference_band=band).et_reference, xy)
    assert abs(output['et_reference'] - expected) <= tol


# def test_Image_et_fraction_default_values(et_reference_source=10,
#                                           expected=2.5799/10, tol=0.0001):
#     # Check that ETf = ET / ETr
#     output = utils.point_image_value(default_image_obj(
#         et_reference_source=et_reference_source).et_fraction, TEST_POINT)
#     assert abs(output['et_fraction'] - expected) <= tol
#
#
# def test_Image_et_fraction_properties():
#     """Test if properties are set on the ETf image"""
#     output = utils.getinfo(default_image_obj().et_fraction)
#     assert output['bands'][0]['id'] == 'et_fraction'
#     assert output['properties']['system:index'] == SCENE_ID
#     assert output['properties']['system:time_start'] == SCENE_TIME


def test_Image_mask_values():
    # Mask is 1 for active pixels and nodata for cloudy pixels (cfmask >= 1)
    output = utils.point_image_value(
        default_image_obj(cfmask=1).mask, TEST_POINT)
    assert output['mask'] == None
    output = utils.point_image_value(
        default_image_obj(cfmask=0).mask, TEST_POINT)
    assert output['mask'] == 1


def test_Image_mask_properties():
    """Test if properties are set on the mask image"""
    output = utils.getinfo(default_image_obj().mask)
    assert output['bands'][0]['id'] == 'mask'
    assert output['properties']['system:index'] == SCENE_ID
    assert output['properties']['system:time_start'] == SCENE_TIME
    assert output['properties']['image_id'] == COLL_ID + SCENE_ID


def test_Image_time_point_values():
    output = utils.point_image_value(default_image_obj().time, TEST_POINT)
    assert output['time'] == utils.millis(SCENE_0UTC_DT)


def test_Image_time_properties():
    """Test if properties are set on the time image"""
    output = utils.getinfo(default_image_obj().time)
    assert output['bands'][0]['id'] == 'time'
    assert output['properties']['system:index'] == SCENE_ID
    assert output['properties']['system:time_start'] == SCENE_TIME
    assert output['properties']['image_id'] == COLL_ID + SCENE_ID


def test_Image_calculate_variables_default():
    output = utils.getinfo(default_image_obj().calculate())
    assert set([x['id'] for x in output['bands']]) == set(['et'])


def test_Image_calculate_variables_custom():
    variables = set(['ndvi'])
    output = utils.getinfo(default_image_obj().calculate(variables))
    assert set([x['id'] for x in output['bands']]) == variables


def test_Image_calculate_variables_all():
    variables = set(['et', 'mask', 'ndvi', 'time'])
    # variables = set(['et', 'etf', 'etr', 'mask', 'ndvi', 'time'])
    output = utils.getinfo(
        default_image_obj().calculate(variables=list(variables)))
    assert set([x['id'] for x in output['bands']]) == variables


def test_Image_calculate_properties():
    """Test if properties are set on the output image"""
    output =  utils.getinfo(default_image_obj().calculate(['ndvi']))
    assert output['properties']['system:index'] == SCENE_ID
    assert output['properties']['system:time_start'] == SCENE_TIME
    assert output['properties']['image_id'] == COLL_ID + SCENE_ID


def test_Image_calculate_values(tol=0.0001):
    """Test if the calculate method returns values"""
    output_img = disalexi.Image(
            default_image(albedo=0.125, cfmask=0, lai=4.7, lst=306, ndvi=0.875),
            ta_source=300, alexi_source=10, elevation_source=10,
            landcover_source=82, rs_daily_source=8600, rs_hourly_source=950,
            windspeed_source=3.25, stability_iterations=6, albedo_iterations=3)\
        .calculate(['et'])
    #     .calculate(['et', 'et_reference', 'et_fraction'])
    output = utils.point_image_value(output_img, TEST_POINT)
    assert abs(output['et'] - 7.080456256866455) <= tol
    # assert abs(output['et_reference'] - 10) <= tol
    # assert abs(output['et_fraction'] - 0.58) <= tol


def test_Image_calculate_variables_valueerror():
    """Test if calculate method raises a valueerror for invalid variables"""
    with pytest.raises(ValueError):
        utils.getinfo(default_image_obj().calculate(['FOO']))


# Landsat Collection 1 SR
def test_Image_from_landsat_c1_sr_default_image():
    """Test that the classmethod is returning a class object"""
    output = disalexi.Image.from_landsat_c1_sr(
        ee.Image(COLL_ID + SCENE_ID))
    assert type(output) == type(disalexi.Image(default_image()))


@pytest.mark.parametrize(
    'image_id',
    [
        'LANDSAT/LC08/C01/T1_SR/LC08_044033_20170716',
        'LANDSAT/LE07/C01/T1_SR/LE07_044033_20170708',
        'LANDSAT/LT05/C01/T1_SR/LT05_044033_20110716',
    ]
)
def test_Image_from_landsat_c1_sr_image_id(image_id):
    """Test instantiating the class from a Landsat image ID"""
    output = utils.getinfo(disalexi.Image.from_landsat_c1_sr(image_id).ndvi)
    assert output['properties']['system:index'] == image_id.split('/')[-1]


def test_Image_from_landsat_c1_sr_image():
    """Test instantiating the class from a Landsat ee.Image"""
    image_id = 'LANDSAT/LC08/C01/T1_SR/LC08_044033_20170716'
    output = utils.getinfo(disalexi.Image.from_landsat_c1_sr(
        ee.Image(image_id)).ndvi)
    assert output['properties']['system:index'] == image_id.split('/')[-1]


# def test_Image_from_landsat_c1_sr_et_fraction():
#     """Test if ETf can be built for a Landsat images"""
#     image_id = 'LANDSAT/LC08/C01/T1_SR/LC08_044033_20170716'
#     output = utils.getinfo(disalexi.Image.from_landsat_c1_sr(image_id).et_fraction)
#     assert output['properties']['system:index'] == image_id.split('/')[-1]


def test_Image_from_landsat_c1_sr_et():
    """Test if ET can be built for a Landsat images"""
    image_id = 'LANDSAT/LC08/C01/T1_SR/LC08_044033_20170716'
    # Using the collection 2 inputs since the c01 inputs no longer exist
    output = utils.getinfo(disalexi.Image.from_landsat_c1_sr(
        image_id, ta_source='CONUS_V003',
        lai_source='projects/earthengine-legacy/assets/projects/openet/lai/landsat/c02',
        tir_source='projects/earthengine-legacy/assets/projects/openet/tir/landsat/c02',
    ).et)
    assert output['properties']['system:index'] == image_id.split('/')[-1]


def test_Image_from_landsat_c1_sr_exception():
    """Test instantiating the class for an invalid image ID"""
    with pytest.raises(Exception):
        utils.getinfo(disalexi.Image.from_landsat_c1_sr(ee.Image('FOO')).ndvi)


# Landsat Collection 2 SR
def test_Image_from_landsat_c2_sr_default_image():
    """Test that the classmethod is returning a class object"""
    output = disalexi.Image.from_landsat_c2_sr(
        ee.Image(COLL_ID + SCENE_ID))
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
def test_Image_from_landsat_c2_sr_image_id(image_id):
    """Test instantiating the class from a Landsat image ID"""
    output = utils.getinfo(disalexi.Image.from_landsat_c2_sr(image_id).ndvi)
    assert output['properties']['system:index'] == image_id.split('/')[-1]


def test_Image_from_landsat_c2_sr_image():
    """Test instantiating the class from a Landsat ee.Image"""
    image_id = 'LANDSAT/LC08/C02/T1_L2/LC08_044033_20170716'
    output = utils.getinfo(disalexi.Image.from_landsat_c2_sr(
        ee.Image(image_id)).ndvi)
    assert output['properties']['system:index'] == image_id.split('/')[-1]


# def test_Image_from_landsat_c2_sr_et_fraction():
#     """Test if ETf can be built for a Landsat images"""
#     image_id = 'LANDSAT/LC08/C02/T1_L2/LC08_044033_20170716'
#     output = utils.getinfo(disalexi.Image.from_landsat_c1_sr(image_id).et_fraction)
#     assert output['properties']['system:index'] == image_id.split('/')[-1]


def test_Image_from_landsat_c2_sr_et():
    """Test if ET can be built for a Landsat images"""
    image_id = 'LANDSAT/LC08/C02/T1_L2/LC08_044033_20170716'
    output = utils.getinfo(disalexi.Image.from_landsat_c2_sr(image_id).et)
    assert output['properties']['system:index'] == image_id.split('/')[-1]


def test_Image_from_landsat_c2_sr_exception():
    """Test instantiating the class for an invalid image ID"""
    with pytest.raises(Exception):
        utils.getinfo(disalexi.Image.from_landsat_c2_sr(ee.Image('FOO')).ndvi)


@pytest.mark.parametrize(
    'image_id',
    [
        'LANDSAT/LC08/C01/T1_SR/LC08_044033_20170716',
        'LANDSAT/LC08/C02/T1_L2/LC08_044033_20170716',
    ]
)
def test_Image_from_image_id(image_id):
    """Test instantiating the class using the from_image_id method"""
    output = utils.getinfo(disalexi.Image.from_image_id(image_id).ndvi)
    assert output['properties']['system:index'] == image_id.split('/')[-1]
    # assert output['properties']['image_id'] == image_id


def test_Image_from_method_kwargs():
    """Test that the init parameters can be passed through the helper methods"""
    # assert disalexi.Image.from_landsat_c1_toa(
    #     'LANDSAT/LC08/C01/T1_TOA/LC08_042035_20150713',
    #     elevation_source='DEADBEEF').elevation_source == 'DEADBEEF'
    assert disalexi.Image.from_landsat_c1_sr(
        'LANDSAT/LC08/C01/T1_SR/LC08_042035_20150713',
        elevation_source='DEADBEEF').elevation_source == 'DEADBEEF'


# # CGM - Track the default input values to see if they change
# @pytest.mark.parametrize(
#     'xy, expected',
#     [
#         [TEST_POINT, 0.1259961],
#         [[-121.50822, 38.71776], 0.1716302],
#     ]
# )
# def test_Image_albedo_values(xy, expected, tol=0.0001):
#     output = utils.point_image_value(disalexi.Image.from_image_id(
#         'LANDSAT/LC08/C01/T1_SR/LC08_044033_20170716').albedo, xy)
#     assert abs(output['albedo'] - expected) <= tol
#
#
# @pytest.mark.parametrize(
#     'xy, expected',
#     [
#         [TEST_POINT, 0.8743930074457752],
#         [[-121.50822, 38.71776], 0.16195230171935662],
#     ]
# )
# def test_Image_ndvi_values(xy, expected, tol=0.0001):
#     output = utils.point_image_value(disalexi.Image.from_image_id(
#         'LANDSAT/LC08/C01/T1_SR/LC08_044033_20170716').ndvi, xy)
#     assert abs(output['ndvi'] - expected) <= tol
#
#
# @pytest.mark.parametrize(
#     'xy, expected',
#     [
#         [TEST_POINT, 4.234],
#         [[-121.50822, 38.71776], 0.5791],
#     ]
# )
# def test_Image_lai_values(xy, expected, tol=0.0001):
#     output = utils.point_image_value(disalexi.Image(default_image()).lai, xy)
#     assert abs(output['lai'] - expected) <= tol
#
#
# @pytest.mark.parametrize(
#     ' xy, expected',
#     [
#         [TEST_POINT, 306.5],
#         [[-121.50822, 38.71776], 323.6],
#     ]
# )
# def test_Image_tir_values(xy, expected, tol=0.0001):
#     output = utils.point_image_value(disalexi.Image(default_image()).tir, xy)
#     assert abs(output['tir'] - expected) <= tol



