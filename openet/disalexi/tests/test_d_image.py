import datetime

import ee
import pytest

import openet.disalexi as disalexi
import openet.disalexi.utils as utils
# TODO: import utils from openet.core
# import openet.core.utils as utils


COLL_ID = 'LANDSAT/LC08/C01/T1_SR/'
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
def default_image(albedo=0.125, cfmask=0, lai=4.7, lst=306, ndvi=0.875):
    mask_img = ee.Image(COLL_ID + SCENE_ID).select(['B1']).multiply(0)
    # mask_img = ee.Image('projects/disalexi/ta/CONUS_V001/{}_0'.format(SCENE_ID))\
    #     .select(0).multiply(0)
    return ee.Image([mask_img.add(albedo), mask_img.add(cfmask),
                     mask_img.add(lai), mask_img.add(lst),
                     mask_img.add(ndvi)])\
        .rename(['albedo', 'cfmask', 'lai', 'lst', 'ndvi'])\
        .set({
            'system:index': SCENE_ID,
            'system:time_start': SCENE_TIME,
            # 'system:time_start': ee.Date(SCENE_DATE).millis(),
            'system:id': COLL_ID + SCENE_ID,
        })

def default_image_args(albedo=0.125, cfmask=0, lai=4.7, lst=306, ndvi=0.875,
                       etr_source=10, etr_band=None):
    return {
        'image': default_image(albedo=albedo, cfmask=cfmask, lai=lai, lst=lst,
                               ndvi=ndvi),
        'etr_source': etr_source, 'etr_band': etr_band
    }

def default_image_obj(albedo=0.125, cfmask=0, lai=4.7, lst=306, ndvi=0.875,
                      etr_source=10, etr_band=None):
    return disalexi.Image(**default_image_args(
        albedo=albedo, cfmask=cfmask, lai=lai, lst=lst, ndvi=ndvi,
        etr_source=etr_source, etr_band=etr_band,
    ))


def test_Image_init_default_parameters():
    m = disalexi.Image(default_image())
    assert m.ta_source == 'CONUS_V001'
    assert m.alexi_source == 'CONUS_V001'
    assert m.elevation_source == 'USGS/SRTMGL1_003'
    assert m.landcover_source == 'NLCD2011'
    assert m.rs_daily_source == 'MERRA2'
    assert m.rs_hourly_source == 'MERRA2'
    assert m.windspeed_source == 'CFSV2'
    assert m.stabil_iter == 36
    assert m.albedo_iter == 10
    assert m.rs_interp_flag == True
    # assert m.ta_interp_flag == True
    assert m.ta_smooth_flag == True
    assert m.etr_source == None
    assert m.etr_band == None
    assert m.etr_factor == 1.0


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
        ['CONUS_V001', TEST_POINT, 289.10],
        [ee.Image('USGS/SRTMGL1_003').multiply(0).add(10), TEST_POINT, 10],
        ['294.8', TEST_POINT, 294.8],  # Check constant values
        [294.8, TEST_POINT, 294.8],    # Check constant values
    ]
)
def test_Image_ta_sources(source, xy, expected, tol=0.01):
    m = disalexi.Image(default_image(), ta_source=source, ta_smooth_flag=False)
    output = utils.point_image_value(ee.Image(m.ta), xy)
    assert abs(output['ta'] - expected) <= tol


def test_Image_ta_smooth_flag(tol=0.01):
    m = disalexi.Image(default_image(), ta_source='CONUS_V001',
                       ta_smooth_flag=True)
    output = utils.point_image_value(ee.Image(m.ta), TEST_POINT)
    assert abs(output['ta'] - 292.4093) <= tol


def test_Image_ta_sources_exception():
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


@pytest.mark.parametrize(
    'source, xy, expected',
    [
        # ALEXI ET is currently in MJ m-2 d-1
        ['CONUS_V001', TEST_POINT, 10.382039 * 0.408],
        [ee.Image('USGS/SRTMGL1_003').multiply(0).add(10), TEST_POINT, 10],
        ['10.382039', TEST_POINT, 10.382039],
        [10.382039, TEST_POINT, 10.382039],
    ]
)
def test_Image_alexi_sources(source, xy, expected, tol=0.0001):
    output = utils.point_image_value(disalexi.Image(
        default_image(), alexi_source=source).et_alexi, xy)
    assert abs(output['et_alexi'] - expected) <= tol


def test_Image_alexi_sources_exception():
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
def test_Image_elevation_sources(source, xy, expected, tol=0.001):
    output = utils.point_image_value(disalexi.Image(
        default_image(), elevation_source=source).elevation, xy)
    assert abs(output['elevation'] - expected) <= tol


# def test_Image_elevation_sources_exception():
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
        ['NLCD2011', TEST_POINT, 82],
        ['NLCD2006', TEST_POINT, 82],
        ['GLOBELAND30', TEST_POINT, 10],
        [ee.Image('USGS/SRTMGL1_003').multiply(0).add(82), TEST_POINT, 82],
        ['82', TEST_POINT, 82],
        [82, TEST_POINT, 82],
    ]
)
def test_Image_landcover_sources(source, xy, expected, tol=0.001):
    output = utils.point_image_value(disalexi.Image(
        default_image(), landcover_source=source).lc_source, xy)
    assert abs(output['landcover'] - expected) <= tol


def test_Image_landcover_sources_exception():
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
        ['MERRA2', TEST_POINT, 8587.4091796875],
        [ee.Image('USGS/SRTMGL1_003').multiply(0).add(10), TEST_POINT, 10],
        ['8587.4091796875', TEST_POINT, 8587.4091796875],
        [8587.4091796875, TEST_POINT, 8587.4091796875],
    ]
)
def test_Image_rs_daily_sources(source, xy, expected, tol=0.0001):
    output = utils.point_image_value(disalexi.Image(
        default_image(), rs_daily_source=source).rs24, xy)
    assert abs(output['rs'] - expected) <= tol


def test_Image_rs_daily_sources_exception():
    with pytest.raises(ValueError):
        utils.getinfo(disalexi.Image(
            default_image(), rs_daily_source='').rs24)


def test_Image_rs_daily_band_name():
    output = utils.getinfo(disalexi.Image(default_image()).rs24)['bands'][0]['id']
    assert output == 'rs'


@pytest.mark.parametrize(
    'source, xy, expected',
    [
        ['MERRA2', TEST_POINT, 946.6906],
        [ee.Image('USGS/SRTMGL1_003').multiply(0).add(10), TEST_POINT, 10],
        ['946.6906', TEST_POINT, 946.6906],
        [946.6906, TEST_POINT, 946.6906],
    ]
)
def test_Image_rs_hourly_sources(source, xy, expected, tol=0.0001):
    output = utils.point_image_value(disalexi.Image(
        default_image(), rs_hourly_source=source).rs1, xy)
    assert abs(output['rs'] - expected) <= tol


def test_Image_rs_hourly_interp(tol=0.001):
    output = utils.point_image_value(disalexi.Image(
        default_image(), rs_hourly_source='MERRA2',
        rs_interp_flag=True).rs1, TEST_POINT)
    assert abs(output['rs'] - 946.6906) <= tol


def test_Image_rs_hourly_no_interp(tol=0.001):
    output = utils.point_image_value(disalexi.Image(
        default_image(), rs_hourly_source='MERRA2',
        rs_interp_flag=False).rs1, TEST_POINT)
    assert abs(output['rs'] - 929.75) <= tol


def test_Image_rs_hourly_sources_exception():
    with pytest.raises(ValueError):
        utils.getinfo(disalexi.Image(default_image(), rs_hourly_source='').rs1)


def test_Image_rs_hourly_band_name():
    output = utils.getinfo(disalexi.Image(default_image()).rs1)['bands'][0]['id']
    assert output == 'rs'


@pytest.mark.parametrize(
    'source, xy, expected',
    [
        ['CFSV2', TEST_POINT, 2.001476],
        [ee.Image('USGS/SRTMGL1_003').multiply(0).add(10), TEST_POINT, 10],
        ['2.001476', TEST_POINT, 2.001476],
        [2.001476, TEST_POINT, 2.001476],
    ]
)
def test_Image_windspeed_sources(source, xy, expected, tol=0.0001):
    output = utils.point_image_value(disalexi.Image(
        default_image(), windspeed_source=source).windspeed, xy)
    assert abs(output['windspeed'] - expected) <= tol


def test_Image_windspeed_sources_exception():
    with pytest.raises(ValueError):
        utils.getinfo(disalexi.Image(
            default_image(), windspeed_source='').windspeed)


def test_Image_windspeed_band_name():
    output = utils.getinfo(disalexi.Image(default_image()).windspeed)['bands'][0]['id']
    assert output == 'windspeed'


def test_Image_et_default_values(expected=2.5799, tol=0.0001):
    output = utils.point_image_value(default_image_obj().et, TEST_POINT)
    assert abs(output['et'] - expected) <= tol


@pytest.mark.parametrize(
    'ta, expected',
    [
        [289.10, 0.3314],  # Unsmoothed Ta at TEST_POINT
        [292.40, 2.5761],  # Smoothed Ta at TEST_POINT
    ]
)
def test_Image_et_fixed_sources(ta, expected, tol=0.0001):
    m = disalexi.Image(default_image(), ta_source=ta)
    output = utils.point_image_value(m.et, TEST_POINT)
    assert abs(output['et'] - expected) <= tol


@pytest.mark.parametrize(
    'albedo, cfmask, lai, lst, ndvi, ta, '
    'alexi, elevation, landcover, rs_daily, rs_hourly, windspeed, '
    'stabil_iter, albedo_iter, lat, lon, expected',
    [
        # Rounded values based on high NDVI site (below)
        [0.125, 0, 4.7, 306, 0.875, 300,
         3.25, 10.0, 82, 8600, 950, 3.25,
         36, 10, 38.7399, -121.5265, 6.94218],
        # Same as above but with fewer iterations
        [0.125, 0, 4.7, 306, 0.875, 300,
         3.25, 10.0, 82, 8600, 950, 3.25,
         6, 3, 38.7399, -121.5265, 6.94414],
        # High NDVI site in LC08_044033_20170716
        [0.1259961, 0, 4.6797005913579, 305.92253850611, 0.87439300744578, 300,
         3.35975885, 3.0, 82, 8603.212890625, 946.69066527778, 3.2665367230039,
         36, 10, 38.7399, -121.5265, 6.9511],
        # Low NDVI site in LC08_044033_20170716
        [0.17163020, 0, 0.029734416071998, 323.59893135545, 0.16195230171936, 300,
         3.35975885, 4.0, 82, 8603.212890625, 946.69066527778, 3.2665367230039,
         36, 10, 38.71776, -121.50822, 4.1705],
    ]
)
def test_Image_et_values(albedo, cfmask, lai, lst, ndvi, ta, alexi, elevation,
                         landcover, rs_daily, rs_hourly, windspeed, stabil_iter,
                         albedo_iter, lat, lon, expected, tol=0.0001):
    output_img = disalexi.Image(
        default_image(albedo=albedo, cfmask=cfmask, lai=lai, lst=lst, ndvi=ndvi),
        ta_source=ta, alexi_source=alexi, elevation_source=elevation,
        landcover_source=landcover, rs_daily_source=rs_daily,
        rs_hourly_source=rs_hourly, windspeed_source=windspeed, lat=lat, lon=lon,
        stabil_iterations=stabil_iter, albedo_iterations=albedo_iter).et
    output = utils.point_image_value(output_img, TEST_POINT)
    assert abs(output['et'] - expected) <= tol


def test_Image_et_properties():
    """Test if properties are set on the ET image"""
    output =  utils.getinfo(default_image_obj().et)
    assert output['bands'][0]['id'] == 'et'
    assert output['properties']['system:index'] == SCENE_ID
    assert output['properties']['system:time_start'] == SCENE_TIME
    assert output['properties']['image_id'] == COLL_ID + SCENE_ID
    # assert output['properties']['ta_iteration'] > 0


def test_Image_etr_properties():
    """Test if properties are set on the ETr image"""
    output =  utils.getinfo(default_image_obj().etr)
    assert output['bands'][0]['id'] == 'etr'
    assert output['properties']['system:index'] == SCENE_ID
    assert output['properties']['system:time_start'] == SCENE_TIME
    assert output['properties']['image_id'] == COLL_ID + SCENE_ID


def test_Image_etr_default_values(expected=10.0, tol=0.0001):
    output = utils.point_image_value(
        default_image_obj(etr_source=expected).etr, TEST_POINT)
    assert abs(output['etr'] - expected) <= tol


@pytest.mark.parametrize(
    'source, band, xy, expected',
    [
        ['IDAHO_EPSCOR/GRIDMET', 'etr', TEST_POINT, 10.8985],
        ['projects/climate-engine/cimis/daily', 'ETr_ASCE', TEST_POINT, 10.124],
        ['10.8985', None, TEST_POINT, 10.8985],
        [10.8985, None, TEST_POINT, 10.8985],
    ]
)
def test_Image_etr_sources(source, band, xy, expected, tol=0.001):
    output = utils.point_image_value(
        default_image_obj(etr_source=source, etr_band=band).etr, xy)
    assert abs(output['etr'] - expected) <= tol


def test_Image_etf_default_values(etr_source=10, expected=2.5799/10, tol=0.0001):
    # Check that ETf = ET / ETr
    output_img = default_image_obj(etr_source=etr_source).etf
    output = utils.point_image_value(output_img, TEST_POINT)
    assert abs(output['etf'] - expected) <= tol


def test_Image_etf_properties():
    """Test if properties are set on the ETf image"""
    output = utils.getinfo(default_image_obj().etf)
    assert output['bands'][0]['id'] == 'etf'
    assert output['properties']['system:index'] == SCENE_ID
    assert output['properties']['system:time_start'] == SCENE_TIME


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
    assert set([x['id'] for x in output['bands']]) == set(['et', 'etr', 'etf'])


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
            windspeed_source=3.25, stabil_iterations=6, albedo_iterations=3)\
        .calculate(['et'])
    #     .calculate(['et', 'etr', 'etf'])
    output = utils.point_image_value(output_img, TEST_POINT)
    assert abs(output['et'] - 6.94414) <= tol
    # assert abs(output['etr'] - 10) <= tol
    # assert abs(output['etf'] - 0.58) <= tol


def test_Image_calculate_variables_valueerror():
    """Test if calculate method raises a valueerror for invalid variables"""
    with pytest.raises(ValueError):
        utils.getinfo(default_image_obj().calculate(['FOO']))


# # How should these @classmethods be tested?
# def test_Image_from_landsat_c1_toa_default_image():
#     """Test that the classmethod is returning a class object"""
#     output = disalexi.Image.from_landsat_c1_toa(ee.Image(COLL_ID + SCENE_ID))
#     assert type(output) == type(default_image_obj())
#
#
# @pytest.mark.parametrize(
#     'image_id',
#     [
#         'LANDSAT/LC08/C01/T1_RT_TOA/LC08_044033_20170716',
#         'LANDSAT/LC08/C01/T1_TOA/LC08_044033_20170716',
#         'LANDSAT/LE07/C01/T1_RT_TOA/LE07_044033_20170708',
#         'LANDSAT/LE07/C01/T1_TOA/LE07_044033_20170708',
#         'LANDSAT/LT05/C01/T1_TOA/LT05_044033_20110716',
#     ]
# )
# def test_Image_from_landsat_c1_toa_image_id(image_id):
#     """Test instantiating the class from a Landsat image ID"""
#     output = utils.getinfo(disalexi.Image.from_landsat_c1_toa(image_id).ndvi)
#     assert output['properties']['system:index'] == image_id.split('/')[-1]
#
#
# def test_Image_from_landsat_c1_toa_image():
#     """Test instantiating the class from a Landsat ee.Image"""
#     image_id = 'LANDSAT/LC08/C01/T1_TOA/LC08_044033_20170716'
#     output = utils.getinfo(
#         disalexi.Image.from_landsat_c1_toa(ee.Image(image_id)).ndvi)
#     assert output['properties']['system:index'] == image_id.split('/')[-1]
#
#
# def test_Image_from_landsat_c1_toa_etf():
#     """Test if ETf can be built for a Landsat images"""
#     image_id = 'LANDSAT/LC08/C01/T1_TOA/LC08_044033_20170716'
#     output = utils.getinfo(disalexi.Image.from_landsat_c1_toa(image_id).etf)
#     assert output['properties']['system:index'] == image_id.split('/')[-1]
#
#
# def test_Image_from_landsat_c1_toa_et():
#     """Test if ET can be built for a Landsat images"""
#     image_id = 'LANDSAT/LC08/C01/T1_TOA/LC08_044033_20170716'
#     output = utils.getinfo(disalexi.Image.from_landsat_c1_toa(
#         image_id, etr_source='IDAHO_EPSCOR/GRIDMET', etr_band='etr').et)
#     assert output['properties']['system:index'] == image_id.split('/')[-1]
#
#
# def test_Image_from_landsat_c1_toa_exception():
#     with pytest.raises(Exception):
#         utils.getinfo(disalexi.Image.from_landsat_c1_toa(ee.Image('FOO')).ndvi)


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
    output = utils.getinfo(
        disalexi.Image.from_landsat_c1_sr(ee.Image(image_id)).ndvi)
    assert output['properties']['system:index'] == image_id.split('/')[-1]


# def test_Image_from_landsat_c1_sr_etf():
#     """Test if ETf can be built for a Landsat images"""
#     image_id = 'LANDSAT/LC08/C01/T1_SR/LC08_044033_20170716'
#     output = utils.getinfo(disalexi.Image.from_landsat_c1_sr(image_id).etf)
#     assert output['properties']['system:index'] == image_id.split('/')[-1]


def test_Image_from_landsat_c1_sr_et():
    """Test if ET can be built for a Landsat images"""
    image_id = 'LANDSAT/LC08/C01/T1_SR/LC08_044033_20170716'
    output = utils.getinfo(disalexi.Image.from_landsat_c1_sr(image_id).et)
    assert output['properties']['system:index'] == image_id.split('/')[-1]


def test_Image_from_landsat_c1_sr_exception():
    """Test instantiating the class for an invalid image ID"""
    with pytest.raises(Exception):
        utils.getinfo(disalexi.Image.from_landsat_c1_sr(ee.Image('FOO')).ndvi)


# @pytest.mark.parametrize(
#     'image_id',
#     [
#         'LANDSAT/LC08/C01/T1_TOA/LC08_044033_20170716',
#         'LANDSAT/LC08/C01/T1_SR/LC08_044033_20170716',
#     ]
# )
# def test_Image_from_image_id(image_id):
#     """Test instantiating the class using the from_image_id method"""
#     output = utils.getinfo(ssebop.Image.from_image_id(image_id).ndvi)
#     assert output['properties']['system:index'] == image_id.split('/')[-1]
#     assert output['properties']['image_id'] == image_id


def test_Image_from_method_kwargs():
    """Test that the init parameters can be passed through the helper methods"""
    # assert disalexi.Image.from_landsat_c1_toa(
    #     'LANDSAT/LC08/C01/T1_TOA/LC08_042035_20150713',
    #     elevation_source='DEADBEEF').elevation_source == 'DEADBEEF'
    assert disalexi.Image.from_landsat_c1_sr(
        'LANDSAT/LC08/C01/T1_SR/LC08_042035_20150713',
        elevation_source='DEADBEEF').elevation_source == 'DEADBEEF'
