import datetime

import ee
import pytest

import openet.disalexi as model
import openet.disalexi.utils as utils
# TODO: import utils from openet.core
# import openet.core.utils as utils


COLL_ID = 'LANDSAT/LC08/C01/T1_SR/'
SCENE_ID = 'LC08_044033_20170716'
SCENE_TIME = 1500230731090
# SCENE_ID = 'LC08_042035_20150713'
# SCENE_TIME = 1436812419150
SCENE_DT = datetime.datetime.utcfromtimestamp(SCENE_TIME / 1000.0)
print(SCENE_DT)
# SCENE_DT = datetime.datetime.strptime(SCENE_ID[-8:], '%Y%m%d')
SCENE_DATE = SCENE_DT.strftime('%Y-%m-%d')
SCENE_DOY = int(SCENE_DT.strftime('%j'))
SCENE_HOUR = (SCENE_DT.hour + SCENE_DT.minute / 60.0 +
              (SCENE_DT.second + SCENE_DT.microsecond / 1000000.0) / 3600.0)
# SCENE_HOUR = 18 + 33.0/60 + 39.15/3600
# SCENE_TIME = utils.millis(SCENE_DT)
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
def default_image_args(albedo=0.125, cfmask=0, lai=4.7, lst=306, ndvi=0.875):
    return {
        'image': default_image(albedo=albedo, cfmask=cfmask, lai=lai, lst=lst,
                               ndvi=ndvi),
    }
def default_image_obj(albedo=0.125, cfmask=0, lai=4.7, lst=306, ndvi=0.875):
    return model.Image(**default_image_args(
        albedo=albedo, cfmask=cfmask, lai=lai, lst=lst, ndvi=ndvi,
    ))


def test_ee_init():
    """Check that Earth Engine was initialized"""
    assert ee.Number(1).getInfo() == 1


def test_Image_init_default_parameters():
    m = model.Image(default_image())
    assert m.ta_source == 'CONUS_V001'
    assert m.alexi_source == 'CONUS_V001'
    assert m.elevation_source == 'USGS/SRTMGL1_003'
    assert m.landcover_source == 'NLCD2011'
    assert m.rs_daily_source == 'MERRA2'
    assert m.rs_hourly_source == 'MERRA2'
    assert m.windspeed_source == 'CFSV2'
    assert m.stabil_iter == 10
    assert m.albedo_iter == 3
    assert m.ta_interp_flag == True
    assert m.rs_interp_flag == True
    # assert m.etr_source == None
    # assert m.etr_band == None
    # assert m.etr_factor == 1.0


# Todo: Break these up into separate functions?
def test_Image_init_calculated_properties():
    d = default_image_obj()
    assert utils.getinfo(d.time_start) == SCENE_TIME
    # assert utils.getinfo(d.scene_id) == SCENE_ID
    # assert utils.getinfo(d.wrs2_tile) == 'p{}r{}'.format(
    #     SCENE_ID.split('_')[1][:3], SCENE_ID.split('_')[1][3:])


def test_Image_init_date_properties():
    d = default_image_obj()
    assert utils.getinfo(d.date)['value'] == SCENE_TIME
    # assert utils.getinfo(d.year) == int(SCENE_DATE.split('-')[0])
    # assert utils.getinfo(d.month) == int(SCENE_DATE.split('-')[1])
    # assert utils.getinfo(d.start_date)['value'] == SCENE_TIME
    # assert utils.getinfo(d.end_date)['value'] == utils.millis(
    #     SCENE_DT + datetime.timedelta(days=1))
    # assert utils.getinfo(d.doy) == SCENE_DOY
    # assert utils.getinfo(d.cycle_day) == int(
    #     (SCENE_DT - datetime.datetime(1970, 1, 3)).days % 8 + 1)
    assert utils.getinfo(d.hour) == SCENE_HOUR
    assert utils.getinfo(d.hour_int) == int(SCENE_HOUR)


# def test_Image_init_scene_id_property():
#     """Test that the system:index from a merged collection is parsed"""
#     input_img = default_image()
#     m = model.Image(input_img.set('system:index', '1_2_' + SCENE_ID))
#     assert utils.getinfo(m.scene_id) == SCENE_ID


# CGM - These aren't lazy properties and don't have properties being set on them
# def test_Image_ndvi_properties():
#     """Test if properties are set on the NDVI image"""
#     output = utils.getinfo(model.Image(default_image()).ndvi)
#     assert output['bands'][0]['id'] == 'ndvi'
#     assert output['properties']['system:index'] == SCENE_ID
#     assert output['properties']['system:time_start'] == SCENE_TIME
#     assert output['properties']['image_id'] == COLL_ID + SCENE_ID
#
#
# def test_Image_lst_properties():
#     """Test if properties are set on the LST image"""
#     output = utils.getinfo(model.Image(default_image()).lst)
#     assert output['bands'][0]['id'] == 'lst'
#     assert output['properties']['system:index'] == SCENE_ID
#     assert output['properties']['system:time_start'] == SCENE_TIME
#     assert output['properties']['image_id'] == COLL_ID + SCENE_ID


@pytest.mark.parametrize(
    'source, xy, expected',
    [
        ['CONUS_V001', TEST_POINT, 291.10235],
        [ee.Image('USGS/SRTMGL1_003').multiply(0).add(10), TEST_POINT, 10],
        ['294.8', TEST_POINT, 294.8],  # Check constant values
        [294.8, TEST_POINT, 294.8],    # Check constant values
    ]
)
def test_Image_ta_sources(source, xy, expected, tol=0.01):
    m = model.Image(default_image(), ta_source=source)
    output = utils.point_image_value(ee.Image(m.ta), xy)
    assert abs(output['ta'] - expected) <= tol


def test_Image_ta_sources_exception():
    with pytest.raises(ValueError):
        utils.getinfo(model.Image(default_image(), ta_source='').ta)


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
        ['CONUS_V001', TEST_POINT, 10.3344 * 0.408],
        # ['CONUS_V001', TEST_POINT, 10.3344],  # This is the MJ m-2 d-1 value
        [ee.Image('USGS/SRTMGL1_003').multiply(0).add(10), TEST_POINT, 10],
        ['10.3344', TEST_POINT, 10.3344],
        [10.3344, TEST_POINT, 10.3344],
    ]
)
def test_Image_alexi_sources(source, xy, expected, tol=0.001):
    output = utils.point_image_value(model.Image(
        default_image(), alexi_source=source).et_alexi, xy)
    assert abs(output['et_alexi'] - expected) <= tol


def test_Image_alexi_sources_exception():
    with pytest.raises(ValueError):
        utils.getinfo(model.Image(default_image(), alexi_source='').et_alexi)


def test_Image_alexi_band_name():
    output = utils.getinfo(
        model.Image(default_image()).et_alexi)['bands'][0]['id']
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
    output = utils.point_image_value(model.Image(
        default_image(), elevation_source=source).elevation, xy)
    assert abs(output['elevation'] - expected) <= tol


# def test_Image_elevation_sources_exception():
#     with pytest.raises(ValueError):
#         utils.getinfo(model.Image(
#             default_image(), elevation_source='').elevation)


def test_Image_elevation_band_name():
    output = utils.getinfo(
        model.Image(default_image()).elevation)['bands'][0]['id']
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
    output = utils.point_image_value(model.Image(
        default_image(), landcover_source=source).lc_source, xy)
    assert abs(output['landcover'] - expected) <= tol


def test_Image_landcover_sources_exception():
    with pytest.raises(ValueError):
        utils.getinfo(model.Image(
            default_image(), landcover_source='').landcover)


def test_Image_landcover_band_name():
    output = utils.getinfo(
        model.Image(default_image()).lc_source)['bands'][0]['id']
    assert output == 'landcover'


@pytest.mark.parametrize(
    'source, xy, expected',
    [
        ['MERRA2', TEST_POINT, 8603.2129],
        [ee.Image('USGS/SRTMGL1_003').multiply(0).add(10), TEST_POINT, 10],
        ['8603.212890625', TEST_POINT, 8603.2129],
        [8603.2129, TEST_POINT, 8603.2129],
    ]
)
def test_Image_rs_daily_sources(source, xy, expected, tol=0.0001):
    output = utils.point_image_value(model.Image(
        default_image(), rs_daily_source=source).rs24, xy)
    assert abs(output['rs'] - expected) <= tol


def test_Image_rs_daily_sources_exception():
    with pytest.raises(ValueError):
        utils.getinfo(model.Image(
            default_image(), rs_daily_source='').rs24)


def test_Image_rs_daily_band_name():
    output = utils.getinfo(model.Image(default_image()).rs24)['bands'][0]['id']
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
def test_Image_rs_hourly_sources(source, xy, expected, tol=0.001):
    output = utils.point_image_value(model.Image(
        default_image(), rs_hourly_source=source).rs1, xy)
    assert abs(output['rs'] - expected) <= tol


def test_Image_rs_hourly_interp(tol=0.001):
    output = utils.point_image_value(model.Image(
        default_image(), rs_hourly_source='MERRA2',
        rs_interp_flag=True).rs1, TEST_POINT)
    assert abs(output['rs'] - 946.6906) <= tol


def test_Image_rs_hourly_no_interp(tol=0.001):
    output = utils.point_image_value(model.Image(
        default_image(), rs_hourly_source='MERRA2',
        rs_interp_flag=False).rs1, TEST_POINT)
    assert abs(output['rs'] - 930.25) <= tol


def test_Image_rs_hourly_sources_exception():
    with pytest.raises(ValueError):
        utils.getinfo(model.Image(default_image(), rs_hourly_source='').rs1)


def test_Image_rs_hourly_band_name():
    output = utils.getinfo(model.Image(default_image()).rs1)['bands'][0]['id']
    assert output == 'rs'


@pytest.mark.parametrize(
    'source, xy, expected',
    [
        ['CFSV2', TEST_POINT, 3.2665],
        [ee.Image('USGS/SRTMGL1_003').multiply(0).add(10), TEST_POINT, 10],
        ['3.2665', TEST_POINT, 3.2665],
        [3.2665, TEST_POINT, 3.2665],
    ]
)
def test_Image_windspeed_sources(source, xy, expected, tol=0.001):
    output = utils.point_image_value(model.Image(
        default_image(), windspeed_source=source).windspeed, xy)
    assert abs(output['windspeed'] - expected) <= tol


def test_Image_windspeed_sources_exception():
    with pytest.raises(ValueError):
        utils.getinfo(model.Image(
            default_image(), windspeed_source='').windspeed)


def test_Image_windspeed_band_name():
    output = utils.getinfo(model.Image(default_image()).windspeed)['bands'][0]['id']
    assert output == 'windspeed'


def test_Image_et_values(tol=0.0001):
    output_img = model.Image(
        default_image(albedo=0.125, cfmask=0, lai=4.7, lst=306, ndvi=0.875),
        ta_source=295, alexi_source=10,elevation_source=10,
        landcover_source=82, rs_daily_source=8600, rs_hourly_source=950,
        windspeed_source=3.25, stabil_iterations=6, albedo_iterations=3).et
    output = utils.point_image_value(output_img, TEST_POINT)
    assert abs(output['et'] - 12.00397) <= tol


def test_Image_et_properties():
    """Test if properties are set on the ET image"""
    output =  utils.getinfo(default_image_obj().et)
    assert output['bands'][0]['id'] == 'et'
    assert output['properties']['system:index'] == SCENE_ID
    assert output['properties']['system:time_start'] == SCENE_TIME
    assert output['properties']['image_id'] == COLL_ID + SCENE_ID
    # assert output['properties']['ta_iteration'] > 0


# @pytest.mark.parametrize(
#     # Note: These are made up values
#     'lst, ndvi, elev, expected',
#     [
#         # Basic ETf test
#         [308, 0.50, 10, 50, 0.98, 310, 0.58],
#     ]
# )
# def test_Image_etf_values(lst, ndvi, elev, expected,
#                           tol=0.0001):
#     output_img = model.Image(
#             default_image(lst=lst, ndvi=ndvi), elev_source=elev).etf
#     output = utils.point_image_value(ee.Image(output_img))
#     assert abs(output['etf'] - expected) <= tol
#
#
# def test_Image_etf_properties():
#     """Test if properties are set on the ETf image"""
#     output = utils.getinfo(model.Image(default_image()).etf)
#     assert output['bands'][0]['id'] == 'etf'
#     assert output['properties']['system:index'] == SCENE_ID
#     assert output['properties']['system:time_start'] == SCENE_TIME
#
#
# def test_Image_etr_properties():
#     """Test if properties are set on the ETr image"""
#     output =  utils.getinfo(default_image_obj().etr)
#     assert output['bands'][0]['id'] == 'etr'
#     assert output['properties']['system:index'] == SCENE_ID
#     assert output['properties']['system:time_start'] == SCENE_TIME
#     assert output['properties']['image_id'] == COLL_ID + SCENE_ID
#
#
# def test_Image_etr_constant_values(etr=10.0, expected=8.5, tol=0.0001):
#     output = utils.point_image_value(default_image_obj(etr_source=etr).etr)
#     assert abs(output['etr'] - expected) <= tol


def test_Image_mask_values():
    output_img = model.Image(
        default_image(albedo=0.125, cfmask=0, lai=4.7, lst=306, ndvi=0.875)).mask
    output = utils.point_image_value(output_img, TEST_POINT)
    assert output['mask'] == 1


def test_Image_mask_band_name():
    output = utils.getinfo(default_image_obj().mask)
    assert output['bands'][0]['id'] == 'mask'


# CGM - Mask doesn't have properties since it isn't a lazy property
# def test_Image_mask_properties():
#     """Test if properties are set on the time image"""
#     output = utils.getinfo(default_image_obj().mask)['properties']
#     assert output['system:index'] == SCENE_ID
#     assert output['system:time_start'] == SCENE_TIME
#     assert output['image_id'] == COLL_ID + SCENE_ID


# def test_Image_time_values():
#     # The time image is currently being built from the etf image, so all the
#     #   ancillary values must be set for the constant_image_value to work.
#     output = utils.constant_image_value(model.Image(
#         default_image(ndvi=0.5, lst=308), dt_source=10, elev_source=50,
#         tcorr_source=0.98, tmax_source=310).time)
#     assert output['time'] == SCENE_TIME
#
#
# def test_Image_time_band_name():
#     output = utils.getinfo(default_image_obj().time)
#     assert output['bands'][0]['id'] == 'time'
#
#
# def test_Image_time_properties():
#     """Test if properties are set on the time image"""
#     output = utils.getinfo(default_image_obj().time)['properties']
#     assert output['system:index'] == SCENE_ID
#     assert output['system:time_start'] == SCENE_TIME
#     assert output['image_id'] == COLL_ID + SCENE_ID


def test_Image_calculate_variables_default():
    output = utils.getinfo(default_image_obj().calculate())
    assert set([x['id'] for x in output['bands']]) == set(['et'])
    # assert set([x['id'] for x in output['bands']]) == set(['et', 'etr', 'etf'])


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
    output_img = model.Image(
            default_image(albedo=0.125, cfmask=0, lai=4.7, lst=306, ndvi=0.875),
            ta_source=295, alexi_source=10, elevation_source=10,
            landcover_source=82, rs_daily_source=8600, rs_hourly_source=950,
            windspeed_source=3.25, stabil_iterations=6, albedo_iterations=3)\
        .calculate(['et'])
    #     .calculate(['et', 'etr', 'etf'])
    output = utils.point_image_value(output_img, TEST_POINT)
    assert abs(output['et'] - 12.00397) <= tol
    # assert abs(output['etr'] - 10) <= tol
    # assert abs(output['etf'] - 0.58) <= tol


def test_Image_calculate_variables_valueerror():
    """Test if calculate method raises a valueerror for invalid variables"""
    with pytest.raises(ValueError):
        utils.getinfo(default_image_obj().calculate(['FOO']))


# # How should these @classmethods be tested?
# def test_Image_from_landsat_c1_toa_default_image():
#     """Test that the classmethod is returning a class object"""
#     output = model.Image.from_landsat_c1_toa(ee.Image(COLL_ID + SCENE_ID))
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
#     output = utils.getinfo(model.Image.from_landsat_c1_toa(image_id).ndvi)
#     assert output['properties']['system:index'] == image_id.split('/')[-1]
#
#
# def test_Image_from_landsat_c1_toa_image():
#     """Test instantiating the class from a Landsat ee.Image"""
#     image_id = 'LANDSAT/LC08/C01/T1_TOA/LC08_044033_20170716'
#     output = utils.getinfo(
#         model.Image.from_landsat_c1_toa(ee.Image(image_id)).ndvi)
#     assert output['properties']['system:index'] == image_id.split('/')[-1]
#
#
# def test_Image_from_landsat_c1_toa_etf():
#     """Test if ETf can be built for a Landsat images"""
#     image_id = 'LANDSAT/LC08/C01/T1_TOA/LC08_044033_20170716'
#     output = utils.getinfo(model.Image.from_landsat_c1_toa(image_id).etf)
#     assert output['properties']['system:index'] == image_id.split('/')[-1]
#
#
# def test_Image_from_landsat_c1_toa_et():
#     """Test if ET can be built for a Landsat images"""
#     image_id = 'LANDSAT/LC08/C01/T1_TOA/LC08_044033_20170716'
#     output = utils.getinfo(model.Image.from_landsat_c1_toa(
#         image_id, etr_source='IDAHO_EPSCOR/GRIDMET', etr_band='etr').et)
#     assert output['properties']['system:index'] == image_id.split('/')[-1]
#
#
# def test_Image_from_landsat_c1_toa_exception():
#     with pytest.raises(Exception):
#         utils.getinfo(model.Image.from_landsat_c1_toa(ee.Image('FOO')).ndvi)


def test_Image_from_landsat_c1_sr_default_image():
    """Test that the classmethod is returning a class object"""
    output = model.Image.from_landsat_c1_sr(
        ee.Image(COLL_ID + SCENE_ID))
    assert type(output) == type(model.Image(default_image()))


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
    output = utils.getinfo(model.Image.from_landsat_c1_sr(image_id).ndvi)
    assert output['properties']['system:index'] == image_id.split('/')[-1]


def test_Image_from_landsat_c1_sr_image():
    """Test instantiating the class from a Landsat ee.Image"""
    image_id = 'LANDSAT/LC08/C01/T1_SR/LC08_044033_20170716'
    output = utils.getinfo(
        model.Image.from_landsat_c1_sr(ee.Image(image_id)).ndvi)
    assert output['properties']['system:index'] == image_id.split('/')[-1]


# def test_Image_from_landsat_c1_sr_etf():
#     """Test if ETf can be built for a Landsat images"""
#     image_id = 'LANDSAT/LC08/C01/T1_SR/LC08_044033_20170716'
#     output = utils.getinfo(model.Image.from_landsat_c1_sr(image_id).etf)
#     assert output['properties']['system:index'] == image_id.split('/')[-1]


def test_Image_from_landsat_c1_sr_et():
    """Test if ET can be built for a Landsat images"""
    image_id = 'LANDSAT/LC08/C01/T1_SR/LC08_044033_20170716'
    output = utils.getinfo(model.Image.from_landsat_c1_sr(image_id).et)
    assert output['properties']['system:index'] == image_id.split('/')[-1]


def test_Image_from_landsat_c1_sr_exception():
    """Test instantiating the class for an invalid image ID"""
    with pytest.raises(Exception):
        utils.getinfo(model.Image.from_landsat_c1_sr(ee.Image('FOO')).ndvi)


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
    # assert model.Image.from_landsat_c1_toa(
    #     'LANDSAT/LC08/C01/T1_TOA/LC08_042035_20150713',
    #     elevation_source='DEADBEEF').elevation_source == 'DEADBEEF'
    assert model.Image.from_landsat_c1_sr(
        'LANDSAT/LC08/C01/T1_SR/LC08_042035_20150713',
        elevation_source='DEADBEEF').elevation_source == 'DEADBEEF'
