import pprint

import ee
import pytest

import openet.disalexi as disalexi
import openet.disalexi.utils as utils
# TODO: import utils from openet.core
# import openet.core.utils as utils


COLLECTIONS = ['LANDSAT/LC08/C01/T1_SR']
SCENE_ID_LIST = sorted(['LC08_044033_20170716'])
# CGM - Need to build air temperature assets for these Landsat 7 images
# COLLECTIONS = ['LANDSAT/LC08/C01/T1_SR', 'LANDSAT/LE07/C01/T1_SR']
# SCENE_ID_LIST = sorted(['LC08_044033_20170716', 'LE07_044033_20170708',
#                         'LE07_044033_20170724'])
START_DATE = '2017-07-01'
END_DATE = '2017-08-01'
SCENE_POINT = (-121.9, 39)
VARIABLES = sorted(['et', 'etf', 'etr'])
TEST_POINT = (-121.5265, 38.7399)


def default_coll_args(collections=COLLECTIONS, start_date=START_DATE,
                      end_date=END_DATE, variables=VARIABLES,
                      cloud_cover_max=70, etr_source='IDAHO_EPSCOR/GRIDMET',
                      etr_band='etr', etr_factor=0.85):
    # Defining inside a function since this uses an ee.Geometry(),
    # but ee.Initialize() isn't called until after all tests are collected.
    return {
        'collections': collections,
        'start_date': start_date,
        'end_date': end_date,
        'geometry': ee.Geometry.Point(SCENE_POINT),
        'variables': variables,
        'cloud_cover_max': cloud_cover_max,
        'etr_source': etr_source,
        'etr_band': etr_band,
        'etr_factor': etr_factor,
    }


def default_coll_obj(collections=COLLECTIONS, start_date=START_DATE,
                     end_date=END_DATE, variables=VARIABLES,
                     cloud_cover_max=70, etr_source='IDAHO_EPSCOR/GRIDMET',
                     etr_band='etr', etr_factor=0.85):
    return disalexi.Collection(**default_coll_args(
        collections=collections, start_date=start_date, end_date=end_date,
        variables=variables, cloud_cover_max=cloud_cover_max,
        etr_source=etr_source, etr_band=etr_band, etr_factor=etr_factor,
    ))


# CGM - Should this be a fixture?
def parse_scene_id(output_info):
    output = [x['properties']['system:index'] for x in output_info['features']]
    # Strip merge indices (this works for Landsat image IDs
    return sorted(['_'.join(x.split('_')[-3:]) for x in output])


def test_Collection_init_default_parameters():
    """Test if init sets default parameters"""
    args = default_coll_args()

    # These values are being set above but have defaults that need to be checked
    del args['etr_source']
    del args['etr_band']
    del args['etr_factor']
    del args['variables']

    m = disalexi.Collection(**args)

    assert m.variables == None
    assert m.etr_source == None
    assert m.etr_band == None
    assert m.etr_factor == 1.0
    assert m.cloud_cover_max == 70
    assert m.model_args == {}
    assert m.filter_args == {}
    assert m._interp_vars == ['ndvi', 'etf']


def test_Collection_init_collection_str(coll_id='LANDSAT/LC08/C01/T1_SR'):
    """Test if a single coll_id str is converted to a single item list"""
    assert default_coll_obj(collections=coll_id).collections == [coll_id]


def test_Collection_init_cloud_cover_max_str():
    """Test if cloud_cover_max strings are converted to float"""
    assert default_coll_obj(cloud_cover_max='70').cloud_cover_max == 70


@pytest.mark.parametrize(
    'coll_id, start_date, end_date',
    [
        # ['LANDSAT/LT04/C01/T1_TOA', '1981-01-01', '1982-01-01'],
        # ['LANDSAT/LT04/C01/T1_TOA', '1994-01-01', '1995-01-01'],
        ['LANDSAT/LT05/C01/T1_SR', '1983-01-01', '1984-01-01'],
        ['LANDSAT/LT05/C01/T1_SR', '2012-01-01', '2013-01-01'],
        ['LANDSAT/LE07/C01/T1_SR', '1998-01-01', '1999-01-01'],
        ['LANDSAT/LC08/C01/T1_SR', '2012-01-01', '2013-01-01'],
    ]
)
def test_Collection_init_collection_filter(coll_id, start_date, end_date):
    """Test that collection IDs are filtered based on start/end dates"""
    # The target collection ID should be removed from the collections lists
    assert default_coll_obj(collections=coll_id, start_date=start_date,
                            end_date=end_date).collections == []


def test_Collection_init_startdate_exception():
    """Test if Exception is raised for invalid start date formats"""
    with pytest.raises(ValueError):
        default_coll_obj(start_date='1/1/2000', end_date='2000-01-02')


def test_Collection_init_enddate_exception():
    """Test if Exception is raised for invalid end date formats"""
    with pytest.raises(ValueError):
        default_coll_obj(start_date='2000-01-01', end_date='1/2/2000')


def test_Collection_init_swapped_date_exception():
    """Test if Exception is raised when start_date == end_date"""
    with pytest.raises(ValueError):
        default_coll_obj(start_date='2017-01-01', end_date='2017-01-01')


def test_Collection_init_invalid_collections_exception():
    """Test if Exception is raised for an invalid collection ID"""
    with pytest.raises(ValueError):
        default_coll_obj(collections='FOO')


# def test_Collection_init_duplicate_collections_exception():
#     """Test if Exception is raised for duplicate Landsat types"""
#     args = default_coll_args()
#     args['collections'] = ['LANDSAT/LC08/C01/T1_RT_TOA',
#                            'LANDSAT/LC08/C01/T1_TOA']
#     with pytest.raises(ValueError):
#         disalexi.Collection(**args)
#     args['collections'] = ['LANDSAT/LC08/C01/T1_SR', 'LANDSAT/LC08/C01/T1_TOA']
#     with pytest.raises(ValueError):
#         disalexi.Collection(**args)


def test_Collection_init_cloud_cover_exception():
    """Test if Exception is raised for an invalid cloud_cover_max"""
    with pytest.raises(TypeError):
        default_coll_obj(cloud_cover_max='A')
    with pytest.raises(ValueError):
        default_coll_obj(cloud_cover_max=-1)
    with pytest.raises(ValueError):
        default_coll_obj(cloud_cover_max=101)


# # TODO: Test for Error if geometry is not ee.Geometry
# def test_Collection_init_geometry_exception():
#     """Test if Exception is raised for an invalid geometry"""
#     args = default_coll_args()
#     args['geometry'] = 'DEADBEEF'
#     m = disalexi.Collection(**args)
#     assert utils.getinfo(m.geometry) ==


# TODO: Test if a geojson string can be passed for the geometry
# def test_Collection_init_geometry_geojson():
#     assert False


def test_Collection_build_default():
    output = utils.getinfo(default_coll_obj()._build())
    assert output['type'] == 'ImageCollection'
    assert parse_scene_id(output) == SCENE_ID_LIST
    # For the default build, check that the target variables are returned also
    assert VARIABLES == sorted(list(set([
        y['id'] for x in output['features'] for y in x['bands']])))


def test_Collection_build_variables():
    output = utils.getinfo(default_coll_obj()._build(variables=['ndvi']))
    assert ['ndvi'] == sorted(list(set([
        y['id'] for x in output['features'] for y in x['bands']])))


def test_Collection_build_dates():
    """Check that dates passed to build function override Class dates"""
    coll_obj = default_coll_obj(start_date='2017-08-01', end_date='2017-09-01')
    output = utils.getinfo(coll_obj._build(
        start_date='2017-07-16', end_date='2017-07-17'))
    assert parse_scene_id(output) == ['LC08_044033_20170716']


# def test_Collection_build_landsat_toa():
#     """Test if the Landsat TOA (non RT) collections can be built"""
#     coll_obj = default_coll_obj(
#         collections=['LANDSAT/LC08/C01/T1_TOA', 'LANDSAT/LE07/C01/T1_TOA'])
#     output = utils.getinfo(coll_obj._build())
#     assert parse_scene_id(output) == SCENE_ID_LIST
#     assert VARIABLES == sorted(list(set([
#         y['id'] for x in output['features'] for y in x['bands']])))


def test_Collection_build_landsat_sr():
    """Test if the Landsat SR collections can be built"""
    coll_obj = default_coll_obj(
        collections=['LANDSAT/LC08/C01/T1_SR', 'LANDSAT/LE07/C01/T1_SR'])
    output = utils.getinfo(coll_obj._build())
    assert parse_scene_id(output) == SCENE_ID_LIST
    assert VARIABLES == sorted(list(set([
        y['id'] for x in output['features'] for y in x['bands']])))


def test_Collection_build_exclusive_enddate():
    """Test if the end_date is exclusive"""
    output = utils.getinfo(default_coll_obj(end_date='2017-07-24')._build())
    assert [x for x in parse_scene_id(output) if int(x[-8:]) >= 20170724] == []


def test_Collection_build_cloud_cover():
    """Test if the cloud cover max parameter is being applied"""
    # CGM - The filtered images should probably be looked up programmatically
    output = utils.getinfo(default_coll_obj(cloud_cover_max=0.5)._build(
        variables=['et']))
    assert 'LE07_044033_20170724' not in parse_scene_id(output)


# DEADBEEF - I'm not sure what this test was supposed to do
#   It appears to be identical to the exclusive_enddate test above
# def test_Collection_build_model_args():
#     """Test if the end_date is exclusive"""
#     output = utils.getinfo(default_coll_obj(end_date='2017-07-24')._build(
#         variables=['et']))
#     assert [x for x in parse_scene_id(output) if int(x[-8:]) >= 20170724] == []


def test_Collection_build_filter_args():
    args = default_coll_args()
    coll_id = 'LANDSAT/LC08/C01/T1_SR'
    args['collections'] = [coll_id]
    args['geometry'] = ee.Geometry.Rectangle(-125, 35, -120, 40)
    args['filter_args'] = {coll_id: [
        {'type': 'equals', 'leftField': 'WRS_PATH', 'rightValue': 44},
        {'type': 'equals', 'leftField': 'WRS_ROW', 'rightValue': 33}]}
    output = utils.getinfo(disalexi.Collection(**args)._build(variables=['et']))
    assert set([x[5:11] for x in parse_scene_id(output)]) == set(['044033'])


def test_Collection_build_variable_exception():
    """Test if Exception is raised for an invalid variable"""
    with pytest.raises(ValueError):
        utils.getinfo(default_coll_obj()._build(variables=['FOO']))


def test_Collection_overpass_default():
    """Test overpass method with default values (variables from Class init)"""
    output = utils.getinfo(default_coll_obj().overpass())
    assert VARIABLES == sorted(list(set([
        y['id'] for x in output['features'] for y in x['bands']])))
    assert parse_scene_id(output) == SCENE_ID_LIST


def test_Collection_overpass_class_variables():
    """Test that custom class variables are passed through to build function"""
    output = utils.getinfo(default_coll_obj(variables=['et']).overpass())
    output = set([y['id'] for x in output['features'] for y in x['bands']])
    assert output == set(['et'])


def test_Collection_overpass_method_variables():
    """Test that custom method variables are passed through to build function"""
    output = utils.getinfo(default_coll_obj().overpass(variables=['et']))
    output = set([y['id'] for x in output['features'] for y in x['bands']])
    assert output == set(['et'])


def test_Collection_overpass_no_variables_exception():
    """Test if Exception is raised if variables is not set in init or method"""
    args = default_coll_args()
    del args['variables']
    with pytest.raises(ValueError):
        disalexi.Collection(**args).overpass().getInfo()


def test_Collection_interpolate_default():
    """Default t_interval should be custom"""
    output = utils.getinfo(default_coll_obj().interpolate())
    assert output['type'] == 'ImageCollection'
    assert parse_scene_id(output) == ['20170701']
    assert VARIABLES == sorted(list(set([
        y['id'] for x in output['features'] for y in x['bands']])))


def test_Collection_interpolate_variables_custom():
    output = utils.getinfo(default_coll_obj().interpolate(variables=['et']))
    assert [y['id'] for x in output['features'] for y in x['bands']] == ['et']


def test_Collection_interpolate_t_interval_daily():
    """Test if the daily time interval parameter works

    Since end_date is exclusive last image date will be one day earlier
    """
    coll_obj = default_coll_obj(start_date='2017-07-01', end_date='2017-07-05')
    output = utils.getinfo(coll_obj.interpolate(t_interval='daily'))
    assert output['type'] == 'ImageCollection'
    assert parse_scene_id(output)[0] == '20170701'
    assert parse_scene_id(output)[-1] == '20170704'
    assert VARIABLES == sorted(list(set([
        y['id'] for x in output['features'] for y in x['bands']])))


def test_Collection_interpolate_t_interval_monthly():
    """Test if the monthly time interval parameter works"""
    output = utils.getinfo(disalexi.Collection(**default_coll_args())
        .interpolate(t_interval='monthly'))
    assert output['type'] == 'ImageCollection'
    assert parse_scene_id(output) == ['201707']
    assert VARIABLES == sorted(list(set([
        y['id'] for x in output['features'] for y in x['bands']])))


# CGM - Commenting out since it takes a really long time to run
#   This function could probably be be tested for a shorter time period
# def test_Collection_interpolate_t_interval_annual():
#     """Test if the annual time interval parameter works"""
#     args = default_coll_args(start_date='2017-01-01', end_date='2018-01-01')
#     output = utils.getinfo(disalexi.Collection(**args)
#         .interpolate(t_interval='annual'))
#     assert output['type'] == 'ImageCollection'
#     assert parse_scene_id(output) == ['2017']
#     assert VARIABLES == sorted(list(set([
#         y['id'] for x in output['features'] for y in x['bands']])))


def test_Collection_interpolate_t_interval_custom():
    """Test if the custom time interval parameter works"""
    output = utils.getinfo(default_coll_obj().interpolate(t_interval='custom'))
    assert output['type'] == 'ImageCollection'
    assert parse_scene_id(output) == ['20170701']
    assert VARIABLES == sorted(list(set([
        y['id'] for x in output['features'] for y in x['bands']])))


# TODO: Write test for annual interpolation with a date range that is too short


# def test_Collection_interpolate_interp_days():
#     """Test if the interpolate interp_days parameter works"""
#     # Is there any way to test this without pulling values at a point?


# This is already being tested by test_Collection_interpolate_default() above
# def test_Collection_interpolate_etr_source_init():
#     """Test setting etr_source in the class init"""
#     args = default_coll_args()
#     args.update({'etr_source': 'IDAHO_EPSCOR/GRIDMET', 'etr_band': 'etr'})
#     output = utils.getinfo(disalexi.Collection(**args).interpolate())
#     assert VARIABLES == sorted(list(set([
#         y['id'] for x in output['features'] for y in x['bands']])))


def test_Collection_interpolate_etr_source_model_args():
    """Test setting etr_source in the model_args"""
    args = default_coll_args()
    del args['etr_source']
    del args['etr_band']
    args['model_args'] = {'etr_source': 'IDAHO_EPSCOR/GRIDMET', 'etr_band': 'etr'}
    output = utils.getinfo(disalexi.Collection(**args).interpolate())
    assert VARIABLES == sorted(list(set([
        y['id'] for x in output['features'] for y in x['bands']])))


def test_Collection_interpolate_etr_source_method():
    """Test setting etr_source in the interpolate call"""
    args = default_coll_args()
    del args['etr_source']
    del args['etr_band']
    etr_kwargs = {'etr_source': 'IDAHO_EPSCOR/GRIDMET', 'etr_band': 'etr'}
    output = utils.getinfo(disalexi.Collection(**args).interpolate(**etr_kwargs))
    assert VARIABLES == sorted(list(set([
        y['id'] for x in output['features'] for y in x['bands']])))


def test_Collection_interpolate_etr_source_not_set():
    """Test if Exception is raised if etr_source is not set"""
    args = default_coll_args()
    del args['etr_source']
    del args['etr_band']
    with pytest.raises(ValueError):
        utils.getinfo(disalexi.Collection(**args).interpolate())


# def test_Collection_interpolate_etr_source_exception():
#     """Test if Exception is raised if etr_source is invalid"""
#     args = default_coll_args()
#     args['model_args'] = {'etr_source': 'DEADBEEF', 'etr_band': 'etr'}
#     with pytest.raises(ValueError):
#         utils.getinfo(disalexi.Collection(**args).interpolate())


# def test_Collection_interpolate_etr_band_exception():
#     """Test if Exception is raised if etr_band is invalid"""
#     args = default_coll_args()
#     args['model_args'] = {'etr_source': 'IDAHO_EPSCOR/GRIDMET',
#                           'etr_band': 'DEADBEEF'}
#     with pytest.raises(ValueError):
#         utils.getinfo(disalexi.Collection(**args).interpolate())


def test_Collection_interpolate_t_interval_exception():
    """Test if Exception is raised for an invalid t_interval parameter"""
    with pytest.raises(ValueError):
        utils.getinfo(default_coll_obj().interpolate(t_interval='DEADBEEF'))


def test_Collection_interpolate_interp_method_exception():
    """Test if Exception is raised for an invalid interp_method parameter"""
    with pytest.raises(ValueError):
        utils.getinfo(default_coll_obj().interpolate(interp_method='DEADBEEF'))


def test_Collection_interpolate_interp_days_exception():
    """Test if Exception is raised for an invalid interp_days parameter"""
    with pytest.raises(ValueError):
        utils.getinfo(default_coll_obj().interpolate(interp_days=0))


def test_Collection_interpolate_no_variables_exception():
    """Test if Exception is raised if variables is not set in init or method"""
    args = default_coll_args()
    del args['variables']
    with pytest.raises(ValueError):
        utils.getinfo(disalexi.Collection(**args).interpolate())


def test_Collection_interpolate_output_type_default():
    """Test if output_type parameter is defaulting to float"""
    output = utils.getinfo(default_coll_obj(
        variables=['et', 'etr', 'etf', 'ndvi', 'count']).interpolate())
    output = output['features'][0]['bands']
    bands = {info['id']: i for i, info in enumerate(output)}
    assert(output[bands['et']]['data_type']['precision'] == 'float')
    assert(output[bands['etr']]['data_type']['precision'] == 'float')
    assert(output[bands['etf']]['data_type']['precision'] == 'float')
    assert(output[bands['ndvi']]['data_type']['precision'] == 'float')
    assert(output[bands['count']]['data_type']['precision'] == 'int')


@pytest.mark.parametrize(
    'output_type, precision, max_value',
    [
        ['int8', 'int', 127],
        ['uint8', 'int', 255],
        ['int16', 'int', 32767],
        ['uint16', 'int', 65535],
        ['float', 'float', None],
        ['double', 'double', None],
    ]
)
def test_Collection_interpolate_output_type_parameter(output_type, precision,
                                                      max_value):
    """Test if changing the output_type parameter works"""
    output = utils.getinfo(default_coll_obj().interpolate(output_type=output_type))
    output = output['features'][0]['bands']
    bands = {info['id']: i for i, info in enumerate(output)}

    assert(output[bands['et']]['data_type']['precision'] == precision)
    assert(output[bands['etr']]['data_type']['precision'] == precision)

    if max_value is not None:
        assert(output[bands['et']]['data_type']['max'] == max_value)
        assert(output[bands['etr']]['data_type']['max'] == max_value)
        