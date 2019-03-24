#--------------------------------
# Name:         ta_export_wrs2_image.py
# Purpose:      Compute/Export WRS2 Ta images
#--------------------------------

import argparse
from builtins import input
import datetime
import json
import logging
import math
import os
import pprint
import random
import sys

import ee
# from osgeo import ogr, osr

import openet.disalexi as disalexi
# from . import utils
import utils


def main(ini_path=None, overwrite_flag=False, delay=0, key=None):
    """Compute WRS2 Ta images

    Parameters
    ----------
    ini_path : str
        Input file path.
    overwrite_flag : bool, optional
        If True, overwrite existing files (the default is False).
    delay : float, optional
        Delay time between each export task (the default is 0).
    key : str, optional
        File path to an Earth Engine json key file (the default is None).

    """
    logging.info('\nCompute WRS2 Ta images')

    ini = utils.read_ini(ini_path)

    # if (ini['DISALEXI']['alexi_source'] == 'projects/disalexi/alexi/CONUS_V001' and
    #         ini['INPUTS']['end_date'] < '2003-10-01'):
    #     logging.error(
    #         '\nCONUS ALEXI is not currently available before 2003-10-01, exiting\n')
    #     sys.exit()

    model_name = 'DISALEXI'
    # model_name = ini['INPUTS']['ET_MODEL'].upper()
    model_args = {
        k.lower(): float(v) if utils.is_number(v) else v
        for k, v in dict(ini[model_name]).items()}
    # This seems unnecessary
    tair_args = {}
    for k, v in dict(ini['TAIR']).items():
        if utils.is_number(v):
            if v.isdigit():
                tair_args[k.lower()] = int(v)
            else:
                tair_args[k.lower()] = float(v)
        else:
            tair_args[k.lower()] = v
    # tair_args = {
    #     k.lower(): int(v) if utils.is_number(v) else v
    #     for k, v in dict(ini['TAIR']).items()}

    try:
        wrs2_tiles = str(ini['INPUTS']['wrs2_tiles'])
        wrs2_tiles = sorted([x.strip() for x in wrs2_tiles.split(',')])
        wrs2_tiles = [x.replace('p', '').replace('r', '') for x in wrs2_tiles]
    except KeyError:
        logging.warning('\nwrs2_tiles parameter was net set, default to []')
        wrs2_tiles = []
        # raise ValueError('"wrs2_tiles" parameter was not set in INI')
    except Exception as e:
        raise e

    if 'cell_size' not in tair_args.keys():
        tair_args['cell_size'] = 30
    if 'retile' not in tair_args.keys():
        tair_args['retile'] = 0
    if ('source_coll' not in tair_args.keys() or
            tair_args['source_coll'].lower() == 'none'):
        tair_args['source_coll'] = None
    if tair_args['source_coll'] is None and 'ta_seed' in tair_args.keys():
        ta_seed = tair_args['ta_seed']
        logging.debug('  Seeding QM with Ta={}'.format(ta_seed))

    logging.info('\nTAIR Parameters')
    logging.info('  Source:     {}'.format(tair_args['source_coll']))
    logging.info('  Cell size:  {}'.format(tair_args['cell_size']))
    logging.info('  Retile:     {}'.format(tair_args['retile']))
    logging.info('  Step Size:  {}'.format(tair_args['step_size']))
    logging.info('  Step Count: {}'.format(tair_args['step_count']))

    # Output Ta daily image collection
    ta_wrs2_coll_id = '{}'.format(ini['EXPORT']['export_coll'])

    # ta_values = list(range(tair_args['ta_start'], tair_args['ta_stop'],
    #                        tair_args['ta_step']))
    # logging.info('  Start: {}'.format(tair_args['ta_start']))
    # logging.info('  Stop:  {}'.format(tair_args['ta_stop']))
    # logging.info('  Step:  {}'.format(tair_args['ta_step']))
    # if 'ta_iterations' not in tair_args.keys():
    #     logging.info('  Iterations: {}'.format(tair_args['ta_iterations']))

    logging.info('\nInitializing Earth Engine')
    if key:
        logging.info('  Using service account key file: {}'.format(key))
        # The "EE_ACCOUNT" parameter is not used if the key file is valid
        ee.Initialize(ee.ServiceAccountCredentials('deadbeef', key_file=key))
    else:
        ee.Initialize()

    # Get an ET image to set the Ta values to
    logging.debug('\nALEXI ET properties')
    alexi_coll_id = ini['DISALEXI']['alexi_et_source']
    if alexi_coll_id.upper() == 'CONUS_V001':
        alexi_coll_id = 'projects/disalexi/alexi/CONUS_V001'
        alexi_mask = ee.Image('projects/disalexi/alexi/conus_v001_mask')
    else:
        raise ValueError('unsupported ALEXI source')
    # alexi_coll = ee.ImageCollection(alexi_coll_id)
    # alexi_proj = alexi_mask.projection().getInfo()
    # alexi_geo = alexi_proj['transform']
    # alexi_crs = alexi_proj['crs']
    alexi_crs = 'EPSG:4326'
    alexi_geo = [0.04, 0.0, -125.04, 0.0, -0.04, 49.8]
    alexi_cs = 0.04
    alexi_x, alexi_y = -125.04, 49.8
    logging.debug('  Collection: {}'.format(alexi_coll_id))


    # if 'study_area_path' in ini['INPUTS'].keys():
    #     logging.info('\nReading study area shapefile')
    #     logging.info('  {}'.format(ini['INPUTS']['study_area_path']))
    #     study_area_ds = ogr.Open(ini['INPUTS']['study_area_path'], 0)
    #     study_area_lyr = study_area_ds.GetLayer()
    #     study_area_osr = study_area_lyr.GetSpatialRef()
    #     study_area_crs = str(study_area_osr.ExportToWkt())
    #     # study_area_proj4 = study_area_osr.ExportToProj4()
    #     logging.debug('  Study area projection: {}'.format(study_area_crs))
    #
    #     # Get the dissolved/unioned geometry of the study area
    #     output_geom = ogr.Geometry(ogr.wkbMultiPolygon)
    #     for study_area_ftr in study_area_lyr:
    #         study_area_geom = output_geom.Union(study_area_ftr.GetGeometryRef())
    #     study_area_ds = None
    #
    #     simplify_buffer = 0
    #     if simplify_buffer:
    #         study_area_geom = output_geom\
    #             .SimplifyPreserveTopology(simplify_buffer)\
    #             .buffer(simplify_buffer)
    #     else:
    #         # Added flatten call to change clockwise geometies to counter cw
    #         study_area_geom.FlattenTo2D()
    #
    #     export_geom = ee.Geometry(
    #         json.loads(study_area_geom.ExportToJson()), study_area_crs, False)
    #
    # else:
    export_geom = alexi_mask.geometry()


    # Get current asset list
    logging.debug('\nGetting asset list')
    # DEADBEEF - daily is hardcoded in the asset_id for now
    asset_list = utils.get_ee_assets(ta_wrs2_coll_id)

    # Get current running tasks
    tasks = utils.get_ee_tasks()
    if logging.getLogger().getEffectiveLevel() == logging.DEBUG:
        logging.debug('  Tasks: {}\n'.format(len(tasks)))
        input('ENTER')

    # Limit by year and month
    month_list = list(range(1, 13))
    year_list = []
    # try:
    #     month_list = sorted(list(utils.parse_int_set(ini['TCORR']['months'])))
    # except:
    #     logging.info('\nTCORR "months" parameter not set in the INI,'
    #                  '\n  Defaulting to all months (1-12)\n')
    #     month_list = list(range(1, 13))
    # try:
    #     year_list = sorted(list(utils.parse_int_set(ini['TCORR']['years'])))
    # except:
    #     logging.info('\nTCORR "years" parameter not set in the INI,'
    #                  '\n  Defaulting to all available years\n')
    #     year_list = []

    # Key is cycle day, value is a reference date on that cycle
    # Data from: https://landsat.usgs.gov/landsat_acq
    # I only need to use 8 cycle days because of 5/7 and 7/8 are offset
    cycle_dates = {
        7: '1970-01-01',
        8: '1970-01-02',
        1: '1970-01-03',
        2: '1970-01-04',
        3: '1970-01-05',
        4: '1970-01-06',
        5: '1970-01-07',
        6: '1970-01-08',
    }
    # cycle_dates = {
    #     1:  '2000-01-06',
    #     2:  '2000-01-07',
    #     3:  '2000-01-08',
    #     4:  '2000-01-09',
    #     5:  '2000-01-10',
    #     6:  '2000-01-11',
    #     7:  '2000-01-12',
    #     8:  '2000-01-13',
    #     # 9:  '2000-01-14',
    #     # 10: '2000-01-15',
    #     # 11: '2000-01-16',
    #     # 12: '2000-01-01',
    #     # 13: '2000-01-02',
    #     # 14: '2000-01-03',
    #     # 15: '2000-01-04',
    #     # 16: '2000-01-05',
    # }
    cycle_base_dt = datetime.datetime.strptime(cycle_dates[1], '%Y-%m-%d')

    iter_start_dt = datetime.datetime.strptime(
        ini['INPUTS']['start_date'], '%Y-%m-%d')
    iter_end_dt = datetime.datetime.strptime(
        ini['INPUTS']['end_date'], '%Y-%m-%d')

    # DEADBEEF - For testing, process the dates in random order
    date_list = list(utils.date_range(iter_start_dt, iter_end_dt))
    random.shuffle(date_list)

    # Iterate over date ranges
    for export_dt in date_list:
        export_date = export_dt.strftime('%Y-%m-%d')
        if ((month_list and export_dt.month not in month_list) or
                ( year_list and export_dt.year not in year_list)):
            logging.debug('Date: {} - skipping'.format(export_date))
            continue
        logging.info('Date: {}'.format(export_date))

        if export_date >= datetime.datetime.today().strftime('%Y-%m-%d'):
            logging.info('  Unsupported date, skipping')
            continue
        elif export_date < '1984-03-23':
            logging.info('  No Landsat 5+ images before 1984-03-16, skipping')
            continue
        elif export_date < '2013-01-01':
            logging.warning(
                '  Script is currently hardcoded for Landsat 7 and 8 only')
            sys.exit()

        # # Eventually build this using the model Collection class
        # Build and merge the Landsat collections
        # Time filters are to remove bad (L5) and pre-op (L8) images
        l8_coll = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR')\
            .filterDate(export_dt, export_dt + datetime.timedelta(days=1))\
            .filterBounds(export_geom)\
            .filter(ee.Filter.gt('system:time_start',
                                 ee.Date('2013-03-24').millis()))\
            .filterMetadata('CLOUD_COVER_LAND', 'less_than',
                            float(ini['INPUTS']['cloud_cover']))
        #    .filterMetadata('DATA_TYPE', 'equals', 'L1TP')\
        l7_coll = ee.ImageCollection('LANDSAT/LE07/C01/T1_SR')\
            .filterDate(export_dt, export_dt + datetime.timedelta(days=1))\
            .filterBounds(export_geom)\
            .filterMetadata('CLOUD_COVER_LAND', 'less_than',
                            float(ini['INPUTS']['cloud_cover']))
        #    .filterMetadata('DATA_TYPE', 'equals', 'L1TP')
        # l5_coll = ee.ImageCollection('LANDSAT/LT05/C01/T1_SR')\
        #     .filterDate(export_dt, export_dt + datetime.timedelta(days=1))\
        #     .filterBounds(export_geom)\
        #     .filterMetadata('CLOUD_COVER_LAND', 'less_than',
        #                     float(ini['INPUTS']['cloud_cover']))\
        #     .filter(ee.Filter.lt('system:time_start',
        #                          ee.Date('2011-12-31').millis()))
        # #     .filterMetadata('DATA_TYPE', 'equals', 'L1TP')
        # l4_coll = ee.ImageCollection('LANDSAT/LT04/C01/T1_SR')\
        #     .filterDate(export_dt, export_dt + datetime.timedelta(days=1))\
        #     .filterBounds(export_geom)\
        #     .filterMetadata('CLOUD_COVER_LAND', 'less_than',
        #                     float(ini['INPUTS']['cloud_cover']))
        # #     .filterMetadata('DATA_TYPE', 'equals', 'L1TP')

        # if export_date <= '1993-12-31':
        #     landsat_coll = ee.ImageCollection(l5_coll.merge(l4_coll))
        # if export_date < '1999-01-01':
        #     landsat_coll = l5_coll
        # elif export_date <= '2011-12-31':
        #     landsat_coll = ee.ImageCollection(l7_coll.merge(l5_coll))
        # elif export_date <= '2013-03-24':
        #     landsat_coll = l7_coll
        # else:
        #     landsat_coll = ee.ImageCollection(l8_coll.merge(l7_coll))
        landsat_coll = ee.ImageCollection(l8_coll.merge(l7_coll))
        # print(landsat_coll.getInfo())

        image_id_list = landsat_coll.aggregate_array('system:id').getInfo()
        # print(image_id_list)
        # input('ENTER')

        for image_id in image_id_list:
            scene_id = image_id.split('/')[-1]
            landsat, path, row, year, month, day = parse_landsat_id(scene_id)
            image_dt = datetime.datetime.strptime(
                '{:04d}{:02d}{:02d}'.format(year, month, day), '%Y%m%d')
            image_date = image_dt.date().isoformat()
            # logging.debug('  Date: {}'.format(image_date))
            # logging.debug('  DOY: {}'.format(doy))

            wrs2_tile = '{:03d}{:03d}'.format(path, row)
            if wrs2_tiles and wrs2_tile not in wrs2_tiles:
                logging.debug('{} - wrs2 not in list, skipping'.format(image_id))
                continue
            else:
                logging.info('{}'.format(image_id))

            export_id = ini['EXPORT']['export_id_fmt'] \
                .format(index=scene_id.lower())
            # export_id = '{}_qm_{}_{}_{}_{}_{}'.format(
            #     export_id,
            #     ini['DISALEXI']['stabil_iterations'],
            #     ini['DISALEXI']['albedo_iterations'],
            #     tair_args['step_size'],
            #     tair_args['retile'],
            #     tair_args['cell_size'],
            # )
            logging.debug('  Export ID: {}'.format(export_id))

            asset_id = '{}/{}'.format(ta_wrs2_coll_id, scene_id)
            # asset_id = '{}/{}_qm_{}_{}_{}_{}_{}'.format(
            #     ta_wrs2_coll_id, scene_id,
            #     ini['DISALEXI']['stabil_iterations'],
            #     ini['DISALEXI']['albedo_iterations'],
            #     tair_args['step_size'],
            #     tair_args['retile'],
            #     tair_args['cell_size'],
            # )
            logging.debug('  Asset ID: {}'.format(asset_id))

            if overwrite_flag:
                if export_id in tasks.keys():
                    logging.info('  Task already submitted, cancelling')
                    ee.data.cancelTask(tasks[export_id])
                # This is intentionally not an "elif" so that a task can be
                # cancelled and an existing image/file/asset can be removed
                if asset_id in asset_list:
                    logging.info('  Asset already exists, removing')
                    ee.data.deleteAsset(asset_id)
            else:
                if export_id in tasks.keys():
                    logging.info('  Task already submitted, exiting')
                    continue
                elif asset_id in asset_list:
                    logging.info('  Asset already exists, skipping')
                    continue

            properties = {
                'system:time_start': utils.millis(export_dt),
                'date_ingested': datetime.datetime.today().strftime('%Y-%m-%d'),
                'landsat': landsat,
                'wrs2_tile': wrs2_tile,
                'date': export_dt.strftime('%Y-%m-%d'),
                'year': int(export_dt.year),
                'month': int(export_dt.month),
                'day': int(export_dt.day),
                'doy': int(export_dt.strftime('%j')),
                'cycle_day': ((export_dt - cycle_base_dt).days % 8) + 1,
                'model_name': 'DISALEXI',
                'model_version': disalexi.__version__,
                'cloud_cover_max': float(ini['INPUTS']['cloud_cover']),
            }
            properties.update(model_args)
            properties.update(tair_args)

            # Either Ta is being read from an existing collection
            #   or it will be computed f
            if tair_args['source_coll'] is not None:
                ta_source_coll = ee.ImageCollection(tair_args['source_coll'])\
                    .filterMetadata('id', 'equals', image_id)
                if ta_source_coll.size().getInfo() == 0:
                    logging.info('  No images in Ta source coll, skipping')
                    continue
                ta_source_img = ee.Image(ta_source_coll.first()).select(['ta'])
            else:
                ta_source_img = alexi_mask.add(ta_seed)

            landsat_img = ee.Image(image_id)
            d_obj = disalexi.Image(
                disalexi.LandsatSR(landsat_img).prep(), **model_args)
            ta_img = d_obj.ta_qm(ta_img=ta_source_img,
                                 step_size=tair_args['step_size'],
                                 step_count=tair_args['step_count'],
                                 cell_size=tair_args['cell_size'])
            # ta = d_obj.ta_iter(cellsize=tair_args['ta_cellsize'],
            #                    ta_init=[tair_args['ta_start'],
            #                             tair_args['ta_stop']],
            #                    iterations=tair_args['ta_iterations'],
            #                    method='golden')
            export_img = ta_img.select(['ta', 'bias']).float()\
                .set({
                    'id': landsat_img.get('system:id'),
                    'system:index': landsat_img.get('system:index'),
                    'system:time_start': landsat_img.get('system:time_start'),
                    'spacecraft_id': landsat_img.get('SATELLITE'),
                })\
                .set(properties)

            if tair_args['retile'] and tair_args['retile'] > 0:
                export_img = export_img.retile(tair_args['retile'])

            # pprint.pprint(ee.ImageCollection(ta_source_img).getRegion(
            #     ee.Geometry.Point(-121.3022, 39.5785), scale=1000).getInfo())
            # pprint.pprint(ee.ImageCollection(ta_source_img).getRegion(
            #     ee.Geometry.Point(-121.2624, 39.5806), scale=1000).getInfo())
            # input('ENTER')
            #
            # pprint.pprint(ee.ImageCollection(d_obj.ta_single(ee.Image.constant(280))).getRegion(
            #     ee.Geometry.Point(-121.3022, 39.5785), scale=1000).getInfo())
            # pprint.pprint(ee.ImageCollection(d_obj.ta_single(ee.Image.constant(280))).getRegion(
            #     ee.Geometry.Point(-121.2624, 39.5806), scale=1000).getInfo())
            # input('ENTER')
            #
            # pprint.pprint(ee.ImageCollection(d_obj.ta_qm(ee.Image.constant(280))).getRegion(
            #     ee.Geometry.Point(-121.3022, 39.5785), scale=1000).getInfo())
            # pprint.pprint(ee.ImageCollection(d_obj.ta_qm(ee.Image.constant(280))).getRegion(
            #     ee.Geometry.Point(-121.2624, 39.5806), scale=1000).getInfo())
            # input('ENTER')

            # Build the export transform and shape from the Landsat image
            image_xy = landsat_img.geometry().bounds(1, 'EPSG:4326')\
                .coordinates().get(0).getInfo()
            export_extent = [min([xy[0] for xy in image_xy]),
                             min([xy[1] for xy in image_xy]),
                             max([xy[0] for xy in image_xy]),
                             max([xy[1] for xy in image_xy])]
            # Adjust extent to the cell size
            export_extent[0] = math.floor((
                export_extent[0] - alexi_x) / alexi_cs) * alexi_cs + alexi_x
            export_extent[1] = math.floor((
                export_extent[1] - alexi_y) / alexi_cs) * alexi_cs + alexi_y
            export_extent[2] = math.ceil((
                export_extent[2] - alexi_x) / alexi_cs) * alexi_cs + alexi_x
            export_extent[3] = math.ceil((
                export_extent[3] - alexi_y) / alexi_cs) * alexi_cs + alexi_y
            export_geo = [alexi_cs, 0, export_extent[0], 0,
                          -alexi_cs, export_extent[3]]
            export_shape = [
                int(abs(export_extent[2] - export_extent[0]) / alexi_cs),
                int(abs(export_extent[3] - export_extent[1]) / alexi_cs)]
            logging.debug('  CRS: {}'.format(alexi_crs))
            logging.debug('  Extent: {}'.format(export_extent))
            logging.debug('  Geo: {}'.format(export_geo))
            logging.debug('  Shape: {}'.format(export_shape))

            # Build export tasks
            logging.debug('  Building export task')
            task = ee.batch.Export.image.toAsset(
                image=export_img,
                description=export_id,
                assetId=asset_id,
                crs=alexi_crs,
                crsTransform='[' + ','.join(list(map(str, export_geo))) + ']',
                dimensions='{0}x{1}'.format(*export_shape),
            )
            logging.info('  Starting export task')
            utils.ee_task_start(task)

            # Pause before starting next task
            utils.delay_task(delay)
            logging.debug('')


def parse_landsat_id(system_index):
    """Return the components of an EE Landsat Collection 1 system:index

    Parameters
    ----------
    system_index : str

    Notes
    -----
    LT05_PPPRRR_YYYYMMDD

    """
    sensor = system_index[0:4]
    path = int(system_index[5:8])
    row = int(system_index[8:11])
    year = int(system_index[12:16])
    month = int(system_index[16:18])
    day = int(system_index[18:20])
    return sensor, path, row, year, month, day


def arg_parse():
    """"""
    parser = argparse.ArgumentParser(
        description='Compute/export WRS2 Ta images',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '-i', '--ini', type=utils.arg_valid_file,
        help='Input file', metavar='FILE')
    parser.add_argument(
        '--delay', default=0, type=float,
        help='Delay (in seconds) between each export tasks')
    parser.add_argument(
        '--key', type=utils.arg_valid_file, metavar='FILE',
        help='JSON key file')
    parser.add_argument(
        '-o', '--overwrite', default=False, action='store_true',
        help='Force overwrite of existing files')
    parser.add_argument(
        '-d', '--debug', default=logging.INFO, const=logging.DEBUG,
        help='Debug level logging', action='store_const', dest='loglevel')
    args = parser.parse_args()

    # Prompt user to select an INI file if not set at command line
    # if not args.ini:
    #     args.ini = utils.get_ini_path(os.getcwd())

    return args


if __name__ == "__main__":
    args = arg_parse()

    logging.basicConfig(level=args.loglevel, format='%(message)s')
    logging.info('\n{0}'.format('#' * 80))
    logging.info('{0:<20s} {1}'.format(
        'Run Time Stamp:', datetime.datetime.now().isoformat(' ')))
    logging.info('{0:<20s} {1}'.format('Current Directory:', os.getcwd()))
    logging.info('{0:<20s} {1}'.format(
        'Script:', os.path.basename(sys.argv[0])))

    main(ini_path=args.ini, overwrite_flag=args.overwrite, delay=args.delay,
         key=args.key)
