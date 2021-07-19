import argparse
from builtins import input
import configparser
import datetime
import json
import logging
import math
import os
import pprint
import re
import time

import ee

import openet.disalexi
import openet.core
import openet.core.utils as utils

TOOL_NAME = 'tair_image_wrs2_export'
TOOL_VERSION = '0.1.7'


def main(ini_path=None, overwrite_flag=False, delay_time=0, gee_key_file=None,
         ready_task_max=-1, reverse_flag=False, tiles=None, update_flag=False,
         log_tasks=True, recent_days=0, start_dt=None, end_dt=None):
    """Compute WRS2 Ta images

    Parameters
    ----------
    ini_path : str
        Input file path.
    overwrite_flag : bool, optional
        If True, overwrite existing files (the default is False).
    delay_time : float, optional
        Delay time in seconds between starting export tasks (or checking the
        number of queued tasks, see "max_ready" parameter).  The default is 0.
    gee_key_file : str, None, optional
        Earth Engine service account JSON key file (the default is None).
    ready_task_max: int, optional
        Maximum number of queued "READY" tasks.  The default is -1 which is
        implies no limit to the number of tasks that will be submitted.
    reverse_flag : bool, optional
        If True, process WRS2 tiles in reverse order (the default is False).
    tiles : str, None, optional
        List of MGRS tiles to process (the default is None).
    update_flag : bool, optional
        If True, only overwrite scenes with an older model version.
    log_tasks : bool, optional
        If True, log task information to the datastore (the default is True).
    recent_days : int, optional
        Limit start/end date range to this many days before the current date
        (the default is 0 which is equivalent to not setting the parameter and
         will use the INI start/end date directly).
    start_dt : datetime, optional
        Override the start date in the INI file
        (the default is None which will use the INI start date).
    end_dt : datetime, optional
        Override the (inclusive) end date in the INI file
        (the default is None which will use the INI end date).

    """
    logging.info('\nCompute WRS2 Ta images')

    # CGM - Which format should we use for the WRS2 tile?
    wrs2_tile_fmt = 'p{:03d}r{:03d}'
    # wrs2_tile_fmt = '{:03d}{:03d}'
    wrs2_tile_re = re.compile('p?(\\d{1,3})r?(\\d{1,3})')

    # List of path/rows to skip
    wrs2_skip_list = [
        'p049r026',  # Vancouver Island, Canada
        # 'p047r031', # North California coast
        'p042r037',  # San Nicholas Island, California
        # 'p041r037', # South California coast
        'p040r038', 'p039r038', 'p038r038',  # Mexico (by California)
        'p037r039', 'p036r039', 'p035r039',  # Mexico (by Arizona)
        'p034r039', 'p033r039',  # Mexico (by New Mexico)
        'p032r040',  # Mexico (West Texas)
        'p029r041', 'p028r042', 'p027r043', 'p026r043',  # Mexico (South Texas)
        'p019r040',  # West Florida coast
        'p016r043', 'p015r043',  # South Florida coast
        'p014r041', 'p014r042', 'p014r043',  # East Florida coast
        'p013r035', 'p013r036',  # North Carolina Outer Banks
        'p013r026', 'p012r026',  # Canada (by Maine)
        'p011r032',  # Rhode Island coast
    ]
    wrs2_path_skip_list = [9, 49]
    wrs2_row_skip_list = [25, 24, 43]

    #date_skip_list = [
    #    '2003-12-15', '2004-12-12', '2004-12-31', '2008-12-31',
    #    '2009-03-20', '2009-03-21', '2010-04-10', '2011-04-10',
    #    '2012-04-09', '2012-12-30', '2012-12-31',
    #    '2013-04-10', '2016-03-28', '2016-12-31',
    #    '2017-08-02', '2017-10-11', '2017-10-12', '2017-12-12',
    #    '2017-12-13', '2017-12-14', '2017-12-15', '2017-12-16',
    #    '2017-12-17', '2017-12-30', '2017-12-31',
    #    '2018-05-25', '2018-05-26', '2018-05-27', '2018-06-30', '2018-07-01',
    #    '2018-10-20', '2018-10-21', '2018-10-22', '2018-10-23', '2018-12-22',
    #    '2018-12-23', '2018-12-24', '2018-12-25', '2018-12-30', '2018-12-31',
    #    '2019-02-23', '2019-02-24', '2019-04-10', '2019-04-11', '2019-04-25',
    #    '2019-04-26', '2019-04-27', '2019-10-17', '2019-10-18',
    #    '2019-10-26', '2019-10-27',
    #]
    date_skip_list = []

    mgrs_skip_list = []

    export_id_fmt = '{model}_{index}'

    # TODO: Move to INI file
    # clip_ocean_flag = True

    # Read config file
    ini = read_ini(ini_path)
    # ini = configparser.ConfigParser(interpolation=None)
    # ini.read_file(open(ini_path, 'r'))
    # # Force conversion of unicode to strings
    # for section in ini.sections():
    #     ini[str(section)] = {}
    #     for k, v in ini[section].items():
    #         ini[str(section)][str(k)] = v

    # TODO: Move to INI parsing function or module
    # Required parameters
    try:
        model_name = str(ini['INPUTS']['et_model']).upper()
    except KeyError:
        raise ValueError('"et_model" parameter was not set in INI')
    except Exception as e:
        raise e
    logging.info('  ET Model: {}'.format(model_name))

    try:
        study_area_coll_id = str(ini['INPUTS']['study_area_coll'])
    except KeyError:
        raise ValueError('"study_area_coll" parameter was not set in INI')
    except Exception as e:
        raise e

    try:
        start_date = str(ini['INPUTS']['start_date'])
    except KeyError:
        raise ValueError('"start_date" parameter was not set in INI')
    except Exception as e:
        raise e

    try:
        end_date = str(ini['INPUTS']['end_date'])
    except KeyError:
        raise ValueError('"end_date" parameter was not set in INI')
    except Exception as e:
        raise e

    try:
        collections = str(ini['INPUTS']['collections'])
        collections = sorted([x.strip() for x in collections.split(',')])
    except KeyError:
        # CGM - I don't think we want to mix the collections here
        logging.info('\nINPUTS collections parameter was net set, '
                        'default to Landsat 5/7/8 C02 L2 collections')
        collections = ['LANDSAT/LC08/C02/T1_L2', 'LANDSAT/LE07/C02/T1_L2',
                       'LANDSAT/LT05/C02/T1_L2']
        # logging.info('\nINPUTS collections parameter was net set, '
        #                 'default to Landsat 5/7/8 C02 L2 and C01 SR collections')
        # collections = ['LANDSAT/LC08/C02/T1_L2', 'LANDSAT/LE07/C02/T1_L2',
        #                'LANDSAT/LT05/C02/T1_L2', 'LANDSAT/LC08/C01/T1_SR',
        #                'LANDSAT/LE07/C01/T1_SR', 'LANDSAT/LT05/C01/T1_SR']
    except Exception as e:
        raise e

    export_coll_id = '{}'.format(ini['EXPORT']['export_coll'])

    try:
        mgrs_ftr_coll_id = str(ini['EXPORT']['mgrs_ftr_coll'])
    except KeyError:
        raise ValueError('"mgrs_ftr_coll" parameter was not set in INI')
    except Exception as e:
        raise e

    # Optional parameters
    try:
        study_area_property = str(ini['INPUTS']['study_area_property'])
    except KeyError:
        study_area_property = None
        logging.debug('  study_area_property: not set in INI, defaulting to None')
    except Exception as e:
        raise e

    try:
        study_area_features = str(ini['INPUTS']['study_area_features'])
        study_area_features = sorted([
            x.strip() for x in study_area_features.split(',')])
    except KeyError:
        study_area_features = []
        logging.debug('  study_area_features: not set in INI, defaulting to []')
    except Exception as e:
        raise e

    try:
        wrs2_tiles = str(ini['INPUTS']['wrs2_tiles'])
        wrs2_tiles = sorted([x.strip() for x in wrs2_tiles.split(',')])
        # wrs2_tiles = [x.replace('p', '').replace('r', '') for x in wrs2_tiles]
    except KeyError:
        logging.debug('  wrs2_tiles: not set in INI, defaulting to []')
        wrs2_tiles = []
    except Exception as e:
        raise e

    try:
        mgrs_tiles = str(ini['EXPORT']['mgrs_tiles'])
        mgrs_tiles = sorted([x.strip() for x in mgrs_tiles.split(',')])
        # CGM - Remove empty strings caused by trailing or extra commas
        mgrs_tiles = [x.upper() for x in mgrs_tiles if x]
        logging.debug('  mgrs_tiles: {}'.format(mgrs_tiles))
    except KeyError:
        mgrs_tiles = []
        logging.debug('  mgrs_tiles: not set in INI, defaulting to []')
    except Exception as e:
        raise e

    try:
        utm_zones = str(ini['EXPORT']['utm_zones'])
        utm_zones = sorted([int(x.strip()) for x in utm_zones.split(',')])
        logging.debug('  utm_zones: {}'.format(utm_zones))
    except KeyError:
        utm_zones = []
        logging.debug('  utm_zones: not set in INI, defaulting to []')
    except Exception as e:
        raise e

    # try:
    #     output_type = str(ini['EXPORT']['output_type'])
    # except KeyError:
    #     output_type = 'float'
    #     # output_type = 'int16'
    #     logging.debug('  output_type: not set in INI, '
    #                   'defaulting to {}'.format(output_type))
    # except Exception as e:
    #     raise e
    #
    # try:
    #     scale_factor = int(ini['EXPORT']['scale_factor'])
    # except KeyError:
    #     scale_factor = 1
    #     # scale_factor = 10000
    #     logging.debug('  scale_factor: not set in INI, '
    #                   'defaulting to {}'.format(scale_factor))
    # except Exception as e:
    #     raise e

    try:
        export_id_name = '_' + str(ini['EXPORT']['export_id_name'])
    except KeyError:
        export_id_name = ''
        logging.debug('  export_id_name: not set in INI, defaulting to ""')
    except Exception as e:
        raise e


    # If the user set the tiles argument, use these instead of the INI values
    if tiles:
        logging.info('\nOverriding INI mgrs_tiles and utm_zones parameters')
        logging.info('  user tiles: {}'.format(tiles))
        mgrs_tiles = sorted([y.strip() for x in tiles for y in x.split(',')])
        mgrs_tiles = [x.upper() for x in mgrs_tiles if x]
        logging.info('  mgrs_tiles: {}'.format(', '.join(mgrs_tiles)))
        utm_zones = sorted(list(set([int(x[:2]) for x in mgrs_tiles])))
        logging.info('  utm_zones:  {}'.format(', '.join(map(str, utm_zones))))

    today_dt = datetime.datetime.now()
    today_dt = today_dt.replace(hour=0, minute=0, second=0, microsecond=0)
    if recent_days:
        logging.info('\nOverriding INI "start_date" and "end_date" parameters')
        logging.info('  Recent days: {}'.format(recent_days))
        end_dt = today_dt - datetime.timedelta(days=1)
        start_dt = today_dt - datetime.timedelta(days=recent_days)
        start_date = start_dt.strftime('%Y-%m-%d')
        end_date = end_dt.strftime('%Y-%m-%d')
    elif start_dt and end_dt:
        # Attempt to use the function start/end dates
        logging.info('\nOverriding INI "start_date" and "end_date" parameters')
        logging.info('  Custom date range')
        start_date = start_dt.strftime('%Y-%m-%d')
        end_date = end_dt.strftime('%Y-%m-%d')
    else:
        # Parse the INI start/end dates
        logging.info('\nINI date range')
        try:
            start_dt = datetime.datetime.strptime(start_date, '%Y-%m-%d')
            end_dt = datetime.datetime.strptime(end_date, '%Y-%m-%d')
        except Exception as e:
            raise e
    logging.info('  Start: {}'.format(start_date))
    logging.info('  End:   {}'.format(end_date))

    # TODO: Add a few more checks on the dates
    if end_dt < start_dt:
        raise ValueError('end date can not be before start date')

    logging.info('\nFilter date range')
    iter_start_dt = start_dt
    iter_end_dt = end_dt + datetime.timedelta(days=1)
    logging.info('  Start: {}'.format(iter_start_dt.strftime('%Y-%m-%d')))
    logging.info('  End:   {}'.format(iter_end_dt.strftime('%Y-%m-%d')))

    model_args = {
        k.lower(): float(v) if utils.is_number(v) else v
        for k, v in dict(ini['DISALEXI']).items()}

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

    if 'cell_size' not in tair_args.keys():
        tair_args['cell_size'] = 30
    if 'retile' not in tair_args.keys():
        tair_args['retile'] = 0
    if ('source_coll' not in tair_args.keys() or
            tair_args['source_coll'].lower() == 'none'):
        tair_args['source_coll'] = None
    # Clear Ta start value if is source is set
    if tair_args['source_coll'] is not None:
        tair_args['ta_start'] = None

    logging.info('\nDISALEXI Parameters')
    if 'stability_iterations' in model_args.keys():
        logging.info('  Stabil iter: {}'.format(int(model_args['stability_iterations'])))
    logging.info('  Albedo iter: {}'.format(int(model_args['albedo_iterations'])))

    logging.info('\nTAIR Parameters')
    logging.info('  Source:     {}'.format(tair_args['source_coll']))
    logging.info('  Ta Start:   {}'.format(tair_args['ta_start']))
    logging.info('  Cell size:  {}'.format(tair_args['cell_size']))
    logging.info('  Retile:     {}'.format(tair_args['retile']))
    logging.info('  Step Size:  {}'.format(tair_args['step_size']))
    logging.info('  Step Count: {}'.format(tair_args['step_count']))


    logging.info('\nInitializing Earth Engine')
    if gee_key_file:
        logging.info('  Using service account key file: {}'.format(gee_key_file))
        # The "EE_ACCOUNT" parameter is not used if the key file is valid
        ee.Initialize(ee.ServiceAccountCredentials(
            'deadbeef', key_file=gee_key_file))
    else:
        ee.Initialize()


    # # DEADBEEF
    # # TODO: set datastore key file as a parameter?
    # datastore_key_file = 'openet-dri-datastore.json'
    # if log_tasks and not os.path.isfile(datastore_key_file):
    #     logging.info('\nTask logging disabled, datastore key does not exist')
    #     log_tasks = False
    #     # input('ENTER')
    # if log_tasks:
    #     logging.info('\nInitializing task datastore client')
    #     # TODO: Move to top and add to requirements.txt and environment.yaml
    #     from google.cloud import datastore
    #     try:
    #         datastore_client = datastore.Client.from_service_account_json(
    #             datastore_key_file)
    #     except Exception as e:
    #         logging.error('{}'.format(e))
    #         return False


    # Get current running tasks
    if ready_task_max == -9999:
        # CGM - Getting the task list can take awhile so set ready tasks to
        #   -9999 to skip requesting it.  Only run this if you are going to
        #   manually avoid running existing tasks.
        # TODO: Check if this should disable delay_task() or set the
        #   ready_task_max to a large value
        tasks = {}
        ready_task_count = 0
    else:
        logging.info('\nRequesting Task List')
        tasks = utils.get_ee_tasks()
        if logging.getLogger().getEffectiveLevel() == logging.DEBUG:
            utils.print_ee_tasks(tasks)
            input('ENTER')
        ready_task_count = len(tasks.keys())
        logging.info(f'  Tasks: {ready_task_count}')
        # CGM - I'm still not sure if it makes sense to hold here or after the
        #   first task is started.
        ready_task_count = delay_task(
            delay_time=0, task_max=ready_task_max, task_count=ready_task_count)


    if not ee.data.getInfo(export_coll_id.rsplit('/', 1)[0]):
        logging.debug('\nFolder does not exist and will be built'
                      '\n  {}'.format(export_coll_id.rsplit('/', 1)[0]))
        input('Press ENTER to continue')
        ee.data.createAsset({'type': 'FOLDER'}, export_coll_id.rsplit('/', 1)[0])
    if not ee.data.getInfo(export_coll_id):
        logging.info('\nExport collection does not exist and will be built'
                     '\n  {}'.format(export_coll_id))
        input('Press ENTER to continue')
        ee.data.createAsset({'type': 'IMAGE_COLLECTION'}, export_coll_id)


    # Get an ET image to set the Ta values to
    logging.debug('\nALEXI ET properties')
    alexi_coll_id = model_args['alexi_source']
    if alexi_coll_id.upper() == 'CONUS_V002':
        alexi_coll_id = 'projects/earthengine-legacy/assets/' \
                        'projects/disalexi/alexi/CONUS_V002'
        alexi_mask = ee.Image('projects/earthengine-legacy/assets/'
                              'projects/disalexi/alexi/conus_v002_mask')\
            .double().multiply(0)
    elif alexi_coll_id.upper() == 'CONUS_V003':
        alexi_coll_id = 'projects/earthengine-legacy/assets/' \
                        'projects/disalexi/alexi/CONUS_V003'
        alexi_mask = ee.Image('projects/earthengine-legacy/assets/'
                              'projects/disalexi/alexi/conus_v002_mask')\
            .double().multiply(0)
    else:
        raise ValueError(f'unsupported ALEXI source: {alexi_coll_id}')
    # alexi_coll = ee.ImageCollection(alexi_coll_id)
    # alexi_proj = alexi_mask.projection().getInfo()
    # alexi_geo = alexi_proj['transform']
    # alexi_crs = alexi_proj['crs']
    alexi_crs = 'EPSG:4326'
    alexi_geo = [0.04, 0.0, -125.04, 0.0, -0.04, 49.8]
    alexi_cs = 0.04
    alexi_x, alexi_y = -125.04, 49.8
    logging.debug('  Collection: {}'.format(alexi_coll_id))


    # Get list of MGRS tiles that intersect the study area
    logging.debug('\nMGRS Tiles/Zones')
    export_list = mgrs_export_tiles(
        study_area_coll_id=study_area_coll_id,
        mgrs_coll_id=mgrs_ftr_coll_id,
        study_area_property=study_area_property,
        study_area_features=study_area_features,
        mgrs_tiles=mgrs_tiles,
        mgrs_skip_list=mgrs_skip_list,
        utm_zones=utm_zones,
        wrs2_tiles=wrs2_tiles,
    )
    if not export_list:
        logging.error('\nEmpty export list, exiting')
        return False
    # pprint.pprint(export_list)
    # input('ENTER')


    # Process each WRS2 tile separately
    logging.info('\nImage Exports')
    wrs2_tiles = []
    for export_info in sorted(export_list, key=lambda i: i['index'],
                              reverse=reverse_flag):
        logging.info('{}'.format(export_info['index']))
        # logging.info('  {} - {}'.format(
        #     export_info['index'], ', '.join(export_info['wrs2_tiles'])))
        tile_count = len(export_info['wrs2_tiles'])
        tile_list = sorted(export_info['wrs2_tiles'], reverse=not(reverse_flag))

        for export_n, wrs2_tile in enumerate(tile_list):
            path, row = map(int, wrs2_tile_re.findall(wrs2_tile)[0])
            if wrs2_tile in wrs2_tiles:
                logging.info('{} {} ({}/{}) - already processed'.format(
                    export_info['index'], wrs2_tile, export_n + 1, tile_count))
                continue
            elif wrs2_skip_list and wrs2_tile in wrs2_skip_list:
                logging.info('{} {} ({}/{}) - in wrs2 skip list'.format(
                    export_info['index'], wrs2_tile, export_n + 1, tile_count))
                continue
            elif wrs2_row_skip_list and row in wrs2_row_skip_list:
                logging.info('{} {} ({}/{}) - in wrs2 row skip list'.format(
                    export_info['index'], wrs2_tile, export_n + 1, tile_count))
                continue
            elif wrs2_path_skip_list and path in wrs2_path_skip_list:
                logging.info('{} {} ({}/{}) - in wrs2 path skip list'.format(
                    export_info['index'], wrs2_tile, export_n + 1, tile_count))
                continue
            else:
                logging.info('{} {} ({}/{})'.format(
                    export_info['index'], wrs2_tile, export_n + 1, tile_count))
            wrs2_tiles.append(wrs2_tile)

            # path, row = map(int, wrs2_tile_re.findall(export_info['index'])[0])
            # logging.info('WRS2 tile: {}  ({}/{})'.format(
            #     export_info['index'], export_n + 1, len(export_list)))
            #
            # logging.debug('  Shape:     {}'.format(export_info['shape']))
            # logging.debug('  Transform: {}'.format(export_info['geo_str']))
            # logging.debug('  Extent:    {}'.format(export_info['extent']))
            # logging.debug('  MaxPixels: {}'.format(export_info['maxpixels']))

            filter_args = {}
            for coll_id in collections:
                filter_args[coll_id] = [
                    {'type': 'equals', 'leftField': 'WRS_PATH', 'rightValue': path},
                    {'type': 'equals', 'leftField': 'WRS_ROW', 'rightValue': row}]
            logging.debug(f'  Filter Args: {filter_args}')

            # Get the full Landsat collection
            # Collection end date is exclusive
            model_obj = openet.disalexi.Collection(
                collections=collections,
                cloud_cover_max=float(ini['INPUTS']['cloud_cover']),
                start_date=iter_start_dt.strftime('%Y-%m-%d'),
                end_date=iter_end_dt.strftime('%Y-%m-%d'),
                geometry=ee.Geometry.Point(openet.core.wrs2.centroids[wrs2_tile]),
                model_args=model_args,
                filter_args=filter_args,
            )
            image_id_list = utils.get_info(ee.List(model_obj.overpass(
                variables=['ndvi']).aggregate_array('image_id')), max_retries=10)

            if not image_id_list:
                logging.info('  Empty image ID list, skipping tile')
                # logging.debug('  Empty image ID list, exiting')
                # return False

            # Get list of existing images for the target tile
            logging.debug('  Getting GEE asset list')
            asset_coll = ee.ImageCollection(export_coll_id) \
                .filterDate(iter_start_dt.strftime('%Y-%m-%d'),
                            iter_end_dt.strftime('%Y-%m-%d')) \
                .filterMetadata('wrs2_tile', 'equals',
                                wrs2_tile_fmt.format(path, row))
            asset_props = {f'{export_coll_id}/{x["properties"]["system:index"]}':
                               x['properties']
                           for x in utils.get_info(asset_coll)['features']}
            # asset_props = {x['id']: x['properties']
            #                for x in assets_info['features']}

            # # Get list of band types for checking to see if any bands are floats
            # asset_types = {
            #     f['id']: {b['id']: b['data_type']['precision'] for b in
            #               f['bands']}
            #     for f in assets_info['features']}

            # Sort image ID list by date
            image_id_list = sorted(
                image_id_list, key=lambda k: k.split('/')[-1].split('_')[-1],
                reverse=reverse_flag)
            # pprint.pprint(image_id_list)
            # input('ENTER')

            # Process each image in the collection by date
            # image_id is the full Earth Engine ID to the asset
            for image_id in image_id_list:
                logging.info('  {}'.format(image_id))
                coll_id, scene_id = image_id.rsplit('/', 1)
                l, p, r, year, month, day = parse_landsat_id(scene_id)
                image_dt = datetime.datetime.strptime(
                    '{:04d}{:02d}{:02d}'.format(year, month, day), '%Y%m%d')
                image_date = image_dt.strftime('%Y-%m-%d')
                next_date = (image_dt + datetime.timedelta(days=1)).strftime('%Y-%m-%d')
                logging.debug('    Date: {}'.format(image_date))
                # logging.debug('    DOY: {}'.format(doy))

                if date_skip_list and image_date in date_skip_list:
                    logging.info('    Date in skip list, skipping')
                    continue

                export_id = export_id_fmt.format(
                    model=ini['INPUTS']['et_model'].lower(),
                    index=image_id.lower().replace('/', '_'))
                export_id = export_id.replace('-', '')
                export_id += export_id_name
                logging.debug('    Export ID:  {}'.format(export_id))

                asset_id = '{}/{}'.format(export_coll_id, scene_id.lower())
                logging.debug('    Collection: {}'.format(
                    os.path.dirname(asset_id)))
                logging.debug('    Image ID:   {}'.format(
                    os.path.basename(asset_id)))

                if update_flag:
                    def version_number(version_str):
                        return list(map(int, version_str.split('.')))

                    if export_id in tasks.keys():
                        logging.info('    Task already submitted, skipping')
                        continue
                    # In update mode only overwrite if the version is old
                    if asset_props and asset_id in asset_props.keys():
                        model_ver = version_number(openet.disalexi.__version__)
                        asset_ver = version_number(
                            asset_props[asset_id]['model_version'])
                        # asset_flt = [
                        #     t == 'float' for b, t in asset_types.items()
                        #     if b in ['et', 'et_reference']]

                        if asset_ver < model_ver:
                            logging.info('    Existing asset model version is old, '
                                         'removing')
                            logging.debug(f'    asset: {asset_ver}\n'
                                          f'    model: {model_ver}')
                            try:
                                ee.data.deleteAsset(asset_id)
                            except:
                                logging.info('    Error removing asset, skipping')
                                continue
                        elif (asset_props[asset_id]['alexi_source'] <
                              model_args['alexi_source']):
                            logging.info('    ALEXI source is old, removing')
                            # input('ENTER')
                            try:
                                ee.data.deleteAsset(asset_id)
                            except:
                                logging.info('    Error removing asset, skipping')
                                continue
                        # elif (asset_props[asset_id]['date_ingested'] <= '2020-04-27'):
                        #     logging.info('    date_ingested is old, removing')
                        #     # input('ENTER')
                        #     try:
                        #         ee.data.deleteAsset(asset_id)
                        #     except:
                        #         logging.info('    Error removing asset, skipping')
                        #         continue
                        # elif ((('T1_RT_TOA' in asset_props[asset_id]['coll_id']) and
                        #        ('T1_RT_TOA' not in image_id)) or
                        #       (('T1_RT' in asset_props[asset_id]['coll_id']) and
                        #        ('T1_RT' not in image_id))):
                        #     logging.info(
                        #         '    Existing asset is from realtime Landsat '
                        #         'collection, removing')
                        #     try:
                        #         ee.data.deleteAsset(asset_id)
                        #     except:
                        #         logging.info('    Error removing asset, skipping')
                        #         continue
                        # elif (version_number(asset_props[asset_id]['tool_version']) <
                        #       version_number(TOOL_VERSION)):
                        #     logging.info('    Asset tool version is old, removing')
                        #     try:
                        #         ee.data.deleteAsset(asset_id)
                        #     except:
                        #         logging.info('    Error removing asset, skipping')
                        #         continue
                        # elif any(asset_flt):
                        #     logging.info(
                        #         '    Asset ET types are float, removing')
                        #     ee.data.deleteAsset(asset_id)
                        # elif 'tool_version' not in asset_props[asset_id].keys():
                        #     logging.info('    TOOL_VERSION property was not set, removing')
                        #     ee.data.deleteAsset(asset_id)

                        # elif asset_props[asset_id]['images'] == '':
                        #     logging.info('    Images property was not set, removing')
                        #     input('ENTER')
                        #     ee.data.deleteAsset(asset_id)
                        else:
                            logging.info('    Asset is up to date, skipping')
                            continue
                elif overwrite_flag:
                    if export_id in tasks.keys():
                        logging.info('    Task already submitted, cancelling')
                        ee.data.cancelTask(tasks[export_id]['id'])
                        # ee.data.cancelOperation(tasks[export_id]['id'])
                    # This is intentionally not an "elif" so that a task can be
                    # cancelled and an existing image/file/asset can be removed
                    if asset_props and asset_id in asset_props.keys():
                        logging.info('    Asset already exists, removing')
                        ee.data.deleteAsset(asset_id)
                else:
                    if export_id in tasks.keys():
                        logging.info('    Task already submitted, skipping')
                        continue
                    elif asset_props and asset_id in asset_props.keys():
                        logging.info('    Asset already exists, skipping')
                        continue

                if tair_args['source_coll'] is None:
                    logging.debug('    Tair source: {}'.format(tair_args['ta_start']))
                    ta_source_img = alexi_mask.add(float(tair_args['ta_start']))\
                        .rename(['ta'])
                elif tair_args['source_coll'] == 'NLDAS':
                    logging.debug('    Tair source: NLDAS')
                    ta_source_coll = ee.ImageCollection('NASA/NLDAS/FORA0125_H002')\
                        .filterDate(image_date, next_date)\
                        .select(['temperature'])
                    # Pulling maximum air temperature instead of 0 UTC
                    # input_image = ee.Image(ta_source_coll.first())\
                    input_image = ee.Image(ta_source_coll.reduce(ee.Reducer.max()))\
                        .add(273.15).floor()
                    ta_source_img = alexi_mask.add(input_image).rename(['ta'])
                else:
                    logging.debug('    Tair source: {}'.format(tair_args['source_coll']))
                    ta_source_coll = ee.ImageCollection(tair_args['source_coll'])\
                        .filterMetadata('image_id', 'equals', image_id)
                    if ta_source_coll.size().getInfo() == 0:
                        logging.info('    No Tair image in source coll, skipping')
                        # input('ENTER')
                        continue
                    ta_source_img = ta_min_bias(ee.Image(ta_source_coll.first()))

                # Manually check if the source LAI and TIR images are present
                # Eventually this should/could be done inside the model instead
                if ('tir_source' in model_args.keys() and
                        type(model_args['tir_source']) is str):
                    # Assumptions: string tir_source is an image collection ID
                    tir_coll = ee.ImageCollection(model_args['tir_source']) \
                        .filterMetadata('scene_id', 'equals', scene_id)
                    tir_img = ee.Image(tir_coll.first())
                    tir_info = tir_img.getInfo()
                    try:
                        sharpen_version = tir_info['properties']['sharpen_version']
                    except:
                        logging.info('    No TIR image in source, skipping')
                        input('ENTER')
                        continue

                if ('lai_source' in model_args.keys() and
                        type(model_args['lai_source']) is str):
                    # Assumptions: string lai_source is an image collection ID
                    lai_coll = ee.ImageCollection(model_args['lai_source']) \
                        .filterMetadata('scene_id', 'equals', scene_id)
                    lai_img = ee.Image(lai_coll.first())
                    lai_info = lai_img.getInfo()
                    try:
                        landsat_lai_version = lai_info['properties']['landsat_lai_version']
                    except:
                        logging.info('    No LAI image in source, skipping')
                        input('ENTER')
                        continue

                if ('alexi_source' in model_args.keys() and
                        type(model_args['alexi_source']) is str and
                        model_args['alexi_source'].upper() == 'CONUS_V003'):
                    alexi_coll_id = 'projects/disalexi/alexi/CONUS_V003'
                    alexi_coll = ee.ImageCollection(alexi_coll_id) \
                        .filterDate(image_date, next_date)
                    if alexi_coll.size().getInfo() == 0:
                        logging.info('    No ALEXI image in source, skipping')
                        input('ENTER')
                        continue

                # CGM: We could pre-compute (or compute once and then save)
                #   the crs, transform, and shape since they should (will?) be
                #   the same for each wrs2 tile
                landsat_img = ee.Image(image_id)
                output_info = utils.get_info(landsat_img.select([1]))
                # output_info = utils.get_info(landsat_img.select(['SR_B2']))

                d_obj = openet.disalexi.Image.from_image_id(image_id, **model_args)
                export_img = d_obj.ta_mosaic(
                    ta_img=ta_source_img,
                    step_size=tair_args['step_size'],
                    step_count=tair_args['step_count'])
                # pprint.pprint(export_img.getInfo())
                # input('ENTER')

                if tair_args['retile'] and tair_args['retile'] > 1:
                    export_img = export_img.retile(tair_args['retile'])

                properties = {
                    # Custom properties
                    'coll_id': coll_id,
                    'core_version': openet.core.__version__,
                    'date_ingested': datetime.datetime.today().strftime('%Y-%m-%d'),
                    'image_id': image_id,
                    'landsat_lai_version': landsat_lai_version,
                    'model_name': model_name,
                    'model_version': openet.disalexi.__version__,
                    'scene_id': scene_id,
                    'sharpen_version': sharpen_version,
                    'tool_name': TOOL_NAME,
                    'tool_version': TOOL_VERSION,
                    'wrs2_tile': wrs2_tile_fmt.format(p, r),
                    # Source properties
                    'CLOUD_COVER': output_info['properties']['CLOUD_COVER'],
                    'CLOUD_COVER_LAND': output_info['properties']['CLOUD_COVER_LAND'],
                    # CGM - Should we use the Landsat time or the ALEXI time?
                    'system:time_start': output_info['properties']['system:time_start'],
                    # 'system:time_start': utils.millis(image_dt),
                    # Other poperties
                    # 'spacecraft_id': landsat_img.get('SATELLITE'),
                    # 'landsat': landsat,
                    # 'date': image_dt.strftime('%Y-%m-%d'),
                    # 'year': int(image_dt.year),
                    # 'month': int(image_dt.month),
                    # 'day': int(image_dt.day),
                    # 'doy': int(image_dt.strftime('%j')),
                }
                properties.update(model_args)
                properties.update(tair_args)
                export_img = export_img.set(properties)

                # CGM: We could pre-compute (or compute once and then save)
                #   the crs, transform, and shape since they should (will?) be
                #   the same for each wrs2 tile
                # Build the export transform and shape from the Landsat image
                image_xy = landsat_img.geometry().bounds(1, 'EPSG:4326')\
                    .coordinates().get(0).getInfo()
                export_extent = [min([xy[0] for xy in image_xy]),
                                 min([xy[1] for xy in image_xy]),
                                 max([xy[0] for xy in image_xy]),
                                 max([xy[1] for xy in image_xy])]

                # Adjust extent to the cell size
                export_extent[0] = round(math.floor((
                    export_extent[0] - alexi_x) / alexi_cs) * alexi_cs + alexi_x, 8)
                export_extent[1] = round(math.floor((
                    export_extent[1] - alexi_y) / alexi_cs) * alexi_cs + alexi_y, 8)
                export_extent[2] = round(math.ceil((
                    export_extent[2] - alexi_x) / alexi_cs) * alexi_cs + alexi_x, 8)
                export_extent[3] = round(math.ceil((
                    export_extent[3] - alexi_y) / alexi_cs) * alexi_cs + alexi_y, 8)
                export_geo = [alexi_cs, 0, export_extent[0], 0,
                              -alexi_cs, export_extent[3]]
                export_shape = [
                    int(abs(export_extent[2] - export_extent[0]) / alexi_cs),
                    int(abs(export_extent[3] - export_extent[1]) / alexi_cs)]
                logging.debug('    CRS: {}'.format(alexi_crs))
                logging.debug('    Extent: {}'.format(export_extent))
                logging.debug('    Geo: {}'.format(export_geo))
                logging.debug('    Shape: {}'.format(export_shape))

                # Build export tasks
                max_retries = 4
                logging.debug('    Building export task')
                # TODO: Move into try/except if getting EEException
                # for i in range(1, max_retries):
                #     try:
                task = ee.batch.Export.image.toAsset(
                    image=export_img,
                    description=export_id,
                    assetId=asset_id,
                    crs=alexi_crs,
                    crsTransform='[' + ','.join(list(map(str, export_geo))) + ']',
                    dimensions='{0}x{1}'.format(*export_shape),
                )
                #     # except ee.ee_exception.EEException as e:
                #     except Exception as e:
                #         if ('Earth Engine memory capacity exceeded' in str(e) or
                #                 'Earth Engine capacity exceeded' in str(e)):
                #             logging.info('    Rebuilding task ({}/{})'.format(
                #                 i, max_retries))
                #             logging.debug('    {}'.format(e))
                #             time.sleep(i ** 2)
                #         else:
                #             logging.warning('Unhandled exception\n{}'.format(e))
                #             break
                #             raise e

                if not task:
                    logging.warning('    Export task was not built, skipping')
                    continue

                logging.info('    Starting export task')
                for i in range(1, max_retries):
                    try:
                        task.start()
                        break
                    except Exception as e:
                        logging.info('    Resending query ({}/{})'.format(
                            i, max_retries))
                        logging.debug('    {}'.format(e))
                        time.sleep(i ** 2)
                # # Not using ee_task_start since it doesn't return the task object
                # utils.ee_task_start(task)

                # DEADBEEF
                # # Write the export task info the openet-dri project datastore
                # if log_tasks:
                #     logging.debug('    Writing datastore entity')
                #     try:
                #         task_obj = datastore.Entity(key=datastore_client.key(
                #             'Task', task.status()['id']),
                #             exclude_from_indexes=['properties'])
                #         for k, v in task.status().items():
                #             task_obj[k] = v
                #         # task_obj['date'] = datetime.datetime.today() \
                #         #     .strftime('%Y-%m-%d')
                #         task_obj['index'] = properties.pop('wrs2_tile')
                #         # task_obj['wrs2_tile'] = properties.pop('wrs2_tile')
                #         task_obj['model_name'] = properties.pop('model_name')
                #         # task_obj['model_version'] = properties.pop('model_version')
                #         task_obj['runtime'] = 0
                #         task_obj['start_timestamp_ms'] = 0
                #         task_obj['tool_name'] = properties.pop('tool_name')
                #         task_obj['properties'] = json.dumps(properties)
                #         datastore_client.put(task_obj)
                #     except Exception as e:
                #         # CGM - The message/handling will probably need to be updated
                #         #   We may want separate try/excepts on the create and the put
                #         logging.warning('\nDatastore entity was not written')
                #         logging.warning('{}\n'.format(e))

                # Pause before starting the next export task
                ready_task_count += 1
                ready_task_count = delay_task(
                    delay_time=delay_time, task_max=ready_task_max,
                    task_count=ready_task_count)

                logging.debug('')


    # # DEADBEEF - Old code for iteratively computing a new Tair image
    # if iteration <= 0:
    #      # For initial iteration compute bias at 250 and 350 K
    #      a_img = d_obj.ta_coarse(ta_img=ee.Image.constant(250)) \
    #          .select(['ta', 'bias'], ['ta_a', 'bias_a'])
    #      b_img = d_obj.ta_coarse(ta_img=ee.Image.constant(350)) \
    #          .select(['ta', 'bias'], ['ta_b', 'bias_b'])
    #      # Applying both bias masks to the output
    #      # This shouldn't be necessary but ta_coarse was returning
    #      #   ta and bias images with different masks
    #      export_img = ee.Image([a_img, b_img])\
    #          .updateMask(a_img.select(['bias_a']).And(
    #              b_img.select(['bias_b'])))\
    #          .double()
    #  # elif iteration <= 4:
    #  #     # Interpolate new Ta from the bias and test directly
    #  #     # Roughly equivalent to false position method
    #  #     ta_img = ee.Image('{}/{}_{:02d}'.format(
    #  #         ta_wrs2_coll_id, scene_id, iteration - 1))
    #  #     # ta_img = ee.Image(ta_coll
    #  #     #     .filterMetadata('date', 'equals', export_date)
    #  #     #     .filterMetadata('iteration', 'equals', iteration - 1)
    #  #     #     .first())
    #  #     ta_a = ta_img.select(['ta_a'])
    #  #     ta_b = ta_img.select(['ta_b'])
    #  #     bias_a = ta_img.select(['bias_a'])
    #  #     bias_b = ta_img.select(['bias_b'])
    #  #
    #  #     ta_x = bias_a.multiply(ta_b).subtract(bias_b.multiply(ta_a))\
    #  #         .divide(bias_a.subtract(bias_b))
    #  #     bias_x = d_obj.et_bias(d_obj.et_coarse(ta_x))
    #  #
    #  #     # Use the new value if it minimizes the bias and brackets 0
    #  #     mask1 = bias_x.lt(bias_b).And(bias_x.gt(0))
    #  #     mask2 = bias_x.gt(bias_a).And(bias_x.lt(0))
    #  #     ta_b = ta_b.where(mask1, ta_x)
    #  #     bias_b = bias_b.where(mask1, bias_x)
    #  #     ta_a = ta_a.where(mask2, ta_x)
    #  #     bias_a = bias_a.where(mask2, bias_x)
    #  #     export_img = ee.Image([ta_a, bias_a, ta_b, bias_b])
    #  else:
    #      # Generate test Ta randomly from a triangular distribution
    #      # Center the distribution on the interpolated zero bias Ta
    #      ta_img = ee.Image('{}/{}_{:02d}'.format(
    #          ta_wrs2_coll_id, scene_id, iteration - 1))
    #      # ta_img = ee.Image(ta_coll
    #      #     .filterMetadata('date', 'equals', export_date)
    #      #     .filterMetadata('iteration', 'equals', iteration - 1)
    #      #     .first())
    #      ta_a = ta_img.select(['ta_a'])
    #      ta_b = ta_img.select(['ta_b'])
    #      bias_a = ta_img.select(['bias_a'])
    #      bias_b = ta_img.select(['bias_b'])
    #
    #      ta_c = bias_a.multiply(ta_b).subtract(bias_b.multiply(ta_a))\
    #          .divide(bias_a.subtract(bias_b))
    #      # ta_c = ta_a.add(ta_b).multiply(0.5)
    #
    #      # For now use a single random number for the whole scene
    #      # Need to check if ee.Image.random() will work though
    #      u = ta_b.multiply(0).add(random.random())
    #      # u = ta_b.multiply(0).add(ee.Image.random(0))
    #
    #      a = u.multiply(ta_b.subtract(ta_a))\
    #          .multiply(ta_c.subtract(ta_a)).sqrt().add(ta_a)
    #      b = u.multiply(-1).add(1)\
    #          .multiply(ta_b.subtract(ta_a))\
    #          .multiply(ta_b.subtract(ta_c))\
    #          .sqrt().multiply(-1).add(ta_b)
    #      fc = ta_c.subtract(ta_a).divide(ta_b.subtract(ta_a))
    #      ta_x = a.where(u.gt(fc), b)
    #      bias_x = d_obj.et_bias(d_obj.et_coarse(ta_x))
    #
    #      # Use the new value if it minimizes the bias and brackets 0
    #      mask1 = bias_x.lt(bias_b).And(bias_x.gt(0))
    #      mask2 = bias_x.gt(bias_a).And(bias_x.lt(0))
    #      ta_b = ta_b.where(mask1, ta_x)
    #      bias_b = bias_b.where(mask1, bias_x)
    #      ta_a = ta_a.where(mask2, ta_x)
    #      bias_a = bias_a.where(mask2, bias_x)
    #      export_img = ee.Image([ta_a, bias_a, ta_b, bias_b])
    #
    #  # else:
    #  #     ta_img = ee.Image('{}/{}_{}'.format(
    #  #         ta_wrs2_coll_id, scene_id, iteration - 1))
    #  #     # ta_img = ee.Image(ta_coll
    #  #     #     .filterMetadata('date', 'equals', export_date)
    #  #     #     .filterMetadata('iteration', 'equals', iteration - 1)
    #  #     #     .first())
    #  #     ta_a = ta_img.select(['ta_a'])
    #  #     ta_b = ta_img.select(['ta_b'])
    #  #     bias_a = ta_img.select(['bias_a'])
    #  #     bias_b = ta_img.select(['bias_b'])
    #  #     # abs_a = ta_img.select(['bias_a']).abs()
    #  #     # abs_b = ta_img.select(['bias_b']).abs()
    #  #
    #  #     # Compute new test Ta and biases
    #  #     ta_c = ta_b.subtract(ta_b.subtract(ta_a).multiply(0.618034))
    #  #     ta_d = ta_a.add(ta_b.subtract(ta_a).multiply(0.618034))
    #  #     bias_c = d_obj.et_bias(d_obj.et_coarse(ta_c))
    #  #     bias_d = d_obj.et_bias(d_obj.et_coarse(ta_d))
    #  #     abs_c = bias_c.abs()
    #  #     abs_d = bias_d.abs()
    #  #
    #  #     # Use the new values if they minimize the bias
    #  #     # If f(c) < f(d): move the data from d to b and c to d
    #  #     mask1 = abs_c.lt(abs_d)
    #  #     # If f(c) > f(d): move the data from c to a and d to c
    #  #     mask2 = abs_c.gte(abs_d)
    #  #     ta_b = ta_b.where(mask1, ta_d)
    #  #     bias_b = bias_b.where(mask1, bias_d)
    #  #     ta_a = ta_a.where(mask2, ta_c)
    #  #     bias_a = bias_a.where(mask2, bias_c)
    #  #
    #  #     export_img = ee.Image([ta_a, bias_a, ta_b, bias_b])


# TODO: Move this function into the model code so it can be tested
def ta_min_bias(input_img):
    """

    Parameters
    ----------
    input_img

    Returns
    -------

    """
    input_img = ee.Image(input_img)

    # Reverse the band order so that we can find the last transition
    #   from decreasing to increasing with a positive bias
    ta_bands = input_img.select('step_\\d+_ta').bandNames().reverse()
    bias_bands = input_img.select('step_\\d+_bias').bandNames().reverse()
    ta_array = input_img.select(ta_bands).toArray()
    bias_array = input_img.select(bias_bands).toArray()
    #Assign the bias that are very similar a very large value that they will not be selected
    diff = bias_array.arraySlice(0,1).subtract(bias_array.arraySlice(0,0,-1))
    bias_array_mask = diff.abs().lt(0.001)
    #repeat the last value to make the array the same length. array is reversed order.
    bias_array_mask = bias_array_mask.arrayCat(bias_array_mask.arraySlice(0,-1),0)
    bias_array = bias_array.add(bias_array_mask.multiply(99))

    # Identify the "last" transition from a negative to positive bias
    # CGM - Having problems with .ceil() limiting result to the image data range
    #   Multiplying by a big number seemed to fix the issue but this could still
    #     be a problem with the bias ranges get really small
    sign_array = bias_array.multiply(1000).ceil().max(0).min(1).int()
    transition_array = sign_array.arraySlice(0, 0, -1)\
        .subtract(sign_array.arraySlice(0, 1))
    # Insert an extra value at the beginning (of reverse, so actually at end)
    #   of the transition array so the indexing lines up for all steps
    transition_array = bias_array.arraySlice(0, 0, 1).multiply(0)\
        .arrayCat(transition_array, 0)
    transition_index = transition_array.arrayArgmax().arrayFlatten([['index']])
    # Get the max transition value in order to know if there was a transition
    transition_max = transition_array\
        .arrayReduce(ee.Reducer.max(), [0]).arrayFlatten([['max']])

    # Identify the position of minimum absolute bias
    min_bias_index = bias_array.abs().multiply(-1).arrayArgmax()\
        .arrayFlatten([['index']])

    # Identify the "bracketing" Ta and bias values
    # If there is a transition, use the "last" transition
    # If there is not a transition, use the minimum absolute bias for both
    # Note, the index is for the reversed arrays
    index_b = transition_index.subtract(1).max(0)\
        .where(transition_max.eq(0), min_bias_index)
    index_a = transition_index.min(ta_bands.size().subtract(1))\
        .where(transition_max.eq(0), min_bias_index)
    ta_b = ta_array.arrayGet(index_b)
    ta_a = ta_array.arrayGet(index_a)
    bias_b = bias_array.arrayGet(index_b)
    bias_a = bias_array.arrayGet(index_a)

    # For now, compute the target Ta as the average of the bracketing Ta values
    # Eventually Ta could be linearly interpolated or computed
    #   as some sort of weighted average (based on the biases)
    ta_source_img = ta_a.add(ta_b).multiply(0.5).rename(['ta'])

    # Mask out Ta cells with all negative biases
    ta_source_img = ta_source_img\
        .updateMask(bias_b.lt(0).And(bias_a.lt(0)).Not())

    return ta_source_img


# import pytest
# import pprint
# import ee
# ee.Initialize()
# @pytest.mark.parametrize(
#     'ta_list, bias_list, expected',
#     [
#         # Normal bias profile, select average of bracketing Ta values
#         [[257, 267, 277, 287, 297, 307, 317, 327, 337, 347, 357],
#          [-0.2, -0.1, 0.1, 0.6, 2.4, 4.6, 5.7, 6.2, 6.5, 6.8, 7.0],
#          272],
#         # Normal bias profile, crossing at top interval
#         [[257, 267, 277, 287, 297, 307, 317, 327, 337, 347, 357],
#          [-1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, 0.1],
#          352],
#         # Normal bias profile, crossing at bottom interval
#         [[257, 267, 277, 287, 297, 307, 317, 327, 337, 347, 357],
#          [-0.2, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
#          262],
#         # Increasing then decreasing then increasing biases
#         # Last transition should be selected
#         [[257, 267, 277, 287, 297, 307, 317, 327, 337, 347, 357],
#          [-0.2, 0.1, -0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9],
#          282],
#         # Increasing then decreasing all positive biases
#         [[257, 267, 277, 287, 297, 307, 317, 327, 337, 347, 357],
#          [0.2, 0.3, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9],
#          277],
#         # All positive biases, none equal
#         [[257, 267, 277, 287, 297, 307, 317, 327, 337, 347, 357],
#          [0.1, 0.2, 0.3, 0.6, 2.4, 4.6, 5.7, 6.2, 6.5, 6.8, 7.0],
#          257],
#         # All positive biases, first two equal
#         [[257, 267, 277, 287, 297, 307, 317, 327, 337, 347, 357],
#          [0.1, 0.1, 0.3, 0.6, 2.4, 4.6, 5.7, 6.2, 6.5, 6.8, 7.0],
#          267],
#         # All positive biases, first three equal
#         [[257, 267, 277, 287, 297, 307, 317, 327, 337, 347, 357],
#          [0.2, 0.2, 0.2, 0.6, 2.4, 4.6, 5.7, 6.2, 6.5, 6.8, 7.0],
#          277],
#         # All negative biases will return a masked out pixel
#         [[257, 267, 277, 287, 297, 307, 317, 327, 337, 347, 357],
#          [-1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1],
#          None],
#         # Normal bias profile, decreasing bias at high end
#         [[257, 267, 277, 287, 297, 307, 317, 327, 337, 347, 357],
#          [-0.2, -0.1, 0.1, 0.6, 2.4, 4.6, 5.7, 6.2, 6.5, 6.8, 6.0],
#          272],
#         # All positive biases, then decreasing bias at high end
#         [[257, 267, 277, 287, 297, 307, 317, 327, 337, 347, 357],
#          [0.1, 0.2, 0.3, 0.6, 2.4, 4.6, 5.7, 6.2, 6.5, 6.8, 6.0],
#          257],
#         # False/early transition with smaller bias than main transition
#         [[257, 267, 277, 287, 297, 307, 317, 327, 337, 347, 357],
#          [-0.1, 0.1, -0.2, -0.3, 0.6, 2.4, 4.6, 5.7, 6.2, 6.5, 6.8],
#          292],
#     ]
# )
# def test_ta_min_bias(ta_list, bias_list, expected, tol=0.0001):
#     ta_image_list = [
#         ee.Image.constant(ta).rename(['step_{:02d}_ta'.format(i+1)])
#         for i, ta in enumerate(ta_list)]
#     bias_image_list = [
#         ee.Image.constant(bias).rename(['step_{:02d}_bias'.format(i+1)])
#         for i, bias in enumerate(bias_list)]
#     input_img = ee.Image(ta_image_list + bias_image_list)
#     output = utils.constant_image_value(ta_min_bias(input_img))['ta']
#     if expected is None:
#         assert output is None
#     else:
#         assert abs(output - expected) <= tol


def mgrs_export_tiles(study_area_coll_id, mgrs_coll_id,
                      study_area_property=None, study_area_features=[],
                      mgrs_tiles=[], mgrs_skip_list=[],
                      utm_zones=[], wrs2_tiles=[],
                      mgrs_property='mgrs', utm_property='utm',
                      wrs2_property='wrs2'):
    """Select MGRS tiles and metadata that intersect the study area geometry

    Parameters
    ----------
    study_area_coll_id : str
        Study area feature collection asset ID.
    mgrs_coll_id : str
        MGRS feature collection asset ID.
    study_area_property : str, optional
        Property name to use for inList() filter call of study area collection.
        Filter will only be applied if both 'study_area_property' and
        'study_area_features' parameters are both set.
    study_area_features : list, optional
        List of study area feature property values to filter on.
    mgrs_tiles : list, optional
        User defined MGRS tile subset.
    mgrs_skip_list : list, optional
        User defined list MGRS tiles to skip.
    utm_zones : list, optional
        User defined UTM zone subset.
    wrs2_tiles : list, optional
        User defined WRS2 tile subset.
    mgrs_property : str, optional
        MGRS property in the MGRS feature collection (the default is 'mgrs').
    utm_property : str, optional
        UTM zone property in the MGRS feature collection (the default is 'utm').
    wrs2_property : str, optional
        WRS2 property in the MGRS feature collection (the default is 'wrs2').

    Returns
    ------
    list of dicts: export information

    """
    # Build and filter the study area feature collection
    logging.debug('Building study area collection')
    logging.debug('  {}'.format(study_area_coll_id))
    study_area_coll = ee.FeatureCollection(study_area_coll_id)
    if (study_area_property == 'STUSPS' and
            'CONUS' in [x.upper() for x in study_area_features]):
        # Exclude AK, HI, AS, GU, PR, MP, VI, (but keep DC)
        study_area_features = [
            'AL', 'AR', 'AZ', 'CA', 'CO', 'CT', 'DC', 'DE', 'FL', 'GA',
            'IA', 'ID', 'IL', 'IN', 'KS', 'KY', 'LA', 'MA', 'MD', 'ME',
            'MI', 'MN', 'MO', 'MS', 'MT', 'NC', 'ND', 'NE', 'NH', 'NJ',
            'NM', 'NV', 'NY', 'OH', 'OK', 'OR', 'PA', 'RI', 'SC', 'SD',
            'TN', 'TX', 'UT', 'VA', 'VT', 'WA', 'WI', 'WV', 'WY']
    # elif (study_area_property == 'STUSPS' and
    #         'WESTERN11' in [x.upper() for x in study_area_features]):
    #     study_area_features = [
    #         'AZ', 'CA', 'CO', 'ID', 'MT', 'NM', 'NV', 'OR', 'UT', 'WA', 'WY']
    study_area_features = sorted(list(set(study_area_features)))

    if study_area_property and study_area_features:
        logging.debug('  Filtering study area collection')
        logging.debug('  Property: {}'.format(study_area_property))
        logging.debug('  Features: {}'.format(','.join(study_area_features)))
        study_area_coll = study_area_coll.filter(
            ee.Filter.inList(study_area_property, study_area_features))

    logging.debug('Building MGRS tile list')
    tiles_coll = ee.FeatureCollection(mgrs_coll_id) \
        .filterBounds(study_area_coll.geometry())

    # Filter collection by user defined lists
    if utm_zones:
        logging.debug('  Filter user UTM Zones:    {}'.format(utm_zones))
        tiles_coll = tiles_coll.filter(ee.Filter.inList(utm_property, utm_zones))
    if mgrs_skip_list:
        logging.debug('  Filter MGRS skip list:    {}'.format(mgrs_skip_list))
        tiles_coll = tiles_coll.filter(
            ee.Filter.inList(mgrs_property, mgrs_skip_list).Not())
    if mgrs_tiles:
        logging.debug('  Filter MGRS tiles/zones:  {}'.format(mgrs_tiles))
        # Allow MGRS tiles to be subsets of the full tile code
        #   i.e. mgrs_tiles = 10TE, 10TF
        mgrs_filters = [
            ee.Filter.stringStartsWith(mgrs_property, mgrs_id.upper())
            for mgrs_id in mgrs_tiles]
        tiles_coll = tiles_coll.filter(ee.call('Filter.or', mgrs_filters))

    def drop_geometry(ftr):
        return ee.Feature(None).copyProperties(ftr)

    logging.debug('  Requesting tile/zone info')
    tiles_info = utils.get_info(tiles_coll.map(drop_geometry))

    # Constructed as a list of dicts to mimic other interpolation/export tools
    tiles_list = []
    for tile_ftr in tiles_info['features']:
        tiles_list.append({
            'index': tile_ftr['properties']['mgrs'].upper(),
            'wrs2_tiles': sorted(utils.wrs2_str_2_set(
                tile_ftr['properties'][wrs2_property])),
        })

    # Apply the user defined WRS2 tile list
    if wrs2_tiles:
        logging.debug('  Filter WRS2 tiles: {}'.format(wrs2_tiles))
        for tile in tiles_list:
            tile['wrs2_tiles'] = sorted(list(
                set(tile['wrs2_tiles']) & set(wrs2_tiles)))

    # Only return export tiles that have intersecting WRS2 tiles
    export_list = [
        tile for tile in sorted(tiles_list, key=lambda k: k['index'])
        if tile['wrs2_tiles']]

    return export_list


# TODO: Move to openet.core.utils?
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


# DEADBEEF - Was in utils.py
def read_ini(ini_path):
    logging.debug('\nReading Input File')
    # Open config file
    config = configparser.ConfigParser()
    try:
        config.read(ini_path)
    except Exception as e:
        logging.error(
            '\nERROR: Input file could not be read, '
            'is not an input file, or does not exist\n'
            '  ini_path={}\n\nException: {}'.format(ini_path, e))
        import sys
        sys.exit()

    # Force conversion of unicode to strings
    ini = dict()
    for section in config.keys():
        ini[str(section)] = {}
        for k, v in config[section].items():
            ini[str(section)][str(k)] = v
    return ini


# CGM - This is a modified copy of openet.utils.delay_task()
#   It was changed to take and return the number of ready tasks
#   This change may eventually be pushed to openet.utils.delay_task()
def delay_task(delay_time=0, task_max=-1, task_count=0):
    """Delay script execution based on number of READY tasks

    Parameters
    ----------
    delay_time : float, int
        Delay time in seconds between starting export tasks or checking the
        number of queued tasks if "ready_task_max" is > 0.  The default is 0.
        The delay time will be set to a minimum of 10 seconds if
        ready_task_max > 0.
    task_max : int, optional
        Maximum number of queued "READY" tasks.
    task_count : int
        The current/previous/assumed number of ready tasks.
        Value will only be updated if greater than or equal to ready_task_max.

    Returns
    -------
    int : ready_task_count

    """
    if task_max > 3000:
        raise ValueError('The maximum number of queued tasks must be less than 3000')

    # Force delay time to be a positive value since the parameter used to
    #   support negative values
    if delay_time < 0:
        delay_time = abs(delay_time)

    if ((task_max is None or task_max <= 0) and (delay_time >= 0)):
        # Assume task_max was not set and just wait the delay time
        logging.debug(f'  Pausing {delay_time} seconds, not checking task list')
        time.sleep(delay_time)
        return 0
    elif task_max and (task_count < task_max):
        # Skip waiting or checking tasks if a maximum number of tasks was set
        #   and the current task count is below the max
        logging.debug(f'  Ready tasks: {task_count}')
        return task_count

    # If checking tasks, force delay_time to be at least 10 seconds if
    #   ready_task_max is set to avoid excessive EE calls
    delay_time = max(delay_time, 10)

    # Make an initial pause before checking tasks lists to allow
    #   for previous export to start up
    # CGM - I'm not sure what a good default first pause time should be,
    #   but capping it at 30 seconds is probably fine for now
    logging.debug(f'  Pausing {min(delay_time, 30)} seconds for tasks to start')
    time.sleep(delay_time)

    # If checking tasks, don't continue to the next export until the number
    #   of READY tasks is greater than or equal to "ready_task_max"
    while True:
        ready_task_count = len(utils.get_ee_tasks(states=['READY']).keys())
        logging.debug(f'  Ready tasks: {ready_task_count}')
        if ready_task_count >= task_max:
            logging.debug(f'  Pausing {delay_time} seconds')
            time.sleep(delay_time)
        else:
            logging.debug(f'  {task_max - ready_task_count} open task '
                          f'slots, continuing processing')
            break

    return ready_task_count


def arg_parse():
    """"""
    parser = argparse.ArgumentParser(
        description='Compute/export WRS2 Ta images',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '-i', '--ini', type=utils.arg_valid_file,
        help='Input file', metavar='FILE')
    parser.add_argument(
        '--debug', default=logging.INFO, const=logging.DEBUG,
        help='Debug level logging', action='store_const', dest='loglevel')
    parser.add_argument(
        '--delay', default=0, type=float,
        help='Delay (in seconds) between each export tasks')
    parser.add_argument(
        '--key', type=utils.arg_valid_file, metavar='FILE',
        help='JSON key file')
    parser.add_argument(
        '--overwrite', default=False, action='store_true',
        help='Force overwrite of existing files')
    parser.add_argument(
        '--ready', default=-1, type=int,
        help='Maximum number of queued READY tasks')
    parser.add_argument(
        '--recent', default=0, type=int,
        help='Number of days to process before current date '
             '(ignore INI start_date and end_date')
    parser.add_argument(
        '--reverse', default=False, action='store_true',
        help='Process WRS2 tiles in reverse order')
    parser.add_argument(
        '--tiles', default='', nargs='+',
        help='Comma/space separated list of tiles to process')
    parser.add_argument(
        '--update', default=False, action='store_true',
        help='Update images with older model version numbers')
    parser.add_argument(
        '--start', type=utils.arg_valid_date, metavar='DATE', default=None,
        help='Start date (format YYYY-MM-DD)')
    parser.add_argument(
        '--end', type=utils.arg_valid_date, metavar='DATE', default=None,
        help='End date (format YYYY-MM-DD)')
    args = parser.parse_args()

    return args


if __name__ == "__main__":
    args = arg_parse()

    logging.basicConfig(level=args.loglevel, format='%(message)s')
    logging.getLogger('googleapiclient').setLevel(logging.ERROR)

    main(ini_path=args.ini, overwrite_flag=args.overwrite,
         delay_time=args.delay, gee_key_file=args.key, ready_task_max=args.ready,
         reverse_flag=args.reverse, tiles=args.tiles, update_flag=args.update,
         recent_days=args.recent, start_dt=args.start, end_dt=args.end,
    )
