import argparse
from builtins import input
from collections import defaultdict
import datetime
import json
import logging
import math
import os
import pprint
import random
import re
import sys
import time

import ee
# from google.cloud import datastore
from osgeo import ogr, osr

import openet.disalexi
import openet.core
import utils
# from . import utils


def main(ini_path=None, overwrite_flag=False, delay_time=0, gee_key_file=None,
         max_ready=-1, reverse_flag=False, tiles=None, update_flag=False,
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
    max_ready: int, optional
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
    wrs2_tile_re = re.compile('p?(\d{1,3})r?(\d{1,3})')

    # List of path/rows to skip
    wrs2_skip_list = [
        'p038r038', 'p039r038', 'p040r038',  # Mexico
        'p042r037',  # San Nicholas Island
        'p049r026',  # Vancouver Island
        # 'p041r037', 'p042r037', 'p047r031',  # CA Coast
    ]

    mgrs_skip_list = []

    export_id_fmt = '{model}_{index}'

    # TODO: Move to INI file
    clip_ocean_flag = True

    # Read config file
    ini = utils.read_ini(ini_path)
    # ini = configparser.ConfigParser(interpolation=None)
    # ini.read_file(open(ini_path, 'r'))

    # # Force conversion of unicode to strings
    # for section in ini.sections():
    #     ini[str(section)] = {}
    #     for k, v in ini[section].items():
    #         ini[str(section)][str(k)] = v

    model_name = 'DISALEXI'
    # model_name = ini['INPUTS']['ET_MODEL'].upper()

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
        logging.info('\nINPUTS collections parameter was net set, '
                        'default to Landsat 5/7/8 C01 SR collections')
        collections = ['LANDSAT/LC08/C01/T1_SR', 'LANDSAT/LE07/C01/T1_SR',
                       'LANDSAT/LT05/C01/T1_SR']
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
        wrs2_tiles = str(ini['INPUTS']['wrs2_tiles'])
        wrs2_tiles = sorted([x.strip() for x in wrs2_tiles.split(',')])
        wrs2_tiles = [x.replace('p', '').replace('r', '') for x in wrs2_tiles]
    except KeyError:
        logging.info('\nINPUTS wrs2_tiles parameter was net set, default to []')
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

    try:
        output_type = str(ini['EXPORT']['output_type'])
    except KeyError:
        output_type = 'float'
        # output_type = 'int16'
        logging.debug('  output_type: not set in INI, '
                      'defaulting to {}'.format(output_type))
    except Exception as e:
        raise e

    try:
        scale_factor = int(ini['EXPORT']['scale_factor'])
    except KeyError:
        scale_factor = 1
        # scale_factor = 10000
        logging.debug('  scale_factor: not set in INI, '
                      'defaulting to {}'.format(scale_factor))
    except Exception as e:
        raise e

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
        for k, v in dict(ini[model_name]).items()}
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

    logging.info('\nDISALEXI Parameters')
    logging.info('  Stabil iter: {}'.format(int(model_args['stabil_iterations'])))
    logging.info('  Albedo iter: {}'.format(int(model_args['albedo_iterations'])))

    logging.info('\nTAIR Parameters')
    logging.info('  Cell size:  {}'.format(tair_args['cell_size']))
    logging.info('  Retile:     {}'.format(tair_args['retile']))
    logging.info('  Min Iter:   {}'.format(tair_args['min_iteration']))
    logging.info('  Max Iter:   {}'.format(tair_args['max_iteration']))
    # logging.info('  Iteration:  {}'.format(tair_args['iteration']))

    # ta_values = list(range(tair_args['ta_start'], tair_args['ta_stop'],
    #                        tair_args['ta_step']))
    # logging.info('  Start: {}'.format(tair_args['ta_start']))
    # logging.info('  Stop:  {}'.format(tair_args['ta_stop']))
    # logging.info('  Step:  {}'.format(tair_args['ta_step']))
    # if 'ta_iterations' not in tair_args.keys():
    #     logging.info('  Iterations: {}'.format(tair_args['ta_iterations']))


    logging.info('\nInitializing Earth Engine')
    if gee_key_file:
        logging.info('  Using service account key file: {}'.format(gee_key_file))
        # The "EE_ACCOUNT" parameter is not used if the key file is valid
        ee.Initialize(ee.ServiceAccountCredentials(
            'deadbeef', key_file=gee_key_file))
    else:
        ee.Initialize()


    # # TODO: set datastore key file as a parameter?
    # datastore_key_file = 'openet-dri-datastore.json'
    # if log_tasks and not os.path.isfile(datastore_key_file):
    #     logging.info('Task logging disabled, datastore key does not exist')
    #     log_tasks = False
    #     # input('ENTER')
    # if log_tasks:
    #     logging.info('\nInitializing task datastore client')
    #     try:
    #         datastore_client = datastore.Client.from_service_account_json(
    #             datastore_key_file)
    #     except Exception as e:
    #         logging.error('{}'.format(e))
    #         return False


    # Get current running tasks
    tasks = utils.get_ee_tasks()
    if logging.getLogger().getEffectiveLevel() == logging.DEBUG:
        utils.print_ee_tasks()
        input('ENTER')


    # Build output collection and folder if necessary
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
    alexi_coll_id = ini['DISALEXI']['alexi_source']
    if alexi_coll_id.upper() == 'CONUS_V002':
        alexi_coll_id = 'projects/earthengine-legacy/assets/' \
                        'projects/disalexi/alexi/CONUS_V002'
        alexi_mask = ee.Image('projects/earthengine-legacy/assets/'
                              'projects/disalexi/alexi/conus_v002_mask')\
            .double().multiply(0)
    elif alexi_coll_id.upper() == 'CONUS_V001':
        alexi_coll_id = 'projects/earthengine-legacy/assets/' \
                        'projects/disalexi/alexi/CONUS_V001'
        alexi_mask = ee.Image('projects/earthengine-legacy/assets/'
                              'projects/disalexi/alexi/conus_v001_mask')\
            .double().multiply(0)
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


    # Get list of MGRS tiles that intersect the study area
    logging.debug('\nMGRS Tiles/Zones')
    export_list = mgrs_export_tiles(
        ini['INPUTS']['study_area_path'],
        mgrs_coll_id=mgrs_ftr_coll_id,
        mgrs_tiles=mgrs_tiles,
        mgrs_skip_list=mgrs_skip_list,
        utm_zones=utm_zones,
        wrs2_tiles=wrs2_tiles,
        # simplify_buffer=simplify_buffer,
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

            # filter_args = {}
            # for coll_id in collections:
            #     filter_args[coll_id] = [
            #         {'type': 'equals', 'leftField': 'WRS_PATH', 'rightValue': path},
            #         {'type': 'equals', 'leftField': 'WRS_ROW', 'rightValue': row}]
            # # logging.debug('  Filter Args: {}'.format(filter_args))

            # TODO: Switch to call to model.Collection()
            landsat_coll = ee.ImageCollection([])
            export_geom = ee.Geometry.Point(openet.core.wrs2.centroids[wrs2_tile])
            if ('LANDSAT/LC08/C01/T1_SR' in collections and
                    ini['INPUTS']['end_date'] > '2013-03-24'):
                l8_coll = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR')\
                    .filterDate(iter_start_dt.strftime('%Y-%m-%d'),
                                iter_end_dt.strftime('%Y-%m-%d'))\
                    .filterMetadata('WRS_PATH', 'equals', path)\
                    .filterMetadata('WRS_ROW', 'equals', row)\
                    .filter(ee.Filter.gt('system:time_start',
                                         ee.Date('2013-03-24').millis()))\
                    .filterMetadata('CLOUD_COVER_LAND', 'less_than',
                                    float(ini['INPUTS']['cloud_cover']))
                #     .filterBounds(export_geom)\
                landsat_coll = ee.ImageCollection(landsat_coll.merge(l8_coll))
            if ('LANDSAT/LE07/C01/T1_SR' in collections and
                    ini['INPUTS']['end_date'] >= '1999-01-01'):
                l7_coll = ee.ImageCollection('LANDSAT/LE07/C01/T1_SR')\
                    .filterDate(iter_start_dt.strftime('%Y-%m-%d'),
                                iter_end_dt.strftime('%Y-%m-%d'))\
                    .filterMetadata('WRS_PATH', 'equals', path)\
                    .filterMetadata('WRS_ROW', 'equals', row)\
                    .filterMetadata('CLOUD_COVER_LAND', 'less_than',
                                    float(ini['INPUTS']['cloud_cover']))
                #     .filterBounds(export_geom)\
                landsat_coll = ee.ImageCollection(landsat_coll.merge(l7_coll))
            if ('LANDSAT/LT05/C01/T1_SR' in collections and
                    ini['INPUTS']['start_date'] <= '2011-12-31'):
                l5_coll = ee.ImageCollection('LANDSAT/LT05/C01/T1_SR')\
                    .filterDate(iter_start_dt.strftime('%Y-%m-%d'),
                                iter_end_dt.strftime('%Y-%m-%d'))\
                    .filterMetadata('WRS_PATH', 'equals', path)\
                    .filterMetadata('WRS_ROW', 'equals', row)\
                    .filter(ee.Filter.lt('system:time_start',
                                         ee.Date('2011-12-31').millis()))\
                    .filterMetadata('CLOUD_COVER_LAND', 'less_than',
                                    float(ini['INPUTS']['cloud_cover']))
                #     .filterBounds(export_geom)\
                landsat_coll = ee.ImageCollection(landsat_coll.merge(l5_coll))
            if ('LANDSAT/LT04/C01/T1_SR' in collections and
                    ini['INPUTS']['start_date'] <= '1993-12-01'):
                l4_coll = ee.ImageCollection('LANDSAT/LT04/C01/T1_SR')\
                    .filterDate(iter_start_dt.strftime('%Y-%m-%d'),
                                iter_end_dt.strftime('%Y-%m-%d'))\
                    .filterMetadata('WRS_PATH', 'equals', path)\
                    .filterMetadata('WRS_ROW', 'equals', row)\
                    .filterMetadata('CLOUD_COVER_LAND', 'less_than',
                                    float(ini['INPUTS']['cloud_cover']))
                #     .filterBounds(export_geom)\
                landsat_coll = ee.ImageCollection(landsat_coll.merge(l4_coll))
            # DATA_TYPE filter only needed if using TOA collections
            #     .filterMetadata('DATA_TYPE', 'equals', 'L1TP')

            image_id_list = sorted(list(set(
                landsat_coll.aggregate_array('system:id').getInfo())))
            if not image_id_list:
                logging.debug('  Empty image ID list, exiting')
                return False

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

                # if update_flag:
                #     def version_number(version_str):
                #         return list(map(int, version_str.split('.')))
                #
                #     if export_id in tasks.keys():
                #         logging.info('    Task already submitted, skipping')
                #         continue
                #     # In update mode only overwrite if the version is old
                #     if asset_props and asset_id in asset_props.keys():
                #         model_ver = version_number(model.__version__)
                #         asset_ver = version_number(
                #             asset_props[asset_id]['model_version'])
                #         # asset_flt = [
                #         #     t == 'float' for b, t in asset_types.items()
                #         #     if b in ['et', 'et_reference']]
                #
                #         if asset_ver < model_ver:
                #             logging.info('    Existing asset model version is old, '
                #                          'removing')
                #             logging.debug(f'    asset: {asset_ver}\n'
                #                           f'    model: {model_ver}')
                #             try:
                #                 ee.data.deleteAsset(asset_id)
                #             except:
                #                 logging.info('    Error removing asset, skipping')
                #                 continue
                #         elif ((('T1_RT_TOA' in asset_props[asset_id]['coll_id']) and
                #                ('T1_RT_TOA' not in image_id)) or
                #               (('T1_RT' in asset_props[asset_id]['coll_id']) and
                #                ('T1_RT' not in image_id))):
                #             logging.info(
                #                 '    Existing asset is from realtime Landsat '
                #                 'collection, removing')
                #             try:
                #                 ee.data.deleteAsset(asset_id)
                #             except:
                #                 logging.info('    Error removing asset, skipping')
                #                 continue
                #         # elif (version_number(asset_props[asset_id]['tool_version']) <
                #         #       version_number(TOOL_VERSION)):
                #         #     logging.info('    Asset tool version is old, removing')
                #         #     try:
                #         #         ee.data.deleteAsset(asset_id)
                #         #     except:
                #         #         logging.info('    Error removing asset, skipping')
                #         #         continue
                #         # elif any(asset_flt):
                #         #     logging.info(
                #         #         '    Asset ET types are float, removing')
                #         #     ee.data.deleteAsset(asset_id)
                #         # elif 'tool_version' not in asset_props[asset_id].keys():
                #         #     logging.info('    TOOL_VERSION property was not set, removing')
                #         #     ee.data.deleteAsset(asset_id)
                #
                #         # elif asset_props[asset_id]['images'] == '':
                #         #     logging.info('    Images property was not set, removing')
                #         #     input('ENTER')
                #         #     ee.data.deleteAsset(asset_id)
                #         else:
                #             logging.info('    Asset is up to date, skipping')
                #             continue
                if overwrite_flag:
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
                    if asset_props and asset_id in asset_props.keys():
                        logging.info('    Asset already exists, skipping')
                        continue

                # CGM: We could pre-compute (or compute once and then save)
                #   the crs, transform, and shape since they should (will?) be
                #   the same for each wrs2 tile
                output_info = utils.get_info(ee.Image(image_id).select(['B2']))
                transform = '[{}]'.format(
                    ','.join(map(str, output_info['bands'][0]['crs_transform'])))






def mgrs_export_tiles(study_area_path, mgrs_coll_id, mgrs_tiles=[],
                      mgrs_skip_list=[], utm_zones=[], wrs2_tiles=[],
                      mgrs_property='mgrs', utm_property='utm',
                      wrs2_property='wrs2', simplify_buffer=0):
    """Select MGRS tiles and metadata that intersect the study area geometry

    Parameters
    ----------
    study_area_path : str
        File path of the study area shapefile.
    mgrs_coll_id : str
        MGRS feature collection asset ID.
    mgrs_tiles : list
        User defined MGRS tile subset.
    mgrs_skip_list : list
        User defined list MGRS tiles to skip.
    utm_zones : list
        User defined UTM zone subset.
    wrs2_tiles : list
        User defined WRS2 tile subset.
    mgrs_property : str, optional
        MGRS property in the MGRS feature collection (the default is 'mgrs').
    utm_property : str, optional
        UTM zone property in the MGRS feature collection (the default is 'wrs2').
    wrs2_property : str, optional
        WRS2 property in the MGRS feature collection (the default is 'wrs2').
    simplify_buffer : float, optional
        Study area simplify tolerance (the default is 0).
        Note, this distance is in the units of the study area shapefile.

    Returns
    ------
    list of dicts: export information

    """
    logging.info('\nReading study area shapefile')
    logging.info('  {}'.format(study_area_path))
    study_area_ds = ogr.Open(study_area_path, 0)
    study_area_lyr = study_area_ds.GetLayer()
    study_area_osr = study_area_lyr.GetSpatialRef()
    study_area_crs = str(study_area_osr.ExportToWkt())
    # study_area_proj4 = study_area_osr.ExportToProj4()
    logging.debug('  Study area projection: {}'.format(study_area_crs))

    # Get the dissolved/unioned geometry of the study area
    output_geom = ogr.Geometry(ogr.wkbMultiPolygon)
    for study_area_ftr in study_area_lyr:
        output_geom = output_geom.Union(study_area_ftr.GetGeometryRef())
    study_area_ds = None

    # # Project the study area geometry to the EPSG:3857
    # #   so units will be meters for buffering and simplifying
    # temp_crs = 'EPSG:3857'
    # temp_osr = osr.SpatialReference()
    # temp_osr.ImportFromEPSG(3857)
    # output_tx = osr.CoordinateTransformation(study_area_osr, temp_osr)
    # output_geom.Transform(output_tx)

    if simplify_buffer:
        output_geom = output_geom.SimplifyPreserveTopology(simplify_buffer) \
            .buffer(simplify_buffer)
    elif study_area_osr.IsGeographic():
        tol = 0.0000001
        logging.debug('  Simplifying study area geometry (tol={})'.format(tol))
        output_geom = output_geom.SimplifyPreserveTopology(tol)
        # output_geom = output_geom.SimplifyPreserveTopology(tol).buffer(tol)
    # else:
    # Added flatten call to change clockwise geometries to counter cw
    output_geom.FlattenTo2D()

    logging.debug('  Building GeoJSON')
    output_geojson = json.loads(output_geom.ExportToJson())

    logging.debug('  Building EE geometry')
    output_ee_geom = ee.Geometry(output_geojson, study_area_crs, False)

    logging.info('Building MGRS tile list')
    tiles_coll = ee.FeatureCollection(mgrs_coll_id) \
        .filterBounds(output_ee_geom)

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

    # pprint.pprint(export_list)
    # input('ENTER')

    return export_list


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
    # parser.add_argument(
    #     '--recent', default=0, type=int,
    #     help='Number of days to process before current date '
    #          '(ignore INI start_date and end_date')
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
        '-s', '--start', type=utils.arg_valid_date, metavar='DATE', default=None,
        help='Start date (format YYYY-MM-DD)')
    parser.add_argument(
        '-e', '--end', type=utils.arg_valid_date, metavar='DATE', default=None,
        help='End date (format YYYY-MM-DD)')
    args = parser.parse_args()

    return args


if __name__ == "__main__":
    args = arg_parse()

    logging.basicConfig(level=args.loglevel, format='%(message)s')
    logging.getLogger('googleapiclient').setLevel(logging.ERROR)

    main(ini_path=args.ini, overwrite_flag=args.overwrite,
         delay_time=args.delay, gee_key_file=args.key, max_ready=args.ready,
         reverse_flag=args.reverse, tiles=args.tiles, update_flag=args.update,
         recent_days=args.recent, start_dt=args.start, end_dt=args.end,
    )
