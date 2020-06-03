import argparse
from builtins import input
import configparser
import datetime
import json
import logging
import math
import os
import pprint
# import random
import re
import time

import ee
from osgeo import ogr, osr

import openet.disalexi
import openet.core
import openet.core.utils as utils

TOOL_NAME = 'tair_image_wrs2_export'
TOOL_VERSION = '0.1.5'


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
    try:
        model_name = str(ini['INPUTS']['et_model']).upper()
    except KeyError:
        raise ValueError('"et_model" parameter was not set in INI')
    except Exception as e:
        raise e
    logging.info('  ET Model: {}'.format(model_name))

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
        # wrs2_tiles = [x.replace('p', '').replace('r', '') for x in wrs2_tiles]
    except KeyError:
        logging.debug('  wrs2_tiles: not set in INI, defaulting to []')
        wrs2_tiles = []
    except Exception as e:
        raise e

    # try:
    #     simplify_buffer = float(ini['INPUTS']['simplify_buffer'])
    # except KeyError:
    #     simplify_buffer = 0
    #     logging.debug('  simplify_buffer: not set in INI, '
    #                   'defaulting to {}'.format(simplify_buffer))
    # except Exception as e:
    #     raise e

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
    logging.info('  Stabil iter: {}'.format(int(model_args['stabil_iterations'])))
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


    # TODO: set datastore key file as a parameter?
    datastore_key_file = 'openet-dri-datastore.json'
    if log_tasks and not os.path.isfile(datastore_key_file):
        logging.info('Task logging disabled, datastore key does not exist')
        log_tasks = False
        # input('ENTER')
    if log_tasks:
        logging.info('\nInitializing task datastore client')
        # TODO: Move to top and add to requirements.txt and environment.yaml
        from google.cloud import datastore
        try:
            datastore_client = datastore.Client.from_service_account_json(
                datastore_key_file)
        except Exception as e:
            logging.error('{}'.format(e))
            return False


    # Get current running tasks
    tasks = utils.get_ee_tasks()
    if logging.getLogger().getEffectiveLevel() == logging.DEBUG:
        logging.debug('  Tasks: {}'.format(len(tasks)))
        input('ENTER')


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
    elif alexi_coll_id.upper() == 'CONUS_V001':
        alexi_coll_id = 'projects/earthengine-legacy/assets/' \
                        'projects/disalexi/alexi/CONUS_V001'
        alexi_mask = ee.Image('projects/earthengine-legacy/assets/'
                              'projects/disalexi/alexi/conus_v001_mask')\
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

            # TODO: Switch to call to openet.disalexi.Collection()
            landsat_coll = ee.ImageCollection([])
            export_geom = ee.Geometry.Point(openet.core.wrs2.centroids[wrs2_tile])
            if ('LANDSAT/LC08/C01/T1_SR' in collections and
                    iter_end_dt.strftime('%Y-%m-%d') > '2013-03-24'):
                l8_coll = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR')\
                    .filterDate(iter_start_dt.strftime('%Y-%m-%d'),
                                iter_end_dt.strftime('%Y-%m-%d'))\
                    .filterMetadata('WRS_PATH', 'equals', path)\
                    .filterMetadata('WRS_ROW', 'equals', row)\
                    .filter(ee.Filter.gt('system:time_start',
                                         ee.Date('2013-03-24').millis()))\
                    .filterMetadata('CLOUD_COVER_LAND', 'less_than',
                                    float(ini['INPUTS']['cloud_cover']))\
                    .filterBounds(export_geom)
                landsat_coll = ee.ImageCollection(landsat_coll.merge(l8_coll))
            if ('LANDSAT/LE07/C01/T1_SR' in collections and
                    iter_end_dt.strftime('%Y-%m-%d') >= '1999-01-01'):
                l7_coll = ee.ImageCollection('LANDSAT/LE07/C01/T1_SR')\
                    .filterDate(iter_start_dt.strftime('%Y-%m-%d'),
                                iter_end_dt.strftime('%Y-%m-%d'))\
                    .filterMetadata('WRS_PATH', 'equals', path)\
                    .filterMetadata('WRS_ROW', 'equals', row)\
                    .filterMetadata('CLOUD_COVER_LAND', 'less_than',
                                    float(ini['INPUTS']['cloud_cover']))\
                    .filterBounds(export_geom)
                landsat_coll = ee.ImageCollection(landsat_coll.merge(l7_coll))
            if ('LANDSAT/LT05/C01/T1_SR' in collections and
                    iter_start_dt.strftime('%Y-%m-%d') <= '2011-12-31'):
                l5_coll = ee.ImageCollection('LANDSAT/LT05/C01/T1_SR')\
                    .filterDate(iter_start_dt.strftime('%Y-%m-%d'),
                                iter_end_dt.strftime('%Y-%m-%d'))\
                    .filterMetadata('WRS_PATH', 'equals', path)\
                    .filterMetadata('WRS_ROW', 'equals', row)\
                    .filter(ee.Filter.lt('system:time_start',
                                         ee.Date('2011-12-31').millis()))\
                    .filterMetadata('CLOUD_COVER_LAND', 'less_than',
                                    float(ini['INPUTS']['cloud_cover']))\
                    .filterBounds(export_geom)
                landsat_coll = ee.ImageCollection(landsat_coll.merge(l5_coll))
            if ('LANDSAT/LT04/C01/T1_SR' in collections and
                    iter_start_dt.strftime('%Y-%m-%d') <= '1993-12-01'):
                l4_coll = ee.ImageCollection('LANDSAT/LT04/C01/T1_SR')\
                    .filterDate(iter_start_dt.strftime('%Y-%m-%d'),
                                iter_end_dt.strftime('%Y-%m-%d'))\
                    .filterMetadata('WRS_PATH', 'equals', path)\
                    .filterMetadata('WRS_ROW', 'equals', row)\
                    .filterMetadata('CLOUD_COVER_LAND', 'less_than',
                                    float(ini['INPUTS']['cloud_cover']))\
                    .filterBounds(export_geom)
                landsat_coll = ee.ImageCollection(landsat_coll.merge(l4_coll))
            # DATA_TYPE filter only needed if using TOA collections
            #     .filterMetadata('DATA_TYPE', 'equals', 'L1TP')


            image_id_list = sorted(list(set(
                landsat_coll.aggregate_array('system:id').getInfo())))
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
                    # TODO - Check if selecting the 0 UTC time is intentional
                    ta_source_coll = ee.ImageCollection('NASA/NLDAS/FORA0125_H002')\
                        .filterDate(image_date, next_date)\
                        .select(['temperature'])
                    input_image = ee.Image(ta_source_coll.first())\
                        .add(273.15).subtract(40).floor()
                    ta_source_img = alexi_mask.add(input_image).rename(['ta'])
                else:
                    logging.debug('    Tair source: {}'.format(tair_args['source_coll']))
                    ta_source_coll = ee.ImageCollection(tair_args['source_coll'])\
                        .filterMetadata('image_id', 'equals', image_id)
                    if ta_source_coll.size().getInfo() == 0:
                        logging.info('  No Tair image in source coll, skipping')
                        # input('ENTER')
                        continue

                    # A lot code to figure out the starting Ta value
                    # This identifies the first Ta that has a positive bias and
                    #   a bias that is larger than the previous bias
                    # It then selects the Ta for the previous step
                    # This should bracket a bias of zero but it is not guaranteed
                    input_img = ee.Image(ta_source_coll.first())
                    ta_array = input_img.select('step_\\d+_ta').toArray()
                    bias_array = input_img.select('step_\\d+_bias').toArray()
                    diff = bias_array.arraySlice(0, 1)\
                        .subtract(bias_array.arraySlice(0, 0, -1))
                    index = diff.gt(0).And(bias_array.arraySlice(0, 1).gt(0))
                    # Intentionally use 0,0,-1 slice here (instead of 0,1)
                    #   to get Ta before bias goes positive
                    ta_source_img = ta_array.arraySlice(0, 0, -1).arrayMask(index)\
                        .arraySlice(0, 0, 1).arrayFlatten([['array']])\
                        .rename(['ta'])


                # Manually check if the source LAI and TIR images are present
                # Eventually this should/could be done inside the model instead
                if ('tir_source' in model_args.keys() and \
                        type(model_args['tir_source']) is str):
                    # Assumptions: string tir_source is an image collection ID
                    tir_coll = ee.ImageCollection(model_args['tir_source']) \
                        .filterMetadata('scene_id', 'equals', scene_id)
                    if tir_coll.size().getInfo() == 0:
                        logging.info('  No TIR image in source, skipping')
                        input('ENTER')
                        continue
                if ('lai_source' in model_args.keys() and \
                        type(model_args['lai_source']) is str):
                    # Assumptions: string lai_source is an image collection ID
                    lai_coll = ee.ImageCollection(model_args['lai_source']) \
                        .filterMetadata('scene_id', 'equals', scene_id)
                    if lai_coll.size().getInfo() == 0:
                        logging.info('  No LAI image in source, skipping')
                        input('ENTER')
                        continue


                landsat_img = ee.Image(image_id)
                # CGM: We could pre-compute (or compute once and then save)
                #   the crs, transform, and shape since they should (will?) be
                #   the same for each wrs2 tile
                output_info = utils.get_info(landsat_img.select(['B2']))

                d_obj = openet.disalexi.Image(
                    openet.disalexi.LandsatSR(landsat_img).prep(), **model_args)
                export_img = d_obj.ta_mosaic(
                    ta_img=ta_source_img,
                    step_size=tair_args['step_size'],
                    step_count=tair_args['step_count'])

                # pprint.pprint(export_img.getInfo())
                # input('ENTER')

                if tair_args['retile'] and tair_args['retile'] > 0:
                    export_img = export_img.retile(tair_args['retile'])

                properties = {
                    # Custom properties
                    'coll_id': coll_id,
                    'core_version': openet.core.__version__,
                    'date_ingested': datetime.datetime.today().strftime('%Y-%m-%d'),
                    'image_id': image_id,
                    'model_name': model_name,
                    'model_version': openet.disalexi.__version__,
                    'scene_id': scene_id,
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

                # Write the export task info the openet-dri project datastore
                if log_tasks:
                    logging.debug('    Writing datastore entity')
                    try:
                        task_obj = datastore.Entity(key=datastore_client.key(
                            'Task', task.status()['id']),
                            exclude_from_indexes=['properties'])
                        for k, v in task.status().items():
                            task_obj[k] = v
                        # task_obj['date'] = datetime.datetime.today() \
                        #     .strftime('%Y-%m-%d')
                        task_obj['index'] = properties.pop('wrs2_tile')
                        # task_obj['wrs2_tile'] = properties.pop('wrs2_tile')
                        task_obj['model_name'] = properties.pop('model_name')
                        # task_obj['model_version'] = properties.pop('model_version')
                        task_obj['runtime'] = 0
                        task_obj['start_timestamp_ms'] = 0
                        task_obj['tool_name'] = properties.pop('tool_name')
                        task_obj['properties'] = json.dumps(properties)
                        datastore_client.put(task_obj)
                    except Exception as e:
                        # CGM - The message/handling will probably need to be updated
                        #   We may want separate try/excepts on the create and the put
                        logging.warning('\nDatastore entity was not written')
                        logging.warning('{}\n'.format(e))

                # Pause before starting the next export task
                utils.delay_task(delay_time, max_ready)

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
         recent_days=0, start_dt=args.start, end_dt=args.end,
    )
