import argparse
from builtins import input
from collections import defaultdict
import configparser
from datetime import datetime, timedelta, timezone
from importlib import metadata
import json
# import logging
import math
import os
import pprint
import re
import time

import ee
import pandas as pd

import openet.disalexi
import openet.core
import openet.core.utils as utils

TOOL_NAME = 'tair_image_wrs2_export'
TOOL_VERSION = '0.3.0'

if 'FUNCTION_REGION' in os.environ:
    # Logging is not working correctly in cloud functions for Python 3.8+
    # Following workflow suggested in this issue:
    # https://issuetracker.google.com/issues/124403972
    import google.cloud.logging
    log_client = google.cloud.logging.Client(project='openet')
    log_client.setup_logging(log_level=20)
    import logging
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
else:
    import logging
    # logging.basicConfig(level=logging.INFO, format='%(message)s')
    logging.getLogger('earthengine-api').setLevel(logging.INFO)
    logging.getLogger('googleapiclient').setLevel(logging.INFO)
    logging.getLogger('requests').setLevel(logging.INFO)
    logging.getLogger('urllib3').setLevel(logging.INFO)


def main(
        ini_path=None,
        overwrite_flag=False,
        delay_time=0,
        gee_key_file=None,
        project_id=None,
        ready_task_max=-1,
        reverse_flag=False,
        tiles=None,
        update_flag=False,
        log_tasks=False,
        recent_days=None,
        start_dt=None,
        end_dt=None,
):
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
    project_id : str, optional
        Google cloud project ID to use for GEE authentication.
        This will be checked after the gee_key_file and before the user Initialize.
        The default is None.
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
        If True, log task information to the datastore (the default is False).
    recent_days : int, str, optional
        Limit start/end date range to this many days before the current date.
        The default is None which will use the INI start/end date directly.
    start_dt : datetime, optional
        Override the start date in the INI file
        (the default is None which will use the INI start date).
    end_dt : datetime, optional
        Override the (inclusive) end date in the INI file
        (the default is None which will use the INI end date).

    """
    logging.info(f'\n{"#" * 80}')
    logging.info('Compute WRS2 Ta images')

    ee.data.setWorkloadTag('disalexi-tair-scene-export')

    wrs2_tile_fmt = 'p{:03d}r{:03d}'
    wrs2_tile_re = re.compile('p?(\\d{1,3})r?(\\d{1,3})')

    # List of path/rows to skip
    wrs2_skip_list = [
        'p049r026',  # Vancouver Island, Canada
        'p048r028',  # OR/WA Coast
        # 'p047r031',  # North California coast
        'p042r037',  # San Nicholas Island, California
        'p040r038', 'p039r038', 'p038r038',  # Mexico (by California)
        'p037r039', 'p036r039', 'p035r039',  # Mexico (by Arizona)
        'p034r039', 'p033r039',  # Mexico (by New Mexico)
        'p032r040',  # Mexico (West Texas)
        'p029r041', 'p028r042', 'p027r043', 'p026r043',  # Mexico (South Texas)
        'p019r040', 'p018r040',  # West Florida coast
        'p016r043', 'p015r043',  # South Florida coast
        'p014r041', 'p014r042', 'p014r043',  # East Florida coast
        'p013r035', 'p013r036',  # North Carolina Outer Banks
        'p013r026', 'p012r026',  # Canada (by Maine)
        'p011r032',  # Rhode Island coast
        'p010r030',  # Maine
    ]
    wrs2_path_skip_list = [9, 49]
    wrs2_row_skip_list = [25, 24, 43]
    mgrs_skip_list = []
    date_skip_list = []

    export_id_fmt = '{model}_{index}'

    # TODO: Move to INI file
    # mask_ocean_flag = True
    # clip_ocean_flag = False

    # Read config file
    logging.info(f'  {os.path.basename(ini_path)}')
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
    logging.info(f'  ET Model: {model_name}')

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
        logging.info('\nINPUTS collections parameter was net set, '
                        'default to Landsat 5/7/8/9 C02 L2 collections')
        collections = ['LANDSAT/LC09/C02/T1_L2', 'LANDSAT/LC08/C02/T1_L2',
                       'LANDSAT/LE07/C02/T1_L2', 'LANDSAT/LT05/C02/T1_L2']
    except Exception as e:
        raise e

    try:
        export_coll_id = ini['EXPORT']['export_coll']
    except KeyError:
        raise ValueError('"export_coll" was not set in INI')
    except Exception as e:
        raise e

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
        scene_id_skip_path = str(ini['INPUTS']['scene_skip_list'])
    except KeyError:
        scene_id_skip_path = None
        logging.debug('  scene_skip_list: not set in INI, defaulting to None')
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
        logging.debug(f'  mgrs_tiles: {mgrs_tiles}')
    except KeyError:
        mgrs_tiles = []
        logging.debug('  mgrs_tiles: not set in INI, defaulting to []')
    except Exception as e:
        raise e

    try:
        utm_zones = str(ini['EXPORT']['utm_zones'])
        utm_zones = sorted([int(x.strip()) for x in utm_zones.split(',')])
        logging.debug(f'  utm_zones: {utm_zones}')
    except KeyError:
        utm_zones = []
        logging.debug('  utm_zones: not set in INI, defaulting to []')
    except Exception as e:
        raise e

    try:
        export_id_name = '_' + str(ini['EXPORT']['export_id_name'])
    except KeyError:
        export_id_name = ''
        logging.debug('  export_id_name: not set in INI, defaulting to ""')
    except Exception as e:
        raise e

    try:
        retile = int(ini['TAIR']['retile'])
    except KeyError:
        retile = 4
        logging.debug(f'  TAIR retile: not set in INI, defaulting to {retile}')
    except Exception as e:
        raise e

    # Model metadata version
    # Newer versions of the models will have the version set in the metadata
    # but try reading from the model version attribute for older model versions
    model_package_name = openet.disalexi.__name__.replace('.', '-')
    model_metadata = metadata.metadata(model_package_name)
    if 'Name' not in model_metadata.keys():
        model_metadata['Name'] = model_name
    if 'Version' not in model_metadata.keys():
        try:
            model_metadata['Version'] = openet.disalexi.__version__
        except AttributeError:
            raise Exception(f'Could not determine the model version')
    logging.info(f'\nModel name: {model_metadata["Name"]}')
    logging.info(f'Model ver.: {model_metadata["Version"]}')

    # If the user set the tiles argument, use these instead of the INI values
    if tiles:
        logging.info('\nOverriding INI mgrs_tiles and utm_zones parameters')
        logging.info(f'  user tiles: {tiles}')
        mgrs_tiles = sorted([y.strip() for x in tiles for y in x.split(',')])
        mgrs_tiles = [x.upper() for x in mgrs_tiles if x]
        logging.info(f'  mgrs_tiles: {", ".join(mgrs_tiles)}')
        utm_zones = sorted(list(set([int(x[:2]) for x in mgrs_tiles])))
        logging.info(f'  utm_zones:  {", ".join(map(str, utm_zones))}')

    today_dt = datetime.now()
    today_dt = today_dt.replace(hour=0, minute=0, second=0, microsecond=0)
    if start_dt and end_dt:
        # Attempt to use the function start/end dates
        logging.info('\nOverriding INI "start_date" and "end_date" parameters')
        logging.info('  Custom date range')
        start_date = start_dt.strftime('%Y-%m-%d')
        end_date = end_dt.strftime('%Y-%m-%d')
    elif recent_days:
        logging.debug('\nOverriding INI "start_date" and "end_date" parameters')
        logging.debug(f'  Recent days: {recent_days}')
        recent_days = list(sorted(utils.parse_int_set(recent_days)))
        # Assume that a single day value should actually be a range?
        if len(recent_days) == 1:
            recent_days = list(range(1, recent_days[0]))
        end_dt = today_dt - timedelta(days=recent_days[0])
        start_dt = today_dt - timedelta(days=recent_days[-1])
        start_date = start_dt.strftime('%Y-%m-%d')
        end_date = end_dt.strftime('%Y-%m-%d')
    else:
        # Parse the INI start/end dates
        logging.info('\nINI date range')
        try:
            start_dt = datetime.strptime(start_date, '%Y-%m-%d')
            end_dt = datetime.strptime(end_date, '%Y-%m-%d')
        except Exception as e:
            raise e
    logging.info(f'  Start: {start_date}')
    logging.info(f'  End:   {end_date}')

    # TODO: Add a few more checks on the dates
    if end_dt < start_dt:
        raise ValueError('end date can not be before start date')

    logging.debug('\nFilter date range')
    iter_start_dt = start_dt
    iter_end_dt = end_dt + timedelta(days=1)
    logging.debug(f'  Start: {iter_start_dt.strftime("%Y-%m-%d")}')
    logging.debug(f'  End:   {iter_end_dt.strftime("%Y-%m-%d")}')

    model_args = {
        k.lower(): float(v) if utils.is_number(v) else v
        for k, v in dict(ini['DISALEXI']).items()
    }

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
    #     for k, v in dict(ini['TAIR']).items()
    # }

    logging.info('\nDISALEXI Parameters')
    if 'stability_iterations' in model_args.keys():
        logging.info(f'  Stabil iter: {int(model_args["stability_iterations"])}')
    logging.info(f'  Albedo iter: {int(model_args["albedo_iterations"])}')

    logging.info('\nTAIR Parameters')
    logging.info(f'  Offsets:  {tair_args["offsets"]}')
    logging.debug(f'  Retile:   {retile}')

    # Read the scene ID skip list
    if (not scene_id_skip_path) or scene_id_skip_path.lower() in ['none', '']:
        logging.info(f'\nScene ID skip list not set')
        scene_id_skip_list = []
    elif scene_id_skip_path and scene_id_skip_path.startswith('https://'):
        logging.info(f'\nScene ID skip URL: {scene_id_skip_path}')
        scene_id_skip_list = {
            scene_id.upper() for scene_id in
            pd.read_csv(scene_id_skip_path)['SCENE_ID'].values
            if re.match('L[TEC]0[45789]_\d{3}\d{3}_\d{8}', scene_id)
        }
        logging.info(f'  Skip list count: {len(scene_id_skip_list)}')
        # print(list(scene_id_skip_list)[:10])
    elif (scene_id_skip_path and
          scene_id_skip_path.endswith('.csv') and
          os.path.isfile(scene_id_skip_path)
    ):
        logging.info(f'\nScene ID skip path: {scene_id_skip_path}')
        scene_id_skip_list = {
            scene_id.upper() for scene_id in
            pd.read_csv(scene_id_skip_path)['SCENE_ID'].values
            if re.match('L[TEC]0[45789]_\d{3}\d{3}_\d{8}', scene_id)
        }
        logging.info(f'  Skip list count: {len(scene_id_skip_list)}')
        # print(list(scene_id_skip_list)[:10])
    else:
        raise Exception(f'Unsupported scene_skip_list parameter: {scene_id_skip_path}')

    # Setup datastore task logging
    if log_tasks:
        # Assume function is being run deployed as a cloud function
        #   and use the default credentials (should be the SA credentials)
        logging.debug('\nInitializing task datastore client')
        try:
            from google.cloud import datastore
            datastore_client = datastore.Client(project='openet-dri')
        except Exception as e:
            logging.info('  Task logging disabled, error setting up datastore client')
            log_tasks = False

    # Initialize Earth Engine
    if gee_key_file:
        logging.info(f'\nInitializing GEE using user key file: {gee_key_file}')
        try:
            ee.Initialize(ee.ServiceAccountCredentials('_', key_file=gee_key_file))
        except ee.ee_exception.EEException:
            logging.warning('Unable to initialize GEE using user key file')
            return False
    elif 'FUNCTION_REGION' in os.environ:
        # Assume code is deployed to a cloud function
        logging.debug(f'\nInitializing GEE using application default credentials')
        import google.auth
        credentials, project_id = google.auth.default(
            default_scopes=['https://www.googleapis.com/auth/earthengine']
        )
        ee.Initialize(credentials)
    elif project_id is not None:
        logging.info(f'\nInitializing Earth Engine using project credentials'
                     f'\n  Project ID: {project_id}')
        try:
            ee.Initialize(project=project_id)
        except Exception as e:
            logging.warning(f'\nUnable to initialize GEE using project ID\n  {e}')
            return False
    else:
        logging.info('\nInitializing Earth Engine using user credentials')
        ee.Initialize()
    # elif 'GOOGLE_APPLICATION_CREDENTIALS' in os.environ:
    #     logging.info(f'\nInitializing GEE using GOOGLE_APPLICATION_CREDENTIALS key')
    #     try:
    #         ee.Initialize(ee.ServiceAccountCredentials(
    #             "_", key_file=os.environ.get('GOOGLE_APPLICATION_CREDENTIALS')))
    #     except Exception as e:
    #         logging.warning('Unable to initialize GEE using '
    #                         'GOOGLE_APPLICATION_CREDENTIALS key file')
    #         return False

    # Build output collection and folder if necessary
    logging.debug(f'\nExport Collection: {export_coll_id}')
    if not ee.data.getInfo(export_coll_id.rsplit('/', 1)[0]):
        logging.debug('\nFolder does not exist and will be built'
                      '\n  {}'.format(export_coll_id.rsplit('/', 1)[0]))
        input('Press ENTER to continue')
        ee.data.createAsset({'type': 'FOLDER'}, export_coll_id.rsplit('/', 1)[0])

    if not ee.data.getInfo(export_coll_id):
        logging.info(f'\nExport collection does not exist and will be built\n  {export_coll_id}')
        input('Press ENTER to continue')
        ee.data.createAsset({'type': 'IMAGE_COLLECTION'}, export_coll_id)

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
        logging.info('\nChecking Task List')
        tasks = utils.get_ee_tasks()
        if logging.getLogger().getEffectiveLevel() == logging.DEBUG:
            utils.print_ee_tasks(tasks)
            input('ENTER')
        running_task_count = sum([1 for v in tasks.values() if v['state'] == 'RUNNING'])
        ready_task_count = sum([1 for v in tasks.values() if v['state'] == 'READY'])
        logging.info(f'  Running Tasks: {running_task_count}')
        logging.info(f'  Ready Tasks:   {ready_task_count}')

        # Hold the job here if the ready task count is already over the max
        ready_task_count = delay_task(
            delay_time=0, task_max=ready_task_max, task_count=ready_task_count
        )

    # Get an ET image to set the Ta values to
    logging.debug('\nALEXI ET properties')
    alexi_coll_id = model_args['alexi_source']
    alexi_crs = 'EPSG:4326'
    if ((alexi_coll_id.upper() == 'CONUS_V006') or
            alexi_coll_id.endswith('projects/ee-tulipyangyun-2/assets/alexi/ALEXI_V006')):
        alexi_coll_id = 'projects/ee-tulipyangyun-2/assets/alexi/ALEXI_V006'
        alexi_cs = 0.04
        alexi_x, alexi_y = -125.02, 49.78
        # alexi_geo = [0.04, 0.0, -125.02, 0.0, -0.04, 49.78]
    else:
        raise ValueError(f'unsupported ALEXI source: {alexi_coll_id}')

    logging.debug(f'  Collection: {alexi_coll_id}')


    # Get list of ALEXI dates (not using the system:index since it is YYYYJJJ)
    logging.info('\nChecking available ALEXI data')
    alexi_coll = (
        ee.ImageCollection(alexi_coll_id)
        .filterDate(iter_start_dt.strftime('%Y-%m-%d'), iter_end_dt.strftime('%Y-%m-%d'))
    )
    def set_date(x):
        return x.set('date', ee.Date(x.get('system:time_start')).format('yyyy-MM-dd'))
    alexi_dates = set(utils.get_info(alexi_coll.map(set_date).aggregate_array('date')))


    # Get list of dates when there is no windspeed to help pre-screen dates to skip
    # Using windspeed as a proxy for all of the meteorology variables
    if ('windspeed_source' in model_args.keys()) and (model_args['windspeed_source'] == 'CFSR'):
        logging.info('\nChecking available meteorology data')
        meteo_coll_id = 'projects/disalexi/meteo_data/windspeed/global_v001_3hour'
        meteo_coll = (
            ee.ImageCollection(meteo_coll_id)
            .filterDate(iter_start_dt.strftime('%Y-%m-%d'), iter_end_dt.strftime('%Y-%m-%d'))
        )
        def set_date(x):
            return x.set('date', ee.Date(x.get('system:time_start')).format('yyyy-MM-dd'))
        meteo_dates = set(utils.get_info(meteo_coll.map(set_date).aggregate_histogram('date').keys()))
    else:
        meteo_dates = {}


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


    # Process each MGRS grid zone tile separately
    logging.info('\nImage Exports')
    processed_image_ids = set()
    for export_info in sorted(export_list, key=lambda i: i['index'], reverse=reverse_flag):
        logging.info(f'{export_info["index"]}')

        tile_count = len(export_info['wrs2_tiles'])
        tile_list = sorted(export_info['wrs2_tiles'], reverse=not(reverse_flag))
        tile_geom = ee.Geometry.Rectangle(export_info['extent'], export_info['crs'], False)

        # Get the available image ID list for the zone
        # Get list of existing image assets and their properties for the zone
        # Process date range by years to avoid requesting more than 3000 images
        logging.debug('  Getting list of available model images and existing assets')
        export_image_id_list = []
        asset_props = {}
        for year_start_dt, year_end_dt in date_range_by_year(start_dt, end_dt, exclusive_end_dates=True):
            year_start_date = year_start_dt.strftime("%Y-%m-%d")
            year_end_date = year_end_dt.strftime("%Y-%m-%d")
            logging.debug(f'  {year_start_date} {year_end_date}')

            # Get the image and asset lists for the full zone
            # Collection end date is exclusive
            model_obj = openet.disalexi.Collection(
                collections=collections,
                cloud_cover_max=float(ini['INPUTS']['cloud_cover']),
                start_date=year_start_date,
                end_date=year_end_date,
                geometry=tile_geom.buffer(1000),
            )
            year_image_id_list = utils.get_info(ee.List(
                model_obj.overpass(variables=['ndvi']).aggregate_array('image_id')
            ))

            # Filter to the wrs2_tile list
            # The WRS2 tile filtering should be done in the Collection call above,
            #   but the DisALEXI model does not currently support this
            year_image_id_list = [
                x for x in year_image_id_list
                if 'p{}r{}'.format(*re.findall('_(\d{3})(\d{3})_', x)[0]) in tile_list
            ]

            # Filter image_ids that have already been processed as part of a
            #   different MGRS tile (might be faster with sets)
            year_image_id_list = [x for x in year_image_id_list if x not in processed_image_ids]
            # Keep track of all the image_ids that have been processed
            processed_image_ids.update(year_image_id_list)

            export_image_id_list.extend(year_image_id_list)

            # Get list of existing image assets and their properties
            logging.debug('  Getting GEE asset list')
            asset_coll = (
                ee.ImageCollection(export_coll_id)
                .filterDate(year_start_date, year_end_date)
                .filter(ee.Filter.inList('wrs2_tile', tile_list))
                .filterBounds(tile_geom)
            )
            year_asset_props = {
                f'{export_coll_id}/{x["properties"]["system:index"]}': x['properties']
                for x in utils.get_info(asset_coll)['features']
            }
            asset_props.update(year_asset_props)

        if not export_image_id_list:
            logging.info('  No Landsat images in date range, skipping zone')
            continue
        export_image_id_list = sorted(
            export_image_id_list,
            key=lambda k: k.split('/')[-1].split('_')[-1],
            reverse=reverse_flag
        )

        # Group images by wrs2 tile
        image_id_lists = defaultdict(list)
        for image_id in export_image_id_list:
            wrs2_tile = 'p{}r{}'.format(*wrs2_tile_re.findall(image_id.split('/')[-1].split('_')[1])[0])
            if wrs2_tile not in tile_list:
                continue
            image_id_lists[wrs2_tile].append(image_id)

        # Process exports by wrs2 tile
        for export_n, wrs2_tile in enumerate(tile_list):
            path, row = map(int, wrs2_tile_re.findall(wrs2_tile)[0])

            if wrs2_skip_list and (wrs2_tile in wrs2_skip_list):
                logging.debug('{} {} ({}/{}) - in wrs2 skip list'.format(
                    export_info['index'], wrs2_tile, export_n + 1, tile_count
                ))
                continue
            elif wrs2_row_skip_list and (row in wrs2_row_skip_list):
                logging.debug('{} {} ({}/{}) - in wrs2 row skip list'.format(
                    export_info['index'], wrs2_tile, export_n + 1, tile_count
                ))
                continue
            elif wrs2_path_skip_list and (path in wrs2_path_skip_list):
                logging.debug('{} {} ({}/{}) - in wrs2 path skip list'.format(
                    export_info['index'], wrs2_tile, export_n + 1, tile_count
                ))
                continue
            else:
                logging.debug('{} {} ({}/{})'.format(
                    export_info['index'], wrs2_tile, export_n + 1, tile_count
                ))

            # Subset the image ID list to the WRS2 tile
            try:
                image_id_list = image_id_lists[wrs2_tile]
            except KeyError:
                image_id_list = []
            if not image_id_list:
                logging.debug('  No Landsat images in date range, skipping tile')
                continue

            # Process each image in the collection by date
            # image_id is the full Earth Engine ID to the asset
            for image_id in image_id_list:
                coll_id, scene_id = image_id.rsplit('/', 1)
                l, p, r, year, month, day = parse_landsat_id(scene_id)
                image_dt = datetime.strptime('{:04d}{:02d}{:02d}'.format(year, month, day), '%Y%m%d')
                image_date = image_dt.strftime('%Y-%m-%d')

                landsat_img = ee.Image(image_id)
                landsat_dt = ee.Date(landsat_img.get('system:time_start'))

                export_id = export_id_fmt.format(
                    model=ini['INPUTS']['et_model'].lower(),
                    index=image_id.lower().replace('/', '_'),
                )
                export_id = export_id.replace('-', '')
                export_id += export_id_name
                asset_id = f'{export_coll_id}/{scene_id.lower()}'

                if date_skip_list and (image_date in date_skip_list):
                    logging.info(f'  {scene_id} - Date in skip list, skipping')
                    continue
                elif scene_id_skip_list and (scene_id.upper() in scene_id_skip_list):
                    logging.debug(f'  {scene_id} - Scene ID in skip list, skipping')
                    continue
                elif image_date not in alexi_dates:
                    logging.info(f'  {scene_id} - No ALEXI image in source, skipping')
                    continue
                elif meteo_dates and (image_date not in meteo_dates):
                    logging.info(f'  {scene_id} - No windspeed images for date, skipping')
                    continue

                if update_flag:
                    if export_id in tasks.keys():
                        logging.info(f'  {scene_id} - Task already submitted, skipping')
                        continue
                    elif asset_props and (asset_id in asset_props.keys()):
                        # In update mode only overwrite if the version is old
                        model_ver = version_number(model_metadata["Version"])
                        # model_ver = version_number(openet.disalexi.__version__)
                        asset_ver = version_number(asset_props[asset_id]['model_version'])

                        if asset_ver < model_ver:
                            logging.info(f'  {scene_id} - Existing asset model version is old, removing')
                            logging.debug(f'  asset: {asset_ver}\n  model: {model_ver}')
                            try:
                                ee.data.deleteAsset(asset_id)
                            except:
                                logging.info(f'  {scene_id} - Error removing asset, skipping')
                                continue
                        # elif (asset_props[asset_id]['alexi_source'] < model_args['alexi_source']):
                        #     logging.info('  ALEXI source is old, removing')
                        #     # input('ENTER')
                        #     try:
                        #         ee.data.deleteAsset(asset_id)
                        #     except:
                        #         logging.info('  Error removing asset, skipping')
                        #         continue
                        # elif (asset_props[asset_id]['build_date'] <= '2020-04-27'):
                        #     logging.info('  build_date is old, removing')
                        #     # input('ENTER')
                        #     try:
                        #         ee.data.deleteAsset(asset_id)
                        #     except:
                        #         logging.info('  Error removing asset, skipping')
                        #         continue
                        # elif (version_number(asset_props[asset_id]['tool_version']) <
                        #       version_number(TOOL_VERSION)):
                        #     logging.info('  Asset tool version is old, removing')
                        #     try:
                        #         ee.data.deleteAsset(asset_id)
                        #     except:
                        #         logging.info('  Error removing asset, skipping')
                        #         continue
                        else:
                            logging.info(f'  {scene_id} - Asset is up to date, skipping')
                            continue
                elif overwrite_flag:
                    if export_id in tasks.keys():
                        logging.info(f'  {scene_id} - Task already submitted, cancelling')
                        ee.data.cancelTask(tasks[export_id]['id'])
                        # ee.data.cancelOperation(tasks[export_id]['id'])
                    # This is intentionally not an "elif" so that a task can be
                    # cancelled and an existing image/file/asset can be removed
                    if asset_props and (asset_id in asset_props.keys()):
                        logging.info(f'  {scene_id} - Asset already exists, removing')
                        try:
                            ee.data.deleteAsset(asset_id)
                        except:
                            logging.info('  Error removing asset, skipping')
                            continue
                else:
                    if export_id in tasks.keys():
                        logging.debug(f'  {scene_id} - Task already submitted, skipping')
                        continue
                    elif asset_props and (asset_id in asset_props.keys()):
                        logging.debug(f'  {scene_id} - Asset already exists, skipping')
                        continue

                logging.debug(f'  Source: {image_id}')
                logging.debug(f'  Date:   {image_date}')
                # logging.debug(f'  DOY:    {doy}')
                logging.debug(f'  Export ID:  {export_id}')
                logging.debug(f'  Collection: {os.path.dirname(asset_id)}')
                logging.debug(f'  Image ID:   {os.path.basename(asset_id)}')

                # Manually check if the source LAI and LST images are present
                # Eventually this should/could be done inside the model instead
                if ('lai_source' in model_args.keys()) and (type(model_args['lai_source']) is str):
                    if model_args['lai_source'].lower() in ['openet-landsat-lai', 'openet-lai']:
                        landsat_lai_version = metadata.metadata('openet-landsat-lai')['Version']
                    else:
                        # Assumptions: string lai_source is an image collection ID
                        lai_coll = (
                            ee.ImageCollection(model_args['lai_source'])
                            .filterMetadata('scene_id', 'equals', scene_id)
                        )
                        if utils.get_info(lai_coll.size()) == 0:
                            logging.info(f'  {scene_id} - No LAI image in source, skipping')
                            continue
                        # TODO: If not exporting to COG, LAI version value could stay server size
                        lai_info = utils.get_info(ee.Image(lai_coll.first()))
                        landsat_lai_version = None
                        try:
                            landsat_lai_version = lai_info['properties']['landsat_lai_version']
                        except:
                            logging.info(f'  {scene_id} - Could not get LAI properties, skipping')
                            continue

                if ('lst_source' in model_args.keys()) and (type(model_args['lst_source']) is str):
                    # Assumptions: string lst_source is an image collection ID
                    lst_coll = (
                        ee.ImageCollection(model_args['lst_source'])
                        .filterMetadata('scene_id', 'equals', scene_id)
                    )
                    if utils.get_info(lst_coll.size()) == 0:
                        logging.info(f'  {scene_id} - No LST image in source, skipping')
                        continue
                    # TODO: If not exporting to COG, LST version value could stay server size
                    lst_info = utils.get_info(ee.Image(lst_coll.first()))
                    landsat_lst_version = None
                    try:
                        landsat_lst_version = lst_info['properties']['landsat_lst_version']
                    except:
                        logging.info(f'  {scene_id} - Could not get LST properties, skipping')
                        continue

                if (('windspeed_source' in model_args.keys()) and
                       (model_args['windspeed_source'] == 'CFSR')):
                    windspeed_coll_id = 'projects/disalexi/meteo_data/windspeed/global_v001_3hour'
                    windspeed_coll = (
                        ee.ImageCollection(windspeed_coll_id)
                        .filterDate(landsat_dt.advance(-3, 'hour'), landsat_dt.advance(3, 'hour'))
                        .select(['windspeed'])
                    )
                    try:
                        windspeed_count = utils.get_info(windspeed_coll.size())
                    except:
                        windspeed_count = -1
                    if windspeed_count <= 0:
                        logging.info(f'  {scene_id} - No windspeed images, skipping')
                        continue

                d_obj = openet.disalexi.Image.from_image_id(image_id, **model_args)

                # Compute the Ta mosaic and final Ta values from the initial Ta image
                #   computed using the direct approach
                export_img = d_obj.ta_coarse_mosaic(
                    ta_initial_img=d_obj.ta_coarse_initial(),
                    offsets=[int(t.strip()) for t in tair_args["offsets"].split(',')]
                )

                properties = {
                    # Custom properties
                    'build_date': datetime.today().strftime('%Y-%m-%d'),
                    'coll_id': coll_id,
                    'core_version': openet.core.__version__,
                    'image_id': image_id,
                    'landsat_lai_version': landsat_lai_version,
                    'landsat_lst_version': landsat_lst_version,
                    'model_name': model_metadata['Name'],
                    'model_version': model_metadata["Version"],
                    # 'model_name': model_name,
                    # 'model_version': openet.disalexi.__version__,
                    'scene_id': scene_id,
                    'tool_name': TOOL_NAME,
                    'tool_version': TOOL_VERSION,
                    'wrs2_tile': wrs2_tile_fmt.format(p, r),
                    # Source properties
                    # CGM - Note, setting the properties as server side objects
                    #   won't work for COG exports
                    'CLOUD_COVER': landsat_img.get('CLOUD_COVER'),
                    'CLOUD_COVER_LAND': landsat_img.get('CLOUD_COVER_LAND'),
                    # 'CLOUD_COVER': landsat_info['properties']['CLOUD_COVER'],
                    # 'CLOUD_COVER_LAND': landsat_info['properties']['CLOUD_COVER_LAND'],
                    # CGM - Should we use the Landsat time or the ALEXI time?
                    'system:time_start': landsat_img.get('system:time_start'),
                    # 'system:time_start': landsat_info['properties']['system:time_start'],
                    # 'system:time_start': utils.millis(image_dt),
                    # Other properties
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

                if retile and (retile in [2, 4, 8, 16, 32, 64, 128]):
                    export_img = export_img.retile(retile)

                # CGM: We could pre-compute (or compute once and then save)
                #   the crs, transform, and shape since they should (will?) be
                #   the same for each wrs2 tile
                # Build the export transform and shape from the Landsat image
                image_xy = utils.get_info(
                    landsat_img.geometry().bounds(1, 'EPSG:4326').coordinates().get(0)
                )
                if image_xy is None:
                    logging.info('  Could not get image extent, skipping')
                    continue
                extent = [
                    min([xy[0] for xy in image_xy]),
                    min([xy[1] for xy in image_xy]),
                    max([xy[0] for xy in image_xy]),
                    max([xy[1] for xy in image_xy])
                ]

                # Adjust extent to the cell size
                export_extent = [
                    round(math.floor((extent[0] - alexi_x) / alexi_cs) * alexi_cs + alexi_x, 8),
                    round(math.floor((extent[1] - alexi_y) / alexi_cs) * alexi_cs + alexi_y, 8),
                    round(math.ceil((extent[2] - alexi_x) / alexi_cs) * alexi_cs + alexi_x, 8),
                    round(math.ceil((extent[3] - alexi_y) / alexi_cs) * alexi_cs + alexi_y, 8)
                ]
                export_geo = [alexi_cs, 0, export_extent[0], 0, -alexi_cs, export_extent[3]]
                export_shape = [
                    int(abs(export_extent[2] - export_extent[0]) / alexi_cs),
                    int(abs(export_extent[3] - export_extent[1]) / alexi_cs)
                ]
                logging.debug(f'  CRS:    {alexi_crs}')
                logging.debug(f'  Extent: {export_extent}')
                logging.debug(f'  Geo:    {export_geo}')
                logging.debug(f'  Shape:  {export_shape}')

                # Build export tasks
                # max_retries = 4
                logging.debug('  Building export task')
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
                    #shardSize=64,
                )
                #     # except ee.ee_exception.EEException as e:
                #     except Exception as e:
                #         if ('Earth Engine memory capacity exceeded' in str(e) or
                #                 'Earth Engine capacity exceeded' in str(e)):
                #             logging.info(f'  Rebuilding task ({i}/{max_retries})')
                #             logging.debug(f'  {e}')
                #             time.sleep(i ** 2)
                #         else:
                #             logging.warning(f'Unhandled exception\n{e}')
                #             break
                #             raise e

                if not task:
                    logging.warning(f'  {scene_id} - Export task was not built, skipping')
                    continue

                logging.info(f'  {scene_id} - Starting export task')
                # TODO: We should not blindly keeping starting the task
                #   Need to only retry on specific errors, otherwise exit
                try:
                    task.start()
                except Exception as e:
                    logging.warning(f'  {scene_id} - Export task was not started, skipping')
                    logging.warning(f'  {scene_id} - {e}')
                    continue
                # for i in range(1, max_retries):
                #     try:
                #         task.start()
                #         break
                #     except Exception as e:
                #         logging.info(f'  Resending query ({i}/{max_retries})')
                #         logging.debug(f'  {e}')
                #         time.sleep(i ** 2)

                # # Not using ee_task_start since it doesn't return the task object
                # utils.ee_task_start(task)

                try:
                    task_id = task.status()['id']
                    logging.info(f'  {scene_id} - {task_id}')
                except:
                    logging.warning(f'  {scene_id} - No task ID')
                    continue

                # Write the export task info the openet-dri project datastore
                if log_tasks and task_id:
                    logging.debug('  Writing datastore entity')
                    try:
                        task_obj = datastore.Entity(
                            key=datastore_client.key('Task', task_id),
                            exclude_from_indexes=['properties'],
                        )
                        for k, v in task.status().items():
                            task_obj[k] = v
                        # task_obj['date'] = datetime.today().strftime('%Y-%m-%d')
                        # task_obj['eecu_seconds'] = 0
                        task_obj['eecu_hours'] = 0
                        # task_obj['eecu'] = 0
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
                        logging.warning(f'\nDatastore entity was not written\n{e}\n')

                # Pause before starting the next export task
                ready_task_count += 1
                ready_task_count = delay_task(
                    delay_time=delay_time,
                    task_max=ready_task_max,
                    task_count=ready_task_count,
                )

                logging.debug('')


def mgrs_export_tiles(
        study_area_coll_id,
        mgrs_coll_id,
        study_area_property=None,
        study_area_features=[],
        mgrs_tiles=[],
        mgrs_skip_list=[],
        utm_zones=[],
        wrs2_tiles=[],
        mgrs_property='mgrs',
        utm_property='utm',
        wrs2_property='wrs2',
):
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
    logging.debug(f'  {study_area_coll_id}')
    study_area_coll = ee.FeatureCollection(study_area_coll_id)
    if (study_area_property == 'STUSPS') and ('CONUS' in [x.upper() for x in study_area_features]):
        # Exclude AK, HI, AS, GU, PR, MP, VI, (but keep DC)
        study_area_features = [
            'AL', 'AR', 'AZ', 'CA', 'CO', 'CT', 'DC', 'DE', 'FL', 'GA',
            'IA', 'ID', 'IL', 'IN', 'KS', 'KY', 'LA', 'MA', 'MD', 'ME',
            'MI', 'MN', 'MO', 'MS', 'MT', 'NC', 'ND', 'NE', 'NH', 'NJ',
            'NM', 'NV', 'NY', 'OH', 'OK', 'OR', 'PA', 'RI', 'SC', 'SD',
            'TN', 'TX', 'UT', 'VA', 'VT', 'WA', 'WI', 'WV', 'WY',
        ]
    study_area_features = sorted(list(set(study_area_features)))

    if study_area_property and study_area_features:
        logging.debug('  Filtering study area collection')
        logging.debug(f'  Property: {study_area_property}')
        logging.debug(f'  Features: {",".join(study_area_features)}')
        study_area_coll = study_area_coll.filter(
            ee.Filter.inList(study_area_property, study_area_features)
        )

    logging.debug('Building MGRS tile list')
    tiles_coll = ee.FeatureCollection(mgrs_coll_id).filterBounds(study_area_coll.geometry())

    # Filter collection by user defined lists
    if utm_zones:
        logging.debug(f'  Filter user UTM Zones:    {utm_zones}')
        tiles_coll = tiles_coll.filter(ee.Filter.inList(utm_property, utm_zones))
    if mgrs_skip_list:
        logging.debug(f'  Filter MGRS skip list:    {mgrs_skip_list}')
        tiles_coll = tiles_coll.filter(ee.Filter.inList(mgrs_property, mgrs_skip_list).Not())
    if mgrs_tiles:
        logging.debug(f'  Filter MGRS tiles/zones:  {mgrs_tiles}')
        # Allow MGRS tiles to be subsets of the full tile code
        #   i.e. mgrs_tiles = 10TE, 10TF
        mgrs_filters = [
            ee.Filter.stringStartsWith(mgrs_property, mgrs_id.upper())
            for mgrs_id in mgrs_tiles
        ]
        tiles_coll = tiles_coll.filter(ee.call('Filter.or', mgrs_filters))

    def drop_geometry(ftr):
        return ee.Feature(None).copyProperties(ftr)

    logging.debug('  Requesting tile/zone info')
    tiles_info = utils.get_info(tiles_coll.map(drop_geometry))

    # Constructed as a list of dicts to mimic other interpolation/export tools
    tiles_list = []
    for tile_ftr in tiles_info['features']:
        tiles_list.append({
            'crs': 'EPSG:{:d}'.format(int(tile_ftr['properties']['epsg'])),
            'extent': [int(tile_ftr['properties']['xmin']),
                       int(tile_ftr['properties']['ymin']),
                       int(tile_ftr['properties']['xmax']),
                       int(tile_ftr['properties']['ymax'])],
            'index': tile_ftr['properties']['mgrs'].upper(),
            'wrs2_tiles': sorted(utils.wrs2_str_2_set(tile_ftr['properties'][wrs2_property])),
        })

    # Apply the user defined WRS2 tile list
    if wrs2_tiles:
        logging.debug(f'  Filter WRS2 tiles: {wrs2_tiles}')
        for tile in tiles_list:
            tile['wrs2_tiles'] = sorted(list(set(tile['wrs2_tiles']) & set(wrs2_tiles)))

    # Only return export tiles that have intersecting WRS2 tiles
    export_list = [tile for tile in sorted(tiles_list, key=lambda k: k['index']) if tile['wrs2_tiles']]

    return export_list


# TODO: Move to openet.core.utils?
def date_range_by_year(start_dt, end_dt, exclusive_end_dates=False):
    """

    Parameters
    ----------
    start_dt : datetime
    end_dt : datetime
    exclusive_end_dates : bool, optional
        If True, set the end dates for each iteration range to be exclusive.

    Returns
    -------
    list of start and end datetimes split by year

    """
    if (end_dt - start_dt).days > 366:
        for year in range(start_dt.year, end_dt.year+1):
            year_start_dt = max(datetime(year, 1, 1), start_dt)
            year_end_dt = datetime(year+1, 1, 1) - timedelta(days=1)
            year_end_dt = min(year_end_dt, end_dt)
            if exclusive_end_dates:
                year_end_dt = year_end_dt + timedelta(days=1)
            yield year_start_dt, year_end_dt
    else:
        if exclusive_end_dates:
            year_end_dt = end_dt + timedelta(days=1)
        yield start_dt, year_end_dt


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

    if ((task_max is None) or (task_max <= 0)) and (delay_time >= 0):
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


def version_number(version_str):
    return list(map(int, version_str.split('.')))


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
        '--log_tasks', default=False, action='store_true',
        help='Log tasks to the datastore')
    parser.add_argument(
        '--overwrite', default=False, action='store_true',
        help='Force overwrite of existing files')
    parser.add_argument(
        '--project', default=None,
        help='Google cloud project ID to use for GEE authentication')
    parser.add_argument(
        '--ready', default=-1, type=int,
        help='Maximum number of queued READY tasks')
    parser.add_argument(
        '--recent', default='',
        help='Day range (or number of days) to process before current date '
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

    main(
        ini_path=args.ini,
        overwrite_flag=args.overwrite,
        gee_key_file=args.key,
        project_id=args.project,
        delay_time=args.delay,
        ready_task_max=args.ready,
        reverse_flag=args.reverse,
        tiles=args.tiles,
        update_flag=args.update,
        log_tasks=args.log_tasks,
        recent_days=args.recent,
        start_dt=args.start,
        end_dt=args.end,
    )
