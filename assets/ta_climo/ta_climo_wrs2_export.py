import argparse
from builtins import input
import calendar
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

import openet.disalexi
import openet.core
import openet.core.utils as utils
import openet.lai

TOOL_NAME = 'tair_climo_wrs2_export'
TOOL_VERSION = '0.1.0'

if 'FUNCTION_REGION' in os.environ:
    # Logging is not working correctly in cloud functions for Python 3.8+
    # Following workflow suggested in this issue:
    # https://issuetracker.google.com/issues/124403972
    import google.cloud.logging
    log_client = google.cloud.logging.Client(project='openet')
    log_client.setup_logging(log_level=20)
    import logging
    # CGM - Not sure if these lines are needed or not
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
    start_dt : datetime, optional
        Override the start date in the INI file
        (the default is None which will use the INI start date).
    end_dt : datetime, optional
        Override the (inclusive) end date in the INI file
        (the default is None which will use the INI end date).

    """
    logging.info('\nCompute WRS2 Ta climatology images')

    ee.data.setWorkloadTag('disalexi-tair-climo-export')

    wrs2_tile_fmt = 'p{:03d}r{:03d}'
    wrs2_tile_re = re.compile('p?(\\d{1,3})r?(\\d{1,3})')

    # List of path/rows to skip
    wrs2_skip_list = [
        'p049r026',  # Vancouver Island, Canada
        # 'p047r031',  # North California coast
        'p042r037',  # San Nicholas Island, California
        # 'p041r037',  # South California coast
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
    ]
    wrs2_path_skip_list = [9, 49]
    wrs2_row_skip_list = [25, 24, 43]
    mgrs_skip_list = []
    date_skip_list = []

    export_id_fmt = '{model}_{index}'

    # TODO: Move to INI file
    # clip_ocean_flag = True

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

    export_coll_id = f'{ini["EXPORT"]["export_coll"]}'

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

    # try:
    #     output_type = str(ini['EXPORT']['output_type'])
    # except KeyError:
    #     output_type = 'float'
    #     # output_type = 'int16'
    #     logging.debug(f'  output_type: not set in INI, defaulting to {output_type}')
    # except Exception as e:
    #     raise e
    #
    # try:
    #     scale_factor = int(ini['EXPORT']['scale_factor'])
    # except KeyError:
    #     scale_factor = 1
    #     # scale_factor = 10000
    #     logging.debug(f'  scale_factor: not set in INI, defaulting to {scale_factor}')
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
    #     for k, v in dict(ini['TAIR']).items()}

    # DEADBEEF - Tair cell size parameter is not used anywhere
    # if 'cell_size' not in tair_args.keys():
    #     tair_args['cell_size'] = 30
    if 'retile' not in tair_args.keys():
        tair_args['retile'] = 0

    logging.info('\nTAIR Parameters')
    # DEADBEEF - Source and start are not used in Ta direct calculation
    logging.info(f'  Source:     {tair_args["source_coll"]}')
    logging.debug(f'  Retile:     {tair_args["retile"]}')
    logging.info(f'  Step Size:  {tair_args["step_size"]}')
    logging.debug(f'  Step Count: {tair_args["step_count"]}')


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
    if (alexi_coll_id.upper() == 'CONUS_V006' or
            alexi_coll_id.endswith('projects/ee-tulipyangyun-2/assets/alexi/ALEXI_V006')):
        # alexi_coll_id = 'projects/ee-tulipyangyun-2/assets/alexi/ALEXI_V006'
        # alexi_mask_id = 'projects/earthengine-legacy/assets/projects/disalexi/alexi/conus_v004_mask'
        # alexi_mask = ee.Image(alexi_mask_id).double().multiply(0)
        # alexi_geo = [0.04, 0.0, -125.02, 0.0, -0.04, 49.78]
        alexi_cs = 0.04
        alexi_x, alexi_y = -125.02, 49.78
    else:
        raise ValueError(f'unsupported ALEXI source: {alexi_coll_id}')
    alexi_crs = 'EPSG:4326'
    logging.debug(f'  Collection: {alexi_coll_id}')


    # Get list of existing images for the target tile
    logging.debug('  Getting GEE asset list')
    asset_coll = ee.ImageCollection(export_coll_id)\
        .filterDate(iter_start_dt.strftime('%Y-%m-%d'), iter_end_dt.strftime('%Y-%m-%d'))
    asset_props = {f'{export_coll_id}/{x["properties"]["system:index"]}': x['properties']
                   for x in utils.get_info(asset_coll)['features']}
    # asset_props = {x['id']: x['properties'] for x in assets_info['features']}


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


    # Process each WRS2 tile separately
    logging.info('\nImage Exports')
    # wrs2_tiles = []
    processed_wrs2_tiles = set()
    for export_info in sorted(export_list, key=lambda i: i['index'], reverse=reverse_flag):
        logging.info(f'{export_info["index"]}')
        # logging.debug('  {} - {}'.format(
        #     export_info['index'], ', '.join(export_info['wrs2_tiles'])))
        tile_count = len(export_info['wrs2_tiles'])
        tile_list = sorted(export_info['wrs2_tiles'], reverse=not(reverse_flag))

        for export_n, wrs2_tile in enumerate(tile_list):
            path, row = map(int, wrs2_tile_re.findall(wrs2_tile)[0])
            if wrs2_tile in processed_wrs2_tiles:
                logging.debug('{} {} ({}/{}) - already processed'.format(
                    export_info['index'], wrs2_tile, export_n + 1, tile_count))
                continue
            processed_wrs2_tiles.add(wrs2_tile)

            if wrs2_skip_list and wrs2_tile in wrs2_skip_list:
                logging.debug('{} {} ({}/{}) - in wrs2 skip list'.format(
                    export_info['index'], wrs2_tile, export_n + 1, tile_count
                ))
                continue
            elif wrs2_row_skip_list and row in wrs2_row_skip_list:
                logging.debug('{} {} ({}/{}) - in wrs2 row skip list'.format(
                    export_info['index'], wrs2_tile, export_n + 1, tile_count
                ))
                continue
            elif wrs2_path_skip_list and path in wrs2_path_skip_list:
                logging.debug('{} {} ({}/{}) - in wrs2 path skip list'.format(
                    export_info['index'], wrs2_tile, export_n + 1, tile_count
                ))
                continue
            else:
                logging.debug('{} {} ({}/{})'.format(
                    export_info['index'], wrs2_tile, export_n + 1, tile_count
                ))

            # path, row = map(int, wrs2_tile_re.findall(export_info['index'])[0])
            # logging.info('WRS2 tile: {}  ({}/{})'.format(
            #     export_info['index'], export_n + 1, len(export_list)))

            # Adjust and snap WRS2 extent to the ALEXI grid
            wrs2_coll_id = 'projects/openet/featureCollections/wrs2/custom'
            wrs2_geom = (
                ee.FeatureCollection(wrs2_coll_id)
                .filterMetadata('wrs2_tile', 'equals', wrs2_tile)
                .first().geometry()
            )
            image_xy = utils.get_info(wrs2_geom.bounds(1, 'EPSG:4326').coordinates().get(0))
            if image_xy is None:
                logging.info('  Could not get image extent, skipping')
                continue
            extent = [
                min([xy[0] for xy in image_xy]),
                min([xy[1] for xy in image_xy]),
                max([xy[0] for xy in image_xy]),
                max([xy[1] for xy in image_xy])
            ]
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


            scene_id = f'{wrs2_tile}'
            export_id = export_id_fmt.format(
                model=ini['INPUTS']['et_model'].lower(),
                index=wrs2_tile,
            )
            export_id = export_id.replace('-', '')
            export_id += export_id_name
            asset_id = f'{export_coll_id}/{scene_id.lower()}'

            if overwrite_flag:
                if export_id in tasks.keys():
                    logging.info(f'  {scene_id} - Task already submitted, cancelling')
                    ee.data.cancelTask(tasks[export_id]['id'])
                    # ee.data.cancelOperation(tasks[export_id]['id'])
                # This is intentionally not an "elif" so that a task can be
                # cancelled and an existing image/file/asset can be removed
                if asset_props and asset_id in asset_props.keys():
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
                elif asset_props and asset_id in asset_props.keys():
                    logging.debug(f'  {scene_id} - Asset already exists, skipping')
                    continue

            logging.debug(f'  Export ID:  {export_id}')
            logging.debug(f'  Collection: {os.path.dirname(asset_id)}')
            logging.debug(f'  Image ID:   {os.path.basename(asset_id)}')

            export_image_list = []
            for month_i, month in enumerate(range(1, 13)):
                ta_coll = (
                    ee.ImageCollection(tair_args['source_coll'])
                    .filterDate(iter_start_dt.strftime('%Y-%m-%d'),
                                iter_end_dt.strftime('%Y-%m-%d'))
                    .filter(ee.Filter.calendarRange(month, month, 'month'))
                    .filterMetadata('wrs2_tile', 'equals', wrs2_tile)
                )

                def ta_mosaic_interpolate(ta_mosaic_img):
                    # Reverse the band order so that we can find the last transition
                    #   from decreasing to increasing with a positive bias
                    ta_bands = ta_mosaic_img.select('step_\\d+_ta').bandNames().reverse()
                    bias_bands = ta_mosaic_img.select('step_\\d+_bias').bandNames().reverse()
                    ta_array = ta_mosaic_img.select(ta_bands).toArray()
                    bias_array = ta_mosaic_img.select(bias_bands).toArray()

                    # Assign the bias that are very similar a very large value so that they will not be selected
                    diff_array = bias_array.arraySlice(0, 1).subtract(bias_array.arraySlice(0, 0, -1))
                    adj_bias_mask = diff_array.abs().lt(0.001)
                    # Repeat the last value to make the array the same length.array is reversed order.
                    adj_bias_mask = adj_bias_mask.arrayCat(adj_bias_mask.arraySlice(0, -1), 0)
                    adj_bias_array = bias_array.add(adj_bias_mask.multiply(99))

                    # Identify the "last" transition from a negative to positive bias
                    # CGM - Having problems with.ceil() limiting result to the image data range
                    # Multiplying by a big number seemed to fix the issue but this could still
                    #   be a problem with the bias ranges get really small
                    sign_array = bias_array.multiply(1000).ceil().max(0).min(1).int()
                    transition_array = sign_array.arraySlice(0, 0, -1).subtract(sign_array.arraySlice(0, 1))
                    # Insert an extra value at the beginning (of reverse, so actually at end)
                    # of the transition array so the indexing lines up for all steps
                    transition_array = bias_array.arraySlice(0, 0, 1).multiply(0).arrayCat(transition_array, 0)
                    transition_index = transition_array.arrayArgmax().arrayFlatten([['index']])
                    # Get the max transition value in order to know if there was a transition
                    transition_max = transition_array.arrayReduce(ee.Reducer.max(), [0]).arrayFlatten([['max']])

                    # Identify the position of minimum absolute bias
                    min_bias_index = adj_bias_array.abs().multiply(-1).arrayArgmax().arrayFlatten([['index']])

                    # Identify the "bracketing" Ta and bias values
                    # If there is a transition, use the "last" transition
                    # If there is not a transition, use the minimum absolute bias for both
                    # Note, the index is for the reversed arrays
                    # B is the "high" value, A is the "low value"
                    index_b = (
                        transition_index.subtract(1).max(0)
                        .where(transition_max.eq(0), min_bias_index)
                    )
                    index_a = (
                        transition_index.min(ta_bands.size().subtract(1))
                        .where(transition_max.eq(0), min_bias_index)
                    )
                    ta_b = ta_array.arrayGet(index_b)
                    ta_a = ta_array.arrayGet(index_a)
                    bias_b = bias_array.arrayGet(index_b)
                    bias_a = bias_array.arrayGet(index_a)

                    # Linearly interpolate Ta
                    # Limit the interpolated value to the bracketing values (don't extrapolate)
                    ta_img = (
                        ta_b.subtract(ta_a).divide(bias_b.subtract(bias_a))
                        .multiply(bias_a.multiply(-1)).add(ta_a)
                        .max(ta_a).min(ta_b)
                    )

                    # # Compute the target Ta as the average of the bracketing Ta values
                    # # instead of interpolating
                    # ta_img = ta_a.add(ta_b).multiply(0.5)

                    # # CGM - This is mostly needed at the 10k step size
                    # # Commenting out for now
                    # # Mask out Ta cells with all negative biases
                    # ta_img = ta_img.updateMask(bias_b.lt(0).And(bias_a.lt(0)).Not())

                    # Round to the nearest tenth (should it be hundredth?)
                    return (
                        ta_img.multiply(10).round().divide(10)
                        .rename(['ta_interp'])
                        # .addBands([ta_a, bias_a, ta_b, bias_b])
                        # .rename(['ta_interp', 'ta_a', 'bias_a', 'ta_b', 'bias_b'])
                        .set({'system:time_start': ta_mosaic_img.get('system:time_start')})
                    )

                # Should the bands be named by month abbreviation or number?
                # Should we return the min and max or other stats?
                ta_interp_img = ta_coll.map(ta_mosaic_interpolate).select('ta_interp') \
                    .mean().rename(calendar.month_abbr[month])
                export_image_list.append(ta_interp_img)

            export_img = ee.Image(export_image_list)

            # TODO: Test if the retile is still needed
            if tair_args['retile'] and tair_args['retile'] > 1:
                export_img = export_img.retile(tair_args['retile'])

            properties = {
                'core_version': openet.core.__version__,
                'date_ingested': datetime.today().strftime('%Y-%m-%d'),
                # 'model_name': model_name,
                # 'model_version': openet.disalexi.__version__,
                'start_date': start_date,
                'end_date': end_date,
                # 'start_year': start_dt.year,
                # 'end_year': end_dt.year,
                'tool_name': TOOL_NAME,
                'tool_version': TOOL_VERSION,
                'wrs2_tile': wrs2_tile,
                # 'system:time_start': image_date.millis(),
            }
            properties.update(model_args)
            properties.update(tair_args)
            export_img = export_img.set(properties)

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
            except:
                logging.warning(f'  {scene_id} - Export task was not started, skipping')
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
        wrs2_property='wrs2'
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
    if (study_area_property == 'STUSPS' and
            'CONUS' in [x.upper() for x in study_area_features]):
        # Exclude AK, HI, AS, GU, PR, MP, VI, (but keep DC)
        study_area_features = [
            'AL', 'AR', 'AZ', 'CA', 'CO', 'CT', 'DC', 'DE', 'FL', 'GA',
            'IA', 'ID', 'IL', 'IN', 'KS', 'KY', 'LA', 'MA', 'MD', 'ME',
            'MI', 'MN', 'MO', 'MS', 'MT', 'NC', 'ND', 'NE', 'NH', 'NJ',
            'NM', 'NV', 'NY', 'OH', 'OK', 'OR', 'PA', 'RI', 'SC', 'SD',
            'TN', 'TX', 'UT', 'VA', 'VT', 'WA', 'WI', 'WV', 'WY',
        ]
    # elif (study_area_property == 'STUSPS' and
    #         'WESTERN11' in [x.upper() for x in study_area_features]):
    #     study_area_features = [
    #         'AZ', 'CA', 'CO', 'ID', 'MT', 'NM', 'NV', 'OR', 'UT', 'WA', 'WY']
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
        tiles_coll = tiles_coll.filter(
            ee.Filter.inList(mgrs_property, mgrs_skip_list).Not()
        )
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
            'wrs2_tiles': sorted(utils.wrs2_str_2_set(
                tile_ftr['properties'][wrs2_property])),
        })

    # Apply the user defined WRS2 tile list
    if wrs2_tiles:
        logging.debug(f'  Filter WRS2 tiles: {wrs2_tiles}')
        for tile in tiles_list:
            tile['wrs2_tiles'] = sorted(list(set(tile['wrs2_tiles']) & set(wrs2_tiles)))

    # Only return export tiles that have intersecting WRS2 tiles
    export_list = [
        tile for tile in sorted(tiles_list, key=lambda k: k['index'])
        if tile['wrs2_tiles']
    ]

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
        description='Compute/export WRS2 Ta climatology images',
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
        '--project', default=None,
        help='Google cloud project ID to use for GEE authentication')
    parser.add_argument(
        '--ready', default=-1, type=int,
        help='Maximum number of queued READY tasks')
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
        start_dt=args.start,
        end_dt=args.end,
    )
