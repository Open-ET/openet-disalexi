#--------------------------------
# Name:         et_export_wrs2_image.py
# Purpose:      Compute/Export WRS2 ET images
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
from osgeo import ogr, osr

import openet.disalexi as disalexi
import openet.core.utils as utils


def main(ini_path=None, overwrite_flag=False, delay=0, key=None,
         random_flag=False):
    """Compute WRS2 ET images

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
    random_flag : bool, optional
        If True, process dates and tiles in random order (the default is False).

    """
    logging.info('\nCompute WRS2 ET images')

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
    # tair_args = {}
    # for k, v in dict(ini['TAIR']).items():
    #     if utils.is_number(v):
    #         if v.isdigit():
    #             tair_args[k.lower()] = int(v)
    #         else:
    #             tair_args[k.lower()] = float(v)
    #     else:
    #         tair_args[k.lower()] = v
    # # tair_args = {
    # #     k.lower(): int(v) if utils.is_number(v) else v
    # #     for k, v in dict(ini['TAIR']).items()}

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

    try:
        wrs2_tiles = str(ini['INPUTS']['wrs2_tiles'])
        wrs2_tiles = sorted([x.strip() for x in wrs2_tiles.split(',')])
        wrs2_tiles = [x.replace('p', '').replace('r', '') for x in wrs2_tiles]
    except KeyError:
        logging.info('\nINPUTS wrs2_tiles parameter was net set, default to []')
        wrs2_tiles = []
    except Exception as e:
        raise e

    logging.info('\nDISALEXI Parameters')
    logging.info('  Stabil iter: {}'.format(int(model_args['stabil_iterations'])))
    logging.info('  Albedo iter: {}'.format(int(model_args['albedo_iterations'])))

    # Output Ta daily image collection
    et_wrs2_coll_id = '{}'.format(ini['EXPORT']['export_coll'])

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


    if 'study_area_path' in ini['INPUTS'].keys():
        logging.info('\nReading study area shapefile')
        logging.info('  {}'.format(ini['INPUTS']['study_area_path']))
        study_area_ds = ogr.Open(ini['INPUTS']['study_area_path'], 0)
        study_area_lyr = study_area_ds.GetLayer()
        study_area_osr = study_area_lyr.GetSpatialRef()
        study_area_crs = str(study_area_osr.ExportToWkt())
        # study_area_proj4 = study_area_osr.ExportToProj4()
        logging.debug('  Study area projection: {}'.format(study_area_crs))

        # Get the dissolved/unioned geometry of the study area
        output_geom = ogr.Geometry(ogr.wkbMultiPolygon)
        for study_area_ftr in study_area_lyr:
            study_area_geom = output_geom.Union(study_area_ftr.GetGeometryRef())
        study_area_ds = None

        simplify_buffer = 0
        if simplify_buffer:
            study_area_geom = output_geom\
                .SimplifyPreserveTopology(simplify_buffer)\
                .buffer(simplify_buffer)
        else:
            # Added flatten call to change clockwise geometies to counter cw
            study_area_geom.FlattenTo2D()

        export_geom = ee.Geometry(
            json.loads(study_area_geom.ExportToJson()), study_area_crs, False)

    else:
        export_geom = alexi_mask.geometry()


    # Get current asset list
    logging.debug('\nGetting asset list')
    # DEADBEEF - daily is hardcoded in the asset_id for now
    asset_list = utils.get_ee_assets(et_wrs2_coll_id)

    # Get current running tasks
    tasks = utils.get_ee_tasks()
    if logging.getLogger().getEffectiveLevel() == logging.DEBUG:
        logging.debug('  Tasks: {}\n'.format(len(tasks)))
        input('ENTER')

    # Limit by year and month
    try:
        month_list = sorted(list(utils.parse_int_set(ini['INPUTS']['months'])))
    except:
        logging.debug('\nINPUTS "months" parameter not set in the INI,'
                      '\n  Defaulting to all months (1-12)\n')
        month_list = list(range(1, 13))
    try:
        year_list = sorted(list(utils.parse_int_set(ini['INPUTS']['years'])))
    except:
        logging.debug('\nINPUTS "years" parameter not set in the INI,'
                      '\n  Defaulting to all available years\n')
        year_list = []
    # month_list = list(range(1, 13))
    # year_list = []

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
    if random_flag:
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

        # Build and merge the Landsat collections
        # Eventually build this using the model Collection class
        # Time filters are to remove bad (L5) and pre-op (L8) images
        landsat_coll = ee.ImageCollection([])
        if ('LANDSAT/LC08/C01/T1_SR' in collections and
                export_date > '2013-03-24'):
            l8_coll = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR')\
                .filterDate(export_dt, export_dt + datetime.timedelta(days=1))\
                .filterBounds(export_geom)\
                .filter(ee.Filter.gt('system:time_start',
                                     ee.Date('2013-03-24').millis()))\
                .filterMetadata('CLOUD_COVER_LAND', 'less_than',
                                float(ini['INPUTS']['cloud_cover']))
            landsat_coll = ee.ImageCollection(landsat_coll.merge(l8_coll))
        if ('LANDSAT/LE07/C01/T1_SR' in collections and
                export_date >= '1999-01-01'):
            l7_coll = ee.ImageCollection('LANDSAT/LE07/C01/T1_SR')\
                .filterDate(export_dt, export_dt + datetime.timedelta(days=1))\
                .filterBounds(export_geom)\
                .filterMetadata('CLOUD_COVER_LAND', 'less_than',
                                float(ini['INPUTS']['cloud_cover']))
            landsat_coll = ee.ImageCollection(landsat_coll.merge(l7_coll))
        if ('LANDSAT/LT05/C01/T1_SR' in collections and
                export_date <= '2011-12-31'):
            l5_coll = ee.ImageCollection('LANDSAT/LT05/C01/T1_SR')\
                .filterDate(export_dt, export_dt + datetime.timedelta(days=1))\
                .filterBounds(export_geom)\
                .filter(ee.Filter.lt('system:time_start',
                                     ee.Date('2011-12-31').millis()))\
                .filterMetadata('CLOUD_COVER_LAND', 'less_than',
                                float(ini['INPUTS']['cloud_cover']))
            landsat_coll = ee.ImageCollection(landsat_coll.merge(l5_coll))
        if ('LANDSAT/LT04/C01/T1_SR' in collections and
                export_date <= '1993-12-01'):
            l4_coll = ee.ImageCollection('LANDSAT/LT04/C01/T1_SR')\
                .filterDate(export_dt, export_dt + datetime.timedelta(days=1))\
                .filterBounds(export_geom)\
                .filterMetadata('CLOUD_COVER_LAND', 'less_than',
                                float(ini['INPUTS']['cloud_cover']))
            landsat_coll = ee.ImageCollection(landsat_coll.merge(l4_coll))

        # DATA_TYPE filter only needed if using TOA collections
        #     .filterMetadata('DATA_TYPE', 'equals', 'L1TP')

        image_id_list = sorted(list(set(
            landsat_coll.aggregate_array('system:id').getInfo())))
        if not image_id_list:
            logging.debug('  Empty image ID list, skipping')
            continue
        # if random_flag:
        #     random.shuffle(image_id_list)

        for image_id in image_id_list:
            scene_id = image_id.split('/')[-1]
            landsat, path, row, year, month, day = parse_landsat_id(scene_id)
            image_dt = datetime.datetime.strptime(
                '{:04d}{:02d}{:02d}'.format(year, month, day), '%Y%m%d')
            # image_date = image_dt.date().isoformat()
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
            logging.debug('  Export ID: {}'.format(export_id))

            asset_id = '{}/{}'.format(et_wrs2_coll_id, scene_id)
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
                    logging.info('  Task already submitted, skipping')
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
                # 'cloud_cover_max': float(ini['INPUTS']['cloud_cover']),
            }
            properties.update(model_args)
            # properties.update(tair_args)

            # if tair_args['iteration'] <= 0:
            landsat_img = ee.Image(image_id)
            d_obj = disalexi.Image(disalexi.LandsatSR(landsat_img).prep(),
                                   **model_args)
            et_a = d_obj.et(ta_img=ee.Image.constant(250)) \
                .select(['et'], ['et_a']).float()
            # et_d = d_obj.et(ta_img=ee.Image.constant(350)) \
            #     .select(['et'], ['et_d']).float()

            # export_img = ee.Image([et_a, et_d])\
            export_img = et_a\
                .set({
                    'id': landsat_img.get('system:id'),
                    'system:index': landsat_img.get('system:index'),
                    'system:time_start': landsat_img.get('system:time_start'),
                    'spacecraft_id': landsat_img.get('SATELLITE'),
                    'cloud_cover_land': landsat_img.get('CLOUD_COVER_LAND'),
                })\
                .set(properties)

            # DEADBEEF - Trying to see if retile will help
            # if tair_args['retile'] and tair_args['retile'] > 0:
            #     export_img = export_img.retile(tair_args['retile'])
            export_img = export_img.retile(8)

            landsat_info = landsat_img.select([0]).getInfo()
            export_crs = landsat_info['bands'][0]['crs']
            export_geo = landsat_info['bands'][0]['crs_transform']
            export_shape = landsat_info['bands'][0]['dimensions']
            logging.debug('  CRS: {}'.format(export_crs))
            logging.debug('  Geo: {}'.format(export_geo))
            logging.debug('  Shape: {}'.format(export_shape))

            # Build export tasks
            logging.debug('  Building export task')
            task = ee.batch.Export.image.toAsset(
                image=export_img,
                description=export_id,
                assetId=asset_id,
                crs=export_crs,
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
        description='Compute/export WRS2 ET images',
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
        '--random', default=False, action='store_true',
        help='Process dates and tiles in random order')
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
         key=args.key, random_flag=args.random)
