#--------------------------------
# Name:         ta_export_wrs2_image.py
# Purpose:      Compute/Export WRS2 Ta images
#--------------------------------

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
import sys
import time

import ee
from osgeo import ogr, osr

import openet.disalexi
import openet.sharpen
import utils
# from . import utils


def main(ini_path=None, overwrite_flag=False, delay_time=0, gee_key_file=None,
         random_flag=False, max_ready=-1, reverse_flag=False):
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
    random_flag : bool, optional
        If True, process dates and tiles in random order (the default is False).
    max_ready: int, optional
        Maximum number of queued "READY" tasks.  The default is -1 which is
        implies no limit to the number of tasks that will be submitted.
    reverse_flag : bool, optional
        If True, process WRS2 tiles in reverse order (the default is False).

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
    if gee_key_file:
        logging.info('  Using service account key file: {}'.format(gee_key_file))
        # The "EE_ACCOUNT" parameter is not used if the key file is valid
        ee.Initialize(ee.ServiceAccountCredentials(
            'deadbeef', key_file=gee_key_file))
    else:
        ee.Initialize()

    if not ee.data.getInfo(ta_wrs2_coll_id.rsplit('/', 1)[0]):
        logging.debug('\nFolder does not exist and will be built'
                      '\n  {}'.format(ta_wrs2_coll_id.rsplit('/', 1)[0]))
        input('Press ENTER to continue')
        ee.data.createAsset({'type': 'FOLDER'},
                            ta_wrs2_coll_id.rsplit('/', 1)[0])
    if not ee.data.getInfo(ta_wrs2_coll_id):
        logging.info('\nExport collection does not exist and will be built'
                     '\n  {}'.format(ta_wrs2_coll_id))
        input('Press ENTER to continue')
        ee.data.createAsset({'type': 'IMAGE_COLLECTION'}, ta_wrs2_coll_id)

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
    excl_end_dt = iter_end_dt + datetime.timedelta(days=1)
    excl_end_date = excl_end_dt.strftime('%Y-%m-%d')

    # DEADBEEF - For testing, process the dates in random order
    date_list = list(utils.date_range(iter_start_dt, iter_end_dt))
    if random_flag:
        random.shuffle(date_list)

    logging.info('\nGetting image ID list')
    # Build and merge the Landsat collections
    # Eventually build this using the model Collection class
    # Time filters are to remove bad (L5) and pre-op (L8) images
    landsat_coll = ee.ImageCollection([])
    if ('LANDSAT/LC08/C01/T1_SR' in collections and
            ini['INPUTS']['end_date'] > '2013-03-24'):
        l8_coll = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR')\
            .filterDate(ini['INPUTS']['start_date'], excl_end_date)\
            .filterBounds(export_geom)\
            .filter(ee.Filter.gt('system:time_start',
                                 ee.Date('2013-03-24').millis()))\
            .filterMetadata('CLOUD_COVER_LAND', 'less_than',
                            float(ini['INPUTS']['cloud_cover']))
        landsat_coll = ee.ImageCollection(landsat_coll.merge(l8_coll))
    if ('LANDSAT/LE07/C01/T1_SR' in collections and
            ini['INPUTS']['end_date'] >= '1999-01-01'):
        l7_coll = ee.ImageCollection('LANDSAT/LE07/C01/T1_SR')\
            .filterDate(ini['INPUTS']['start_date'], excl_end_date)\
            .filterBounds(export_geom)\
            .filterMetadata('CLOUD_COVER_LAND', 'less_than',
                            float(ini['INPUTS']['cloud_cover']))
        landsat_coll = ee.ImageCollection(landsat_coll.merge(l7_coll))
    if ('LANDSAT/LT05/C01/T1_SR' in collections and
            ini['INPUTS']['start_date'] <= '2011-12-31'):
        l5_coll = ee.ImageCollection('LANDSAT/LT05/C01/T1_SR')\
            .filterDate(ini['INPUTS']['start_date'], excl_end_date)\
            .filterBounds(export_geom)\
            .filter(ee.Filter.lt('system:time_start',
                                 ee.Date('2011-12-31').millis()))\
            .filterMetadata('CLOUD_COVER_LAND', 'less_than',
                            float(ini['INPUTS']['cloud_cover']))
        landsat_coll = ee.ImageCollection(landsat_coll.merge(l5_coll))
    if ('LANDSAT/LT04/C01/T1_SR' in collections and
            ini['INPUTS']['start_date'] <= '1993-12-01'):
        l4_coll = ee.ImageCollection('LANDSAT/LT04/C01/T1_SR')\
            .filterDate(ini['INPUTS']['start_date'], excl_end_date)\
            .filterBounds(export_geom)\
            .filterMetadata('CLOUD_COVER_LAND', 'less_than',
                            float(ini['INPUTS']['cloud_cover']))
        landsat_coll = ee.ImageCollection(landsat_coll.merge(l4_coll))
    # DATA_TYPE filter only needed if using TOA collections
    #     .filterMetadata('DATA_TYPE', 'equals', 'L1TP')

    image_id_list = sorted(list(set(
        landsat_coll.aggregate_array('system:id').getInfo())))
    if not image_id_list:
        logging.debug('  Empty image ID list, exiting')
        return False

    for iteration in range(tair_args['min_iteration'], tair_args['max_iteration'] + 1):
        logging.info('\nIteration {}'.format(iteration))

        # Get current asset list
        logging.debug('\nGetting asset list')
        asset_list = utils.get_ee_assets(ta_wrs2_coll_id)

        # Get current running tasks
        tasks = utils.get_ee_tasks()
        if logging.getLogger().getEffectiveLevel() == logging.DEBUG:
            logging.debug('  Tasks: {}\n'.format(len(tasks)))
            input('ENTER')

        # CGM - Rebuild the iteration dictionary each time to catch new images
        ta_coll = ee.ImageCollection(ta_wrs2_coll_id)
        ta_iter_dict = defaultdict(list)
        for img in ta_coll.getInfo()['features']:
            ta_iter_dict[img['properties']['id']].append(
                int(img['properties']['iteration']))
        # pprint.pprint(ta_iter_dict)
        # input('ENTER')

        # Iterate over date ranges
        for export_dt in sorted(date_list, reverse=reverse_flag):
            export_date = export_dt.strftime('%Y-%m-%d')
            if ((month_list and export_dt.month not in month_list) or
                    ( year_list and export_dt.year not in year_list)):
                logging.debug('Date: {} - skipping'.format(export_date))
                continue
            logging.debug('Date: {}'.format(export_date))

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

            # Checking date with underscore just in case there are other numbers
            export_id_list = [
                image_id for image_id in image_id_list
                if export_dt.strftime('_%Y%m%d') in image_id]

            for image_id in export_id_list:
                scene_id = image_id.split('/')[-1]
                landsat, path, row, year, month, day = parse_landsat_id(scene_id)
                # image_dt = datetime.datetime.strptime(
                #     '{:04d}{:02d}{:02d}'.format(year, month, day), '%Y%m%d')
                # image_date = image_dt.date().isoformat()
                # logging.debug('  Date: {}'.format(image_date))
                # logging.debug('  DOY: {}'.format(doy))

                wrs2_tile = '{:03d}{:03d}'.format(path, row)
                if wrs2_tiles and wrs2_tile not in wrs2_tiles:
                    logging.debug('{} - wrs2 not in list, skipping'.format(image_id))
                    continue
                else:
                    logging.info('{}'.format(image_id))

                # Check if there is a previous iteration
                if image_id not in ta_iter_dict.keys():
                    ta_iter_dict[image_id] = []
                elif iteration in ta_iter_dict[image_id] and not overwrite_flag:
                    logging.info('  Iteration {} already computed, '
                                 'skipping'.format(iteration))
                    continue
                elif iteration == 0:
                    pass
                elif iteration - 1 not in ta_iter_dict[image_id]:
                    logging.info('  Missing previous iteration, skipping')
                    continue

                export_id = ini['EXPORT']['export_id_fmt'] \
                    .format(index=scene_id.lower(), iter=iteration)
                logging.debug('  Export ID: {}'.format(export_id))

                asset_id = '{}/{}_{:02d}'.format(
                    ta_wrs2_coll_id, scene_id, iteration)
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
                    'model_name': model_name,
                    'model_version': openet.disalexi.__version__,
                    # 'cloud_cover_max': float(ini['INPUTS']['cloud_cover']),
                    'iteration': int(iteration),
                }
                properties.update(model_args)

                landsat_img = ee.Image(image_id)
                d_obj = openet.disalexi.Image(
                    openet.disalexi.LandsatSR(landsat_img).prep(), **model_args)

                if iteration <= 0:
                    # For initial iteration compute bias at 250 and 350 K
                    a_img = d_obj.ta_coarse(ta_img=ee.Image.constant(250)) \
                        .select(['ta', 'bias'], ['ta_a', 'bias_a'])
                    b_img = d_obj.ta_coarse(ta_img=ee.Image.constant(350)) \
                        .select(['ta', 'bias'], ['ta_b', 'bias_b'])
                    # Applying both bias masks to the output
                    # This shouldn't be necessary but ta_coarse was returning
                    #   ta and bias images with different masks
                    export_img = ee.Image([a_img, b_img])\
                        .updateMask(a_img.select(['bias_a']).And(
                            b_img.select(['bias_b'])))\
                        .double()
                # elif iteration <= 4:
                #     # Interpolate new Ta from the bias and test directly
                #     # Roughly equivalent to false position method
                #     ta_img = ee.Image('{}/{}_{:02d}'.format(
                #         ta_wrs2_coll_id, scene_id, iteration - 1))
                #     # ta_img = ee.Image(ta_coll
                #     #     .filterMetadata('date', 'equals', export_date)
                #     #     .filterMetadata('iteration', 'equals', iteration - 1)
                #     #     .first())
                #     ta_a = ta_img.select(['ta_a'])
                #     ta_b = ta_img.select(['ta_b'])
                #     bias_a = ta_img.select(['bias_a'])
                #     bias_b = ta_img.select(['bias_b'])
                #
                #     ta_x = bias_a.multiply(ta_b).subtract(bias_b.multiply(ta_a))\
                #         .divide(bias_a.subtract(bias_b))
                #     bias_x = d_obj.et_bias(d_obj.et_coarse(ta_x))
                #
                #     # Use the new value if it minimizes the bias and brackets 0
                #     mask1 = bias_x.lt(bias_b).And(bias_x.gt(0))
                #     mask2 = bias_x.gt(bias_a).And(bias_x.lt(0))
                #     ta_b = ta_b.where(mask1, ta_x)
                #     bias_b = bias_b.where(mask1, bias_x)
                #     ta_a = ta_a.where(mask2, ta_x)
                #     bias_a = bias_a.where(mask2, bias_x)
                #     export_img = ee.Image([ta_a, bias_a, ta_b, bias_b])
                else:
                    # Generate test Ta randomly from a triangular distribution
                    # Center the distribution on the interpolated zero bias Ta
                    ta_img = ee.Image('{}/{}_{:02d}'.format(
                        ta_wrs2_coll_id, scene_id, iteration - 1))
                    # ta_img = ee.Image(ta_coll
                    #     .filterMetadata('date', 'equals', export_date)
                    #     .filterMetadata('iteration', 'equals', iteration - 1)
                    #     .first())
                    ta_a = ta_img.select(['ta_a'])
                    ta_b = ta_img.select(['ta_b'])
                    bias_a = ta_img.select(['bias_a'])
                    bias_b = ta_img.select(['bias_b'])

                    ta_c = bias_a.multiply(ta_b).subtract(bias_b.multiply(ta_a))\
                        .divide(bias_a.subtract(bias_b))
                    # ta_c = ta_a.add(ta_b).multiply(0.5)

                    # For now use a single random number for the whole scene
                    # Need to check if ee.Image.random() will work though
                    u = ta_b.multiply(0).add(random.random())
                    # u = ta_b.multiply(0).add(ee.Image.random(0))

                    a = u.multiply(ta_b.subtract(ta_a))\
                        .multiply(ta_c.subtract(ta_a)).sqrt().add(ta_a)
                    b = u.multiply(-1).add(1)\
                        .multiply(ta_b.subtract(ta_a))\
                        .multiply(ta_b.subtract(ta_c))\
                        .sqrt().multiply(-1).add(ta_b)
                    fc = ta_c.subtract(ta_a).divide(ta_b.subtract(ta_a))
                    ta_x = a.where(u.gt(fc), b)
                    bias_x = d_obj.et_bias(d_obj.et_coarse(ta_x))

                    # Use the new value if it minimizes the bias and brackets 0
                    mask1 = bias_x.lt(bias_b).And(bias_x.gt(0))
                    mask2 = bias_x.gt(bias_a).And(bias_x.lt(0))
                    ta_b = ta_b.where(mask1, ta_x)
                    bias_b = bias_b.where(mask1, bias_x)
                    ta_a = ta_a.where(mask2, ta_x)
                    bias_a = bias_a.where(mask2, bias_x)
                    export_img = ee.Image([ta_a, bias_a, ta_b, bias_b])

                # else:
                #     ta_img = ee.Image('{}/{}_{}'.format(
                #         ta_wrs2_coll_id, scene_id, iteration - 1))
                #     # ta_img = ee.Image(ta_coll
                #     #     .filterMetadata('date', 'equals', export_date)
                #     #     .filterMetadata('iteration', 'equals', iteration - 1)
                #     #     .first())
                #     ta_a = ta_img.select(['ta_a'])
                #     ta_b = ta_img.select(['ta_b'])
                #     bias_a = ta_img.select(['bias_a'])
                #     bias_b = ta_img.select(['bias_b'])
                #     # abs_a = ta_img.select(['bias_a']).abs()
                #     # abs_b = ta_img.select(['bias_b']).abs()
                #
                #     # Compute new test Ta and biases
                #     ta_c = ta_b.subtract(ta_b.subtract(ta_a).multiply(0.618034))
                #     ta_d = ta_a.add(ta_b.subtract(ta_a).multiply(0.618034))
                #     bias_c = d_obj.et_bias(d_obj.et_coarse(ta_c))
                #     bias_d = d_obj.et_bias(d_obj.et_coarse(ta_d))
                #     abs_c = bias_c.abs()
                #     abs_d = bias_d.abs()
                #
                #     # Use the new values if they minimize the bias
                #     # If f(c) < f(d): move the data from d to b and c to d
                #     mask1 = abs_c.lt(abs_d)
                #     # If f(c) > f(d): move the data from c to a and d to c
                #     mask2 = abs_c.gte(abs_d)
                #     ta_b = ta_b.where(mask1, ta_d)
                #     bias_b = bias_b.where(mask1, bias_d)
                #     ta_a = ta_a.where(mask2, ta_c)
                #     bias_a = bias_a.where(mask2, bias_c)
                #
                #     export_img = ee.Image([ta_a, bias_a, ta_b, bias_b])

                export_img = export_img\
                    .set({
                        'id': landsat_img.get('system:id'),
                        'system:index': landsat_img.get('system:index'),
                        'system:time_start': landsat_img.get('system:time_start'),
                        'spacecraft_id': landsat_img.get('SATELLITE'),
                        'cloud_cover_land': landsat_img.get('CLOUD_COVER_LAND'),
                    })\
                    .set(properties)

                if tair_args['retile'] and tair_args['retile'] > 0:
                    export_img = export_img.retile(tair_args['retile'])

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
                logging.debug('    {}'.format(task.id))

                # Pause before starting next task
                utils.delay_task(delay)
                logging.debug('')

        # CGM - Add check here to not go on until all tasks have finished
        while True:
            tasks = {k: v for k, v in utils.get_ee_tasks(verbose=False).items()
                     if k.startswith('disalexi_tair_')}
            if len(tasks.keys()) > 0:
                logging.info('Pausing until all tasks are complete')
                time.sleep(60)
            else:
                break


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
        '--random', default=False, action='store_true',
        help='Process dates and tiles in random order')
    parser.add_argument(
        '--ready', default=-1, type=int,
        help='Maximum number of queued READY tasks')
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
    logging.getLogger('googleapiclient').setLevel(logging.ERROR)

    main(ini_path=args.ini, overwrite_flag=args.overwrite,
         delay_time=args.delay, gee_key_file=args.key, random_flag=args.random,
         max_ready=args.ready,
         )
