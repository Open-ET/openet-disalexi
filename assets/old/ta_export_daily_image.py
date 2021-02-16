#--------------------------------
# Name:         ta_export_daily_image.py
# Purpose:      Compute/Export daily Ta images
#--------------------------------

import argparse
from builtins import input
import datetime
import logging
import math
import os
import pprint
import sys

import ee

import openet.disalexi as disalexi
# from . import utils
import utils


def main(ini_path=None, overwrite_flag=False, delay=0, key=None):
    """Compute daily Ta images

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
    logging.info('\nCompute daily Ta images')

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
    tair_args = {
        k.lower(): int(v) if utils.is_number(v) else v
        for k, v in dict(ini['TAIR']).items()}

    if 'ta_retile' not in tair_args.keys():
        tair_args['ta_retile'] = 0
    ta_values = list(range(
        tair_args['ta_start'], tair_args['ta_stop'], tair_args['ta_step']))
    # if 'ta_iterations' not in tair_args.keys():
    #     logging.info('  Iterations: {}'.format(tair_args['ta_iterations']))

    logging.info('  Start: {}'.format(tair_args['ta_start']))
    logging.info('  Stop:  {}'.format(tair_args['ta_stop']))
    logging.info('  Step:  {}'.format(tair_args['ta_step']))
    logging.info('  Cellsize: {}'.format(tair_args['ta_cellsize']))
    logging.info('  Retile:   {}'.format(tair_args['ta_retile']))

    # Output Ta daily image collection
    ta_daily_coll_id = '{}'.format(ini['EXPORT']['export_coll'])

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

    logging.debug('\nExport properties')
    export_geo = ee.Image(alexi_mask).projection().getInfo()['transform']
    export_shape = ee.Image(alexi_mask).getInfo()['bands'][0]['dimensions']
    export_extent = [
        export_geo[2], export_geo[5] + export_shape[1] * export_geo[4],
        export_geo[2] + export_shape[0] * export_geo[0], export_geo[5]]
    logging.debug('  CRS: {}'.format(alexi_crs))
    logging.debug('  Extent: {}'.format(export_extent))
    logging.debug('  Geo: {}'.format(export_geo))
    logging.debug('  Shape: {}'.format(export_shape))

    # # Limit export to a user defined study area or geometry?
    # export_geom = ee.Geometry.Rectangle(
    #     [-125, 24, -65, 50], proj='EPSG:4326', geodesic=False)  # CONUS
    # export_geom = ee.Geometry.Rectangle(
    #     [-124, 35, -119, 42], proj='EPSG:4326', geodesic=False)  # California

    # If cell_size parameter is set in the INI,
    # adjust the output` cellsize and recompute the transform and shape
    try:
        export_cs = float(ini['EXPORT']['cell_size'])
        export_shape = [
            int(math.ceil(abs((export_shape[0] * export_geo[0]) / export_cs))),
            int(math.ceil(abs((export_shape[1] * export_geo[4]) / export_cs)))]
        export_geo = [export_cs, 0.0, export_geo[2], 0.0, -export_cs, export_geo[5]]
        logging.debug('  Custom export cell size: {}'.format(export_cs))
        logging.debug('  Geo: {}'.format(export_geo))
        logging.debug('  Shape: {}'.format(export_shape))
    except KeyError:
        pass

    # Get current asset list
    logging.debug('\nGetting asset list')
    # DEADBEEF - daily is hardcoded in the asset_id for now
    asset_list = utils.get_ee_assets(ta_daily_coll_id)

    # Get current running tasks
    tasks = utils.get_ee_tasks()
    if logging.getLogger().getEffectiveLevel() == logging.DEBUG:
        utils.print_ee_tasks()
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

    # Iterate over date ranges
    for export_dt in utils.date_range(iter_start_dt, iter_end_dt):
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

        export_id = ini['EXPORT']['export_id_fmt'] \
            .format(date=export_dt.strftime('%Y%m%d'))
        # export_id = '{}_coarse_{}_{}_{}_{}'.format(
        export_id = '{}_qm_{}_{}_{}_{}_{}'.format(
            export_id,
            ini['DISALEXI']['stabil_iterations'],
            ini['DISALEXI']['albedo_iterations'],
            tair_args['ta_step'],
            tair_args['ta_retile'],
            tair_args['ta_cellsize'],
        )
        logging.debug('  Export ID: {}'.format(export_id))

        asset_id = '{}/{}_qm_{}_{}_{}_{}_{}'.format(
            ta_daily_coll_id, export_dt.strftime('%Y%m%d'),
            ini['DISALEXI']['stabil_iterations'],
            ini['DISALEXI']['albedo_iterations'],
            tair_args['ta_step'],
            tair_args['ta_retile'],
            tair_args['ta_cellsize'],
        )
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

        # Build and merge the Landsat collections
        # Time filters are to remove bad (L5) and pre-op (L8) images
        #     .filterBounds(export_geom) \
        l8_coll = ee.ImageCollection('LANDSAT/LC08/C01/T1_RT_TOA') \
            .filterDate(export_dt, export_dt + datetime.timedelta(days=1)) \
            .filterBounds(alexi_mask.geometry()) \
            .filterMetadata('CLOUD_COVER_LAND', 'less_than',
                            float(ini['INPUTS']['cloud_cover'])) \
            .filterMetadata('DATA_TYPE', 'equals', 'L1TP') \
            .filter(ee.Filter.gt('system:time_start',
                                 ee.Date('2013-03-24').millis()))
        l7_coll = ee.ImageCollection('LANDSAT/LE07/C01/T1_RT_TOA') \
            .filterDate(export_dt, export_dt + datetime.timedelta(days=1)) \
            .filterBounds(alexi_mask.geometry()) \
            .filterMetadata('CLOUD_COVER_LAND', 'less_than',
                            float(ini['INPUTS']['cloud_cover'])) \
            .filterMetadata('DATA_TYPE', 'equals', 'L1TP')
        l5_coll = ee.ImageCollection('LANDSAT/LT05/C01/T1_TOA') \
            .filterDate(export_dt, export_dt + datetime.timedelta(days=1)) \
            .filterBounds(alexi_mask.geometry()) \
            .filterMetadata('CLOUD_COVER_LAND', 'less_than',
                            float(ini['INPUTS']['cloud_cover'])) \
            .filterMetadata('DATA_TYPE', 'equals', 'L1TP') \
            .filter(ee.Filter.lt('system:time_start',
                                 ee.Date('2011-12-31').millis()))
        # l4_coll = ee.ImageCollection('LANDSAT/LT04/C01/T1_TOA') \
        #     .filterDate(export_dt, export_dt + datetime.timedelta(days=1)) \
        #     .filterBounds(alexi_mask.geometry()) \
        #     .filterMetadata('CLOUD_COVER_LAND', 'less_than',
        #                     float(ini['INPUTS']['cloud_cover'])) \
        #     .filterMetadata('DATA_TYPE', 'equals', 'L1TP')

        # if export_date <= '1993-12-31':
        #     landsat_coll = ee.ImageCollection(l5_coll.merge(l4_coll))
        if export_date < '1999-01-01':
            landsat_coll = l5_coll
        elif export_date <= '2011-12-31':
            landsat_coll = ee.ImageCollection(l7_coll.merge(l5_coll))
        elif export_date <= '2013-03-24':
            landsat_coll = l7_coll
        else:
            landsat_coll = ee.ImageCollection(l8_coll.merge(l7_coll))



        # # # image_id_list = landsat_coll.aggregate_array('system:id').getInfo()
        # # wrs2_tile_list = sorted(list(set([
        # #     x.split('/')[-1][5:11]
        # #     for x in landsat_coll.aggregate_array('system:id').getInfo()])))
        # # landsat_list = sorted(list(set(
        # #     landsat_coll.aggregate_array('SPACECRAFT_ID').getInfo())))
        # # print(wrs2_tile_list)
        # # print(landsat_list)
        # # input('ENTER')
        #
        # def landsat_prep_func(image):
        #     return disalexi.LandsatTOA(image).prep()
        #
        # input_img = ee.ImageCollection(landsat_coll.map(landsat_prep_func))\
        #     .mean()\
        #     .set({'system:time_start': ee.Date(export_date).millis(),
        #           'system:index': export_date})
        #
        # properties = {
        #     'system:time_start': utils.millis(export_dt),
        #     'date_ingested': datetime.datetime.today().strftime('%Y-%m-%d'),
        #     # 'landsat': ','.join(landsat_list),
        #     # 'wrs2_tiles': ','.join(wrs2_tile_list),
        #     'date': export_dt.strftime('%Y-%m-%d'),
        #     'year': int(export_dt.year),
        #     'month': int(export_dt.month),
        #     'day': int(export_dt.day),
        #     'doy': int(export_dt.strftime('%j')),
        #     'cycle_day': ((export_dt - cycle_base_dt).days % 8) + 1,
        #     'model_name': 'DISALEXI',
        #     'model_version': disalexi.__version__,
        #     'cloud_cover_max': float(ini['INPUTS']['cloud_cover']),
        #     # 'collections': ', '.join(collections),
        # }
        # properties.update(model_args)
        # properties.update(tair_args)
        #
        # d_obj = disalexi.Image(input_img, **model_args)
        # ta = d_obj.ta_qm(cellsize=tair_args['ta_cellsize'],
        #                  ta_values=ta_values)
        # # ta = d_obj.ta_iter(cellsize=tair_args['ta_cellsize'],
        # #                    ta_init=[tair_args['ta_start'],
        # #                             tair_args['ta_stop']],
        # #                    iterations=tair_args['ta_iterations'],
        # #                    method='golden')
        # export_img = ta.select(['ta', 'bias']) \
        #     .float().set(properties)
        #     # .clip(ee.Image(image).geometry())
        #
        # if tair_args['ta_retile'] and tair_args['ta_retile']  > 0:
        #     export_img.retile(tair_args['ta_retile'])
        #
        # print(export_img.getInfo())
        #
        # # Build export tasks
        # logging.debug('  Building export task')
        # task = ee.batch.Export.image.toAsset(
        #     image=export_img,
        #     description=export_id,
        #     assetId=asset_id,
        #     crs=export_crs,
        #     crsTransform='[' + ','.join(list(map(str, export_geo))) + ']',
        #     dimensions='{0}x{1}'.format(*export_shape),
        # )
        # logging.info('  Starting export task')
        # utils.ee_task_start(task)
        #
        # # Pause before starting next task
        # utils.delay_task(delay)
        # logging.debug('')



        def ta_img_func(image):
            d_obj = disalexi.Image(
                disalexi.LandsatTOA(image).prep(), **model_args)

            # Remove the merged collection indices from the system:index
            scene_id = ee.List(
                ee.String(image.get('system:index')).split('_')).slice(-3)
            scene_id = ee.String(scene_id.get(0)).cat('_') \
                .cat(ee.String(scene_id.get(1))).cat('_') \
                .cat(ee.String(scene_id.get(2)))

            ta = d_obj.ta_qm(cellsize=tair_args['ta_cellsize'],
                             ta_values=ta_values)
            # ta = d_obj.ta_iter(cellsize=tair_args['ta_cellsize'],
            #                    ta_init=[tair_args['ta_start'],
            #                             tair_args['ta_stop']],
            #                    iterations=tair_args['ta_iterations'],
            #                    method='golden')

            return alexi_mask.add(ta)\
                .select(['ta', 'bias'])\
                .set({
                    'system:time_start': image.get('system:time_start'),
                    'wrs2_tile': scene_id.slice(5, 11),
                    'spacecraft_id': image.get('SPACECRAFT_ID'),
                })
                # .clip(ee.Image(image).geometry())\

        ta_img_coll = ee.ImageCollection(landsat_coll.map(ta_img_func)).mean()

        def unique_properties(coll, property):
            return ee.String(ee.List(ee.Dictionary(
                coll.aggregate_histogram(property)).keys()).join(','))
        wrs2_tile_list = ee.String('').cat(unique_properties(
            ta_img_coll, 'wrs2_tile'))
        landsat_list = ee.String('').cat(unique_properties(
            ta_img_coll, 'spacecraft_id'))

        properties = {
            'system:time_start': utils.millis(export_dt),
            'date_ingested': datetime.datetime.today().strftime('%Y-%m-%d'),
            'landsat': landsat_list,
            'wrs2_tiles': wrs2_tile_list,
            'date': export_dt.strftime('%Y-%m-%d'),
            'year': int(export_dt.year),
            'month': int(export_dt.month),
            'day': int(export_dt.day),
            'doy': int(export_dt.strftime('%j')),
            'cycle_day': ((export_dt - cycle_base_dt).days % 8) + 1,
            'model_name': 'DISALEXI',
            'model_version': disalexi.__version__,
            'cloud_cover_max': float(ini['INPUTS']['cloud_cover']),
            # 'collections': ', '.join(collections),
        }
        properties.update(model_args)
        properties.update(tair_args)

        # Cast to float and set properties
        ta_img = ta_img.rename(['ta', 'bias']).float().set(properties)
        # ta_img = ta_img.rename(['ta']).float().set(properties)

        if tair_args['ta_retile'] and tair_args['ta_retile']  > 0:
            ta_img.retile(tair_args['ta_retile'])

        # Build export tasks
        logging.debug('  Building export task')
        task = ee.batch.Export.image.toAsset(
            image=ta_img,
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


def arg_parse():
    """"""
    parser = argparse.ArgumentParser(
        description='Compute/export daily Ta images',
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
