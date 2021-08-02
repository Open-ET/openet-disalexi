import argparse
import datetime
import logging
import pprint

import ee

import openet.core.utils as utils


def main(overwrite_flag=False, gee_key_file=None, reverse_flag=False,
         recent_days=0, start_dt=None, end_dt=None,
         daily_flag=True, hourly_flag=False):
    """Unpack multiband insol assets

    Parameters
    ----------
    overwrite_flag : bool, optional
        If True, overwrite existing files (the default is False).
    gee_key_file : str, None, optional
        Earth Engine service account JSON key file (the default is None).
    reverse_flag : bool, optional
        If True, process dates in reverse order (the default is False).
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
    daily_flag : bool, optional
        If True, build daily images (the default is True).
    hourly_flag : bool, optional
        If True, build hourly images (the default is False).

    """
    logging.info('\nUnpack multiband insol assets')

    input_coll_id = 'projects/earthengine-legacy/assets/' \
                    'projects/disalexi/insol_data/GLOBAL_V001'
    # input_coll_id = 'projects/disalexi/insol_data/GLOBAL_V001'

    daily_coll_id = 'projects/earthengine-legacy/assets/' \
                    'projects/disalexi/insol_data/global_v001_daily_conus'
    hourly_coll_id = 'projects/earthengine-legacy/assets/' \
                     'projects/disalexi/insol_data/global_v001_hourly_conus'

    resample_method = 'bicubic'
    utc_offset = 6

    # Project to ALEXI projection/shape/domain
    export_crs = 'EPSG:4326'
    export_transform = [0.04, 0, -125.04, 0, -0.04, 49.8]
    export_shape = '1456x625'

    # export_crs = 'EPSG:4326'
    # export_transform = [0.25, 0, -180, 0, -0.25, 90]
    # export_shape = '1440x720'
    # # export_info = ee.Image(ee.ImageCollection(input_coll_id).first())\
    # #     .select([0]).getInfo()
    # # export_crs = export_info['bands][0]['crs']
    # # export_transform = export_info['bands][0]['crs_transform']
    # # export_shape = '{0}x{1}'.format(*export_info['bands][0]['dimensions'])

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
        logging.info('\nUser defined date range')
        start_date = start_dt.strftime('%Y-%m-%d')
        # Make end date exclusive for filterDate call
        end_dt = end_dt + datetime.timedelta(days=1)
        end_date = end_dt.strftime('%Y-%m-%d')
    else:
        logging.info('\nDefault date range')
        start_date = '2001-01-01'
        end_date = '2019-11-01'
        # start_date = '2001-01-01'
        # end_date = '2019-11-01'
    # TODO: Add a few more checks on the dates
    if end_dt < start_dt:
        raise ValueError('end date can not be before start date')
    logging.info('  Start: {}'.format(start_date))
    logging.info('  End:   {}'.format(end_date))


    logging.info('\nInitializing Earth Engine')
    if gee_key_file:
        logging.info('  Using service account key file: {}'.format(gee_key_file))
        # The "EE_ACCOUNT" parameter is not used if the key file is valid
        ee.Initialize(ee.ServiceAccountCredentials(
            'deadbeef', key_file=gee_key_file))
    else:
        ee.Initialize()


    # if not ee.data.getInfo(daily_coll_id.rsplit('/', 1)[0]):
    #     logging.debug('\nFolder does not exist and will be built'
    #                   '\n  {}'.format(daily_coll_id.rsplit('/', 1)[0]))
    #     input('Press ENTER to continue')
    #     ee.data.createAsset({'type': 'FOLDER'}, daily_coll_id.rsplit('/', 1)[0])
    if daily_flag and not ee.data.getInfo(daily_coll_id):
        logging.info('\nExport collection does not exist and will be built'
                     '\n  {}'.format(daily_coll_id))
        input('Press ENTER to continue')
        ee.data.createAsset({'type': 'IMAGE_COLLECTION'}, daily_coll_id)
    if hourly_flag and not ee.data.getInfo(hourly_coll_id):
        logging.info('\nExport collection does not exist and will be built'
                     '\n  {}'.format(hourly_coll_id))
        input('Press ENTER to continue')
        ee.data.createAsset({'type': 'IMAGE_COLLECTION'}, hourly_coll_id)


    # Build the input image ID list
    input_coll = ee.ImageCollection(input_coll_id)
    if start_date and end_date:
        input_coll = input_coll.filterDate(start_date, end_date)
    input_id_list = input_coll.aggregate_array('system:index').getInfo()
    # pprint.pprint(input_id_list)


    # Get current running tasks
    logging.info('\nRequesting Task List')
    tasks = utils.get_ee_tasks()
    if logging.getLogger().getEffectiveLevel() == logging.DEBUG:
        utils.print_ee_tasks(tasks)
        input('ENTER')


    if daily_flag:
        logging.debug('\nGetting daily asset list')
        asset_coll = ee.ImageCollection(daily_coll_id)\
            .filterDate(start_date, end_date)
        asset_props = {
            f'{daily_coll_id}/{x["properties"]["system:index"]}': x['properties']
            for x in utils.get_info(asset_coll)['features']}

        logging.debug('\nExporting daily assets')
        for input_id in sorted(input_id_list, reverse=reverse_flag):
            logging.info(f'{input_id}')
            export_dt = datetime.datetime.strptime(input_id, '%Y%j')
            logging.debug(f'  {export_dt}')

            input_img_id = f'{input_coll_id}/{input_id}'
            input_img = ee.Image(input_img_id)
            logging.debug(f'  {input_img_id}')

            next_id = (export_dt + datetime.timedelta(days=1)).strftime('%Y%j')
            next_img_id = f'{input_coll_id}/{next_id}'
            next_img = ee.Image(next_img_id)
            logging.debug(f'  {next_id}')

            image_id = f'{daily_coll_id}/{export_dt.strftime("%Y%m%d")}'
            logging.debug(f'  {image_id}')

            export_id = image_id \
                .replace('projects/earthengine-legacy/assets/', '') \
                .replace('projects/', '') \
                .replace('/', '_')
            logging.debug(f'  {export_id}')

            if overwrite_flag:
                if export_id in tasks.keys():
                    logging.info('  Task already submitted, cancelling')
                    ee.data.cancelTask(tasks[export_id]['id'])
                    # ee.data.cancelOperation(tasks[export_id]['id'])
                # This is intentionally not an "elif" so that a task can be
                # cancelled and an existing image/file/asset can be removed
                if asset_props and image_id in asset_props.keys():
                    logging.info('  Asset already exists, removing')
                    ee.data.deleteAsset(image_id)
            else:
                if export_id in tasks.keys():
                    logging.info('    Task already submitted, skipping')
                    continue
                elif asset_props and image_id in asset_props.keys():
                    logging.info('  Asset already exists, skipping')
                    continue

            if utc_offset > 0:
                input_img = input_img\
                    .select(list(range(utc_offset, 24)))\
                    .addBands(next_img.select(list(range(utc_offset))))
                # print(input_img.bandNames().getInfo())
                # input('ENTER')
            elif utc_offset == 0:
                pass
            else:
                raise ValueError('Negative utc_offsets are not supported')

            # Sum the hourly bands to daily
            output_img = input_img.reduce(ee.Reducer.sum()) \
                .rename(['rs']).toInt16()

            if resample_method != 'nearest':
                output_img = output_img.resample('bicubic')
                #     .reproject(crs=export_crs, crsTransform=export_transform)

            output_img = output_img.set({
                'system:time_start': utils.millis(export_dt),
                'date_ingested': input_img.get('DATE_INGESTED'),
                'doy': export_dt.strftime('%j'),
                'insolation_version': input_img.get('INSOLATION_VERSION'),
                'resample_method': resample_method,
                'units': 'W m-2',
                'utc_offset': utc_offset,
            })

            task = ee.batch.Export.image.toAsset(
                image=output_img,
                description=export_id,
                assetId=image_id,
                crs=export_crs,
                crsTransform=export_transform,
                dimensions=export_shape,
            )
            # logging.debug(task.status())

            logging.debug('  Starting export task')
            try:
                task.start()
            except Exception as e:
                logging.error(e)
                input('Press ENTER to continue')
                continue

            # logging.debug(task.status())
            logging.debug('')


    if hourly_flag:
        logging.debug('\nGetting hourly asset list')
        asset_coll = ee.ImageCollection(hourly_coll_id)\
            .filterDate(start_date, end_date)
        asset_props = {
            f'{hourly_coll_id}/{x["properties"]["system:index"]}': x['properties']
            for x in utils.get_info(asset_coll)['features']}

        logging.debug('\nExporting hourly assets')
        for input_id in sorted(input_id_list, reverse=reverse_flag):
            logging.info(f'{input_id}')
            input_img_id = f'{input_coll_id}/{input_id}'
            input_img = ee.Image(input_img_id)
            logging.debug(f'  {input_img_id}')

            input_dt = datetime.datetime.strptime(input_id, '%Y%j')
            logging.debug(f'  {input_dt}')

            # Process hourly images
            for hour in range(0, 24):
                logging.info(f'Hour: {hour}')

                export_dt = datetime.datetime.strptime(
                    '{}_{:02d}'.format(input_id, hour), '%Y%j_%H')
                logging.debug(f'  {export_dt}')

                image_id = '{}/{}{:02d}'.format(
                    hourly_coll_id, input_dt.strftime("%Y%m%d"), hour)
                logging.debug(f'  {image_id}')

                export_id = image_id \
                    .replace('projects/earthengine-legacy/assets/', '') \
                    .replace('projects/', '') \
                    .replace('/', '_')

                if overwrite_flag:
                    if export_id in tasks.keys():
                        logging.info('  Task already submitted, cancelling')
                        ee.data.cancelTask(tasks[export_id]['id'])
                        # ee.data.cancelOperation(tasks[export_id]['id'])
                    # This is intentionally not an "elif" so that a task can be
                    # cancelled and an existing image/file/asset can be removed
                    if asset_props and image_id in asset_props.keys():
                        logging.info('  Asset already exists, removing')
                        ee.data.deleteAsset(image_id)
                else:
                    if export_id in tasks.keys():
                        logging.info('  Task already submitted, skipping')
                        continue
                    elif asset_props and image_id in asset_props.keys():
                        logging.info('  Asset already exists, skipping')
                        continue

                output_img = ee.Image(input_img_id).select([hour], ['rs']) \
                    .multiply(1) \
                    .toInt16() \
                    .set({'system:time_start': utils.millis(export_dt),
                          'date_ingested': input_img.get('DATE_INGESTED'),
                          'insolation_version': input_img.get('INSOLATION_VERSION'),
                          'doy': export_dt.strftime('%j'),
                          'units': 'W m-2',
                          })

                task = ee.batch.Export.image.toAsset(
                    image=output_img,
                    description=export_id,
                    assetId=image_id,
                    crs=export_crs,
                    crsTransform=export_transform,
                    dimensions=export_shape,
                )

                logging.debug('  Starting export task')
                try:
                    task.start()
                except Exception as e:
                    logging.error(e)
                    input('Press ENTER to continue')
                    continue

                # logging.debug(task.status())
            logging.debug('')


def arg_parse():
    """"""
    parser = argparse.ArgumentParser(
        description='Unpack multiband insol assets',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '--debug', default=logging.INFO, const=logging.DEBUG,
        help='Debug level logging', action='store_const', dest='loglevel')
    parser.add_argument(
        '--key', type=utils.arg_valid_file, metavar='FILE',
        help='JSON key file')
    parser.add_argument(
        '--overwrite', default=False, action='store_true',
        help='Force overwrite of existing files')
    parser.add_argument(
        '--recent', default=0, type=int,
        help='Number of days to process before current date '
             '(ignore INI start_date and end_date')
    parser.add_argument(
        '--reverse', default=False, action='store_true',
        help='Process dates in reverse order')
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

    main(overwrite_flag=args.overwrite, gee_key_file=args.key,
         reverse_flag=args.reverse, recent_days=args.recent,
         start_dt=args.start, end_dt=args.end,
    )