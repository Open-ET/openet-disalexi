import argparse
from calendar import monthrange
import datetime
import logging
import os
import re
import time

from dateutil.relativedelta import relativedelta
import ee
from flask import abort, Response
import google.auth
from google.auth.transport.requests import AuthorizedSession

# import openet.core.utils as utils

SCOPES = [
    'https://www.googleapis.com/auth/earthengine',
    # 'https://www.googleapis.com/auth/devstorage.full_control'
]
CREDENTIALS, project_id = google.auth.default(default_scopes=SCOPES)

ASSET_COLL_ID = 'projects/earthengine-legacy/assets/' \
                'projects/disalexi/insol_data/global_v001_daily_conus'
ASSET_DT_FMT = '%Y%m%d'
SOURCE_COLL_ID = 'projects/earthengine-legacy/assets/' \
                 'projects/disalexi/insol_data/GLOBAL_V001'
RESAMPLE_METHOD = 'bicubic'
UTC_OFFSET = 6
START_DAY_OFFSET = 60
END_DAY_OFFSET = 10


def cfsr_daily_ingest(tgt_dt, overwrite_flag=False,
                      user_credentials_flag=False):
    """

    Parameters
    ----------
    tgt_dt : datetime
    overwrite_flag : bool, optional
    user_credentials_flag : bool, optional
        If True, the GEE key argument will not be set and the export tools will
        attempt to use the user's credentials.
        If False, the tool will attempt to use the project credentials
        (the default is False).

    Returns
    -------
    str : response string

    """

    # tgt_date = tgt_dt.strftime('%Y%m%d')
    logging.info(f'CFSR Daily Insolation - {tgt_dt.strftime("%Y-%m-%d")}')
    # response = f'CFSR Daily Insolation - {tgt_dt.strftime("%Y-%m")}'

    asset_id = f'{ASSET_COLL_ID}/{tgt_dt.strftime(ASSET_DT_FMT)}'
    export_name = f'cfsr_daily_insolation_{tgt_dt.strftime("%Y%m%d")}'

    logging.debug(f'  {SOURCE_COLL_ID}')
    logging.debug(f'  {asset_id}')
    logging.debug(f'  {export_name}')

    # TODO: Move to config.py
    # TODO: Pull credentials automatically when deployed
    logging.debug('\nInitializing Earth Engine')
    if user_credentials_flag:
        logging.debug('  Using user GEE credentials')
        ee.Initialize()
    else:
        ee.Initialize(CREDENTIALS)
    # elif os.path.isfile(GEE_KEY_FILE):
    #     logging.debug(f'  Using service account key file: {GEE_KEY_FILE}')
    #     ee.Initialize(ee.ServiceAccountCredentials('', key_file=GEE_KEY_FILE))
    # else:
    #     raise Exception('EE not initialized')

    if ee.data.getInfo(asset_id):
        if overwrite_flag:
            try:
                ee.data.deleteAsset(asset_id)
            except Exception as e:
                return f'{export_name} - An error occured while trying to '\
                       f'delete the existing asset, skipping\n{e}\n'
        else:
            return f'{export_name} - The asset already exists and overwrite '\
                   f'is False, skipping\n'

    # Project to ALEXI projection/shape/domain
    asset_crs = 'EPSG:4326'
    asset_transform = [0.04, 0, -125.04, 0, -0.04, 49.8]
    asset_shape = '1456x625'

    input_img_id = f'{SOURCE_COLL_ID}/{tgt_dt.strftime("%Y%j")}'
    input_img = ee.Image(input_img_id)
    logging.debug(f'  {input_img_id}')

    if UTC_OFFSET > 0:
        logging.debug('  Positive UTC offset, adding bands for next day')
        next_id = (tgt_dt + datetime.timedelta(days=1)).strftime('%Y%j')
        next_img_id = f'{SOURCE_COLL_ID}/{next_id}'
        logging.debug(f'  {next_img_id}')
        input_img = input_img\
            .select(list(range(UTC_OFFSET, 24)))\
            .addBands(ee.Image(next_img_id).select(list(range(UTC_OFFSET))))
    elif UTC_OFFSET == 0:
        pass
    else:
        raise ValueError('Negative utc_offsets are not supported')

    # TODO: Decide if the scheduler should be responsible for checking if there
    #   are enough source images
    # TODO: Wrap getInfo call in a try/except loop
    try:
        source_count = input_img.bandNames().size().getInfo()
    except:
        source_count = 0
    if source_count == 0:
        return f'{export_name} - source image does not exist\n'
    elif source_count < 24:
        return f'{export_name} - too few source images ({source_count}) for day\n'

    # Sum the hourly bands to daily
    output_img = input_img.reduce(ee.Reducer.sum()) \
        .rename(['rs']).toInt16()

    if RESAMPLE_METHOD != 'nearest':
        output_img = output_img.resample('bicubic')
        #     .reproject(crs=export_crs, crsTransform=export_transform)

    output_img = output_img.set({
        'system:time_start': millis(tgt_dt),
        'date_ingested': input_img.get('DATE_INGESTED'),
        'doy': tgt_dt.strftime('%j'),
        'insolation_version': input_img.get('INSOLATION_VERSION'),
        'resample_method': RESAMPLE_METHOD,
        'units': 'W m-2',
        'utc_offset': UTC_OFFSET,
    })

    task = ee.batch.Export.image.toAsset(
        image=output_img,
        description=export_name,
        assetId=asset_id,
        crs=asset_crs,
        crsTransform=asset_transform,
        dimensions=asset_shape,
    )

    # Start the export task
    task.start()

    # # Try to start the task a couple of times
    # for i in range(1, 6):
    #     try:
    #         task.start()
    #         break
    #     except ee.ee_exception.EEException as e:
    #         logging.warning('EE Exception, retry {}\n{}'.format(i, e))
    #     except Exception as e:
    #         logging.warning('Unhandled Exception: {}'.format(e))
    #         return 'Unhandled Exception: {}'.format(e)
    #     time.sleep(i ** 2)

    return f'{export_name} - {task.id}\n'
    # return f'{export_name}\n'


def cron_scheduler(request):
    """Responds to any HTTP request.

    Parameters
    ----------
    request (flask.Request): HTTP request object.

    Returns
    -------
    The response text or any set of values that can be turned into a
    Response object using
    `make_response <http://flask.pocoo.org/docs/1.0/api/#flask.Flask.make_response>`.

    """
    response = 'Generate Daily CFSR Insolation\n'

    request_json = request.get_json(silent=True)
    request_args = request.args

    # Default start and end date to None if not set
    if request_json and 'start' in request_json:
        start_date = request_json['start']
    elif request_args and 'start' in request_args:
        start_date = request_args['start']
    else:
        start_date = None

    if request_json and 'end' in request_json:
        end_date = request_json['end']
    elif request_args and 'end' in request_args:
        end_date = request_args['end']
    else:
        end_date = None

    if not start_date and not end_date:
        today = datetime.datetime.today()
        start_dt = (datetime.datetime(today.year, today.month, today.day) -
                    relativedelta(days=START_DAY_OFFSET))
        end_dt = (datetime.datetime(today.year, today.month, today.day) -
                  relativedelta(days=END_DAY_OFFSET))
    elif start_date and end_date:
        # Only process custom range if start and end are both set
        # Limit the end date to the last full month date
        try:
            start_dt = datetime.datetime.strptime(start_date, '%Y-%m-%d')
            end_dt = datetime.datetime.strptime(end_date, '%Y-%m-%d')
        except ValueError as e:
            response = 'Error parsing start and/or end date\n'
            response += str(e)
            abort(404, description=response)

        # Force end date to be last day of previous month
        # end_dt = min(end_dt,
        #              datetime.datetime.today() - datetime.timedelta(days=1))

        # TODO: Force start date to be at least one month before end
        # start_dt = min(
        #     start_dt,
        #     end_dt - relativedelta(months=1) + relativedelta(days=1))

        if start_dt > end_dt:
            abort(404, description='Start date must be before end date')
        # elif (end_dt - start_dt) > datetime.timedelta(days=400):
        #     abort(404, description='No more than 1 year can be processed in a single request')
        # if start_dt < datetime.datetime(1980, 1, 1):
        #     logging.debug('Start Date: {} - no CIMIS images before '
        #                   '1980-01-01'.format(start_dt.strftime('%Y-%m-%d')))
        #     start_dt = datetime.datetime(1980, 1, 1)
    else:
        abort(404, description='Both start and end date must be specified')
    response += 'Start Date: {}\n'.format(start_dt.strftime('%Y-%m-%d'))
    response += 'End Date:   {}\n'.format(end_dt.strftime('%Y-%m-%d'))

    args = {
        'start_dt': start_dt,
        'end_dt': end_dt,
    }

    for tgt_dt in cfsr_daily_dates(**args):
        logging.info(f'Date: {tgt_dt.strftime("%Y-%m-%d")}')
        # response += 'Date: {}\n'.format(tgt_dt.strftime('%Y-%m-%d'))
        response += cfsr_daily_ingest(tgt_dt, overwrite_flag=True)

    return Response(response, mimetype='text/plain')


def cfsr_daily_dates(start_dt, end_dt, overwrite_flag=False,
                     user_credentials_flag=False):
    """"""
    logging.debug('\nBuilding CFSR daily date list')

    logging.debug('{}'.format(start_dt.strftime('%Y-%m-%d')))
    logging.debug('{}'.format(end_dt.strftime('%Y-%m-%d')))

    # TODO: Move to config.py
    logging.debug('\nInitializing Earth Engine')
    if user_credentials_flag:
        logging.debug('  Using user GEE credentials')
        ee.Initialize()
    else:
        ee.Initialize(CREDENTIALS)
    # elif os.path.isfile(GEE_KEY_FILE):
    #     logging.debug(f'  Using service account key file: {GEE_KEY_FILE}')
    #     ee.Initialize(ee.ServiceAccountCredentials('', key_file=GEE_KEY_FILE))
    # else:
    #     raise Exception('EE not initialized')

    task_id_re = re.compile('cfsr_daily_insolation_(?P<date>\d{8})')
    # asset_id_re = re.compile(
    #     ASSET_COLL_ID.split('projects/')[-1] + '/(?P<date>\d{8})$')

    # Figure out which asset dates need to be ingested
    # Start with a list of dates to check
    # logging.debug('\nBuilding Date List')
    test_dt_list = list(date_range(start_dt, end_dt, skip_leap_days=False))
    if not test_dt_list:
        logging.info('Empty date range')
        return []
    # logging.info('\nTest dates: {}'.format(
    #     ', '.join(map(lambda x: x.strftime('%Y-%m-%d'), test_dt_list))))

    # Check if any of the needed dates are currently being ingested
    # Check task list before checking asset list in case a task switches
    #   from running to done before the asset list is retrieved.
    task_id_list = [
        desc.replace('\nAsset ingestion: ', '')
        for desc in get_ee_tasks(states=['RUNNING', 'READY']).keys()]
    task_date_list = [
        datetime.datetime.strptime(m.group('date'), '%Y%m%d').strftime('%Y-%m-%d')
        for task_id in task_id_list
        for m in [task_id_re.search(task_id)] if m]
    # logging.info('Task dates: {}'.format(', '.join(task_date_list)))

    # Switch date list to be dates that are missing
    test_dt_list = [
        dt for dt in test_dt_list
        if overwrite_flag or dt.strftime('%Y-%m-%d') not in task_date_list]
    if not test_dt_list:
        logging.info('All dates are queued for export')
        return []
    # else:
    #     logging.info('\nMissing asset dates: {}'.format(', '.join(
    #         map(lambda x: x.strftime('%Y-%m-%d'), test_dt_list))))

    # Check if the assets already exist
    # For now, assume the collection exists
    # Bump end date for filterDate() calls
    logging.debug('\nChecking existing assets')
    filter_end_dt = end_dt + datetime.timedelta(days=1)
    asset_date_coll = ee.ImageCollection(ASSET_COLL_ID) \
            .filterDate(start_dt.strftime('%Y-%m-%d'),
                        filter_end_dt.strftime('%Y-%m-%d'))
    asset_date_list = asset_date_coll \
        .aggregate_array('system:index').getInfo()
    # asset_id_list = get_ee_assets(
    #     ASSET_COLL_ID, start_dt, end_dt + datetime.timedelta(days=1))
    # asset_date_list = [
    #     datetime.datetime.strptime(m.group('date'), ASSET_DT_FMT)
    #         .strftime('%Y-%m-%d')
    #     for asset_id in asset_id_list
    #     for m in [asset_id_re.search(asset_id)] if m]
    logging.debug(f'\nAsset dates: {", ".join(asset_date_list)}')

    # Switch date list to be dates that are missing
    test_dt_list = [
        dt for dt in test_dt_list
        if overwrite_flag or dt.strftime(ASSET_DT_FMT) not in asset_date_list]
    if not test_dt_list:
        logging.info('No dates to process after filtering existing assets')
        return []
    logging.debug('\nDates (after filtering existing assets): {}'.format(
        ', '.join(map(lambda x: x.strftime('%Y-%m-%d'), test_dt_list))))

    # TODO: Should the source collection be checked here to see if there are
    #   enough images?

    return test_dt_list


def date_range(start_dt, end_dt, days=1, skip_leap_days=False):
    """Generate dates within a range (inclusive)

    Parameters
    ----------
    start_dt : datetime
        Start date.
    end_dt : datetime
        End date.
    days : int, optional
        Step size (the default is 1).
    skip_leap_days : bool, optional
        If True, skip leap days while incrementing (the default is True).

    Yields
    ------
    datetime

    """
    import copy
    curr_dt = copy.copy(start_dt)
    while curr_dt <= end_dt:
        if not skip_leap_days or curr_dt.month != 2 or curr_dt.day != 29:
            yield curr_dt
        curr_dt += datetime.timedelta(days=days)


def get_ee_tasks(states=['RUNNING', 'READY'], verbose=False, retries=6):
    """Return current active tasks

    Parameters
    ----------
    states : list, optional
        List of task states to check (the default is ['RUNNING', 'READY']).
    verbose : bool, optional
        This parameter is deprecated and is no longer being used.
        To get verbose logging of the active tasks use utils.print_ee_tasks().
    retries : int, optional
        The number of times to retry getting the task list if there is an error.

    Returns
    -------
    dict : task descriptions (key) and full task info dictionary (value)

    """
    logging.debug('\nRequesting Task List')
    task_list = None
    for i in range(retries):
        try:
            # TODO: getTaskList() is deprecated, switch to listOperations()
            task_list = ee.data.getTaskList()
            # task_list = ee.data.listOperations()
            break
        except Exception as e:
            logging.warning(
                f'  Error getting task list, retrying ({i}/{retries})\n  {e}')
            time.sleep((i+1) ** 2)
    if task_list is None:
        raise Exception('\nUnable to retrieve task list, exiting')

    task_list = sorted(
        [task for task in task_list if task['state'] in states],
        key=lambda t: (t['state'], t['description'], t['id']))
    # task_list = sorted([
    #     [t['state'], t['description'], t['id']] for t in task_list
    #     if t['state'] in states])

    # Convert the task list to a dictionary with the task name as the key
    return {task['description']: task for task in task_list}


# TODO: Pull from openet.core.utils
def millis(input_dt):
    """Convert datetime to milliseconds since epoch

    Parameters
    ----------
    input_dt : datetime

    Returns
    -------
    int

    """
    import calendar
    return 1000 * int(calendar.timegm(input_dt.timetuple()))


# TODO: Pull from openet.core.utils
def arg_valid_date(input_date):
    """Check that a date string is ISO format (YYYY-MM-DD)

    This function is used to check the format of dates entered as command
      line arguments.
    DEADBEEF - It would probably make more sense to have this function
      parse the date using dateutil parser (http://labix.org/python-dateutil)
      and return the ISO format string

    Parameters
    ----------
    input_date : string

    Returns
    -------
    datetime

    Raises
    ------
    ArgParse ArgumentTypeError

    """
    try:
        return datetime.datetime.strptime(input_date, "%Y-%m-%d")
    except ValueError:
        msg = "Not a valid date: '{}'.".format(input_date)
        raise argparse.ArgumentTypeError(msg)


def arg_parse():
    """"""
    today = datetime.date.today()

    parser = argparse.ArgumentParser(
        description='Generate daily CFSR insolation assets',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '--start', type=arg_valid_date, metavar='DATE',
        default=(datetime.datetime(today.year, today.month, today.day) -
                 relativedelta(days=START_DAY_OFFSET)).strftime('%Y-%m-%d'),
        help='Start date (format YYYY-MM-DD)')
    parser.add_argument(
        '--end', type=arg_valid_date, metavar='DATE',
        default=(datetime.datetime(today.year, today.month, today.day) -
                 relativedelta(days=END_DAY_OFFSET)).strftime('%Y-%m-%d'),
        help='End date (format YYYY-MM-DD)')
    parser.add_argument(
        '--overwrite', default=False, action='store_true',
        help='Force overwrite of existing files')
    parser.add_argument(
        '--reverse', default=False, action='store_true',
        help='Process dates in reverse order')
    parser.add_argument(
        '--user_credentials', default=False, action='store_true',
        help='Use the user\'s credentials (instead of the default service '
             'account key file)')
    parser.add_argument(
        '--debug', default=logging.INFO, const=logging.DEBUG,
        help='Debug level logging', action='store_const', dest='loglevel')
    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = arg_parse()
    logging.basicConfig(level=args.loglevel, format='%(message)s')

    # Build the image collection if it doesn't exist
    logging.debug('Image Collection: {}'.format(ASSET_COLL_ID))
    ee.Initialize()
    if not ee.data.getInfo(ASSET_COLL_ID):
        logging.info('\nImage collection does not exist and will be built'
                     '\n  {}'.format(ASSET_COLL_ID))
        input('Press ENTER to continue')
        ee.data.createAsset({'type': 'IMAGE_COLLECTION'}, ASSET_COLL_ID)

    ingest_dt_list = cfsr_daily_dates(
        args.start, args.end, overwrite_flag=args.overwrite,
        user_credentials_flag=args.user_credentials)

    for ingest_dt in sorted(ingest_dt_list, reverse=args.reverse):
        # logging.info(f'Date: {ingest_dt.strftime("%Y-%m-%d")}')
        response = cfsr_daily_ingest(
            ingest_dt, overwrite_flag=args.overwrite,
            user_credentials_flag=args.user_credentials)
        logging.info(f'  {response}')
