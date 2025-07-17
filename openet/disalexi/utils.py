import calendar
from datetime import datetime, UTC
import logging
import pprint
import time

import ee

# Use US-NE1 as default sample point
default_xy = (-96.47672812080845, 41.16506126041818)
default_crs = 'EPSG:32614'
default_geo = [30, 0, 632685, 0, -30, 4742715]
default_scale = 0.1


def getinfo(ee_obj, n=4):
    """Make an exponential back off getInfo call on an Earth Engine object"""
    output = None
    for i in range(1, n):
        try:
            output = ee_obj.getInfo()
        except ee.ee_exception.EEException as e:
            if 'Earth Engine memory capacity exceeded' in str(e):
                logging.info(f'    Resending query ({i}/{n})')
                logging.debug(f'    {e}')
                time.sleep(i ** 2)
            else:
                raise e

        if output:
            break

    return output


# TODO: Import from common.utils
# Should these be test fixtures instead?
# I'm not sure how to make them fixtures and allow input parameters
def constant_image_value(image, crs='EPSG:32613', scale=1):
    """Extract the output value from a calculation done with constant images"""
    return getinfo(ee.Image(image).reduceRegion(
        reducer=ee.Reducer.first(), scale=scale,
        geometry=ee.Geometry.Rectangle([0, 0, 10, 10], crs, False))
    )


def point_image_value(image, xy, scale=1):
    """Extract the output value from a calculation at a point"""
    return getinfo(ee.Image(image).reduceRegion(
        reducer=ee.Reducer.first(), geometry=ee.Geometry.Point(xy), scale=scale)
    )


def point_coll_value(coll, xy, scale=1):
    """Extract the output value from a calculation at a point"""
    output = getinfo(coll.getRegion(ee.Geometry.Point(xy), scale=scale))

    # Structure output to easily be converted to a Pandas dataframe
    # First key is band name, second key is the date string
    col_dict = {}
    info_dict = {}
    for i, k in enumerate(output[0][4:]):
        col_dict[k] = i + 4
        info_dict[k] = {}
    for row in output[1:]:
        date = datetime.fromtimestamp(row[3] / 1000.0, UTC).strftime('%Y-%m-%d')
        for k, v in col_dict.items():
            info_dict[k][date] = row[col_dict[k]]

    return info_dict


def boolean(x):
    """Convert boolean like objects to boolean"""
    if type(x) is str:
        if x.upper() in ['TRUE', 'T']:
            return True
        elif x.upper() in ['FALSE', 'F']:
            return False
        else:
            raise ValueError('"{}" could not be interpreted as bool'.format(x))
    elif type(x) is bool:
        return x
    else:
        raise ValueError('"{}" could not be interpreted as bool'.format(x))


def date_to_time_0utc(date):
    """Get the 0 UTC time_start for a date

    Parameters
    ----------
    date : ee.Date

    Returns
    -------
    ee.Number

    """
    return ee.Date.fromYMD(date.get('year'), date.get('month'), date.get('day')).millis()


def interpolate(coll, interp_dt, timestep=3, offset=0):
    """Temporally interpolate an image collection

    Parameters
    ----------
    coll : ee.ImageCollection
        Single band image collection.
    interp_dt : ee.Date
        The datetime to interpolate to.
    timestep : int
        Image collection time step in hours.
    offset : float
        Not currently implemented.  This parameter would be useful if the values
        in the image collection where the average over timestep instead of the
        instantaneous value at the time start.

    Returns
    -------
    ee.Image

    Notes
    -----
    This particular implementation does not work correctly if the interpolation
    time is exactly the time start of one of the images.

    """
    # interp_dt = interp_dt.advance(offset, 'hour')
    a_img = ee.Image(coll.filterDate(interp_dt.advance(-timestep, 'hour'), interp_dt).first())
    b_img = ee.Image(coll.filterDate(interp_dt, interp_dt.advance(timestep, 'hour')).first())
    a_time = ee.Number(a_img.get('system:time_start'))
    b_time = ee.Number(b_img.get('system:time_start'))
    return (
        b_img.subtract(a_img)
        .multiply(interp_dt.millis().subtract(a_time).divide(b_time.subtract(a_time)))
        .add(a_img)
    )


def is_number(x):
    try:
        float(x)
        return True
    except:
        return False


def millis(input_dt):
    """Convert datetime to milliseconds since epoch

    Parameters
    ----------
    input_dt : datetime

    Returns
    -------
    int

    """
    return 1000 * int(calendar.timegm(input_dt.timetuple()))


def valid_date(date_str, date_fmt='%Y-%m-%d'):
    """Check if a datetime can be built from a date string and if it is valid"""
    try:
        datetime.strptime(date_str, date_fmt)
        return True
    except Exception as e:
        return False


# def point_image_value(image, xy, scale=None, transform=None, crs=None,):
#     """Extract the output value from a calculation at a point"""
#     rr_args = {
#         'reducer': ee.Reducer.first(),
#         'geometry': ee.Geometry.Point(xy),
#         'crs': 'EPSG:4326',
#         'scale': 1,
#     }
#     if transform:
#         rr_args.update({'crsTransform': transform})
#         del rr_args['scale']
#     if scale:
#         rr_args.update({'scale': scale})
#         del rr_args['crsTransform']
#     if crs:
#         rr_args.update({'crs': crs})
#
#     return ee.Image(image).select([0], ['value']).reduceRegion(**rr_args).get('value').getInfo()


# def point_image_values(image, band_name):
#     """Extract the output value from a calculation at a point"""
#     # print(f'\n{band_name}')
#     points = [
#         # Bad XY
#         [-106.60, 39.72],
#         # # Good XYs
#         # [-106.56, 39.72],
#         # [-106.60, 39.76],
#         # [-106.64, 39.72],
#         # [-106.60, 39.68],
#     ]
#     for xy in points:
#         print(f'{band_name}: {point_image_value(image, xy, transform=[0.04, 0, -125.02, 0, -0.04, 49.78])}')
