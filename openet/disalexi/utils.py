import logging
import time

import ee

# Use US-NE1 as default sample point
default_xy = (-96.47672812080845, 41.16506126041818)
default_crs = 'EPSG:32614'
default_geo = [30, 0, 632685, 0, -30, 4742715]
default_scale = 0.1


def constant_image_value(image, crs=default_crs):
    value = image.reduceRegion(
        reducer=ee.Reducer.first(),
        geometry=ee.Geometry.Rectangle([0, 0, 10, 10], crs, False),
        scale=1)

    return ee_getinfo(value)
    # return value.getInfo()


def image_value(image, xy=default_xy, scale=None,
                crs=default_crs, crsTransform=default_geo, tile_scale=1):
    # Default to using Landsat crsTransform unless scale is set

    if scale is not None:
        value = image.reduceRegion(
            reducer=ee.Reducer.first(), geometry=ee.Geometry.Point(xy),
            scale=default_scale, tileScale=tile_scale)
    else:
        value = image.reduceRegion(
            reducer=ee.Reducer.first(), geometry=ee.Geometry.Point(xy),
            crs=crs, crsTransform=crsTransform, tileScale=tile_scale)

    return ee_getinfo(value)
    # return value.getInfo()


def coll_value(coll, xy=default_xy, crs=default_crs,
               crsTransform=default_geo):
    # values = coll.getRegion(
    #     geometry=ee.Geometry.Point(xy),
    #     crs=crs, crsTransform=crsTransform).getInfo()
    values = coll.getRegion(
        geometry=ee.Geometry.Point(xy),
        scale=default_scale).getInfo()
    return [item for item in values]


def ee_getinfo(ee_obj, n=4):
    """Make an exponential backoff getInfo call on the Earth Engine object

    getInfo() does this internally also, so the purpose is primarily to
    attempt to catch other exceptions and timeout issues.

    """
    output = None
    for i in range(1, n):
        try:
            output = ee_obj.getInfo()
        except Exception as e:
            logging.info('    Resending query ({}/10)'.format(i))
            logging.debug('    {}'.format(e))
            time.sleep(i ** 2)
        if output:
            break
    return output
