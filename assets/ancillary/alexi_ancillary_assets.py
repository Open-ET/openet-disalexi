import argparse
import logging
import pprint

import ee

logging.getLogger('earthengine-api').setLevel(logging.INFO)
#logging.getLogger('googleapiclient').setLevel(logging.INFO)
#logging.getLogger('requests').setLevel(logging.INFO)
#logging.getLogger('urllib3').setLevel(logging.INFO)


def main(project_id, overwrite_flag=False):
    """

    """

    ee.Initialize(project=project_id)

    output_folder = 'projects/openet/assets/alexi/ancillary/conus/v006'

    elevation_source = 'USGS/SRTMGL1_003'
    elevation_img = ee.Image(elevation_source)

    # Hardcoding export parameters for now instead of reading from the image
    alexi_shape = [1456, 625]
    alexi_geo = [0.04, 0, -125.02, 0, -0.04, 49.78]
    alexi_crs = 'EPSG:4326'
    # alexi_source = 'projects/ee-tulipyangyun-2/assets/alexi/ALEXI_V006'
    # alexi_img = ee.ImageCollection(alexi_source).first().multiply(0).unmask(0, True)
    # alexi_crs = alexi_img.projection().crs()
    # alexi_geo = ee.List(ee.Dictionary(ee.Algorithms.Describe(alexi_img.projection())).get('transform'))

    elev_coarse = (
        elevation_img.multiply(1)
        .reduceResolution(reducer=ee.Reducer.mean().unweighted(), maxPixels=30000, bestEffort=False)
        .reproject(crs=alexi_crs, crsTransform=alexi_geo)
    )
    lat_coarse = elev_coarse.multiply(0).add(ee.Image.pixelLonLat().select('latitude'))
    lon_coarse = elev_coarse.multiply(0).add(ee.Image.pixelLonLat().select('longitude'))

    elev_asset_id = f'{output_folder}/elevation'
    lat_asset_id = f'{output_folder}/latitude'
    lon_asset_id = f'{output_folder}/longitude'

    if overwrite_flag:
        if ee.data.getInfo(elev_asset_id):
            ee.data.deleteAsset(elev_asset_id)
        if ee.data.getInfo(lat_asset_id):
            ee.data.deleteAsset(lat_asset_id)
        if ee.data.getInfo(lon_asset_id):
            ee.data.deleteAsset(lon_asset_id)

    if not ee.data.getInfo(elev_asset_id):
        logging.info('Starting elevation asset export')
        task = ee.batch.Export.image.toAsset(
            image=elev_coarse.rename('elevation'),
            description='alexi_elevation_asset',
            assetId=f'{output_folder}/elevation',
            crs=alexi_crs,
            crsTransform='[' + ','.join(list(map(str, alexi_geo))) + ']',
            dimensions='{0}x{1}'.format(*alexi_shape),
        )
        task.start()

    if not ee.data.getInfo(lat_asset_id):
        logging.info('Starting latitude asset export')
        task = ee.batch.Export.image.toAsset(
            image=lat_coarse.rename('latitude'),
            description='alexi_latitude_asset',
            assetId=lat_asset_id,
            crs=alexi_crs,
            crsTransform='[' + ','.join(list(map(str, alexi_geo))) + ']',
            dimensions='{0}x{1}'.format(*alexi_shape),
        )
        task.start()

    if not ee.data.getInfo(lon_asset_id):
        logging.info('Starting longitude asset export')
        task = ee.batch.Export.image.toAsset(
            image=lon_coarse.rename('longitude'),
            description='alexi_longitude_asset',
            assetId=lon_asset_id,
            crs=alexi_crs,
            crsTransform='[' + ','.join(list(map(str, alexi_geo))) + ']',
            dimensions='{0}x{1}'.format(*alexi_shape),
        )
        task.start()

    # TODO: Add code to set a public acl


def arg_parse():
    """"""
    parser = argparse.ArgumentParser(
        description='Compute/export WRS2 Ta images',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '--project', default=None,
        help='Google cloud project ID to use for GEE authentication')
    parser.add_argument(
        '--overwrite', default=False, action='store_true',
        help='Force overwrite of existing files')
    parser.add_argument(
        '--debug', default=logging.INFO, const=logging.DEBUG,
        help='Debug level logging', action='store_const', dest='loglevel')
    args = parser.parse_args()

    return args


if __name__ == "__main__":
    args = arg_parse()

    logging.basicConfig(level=args.loglevel, format='%(message)s')

    main(project_id=args.project, overwrite_flag=args.overwrite)
