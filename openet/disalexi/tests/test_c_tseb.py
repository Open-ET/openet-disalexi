import logging

import ee
import pytest

import openet.disalexi.tseb as tseb
import openet.disalexi.utils as utils

# AmeriFlux sites adjusted to nearest Landsat cell centroid
ne1_xy = [-96.47672812080845, 41.16506126041818]
ne2_xy = [-96.46994024736414, 41.16491226772292]
ne3_xy = [-96.43968912903934, 41.17964494123755]


@pytest.mark.parametrize(
    't_air, t_rad, t_air0, '
    'e_air, u, p, z, '
    'rs_1, rs24, vza,'
    'aleafv, aleafn, aleafl, adeadv, adeadn, adeadl,'
    'albedo, ndvi, lai, clump,'
    'hc_min, hc_max, leaf_width, datetime, lon, lat, a_pt_in, '
    'stabil_iter, albedo_iter, expected, tol',
    [
        # High NDVI site in LC08_044033_20170716
        [300.0, 306.5, 306.3861390180871,
         0.84104, 3.2665367230039, 101.264543111959, 3.0,
         946.69066527778, 8603.212890625, 0,
         0.83, 0.35, 0.95, 0.49, 0.13, 0.95,
         0.1259961, 0.8743930074457752, 4.234, 0.83,
         0.0, 0.6, 0.05, 1500230731090, -121.5265, 38.7399, 1.32,
         25, 10, 6.713934963976832, 0.1],
        # Low NDVI site in LC08_044033_20170716
        [300.0, 323.6, 306.3840757488451,
         0.84104, 3.2665367230039, 101.252726383124, 4.0,
         946.69066527778, 8603.212890625, 0,
         0.83, 0.35, 0.95, 0.49, 0.13, 0.95,
         0.1716302, 0.16195230171935662, 0.5791, 0.83,
         0.0, 0.6, 0.05, 1500230731090, -121.50822, 38.71776, 1.32,
         25, 10, 2.2609140660778926, 0.1],
        # High NDVI site in LC08_044033_20170716
        [298.4269785619036, 306.5, 306.3861390180871,
         0.84104, 3.2665367230039, 101.264543111959, 3.0,
         946.69066527778, 8603.212890625, 0,
         0.83, 0.35, 0.95, 0.49, 0.13, 0.95,
         0.1259961, 0.8743930074457752, 4.234, 0.83,
         0.0, 0.6, 0.05, 1500230731090, -121.5265, 38.7399, 1.32,
         25, 10, 5.84718054547899, 0.1],
        # Low NDVI site in LC08_044033_20170716
        [300.36076433814867, 323.6, 306.3840757488451,
         0.84104, 3.2665367230039, 101.252726383124, 4.0,
         946.69066527778, 8603.212890625, 0,
         0.83, 0.35, 0.95, 0.49, 0.13, 0.95,
         0.1716302, 0.16195230171935662, 0.5791, 0.83,
         0.0, 0.6, 0.05, 1500230731090, -121.50822, 38.71776, 1.32,
         25, 10, 2.3879099490407283, 0.1],
    ]
)
def test_tseb_pt(
        t_air, t_rad, t_air0, e_air, u, p, z, rs_1, rs24, vza,
        aleafv, aleafn, aleafl, adeadv, adeadn, adeadl,
        albedo, ndvi, lai, clump,
        hc_min, hc_max, leaf_width, datetime, lon, lat, a_pt_in,
        stabil_iter, albedo_iter, expected, tol,
        ):

    # study_area = ee.Geometry.Rectangle(-122.00, 38.60, -121.00, 39.0)
    # mask = ee.Image('LANDSAT/LC08/C02/T1_L2/LC08_044033_20170716')\
    #     .select(['SR_B3']).float().multiply(0)\
    #     .clip(study_area)

    output_image = tseb.tseb_pt(
        t_air=ee.Image.constant(t_air), t_rad=ee.Image.constant(t_rad),
        t_air0=ee.Image.constant(t_air0),
        # t_air0=mask.add(t_air0),
        e_air=ee.Image.constant(e_air), u=ee.Image.constant(u),
        p=ee.Image.constant(p), z=ee.Image.constant(z),
        rs_1=ee.Image.constant(rs_1), rs24=ee.Image.constant(rs24),
        vza=ee.Image.constant(vza),
        aleafv=ee.Image.constant(aleafv), aleafn=ee.Image.constant(aleafn),
        aleafl=ee.Image.constant(aleafl), adeadv=ee.Image.constant(adeadv),
        adeadn=ee.Image.constant(adeadn), adeadl=ee.Image.constant(adeadl),
        albedo=ee.Image.constant(albedo), ndvi=ee.Image.constant(ndvi),
        lai=ee.Image.constant(lai),
        # lai=mask.add(lai),
        clump=ee.Image.constant(clump), leaf_width=ee.Image.constant(leaf_width),
        hc_min=ee.Image.constant(hc_min), hc_max=ee.Image.constant(hc_max),
        lon=ee.Image.constant(lon), lat=ee.Image.constant(lat),
        datetime=ee.Date(datetime), a_pt_in=ee.Image.constant(a_pt_in),
        stabil_iter=stabil_iter, albedo_iter=albedo_iter,
    )

    output = list(utils.constant_image_value(output_image).values())[0]
    # output = list(utils.point_image_value(output_image, [lon, lat], 30).values())[0]
    logging.debug('\n  Target values: {}'.format(expected))
    logging.debug('  Output values: {}'.format(output))
    assert abs(output - expected) <= tol
