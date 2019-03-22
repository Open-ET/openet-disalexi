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
    'T_air, T_rad, u,'
    'p, z, Rs_1, '
    'Rs24, vza,'
    'aleafv, aleafn, aleafl, adeadv, adeadn, adeadl,'
    'albedo, ndvi, lai, clump,'
    'hc, leaf_width, datetime, lon, lat, a_PT_in, '
    'stabil_iter, albedo_iter, expected',
    [
        # US-NE1 - IDL iteration 1
        [296.86259968439742, 302.49074951171878, 7.02662301063538,
         97.23062487251868, 350, 917.87845865885413,
         30.81263826599121 / (0.0864 / 24), 0,
         0.83, 0.35, 0.95, 0.49, 0.13, 0.95,
         0.19908118247986, 0.84600001573563, 2.82700014114380, 0.83,
         0.44531310527793, 0.05, 1404839150550, ne1_xy[0], ne1_xy[1], 1.32,
         36, 10,5.80784207568089 / 0.408],
        # US-NE1 - IDL iteration 2
        [298.76438086370422, 302.49074951171878, 7.02662301063538,
         97.23062487251868, 350, 917.87845865885413,
         30.81263826599121 / (0.0864 / 24), 0,
         0.83, 0.35, 0.95, 0.49, 0.13, 0.95,
         0.19908118247986, 0.84600001573563, 2.82700014114380, 0.83,
         0.44531310527793, 0.05, 1404839150550, ne1_xy[0], ne1_xy[1], 1.32,
         36, 10, 6.78850125182576 / 0.408],
        # US-NE1 - IDL iteration 3
        [298.47296436495526, 302.49074951171878, 7.02662301063538,
         97.23062487251868, 350, 917.87845865885413,
         30.81263826599121 / (0.0864 / 24), 0,
         0.83, 0.35, 0.95, 0.49, 0.13, 0.95,
         0.19908118247986, 0.84600001573563, 2.82700014114380, 0.83,
         0.44531310527793, 0.05, 1404839150550, ne1_xy[0], ne1_xy[1], 1.32,
         36, 10, 6.69752524158102 / 0.408],
    ]
)
def test_tseb_pt(T_air, T_rad, u, p, z, Rs_1, Rs24, vza, aleafv, aleafn,
                 aleafl, adeadv, adeadn, adeadl, albedo, ndvi, lai, clump, hc,
                 leaf_width, datetime, lon, lat, a_PT_in,
                 stabil_iter, albedo_iter, expected, tol=1E-8):
    output_image = tseb.tseb_pt(
        ee.Image.constant(T_air), ee.Image.constant(T_rad),
        ee.Image.constant(u), ee.Image.constant(p), ee.Image.constant(z),
        ee.Image.constant(Rs_1), ee.Image.constant(Rs24), ee.Image.constant(vza),
        ee.Image.constant(aleafv), ee.Image.constant(aleafn),
        ee.Image.constant(aleafl), ee.Image.constant(adeadv),
        ee.Image.constant(adeadn), ee.Image.constant(adeadl),
        ee.Image.constant(albedo), ee.Image.constant(ndvi), ee.Image.constant(lai),
        ee.Image.constant(clump), ee.Image.constant(leaf_width),
        ee.Image.constant(hc), ee.Image.constant(hc),
        ee.Date(datetime), ee.Image.constant(lon), ee.Image.constant(lat),
        ee.Image.constant(a_PT_in), stabil_iter, albedo_iter)

    output = list(utils.constant_image_value(output_image).values())[0]
    logging.debug('\n  Target values: {}'.format(expected))
    logging.debug('  Output values: {}'.format(output))
    assert abs(output - expected) <= tol
