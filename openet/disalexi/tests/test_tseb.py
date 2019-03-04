import logging

import ee
import pytest

import openet.disalexi.tseb as tseb
import openet.disalexi.utils as utils


@pytest.mark.parametrize(
    'T_air, T_rad, u,'
    'p, z, Rs_1, Rs24, vza,'
    'zs, aleafv, aleafn, aleafl, adeadv, adeadn, adeadl,'
    'albedo, ndvi, lai, clump,'
    'hc, time, t_rise,'
    't_end, leaf_width, a_PT_in, iterations,expected',
    [
        # US-NE1 - IDL iteration 1
        [296.86259968439742, 302.49074951171878, 7.02662301063538,
         97.23062487251868, 350, 917.87845865885413,
         30.81263826599121 / (0.0864 / 24), 0, 0.45641128977509,
         0.83, 0.35, 0.95, 0.49, 0.13, 0.95,
         0.19908118247986, 0.84600001573563, 2.82700014114380, 0.83,
         0.44531310527793, 17 + 5.0 / 60, 11.02448526443880,
         26.01087501850882, 0.05, 1.32, 36, 5.80784207568089],
        # US-NE1 - IDL iteration 2
        [298.76438086370422, 302.49074951171878, 7.02662301063538,
         97.23062487251868, 350, 917.87845865885413,
         30.81263826599121 / (0.0864 / 24), 0, 0.45641128977509,
         0.83, 0.35, 0.95, 0.49, 0.13, 0.95,
         0.19908118247986, 0.84600001573563, 2.82700014114380, 0.83,
         0.44531310527793, 17 + 5.0 / 60, 11.02448526443880,
         26.01087501850882, 0.05, 1.32, 36, 6.78850125182576],
        # US-NE1 - IDL iteration 3
        [298.47296436495526, 302.49074951171878, 7.02662301063538,
         97.23062487251868, 350, 917.87845865885413,
         30.81263826599121 / (0.0864 / 24), 0, 0.45641128977509,
         0.83, 0.35, 0.95, 0.49, 0.13, 0.95,
         0.19908118247986, 0.84600001573563, 2.82700014114380, 0.83,
         0.44531310527793, 17 + 5.0 / 60, 11.02448526443880,
         26.01087501850882, 0.05, 1.32, 36, 6.69752524158102],
    ]
)
def test_tseb_pt(T_air, T_rad, u, p, z, Rs_1, Rs24, vza, zs, aleafv, aleafn,
                 aleafl, adeadv, adeadn, adeadl, albedo, ndvi, lai, clump, hc,
                 time, t_rise, t_end, leaf_width, a_PT_in, iterations,
                 expected, tol=1E-8):
    output_image = tseb.tseb_pt(
        ee.Image.constant(T_air), ee.Image.constant(T_rad),
        ee.Image.constant(u), ee.Image.constant(p), ee.Image.constant(z),
        ee.Image.constant(Rs_1), ee.Image.constant(Rs24),
        ee.Image.constant(vza), ee.Image.constant(zs),
        ee.Image.constant(aleafv), ee.Image.constant(aleafn),
        ee.Image.constant(aleafl), ee.Image.constant(adeadv),
        ee.Image.constant(adeadn), ee.Image.constant(adeadl),
        ee.Image.constant(albedo), ee.Image.constant(ndvi),
        ee.Image.constant(lai), ee.Image.constant(clump), ee.Image.constant(hc),
        ee.Image.constant(time), ee.Image.constant(t_rise),
        ee.Image.constant(t_end), ee.Image.constant(leaf_width),
        ee.Image.constant(a_PT_in),
        iterations)

    output = list(utils.constant_image_value(output_image).values())[0]
    logging.debug('\n  Target values: {}'.format(expected))
    logging.debug('  Output values: {}'.format(output))
    assert abs(output - expected) <= tol
