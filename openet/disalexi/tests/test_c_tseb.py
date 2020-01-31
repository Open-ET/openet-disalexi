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
    'T_air, T_rad, e_air, u,'
    'p, z, Rs_1, '
    'Rs24, vza,'
    'aleafv, aleafn, aleafl, adeadv, adeadn, adeadl,'
    'albedo, ndvi, lai, clump,'
    'hc_min, hc_max, leaf_width, datetime, lon, lat, a_PT_in, '
    'stabil_iter, albedo_iter, expected, tol',
    [
        # # US-NE1 - IDL iteration 1
        # [296.86259968439742, 302.49074951171878, e_air, 7.02662301063538,
        #  97.23062487251868, 350, 917.87845865885413,
        #  30.81263826599121 / (0.0864 / 24), 0,
        #  0.83, 0.35, 0.95, 0.49, 0.13, 0.95,
        #  0.19908118247986, 0.84600001573563, 2.82700014114380, 0.83,
        #  0.44531310527793, 0.05, 1404839150550, ne1_xy[0], ne1_xy[1], 1.32,
        #  36, 10,5.80784207568089 / 0.408, 1E-8],
        # # US-NE1 - IDL iteration 2
        # [298.76438086370422, 302.49074951171878, e_air, 7.02662301063538,
        #  97.23062487251868, 350, 917.87845865885413,
        #  30.81263826599121 / (0.0864 / 24), 0,
        #  0.83, 0.35, 0.95, 0.49, 0.13, 0.95,
        #  0.19908118247986, 0.84600001573563, 2.82700014114380, 0.83,
        #  0.44531310527793, 0.05, 1404839150550, ne1_xy[0], ne1_xy[1], 1.32,
        #  36, 10, 6.78850125182576 / 0.408, 1E-8],
        # # US-NE1 - IDL iteration 3
        # [298.47296436495526, 302.49074951171878, e_air, 7.02662301063538,
        #  97.23062487251868, 350, 917.87845865885413,
        #  30.81263826599121 / (0.0864 / 24), 0,
        #  0.83, 0.35, 0.95, 0.49, 0.13, 0.95,
        #  0.19908118247986, 0.84600001573563, 2.82700014114380, 0.83,
        #  0.44531310527793, 0.05, 1404839150550, ne1_xy[0], ne1_xy[1], 1.32,
        #  36, 10, 6.69752524158102 / 0.408, 1E-8],

        # High NDVI site in LC08_044033_20170716
        [300.0, 305.92253850611, 0.84104, 3.2665367230039,
         101.264543111959, 3.0, 946.69066527778,
         8603.212890625, 0,
         0.83, 0.35, 0.95, 0.49, 0.13, 0.95,
         0.1259961, 0.87439300744578, 4.6797005913579, 0.83,
         0.0, 0.6, 0.05, 1500230731090, -121.5265, 38.7399, 1.32,
         36, 10, 6.9511, 1E-3],
        # Low NDVI site in LC08_044033_20170716
        [300.0, 323.59893135545, 0.84104, 3.2665367230039,
         101.252726383124, 4.0, 946.69066527778,
         8603.212890625, 0,
         0.83, 0.35, 0.95, 0.49, 0.13, 0.95,
         0.1716302, 0.16195230171936, 0.029734416071998, 0.83,
         0.0, 0.6, 0.05, 1500230731090, -121.50822, 38.71776, 1.32,
         36, 10, 4.1705, 1E-3],
    ]
)
def test_tseb_pt(T_air, T_rad, e_air, u, p, z, Rs_1, Rs24, vza, aleafv, aleafn,
                 aleafl, adeadv, adeadn, adeadl, albedo, ndvi, lai, clump,
                 hc_min, hc_max, leaf_width, datetime, lon, lat, a_PT_in,
                 stabil_iter, albedo_iter, expected, tol):
    output_image = tseb.tseb_pt(
        t_air=ee.Image.constant(T_air), t_rad=ee.Image.constant(T_rad),
        e_air=ee.Image.constant(e_air), u=ee.Image.constant(u),
        p=ee.Image.constant(p), z=ee.Image.constant(z),
        rs_1=ee.Image.constant(Rs_1), rs24=ee.Image.constant(Rs24),
        vza=ee.Image.constant(vza),
        aleafv=ee.Image.constant(aleafv), aleafn=ee.Image.constant(aleafn),
        aleafl=ee.Image.constant(aleafl), adeadv=ee.Image.constant(adeadv),
        adeadn=ee.Image.constant(adeadn), adeadl=ee.Image.constant(adeadl),
        albedo=ee.Image.constant(albedo), ndvi=ee.Image.constant(ndvi),
        lai=ee.Image.constant(lai),
        clump=ee.Image.constant(clump), leaf_width=ee.Image.constant(leaf_width),
        hc_min=ee.Image.constant(hc_min), hc_max=ee.Image.constant(hc_max),
        lon=ee.Image.constant(lon), lat=ee.Image.constant(lat),
        datetime=ee.Date(datetime), a_pt_in=ee.Image.constant(a_PT_in),
        stabil_iter=stabil_iter, albedo_iter=albedo_iter)

    output = list(utils.constant_image_value(output_image).values())[0]
    logging.debug('\n  Target values: {}'.format(expected))
    logging.debug('  Output values: {}'.format(output))
    assert abs(output - expected) <= tol
