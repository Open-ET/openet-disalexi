import logging
import math

import ee
import pytest

import openet.disalexi.tseb_utils as tseb_utils
import openet.disalexi.utils as utils
from .idl_values import ne1, ne2, ne3


# AmeriFlux sites adjusted to nearest Landsat cell centroid
# ne1_xy = [-96.47672812080845, 41.16506126041818]
# ne2_xy = [-96.46994024736414, 41.16491226772292]
# ne3_xy = [-96.43968912903934, 41.17964494123755]


@pytest.mark.parametrize(
    'timestamp, expected',
    [
        # IDL only computes integer JD
        [ne1['timestamp'], ne1['jd']],
        [946728000000, 2451545],
    ]
)
def test_to_jd(timestamp, expected, tol=0.001):
    """"""
    date = ee.Date(timestamp)
    output = tseb_utils._to_jd(date).getInfo()

    logging.debug('\n  Target values: {}'.format(expected))
    logging.debug('  Output values: {}'.format(output))
    assert abs(output - expected) <= tol



@pytest.mark.parametrize(
    'timestamp, xy, expected',
    [
        [ne1['timestamp'], ne1['xy'], 0.5 * (ne1['t_rise'] + ne1['t_end'])],
        [ne2['timestamp'], ne2['xy'], 0.5 * (ne2['t_rise'] + ne2['t_end'])],
        [ne2['timestamp'], ne3['xy'], 0.5 * (ne3['t_rise'] + ne3['t_end'])],
    ]
)
def test_solar_noon_image(timestamp, xy, expected, tol=1E-6):
    """Check that the sunset_sunrise function works for real images"""
    output_images = tseb_utils.solar_noon(
        datetime=ee.Date(timestamp),
        # date=ee.Date(datetime.datetime.utcfromtimestamp(timestamp / 1000)),
        lon=ee.Image.pixelLonLat().select(['longitude']).multiply(math.pi / 180))
    output = utils.point_image_value(
        ee.Image(output_images).rename(['t_noon']), xy=xy)['t_noon']

    logging.debug('  Target values: {:.12f}'.format(expected))
    logging.debug('  Output values: {:.12f}'.format(output))
    assert abs(output - expected) <= tol


@pytest.mark.parametrize(
    'timestamp, xy, expected',
    [
        [ne1['timestamp'], ne1['xy'], 0.5 * (ne1['t_rise'] + ne1['t_end'])],
        [ne2['timestamp'], ne2['xy'], 0.5 * (ne2['t_rise'] + ne2['t_end'])],
        [ne2['timestamp'], ne3['xy'], 0.5 * (ne3['t_rise'] + ne3['t_end'])],
    ]
)
def test_solar_noon_constant(timestamp, xy, expected, tol=1E-6):
    """Check that the sunset_sunrise function works for constant images"""
    output_images = tseb_utils.solar_noon(
        datetime=ee.Date(timestamp),
        lon=ee.Image.constant(xy[0]).multiply(math.pi / 180))
    output = utils.constant_image_value(ee.Image(output_images))['t_noon']

    logging.debug('  Target values: {:.12f}'.format(expected))
    logging.debug('  Output values: {:.12f}'.format(output))
    assert abs(output - expected) <= tol


@pytest.mark.parametrize(
    'timestamp, xy, expected',
    [
        [ne1['timestamp'], ne1['xy'], ne1['zs']],
        [ne2['timestamp'], ne2['xy'], ne2['zs']],
        [ne3['timestamp'], ne3['xy'], ne3['zs']],
    ]
)
def test_solar_zenith_image(timestamp, xy, expected, tol=1E-6):
    """Check that the solar zenith function works for real images"""
    output_images = tseb_utils.solar_zenith(
        datetime=ee.Date(timestamp),
        lon=ee.Image.pixelLonLat().select(['longitude']).multiply(math.pi / 180),
        lat=ee.Image.pixelLonLat().select(['latitude']).multiply(math.pi / 180))
    output = utils.point_image_value(
        ee.Image(output_images).rename(['vs']), xy=xy)['vs']

    logging.debug('\n  Target values: {}'.format(expected))
    logging.debug('  Output values: {}'.format(output))
    assert abs(output - expected) <= tol


@pytest.mark.parametrize(
    'timestamp, xy, expected',
    [
        [ne1['timestamp'], ne1['xy'], ne1['zs']],
        [ne2['timestamp'], ne2['xy'], ne2['zs']],
        [ne3['timestamp'], ne3['xy'], ne3['zs']],
    ]
)
def test_solar_zenith_constant(timestamp, xy, expected, tol=1E-10):
    """Check that the solar_zenith function also works for constant images"""
    output_images = tseb_utils.solar_zenith(
        datetime=ee.Date(timestamp),
        lon=ee.Image.constant(xy[0]).multiply(math.pi / 180),
        lat=ee.Image.constant(xy[1]).multiply(math.pi / 180))
    output = utils.constant_image_value(
        ee.Image(output_images).rename(['vs']))['vs']

    logging.debug('\n  Target values: {}'.format(expected))
    logging.debug('  Output values: {}'.format(output))
    assert abs(output - expected) <= tol


@pytest.mark.parametrize(
    'T_air, expected',
    [
        # US-NE1
        [296.86259968439742, 0.76937715366883],
        # [298.76438086370422, 0.751758652383],
        # [291.9443, 0.751758652383],
        # [273.9443, 0.718960910317]
    ]
)
def test_emissivity(T_air, expected, tol=1E-8):
    output_image = tseb_utils.emissivity(ee.Image.constant(T_air))

    output = list(utils.constant_image_value(output_image).values())[0]
    logging.debug('\n  Target values: {}'.format(expected))
    logging.debug('  Output values: {}'.format(output))
    assert abs(output - expected) <= tol


@pytest.mark.parametrize(
    'albedo, Rs_1, F, fc, aleafv, aleafn, aleafl, adeadv, adeadn, adeadl, '
    'zs, albedo_iter, expected',
    [
        # US-NE1
        [0.19908118247986, 917.87845865885413, 2.34641011714935,
         0.69062621055586, 0.83, 0.35, 0.95, 0.49, 0.13, 0.95,
         ne1['zs'], 10,
         {'Rs_c': 563.86830655658594, 'Rs_s': 354.01015210226819,
          'albedo_c': 0.14601381346619, 'albedo_s': 0.30626749091907}],
        # US-NE2
        [0.21406339108944, 917.87921142578125, 0.85905007123947,
         0.34918186324148, 0.83, 0.35, 0.95, 0.49, 0.13, 0.95,
         ne2['zs'], 10,
         {'Rs_c': 278.58009845129851, 'Rs_s': 639.29911297448280,
          'albedo_c': 0.17126935827424, 'albedo_s': 0.22970060954277}],
        # US-NE3
        [0.21599538624287, 917.72406005859375, 0.92213005065918,
         0.36938832966689, 0.83, 0.35, 0.95, 0.49, 0.13, 0.95,
         ne3['zs'], 10,
         {'Rs_c': 294.95646120840320, 'Rs_s': 622.76759885019055,
          'albedo_c': 0.16847370243644, 'albedo_s': 0.22970059209214}],
    ]
)
def test_albedo_separation(albedo, Rs_1, F, fc, aleafv, aleafn, aleafl, adeadv,
                           adeadn, adeadl, zs, albedo_iter, expected,
                           tol=1E-10):
    output_images = tseb_utils.albedo_separation(
        ee.Image.constant(albedo), ee.Image.constant(Rs_1),
        ee.Image.constant(F), ee.Image.constant(fc), ee.Image.constant(aleafv),
        ee.Image.constant(aleafn), ee.Image.constant(aleafl),
        ee.Image.constant(adeadv), ee.Image.constant(adeadn),
        ee.Image.constant(adeadl), ee.Image.constant(zs), albedo_iter)
    output = utils.constant_image_value(ee.Image(output_images).rename(
        ['Rs_c', 'Rs_s', 'albedo_c', 'albedo_s']))

    for k in expected.keys():
        logging.debug('\n  {}'.format(k))
        logging.debug('  Target values: {}'.format(expected[k]))
        logging.debug('  Output values: {}'.format(output[k]))
        assert abs(output[k] - expected[k]) <= tol


@pytest.mark.parametrize(
    'u, d0, z0m, z_u, fm, expected',
    [
        # Test neutral conditions first
        # US-NE1
        [7.02662301063538, 0.29687540351862, 0.05477351194919, 50.0,
         0, 0.42300362871195],
        # US-NE2
        [7.02415990829468, 0.18306062108049, 0.03377468458935, 50.0,
         0, 0.39470232345209],
        # US-NE3
        [7.01548242568970, 0.18979610988896, 0.03501738227451, 50.0,
         0, 0.39618403124601],
    ]
)
def test_compute_u_attr(u, d0, z0m, z_u, fm, expected, tol=1E-10):
    output_image = tseb_utils.compute_u_attr(
        ee.Image.constant(u), ee.Image.constant(d0), ee.Image.constant(z0m),
        ee.Image.constant(z_u), ee.Image.constant(fm))

    output = list(utils.constant_image_value(output_image).values())[0]
    logging.debug('\n  Target values: {}'.format(expected))
    logging.debug('  Output values: {}'.format(output))
    assert abs(output - expected) <= tol


@pytest.mark.parametrize(
    'u_attr, d0, z0h, z_t, fh, expected',
    [
        # Test neutral conditions first
        # US-NE1
        [0.42300362871195, 0.29687540351862, 0.05477351194919, 50.0,
         0, 39.26977994736171],
        # US-NE2
        [0.39470232345209, 0.18306062108049, 0.03377468458935, 50.0,
         0, 45.08738255788060],
        # US-NE3
        [0.39618403124601, 0.18979610988896, 0.03501738227451, 50.0,
         0, 44.69548019943959],
        # Test conditionals
        [1, 0, 1, 1, 0, 500],
        [1, 0, 1, 1.5, 0, 1],
        # Made up values
        [0.5, 0.2, 0.05, 50, 0, 33.67681589065658],
    ]
)
def test_compute_r_ah(u_attr, d0, z0h, z_t, fh, expected, tol=1E-10):
    output_image = tseb_utils.compute_r_ah(
        ee.Image.constant(u_attr), ee.Image.constant(d0),
        ee.Image.constant(z0h), ee.Image.constant(z_t), ee.Image.constant(fh))

    output = list(utils.constant_image_value(output_image).values())[0]
    logging.debug('\n  Target values: {}'.format(expected))
    logging.debug('  Output values: {}'.format(output))
    assert abs(output - expected) <= tol


@pytest.mark.parametrize(
    'u_attr, T_s, T_c, hc, F, d0, z0m, leaf, leaf_s, fm_h, expected',
    [
        # Intentionally using LAI value for F to match Python code
        # Test neutral conditions first (use T_air for T_s and T_c)
        # US-NE1
        [0.42300362871195, 296.86259968439742, 296.86259968439742,
         0.44531310527793, 2.34641011714935, 0.29687540351862, 0.05477351194919,
         1.02484268608471, 0.12504030461746, 0, 111.48941465802088],
        # US-NE2
        [0.39470232345209, 296.86320701599124, 296.86320701599124,
         0.27459093162074, 1.03500008583069, 0.18306062108049, 0.03377468458935,
         0.44641597407594, 0.10642842315602, 0, 83.37415489800380],
        # US-NE3
        [0.39618403124601, 296.85014930725100, 296.85014930725100,
         0.28469416483344, 1.11100006103516, 0.18979610988896, 0.03501738227451,
         0.47368143806798, 0.10771802129930, 0, 84.57955576951095],
    ]
)
def test_compute_r_s(u_attr, T_s, T_c, hc, F, d0, z0m, leaf, leaf_s, fm_h,
                     expected, tol=1E-10):
    output_image = tseb_utils.compute_r_s(
        ee.Image.constant(u_attr), ee.Image.constant(T_s),
        ee.Image.constant(T_c), ee.Image.constant(hc), ee.Image.constant(F),
        ee.Image.constant(d0), ee.Image.constant(z0m), ee.Image.constant(leaf),
        ee.Image.constant(leaf_s), ee.Image.constant(fm_h))

    output = list(utils.constant_image_value(output_image).values())[0]
    logging.debug('\n  Target values: {}'.format(expected))
    logging.debug('  Output values: {}'.format(output))
    assert abs(output - expected) <= tol


@pytest.mark.parametrize(
    'u_attr, hc, F, d0, z0m, xl, leaf_c, fm_h, expected',
    [
        # Intentionally using LAI value for F to match Python code
        # Test neutral conditions first
        # US-NE1
        [0.42300362871195, 0.44531310527793, 2.82700014114380,
         0.29687540351862, 0.05477351194919, 0.05, 1.48517368527037,
         0, 15.95553315272425],
        # US-NE2
        [0.39470232345209, 0.27459093162074, 1.03500008583069,
         0.18306062108049, 0.03377468458935, 0.05, 1.01934611604499,
         0, 42.95938805710947],
        # US-NE3
        [0.39618403124601, 0.28469416483344, 1.11100006103516,
         0.18979610988896, 0.03501738227451, 0.05, 1.04179092390976,
         0, 40.04016661292570],
    ]
)
def test_compute_r_x(u_attr, hc, F, d0, z0m, xl, leaf_c, fm_h,
                     expected, tol=1E-10):
    output_image = tseb_utils.compute_r_x(
        ee.Image.constant(u_attr), ee.Image.constant(hc), ee.Image.constant(F),
        ee.Image.constant(d0), ee.Image.constant(z0m), ee.Image.constant(xl),
        ee.Image.constant(leaf_c), ee.Image.constant(fm_h))

    output = list(utils.constant_image_value(output_image).values())[0]
    logging.debug('\n  Target values: {}'.format(expected))
    logging.debug('  Output values: {}'.format(output))
    assert abs(output - expected) <= tol


@pytest.mark.parametrize(
    'albedo_c, T_air, T_c, T_s, e_atm, Rs_c, F, expected',
    [
        # US-NE1
        [0.14601381353606, 296.86260047912600, 296.86260047912600,
         315.05466853603883, 0.76937715672980, 563.86831631702330,
         2.34641011714935, 474.40464299790841],
        # US-NE2
        [0.17126935827424, 296.86320701599124, 296.86320701599124,
         310.81055724970071, 0.76937949288406, 278.58009845129851,
         0.85905007123947, 210.93979017388850],
        # US-NE3
        [0.16847370243644, 296.85014930725100, 296.85014930725100,
         312.65231892496638, 0.76932920743274, 294.95646120840320,
         0.92213005065918, 231.40032270633148],
    ]
)
def test_compute_Rn_c(albedo_c, T_air, T_c, T_s, e_atm, Rs_c, F,
                      expected, tol=1E-10):
    output_image = tseb_utils.compute_Rn_c(
        ee.Image.constant(albedo_c), ee.Image.constant(T_air),
        ee.Image.constant(T_c), ee.Image.constant(T_s),
        ee.Image.constant(e_atm), ee.Image.constant(Rs_c),
        ee.Image.constant(F))

    output = list(utils.constant_image_value(output_image).values())[0]
    logging.debug('\n  Target values: {}'.format(expected))
    logging.debug('  Output values: {}'.format(output))
    assert abs(output - expected) <= tol


@pytest.mark.parametrize(
    'albedo_s, T_air, T_c, T_s, e_atm, Rs_s, F, expected',
    [
        # US-NE1
        [0.30626749091907, 296.86260047912600, 296.86260047912600,
         315.05466853603883, 0.76937715672980, 354.01016268688295,
         2.34641011714935, 145.97048553078571],
        # US-NE2
        [0.22970060954277, 296.86320701599124, 296.86320701599124,
         310.81055724970071, 0.76937949288406, 639.29911297448280,
         0.85905007123947, 388.06596376956912],
        # US-NE3
        [0.22970059209214, 296.85014930725100, 296.85014930725100,
         312.65231892496638, 0.76932920743274, 622.76759885019055,
         0.92213005065918, 365.85694730307142],
    ]
)
def test_compute_Rn_s(albedo_s, T_air, T_c, T_s, e_atm, Rs_s, F,
                      expected, tol=1E-10):
    output_image = tseb_utils.compute_Rn_s(
        ee.Image.constant(albedo_s), ee.Image.constant(T_air),
        ee.Image.constant(T_c), ee.Image.constant(T_s),
        ee.Image.constant(e_atm), ee.Image.constant(Rs_s),
        ee.Image.constant(F))

    output = list(utils.constant_image_value(output_image).values())[0]
    logging.debug('\n  Target values: {}'.format(expected))
    logging.debug('  Output values: {}'.format(output))
    assert abs(output - expected) <= tol


@pytest.mark.parametrize(
    'Rn, Rn_s, albedo, ndvi, t_rise, t_end, time, EF_s, expected',
    [
        # US-NE1
        [620.37512852869418, 145.97048553078571, 0.19908118247986,
         0.84600001573563, ne1['t_rise'], ne1['t_end'], 17 + 5.0/60,
         0.0, 47.91926135598481],
        # US-NE2
        [599.00575394345765, 388.06596376956912, 0.21406339108944,
         0.51700001955032, ne2['t_rise'], ne2['t_end'], 17 + 5.0/60,
         0.0, 127.38965973286223],
        # US-NE3
        [597.25727000940287, 365.85694730307142, 0.21599538624287,
         0.54000002145767, ne3['t_rise'], ne3['t_end'], 17 + 5.0/60,
         0.0, 120.07887402772921],
    ]
)
def test_compute_G0(Rn, Rn_s, albedo, ndvi, t_rise, t_end, time, EF_s,
                    expected, tol=1E-10):
    output_image = tseb_utils.compute_G0(
        ee.Image.constant(Rn), ee.Image.constant(Rn_s),
        ee.Image.constant(albedo), ee.Image.constant(ndvi),
        ee.Image.constant(0.5 * (t_rise + t_end)),
        ee.Image.constant(time), ee.Image.constant(EF_s))

    output = list(utils.constant_image_value(output_image).values())[0]
    logging.debug('\n  Target values: {}'.format(expected))
    logging.debug('  Output values: {}'.format(output))
    assert abs(output - expected) <= tol


@pytest.mark.parametrize(
    'H_c, fc_q, T_air, t0, r_ah, r_s, r_x, r_air, expected',
    [
        # H_c is intermediate H_c
        # US-NE1
        [15.38257790294108, 0.69062621055586, 296.86260047912600,
         302.49074951171878, 39.26977861490872, 111.48941256212345,
         15.95553288203267, 1.13051648163638, 300.21526728574315],
        # US-NE2
        [6.83802281209108, 0.34918186324148, 296.86320701599124,
         305.94039550781253, 45.08738255788060, 83.37415489800380,
         42.95938805710947, 1.13051426564594, 301.32407278194165],
        # US-NE3
        [7.54152132085690, 0.36938832966689, 296.85014930725100,
         306.81518188476565, 44.69548019943959, 84.57955576951095,
         40.04016661292570, 1.13056197407484, 301.74371532209193],
        # Test fc conditionals
        [10, 0.01, 295, 300, 1, 100, 20, 1, 300],
        [10, 0.99, 295, 300, 50, 100, 20, 1, 300],
        # Test Tair conditionals (using fc conditional to set T_c)
        [10, 0.01, 295, 280, 1, 100, 20, 1, 295-10],
        [10, 0.99, 295, 400, 50, 100, 20, 1, 295+50],
    ]
)
def test_temp_separation_tc(H_c, fc_q, T_air, t0, r_ah, r_s, r_x, r_air,
                            expected, tol=1E-10):
    output_image = tseb_utils.temp_separation_tc(
        ee.Image.constant(H_c), ee.Image.constant(fc_q), ee.Image.constant(T_air),
        ee.Image.constant(t0), ee.Image.constant(r_ah), ee.Image.constant(r_s),
        ee.Image.constant(r_x), ee.Image.constant(r_air))

    output = list(utils.constant_image_value(output_image).values())[0]
    logging.debug('\n  Target values: {}'.format(expected))
    logging.debug('  Output values: {}'.format(output))
    assert abs(output - expected) <= tol


@pytest.mark.parametrize(
    'T_c, fc_q, T_air, t0, expected',
    [
        # US-NE1
        [300.21526728574315, 0.69062621055586, 296.86260047912600,
         302.49074951171878, 307.39290099869856],
        # US-NE2
        [301.32407278194165, 0.34918186324148, 296.86320701599124,
         305.94039550781253, 308.33345898581922],
        # US-NE3
        [301.74371532209193, 0.36938832966689, 296.85014930725100,
         306.81518188476565, 309.67283488882907],
        # Test fc conditionals
        [10, 0.01, 295, 300, 300],
        [10, 0.99, 295, 300, 300],
        # Test Tair conditionals (using fc conditional to set T_s)
        [10, 0.01, 295, 280, 295-10],
        [10, 0.99, 295, 400, 295+50],
    ]
)
def test_temp_separation_ts(T_c, fc_q, T_air, t0,
                            expected, tol=1E-10):
    output_image = tseb_utils.temp_separation_ts(
        ee.Image.constant(T_c), ee.Image.constant(fc_q),
        ee.Image.constant(T_air), ee.Image.constant(t0))

    output = list(utils.constant_image_value(output_image).values())[0]
    logging.debug('\n  Target values: {}'.format(expected))
    logging.debug('  Output values: {}'.format(output))
    assert abs(output - expected) <= tol


@pytest.mark.parametrize(
    'T_c, T_s, fc_q, T_air, r_ah, r_s, r_x, expected',
    [
        # US-NE1
        [300.21526728574315, 307.39290099869856, 0.68492516637839,
         296.86260047912600, 39.26977861490872, 111.48941256212345,
         15.95553288203267, 299.99905833599252],
        # US-NE2
        [301.32407278194165, 308.33345898581922, 0.34918186324148,
         296.86320701599124, 45.08738255788060, 83.37415489800380,
         42.95938805710947, 301.06529996155803],
        # US-NE3
        [301.74371532209193, 309.67283488882907, 0.36938832966689,
         296.85014930725100, 44.69548019943959, 84.57955576951095,
         40.04016661292570, 301.47772096219484],
    ]
)
def test_temp_separation_tac(T_c, T_s, fc_q, T_air, r_ah, r_s, r_x,
                             expected, tol=1E-6):
    output_image = tseb_utils.temp_separation_tac(
        ee.Image.constant(T_c), ee.Image.constant(T_s), ee.Image.constant(fc_q),
        ee.Image.constant(T_air), ee.Image.constant(r_ah),
        ee.Image.constant(r_s), ee.Image.constant(r_x))

    output = list(utils.constant_image_value(output_image).values())[0]
    logging.debug('\n  Target values: {}'.format(expected))
    logging.debug('  Output values: {}'.format(output))
    assert abs(output - expected) <= tol


@pytest.mark.parametrize(
    'H, t0, u_attr, r_air, z_t, d0, expected',
    [
        # US-NE1
        [90.66942780971830, 302.49074951171878, 0.42300362871195,
         1.13051648453994, 50.0, 0.29687540351862, 1.61571998736168],
        [118.59153871239729, 302.49074951171878, 0.49080214694589,
         1.13051648453994, 50.0, 0.29687540351862, 1.49124364135338],
        [117.97091195540997, 302.49074951171878, 0.48416306927365,
         1.13051648453994, 50.0, 0.29687540351862, 1.51586847161076],
        # Test L < -100
        [57.64681624293684, 302.49074951171878, 0.42300362871195,
         1.13051648453994, 50.0, 0.29687540351862, 0.00000000000000],
        # # US-NE2
        # [6.83816439594116 + 98.96279228312359, 1.13051426564594, 305.94039550781253,
        #  0.39470232345209, 50.0, 0.18306062108049, 1.87533132192393],
        # # US-NE3
        # [7.54177972507236 + 109.99853149146979, 1.13056197407484, 306.81518188476565,
        #  0.39618403124601, 50.0, 0.18979610988896, 1.94460807517059],
    ]
)
def test_compute_stability_fh(H, t0, u_attr, r_air, z_t, d0, expected,
                              tol=1E-10):
    output_image = tseb_utils.compute_stability_fh(
        ee.Image.constant(H), ee.Image.constant(t0), ee.Image.constant(u_attr),
        ee.Image.constant(r_air), ee.Image.constant(z_t), ee.Image.constant(d0))

    output = list(utils.constant_image_value(output_image).values())[0]
    logging.debug('\n  Target values: {}'.format(expected))
    logging.debug('  Output values: {}'.format(output))
    assert abs(output - expected) <= tol


@pytest.mark.parametrize(
    'H, t0, u_attr, r_air, z_u, d0, z0m, expected',
    [
        # US-NE1
        [90.66942780971830, 302.49074951171878, 0.42300362871195,
         1.13051648453994, 50.0, 0.29687540351862, 0.05477351194919,
         0.94080618919816],
        [118.59153871239729, 302.49074951171878, 0.49080214694589,
         1.13051648453994, 50.0, 0.29687540351862, 0.05477351194919,
         0.86031651924326],
        [117.97091195540997, 302.49074951171878, 0.48416306927365,
         1.13051648453994, 50.0, 0.29687540351862, 0.05477351194919,
         0.87614809449582],
        # Test L < -100
        [57.64681624293684, 302.49074951171878, 0.42300362871195,
         1.13051648453994, 50.0, 0.29687540351862, 0.05477351194919,
         0.00000000000000],
        # # US-NE2
        # [6.83816439594116 + 98.96279228312359, 305.94039550781253, 0.39470232345209,
        #  1.13051426564594, 50.0, 0.18306062108049, 0.03377468458935,
        #  1.11228299509438],
        # # US-NE3
        # [7.54177972507236 + 109.99853149146979, 306.81518188476565, 0.39618403124601,
        #  1.13056197407484, 50.0, 0.18979610988896, 0.03501738227451,
        #  1.15883632516281],
    ]
)
def test_compute_stability_fm(H, t0, u_attr, r_air, z_u, d0, z0m, expected,
                              tol=1E-10):
    output_image = tseb_utils.compute_stability_fm(
        ee.Image.constant(H), ee.Image.constant(t0), ee.Image.constant(u_attr),
        ee.Image.constant(r_air), ee.Image.constant(z_u),
        ee.Image.constant(d0), ee.Image.constant(z0m))

    output = list(utils.constant_image_value(output_image).values())[0]
    logging.debug('\n  Target values: {}'.format(expected))
    logging.debug('  Output values: {}'.format(output))
    assert abs(output - expected) <= tol


@pytest.mark.parametrize(
    'H, t0, u_attr, r_air, hc, d0, z0m, expected',
    [
        # US-NE1
        [90.66942780971830, 302.49074951171878, 0.42300362871195,
         1.13051648453994, 0.44531310527793, 0.29687540351862,
         0.05477351194919, 0.00824226853658],
        [118.59153871239729, 302.49074951171878, 0.49080214694589,
         1.13051648453994, 0.44531310527793, 0.29687540351862,
         0.05477351194919, 0.00691309800770],
        [117.97091195540997, 302.49074951171878, 0.48416306927365,
         1.13051648453994,
         0.44531310527793, 0.29687540351862,
         0.05477351194919, 0.00716149541408],
        # Test L < -100
        [57.64681624293684, 302.49074951171878, 0.42300362871195,
         1.13051648453994, 0.44531310527793, 0.29687540351862,
         0.05477351194919, 0.00000000000000],
        # # US-NE2
        # [6.83816439594116 + 98.96279228312359, 305.94039550781253, 0.39470232345209,
        #  1.13051426564594, 0.27459093162074, 0.18306062108049,
        #  0.0079049231986360, 0.00722679542184],
        # # US-NE3
        # [7.54177972507236 + 109.99853149146979, 306.81518188476565, 0.39618403124601,
        #  1.13056197407484, 0.28469416483344, 0.18979610988896,
        #  0.03501738227451, 0.00819727197171],
    ]
)
def test_compute_stability_fm_h(H, t0, u_attr, r_air, hc, d0, z0m, expected,
                                tol=1E-10):
    output_image = tseb_utils.compute_stability_fm_h(
        ee.Image.constant(H), ee.Image.constant(t0), ee.Image.constant(u_attr),
        ee.Image.constant(r_air), ee.Image.constant(hc), ee.Image.constant(d0),
        ee.Image.constant(z0m))

    output = list(utils.constant_image_value(output_image).values())[0]
    logging.debug('\n  Target values: {}'.format(expected))
    logging.debug('  Output values: {}'.format(output))
    assert abs(output - expected) <= tol
