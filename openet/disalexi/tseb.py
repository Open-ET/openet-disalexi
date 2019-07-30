import math

import ee

from . import tseb_utils
from . import utils

# DEADBEEF
# test_xy = [-121.50822, 38.71776] # Low NDVI
# test_xy = [-121.5265, 38.7399] # High NDVI A
test_xy = [-121.51465, 38.71846] # High NDVI B
def debug(x_var, x_str, xy=test_xy, scale=1):
    print('{:12s} {:>20.14}'.format(x_str+':', float(utils.point_image_value(
        x_var.rename(['test']), xy=xy, scale=scale)['test'])))


def tseb_pt(t_air, t_rad, u, p, z, rs_1, rs24, vza,
            aleafv, aleafn, aleafl, adeadv, adeadn, adeadl,
            albedo, ndvi, lai, clump, leaf_width, hc_min, hc_max,
            datetime, lon=None, lat=None, a_pt_in=1.32,
            stabil_iter=36, albedo_iter=10):
    """Priestley-Taylor TSEB

    Calculates the Priestley Taylor TSEB fluxes using a single observation of
    composite radiometric temperature and using resistances in series.

    Parameters
    ----------
    t_air : ee.Image
        Air temperature [K].
    t_rad : ee.Image
        Radiometric composite temperature [K].
    u : ee.Image
        Wind speed above the canopy [m s-1].
    p : ee.Image
        Atmospheric pressure [kPa].
    z : ee.Image
        Elevation [m].
    rs_1 : ee.Image
        Overpass insolation [w m-2].
    rs24 : ee.Image
        Daily insolation [w m-2].
    vza : float
        View zenith angle [radians].
    aleafv : ee.Image

    aleafn : ee.Image

    aleafl : ee.Image

    adeadv : ee.Image

    adeadn : ee.Image

    adeadl : ee.Image

    albedo : ee.Image

    ndvi : ee.Image
        Normalized Difference Vegetation Index
    lai : ee.Image
        Effective Leaf Area Index [m2 m-2].
    clump : ee.Image

    leaf_width : ee.Image
        Average/effective leaf width [m].
    hc_min : ee.Image
        Canopy height [m].
    hc_max : ee.Image
        Canopy height [m].
    datetime : ee.Date
        Image datetime.
    lat : ee.Image
        Latitude [deg].  If not set will default to ee.Image.pixelLonLat().
    lon : ee.Image
        Longitude [deg].  If not set will default to ee.Image.pixelLonLat().
    a_pt_in : float, optional
        Priestley Taylor coefficient for canopy potential transpiration
        (the default is 1.32).
    stabil_iter: int, optional
        Number of iterations of stability calculation (the default is 36).
    albedo_iter: int, optional
        Number of iterations of albedo separation calculation
        (the default is 10).

    Returns
    -------
    ET : ee.Image
        Evapotranspiration [mm].

    References
    ----------
    .. [Norman1995] J.M. Norman, W.P. Kustas, & K.S. Humes (1995),
        Source approach for estimating soil and vegetation energy fluxes in
        observations of directional radiometric surface temperature,
        Agricultural and Forest Meteorology,
        Volume 77, Issues 3-4, Pages 263-293,
        http://dx.doi.org/10.1016/0168-1923(95)02265-Y.
    .. [Kustas1999] W.P. Kustas, & J.M. Norman (1999), Evaluation of soil
        and vegetation heat flux predictions using a simple two-source
        model with radiometric temperatures for partial canopy cover,
        Agricultural and Forest Meteorology, Volume 94, Issue 1, Pages 13-29,
        http://dx.doi.org/10.1016/S0168-1923(99)00005-2.

    """

    # # DEADBEEF
    # print('\nInputs')
    # debug(t_air, 't_air')
    # debug(t_rad, 't_rad')
    # debug(u, 'u')
    # debug(p, 'p')
    # debug(z, 'z')
    # debug(rs_1, 'Rs_1')
    # debug(rs24, 'Rs24')
    # debug(albedo, 'albedo')
    # debug(ndvi, 'ndvi')
    # debug(lai, 'lai')
    # debug(aleafv, 'aleafv')
    # debug(aleafn, 'aleafn')
    # debug(aleafl, 'aleafl')
    # debug(adeadv, 'adeadv')
    # debug(adeadn, 'adeadn')
    # debug(adeadl, 'adeadl')
    # debug(clump, 'clump')
    # debug(leaf_width, 'leaf_width')
    # debug(hc_min, 'hc_min')
    # debug(hc_max, 'hc_max')

    # ************************************************************************
    # Apply met bands directly to Landsat image
    # CGM - This can probably be removed if Rs is resampled/smoothed in disalexi.py
    rs_1 = lai.multiply(0).add(rs_1).rename(['rs'])
    rs24 = lai.multiply(0).add(rs24).rename(['rs'])
    # t_air = lai.multiply(0).add(t_air).rename(['ta'])
    # u = lai.multiply(0).add(u).rename(['windspeed'])

    # ************************************************************************
    # CGM - Moved from disalexi.py to here since these are not parameters
    #   you would typically vary and could be functions of the input images.
    # CGM - time and t_noon could probably be moved into compute_G0 since that
    #   is the only place they are used.
    time = ee.Date(datetime).get('hour')\
        .add(ee.Date(datetime).get('minute').divide(60))
    if lat is None:
        lat = ee.Image.pixelLonLat().select(['latitude'])
    if lon is None:
        lon = ee.Image.pixelLonLat().select(['longitude'])
    t_noon = tseb_utils.solar_noon(
        datetime=datetime, lon=lon.multiply(math.pi / 180))
    zs = tseb_utils.solar_zenith(
        datetime=datetime, lon=lon.multiply(math.pi / 180),
        lat=lat.multiply(math.pi / 180))

    # ************************************************************************
    # Correct Clumping Factor
    f_green = 1.

    # LAI for leaf spherical distribution
    F = lai.multiply(clump)
    # F = lai.expression('lai * clump', {'lai': lai, 'clump': clump})

    # Fraction cover at nadir (view=0)
    fc = F.multiply(-0.5).exp().multiply(-1).add(1.0).clamp(0.01, 0.9)
    # fc = lai.expression('1.0 - exp(-0.5 * F)', {'F': F}).clamp(0.01, 0.9)

    # Compute canopy height and roughness parameters
    # CGM - Moved from _set_landcover_vars()
    hc = hc_max.subtract(hc_min).multiply(fc).add(hc_min)
    # hc = lai.expression(
    #     'hc_min + ((hc_max - hc_min) * fc)',
    #     {'hc_min': hc_min, 'hc_max': hc_max, 'fc': fc})

    # LAI relative to canopy projection only
    lai_c = lai.divide(fc)
    # lai_c = lai.expression('lai / fc', {'lai': lai, 'fc': fc})

    # Houborg modification (according to Anderson et al. 2005)
    fc_q = lai \
        .expression('1 - (exp(-0.5 * F / cos(vza)))', {'F': F, 'vza': vza}) \
        .clamp(0.05, 0.90)

    # Brutsaert (1982)
    z0m = hc.multiply(0.123)
    # z0m = hc.expression('hc * 0.123', {'hc': hc})
    # CGM - add(0) is to mimic numpy copy, check if needed
    z0h = z0m.add(0)
    d0 = hc.multiply(2.0 / 3.0)
    # d0 = hc.expression('hc * (2.0 / 3.0)', {'hc': hc})

    # Correction of roughness parameters for bare soils (F < 0.1)
    d0 = d0.where(F.lte(0.1), 0.00001)
    z0m = z0m.where(F.lte(0.1), 0.01)
    z0h = z0h.where(F.lte(0.1), 0.0001)

    # Correction of roughness parameters for water bodies
    # (NDVI < 0 and albedo < 0.05)
    water_mask = ndvi.lte(0).And(albedo.lte(0.05))
    d0 = d0.where(water_mask, 0.00001)
    z0m = z0m.where(water_mask, 0.00035)
    z0h = z0h.where(water_mask, 0.00035)

    # Check to avoid division by 0 in the next computations
    z0h = z0h.where(z0h.eq(0), 0.001)
    z0m = z0m.where(z0m.eq(0), 0.01)

    # DEADBEEF
    # z_u = ee.Number(50.0)
    # z_t = ee.Number(50.0)
    z_u = ee.Image.constant(50.0)
    z_t = ee.Image.constant(50.0)
    # z_u = lai.multiply(0).add(50)
    # z_t = lai.multiply(0).add(50)

    # Parameters for In-Canopy Wind Speed Extinction
    leaf = lai.expression(
        '(0.28 * (F ** (0.66667)) * (hc ** (0.33333)) * '
        '(leaf_width ** (-0.33333)))',
        {'F': F, 'hc': hc, 'leaf_width': leaf_width})
    leaf_c = lai.expression(
        '(0.28 * (lai_c ** (0.66667)) * (hc ** (0.33333)) * '
        '(leaf_width ** (-0.33333)))',
        {'lai_c': lai_c, 'hc': hc, 'leaf_width': leaf_width})
    leaf_s = lai.expression(
        '(0.28 * (0.1 ** (0.66667)) * (hc ** (0.33333)) * '
        '(leaf_width ** (-0.33333)))',
        {'hc': hc, 'leaf_width': leaf_width})

    # ************************************************************************
    # Atmospheric Parameters
    # Saturation vapour pressure [kPa] (FAO56 3-8)
    e_s = t_air.expression(
        '0.6108 * exp((17.27 * (t_air - 273.16)) / ((t_air - 273.16) + 237.3))',
        {'t_air': t_air})
    # Slope of the saturation vapor pressure [kPa] (FAO56 3-9)
    Ss = t_air.expression(
        '4098. * e_s / (((t_air - 273.16) + 237.3) ** 2)',
        {'e_s': e_s, 't_air': t_air})
    # Latent heat of vaporization (~2.45 at 20 C) [MJ kg-1] (FAO56 3-1)
    lambda1 = t_air.expression(
        '(2.501 - (2.361e-3 * (t_air - 273.16)))',
        {'t_air': t_air})
    # Psychrometric constant [kPa C-1] (FAO56 3-10)
    g = p.expression('1.615E-3 * p / lambda1', {'p': p, 'lambda1': lambda1})

    # ************************************************************************
    # Initialization of
    a_pt = albedo.multiply(0).add(a_pt_in)
    # a_pt = ee.Image.constant(a_pt_in)
    # a_pt = mask.multiply(a_pt)

    # CGM - This was also being computed inside albedo_separation function below
    # Commented out from here for now.
    # e_atm = t_air.expression(
    #     '1.0 - (0.2811 * (exp(-0.0003523 * ((t_air - 273.16) ** 2))))',
    #     {'t_air': t_air})

    Rs_c, Rs_s, albedo_c, albedo_s = tseb_utils.albedo_separation(
        albedo, rs_1, F, fc, aleafv, aleafn, aleafl, adeadv, adeadn, adeadl,
        zs, albedo_iter)

    # CGM - Moved emissivity calculation to separate function.
    #   I removed the Rs0 check.
    e_atm = tseb_utils.emissivity(t_air)
    # p = t_air.expression(
    #     '101.3 * (((t_air - (0.0065 * z)) / t_air) ** 5.26)',
    #     {'t_air': t_air, 'z': z})
    # Density of air? (kg m-3)
    r_air = t_air.expression(
        '101.3 * (((t_air - (0.0065 * z)) / t_air) ** 5.26) / 1.01 / t_air / 0.287',
        {'t_air': t_air, 'z': z})
    cp = ee.Number(1004.16)
    # cp = ee.Image.constant(1004.16)

    # Assume neutral conditions on first iteration (use t_air for Ts and Tc)
    # CGM - Using lai for F to match Python code
    u_attr = tseb_utils.compute_u_attr(
        u=u, d0=d0, z0m=z0m, z_u=z_u, fm=0)
    r_ah = tseb_utils.compute_r_ah(
        u_attr=u_attr, d0=d0, z0h=z0h, z_t=z_t, fh=0)
    # CGM - Why is this function is passing "lai" to "F"?
    r_s = tseb_utils.compute_r_s(
        u_attr=u_attr, T_s=t_air, T_c=t_air, hc=hc, F=lai, d0=d0, z0m=z0m,
        leaf=leaf, leaf_s=leaf_s, fm_h=0)
    r_x = tseb_utils.compute_r_x(
        u_attr=u_attr, hc=hc, F=lai, d0=d0, z0m=z0m, xl=leaf_width,
        leaf_c=leaf_c, fm_h=0)
    # r_ah, r_s, r_x, u_attr = tseb_utils.compute_resistance(
    #     u, t_air, t_air, hc, lai, d0, z0m, z0h, z_u, z_t, leaf_width, leaf,
    #     leaf_s, leaf_c, 0, 0, 0)

    T_c = t_air.multiply(1)
    # DEADBEEF - In IDL, this calculation is in C, not K?
    T_s = lai.expression(
        '(((t_rad - 273.16) - (fc_q * (T_c - 273.16))) / (1 - fc_q)) + 273.16',
        {'t_rad': t_rad, 'T_c': T_c, 'fc_q': fc_q})
    # T_s = lai.expression(
    #     '(t_rad - (fc_q * T_c)) / (1 - fc_q)',
    #     {'t_rad': t_rad, 'T_c': T_c, 'fc_q': fc_q})

    # CGM - Initialize to match t_air shape
    # This doesn't seem to do anything, commenting out for now
    # H_iter = t_air.multiply(0).add(200.16)
    EF_s = t_air.multiply(0)

    # # DEADBEEF
    # debug(t_noon, 't_noon')
    # print('{:12s} {:>20.14}'.format('time:', time.getInfo()))
    # debug(zs, 'zs')
    # debug(F, 'F')
    # debug(fc, 'fc')
    # debug(hc, 'hc')
    # debug(lai_c, 'lai_c')
    # debug(fc_q, 'fc_q')
    # debug(d0, 'd0')
    # debug(z0h, 'z0h')
    # debug(z0m, 'z0m')
    # debug(leaf, 'leaf')
    # debug(leaf_c, 'leaf_c')
    # debug(leaf_s, 'leaf_s')
    # debug(e_s, 'e_s')
    # debug(Ss, 'Ss')
    # debug(lambda1, 'lambda1')
    # debug(g, 'g')
    # debug(Rs_c, 'Rs_c')
    # debug(Rs_s, 'Rs_s')
    # debug(albedo_c, 'albedo_c')
    # debug(albedo_s, 'albedo_s')
    # debug(e_atm, 'e_atm')
    # debug(r_air, 'r_air')
    # debug(u_attr, 'u_attr')
    # debug(r_ah, 'r_ah')
    # debug(r_s, 'r_s')
    # debug(r_x, 'r_x')
    # debug(T_c, 'T_c')
    # debug(T_s, 'T_s')
    # debug(EF_s, 'EF_s')

    # ************************************************************************
    # Start Loop for Stability Correction and Water Stress
    def iter_func(n, prev):
        # Extract inputs from previous iteration
        a_pt_iter = ee.Image(ee.Dictionary(prev).get('a_pt'))
        EF_s_iter = ee.Image(ee.Dictionary(prev).get('EF_s'))
        r_ah_iter = ee.Image(ee.Dictionary(prev).get('r_ah'))
        r_s_iter = ee.Image(ee.Dictionary(prev).get('r_s'))
        r_x_iter = ee.Image(ee.Dictionary(prev).get('r_x'))
        T_c_iter = ee.Image(ee.Dictionary(prev).get('T_c'))
        T_s_iter = ee.Image(ee.Dictionary(prev).get('T_s'))
        u_attr_iter = ee.Image(ee.Dictionary(prev).get('u_attr'))

        Rn_c = tseb_utils.compute_Rn_c(
            albedo_c, t_air, T_c_iter, T_s_iter, e_atm, Rs_c, F)
        Rn_s = tseb_utils.compute_Rn_s(
            albedo_s, t_air, T_c_iter, T_s_iter, e_atm, Rs_s, F)
        Rn = Rn_c.add(Rn_s)
        # Rn_s, Rn_c, Rn = tseb_utils.compute_Rn(
        #     albedo_c, albedo_s, t_air, T_c_iter, T_s_iter, e_atm, Rs_c, Rs_s, F)

        G = tseb_utils.compute_G0(
            Rn, Rn_s, albedo, ndvi, t_noon, time, EF_s_iter)

        LE_c_init = albedo \
            .expression(
                'f_green * (a_pt * Ss / (Ss + g)) * Rn_c',
                {'f_green': f_green, 'a_pt': a_pt_iter, 'Ss': Ss, 'g': g,
                 'Rn_c': Rn_c}) \
            .max(0)
        H_c_init = Rn_c.subtract(LE_c_init)
        # H_c = albedo.expression(
        #     'Rn_c - LE_c', {'Rn_c': Rn_c, 'LE_c': LE_c})

        T_c_iter = tseb_utils.temp_separation_tc(
            H_c_init, fc_q, t_air, t_rad, r_ah_iter, r_s_iter, r_x_iter, r_air, cp)
        T_s_iter = tseb_utils.temp_separation_ts(T_c_iter, fc_q, t_air, t_rad)
        T_ac = tseb_utils.temp_separation_tac(
            T_c_iter, T_s_iter, fc_q, t_air, r_ah_iter, r_s_iter, r_x_iter)
        # T_c_iter, T_s_iter, T_ac = tseb_utils.temp_separation(
        #     H_c, fc_q, t_air, t_rad, r_ah_iter, r_s_iter, r_x_iter, r_air, cp)

        H_s = albedo.expression(
            'r_air * cp * (T_s - T_ac) / r_s',
            {'r_air': r_air, 'cp': cp, 'T_s': T_s_iter, 'T_ac': T_ac,
             'r_s': r_s_iter})
        H_c = albedo.expression(
            'r_air * cp * (T_c - T_ac) / r_x',
            {'r_air': r_air, 'cp': cp, 'T_c': T_c_iter, 'T_ac': T_ac,
             'r_x': r_x_iter})
        H = albedo.expression('H_s + H_c', {'H_s': H_s, 'H_c': H_c})

        LE_s = Rn_s.subtract(G).subtract(H_s)
        LE_c = Rn_c.subtract(H_c)
        # LE_s = albedo.expression(
        #     'Rn_s - G - H_s', {'Rn_s': Rn_s, 'G': G, 'H_s': H_s})
        # LE_c = albedo.expression('Rn_c - H_c', {'Rn_c': Rn_c, 'H_c': H_c})

        # CGM - Is there a reason this isn't up with the H calculation?
        H = H.where(H.eq(0), 10.0)

        # CGM - This wont doing anything at this position in the code.
        #   Commenting out for now.
        # r_ah_iter = r_ah_iter.where(r_ah_iter.eq(0), 10.0)

        # CGM - This doesn't seem to do anything, commenting out for now
        # mask_iter = H_iter.divide(H).lte(1.05).And(H_iter.divide(H).gte(0.95))
        # chk_iter = np.sum(mask_iter) / np.size(mask_iter)

        fh = tseb_utils.compute_stability_fh(
            H, t_rad, u_attr_iter, r_air, z_t, d0, cp)
        fm = tseb_utils.compute_stability_fm(
            H, t_rad, u_attr_iter, r_air, z_u, d0, z0m, cp)
        fm_h = tseb_utils.compute_stability_fm_h(
            H, t_rad, u_attr_iter, r_air, hc, d0, z0m, cp)
        # CGM - z0h is not used in this function, should it be?
        # fm, fh, fm_h = tseb_utils.compute_stability(
        #     H, t_rad, r_air, cp, u_attr, z_u, z_t, hc, d0, z0m, z0h)

        u_attr_iter = tseb_utils.compute_u_attr(
            u=u, d0=d0, z0m=z0m, z_u=z_u, fm=fm)
        r_ah_iter = tseb_utils.compute_r_ah(
            u_attr=u_attr_iter, d0=d0, z0h=z0h, z_t=z_t, fh=fh)
        r_s_iter = tseb_utils.compute_r_s(
            u_attr=u_attr_iter, T_s=T_s_iter, T_c=T_c_iter, hc=hc, F=lai,
            d0=d0, z0m=z0m, leaf=leaf, leaf_s=leaf_s, fm_h=fm_h)
        # CGM - Why is this function is passing "lai" to "F"?
        r_x_iter = tseb_utils.compute_r_x(
            u_attr=u_attr_iter, hc=hc, F=lai, d0=d0, z0m=z0m, xl=leaf_width,
            leaf_c=leaf_c, fm_h=fm_h)
        # r_ah_iter, r_s_iter, r_x_iter, u_attr_iter = tseb_utils.compute_resistance(
        #     u, T_s_iter, T_c_iter, hc, lai, d0, z0m, z0h, z_u, z_t,
        #     leaf_width, leaf, leaf_s, leaf_c, fm, fh, fm_h)

        a_pt_iter = a_pt_iter \
            .where(LE_s.lte(0), a_pt_iter.subtract(0.05)) \
            .where(a_pt_iter.lte(0), 0.01)

        den_s = albedo.expression('Rn_s - G', {'Rn_s': Rn_s, 'G': G})
        den_s = den_s.updateMask(den_s.neq(0))
        # den_s[den_s == 0.] = np.nan

        EF_s_iter = LE_s.divide(den_s)
        # EF_s_iter = albedo.expression(
        #     'LE_s / den_s', {'LE_s': LE_s, 'den_s': den_s})

        return ee.Dictionary({
            'a_pt': a_pt_iter, 'EF_s': EF_s_iter, 'G': G,
            'H_c': H_c, 'H_s': H_s, 'LE_c': LE_c, 'LE_s': LE_s,
            'Rn_c': Rn_c, 'Rn_s': Rn_s,
            'r_ah': r_ah_iter, 'r_s': r_s_iter, 'r_x': r_x_iter,
            'T_ac': T_ac, 'T_c': T_c_iter, 'T_s': T_s_iter,
            'u_attr': u_attr_iter})

    # Iterate the function n times
    # CGM - Iteration count is an input to the function
    input_images = ee.Dictionary({
        'a_pt': a_pt, 'EF_s': EF_s, 'G': ee.Image(0),
        'H_c': ee.Image(0), 'H_s': ee.Image(0),
        'LE_c': ee.Image(0), 'LE_s': ee.Image(0),
        'Rn_c': ee.Image(0), 'Rn_s': ee.Image(0),
        'r_ah': r_ah, 'r_s': r_s, 'r_x': r_x,
        'T_ac': ee.Image(0), 'T_c': T_c, 'T_s': T_s, 'u_attr': u_attr
    })
    iter_output = ee.Dictionary(
        ee.List.sequence(1, stabil_iter).iterate(iter_func, input_images))

    # Unpack the iteration output
    a_pt = ee.Image(iter_output.get('a_pt'))
    Rn_c = ee.Image(iter_output.get('Rn_c'))
    Rn_s = ee.Image(iter_output.get('Rn_s'))
    G = ee.Image(iter_output.get('G'))
    H_c = ee.Image(iter_output.get('H_c'))
    H_s = ee.Image(iter_output.get('H_s'))
    LE_c = ee.Image(iter_output.get('LE_c'))
    LE_s = ee.Image(iter_output.get('LE_s'))

    # # DEADBEEF
    # print('\nAfter Stability Iteration')
    # # T_ac = ee.Image(iter_output.get('T_ac'))
    # # T_c = ee.Image(iter_output.get('T_c'))
    # # T_s = ee.Image(iter_output.get('T_s'))
    # # r_ah = ee.Image(iter_output.get('r_ah'))
    # # r_s = ee.Image(iter_output.get('r_s'))
    # # r_x = ee.Image(iter_output.get('r_x'))
    # # debug(T_ac, 'T_ac')
    # # debug(T_c, 'T_c')
    # # debug(T_s, 'T_s')
    # # debug(r_ah, 'r_ah')
    # # debug(r_s, 'r_s')
    # # debug(r_x, 'r_x')
    # debug(a_pt, 'a_pt')
    # debug(Rn_c, 'Rn_c')
    # debug(Rn_s, 'Rn_s')
    # debug(Rn_c.add(Rn_s), 'Rn')
    # debug(G, 'G')
    # debug(H_c, 'H_c')
    # debug(H_s, 'H_s')
    # debug(H_c.add(H_s), 'H')
    # debug(LE_c, 'LE_c')
    # debug(LE_s, 'LE_s')
    # debug(LE_c.add(LE_s), 'LE')

    # ************************************************************************
    # Check Energy Balance Closure
    ind = a_pt.lte(0.01)
    LE_s = LE_s.where(ind, 1.0)
    LE_c = LE_c.where(ind, 1.0)
    G = G.where(ind, Rn_s.subtract(H_s))

    ind = LE_s.gt(Rn_s)
    LE_s = LE_s.where(ind,  Rn_s)
    H_s = H_s.where(ind,  Rn_s.subtract(G).subtract(LE_s))

    # CGM - Check order of operations
    ind = LE_c.gt(Rn_c.add(100))
    # CGM - Not used below since LE_c is recomputed
    # LE_c = LE_c.where(ind, Rn_c.add(100))
    H_c = H_c.where(ind, -100)

    LE_s = Rn_s.subtract(G).subtract(H_s)
    LE_c = Rn_c.subtract(H_c)

    # The latent heat of vaporization is 2.45 MJ kg-1
    # Assume rs24 is still in W m-2 day-1 and convert to MJ kg-1
    # CGM - Leaving out scaling value for now
    ET = albedo \
        .expression(
            '((LE_c + LE_s) / rs_1) * (rs24 / 2.45) * scaling',
            {'LE_c': LE_c, 'LE_s': LE_s, 'rs_1': rs_1,
             'rs24': rs24.multiply(0.0864 / 24.0), 'scaling': 1}) \
        .max(0.01) \

    # # DEADBEEF
    # print('\nAfter Checking EBC')
    # debug(G, 'G')
    # debug(H_c, 'H_c')
    # debug(H_s, 'H_s')
    # debug(LE_c, 'LE_c')
    # debug(LE_s, 'LE_s')
    # debug(ET, 'ET')

    return ET.rename(['et'])
