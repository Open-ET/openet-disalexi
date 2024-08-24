import math

import ee

from . import tseb_utils
from . import utils

def tseb_pt(
        t_air, t_rad, t_air0, e_air, u, p, z, rs_1, rs24, vza,
        aleafv, aleafn, aleafl, adeadv, adeadn, adeadl,
        albedo, ndvi, lai, clump, leaf_width, hc_min, hc_max,
        datetime, lon=None, lat=None, a_pt_in=1.32,
        stabil_iter=None, albedo_iter=10, et_min=0.01
):
    """Priestley-Taylor TSEB

    Calculates the Priestley Taylor TSEB fluxes using a single observation of
    composite radiometric temperature and using resistances in series.

    Parameters
    ----------
    t_air : ee.Image
        Air temperature [K].
    t_rad : ee.Image
        Radiometric composite temperature [K].
    t_air0 : ee.Image
        Measured Air Temperature [K]
    e_air : ee.Image
        Vapour pressure [kPa]
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
        Albedo
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
    stabil_iter : int, optional
        Number of iterations of stability calculation.  If not set the number
        of iterations will be computed dynamically.
    albedo_iter : int, optional
        Number of iterations of albedo separation calculation
        (the default is 10).
    et_min : float, optinal
        Minimum output ET value (the default is 0.01).

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
    mask = lai.double().multiply(0).rename(['mask'])

    # ************************************************************************
    # Apply met bands directly to Landsat image
    # CGM - This can probably be removed if Rs is resampled/smoothed in disalexi.py
    rs_1 = mask.add(rs_1).rename(['rs'])
    rs24 = mask.add(rs24).rename(['rs'])
    # t_air = mask.add(t_air).rename(['ta'])
    # u = mask.add(u).rename(['windspeed'])

    # ************************************************************************
    # time and t_noon could be moved into compute_G0 since that is the only place they are used
    time = ee.Date(datetime).get('hour').add(ee.Date(datetime).get('minute').divide(60))
    if lat is None:
        lat = ee.Image.pixelLonLat().select(['latitude'])
    if lon is None:
        lon = ee.Image.pixelLonLat().select(['longitude'])
    t_noon = tseb_utils.solar_noon(datetime=datetime, lon=lon.multiply(math.pi / 180))
    zs = tseb_utils.solar_zenith(
        datetime=datetime,
        lon=lon.multiply(math.pi / 180),
        lat=lat.multiply(math.pi / 180),
    )

    # ************************************************************************
    # Correct Clumping Factor
    f_green = 1.0

    # LAI for leaf spherical distribution
    F = lai.multiply(clump)

    # Fraction cover at nadir (view=0)
    # fc = 1.0 - exp(-0.5 * F)
    fc = F.multiply(-0.5).exp().multiply(-1).add(1.0).clamp(0.01, 0.9)

    # Compute canopy height and roughness parameters
    # hc = hc_min + ((hc_max - hc_min) * fc)
    hc = hc_max.subtract(hc_min).multiply(fc).add(hc_min)

    # LAI relative to canopy projection only
    lai_c = lai.divide(fc)

    # Houborg modification (according to Anderson et al. 2005)
    fc_q = (
        lai
        .expression('1 - (exp(-0.5 * F / cos(vza)))', {'F': F, 'vza': vza})
        .clamp(0.05, 0.90)
    )

    # Brutsaert (1982)
    z0m = hc.multiply(0.123)
    # CGM - add(0) is to mimic numpy copy, check if needed
    z0h = z0m.add(0)
    d0 = hc.multiply(2.0 / 3.0)

    # Correction of roughness parameters for water bodies
    # (NDVI < 0 and albedo < 0.05)
    water_mask = ndvi.lte(0).And(albedo.lte(0.05))
    d0 = d0.where(water_mask, 0.00001)
    z0m = z0m.where(water_mask, 0.00035)
    z0h = z0h.where(water_mask, 0.00035)

    # Check to avoid division by 0 in the next computations
    z0h = z0h.where(z0h.eq(0), 0.001)
    z0m = z0m.where(z0m.eq(0), 0.01)

    z_u = 30.0
    z_t = 30.0

    # Parameters for In-Canopy Wind Speed Extinction
    leaf = lai.expression(
        '(0.28 * (F ** (0.66667)) * (hc ** (0.33333)) * (leaf_width ** (-0.33333)))',
        {'F': F, 'hc': hc, 'leaf_width': leaf_width}
    )
    leaf_c = lai.expression(
        '(0.28 * (lai_c ** (0.66667)) * (hc ** (0.33333)) * (leaf_width ** (-0.33333)))',
        {'lai_c': lai_c, 'hc': hc, 'leaf_width': leaf_width}
    )
    leaf_s = lai.expression(
        '(0.28 * (0.1 ** (0.66667)) * (hc ** (0.33333)) * (leaf_width ** (-0.33333)))',
        {'hc': hc, 'leaf_width': leaf_width}
    )

    # ************************************************************************
    # Atmospheric Parameters
    # Saturation vapour pressure [kPa] (FAO56 3-8)
    #Yun modified to use METEO air temperature
    e_s0 = t_air0.expression(
        '0.6108 * exp((17.27 * (t_air - 273.16)) / ((t_air - 273.16) + 237.3))',
        {'t_air': t_air0}
    )
    vpd = e_s0.subtract(e_air)

    # Saturation vapor pressure [kpa] using iterated air temperature
    e_s = t_air.expression(
        '0.6108 * exp((17.27 * (t_air - 273.16)) / ((t_air - 273.16) + 237.3))',
        {'t_air': t_air}
    )

    # Slope of the saturation vapor pressure [kPa] (FAO56 3-9)
    Ss = t_air.expression(
        '4098. * e_s / (((t_air - 273.16) + 237.3) ** 2)',
        {'e_s': e_s, 't_air': t_air}
    )

    # Latent heat of vaporization (~2.45 at 20 C) [MJ kg-1] (FAO56 3-1)
    # lambda1 = (2.501 - (2.361e-3 * (t_air - 273.16)))
    lambda1 = t_air.subtract(273.16).multiply(2.361e-3).multiply(-1).add(2.501)

    # Psychrometric constant [kPa C-1] (FAO56 3-10)
    g = p.multiply(1.615E-3).divide(lambda1)

    # ************************************************************************
    a_pt = mask.add(a_pt_in)
    a_pt = a_pt.add(vpd.subtract(2.0).multiply(0.4)).max(a_pt_in).min(2.5).rename('a_pt')

    if stabil_iter is None:
        # Compute the number of stability iterations dynamically
        a_pt_max = ee.Number(
            a_pt.reduceRegion(reducer=ee.Reducer.max(), scale=4000, maxPixels=1E10)
            .get('a_pt')
        )
        stabil_iter = a_pt_max.divide(0.05).ceil().max(25).min(40)

    # ************************************************************************
    Rs_c, Rs_s, albedo_c, albedo_s, taudl, tausolar = tseb_utils.albedo_separation(
        albedo, rs_1, F, fc, aleafv, aleafn, aleafl, adeadv, adeadn, adeadl, zs, albedo_iter
    )

    # CGM - Moved emissivity calculation to separate function.
    #   I removed the Rs0 check.
    e_atm = tseb_utils.emissivity(t_air)
    # p = t_air.expression(
    #     '101.3 * (((t_air - (0.0065 * z)) / t_air) ** 5.26)', {'t_air': t_air, 'z': z}
    # )
    # Density of air? (kg m-3)
    r_air = t_air.expression(
        '101.3 * (((t_air - (0.0065 * z)) / t_air) ** 5.26) / 1.01 / t_air / 0.287',
        {'t_air': t_air, 'z': z}
    )
    cp = 1004.16

    # Assume neutral conditions on first iteration (use t_air for Ts and Tc)
    u_attr = tseb_utils.compute_u_attr(u=u, d0=d0, z0m=z0m, z_u=z_u, fm=0)
    r_ah = tseb_utils.compute_r_ah(u_attr=u_attr, d0=d0, z0h=z0h, z_t=z_t, fh=0)
    # CGM - Why is this function passing "lai" to "F"?
    r_s = tseb_utils.compute_r_s(
        u_attr=u_attr, T_s=t_air, T_c=t_air, hc=hc, F=lai, d0=d0, z0m=z0m,
        leaf=leaf, leaf_s=leaf_s, fm_h=0
    )
    r_x = tseb_utils.compute_r_x(
        u_attr=u_attr, hc=hc, F=lai, d0=d0, z0m=z0m, xl=leaf_width, leaf_c=leaf_c, fm_h=0
    )

    T_c = t_air.multiply(1)
    T_s = lai.expression(
        '(((t_rad - 273.16) - (fc_q * (T_c - 273.16))) / (1 - fc_q)) + 273.16',
        {'t_rad': t_rad, 'T_c': T_c, 'fc_q': fc_q}
    )
    # T_s = lai.expression(
    #     '(t_rad - (fc_q * T_c)) / (1 - fc_q)',
    #     {'t_rad': t_rad, 'T_c': T_c, 'fc_q': fc_q}
    # )

    # CGM - This doesn't seem to do anything, commenting out for now
    # H_iter = t_air.multiply(0).add(200.16)

    # CGM - Initialize to match t_air shape
    EF_s = t_air.multiply(0)

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

        Rn_c = tseb_utils.compute_Rn_c(albedo_c, t_air, T_c_iter, T_s_iter, e_atm, Rs_c, F)
        Rn_s = tseb_utils.compute_Rn_s(albedo_s, t_air, T_c_iter, T_s_iter, e_atm, Rs_s, F)
        Rn = Rn_c.add(Rn_s)

        G = tseb_utils.compute_G0(Rn, Rn_s, albedo, ndvi, t_noon, time, EF_s_iter)

        LE_c_init = (
            albedo
            .expression(
                'f_green * (a_pt * Ss / (Ss + g)) * Rn_c',
                {'f_green': f_green, 'a_pt': a_pt_iter, 'Ss': Ss, 'g': g, 'Rn_c': Rn_c}
            )
            .max(0)
        )

        H_c_init = Rn_c.subtract(LE_c_init)

        T_c_iter = tseb_utils.temp_separation_tc(
            H_c_init, fc_q, t_air, t_rad, r_ah_iter, r_s_iter, r_x_iter, r_air, cp
        )
        T_s_iter = tseb_utils.temp_separation_ts(T_c_iter, fc_q, t_air, t_rad)
        T_ac = tseb_utils.temp_separation_tac(
            T_c_iter, T_s_iter, fc_q, t_air, r_ah_iter, r_s_iter, r_x_iter
        )

        H_s = albedo.expression(
            'r_air * cp * (T_s - T_ac) / r_s',
            {'r_air': r_air, 'cp': cp, 'T_s': T_s_iter, 'T_ac': T_ac, 'r_s': r_s_iter}
        )
        H_c = albedo.expression(
            'r_air * cp * (T_c - T_ac) / r_x',
            {'r_air': r_air, 'cp': cp, 'T_c': T_c_iter, 'T_ac': T_ac, 'r_x': r_x_iter}
        )
        H = albedo.expression('H_s + H_c', {'H_s': H_s, 'H_c': H_c})
        H = H.where(H.eq(0), 10.0)

        LE_s = Rn_s.subtract(G).subtract(H_s)
        LE_c = Rn_c.subtract(H_c)
        # LE_s = albedo.expression('Rn_s - G - H_s', {'Rn_s': Rn_s, 'G': G, 'H_s': H_s})
        # LE_c = albedo.expression('Rn_c - H_c', {'Rn_c': Rn_c, 'H_c': H_c})

        # CGM - This won't do anything at this position in the code.
        #   Commenting out for now.
        # r_ah_iter = r_ah_iter.where(r_ah_iter.eq(0), 10.0)

        # CGM - This doesn't seem to do anything, commenting out for now
        # mask_iter = H_iter.divide(H).lte(1.05).And(H_iter.divide(H).gte(0.95))
        # chk_iter = np.sum(mask_iter) / np.size(mask_iter)

        fh = tseb_utils.compute_stability_fh(H, t_rad, u_attr_iter, r_air, z_t, d0, cp)
        fm = tseb_utils.compute_stability_fm(H, t_rad, u_attr_iter, r_air, z_u, d0, z0m, cp)
        fm_h = tseb_utils.compute_stability_fm_h(H, t_rad, u_attr_iter, r_air, hc, d0, z0m, cp)

        u_attr_iter = tseb_utils.compute_u_attr(u=u, d0=d0, z0m=z0m, z_u=z_u, fm=fm)
        r_ah_iter = tseb_utils.compute_r_ah(u_attr=u_attr_iter, d0=d0, z0h=z0h, z_t=z_t, fh=fh)
        r_s_iter = tseb_utils.compute_r_s(
            u_attr=u_attr_iter, T_s=T_s_iter, T_c=T_c_iter, hc=hc, F=lai,
            d0=d0, z0m=z0m, leaf=leaf, leaf_s=leaf_s, fm_h=fm_h
        )
        # CGM - Why is this function is passing "lai" to "F"?
        r_x_iter = tseb_utils.compute_r_x(
            u_attr=u_attr_iter, hc=hc, F=lai, d0=d0, z0m=z0m,
            xl=leaf_width, leaf_c=leaf_c, fm_h=fm_h
        )

        a_pt_iter = (
            a_pt_iter
            .where(LE_s.lte(0), a_pt_iter.subtract(0.05))
            .where(a_pt_iter.lte(0), 0.01)
        )

        den_s = Rn_s.subtract(G)
        den_s = den_s.updateMask(den_s.neq(0))
        # den_s[den_s == 0.] = np.nan

        EF_s_iter = LE_s.divide(den_s)

        return ee.Dictionary({
            'a_pt': a_pt_iter, 'EF_s': EF_s_iter, 'G': G,
            'H_c': H_c, 'H_s': H_s, 'LE_c': LE_c, 'LE_s': LE_s,
            'Rn_c': Rn_c, 'Rn_s': Rn_s,
            'r_ah': r_ah_iter, 'r_s': r_s_iter, 'r_x': r_x_iter,
            'T_ac': T_ac, 'T_c': T_c_iter, 'T_s': T_s_iter,
            'u_attr': u_attr_iter,
        })

    # Iterate the function n times
    # CGM - Iteration count is an input to the function
    input_images = ee.Dictionary({
        'a_pt': a_pt, 'EF_s': EF_s, 'G': ee.Image(0),
        'H_c': ee.Image(0), 'H_s': ee.Image(0),
        'LE_c': ee.Image(0), 'LE_s': ee.Image(0),
        'Rn_c': ee.Image(0), 'Rn_s': ee.Image(0),
        'r_ah': r_ah, 'r_s': r_s, 'r_x': r_x,
        'T_ac': ee.Image(0), 'T_c': T_c, 'T_s': T_s,
        'u_attr': u_attr,
    })
    iter_output = ee.Dictionary(
        ee.List.sequence(1, stabil_iter).iterate(iter_func, input_images)
    )

    # Unpack the iteration output
    a_pt = ee.Image(iter_output.get('a_pt'))
    Rn_c = ee.Image(iter_output.get('Rn_c'))
    Rn_s = ee.Image(iter_output.get('Rn_s'))
    G = ee.Image(iter_output.get('G'))
    H_c = ee.Image(iter_output.get('H_c'))
    H_s = ee.Image(iter_output.get('H_s'))
    LE_c = ee.Image(iter_output.get('LE_c'))
    LE_s = ee.Image(iter_output.get('LE_s'))

    # ************************************************************************
    # Check Energy Balance Closure
    ind = a_pt.lte(0.01)
    LE_s = LE_s.where(ind, 1.0)
    LE_c = LE_c.where(ind, 1.0)
    G = G.where(ind, Rn_s.subtract(H_s))

    ind = LE_s.gt(Rn_s)
    LE_s = LE_s.where(ind,  Rn_s)
    # H_s = H_s.where(ind,  Rn_s.subtract(G).subtract(LE_s))
    H_s = Rn_s.subtract(G).subtract(LE_s)
    # CGM - Check order of operations
    ind = LE_c.gt(Rn_c.add(100))
    # CGM - Not used below since LE_c is recomputed
    LE_c = LE_c.where(ind, Rn_c.add(100))
    # H_c = H_c.where(ind, -100)
    H_c = Rn_c.subtract(LE_c)

    #LE_s = Rn_s.subtract(G).subtract(H_s)
    #LE_c = Rn_c.subtract(H_c)

    # The latent heat of vaporization is 2.45 MJ kg-1
    # Assume rs24 is still in W m-2 day-1 and convert to MJ kg-1
    # CGM - Leaving out scaling value for now
    ET = (
        albedo
        .expression(
            '((LE_c + LE_s) / rs_1) * (rs24 / 2.45) * scaling',
            {
                'LE_c': LE_c, 'LE_s': LE_s, 'rs_1': rs_1,
                'rs24': rs24.multiply(0.0864 / 24.0), 'scaling': 1,
            }
        )
        .max(et_min)
    )

    return ET.rename(['et'])


def tseb_invert(
        et_alexi, t_air, t_rad, t_air0, e_air, u, p, z, rs_1, rs24, vza,
        aleafv, aleafn, aleafl, adeadv, adeadn, adeadl,
        albedo, ndvi, lai, clump, leaf_width, hc_min, hc_max,
        datetime, lon=None, lat=None, a_pt_in=1.32,
        stabil_iter=None, albedo_iter=10,
):
    """Priestley-Taylor TSEB

    Calculates the Priestley Taylor TSEB fluxes using a single observation of
    composite radiometric temperature and using resistances in series.

    Parameters
    ----------
    et_alexi : ee.Image
        ALEXI scale ET [mm].
    t_air : ee.Image
        Iterated air temperature [K].
    t_rad : ee.Image
        Radiometric composite temperature [K].
    t_air0 : ee.Image
        Measured air Temperature [K]
    e_air : ee.Image
        Vapour pressure [kPa]
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
    stabil_iter : int, optional
        Number of iterations of stability calculation.  If not set the number
        of iterations will be computed dynamically.
    albedo_iter : int, optional
        Number of iterations of albedo separation calculation
        (the default is 10).

    Returns
    -------
    ta : ee.Image
        Air temperature [K?].

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

    mask = lai.double().multiply(0).rename(['mask'])

    # ************************************************************************
    # Apply met bands directly to Landsat image
    # CGM - This can probably be removed if Rs is resampled/smoothed in disalexi.py
    rs_1 = mask.add(rs_1).rename(['rs'])
    rs24 = mask.add(rs24).rename(['rs'])
    # t_air = mask.add(t_air).rename(['ta'])
    # u = mask.add(u).rename(['windspeed'])

    # Estimate total daily latent heat (not in tseb_pt())
    LEtot = et_alexi.expression(
        '(et_alexi * 2.45 / (Rsd * 0.0864 / 24)) * Rs',
        {'et_alexi': et_alexi, 'Rsd': rs24, 'Rs': rs_1}
    )

    # ************************************************************************
    # time and t_noon could be moved into compute_G0 since that is the only place they are used
    time = ee.Date(datetime).get('hour').add(ee.Date(datetime).get('minute').divide(60))
    if lat is None:
        lat = ee.Image.pixelLonLat().select(['latitude'])
    if lon is None:
        lon = ee.Image.pixelLonLat().select(['longitude'])
    t_noon = tseb_utils.solar_noon(datetime=datetime, lon=lon.multiply(math.pi / 180))
    zs = tseb_utils.solar_zenith(
        datetime=datetime,
        lon=lon.multiply(math.pi / 180),
        lat=lat.multiply(math.pi / 180),
    )

    # ************************************************************************
    # Correct Clumping Factor
    f_green = 1.0

    # LAI for leaf spherical distribution
    F = lai.multiply(clump)

    # Fraction cover at nadir (view=0)
    # fc = 1.0 - exp(-0.5 * F)
    fc = F.multiply(-0.5).exp().multiply(-1).add(1.0).clamp(0.01, 0.9)

    # Compute canopy height and roughness parameters
    # hc = hc_min + ((hc_max - hc_min) * fc)
    hc = hc_max.subtract(hc_min).multiply(fc).add(hc_min)

    # LAI relative to canopy projection only
    lai_c = lai.divide(fc)

    # Houborg modification (according to Anderson et al. 2005)
    fc_q = (
        et_alexi
        .expression('1 - (exp(-0.5 * F / cos(vza)))', {'F': F, 'vza': vza})
        .clamp(0.05, 0.90)
    )

    # Brutsaert (1982)
    z0m = hc.multiply(0.123)
    # CGM - add(0) is to mimic numpy copy, check if needed
    z0h = z0m.add(0)
    d0 = hc.multiply(2.0 / 3.0)

    # Correction of roughness parameters for water bodies
    # (NDVI < 0 and albedo < 0.05)
    water_mask = ndvi.lte(0).And(albedo.lte(0.05))
    d0 = d0.where(water_mask, 0.00001)
    z0m = z0m.where(water_mask, 0.00035)
    z0h = z0h.where(water_mask, 0.00035)

    # Check to avoid division by 0 in the next computations
    z0h = z0h.where(z0h.eq(0), 0.001)
    z0m = z0m.where(z0m.eq(0), 0.01)

    z_u = 30.0
    z_t = 30.0

    # Modify z0m when using it at alexi scale (not in tseb_pt())
    z0m = et_alexi.expression(
        '1.0 / ((log((z_U - D0) / z0m)) * (log((z_U - D0) / z0m)))',
        {'z_U': z_u, 'D0': d0, 'z0m': z0m}
    )
    # Further for z0m
    z0m = et_alexi.expression(
        '(z_U - D0) / (exp(1.0 / sqrt(z0m)))', {'z_U': z_u, 'D0': d0, 'z0m': z0m}
    )

    # Redefine z0h
    z0h = z0m.add(0)

    # Parameters for In-Canopy Wind Speed Extinction
    leaf = et_alexi.expression(
        '(0.28 * (F ** (0.66667)) * (hc ** (0.33333)) * (leaf_width ** (-0.33333)))',
        {'F': F, 'hc': hc, 'leaf_width': leaf_width}
    )
    leaf_c = et_alexi.expression(
        '(0.28 * (lai_c ** (0.66667)) * (hc ** (0.33333)) * (leaf_width ** (-0.33333)))',
        {'lai_c': lai_c, 'hc': hc, 'leaf_width': leaf_width}
    )
    leaf_s = et_alexi.expression(
        '(0.28 * (0.1 ** (0.66667)) * (hc ** (0.33333)) * (leaf_width ** (-0.33333)))',
        {'hc': hc, 'leaf_width': leaf_width}
    )

    # ************************************************************************
    # Atmospheric Parameters
    # Saturation vapour pressure [kPa] (FAO56 3-8)
    # Yun modified to use METEO air temperature
    e_s0 = t_air0.expression(
        '0.6108 * exp((17.27 * (t_air - 273.16)) / ((t_air - 273.16) + 237.3))',
        {'t_air': t_air0}
    )
    vpd = e_s0.subtract(e_air)

    # Saturation vapor pressure [kpa] using iterated air temperature
    e_s = t_air.expression(
        '0.6108 * exp((17.27 * (t_air - 273.16)) / ((t_air - 273.16) + 237.3))',
        {'t_air': t_air}
    )

    # Slope of the saturation vapor pressure [kPa] (FAO56 3-9)
    Ss = t_air.expression(
        '4098. * e_s / (((t_air - 273.16) + 237.3) ** 2)',
        {'e_s': e_s, 't_air': t_air}
    )

    # Latent heat of vaporization (~2.45 at 20 C) [MJ kg-1] (FAO56 3-1)
    # lambda1 = (2.501 - (2.361e-3 * (t_air - 273.16)))
    lambda1 = t_air.subtract(273.16).multiply(2.361e-3).multiply(-1).add(2.501)

    # Psychrometric constant [kPa C-1] (FAO56 3-10)
    g = p.multiply(1.615E-3).divide(lambda1)

    # ************************************************************************
    a_pt = mask.add(a_pt_in)
    a_pt = (
        a_pt.add(vpd.subtract(2.0).multiply(0.4))
        .max(a_pt_in).min(2.5)
        .rename('a_pt')
    )

    if stabil_iter is None:
        # Compute the number of stability iterations dynamically
        a_pt_max = ee.Number(
            a_pt.reduceRegion(reducer=ee.Reducer.max(), scale=4000, maxPixels=1E10)
            .get('a_pt')
        )
        stabil_iter = a_pt_max.divide(0.05).ceil().max(25).min(40)

    # ************************************************************************
    Rs_c, Rs_s, albedo_c, albedo_s, taudl, tausolar = tseb_utils.albedo_separation(
        albedo, rs_1, F, fc, aleafv, aleafn, aleafl, adeadv, adeadn, adeadl, zs, albedo_iter
    )

    # CGM - Moved emissivity calculation to separate function.
    #   I removed the Rs0 check.
    e_atm = tseb_utils.emissivity(t_air)
    # p = t_air.expression(
    #     '101.3 * (((t_air - (0.0065 * z)) / t_air) ** 5.26)',
    #     {'t_air': t_air, 'z': z}
    # )

    # Density of air? (kg m-3)
    r_air = t_air.expression(
        '101.3 * (((t_air - (0.0065 * z)) / t_air) ** 5.26) / 1.01 / t_air / 0.287',
        {'t_air': t_air, 'z': z}
    )
    cp = 1004.16

    # Assume neutral conditions on first iteration (use t_air for Ts and Tc)
    u_attr = tseb_utils.compute_u_attr(u=u, d0=d0, z0m=z0m, z_u=z_u, fm=0)
    r_ah = tseb_utils.compute_r_ah(u_attr=u_attr, d0=d0, z0h=z0h, z_t=z_t, fh=0)
    # CGM - Why is this function passing "lai" to "F"?
    r_s = tseb_utils.compute_r_s(
        u_attr=u_attr, T_s=t_air, T_c=t_air, hc=hc, F=lai, d0=d0, z0m=z0m,
        leaf=leaf, leaf_s=leaf_s, fm_h=0
    )
    r_x = tseb_utils.compute_r_x(
        u_attr=u_attr, hc=hc, F=lai, d0=d0, z0m=z0m, xl=leaf_width, leaf_c=leaf_c, fm_h=0
    )

    T_c = t_air.multiply(1)
    # Modified from tseb_pt() approach
    T_s = et_alexi.expression(
        '(((t_rad - 273.16) - (fc_q * (T_c - 273.16))) / (1 - fc_q)) + 273.16',
        {'t_rad': t_rad, 'T_c': T_c, 'fc_q': fc_q}
    )

    # CGM - This doesn't seem to do anything, commenting out for now
    # H_iter = t_air.multiply(0).add(200.16)

    # CGM - Initialize to match t_air shape
    EF_s = t_air.multiply(0)

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
        ta_iter = ee.Image(ee.Dictionary(prev).get('ta'))
        #ta_iter=t_air # test tair not change

        # Yun whether tair changes in each iteration???
        Rn_c = tseb_utils.compute_Rn_c(albedo_c, ta_iter, T_c_iter, T_s_iter, e_atm, Rs_c, F)
        Rn_s = tseb_utils.compute_Rn_s(albedo_s, ta_iter, T_c_iter, T_s_iter, e_atm, Rs_s, F)
        Rn = Rn_c.add(Rn_s)

        G = tseb_utils.compute_G0(Rn, Rn_s, albedo, ndvi, t_noon, time, EF_s_iter)

        # Yun - use the ALEXI ET to get daily LE and then to get H
        H = Rn.subtract(LEtot).subtract(G)

        LE_c_init = (
            et_alexi
            .expression(
                'f_green * (a_pt * Ss / (Ss + g)) * Rn_c',
                {'f_green': f_green, 'a_pt': a_pt_iter, 'Ss': Ss, 'g': g, 'Rn_c': Rn_c})
            .max(0)
        )
        H_c_init = Rn_c.subtract(LE_c_init)
        H_s = H.subtract(H_c_init)

        LE_s = Rn_s.subtract(G).subtract(H_s)

        # CGM - This won't do anything at this position in the code. Commenting out for now.
        # r_ah_iter = r_ah_iter.where(r_ah_iter.eq(0), 10.0)

        # CGM - This doesn't seem to do anything, commenting out for now
        # mask_iter = H_iter.divide(H).lte(1.05).And(H_iter.divide(H).gte(0.95))
        # chk_iter = np.sum(mask_iter) / np.size(mask_iter)

        fh = tseb_utils.compute_stability_fh(H, t_rad, u_attr_iter, r_air, z_t, d0, cp)
        fm = tseb_utils.compute_stability_fm(H, t_rad, u_attr_iter, r_air, z_u, d0, z0m, cp)
        fm_h = tseb_utils.compute_stability_fm_h(H, t_rad, u_attr_iter, r_air, hc, d0, z0m, cp)

        u_attr_iter = tseb_utils.compute_u_attr(u=u, d0=d0, z0m=z0m, z_u=z_u, fm=fm)
        r_ah_iter = tseb_utils.compute_r_ah(u_attr=u_attr_iter, d0=d0, z0h=z0h, z_t=z_t, fh=fh)
        r_s_iter = tseb_utils.compute_r_s(
            u_attr=u_attr_iter, T_s=T_s_iter, T_c=T_c_iter, hc=hc, F=lai,
            d0=d0, z0m=z0m, leaf=leaf, leaf_s=leaf_s, fm_h=fm_h
        )
        # CGM - Why is this function passing "lai" to "F"?
        r_x_iter = tseb_utils.compute_r_x(
            u_attr=u_attr_iter, hc=hc, F=lai, d0=d0, z0m=z0m,
            xl=leaf_width, leaf_c=leaf_c, fm_h=fm_h
        )

        T_c_iter = tseb_utils.temp_separation_tc(
            H_c_init, fc_q, ta_iter, t_rad, r_ah_iter, r_s_iter, r_x_iter, r_air, cp
        )
        T_s_iter = tseb_utils.temp_separation_ts(T_c_iter, fc_q, ta_iter, t_rad)
        T_ac = tseb_utils.temp_separation_tac(
            T_c_iter, T_s_iter, fc_q, ta_iter, r_ah_iter, r_s_iter, r_x_iter
        )
        ta_iter = T_ac.subtract(H.multiply(r_ah_iter).divide(r_air.multiply(cp)))

        a_pt_iter = (
            a_pt_iter
            .where(LE_s.lte(0), a_pt_iter.subtract(0.05))
            .where(a_pt_iter.lte(0), 0.01)
        )

        den_s = Rn_s.subtract(G)
        den_s = den_s.updateMask(den_s.neq(0))
        # den_s[den_s == 0.] = np.nan

        EF_s_iter = LE_s.divide(den_s)

        return ee.Dictionary({
            'a_pt': a_pt_iter, 'EF_s': EF_s_iter, 'G': G,
            'Rn_c': Rn_c, 'Rn_s': Rn_s,
            'r_ah': r_ah_iter, 'r_s': r_s_iter, 'r_x': r_x_iter,
            'T_ac': T_ac, 'T_c': T_c_iter, 'T_s': T_s_iter,
            'u_attr': u_attr_iter, 'ta': ta_iter,
        })

    # Iterate the function n times
    # CGM - Iteration count is an input to the function
    input_images = ee.Dictionary({
        'a_pt': a_pt, 'EF_s': EF_s, 'G': ee.Image(0),
        'Rn_c': ee.Image(0), 'Rn_s': ee.Image(0),
        'r_ah': r_ah, 'r_s': r_s, 'r_x': r_x, 'rs_1': rs_1,
        'T_ac': ee.Image(0), 'T_c': T_c, 'T_s': T_s,
        'u_attr': u_attr, 'ta': t_air,
    })
    iter_output = ee.Dictionary(
        ee.List.sequence(1, stabil_iter).iterate(iter_func, input_images)
    )

    # Extract the Ta from the iteration output
    ta = ee.Image(iter_output.get('ta'))
    ta = ta.updateMask(ta.gte(173.16))
    ta = ta.updateMask(ta.lte(473.16))

    # # Uncomment to extract the other iteration variables
    # a_pt = ee.Image(iter_output.get('a_pt'))
    # EF_s = ee.Image(iter_output.get('EF_s'))
    # G = ee.Image(iter_output.get('G'))
    # Rn_c = ee.Image(iter_output.get('Rn_c'))
    # Rn_s = ee.Image(iter_output.get('Rn_s'))
    # r_ah = ee.Image(iter_output.get('r_ah'))
    # r_s = ee.Image(iter_output.get('r_s'))
    # r_x = ee.Image(iter_output.get('r_x'))
    # T_ac = ee.Image(iter_output.get('T_ac'))
    # T_c = ee.Image(iter_output.get('T_c'))
    # T_s = ee.Image(iter_output.get('T_s'))
    # u_attr = ee.Image(iter_output.get('u_attr'))

    return ta.rename(['ta'])
