import math

import ee

deg2rad = math.pi / 180.0
rad2deg = 180.0 / math.pi


def solar_noon(datetime, lon):
    """Computes sunrise/sunset times

    Parameters
    ----------
    datetime : ee.Date
    lon : ee.Image
        Longitude (radians)
    lat : ee.Image
        Latitude (radians)

    Returns
    -------
    t_noon : ee.Image

    """
    # Adjust image datetime to start of day
    d, eq_t = _solar_time(ee.Date(datetime.format('yyyy-MM-dd')))

    # TODO: Check if the order of operations is correct
    #   Should the 1440 * 24.0 be in parenthesis?
    # t_noon = (720.0 - 4 * lon - eq_t) / 1440 * 24.0
    # t_noon = ((720.0 - (4 * lon) - eq_t) / 1440) * 24.0
    t_noon = lon.multiply(rad2deg).multiply(-4).subtract(eq_t).add(720).divide(1440).multiply(24.0)

    return t_noon.rename(['t_noon'])


def solar_zenith(datetime, lon, lat):
    """Computes zenith angle

    Parameters
    ----------
    datetime : ee.Date
    lon : ee.Image
        Longitude (radians)
    lat : ee.Image
        Latitude (radians)

    Returns
    -------
    zs : ee.Image

    """
    # IDL is computing time_t as hours and fractional minutes (no seconds)
    time_t = ee.Date(datetime).get('hour').add(ee.Date(datetime).get('minute').divide(60))
    # # This will return the hour floating point value
    # time_t = ee.Date(date).get('hour').add(ee.Date(date).getFraction('hour'))

    d, eq_t = _solar_time(datetime)

    ts_time = lon.expression(
        '((time_t / 24.0) * 1440.0 + eq_t + 4.0 * lon) % 1440.0',
        {'lon': lon.multiply(rad2deg), 'time_t': time_t, 'eq_t': eq_t}
    )
    ts_time = ts_time.where(ts_time.gt(1440), ts_time.subtract(1440))

    w = ts_time.divide(4).add(180)
    w = w.where(ts_time.divide(4).gt(0), ts_time.divide(4).subtract(180))

    zs = lat.expression(
        'acos( (sin(lat) * sin(d)) + (cos(lat) * cos(d) * cos(w)) )',
        {'lat': lat, 'd': d, 'w': w.multiply(deg2rad)}
    )

    return zs.rename(['zs'])


def _solar_time(datetime):
    """Computes solar time variables following Campbell & Norman 1998

    Parameters
    ----------
    datetime : ee.Date
    lon : ee.Image
        Longitude (radians)
    lat : ee.Image
        Latitude (radians)

    Returns
    -------
    eq_t : ee.Number

    """
    # IDL is computing time_t as hours and fractional minutes
    time_t = ee.Date(datetime).get('hour').add(ee.Date(datetime).get('minute').divide(60))
    # This will return the hour floating point value
    # time_t = ee.Date(date).get('hour').add(ee.Date(date).getFraction('hour'))

    # CGM - DOY and hour could be images in order to make these expressions
    julian = _to_jd(datetime)

    # Sunrise time
    julian_ = time_t.divide(24.0).add(julian)
    j_cen = julian_.add(0.5 - 2451545.0).divide(36525.0)
    # CGM - Does the mod happen before or after the multiply
    lon_sun = (
        j_cen.multiply(0.0003032).add(36000.76983)
        .multiply(j_cen).mod(360.0).add(280.46646).subtract(360.0)
    )
    an_sun = j_cen.multiply(-0.0001537).add(35999.05029).multiply(j_cen).add(357.52911)
    ecc = (
        j_cen.multiply(0.0000001267).add(0.000042037)
        .multiply(j_cen).multiply(-1).add(0.016708634)
    )
    ob_ecl = (
        j_cen.multiply(-0.001813).add(0.00059)
        .multiply(j_cen).add(46.815)
        .multiply(j_cen).multiply(-1).add(21.448)
        .divide(60.0).add(26).divide(60).add(23)
    )
    ob_corr = (
        j_cen.multiply(-1934.136).add(125.04).multiply(deg2rad).cos()
        .multiply(0.00256).add(ob_ecl)
    )
    var_y = (
        ob_corr.divide(2.0).multiply(deg2rad).tan()
        .multiply(ob_corr.divide(2.0).multiply(deg2rad).tan())
    )
    eq_t = (
        lon_sun.multiply(2.0).multiply(deg2rad).sin().multiply(var_y)
        .subtract(an_sun.multiply(deg2rad).sin().multiply(ecc).multiply(2.0))
        .add(
            an_sun.multiply(deg2rad).sin()
            .multiply(lon_sun.multiply(2.0).multiply(deg2rad).cos())
            .multiply(var_y).multiply(ecc).multiply(4.0)
        )
        .subtract(
            lon_sun.multiply(4.0).multiply(deg2rad).sin()
            .multiply(var_y).multiply(var_y).multiply(0.5)
        )
        .subtract(
            an_sun.multiply(2.0).multiply(deg2rad).sin()
            .multiply(ecc).multiply(ecc).multiply(1.25)
        )
        .multiply(4.0).multiply(rad2deg)
    )
    sun_eq = (
        an_sun.multiply(deg2rad).sin()
        .multiply(
            j_cen.multiply(0.000014).add(0.004817)
            .multiply(j_cen).multiply(-1).add(1.914602)
        )
        .add(
            an_sun.multiply(2.0).multiply(deg2rad).sin()
            .multiply(j_cen.multiply(-0.000101).add(0.019993))
        )
        .add(an_sun.multiply(3.0).multiply(deg2rad).sin().multiply(0.000289))
    )
    sun_true = sun_eq.add(lon_sun)
    sun_app = (
        j_cen.multiply(-1934.136).add(125.04).multiply(deg2rad).sin()
        .multiply(-0.00478).subtract(0.00569).add(sun_true)
    )

    # CGM - Intentionally not converting back to degrees
    d = ob_corr.multiply(deg2rad).sin().multiply(sun_app.multiply(deg2rad).sin()).asin()

    return d, eq_t


def _to_jd(date):
    """Computes the Julian Day

    Follows equations in https://en.wikipedia.org/wiki/Julian_day
    IDL code is only computing Julian day from year and DOY which is equivalent
        to the JDN term.
    Python code is performing integer division (assuming Python 2.7)
        in JD calculation when it should be floating point division.

    Parameters
    ----------
    date : ee.Date

    Returns
    -------
    julian day : float

    """
    a = ee.Date(date).get('month').multiply(-1).add(14).divide(12).floor()
    y = ee.Date(date).get('year').add(4800).subtract(a)
    m = ee.Date(date).get('month').add(a.multiply(12)).subtract(3)

    jdn = (
        ee.Date(date).get('day')
        .add(m.multiply(153).add(2).divide(5).floor())
        .add(y.multiply(365))
        .add(y.divide(4.0).floor())
        .subtract(y.divide(100.0).floor())
        .add(y.divide(400.0).floor())
        .subtract(32045)
    )

    # To match IDL, only return JDN term
    return jdn

    # CGM - Forcing to float types to match wikipedia equations.
    #   Original equation was doing integer division in Python 2.7.
    #   Microseconds are not in the original equation.
    # jd = (
    #     jdn
    #     .add(ee.Date(date).get('hour').float().subtract(12).divide(24))
    #     .add(ee.Date(date).get('minute').float().divide(1440))
    #     .add(ee.Date(date).get('second').float().divide(86400))
    #     # .add(date.get('microsecond') / 86400000000)
    # )
    # return jd


def emissivity(t_air):
    """Apparent atmospheric emissivity

    Apparent emissivity (Sedlar and Hock, 2009: Cryosphere 3:75-84)
    Atmospheric emissivity (clear-sky) Idso and Jackson (1969)

    Parameters
    ----------
    t_air : ee.Image
        Air temperature (Kelvin)

    Returns
    -------
    e_atm : ee.Image

    """
    e_atm = t_air.subtract(273.16).pow(2.0).multiply(-0.0003523).exp().multiply(-0.2811).add(1)

    return e_atm


def albedo_separation(
        albedo, rs_1, f, fc, aleafv, aleafn, aleafl, adeadv, adeadn, adeadl, zs, iterations=10
):
    """Compute Solar Components and atmospheric properties ([Campbell1998]_)

    Parameters
    ----------
    albedo : ee.Image
    rs_1 : ee.Image
    f : ee.Image
    fc : ee.Image
    aleafv : ee.Image
    aleafn : ee.Image
    aleafl : ee.Image
    adeadv : ee.Image
    adeadn : ee.Image
    adeadl : ee.Image
    zs : ee.Image
    iterations: int

    Returns
    -------
    Rs_c : ee.Image
    Rs_s : ee.Image
    albedo_c : ee.Image
    albedo_s : ee.Image
    taudl : ee.Image
    tausolar : ee.Image

    References
    ----------
    .. [Campbell1998] Campbell, G. S. & Norman, J. M. (1998),
        An introduction to environmental biophysics. Springer, New York
        https://archive.org/details/AnIntroductionToEnvironmentalBiophysics.

    """
    rad2deg = 180.0 / math.pi

    # Correct for curvature of atmos in airmas
    airmas = zs.cos().pow(2).add(0.0025).sqrt().subtract(zs.cos()).divide(0.00125)

    # Correct for refraction(good up to 89.5 deg.)
    airmas = airmas.where(
        zs.multiply(rad2deg).lt(89.5),
        zs.multiply(rad2deg).multiply(-1).add(90).pow(-2).multiply(-2.8).add(airmas)
    )

    potbm1 = airmas.multiply(-0.160).exp().multiply(600)

    potvis = potbm1.multiply(-1).add(600).multiply(0.4).add(potbm1).multiply(zs.cos())

    # CGM - Applying log10 call to this variable instead of in a calculation below
    uu_log10 = zs.cos().pow(-1).max(0.01).log10()

    a_temp = uu_log10.multiply(-0.0345).add(0.4459).multiply(uu_log10).add(-1.195)
    a = zs.expression('10 ** a_temp', {'a_temp': a_temp})

    watabs = a.multiply(1320.0)

    potbm2 = airmas.multiply(-0.05).exp().multiply(720).subtract(watabs)

    # CGM - Reordered equation
    evaL = watabs.add(potbm2).subtract(720).multiply(-0.54).multiply(zs.cos())

    # potnir = evaL + potbm2 * cos(zs)
    potnir = zs.cos().multiply(potbm2).add(evaL)

    fclear = potvis.add(potnir).pow(-1).multiply(rs_1).clamp(0.01, 1.0).where(zs.cos().lte(0.01), 1)
  
    # Partition SDN into VIS and NIR
    fvis = potvis.add(potnir).pow(-1).multiply(potvis)
    fnir = potvis.add(potnir).pow(-1).multiply(potnir)
  
    # Estimate direct beam and diffuse fraction in VIS and NIR wavebands
    fb1 = zs.cos().multiply(potbm1).divide(potvis)
    # CGM - fb2 not used
    # fb2 = zs.cos().multiply(potbm2).divide(potnir)

    dirvis = (
        fclear.min(0.9).multiply(-1).add(0.9).divide(0.7).pow(0.6667).multiply(-1).add(1)
        .multiply(fb1).min(fb1)
    )
    dirnir = (
        fclear.min(0.88).multiply(-1).add(0.88).divide(0.68).pow(0.6667).multiply(-1).add(1)
        .multiply(fb1).min(fb1)
    )

    dirvis = dirvis.where(dirvis.lt(0.01).And(dirnir.gt(0.01)), 0.011)
    dirnir = dirnir.where(dirnir.lt(0.01).And(dirvis.gt(0.01)), 0.011)

    # difvis = 1.0 - dirvis
    # difnir = 1.0 - dirnir
    difvis = dirvis.multiply(-1).add(1)
    difnir = dirnir.multiply(-1).add(1)

    # Correction for NIGHTIME
    ind = zs.cos().lte(0.01)
    fvis = fvis.where(ind, 0.5)
    fnir = fnir.where(ind, 0.5)
    difvis = difvis.where(ind, 0.0)
    difnir = difnir.where(ind, 0.0)

    #**********************************************
    # Compute Albedo
    ratio_soil = 2.0

    # CGM - Initialize rsoilv and fg from f and albedo
    rsoilv = f.multiply(0).add(0.12)
    fg = albedo.multiply(0).add(1)

    # CGM - Switched to an iterate call
    def iter_func(n, prev):
        # Extract inputs from previous iteration
        fg_iter = ee.Image(ee.Dictionary(prev).get('fg'))
        rsoilv_iter = ee.Image(ee.Dictionary(prev).get('rsoilv'))
        # CGM - Variables that are commented out only need to be returned
        # akb = ee.Image(ee.Dictionary(prev).get('akb'))
        # albedo_c = ee.Image(ee.Dictionary(prev).get('albedo_c'))
        # albedo_s = ee.Image(ee.Dictionary(prev).get('albedo_s'))
        # ameann = ee.Image(ee.Dictionary(prev).get('ameann'))
        # ameanv = ee.Image(ee.Dictionary(prev).get('ameanv'))
        # diff = ee.Image(ee.Dictionary(prev).get('diff'))
        # rbcpyn = ee.Image(ee.Dictionary(prev).get('rbcpyn'))
        # rbcpyv = ee.Image(ee.Dictionary(prev).get('rbcpyv'))
        # taudn = ee.Image(ee.Dictionary(prev).get('taudn'))
        # taudv = ee.Image(ee.Dictionary(prev).get('taudv'))

        rsoiln = rsoilv_iter.multiply(ratio_soil)

        # Weighted live/dead leaf average properties
        ameanv = aleafv.expression(
            'aleafv * fg + adeadv * (1.0 - fg)',
            {'adeadv': adeadv, 'aleafv': aleafv, 'fg': fg_iter}
        )
        ameann = aleafn.expression(
            'aleafn * fg + adeadn * (1.0 - fg)',
            {'adeadn': adeadn, 'aleafn': aleafn, 'fg': fg_iter}
        )
        ameanl = aleafl.expression(
            'aleafl * fg + adeadl * (1.0 - fg)',
            {'adeadl': adeadl, 'aleafl': aleafl, 'fg': fg_iter}
        )

        # DIFFUSE COMPONENT
        #*******************************
        # Canopy reflection (deep canopy)
        # Fit to Fig 15.4 for x=1
        akd = f.log().multiply(-0.0683).add(0.804)

        # Eq 15.7
        rcpyn = ameann.expression('(1.0 - sqrt(ameann)) / (1.0 + sqrt(ameann))', {'ameann': ameann})
        rcpyv = ameanv.expression('(1.0 - sqrt(ameanv)) / (1.0 + sqrt(ameanv))', {'ameanv': ameanv})

        # Eq 15.8
        rdcpyn = akd.add(1).pow(-1).multiply(akd).multiply(rcpyn).multiply(2)
        rdcpyv = akd.add(1).pow(-1).multiply(akd).multiply(rcpyv).multiply(2)

        # Eq 15.11 - Canopy transmission (VIS)
        expfac = ameanv.sqrt().multiply(akd).multiply(f).max(0.001)
        xnum = f.expression('(rdcpyv * rdcpyv - 1.0) * exp(-expfac)', {'rdcpyv': rdcpyv, 'expfac': expfac})
        xden = f.expression(
            '(rdcpyv * rsoilv - 1.0) + rdcpyv * (rdcpyv - rsoilv) * exp(-2.0 * expfac)',
            {'expfac': expfac, 'rdcpyv': rdcpyv, 'rsoilv': rsoilv_iter}
        )
        taudv = xnum.divide(xden)

        # Eq 15.11 - Canopy transmission (NIR)
        expfac = ameann.sqrt().multiply(akd).multiply(f).max(0.001)
        xnum = f.expression('(rdcpyn * rdcpyn - 1.0) * exp(-expfac)', {'expfac': expfac, 'rdcpyn': rdcpyn})
        xden = f.expression(
            '(rdcpyn * rsoiln - 1.0) + rdcpyn * (rdcpyn - rsoiln) * exp(-2.0 * expfac)',
            {'expfac': expfac, 'rdcpyn': rdcpyn, 'rsoiln': rsoiln}
        )
        taudn = xnum.divide(xden)

        # Canopy transmission (LW)
        taudl = ameanl.sqrt().multiply(-1).multiply(akd).multiply(f)

        # Eq 15.9 - Diffuse albedo for generic canopy
        fact = f.expression(
            '((rdcpyn - rsoiln) / (rdcpyn * rsoiln - 1.0)) * exp(-2.0 * sqrt(ameann) * akd * f)',
            {'akd': akd, 'ameann': ameann, 'f': f, 'rdcpyn': rdcpyn, 'rsoiln': rsoiln}
        )
        albdn = f.expression('(rdcpyn + fact) / (1.0 + rdcpyn * fact)', {'fact': fact, 'rdcpyn': rdcpyn})

        # Eq 15.9
        fact = f.expression(
            '((rdcpyv - rsoilv) / (rdcpyv * rsoilv - 1.0)) * exp(-2.0 * sqrt(ameanv) * akd * f)',
            {'akd': akd, 'ameanv': ameanv, 'f': f, 'rdcpyv': rdcpyv, 'rsoilv': rsoilv_iter}
        )
        albdv = f.expression('(rdcpyv + fact) / (1.0 + rdcpyv * fact)', {'fact': fact, 'rdcpyv': rdcpyv})

        # BEAM COMPONENT
        #*******************************
        # Canopy reflection (deep canopy)
        akb = zs.cos().pow(-1).multiply(0.5).where(zs.cos().lte(0.01), 0.5)

        # Eq 15.7
        rcpyn = ameann.expression('(1.0 - sqrt(ameann)) / (1.0 + sqrt(ameann))', {'ameann': ameann})
        rcpyv = ameanv.expression('(1.0 - sqrt(ameanv)) / (1.0 + sqrt(ameanv))', {'ameanv': ameanv})

        # Eq 15.8
        rbcpyn = akb.add(1).pow(-1).multiply(rcpyn).multiply(akb).multiply(2.0)
        rbcpyv = akb.add(1).pow(-1).multiply(rcpyv).multiply(akb).multiply(2.0)

        # Eq 15.9 - Beam albedo for generic canopy
        fact = f.expression(
            '((rbcpyn - rsoiln) / (rbcpyn * rsoiln - 1.0)) * exp(-2.0 * sqrt(ameann) * akb * f)',
            {'akb': akb, 'ameann': ameann, 'f': f, 'rbcpyn': rbcpyn, 'rsoiln': rsoiln}
        )
        albbn = f.expression('(rbcpyn + fact) / (1.0 + rbcpyn * fact)', {'fact': fact, 'rbcpyn': rbcpyn})

        # Eq 15.9
        fact = f.expression(
            '((rbcpyv - rsoilv) / (rbcpyv * rsoilv - 1.0)) * exp(-2.0 * sqrt(ameanv) * akb * f)',
            {'akb': akb, 'ameanv': ameanv, 'f': f, 'rbcpyv': rbcpyv, 'rsoilv': rsoilv_iter}
        )
        albbv = f.expression('(rbcpyv + fact) / (1.0 + rbcpyv * fact)', {'fact': fact, 'rbcpyv': rbcpyv})

        # Weighted albedo (canopy)
        albedo_c = f.expression(
            'fvis * (dirvis * albbv + difvis * albdv) + fnir * (dirnir * albbn + difnir * albdn)',
            {
                'albbn': albbn, 'albbv': albbv, 'albdn': albdn, 'albdv': albdv,
                'difnir': difnir, 'difvis': difvis, 'dirvis': dirvis,
                'dirnir': dirnir, 'fnir': fnir, 'fvis': fvis,
            }
        )
        albedo_c = albedo_c.where(
            zs.cos().lte(0.01),
            f.expression(
                'fvis * (difvis * albdv) + fnir * (difnir * albdn)',
                {
                    'albdn': albdn, 'albdv': albdv, 'difnir': difnir,
                    'difvis': difvis, 'fnir': fnir, 'fvis': fvis,
                }
            )
        )

        albedo_s = rsoilv.expression(
            'fvis * rsoilv + fnir * rsoiln',
            {'fnir': fnir, 'fvis': fvis, 'rsoiln': rsoiln, 'rsoilv': rsoilv_iter}
        )

        albedo_avg = fc.expression(
            '(fc * albedo_c) + ((1 - fc) * albedo_s)',
            {'albedo_c': albedo_c, 'albedo_s': albedo_s, 'fc': fc}
        )
        diff = albedo_avg.subtract(albedo)

        # CGM - Check what this is doing
        # Extra select call is needed if LAI is multiband
        # Added fc_mask call
        fc_mask = fc.select([0]).lt(0.75)
        rsoilv_iter = (
            rsoilv_iter
            .where(fc_mask.And(diff.lte(-0.01)), rsoilv_iter.add(0.01))
            .where(fc_mask.And(diff.gt(0.01)), rsoilv_iter.add(-0.01))
        )

        # CGM - Extra select call is needed since fc is multiband
        fc_mask = fc.select([0]).gte(0.75)
        fg_iter = (
            fg_iter
            .where(fc_mask.And(diff.lte(-0.01)), fg_iter.subtract(0.05))
            .where(fc_mask.And(diff.gt(0.01)), fg_iter.add(0.05))
            .clamp(0.01, 1)
        )

        return ee.Dictionary({
            'akb': akb,
            'albedo_c': albedo_c,
            'albedo_s': albedo_s,
            'ameann': ameann,
            'ameanv': ameanv,
            'diff': diff,
            'fg': fg_iter,
            'rbcpyn': rbcpyn,
            'rbcpyv': rbcpyv,
            'rsoiln': rsoiln,
            'rsoilv': rsoilv_iter,
            'taudn': taudn,
            'taudv': taudv,
            'taudl': taudl,
        })

    # Iterate the function n times
    input_images = ee.Dictionary({
        'akb': None,
        'albedo_c': None,
        'albedo_s': None,
        'ameann': None,
        'ameanv': None,
        'diff': None,
        'fg': fg,
        'rbcpyn': None,
        'rbcpyv': None,
        'rsoiln': None,
        'rsoilv': rsoilv,
        'taudn': None,
        'taudv': None,
        'taudl': None,
    })
    iter_output = ee.Dictionary(
        ee.List.repeat(input_images, iterations).iterate(iter_func, input_images)
        # ee.List.sequence(1, iterations).iterate(iter_func, input_images)
    )

    # Unpack the iteration output
    akb = ee.Image(iter_output.get('akb'))
    albedo_c = ee.Image(iter_output.get('albedo_c'))
    albedo_s = ee.Image(iter_output.get('albedo_s'))
    ameann = ee.Image(iter_output.get('ameann'))
    ameanv = ee.Image(iter_output.get('ameanv'))
    diff = ee.Image(iter_output.get('diff'))
    rbcpyn = ee.Image(iter_output.get('rbcpyn'))
    rbcpyv = ee.Image(iter_output.get('rbcpyv'))
    rsoilv = ee.Image(iter_output.get('rsoilv'))
    rsoiln = ee.Image(iter_output.get('rsoiln'))
    # rsoiln = rsoilv.multiply(ratio_soil)
    taudn = ee.Image(iter_output.get('taudn'))
    taudv = ee.Image(iter_output.get('taudv'))
    taudl = ee.Image(iter_output.get('taudl'))

    # if a solution is not reached, alb_c=alb_s=alb
    albedo_c = albedo_c.where(diff.abs().gt(0.05), albedo)
    albedo_s = albedo_s.where(diff.abs().gt(0.05), albedo)

    # Direct beam+scattered canopy transmission coefficient (visible)
    expfac = ameanv.sqrt().multiply(akb).multiply(f)
    xnum = f.expression('(rbcpyv * rbcpyv - 1.0) * exp(-expfac)', {'rbcpyv': rbcpyv, 'expfac': expfac})
    xden = f.expression(
        '(rbcpyv * rsoilv - 1.0) + rbcpyv * (rbcpyv - rsoilv) * exp(-2.0 * expfac)',
        {'rbcpyv': rbcpyv, 'rsoilv': rsoilv, 'expfac': expfac}
    )
    # Eq 15.11
    taubtv = xnum.divide(xden)

    # Direct beam+scattered canopy transmission coefficient (NIR)
    expfac = ameann.sqrt().multiply(akb).multiply(f)
    xnum = f.expression('(rbcpyn * rbcpyn - 1.0) * exp(-expfac)', {'rbcpyn': rbcpyn, 'expfac': expfac})
    xden = f.expression(
        '(rbcpyn * rsoiln - 1.0) + rbcpyn * (rbcpyn - rsoiln) * exp(-2.0 * expfac)',
        {'rbcpyn': rbcpyn, 'rsoiln': rsoiln, 'expfac': expfac}
    )
    # Eq 15.11
    taubtn = xnum.divide(xden)

    # Shortwave radiation components
    tausolar = f.expression(
        'fvis * (difvis * taudv + dirvis * taubtv) + fnir * (difnir * taudn + dirnir * taubtn)',
        {
            'difnir': difnir, 'difvis': difvis, 'dirnir': dirnir, 'dirvis': dirvis,
            'fnir': fnir, 'fvis': fvis,
            'taubtn': taubtn, 'taubtv': taubtv, 'taudn': taudn, 'taudv': taudv,
        })

    rs_c = tausolar.multiply(-1).add(1).multiply(rs_1)
    rs_s = rs_1.multiply(tausolar)

    return rs_c, rs_s, albedo_c, albedo_s, taudl, tausolar


def compute_G0(rn, rn_s, ef_s, water_mask, lon, datetime):
    """

    Parameters
    ----------
    rn : ee.Image
    rn_s : ee.Image
    ef_s : ee.Image
        Evaporative fraction from the soil?
    water_mask : ee.Image
    lon : ee.Image
    datetime: datetime

    Returns
    -------
    g0 : ee.Image

    """
    w = ef_s.divide(0.5).pow(8).add(1).pow(-1)

    # Maximum fraction of Rn,s that become G0
    # (0.35 for dry soil and 0.31 for wet soil)
    c_g = w.expression('(w * 0.35) + ((1 - w) * 0.31)', {'w': w})
    t_g = w.expression('(w * 100000.0) + ((1 - w) * 74000.0)', {'w': w})

    time = ee.Date(datetime).get('hour').add(ee.Date(datetime).get('minute').divide(60))
    t_noon = solar_noon(datetime=datetime, lon=lon.multiply(math.pi / 180))
    t_g0 = t_noon.multiply(-1).add(time).multiply(3600)

    g0_temp = t_g0.add(10800).multiply(2 * math.pi).divide(t_g).cos()
    g0 = rn_s.multiply(c_g).multiply(g0_temp)
    g0 = g0.where(water_mask, rn.multiply(0.5))

    return g0


def compute_u_attr(u, d0, z0m, z_u, fm):
    """Friction Velocity

    Parameters
    ----------
    u : ee.Image
    d0 : ee.Image
    z0m : ee.Image
    z_u : ee.Image
    fm : ee.Image

    Returns
    -------
    u_attr : ee.Image

    """
    u_attr = d0.multiply(-1).add(z_u).divide(z0m).log().subtract(fm).pow(-1).multiply(u).multiply(0.41)
    u_attr = u_attr.where(u_attr.eq(0), 10)
    u_attr = u_attr.where(u_attr.lte(0), 0.01)

    return u_attr


def compute_r_ah(u_attr, d0, z0h, z_t, fh):
    """

    Parameters
    ----------
    u_attr : ee.Image
    d0 : ee.Image
    z0h : ee.Image
    z_t : ee.Image
    fh : ee.Image

    Returns
    -------
    r_ah : ee.Image

    """
    r_ah = d0.multiply(-1).add(z_t).divide(z0h).log().subtract(fh).divide(u_attr).divide(0.41)
    r_ah = r_ah.where(r_ah.eq(0), 500)
    r_ah = r_ah.max(1)

    return r_ah


def compute_r_s(u_attr, t_s, t_c, hc, f, d0, z0m, leaf, leaf_s, fm_h):
    """

    Parameters
    ----------
    u_attr : ee.Image
    t_s : ee.Image
        Soil temperature (Kelvin).
    t_c : ee.Image
        Canopy temperature (Kelvin).
    hc : ee.Image
    f : ee.Image
        Input is LAI?
    d0 : ee.Image
    z0m : ee.Image
    leaf : ee.Image
    leaf_s : ee.Image
    fm_h : ee.Image

    Returns
    -------
    r_s : ee.Image

    """
    # Free convective velocity constant for r_s modelling
    c_a = 0.004
    # Empirical constant for r_s modelling
    c_b = 0.012
    # Empirical constant for r_s modelling
    # (new formulation Kustas and Norman, 1999)
    c_c = 0.0025

    # Computation of the resistance of the air between soil and canopy space
    u_c = d0.multiply(-1).add(hc).divide(z0m).log().subtract(fm_h).multiply(u_attr).divide(0.41)
    u_c = u_c.where(u_c.lte(0), 0.1)

    u_s = hc.pow(-1).multiply(0.05).add(-1).multiply(leaf).exp().multiply(u_c)

    r_ss = (
        hc.pow(-1).multiply(-0.05).add(1).multiply(leaf_s).multiply(-1).exp()
        .multiply(c_b).add(c_a).pow(-1)
    )
    r_s1 = t_s.subtract(t_c).abs().pow(1.0 / 3).multiply(c_c).add(u_s.multiply(c_b)).pow(-1)
    r_s2 = u_s.multiply(c_b).add(c_a).pow(-1)
    r_s = f.subtract(0.01).multiply(0.09).pow(-1).multiply(r_ss.subtract(1)).add(1)

    # Linear function between 0 (bare soil) and the value at f=0.1
    r_s = r_s.where(f.gt(0.1), r_s1)
    r_s = r_s.where(t_s.subtract(t_c).abs().lt(1), r_s2)

    # Use "new" formula only for high DT values
    # Use "new" formula only for partial coverage (lai<3)
    r_s = r_s.where(f.gt(3), r_s2)

    return r_s


def compute_r_x(u_attr, hc, f, d0, z0m, xl, leaf_c, fm_h):
    """

    Parameters
    ----------
    u_attr : ee.Image
    hc : ee.Image
    f : ee.Image
        Input is LAI?
    d0 : ee.Image
    z0m : ee.Image
    xl : ee.Image
        Leaf width
    leaf_c : ee.Image
    fm_h : ee.Image

    Returns
    -------
    r_x : ee.Image

    """

    # Parameter for canopy boundary-layer resistance
    # (C=90 Grace '81, C=175 Cheubouni 2001, 144 Li '98)
    c = 175.0

    # Computation of the resistance of the air between soil and canopy space
    u_c = d0.multiply(-1).add(hc).divide(z0m).log().subtract(fm_h).multiply(u_attr).divide(0.41)
    u_c = u_c.where(u_c.lte(0), 0.1)

    # Computation of the canopy boundary layer resistance
    u_d = d0.add(z0m).divide(hc).add(-1).multiply(leaf_c).exp().multiply(u_c)
    u_d = u_d.where(u_d.lte(0), 100)

    r_x = xl.divide(u_d).pow(0.5).multiply(c).divide(f).where(u_d.eq(100), 0.1)

    return r_x


def compute_Rn_c(albedo_c, t_air, t_c, t_s, e_atm, rs_c, f):
    """Compute Canopy Net Radiation

    Parameters
    ----------
    albedo_c : ee.Image
    t_air : ee.Image
        Air temperature (Kelvin).
    t_c : ee.Image
        Canopy temperature (Kelvin).
    t_s : ee.Image
        Soil temperature (Kelvin).
    e_atm : ee.Image
    rs_c : ee.Image
    f : ee.Image

    Returns
    -------
    rn_c : ee.Image

    """
    # Long-wave extinction coefficient [-]
    kl = 0.95
    # Soil Emissivity [-]
    eps_s = 0.94
    # Canopy emissivity [-]
    eps_c = 0.99

    # Stephan Boltzmann constant (W m-2 K-4)
    # CGM - Original DisALEXI code used shortened version
    sb = 5.67E-8
    # sb = 5.670373e-8

    lc = t_c.pow(4).multiply(sb).multiply(eps_c)
    ls = t_s.pow(4).multiply(sb).multiply(eps_s)
    rle = t_air.pow(4).multiply(sb).multiply(e_atm)

    rn_c = albedo_c.expression(
        '((1 - albedo_c) * rs_c) + ((1 - exp(-kl * f)) * (rle + ls - 2 * lc))',
        {'albedo_c': albedo_c, 'f': f, 'kl': kl, 'lc': lc, 'ls': ls, 'rle': rle, 'rs_c': rs_c}
    )

    return rn_c


def compute_Rn_s(albedo_s, t_air, t_c, t_s, e_atm, rs_s, f):
    """Compute Soil Net Radiation

    Parameters
    ----------
    albedo_s : ee.Image
    t_air : ee.Image
        Air temperature (Kelvin).
    t_c : ee.Image
        Canopy temperature (Kelvin).
    t_s : ee.Image
        Soil temperature (Kelvin).
    e_atm : ee.Image
    rs_c : ee.Image
    f : ee.Image

    Returns
    -------
    rn_s : ee.Image

    """
    # Long-wave extinction coefficient [-]
    kl = 0.95
    # Soil Emissivity [-]
    eps_s = 0.94
    # Canopy emissivity [-]
    eps_c = 0.99

    # Stephan Boltzmann constant (W m-2 K-4)
    # CGM - Original DisALEXI code used shortened version version
    sb = 5.67E-8
    # sb = 5.670373e-8

    lc = t_c.pow(4).multiply(sb).multiply(eps_c)
    ls = t_s.pow(4).multiply(sb).multiply(eps_s)
    rle = t_air.pow(4).multiply(sb).multiply(e_atm)

    rn_s = albedo_s.expression(
        '((1 - albedo_s) * rs_s) + ((exp(-kl * f)) * rle) + ((1 - exp(-kl * f)) * lc) - ls',
        {'albedo_s': albedo_s, 'f': f, 'kl': kl, 'lc': lc, 'ls': ls, 'rle': rle, 'rs_s': rs_s}
    )

    return rn_s


def temp_separation_tc(h_c, fc, t_air, t0, r_ah, r_s, r_x, r_air, cp=1004.16):
    """Compute canopy temperature

    Parameters
    ----------
    h_c : ee.Image
    fc : ee.Image
    t_air : ee.Image
        Air temperature (Kelvin).
    t0 : ee.Image
    r_ah : ee.Image
    r_s : ee.Image
        Soil aerodynamic resistance to heat transport (s m-1).
    r_x : ee.Image
        Bulk canopy aerodynamic resistance to heat transport (s m-1).
    r_air : ee.Image
    cp : float

    Returns
    -------
    t_c : ee.Image
        Canopy temperature (Kelvin).

    """
    t_c_lin = fc.expression(
        '((t_air / r_ah) + (t0 / r_s / (1 - fc)) + '
        ' (h_c * r_x / r_air / cp * ((1 / r_ah) + (1 / r_s) + (1 / r_x)))) / '
        '((1 / r_ah) + (1 / r_s) + (fc / r_s / (1 - fc)))',
        {
            'cp': cp, 'fc': fc, 'h_c': h_c, 'r_ah': r_ah, 'r_air': r_air,
            'r_s': r_s, 'r_x': r_x, 't0': t0, 't_air': t_air,
        }
    )

    Td = fc.expression(
        '(t_c_lin * (1 + (r_s / r_ah))) - '
        '(h_c * r_x / r_air / cp * (1 + (r_s / r_x) + (r_s / r_ah))) - '
        '(t_air * r_s / r_ah)',
        {
            'cp': cp, 'h_c': h_c, 'r_ah': r_ah, 'r_air': r_air, 'r_s': r_s,
            'r_x': r_x, 't_air': t_air, 't_c_lin': t_c_lin,
        }
    )

    delta_t_c = fc.expression(
        '((t0 ** 4) - (fc * (t_c_lin ** 4)) - ((1 - fc) * (Td ** 4))) / '
        '((4 * (1 - fc) * (Td ** 3) * (1 + (r_s / r_ah))) + (4 * fc * (t_c_lin ** 3)))',
        {'fc': fc, 'r_ah': r_ah, 'r_s': r_s, 't0': t0, 'Td': Td, 't_c_lin': t_c_lin}
    )

    t_c = (
        t_c_lin.add(delta_t_c)
        .where(fc.lt(0.10), t0)
        .where(fc.gt(0.90), t0)
        .max(t_air.subtract(10.0))
        .min(t_air.add(50.0))
    )

    return t_c


def temp_separation_ts(t_c, fc, t_air, t0):
    """Compute soil temperature

    Parameters
    ----------
    t_c : ee.Image
        Canopy temperature (Kelvin).
    fc : ee.Image
    t_air : ee.Image
        Air temperature (Kelvin).
    t0 : ee.Image

    Returns
    -------
    t_s : ee.Image

    """
    delta = t_c.pow(4).multiply(fc).multiply(-1).add(t0.pow(4))
    delta = delta.where(delta.lte(0), 10)

    # CGM - This could probably be simplified further
    # t_s = (delta / (1 - fc)) ** 0.25
    # where:
    #     ((t0 ** 4) - (fc * t_c ** 4)) <= 0
    #     t_s = (t0 - (fc * t_c)) / (1 - fc)
    t_s = (
        fc.multiply(-1).add(1).pow(-1).multiply(delta).pow(0.25)
        .where(
            t0.pow(4).lte(t_c.pow(4).multiply(fc)),
            fc.multiply(t_c).multiply(-1).add(t0).divide(fc.multiply(-1).add(1))
        )
        .where(fc.lt(0.1), t0)
        .where(fc.gt(0.9), t0)
        .max(t_air.subtract(10.0))
        .min(t_air.add(50.0))
    )

    return t_s


def temp_separation_tac(t_c, t_s, fc, t_air, r_ah, r_s, r_x):
    """Compute air temperature at the canopy interface

    Parameters
    ----------
    t_c : ee.Image
        Canopy temperature (Kelvin).
    t_s : ee.Image
        Soil temperature (Kelvin).
    fc : ee.Image
    t_air : ee.Image
    r_ah : ee.Image
    r_s : ee.Image
    r_x : ee.Image

    Returns
    -------
    t_ac : ee.Image
        Air temperature at the canopy interface (Kelvin).

    """
    t_ac = fc.expression(
        '((t_air / r_ah) + (t_s / r_s) + (t_c / r_x)) / ((1 / r_ah) + (1 / r_s) + (1 / r_x))',
        {'r_ah': r_ah, 'r_s': r_s, 'r_x': r_x, 't_c': t_c, 't_s': t_s, 't_air': t_air}
    )

    return t_ac


def compute_stability_fh(h, t0, u_attr, r_air, z_t, d0, cp=1004.16):
    """

    To Match IDL, fh defaults to 0 for L values outside the range -100 to 100
    It is unclear why the range -100 to 100 is used since L is forced to be < 0

    Parameters
    ----------
    h : ee.Image
    t0 : ee.Image
    u_attr : ee.Image
    r_air : ee.Image
    z_t : float
    d0 : ee.Image
    cp : float

    Returns
    -------
    fh : ee.Image

    """
    l_ob = u_attr.pow(3).multiply(t0).multiply(r_air).multiply(cp / -0.41 / 9.806).divide(h)
    l_ob = l_ob.where(l_ob.gte(0), -99)

    mh = d0.subtract(z_t).multiply(16).divide(l_ob).add(1).pow(0.25).where(l_ob.eq(-99), 0.0)
    fh = mh.pow(2).add(1).divide(2).log().multiply(2).where(l_ob.lte(-100).Or(l_ob.gte(100)), 0)

    return fh


def compute_stability_fm(h, t0, u_attr, r_air, z_u, d0, z0m, cp=1004.16):
    """

    To Match IDL, fh defaults to 0 for L values outside the range -100 to 100
    It is unclear why the range -100 to 100 is used since L is forced to be < 0

    Parameters
    ----------
    h : ee.Image
    t0 : ee.Image
    u_attr : ee.Image
    r_air : ee.Image
    z_u : float
    d0 : ee.Image
    z0m : ee.Image
    cp : float

    Returns
    -------
    fm : ee.Image

    """

    l_ob = u_attr.pow(3).multiply(t0).multiply(r_air).multiply(cp / -0.41 / 9.806).divide(h)
    l_ob = l_ob.where(l_ob.gte(0), -99)

    mh = d0.subtract(z_u).multiply(16).divide(l_ob).add(1).pow(0.25).where(l_ob.eq(-99), 0.0)
    fm = (
        h.expression(
            '2.0 * log((1.0 + mh) / 2.0) + log((1.0 + (mh ** 2)) / 2.0) - 2.0 * atan(mh) + (pi / 2)',
            {'mh': mh, 'pi': math.pi}
        )
        .where(l_ob.lte(-100).Or(l_ob.gte(100)), 0)
    )

    # CGM - Swapped order of calc since d0 is an image computed from hc and
    #   z_u is being set as a constant number (for now).
    fm = fm.where(fm.eq(d0.multiply(-1).add(z_u).divide(z0m).log()), fm.add(1.0))

    return fm


def compute_stability_fm_h(h, t0, u_attr, r_air, hc, d0, z0m, cp=1004.16):
    """

    To Match IDL, fh defaults to 0 for L values outside the range -100 to 100
    It is unclear why the range -100 to 100 is used since L is forced to be < 0

    Parameters
    ----------
    h : ee.Image
    t0 : ee.Image
    u_attr : ee.Image
    r_air : ee.Image
    hc : ee.Image
    d0 : ee.Image
    z0m : ee.Image
    cp : float

    Returns
    -------
    fm_h : ee.Image

    """
    l_ob = u_attr.pow(3).multiply(t0).multiply(r_air).multiply(cp / -0.41 / 9.806).divide(h)
    l_ob = l_ob.where(l_ob.gte(0), -99)

    mm_h = d0.subtract(hc).multiply(16).divide(l_ob).add(1).pow(0.25).where(l_ob.eq(-99), 0.0)
    fm_h = (
        h.expression(
            '2.0 * log((1.0 + mm_h) / 2.0) + log((1.0 + (mm_h ** 2)) / 2.0) - 2.0 * atan(mm_h) + (pi / 2)',
            {'mm_h': mm_h, 'pi': math.pi}
        )
        .where(l_ob.lte(-100).Or(l_ob.gte(100)), 0)
    )

    # CGM - Swapped order of calc since d0 is an image computed from hc and
    #   z_u is being set as a constant number (for now).
    fm_h = fm_h.where(fm_h.eq(hc.subtract(d0).divide(z0m).log()), fm_h.add(1.0))

    return fm_h
