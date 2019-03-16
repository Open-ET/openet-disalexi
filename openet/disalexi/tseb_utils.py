import math

import ee

from . import utils

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

    t_noon = lon.expression(
        '(720.0 - 4 * lon - eq_t) / 1440 * 24.0',
        {'lon': lon.multiply(rad2deg), 'eq_t': eq_t})
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
    time_t = ee.Date(datetime).get('hour').add(
        ee.Date(datetime).get('minute').divide(60))
    # # This will return the hour floating point value
    # time_t = ee.Date(date).get('hour').add(ee.Date(date).getFraction('hour'))

    d, eq_t = _solar_time(datetime)

    ts_time = lon.expression(
        '(time_t / 24.0 * 1440 + eq_t + 4.0 * lon) % 1440.',
        {'lon': lon.multiply(rad2deg), 'time_t': time_t, 'eq_t': eq_t})
    ts_time = ts_time.where(ts_time.gt(1440), ts_time.subtract(1440))
    # ts_time[ts_time > 1440.] = ts_time[ts_time > 1440.] - 1440.

    w = lat.expression('ts_time / 4.0 + 180.0', {'ts_time': ts_time})
    w = w.where(ts_time.divide(4).gt(0), ts_time.divide(4).subtract(180))
    # w[ts_time/4.0 >= 0] = ts_time[ts_time/4.0 >= 0.] / 4.-180.

    zs = lat.expression(
        'acos( (sin(lat) * sin(d)) + (cos(lat) * cos(d) * cos(w)) )',
        {'lat': lat, 'd': d, 'w': w.multiply(deg2rad)})

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
    time_t = ee.Date(datetime).get('hour').add(
        ee.Date(datetime).get('minute').divide(60))

    # This will return the hour floating point value
    # time_t = ee.Date(date).get('hour').add(ee.Date(date).getFraction('hour'))

    # CGM - DOY and hour could be images in order to make these expressions
    julian = _to_jd(datetime)

    # Sunrise time
    julian_ = time_t.divide(24.0).add(julian)
    j_cen = julian_.add(0.5 - 2451545.0).divide(36525.0)
    # CGM - Does the mod happen before or after the multiply
    lon_sun = j_cen.multiply(0.0003032).add(36000.76983) \
        .multiply(j_cen).mod(360.0).add(280.46646).subtract(360)
    an_sun =  j_cen.multiply(-0.0001537).add(35999.05029) \
        .multiply(j_cen).add(357.52911)
    ecc = j_cen.multiply(0.0000001267).add(0.000042037) \
        .multiply(j_cen).multiply(-1).add(0.016708634)
    ob_ecl = j_cen.multiply(-0.001813).add(0.00059) \
        .multiply(j_cen).add(46.815) \
        .multiply(j_cen).multiply(-1).add(21.448) \
        .divide(60.0).add(26).divide(60).add(23)
    ob_corr = j_cen.multiply(-1934.136).add(125.04).multiply(deg2rad).cos() \
        .multiply(0.00256).add(ob_ecl)
    var_y = ob_corr.divide(2.0).multiply(deg2rad).tan().multiply(
        ob_corr.divide(2.0).multiply(deg2rad).tan())
    eq_t = (
        lon_sun.multiply(2.0).multiply(deg2rad).sin().multiply(var_y)
            .subtract(an_sun.multiply(deg2rad).sin().multiply(ecc).multiply(2.0))
            .add(an_sun.multiply(deg2rad).sin()
                .multiply(lon_sun.multiply(2.0).multiply(deg2rad).cos())
                .multiply(var_y).multiply(ecc).multiply(4.0))
            .subtract(lon_sun.multiply(4.0).multiply(deg2rad).sin()
                      .multiply(var_y).multiply(var_y).multiply(0.5))
            .subtract(an_sun.multiply(2.0).multiply(deg2rad).sin()
                      .multiply(ecc).multiply(ecc).multiply(1.25))
        .multiply(4.0).multiply(rad2deg))
    sun_eq = (
        an_sun.multiply(deg2rad).sin().multiply(
            j_cen.multiply(0.000014).add(0.004817)
                .multiply(j_cen).multiply(-1).add(1.914602))
        .add(an_sun.multiply(2.0).multiply(deg2rad).sin()
             .multiply(j_cen.multiply(-0.000101).add(0.019993)))
        .add(an_sun.multiply(3.0).multiply(deg2rad).sin().multiply(0.000289)))
    sun_true = sun_eq.add(lon_sun)
    sun_app = j_cen.multiply(-1934.136).add(125.04).multiply(deg2rad).sin() \
        .multiply(-0.00478).subtract(0.00569).add(sun_true)

    # CGM - Intentionally not converting back to degrees
    d = ob_corr.multiply(deg2rad).sin() \
        .multiply(sun_app.multiply(deg2rad).sin()) \
        .asin()

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

    jdn = ee.Date(date).get('day') \
        .add(m.multiply(153).add(2).divide(5).floor()) \
        .add(y.multiply(365)) \
        .add(y.divide(4.0).floor()) \
        .subtract(y.divide(100.0).floor()) \
        .add(y.divide(400.0).floor()) \
        .subtract(32045)

    # To match IDL, only return JDN term
    return jdn

    # CGM - Forcing to float types to match wikipedia equations.
    #   Original equation was doing integer division in Python 2.7.
    #   Microseconds are not in the original equation.
    # jd = jdn \
    #     .add(ee.Date(date).get('hour').float().subtract(12).divide(24)) \
    #     .add(ee.Date(date).get('minute').float().divide(1440)) \
    #     .add(ee.Date(date).get('second').float().divide(86400))
    #     # .add(date.get('microsecond') / 86400000000)
    # return jd


def emissivity(T_air):
    """Apparent atmospheric emissivity

    Apparent emissivity (Sedlar and Hock, 2009: Cryosphere 3:75-84)
    Atmospheric emissivity (clear-sky) Idso and Jackson (1969)

    Parameters
    ----------
    T_air : ee.Image
        Air temperature (Kelvin)

    Returns
    -------
    e_atm : ee.Image

    """
    e_atm = T_air \
        .expression(
            '1.0 - (0.2811 * (exp(-0.0003523 * ((T_air - 273.16) ** 2.0))))',
            {'T_air': T_air})
        # .where(Rs0.lte(50.0), 1.0)
    return e_atm


# CGM - Not including parameter z since is not used in this function
def albedo_separation(albedo, Rs_1, F, fc, aleafv, aleafn, aleafl, adeadv,
                      adeadn, adeadl, zs, iterations=10):
    """Compute Solar Components and atmospheric properties ([Campbell1998]_)

    Parameters
    ----------
    albedo : ee.Image
    Rs_1 : ee.Image
    F : ee.Image
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

    References
    ----------
    .. [Campbell1998] Campbell, G. S. & Norman, J. M. (1998),
        An introduction to environmental biophysics. Springer, New York
        https://archive.org/details/AnIntroductionToEnvironmentalBiophysics.

    """
    # DAYTIME
    # Calculate potential (clear-sky) VIS and NIR solar components

    rad2deg = 180.0 / math.pi
    # deg2rad = math.pi / 180.0

    # Correct for curvature of atmos in airmas
    airmas = zs.expression(
        '(sqrt(cos(zs) ** 2 + 0.0025) - cos(zs)) / 0.00125', {'zs': zs})

    # Correct for refraction(good up to 89.5 deg.)
    airmas = airmas.where(
        zs.multiply(rad2deg).lt(89.5),
        zs.expression(
            'airmas - (2.8 / (90.0 - zs_temp) ** 2)',
            {'airmas': airmas, 'zs_temp': zs.multiply(rad2deg)}))

    potbm1 = zs.expression(
        '600.0 * exp(-0.160 * airmas)', {'airmas': airmas})
    potvis = zs.expression(
        '(potbm1 + (600.0 - potbm1) * 0.4) * cos(zs)',
        {'potbm1': potbm1, 'zs': zs})
    # CGM - Not used
    # potdif = zs.expression(
    #     '(600.0 - potbm1) * 0.4 * cos(zs)', {'potbm1': potbm1, 'zs': zs})
    uu = zs.expression('1.0 / cos(zs)', {'zs': zs}).max(0.01)
    a = zs.expression(
        '10 ** (-1.195 + 0.4459 * axlog - 0.0345 * axlog * axlog)',
        {'axlog': uu.log10()})
    watabs = a.multiply(1320.0)
    potbm2 = zs.expression(
        '720.0 * exp(-0.05 * airmas) - watabs',
        {'airmas': airmas, 'watabs': watabs})
    evaL = zs.expression(
        '(720.0 - potbm2 - watabs) * 0.54 * cos(zs)',
        {'potbm2': potbm2, 'watabs': watabs, 'zs': zs})
    potnir = zs.expression(
        'evaL + potbm2 * cos(zs)',
        {'evaL': evaL, 'potbm2': potbm2, 'zs': zs})

    fclear = zs \
        .expression(
            'Rs_1 / (potvis + potnir)',
            {'potvis': potvis, 'potnir': potnir, 'Rs_1': Rs_1}) \
        .clamp(0.01, 1.0) \
        .where(zs.cos().lte(0.01), 1)
  
    # Partition SDN into VIS and NIR
    fvis = zs.expression(
        'potvis / (potvis + potnir)', {'potvis': potvis, 'potnir': potnir})
    fnir = zs.expression(
        'potnir / (potvis + potnir)', {'potvis': potvis, 'potnir': potnir})
  
    # Estimate direct beam and diffuse fraction in VIS and NIR wavebands
    fb1 = zs.expression(
        'potbm1 * cos(zs) / potvis',
        {'potbm1': potbm1, 'potvis': potvis, 'zs': zs})
    # CGM - Not used
    # fb2 = zs.expression(
    #     'potbm2 * cos(zs) / potnir',
    #     {'potbm2': potbm2, 'potnir': potnir, 'zs': zs})

    dirvis = zs \
        .expression(
            'fb1 * (1.0 - ((0.9 - ratiox) / 0.7) ** 0.6667)',
            {'fb1': fb1, 'ratiox': fclear.min(0.9)}) \
        .min(fb1)
    dirnir = zs \
        .expression(
            'fb1 * (1.0 - ((0.88 - ratiox) / 0.68) ** 0.6667)',
            {'fb1': fb1, 'ratiox': fclear.min(0.88)}) \
        .min(fb1)

    dirvis = dirvis.where(dirvis.lt(0.01).And(dirnir.gt(0.01)), 0.011)
    dirnir = dirnir.where(dirnir.lt(0.01).And(dirvis.gt(0.01)), 0.011)

    difvis = zs.expression('1.0 - dirvis', {'dirvis': dirvis})
    difnir = zs.expression('1.0 - dirnir', {'dirnir': dirnir})
  
    # Correction for NIGHTIME
    ind = zs.cos().lte(0.01)
    fvis = fvis.where(ind, 0.5)
    fnir = fnir.where(ind, 0.5)
    difvis = difvis.where(ind, 0.0)
    difnir = difnir.where(ind, 0.0)


    # CGM - Not used anymore in function since e_atm is not computed
    # Rs0 = zs \
    #     .expression('potvis + potnir', {'potnir': potnir, 'potvis': potvis}) \
    #     .where(zs.cos().lte(0.01), 0.0)

    #**********************************************
    # Compute Albedo
    ratio_soil = 2.

    # CGM - Initialize rsoilv and fg from F and albedo
    rsoilv = F.multiply(0).add(0.12)
    fg = albedo.multiply(0).add(1)
    # rsoilv = ee.Image.constant(0.12)
    # fg = ee.Image.constant(1.0)

    # CGM - Switched to an iterate call
    def iter_func(n, prev):
        # Extract inputs from previous iteration
        # CGM - Variables that are commented out only need to be returned
        # akb = ee.Image(ee.Dictionary(prev).get('akb'))
        # albedo_c = ee.Image(ee.Dictionary(prev).get('albedo_c'))
        # albedo_s = ee.Image(ee.Dictionary(prev).get('albedo_s'))
        # ameann = ee.Image(ee.Dictionary(prev).get('ameann'))
        # ameanv = ee.Image(ee.Dictionary(prev).get('ameanv'))
        # diff = ee.Image(ee.Dictionary(prev).get('diff'))
        fg_iter = ee.Image(ee.Dictionary(prev).get('fg'))
        # rbcpyn = ee.Image(ee.Dictionary(prev).get('rbcpyn'))
        # rbcpyv = ee.Image(ee.Dictionary(prev).get('rbcpyv'))
        rsoilv_iter = ee.Image(ee.Dictionary(prev).get('rsoilv'))
        # taudn = ee.Image(ee.Dictionary(prev).get('taudn'))
        # taudv = ee.Image(ee.Dictionary(prev).get('taudv'))

        rsoiln = rsoilv_iter.multiply(ratio_soil)
        # rsoiln = .expression(
        #     'rsoilv * ratio_soil',
        #     {'rsoilv': rsoilv, 'ratio_soil': ratio_soil})

        # Weighted live/dead leaf average properties
        ameanv = aleafv.expression(
            'aleafv * fg + adeadv * (1.0 - fg)',
            {'adeadv': adeadv, 'aleafv': aleafv, 'fg': fg_iter})
        ameann = aleafn.expression(
            'aleafn * fg + adeadn * (1.0 - fg)',
            {'adeadn': adeadn, 'aleafn': aleafn, 'fg': fg_iter})
        ameanl = aleafl.expression(
            'aleafl * fg + adeadl * (1.0 - fg)',
            {'adeadl': adeadl, 'aleafl': aleafl, 'fg': fg_iter})

        # DIFFUSE COMPONENT
        #*******************************
        # Canopy reflection (deep canopy)
        # Fit to Fig 15.4 for x=1
        akd = F.expression('-0.0683 * log(F) + 0.804', {'F': F})

        # Eq 15.7
        rcpyn = ameann.expression(
            '(1.0 - sqrt(ameann)) / (1.0 + sqrt(ameann))', {'ameann': ameann})
        rcpyv = ameanv.expression(
            '(1.0 - sqrt(ameanv)) / (1.0 + sqrt(ameanv))', {'ameanv': ameanv})
        # rcpyl = ameanl.expression(
        #     '(1.0 - sqrt(ameanl)) / (1.0 + sqrt(ameanl))', {'ameanl': ameanl})

        # Eq 15.8
        rdcpyn = akd.expression(
            '2.0 * akd * rcpyn / (akd + 1.0)', {'akd': akd, 'rcpyn': rcpyn})
        rdcpyv = akd.expression(
            '2.0 * akd * rcpyv / (akd + 1.0)', {'akd': akd, 'rcpyv': rcpyv})
        # rdcpyl = akd.expression(
        #     '2.0 * akd * rcpyl / (akd + 1.0)', {'akd': akd, 'rcpyl': rcpyl})

        # Canopy transmission (VIS)
        expfac = F.expression(
            'sqrt(ameanv) * akd * F', {'akd': akd, 'ameanv': ameanv, 'F': F})
        expfac = expfac.max(0.001)
        # expfac = expfac.where(expfac.lt(0.001), 0.001)
        xnum = F.expression(
            '(rdcpyv * rdcpyv - 1.0) * exp(-expfac)',
            {'rdcpyv': rdcpyv, 'expfac': expfac})
        xden = F.expression(
            '(rdcpyv * rsoilv - 1.0) + '
            'rdcpyv * (rdcpyv - rsoilv) * exp(-2.0 * expfac)',
            {'expfac': expfac, 'rdcpyv': rdcpyv, 'rsoilv': rsoilv_iter})
        # Eq 15.11
        taudv = F.expression('xnum / xden', {'xden': xden, 'xnum': xnum})
        # taudv = xnum.divide(xden)

        # Canopy transmission (NIR)
        expfac = F.expression(
            'sqrt(ameann) * akd * F', {'akd': akd, 'ameann': ameann, 'F': F})
        expfac = expfac.max(0.001)
        # expfac = expfac.where(expfac.lt(0.001), 0.001)
        xnum = F.expression(
            '(rdcpyn * rdcpyn - 1.0) * exp(-expfac)',
            {'expfac': expfac, 'rdcpyn': rdcpyn})
        xden = F.expression(
            '(rdcpyn * rsoiln - 1.0) + '
            'rdcpyn * (rdcpyn - rsoiln) * exp(-2.0 * expfac)',
            {'expfac': expfac, 'rdcpyn': rdcpyn, 'rsoiln': rsoiln})
        # Eq 15.11
        taudn = F.expression('xnum / xden', {'xden': xden, 'xnum': xnum})
        # taudn = xnum.divide(nden)

        # Canopy transmission (LW)
        taudl = F.expression(
            'exp(-sqrt(ameanl) * akd * F)',
            {'akd': akd, 'ameanl': ameanl, 'F': F})

        # Diffuse albedo for generic canopy
        # Eq 15.9
        fact = F.expression(
            '((rdcpyn - rsoiln) / (rdcpyn * rsoiln - 1.0)) * '
            'exp(-2.0 * sqrt(ameann) * akd * F)',
            {'akd': akd, 'ameann': ameann, 'F': F, 'rdcpyn': rdcpyn,
             'rsoiln': rsoiln})
        albdn = F.expression(
            '(rdcpyn + fact) / (1.0 + rdcpyn * fact)',
            {'fact': fact, 'rdcpyn': rdcpyn})

        # Eq 15.9
        fact = F.expression(
            '((rdcpyv - rsoilv) / (rdcpyv * rsoilv - 1.0)) * '
            'exp(-2.0 * sqrt(ameanv) * akd * F)',
            {'akd': akd, 'ameanv': ameanv, 'F': F, 'rdcpyv': rdcpyv,
             'rsoilv': rsoilv_iter})
        albdv = F.expression(
            '(rdcpyv + fact) / (1.0 + rdcpyv * fact)',
            {'fact': fact, 'rdcpyv': rdcpyv})

        # BEAM COMPONENT
        #*******************************
        # Canopy reflection (deep canopy)
        akb = zs.expression('0.5 / cos(zs)', {'zs': zs})
        akb = akb.where(zs.cos().lte(0.01), 0.5)

        # Eq 15.7
        rcpyn = ameann.expression(
            '(1.0 - sqrt(ameann)) / (1.0 + sqrt(ameann))',
            {'ameann': ameann})
        rcpyv = ameanv.expression(
            '(1.0 - sqrt(ameanv)) / (1.0 + sqrt(ameanv))',
            {'ameanv': ameanv})

        # Eq 15.8
        rbcpyn = rcpyn.expression(
            '2.0 * akb * rcpyn / (akb + 1.0)', {'akb': akb, 'rcpyn': rcpyn})
        rbcpyv = rcpyv.expression(
            '2.0 * akb * rcpyv / (akb + 1.0)', {'akb': akb, 'rcpyv': rcpyv})

        # Beam albedo for generic canopy
        # Eq 15.9
        fact = F.expression(
            '((rbcpyn - rsoiln) / (rbcpyn * rsoiln - 1.0)) * '
            'exp(-2.0 * sqrt(ameann) * akb * F)',
            {'akb': akb, 'ameann': ameann, 'F': F, 'rbcpyn': rbcpyn,
             'rsoiln': rsoiln})
        albbn = F.expression(
            '(rbcpyn + fact) / (1.0 + rbcpyn * fact)',
            {'fact': fact, 'rbcpyn': rbcpyn})

        # Eq 15.9
        fact = F.expression(
            '((rbcpyv - rsoilv) / (rbcpyv * rsoilv - 1.0)) * '
            'exp(-2.0 * sqrt(ameanv) * akb * F)',
            {'akb': akb, 'ameanv': ameanv, 'F': F, 'rbcpyv': rbcpyv,
             'rsoilv': rsoilv_iter})
        albbv = F.expression(
            '(rbcpyv + fact) / (1.0 + rbcpyv * fact)',
            {'fact': fact, 'rbcpyv': rbcpyv})

        # CGM - finish
        # Weighted albedo (canopy)
        albedo_c = F.expression(
            'fvis * (dirvis * albbv + difvis * albdv) + '
            'fnir * (dirnir * albbn + difnir * albdn)',
            {'albbn': albbn, 'albbv': albbv, 'albdn': albdn, 'albdv': albdv,
             'difnir': difnir, 'difvis': difvis, 'dirvis': dirvis,
             'dirnir': dirnir, 'fnir': fnir, 'fvis': fvis, })
        albedo_c = albedo_c.where(
            zs.cos().lte(0.01),
            F.expression(
                'fvis * (difvis * albdv) + fnir * (difnir * albdn)',
                {'albdn': albdn, 'albdv': albdv, 'difnir': difnir,
                 'difvis': difvis, 'fnir': fnir, 'fvis': fvis}))

        albedo_s = rsoilv.expression(
            'fvis * rsoilv + fnir * rsoiln',
            {'fnir': fnir, 'fvis': fvis, 'rsoiln': rsoiln, 'rsoilv': rsoilv_iter})

        albedo_avg = fc.expression(
            '(fc * albedo_c) + ((1 - fc) * albedo_s)',
            {'albedo_c': albedo_c, 'albedo_s': albedo_s, 'fc': fc})
        diff = albedo_avg.subtract(albedo)
        # diff = albedo_avg.expression(
        #     'albedo_avg - albedo',
        #     {'albedo_avg': albedo_avg, 'albedo': albedo})

        # CGM - Check what this is doing
        # Extra select call is needed if LAI is multiband
        # Added fc_mask call
        fc_mask = fc.select([0]).lt(0.75)
        rsoilv_iter = rsoilv_iter \
            .where(fc_mask.And(diff.lte(-0.01)), rsoilv_iter.add(0.01)) \
            .where(fc_mask.And(diff.gt(0.01)), rsoilv_iter.add(-0.01))
        # # CGM - IDL function
        # rsoilv = ((fc lt 0.75) * (
        #             ((abs(diff) le 0.01) * rsoilv) +
        #             ((diff le -0.01)*(rsoilv + 0.01)) +
        #             ((diff gt 0.01)*(rsoilv - 0.01))))+
        #          ((fc ge 0.75) * rsoilv)

        # CGM - Extra select call is needed since fc is multiband
        fc_mask = fc.select([0]).gte(0.75)
        fg_iter = fg_iter \
            .where(fc_mask.And(diff.lte(-0.01)), fg_iter.subtract(0.05)) \
            .where(fc_mask.And(diff.gt(0.01)), fg_iter.add(0.05)) \
            .clamp(0.01, 1)
        # # CGM - IDL function
        # fg = ((fc ge 0.75) * (
        #          ((abs(diff) le 0.01)*fg) +
        #          ((diff le -0.01) * (fg - 0.05d0)) +
        #          ((diff gt 0.01) * (fg + 0.05d0)))) +
        #      ((fc lt 0.75) * fg)

        return ee.Dictionary({
            'akb': akb, 'albedo_c': albedo_c, 'albedo_s': albedo_s,
            'ameann': ameann, 'ameanv': ameanv, 'diff': diff, 'fg': fg_iter,
            'rbcpyn': rbcpyn, 'rbcpyv': rbcpyv,
            'rsoiln': rsoiln, 'rsoilv': rsoilv_iter,
            'taudn': taudn, 'taudv': taudv
        })

    # Iterate the function n times
    input_images = ee.Dictionary({
        'akb': None, 'albedo_c': None, 'albedo_s': None,
        'ameann': None, 'ameanv': None, 'diff': None, 'fg': fg,
        'rbcpyn': None, 'rbcpyv': None,
        'rsoiln': None, 'rsoilv': rsoilv,
        'taudn': None, 'taudv': None
    })
    iter_output = ee.Dictionary(
        # ee.List.sequence(1, iterations) \
        ee.List.repeat(input_images, iterations) \
            .iterate(iter_func, input_images))

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

    # if a solution is not reached, alb_c=alb_s=alb
    albedo_c = albedo_c.where(diff.abs().gt(0.05), albedo)
    albedo_s = albedo_s.where(diff.abs().gt(0.05), albedo)

    # Direct beam+scattered canopy transmission coefficient (visible)
    expfac = F.expression(
        'sqrt(ameanv) * akb * F',
        {'ameanv': ameanv, 'akb': akb, 'F': F})
    xnum = F.expression(
        '(rbcpyv * rbcpyv - 1.0) * exp(-expfac)',
        {'rbcpyv': rbcpyv, 'expfac': expfac})
    xden = F.expression(
        '(rbcpyv * rsoilv - 1.0) + '
        'rbcpyv * (rbcpyv - rsoilv) * exp(-2.0 * expfac)',
        {'rbcpyv': rbcpyv, 'rsoilv': rsoilv, 'expfac': expfac})
    # Eq 15.11
    taubtv = F.expression('xnum / xden', {'xnum': xnum, 'xden': xden})

    # Direct beam+scattered canopy transmission coefficient (NIR)
    expfac = F.expression(
        'sqrt(ameann) * akb * F',
        {'ameann': ameann, 'akb': akb, 'F': F})
    xnum = F.expression(
        '(rbcpyn * rbcpyn - 1.0) * exp(-expfac)',
        {'rbcpyn': rbcpyn, 'expfac': expfac})
    xden = F.expression(
        '(rbcpyn * rsoiln - 1.0) + '
        'rbcpyn * (rbcpyn - rsoiln) * exp(-2.0 * expfac)',
        {'rbcpyn': rbcpyn, 'rsoiln': rsoiln, 'expfac': expfac})
    # Eq 15.11
    taubtn = F.expression('xnum / xden', {'xnum': xnum, 'xden': xden})

    # Shortwave radiation components
    tausolar = F.expression(
        'fvis * (difvis * taudv + dirvis * taubtv) + '
        'fnir * (difnir * taudn + dirnir * taubtn)',
        {'difnir': difnir, 'difvis': difvis,
         'dirnir': dirnir, 'dirvis': dirvis,
         'fnir': fnir, 'fvis': fvis,
         'taubtn': taubtn, 'taubtv': taubtv,
         'taudn': taudn, 'taudv': taudv})

    Rs_c = Rs_1.expression(
        'Rs_1 * (1.0 - tausolar)', {'Rs_1': Rs_1, 'tausolar': tausolar})
    Rs_s = Rs_1.expression(
        'Rs_1 * tausolar', {'Rs_1': Rs_1, 'tausolar': tausolar})

    return Rs_c, Rs_s, albedo_c, albedo_s


def compute_G0(Rn, Rn_s, albedo, ndvi, t_noon, time, EF_s):
    """

    Parameters
    ----------
    Rn : ee.Image
    Rn_s : ee.Image
    albedo : ee.Image
    ndvi : ee.Image
    t_noon: ee.Image
    time :
    EF_s :

    Returns
    -------
    G0 : ee.Image

    """
    w = EF_s.expression('1 / (1 + (EF_s / 0.5) ** 8.0)', {'EF_s': EF_s})

    # Maximum fraction of Rn,s that become G0
    # (0.35 for dry soil and 0.31 for wet soil)
    c_g = w.expression('(w * 0.35) + ((1 - w) * 0.31)', {'w': w})
    t_g = w.expression('(w * 100000.0) + ((1 - w) * 74000.0)', {'w': w})

    t_g0 = t_noon.expression(
        '(time - t_noon) * 3600.0', {'time': time, 't_noon': t_noon})
      
    G0 = Rn_s.expression(
        'c_g * cos(2 * pi * (t_g0 + 10800.0) / t_g) * Rn_s',
        {'c_g': c_g, 'pi': math.pi, 'Rn_s': Rn_s, 't_g': t_g, 't_g0': t_g0})

    water_mask = ndvi.lte(0).And(albedo.lte(0.05))
    G0 = G0.where(water_mask, Rn.multiply(0.5))

    return G0


def compute_resistance(u, T_s, T_c, hc, F, d0, z0m, z0h, z_u, z_t, xl,
                       leaf, leaf_s, leaf_c, fm, fh, fm_h):
    """

    Parameters
    ----------
    u : ee.Image
    T_s : ee.Image
    T_c : ee.Image
    hc : ee.Image
    F : ee.Image
        Input is LAI?
    d0
    z0m
    z0h
    z_u
    z_t
    xl
    leaf
    leaf_s
    leaf_c
    fm
    fh
    fm_h

    Returns
    -------
    r_ah
    r_s
    r_x
    u_attr

    """
    # Free convective velocity constant for r_s modelling
    c_a = 0.004
    # Empirical constant for r_s modelling
    c_b = 0.012
    # Empirical constant for r_s modelling
    # (new formulation Kustas and Norman, 1999)
    c_c = 0.0025

    # Parameter for canopy boundary-layer resistance
    # (C=90 Grace '81, C=175 Cheubouni 2001, 144 Li '98)
    C = 175.

    # Computation of friction velocity and aerodynamic resistance
    u_attr = u \
        .expression(
            '0.41 * u / ((log((z_u - d0) / z0m)) - fm)',
            {'d0': d0, 'fm': fm, 'u': u, 'z0m': z0m, 'z_u': z_u})
    u_attr = u_attr.where(u_attr.eq(0), 10)
    u_attr = u_attr.where(u_attr.lte(0), 0.01)

    r_ah = u.expression(
        '((log((z_t - d0) / z0h)) - fh) / u_attr / 0.41',
        {'d0': d0, 'fh': fh, 'u_attr': u_attr, 'z0h': z0h, 'z_t': z_t})
    # CGM - The second conditional will overwrite the first one?
    r_ah = r_ah.where(r_ah.eq(0), 500)
    r_ah = r_ah.where(r_ah.lte(1.0), 1.0)
    # DEADBEEF
    # r_ah[r_ah == 0] = 500.
    # r_ah[r_ah <= 1.] = 1.

    # Computation of the resistance of the air between soil and canopy space
    u_c = u.expression(
        'u_attr / 0.41 * ((log((hc - d0) / z0m)) - fm_h)',
        {'d0': d0, 'fm_h': fm_h, 'hc': hc, 'u_attr': u_attr, 'z0m': z0m})
    u_c = u_c.where(u_c.lte(0), 0.1)
    u_s = u.expression(
        'u_c * exp(-leaf * (1 - (0.05 / hc)))',
        {'hc': hc, 'leaf': leaf, 'u_c': u_c})

    r_ss = u.expression(
        '1.0 / (c_a + (c_b * (u_c * exp(-leaf_s * (1.0 - (0.05 / hc))))))',
        {'c_a': c_a, 'c_b': c_b, 'hc': hc, 'leaf_s': leaf_s, 'u_c': u_c})
    r_s1 = T_s.expression(
        '1.0 / ((((abs(T_s - T_c)) ** (1.0 / 3.0)) * c_c) + (c_b * u_s))',
        {'c_b': c_b, 'c_c': c_c, 'T_c': T_c, 'T_s': T_s, 'u_s': u_s})
    r_s2 = u.expression(
        '1.0 / (c_a + (c_b * u_s))', {'c_a': c_a, 'c_b': c_b, 'u_s': u_s})
    r_s = u.expression(
        '(((r_ss - 1.0) / 0.09 * (F - 0.01)) + 1.0)', {'F': F, 'r_ss': r_ss})

    # Linear function between 0 (bare soil) and the value at F=0.1
    r_s = r_s.where(F.gt(0.1), r_s1)
    r_s = r_s.where(T_s.subtract(T_c).abs().lt(1), r_s2)

    # Use "new" formula only for high DT values
    # Use "new" formula only for partial coverage (lai<3)
    r_s = r_s.where(F.gt(3), r_s2)

    # Computation of the canopy boundary layer resistance
    u_d = u.expression(
        'u_c * exp(-leaf_c * (1 - ((d0 + z0m) / hc)))',
        {'d0': d0, 'hc': hc, 'leaf_c': leaf_c, 'u_c': u_c, 'z0m': z0m})
    u_d = u_d.where(u_d.lte(0), 100)

    r_x = u.expression(
        'C / F * ((xl / u_d) ** 0.5)', {'C': C, 'F': F, 'u_d': u_d, 'xl': xl})
    r_x = r_x.where(u_d.eq(100), 0.1)

    return r_ah, r_s, r_x, u_attr


def compute_u_attr(u, d0, z0m, z_u, fm):
    """Friction Velocity

    Parameters
    ----------
    u : ee.Image
    d0
    z0m
    z_u
    fm

    Returns
    -------
    u_attr

    """
    u_attr = u.expression(
        '0.41 * u / ((log((z_u - d0) / z0m)) - fm)',
        {'d0': d0, 'fm': fm, 'u': u, 'z0m': z0m, 'z_u': z_u})
    u_attr = u_attr.where(u_attr.eq(0), 10)
    u_attr = u_attr.where(u_attr.lte(0), 0.01)
    return u_attr


def compute_r_ah(u_attr, d0, z0h, z_t, fh):
    """

    Parameters
    ----------
    u_attr : ee.Image
    d0
    z0h
    z_t
    fh

    Returns
    -------
    r_ah

    """
    r_ah = u_attr.expression(
        '(log((z_t - d0) / z0h) - fh) / u_attr / 0.41',
        {'d0': d0, 'fh': fh, 'u_attr': u_attr, 'z0h': z0h, 'z_t': z_t})
    r_ah = r_ah.where(r_ah.eq(0), 500)
    # r_ah = r_ah.max(1)
    r_ah = r_ah.where(r_ah.lte(1.0), 1.0)
    return r_ah


def compute_r_s(u_attr, T_s, T_c, hc, F, d0, z0m, leaf, leaf_s, fm_h):
    """

    Parameters
    ----------
    u_attr : ee.Image
    T_s : ee.Image
        Soil temperature (Kelvin).
    T_c : ee.Image
        Canopy temperature (Kelvin).
    hc : ee.Image
    F : ee.Image
        Input is LAI?
    d0
    z0m
    leaf
    leaf_s
    fm_h

    Returns
    -------
    r_s

    """
    # Free convective velocity constant for r_s modelling
    c_a = 0.004
    # Empirical constant for r_s modelling
    c_b = 0.012
    # Empirical constant for r_s modelling
    # (new formulation Kustas and Norman, 1999)
    c_c = 0.0025

    # Computation of the resistance of the air between soil and canopy space
    u_c = u_attr.expression(
        'u_attr / 0.41 * ((log((hc - d0) / z0m)) - fm_h)',
        {'d0': d0, 'fm_h': fm_h, 'hc': hc, 'u_attr': u_attr, 'z0m': z0m})
    u_c = u_c.where(u_c.lte(0), 0.1)
    u_s = u_attr.expression(
        'u_c * exp(-leaf * (1 - (0.05 / hc)))',
        {'hc': hc, 'leaf': leaf, 'u_c': u_c})

    r_ss = u_attr.expression(
        '1.0 / (c_a + (c_b * (u_c * exp(-leaf_s * (1.0 - (0.05 / hc))))))',
        {'c_a': c_a, 'c_b': c_b, 'hc': hc, 'leaf_s': leaf_s, 'u_c': u_c})
    r_s1 = T_s.expression(
        '1.0 / ((((abs(T_s - T_c)) ** (1.0 / 3.0)) * c_c) + (c_b * Us))',
        {'c_b': c_b, 'c_c': c_c, 'T_c': T_c, 'T_s': T_s, 'Us': u_s})
    r_s2 = u_attr.expression(
        '1.0 / (c_a + (c_b * Us))', {'c_a': c_a, 'c_b': c_b, 'Us': u_s})
    r_s = u_attr.expression(
        '(((r_ss - 1.0) / 0.09 * (F - 0.01)) + 1.0)', {'F': F, 'r_ss': r_ss})

    # Linear function between 0 (bare soil) and the value at F=0.1
    r_s = r_s.where(F.gt(0.1), r_s1)
    r_s = r_s.where(T_s.subtract(T_c).abs().lt(1), r_s2)

    # Use "new" formula only for high DT values
    # Use "new" formula only for partial coverage (lai<3)
    r_s = r_s.where(F.gt(3), r_s2)
    return r_s


def compute_r_x(u_attr, hc, F, d0, z0m, xl, leaf_c, fm_h):
    """

    Parameters
    ----------
    u_attr : ee.Image
    hc : ee.Image
    F : ee.Image
        Input is LAI?
    d0
    z0m
    xl
    leaf_c
    fm_h

    Returns
    -------
    r_x

    """

    # Parameter for canopy boundary-layer resistance
    # (C=90 Grace '81, C=175 Cheubouni 2001, 144 Li '98)
    C = 175.0

    # Computation of the resistance of the air between soil and canopy space
    u_c = u_attr.expression(
        'u_attr / 0.41 * ((log((hc - d0) / z0m)) - fm_h)',
        {'d0': d0, 'fm_h': fm_h, 'hc': hc, 'u_attr': u_attr, 'z0m': z0m})
    u_c = u_c.where(u_c.lte(0), 0.1)

    # Computation of the canopy boundary layer resistance
    u_d = u_attr.expression(
        'u_c * exp(-leaf_c * (1 - ((d0 + z0m) / hc)))',
        {'d0': d0, 'hc': hc, 'leaf_c': leaf_c, 'u_c': u_c, 'z0m': z0m})
    u_d = u_d.where(u_d.lte(0), 100)

    r_x = u_attr.expression(
        'C / F * ((xl / u_d) ** 0.5)', {'C': C, 'F': F, 'u_d': u_d, 'xl': xl})
    r_x = r_x.where(u_d.eq(100), 0.1)
    return r_x


def compute_Rn(albedo_c, albedo_s, T_air, T_c, T_s, e_atm, Rs_c, Rs_s, F):
    """Compute Soil and Canopy Net Radiation

    Parameters
    ----------
    albedo_c : ee.Image
    albedo_s : ee.Image
    T_air : ee.Image
        Air temperature (Kelvin).
    T_c : ee.Image
        Canopy temperature (Kelvin).
    T_s : ee.Image
        Soil temperature (Kelvin).
    e_atm : ee.Image
    Rs_c : ee.Image
    Rs_s : ee.Image
    F : ee.Image

    Returns
    -------
    Rn_s : ee.Image
        Soil net radiation (W m-2)
    Rn_c : ee.Image
        Canopy net radiation (W m-2)
    Rn : ee.Image
        Net radiation (W m-2)

    """
    # Long-wave extinction coefficient [-]
    kL = 0.95
    # Soil Emissivity [-]
    eps_s = 0.94
    # Canopy emissivity [-]
    eps_c = 0.99
      
    L_c = T_c.expression(
        'eps_c * 0.0000000567 * (T_c ** 4)', {'eps_c': eps_c, 'T_c': T_c})
    L_s = T_s.expression(
        'eps_s * 0.0000000567 * (T_s ** 4)', {'eps_s': eps_s, 'T_s': T_s})
    Rle = T_air.expression(
        'e_atm * 0.0000000567 * (T_air ** 4)',
        {'e_atm': e_atm, 'T_air': T_air})
    Rn_c = albedo_c.expression(
        '((1 - albedo_c) * Rs_c) + '
        '((1 - exp(-kL * F)) * (Rle + Ls - 2 * L_c))',
        {'albedo_c': albedo_c, 'F': F, 'kL': kL, 'L_c': L_c, 'L_s': L_s,
         'Rle': Rle, 'Rs_c': Rs_c})
    Rn_s = albedo_s.expression(
        '((1 - albedo_s) * Rs_s) + '
        '((exp(-kL * F)) * Rle) + ((1 - exp(-kL * F)) * L_c) - L_s',
        {'albedo_s': albedo_s, 'F': F, 'kL': kL, 'Lc': L_c, 'Ls': L_s,
         'Rle': Rle, 'Rs_s': Rs_s})
    Rn = Rn_s.expression(
        'Rn_s + Rn_c', {'Rn_s': Rn_s, 'Rn_c': Rn_c})

    return Rn_s, Rn_c, Rn


def compute_Rn_c(albedo_c, T_air, T_c, T_s, e_atm, Rs_c, F):
    """Compute Canopy Net Radiation

    Parameters
    ----------
    albedo_c : ee.Image
    T_air : ee.Image
        Air temperature (Kelvin).
    T_c : ee.Image
        Canopy temperature (Kelvin).
    T_s : ee.Image
        Soil temperature (Kelvin).
    e_atm : ee.Image
    Rs_c : ee.Image
    F : ee.Image

    Returns
    -------
    Rn_c : ee.Image

    """
    # Long-wave extinction coefficient [-]
    kL = 0.95
    # Soil Emissivity [-]
    eps_s = 0.94
    # Canopy emissivity [-]
    eps_c = 0.99

    # Stephan Boltzmann constant (W m-2 K-4)
    # sb = 5.670373e-8
    Lc = T_c.expression(
        'eps_c * 5.67E-8 * (T_c ** 4)', {'eps_c': eps_c, 'T_c': T_c})
    Ls = T_s.expression(
        'eps_s * 5.67E-8 * (T_s ** 4)', {'eps_s': eps_s, 'T_s': T_s})
    Rle = T_air.expression(
        'e_atm * 5.67E-8 * (T_air ** 4)',
        {'e_atm': e_atm, 'T_air': T_air})
    Rn_c = albedo_c.expression(
        '((1 - albedo_c) * Rs_c) + '
        '((1 - exp(-kL * F)) * (Rle + Ls - 2 * Lc))',
        {'albedo_c': albedo_c, 'F': F, 'kL': kL, 'Lc': Lc, 'Ls': Ls,
         'Rle': Rle, 'Rs_c': Rs_c})
    return Rn_c


def compute_Rn_s(albedo_s, T_air, T_c, T_s, e_atm, Rs_s, F):
    """Compute Soil Net Radiation

    Parameters
    ----------
    albedo_s : ee.Image
    T_air : ee.Image
        Air temperature (Kelvin).
    T_c : ee.Image
        Canopy temperature (Kelvin).
    T_s : ee.Image
        Soil temperature (Kelvin).
    e_atm : ee.Image
    Rs_c : ee.Image
    F : ee.Image

    Returns
    -------
    Rn_s : ee.Image

    """
    # Long-wave extinction coefficient [-]
    kL = 0.95
    # Soil Emissivity [-]
    eps_s = 0.94
    # Canopy emissivity [-]
    eps_c = 0.99

    L_c = T_c.expression(
        'eps_c * 0.0000000567 * (T_c ** 4)', {'eps_c': eps_c, 'T_c': T_c})
    L_s = T_s.expression(
        'eps_s * 0.0000000567 * (T_s ** 4)', {'eps_s': eps_s, 'T_s': T_s})
    Rle = T_air.expression(
        'e_atm * 0.0000000567 * (T_air ** 4)',
        {'e_atm': e_atm, 'T_air': T_air})
    Rn_s = albedo_s.expression(
        '((1 - albedo_s) * Rs_s) + '
        '((exp(-kL * F)) * Rle) + ((1 - exp(-kL * F)) * L_c) - L_s',
        {'albedo_s': albedo_s, 'F': F, 'kL': kL, 'L_c': L_c, 'L_s': L_s,
         'Rle': Rle, 'Rs_s': Rs_s})
    return Rn_s


# def temp_separation(H_c, fc, T_air, t0, r_ah, r_s, r_x, r_air, cp=1004.16):
#     """
#
#     Parameters
#     ----------
#     H_c : ee.Image
#     fc : ee.Image
#     T_air : ee.Image
#         Air temperature (Kelvin).
#     t0 :
#     r_ah : ee.Image
#     r_s : ee.Image
#         Soil aerodynamic resistance to heat transport (s m-1).
#     r_x : ee.Image
#         Bulk canopy aerodynamic resistance to heat transport (s m-1).
#     r_air : ee.Image
#     cp :
#
#     Returns
#     -------
#     T_c
#     T_s
#     Tac
#
#     """
#
#     T_c_lin = fc.expression(
#         '((T_air / r_ah) + '
#         ' (t0 / r_s / (1 - fc)) + '
#         ' (H_c * r_x / r_air / cp * ((1 / r_ah) + (1 / r_s) + (1 / r_x)))) / '
#         '((1 / r_ah) + (1 / r_s) + (fc / r_s / (1 - fc)))',
#         {'cp': cp, 'fc': fc, 'H_c': H_c, 'r_ah': r_ah, 'r_air': r_air,
#          'r_s': r_s, 'r_x': r_x, 't0': t0, 'T_air': T_air})
#
#     Td = fc.expression(
#         '(T_c_lin * (1 + (r_s / r_ah))) - '
#         '(H_c * r_x / r_air / cp * (1 + (r_s / r_x) + (r_s / r_ah))) - '
#         '(T_air * r_s / r_ah)',
#         {'cp': cp, 'H_c': H_c, 'r_ah': r_ah, 'r_air': r_air, 'r_s': r_s,
#          'r_x': r_x, 'T_air': T_air, 'T_c_lin': T_c_lin})
#
#     delta_T_c = fc.expression(
#         '((t0 ** 4) - (fc * (T_c_lin ** 4)) - ((1 - fc) * (Td ** 4))) / '
#         '((4 * (1 - fc) * (Td ** 3) * (1 + (r_s / r_ah))) + (4 * fc * (T_c_lin ** 3)))',
#         {'fc': fc, 'r_ah': r_ah, 'r_s': r_s, 't0': t0, 'Td': Td,
#          'T_c_lin': T_c_lin})
#
#     T_c = fc \
#         .expression(
#             'T_c_lin + delta_T_c', {'T_c_lin': T_c_lin, 'delta_T_c': delta_T_c}) \
#         .where(fc.lt(0.10), t0) \
#         .where(fc.gt(0.90), t0)
#
#     # Get T_s
#     Delta = fc.expression(
#         '(t0 ** 4) - (fc * (T_c ** 4))', {'fc': fc, 't0': t0, 'T_c': T_c})
#     Delta = Delta.where(Delta.lte(0), 10)
#
#     # CGM - This could probably be simplified
#     T_s = fc \
#         .expression('(Delta / (1 - fc)) ** 0.25', {'Delta': Delta, 'fc': fc}) \
#         .where(
#             fc.expression(
#                 '((t0 ** 4) - (fc * T_c ** 4)) <= 0.',
#                 {'fc': fc, 't0': t0, 'T_c': T_c}),
#             fc.expression(
#                 '(t0 - (fc * T_c)) / (1 - fc)',
#                 {'fc': fc, 't0': t0, 'T_c': T_c})) \
#         .where(fc.lt(0.1), t0) \
#         .where(fc.gt(0.9), t0)
#
#     T_c = T_c.where(T_c.lte(T_air.subtract(10.0)), T_air.subtract(10.0))
#     T_c = T_c.where(T_c.gte(T_air.add(50.0)), T_air.add(50.0))
#     T_s = T_s.where(T_s.lte(T_air.subtract(10.0)), T_air.subtract(10.0))
#     T_s = T_s.where(T_s.gte(T_air.add(50.0)), T_air.add(50.0))
#
#     T_ac = fc.expression(
#         '((T_air / r_ah) + (T_s / r_s) + (T_c / r_x)) / '
#         '((1 / r_ah) + (1 / r_s) + (1 / r_x))',
#         {'r_ah': r_ah, 'r_s': r_s, 'r_x': r_x, 'T_c': T_c, 'T_s': T_s,
#          'T_air': T_air})
#
#     return T_c, T_s, T_ac


def temp_separation_tc(H_c, fc, T_air, t0, r_ah, r_s, r_x, r_air, cp=1004.16):
    """Compute canopy temperature

    Parameters
    ----------
    H_c : ee.Image
    fc : ee.Image
    T_air : ee.Image
        Air temperature (Kelvin).
    t0 :
    r_ah : ee.Image
    r_s : ee.Image
        Soil aerodynamic resistance to heat transport (s m-1).
    r_x : ee.Image
        Bulk canopy aerodynamic resistance to heat transport (s m-1).
    r_air : ee.Image
    cp :

    Returns
    -------
    T_c : ee.Image
        Canopy temperature (Kelvin).

    """
    T_c_lin = fc.expression(
        '((T_air / r_ah) + '
        ' (t0 / r_s / (1 - fc)) + '
        ' (H_c * r_x / r_air / cp * ((1 / r_ah) + (1 / r_s) + (1 / r_x)))) / '
        '((1 / r_ah) + (1 / r_s) + (fc / r_s / (1 - fc)))',
        {'cp': cp, 'fc': fc, 'H_c': H_c, 'r_ah': r_ah, 'r_air': r_air,
         'r_s': r_s, 'r_x': r_x, 't0': t0, 'T_air': T_air})
    Td = fc.expression(
        '(T_c_lin * (1 + (r_s / r_ah))) - '
        '(H_c * r_x / r_air / cp * (1 + (r_s / r_x) + (r_s / r_ah))) - '
        '(T_air * r_s / r_ah)',
        {'cp': cp, 'H_c': H_c, 'r_ah': r_ah, 'r_air': r_air, 'r_s': r_s,
         'r_x': r_x, 'T_air': T_air, 'T_c_lin': T_c_lin})
    delta_T_c = fc.expression(
        '((t0 ** 4) - (fc * (T_c_lin ** 4)) - ((1 - fc) * (Td ** 4))) / '
        '((4 * (1 - fc) * (Td ** 3) * (1 + (r_s / r_ah))) + (4 * fc * (T_c_lin ** 3)))',
        {'fc': fc, 'r_ah': r_ah, 'r_s': r_s, 't0': t0, 'Td': Td,
         'T_c_lin': T_c_lin})
    T_c = fc \
        .expression(
        'T_c_lin + delta_T_c', {'T_c_lin': T_c_lin, 'delta_T_c': delta_T_c}) \
        .where(fc.lt(0.10), t0) \
        .where(fc.gt(0.90), t0)
    return T_c.max(T_air.subtract(10.0)).min(T_air.add(50.0))
    # T_c = T_c.where(T_c.lte(T_air.subtract(10.0)), T_air.subtract(10.0))
    # T_c = T_c.where(T_c.gte(T_air.add(50.0)), T_air.add(50.0))
    # return T_c


def temp_separation_ts(T_c, fc, T_air, t0):
    """Compute soil temperature

    Parameters
    ----------
    T_c : ee.Image
        Canopy temperature (Kelvin).
    fc : ee.Image
    T_air : ee.Image
        Air temperature (Kelvin).
    t0

    Returns
    -------
    T_s : ee.Image

    """
    Delta = fc.expression(
        '(t0 ** 4) - (fc * (T_c ** 4))', {'fc': fc, 't0': t0, 'T_c': T_c})
    Delta = Delta.where(Delta.lte(0), 10)

    # CGM - This could probably be simplified
    T_s = fc \
        .expression('(Delta / (1 - fc)) ** 0.25', {'Delta': Delta, 'fc': fc}) \
        .where(
            fc.expression(
                '((t0 ** 4) - (fc * T_c ** 4)) <= 0.',
                {'fc': fc, 't0': t0, 'T_c': T_c}),
            fc.expression(
                '(t0 - (fc * T_c)) / (1 - fc)',
                {'fc': fc, 't0': t0, 'T_c': T_c})) \
        .where(fc.lt(0.1), t0) \
        .where(fc.gt(0.9), t0)
    return T_s.max(T_air.subtract(10.0)).min(T_air.add(50.0))
    # T_s = T_s.where(T_s.lte(T_air.subtract(10.0)), T_air.subtract(10.0))
    # T_s = T_s.where(T_s.gte(T_air.add(50.0)), T_air.add(50.0))
    # return T_s


def temp_separation_tac(T_c, T_s, fc, T_air, r_ah, r_s, r_x):
    """Compute air temperature at the canopy interface

    Parameters
    ----------
    T_c : ee.Image
        Canopy temperature (Kelvin).
    T_s : ee.Image
        Soil temperature (Kelvin).
    fc
    T_air : ee.Image
    r_ah : ee.Image
    r_s : ee.Image
    r_x : ee.Image

    Returns
    -------
    T_ac : ee.Image
        Air temperature at the canopy interface (Kelvin).

    """
    T_ac = fc.expression(
        '((T_air / r_ah) + (T_s / r_s) + (T_c / r_x)) / '
        '((1 / r_ah) + (1 / r_s) + (1 / r_x))',
        {'r_ah': r_ah, 'r_s': r_s, 'r_x': r_x, 'T_c': T_c, 'T_s': T_s,
         'T_air': T_air})
    return T_ac


# CGM - z0h is not used in this function, should it be?
def compute_stability(H, t0, r_air, u_attr, z_u, z_t, hc, d0, z0m, z0h,
                      cp=1004.16):
    """

    Parameters
    ----------
    H
    t0
    r_air
    u_attr
    z_u
    z_t
    hc
    d0
    z0m
    z0h
    cp

    Returns
    -------
    fm
    fh
    fm_h

    """
    t0 = t0.where(t0.eq(0), 100)
    L_ob = H \
        .expression(
            '-(r_air * cp * t0 * (u_attr ** 3.0) / 0.41 / 9.806 / H)',
            {'cp': cp, 'H': H, 'r_air': r_air, 't0': t0, 'u_attr': u_attr})
    L_ob = L_ob.where(L_ob.gte(0), -99.0)

    L_mask = L_ob.eq(-99.0)
    mm = H \
        .expression(
            '((1 - (16.0 * (z_u - d0) / L_ob)) ** 0.25)',
            {'d0': d0, 'L_ob': L_ob, 'z_u': z_u}) \
        .where(L_mask, 0.0)
    mm_h = H \
        .expression(
            '((1 - (16.0 * (hc - d0) / L_ob)) ** 0.25)',
            {'d0': d0, 'hc': hc, 'L_ob': L_ob}) \
        .where(L_mask, 0.0)
    mh = H \
        .expression(
            '((1 - (16.0 * (z_t - d0) / L_ob)) ** 0.25)',
            {'d0': d0, 'L_ob': L_ob, 'z_t': z_t}) \
        .where(L_mask, 0.0)

    L_mask = L_ob.gt(-100).And(L_ob.lt(100))
    # CGM - Intentionally using mh here (instead of mm_h) to match Python code
    fm = mh.where(
        L_mask,
        H.expression(
            '2.0 * log((1.0 + mm) / 2.0) + log((1.0 + (mm ** 2)) / 2.0) - '
            '2.0 * atan(mm) + (pi / 2)',
            {'mm': mm, 'pi': math.pi}))
    fm_h = mh.where(
        L_mask,
        H.expression(
            '2.0 * log((1.0 + mm_h) / 2.0) + log((1.0 + (mm_h ** 2)) / 2.0) - '
            '2.0 * atan(mm_h) + (pi / 2)',
            {'mm_h': mm_h, 'pi': math.pi}))
    fh = mh.where(
        L_mask, H.expression(
            '(2.0 * log((1.0 + (mh ** 2.0)) / 2.0))', {'mh': mh}))

    # CGM - Swapped order of calc since d0 is an image compute from hc and
    #   z_u is being set as a constant number (for now).
    fm = fm.where(fm.eq(d0.multiply(-1).add(z_u).divide(z0m).log()), fm.add(1.0))
    # fm = fm.where(fm.eq(z_u.subtract(d0).divide(z0m).log()), fm.add(1.0))
    fm_h = fm_h.where(fm_h.eq(hc.subtract(d0).divide(z0m).log()), fm_h.add(1.0))
    # fm_h = fm_h.where(fm_h.eq(hc.subtract(d0).divide(z0m).log()), fm_h.add(1.0))

    return fm, fh, fm_h


def compute_stability_fh(H, t0, u_attr, r_air, z_t, d0, cp=1004.16):
    """

    To Match IDL, fh defaults to 0 for L values outside the range -100 to 100
    It is unclear why the range -100 to 100 is used since L is forced to be < 0

    Parameters
    ----------
    H : ee.Image
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
    L_ob = H .expression(
        '-(r_air * cp * t0 * (u_attr ** 3.0) / 0.41 / 9.806 / H)',
        {'cp': cp, 'H': H, 'r_air': r_air, 't0': t0, 'u_attr': u_attr})
    L_ob = L_ob.where(L_ob.gte(0), -99)
    mh = H \
        .expression(
            '((1 - (16.0 * (z_t - d0) / L_ob)) ** 0.25)',
            {'d0': d0, 'L_ob': L_ob, 'z_t': z_t}) \
        .where(L_ob.eq(-99), 0.0)
    fh = H \
        .expression('(2.0 * log((1.0 + (mh ** 2.0)) / 2.0))', {'mh': mh}) \
        .where(L_ob.lte(-100).Or(L_ob.gte(100)), 0)

    return fh


def compute_stability_fm(H, t0, u_attr, r_air, z_u, d0, z0m, cp=1004.16):
    """

    To Match IDL, fh defaults to 0 for L values outside the range -100 to 100
    It is unclear why the range -100 to 100 is used since L is forced to be < 0

    Parameters
    ----------
    H : ee.Image
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
    L_ob = H.expression(
        '-(r_air * cp * t0 * (u_attr ** 3.0) / 0.41 / 9.806 / H)',
        {'cp': cp, 'H': H, 'r_air': r_air, 't0': t0, 'u_attr': u_attr})
    L_ob = L_ob.where(L_ob.gte(0), -99.0)
    mh = H \
        .expression(
            '((1 - (16.0 * (z_u - d0) / L_ob)) ** 0.25)',
            {'d0': d0, 'L_ob': L_ob, 'z_u': z_u}) \
        .where(L_ob.eq(-99.0), 0.0)
    fm = H \
        .expression(
            '2.0 * log((1.0 + mh) / 2.0) + log((1.0 + (mh ** 2)) / 2.0) - '
            '2.0 * atan(mh) + (pi / 2)',
            {'mh': mh, 'pi': math.pi}) \
        .where(L_ob.lte(-100).Or(L_ob.gte(100)), 0)

    # CGM - Swapped order of calc since d0 is an image compute from hc and
    #   z_u is being set as a constant number (for now).
    fm = fm.where(fm.eq(d0.multiply(-1).add(z_u).divide(z0m).log()), fm.add(1.0))
    # fm = fm.where(fm.eq(z_u.subtract(d0).divide(z0m).log()), fm.add(1.0))
    return fm


def compute_stability_fm_h(H, t0, u_attr, r_air, hc, d0, z0m, cp=1004.16):
    """

    To Match IDL, fh defaults to 0 for L values outside the range -100 to 100
    It is unclear why the range -100 to 100 is used since L is forced to be < 0

    Parameters
    ----------
    H : ee.Image
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
    L_ob = H.expression(
        '-(r_air * cp * t0 * (u_attr ** 3.0) / 0.41 / 9.806 / H)',
        {'cp': cp, 'H': H, 'r_air': r_air, 't0': t0, 'u_attr': u_attr})
    L_ob = L_ob.where(L_ob.gte(0), -99.0)
    mm_h = H \
        .expression(
            '((1 - (16.0 * (hc - d0) / L_ob)) ** 0.25)',
            {'d0': d0, 'hc': hc, 'L_ob': L_ob}) \
        .where(L_ob.eq(-99.0), 0.0)
    fm_h = H \
        .expression(
            '2.0 * log((1.0 + mm_h) / 2.0) + log((1.0 + (mm_h ** 2)) / 2.0) - '
            '2.0 * atan(mm_h) + (pi / 2)',
            {'mm_h': mm_h, 'pi': math.pi}) \
        .where(L_ob.lte(-100).Or(L_ob.gte(100)), 0)

    # CGM - Swapped order of calc since d0 is an image compute from hc and
    #   z_u is being set as a constant number (for now).
    fm_h = fm_h.where(fm_h.eq(hc.subtract(d0).divide(z0m).log()), fm_h.add(1.0))
    # fm_h = fm_h.where(fm_h.eq(hc.subtract(d0).divide(z0m).log()), fm_h.add(1.0))
    return fm_h
