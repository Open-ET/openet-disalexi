// DisALEXI - Test building ALEXI scale air temperature assets

var pi = Math.PI;
var rad2deg = 180.0 / Math.PI;
var deg2rad = Math.PI / 180.0;


// DisALEXI Image class
function disalexi_image(image) {

  // Variables from Image Class
  var iterations = 10;
  var t_air_values = range(273, 321, 1);

  // Set server side date/time properties using the 'system:time_start'
  var image_dt = ee.Date(ee.Image(image).get('system:time_start'));
  var image_date = ee.Date(image_dt.format('yyyy-MM-dd'));
  var doy = ee.Number(image_dt.getRelative('day', 'year')).add(1).double();
  var hour = ee.Number(image_dt.getFraction('day')).multiply(24);
  var hour_int = hour.floor();
  // Time used in IDL is hours and fractional minutes (no seconds)
  var time = ee.Date(image_dt).get('hour').add(
    ee.Date(image_dt).get('minute').divide(60));

  // CGM - Applying cloud mask directly to input image
  //   instead of to a_pt in main TSEB function
  var cfmask = ee.Image(image).select('cfmask');
  var mask = cfmask.eq(0);
  var input_image = ee.Image(image).updateMask(mask);

  // Get input bands from the image
  var albedo = input_image.select('albedo');
  var lai = input_image.select('lai');
  // lai = lai.where(lai.mask(), 0.01);
  var lst = input_image.select('lst');
  var ndvi = input_image.select('ndvi');


  // Elevation [m]
  var elevation = ee.Image('USGS/NED').rename(['elevation']);
  // elevation = ee.Image('USGS/SRTMGL1_003').rename(['elevation']);
  // elevation = ee.Image('USGS/GTOPO30').rename(['elevation']);
  // elevation = ee.Image('USGS/GMTED2010').rename(['elevation']);
  // elevation = ee.Image('CGIAR/SRTM90_V4').rename(['elevation']);
  // elevation = ee.Image.constant(350.0).rename(['elevation']);


  // Set default land cover image and type
  // For now default to CONUS and use default if image and type were not set
  // GlobeLand30 values need to be set to the lowest even multiple of 10,
  //   since that is currently what is in the landcover.xlsx file.
  // http://www.globallandcover.com/GLC30Download/index.aspx
  // Using NLCD as default land cover and type
  var landcover = ee.Image('USGS/NLCD/NLCD2011').select(['landcover']);
  var lc_type = 'NLCD';


  // ALEXI ET - CONUS
  var et_coll = ee.ImageCollection('projects/disalexi/alexi/CONUS');
  var et_transform = [0.04, 0, -125.04, 0, -0.04, 49.82];
  var et_crs = 'EPSG:4326';

  // Hard coding using CFSR for wind speed
  var windspeed_coll = ee.ImageCollection('NOAA/CFSV2/FOR6H')
    .select([
      'u-component_of_wind_height_above_ground',
      'v-component_of_wind_height_above_ground']);

  // Hard coding using MERRA2 solar insolation [W m-2]
  var rs_hourly_coll = ee.ImageCollection('projects/climate-engine/merra2/rs_hourly');
  var rs_daily_coll = ee.ImageCollection('projects/climate-engine/merra2/rs_daily');


  // def _set_alexi_et_vars(self):
  // Extract ALEXI ET image for the target image time
  var alexi_et = ee.Image(ee.ImageCollection(et_coll)
    .filterDate(image_date, image_date.advance(1, 'day')).first());
  alexi_et = alexi_et.rename(['alexi_et']);


  // def _set_elevation_vars(self):
  // Compute elevation derived variables
  var pressure = elevation
    .expression(
      '101.3 * (((293.0 - 0.0065 * z) / 293.0) ** 5.26)',
      {'z': elevation})
    .rename(['pressure']);

  var remaps = {
    'NLCD': {
      'aleafv': {
        11: 0.82, 12: 0.82, 21: 0.84, 22: 0.84, 23: 0.84, 24: 0.84,
        31: 0.82, 32: 0.82, 41: 0.86, 42: 0.89, 43: 0.87,
        51: 0.83, 52: 0.83, 71: 0.82, 72: 0.82, 73: 0.82, 74: 0.82,
        81: 0.82, 82: 0.83, 90: 0.85, 91: 0.85, 92: 0.85, 93: 0.85,
        94: 0.85, 95: 0.85, 96: 0.85, 97: 0.85, 98: 0.85, 99: 0.85},
      'aleafn': {
        11: 0.28, 12: 0.28, 21: 0.37, 22: 0.37, 23: 0.37, 24: 0.37,
        31: 0.57, 32: 0.57, 41: 0.37, 42: 0.60, 43: 0.48,
        51: 0.35, 52: 0.35, 71: 0.28, 72: 0.28, 73: 0.28, 74: 0.28,
        81: 0.28, 82: 0.35, 90: 0.36, 91: 0.36, 92: 0.36, 93: 0.36,
        94: 0.36, 95: 0.36, 96: 0.36, 97: 0.36, 98: 0.36, 99: 0.36},
      'aleafl': {
        11: 0.95, 12: 0.95, 21: 0.95, 22: 0.95, 23: 0.95, 24: 0.95,
        31: 0.95, 32: 0.95, 41: 0.95, 42: 0.95, 43: 0.95,
        51: 0.95, 52: 0.95, 71: 0.95, 72: 0.95, 73: 0.95, 74: 0.95,
        81: 0.95, 82: 0.95, 90: 0.95, 91: 0.95, 92: 0.95, 93: 0.95,
        94: 0.95, 95: 0.95, 96: 0.95, 97: 0.95, 98: 0.95, 99: 0.95},
      'adeadv': {
        11: 0.42, 12: 0.42, 21: 0.58, 22: 0.58, 23: 0.58, 24: 0.58,
        31: 0.92, 32: 0.92, 41: 0.84, 42: 0.84, 43: 0.84,
        51: 0.77, 52: 0.77, 71: 0.42, 72: 0.42, 73: 0.42, 74: 0.42,
        81: 0.42, 82: 0.49, 90: 0.58, 91: 0.58, 92: 0.58, 93: 0.58,
        94: 0.58, 95: 0.58, 96: 0.58, 97: 0.58, 98: 0.58, 99: 0.58},
      'adeadn': {
        11: 0.04, 12: 0.04, 21: 0.26, 22: 0.26, 23: 0.26, 24: 0.26,
        31: 0.80, 32: 0.80, 41: 0.61, 42: 0.61, 43: 0.61,
        51: 0.52, 52: 0.52, 71: 0.04, 72: 0.04, 73: 0.04, 74: 0.04,
        81: 0.04, 82: 0.13, 90: 0.26, 91: 0.26, 92: 0.26, 93: 0.26,
        94: 0.26, 95: 0.26, 96: 0.26, 97: 0.26, 98: 0.26, 99: 0.26},
      'adeadl': {
        11: 0.95, 12: 0.95, 21: 0.95, 22: 0.95, 23: 0.95, 24: 0.95,
        31: 0.95, 32: 0.95, 41: 0.95, 42: 0.95, 43: 0.95,
        51: 0.95, 52: 0.95, 71: 0.95, 72: 0.95, 73: 0.95, 74: 0.95,
        81: 0.95, 82: 0.95, 90: 0.95, 91: 0.95, 92: 0.95, 93: 0.95,
        94: 0.95, 95: 0.95, 96: 0.95, 97: 0.95, 98: 0.95, 99: 0.95},
      'hmin': {
        11: 0.00, 12: 0.00, 21: 0.00, 22: 0.00, 23: 1.00, 24: 6.00,
        31: 0.00, 32: 0.00, 41: 10.0, 42: 15.0, 43: 12.0,
        51: 0.20, 52: 1.00, 71: 0.10, 72: 0.10, 73: 0.10, 74: 0.10,
        81: 0.10, 82: 0.00, 90: 5.00, 91: 1.00, 92: 1.00, 93: 1.00,
        94: 1.00, 95: 1.00, 96: 1.00, 97: 1.00, 98: 1.00, 99: 1.00},
      'hmax': {
        11: 0.60, 12: 0.60, 21: 0.60, 22: 0.60, 23: 1.00, 24: 6.00,
        31: 0.20, 32: 0.20, 41: 10.0, 42: 15.0, 43: 12.0,
        51: 0.20, 52: 1.00, 71: 0.60, 72: 0.60, 73: 0.10, 74: 0.10,
        81: 0.60, 82: 0.60, 90: 5.00, 91: 2.50, 92: 2.50, 93: 2.50,
        94: 2.50, 95: 2.50, 96: 2.50, 97: 2.50, 98: 2.50, 99: 2.50},
      'xl': {
        11: 0.02, 12: 0.02, 21: 0.02, 22: 0.02, 23: 0.02, 24: 0.02,
        31: 0.02, 32: 0.02, 41: 0.10, 42: 0.05, 43: 0.08,
        51: 0.02, 52: 0.02, 71: 0.02, 72: 0.02, 73: 0.02, 74: 0.02,
        81: 0.02, 82: 0.05, 90: 0.05, 91: 0.05, 92: 0.05, 93: 0.05,
        94: 0.05, 95: 0.05, 96: 0.05, 97: 0.05, 98: 0.05, 99: 0.05},
      'omega': {
        11: 0.99, 12: 0.99, 21: 0.99, 22: 0.99, 23: 0.99, 24: 0.99,
        31: 0.99, 32: 0.99, 41: 0.78, 42: 0.68, 43: 0.75,
        51: 0.84, 52: 0.84, 71: 0.83, 72: 0.83, 73: 0.83, 74: 0.83,
        81: 0.83, 82: 0.83, 90: 0.86, 91: 0.86, 92: 0.86, 93: 0.86,
        94: 0.86, 95: 0.86, 96: 0.86, 97: 0.86, 98: 0.86, 99: 0.86}
      },
  };

  // def _set_landcover_vars(self):
  // Compute Land Cover / LAI derived variables
  function lc_remap(landcover, lc_type, lc_var) {
    // Remap land cover values to target values
    // Scale output values by 100 since remap values must be integer
    var lc_items = remaps[lc_type][lc_var];
    var input_values = [];
    var output_values = [];
    for (var key in lc_items) {
      if (lc_items.hasOwnProperty(key)) {
        input_values.push(parseInt(key, 10));
        output_values.push(parseInt(lc_items[key] * 100, 10));
      }
    }
    return ee.Image(landcover)
      .remap(input_values, output_values)
      .divide(100)
      .rename([lc_var]);
  }

  // Get LC based variables
  var aleafv = lc_remap(landcover, lc_type, 'aleafv');
  var aleafn = lc_remap(landcover, lc_type, 'aleafn');
  var aleafl = lc_remap(landcover, lc_type, 'aleafl');
  var adeadv = lc_remap(landcover, lc_type, 'adeadv');
  var adeadn = lc_remap(landcover, lc_type, 'adeadn');
  var adeadl = lc_remap(landcover, lc_type, 'adeadl');
  var hc_min = lc_remap(landcover, lc_type, 'hmin');
  var hc_max = lc_remap(landcover, lc_type, 'hmax');
  var leaf_width = lc_remap(landcover, lc_type, 'xl');
  var clump = lc_remap(landcover, lc_type, 'omega');

  // LAI for leafs spherical distribution
  var F = lai.multiply(clump).rename(['F']);

  // Fraction cover at nadir (view=0)
  var f_c = lai.expression('1.0 - exp(-0.5 * F)', {'F': F})
    .clamp(0.01, 0.9)
    .rename(['f_c']);

  // ======================================================================
  // Compute canopy height and roughness parameters
  var hc = lai
    .expression(
      'hc_min + ((hc_max - hc_min) * f_c)',
      {'hc_min': hc_min, 'hc_max': hc_max, 'f_c': f_c})
    .rename(['hc']);


  // function _set_solar_vars(self, interpolate_flag=True)
  // Extract MERRA2 solar images for the target image time

  // Interpolate rs hourly image at image time
  // Hourly Rs is time average so time starts are 30 minutes early
  // Move image time 30 minutes earlier to simplify filtering/interpolation
  // This interpolation scheme will only work for hourly data
  // if (interpolate_flag) {
  var interp_dt = image_dt.advance(-0.5, 'hour');
  var rs_a_img = ee.Image(rs_hourly_coll
    .filterDate(interp_dt.advance(-1, 'hour'), interp_dt).first());
  var rs_b_img = ee.Image(rs_hourly_coll
    .filterDate(interp_dt, interp_dt.advance(1, 'hour')).first());
  var t_a = ee.Number(rs_a_img.get('system:time_start'));
  var t_b = ee.Number(rs_b_img.get('system:time_start'));
  var rs1 = rs_b_img.subtract(rs_a_img)
    .multiply(interp_dt.millis().subtract(t_a).divide(t_b.subtract(t_a)))
    .add(rs_a_img)
    .rename(['rs']);
  // } else {
  //   var rs1 = ee.Image(
  //     ee.ImageCollection(rs_hourly_coll)
  //       .filterDate(image_date, image_date.advance(1, 'day'))
  //       .filter(ee.Filter.calendarRange(hour_int, hour_int, 'hour'))
  //       .first())
  //     .rename(['rs']);
  // }
  var rs24 = ee.Image(
    ee.ImageCollection(rs_daily_coll)
      .filterDate(image_date, image_date.advance(1, 'day'))
      .first())
    .rename(['rs']);


  // def _set_weather_vars(self):
  // Compute weather derived variables (only wind from CFSv2 for now)
  var windspeed_img = ee.Image(
    windspeed_coll.filterDate(image_date, image_date.advance(1, 'day')).mean());
  var windspeed = windspeed_img
    .expression('sqrt(b(0) ** 2 + b(1) ** 2)')
    .rename(['windspeed']);


  // def _set_time_vars(self):
  // Sub functions in tseb_utils.py
  var lon = ee.Image.pixelLonLat().select(['longitude']).multiply(pi / 180);
  var lat = ee.Image.pixelLonLat().select(['latitude']).multiply(pi / 180);

  // function _solar_time(date, lon, lat):
  // Computes solar time variables following Campbell & Norman 1998
  // IDL is computing time_t as hours and fractional minutes
  var time_t = image_dt.get('hour').add(image_dt.get('minute').divide(60));

  // This will return the hour floating point value
  // time_t = image_date.get('hour').add(image_date.getFraction('hour'))

  // CGM - DOY and hour could be images in order to make these expressions
  var a = image_date.get('month').multiply(-1).add(14).divide(12).floor();
  var y = image_date.get('year').add(4800).subtract(a);
  var m = image_date.get('month').add(a.multiply(12)).subtract(3);
  var julian = image_date.get('day')
    .add(m.multiply(153).add(2).divide(5).floor())
    .add(y.multiply(365))
    .add(y.divide(4.0).floor())
    .subtract(y.divide(100.0).floor())
    .add(y.divide(400.0).floor())
    .subtract(32045);

  // Sunrise time
  var julian_ = time_t.divide(24.0).add(julian);
  var j_cen = julian_.add(0.5 - 2451545.0).divide(36525.0);
  // CGM - Does the mod happen before or after the multiply
  var lon_sun = j_cen.multiply(0.0003032).add(36000.76983)
    .multiply(j_cen).mod(360.0).add(280.46646).subtract(360);
  var an_sun = j_cen.multiply(-0.0001537).add(35999.05029)
    .multiply(j_cen).add(357.52911);
  var ecc = j_cen.multiply(0.0000001267).add(0.000042037)
    .multiply(j_cen).multiply(-1).add(0.016708634);
  var ob_ecl = j_cen.multiply(-0.001813).add(0.00059)
    .multiply(j_cen).add(46.815)
    .multiply(j_cen).multiply(-1).add(21.448)
    .divide(60.0).add(26).divide(60).add(23);
  var ob_corr = j_cen.multiply(-1934.136).add(125.04).multiply(deg2rad).cos()
    .multiply(0.00256).add(ob_ecl);
  var var_y = ob_corr.divide(2.0).multiply(deg2rad).tan()
    .multiply(ob_corr.divide(2.0).multiply(deg2rad).tan());
  var eq_t = (
    lon_sun.multiply(2.0).multiply(deg2rad).sin().multiply(var_y)
      .subtract(an_sun.multiply(deg2rad).sin().multiply(ecc).multiply(2.0))
      .add(an_sun.multiply(deg2rad).sin()
        .multiply(lon_sun.multiply(2.0).multiply(deg2rad).cos())
        .multiply(var_y).multiply(ecc).multiply(4.0))
      .subtract(lon_sun.multiply(4.0).multiply(deg2rad).sin().multiply(var_y).multiply(var_y).multiply(0.5))
      .subtract(an_sun.multiply(2.0).multiply(deg2rad).sin().multiply(ecc).multiply(ecc).multiply(1.25))
    .multiply(4.0).multiply(rad2deg));
  var sun_eq = an_sun.multiply(deg2rad).sin()
    .multiply(j_cen.multiply(0.000014).add(0.004817).multiply(j_cen).multiply(-1).add(1.914602))
    .add(an_sun.multiply(2.0).multiply(deg2rad).sin().multiply(j_cen.multiply(-0.000101).add(0.019993)))
    .add(an_sun.multiply(3.0).multiply(deg2rad).sin().multiply(0.000289));
  var sun_true = sun_eq.add(lon_sun);
  var sun_app = j_cen.multiply(-1934.136).add(125.04).multiply(deg2rad).sin()
    .multiply(-0.00478).subtract(0.00569).add(sun_true);

  // CGM - Intentionally not converting back to degrees
  var d = ob_corr.multiply(deg2rad).sin()
    .multiply(sun_app.multiply(deg2rad).sin())
    .asin();

  // CGM - Functions below are lat/lon dependent and can be written as ee.Image expressions
  // CGM - d is still in radians, not converting
  var ha_t = lat.expression(
    'acos((cos(90.833 * pi / 180) / (cos(lat) * cos(d))) - tan(lat) * tan(d)) * (180 / pi)',
    {'lat': lat, 'd': d, 'pi': pi});

  // def sunrise_sunset(date, lon, lay):
  // Computes sunrise/sunset times
  // Adjust image date time to start of day
  var t_noon = lon.expression(
    '(720.0 - 4 * lon - eq_t) / 1440 * 24.0',
    {'lon': lon.multiply(rad2deg), 'eq_t': eq_t});
  var t_rise = lat
    .expression(
      '((t_noon / 24.0) - (ha_t * 4.0 / 1440)) * 24.0',
      {'t_noon': t_noon, 'ha_t': ha_t})
    .rename(['t_rise']);
  var t_end = lat
    .expression(
      '((t_noon / 24.0) + (ha_t * 4.0 / 1440)) * 24.0',
      {'t_noon': t_noon, 'ha_t': ha_t})
    .rename(['t_end']);

  // Compute time and position related variables
  // CGM - The zs returned by this function is not used in the original Python code.
  // The hour in the original call was the integer hour, even though it
  //   seems like it should be the float hour.

  // def solar_zenith(datetime, lon, lat):
  // Computes zenith angle
  // IDL is computing time_t as hours and fractional minutes (no seconds)
  time_t = image_dt.get('hour').add(image_dt.get('minute').divide(60));
  // # This will return the hour floating point value
  // time_t = image_dt.get('hour').add(image_dt.getFraction('hour'))

  var ts_time = lon.expression(
    '(time_t / 24.0 * 1440 + eq_t + 4.0 * lon) % 1440.',
    {'lon': lon.multiply(rad2deg), 'time_t': time_t, 'eq_t': eq_t});
  ts_time = ts_time.where(ts_time.gt(1440), ts_time.subtract(1440));
  // ts_time[ts_time > 1440.] = ts_time[ts_time > 1440.] - 1440.

  var w = lon.expression('ts_time / 4.0 + 180.0', {'ts_time': ts_time});
  w = w.where(ts_time.divide(4).gt(0), ts_time.divide(4).subtract(180));
  // w[ts_time/4.0 >= 0] = ts_time[ts_time/4.0 >= 0.] / 4.-180.

  var sol_zenith = lat
    .expression(
      'acos( (sin(lat) * sin(d)) + (cos(lat) * cos(d) * cos(w)) )',
      {'lat': lat, 'd': d, 'w': w.multiply(deg2rad)})
    .rename(['zs']);


  function t_air_func (t_air) {
    // Compute TSEB ET for each T_air value
    // Assume the function is being mapped over a FC with a "t_air" property
    // Return bias (ALEXI ET - TSEB ET) as first band for quality mosaic

    // Apply T_air values to Landsat pixels
    var t_air_img = lst.multiply(0).add(ee.Number(t_air.get('t_air')));
    // t_air_img = ee.Image.constant(ee.Number(t_air.get('t_air'))).double();

    var et = tseb_pt(
      t_air_img, lst, windspeed, pressure, elevation, rs1, rs24, 0, sol_zenith,
      aleafv, aleafn, aleafl, adeadv, adeadn, adeadl, albedo, ndvi, lai, clump, hc,
      time, t_rise, t_end, leaf_width, 1.32, iterations);

    // Invert the abs(bias) since quality mosaic sorts descending
    var bias = ee.Image(et).subtract(alexi_et).abs().multiply(-1);
    // bias = et_img.subtract(ET).abs().multiply(-1)
    return ee.Image([bias, ee.Image(et), t_air_img])
      .rename(['bias', 'et', 't_air']);
  }


  // Get output for a range of Tair values
  var t_air_ftr_list = [];
  for (var i = 0; i < t_air_values.length; i++) {
    t_air_ftr_list.push(ee.Feature(null, {'t_air': t_air_values[i]}));
  }
  var t_air_ftr_coll = ee.FeatureCollection(t_air_ftr_list);

  // Mapping over the list seemed a little slower than the FC
  // t_air_ftr_coll = ee.List.sequence(275, 335, 5)

  var output_coll = ee.ImageCollection(t_air_ftr_coll.map(t_air_func));

  // Return the air temperature associated with the lowest difference
  //   in ET from the ALEXI ET value
  var t_air = ee.Image(output_coll.qualityMosaic('bias'))
    .select(['t_air']);

  // Project the output back to the landsat scene
  //     .reproject(crs=, crsTransform=) \
  //     .rename(['t_air'])

  return t_air;
}




function tseb_pt(T_air, T_rad, u, p, z, Rs_1, Rs24, vza, zs,
                 aleafv, aleafn, aleafl, adeadv, adeadn, adeadl,
                 albedo, ndvi, lai, clump, hc, time, t_rise, t_end,
                 leaf_width, a_PT_in, iterations) {
  // Correct Clumping Factor
  var f_green = 1.0;

  // LAI for leaf spherical distribution
  var F = lai.expression('lai * clump', {'lai': lai, 'clump': clump});

  // Fraction cover at nadir (view=0)
  var fc = F.expression('1.0 - exp(-0.5 * F)', {'F': F})
    .clamp(0.01, 0.9);

  // LAI relative to canopy projection only
  var lai_c = lai.expression('lai / fc', {'lai': lai, 'fc': fc});

  // Houborg modification (according to Anderson et al. 2005)
  var fc_q = lai
    .expression('1 - (exp(-0.5 * F / cos(vza)))', {'F': F, 'vza': vza})
    .clamp(0.05, 0.90);

  // Brutsaert (1982)
  var z0m = hc.expression('hc * 0.123', {'hc': hc});
  // CGM - add(0) is to mimic numpy copy, check if needed
  var z0h = z0m.add(0);
  var d_0 = hc.expression('hc * (2.0 / 3.0)', {'hc': hc});

  // Correction of roughness parameters for bare soils (F < 0.1)
  d_0 = d_0.where(F.lte(0.1), 0.00001);
  z0m = z0m.where(F.lte(0.1), 0.01);
  z0h = z0h.where(F.lte(0.1), 0.0001);

  // Correction of roughness parameters for water bodies
  // (NDVI < 0 and albedo < 0.05)
  var water_mask = ndvi.lte(0).and(albedo.lte(0.05));
  d_0 = d_0.where(water_mask, 0.00001);
  z0m = z0m.where(water_mask, 0.00035);
  z0h = z0h.where(water_mask, 0.00035);

  // Check to avoid division by 0 in the next computations
  z0h = z0h.where(z0h.eq(0), 0.001);
  z0m = z0m.where(z0m.eq(0), 0.01);

  // DEADBEEF
  // z_u = ee.Number(50.0)
  // z_t = ee.Number(50.0)
  var z_u = ee.Image.constant(50.0);
  var z_t = ee.Image.constant(50.0);
  // z_u = lai.multiply(0).add(50)
  // z_t = lai.multiply(0).add(50)

  // Parameters for In-Canopy Wind Speed Extinction
  var leaf = lai.expression(
    '(0.28 * (F ** (0.66667)) * (hc ** (0.33333)) * (leaf_width ** (-0.33333)))',
    {'F': F, 'hc': hc, 'leaf_width': leaf_width});
  var leaf_c = lai.expression(
    '(0.28 * (lai_c ** (0.66667)) * (hc ** (0.33333)) * (leaf_width ** (-0.33333)))',
    {'lai_c': lai_c, 'hc': hc, 'leaf_width': leaf_width});
  var leaf_s = lai.expression(
    '(0.28 * (0.1 ** (0.66667)) * (hc ** (0.33333)) * (leaf_width ** (-0.33333)))',
    {'hc': hc, 'leaf_width': leaf_width});

  // ************************************************************************
  // Atmospheric Parameters
  // Saturation vapour pressure [kPa] (FAO56 3-8)
  var e_s = T_air.expression(
    '0.6108 * exp((17.27 * (T_air - 273.16)) / ((T_air - 273.16) + 237.3))',
    {'T_air': T_air});
  // Slope of the saturation vapor pressure [kPa] (FAO56 3-9)
  var Ss = T_air.expression(
    '4098. * e_s / (((T_air - 273.16) + 237.3) ** 2)',
    {'e_s': e_s, 'T_air': T_air});
  // Latent heat of vaporization (~2.45 at 20 C) [MJ kg-1] (FAO56 3-1)
  var lambda1 = T_air.expression(
    '(2.501 - (2.361e-3 * (T_air - 273.16)))',
    {'T_air': T_air});
  // Psychrometric constant [kPa C-1] (FAO56 3-10)
  var g = p.expression('1.615E-3 * p / lambda1', {'p': p, 'lambda1': lambda1});

  // ************************************************************************
  // Initialization of
  a_PT = albedo.multiply(0).add(a_PT_in);
  // a_PT = ee.Image.constant(a_PT_in)
  // a_PT = mask.multiply(a_PT)

  // CGM - This was also being computed inside albedo_separation function below
  // Commented out from here for now.
  // e_atm = T_air.expression(
  //     '1.0 - (0.2811 * (exp(-0.0003523 * ((T_air - 273.16) ** 2))))',
  //     {'T_air': T_air})


  // #####################################################################
  // def albedo_separation(albedo, Rs_1, F, fc, aleafv, aleafn, aleafl, adeadv,
  //                       adeadn, adeadl, zs, iterations=10):
  // Compute Solar Components and atmospheric properties ([Campbell1998]_)

  // Correct for curvature of atmos in airmas
  var airmas = zs.expression(
    '(sqrt(cos(zs) ** 2 + 0.0025) - cos(zs)) / 0.00125', {'zs': zs});

  // Correct for refraction(good up to 89.5 deg.)
  airmas = airmas.where(
    zs.multiply(rad2deg).lt(89.5),
    zs.expression(
      'airmas - (2.8 / (90.0 - zs_temp) ** 2)',
      {'airmas': airmas, 'zs_temp': zs.multiply(rad2deg)}));

  var potbm1 = zs.expression(
    '600.0 * exp(-0.160 * airmas)', {'airmas': airmas});
  var potvis = zs.expression(
    '(potbm1 + (600.0 - potbm1) * 0.4) * cos(zs)',
    {'potbm1': potbm1, 'zs': zs});
  // CGM - Not used
  var potdif = zs.expression(
    '(600.0 - potbm1) * 0.4 * cos(zs)', {'potbm1': potbm1, 'zs': zs});
  var uu = zs.expression('1.0 / cos(zs)', {'zs': zs})
    .max(0.01);
  var a = zs.expression(
    '10 ** (-1.195 + 0.4459 * axlog - 0.0345 * axlog * axlog)',
    {'axlog': uu.log10()});
  var watabs = zs.expression('1320.0 * a', {'a': a});
  var potbm2 = zs.expression(
    '720.0 * exp(-0.05 * airmas) - watabs',
    {'airmas': airmas, 'watabs': watabs});
  var evaL = zs.expression(
    '(720.0 - potbm2 - watabs) * 0.54 * cos(zs)',
    {'potbm2': potbm2, 'watabs': watabs, 'zs': zs});
  var potnir = zs.expression(
    'evaL + potbm2 * cos(zs)',
    {'evaL': evaL, 'potbm2': potbm2, 'zs': zs});

  var fclear = zs
    .expression(
      'Rs_1 / (potvis + potnir)',
      {'potvis': potvis, 'potnir': potnir, 'Rs_1': Rs_1})
    .clamp(0.01, 1.0)
    .where(zs.cos().lte(0.01), 1);

  // Partition SDN into VIS and NIR
  var fvis = zs.expression(
    'potvis / (potvis + potnir)',
      {'potvis': potvis, 'potnir': potnir});
  var fnir = zs.expression(
    'potnir / (potvis + potnir)',
      {'potvis': potvis, 'potnir': potnir});

  // Estimate direct beam and diffuse fraction in VIS and NIR wavebands
  var fb1 = zs.expression(
    'potbm1 * cos(zs) / potvis',
    {'potbm1': potbm1, 'potvis': potvis, 'zs': zs});
  var fb2 = zs.expression(
    'potbm2 * cos(zs) / potnir',
    {'potbm2': potbm2, 'potnir': potnir, 'zs': zs});

  var dirvis = zs
    .expression(
      'fb1 * (1.0 - ((0.9 - ratiox) / 0.7) ** 0.6667)',
      {'fb1': fb1, 'ratiox': fclear.min(0.9)})
    .min(fb1);
  var dirnir = zs
    .expression(
      'fb1 * (1.0 - ((0.88 - ratiox) / 0.68) ** 0.6667)',
      {'fb1': fb1, 'ratiox': fclear.min(0.88)})
    .min(fb1);

  dirvis = dirvis.where(dirvis.lt(0.01).and(dirnir.gt(0.01)), 0.011);
  dirnir = dirnir.where(dirnir.lt(0.01).and(dirvis.gt(0.01)), 0.011);

  var difvis = zs.expression('1.0 - dirvis', {'dirvis': dirvis});
  var difnir = zs.expression('1.0 - dirnir', {'dirnir': dirnir});

  // Correction for NIGHTIME
  ind = zs.cos().lte(0.01);
  fvis = fvis.where(ind, 0.5);
  fnir = fnir.where(ind, 0.5);
  difvis = difvis.where(ind, 0.0);
  difnir = difnir.where(ind, 0.0);

  // CGM - Not used anymore in function since e_atm is not computed
  // Rs0 = zs \
  //     .expression('potvis + potnir', {'potnir': potnir, 'potvis': potvis}) \
  //     .where(zs.cos().lte(0.01), 0.0)

  //**********************************************
  // Compute Albedo
  var ratio_soil = 2.0;

  // CGM - Initialize rsoilv and fg from F and albedo
  rsoilv = F.multiply(0).add(0.12);
  var fg = albedo.multiply(0).add(1);
  // rsoilv = ee.Image.constant(0.12)
  // fg = ee.Image.constant(1.0)

  // CGM - Switched to an iterate call
  function albedo_iter_func(n, prev) {
    // Extract inputs from previous iteration
    // CGM - Variables that are commented out only need to be returned
    // akb = ee.Image(ee.Dictionary(prev).get('akb'));
    // albedo_c = ee.Image(ee.Dictionary(prev).get('albedo_c'));
    // albedo_s = ee.Image(ee.Dictionary(prev).get('albedo_s'));
    // ameann = ee.Image(ee.Dictionary(prev).get('ameann'));
    // ameanv = ee.Image(ee.Dictionary(prev).get('ameanv'));
    // diff = ee.Image(ee.Dictionary(prev).get('diff'));
    var fg_iter = ee.Image(ee.Dictionary(prev).get('fg'));
    // rbcpyn = ee.Image(ee.Dictionary(prev).get('rbcpyn'));
    // rbcpyv = ee.Image(ee.Dictionary(prev).get('rbcpyv'));
    var rsoilv_iter = ee.Image(ee.Dictionary(prev).get('rsoilv'));
    // taudn = ee.Image(ee.Dictionary(prev).get('taudn'));
    // taudv = ee.Image(ee.Dictionary(prev).get('taudv'));

    var rsoiln = rsoilv_iter.multiply(ratio_soil);
    // rsoiln = .expression(
    //   'rsoilv * ratio_soil',
    //   {'rsoilv': rsoilv, 'ratio_soil': ratio_soil})

    // Weighted live/dead leaf average properties
    var ameanv = aleafv.expression(
      'aleafv * fg + adeadv * (1.0 - fg)',
      {'adeadv': adeadv, 'aleafv': aleafv, 'fg': fg_iter});
    var ameann = aleafn.expression(
      'aleafn * fg + adeadn * (1.0 - fg)',
      {'adeadn': adeadn, 'aleafn': aleafn, 'fg': fg_iter});
    var ameanl = aleafl.expression(
      'aleafl * fg + adeadl * (1.0 - fg)',
      {'adeadl': adeadl, 'aleafl': aleafl, 'fg': fg_iter});

    // DIFFUSE COMPONENT
    //*******************************
    // Canopy reflection (deep canopy)
    // Fit to Fig 15.4 for x=1
    var akd = F.expression('-0.0683 * log(F) + 0.804', {'F': F});

    // Eq 15.7
    var rcpyn = ameann.expression(
      '(1.0 - sqrt(ameann)) / (1.0 + sqrt(ameann))',
      {'ameann': ameann});
    var rcpyv = ameanv.expression(
      '(1.0 - sqrt(ameanv)) / (1.0 + sqrt(ameanv))',
      {'ameanv': ameanv});
    // rcpyl = ameanl.expression(
    //   '(1.0 - sqrt(ameanl)) / (1.0 + sqrt(ameanl))',
    //   {'ameanl': ameanl})

    // Eq 15.8
    var rdcpyn = akd.expression(
      '2.0 * akd * rcpyn / (akd + 1.0)',
      {'akd': akd, 'rcpyn': rcpyn});
    var rdcpyv = akd.expression(
      '2.0 * akd * rcpyv / (akd + 1.0)',
      {'akd': akd, 'rcpyv': rcpyv});
    // rdcpyl = akd.expression(
    //   '2.0 * akd * rcpyl / (akd + 1.0)', {'akd': akd, 'rcpyl': rcpyl})

    // Canopy transmission (VIS)
    var expfac = F.expression(
      'sqrt(ameanv) * akd * F', {'akd': akd, 'ameanv': ameanv, 'F': F});
    expfac = expfac.max(0.001);
    // expfac = expfac.where(expfac.lt(0.001), 0.001)
    var xnum = F.expression(
      '(rdcpyv * rdcpyv - 1.0) * exp(-expfac)',
      {'rdcpyv': rdcpyv, 'expfac': expfac});
    var xden = F.expression(
      '(rdcpyv * rsoilv - 1.0) + rdcpyv * (rdcpyv - rsoilv) * exp(-2.0 * expfac)',
      {'expfac': expfac, 'rdcpyv': rdcpyv, 'rsoilv': rsoilv_iter});
    // Eq 15.11
    var taudv = F.expression('xnum / xden', {'xden': xden, 'xnum': xnum});
    // taudv = xnum.divide(xden)

    // Canopy transmission (NIR)
    expfac = F.expression(
      'sqrt(ameann) * akd * F', {'akd': akd, 'ameann': ameann, 'F': F});
    expfac = expfac.max(0.001);
    // expfac = expfac.where(expfac.lt(0.001), 0.001)
    xnum = F.expression(
      '(rdcpyn * rdcpyn - 1.0) * exp(-expfac)',
      {'expfac': expfac, 'rdcpyn': rdcpyn});
    xden = F.expression(
      '(rdcpyn * rsoiln - 1.0) + rdcpyn * (rdcpyn - rsoiln) * exp(-2.0 * expfac)',
      {'expfac': expfac, 'rdcpyn': rdcpyn, 'rsoiln': rsoiln});
    // Eq 15.11
    var taudn = F.expression('xnum / xden', {'xden': xden, 'xnum': xnum});
    // taudn = xnum.divide(nden);

    // Canopy transmission (LW)
    var taudl = F.expression(
      'exp(-sqrt(ameanl) * akd * F)',
      {'akd': akd, 'ameanl': ameanl, 'F': F});

    // Diffuse albedo for generic canopy
    // Eq 15.9
    var fact = F.expression(
      '((rdcpyn - rsoiln) / (rdcpyn * rsoiln - 1.0)) * '+
      'exp(-2.0 * sqrt(ameann) * akd * F)',
      {'akd': akd, 'ameann': ameann, 'F': F, 'rdcpyn': rdcpyn,
      'rsoiln': rsoiln});
    var albdn = F.expression(
      '(rdcpyn + fact) / (1.0 + rdcpyn * fact)',
      {'fact': fact, 'rdcpyn': rdcpyn});

    // Eq 15.9
    fact = F.expression(
      '((rdcpyv - rsoilv) / (rdcpyv * rsoilv - 1.0)) * exp(-2.0 * sqrt(ameanv) * akd * F)',
      {'akd': akd, 'ameanv': ameanv, 'F': F, 'rdcpyv': rdcpyv, 'rsoilv': rsoilv_iter});
    var albdv = F.expression(
      '(rdcpyv + fact) / (1.0 + rdcpyv * fact)',
      {'fact': fact, 'rdcpyv': rdcpyv});

    // BEAM COMPONENT
    //*******************************
    // Canopy reflection (deep canopy)
    var akb = zs.expression('0.5 / cos(zs)', {'zs': zs});
    akb = akb.where(zs.cos().lte(0.01), 0.5);

    // Eq 15.7
    rcpyn = ameann.expression(
      '(1.0 - sqrt(ameann)) / (1.0 + sqrt(ameann))',
      {'ameann': ameann});
    rcpyv = ameanv.expression(
      '(1.0 - sqrt(ameanv)) / (1.0 + sqrt(ameanv))',
      {'ameanv': ameanv});

    // Eq 15.8
    var rbcpyn = rcpyn.expression(
      '2.0 * akb * rcpyn / (akb + 1.0)',
      {'akb': akb, 'rcpyn': rcpyn});
    var rbcpyv = rcpyv.expression(
      '2.0 * akb * rcpyv / (akb + 1.0)',
      {'akb': akb, 'rcpyv': rcpyv});

    // Beam albedo for generic canopy
    // Eq 15.9
    fact = F.expression(
      '((rbcpyn - rsoiln) / (rbcpyn * rsoiln - 1.0)) * exp(-2.0 * sqrt(ameann) * akb * F)',
      {'akb': akb, 'ameann': ameann, 'F': F, 'rbcpyn': rbcpyn, 'rsoiln': rsoiln});
    var albbn = F.expression(
      '(rbcpyn + fact) / (1.0 + rbcpyn * fact)',
      {'fact': fact, 'rbcpyn': rbcpyn});

    // Eq 15.9
    fact = F.expression(
      '((rbcpyv - rsoilv) / (rbcpyv * rsoilv - 1.0)) * exp(-2.0 * sqrt(ameanv) * akb * F)',
      {'akb': akb, 'ameanv': ameanv, 'F': F, 'rbcpyv': rbcpyv, 'rsoilv': rsoilv_iter});
    var albbv = F.expression(
      '(rbcpyv + fact) / (1.0 + rbcpyv * fact)',
      {'fact': fact, 'rbcpyv': rbcpyv});

    // CGM - finish
    // Weighted albedo (canopy)
    var albedo_c = F.expression(
      'fvis * (dirvis * albbv + difvis * albdv) + fnir * (dirnir * albbn + difnir * albdn)',
      {'albbn': albbn, 'albbv': albbv, 'albdn': albdn, 'albdv': albdv,
      'difnir': difnir, 'difvis': difvis, 'dirvis': dirvis,
      'dirnir': dirnir, 'fnir': fnir, 'fvis': fvis, });
    albedo_c = albedo_c.where(
      zs.cos().lte(0.01),
      F.expression(
        'fvis * (difvis * albdv) + fnir * (difnir * albdn)',
        {'albdn': albdn, 'albdv': albdv, 'difnir': difnir,
        'difvis': difvis, 'fnir': fnir, 'fvis': fvis}));

    var albedo_s = rsoilv.expression(
      'fvis * rsoilv + fnir * rsoiln',
      {'fnir': fnir, 'fvis': fvis, 'rsoiln': rsoiln, 'rsoilv': rsoilv_iter});

    var albedo_avg = fc.expression(
      '(fc * albedo_c) + ((1 - fc) * albedo_s)',
      {'albedo_c': albedo_c, 'albedo_s': albedo_s, 'fc': fc});
    var diff = albedo_avg.subtract(albedo);
    // diff = albedo_avg.expression(
    //   'albedo_avg - albedo',
    //   {'albedo_avg': albedo_avg, 'albedo': albedo})

    // CGM - Check what this is doing
    // Extra select call is needed if LAI is multiband
    // Added fc_mask call
    var fc_mask = fc.select([0]).lt(0.75);
    rsoilv_iter = rsoilv_iter
      .where(fc_mask.and(diff.lte(-0.01)), rsoilv_iter.add(0.01))
      .where(fc_mask.and(diff.gt(0.01)), rsoilv_iter.add(-0.01));
    // # CGM - IDL function
    // rsoilv = ((fc lt 0.75) * (
    //             ((abs(diff) le 0.01) * rsoilv) +
    //             ((diff le -0.01)*(rsoilv + 0.01)) +
    //             ((diff gt 0.01)*(rsoilv - 0.01))))+
    //          ((fc ge 0.75) * rsoilv)

    // CGM - Extra select call is needed since fc is multiband
    fc_mask = fc.select([0]).gte(0.75);
    fg_iter = fg_iter
      .where(fc_mask.and(diff.lte(-0.01)), fg_iter.subtract(0.05))
      .where(fc_mask.and(diff.gt(0.01)), fg_iter.add(0.05))
      .clamp(0.01, 1);
    // # CGM - IDL function
    // fg = ((fc ge 0.75) * (
    //          ((abs(diff) le 0.01)*fg) +
    //          ((diff le -0.01) * (fg - 0.05d0)) +
    //          ((diff gt 0.01) * (fg + 0.05d0)))) +
    //      ((fc lt 0.75) * fg)

    return ee.Dictionary({
      'akb': akb, 'albedo_c': albedo_c, 'albedo_s': albedo_s,
      'ameann': ameann, 'ameanv': ameanv, 'diff': diff, 'fg': fg_iter,
      'rbcpyn': rbcpyn, 'rbcpyv': rbcpyv,
      'rsoiln': rsoiln, 'rsoilv': rsoilv_iter,
      'taudn': taudn, 'taudv': taudv
    });
  }

  // Iterate the function n times
  var albedo_input_images = ee.Dictionary({
    'akb': null, 'albedo_c': null, 'albedo_s': null,
    'ameann': null, 'ameanv': null, 'diff': null, 'fg': fg,
    'rbcpyn': null, 'rbcpyv': null,
    'rsoiln': null, 'rsoilv': rsoilv,
    'taudn': null, 'taudv': null
  });

  var albedo_iter_output = ee.Dictionary(
    ee.List.repeat(albedo_input_images, iterations)
      .iterate(albedo_iter_func, albedo_input_images));
    // ee.List.sequence(1, iterations)

  // Unpack the iteration output
  var akb = ee.Image(albedo_iter_output.get('akb'));
  var albedo_c = ee.Image(albedo_iter_output.get('albedo_c'));
  var albedo_s = ee.Image(albedo_iter_output.get('albedo_s'));
  var ameann = ee.Image(albedo_iter_output.get('ameann'));
  var ameanv = ee.Image(albedo_iter_output.get('ameanv'));
  var diff = ee.Image(albedo_iter_output.get('diff'));
  var rbcpyn = ee.Image(albedo_iter_output.get('rbcpyn'));
  var rbcpyv = ee.Image(albedo_iter_output.get('rbcpyv'));
  var rsoilv = ee.Image(albedo_iter_output.get('rsoilv'));
  var rsoiln = ee.Image(albedo_iter_output.get('rsoiln'));
  // var rsoiln = rsoilv.multiply(ratio_soil)
  var taudn = ee.Image(albedo_iter_output.get('taudn'));
  var taudv = ee.Image(albedo_iter_output.get('taudv'));

  // if a solution is not reached, alb_c=alb_s=alb
  albedo_c = albedo_c.where(diff.abs().gt(0.05), albedo);
  albedo_s = albedo_s.where(diff.abs().gt(0.05), albedo);

  // Direct beam+scattered canopy transmission coefficient (visible)
  var expfac = F.expression(
    'sqrt(ameanv) * akb * F',
    {'ameanv': ameanv, 'akb': akb, 'F': F});
  var xnum = F.expression(
    '(rbcpyv * rbcpyv - 1.0) * exp(-expfac)',
    {'rbcpyv': rbcpyv, 'expfac': expfac});
  var xden = F.expression(
    '(rbcpyv * rsoilv - 1.0) + rbcpyv * (rbcpyv - rsoilv) * exp(-2.0 * expfac)',
    {'rbcpyv': rbcpyv, 'rsoilv': rsoilv, 'expfac': expfac});
  // Eq 15.11
  var taubtv = F.expression('xnum / xden', {'xnum': xnum, 'xden': xden});

  // Direct beam+scattered canopy transmission coefficient (NIR)
  expfac = F.expression(
    'sqrt(ameann) * akb * F',
    {'ameann': ameann, 'akb': akb, 'F': F});
  xnum = F.expression(
    '(rbcpyn * rbcpyn - 1.0) * exp(-expfac)',
    {'rbcpyn': rbcpyn, 'expfac': expfac});
  xden = F.expression(
    '(rbcpyn * rsoiln - 1.0) + rbcpyn * (rbcpyn - rsoiln) * exp(-2.0 * expfac)',
    {'rbcpyn': rbcpyn, 'rsoiln': rsoiln, 'expfac': expfac});
  // Eq 15.11
  var taubtn = F.expression('xnum / xden', {'xnum': xnum, 'xden': xden});

  // Shortwave radiation components
  var tausolar = F.expression(
    'fvis * (difvis * taudv + dirvis * taubtv) + fnir * (difnir * taudn + dirnir * taubtn)',
    {'difnir': difnir, 'difvis': difvis, 'dirnir': dirnir, 'dirvis': dirvis,
    'fnir': fnir, 'fvis': fvis,
    'taubtn': taubtn, 'taubtv': taubtv, 'taudn': taudn, 'taudv': taudv});
  var Rs_c = Rs_1.expression(
    'Rs_1 * (1.0 - tausolar)',
    {'Rs_1': Rs_1, 'tausolar': tausolar});
  var Rs_s = Rs_1.expression(
    'Rs_1 * tausolar',
    {'Rs_1': Rs_1, 'tausolar': tausolar});

  // CGM - Moved emissivity calculation to separate function.
  //   I removed the Rs0 check.
  var e_atm = emissivity(T_air);
  // p = T_air.expression(
  //   '101.3 * (((T_air - (0.0065 * z)) / T_air) ** 5.26)',
  //   {'T_air': T_air, 'z': z})
  // Density of air? (kg m-3)
  var r_air = T_air.expression(
    '101.3 * (((T_air - (0.0065 * z)) / T_air) ** 5.26) / 1.01 / T_air / 0.287',
    {'T_air': T_air, 'z': z});
  var cp = ee.Number(1004.16);
  // cp = ee.Image.constant(1004.16)

  // Assume neutral conditions on first iteration (use T_air for Ts and Tc)
  // CGM - Using lai for F to match Python code
  var u_attr = compute_u_attr(u, d_0, z0m, z_u, 0);
  var r_ah = compute_r_ah(u_attr, d_0, z0h, z_t, 0);
  // CGM - Why is this function is passing "lai" to "F"?
  var r_s = compute_r_s(u_attr, T_air, T_air, hc, lai, d_0, z0m, leaf, leaf_s, 0);
  var r_x = compute_r_x(u_attr, hc, lai, d_0, z0m, leaf_width, leaf_c, 0);

  var T_c = T_air;
  // DEADBEEF - In IDL, this calculation is in C, not K?
  var T_s = lai.expression(
    '((T_rad - 273.16) - (fc_q * (T_c - 273.16))) / (1 - fc_q) + 273.16',
    {'T_rad': T_rad, 'T_c': T_c, 'fc_q': fc_q});
  // T_s = lai.expression(
  //   '(T_rad - (fc_q * T_c)) / (1 - fc_q)',
  //   {'T_rad': T_rad, 'T_c': T_c, 'fc_q': fc_q})

  // CGM - Initialize to match T_air shape
  // This doesn't seem to do anything, commenting out for now
  // H_iter = T_air.multiply(0).add(200.16)
  var EF_s = T_air.multiply(0);

  // ************************************************************************
  // Start Loop for Stability Correction and Water Stress
  function iter_func(n, prev) {
    // Extract inputs from previous iteration
    var a_PT_iter = ee.Image(ee.Dictionary(prev).get('a_PT'));
    var EF_s_iter = ee.Image(ee.Dictionary(prev).get('EF_s'));
    var r_ah_iter = ee.Image(ee.Dictionary(prev).get('r_ah'));
    var r_s_iter = ee.Image(ee.Dictionary(prev).get('r_s'));
    var r_x_iter = ee.Image(ee.Dictionary(prev).get('r_x'));
    var T_c_iter = ee.Image(ee.Dictionary(prev).get('T_c'));
    var T_s_iter = ee.Image(ee.Dictionary(prev).get('T_s'));
    var u_attr_iter = ee.Image(ee.Dictionary(prev).get('u_attr'));

    var Rn_c = compute_Rn_c(albedo_c, T_air, T_c_iter, T_s_iter, e_atm, Rs_c, F);
    var Rn_s = compute_Rn_s(albedo_s, T_air, T_c_iter, T_s_iter, e_atm, Rs_s, F);
    var Rn = Rn_c.add(Rn_s);

    var G = compute_G0(Rn, Rn_s, albedo, ndvi, t_rise, t_end, time, EF_s_iter);

    var LE_c = albedo
      .expression(
        'f_green * (a_PT * Ss / (Ss + g)) * Rn_c',
        {'f_green': f_green, 'a_PT': a_PT_iter, 'Ss': Ss, 'g': g, 'Rn_c': Rn_c})
      .max(0);
    var H_c = albedo.expression('Rn_c - LE_c', {'Rn_c': Rn_c, 'LE_c': LE_c});

    T_c_iter = temp_separation_tc(
      H_c, fc_q, T_air, T_rad, r_ah_iter, r_s_iter, r_x_iter, r_air, cp);
    T_s_iter = temp_separation_ts(T_c_iter, fc_q, T_air, T_rad);
    var T_ac = temp_separation_tac(
      T_c_iter, T_s_iter, fc_q, T_air, r_ah_iter, r_s_iter, r_x_iter);

    H_s = albedo.expression(
      'r_air * cp * (T_s - T_ac) / r_s',
      {'r_air': r_air, 'cp': cp, 'T_s': T_s_iter, 'T_ac': T_ac, 'r_s': r_s_iter});
    H_c = albedo.expression(
      'r_air * cp * (T_c - T_ac) / r_x',
      {'r_air': r_air, 'cp': cp, 'T_c': T_c_iter, 'T_ac': T_ac, 'r_x': r_x_iter});
    var H = albedo.expression('H_s + H_c', {'H_s': H_s, 'H_c': H_c});

    LE_s = albedo.expression(
      'Rn_s - G - H_s',
      {'Rn_s': Rn_s, 'G': G, 'H_s': H_s});
    LE_c = albedo.expression('Rn_c - H_c', {'Rn_c': Rn_c, 'H_c': H_c});

    // CGM - Is there a reason this isn't up with the H calculation?
    H = H.where(H.eq(0), 10.0);

    // CGM - This wont doing anything at this position in the code.
    //   Commenting out for now.
    // r_ah_iter = r_ah_iter.where(r_ah_iter.eq(0), 10.0)

    // CGM - This doesn't seem to do anything, commenting out for now
    // mask_iter = H_iter.divide(H).lte(1.05).and(H_iter.divide(H).gte(0.95))
    // chk_iter = np.sum(mask_iter) / np.size(mask_iter)

    var fh = compute_stability_fh(H, T_rad, u_attr_iter, r_air, z_t, d_0, cp);
    var fm = compute_stability_fm(H, T_rad, u_attr_iter, r_air, z_u, d_0, z0m, cp);
    var fm_h = compute_stability_fm_h(H, T_rad, u_attr_iter, r_air, hc, d_0, z0m, cp);

    u_attr_iter = compute_u_attr(u, d_0, z0m, z_u, fm);
    r_ah_iter = compute_r_ah(u_attr_iter, d_0, z0h, z_t, fh);
    r_s_iter = compute_r_s(u_attr_iter, T_s_iter, T_c_iter, hc, lai, d_0, z0m, leaf, leaf_s, fm_h);
    // CGM - Why is this function is passing "lai" to "F"?
    r_x_iter = compute_r_x(u_attr_iter, hc, lai, d_0, z0m, leaf_width, leaf_c, fm_h);

    a_PT_iter = a_PT_iter
      .where(LE_s.lte(0), a_PT_iter.subtract(0.05))
      .where(a_PT_iter.lte(0), 0.01);

    var den_s = albedo.expression('Rn_s - G', {'Rn_s': Rn_s, 'G': G});
    den_s = den_s.updateMask(den_s.neq(0));
    // den_s[den_s == 0.] = np.nan

    EF_s_iter = albedo.expression('LE_s / den_s', {'LE_s': LE_s, 'den_s': den_s});

    return ee.Dictionary({
      'a_PT': a_PT_iter, 'EF_s': EF_s_iter, 'G': G,
      'H_c': H_c, 'H_s': H_s, 'LE_c': LE_c, 'LE_s': LE_s,
      'Rn_c': Rn_c, 'Rn_s': Rn_s,
      'r_ah': r_ah_iter, 'r_s': r_s_iter, 'r_x': r_x_iter,
      'T_ac': T_ac, 'T_c': T_c_iter, 'T_s': T_s_iter,
      'u_attr': u_attr_iter
    });
  }

  // Iterate the function n times
  // CGM - Iteration count is an input to the function
  var input_images = ee.Dictionary({
    'a_PT': a_PT, 'EF_s': EF_s, 'G': ee.Image(0),
    'H_c': ee.Image(0), 'H_s': ee.Image(0),
    'LE_c': ee.Image(0), 'LE_s': ee.Image(0),
    'Rn_c': ee.Image(0), 'Rn_s': ee.Image(0),
    'r_ah': r_ah, 'r_s': r_s, 'r_x': r_x,
    'T_ac': ee.Image(0), 'T_c': T_c, 'T_s': T_s, 'u_attr': u_attr
  });
  var iter_output = ee.Dictionary(
    ee.List.sequence(1, iterations).iterate(iter_func, input_images));

  // Unpack the iteration output
  var a_PT = ee.Image(iter_output.get('a_PT'));
  var Rn_c = ee.Image(iter_output.get('Rn_c'));
  var Rn_s = ee.Image(iter_output.get('Rn_s'));
  var G = ee.Image(iter_output.get('G'));
  var H_c = ee.Image(iter_output.get('H_c'));
  var H_s = ee.Image(iter_output.get('H_s'));
  var LE_c = ee.Image(iter_output.get('LE_c'));
  var LE_s = ee.Image(iter_output.get('LE_s'));
  // T_ac = ee.Image(iter_output.get('T_ac'));
  // T_c = ee.Image(iter_output.get('T_c'));
  // T_s = ee.Image(iter_output.get('T_s'));
  // r_ah = ee.Image(iter_output.get('r_ah'));
  // r_s = ee.Image(iter_output.get('r_s'));
  // r_x = ee.Image(iter_output.get('r_x'));

  // ************************************************************************
  // Check Energy Balance Closure
  var ind = a_PT.lte(0.01);
  LE_s = LE_s.where(ind, 1.0);
  LE_c = LE_c.where(ind, 1.0);
  G = G.where(ind, Rn_s.subtract(H_s));

  ind = LE_s.gt(Rn_s);
  LE_s = LE_s.where(ind,  Rn_s);
  H_s = H_s.where(ind,  Rn_s.subtract(G).subtract(LE_s));

  // CGM - Check order of operations
  ind = LE_c.gt(Rn_c.add(100));
  // CGM - Not used below since LE_c is recomputed
  LE_c = LE_c.where(ind, Rn_c.add(100));
  H_c = H_c.where(ind, -100);

  LE_s = albedo.expression('Rn_s - G - H_s', {'Rn_s': Rn_s, 'G': G, 'H_s': H_s});
  LE_c = albedo.expression('Rn_c - H_c', {'Rn_c': Rn_c, 'H_c': H_c});

  // The latent heat of vaporization is 2.45 MJ kg-1
  // Assume Rs24 is still in W m-2 day-1 and convert to MJ kg-1
  // CGM - Leaving out scaling value for now
  var ET = albedo
    .expression(
      '((LE_c + LE_s) / Rs_1) * (Rs24 / 2.45) * scaling',
      {'LE_c': LE_c, 'LE_s': LE_s, 'Rs_1': Rs_1,
      'Rs24': Rs24.multiply(0.0864 / 24.0), 'scaling': 1})
    .max(0.01);

  return ET;
}


function emissivity(T_air) {
  // Apparent atmospheric emissivity
  var e_atm = T_air
    .expression(
      '1.0 - (0.2811 * (exp(-0.0003523 * ((T_air - 273.16) ** 2.0))))',
      {'T_air': T_air});
    // .where(Rs0.lte(50.0), 1.0);
  return e_atm;
}

function compute_G0(Rn, Rn_s, albedo, ndvi, t_rise, t_end, time, EF_s) {
  var w = EF_s.expression('1 / (1 + (EF_s / 0.5) ** 8.0)', {'EF_s': EF_s});

  // Maximum fraction of Rn,s that become G0
  // (0.35 for dry soil and 0.31 for wet soil)
  var c_g = w.expression('(w * 0.35) + ((1 - w) * 0.31)', {'w': w});
  var t_g = w.expression('(w * 100000.0) + ((1 - w) * 74000.0)', {'w': w});

  var t_noon = t_rise.expression(
    '0.5 * (t_rise + t_end)', {'t_rise': t_rise, 't_end': t_end});
  var t_g0 = t_noon.expression(
    '(time - t_noon) * 3600.0', {'time': time, 't_noon': t_noon});

  var G0 = Rn_s.expression(
    'c_g * cos(2 * pi * (t_g0 + 10800.0) / t_g) * Rn_s',
    {'c_g': c_g, 'pi': pi, 'Rn_s': Rn_s, 't_g': t_g, 't_g0': t_g0});

  var water_mask = ndvi.lte(0).and(albedo.lte(0.05));
  G0 = G0.where(water_mask, Rn.multiply(0.5));
  return G0;
}

function compute_u_attr(u, d0, z0m, z_u, fm) {
  // Friction Velocity
  var u_attr = u.expression(
    '0.41 * u / ((log((z_u - d0) / z0m)) - fm)',
    {'d0': d0, 'fm': fm, 'u': u, 'z0m': z0m, 'z_u': z_u});
  u_attr = u_attr.where(u_attr.eq(0), 10);
  u_attr = u_attr.where(u_attr.lte(0), 0.01);
  return u_attr;
}

function compute_r_ah(u_attr, d0, z0h, z_t, fh) {
  var r_ah = u_attr.expression(
    '((log((z_t - d0) / z0h)) - fh) / u_attr / 0.41',
    {'d0': d0, 'fh': fh, 'u_attr': u_attr, 'z0h': z0h, 'z_t': z_t});
  // CGM - The second conditional will overwrite the first one?
  r_ah = r_ah.where(r_ah.eq(0), 500);
  r_ah = r_ah.where(r_ah.lte(1.0), 1.0);
  return r_ah;
}

function compute_r_s(u_attr, T_s, T_c, hc, F, d0, z0m, leaf, leaf_s, fm_h) {
  // Free convective velocity constant for r_s modelling
  var c_a = 0.004;
  // Empirical constant for r_s modelling
  var c_b = 0.012;
  // Empirical constant for r_s modelling
  // (new formulation Kustas and Norman, 1999)
  var c_c = 0.0025;

  // Computation of the resistance of the air between soil and canopy space
  var u_c = u_attr.expression(
    'u_attr / 0.41 * ((log((hc - d0) / z0m)) - fm_h)',
    {'d0': d0, 'fm_h': fm_h, 'hc': hc, 'u_attr': u_attr, 'z0m': z0m});
  u_c = u_c.where(u_c.lte(0), 0.1);
  var u_s = u_attr.expression(
    'u_c * exp(-leaf * (1 - (0.05 / hc)))',
    {'hc': hc, 'leaf': leaf, 'u_c': u_c});

  var r_ss = u_attr.expression(
    '1.0 / (c_a + (c_b * (u_c * exp(-leaf_s * (1.0 - (0.05 / hc))))))',
    {'c_a': c_a, 'c_b': c_b, 'hc': hc, 'leaf_s': leaf_s, 'u_c': u_c});
  var r_s1 = T_s.expression(
    '1.0 / ((((abs(T_s - T_c)) ** (1.0 / 3.0)) * c_c) + (c_b * Us))',
    {'c_b': c_b, 'c_c': c_c, 'T_c': T_c, 'T_s': T_s, 'Us': u_s});
  var r_s2 = u_attr.expression(
    '1.0 / (c_a + (c_b * Us))', {'c_a': c_a, 'c_b': c_b, 'Us': u_s});
  var r_s = u_attr.expression(
    '(((r_ss - 1.0) / 0.09 * (F - 0.01)) + 1.0)', {'F': F, 'r_ss': r_ss});

  // Linear function between 0 (bare soil) and the value at F=0.1
  r_s = r_s.where(F.gt(0.1), r_s1);
  r_s = r_s.where(T_s.subtract(T_c).abs().lt(1), r_s2);

  // Use "new" formula only for high DT values
  // Use "new" formula only for partial coverage (lai<3)
  r_s = r_s.where(F.gt(3), r_s2);
  return r_s;
}

function compute_r_x(u_attr, hc, F, d0, z0m, xl, leaf_c, fm_h) {
  // Parameter for canopy boundary-layer resistance
  // (C=90 Grace '81, C=175 Cheubouni 2001, 144 Li '98)
  var C = 175.0;

  // Computation of the resistance of the air between soil and canopy space
  var u_c = u_attr.expression(
    'u_attr / 0.41 * ((log((hc - d0) / z0m)) - fm_h)',
    {'d0': d0, 'fm_h': fm_h, 'hc': hc, 'u_attr': u_attr, 'z0m': z0m});
  u_c = u_c.where(u_c.lte(0), 0.1);

  // Computation of the canopy boundary layer resistance
  var u_d = u_attr.expression(
    'u_c * exp(-leaf_c * (1 - ((d0 + z0m) / hc)))',
    {'d0': d0, 'hc': hc, 'leaf_c': leaf_c, 'u_c': u_c, 'z0m': z0m});
  u_d = u_d.where(u_d.lte(0), 100);

  var r_x = u_attr.expression(
    'C / F * ((xl / u_d) ** 0.5)', {'C': C, 'F': F, 'u_d': u_d, 'xl': xl});
  r_x = r_x.where(u_d.eq(100), 0.1);
  return r_x;
}

function compute_Rn_c(albedo_c, T_air, T_c, T_s, e_atm, Rs_c, F) {
  // Compute Canopy Net Radiation

  // Long-wave extinction coefficient [-]
  var kL = 0.95;
  // Soil Emissivity [-]
  var eps_s = 0.94;
  // Canopy emissivity [-]
  var eps_c = 0.99;

  // Stephan Boltzmann constant (W m-2 K-4)
  // sb = 5.670373e-8
  var Lc = T_c.expression(
    'eps_c * 5.67E-8 * (T_c ** 4)', {'eps_c': eps_c, 'T_c': T_c});
  var Ls = T_s.expression(
    'eps_s * 5.67E-8 * (T_s ** 4)', {'eps_s': eps_s, 'T_s': T_s});
  var Rle = T_air.expression(
    'e_atm * 5.67E-8 * (T_air ** 4)',
    {'e_atm': e_atm, 'T_air': T_air});
  var Rn_c = albedo_c.expression(
    '((1 - albedo_c) * Rs_c) + ' +
    '((1 - exp(-kL * F)) * (Rle + Ls - 2 * Lc))',
    {'albedo_c': albedo_c, 'F': F, 'kL': kL, 'Lc': Lc, 'Ls': Ls,
     'Rle': Rle, 'Rs_c': Rs_c});
  return Rn_c;
}

function compute_Rn_s(albedo_s, T_air, T_c, T_s, e_atm, Rs_s, F) {
  // Compute Soil Net Radiation

  // Long-wave extinction coefficient [-]
  var kL = 0.95;
  // Soil Emissivity [-]
  var eps_s = 0.94;
  // Canopy emissivity [-]
  var eps_c = 0.99;

  var L_c = T_c.expression(
    'eps_c * 0.0000000567 * (T_c ** 4)', {'eps_c': eps_c, 'T_c': T_c});
  var L_s = T_s.expression(
    'eps_s * 0.0000000567 * (T_s ** 4)', {'eps_s': eps_s, 'T_s': T_s});
  var Rle = T_air.expression(
    'e_atm * 0.0000000567 * (T_air ** 4)',
    {'e_atm': e_atm, 'T_air': T_air});
  var Rn_s = albedo_s.expression(
    '((1 - albedo_s) * Rs_s) + ' +
    '((exp(-kL * F)) * Rle) + ((1 - exp(-kL * F)) * L_c) - L_s',
    {'albedo_s': albedo_s, 'F': F, 'kL': kL, 'L_c': L_c, 'L_s': L_s,
     'Rle': Rle, 'Rs_s': Rs_s});
  return Rn_s;
}

function temp_separation_tc(H_c, fc, T_air, t0, r_ah, r_s, r_x, r_air, cp) {
  // Compute canopy temperature
  var T_c_lin = fc.expression(
    '((T_air / r_ah) + ' +
    ' (t0 / r_s / (1 - fc)) + ' +
    ' (H_c * r_x / r_air / cp * ((1 / r_ah) + (1 / r_s) + (1 / r_x)))) / ' +
    '((1 / r_ah) + (1 / r_s) + (fc / r_s / (1 - fc)))',
    {'cp': cp, 'fc': fc, 'H_c': H_c, 'r_ah': r_ah, 'r_air': r_air,
     'r_s': r_s, 'r_x': r_x, 't0': t0, 'T_air': T_air});
  var Td = fc.expression(
    '(T_c_lin * (1 + (r_s / r_ah))) - ' +
    '(H_c * r_x / r_air / cp * (1 + (r_s / r_x) + (r_s / r_ah))) - ' +
    '(T_air * r_s / r_ah)',
    {'cp': cp, 'H_c': H_c, 'r_ah': r_ah, 'r_air': r_air, 'r_s': r_s,
     'r_x': r_x, 'T_air': T_air, 'T_c_lin': T_c_lin});
  var delta_T_c = fc.expression(
    '((t0 ** 4) - (fc * (T_c_lin ** 4)) - ((1 - fc) * (Td ** 4))) / ' +
    '((4 * (1 - fc) * (Td ** 3) * (1 + (r_s / r_ah))) + (4 * fc * (T_c_lin ** 3)))',
    {'fc': fc, 'r_ah': r_ah, 'r_s': r_s, 't0': t0, 'Td': Td,
     'T_c_lin': T_c_lin});
  var T_c = fc
    .expression('T_c_lin + delta_T_c', {'T_c_lin': T_c_lin, 'delta_T_c': delta_T_c})
    .where(fc.lt(0.10), t0)
    .where(fc.gt(0.90), t0);
  T_c = T_c.where(T_c.lte(T_air.subtract(10.0)), T_air.subtract(10.0));
  T_c = T_c.where(T_c.gte(T_air.add(50.0)), T_air.add(50.0));
  return T_c;
}

function temp_separation_ts(T_c, fc, T_air, t0) {
  // Compute soil temperature
  var Delta = fc.expression(
    '(t0 ** 4) - (fc * (T_c ** 4))', {'fc': fc, 't0': t0, 'T_c': T_c});
  Delta = Delta.where(Delta.lte(0), 10);

  // CGM - This could probably be simplified
  var T_s = fc
    .expression('(Delta / (1 - fc)) ** 0.25', {'Delta': Delta, 'fc': fc})
    .where(
      fc.expression(
        '((t0 ** 4) - (fc * T_c ** 4)) <= 0.',
        {'fc': fc, 't0': t0, 'T_c': T_c}),
      fc.expression(
        '(t0 - (fc * T_c)) / (1 - fc)',
        {'fc': fc, 't0': t0, 'T_c': T_c}))
    .where(fc.lt(0.1), t0)
    .where(fc.gt(0.9), t0);
  T_s = T_s.where(T_s.lte(T_air.subtract(10.0)), T_air.subtract(10.0));
  T_s = T_s.where(T_s.gte(T_air.add(50.0)), T_air.add(50.0));
  return T_s;
}

function temp_separation_tac(T_c, T_s, fc, T_air, r_ah, r_s, r_x) {
  // Compute air temperature at the canopy interface
  var T_ac = fc.expression(
    '((T_air / r_ah) + (T_s / r_s) + (T_c / r_x)) / '+
    '((1 / r_ah) + (1 / r_s) + (1 / r_x))',
    {'r_ah': r_ah, 'r_s': r_s, 'r_x': r_x, 'T_c': T_c, 'T_s': T_s,
     'T_air': T_air});
  return T_ac;
}

function compute_stability_fh(H, t0, u_attr, r_air, z_t, d0, cp) {
  var L_ob = H.expression(
    '-(r_air * cp * t0 * (u_attr ** 3.0) / 0.41 / 9.806 / H)',
    {'cp': cp, 'H': H, 'r_air': r_air, 't0': t0, 'u_attr': u_attr});
  L_ob = L_ob.where(L_ob.gte(0), -99);
  var mh = H
    .expression(
      '((1 - (16.0 * (z_t - d0) / L_ob)) ** 0.25)',
      {'d0': d0, 'L_ob': L_ob, 'z_t': z_t})
    .where(L_ob.eq(-99), 0.0);
  var fh = H
    .expression('(2.0 * log((1.0 + (mh ** 2.0)) / 2.0))', {'mh': mh})
    .where(L_ob.lte(-100).or(L_ob.gte(100)), 0);
  return fh;
}

function compute_stability_fm(H, t0, u_attr, r_air, z_u, d0, z0m, cp) {
  var L_ob = H.expression(
    '-(r_air * cp * t0 * (u_attr ** 3.0) / 0.41 / 9.806 / H)',
    {'cp': cp, 'H': H, 'r_air': r_air, 't0': t0, 'u_attr': u_attr});
  L_ob = L_ob.where(L_ob.gte(0), -99.0);
  var mh = H
    .expression(
      '((1 - (16.0 * (z_u - d0) / L_ob)) ** 0.25)',
      {'d0': d0, 'L_ob': L_ob, 'z_u': z_u})
    .where(L_ob.eq(-99.0), 0.0);
  var fm = H
    .expression(
      '2.0 * log((1.0 + mh) / 2.0) + log((1.0 + (mh ** 2)) / 2.0) - ' +
      '2.0 * atan(mh) + (pi / 2)',
      {'mh': mh, 'pi': pi})
    .where(L_ob.lte(-100).or(L_ob.gte(100)), 0);

  // CGM - Swapped order of calc since d0 is an image compute from hc and
  //   z_u is being set as a constant number (for now).
  fm = fm.where(fm.eq(d0.multiply(-1).add(z_u).divide(z0m).log()), fm.add(1.0));
  // fm = fm.where(fm.eq(z_u.subtract(d0).divide(z0m).log()), fm.add(1.0))
  return fm;
}

function compute_stability_fm_h(H, t0, u_attr, r_air, hc, d0, z0m, cp) {
  var L_ob = H.expression(
    '-(r_air * cp * t0 * (u_attr ** 3.0) / 0.41 / 9.806 / H)',
    {'cp': cp, 'H': H, 'r_air': r_air, 't0': t0, 'u_attr': u_attr});
  L_ob = L_ob.where(L_ob.gte(0), -99.0);
  var mm_h = H
    .expression(
      '((1 - (16.0 * (hc - d0) / L_ob)) ** 0.25)',
      {'d0': d0, 'hc': hc, 'L_ob': L_ob})
    .where(L_ob.eq(-99.0), 0.0);
  var fm_h = H
    .expression(
      '2.0 * log((1.0 + mm_h) / 2.0) + log((1.0 + (mm_h ** 2)) / 2.0) - ' +
      '2.0 * atan(mm_h) + (pi / 2)',
      {'mm_h': mm_h, 'pi': pi})
    .where(L_ob.lte(-100).or(L_ob.gte(100)), 0);

  // CGM - Swapped order of calc since d0 is an image compute from hc and
  //   z_u is being set as a constant number (for now).
  fm_h = fm_h.where(fm_h.eq(hc.subtract(d0).divide(z0m).log()), fm_h.add(1.0));
  // fm_h = fm_h.where(fm_h.eq(hc.subtract(d0).divide(z0m).log()), fm_h.add(1.0));
  return fm_h;
}

// Landsat functions
function landsat_init(raw_img) {
  // "Prep" the Landsat image for DisALEXI
  // Rename bands to generic names
  var input_bands = ee.Dictionary({
    'LANDSAT_5': ['B1', 'B2', 'B3', 'B4', 'B5', 'B7', 'B6', 'BQA'],
    'LANDSAT_7': ['B1', 'B2', 'B3', 'B4', 'B5', 'B7', 'B6_VCID_1', 'BQA'],
    'LANDSAT_8': ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B10', 'BQA']});
  var output_bands = ['blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'lst', 'bqa'];
  // Rename thermal band "k" coefficients to generic names
  var input_k1 = ee.Dictionary({
    'LANDSAT_5': 'K1_CONSTANT_BAND_6',
    'LANDSAT_7': 'K1_CONSTANT_BAND_6_VCID_1',
    'LANDSAT_8': 'K1_CONSTANT_BAND_10'});
  var input_k2 = ee.Dictionary({
    'LANDSAT_5': 'K2_CONSTANT_BAND_6',
    'LANDSAT_7': 'K2_CONSTANT_BAND_6_VCID_1',
    'LANDSAT_8': 'K2_CONSTANT_BAND_10'});
  var spacecraft_id = ee.String(raw_img.get('SPACECRAFT_ID'));
  var landsat_img = ee.Image(raw_img)
    .select(input_bands.get(spacecraft_id), output_bands)
    .set('k1_constant', ee.Number(raw_img.get(input_k1.get(spacecraft_id))))
    .set('k2_constant', ee.Number(raw_img.get(input_k2.get(spacecraft_id))));
  return landsat_img;
}

function landsat_albedo(landsat_img) {
  // Compute total shortwave broadband albedo following [Liang2001]
  return ee.Image(landsat_img)
    .select(['blue', 'red', 'nir', 'swir1', 'swir2'])
    .multiply([0.356, 0.130, 0.373, 0.085, 0.072])
    .reduce(ee.Reducer.sum())
    .subtract(0.0018)
    .rename(['albedo']);
}

function landsat_lai(landsat_img) {
  // Compute LAI using METRIC NDVI / LAI empirical equation
  var ndvi = landsat_img.normalizedDifference(['nir', 'red']).rename(['lai']);
  return ndvi.pow(3).multiply(7.0).clamp(0, 6);
}

function landsat_emissivity(landsat_img) {
  // METRIC narrowband emissivity
  var lai = landsat_lai(landsat_img);
  var ndvi = landsat_ndvi(landsat_img);
  // Initial values are for NDVI > 0 and LAI <= 3
  return lai.divide(300).add(0.97)
    .where(ndvi.lte(0), 0.99)
    .where(ndvi.gt(0).and(lai.gt(3)), 0.98);
}

function landsat_lst(landsat_img) {
  // Compute emissivity corrected land surface temperature (LST) from brightness temperature.
  // Get properties from image
  var k1 = ee.Number(landsat_img.get('k1_constant'));
  var k2 = ee.Number(landsat_img.get('k2_constant'));

  var ts_brightness = landsat_img.select(['lst']);
  var emissivity = landsat_emissivity(landsat_img);

  // First back out radiance from brightness temperature
  // Then recalculate emissivity corrected Ts
  var thermal_rad_toa = ts_brightness.expression(
    'k1 / (exp(k2 / ts_brightness) - 1)',
    {'ts_brightness': ts_brightness, 'k1': k1, 'k2': k2});

  // tnb = 0.866   # narrow band transmissivity of air
  // rp = 0.91     # path radiance
  // rsky = 1.32   # narrow band clear sky downward thermal radiation
  var rc = thermal_rad_toa.expression(
    '((thermal_rad_toa - rp) / tnb) - ((1. - emiss) * rsky)',
    {'thermal_rad_toa': thermal_rad_toa, 'emiss': emissivity,
     'rp': 0.91, 'tnb': 0.866, 'rsky': 1.32});
  var lst = rc.expression(
    'k2 / log(emiss * k1 / rc + 1)',
    {'emiss': emissivity, 'rc': rc, 'k1': k1, 'k2': k2});
  return lst.rename(['lst']);
}

function landsat_ndvi(landsat_img) {
  return landsat_img.normalizedDifference(['nir', 'red']).rename(['ndvi']);
}

function landsat_cloud_mask(landsat_img) {
  // Extract CFmask from Landsat Collection 1 BQA band
  var bqa_image = landsat_img.select(['bqa']);

  // Extract the various masks from the QA band
  var fill_mask = bqa_image.rightShift(0).bitwiseAnd(1);
  // Landsat 8 - drop_mask is terrain_mask
  // drop_mask = bqa_image.rightShift(1).bitwiseAnd(1);
  // saturation_mask = bqa_image.rightShift(2).bitwiseAnd(3).gte(2);
  // cloud_mask = bqa_image.rightShift(4).bitwiseAnd(1);
  var cloud_mask = bqa_image.rightShift(5).bitwiseAnd(3).gte(2);
  var shadow_mask = bqa_image.rightShift(7).bitwiseAnd(3).gte(3);
  var snow_mask = bqa_image.rightShift(9).bitwiseAnd(3).gte(3);
  // cirrus_mask = bqa_image.rightShift(11).bitwiseAnd(3).gte(3)  // Landsat 8 only

  // Convert masks to old style Fmask values
  // 0 - Clear land
  // 1 - Clear water
  // 2 - Cloud shadow
  // 3 - Snow
  // 4 - Cloud
  return fill_mask
    .add(shadow_mask.multiply(2))
    .add(snow_mask.multiply(3))
    .add(cloud_mask.multiply(4))
    .rename(['cfmask']);
}

function range(start, stop, step) {
  // JS equivalent to Python range function
  if (typeof stop == 'undefined') {
    // one param defined
    stop = start;
    start = 0;
  }
  if (typeof step == 'undefined') {
    step = 1;
  }
  if ((step > 0 && start >= stop) || (step < 0 && start <= stop)) {
    return [];
  }
  var result = [];
  for (var i = start; step > 0 ? i < stop : i > stop; i += step) {
    result.push(i);
  }
  return result;
}


// Start with a single Landsat image
var image_id = 'LC08_044033_20150711';
var raw_img = ee.Image('LANDSAT/LC08/C01/T1_RT_TOA/' + image_id);
Map.addLayer(raw_img.select([3,2,1]), {min:0, max:0.3}, 'True Color');
Map.centerObject(raw_img);


// Clipped ALEXI grid and spatial reference
var output_crs = 'EPSG:4326';
var output_geo = [0.04,0.0,-123.04,0.0,-0.04,39.98];
var output_cs_str = '0p04';
var output_shape_str = '69x54';

// Landsat grid and spatial reference
var landsat_crs = raw_img.select(['B2']).projection().getInfo()['crs'];
var landsat_geo = raw_img.select(['B2']).projection().getInfo()['transform'];
var landsat_shape = raw_img.select(['B2']).getInfo()['bands'][0]['dimensions'];
var landsat_shape_str = landsat_shape[0] + 'x' + landsat_shape[1];



// Prep the Landsat image for DisALEXI
// Rename Landsat bands to common band names
var landsat_img = landsat_init(raw_img);
//Map.addLayer(landsat_img, {}, 'Prepped Image');


// Build the DisALEXI input image
var input_img = ee.Image([
  landsat_albedo(landsat_img),
  landsat_cloud_mask(landsat_img),
  landsat_lai(landsat_img),
  landsat_lst(landsat_img),
  landsat_ndvi(landsat_img)]);
input_img = ee.Image(input_img.setMulti({
  'SCENE_ID': landsat_img.get('system:index'),
  'system:time_start': landsat_img.get('system:time_start')
}));
//Map.addLayer(input_img, {}, 'Input Image');


// Compute Tair at the Landsat scale
var tair_img = disalexi_image(input_img);
Map.addLayer(tair_img, {min:273, max:325, 'palette': ['FF0000', 'FFFF00', '00FFFF', '0000FF']}, 'Tair Coarse');


// This asset already exists, but uncomment to rebuild it
// // Export the Landsat scale Tair image
// var asset_id = 'users/cgmorton/disalexi/ta/landsat/' + image_id;
// var task_id = 'disalexi_tair_landsat_' + image_id;
// Export.image.toAsset({
//   image: ee.Image(tair_img).toFloat(),
//   description: task_id,
//   assetId: asset_id,
//   crs: landsat_crs,
//   crsTransform: landsat_geo,
//   dimensions: landsat_shape_str,
// });


// Aggregate Landsat scale Tair to the coarse ALEXI grid
var tair_coarse_img = tair_img
  .reproject(landsat_crs, null, 30)
  .reduceResolution(ee.Reducer.mean(), false, 30000)
  .reproject(output_crs, output_geo)
  .updateMask(1);
Map.addLayer(tair_coarse_img, {min:273, max:325, 'palette': ['FF0000', 'FFFF00', '00FFFF', '0000FF']}, 'Tair Coarse');


// Export the coarse scale Tair image
var asset_id = 'users/cgmorton/disalexi/ta/coarse/' + image_id + '_' + output_cs_str;
var task_id = 'disalexi_tair_coarse_' + image_id + '_' + output_cs_str;
Export.image.toAsset({
  image: ee.Image(tair_coarse_img).toFloat(),
  description: task_id,
  assetId: asset_id,
  crs: output_crs,
  crsTransform: output_geo,
  dimensions: output_shape_str,
});

