# Coarse Resolution Air Temperature Assets

The coarse resolution air temperature assets will need to be generated for each Landsat image before DisALEXI can be run.

Note, the provided "ini" file is setup to generate the air temperature assets for the full CONUS.  It is recommended to create a separate INI file for specific study areas or date ranges of interest.

```
python tair_image_wrs2_export.py --project <YOUR_GCP_PROJECT_NAME> -i ta_wrs2_conus_v006.ini
```

# ALEXI Ancillary Assets


