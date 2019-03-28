# Weather Data Processor
#### Converts gridded binary (.grb) files to dataframes for analysis, and hdf5 for storage 

### Modules

datablock.py - This file contains functions and classes for formatting weather data in dataframes "wxblocks". These wxblocks are multi-index frames that can compute and store weather balloon data with pressure (re: geopotential height) levels, or time-series data with datetime levels. It includes functions for using balloon sounding data to compute severe weather predictors.

dataform.py - This file contains functions and classes for downloading and processing gridded binary (.grb) files from the NOAA website

https://nomads.ncdc.noaa.gov/data

### Notebook for validating models with sounding data

sounding_validation.ipynb

### Notebook for running ETL processes

etl_gribdata.ipynb
