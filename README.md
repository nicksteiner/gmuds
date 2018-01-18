# George Mason University Down Scale (GMUDS)
Python library to subset GMUDS dataset for HiMAT.


# Data Format
Description - The current George Mason University Down Scale (GMUDS) product on the NCCS ADAPT server 
 at the NASA Goddard Space Flight Center is spatially subset and temporally concatenated to create 
 daily collections of hourly surface state variables (**Table 1**). File naming, current variables and 
 descriptions of variable coordinates are described in the following.

## File Naming

himat_gmuds_trishuli_daily_YYYYMMDD.nc4  
 * YYYY - 4 digit year
 * MM - 2 digit month of year
 * DD - 2 digit day of month  

*Example*: File, himat_gmuds_trishuli_daily_20071009.nc4, contains hourly surface state variables for 
the October 4th, 2009.
  
## Current Data Variables

**Table 1**

|   | name | units | variable | description |
|---| ---  | ---   | ---      | ---         |
|1  | Surface Air Temperature   | [K]  | Tad | surface air temperature  |
|2. | Surface Pressure          | [Pa] | Pad | surface pressure |
|3. | Surface Relative Humidity | [%]  | RHd | surface relative humidity |
|4. | Surface Net Downward Longwave Flux  | [W/m^2] | Ld | surface net downward longwave flux |
|5. | Surface Net Downward Shortwave Flux | [W/m^2] | Sd | surface net downward shortwave flux |
|6. | Surface Wind Eastward     | [m/s] | U | surface eastward wind  |
|7. | Surface Wind Northward    | [m/s] | V | surface northward wind |

## Current Data Coordinates  

| -- | name | units | variable | description |
| --- | ---- | ----- | -------- | ----------- |
|1.  | Easting | [m]    | easting | meters east of grid origin |
|2.  | Norting | [m]    | northing | meters north of grid origin |
|3.  | Time    | [hour] | time | hour of day |


## Latitude Longitude Grid
Latitude and longitude coordinates, in decimal degrees, for each grid-cell center are contained in 
**latlon_trishuli_gmu_ds.nc4**.

## GMUDS Grid
```python
crs = '+proj=lcc +lat_1=30 +lat_2=62 +lat_0=0 +lon_0=105 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'
grid_size =  (3000, 2000) # row, col
geotransform =  (-3000000, 1000, 0, 5000000, 0, -1000)  # GDAL
```
## Trishuli GMUDS Subset
```python
crs = '+proj=lcc +lat_1=30 +lat_2=62 +lat_0=0 +lon_0=105 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'
grid_size =  (3000, 2000) # row, col
## TO GRID EDGES
geotransform =  (-3000000, 1000, 0, 5000000, 0, -1000)
## TO GRID CENTERS

```

## Python

### Requires
xarray


### Usage

NetCDF files are generated to be compatible with the [xarray](http://xarray.pydata.org/en/stable/) python package.
 
 ```python
import xarray as xr

# open as xarray
ds_file = 'himat_gmuds_trishuli_daily_20080522.nc4'
ds = xr.open_dataset(ds_file)
print(ds)
```

```text
Dimensions:   (time: 24, x: 168, y: 221)
Coordinates:
    easting   (x) float64 ...
    northing  (y) float64 ...
  * time      (time) datetime64[ns] 2008-05-22 2008-05-22T01:00:00 ...
Dimensions without coordinates: x, y
Data variables:
    Tad       (y, x, time) float64 ...
    Pad       (y, x, time) float64 ...
    RHd       (y, x, time) float64 ...
    Ld        (y, x, time) float64 ...
    Sd        (y, x, time) float64 ...
    U         (y, x, time) float64 ...
    V         (y, x, time) float64 ..."""
```

```python
# calculate stats
avg_temp  = ds.Tad.mean(dim='time')
min_temp  = ds.Tad.min(dim='time')
max_temp  = ds.Tad.max(dim='time')

# mean from 10 to 2 pm
tem_midday = ds.sel(time=slice('2008-05-22T10:00:00', '2008-05-22T14:00:00')).mean(dim='time')

# combine days
ds_file_next = 'himat_gmuds_trishuli_daily_20080523.nc4'
ds_next = xr.open_dataset(ds_file_next)
ds_combine = xr.concat([ds, ds_next])

# mean from 10 to 2 am
tem_midngt = ds_combine.sel(time=slice('2008-05-22T22:00:00', '2008-05-23T02:00:00')).mean(dim='time')


```

# Notes
## From Email Chain 9/25/16 - Based on Vinod's Input

### Primary Dynamic
* Precipitation <  
* Air Temperature 1   
* Wind Speed 6+7  
* Relative humidity 3   
* Incoming SWR 5 
* Incoming LWR 4   
* Pressure 2

### Secondary Dynamic
* Albedo, Glacier Type etc.  

### Static
* Elevation   
* Leaf Area Index  
* Vegetation Cover  
* Vegetation Height  
* Canopy Cover Fraction
* Soil Bulk Density  

### For Model Validation
* SWE/SWE Change  
* Snow/Glacier Thickness Change  
* Snow/Glacier Temperature  
* Snow/Glacier Melt  
* Evaporation/Sublimation Measurements


## From Rijan - MPDDM Version 2-Rijan.ppx, Slide 1

### Dynamic

* Temperature
* Precipitation
* Lapse Rate
* Precipitation Gradient
* Critical Temperature
* Degree Day Factor

### Static

* Elevation   
* Landuse Class
* Clean Glacier
* Debris
* In situ measurement of flow

# Matching GMUDS Data Available

| --- | name | units |  collection  | shrtname | folder | variable |
| --- | ---- | ----- |  ----------  | -------- | ------ | -------- | 
|1.  | Surface Air Temperature | [K]  |  M2T1NXFLX  |  TLML   |Ambient| Tad | 
|2.  | Surface Pressure |[Pa] | M2T1NXSLV |PS|  Ambient| Pad |
|3.  | Surface Relative Humidity |[%]| M2T1NXFLX | QLML  |Ambient| RHd | 
|4.  | Surface Net Downward Longwave Flux |[W/m^2] | M2T1NXINT | LWGNET | Ambient | Ld |
|5.  | Surface Net Downward Shortwave Flux |[W/m^2] | M2T1NXINT | SWNETSRF|SW| Sd |
|6.  | Surface Wind Eastward |[m/s] | M2T1NXFLX | ULML |UVWind| U |
|7.  | Surface Wind Northward |[m/s] | M2T1NXFLX | VLML  |UVWind| V |
|8.  | Convective Rainfall |[kg/m^2s] | M2T1NXINT | PRECCU  |N/A|
|9.  | Large Scale Rainfall |[kg/m^2s] | M2T1NXINT | PRECLS  |N/A|
|10. | Snowfall |[kg/m^2s] | M2T1NXINT | PRECSN  |N/A|
|11. | Total Precipitation |[kg/m^2s] | M2T1NXINT | PRECTOT  |N/A|
|12. | Total Precipitation - Bias Corrected|[kg/m^2s] | M2T1NXINT | PRECTOTCORR  |N/A|



