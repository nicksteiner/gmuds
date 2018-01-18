# George Mason University Down Scale (GMUDS)
Python library to subset GMUDS dataset for HiMAT.

## Requires
xarray, pyproj, affine

## Usage



## GMUDS Grid
```
CRS: '+proj=lcc +lat_1=30 +lat_2=62 +lat_0=0 +lon_0=105 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'
Grid Size: (3000, 2000) # row, col
Geotransform: (-3000000, 1000, 0, 5000000, 0, -1000)  # GDAL

```

# Project Needs
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
|3.  | Surface Relative Humidity |[K]| M2T1NXFLX | QLML  |Ambient| RHd | 
|4.  | Surface Net Downward Longwave Flux |[W/m^2] | M2T1NXINT | LWGNET | Ambient | Ld |
|5.  | Surface Net Downward Shortwave Flux |[W/m^2] | M2T1NXINT | SWNETSRF|SW| Sd |
|6.  | Surface Wind Eastward |[m/s] | M2T1NXFLX | ULML |UVWind| U |
|7.  | Surface Wind Northward |[m/s] | M2T1NXFLX | VLML  |UVWind| V |
|8.  | Convective Rainfall |[kg/m^2s] | M2T1NXINT | PRECCU  |N/A|
|9.  | Large Scale Rainfall |[kg/m^2s] | M2T1NXINT | PRECLS  |N/A|
|10. | Snowfall |[kg/m^2s] | M2T1NXINT | PRECSN  |N/A|
|11. | Total Precipitation |[kg/m^2s] | M2T1NXINT | PRECTOT  |N/A|
|12. | Total Precipitation - Bias Corrected|[kg/m^2s] | M2T1NXINT | PRECTOTCORR  |N/A|



