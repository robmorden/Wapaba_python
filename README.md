# WAPABA_python
Python implementation of the rainfall runoff model WAPABA.

This repository includes code for the rainfall runoff model WAPABA (Spatial tool for estimating dam impacts). For further information about the purpose of this code, please refer to the following two papers:

> Morden R, Horne A, Nathan R, Bond N R (2024), XXXXXXXXXXXXXXXX

> Wang, Q.J., Pagano, T.C., Zhou, S.L., Hapuarachchi, H.A.P., Zhang, L., Robertson, D.E., 2011. Monthly versus daily water balance models in simulating monthly runoff. Journal of Hydrology 404, 166–175. https://doi.org/10.1016/J.JHYDROL.2011.04.027

## Installation
I am a complete GitHub newbie as of August 2022. There is no prepared package to install, this is simply a repository of my code. If it looks useful, please download it and reference this repository or my paper Morden et al (2024).

## Overview
Wapaba is a rainfall runoff model which operates on a monthly timestep.

This code has been prepared with a short 'wrapper' function to allow it to load and run sample data for demonstration purposes. The sample data includes basic information for 5 hypothetical catchments.

## Inputs
The model requires basic information about catchments, climate, plus the model parameters. Each input is discussed below.
    
### List of catchments   
See sample data in the file `WAPABA_catlist.csv`.

The list of catchments must include the following fields:
                    
* `catID`      a unique alphanumeric ID for each catchment (should match the field name for the respective climate data, and model parameters)
* `area_km2`   area in km2
    
### Rainfall
Catchment rainfall in mm per month. See sample data in the file `WAPABA_rain.csv`.

### Potential evapotranspiration
Catchment evapotranspiration in mm per month. See sample data in the file `WAPABA_pot_et.csv`.

The format for both rain and evap should be as follows for the purposes of this code:

* `field 1`    date (month, year) in some sort of python interpretable format such as `1895-12-01`.
* `field 2`    site 1 rain/evap (field name should match the unique id in catlist, see above)
* `field 3`    site 2 rain/evap
etc...
    
### Model parameters
See sample data in the file `WAPABA_rr_params.csv`.

The format for the model parameters should be as follows:
                    
* `catID`    a unique alphanumeric ID for each catchment (should match the ID given in 'catlist', refer above)
* `a1`       model parameter a1
* `a2`       model parameter a2
* `B`        model parameter B
* `Smax`     model parameter Smax
* `K`        model parameter K
* `Const`    Special flag. If this field is blank, it is ignored. If this field is set to a positive number, the output will be a monthly timeseries the same length as the input rain, but with all flows set to equal the given value.

The `Const` flag is particularly useful for large batch run situations where catchments are connected/nested. Sometimes small catchments do not contribute any meaningful runoff (const=0) or the flow is spring based and can be set to a steady flow of, for example, 5 ML per month (const=5.0).
