# -*- coding: utf-8 -*-
"""
WAPABA rainfall runoff model for python

Robert Morden
University of Melbourne
February 2024

This code has been published as part of the following journal paper:
    Morden R, Horne A, Nathan R, Bond N R (2024), 
    XXXXXXXXXXXXXXXX
    

"""

import pandas as pd
import numpy as np
import math
from numba import jit                                                          # if you are not using NUMBA, just comment this entire line out

# =======================================================================================================================
def WAPABA_wrapper():

    """
    This subroutine is a wrapper to demonstrate the usage of WAPABA.
    
    It loads some basic information for 5 hypothetical catchments and calls the WAPABA engine. 
    It requires basic information about catchments, climate, plus the model parameters.
    
    Please refer to Wang et al (2011) for a detailed discussion of the model.
    
    Wang, Q.J., Pagano, T.C., Zhou, S.L., Hapuarachchi, H.A.P., Zhang, L., Robertson, D.E., 2011.
    Monthly versus daily water balance models in simulating monthly runoff. Journal of Hydrology
    404, 166–175. https://doi.org/10.1016/J.JHYDROL.2011.04.027
    
    Input files
    -----------
    
    catlist     A list of catchments including the following fields:
                    
                    'catID'    : a unique alphanumeric ID for each catchment (should match the
                                 field name for the respective climate data, and model parameters)
                    'area_km2' : area in km2
    
    rain        Catchment rainfall in mm per month.

    pot_et      Catchment evapotranspiration in mm per month.
                The format for both rain and evap should be as follows for the purposes of this wrapper:

                    field 1  : date (month, year) in some sort of python interpretable format
                    field 2  : site 1  rain/evap (field name should match the unique id in catlist, see above)
                    field 3  : site 2  rain/evap
                    etc...
    
    rr_params   Model parameters for each catchment.
                The format for the parameters should be as follows:
                    
                    'catID'  : a unique alphanumeric ID for each catchment (should match the
                                 ID given in 'catlist', refer above)
                    'a1'     : model parameter a1
                    'a2'     : model parameter a2
                    'B'      : model parameter B
                    'Smax'   : model parameter Smax
                    'K'      : model parameter K
                    'Const'  : Special flag. If this field is blank, it is ignored.
                               If this field is set to a positive number, the output will be a monthly
                               timeseries the same length as the input rain, but with all
                               flows set to equal the given value.
                               This is particularly useful for large batch run situations where
                               catchments are connected. Sometimes small catchments do not
                               contribute any meaningful runoff (const=0) or the flow is
                               spring based and can be set to a steady flow of, for example, 5 ML
                               per month (const=5.0).
     """

    # load basic data
    catlist = pd.read_csv('WAPABA_catlist.csv',dtype={'catID':str})                  # database, string index called catID
    catlist.set_index('catID',inplace=True)
    rain = pd.read_csv('WAPABA_rain.csv',index_col=0,parse_dates=True)               # timeseries, date index in first column
    pot_et = pd.read_csv('WAPABA_pot_et.csv',index_col=0,parse_dates=True)           # timeseries, date index in first column
    rr_params = pd.read_csv('WAPABA_rr_params.csv',dtype={'catID':str})              # database, string index called catID
    rr_params.set_index('catID',inplace=True)
    
    # run wapaba
    flowout = WAPABA_setup(catlist,rain,pot_et,rr_params)
    
    return flowout
                   
# =======================================================================================================================
# generate wapaba flows for each FSR
# pass each catchment to numpy routine
def WAPABA_setup(catlist,rain,pet,rr_params):             #(cat,cat_area,rain_series,evap_series,rr_params):
    
    flowML = pd.DataFrame(index=rain.index,columns=rain.columns)
    
    # for each catchment
    for irow in range(len(catlist)):
        
        cat_name = catlist.index[irow]
        cat_area = catlist.loc[cat_name,'area_km2']
        
        # convert rain and evap
        rain_array = rain[cat_name].to_numpy()
        evap_array = pet[cat_name].to_numpy()
        month_array = rain.index.month.to_numpy()    # need the month number, so that the number of days per month can be determined
        
        # get parameters
        Smax = rr_params.loc[cat_name,'Smax']
        a1 = rr_params.loc[cat_name,'a1']
        a2 = rr_params.loc[cat_name,'a2']
        B = rr_params.loc[cat_name,'B']
        K = rr_params.loc[cat_name,'K']
        const = rr_params.loc[cat_name,'Const']
        
        # run model
        if const>=0.0:           # cases where flow is set to a constant
            flowML[cat_name] = const
        else:
            flowmm_array = WAPABA_engine(rain_array,evap_array,month_array,Smax,a1,a2,B,K)
        
            # convert to ML
            flowML[cat_name] = flowmm_array * cat_area
        
    return flowML

# =======================================================================================================================
# run wapaba model using pure numpy
# references to equations are from Wang et al (2011)

@jit(nopython=True) 
def WAPABA_engine(rain,evap,month,Smax,a1,a2,B,K):

    # set up flowout array same as rainfall but with all nans
    flowout = np.empty_like(rain)
    flowout[:] = np.nan
    
    # start and end of monthly loop
    istart = 0
    iend = len(rain)
    monthdays = {1:31.0, 2:28.25, 3:31.0, 4:30.0, 5:31.0, 6:30.0, 7:31.0, 8:31.0, 9:30.0, 10:31.0, 11:30.0, 12:31.0}        # set up days in month
        
    # initialise soil storage
    Send = Smax / 2
    Gend = 1
    # for each month
    for t in range(istart,iend):
        
        Sstart = Send                                                # intialise timestep
        Gstart = Gend                                                # intialise timestep
        numdays = monthdays[month[t]]
        
        X0 = evap[t] + (Smax - Sstart)                               # eqn (3)
        Fr = 1 + (rain[t]/X0) - ((1 + ((rain[t]/X0)**a1))**(1/a1))   # eqn (1)
        X = X0 * Fr                                                  # eqn (2)
        Y = max(0.0, rain[t] - X)                                    # eqn (4)
        W = Sstart + X                                               # eqn (5)
        Fe = 1 + (W/evap[t]) - ((1 + ((W/evap[t])**a2))**(1/a2))    
        ET = evap[t] * Fe                                            # eqn (6)
        Send = max(0.0,min(Smax, W - ET))                            # eqn (7)
        QsAdd = max(0.0,max(0.0,W-ET)-Smax)                          # additional runoff if soil moisture store is full
        R = B * Y                                                    # eqn (8)
        Qs = max(0.0, Y - R) + QsAdd                                 # eqn (9)
        temp1 = (1-math.exp(-numdays/K))
        temp2 = 1-(K/numdays*temp1)
        Qb = (Gstart * temp1) + (R * temp2)                          # eqn (10)
        Gend = Gstart + R - Qb                                       # eqn (11)
        flowout[t] = Qs + Qb
        
    return flowout

# =======================================================================================================================
flowout = WAPABA_wrapper()