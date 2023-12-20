# Functions to get T and EM info out of the GOES/XRS data
# Original code from Iain Hannah's GitHub Repository: https://github.com/ianan/xrs_example/blob/main/goesxrs_temp.py
# which is based off of https://hesperia.gsfc.nasa.gov/ssw/gen/idl/synoptic/goes/goes_tem_calc.pro
# and specifically the CHIANTI versions via
# https://hesperia.gsfc.nasa.gov/ssw/gen/idl/synoptic/goes/goes_get_chianti_temp.pro
# and
# https://hesperia.gsfc.nasa.gov/ssw/gen/idl/synoptic/goes/goes_get_chianti_em.pro
# 
# Response calculated in sswidl:
#   https://hesperia.gsfc.nasa.gov/ssw/gen/idl/synoptic/goes/buildresponse/goes_chianti_response.pro
# the response fits is always called "goes_chianti_response_latest.fits" now
# Version have here is from 21-05-2022 and v10 CHIANTI:
#   https://hesperia.gsfc.nasa.gov/ssw/gen/idl/synoptic/goes/goes_chianti_response_latest.fits
# Older V9 CHIANTI one also in repo "goes_chianti_resp_20200812.fits", but have wrong factor in GOES13/15 short
#  
# 
# Values returned << 1% out from sswidl versions, probably due to slight difference in the spline interpolation
# 
# 24-10-2021    IGH    Created
# 25-10-2021    IGH    Updated to return newer CHIANTI v9 responses, and option of abundance and sat
# 23-05-2022    IGH    Updated with new response (CHIANTI v10 and G13/15 short fix)
#                      Default loads new response, but old_ver=True gives 20200812 version
# 18-07-2022    IGH    Fixed get_tem() bug that loaded default response only for the EM calc
#
# 16-03-2023    WEO    Added temperature and emission measure uncertainty calculations
# 16-12-2023    WEO    Added in docstring documentation and additional comments
#-----------------------------------------------
#-----------------------------------------------

import numpy as np
from scipy import interpolate
from astropy.io import fits
#-----------------------------------------------
#-----------------------------------------------

def get_resps(sat=15,cor_not_pho=True,old_ver=False):
    """
    Returns the GOES/XRS temperature response functions for long and short channels.

    The response is calculated in sswidl:
    https://hesperia.gsfc.nasa.gov/ssw/gen/idl/synoptic/goes/buildresponse/goes_chianti_response.pro
    The response fits is always called "goes_chianti_response_latest.fits" now.
    Version have here is from 21-05-2022 and v10 CHIANTI:
    https://hesperia.gsfc.nasa.gov/ssw/gen/idl/synoptic/goes/goes_chianti_response_latest.fits

    Parameters:
    sat (int, optional): Which GOES satellite to use? Defaults to 15.
    cor_not_pho (bool, optional): Use coronal not photospheric abundances. Defaults to True.
    old_ver (bool, optional): Use previous responses of 20200812 + wrong G13/15 short. Defaults to False.

    Returns:
    resps (numpy.ndarray): Temperature response for 1-8\AA and 0.5-4\AA in units of 10^{-55} W m^{-2} cm^{3}.
    resptmk (numpy.ndarray): Temperature binning of TR ratio in MK.
    """

    # Check if old version is to be used
    if old_ver:
        # Use old response file
        rfile='Response_Functions/goes_chianti_resp_20200812.fits' # NOTE: Only works if files are located in Repsonse_Functions subdirectory
    else:
        # Use latest response file
        rfile='Response_Functions/goes_chianti_response_latest.fits' # NOTE: Only works if files are located in Repsonse_Functions subdirectory
    
    # Open the response file
    hdulist = fits.open(rfile)
    # Get the data from the file
    respdat=hdulist[1].data
    # Close the file
    hdulist.close()

    # Get the temperature data for the specified satellite
    resptmk=np.array(respdat["TEMP_MK"][sat-1])
    
    # Check if coronal abundances are to be used
    if cor_not_pho:
        abdun="COR"
    else:
        # Use photospheric abundances
        abdun="PHO"  
    
    # Initialize response array
    resps=np.empty((101,2))
    # Get the long and short channel responses
    resps[:,0]=respdat["FLONG_"+abdun][sat-1]
    resps[:,1]=respdat["FSHORT_"+abdun][sat-1]
   
    # Return the responses and temperature data
    return resps, resptmk
#-----------------------------------------------
#-----------------------------------------------

def get_tem(fl,fs,fl_err,fs_err,sat=15,cor_not_pho=True,old_ver=False):
    """
    Returns the temperature and emission measure for ratio of fluxes in short/long of GOES/XRS channels.

    This function is based on the script at https://hesperia.gsfc.nasa.gov/ssw/gen/idl/synoptic/goes/goes_tem_calc.pro
    The response is calculated in sswidl:
    https://hesperia.gsfc.nasa.gov/ssw/gen/idl/synoptic/goes/buildresponse/goes_chianti_response.pro
    The response fits is always called "goes_chianti_response_latest.fits" now.
    Version have here is from 21-05-2022 and v10 CHIANTI:
    https://hesperia.gsfc.nasa.gov/ssw/gen/idl/synoptic/goes/goes_chianti_response_latest.fits

    Parameters:
    fl (float or array): GOES/XRS flux in long channel, 1-8\AA no scaling, Wm^{-2}.
    fs (float or array): GOES/XRS flux in short channel,0.5-4\AA scaling, Wm^{-2}.
    fl_err (float or array): Uncertainty in the long channel flux.
    fs_err (float or array): Uncertainty in the short channel flux.
    sat (int, optional): Which GOES satellite to use? Defaults to 15.
    cor_not_pho (bool, optional): Use coronal not photospheric abundances. Defaults to True.
    old_ver (bool, optional): Use previous responses of 20200812 + wrong G13/15 short. Defaults to False.

    Returns:
    tmk (numpy.ndarray): Temperature in MK.
    em (numpy.ndarray): Emission Measure in cm^{-3}.
    tmk_upper (numpy.ndarray): Upper bound of the temperature in MK.
    tmk_lower (numpy.ndarray): Lower bound of the temperature in MK.
    em_err (numpy.ndarray): Uncertainty in the emission measure.
    """
#   NOTE -  No longer need scaling for GOES16/17 as was correct but GOES13/15 short resp wrong (now fixed)

#   Get the TR to work out the EM
    resps, resptmk=get_resps(sat,cor_not_pho,old_ver)

#   Then calculate the resps ratio
    resprat=resps[:,1]/resps[:,0]
    
#   Use scipy cubic spline interpolation
    rat_func=interpolate.interp1d(resprat, resptmk,kind='cubic')
    
#   Ratio can't be outside of values of TR ratio, so
    gfs=np.array(fs)
    gfl=np.array(fl)
    grat=np.array(gfs/gfl)
    grat[grat < resprat.min()]=np.nan 
    grat[grat > resprat.max()]=np.nan 

# Calculating upper and lower bounds for the ratios based on uncertainties
    gfs_errs = np.array(fs_err)
    gfl_errs = np.array(fl_err)
    grat_err = grat * np.sqrt((gfs_errs/gfs)**2 + (gfl_errs/gfl)**2) # Propagating errors properly
    grat_err_upper = grat + grat_err
    grat_err_lower = grat - grat_err

#   Work out the temperature 
    tmk=np.array(rat_func(grat))
    tmk_upper = np.array(rat_func(grat_err_upper))
    tmk_lower = np.array(rat_func(grat_err_lower))
    
#   Use scipy cubic spline interpolation
    tr18_func=interpolate.interp1d(resptmk,resps[:,0],kind='cubic')
#   Can't use TMK values at/outside the range, so
    tmk[tmk < resptmk.min()]=np.nan
    tmk[tmk > resptmk.max()]=np.nan

#   Work out the Emission Measure
    em=np.array(1e55*gfl/tr18_func(tmk))
    # em[em < 0]=0.
    # em[tmk <= 1.0001]=0.
    em[em < 0]=np.nan
    em[tmk <= 1.0001]=np.nan

# Calculating upper and lower bounds for the EM based on uncertainties
    resp_upper = np.array(tr18_func(tmk_upper)) 
    resp_lower = np.array(tr18_func(tmk_lower))
    resp_avg = (resp_upper - resp_lower)/2 # Averaging upper and lower limits of response to get an uncertainty
    em_err = em*np.sqrt((resp_avg/tr18_func(tmk))**2 + (gfl_errs/gfl)**2) #propagating uncertainty
    
    return tmk, em, tmk_upper, tmk_lower, em_err
#-----------------------------------------------
#-----------------------------------------------

def get_tem_old(fl,fs,sat=15,cor_not_pho=True,old_ver=False):
    
    """
    Returns the temperature and emission measure for ratio of fluxes in short/long of GOES/XRS channels.

    This function is based on the script at https://hesperia.gsfc.nasa.gov/ssw/gen/idl/synoptic/goes/goes_tem_calc.pro
    The response is calculated in sswidl:
    https://hesperia.gsfc.nasa.gov/ssw/gen/idl/synoptic/goes/buildresponse/goes_chianti_response.pro
    The response fits is always called "goes_chianti_response_latest.fits" now.
    Version have here is from 21-05-2022 and v10 CHIANTI:
    https://hesperia.gsfc.nasa.gov/ssw/gen/idl/synoptic/goes/goes_chianti_response_latest.fits

    Parameters:
    fl (float or array): GOES/XRS flux in long channel, 1-8\AA no scaling, Wm^{-2}.
    fs (float or array): GOES/XRS flux in short channel,0.5-4\AA scaling, Wm^{-2}.
    sat (int, optional): Which GOES satellite to use? Defaults to 15.
    cor_not_pho (bool, optional): Use coronal not photospheric abundances. Defaults to True.
    old_ver (bool, optional): Use previous responses of 20200812 + wrong G13/15 short. Defaults to False.

    Returns:
    tmk (numpy.ndarray): Temperature in MK.
    em (numpy.ndarray): Emission Measure in cm^{-3}.
    """
#   NOTE -  No longer need scaling for GOES16/17 as was correct but GOES13/15 short resp wrong (now fixed)
# 

#   Get the TR to work out the EM
    resps, resptmk=get_resps(sat,cor_not_pho,old_ver)

#   Then calculate the resps ratio
    resprat=resps[:,1]/resps[:,0]
    
#   Use scipy cubic spline interpolation
    rat_func=interpolate.interp1d(resprat, resptmk,kind='cubic')
    
#   Ratio can't be outside of values of TR ratio, so
    gfs=np.array(fs)
    gfl=np.array(fl)
    grat=np.array(gfs/gfl)
    grat[grat < resprat.min()]=np.nan #convert to np.nan??
    grat[grat > resprat.max()]=np.nan #convert to np.nan??


#   Work out the temperature 
    tmk=np.array(rat_func(grat))
    
#   Use scipy cubic spline interpolation
    tr18_func=interpolate.interp1d(resptmk,resps[:,0],kind='cubic')
#   Can't use TMK values at/outside the range, so
#   (might this causes issues - better way of doing this??)
    # tmk[tmk < resptmk.min()]=1.0001
    # tmk[tmk >resptmk.max()]=resptmk.max()
    tmk[tmk < resptmk.min()]=np.nan
    tmk[tmk > resptmk.max()]=np.nan

#   Work out the Emission Measure
    em=np.array(1e55*gfl/tr18_func(tmk))
    # em[em < 0]=0.
    # em[tmk <= 1.0001]=0.
    em[em < 0]=np.nan
    em[tmk <= 1.0001]=np.nan

    return tmk, em
#-----------------------------------------------
#-----------------------------------------------

def get_resprat(sat=15,cor_not_pho=True,old_ver=False):
    """
    Returns the ratio of temperature response and temperature binning in MK for GOES/XRS channels.

    The response is calculated in sswidl:
    https://hesperia.gsfc.nasa.gov/ssw/gen/idl/synoptic/goes/buildresponse/goes_chianti_response.pro
    The response fits is always called "goes_chianti_response_latest.fits" now.
    Version have here is from 21-05-2022 and v10 CHIANTI:
    https://hesperia.gsfc.nasa.gov/ssw/gen/idl/synoptic/goes/goes_chianti_response_latest.fits

    Parameters:
    sat (int, optional): Which GOES satellite to use? Defaults to 15.
    cor_not_pho (bool, optional): Use coronal not photospheric abundances. Defaults to True.
    old_ver (bool, optional): Use previous responses of 20200812 + wrong G13/15 short. Defaults to False.

    Returns:
    resprat (numpy.ndarray): Ratio of TR_(0.5-4) / TR_(1-8).
    resptmk (numpy.ndarray): Temperature binning of TR ratio in MK.
    """
    
    if old_ver:
        rfile='goes_chianti_resp_20200812.fits'
    else:
        rfile='goes_chianti_response_latest.fits'
    hdulist = fits.open(rfile)
    respdat=hdulist[1].data
    hdulist.close()
    
    resptmk=np.array(respdat["TEMP_MK"][sat-1])
    
    if cor_not_pho:
        abdun="COR"
    else:
        abdun="PHO"   
    resprat=np.empty((101))
    resprat=respdat["FSHORT_"+abdun][sat-1]/respdat["FLONG_"+abdun][sat-1]
    
    return resprat, resptmk
#-----------------------------------------------
#-----------------------------------------------