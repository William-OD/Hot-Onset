# William O'Donnell
# Date: March 2022
# 
# Functions: 
#   data_download(fl_data)  - fetches the GOES data for a given flare, given the start and end times of the flare
#   
#
#
#
#
#
#
#
# 
#
# ----------------------------------------------------------------------------------------------
import numpy as np
import pandas as pd
from astropy.time import Time, TimeDelta
# ----------------------------------------------------------------------------------------------

import re
from datetime import time
from datetime import datetime

from sunpy.net import attrs as a
from sunpy.net import Fido
from scipy import signal
import goesxrs_temp as gtem

def data_download(fl_data):
    #Imports one row of a table with start, peak, end times..
    # Want data for +/-10 mins of flare, using astopy time to do this (easier sunpy way?)
    gtstart=Time(fl_data["event_starttime"],scale='utc')-TimeDelta(10*60,format='sec')
    gtend=Time(fl_data["event_endtime"],scale='utc')+TimeDelta(10*60,format='sec')

    # Search and get the XRS data
    tflrange=a.Time(gtstart.iso,gtend.iso)
    rg15 = Fido.search(tflrange, a.Instrument("XRS"), a.goes.SatelliteNumber(15))
    fg15 = Fido.fetch(rg15,path=data_dir)


# Function to truncate the data and output truncated data in 2s and 10s sampling
# Modified to have the first window begin at 5 mins after the posted start time.
def bck_truncate(short_data,long_data,fl_startt,startt = 60, endt = 0):
    bckstart=Time(fl_startt,scale='utc')+TimeDelta(300,format='sec')-TimeDelta(startt,format='sec')
    bckend=Time(fl_startt,scale='utc')+TimeDelta(300,format='sec')-TimeDelta(endt, format='sec')
    bcktime=a.Time(bckstart.iso,bckend.iso)
    bcktrunc_054= short_data.truncate(bcktime.start.iso,bcktime.end.iso)
    bcktrunc_18= long_data.truncate(bcktime.start.iso,bcktime.end.iso)
    return bcktime, bcktrunc_054, bcktrunc_18


# Function to take the raw data and determine a suitable 1-minute background based on 
# variance fluctuations
def background(data_054, data_18, start_time, peak_time):
    backg_flag = False

    srch_start=Time(start_time, scale='utc')-TimeDelta(300,format='sec')
    srch_end=Time(peak_time,scale='utc')
    srch_time=a.Time(srch_start.iso, srch_end.iso)
    bcktrunc_054=data_054.truncate(srch_time.start.iso,srch_time.end.iso)
    bcktrunc_18=data_18.truncate(srch_time.start.iso,srch_time.end.iso)
    bck_t = bcktrunc_054.index
    
    vars = []
    for i in range(0,len(bcktrunc_054)):
        var = np.var(bcktrunc_054[i-30:i])
        vars.append(var)

    smooth_df = pd.DataFrame({'Variance':vars, 'Flux':bcktrunc_054})
    smooth_df.set_index(bck_t)
    smooth_df.index = pd.to_datetime(smooth_df.index)
    smooth_df = smooth_df.resample('10S').mean()
    df = pd.DataFrame({'Time': smooth_df.index, 'Variance': smooth_df['Variance'].values, 'Flux': smooth_df['Flux'].values})
    df = df.dropna()
    df.reset_index()

    # Finding Troughs
    trs, _ = signal.find_peaks(np.negative(df['Variance']))
    factor = 4 #Range number for looking for background interval

    # Finding the min trough and any troughs within a range till double the variance of that min trough.
    if len(trs) >= 1:
        tr_time = []
        tr_var = []
        for t in trs:
            if df.iloc[t]['Variance'] <= min(df['Variance'].values)*factor and df.iloc[t]['Flux'] <= df.iloc[min(df['Variance'].index)]['Flux']*2:
                tr_time.append(df.iloc[t]['Time'])
                tr_var.append(df.iloc[t]['Variance'])
            else:
                continue
        if len(tr_time) >= 1:
        # Selecting the interval to be the latest min trough in that range
            bck_start = Time(max(tr_time), scale='utc') - TimeDelta(60,format='sec')
            bck_end = Time(max(tr_time), scale='utc')

        else:
            print("---------------------")
            print("No Suitable Background Detected... Using a 1 min window starting from start time")
            print("---------------------")
            bck_start = (Time(start_time))
            bck_end = (Time(start_time) + TimeDelta(60, format = 'sec'))
            backg_flag = True
    else:
        print("---------------------")
        print("Not enough troughs located!")
        print("---------------------")
        tr_time = np.nan
        tr_var = np.nan
        bck_start = (Time(start_time))
        bck_end = (Time(start_time) + TimeDelta(60, format = 'sec'))
        backg_flag = True

    return bck_start.datetime, bck_end.datetime, smooth_df, tr_time, tr_var, backg_flag


# Function to take the start of a fixed onset (either the end of background or 1/8th into the impulsive phase)
# and calculate the end time of a fixed fraction of the impulsive phase based on n.
def scaled_onset(onset_start, true_peakt, n):   
    onset_end = Time(onset_start) + (Time(true_peakt)-Time(onset_start))/n
    return onset_end.datetime


# Function to determine the average flux of a pre-calculated onset interval for both channels
# Then uses goesxrs.temp.py fucntion to determine the temperature and emission measure of this flux value.
def onset_tem(short_raw, long_raw, short_backsub, long_backsub, onset_start, onset_end, bck_short_std, bck_long_std): 
    trunc_054_flux_bck = short_backsub.truncate(onset_start, onset_end)
    trunc_18_flux_bck = long_backsub.truncate(onset_start, onset_end)
    trunc_054_flux_raw = short_raw.truncate(onset_start, onset_end)
    trunc_18_flux_raw = long_raw.truncate(onset_start, onset_end)

#Uncertainty on the raw unsubbed flux measurement   
    counts_unc_054 = np.sqrt(trunc_054_flux_raw.values)
    counts_unc_18 = np.sqrt(trunc_18_flux_raw.values)
#Uncertainty on the backsubbed flux
    flux_unc_054 = np.sqrt(counts_unc_054**2 + bck_short_std) 
    flux_unc_18 = np.sqrt(counts_unc_18**2 + bck_long_std)

    if len(trunc_054_flux_bck) and len(trunc_18_flux_bck) >= 1:
    #Calculating the weighted mean flux
        weights_054 = 1/(flux_unc_054)**2
        weights_18 = 1/(flux_unc_18)**2
        weighted_mean_054 = np.sum(weights_054 * trunc_054_flux_bck) / np.sum(weights_054)
        weighted_mean_18 = np.sum(weights_18 * trunc_18_flux_bck) / np.sum(weights_18)

    #Calculating the error on each weighted mean flux
        err_054 = np.sqrt(np.sum(weights_054 * (trunc_054_flux_bck - weighted_mean_054)**2) / (np.sum(weights_054) * (len(trunc_054_flux_bck)) - 1))
        err_18 = np.sqrt(np.sum(weights_18 * (trunc_18_flux_bck - weighted_mean_18)**2) / (np.sum(weights_18) * (len(trunc_18_flux_bck)) - 1))

    else:
        weighted_mean_054 = trunc_054_flux_bck
        weighted_mean_18 = trunc_18_flux_bck
        err_054 = flux_unc_054
        err_18 = flux_unc_18


#Converting to TEM (including uncertainties)
    tmk_onset, em_onset, tmk_upper, tmk_lower, em_err = gtem.get_tem(weighted_mean_18, weighted_mean_054, err_18, err_054) 

    return tmk_onset, em_onset, tmk_upper, tmk_lower, em_err


# Fancy Onset Calculation Function
def onset_times(df_short_bcksub, background_end, peak_time, bck_sigma):
   
    onset_flag = False
    
    #Start by truncating short channel data from end of background interval - 1 min to peak time.
    srch_start = Time(background_end) - TimeDelta(60, format = 'sec')
    trange_search = a.Time(srch_start, peak_time)
    g_short_srch = df_short_bcksub.truncate(trange_search.start.iso,trange_search.end.iso)
    g_tims_srch = g_short_srch.index
    
    #Calculating the cumulative variance at each point
    vars = []
    for i in range(0,len(g_short_srch)):
        var = np.var(g_short_srch[i-30:i])
        vars.append(var)

    smooth = pd.Series(vars, index = pd.DatetimeIndex(g_tims_srch))
    smooth = smooth.resample('10S').mean()
    df = pd.DataFrame({'Time': smooth.index, 'Variance': smooth.values})

    #Calculating the turning points through 'find peaks'
    peaks, _ = signal.find_peaks(df['Variance'])
    troughs, _ = signal.find_peaks(np.negative(df['Variance']))
    
    peak_tims = []
    for p in peaks:
        peak_tims.append(df.iloc[p]['Time'])

    if len(troughs) >= 2:
        trough_time = []
        trough_var =[]
        for t in troughs:
            trough_time.append(df.iloc[t]['Time'])
            trough_var.append(df.iloc[t]['Variance'])

         #Find first point above 5-sigma, and if not, choosing first peak in variance
        try:
            onset_find = g_tims_srch[(next(i for i, v in enumerate(g_short_srch) if v >= 5*bck_sigma))]
        except:
            print("Using First Peak as 5sigma intercept had an error!")
            onset_find = df.iloc[peaks[0]]['Time']

        first_peak_time = min(peak_tims, key=lambda sub: abs(sub - onset_find))

        for t in range(0, len(trough_time)):
            if trough_time[t] < first_peak_time:
                continue
            else:
                onset_end = trough_time[t]
                if t >= 1:
                    onset_start = trough_time[t-1]
                else:
                    onset_start = (Time(background_end) + TimeDelta(2, format = 'sec')).datetime
                break

    else:
        print("")
        print("Not enough troughs for located! To calculate onset!")
        print("Going to caculate the onset time to be the first 1/8 of the Impulsive Phase.")
        onset_start = (Time(background_end) + TimeDelta(2, format = 'sec')).datetime
        onset_end = scaled_onset(Time(background_end)+ TimeDelta(2, format = 'sec'), peak_time, 8)
        onset_find = np.nan
        trough_time = np.nan
        trough_var = np.nan

        onset_flag = True

    try:
        if onset_end >= peak_time:
            onset_end = (Time(peak_time) - TimeDelta(2, format = 'sec')).datetime
        return onset_start, onset_end, onset_flag
    except:
        return np.nan, np.nan