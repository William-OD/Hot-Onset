import numpy as np
import re
from datetime import time
from datetime import datetime 
import os

from scipy.io import readsav
from scipy import signal

from astropy.time import Time
from astropy.time import TimeDelta

from sunpy import timeseries as ts
from sunpy.net import attrs as a

import pandas as pd

import Scripts.goesxrs_temp as gtem #from Ian's functions
import Scripts.Onsets_temp as onsets

import warnings #Need to remove for future
warnings.filterwarnings("ignore", module="sunpy.timeseries.timeseriesbase")

# data_dir = r"D:\MastersProj\Data\goes15"         #Directory to read in GOES data from
# hek_dir = r"D:\MastersProj\Data\goes15\HEK_Data" #Directory to read in HEK data from
# onset_length = 60                                #Time of fixed onset in seconds


def print_fl_info(fl_hek_start, fl_class):
    """
    Prints information about a flare event.

    Parameters:
    fl_hek_start (str): The start time of the flare event.
    fl_class (str): The class of the flare event.

    Returns:
    None

    """

    print("------------------------------")
    print(f"{fl_class} Flare With start Time: {fl_hek_start}")


def load_timeseries(data_dir, fl_hek_start, fl_hek_end):
    """
    Loads time series data based on the specified directory and flare start/end times.

    Parameters:
        data_dir (str): The directory path where the data is stored.
        fl_hek_start (str): The start time of the flare in UTC format.
        fl_hek_end (str): The end time of the flare in UTC format.

    Returns:
        g15 (ts.TimeSeries): The loaded time series data.
        flare_time (a.Time): The time range of the flare.

    """
    flare_time=a.Time(Time(fl_hek_start, scale='utc').iso, Time(fl_hek_end,scale='utc').iso)
    day_string = flare_time.start.datetime.date().strftime("%Y%m%d")
    yr_str = flare_time.start.datetime.date().strftime("%Y")
    mnth_str = flare_time.start.datetime.date().strftime("%m")
    date_path = os.path.join(yr_str, mnth_str)
    fl_path = os.path.join(data_dir, date_path,f'*sci_gxrs-l2-irrad_g15_d*{day_string}*.nc').replace("\\", "/")
    print(fl_path)

    if Time(fl_hek_start).datetime.time() <= time(hour=1, minute=0, second=0):
        day_str_prev = (Time(fl_hek_start) - TimeDelta(1, format = 'jd')).datetime.date().strftime("%Y%m%d")
        fl_path_prev = os.path.join(data_dir, date_path,f'*sci_gxrs-l2-irrad_g15_d*{day_str_prev}*.nc')
        g15 = ts.TimeSeries([fl_path, fl_path_prev],concatenate=True)
    elif Time(fl_hek_start).datetime.time() >= time(hour=23, minute=0, second=0):
        day_str_post = (Time(fl_hek_start)+ TimeDelta(1, format = 'jd')).datetime.date().strftime("%Y%m%d")
        fl_path_post = os.path.join(data_dir, date_path,f'*sci_gxrs-l2-irrad_g15_d*{day_str_post}*.nc')
        g15 = ts.TimeSeries([fl_path, fl_path_post],concatenate=True)
    else:
        g15 = ts.TimeSeries(fl_path, concatenate=True)
    
    return g15, flare_time


def pre_process_timeseries(g15):
    """
    Preprocesses the flare data into a time series format.

    Args:
        g15 (TimeSeries): A time series object containing the flare data.
        flare_time (TimeRange): The time range of the flare.

    Returns:
        tuple: A tuple containing the preprocessed short and long wavelength data.
            - df_short (Series): A pandas series object containing the preprocessed short wavelength data.
            - df_long (Series): A pandas series object containing the preprocessed long wavelength data.
    """
# Load flare into time series
    print("Test1")
    g_short = g15.quantity("xrsa").value 
    g_long = g15.quantity("xrsb").value
    g_short_flag = g15.quantity("xrsa_quality").value 
    g_long_flag = g15.quantity("xrsb_quality").value 

# Create a Boolean mask where any flags are located
    print("Test2")
    mask_long = g_long_flag != 0
    mask_short = g_short_flag != 0

# Set the corresponding values in g_short and g_long to NaN
    print("Test3")
    g_short[mask_short] = np.nan
    g_long[mask_long] = np.nan

# Set sub-zero values not picked up by flags to nans
    print("Test4")
    mask_short_mask = g_short <= 0 
    g_short[mask_short_mask] = np.nan 

# Load into dataframes
    print("Test5")
    df_long = pd.Series(g_long, index = pd.DatetimeIndex(g15.data.index))
    df_short = pd.Series(g_short, index = pd.DatetimeIndex(g15.data.index))
    print("Test6")

    return df_short, df_long


def backsub_timeseries(df_short, df_long, fl_hek_start, flare_time):
    """
    Subtract the background from the given time series data.

    Parameters:
    - df_short (pandas.DataFrame): Short cadence time series data.
    - df_long (pandas.DataFrame): Long cadence time series data.
    - fl_hek_start (str): Start time of the flare in HEK format.
    - flare_time (TimeRange): Time range of the flare.

    Returns:
    - short_backsub (pandas.DataFrame): Short cadence data with background subtracted.
    - long_backsub (pandas.DataFrame): Long cadence data with background subtracted.
    - bck_startt (str): Start time of the background.
    - bck_endt (str): End time of the background.
    - bck_short_std (float): Standard deviation of the short cadence background.
    - bck_long_std (float): Standard deviation of the long cadence background.
    - bck_flag (bool): Flag indicating if the background calculation was successful.
    """
    # Calculating time of non-backsubbed peak (for 2s cadence instead of 1-min)
    bcktrunc_18=df_long.truncate(flare_time.start.iso,flare_time.end.iso)
    old_peak = bcktrunc_18.idxmax()

    # Background Time Calculation Function
    bck_startt, bck_endt, _, _, _, bck_flag = onsets.background(df_short, df_long, start_time = Time(fl_hek_start), peak_time = old_peak)
    srch_time = a.Time(Time(bck_startt, scale='utc').iso, Time(bck_endt,scale='utc').iso)
    trunc_054 = df_short.truncate(srch_time.start.iso,srch_time.end.iso)
    trunc_18 = df_long.truncate(srch_time.start.iso,srch_time.end.iso)
    bck_short = np.mean(trunc_054)
    bck_long = np.mean(trunc_18)

    bck_short_std = np.std(trunc_054)
    bck_long_std = np.std(trunc_18)

    # Background Subtracting the Data
    short_backsub = df_short - bck_short
    long_backsub = df_long - bck_long
    
    return short_backsub, long_backsub, bck_startt, bck_endt, bck_short_std, bck_long_std, bck_flag

def calc_peak(backsub_ts, flare_time):
    """
    Calculate the true peak time and flux for a given flare using the backsubbed long channel flux.

    Parameters:
    backsub_ts (TimeSeries): The backsubbed time series containing the long channel flux.
    flare_time (TimeRange): The time range of the flare.

    Returns:
    tuple: A tuple containing the true peak time and the corresponding flux value.

    """
    # Calculating true peak time (using backsubbed long channel flux)
    flare_long = backsub_ts.truncate(flare_time.start.iso, flare_time.end.iso)
    long_peakt = flare_long.idxmax()
    long_peakfl_18 = flare_long[long_peakt]
    return long_peakt, long_peakfl_18

def find_onset_start(short_backsub, long_backsub, bck_endt, true_peak, bck_short_std, bck_long_std, sigma = 4):
    """
    Finds the onset start time based on the given parameters.

    Parameters:
    - short_backsub: The short background subtracted signal.
    - long_backsub: The long background subtracted signal.
    - bck_endt: The end time of the background.
    - true_peak: The true peak time.
    - bck_short_std: The standard deviation of the short background.
    - bck_long_std: The standard deviation of the long background.
    - sigma: The threshold value in terms of standard deviations (default: 4).

    Returns:
    - onset_start: The onset start time.

    Raises:
    - ValueError: If no suitable onset start time is found.
    """
    srch_start = Time(bck_endt) - TimeDelta(60, format = 'sec')
    trange_search = a.Time(srch_start, true_peak)
    g_long_srch = long_backsub.truncate(trange_search.start.iso,trange_search.end.iso)
    g_short_srch = short_backsub.truncate(trange_search.start.iso,trange_search.end.iso)

    indexes = g_long_srch[(g_long_srch >= sigma*bck_long_std) & (g_short_srch >= sigma* bck_short_std)].index ## NOTE: 4 sigma threshold for both channels
    if len(indexes) == 0:
        raise ValueError
    else:
        first_index = indexes[0]

    if first_index > bck_endt:
        onset_start = first_index
    else:
        onset_start = bck_endt

    return onset_start

def calculate_temperature(short_backsub, long_backsub, index):
    """
    Calculate the temperature using the given short background subtraction, long background subtraction, and index.

    Parameters:
    short_backsub (list): The short background subtraction values.
    long_backsub (list): The long background subtraction values.
    index (int): The index of the temperature value to calculate.

    Returns:
    tuple: A tuple containing the calculated temperature (tmk) and the emission (em).
    """
    trunc_054_tem = short_backsub[index]
    trunc_18_tem = long_backsub[index]
    tmk, em = gtem.get_tem_old(trunc_18_tem, trunc_054_tem)
    return tmk, em

def set_nan_values(filenotfound_flag, index_flag_err, flare_data):
    """
    Concatenates a DataFrame with NaN values and the appropriate boolean flags (reason for nan values) to the given `flare_data` DataFrame.

    Parameters:
    - filenotfound_flag (bool): Flag indicating if a file was not found.
    - index_flag_err (bool): Flag indicating if there was an index error.
    - flare_data (pandas.DataFrame): DataFrame containing flare data.

    Returns:
    - pandas.DataFrame: Concatenated DataFrame with NaN values and implemented flags.

    """
    return pd.concat([flare_data, pd.DataFrame({'Peak Time Long': [np.nan], 
                                                'Peak Time Short':[np.nan],
                                                'Peak Flux': [np.nan], 
                                                'Background Start Time':[np.nan], 
                                                'Background End Time': [np.nan], 
                                                'Onset Start Time':[np.nan], 
                                                'Peak Temp Long': [np.nan], 
                                                'Peak EM Long': [np.nan],
                                                'Peak Temp Short': [np.nan], 
                                                'Peak EM Short': [np.nan], 
                                                'Flare Max Temp': [np.nan], 
                                                'Flare Max EM': [np.nan], 
                                                'Temp 1/8': [np.nan],
                                                'Temp 1/8 Upper': [np.nan], 
                                                'Temp 1/8 Lower': [np.nan],
                                                'EM 1/8': [np.nan], 
                                                'EM 1/8 Error': [np.nan],
                                                'Tdelta 1/8': [np.nan],
                                                'Temp 1/6': [np.nan], 
                                                'Temp 1/6 Upper': [np.nan], 
                                                'Temp 1/6 Lower': [np.nan],
                                                'EM 1/6': [np.nan],
                                                'EM 1/6 Error': [np.nan], 
                                                'Tdelta 1/6': [np.nan],
                                                'Temp 1/4': [np.nan], 
                                                'Temp 1/4 Upper': [np.nan], 
                                                'Temp 1/4 Lower': [np.nan],
                                                'EM 1/4': [np.nan], 
                                                'EM 1/4 Error': [np.nan], 
                                                'Tdelta 1/4': [np.nan],
                                                'Temp 1/3': [np.nan], 
                                                'Temp 1/3 Upper': [np.nan], 
                                                'Temp 1/3 Lower': [np.nan], 
                                                'EM 1/3': [np.nan],
                                                'EM 1/3 Error': [np.nan], 
                                                'Tdelta 1/3':[np.nan],
                                                'Temp 1/2': [np.nan], 
                                                'Temp 1/2 Upper': [np.nan], 
                                                'Temp 1/2 Lower': [np.nan], 
                                                'EM 1/2': [np.nan], 
                                                'EM 1/2 Error': [np.nan], 
                                                'Tdelta 1/2': [np.nan],
                                                'Background Flag': [np.nan], 
                                                'FileNotFound Flag':[filenotfound_flag], 
                                                'IndexError Flag': [index_flag_err]})], 
                     ignore_index=True)

def process_data(data_dir, hek_dir, output_dir, onset_length = 60, sigma = 4):
    """
    Process the data for GOES-15 flares. Saves data to a specified csv file.

    Args:
        data_dir (str): The directory path where the timeseries data files are located.
        hek_dir (str): The directory path where the HEK data file is located.
        output_dir (str): The directory path where the processed data will be saved.
        onset_length (int, optional): The length of the flare onset in minutes. Defaults to 60.
        sigma (int, optional): The sigma threshold for onset detection. Defaults to 4.

    Returns:
        None
    """
# Looping through all goes-15 flares in directory
    flares_df = pd.read_csv(hek_dir)
    flare_data = pd.DataFrame(columns = [])
    for index, row in flares_df.iterrows():
        print_fl_info(row['event_starttime'], row['fl_goescls'])
        print("The flare is of class: " + row['fl_goescls'])
        try:
            print("Loading Timeseries Data...")
        # Loading Timeseries Data from file in given data directory
            g15, flare_time = load_timeseries(data_dir, row['event_starttime'], row['event_endtime'])
        # Pre-processing the data - removing flagged data points
            print("Pre-processing Timeseries Data...")
            df_short, df_long = pre_process_timeseries(g15)
        # Background subtracting the data
            print("Background Subtracting Timeseries Data...")
            short_backsub, long_backsub, bck_startt, bck_endt, bck_short_std, bck_long_std, bck_flag = backsub_timeseries(df_short, df_long, row['event_starttime'], flare_time)
        # Calculating peak time and flux
            print("Calculating Peak Time and Flux...")
            true_peak, long_peakfl_18 = calc_peak(long_backsub, flare_time)
            short_peakt, _ = calc_peak(short_backsub, flare_time)
        # Calculating onset start time
            print("Calculating Onset Start Time...")
            onset_start = find_onset_start(short_backsub, long_backsub, bck_endt, true_peak, bck_short_std, bck_long_std, sigma = 4) 
            #NOTE: 4 sigma threshold for both channels, defaulted at 4 but can be altered.

        # Calculating fractional onset end times
            print("Calculating Fractional Onset End Times...")
            fractions = [8, 6, 4, 3, 2, 3/2, 4/3] # NOTE: These are the reciprocal fractions of the onset length
            endt_values = {}
            timedelta_values = {}
            tmk_values = {}
            em_values = {}
            tmk_err_upper = {}
            tmk_err_lower = {}
            em_err = {}

            for fraction in fractions:
                endt_values[fraction] = onsets.scaled_onset(onset_start, true_peak, fraction)
                timedelta_values[fraction] = pd.to_timedelta(endt_values[fraction] - onset_start).total_seconds()

                if onset_start < endt_values[fraction]:
                    tmk_values[fraction], em_values[fraction], tmk_err_upper[fraction], tmk_err_lower[fraction], em_err[fraction] = onsets.onset_tem(df_short, df_long, short_backsub, long_backsub, onset_start, endt_values[fraction], bck_short_std, bck_long_std)
                else:
                    tmk_values[fraction] = np.nan
                    em_values[fraction] = np.nan

            print("Calculated onset times!")
            
            try:
                # Calculating temperature at flare peak (long channel)
                tmk_peakl, em_peakl = calculate_temperature(short_backsub, long_backsub, true_peak)

                # Calculating temperature at flare peak (short_channel)
                tmk_peaks, em_peaks = calculate_temperature(short_backsub, long_backsub, short_peakt)

                # Calculating max flare temperature
                endt_half = endt_values[2]
                trunc_054_tem_max = short_backsub.truncate(endt_half, Time(true_peak).datetime)
                trunc_18_tem_max = long_backsub.truncate(endt_half, Time(true_peak).datetime)
                tmk_maxx, em_maxx = gtem.get_tem_old(trunc_18_tem_max, trunc_054_tem_max)
                
                tmk_max = np.nanmax(tmk_maxx)
                tmk_max_index = np.nanargmax(tmk_maxx) 
                em_max = em_maxx[tmk_max_index]
            except:
                tmk_max = np.nan
                em_max = np.nan

            print("SUCCESS!!!")
            
            filenotfound_flag = False
            index_flag_err = False

            flare_data = pd.concat([flare_data, pd.DataFrame({'Peak Time Long': [true_peak], 
                                                            'Peak Time Short':[short_peakt],
                                                            'Peak Flux': [long_peakfl_18], 
                                                            'Background Start Time':[bck_startt], 
                                                            'Background End Time': [bck_endt], 
                                                            'Onset Start Time':[onset_start], 
                                                            'Peak Temp Long': [tmk_peakl], 
                                                            'Peak EM Long': [em_peakl],
                                                            'Peak Temp Short': [tmk_peaks], 
                                                            'Peak EM Short': [em_peaks], 
                                                            'Flare Max Temp': [tmk_max], 
                                                            'Flare Max EM': [em_max], 
                                                            'Temp 1/8': [tmk_values[8]], 
                                                            'Temp 1/8 Upper': [tmk_err_upper[8]], 
                                                            'Temp 1/8 Lower': [tmk_err_lower[8]],
                                                            'EM 1/8': [em_values[8]], 
                                                            'EM 1/8 Error': [em_err[8]], 
                                                            'Tdelta 1/8': [timedelta_values[8]],
                                                            'Temp 1/6': [tmk_values[6]], 
                                                            'Temp 1/6 Upper': [tmk_err_upper[6]], 
                                                            'Temp 1/6 Lower': [tmk_err_lower[6]], 
                                                            'EM 1/6': [em_values[6]], 
                                                            'EM 1/6 Error': [em_err[6]], 
                                                            'Tdelta 1/6': [timedelta_values[6]],
                                                            'Temp 1/4': [tmk_values[4]], 
                                                            'Temp 1/4 Upper': [tmk_err_upper[4]], 
                                                            'Temp 1/4 Lower': [tmk_err_lower[4]],
                                                            'EM 1/4': [em_values[4]], 
                                                            'EM 1/4 Error': [em_err[4]], 
                                                            'Tdelta 1/4': [timedelta_values[4]],
                                                            'Temp 1/3': [tmk_values[3]], 
                                                            'Temp 1/3 Upper': [tmk_err_upper[3]], 
                                                            'Temp 1/3 Lower': [tmk_err_lower[3]], 
                                                            'EM 1/3': [em_values[3]], 
                                                            'EM 1/3 Error': [em_err[3]],
                                                            'Tdelta 1/3':[timedelta_values[3]],
                                                            'Temp 1/2': [tmk_values[2]], 
                                                            'Temp 1/2 Upper': [tmk_err_upper[2]], 
                                                            'Temp 1/2 Lower': [tmk_err_lower[2]], 
                                                            'EM 1/2': [em_values[2]], 
                                                            'EM 1/2 Error': [em_err[2]], 
                                                            'Tdelta 1/2': [timedelta_values[2]],
                                                            'Background Flag': [bck_flag], 
                                                            'FileNotFound Flag':[filenotfound_flag], 
                                                            'IndexError Flag': [index_flag_err]
                                                            })], ignore_index=True)

        except ValueError:
            filenotfound_flag = True
            index_flag_err = False
            print("NO FILE FOUND")
            print("ERROR: ValueError - NAN ROW PLACED")
            print("---------------------------------------")
            flare_data = set_nan_values(filenotfound_flag, index_flag_err, flare_data)
            continue

        except IndexError:
            filenotfound_flag = False
            index_flag_err = True
            print("ERROR: IndexError - NAN ROW PLACED")
            print("---------------------------------------")
            flare_data = set_nan_values(filenotfound_flag, index_flag_err, flare_data)
            continue
        except:
            filenotfound_flag = False
            index_flag_err = False
            print("ERROR: Unknown Error - NAN ROW PLACED")
            print("---------------------------------------")
            flare_data = set_nan_values(filenotfound_flag, index_flag_err, flare_data)
            continue


    if len(flare_data) == len(flares_df):
        full_data = pd.concat([flares_df, flare_data],axis = 1)
    else:
        print("ERROR, data lengths do not match")

    full_data.to_csv(output_dir)