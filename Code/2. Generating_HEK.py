import os 
import pandas as pd
import numpy as np
import re

from datetime import datetime
from astropy.time import Time
from astropy.time import TimeDelta

from sunpy.net import attrs as a
from sunpy.net import hek
from sunpy.net import Fido

def get_dates_from_files(data_dir):
    """
    Extracts dates from file names in a given directory.

    Parameters:
    data_dir (str): The directory to search for files.

    Returns:
    np.array: An array of dates extracted from file names.
    """
    # Initialize an empty list to store dates
    dates = []
    
    # Walk through the directory
    for root, dirs, files in os.walk(data_dir):
        # Iterate over each file
        for f in files:
            # Check if the file is a .nc file and contains the specified string
            if f.__contains__('.nc') and f.__contains__('sci_gxrs-l2-irrad_g15_d'):
                # Print the file name
                print(f)
                
                # Search for a date in the file name using regular expressions
                match = re.search(r'd(\d{8})', f)
                
                # Extract the date string from the match
                date_str = match.group(1)
                
                # Convert the date string to a datetime object
                date = datetime.strptime(date_str, '%Y%m%d')
                
                # Append the date to the list
                dates.append(date)
    
    # Convert the list of dates to a numpy array and return it
    return np.array(dates, dtype='datetime64[s]')

def search_and_save_flare_data(fl_search, event_type, condition, filename):
    """
    Searches for and saves flare data.

    Parameters:
    fl_search (Time): The time period to search for flares.
    event_type (str): The type of event to search for.
    condition (str): The condition to apply to the search.
    filename (str): The name of the file to save the results to.

    Returns:
    pd.DataFrame: A DataFrame containing the search results.
    """
    # Use Fido to search for flares in the specified time period that meet the condition
    res = Fido.search(fl_search, a.hek.EventType(event_type), condition)

    # Print a test message
    print("Test")

    # Extract the full results from the search
    fullres = res["hek"]

    # Extract the required columns from the results
    srch_res = fullres["event_starttime", "event_peaktime", "event_endtime", "fl_goescls", "ar_noaanum"]

    # Print the total number of flares in the period
    print(f"Total Number of flares in period: {len(srch_res)}\n")

    # Write the results to a CSV file
    srch_res.write(filename, overwrite = True, format="csv")

    # Return the results as a DataFrame
    return pd.read_csv(filename)

def process_data(data, origin):
    """
    Processes flare data by removing duplicates and adding an origin column.

    Parameters:
    data (pd.DataFrame): The data to process. It's expected to have columns 'event_peaktime' and 'fl_goescls'.
    origin (str): The origin of the data. This will be added as a new column 'Origin' in the DataFrame.

    Returns:
    pd.DataFrame: The processed data with duplicates removed and a new 'Origin' column.
    """
    # Remove duplicate rows based on 'event_peaktime' and 'fl_goescls' columns
    data = data.drop_duplicates(['event_peaktime', 'fl_goescls'])
    
    # Add a new column 'Origin' with the provided origin value
    data['Origin'] = origin
    
    # Return the processed DataFrame
    return data

def set_proximity_flag(data):
    """
    Sets a proximity flag for flare events that are close in time.

    Parameters:
    data (pd.DataFrame): The data to set the proximity flag for.

    Returns:
    pd.DataFrame: The data with the proximity flag set.
    """
    # Initialize the 'Proximity Flag' column with False
    data['Proximity Flag'] = False
    # Explicitly set the 'Proximity Flag' for the first row as False
    data.at[0, 'Proximity Flag'] = False

    # Iterate over the rows of the DataFrame, starting from the second row
    for i in range(1, len(data)):
        # Calculate the time difference between the start time of the current event and the end time of the previous event
        time_diff = pd.to_datetime(data.iloc[i]['event_starttime']) - pd.to_datetime(data.iloc[i-1]['event_endtime'])
        # Convert the time difference to minutes
        time_diff_minutes = time_diff.total_seconds() / 60

        # If the time difference is less than or equal to 30 minutes, set the 'Proximity Flag' for the current row as True
        if time_diff_minutes <= 30:
            data.at[i, 'Proximity Flag'] = True

    # Return the DataFrame with the 'Proximity Flag' set
    return data


def process_directory(input_dir, output_file):
    """
    Processes all files in a directory and saves the results to a .csv file.

    Parameters:
    input_dir (str): The directory to process files from.
    output_file (str): The name of the .csv file to save the results to.
    """
    # Get the dates from the file names in the input directory
    datetime_array = get_dates_from_files(input_dir)

    # Define the time period for the flare search
    first_date = min(datetime_array)
    last_date = Time(max(datetime_array)) + TimeDelta(1, format = 'jd') - TimeDelta(1, format ='sec')
    fl_search = a.Time(first_date, last_date)

    # Search for and save flare data from SSW Latest Events and GOES
    data_SSW = search_and_save_flare_data(fl_search, "FL", a.hek.FRM.Name == "SSW Latest Events", "GOES15_HEK_Data_SSW.csv")
    data_GOES = search_and_save_flare_data(fl_search, "FL", a.hek.OBS.Observatory == "GOES", "GOES15_HEK_Data_GOES.csv")

    # Process the data
    data_SSW = process_data(data_SSW, "SSW")
    data_GOES = process_data(data_GOES, "GOES")

    # Combine the data from SSW and GOES
    data_full = pd.concat([data_SSW, data_GOES],axis = 0)
    # Sort the data by event start time and origin, and remove duplicates
    data_full = data_full.sort_values(by = ['event_starttime', 'Origin'])
    data_full = data_full.drop_duplicates(['event_peaktime', 'fl_goescls'], keep = "last")
    data_full = data_full.sort_values(by = ['event_starttime'])
    data_full = data_full.reset_index(drop = True)
    # Save the full data to a CSV file
    data_full.to_csv(output_file)

    # Set the proximity flag for flare events that are close in time
    data = set_proximity_flag(data_full)

    # Remove duplicates based on event peak time and flare class, and reset the index
    data1 = data.drop_duplicates(subset = ['event_peaktime','fl_goescls'])
    data1 = data1.reset_index(drop=True)
    # Print the number of duplicates removed
    print(f"There were: {len(data) - len(data1)} duplicates removed")

    # Save the final data to a CSV file
    data1.to_csv(output_file)


## Generating the flare database
process_directory("D:\MastersProj\Data\Run_1", "Graphs/GOES15_HEK_testN.csv")