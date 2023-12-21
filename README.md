# Hot-Onset
Code and analysis from my MSci Project into the hot onset of Solar Flares. Information on the methodology of the algorithms, as well as results and analysis of data can be found in my Masters Thesis at: {Placeholder}
This project includes a new background-subtraction algorithm which is sensitive to low GOES XRS fluxes, as well as method for determining the hot onset phase of a solar flare.

## Getting Started
Note this guide is for setting up the tool from scratch. If you are not wanting to process and background subtract the raw files yourself, you can do so by skipping to step 5 and use the most recent `full_run` version in `legacy/Processed_Run_Data`.

1. Clone the repository to your local machine.
2. Install the required Python packages using a `pip` virtual environment: `pip install -r requirements.txt`
3. Download GOES data. You can use many methods such a `wget` to fetch this data from https://www.ncei.noaa.gov/data/goes-space-environment-monitor/access/ and store in the Data/ directory. As of the current revision the program will only work if the same directory structure is used e.g. `Data/goes15/year/month/*.nc`
4. To perform operations, run main.py by typing `python main.py`. This program will currently give three options:

    1. To fetch Heliophysics Event Knowledgebase Data (HEK) for all flares in the given data directory. Generates a list of flares and information from the SSW Latest events and GOES databases through HEK. This also performs some simple post-processing to the flare list to remove duplicates and set up a proximity flag for flares that start within 30 minutes of another flare ending.

    2. To use process HEK data and flare data to background subtract data and calculate the following information for each flare:
        * `Peak Time Long` - The time in which the XRS-B channel data is at it's peak flux.
        * `Peak Time Short` - The time in which the XRS-A channel data is at it's peak flux.
        * `Peak Flux` - The maximum flux value of the long channel flux. This can be used for flare re-classification.
        * `Background Start Time` - The start point of the one-minute background window, calculated by my new background-subtraction algorithm.
        * `Background End Time` - The end point of the one-minute background window, calculated by my new background-subtraction algorithm.
        * `Onset Start Time` - The start point of the onset, calculated by my onset-selection algorithm.
        * `Peak Temp Long` - The peak temperature of the flare from the long channel flux.
        * `Peak Temp Short` - The peak temperature of the flare from the short channel flux.
        * `Flare Max Temp` - The Maximum temperature of the flare.
        * `Peak EM Long` - The peak emission measure of the flare from the long channel flux.
        * `Peak EM short` - The peak emission measure of the flare from the short channel flux.
        * `Flare Max EM` - The maximum emission measure of the flare.

        In addition to these, data at given fractions (1/8, 1/6, 1/4, 1/3, 1/2) of time between the 'Onset Start Time' and the 'Peak Time Long' are labelled as:
        * `Temp {frac}` - The average temperature of the fractional onset period.
        * `Temp {frac} Upper` - The upper limit of the temperature uncertainty.
        * `Temp {frac} Lower` - The lower limit of the temperature uncertainty.
        * `EM {frac} Error` - The uncertainty in emission measure.
        * `Tdelta {frac}` - The time between the onset start and the fractional onset end.

        Additional Data Flags:
        * `Background Flag` - Set to true if the background subtraction algorithm failed, otherwise flase.
        * `FileNotFound Flag` - Set to true if the file(s) containing the data for the flare that is being searched don't exist. This is often caused by data gaps in the GOES-15 operations.
        * `IndexError Flag` - Set to true if there exists an IndexError at any point in the data run. This can happen for a variety of reasons, and exists to provide complete running of the batch.

        Original Analysis is saved to 'Full_Run_v{version number}.csv' in `Legacy/Processed_Run_Data`.

    3. To perform both 1 and 2 in the same operation.

5. Use the `Data_Analysis.ipynb` notebook in `Notebooks/` to analyse the onset data. Note that this section is in early development and needs updated for cleanliness. 


## Project Structure
`HEK_Data/`: Contains data such as `GOES15_HEK_Data_flag.csv`.
`Legacy/`: Contains older versions of the code and data.
`Notebooks/`: Contains Jupyter notebooks for data analysis and downloading data.
`Onset_Data/`: Contains Onset_Data.csv.
`Response_Functions/`: Contains response function files.
`Scripts/`: Contains Python scripts for data processing.

## Licences

## Statement of Work
Usage of all content in this repository was created and owned by  William O'Donnell (University of Glasgow), and is freely available for research purposes so long as creditation is given. This is with exception of:
* The 'Code/goesxrs_temp.py' file whose original author was Dr Iain Hannah (University of Glasgow), which was lightly edited to include uncertainty calculations and some data handling (original repository found at: https://github.com/ianan/xrs_example).

## TO DO
 0. Add in functionality to alter output filenames based on a label from user input. (e.g. "GOES15_HEK_Data.csv" -> "GOES15_HEK_Data_2010.csv")
 1. Add in functionality to read data from goes multi-directory format (year/month) and also a single directory format (just files).
 2. Check over general customization of the program. Any values that are non-fixed and could be experimented with.
 3. Add in functionality to not rely on first, last in the data directory when generating HEK file. Instead, create a list of all flare files and specifically search for those (consider flares spanning two files).
 4. Bug checking, assess input validation procedures and make robust.
 5. Add in the option to change the filename of the output HEK file and output onset file (but not the directory).
 6. Update the readme file (especially the setup instructions for downloading data and setting up the repo). Outline customization options.
 7. Add licence.
 8. Add in functionality for multi-satellite support. (GOES 3-14, 16-17).
