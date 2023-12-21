# Hot-Onset
Code and analysis from my MSci Project into the hot onset of Solar Flares. Information on the methodology of the algorithms, as well as results and analysis of data can be found in my Masters Thesis at: {Placeholder}
This project includes a new background-subtraction algorithm which is sensitive to low GOES XRS fluxes, as well as method for determining the hot onset phase of a solar flare.

All directory paths will have to be re-shuffled for use on your local machine. I have re-organsied this repo for readability, however the directory paths between each file may no longer match up. 

Usage of all content in this repository was created and owned by  William O'Donnell (University of Glasgow), and is freely available for research purposes so long as creditation is given. This is with exception of:
* The 'Code/goesxrs_temp.py' file whose original author was Dr Iain Hannah (University of Glasgow), which was lightly edited to include uncertainty calculations and some data handling (original repository found at: https://github.com/ianan/xrs_example).

## Downloading GOES XRS data
Can use 'Bulk_Downloads.ipynb' to use wget to bulk download GOES XRS data from https://www.ncei.noaa.gov/data/goes-space-environment-monitor/access/science/xrs/ to a specified data directory. There are, of course, other ways of downloading this flare data locally.

## Generating a Flare List of Relevant Flare Information.
Use 'Generating_HEK.ipynb' to generate a list of flares and information from the SSW Latest events and GOES databases through HEK.
This also performs some simple post-processing to the flare list too remove duplicates and set up a proximity flag for flares that start within 30 minutes of another flare ending.

## Calculation of Background, Onset Times and Precise Calculation of Start and Peak times.
Use 'Full_run.ipynb' to consolidate downloaded raw data and HEK flare list to calculate the following information for each flare:
* Peak Time Long - The time in which the XRS-B channel data is at it's peak flux.
* Peak Time Short - The time in which the XRS-A channel data is at it's peak flux.
* Peak Flux - The maximum flux value of the long channel flux. This can be used for flare re-classification.
* Background Start Time - The start point of the one-minute background window, calculated by my new background-subtraction algorithm.
* Background End Time - The end point of the one-minute background window, calculated by my new background-subtraction algorithm.
* Onset Start Time - The start point of the onset, calculated by my onset-selection algorithm.
* Peak Temp Long - The peak temperature of the flare from the long channel flux.
* Peak Temp Short - The peak temperature of the flare from the short channel flux.
* Flare Max Temp - The Maxuimum temperature of the flare.
* Peak EM Long - The peak emission measure of the flare from the long channel flux.
* Peak EM short - The peak emission measure of the flare from the short channel flux.
* Flare Max EM - The maximum emission measure of the flare.

In addition to these, data at given fractions (1/8, 1/6, 1/4, 1/3, 1/2) of time between the 'Onset Start Time' and the 'Peak Time Long' are labelled as:
* Temp {frac} - The average temperature of the fractional onset period.
* Temp {frac} Upper - The upper limit of the temperature uncertainty.
* Temp {frac} Lower - The lower limit of the temperature uncertainty.
* EM {frac} Error - The uncertainty in emission measure.
* Tdelta {frac} - The time between the onset start and the fractional onset end.

Additional Data Flags:
* Background Flag - Set to true if the background subtraction algorithm failed, otherwise flase.
* FileNotFound Flag - Set to true if the file(s) containing the data for the flare that is being searched don't exist. This is often caused by data gaps in the GOES-15 operations.
* IndexError Flag - Set to true if there exists an IndexError at any point in the data run. This can happen for a variety of reasons, and exists to provide complete running of the batch.

This data is saved to 'Full_Run_v{version number}.csv' in the Processed_Run_Data folder.

## Data Analysis
Analysis of these gathered data points is processed in 'Data_Analysis.ipynb'. 


Further documentation for this section is incomplete and will be updated in the near future. 


