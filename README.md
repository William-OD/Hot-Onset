# Hot-Onset
Code and analysis from my MSci Project into the hot onset of Solar Flares

## Downloading GOES XRS data
Can use 'Bulk_Downloads.ipynb' to use wget to bulk download GOES XRS data from https://www.ncei.noaa.gov/data/goes-space-environment-monitor/access/science/xrs/ to a specified data directory. There are, of course, other ways of downloading this flare data locally.

## Generating a Flare List of Relevant Flare Information.
Use 'Generating_HEK.ipynb' to generate a list of flares and information from the SSW Latest events and GOES databases through HEK.
This also performs some simple post-processing to the flare list too remove duplicates and set up a proximity flag for flares that start within 30 minutes of another flare ending.

## Calculation of Onset Times and Precise Calculation of Start and Peak times.
Use 'Full_run.ipynb' to consolidate downloaded raw data and HEK flare list to calculate the following information for each flare:
* Peak Time Long - The time in which the XRS-B channel data is at it's peak flux.
* Peak Time Short - The time in which the XRS-A channel data is at it's peak flux.
* Peak Flux - The maximum flux value of the long channel flux. This can be used for flare re-classification.
* Background Start Time - The start point of the one-minute background window, calculated by my new background-subtraction algorithm.
* Background End Time - The end point of the one-minute background window, calculated by my new background-subtraction algorithm.
* Onset Start Time - The start point of the onset, calculated by my onset-selection algorithm.
* Onset End Time - The start point of the onset, calculated by my onset-selection algorithm.
* Peak Temp Long - The peak temperature of the flare from the long channel flux.
* Peak Temp Short - The peak temperature of the flare from the short channel flux.
* Flare Max Temp - The Maxuimum temperature of the flare.
* Peak EM Long - The peak emission measure of the flare from the long channel flux.
* Peak EM short - The peak emission measure of the flare from the short channel flux.
* Flare Max EM - The maximum emission measure of the flare.
