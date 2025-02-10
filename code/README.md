# Code for emperor penguin data cleaning and analysis
### Author: Taylor Azizeh
### Date: 10 February 2025

These folders run through a few different steps of preparing data and running code.

**data_cleaning**: 1) Raw data files are first conatenate with depth data that has been ZOC using the IKNOS program. 2) The depth data are then interpolated from 1Hz to either 50Hz or 100Hz, which allows the files to be further decimated or downsampled without losing depth information. 3) The 100Hz acceleration data from 2019 are then decimated to 50Hz in MATLAB and exported.

**xx**: 
