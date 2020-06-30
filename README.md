# TaoTritonPirataRama_gridding_software

This code combines chipod data and all ancillary data (temperature, salinity, velocity from both ADCPs and current meters, winds, surface fluxes, etc.) into .mat files on evenly spaced temporal and depth grids. Temporal resolution is at 10 minute, hourly, daily and monthly. Temperature, salinity and velocity are on 5m depth grids from 0 to 300m (0:5:300).

 

BASICS:

To run all of the code, use main_driver.m, which will call all of the other routines. Alternately, each routine can be run separately.

Before running any of the the code, be sure to set the location that you want to process in set_location.m

Also, be sure that the time grid encompasses all of the data that you want to process. As of now, it goes through June 1, 2020, but eventually that will need to be pushed back. When changing the time grid, ALL locations will need to be rerun so that they all remain on the same time grid. The common time grid is set in:
taotime.m




HOW TO RUN THE CODE:

First update the raw data by running:
- download.bash in ~/ganges/data/TaoTritonPirataRama/

The modules that need to be updated before processing a new location:
- set_location.m 
- taotime.m (may or may not need to be updated)

The code that will run all of the individual modules:
- main_driver.m

All of the modules that load and save raw data are:
- load_and_save_tao_stratification.m
- load_and_save_tao_winds.m
- load_and_save_enso_oni_index.m           
- load_and_save_tao_adcp_cur.m             
- load_and_save_tao_chipods_old_and_new.m  
- load_and_save_tao_Jq0.m                  
- load_and_save_tao_sst.m

All of the modules that do other calculations (such as calculating Ri) are:
- calculate_eucdepth.m                     
- calculate_Ri.m                           
- calculate_TIW_index.m                    

Plots can be made by running the last section of each individual module or by running:
- plot_all_data.m

Some of the codes need extra info and data for the specific location that you are processing. For instance, in order to grid chipods, an extra module telling the code where to find the chipod data will be needed:
./mfiles/chipod_info/chipod_info_at_0_XXXW.m

If you want to use ADCP or current meter files that aren’t in TaoTritonPirataRama, then you will need to set up an extra module in:
./mfiles/adcp_info/adcp_info_at_0_XXXW.m



TO CREATE MOORING FILES FOR CHIPOD PROCESSING:

To prepare mooring files for chipod processing, run the following scripts:
- set_location.m
- load_and_save_tao_stratification.m
- load_and_save_tao_adcp_cur.m             
- save_mooring_files_at_chipod_depths.m 
- plot_timeseries_at_chipod_depths.m

Then move the mooring files (both processed data and figures) to the input directory in the chipod_gust processing software.



UPDATING THE RAW DATA:

Most of the data that is used by these codes is saved in:
~/ganges/data/TaoTritonPirataRama/

Update the raw data by running:
- download.bash in ~/ganges/data/TaoTritonPirataRama/

