% main_driver.m
%
% This code runs all of the other codes that:
%   1. processes all auxiliary data
%   2. combine processed chipod data
%   3. saves 10 min, hourly, daily and monthly averages of all data
%
% written by Sally Warner


%%%%%%% do you want to update netcdf files of auxiliary data? %%%%%%%
disp('To download updated TAO netcdf files, in a terminal, go to:')
disp('~/ganges/data/TaoTritonPirataRama/ and run download.bash')
disp(' ')

clear

%%%%%%%%%%%%%% set location %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set_location;

%%%%%%%%%%%%%% set processing flags %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% process chiopds? 
% (If you're running the code to make mooring files for chipod processing,
% turn off chipod processing off until after chipods are processed.)
do_chi = 1;

% process chiopds? 
% (If you're running the code to make mooring files for chipod processing,
% turn off chipod processing off until after chipods are processed.)
do_plot = 1;

% process temperature, salinity and stratification?
do_T = 1;

% process velocity?
do_vel = 1;

% Override velocity flag for locations where no velocity data has ever been 
% collected. load_and_save_tao_adcp_cur will crash for these cases.
if lon == 235 % (No velocity measured at 125W)
    do_vel = 0;
end

% save flags so they aren't cleared by 'clear' statements in the
% individual codes
save ./proc_flags/flags.mat do_chi do_plot do_T do_vel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% (1) wind stress %%%%%%%
%%%%%%% no dependencies
disp('****************** processing wind stress...')
load_and_save_tao_winds
 
%%%%%%% (2) Jq0 %%%%%%%
%%%%%%% no dependencies
disp('****************** processing surface fluxes...')
load_and_save_tao_Jq0

%%%%%%% (3) SST %%%%%%%
%%%%%%% no dependencies
disp('****************** processing SST...')
load_and_save_tao_sst  

%%%%%%% (8) temperature and stratification %%%%%%%
%%%%%%% no dependencies
load ./proc_flags/flags.mat
if do_T
    disp('****************** processing temperature and stratification...')
    load_and_save_tao_stratification  
end

%%%%%%% (6) velocity and shear %%%%%%%
%%%%%%% uses winds
%%%%%%% uses mld (which is now create by load_and_save_tao_sst)
load ./proc_flags/flags.mat
if do_vel
    disp('****************** processing TAO velocity from ADCPs and current meters...')
    load_and_save_tao_adcp_cur
end

%%%%%%% (7) EUC depth %%%%%%%
%%%%%%% uses velocity
load ./proc_flags/flags.mat
if do_vel
    disp('****************** processing EUC depth...')
    calculate_eucdepth  
end

%%%%%%% (9) Ri %%%%%%%
%%%%%%% uses vel
%%%%%%% uses T
%%%%%%% uses EUC
load ./proc_flags/flags.mat
if do_vel & do_T
    disp('****************** calculating Richardson number...')
    calculate_Ri
end

%%%%%%% (10) TIW KE %%%%%%%
%%%%%%% uses vel
%%%%%%% uses EUC
load ./proc_flags/flags.mat
if do_vel
    disp('****************** calculating TIW KE...')
    calculate_TIW_index
end

%%%%%%% (11) chipods %%%%%%%
%%%%%%% no dependencies
load ./proc_flags/flags.mat
if do_chi
    disp('****************** loading and saving all chipods...')
    load_and_save_tao_chipods_old_and_new
end


%%%%%%% plot all data %%%%%%%
load ./proc_flags/flags.mat
if do_plot
    disp('****************** plotting all data...')
    plot_all_data
end



%%%%%%%%%%%%% optional sound to let you know when the code has finished
% return
load handel
sound(y,Fs)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% portions of the code that have been removed or replaced

%%%%%% (4) mixed layer depth %%%%%%%
% % % no dependencies (because it loads T from netcdf)
% % % disp('****************** calculating mixed layer depth...')
% % % calculate_mixed_layer_depth
%%% THIS IN NOW INCORPORATED INTO STRATIFICATION

%%%%%%% (5) enso index %%%%%%%
% % % % % uses SST
% % % % % uses ENSO_ONI.mat (which may need to be updated)
% % % % disp('****************** processing enso index...')
% % % % load_and_save_enso_oni_index
%%%% THIS CAN BE LEFT OUT FOR NOW   

