function [chi, euc, Jq0, mld, Ri, sst, t, vel, wind, enso, tiw] = load_all(avg)
%
% [chi, euc, Jq0, mld, Ri, sst, t, vel, wind, enso, tiw] = load_all(avg)
%
% load in all of the data for the tao_enso analysis
% you must specify an averaging period as a string:
%       '10min'
%       'hourly'
%       'daily'
%       'monthly'
   
    


% chipods
load(['../processed/chipods_0_140W_alldepths_' avg '.mat'])
if strcmp(avg,'monthly')
    chi = cmo;
elseif strcmp(avg,'daily')
    chi = cdy;
elseif strcmp(avg,'hourly')
    chi = chr;
elseif strcmp(avg,'10min')
    chi = c;
end

% euc depth and speed
load(['../processed/euc_0_140W_' avg '.mat'])
if strcmp(avg,'monthly')
    euc = eucmo;
elseif strcmp(avg,'daily')
    euc = eucdy;
elseif strcmp(avg,'hourly')
    euc = euchr;
end

% surface fluxes
load(['../processed/Jq0_0_140W_' avg '.mat'])
if strcmp(avg,'monthly')
    Jq0 = Jq0mo;
elseif strcmp(avg,'daily')
    Jq0 = Jq0dy;
elseif strcmp(avg,'hourly')
    Jq0 = Jq0hr;
end

% mixed layer depth
load(['../processed/mld_0_140W_' avg '.mat'])
if strcmp(avg,'monthly')
    mld = mldmo;
elseif strcmp(avg,'daily')
    mld = mlddy;
elseif strcmp(avg,'hourly')
    mld = mldhr;
end

% Ri
load(['../processed/Ri_0_140W_' avg '.mat'])
if strcmp(avg,'monthly')
    Ri = Rimo;
elseif strcmp(avg,'daily')
    Ri = Ridy;
elseif strcmp(avg,'hourly')
    Ri = Rihr;
end

% sst
load(['../processed/sst_0_140W_' avg '.mat'])
if strcmp(avg,'monthly')
    sst = sstmo;
elseif strcmp(avg,'daily')
    sst = sstdy;
elseif strcmp(avg,'hourly')
    sst = ssthr;
end

% temperature and stratification
load(['../processed/T_0_140W_' avg '.mat'])
if strcmp(avg,'monthly')
    t = tmo;
elseif strcmp(avg,'daily')
    t = tdy;
elseif strcmp(avg,'hourly')
    t = thr;
end

% velocity and shear
load(['../processed/vel_0_140W_' avg '.mat'])
if strcmp(avg,'monthly')
    vel = velmo;
elseif strcmp(avg,'daily')
    vel = veldy;
elseif strcmp(avg,'hourly')
    vel = velhr;
end

% winds
load(['../processed/wind_0_140W_' avg '.mat'])
if strcmp(avg,'monthly')
    wind = windmo;
elseif strcmp(avg,'daily')
    wind = winddy;
elseif strcmp(avg,'hourly')
    wind = windhr;
end

% ENSO index
load(['../processed/enso_oni_' avg '.mat'])
if strcmp(avg,'monthly')
    enso = ensomo;
elseif strcmp(avg,'daily')
    enso = ensody;
elseif strcmp(avg,'hourly')
    enso = ensohr;
end

% TIW kinetic energy
load(['../processed/TIW_0_140W_' avg '.mat'])
if strcmp(avg,'monthly')
    tiw = tiwmo;
elseif strcmp(avg,'daily')
    tiw = ensody;
elseif strcmp(avg,'hourly')
    tiw = tiwhr;
end


end