% calculate_Ri
%
% use the velocity to calculate the depth of the EUC. Hourly data.
% 
% The depth of the EUC core is assumed to be the depth of the maximum
% eastward velocity.

clear
addpath('./utilities/')

%% define location
% this creates correct strings and folders

%%%%%%%%% be sure to set the correct location in set_location.m! %%%%%%%%%%

set_location;

%% load common tao time

time = taotime;


%% load

disp('loading...')

% velocity and shear
load([savedir 'vel_' latstrshort '_' lonstr '_10min.mat'])

% temperature and stratification
load([savedir 'T_' latstrshort '_' lonstr '_10min.mat'])



%% calculate Ri

Ri.time = time;
Ri.depth = vel.depth;

% Ri = N^2/S^2
Ri.Ri = t.NTsq./vel.Ssq;

% Rir = S^2 - 4N^2;
Ri.Rir = vel.Ssq - 4*t.NTsq;


%% calculate % of time when Rir is > or < 0
% for the 10 minute data, the answers are either 1 or 0 because S^2 is
% either larger or smaller than 4N^2. As averages are taken, the percentage
% over each time period will emerge.
Ri.Rirpct = double(Ri.Rir >= 0);
Ri.Rirpct(isnan(Ri.Ri)) = NaN;




%% calculate hourly and daily

disp('calculating hourly, daily and monthly averages...')

% hourly average (dt = 10 min, so put 6 points into each hour)
Rihr = bin_average(Ri,6);
Rihr.depth = Ri.depth;

% daily average (dt = 10 min, so put 6*24=144 points into each average)
Ridy = bin_average(Ri,144);
Ridy.depth = Ri.depth;

% calculate monthly averages
Rimo = monthly_average(Ridy);
Rimo.depth = Ri.depth;

% %%%%%%%%%%% calculate medians %%%%%%%%%%%%%%%
% % hourly average (dt = 10 min, so put 6 points into each hour)
% Rimedianhr = bin_median(Ri,6);
% 
% % daily average (dt = 10 min, so put 6*24=144 points into each average)
% Rimediandy = bin_median(Ri,144);
% 
% % calculate monthly averages
% % note, I know this should be done with 10 min data, but it's just too time
% % consuming. The result is NOT the same as finding the median of the 10 min
% % data over a month as opposed to finding the median of the daily data over
% % each month. BUT the difference is not super huge whereas the time to run
% % the code IS huge.
% Rimedianmo = monthly_median(Rimediandy);
% 
% %%%%% put medians into other structures
% Ri.Rimedian = Ri.Ri;
% Ri.Rirmedian = Ri.Rir;
% 
% Rihr.Rimedian = Rimedianhr.Ri;
% Rihr.Rirmedian = Rimedianhr.Rir;
% 
% Ridy.Rimedian = Rimediandy.Ri;
% Ridy.Rirmedian = Rimediandy.Rir;
% 
% Rimo.Rimedian = Rimedianmo.Ri;
% Rimo.Rirmedian = Rimedianmo.Rir;


%% save

disp('saving...')

Ri.readme = strvcat(['Richardson number (Ri) at ' latstrshort ', ' lonstr],...
    'Calculated from 10 min time bins and 5 m depth bins of S^2 and N^2',...
    ' ',...
    'Ri         = Richardson number [ ],',...
    'Rir        = Reduced shear S^2-4N^2 [s^-2]',...
    'Rirpct     = percent of time when S^2 >= 4N^2 (note this works better',...
    '               for data in large time bins, for instance for daily',...
    '               averages, the percenage is of the 144 10-min points in a day',...
    'Rimedian   = median of Richardson numbers over a given averaging period',...
    '               (note, for 10 min data Rimedian = Ri)',...
    'Rirmedian  = median of reduced shear over given averaging period',...
    '',...
    ['created by ' pwd '/calculate_Ri.m'],...
    ['by Sally Warner on ' datestr(now,'mmm dd, yyyy')]);

Rihr.readme = strvcat(Ri.readme,...
    ' ',...
    '10 min data has been binned into hourly averages');

Ridy.readme = strvcat(Ri.readme,...
    ' ',...
    '10 min data has been binned into daily averages');

Rimo.readme = strvcat(Ri.readme,...
    ' ',...
    'daily data has been averaged by month. A month nust have data from at',...
    'least 10 days to be averaged');

save([savedir 'Ri_' latstrshort '_' lonstr '_10min.mat'],'Ri')
save([savedir 'Ri_' latstrshort '_' lonstr '_hourly.mat'],'Rihr')
save([savedir 'Ri_' latstrshort '_' lonstr '_daily.mat'],'Ridy')
save([savedir 'Ri_' latstrshort '_' lonstr '_monthly.mat'],'Rimo')

