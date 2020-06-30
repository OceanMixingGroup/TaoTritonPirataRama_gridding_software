% load_and_save_tao_Jq0
%
% load in radiation netcdf files
% calc hourly, daily, monthly averages
% save
%
% written by Sally Warner


clear
addpath('./utilities/')

%% define location
% this creates correct strings and folders

%%%%%%%%% be sure to set the correct location in set_location.m! %%%%%%%%%%

set_location;

%% load common tao time

time = taotime;


%% load TAO surface fluxes


%%%%% shortwave radiation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename = ['../../TaoTritonPirataRama/high_resolution/2m/rad' ...
    lower(latstr) lower(lonstr) '_2m.cdf'];
if exist(filename,'file')
    raw = loadnc(filename);
else
    disp('WARNING: load_and_save_tao_Jq0 cannot find shortwave radiation file!')
    disp(['looking for: ~/ganges/data/TaoTritonPirataRama/high_resolution/2m/rad' lower(latstr) lower(lonstr) '_2m.cdf'])
    disp('Setting all shortwave radiation (Jq0.swr) to NaN.')
    disp(' ')
    raw.time = time;
    raw.RD = NaN*time;
end

% shorten SWR to just 2005 onward and interpolate to 10 minute time grid
Jq0.time = time;
Jq0.swr = interp1(raw.time,raw.RD,time);




%%%%% longwave radiation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note: There are two versions of longwave radiation:
%       1. 2m/lw0n140w_2m.cdf: Referred to as "LONGWAVE RADIATION" in the
%          attributes of the netcdf file. This is just the upward portion
%          of longwave radiation. Must subtract the downward portion to get
%          the NET longwave radiation which is what we want here. Has an 
%          average value around 400 W/m^2. DO NOT USE here!
%       2. hr/lwnetnc0n140w_hr.cdf: Referred to as "NET LONGWAVE RADIATION"
%          in the attributes of the netcdf file. Average value around 
%          50 W/m^2. This is the difference between the upward (~420 W/m^2)
%          and downward (~370 W/m^2) longwave radiations. Use this version
%          of longwave radiation in this code!

filename = ['../../TaoTritonPirataRama/high_resolution/hr/lwnetnc' ...
    lower(latstr) lower(lonstr) '_hr.cdf'];
if exist(filename,'file')
    raw = loadnc(filename);
else
    disp('WARNING: load_and_save_tao_Jq0 cannot find longwave radiation file!')
    disp(['looking for: ~/ganges/data/TaoTritonPirataRama/high_resolution/hr/lwnetnc' lower(latstr) lower(lonstr) '_hr.cdf'])
    disp('Setting all longwave radiation (Jq0.lwr) to NaN.')
    if lon == 350 | lon == 67
        disp('(Note: This is a known issue at 0, 10W and at 0, 67E where longwave')
        disp('radiation appears not to have ever been measured.)')
    end
    disp(' ')
    raw.time = time;
    raw.LWN = NaN*time;
end

Jq0.lwr = interp1(raw.time,raw.LWN,time);




%%%%% latent heat flux %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note: There are four versions of latent heat flux (on hourly time scale): 
%         1. qlat0n140w_hr.cdf: Fluxes computed from relative wind speeds 
%            (very sparce; starts in 2006)
%         2. qlat_nclw0n140w_hr.cdf: Fluxes computed using absolute wind speeds 
%            and longwave climatology (good coverage; starts in 1998)
%         3. qlat_nlw0n140w_hr.cdf: Fluxes computed using relative wind speeds 
%            and longwave climatology (very sparce; starts in 1998)
%         4. qlatnc0n140w_hr.cdf: Fluxes computed from absolute wind speeds 
%            (good coverage; starts in 2006)
% ** The differences are essentially negligable ±5 W/m^2 (most noticable in
% the relative vs. absolute wind speeds). So using version 2 which has the
% best coverage.

% filename = ['../../TaoTritonPirataRama/high_resolution/hr/qlat_nclw' ...
%     lower(latstr) lower(lonstr) '_hr.cdf'];
filename = ['../../TaoTritonPirataRama/high_resolution/hr/qlatnc' ...
    lower(latstr) lower(lonstr) '_hr.cdf'];
if exist(filename,'file')
    raw = loadnc(filename);
else
    disp('WARNING: load_and_save_tao_Jq0 cannot find latent heat flux file!')
    disp(['looking for: ~/ganges/data/TaoTritonPirataRama/high_resolution/hr/qlat_nclw' lower(latstr) lower(lonstr) '_hr.cdf'])
    disp('Setting all latent heat flux (Jq0.lhf) to NaN.')
    raw.time = time;
    raw.QL = NaN*time;
end

Jq0.lhf = interp1(raw.time,raw.QL,time);




%%%%% sensible heat flux %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: There are four versions of sensible heat flux (on hourly time scale): 
%         1. qsen0n140w_hr.cdf: Fluxes computed from relative wind speeds 
%            (very sparce; starts in 2006)
%         2. qsen_nclw0n140w_hr.cdf: Fluxes computed using absolute wind speeds 
%            and longwave climatology (good coverage; starts in 1998)
%         3. qsen_nlw0n140w_hr.cdf: Fluxes computed using relative wind speeds 
%            and longwave climatology (very sparce; starts in 1998)
%         4. qsennc0n140w_hr.cdf: Fluxes computed from absolute wind speeds 
%            (good coverage; starts in 2006)
% ** The differences are essentially negligable ±2 W/m^2. So I'm using 
% version 2 which has the best coverage.

% filename = ['../../TaoTritonPirataRama/high_resolution/hr/qsen_nclw' ...
%     lower(latstr) lower(lonstr) '_hr.cdf'];
filename = ['../../TaoTritonPirataRama/high_resolution/hr/qsennc' ...
    lower(latstr) lower(lonstr) '_hr.cdf'];
if exist(filename,'file')
    raw = loadnc(filename);
else
    disp('WARNING: load_and_save_tao_Jq0 cannot find sensible heat flux file!')
    disp(['looking for: ~/ganges/data/TaoTritonPirataRama/high_resolution/hr/qsen_nclw' lower(latstr) lower(lonstr) '_hr.cdf'])
    disp('Setting all sensible heat flux (Jq0.shf) to NaN.')
    raw.time = time;
    raw.QS = NaN*time;
end

Jq0.shf = interp1(raw.time,raw.QS,time);




%%%%% net heat flux %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note: There are two versions of net heat flux (on hourly time scale): 
%         1. qnet0n140w_hr.cdf: Fluxes computed from relative wind speeds 
%            (very sparce; starts in 2006)
%         2. qnetnc0n140w_hr.cdf: Fluxes computed from absolute wind speeds 
%            (good coverage; starts in 2006)
% I'm using version 2 which has the best coverage.
% *** The problem is that I don't have LWR before 2006 and I don't know
% what the climatology is... So the pre-calculated net flux only starts in 2006
filename = ['../../TaoTritonPirataRama/high_resolution/hr/qnetnc' ...
    lower(latstr) lower(lonstr) '_hr.cdf'];
if exist(filename,'file')
    raw = loadnc(filename);
else
    disp('WARNING: load_and_save_tao_Jq0 cannot find the net heat flux file!')
    disp(['looking for: ~/ganges/data/TaoTritonPirataRama/high_resolution/hr/qnetnc' lower(latstr) lower(lonstr) '_hr.cdf'])
    disp('Setting all net heat flux (Jq0.net) to NaN.')
    if lon == 350
        disp('(Note: This is a known issue at 0°, 10°W, where there is no file')
        disp('for net heat flux, likely because there is no longwave radiation.)')
    end
    raw.time = time;
    raw.QT = NaN*time;
end

Jq0.net = interp1(raw.time,raw.QT,time);







%% calculate hourly, daily, and monthly averages

% hourly average (Jq0 dt = 10 min, so put 6 points into each hour)
% (this should be equivalent to the original data before interpolation)
Jq0hr = bin_average(Jq0,6);

% daily average (Jq0 dt = 10 min, so put 6*24=144 points into each average)
Jq0dy = bin_average(Jq0,144);

% monthly averages
Jq0mo = monthly_average(Jq0dy);

%% save

disp('saving...')

Jq0.readme = strvcat(['Surface fluxes at ' latstr ', ' lonstr ' from TAO array'],...
    ' ',...
    'swr = downward shortwave radiation (positive = heat INTO the ocean)',...
    'lwr = longwave radiation (positive = heat OUT of the ocean)',...
    'lhf = latent heat flux (positive = heat OUT of the ocean)',...
    'shf = sensible heat flux (positive = heat OUT of the ocean)',...
    'net = net heat flux (positive = heat INTO the ocean)',...
    '      net = swr - lwr - lhf - shf',...
    'all variables in [W m^-2]',...
    ' ',...
    ['created by ' pwd '/load_and_save_tao_Jq0.m'],...
    ['by Sally Warner on ' datestr(now,'mmm dd, yyyy')],...
    ' ',...
    'note: TAO hourly data is used for these calculations and interpolated',...
    'to 10 minute intervals to match chipods. Shortwave and longwave are',...
    'measured at 2 min intervals and met data is saved in 10 min intervals',...
    'so if higher resolution data is needed, it can be calculated from',...
    'raw tao data.');

Jq0hr.readme = strvcat(Jq0.readme,...
    ' ',...
    '10 min data has been binned into hourly averages. Matches original hourly data.');

Jq0dy.readme = strvcat(Jq0.readme,...
    ' ',...
    '10 min data has been binned into daily averages');

Jq0mo.readme = strvcat(Jq0.readme,...
    ' ',...
    'daily data has been averaged by month. A month nust have data from at',...
    'least 10 days to be averaged');

save([savedir 'Jq0_' latstrshort '_' lonstr '_10min.mat'],'Jq0')
save([savedir 'Jq0_' latstrshort '_' lonstr '_hourly.mat'],'Jq0hr')
save([savedir 'Jq0_' latstrshort '_' lonstr '_daily.mat'],'Jq0dy')
save([savedir 'Jq0_' latstrshort '_' lonstr '_monthly.mat'],'Jq0mo')


return
%% plot

dd = dir([savedir 'Jq0*']);

for ii = 1:length(dd)
    load([savedir dd(ii).name])
end

varsrm = fields(Jq0);
vars = varsrm(2:end-1);


figure(265)
clf
    set(gcf,'position',[495          57         766        1048])


for ii = 1:length(vars)
    
    ax(ii) = subplot(length(vars),1,ii);
    
    plot(Jq0.time,Jq0.(vars{ii}))
    hold on
    plot(Jq0hr.time,Jq0hr.(vars{ii}))
    plot(Jq0dy.time,Jq0dy.(vars{ii}))
    plot(Jq0mo.time,Jq0mo.(vars{ii}))
    
    datetick
    
    ylabel(vars{ii})
    set(gca,'tickdir','out')
    
    if ii == 1
        title(['Jq0 at ' latstrshort ', ' lonstr])
        legend('10min','hr','daily','monthly','location','northwest')
        legend boxoff
    end
end

export_fig([figdir 'Jq0.png'],'-r200')

