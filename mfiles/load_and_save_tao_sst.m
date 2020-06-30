% load_and_save_tao_sst
%
% load in sst netcdf file
% calc hourly, daily, monthly averages
% save
% plot if desired
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


%% load in high res data

%%%%% sst %%%%%
filename = ['../../TaoTritonPirataRama/high_resolution/10m/sst' ...
    lower(latstr) lower(lonstr) '_10m.cdf'];
if exist(filename,'file')
    raw = loadnc(filename);
else
    disp('WARNING: load_and_save_tao_sst cannot find raw sst file!')
    disp(['looking for: ' filename])
    disp('Setting all sst data (sst.T) to NaN.')
    disp(' ')
    raw.time = time;
    raw.T = NaN*time;
end

% shorten to just 2005 onward
sst.time = time;
sst.T = interp1(raw.time,raw.T,time);



%% load in daily data (for more temporal coverage)

gooddailydata = 1;

%%%%% sst %%%%%
filename = ['../../TaoTritonPirataRama/' mooringarray '/sst_xyt_dy.cdf'];
if exist(filename,'file')
    rawsstdy = loadnc(filename);
else
    gooddailydata = 0;
    disp('WARNING: load_and_save_tao_sst could not find the')
    disp(['daily sst file: ../../TaoTritonPirataRama/' mooringarray '/sst_xyt_dy.cdf'])
    disp('No daily data will be used to fill gaps in the high resolution data.')        
    disp(' ')
end

%%%%% extrapolate daily data at the correct lat and lon %%%%%
if gooddailydata
    indlat = find(rawsstdy.lat == lat);
    indlon = find(rawsstdy.lon == lon);
    
    if isempty(indlat) | isempty(indlon)
        gooddailydata = 0;
        disp(['WARNING: load_and_save_tao_sst cannot find the desired lat and lon'])
        disp(['in the daily data: ../../TaoTritonPirataRama/' mooringarray '/sst_xyt_dy.cdf'])
        disp('No daily data will be used to fill gaps in the high resolution data')
        disp(' ')
    else
        rawdy.time = rawsstdy.time;
        rawdy.lat = rawsstdy.lat(indlat);
        rawdy.lon = rawsstdy.lon(indlon);
        rawdy.T = squeeze(rawsstdy.T(indlon,indlat,:));
    end
end



%% calculate hourly, daily, and monthly averages

% hourly average (wind dt = 10 min, so put 6 points into each hour)
ssthr = bin_average(sst,6);

% daily average (wind dt = 10 min, so put 6*24=144 points into each average)
sstdy = bin_average(sst,144);


% splice in daily data into high resolution data
if gooddailydata
    
    disp('Filling gaps in high resolution data with daily data')

    rawdyinterp.time = sstdy.time;
    rawdyinterp.T = interp1(rawdy.time,rawdy.T,sstdy.time);

    indnan = isnan(sstdy.T);
    sstdy.T(indnan) = rawdyinterp.T(indnan);
end

% monthly averages
sstmo = monthly_average(sstdy);


%% calculate dSST/dt to double check averaging
% important to do this AFTER averaging temperature. It's a lot less noisy
% and more physically correct to take the gradient of the temperatures at
% a given time step rather than the average of the 10min gradients.

sst.dTdt = gradient(sst.T,1/(24*6));            % °C/day
ssthr.dTdt = gradient(ssthr.T,1/24);       % °C/day
sstdy.dTdt = gradient(sstdy.T,1);          % °C/day
sstmo.dTdt = gradient(sstmo.T,1);            % °C/month

% sstmo.dTdt = sstmo.T*NaN;          
% sstmo.dTdt(1:end-1) = sstmo.T(2:end) - sstmo.T(1:end-1);     % °C/month


%% save


disp('saving...')

sst.readme = strvcat(['Sea surface temperature (SST) at ' latstrshort ', ' lonstr ' from ' mooringarray ' array'],...
    'T [°C]',...
    'dTdt [°C/day] = gradient of 10 min temperature',...
    ' ',...
    ['created by ' pwd '/load_and_save_tao_sst.m'],...
    ['by Sally Warner on ' datestr(now,'mmm dd, yyyy')]);

ssthr.readme = strvcat(['Sea surface temperature (SST) at ' latstrshort ', ' lonstr ' from ' mooringarray ' array'],...
    'T [°C] (10 min data has been binned into hourly averages)',...
    'dTdt [°C/day] = gradient of hourly temperature',...
    ' ',...
    ['created by ' pwd '/load_and_save_tao_sst.m'],...
    ['by Sally Warner on ' datestr(now,'mmm dd, yyyy')]);

sstdy.readme = strvcat(['Sea surface temperature (SST) at ' latstrshort ', ' lonstr ' from ' mooringarray ' array'],...
    'T [°C] (10 min data has been binned into daily averages)',...
    'dTdt [°C/day] = gradient of daily temperature',...
    ' ',...
    ['created by ' pwd '/load_and_save_tao_sst.m'],...
    ['by Sally Warner on ' datestr(now,'mmm dd, yyyy')]);

sstmo.readme = strvcat(['Sea surface temperature (SST) at ' latstrshort ', ' lonstr ' from ' mooringarray ' array'],...
    'T [°C] (daily data has been averaged by month. A month nust have data from at',...
    '    least 10 days to be averaged)',...
    'dTdt [°C/month] = gradient of monthly temperature',...
    '    (divide by 30 to get units of °C/day)',...
    ' ',...
    ['created by ' pwd '/load_and_save_tao_sst.m'],...
    ['by Sally Warner on ' datestr(now,'mmm dd, yyyy')]);

save([savedir 'sst_' latstrshort '_' lonstr '_10min.mat'],'sst')
save([savedir 'sst_' latstrshort '_' lonstr '_hourly.mat'],'ssthr')
save([savedir 'sst_' latstrshort '_' lonstr '_daily.mat'],'sstdy')
save([savedir 'sst_' latstrshort '_' lonstr '_monthly.mat'],'sstmo')


return

%% plot

dd = dir([savedir 'sst*']);

for ii = 1:length(dd)
    load([savedir dd(ii).name])
end

varsrm = fields(sst);
vars = varsrm(2:end-1);


figure(265)
clf
    set(gcf,'position',[495          57         766        1048])


for ii = 1:length(vars)
    
    ax(ii) = subplot(length(vars),1,ii);
    
    plot(sst.time,sst.(vars{ii}))
    hold on
    plot(ssthr.time,ssthr.(vars{ii}))
    plot(sstdy.time,sstdy.(vars{ii}))
    plot(sstmo.time,sstmo.(vars{ii}))
    
    datetick
    
    ylabel(vars{ii})
    set(gca,'tickdir','out')
    
    if ii == 1
        title(['SST at ' latstrshort ', ' lonstr])
        legend('10min','hr','daily','monthly','location','northwest')
        legend boxoff
    end
end

export_fig([figdir 'sst.png'],'-r200')

