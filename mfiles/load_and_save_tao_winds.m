% load_and_save_tao_winds
%
% load in wind netcdf file
% calc wind stresses
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


%% load high res data

%%%%% winds %%%%%

filename = ['../../TaoTritonPirataRama/high_resolution/10m/w' ...
    lower(latstr) lower(lonstr) '_10m.cdf'];
if exist(filename,'file')
    raw = loadnc(filename);
else
    disp('WARNING: load_and_save_tao_winds cannot find raw wind file!')
    disp(['looking for: ' filename])
    disp('Setting all wind (wind.u, wind.v, wind.speed, wind.dir) to NaN.')
    disp(' ')
    raw.time = time;
    raw.WU = NaN*time;
    raw.WV = NaN*time;
    raw.WS = NaN*time;
    raw.WD = NaN*time;
    raw.depth = -4;
end

% shorten winds to just 2005 onward
wind.time = time;
wind.u = interp1(raw.time,raw.WU,time);
wind.v = interp1(raw.time,raw.WV,time);
wind.speed = interp1(raw.time,raw.WS,time);
wind.dir = interp1(raw.time,raw.WD,time);
winddepth = raw.depth;


%%%%% air temp %%%%%

filename = ['../../TaoTritonPirataRama/high_resolution/10m/airt' ...
    lower(latstr) lower(lonstr) '_10m.cdf'];
if exist(filename,'file')
    raw = loadnc(filename);
else
    disp('WARNING: load_and_save_tao_winds cannot find raw air temp file!')
    disp(['looking for: ' filename])
    disp('Setting air temp to NaN, which will lead to NaN for wind stress.')
    disp(' ')
    raw.time = time;
    raw.AT = NaN*time;
end
airt.time = time;
airt.airt = interp1(raw.time,raw.AT,time);


%%%%% humidity %%%%%

filename = ['../../TaoTritonPirataRama/high_resolution/10m/rh' ...
    lower(latstr) lower(lonstr) '_10m.cdf'];
if exist(filename,'file')
    raw = loadnc(filename);
else
    disp('WARNING: load_and_save_tao_winds cannot find raw humidity file!')
    disp(['looking for: ' filename])
    disp('Setting humidity to NaN, which will lead to NaN for wind stress.')
    disp(' ')
    raw.time = time;
    raw.RH = NaN*time;
end
rh.time = time;
rh.rh = interp1(raw.time,raw.RH,time);



%% load in daily data (for more temporal coverage)

gooddailydata = 1;
  
%%%%% wind_u %%%%%
filename = ['../../TaoTritonPirataRama/' mooringarray '/uwnd_xyt_dy.cdf'];
if exist(filename,'file')
    rawudy = loadnc(filename);
else
    gooddailydata = 0;
end

%%%%% wind_v %%%%%
filename = ['../../TaoTritonPirataRama/' mooringarray '/vwnd_xyt_dy.cdf'];
if exist(filename,'file')
    rawvdy = loadnc(filename);
else
    gooddailydata = 0;
end

%%%%% wind_spd %%%%%
filename = ['../../TaoTritonPirataRama/' mooringarray '/wspd_xyt_dy.cdf'];
if exist(filename,'file')
    rawspddy = loadnc(filename);
else
    gooddailydata = 0;
end

%%%%% wind_dir %%%%%
filename = ['../../TaoTritonPirataRama/' mooringarray '/wdir_xyt_dy.cdf'];
if exist(filename,'file')
    rawdirdy = loadnc(filename);
else
    gooddailydata = 0;
end
    
%%%%% humidity %%%%%
filename = ['../../TaoTritonPirataRama/' mooringarray '/rh_xyt_dy.cdf'];
if exist(filename,'file')
    rawrhdy = loadnc(filename);
else
    gooddailydata = 0;
end
       
%%%%% air temp %%%%%
filename = ['../../TaoTritonPirataRama/' mooringarray '/airt_xyt_dy.cdf'];
if exist(filename,'file')
    rawairtdy = loadnc(filename);
else
    gooddailydata = 0;
end

% warning if at least one of the daily files could not be found
if gooddailydata == 0   
    disp('WARNING: load_and_save_tao_winds could not find one or more of the')
    disp(['daily wind files in ../../TaoTritonPirataRama/' mooringarray '/'])
    disp('     uwnd_xyt_dy.cdf, vwnd_xyt_dy.cdf, wspd_xyt_dy.cdf,')
    disp('     wdir_xyt_dy.cdf, rh_xyt_dy.cdf, airt_xyt_dy.cdf, etc.')
    disp('No daily data will be used to fill gaps in the high resolution data.')        
    disp(' ')
end
 
%%%%% extrapolate daily data at the correct lat and lon %%%%%
if gooddailydata
    indlat = find(rawudy.lat == lat);
    indlon = find(rawudy.lon == lon);
    
    if isempty(indlat) | isempty(indlon)
        gooddailydata = 0;
        disp(['WARNING: load_and_save_tao_winds cannot find the desired lat and lon'])
        disp(['in the daily ../../TaoTritonPirataRama/' mooringarray '/ mooring files: '])
        disp('     uwnd_xyt_dy.cdf, vwnd_xyt_dy.cdf, wspd_xyt_dy.cdf,')
        disp('     wdir_xyt_dy.cdf, rh_xyt_dy.cdf, airt_xyt_dy.cdf, etc.')
        disp('No daily data will be used to fill gaps in the high resolution data')
        disp(' ')
    else
        rawdy.time = rawudy.time;
        rawdy.lat = rawudy.lat(indlat);
        rawdy.lon = rawudy.lon(indlon);
        rawdy.u = squeeze(rawudy.WU(indlon,indlat,:));
        rawdy.v = squeeze(rawvdy.WV(indlon,indlat,:));
        rawdy.speed = squeeze(rawspddy.WS(indlon,indlat,:));
        rawdy.dir = squeeze(rawdirdy.WD(indlon,indlat,:));
        rawdy.rh = interp1(rawrhdy.time,squeeze(rawrhdy.RH(indlon,indlat,:)),rawdy.time);
        rawdy.airt = interp1(rawairtdy.time,squeeze(rawairtdy.AT(indlon,indlat,:)),rawdy.time);
    end
end



%% calculate wind stress for high res data

% usually winds are measured at 4m on TaoTritonPirataRama buoys
% calculate u10
u10 = sw_u10(wind.u,-winddepth);
v10 = sw_u10(wind.v,-winddepth);

% calculate Cd
U = sqrt(u10.^2 + v10.^2);
cd = sw_drag(U);

% calculate air density
rhoagappy = sw_airden(airt.airt,1013.25,rh.rh);
% Unfortunately, rhoagappy is super gappy, but varies by very little over a
% range of 1.16 and 1.19. Using a constant air density for rhoagappy won't 
% affect the final wind stress all that much (<3%). Filling gaps with mean.
rhoa = rhoagappy;
rhoa(isnan(rhoagappy)) = nanmean(rhoagappy);

% then calculate tau
wind.tau = cd.*rhoa.*U.^2;
wind.taux = cd.*rhoa.*abs(u10).*u10;
wind.tauy = cd.*rhoa.*abs(v10).*v10;



%% calculate wind stress for daily data

if gooddailydata
    
    % these winds are measured at 4m
    % first calculate u10
    u10dy = sw_u10(rawdy.u,4);
    v10dy = sw_u10(rawdy.v,4);

    % then calculate Cd
    Udy = sqrt(u10dy.^2 + v10dy.^2);
    cddy = sw_drag(Udy);

    % calculate air density
    rhoadygappy = sw_airden(rawdy.airt,1013.25,rawdy.rh);
    % similar to the high resolution data, we fill gaps with the mean value
    % since rhoagappy is so sparse.
    rhoady = rhoadygappy;
    rhoady(isnan(rhoadygappy)) = nanmean(rhoadygappy);

    % then calculate tau
    rawdy.tau = cddy.*rhoady.*Udy.^2;
    rawdy.taux = cddy.*rhoady.*abs(u10dy).*u10dy;
    rawdy.tauy = cddy.*rhoady.*abs(v10dy).*v10dy;
end

%% calculate hourly, daily, and monthly averages

% hourly average (wind dt = 10 min, so put 6 points into each hour)
windhr = bin_average(wind,6);

% daily average (wind dt = 10 min, so put 6*24=144 points into each average)
winddy = bin_average(wind,144);

% splice in daily data into high resolution data
if gooddailydata
    
    disp('Filling gaps in high resolution data with daily data')

    rawdyinterp.time = winddy.time;
    rawdyinterp.u = interp1(rawdy.time,rawdy.u,winddy.time);
    rawdyinterp.v = interp1(rawdy.time,rawdy.v,winddy.time);
    rawdyinterp.speed = interp1(rawdy.time,rawdy.speed,winddy.time);
    rawdyinterp.dir = interp1(rawdy.time,rawdy.dir,winddy.time);
    rawdyinterp.taux = interp1(rawdy.time,rawdy.taux,winddy.time);
    rawdyinterp.tauy = interp1(rawdy.time,rawdy.tauy,winddy.time);
    rawdyinterp.tau = interp1(rawdy.time,rawdy.tau,winddy.time);

    indnan = isnan(winddy.u);
    winddy.u(indnan) = rawdyinterp.u(indnan);
    winddy.v(indnan) = rawdyinterp.v(indnan);
    winddy.speed(indnan) = rawdyinterp.speed(indnan);
    winddy.dir(indnan) = rawdyinterp.dir(indnan);
    indnan = isnan(winddy.tau);
    winddy.tau(indnan) = rawdyinterp.tau(indnan);
    winddy.taux(indnan) = rawdyinterp.taux(indnan);
    winddy.tauy(indnan) = rawdyinterp.tauy(indnan);
end

% monthly averages
windmo = monthly_average(winddy);


%% save

disp('saving...')

wind.readme = strvcat(['Wind stress at ' latstrshort ', ' lonstr ' from ' mooringarray ' array'],...
    ['wind sensor located at ' num2str(-winddepth) 'm above sea level'],...
    'u and v [m/s]',...
    'note that u and v are in oceanographic convention with u being wind to ',...
    '   the east and v being wind to the north',...
    'tau: wind stress = cd * rhoair * (u10.^2 + v10.^2) [N/m^2]',...
    'taux = zonal compondent (u) of windstress (cd*rhoair*abs(u10)*u10)',...
    'tauy = meridional compondent (v) of windstress (cd*rhoair*abs(v10)*v10)',...
    'cd, rhoair, u10 and v10 all calculated with seawater toolbox',...
    ' ',...
    'Note, there are many holes in the high resolution data, especially',...
    'during 2007-2011. This data is NOT spliced into the 10 min or hourly',...
    'averages, but it IS added to the daily and monthly averaged data.',...
    ' ',...
    ['created by ' pwd '/load_and_save_tao_winds.m'],...
    ['by Sally Warner on ' datestr(now,'mmm dd, yyyy')]);

windhr.readme = strvcat(wind.readme,...
    ' ',...
    '10 min data has been binned into hourly averages');

winddy.readme = strvcat(wind.readme,...
    ' ',...
    '10 min data has been binned into daily averages');

windmo.readme = strvcat(wind.readme,...
    ' ',...
    'daily data has been averaged by month. A month nust have data from at',...
    'least 10 days to be averaged');

save([savedir 'wind_' latstrshort '_' lonstr '_10min.mat'],'wind')
save([savedir 'wind_' latstrshort '_' lonstr '_hourly.mat'],'windhr')
save([savedir 'wind_' latstrshort '_' lonstr '_daily.mat'],'winddy')
save([savedir 'wind_' latstrshort '_' lonstr '_monthly.mat'],'windmo')





return

%% make plots


dd = dir([savedir 'wind*']);

for ii = 1:length(dd)
    load([savedir dd(ii).name])
end

vars = {'u';'v';'speed';'dir';'tau';'taux';'tauy'};


figure(265)
clf
    set(gcf,'position',[495          57         766        1048])


for ii = 1:length(vars)
    
    ax(ii) = subplot(length(vars),1,ii);
    
    plot(wind.time,wind.(vars{ii}))
    hold on
    plot(windhr.time,windhr.(vars{ii}))
    plot(winddy.time,winddy.(vars{ii}))
    plot(windmo.time,windmo.(vars{ii}))
    
    datetick
    
    ylabel(vars{ii})
    set(gca,'tickdir','out')
    
    if ii == 1
        title(['Winds at ' latstrshort ', ' lonstr])
        legend('10min','hr','daily','monthly','location','northwest')
        legend boxoff
    end
end

export_fig([figdir 'wind.png'],'-r200')


