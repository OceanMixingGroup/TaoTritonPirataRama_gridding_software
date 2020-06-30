% calculate_TIW_index
%
% use adcp and current meter velocity to calculate the TIW index
%
% steps:
% 1. bandpass meridional velocity between 12 and 33 days
% 2. calculate TIWKE as v^2
% 3. lowpass filter with a cutoff of 20 days



clear
addpath('./utilities/')

%% define location
% this creates correct strings and folders

%%%%%%%%% be sure to set the correct location in set_location.m! %%%%%%%%%%

set_location;

%% load common tao time

time = taotime;
time10min = time;

load([savedir 'euc_' latstrshort '_' lonstr '_hourly.mat'])
timehr = euchr.time;

%% load

% processed velocity data
load([savedir 'vel_' latstrshort '_' lonstr '_daily.mat'])
vel = veldy;
clear veldy

%% calculate average meridional velocity in upper 50 m

ind = vel.depth <= 50;
vgappy = nanmean(vel.v(ind,:),1);
v = fast_fillgap(vgappy,30);

%% band pass between 12 and 33 days

vlow12 = gappy_filt(12,'l1',1,v);
vlow33 = gappy_filt(33,'l1',1,v);
vband = vlow12-vlow33;

%% calculate KE

KE = vband.^2;

%% lowpass filter

tiwke = gappy_filt(20,'l1',1,KE);


%% put on other time grids

% daily
tiwdy.time = vel.time;
tiwdy.ke = tiwke;

% monthly
tiwmo = monthly_average(tiwdy);

% interp to hourly
tiwhr.time = timehr;
tiwhr.ke = interp1(tiwdy.time,tiwdy.ke,tiwhr.time);

% interp to 10 min
tiw.time = time10min;
tiw.ke = interp1(tiwdy.time,tiwdy.ke,tiw.time);


%% calculate annual cycle

for ii = 1:12
    clear ind
    ind = find(str2num(datestr(tiwmo.time,'mm')) == ii);
    tiwclm.month(ii) = ii;
    tiwclm.ke(ii) = nanmean(tiwmo.ke(ind));
end


%% save


tiw.readme = strvcat(['Tropical instability wave (TIW) kinetic energy (KE) at ' latstrshort ', ' lonstr],...
    'calculated as:',...
    '    1 averaged daily meridional velocity from ADCPs and current meters in upper 50 m',...
    '    2 bandpass filtered between 12 and 33 days',...
    '    3 KE assumed to be proportional to v^2 (i.e. KE = v_bandpass.^2)',...
    '    4 lowpass filtered with a 20 day cutoff',...
    ' ',...
    'time = matlab datenum',...
    'tiw.ke = TIW kinetic energy [m^2 s^-2]',...
    ' ',...
    ['created by ' pwd '/calculate_TIW_index.m'],...
    ['by Sally Warner on ' datestr(now,'mmm dd, yyyy')]);

tiwhr.readme = strvcat(tiw.readme,...
    ' ',...
    'daily data has been interpolated into hourly averages');

tiwmo.readme = strvcat(tiw.readme,...
    ' ',...
    'daily data has been averaged by month. A month nust have data from at',...
    'least 10 days to be averaged');

tiw.readme = strvcat(tiw.readme,...
    ' ',...
    'daily data has been interpolated into 10 min averages');

tiwclm.readme = strvcat(tiw.readme,...
    ' ',...
    'climatology is found by averaging monthly data.',...
    'i.e. All Januarys are averaged, all Februarys are averaged, etc.');

save([savedir 'tiw_' latstrshort '_' lonstr '_10min.mat'],'tiw')
save([savedir 'tiw_' latstrshort '_' lonstr '_hourly.mat'],'tiwhr')
save([savedir 'tiw_' latstrshort '_' lonstr '_daily.mat'],'tiwdy')
save([savedir 'tiw_' latstrshort '_' lonstr '_monthly.mat'],'tiwmo')
save([savedir 'tiw_' latstrshort '_' lonstr '_monthly_climatology.mat'],'tiwclm')


return 

%% plot

dd = dir([savedir 'tiw*']);

for ii = 1:length(dd)
    load([savedir dd(ii).name])
end

varsrm = fields(tiw);
vars = varsrm(2:end-1);


figure(265)
clf
    set(gcf,'position',[495          57         766        1048])


for ii = 1:length(vars)
    
    ax(ii) = subplot(length(vars),1,ii);
    
    plot(tiw.time,tiw.(vars{ii}))
    hold on
    plot(tiwhr.time,tiwhr.(vars{ii}))
    plot(tiwdy.time,tiwdy.(vars{ii}))
    plot(tiwmo.time,tiwmo.(vars{ii}))
    
    datetick
    
    ylabel(vars{ii})
    set(gca,'tickdir','out')
    
    if ii == 1
        title(['TIW KE at ' latstrshort ', ' lonstr])
        legend('10min','hr','daily','monthly','location','northwest')
        legend boxoff
    end
end

export_fig([figdir 'tiw.png'],'-r200')


