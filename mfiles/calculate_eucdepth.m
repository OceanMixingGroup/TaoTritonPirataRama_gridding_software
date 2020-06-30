% calculate_eucdepth
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

% processed ADCP data
load([savedir 'vel_' latstrshort '_' lonstr '_10min.mat'])


%% calculate depth of EUC
% the depth of the core of the EUC is assumed to be the depth of the
% maximum zonal velocity

euc.time = time;

% find depth of maximum velocity
[euc.speed,euc.ind] = nanmax(vel.u,[],1);
euc.depth = vel.depth(euc.ind)';

% nan out values that are too shallow
indbad = euc.depth < 40;
euc.depth(indbad) = NaN;
euc.speed(indbad) = NaN;
euc.ind(indbad) = NaN;


%% calculate hourly and daily

% hourly average (dt = 10 min, so put 6 points into each hour)
euchr = bin_average(euc,6);

% daily average (dt = 10 min, so put 6*24=144 points into each average)
eucdy = bin_average(euc,144);


%% calcualte the vertical velocity at the depth of the EUC
% use the daily averaged data

eucdy.w = gradient(-eucdy.depth,1); %[m/day]

% interpolate to faster time scales
euc.w = interp1(eucdy.time,eucdy.w,euc.time);
euchr.w = interp1(eucdy.time,eucdy.w,euchr.time);


%% calculate monthly averages
% do this after calculating vertical velocity from daily averaged depths

eucmo = monthly_average(eucdy);


%% calculate annual cycle

for ii = 1:12
    clear ind
    ind = find(str2num(datestr(eucmo.time,'mm')) == ii);
    eucclm.month(ii) = ii;
    eucclm.depth(ii) = nanmean(eucmo.depth(ind));
    eucclm.speed(ii) = nanmean(eucmo.speed(ind));
%     eucclm.w(ii) = nanmean(eucmo.w(ind));
end
eucclm.w = (-eucclm.depth([2:12 1]) - -eucclm.depth([12 1:11]))/30;


%% save

euc.readme = strvcat(['Equatorial Undercurrent (EUC) depth at ' latstrshort ', ' lonstr],...
    'Calculated as the depth of maximum zonal velocity as measured by ADCPs',...
    'depth      = depth of the EUC [m],',...
    'speed      = zonal velocity at the depth of the euc core [m/s]',...
    'w          = the vertical velocity at the depth of the EUC [m/day]',...
    '             (time derivative of the daily EUC depth)',...
    'ind        = index of the EUC depth corresponding to vel.depth',...
    '',...
    ['created by ' pwd '/calculate_eucdepth'],...
    ['by Sally Warner on ' datestr(now,'mmm dd, yyyy')]);

euchr.readme = strvcat(euc.readme,...
    ' ',...
    '10 min data has been binned into hourly averages');

eucdy.readme = strvcat(euc.readme,...
    ' ',...
    '10 min data has been binned into daily averages');

eucmo.readme = strvcat(euc.readme,...
    ' ',...
    'daily data has been averaged by month. A month nust have data from at',...
    'least 10 days to be averaged');

eucclm.readme = strvcat(euc.readme,...
    ' ',...
    'climatology is found by averaging monthly data.',...
    'i.e. All Januarys are averaged, all Februarys are averaged, etc.');

save([savedir 'euc_' latstrshort '_' lonstr '_10min.mat'],'euc')
save([savedir 'euc_' latstrshort '_' lonstr '_hourly.mat'],'euchr')
save([savedir 'euc_' latstrshort '_' lonstr '_daily.mat'],'eucdy')
save([savedir 'euc_' latstrshort '_' lonstr '_monthly.mat'],'eucmo')
save([savedir 'euc_' latstrshort '_' lonstr '_monthly_climatology.mat'],'eucclm')


return
%% plot

dd = dir([savedir 'vel*']);

for ii = 1:length(dd)
    load([savedir dd(ii).name])
end

varsrm = fields(euc);
vars = varsrm(2:end-1);


figure(265)
clf
    set(gcf,'position',[495          57         766        1048])


for ii = 1:length(vars)
    
    ax(ii) = subplot(length(vars),1,ii);
    
    plot(euc.time,euc.(vars{ii}))
    hold on
    plot(euchr.time,euchr.(vars{ii}))
    plot(eucdy.time,eucdy.(vars{ii}))
    plot(eucmo.time,eucmo.(vars{ii}))
    
    datetick
    
    ylabel(vars{ii})
    set(gca,'tickdir','out')
    
    if ii == 1
        title(['Equatorial undercurrent depth at ' latstrshort ', ' lonstr])
        legend('10min','hr','daily','monthly','location','northwest')
        legend boxoff
    end
end

export_fig([figdir 'euc.png'],'-r200')

