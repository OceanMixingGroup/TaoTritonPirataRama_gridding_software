% load_and_save_enso_oni_index
%
% load in previously processed enso index
% put on monthly time grid used by all other data
% save



clear
addpath('./utilities/')

%% define location
% this creates correct strings and folders

%%%%%%%%% be sure to set the correct location in set_location.m! %%%%%%%%%%

set_location;

%% load common tao time
% rather than using the common time grid which has a 10 minute interval,
% want to use monthly time for this code:

load([savedir 'sst_' latstrshort '_' lonstr '_monthly.mat'])
time = sstmo.time;
clear sstmo


%% load

%%%%% enso oni %%%%%
rawdatadir = '~/Dropbox/data/enso/enso_oni.mat';
load(rawdatadir)


%% interpolate to sstmo time grid

ensomo.time = time;
ensomo.indx = interp1(enso.time,enso.indx,time,'nearest');
ensomo.dindxdt = interp1(enso.time,enso.dindxdt,time,'nearest');



%% interpolate to the faster time grids

moyear = str2num(datestr(ensomo.time,'yyyy'));
momonth = str2num(datestr(ensomo.time,'mm'));

% monthly to daily
disp('interpolating monthly to daily...')
load([savedir 'sst_' latstrshort '_' lonstr '_daily.mat'])
years = str2num(datestr(sstdy.time,'yyyy'));
months = str2num(datestr(sstdy.time,'mm'));

ensody.time = sstdy.time;
ensody.indx = NaN*ensody.time;
ensody.dindxdt = NaN*ensody.time;
for ii = 1:length(ensody.time)
    year = years(ii);
    month = months(ii);
    ind = find(moyear == year & momonth == month);
    ensody.indx(ii) = ensomo.indx(ind);
    ensody.dindxdt(ii) = ensomo.dindxdt(ind);
end

% monthly to hourly
disp('interpolating monthly to hourly...')
load([savedir 'sst_' latstrshort '_' lonstr '_hourly.mat'])
years = str2num(datestr(ssthr.time,'yyyy'));
months = str2num(datestr(ssthr.time,'mm'));

ensohr.time = ssthr.time;
ensohr.indx = NaN*ensohr.time;
ensohr.dindxdt = NaN*ensohr.time;
for ii = 1:length(ensohr.time)
    year = years(ii);
    month = months(ii);
    ind = find(moyear == year & momonth == month);
    ensohr.indx(ii) = ensomo.indx(ind);
    ensohr.dindxdt(ii) = ensomo.dindxdt(ind);
end

% monthly to 10min
disp('interpolating monthly to 10min...')
time = taotime;
years = str2num(datestr(time,'yyyy'));
months = str2num(datestr(time,'mm'));

enso10.time = time;
enso10.indx = NaN*enso10.time;
enso10.dindxdt = NaN*enso10.time;
for ii = 1:length(enso10.time)
    year = years(ii);
    month = months(ii);
    ind = find(moyear == year & momonth == month);
    enso10.indx(ii) = ensomo.indx(ind);
    enso10.dindxdt(ii) = ensomo.dindxdt(ind);
end



%% save

disp('saving...')

enso.readme = strvcat('The Oceanic Nino Index (ONI) is calculated to show the relative',...
    'strength of El Niño-Souther Oscillation (ENSO) events.',...
    'It is not location dependent.',...
    ' ',...
    'More info about the ONI index can be found at:',...
    'https://origin.cpc.ncep.noaa.gov/products/analysis_monitoring/ensostuff/ONI_v5.php',...
    ['original data saved in: ' rawdatadir],...
    ' ',...
    ['created by ' pwd '/load_and_save_enso_oni_index.m'],...
    ['by Sally Warner on ' datestr(now,'mmm dd, yyyy')]);


ensomo.readme = strvcat(enso.readme,...
    ' ',...
    'data has been interpolated to chipod analysis monthly time grid,');
    
ensody.readme = strvcat(enso.readme,...
    ' ',...
    'data has been interpolated to chipod analysis daily time grid,');

ensohr.readme = strvcat(enso.readme,...
    ' ',...
    'data has been interpolated to chipod analysis hourly time grid,');

enso10.readme = strvcat(enso.readme,...
    ' ',...
    'data has been interpolated to chipod analysis 10 minute time grid,');

clear enso
enso = enso10;

save([savedir 'enso_oni_' latstrshort '_' lonstr '_10min.mat'],'enso')
save([savedir 'enso_oni_' latstrshort '_' lonstr '_hourly.mat'],'ensohr')
save([savedir 'enso_oni_' latstrshort '_' lonstr '_daily.mat'],'ensody')
save([savedir 'enso_oni_' latstrshort '_' lonstr '_monthly.mat'],'ensomo')


return

%% plot


dd = dir([savedir 'enso*']);

for ii = 1:length(dd)
    load([savedir dd(ii).name])
end

varsrm = fields(enso);
vars = varsrm(2:end-1);


figure(265)
clf
    set(gcf,'position',[495          57         766        1048])


for ii = 1:length(vars)
    
    ax(ii) = subplot(length(vars),1,ii);
    
    plot(enso.time,enso.(vars{ii}))
    hold on
    plot(ensohr.time,ensohr.(vars{ii}))
    plot(ensody.time,ensody.(vars{ii}))
    plot(ensomo.time,ensomo.(vars{ii}))
    
    datetick
    
    ylabel(vars{ii})
    set(gca,'tickdir','out')
    
    if ii == 1
        title(['ENSO ONI at 0, ' loc])
        legend('10min','hr','daily','monthly','location','northwest')
        legend boxoff
    end
end

export_fig([figdir 'enso.png'],'-r200')

