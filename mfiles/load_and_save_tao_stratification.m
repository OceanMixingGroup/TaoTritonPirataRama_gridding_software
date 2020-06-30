% load_and_save_tao_stratification
% based on TAO140/mfiles/stratification_0_140W.m
%
% The goal of this code is to input all temperature and salinity data from
% 0 140W and interpolate it onto a regular grid (the same grid that's used
% for load_and_save_tao_adcp_cur.m:
%    vel.time = taotime;
%    vel.depth = (5:5:300)';
%    vel.u = NaN*ones(length(vel.depth),length(vel.time));
%
% From the stratification data, calculate the buoyancy frequency
%
% Note that salinity data is sparce both in time and space compared to
% temperature. Compare stratification calculated only with temperature and
% to that calcualted both with salinity and temperature. Bill Smyth says it
% does not really make a big difference if salinity is used or not for 
% equatorial TAO locations. Especially, with the historical data, salinity 
% is very sparce.
%
% written by Sally Warner

clear
tic
addpath('./utilities/')

%% set parameters

% Do you want to save all three interpolation schemes?
%
% (Temperature is interpolated onto the even depth grid using three
% different interpolation schemes. Typically, only the best one is saved,
% but if you specifically want all three interpolation schemes, set
% allthree = 1;)
allthree = 0;


%% define location
% this creates correct strings and folders

%%%%%%%%% be sure to set the correct location in set_location.m! %%%%%%%%%%

set_location;


%% load common tao time

time = taotime;


%% load in temperature

% load both the 10min and the daily temperature data. THe daily data can
% fill holes in the 10min data (daily data is sent back via satellite, so if
% the mooring is vandalized and the instruments are lost, there will still
% be daily data)

%%%%% 10 minute T %%%%%

disp('loading 10 min and daily temperature netcdf files...')

filename = ['../../TaoTritonPirataRama/high_resolution/10m/t' ...
    lower(latstr) lower(lonstr) '_10m.cdf'];
if exist(filename,'file')
    rawd = loadnc(filename);
else
    disp('WARNING: load_and_save_tao_stratification cannot find raw temperature file!')
    disp(['looking for: ' filename])
    disp('Setting all temperature data to NaN.')
    disp(' ')
    rawd.time = time;
    rawd.depth = 0;
    rawd.T = NaN*time;
end
% shorten to just 2005 onward
raw.time = time;
raw.depth = rawd.depth;
for ii = 1:length(raw.depth)
    raw.T(ii,:) = interp1(rawd.time,rawd.T(ii,:),time);
end


%%%%% daily T %%%%%

gooddailydata = 1;

filename = ['../../TaoTritonPirataRama/' mooringarray '/t_xyzt_dy.cdf'];
if exist(filename,'file')
    dummy = loadnc(filename);
else
    gooddailydata = 0;
    disp('WARNING: load_and_save_tao_stratification could not find the')
    disp(['daily tempearture file: ../../TaoTritonPirataRama/' mooringarray '/t_xyzt_dy.cdf'])
    disp('No daily data will be used to fill gaps in the high resolution data.')        
    disp(' ')
end

if gooddailydata
    indlat = find(dummy.lat == lat);
    indlon = find(dummy.lon == lon);
    
    if isempty(indlat) | isempty(indlon)
        gooddailydata = 0;
        disp(['WARNING: load_and_save_tao_stratification cannot find the desired lat and lon'])
        disp(['in the daily ../../TaoTritonPirataRama/' mooringarray '/t_xyzt_dy.cdf. '])
        disp('No daily data will be used to fill gaps in the high resolution data')
        disp(' ')
    else
        rawddy.time = dummy.time;
        rawddy.depth = dummy.depth;
        rawddy.lat = dummy.lat(indlat);
        rawddy.lon = dummy.lon(indlon);
        rawddy.T = squeeze(dummy.T(:,:,indlat,indlon))';
        
        % interpolate to 10 minute time grid
        rawdy.time = time;
        rawdy.depth = rawddy.depth;
        rawdy.lat = rawddy.lat;
        rawdy.lon = rawddy.lon;
        for ii = 1:length(rawdy.depth)
            rawdy.T(ii,:) = interp1(rawddy.time,rawddy.T(ii,:),time);
        end
	end
end



%% load in salinity

% salinity is likely quite sparce


%%%%% hourly S %%%%%

disp('loading hourly and daily salinity netcdf files...')

filename = ['../../TaoTritonPirataRama/high_resolution/hr/s' ...
    lower(latstr) lower(lonstr) '_hr.cdf'];
if exist(filename,'file')
    rawds = loadnc(filename);
else
    disp('WARNING: load_and_save_tao_stratification cannot find raw salinity file!')
    disp(['looking for: ' filename])
    disp('Setting all salinity data to NaN.')
    disp(' ')
    rawds.time = time;
    rawds.depth = 0;
    rawds.S = NaN*time;
end
% shorten to just 2005 onward
raws.time = time;
raws.depth = rawds.depth;
for ii = 1:length(raws.depth)
    raws.S(ii,:) = interp1(rawds.time,rawds.S(ii,:),time);
end


%%%%% daily S %%%%%

gooddailydataS = 1;

filename = ['../../TaoTritonPirataRama/' mooringarray '/s_xyzt_dy.cdf'];
if exist(filename,'file')
    dummy = loadnc(filename);
else
    gooddailydataS = 0;
    disp('WARNING: load_and_save_tao_stratification could not find the')
    disp(['daily salinity file: ../../TaoTritonPirataRama/' mooringarray '/s_xyzt_dy.cdf'])
    disp('No daily data will be used to fill gaps in the high resolution data.')        
    disp(' ')
end

if gooddailydataS
    indlat = find(dummy.lat == lat);
    indlon = find(dummy.lon == lon);
    
    if isempty(indlat) | isempty(indlon)
        gooddailydataS = 0;
        disp(['WARNING: load_and_save_tao_stratification cannot find the desired lat and lon'])
        disp(['in the daily ../../TaoTritonPirataRama/' mooringarray '/s_xyzt_dy.cdf. '])
        disp('No daily data will be used to fill gaps in the high resolution data')
        disp(' ')
    else
        rawsddy.time = dummy.time;
        rawsddy.depth = dummy.depth;
        rawsddy.lat = dummy.lat(indlat);
        rawsddy.lon = dummy.lon(indlon);
        rawsddy.S = squeeze(dummy.S(:,:,indlat,indlon))';
        
        % interpolate to 10 minute time grid
        rawsdy.time = time;
        rawsdy.depth = rawsddy.depth;
        rawsdy.lat = rawsddy.lat;
        rawsdy.lon = rawsddy.lon;
        for ii = 1:length(rawsddy.depth)
            rawsdy.S(ii,:) = interp1(rawsddy.time,rawsddy.S(ii,:),time);
        end
	end
end


%% fill holes in 10min data with daily data

if gooddailydata

    disp('filling holes in 10min temperature data with daily data...')

    % if you only want to fill NaNs in high-temporal-resolution data with daily
    % data below a specific depth, you can set the minimum depth (mninfilldepth)
    % to be greater than 1. 

    minfilldepth = 0;
    indstart = find(raw.depth >= minfilldepth,1,'first'); 

    for ii = indstart:length(raw.depth)
        clear inddepth
        inddepth = find(rawdy.depth == raw.depth(ii));
        if ~isempty(inddepth)
            if raw.depth(ii) == rawdy.depth(inddepth)
                indbad = isnan(raw.T(ii,:));
                raw.T(ii,indbad) = rawdy.T(inddepth,indbad);
            end            
        end    
    end    
end


if gooddailydataS

    disp('filling holes in hourly salinity data with daily data...')

    % if you only want to fill NaNs in high-temporal-resolution data with daily
    % data below a specific depth, you can set the minimum depth (mninfilldepth)
    % to be greater than 1. 

    minfilldepth = 0;
    indstart = find(raws.depth >= minfilldepth,1,'first'); 

    for ii = indstart:length(raws.depth)
        clear inddepth
        inddepth = find(rawsdy.depth == raws.depth(ii));
        if ~isempty(inddepth)
            if raws.depth(ii) == rawsdy.depth(inddepth)
                indbad = isnan(raws.S(ii,:));
                raws.S(ii,indbad) = rawsdy.S(inddepth,indbad);
            end            
        end    
    end    
end


%% fill gaps horizontally (in time) also 

disp('filling small gaps in time...')

maxgap = 6*24; % 6*24 is a gap length of one day

for ii = 1:length(raw.depth)
    raw.Tfill(ii,:) = fast_fillgap(raw.T(ii,:),6*24);
end


for ii = 1:length(raws.depth)
    raws.Sfill(ii,:) = fast_fillgap(raws.S(ii,:),6*24);
end
    

%% fill holes in 1 m data with 5 m data

% % % there is a period from January through March 2015 where the 1m sensor
% % % fails but the 5m sensor is good. For the period before and after, the
% % % temperature at these two depths is nearly equal. Fill gaps in 1m data
% % % with 5m data
% % 
% % nan1 = ~isnan(raw.Tfill(1,:));
% % nan5 = ~isnan(raw.Tfill(2,:));
% % indsub = find(nan1 == 0 & nan5 == 1);
% % raw.Tfill(1,indsub) = raw.Tfill(2,indsub);


%% put 1m data at 0m

% Usually, the raw data from the Tao moorings is measured at 1m, 5m and
% deeper. when interpolating from the raw data to an even depth grid that
% goes from 0:5:300 m, the data from the 1m sensor is often dropped.
% Therefore, put a new row at 0m that copies the 1m raw data. In the next
% step, when the raw data is interpolated to an even depth grid, the 1m
% won't be left out.
%
% (The other alternative of having an uneven depth grid of [1 5:5:300]
% could work well for temperature and stratification, but it really doesn't
% work well for the ADCP and I want to keep the depth grids the same for
% these two instruments.)

raw.depth = [0; raw.depth];
raw.Tfill = [raw.Tfill(1,:); raw.Tfill];

raws.depth = [0; raws.depth];
raws.Sfill = [raws.Sfill(1,:); raws.Sfill];

%% vertical interpolation to even depth grid

disp('interpolating to regularly spaced depth grid...')

% create structure for gridded data
tnan.time = time;

% use the common depth grid set for all ADCP and density data available
[tnan.depth,dz] = taodepth;

tnan.tlinear = NaN*ones(length(tnan.depth),length(tnan.time));
tnan.takima = tnan.tlinear;
tnan.takimabin = tnan.tlinear;

tnan.slinear = tnan.tlinear;
tnan.sakima = tnan.tlinear;
tnan.sakimabin = tnan.tlinear;



for ii = 1:length(tnan.time)
    
    clear depth tt indgood depthdiff indavg depth1 tt1 depthavg ttavg
    
    if mod(ii,10000) == 0
        disp([num2str(ii) '/' num2str(length(tnan.time))])
    end
    
    indgood = ~isnan(raw.Tfill(:,ii));
    if sum(indgood) > 1
        dd = raw.depth(indgood);
        tt = raw.Tfill(indgood,ii);
        tnan.tlinear(:,ii) = interp1(dd,tt,tnan.depth);
    
        % (sjw 7/14/14) Talked with Bill this morning. He suggested using a
        % cubic Akima interpolation rather than linear interpolation. This
        % will make the dataset smoother and therefore better for
        % calculating N^2.
        tnan.takima(:,ii) = cakima(dd,tt,tnan.depth);
        
        % (sjw 7/15/14) In contining my discussion with Bill about the best way to
        % interpolate the data, we decided that places that have multiple sensors
        % close together need to be meaned before the akima spline is used. 
        %              1 2 3  4  5  6  7  8  9  10 11 12 13 14 15  16  17  18  19  20  21
        % raw.depth = [1 5 10 13 20 25 28 40 45 48 60 79 80 83 100 120 123 140 180 300 500];
        % So, any case where good sensors are 3m or less apart, they will
        % be averaged together. Their depth will also be averaged togehter.
        
        depthdiff = diff(dd);
        indavg = find(depthdiff <= 3);
            
        dd1 = dd;
        tt1 = tt;
        for nn = 1:length(indavg)
            if dd(indavg(nn)) > 0 %& dd(indavg(nn)+1) == 1
                % keep the bins at 0m and 1m in the spline to prevent large
                % gradients in the mixed layer created by cakima
%                 dd1(indavg(nn)) = 0;
%                 tt1(indavg(nn)) = tt(indavg(nn));
%                 dd1(indavg(nn)+1) = 1;
%                 tt1(indavg(nn)+1) = tt(indavg(nn)+1);                
%             else
                dd1(indavg(nn)) = nanmean(dd(indavg(nn):indavg(nn)+1));
                dd1(indavg(nn)+1) = NaN;
%                 if dd1(indavg(nn)) == 0.5
%                     dd1(indavg(nn)) = 0; % for 0m and 1m bins, want dd1(1) = 0;
%                 end
                tt1(indavg(nn)) = nanmean(tt(indavg(nn):indavg(nn)+1));
                tt1(indavg(nn)+1) = NaN;
            end
        end
        ddavg = dd1(~isnan(dd1));
        ttavg = tt1(~isnan(tt1));  
        
        tnan.takimabin(:,ii) = cakima(ddavg,ttavg,tnan.depth);
        
    elseif sum(indgood) == 1 % case with only one value
        clear indd1
        indd1 = find(tnan.depth <= raw.depth(indgood),1,'last');
        tnan.slinear(indd1,ii) = raw.Tfill(indgood,ii);
        tnan.sakima(indd1,ii) = raw.Tfill(indgood,ii);
        tnan.sakimabin(indd1,ii) = raw.Tfill(indgood,ii);
    end
    
 %%%%% salinity %%%%%
    clear depths ss indgoods depthdiffs indavgs depth1s ss1 depthavgs ssavg

    indgoods = ~isnan(raws.Sfill(:,ii));
    if sum(indgoods) > 1
        dds = raws.depth(indgoods);
        ss = raws.Sfill(indgoods,ii);
        tnan.slinear(:,ii) = interp1(dds,ss,tnan.depth);

        tnan.sakima(:,ii) = cakima(dds,ss,tnan.depth);

        depthdiffs = diff(dds);
        indavgs = find(depthdiffs <= 3);

        dds1 = dds;
        ss1 = ss;
        for nn = 1:length(indavgs)
            dds1(indavgs(nn)) = nanmean(dds(indavgs(nn):indavgs(nn)+1));
            dds1(indavgs(nn)+1) = NaN;
            if dds1(indavgs(nn)) == 0.5
                dds1(indavgs(nn)) = 0; % for 0m and 1m bins, want dd1(1) = 0;
            end

            ss1(indavgs(nn)) = nanmean(ss(indavgs(nn):indavgs(nn)+1));
            ss1(indavgs(nn)+1) = NaN;
        end
        ddavgs = dds1(~isnan(dds1));
        ssavg = ss1(~isnan(ss1));  

        if length(ddavgs) > 1
            tnan.sakimabin(:,ii) = cakima(ddavgs,ssavg,tnan.depth); 
        else
            clear indd1s
            indd1s = find(tnan.depth <= ddavgs,1,'last');
            tnan.sakimabin(indd1s,ii) = ssavg;
        end
    elseif sum(indgoods) == 1 % case with only one value
        clear indd1s
        indd1s = find(tnan.depth <= raws.depth(indgoods),1,'last');
        tnan.slinear(indd1s,ii) = raws.Sfill(indgoods,ii);
        tnan.sakima(indd1s,ii) = raws.Sfill(indgoods,ii);
        tnan.sakimabin(indd1s,ii) = raws.Sfill(indgoods,ii);
    end
end




%% plot comparision between linear interpolation and akima spline

if allthree == 1
    
    lw1 = 2;
    lw2 = 2;
    col1 = 'k';
    col2 = nicecolor('G');
    col3 = nicecolor('cb');
    col4 = nicecolor('r');
    
    % choose three time points you want to plot
    indp = [70000 0 100000 0 120000];
%     indp = [305000 0 505000 0 705000];
    


    %%%%% temperature %%%%%
    figure(1)
    clf
    set(gcf,'position',[92         561        1710         531])
    
    for ii = 1:2:6
    
    axT(ii) = subplot(1,6,ii);
    plot(tnan.tlinear(:,indp(ii)),-tnan.depth,'color',col1,'linewidth',lw1)
    hold on
    plot(tnan.takima(:,indp(ii)),-tnan.depth,'-','color',col2,'linewidth',lw2)
    plot(tnan.takimabin(:,indp(ii)),-tnan.depth,'-','color',col3,'linewidth',lw2)
    plot(raw.Tfill(:,indp(ii)),-raw.depth,'o','color',col4,'linewidth',lw2)
    title(datestr(tnan.time(indp(ii))))
    if ii == 1
        ylabel('depth [m]')
        legend('linear','akima','akima w/bins','location','southeast')
        legend('boxoff')
    end
    xlabel('intepolated T [°C]')
    ylim([-300 0])
    xlim([10 30])
    
    axT(ii+1) = subplot(1,6,ii+1);
    plot(gradient(tnan.tlinear(:,indp(ii)),-5),-tnan.depth,'color',col1,'linewidth',lw1)
    hold on
    plot(gradient(tnan.takima(:,indp(ii)),-5),-tnan.depth,'-','color',col2,'linewidth',lw2)
    plot(gradient(tnan.takimabin(:,indp(ii)),-5),-tnan.depth,'-','color',col3,'linewidth',lw2)
    title(datestr(tnan.time(indp(ii))))
    xlabel('dT/dz [°C m^{-1}]')
    ylim([-300 0])
    xlim([-0.01 0.21])

    end
    linkaxes(axT,'y')
    
%     export_fig(figdir 'compare_linear_vs_akimaspline_T.pdf'])



    %%%%% salinity %%%%%

    figure(2)
    clf
    set(gcf,'position',[92         561        1710         531])
    
    for ii = 1:2:6
    
    axS(ii) = subplot(1,6,ii);
    plot(tnan.slinear(:,indp(ii)),-tnan.depth,'-','color',col1,'linewidth',lw1)
    hold on
    plot(tnan.sakima(:,indp(ii)),-tnan.depth,'-','color',col2,'linewidth',lw2)
    plot(tnan.sakimabin(:,indp(ii)),-tnan.depth,'-','color',col3,'linewidth',lw2)
    plot(raws.Sfill(:,indp(ii)),-raws.depth,'o','color',col4,'linewidth',lw2)
    title(datestr(tnan.time(indp(ii))))
    if ii == 1
        ylabel('depth [m]')
        legend('linear','akima','akima w/bins','S measurements','location','southeast')
        legend('boxoff')
    end
    xlabel('intepolated S [psu]')
    ylim([-300 0])
    xx = minmax(tnan.sakima(:,indp(ii)));
    if isnan(xx)
        xx = [0 1];
    end
    xlim([xx(1)*0.999 xx(2)*1.001])
    
    axS(ii+1) = subplot(1,6,ii+1);
    plot(gradient(tnan.slinear(:,indp(ii)),-5),-tnan.depth,'color',col1,'linewidth',lw1)
    hold on
    plot(gradient(tnan.sakima(:,indp(ii)),-5),-tnan.depth,'-','color',col2,'linewidth',lw2)
    plot(gradient(tnan.sakimabin(:,indp(ii)),-5),-tnan.depth,'-','color',col3,'linewidth',lw2)
    title(datestr(tnan.time(indp(ii))))
    xlabel('dS/dz [psu m^{-1}]')
    ylim([-300 0])
    xxg = minmax(gradient(tnan.sakimabin(:,indp(ii)),-5));
    if isnan(xxg)
        xxg = [0 1];
    end
    xlim(xxg*1.1)

    end
    linkaxes(axS,'y')
    
%     export_fig([figdir 'compare_linear_vs_akimaspline_S.pdf'])

end



%% calcualte buoyancy frequency

disp('calculating buoyancy frequencies...')

% calculate stratification and buoyancy frequency with just temperature
[~,tnan.dTdzlinear] = gradient(tnan.tlinear,-dz);
tnan.NTsqlinear = 9.81*sw_alpha(30+0*tnan.tlinear,tnan.tlinear,0,'ptmp').*tnan.dTdzlinear;

[~,tnan.dTdzakima] = gradient(tnan.takima,-dz);
tnan.NTsqakima = 9.81*sw_alpha(30+0*tnan.takima,tnan.takima,0,'ptmp').*tnan.dTdzakima;

[~,tnan.dTdzakimabin] = gradient(tnan.takimabin,-dz);
tnan.NTsqakimabin = 9.81*sw_alpha(30+0*tnan.takimabin,tnan.takimabin,0,'ptmp').*tnan.dTdzakimabin;


% calculate stratification and buoyancy frequency with just salinity
[~,tnan.dSdzlinear] = gradient(tnan.slinear,-dz);
tnan.NSsqlinear = -9.81*sw_beta(tnan.slinear,tnan.tlinear,tnan.depth,'ptmp').*tnan.dSdzlinear;

[~,tnan.dSdzakima] = gradient(tnan.sakima,-dz);
tnan.NSsqakima = -9.81*sw_beta(tnan.sakima,tnan.takima,tnan.depth,'ptmp').*tnan.dSdzakima;

[~,tnan.dSdzakimabin] = gradient(tnan.sakimabin,-dz);
tnan.NSsqakimabin = -9.81*sw_beta(tnan.sakimabin,tnan.takimabin,tnan.depth,'ptmp').*tnan.dSdzakimabin;



% calculate stratification and buoyancy frequency with both temperature and salinity
% note that the salinity_akima is often very spurious due to lack of
% measurements. It's likely most correct to use the linear option for
% salinity paired with the akimabin option for temperature.
tnan.rho = sw_dens(tnan.slinear,tnan.takimabin,tnan.depth);
tnan.sigma = sw_pden(tnan.slinear,tnan.takimabin,tnan.depth,0);
[~,tnan.dsigmadz] = gradient(tnan.sigma,-dz);
tnan.Nsq = -9.81./tnan.sigma.*tnan.dsigmadz;


%% calculate depth of 20°C isotherm

disp('calculating depth of 20°C isotherm...')

tic
tnan.depth20 = NaN*tnan.time;
for tt = 1:length(tnan.time)
    if mod(tt,10000) == 0
        disp(num2str(tt))
    end
   clear good ind1 ind2
    % check if there are non-nan values of temperature at that time
    good = ~isnan(tnan.takimabin(:,tt));
    if sum(good) == 0
        tnan.depth20(tt) = NaN;
    else        
        % linearly interpolate the depth of the 20deg isotherm
        ind1 = find(tnan.takimabin(:,tt) > 20,1,'last');
        ind2 = find(tnan.takimabin(:,tt) < 20,1,'first');
        
        if isempty(ind2) | isempty(ind1)
            tnan.depth20(tt) = NaN;
        else
        tnan.depth20(tt) = (tnan.depth(ind2)-tnan.depth(ind1))./...
            (tnan.takimabin(ind2,tt)-tnan.takimabin(ind1,tt))...
            .*(20-tnan.takimabin(ind1,tt)) + tnan.depth(ind1);
        end
    end
end
toc


%% calculate mixed layer depth

%%%%%% interpolate to find the MLD %%%%%%

% deltaT is the change in temperature from the 1m bin to the bottom of the
% mixed layer. This value of 0.015°C comes from the literature, such as 
% Moum et al., 2013. We interpolate to find the depth that is deltaT
% degrees colder than the 1m temperature.
deltaT = 0.015;

mldt = tnan.takimabin(1,:) - deltaT;
tnan.mld = NaN*tnan.time;

disp('calculating first of two mixed layer depths...')
for ii = 1:length(tnan.time)
    
    if mod(ii,10000) == 0
        disp([num2str(ii) '/' num2str(length(tnan.time))])
    end
    
    clear ind1 ind2 d1 d2 t1 t2

    ind1 = find(tnan.takimabin(:,ii) > mldt(ii),1,'last');
    ind2 = find(tnan.takimabin(:,ii) < mldt(ii),1,'first');
    
    if ~isempty(ind1) && ~isempty(ind2)
    
        d1 = tnan.depth(ind1);
        d2 = tnan.depth(ind2);

        t1 = tnan.takimabin(ind1,ii);
        t2 = tnan.takimabin(ind2,ii);

        % linearly interpolate to find depth of mixed layer
        tnan.mld(ii) = d1 - (d1-d2)/(t1-t2)*(t1-mldt(ii));
    
    end  
end






%%%%%% repeat the MLD interpolation with deltaT = 0.04K %%%%%%
% this is for Bill's velocity and shear extrapolation calculation

disp('calculating second of two mixed layer depths...')

deltaTuz = 0.04;

mldtuz = tnan.takimabin(1,:) - deltaTuz;
tnan.mlduz = NaN*tnan.time;

for ii = 1:length(tnan.time)
    
    if mod(ii,10000) == 0
        disp([num2str(ii) '/' num2str(length(tnan.time))])
    end
    
    clear ind1 ind2 d1 d2 t1 t2

    ind1 = find(tnan.takimabin(:,ii) > mldtuz(ii),1,'last');
    ind2 = find(tnan.takimabin(:,ii) < mldtuz(ii),1,'first');
    
    if ~isempty(ind1) && ~isempty(ind2)
    
        d1 = tnan.depth(ind1);
        d2 = tnan.depth(ind2);

        t1 = tnan.takimabin(ind1,ii);
        t2 = tnan.takimabin(ind2,ii);

        tnan.mlduz(ii) = d1 - (d1-d2)/(t1-t2)*(t1-mldtuz(ii));
    
    end  
end




%% create structure that will be saved
% - use akima bin interpolation for temperature becuase there are usuall
% enough temperature measurements to make this nicely smooth
% - use linear interpolation for salinity becuase there just aren't enough
% salinity measurements which results in often very erronious salinity
% profiles
% - use the same combination for calculating rho, sigma and Nsq

t.time  = time;
t.depth = tnan.depth;

t.T     = tnan.takimabin;
t.dTdz  = tnan.dTdzakimabin;
t.NTsq  = tnan.NTsqakimabin;

t.S     = tnan.slinear;
t.dSdz  = tnan.dSdzlinear;
t.NSsq  = tnan.NSsqlinear;

t.rho   = tnan.rho;
t.sigma = tnan.sigma;
t.Nsq   = tnan.Nsq;

t.depth20 = tnan.depth20;
t.mld = tnan.mld;
t.mld_forshear = tnan.mlduz;


%% interpolate to hourly, daily and monthly time grids

disp('calculating hourly, daily, and monthly averages...')

% hourly average (dt = 10 min, so put 6 points into each hour)
thr = bin_average(t,6);
thr.depth = t.depth;

% daily average (wind dt = 10 min, so put 6*24=144 points into each average)
tdy = bin_average(t,144);
tdy.depth = t.depth;

% monthly averages
tmo = monthly_average(tdy);
tmo.depth = t.depth;



%% save best data

disp('saving...')

endtimestr = datestr(time(end),'mmmm yyyy');

t.readme = strvcat(['Temperature and salinity from ' mooringarray ' mooring at' latstrshort ', ' lonstr],...
    '',...
    'Temperature (T, in °C) is measured every 10 minutes.',...
    'At times when there is no high resolution temperature, daily',...
    '     temperature is used to fill gaps.',... 
    ['    (The minimum fill gap depth is ' num2str(minfilldepth) 'm)'],...
    'Temperature is interpolated to 5m depth grid from measurement depths of',...
    ['    ' num2str(rawd.depth')],...
    '      using a cubic Akima spline interpolation scheme (Akima 1970 & Pham et al., 2017)',...
    ' ',...
    'Salinity (S, in psu) is measured hourly.',...
    'At times when there is no high resolution salinity, daily',...
    '     salinity is used to fill gaps.',... 
    ['    (The minimum fill gap depth is ' num2str(minfilldepth) 'm)'],...
    'Salinity is interpolated to 5m depth grid from measurement depths of:',...
    ['    ' num2str(rawds.depth')],...
    '      using a linear interpolation scheme. (There are too few',...
    '      measurements to use a cubic akima spline interpolation.)',...
    ' ',...
    'rho is density [kg m^-3]',...
    'sigma is potential density [kg m^-3]',...
    ' ',...
    'dTdz is the vertical gradient of temperature (dT/dz) [°C/m]',...
    'dSdz is the vertical gradient of salinity (dS/dz) [psu/m]',...
    'dsigmadz is the vertical gradient of potential density (dsigma/dz) [kg m^-4]',...
    'NTsq is the buoyancy frequency calculated only with temperature,',...
    '    NTsq = g*alpha*dT/dz [s^-2]',...
    'NSsq is the buoyancy frequency calculated only with salinity,',...
    '    NSsq = -g*beta*dS/dz [s^-2]',...
    'Nsq is the buoyancy frequency calculated with BOTH salinity and temperature',...
    '    Nsq = -g./sigma*dsigma/dz [s^-2]',...
    '** For many moorings, there will be WAY better temporal and depth',...
    '    coverage using NTsq rather than Nsq***',...
    ' ',...
    'The mixed layer depth (mld) is calculated as the depth at which the',...
    ['    temperature is ' num2str(deltaT) '°C colder than the 1m temperature.'],...
    'A second mixed layer depth (ml_forshear) is calculated as the depth at which the',...
    ['    temperature is ' num2str(deltaTuz) '°C colder than the 1m temperature.'],...
    '    (This version is used for the interpolation of velocity to the ',...
    '    surface following Pham et al 2017).',...
    ' ',...
    'time is matlab datenum',...
    'depth in m (5 m bins)',...
    ' ',...
    ['created by ' pwd '/load_and_save_tao_stratification.m'],...
    ['by Sally Warner on ' datestr(now,'mmm dd, yyyy')]);

thr.readme = strvcat(t.readme,...
    ' ',...
    '10 min data has been binned into hourly averages');

tdy.readme = strvcat(t.readme,...
    ' ',...
    '10 min data has been binned into daily averages');

tmo.readme = strvcat(t.readme,...
    ' ',...
    'daily data has been averaged by month. A month nust have data from at',...
    'least 10 days to be averaged');

save([savedir 'T_' latstrshort '_' lonstr '_10min.mat'],'t','-v7.3')
save([savedir 'T_' latstrshort '_' lonstr '_hourly.mat'],'thr','-v7.3')
save([savedir 'T_' latstrshort '_' lonstr '_daily.mat'],'tdy')
save([savedir 'T_' latstrshort '_' lonstr '_monthly.mat'],'tmo')

%% save three interpolation schemes (at 10minute time interval only)

if allthree == 1
    
    tall.time = tnan.time;
    tall.depth = tnan.depth;
    tall.T_linear = tnan.tlinear;
    tall.T_akima = tnan.takima;
    tall.T_akimabin = tnan.takimabin;
    tall.dTdz_linear = tnan.dTdzlinear;
    tall.dTdz_akima = tnan.dTdzakima;
    tall.dTdz_akimabin = tnan.dTdzakimabin;
    tall.NTsq_linear = tnan.NTsqlinear;
    tall.NTsq_akima = tnan.NTsqakima;
    tall.NTsq_akimabin = tnan.NTsqakimabin;
    tall.S_linear = tnan.slinear;
    tall.S_akima = tnan.sakima;
    tall.S_akimabin = tnan.sakimabin;
    tall.dSdz_linear = tnan.dSdzlinear;
    tall.dSdz_akima = tnan.dSdzakima;
    tall.dSdz_akimabin = tnan.dSdzakimabin;
    tall.NSsq_linear = tnan.NSsqlinear;
    tall.NSsq_akima = tnan.NSsqakima;
    tall.NSsq_akimabin = tnan.NSsqakimabin;
    tall.rho = tnan.rho;
    tall.sigma = tnan.sigma;
    tall.dsigmadz = tnan.dsigmadz;
    tall.Nsq = tnan.Nsq;

    tall.readme = strvcat(['Temperature and salinity from ' mooringarray ' mooring at' latstrshort ', ' lonstr],...
    '',...
    'Temperature (T, in °C) is measured every 10 minutes.',...
    'At times when there is no high resolution temperature, daily',...
    '     temperature is used to fill gaps.',... 
    ['    (The minimum fill gap depth is ' num2str(minfilldepth) 'm)'],...
    'Temperature is interpolated to 5m depth grid from measurement depths of',...
    ['    ' num2str(rawd.depth')],...
    '      using THREE interpolation scheme: linear, cubic Akima spline,',... 
    '      and binned cubic Akima spline (which averages measurements 3m or less apart).',...
    ' ',...
    'Salinity (S, in psu) is measured hourly.',...
    'At times when there is no high resolution salinity, daily',...
    '     salinity is to fill gaps.',... 
    ['    (The minimum fill gap depth is ' num2str(minfilldepth) 'm)'],...
    'Salinity is interpolated to 5m depth grid from measurement depths of:',...
    ['    ' num2str(rawds.depth')],...
    '      using the same three interpolation schemes as for temperature.',...
    '      WARNING: be cautious using any interpolation except for linear',...
    '      for salinity as the sparceness of salinity measurements often',...
    '      leads to spurious results for binned data.',...
    ' ',...
    'rho is density [kg m^-3]',...
    'sigma is potential density [kg m^-3]',...
    '      both rho and sigma are calculated using the seawater toolbox',...
    '      rho = sw_dens(S,T,depth);',...
    '      sigma = sw_pden(S,T,depth);',...
    ' ',...
    'dTdz is the vertical gradient of temperature (dT/dz) [°C/m]',...
    'dSdz is the vertical gradient of salinity (dS/dz) [psu/m]',...
    'dsigmadz is the vertical gradient of potential density (dsigma/dz) [kg m^-4]',...
    'NTsq is the buoyancy frequency calculated only with temperature,',...
    '    NTsq = g*alpha*dT/dz [s^-2]',...
    'NSsq is the buoyancy frequency calculated only with salinity,',...
    '    NSsq = -g*beta*dS/dz [s^-2]',...
    'Nsq is the buoyancy frequency calculated with BOTH salinity and temperature',...
    '    Nsq = -g./sigma*dsigma/dz [s^-2]',...
    '** For many moorings, there will be WAY better temporal and depth',...
    '    coverage using NTsq rather than Nsq***',...
    ' ',...
    'time is matlab datenum',...
    'depth in m (5 m bins)',...
    ' ',...
    'The three different interpolation schemes were used to find the',...
    'best way to interpolate the TAO temperature to a regular grid:',...
    '    - linear interpolation',...
    '    - akima spline interpolation (Akima 1970 & Pham et al., 2017)',...
    '    - akima spline interpolation with sensors less than 3m apart averaged',...
    'Of these three, the best is the akima binned for temperature because',...
    'it has a continuous derivative and no erronious wiggles.',...
    'Of the three, the linear interpolation is probably best for salinity,',...
    'especially on moorings where there were few salinity measurements.',...
    ' ',...
    '*** ALL THREE SCHEMES SAVED IN THIS .mat FILE ***',...
    ['*** USE T_' latstrshort '_' lonstr '_10min.mat FOR BEST DATA ***'],...
    ' ',...
    ['created by ' pwd '/load_and_save_tao_stratification.m'],...
    ['by Sally Warner on ' datestr(now,'mmm dd, yyyy')]);

    t = tall;
    save([savedir 'T_' latstrshort '_' lonstr '_10min_three_interpolation_schemes.mat'],'t','-v7.3')

end


%% save MLD

mld.time = tnan.time;
mld.mld = tnan.mld;
mld.mld_forshear = tnan.mlduz;

mldhr.time = thr.time;
mldhr.mld = thr.mld;
mldhr.mld_forshear = thr.mld_forshear;

mlddy.time = tdy.time;
mlddy.mld = tdy.mld;
mlddy.mld_forshear = tdy.mld_forshear;

mldmo.time = tmo.time;
mldmo.mld = tmo.mld;
mldmo.mld_forshear = tmo.mld_forshear;


mld.readme = strvcat(['Mixed Layer Depth at ' latstrshort ', ' lonstr],...
    'Calculated from 10 min TAO temperature data',...
    ' ',...
    'Bottom of mixed layer is taken to be the depth of temperature',...
    ['   that is ' num2str(deltaT) '°C less than the 1m TAO temperature.'],...
    'MLD in m below the surface',...
    'time in datenum time',...
    '  ',...
    'mld_forshear uses a temperature difference of 0.04°C. This is used in',...
    '    the calculation that extrapolates velocity/shear data to surface.',...
    '  ',...
    ['created by ' pwd '/load_and_save_tao_stratification'],...
    ['by Sally Warner on ' datestr(now,'mmm dd, yyyy')]);

mldhr.readme = strvcat(mld.readme,...
    ' ',...
    '10 min data has been binned into hourly averages');

mlddy.readme = strvcat(mld.readme,...
    ' ',...
    '10 min data has been binned into daily averages');

mldmo.readme = strvcat(mld.readme,...
    ' ',...
    'daily data has been averaged by month. A month nust have data from at',...
    'least 10 days to be averaged');

save([savedir 'mld_' latstrshort '_' lonstr '_10min.mat'],'mld')
save([savedir 'mld_' latstrshort '_' lonstr '_hourly.mat'],'mldhr')
save([savedir 'mld_' latstrshort '_' lonstr '_daily.mat'],'mlddy')
save([savedir 'mld_' latstrshort '_' lonstr '_monthly.mat'],'mldmo')



return

%% plot

load([savedir 'T_' latstrshort '_' lonstr '_daily.mat'])


%%%%%%%%% plot temperature %%%%%%%%%
clear ax

figure(346)
clf
set(gcf,'position',[512          57        1117        1048])


ax(1) = subplot(311);
pcolor(tdy.time,tdy.depth,tdy.T)
shading flat
axis ij
hold on
% plot(tdy.time,tdy.depth20,'k')
colormap(gca,flipud(lbmap(24,'RedBlue')))
datetick
colorbar
title(['Temperature and stratification at ' latstrshort ', ' lonstr])
ylabel({'Temperature [°C]';'depth [m]'})
set(gca,'tickdir','out')

ax(2) = subplot(312);
pcolor(tdy.time,tdy.depth,tdy.dTdz)
shading flat
axis ij
colormap(gca,flipud(lbmap(24,'RedBlue')))
datetick
colorbar
caxis([-1 1]*0.25)
ylabel({'dT/dz [°C m^{-1}]';'depth [m]'})
set(gca,'tickdir','out')


ax(3) = subplot(313);
pcolor(tdy.time,tdy.depth,tdy.NTsq)
shading flat
axis ij
colormap(gca,flipud(lbmap(24,'RedBlue')))
datetick
colorbar
caxis([-1 1]*1e-3)
ylabel({'N_T^2 [s^{-2}]';'depth [m]'})
set(gca,'tickdir','out')

linkaxes

export_fig([figdir 'Temperature.png'],'-r200')


%%%%%%%%% plot salinity %%%%%%%%%
clear ax

figure(347)
clf
set(gcf,'position',[512          57        1117        1048])


ax(1) = subplot(311);
pcolor(tdy.time,tdy.depth,tdy.S)
shading flat
axis ij
hold on
colormap(gca,flipud(lbmap(24,'RedBlue')))
datetick
colorbar
title(['Salinity and stratification at ' latstrshort ', ' lonstr])
ylabel({'Salinity [psu]';'depth [m]'})
set(gca,'tickdir','out')

ax(2) = subplot(312);
pcolor(tdy.time,tdy.depth,tdy.dSdz)
shading flat
axis ij
colormap(gca,flipud(lbmap(24,'RedBlue')))
datetick
colorbar
caxis([-1 1]*0.25)
ylabel({'dS/dz [psu m^{-1}]';'depth [m]'})
set(gca,'tickdir','out')


ax(3) = subplot(313);
pcolor(tdy.time,tdy.depth,tdy.NSsq)
shading flat
axis ij
colormap(gca,flipud(lbmap(24,'RedBlue')))
datetick
colorbar
caxis([-1 1]*1e-3)
ylabel({'N_S^2 [s^{-2}]';'depth [m]'})
set(gca,'tickdir','out')

linkaxes

export_fig([figdir 'Salinity.png'],'-r200')





%%%%%%%%% compare NTsq, NSsq, Nsq %%%%%%%%%
clear ax

figure(348)
clf
set(gcf,'position',[512          57        1117        1048])


ax(1) = subplot(511);
pcolor(tdy.time,tdy.depth,tdy.NTsq)
shading flat
axis ij
hold on
colormap(gca,flipud(lbmap(24,'RedBlue')))
datetick
colorbar
caxis([-1 1]*1e-3)
title(['Compare buoyancy frequencies (NTsq, NSsq, Nsq) at ' latstrshort ', ' lonstr])
ylabel({'NTsq [s^{-2}]';'depth [m]'})
set(gca,'tickdir','out')

ax(2) = subplot(512);
pcolor(tdy.time,tdy.depth,tdy.NSsq)
shading flat
axis ij
colormap(gca,flipud(lbmap(24,'RedBlue')))
datetick
colorbar
caxis([-1 1]*1e-3)
ylabel({'NSsq [s^{-2}]';'depth [m]'})
set(gca,'tickdir','out')


ax(3) = subplot(513);
pcolor(tdy.time,tdy.depth,tdy.Nsq)
shading flat
axis ij
colormap(gca,flipud(lbmap(24,'RedBlue')))
datetick
colorbar
caxis([-1 1]*1e-3)
ylabel({'N^2 [s^{-2}]';'depth [m]'})
set(gca,'tickdir','out')

ax(4) = subplot(514);
pcolor(tdy.time,tdy.depth,tdy.NTsq+tdy.NSsq)
shading flat
axis ij
colormap(gca,flipud(lbmap(24,'RedBlue')))
datetick
colorbar
caxis([-1 1]*1e-3)
ylabel({'N_T^2+N_S^2 [s^{-2}]';'depth [m]'})
set(gca,'tickdir','out')

ax(5) = subplot(515);
pcolor(tdy.time,tdy.depth,(tdy.Nsq - (tdy.NTsq+tdy.NSsq))./tdy.Nsq*100)
shading flat
axis ij
colormap(gca,flipud(lbmap(24,'RedBlue')))
datetick
colorbar
caxis([-1 1]*10)
ylabel({'error [%]';'(N^2-(N_T^2+N_S^2))/N^2*100';'depth [m]'})
set(gca,'tickdir','out')

linkaxes

export_fig([figdir 'Nsq_compare.png'],'-r200')




%%%%%%%%% plot mld and z20 %%%%%%%%%
clear ax

figure(349)
clf
set(gcf,'position',[512          57        1117        1048])


pcolor(tdy.time,tdy.depth,tdy.T)
shading flat
axis ij
hold on
plot(tdy.time,tdy.depth20,'k')
plot(tdy.time,tdy.mld,'k')
colormap(gca,flipud(lbmap(24,'RedBlue')))
datetick
colorbar
title(['Temperature, MLD & z_{20} at ' latstrshort ', ' lonstr])
ylabel({'temp [C]';'depth [m]'})
set(gca,'tickdir','out')

export_fig([figdir 'T_mld_z20.png'],'-r200')





%%%%%%%%% plot mld at different time-spacing %%%%%%%%%

dd = dir([savedir 'mld*']);

for ii = 1:length(dd)
    load([savedir dd(ii).name])
end

varsrm = fields(mld);
vars = varsrm(2:end-1);


figure(265)
clf
    set(gcf,'position',[495          57         766        1048])


for ii = 1:length(vars)
    
    ax(ii) = subplot(length(vars),1,ii);
    
    plot(mld.time,mld.(vars{ii}))
    hold on
    plot(mldhr.time,mldhr.(vars{ii}))
    plot(mlddy.time,mlddy.(vars{ii}))
    plot(mldmo.time,mldmo.(vars{ii}))
    
    datetick
    
    ylabel(vars{ii})
    set(gca,'tickdir','out')
    
    if ii == 1
        title(['Mixed Layer depth at ' latstrshort ', ' lonstr])
        legend('10min','hr','daily','monthly','location','northwest')
        legend boxoff
    end
end

export_fig([figdir 'mld.png'],'-r200')
