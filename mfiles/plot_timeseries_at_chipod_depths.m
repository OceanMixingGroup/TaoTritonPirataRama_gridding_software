% plot_timeseries_at_chipod_depths.m
% 
% plot a timeseries the following variables at chipod depths:
%   T
%   dT/dz
%   Nsq
%   u
%   v
%   du/dz
%   dv/dz
%   Ssq

clear

%% inputs

% longitudinal location of the mooring
loc = '10W';

% depths of chipods
chipod_depths = [21 35 50 65 81];

% start and end times of deployment
% (if left blank, the entire timeseries will be plotted)
tstart = datenum(2017,1,1,0,0,0);
tend = datenum(2020,1,1,0,0,0);

% temporal resolution of the processed data
% (4 resolutions are available: monthly, daily, hourly, and 10 minute (''). 
% Note that data collected at a slower frequency are interpolated the given 
% time grid and data collected at a faster frequency are bin-averaged to 
% the given time grid.
% frq = 'monthly';     % plot monthly averaged data
% frq = 'daily';     % plot daily averaged data
frq = 'hourly';     % plot hourly data
% frq = '10min';       % plot 10 minute data


%% load in processed data files

disp('loading temperature file...')
load(['../processed/T_0_' loc '_' frq '.mat'])

disp('loading velocity file...')
load(['../processed/vel_0_' loc '_' frq '.mat'])

if strcmp(frq(1),'h')
    vel = velhr;
    t = thr;
    clear velhr thr
elseif strcmp(frq(1),'m')
    vel = velmo;
    t = tmo;
    clear velmo tmo
elseif strcmp(frq(1),'d')
    vel = veldy;
    t = tdy;
    clear veldy tdy
end
    

%% plot variables all variables at each depth
% (one plot for each depth)

vars = {'T';'dTdz';'NTsq';'u';'v';'dudz';'dvdz';'Ssq'};
labs = {'[°C]';'[°C m^{-1}]';'[s^{-2}]';'[m s^{-1}]';'[m s^{-1}]';...
    '[s^{-1}]';'[s^{-1}]';'[s^{-2}]'};

depth = t.depth;
time = t.time;

for jj = 1:length(chipod_depths)
    
    clear indd T dTdz NTsq u v dudz dvdz Ssq

    figure(743)
    clf
    set(gcf,'position',[376          57        1157        1048])
    
    % interpolate to the correct depth    
    indd = find(depth == chipod_depths(jj));
    
    if ~isempty(indd)
        A.T = t.T(indd,:);
        A.dTdz = t.dTdz(indd,:);
        A.NTsq = t.NTsq(indd,:);
        A.u = vel.u(indd,:);
        A.v = vel.v(indd,:);
        A.dudz = vel.dudz(indd,:);
        A.dvdz = vel.dvdz(indd,:);
        A.Ssq = vel.Ssq(indd,:);
    else
        indda = find(depth < chipod_depths(jj),1,'last');
        inddb = find(depth > chipod_depths(jj),1,'first');
        
        za = depth(indda);
        zb = depth(inddb);        
        dratio = (chipod_depths(jj)-za)/(zb-za);
        
        A.T    = t.T(inddb,:)*dratio + t.T(indda,:)*(1-dratio);
        A.dTdz = t.dTdz(inddb,:)*dratio + t.dTdz(indda,:)*(1-dratio);
        A.NTsq = t.NTsq(inddb,:)*dratio + t.NTsq(indda,:)*(1-dratio);
        A.u    = vel.u(inddb,:)*dratio + vel.u(indda,:)*(1-dratio);
        A.v    = vel.v(inddb,:)*dratio + vel.v(indda,:)*(1-dratio);
        A.dudz = vel.dudz(inddb,:)*dratio + vel.dudz(indda,:)*(1-dratio);
        A.dvdz = vel.dvdz(inddb,:)*dratio + vel.dvdz(indda,:)*(1-dratio);
        A.Ssq = vel.Ssq(inddb,:)*dratio + vel.Ssq(indda,:)*(1-dratio);
    end
    

    for ii = 1:length(vars)
        
        ax(ii) = subplot(4,2,ii);
        
        plot(time,A.(vars{ii}),'k')
        ylabel({vars{ii};labs{ii}})
        set(gca,'tickdir','out')
        
        if tstart ~= 0
            xlim([tstart tend])
        else
            xlim(time([1 end]))
        end
        datetick('x','keeplimits')
        
        if ii == 1
            title(['0, ' loc '   |   depth = ' num2str(chipod_depths(jj)) ...
                'm   |   ' frq ' data'])
        end
        
        warning off
        linkaxes(ax,'x')
        warning on
    end
    
    export_fig(['../figures/timeseries_' loc '_at_chipod_depth_' ...
        num2str(chipod_depths(jj)) '.png'],'-r200')
    savefig(['../figures/timeseries_' loc '_at_chipod_depth_' ...
        num2str(chipod_depths(jj)) '.fig'])

end



