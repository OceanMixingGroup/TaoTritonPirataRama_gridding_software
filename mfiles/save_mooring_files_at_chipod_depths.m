% save_mooring_files_at_chipod_depths.m
%
% Use the data from this mooring to create input files that can be read by
% the chipod_gust processing software.
%
% These files contain the following structures and variables:
%
%   PROCDIR/input/dTdz_m.mat
%       Tz_m.time
%       Tz_m.Tz
%       Tz_m.N2
%
%   PROCDIR/input/vel_m.mat
%       vel_m.u
%       vel_m.v
%       vel_m.spd
%       vel_m.time
%       vel_m.depth
%       vel_m.U
%       vel_m.comment
%
% Note, you will have to move the files manually to the right location in
% the input folder in the chipod_gust software.

clear
addpath('./utilities/')

%% define location
% this creates correct strings and folders

%%%%%%%%% be sure to set the correct location in set_location.m! %%%%%%%%%%

set_location;

%% inputs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEPLOYMENT NAME
% note, the deployment name is used for picking the start and end times of 
% the deployment and saving and naming the files. The location where the 
% data comes from is determined by set_location.

% dplname = 'Pirata14';
% dplname = 'Pirata15';
% dplname = 'Pirata16';
% dplname = 'Pirata17'; 
dplname = 'Pirata18';  
% dplname = 'tao16_140';
% dplname = 'tao17_140';  
% dplname = 'tao17_110';
% dplname = 'tao18_140';  
% dplname = 'tao18_110';
% dplname = 'tao18_125';
% dplname = 'tao19_140';
% dplname = 'RAMA18_67E';
% dplname = 'RAMA18_80_5E';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEPTHS OF CHIPODS (AND ADDITIONAL USEFUL VELOCITY DEPTHS)
if strcmp(dplname(1:3),'Pir')
    chipod_depths = [21 25 30 35 50 65 81];
elseif strcmp(dplname(end-2:end),'140')
    chipod_depths = [29 30 35 49 69 89 119];    % add in 30m and 35m for velocity of 29m chipod
elseif strcmp(dplname(end-2:end),'110')
    chipod_depths = [9 29 30 35 49];    % add in 30m and 35m for velocity of 29m chipod
elseif strcmp(dplname(end-2:end),'125')
    chipod_depths = [30 50 70];    % no velocity at this location
elseif strcmp(dplname(1:3),'RAM')
    chipod_depths = [15 30 45];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START AND END TIMES OF THE DEPLOYMENT
% (if left blank, the entire timeseries from TaoTritonPirataRama will be 
% used to make the mooring files)
% tstart = datenum(0,0,0,0,0,0);
% tend = datenum(0,0,0,0,0,0);

if strcmp(dplname,'Pirata14')
    tstart = datenum(2014,4,1); % date rande is long enough for both Pirata14_23 and Pirata14_10
    tend = datenum(2015,4,1);
elseif strcmp(dplname,'Pirata15')
    tstart = datenum(2015,3,1);
    tend = datenum(2016,4,1);
elseif strcmp(dplname,'Pirata16')
    tstart = datenum(2016,3,1);
    tend = datenum(2017,4,1);
elseif strcmp(dplname,'Pirata17')
    tstart = datenum(2017,3,1);
    tend = datenum(2018,4,1);
elseif strcmp(dplname,'Pirata18')
    tstart = datenum(2018,3,1);
    tend = datenum(2019,4,1);
elseif strcmp(dplname,'tao15_140')
    tstart = datenum(2015, 03, 23, 03, 00, 00);
    tend = datenum(2016, 02, 17, 11, 00, 00);
elseif strcmp(dplname,'tao16_140')  % all chipods die within about a week of deployment, but includeing a bit more data in mooring files
    tstart = datenum(2016, 2, 20, 00, 00, 00);
    tend = datenum(2016, 4, 01, 00, 00, 00);
elseif strcmp(dplname,'tao17_140')
    tstart = datenum(2017, 1, 25, 00, 00, 00);
    tend = datenum(2018, 09, 17, 00, 00, 00);
elseif strcmp(dplname,'tao17_110')
    tstart = datenum(2017, 6, 1, 00, 00, 00);
    tend = datenum(2017, 7, 1, 00, 00, 00);
elseif strcmp(dplname,'tao18_140')
    tstart = datenum(2018, 09, 17, 0, 00, 00);
    tend = datenum(2019, 10, 09, 00, 00, 00);
elseif strcmp(dplname,'tao18_110')
    tstart = datenum(2018, 04, 01, 0, 00, 00);
    tend = datenum(2019, 07, 31, 00, 00, 00);
elseif strcmp(dplname,'tao18_125')
    tstart = datenum(2018, 09, 30, 0, 00, 00);
    tend = datenum(2019, 10, 25, 00, 00, 00);
elseif strcmp(dplname,'tao19_140')
    tstart = datenum(2019, 10, 01, 0, 00, 00);
    tend = datenum(2020, 04, 01, 00, 00, 00);
elseif strcmp(dplname,'RAMA18_67E')
    tstart = datenum(2018,06,28);
    tend = datenum(2019,10,31);
elseif strcmp(dplname,'RAMA18_80_5E')
    tstart = datenum(2018,11,03);
    tend = datenum(2019,08,19);
end


% Do you want to process just temperature or velocity too?
% For some deployments, no velocity data was collected (such as tao18_125),
% so we need to skip over all the velocity portions of the code.
if lon == (360-125)
    do_vel = 0;
else
    do_vel = 1;
end
	

%% load in processed data files

disp('loading temperature file...')
load([savedir 'T_' latstrshort '_' lonstr '_10min.mat'])

if do_vel
    disp('loading velocity file...')
    load([savedir 'vel_' latstrshort '_' lonstr '_10min.mat'])
end

%% loop through chipod depths and interpolate, save and plot

depth = t.depth;
time = t.time;

% find indices of time limits
if tstart ~= 0
    indt1 = find(time >= tstart,1,'first');
    indt2 = find(time <= tend,1,'last');
else
    indt1 = 1;
    indt2 = length(time);
end
indt = indt1:indt2;

% loop through all chipod depths
for jj = 1:length(chipod_depths)
    
    clear indd Tz_m vel_m indda inddb za dratio
    
    % interpolate data to the chipod depth    
    indd = find(depth == chipod_depths(jj));
    
    if ~isempty(indd)
        Tz_m.time = time(indt);
        Tz_m.T = t.T(indd,indt);
        Tz_m.Tz = t.dTdz(indd,indt);
        Tz_m.N2 = t.NTsq(indd,indt);
        Tz_m.depth = chipod_depths(jj);
        
        if do_vel
            vel_m.time = time(indt);
            vel_m.u = vel.u(indd,indt);
            vel_m.v = vel.v(indd,indt);
            vel_m.spd = (vel_m.u.^2 + vel_m.v.^2).^(0.5);
            vel_m.U = vel_m.u + vel_m.v*1i;
            vel_m.depth = chipod_depths(jj);
        end

    else
        % linearly interpolate between vertical grid points to chipod
        % depths
        indda = find(depth < chipod_depths(jj),1,'last');
        inddb = find(depth > chipod_depths(jj),1,'first');
        
        za = depth(indda);
        zb = depth(inddb);        
        dratio = (chipod_depths(jj)-za)/(zb-za);
        
        Tz_m.time = time(indt);
        Tz_m.T    = t.T(inddb,indt)*dratio + t.T(indda,indt)*(1-dratio);
        Tz_m.Tz = t.dTdz(inddb,indt)*dratio + t.dTdz(indda,indt)*(1-dratio);
        Tz_m.N2   = t.NTsq(inddb,indt)*dratio + t.NTsq(indda,indt)*(1-dratio);
        Tz_m.depth = chipod_depths(jj);
        
        if do_vel
            vel_m.time = time(indt);
            vel_m.u = vel.u(inddb,indt)*dratio + vel.u(indda,indt)*(1-dratio);
            vel_m.v = vel.v(inddb,indt)*dratio + vel.v(indda,indt)*(1-dratio);
            vel_m.spd = (vel_m.u.^2 + vel_m.v.^2).^(0.5);
            vel_m.U = vel_m.u + vel_m.v*1i;
            vel_m.depth = chipod_depths(jj);
        end

    end
    
    Tz_m.readme = strvcat('Temperature mooring file for chipod to use for chipod_gust processing:',...
        ['location: ' latstrshort ', ' lonstr],...
        ['depth: ' num2str(chipod_depths(jj))],...
        ' ',...
        'Tz_m.time = matlab datenum',...
        'Tz_m.T = temperature [°C] profile calculated at every time step ...',...
        '    using a cubic Akima spline to combine irregularly spaced ...',...
        '    temperature sensors on mooring into dz = 5m grid ...',...
        '    (following appendix of Pham et al. 2017). Temperature on ...',...
        '    dz = 5m grid is then interpolated to chipod depths.',...
        'Tz_m.Tz = vertical gradient of temperature profile (dz = 5m)...',...
        '    interpolated to chipod depths.',...
        'Tz_m.N2 = buoyancy frequency [s^{-2}] calculated using ...',...
        '    N_T^2 = -g*alpha*dT/dz (i.e. only temperature is used, not ...',...
        '    salinity) using the temperature gradient profile described above.',...
        'Tz_m.depth = depth of chipod for which this data corresponds',...        
        ' ',...
        'Mooring file created with TaoTritonPirataRama data by',...
        ['~/ganges/data/TaoTritonPirataRama_processed_at_' latstrshort '_' lonstr '/...'],...
        '    mfiles/save_mooring_files_at_chipod_depths.m');
        
    if do_vel
        vel_m.readme = strvcat('Velocity mooring file for chipod to use for chipod_gust processing:',...
            ['location: ' latstrshort ', ' lonstr],...
            ['depth: ' num2str(chipod_depths(jj))],...
            ' ',...
            'vel_m.time = matlab datenum',...
            'vel_m.u and vel_m.v = zonal and meridional velocities from ADCP ...',...
            '    and/or current meters on dz = 5m grid. If current meters are ...',...
            '    available in upper part of the water column, a high-order ...',...
            '    scheme is used to combine ADCP, current meters, and surface ...',...
            '    wind speeds (following Pham et al. 2017 Appendix B). ...',...
            'vel_m.spd = sqrt(u.^2 + v.^2)',...
            'vel_m.U = u + iv',...
            'vel_m.N2 = buoyancy frequency [s^{-2}] calculated using ...',...
            '    N_T^2 = -g*alpha*dT/dz (i.e. only temperature is used, not ...',...
            '    salinity) using the temperature gradient profile described above.',...
            'vel_m.depth = depth of chipod for which this data corresponds',...        
            ' ',...
            'Mooring file created with TaoTritonPirataRama data by',...
            ['~/ganges/data/TaoTritonPirataRama_processed_at_' latstrshort '_' lonstr '/...'],...
            '    mfiles/save_mooring_files_at_chipod_depths.m');    
    end
    
    
    % save files
    disp('saving...')
    mkdir([savedir 'mooring_files_for_chipod_gust_processing_' dplname])
    save([savedir 'mooring_files_for_chipod_gust_processing_' dplname '/dTdz_m_' ...
        latstrshort '_' lonstr '_' num2str(chipod_depths(jj)) 'm.mat'],'Tz_m')
    if do_vel
        save([savedir 'mooring_files_for_chipod_gust_processing_' dplname '/vel_m_' ...
            latstrshort '_' lonstr '_' num2str(chipod_depths(jj)) 'm.mat'],'vel_m')
    end
    


    % plotting
    disp('plotting...')
    
    varsT = {'T';'Tz';'N2'};
    labsT = {'[°C]';'[°C m^{-1}]';'(N_T^2) [s^{-2}]'};

    varsu = {'u';'v';'spd'};
    labsu = {'[m s^{-1}]';'[m s^{-1}]';'[m s^{-1}]'};
    
    figure(567)
    clf
    set(gcf,'position',[229         364        1303         591])
    
    
    for ii = 1:3        
        ax(ii) = subplot(3,2,ii*2-1);        
        plot(Tz_m.time,Tz_m.(varsT{ii}),'k')
        ylabel({varsT{ii};labsT{ii}})
        set(gca,'tickdir','out')
        xlim([tstart tend])
        datetick('x','mmmyy','keeplimits')        
        if ii == 1
            title([latstrshort ', ' lonstr '   |   depth = ' num2str(chipod_depths(jj))  'm'])
        end
    end
    
    if do_vel
        for ii = 1:3        
            ax(ii) = subplot(3,2,ii*2);        
            plot(vel_m.time,vel_m.(varsu{ii}),'k')
            ylabel({varsu{ii};labsu{ii}})
            set(gca,'tickdir','out')
            xlim([tstart tend])
            datetick('x','mmmyy','keeplimits')        
            if ii == 1
                title([latstrshort ', ' lonstr '   |   depth = ' num2str(chipod_depths(jj))  'm'])
            end
        end
    end

    warning off
    linkaxes(ax,'x')
    warning on
        
    mkdir([figdir 'mooring_files_for_chipod_gust_processing_' dplname])
    export_fig([figdir 'mooring_files_for_chipod_gust_processing_' dplname '/mooring_input_files_' ...
        latstrshort '_' lonstr '_' num2str(chipod_depths(jj)) 'm.png'],'-r200')
    savefig([figdir 'mooring_files_for_chipod_gust_processing_' dplname '/mooring_input_files_' ...
        latstrshort '_' lonstr '_' num2str(chipod_depths(jj)) 'm.fig'])
    
end




