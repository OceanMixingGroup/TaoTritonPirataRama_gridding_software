% load_and_save_tao_adcp_cur
% 
% In order to preserve some continuity, this code is based on:
% ~/Dropbox/work/RESEARCH/TAO140/mfiles/velocity_0_140W.m
%
% the goal of this code is to load in the RAW adcp and current meter records
% from the tao array and make a gridded dataset that contains velcity and shear.
%
% changes from old version:
%   - loads in netcdf files from TAO rather than .mat files
%   - interpolates to the standard time grids (10 min, hourly, daily, monthly)


clear
addpath('./utilities/')
addpath('./adcp_info/')


%% define location
% this creates correct strings and folders

%%%%%%%%% be sure to set the correct location in set_location.m! %%%%%%%%%%

set_location;

%% load common tao time

time = taotime;

%% load common tao depth grid

[depth,dz] = taodepth;

%% load in high res adcp data

disp('loading and interpolating raw adcp data from TaoTritonPirataRama...')

%%%%% adcp from TaoTritonPirataRama netcdf files %%%%%
filename = ['../../TaoTritonPirataRama/high_resolution/hr/adcp' ...
    lower(latstr) lower(lonstr) '_hr.cdf'];
if exist(filename,'file')
    raw = loadnc(filename);
    if nanmean(nanmean(raw.u)) > 5
        raw.u = raw.u/100;  % convert cm/s to m/s
        raw.v = raw.v/100;
    end
else
    disp('WARNING: load_and_save_tao_adcp_cur cannot find raw ADCP file')
    disp('in the TaoTritonPirataRama databank.')
    disp(['looking for: ' filename])
    disp(' ')
    raw.time = time;
    raw.depth = 0;
    raw.u = NaN*time;
    raw.v = NaN*time;
end

%%%%% interpolate to common depth and time grid %%%%%
% note, this is much faster and more accurate to do this interpolation
% twice using interp1 rather than using griddata

% predefine adcp structure
adcp.time = time;
adcp.depth = depth;
adcp.u = NaN*ones(length(adcp.depth),length(adcp.time));
adcp.v = adcp.u;

% first put the raw data on an even depth grid
disp('interpolating to correct depth grid...')

% first try a fast way
if raw.depth ~= 0
    [i,j] = size(raw.depth);
    if (i==1 | j==1) & nanmean(diff(raw.depth)) == dz
        cc = 1;
        for ii = 1:length(raw.depth)
            clear inda
            inda = find(adcp.depth == raw.depth(ii));
            if ~isempty(inda)
                inddadcp(cc) = inda;
                cc = cc + 1;
            end
        end
        cc = 1;
        for ii = 1:length(adcp.depth)
            clear inda
            inda = find(raw.depth == adcp.depth(ii));
            if ~isempty(inda)
                inddraw(cc) = inda;
                cc = cc + 1;
            end
        end

        raw.depthinterp = depth;
        raw.uinterpz = NaN*ones(length(raw.depthinterp),length(raw.time));
        raw.vinterpz = raw.uinterpz;

        raw.uinterpz(inddadcp,:) = raw.u(inddraw,:);
        raw.vinterpz(inddadcp,:) = raw.v(inddraw,:);

        clear inda cc indadcp indraw ii 
    else
        % damn, fast way doesn't work. Interpolate in the slow way...
        for ii = 1:length(raw.time)
            if mod(ii,1000) == 0
                disp([num2str(ii) '/' num2str(length(raw.time))])
            end
            raw.uinterpz(ii,:) = interp1(raw.depth,raw.u(:,ii),adcp.depth);
            raw.vinterpz(ii,:) = interp1(raw.depth,raw.v(:,ii),adcp.depth);
        end
    end

    % loop through all depths to interpolate to common time grid
    disp('interpolating to correct time grid...')
    for ii = 1:length(adcp.depth)
        adcp.u(ii,:) = interp1(raw.time,raw.uinterpz(ii,:),adcp.time);
        adcp.v(ii,:) = interp1(raw.time,raw.vinterpz(ii,:),adcp.time);
    end

end


%% load in daily ADCP data
% this only pertains to old ADCP data before the mid-90s

if time(1) < datenum(2005,1,1)

    gooddailydata = 1;

    filename = ['../../TaoTritonPirataRama/' mooringarray '/adcp_xyzt_dy.cdf'];
    if exist(filename,'file')
        dummy = loadnc(filename);
    else
        gooddailydata = 0;
        disp('WARNING: load_and_save_tao_adcp_cur could not find the')
        disp(['daily ADCP file: ../../TaoTritonPirataRama/' mooringarray '/adcp_xyzt_dy.cdf'])
        disp('No daily data will be used to fill gaps in the high resolution data.')        
        disp(' ')
    end

    if gooddailydata
        indlat = find(dummy.lat == lat);
        indlon = find(dummy.lon == lon);

        if isempty(indlat) | isempty(indlon)
            gooddailydata = 0;
            disp(['WARNING: load_and_save_tao_sadcp_cur cannot find the desired lat and lon'])
            disp(['in the daily ../../TaoTritonPirataRama/' mooringarray '/adcp_xyzt_dy.cdf. '])
            disp('No daily data will be used to fill gaps in the high resolution data')
            disp(' ')
        else
            rawddy.time = dummy.time;
            rawddy.depth = dummy.depth;
            rawddy.lat = dummy.lat(indlat);
            rawddy.lon = dummy.lon(indlon);
            rawddy.u = squeeze(dummy.U(:,:,indlat,indlon))'/100;  % convert cm/s to m/s
            rawddy.v = squeeze(dummy.V(:,:,indlat,indlon))'/100;  % convert cm/s to m/s

            % interpolate to 10 minute time grid
            rawdy.time = time;
            rawdy.depth = rawddy.depth;
            rawdy.lat = rawddy.lat;
            rawdy.lon = rawddy.lon;
            for ii = 1:length(rawdy.depth)
                rawdy.u(ii,:) = interp1(rawddy.time,rawddy.u(ii,:),time);
                rawdy.v(ii,:) = interp1(rawddy.time,rawddy.v(ii,:),time);
            end
        end
    end
    
    % fill holes in 10min data with daily data
    if gooddailydata

        disp('filling holes in 10min ADCP data with daily data...')

        for ii = 1:length(adcp.depth)
            clear inddepth
            inddepth = find(rawdy.depth == adcp.depth(ii));
            if ~isempty(inddepth)
                if adcp.depth(ii) == rawdy.depth(inddepth)
                    indbad = isnan(adcp.u(ii,:));
                    adcp.u(ii,indbad) = rawdy.u(inddepth,indbad);
                    adcp.v(ii,indbad) = rawdy.v(inddepth,indbad);
                end            
            end    
        end    
    end 
end



%% load in high res current meter data

yncur = 0; % this flag will change if current meter data is found

%%%%% load current meter data from HOURLY TaoTritonPirataRama netcdf files %%%%%
filename = ['../../TaoTritonPirataRama/high_resolution/hr/cur' ...
    lower(latstr) lower(lonstr) '_hr.cdf'];
if exist(filename,'file')
    yncur = 1;
    disp('loading and interpolating raw hourly current meter data from TaoTritonPirataRama...')
    raw = loadnc(filename);
    if isfield(raw,'U')
        raw.u = raw.U;
        raw.v = raw.V;
    end
    if nanmean(nanmean(raw.u)) > 5
        raw.u = raw.u/100;  % convert cm/s to m/s
        raw.v = raw.v/100;
    end
else
    raw.time = time;
    raw.depth = 0;
    raw.u = NaN*time;
    raw.v = NaN*time;
end

%%%%% interpolate to common time grid %%%%%

cur.time = time;
cur.depth = raw.depth;
if length(cur.depth) == 1
    cur.u(1,:) = interp1(raw.time,raw.u(:),cur.time);
    cur.v(1,:) = interp1(raw.time,raw.v(:),cur.time);
else
    for ii = 1:length(cur.depth)
        cur.u(ii,:) = interp1(raw.time,raw.u(ii,:),cur.time);
        cur.v(ii,:) = interp1(raw.time,raw.v(ii,:),cur.time);
    end
end


%%%%% load current meter data from 20min TaoTritonPirataRama netcdf files %%%%%
filename = ['../../TaoTritonPirataRama/high_resolution/20m/cur' ...
    lower(latstr) lower(lonstr) '_20m.cdf'];
if exist(filename,'file')
    yncur = 1;
    disp('loading and interpolating raw 20 min current meter data from TaoTritonPirataRama...')
    raw = loadnc(filename);
    if isfield(raw,'U')
        raw.u = raw.U;
        raw.v = raw.V;
    end
    if nanmean(nanmean(raw.u)) > 5
        raw.u = raw.u/100;  % convert cm/s to m/s
        raw.v = raw.v/100;
    end
else
    raw.time = time;
    raw.depth = 0;
    raw.u = NaN*time;
    raw.v = NaN*time;
end

%%%%% interpolate to common time grid %%%%%
cur20.time = time;
cur20.depth = raw.depth;
if length(cur20.depth) == 1
    cur20.u(1,:) = interp1(raw.time,raw.u(:),cur.time);
    cur20.v(1,:) = interp1(raw.time,raw.v(:),cur.time);
else
    for ii = 1:length(cur20.depth)
        cur20.u(ii,:) = interp1(raw.time,raw.u(ii,:),cur.time);
        cur20.v(ii,:) = interp1(raw.time,raw.v(ii,:),cur.time);
    end
end



%%%%% load current meter data from 30min TaoTritonPirataRama netcdf files %%%%%
filename = ['../../TaoTritonPirataRama/high_resolution/30m/cur' ...
    lower(latstr) lower(lonstr) '_30m.cdf'];
if exist(filename,'file')
    yncur = 1;
    disp('loading and interpolating raw 30 min current meter data from TaoTritonPirataRama...')
    raw = loadnc(filename);
    if isfield(raw,'U')
        raw.u = raw.U;
        raw.v = raw.V;
    end
    if nanmean(nanmean(raw.u)) > 5
        raw.u = raw.u/100;  % convert cm/s to m/s
        raw.v = raw.v/100;
    end
else
    raw.time = time;
    raw.depth = 0;
    raw.u = NaN*time;
    raw.v = NaN*time;
end

%%%%% interpolate to common time grid %%%%%
cur30.time = time;
cur30.depth = raw.depth;
if length(cur30.depth) == 1
    cur30.u(1,:) = interp1(raw.time,raw.u(:),cur.time);
    cur30.v(1,:) = interp1(raw.time,raw.v(:),cur.time);
else
    for ii = 1:length(cur30.depth)
        cur30.u(ii,:) = interp1(raw.time,raw.u(ii,:),cur.time);
        cur30.v(ii,:) = interp1(raw.time,raw.v(ii,:),cur.time);
    end
end




%%%%%%%% combine hourly and 20 minute current meter records %%%%%%%%
if cur20.depth ~= 0
    for ii = 1:length(cur.depth)
        clear indd
        indd = find(cur20.depth == cur.depth(ii));
        if ~isempty(indd)
            clear nanmask
            nanmask = isnan(cur.u(ii,:));
            cur.u(ii,nanmask) = cur20.u(indd,nanmask);
            clear nanmask
            nanmask = isnan(cur.v(ii,:));
            cur.v(ii,nanmask) = cur20.v(indd,nanmask);
        end
    end
end

% add cur20 at depths that don't already exist in cur.depth
if cur20.depth ~= 0
    for ii = 1:length(cur20.depth)
        clear indd depthorig uorig vorig
        indd = find(cur.depth == cur20.depth(ii));
        if ~isempty(indd)
            % this depth should have already been added to cur
        else
            depthorig = cur.depth;
            uorig = cur.u;
            vorig = cur.v;

            indd = find(cur.depth < cur20.depth(ii),1,'last');
            if indd ~= length(cur.depth)
                cur.depth = [depthorig(1:indd); cur20.depth(ii); ...
                    depthorig(indd+1:end)];
                cur.u = [cur.u(1:indd,:); cur20.u(ii,:); cur.u(indd+1:end,:)];
                cur.v = [cur.v(1:indd,:); cur20.v(ii,:); cur.v(indd+1:end,:)];
            else
                cur.depth = [depthorig(1:indd); cur20.depth(ii)];
                cur.u = [cur.u(1:indd,:); cur20.u(ii,:)];
                cur.v = [cur.v(1:indd,:); cur20.v(ii,:)];
            end
        end
    end
end


%%%%%%%% combine hourly and 30 minute current meter records %%%%%%%%
if cur30.depth ~= 0
    for ii = 1:length(cur.depth)
        clear indd
        indd = find(cur30.depth == cur.depth(ii));
        if ~isempty(indd)
            clear nanmask
            nanmask = isnan(cur.u(ii,:));
            cur.u(ii,nanmask) = cur30.u(indd,nanmask);
            clear nanmask
            nanmask = isnan(cur.v(ii,:));
            cur.v(ii,nanmask) = cur30.v(indd,nanmask);
        end
    end
end

% add cur30 at depths that don't already exist in cur.depth
if cur30.depth ~= 0
    for ii = 1:length(cur30.depth)
        clear indd depthorig uorig vorig
        indd = find(cur.depth == cur30.depth(ii));
        if ~isempty(indd)
            % this depth should have already been added to cur
        else
            depthorig = cur.depth;
            uorig = cur.u;
            vorig = cur.v;

            indd = find(cur.depth < cur30.depth(ii),1,'last');
            if indd ~= length(cur.depth)
                cur.depth = [depthorig(1:indd); cur30.depth(ii); ...
                    depthorig(indd+1:end)];
                cur.u = [cur.u(1:indd,:); cur30.u(ii,:); cur.u(indd+1:end,:)];
                cur.v = [cur.v(1:indd,:); cur30.v(ii,:); cur.v(indd+1:end,:)];
            else
                cur.depth = [depthorig(1:indd); cur30.depth(ii)];
                cur.u = [cur.u(1:indd,:); cur30.u(ii,:)];
                cur.v = [cur.v(1:indd,:); cur30.v(ii,:)];
            end
        end
    end
end


%% load in daily Current Meter data
% some of the RAMA locations have daily current meter data that is not
% saved at any high-resolution time scales


gooddailycurrents = 1;

filename = ['../../TaoTritonPirataRama/' mooringarray '/cur_xyzt_dy.cdf'];
if exist(filename,'file')
    dummy = loadnc(filename);
else
    gooddailycurrents = 0;
    disp('WARNING: load_and_save_tao_adcp_cur could not find the')
    disp(['daily ADCP file: ../../TaoTritonPirataRama/' mooringarray '/adcp_xyzt_dy.cdf'])
    disp('No daily data will be used to fill gaps in the high resolution data.')        
    disp(' ')
end

if gooddailycurrents
    indlat = find(dummy.lat == lat);
    indlon = find(dummy.lon == lon);

    if isempty(indlat) | isempty(indlon)
        gooddailycurrents = 0;
        disp(['WARNING: load_and_save_tao_sadcp_cur cannot find the desired lat and lon'])
        disp(['in the daily ../../TaoTritonPirataRama/' mooringarray '/adcp_xyzt_dy.cdf. '])
        disp('No daily data will be used to fill gaps in the high resolution data')
        disp(' ')
    else
        rawddy.time = dummy.time;
        rawddy.depth = dummy.depth;
        rawddy.lat = dummy.lat(indlat);
        rawddy.lon = dummy.lon(indlon);
        rawddy.u = squeeze(dummy.U(:,:,indlat,indlon))'/100;  % convert cm/s to m/s
        rawddy.v = squeeze(dummy.V(:,:,indlat,indlon))'/100;  % convert cm/s to m/s

        % interpolate to 10 minute time grid
        rawdy.time = time;
        rawdy.depth = rawddy.depth;
        rawdy.lat = rawddy.lat;
        rawdy.lon = rawddy.lon;
        for ii = 1:length(rawdy.depth)
            rawdy.u(ii,:) = interp1(rawddy.time,rawddy.u(ii,:),time);
            rawdy.v(ii,:) = interp1(rawddy.time,rawddy.v(ii,:),time);
        end
    end
end

% fill holes in high-resolution data with daily data
if gooddailycurrents

    disp('filling holes in high-resolution current meter data with daily data...')

    for ii = 1:length(cur.depth)
        clear inddepth
        inddepth = find(rawdy.depth == cur.depth(ii));
        if ~isempty(inddepth)
            clear nanmask
            nanmask = isnan(cur.u(ii,:));
            cur.u(ii,nanmask) = rawdy.u(inddepth,nanmask);
            clear nanmask
            nanmask = isnan(cur.v(ii,:));
            cur.v(ii,nanmask) = rawdy.v(inddepth,nanmask);           
        end    
    end  
    
    for ii = 1:length(rawdy.depth)
        clear indd depthorig uorig vorig
        indd = find(cur.depth == rawdy.depth(ii));
        if ~isempty(indd)
            % this depth should have already been added to cur
        else
            depthorig = cur.depth;
            uorig = cur.u;
            vorig = cur.v;

            indd = find(cur.depth < rawdy.depth(ii),1,'last');
            if indd ~= length(cur.depth)
                cur.depth = [depthorig(1:indd); rawdy.depth(ii); ...
                    depthorig(indd+1:end)];
                cur.u = [cur.u(1:indd,:); rawdy.u(ii,:); cur.u(indd+1:end,:)];
                cur.v = [cur.v(1:indd,:); rawdy.v(ii,:); cur.v(indd+1:end,:)];
            else
                cur.depth = [depthorig(1:indd); rawdy.depth(ii)];
                cur.u = [cur.u(1:indd,:); rawdy.u(ii,:)];
                cur.v = [cur.v(1:indd,:); rawdy.v(ii,:)];
            end
        end
    end
end 




%% load in additonal ADCP and current meter data from random files if they exist

%%%%% parse in additional ADCP files %%%%%
ynadditionalfiles = exist(['./adcp_info/adcp_info_at_' ...
    latstrshortnop '_' lonstrnop '.m'],'file');

if ynadditionalfiles == 2
    % run the function that describes extra adcp and current meter files
    run(['./adcp_info/adcp_info_at_' latstrshortnop '_' lonstrnop '.m']);
end


%% parse additional ADCP data (if it exists) into adcp structure 

if ynadditionalfiles == 2 & exist('ca','var')

    for ii = 1:adcpinfo.numfiles
        % parse new data into adcp.u and adcp.v in regions where adcp.u and
        % adcp.v are NaN
        nanmask = isnan(adcp.u);
        adcp.u(nanmask) = ca(ii).u(nanmask);
        nanmask = isnan(adcp.v);
        adcp.v(nanmask) = ca(ii).v(nanmask);        
    end  
end


clear ca


%% calculate shear in ADCP
% only use adcps... when current meters are used, erronous results are
% calculated for shear

[~, adcp.dudz] = gradient(adcp.u,dz);
[~, adcp.dvdz] = gradient(adcp.v,dz);
adcp.Ssq = adcp.dudz.^2 + adcp.dvdz.^2;
adcp.S = sqrt(adcp.Ssq);


%% load in additonal current meter datafiles if they exist

%%%%% parse in additional current meter files %%%%%

if ynadditionalfiles == 2 & exist('cc','var')
    
    % adcp_info_at_XX_XXX.m will have already been run above
%     run(['./adcp_info/adcp_info_at_' latstrshort '_' lonstr '.m']);
    
    % current meter data should already be concatenated into structure cc,
    % and it should be on the universal time grid. However, they may be on
    % different depth grids.
    
    % Parse new data into cur.u and cur.v in regions where cur.u and
    % cur.v are NaN    
    for ii = 1:length(cur.depth)
        clear indd
        indd = find(cc.depth == cur.depth(ii));
        if ~isempty(indd)
            clear nanmask
            nanmask = isnan(cur.u(ii,:));
            cur.u(ii,nanmask) = cc.u(indd,nanmask);
            clear nanmask
            nanmask = isnan(cur.v(ii,:));
            cur.v(ii,nanmask) = cc.v(indd,nanmask);
        end
    end
    
    % add cc.u and cc.v at depths that don't already exist in cur.depth
    for ii = 1:length(cc.depth)
        clear indd depthorig uorig vorig
        indd = find(cur.depth == cc.depth(ii));
        if ~isempty(indd)
            % this depth should have already been added to cur
        else
            depthorig = cur.depth;
            uorig = cur.u;
            vorig = cur.v;
            
            indd = find(cur.depth < cc.depth(ii),1,'last');
            if indd ~= length(cur.depth)
                cur.depth = [depthorig(1:indd); cc.depth(ii); ...
                    depthorig(indd+1:end)];
                cur.u = [cur.u(1:indd,:); cc.u(ii,:); cur.u(indd+1:end,:)];
                cur.v = [cur.v(1:indd,:); cc.v(ii,:); cur.v(indd+1:end,:)];
            else
                cur.depth = [depthorig(1:indd); cc.depth(ii)];
                cur.u = [cur.u(1:indd,:); cc.u(ii,:)];
                cur.v = [cur.v(1:indd,:); cc.v(ii,:)];
            end
        end
    end

    clear cc
end
    
%% remove extraneous depths in current meter files

DD = length(cur.depth);
indnext = 1;

if DD == 1 & cur.depth == 0
    % no good current meter data. Only dummy structure with NaNs remains
    clear cur
    
else
    % loop through and remove depths with all NaNs
    curgood.time = cur.time;
    for ii = 1:DD
        clear notnan
        notnan = ~isnan(cur.u(ii,:));
        if sum(notnan) == 0
            % this depth has non-NaN data: skip this depth
        else
            curgood.depth(indnext) = cur.depth(ii);
            curgood.u(indnext,:) = cur.u(ii,:);
            curgood.v(indnext,:) = cur.v(ii,:);
            indnext = indnext + 1;
        end 
    end

    clear cur
    cur = curgood;
    clear curgood
end

%% save the ADCP and current meter data in separate files
% this means there are now 3 velocity files being saved:
%   adcp: just ADCP data (saved here)
%   cur:  just current meter data (saved here)
%   vel:  adcp and current meter data combined (saved below)


%% calculate hourly, daily, and monthly averages

% hourly average (dt = 10 min, so put 6 points into each hour)
adcphr = bin_average(adcp,6);
adcphr.depth = adcp.depth;

% daily average (dt = 10 min, so put 6*24=144 points into each average)
adcpdy = bin_average(adcp,144);
adcpdy.depth = adcp.depth;

% monthly averages
adcpmo = monthly_average(adcpdy);
adcpmo.depth = adcp.depth;

% average current meter (if current meter data exists)
if exist('cur','var')
    curhr = bin_average(cur,6);
    curhr.depth = cur.depth;
    curdy = bin_average(cur,144);
    curdy.depth = cur.depth;
    curmo = monthly_average(curdy);
    curmo.depth = cur.depth;
end


% save

disp('saving...')

adcp.readme = strvcat(['ADCP velocity from the TAO mooring at ' latstrshort ',' lonstr],...
    ' ',...
    'ADCP is NOT augmented by current meters at shallow depths',...
    'If you would like the dataset where ADCP and current meters',...
    'have been combined, use vel.',...
    ' ',...
    'u is zonal velocity [m/s]',...
    'v is meridional velocity [m/s]',...
    'depth in meters',...
    'time is in matlab datenum format',...
    'Ssq is shear squared: (du/dz)^2 + (dv/dz)^2 [s^-2]',...
    'S is shear: sqrt(Ssq) [s^-1]',...
    '',...
    ['created by ' pwd '/load_and_save_tao_adcp_cur.m'],...
    ['by Sally Warner on ' datestr(now,'mmm dd, yyyy')]);

adcphr.readme = strvcat(adcp.readme,...
    ' ',...
    '10 min data has been binned into hourly averages');

adcpdy.readme = strvcat(adcp.readme,...
    ' ',...
    '10 min data has been binned into daily averages');

adcpmo.readme = strvcat(adcp.readme,...
    ' ',...
    'daily data has been averaged by month. A month nust have data from at',...
    'least 10 days to be averaged');



save([savedir 'adcp_' latstrshort '_' lonstr '_10min.mat'],'adcp','-v7.3')
save([savedir 'adcp_' latstrshort '_' lonstr '_hourly.mat'],'adcphr')
save([savedir 'adcp_' latstrshort '_' lonstr '_daily.mat'],'adcpdy')
save([savedir 'adcp_' latstrshort '_' lonstr '_monthly.mat'],'adcpmo')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('cur','var')

    cur.readme = strvcat(['Current meter velocity from the TAO mooring at ' latstrshort ',' lonstr],...
        ' ',...
        'Current meters are NOT combined with ADCP.',...
        'If you would like the dataset where ADCP and current meters',...
        'have been combined, use vel.',...
        ' ',...
        'u is zonal velocity [m/s]',...
        'v is meridional velocity [m/s]',...
        'depth in meters',...
        'time is in matlab datenum format',...
        'shear is not calculated with current meters due to sparceness in depth',...
        '',...
        ['created by ' pwd '/load_and_save_tao_adcp_cur.m'],...
        ['by Sally Warner on ' datestr(now,'mmm dd, yyyy')]);

    curhr.readme = strvcat(cur.readme,...
        ' ',...
        '10 min data has been binned into hourly averages');

    curdy.readme = strvcat(cur.readme,...
        ' ',...
        '10 min data has been binned into daily averages');

    curmo.readme = strvcat(cur.readme,...
        ' ',...
        'daily data has been averaged by month. A month nust have data from at',...
        'least 10 days to be averaged');

    save([savedir 'cur_' latstrshort '_' lonstr '_10min.mat'],'cur','-v7.3')
    save([savedir 'cur_' latstrshort '_' lonstr '_hourly.mat'],'curhr')
    save([savedir 'cur_' latstrshort '_' lonstr '_daily.mat'],'curdy')
    save([savedir 'cur_' latstrshort '_' lonstr '_monthly.mat'],'curmo')

end

%% special instructions for 140W
% 
% somehow there is bad data from February 29, 2016 through April 7, 2016 
% the dataset that interplates together the current meter and adcp.
% The current meter data itself doesn't look bad and neither do the winds,
% the badness somehow comes from the interpolating scheme. 
% I will just NaN out the current meter data here and hope to fix the
% problem by not doing any interpolation for this time period and just
% using the straight up adcp data.

if lon == 220
    indbadadcp = find(cur.time >= datenum(2016,02,28) & ...
        cur.time <= datenum(2016,04,08));
    cur.u(1:2,indbadadcp) = NaN; % index 1-8 corresponds to 10m and 25m current meters
    cur.v(1:2,indbadadcp) = NaN;

    clear indbadadcp
end


%% follow Pham et al 2017 Appendix B to calculate velocities to surface
% Bill came up with a method to use the ADCP and current meter data to
% estimate the currents all the way to the surface

% H = mixed layer depth (ie. depth where temp is 0.04°C less than surface temp)
% uw = wind velocity
% cw = constant = 0.02
% uz17 = shear at 17.5 m depth (between current meters at 10 and 25 m)
% uza1 = shear between upper two most adcp bins
% uza2 = shear between second and third adcp bins
% zb = depth of upmost ADCP bin

if exist('cur','var')

    disp('extrapolating velocity to surface using Bill''s method...')

    % load in wind data
    load([savedir 'wind_' latstrshort '_' lonstr '_10min.mat'])
    uw = wind.u;
    vw = wind.v;
    cw = 0.02;

    % load in mixed layer depth data
    load([savedir 'mld_' latstrshort '_' lonstr '_10min.mat'])
    H = mld.mld_forshear;

    % calculate surface shear based on wind and MLD
    uz0 = cw * uw ./ H;
    vz0 = cw * vw ./ H;
    z0 = 0*wind.time;

    %%%%%%%% find shear at top three points where we know shear
    % shear from current meters at 10 and 25 m depths
    ind10 = find(cur.depth == 10);
    ind25 = find(cur.depth == 25);

    if ~isempty(ind10) & ~isempty(ind25)
        z1 = 0.5*(cur.depth(ind10)+cur.depth(ind25)) + 0*cur.time;
        uz1 = (cur.u(ind25,:) - cur.u(ind10,:)) / (cur.depth(ind25)-cur.depth(ind10));
        vz1 = (cur.v(ind25,:) - cur.v(ind10,:)) / (cur.depth(ind25)-cur.depth(ind10));

        % shear from first and second adcp bins
        inda1 = NaN*adcp.time;
        z2 = NaN*adcp.time;
        z3 = z2;
        uz2 = z2;
        uz3 = z2;
        vz2 = z2;
        vz3 = z2;
        for ii = 1:length(adcp.time)

            if ~isnan(cur.u(ind25,ii)) & ~isnan(cur.u(ind10,ii))

                % find depth up uppermost adcp bin with data
                dummy = find(~isnan(adcp.u(:,ii)),1,'first');
                if ~isempty(dummy)
                    inda1(ii) = dummy;

                    % depths of shear
                    z2(ii) = 0.5*(adcp.depth(inda1(ii))+adcp.depth(inda1(ii)+1));
                    z3(ii) = 0.5*(adcp.depth(inda1(ii)+1)+adcp.depth(inda1(ii)+2));

                    % shear
                    uz2(ii) = (adcp.u(inda1(ii)+1,ii) - adcp.u(inda1(ii),ii)) /...
                        (adcp.depth(inda1(ii)+1) - adcp.depth(inda1(ii)));
                    vz2(ii) = (adcp.v(inda1(ii)+1,ii) - adcp.v(inda1(ii),ii)) /...
                        (adcp.depth(inda1(ii)+1) - adcp.depth(inda1(ii)));
                    uz3(ii) = (adcp.u(inda1(ii)+2,ii) - adcp.u(inda1(ii)+1,ii)) /...
                        (adcp.depth(inda1(ii)+2) - adcp.depth(inda1(ii)+1));
                    vz3(ii) = (adcp.v(inda1(ii)+2,ii) - adcp.v(inda1(ii)+1,ii)) /...
                        (adcp.depth(inda1(ii)+2) - adcp.depth(inda1(ii)+1));
                else
                    z2(ii) = 35;
                    z3(ii) = 40;
                    uz2(ii) = NaN;
                    vz2(ii) = NaN;
                    uz3(ii) = NaN;
                    vz3(ii) = NaN;       
                end
            end   
        end


        % calculate constants at each time step
        % calculate shear at all depths from 0 to adcp.depth(inda1(ii))
        % calculate velocity at all depths from 0 to adcp.depth(inda1(ii))
        uz = NaN*adcp.u;
        u = NaN*adcp.u;
        vz = NaN*adcp.u;
        v = NaN*adcp.u;
        zz = adcp.depth;
        % create new structure for the data that will be filled to the surface
        vel = adcp;
        % go through step by step
        for ii = 1:length(adcp.time)

            if mod(ii,10000) == 0
                disp([num2str(ii) '/' num2str(length(adcp.time))])
            end

            if ~isnan(inda1(ii))
                clear abc abcinv coefsu bu cu du zb coefsv bv cv dv surfind
                abc = [z1(ii) z1(ii)^2 z1(ii)^3; z2(ii) z2(ii)^2 z2(ii)^3; z3(ii) z3(ii)^2 z3(ii)^3];
                abcinv = abc^-1;

                zb = adcp.depth(inda1(ii));

                coefsu = abcinv*[uz1(ii)-uz0(ii); uz2(ii)-uz0(ii); uz3(ii)-uz0(ii)];
                bu = coefsu(1);
                cu = coefsu(2);
                du = coefsu(3);    
                uz(:,ii) = uz0(ii) + bu*zz + cu*zz.^2 + du*zz.^3;
                u(:,ii) = adcp.u(inda1(ii),ii) + uz0(ii)*(zz - zb) + bu/2*(zz.^2 - zb^2) + ...
                    cu/3*(zz.^3 - zb^3) + du/4*(zz.^4 - zb^4);

                coefsv = abcinv*[vz1(ii)-vz0(ii); vz2(ii)-vz0(ii); vz3(ii)-vz0(ii)];
                bv = coefsv(1);
                cv = coefsv(2);
                dv = coefsv(3);    
                vz(:,ii) = vz0(ii) + bv*zz + cv*zz.^2 + dv*zz.^3;
                v(:,ii) = adcp.v(inda1(ii),ii) + vz0(ii)*(zz - zb) + bv/2*(zz.^2 - zb^2) + ...
                    cv/3*(zz.^3 - zb^3) + dv/4*(zz.^4 - zb^4);

                % put data into filled structure
                surfind = 1:inda1(ii);
                vel.u(surfind,ii) = u(surfind,ii);
                vel.v(surfind,ii) = v(surfind,ii);
                vel.dudz(surfind,ii) = uz(surfind,ii);
                vel.dvdz(surfind,ii) = vz(surfind,ii); 

            else
                % leave uz, u, vz, and v as NaNs when inda1(ii) is NaN
            end

        end

        % calculate Ssq and S
        vel.Ssq = vel.dudz.^2 + vel.dvdz.^2;
        vel.S = sqrt(vel.Ssq);

    else

        % case with no current meters at 25m and 10m at any time

        % If there are current meters at ADCP depths, parse in that data in
        % regions where ADCP is NaN
        for ii = 1:length(cur.depth)
            clear indd
            indd = find(adcp.depth == cur.depth(ii));
            if ~isempty(indd)
                clear indnan
                indnan = isnan(adcp.u(indd,:));
                adcp.u(indd,indnan) = cur.u(ii,indnan);
                adcp.v(indd,indnan) = cur.v(ii,indnan);
            end
        end
        % don't calculate shear with current meters in this case

        % save as vel (the structure that has both ADCP and current meter data)
        vel = adcp;
    end

else
    vel = adcp;
end



%% calculate hourly, daily, and monthly averages

disp('binning velocity into hourly, daily, and monthly averages...')

% hourly average (dt = 10 min, so put 6 points into each hour)
velhr = bin_average(vel,6);
velhr.depth = vel.depth;

% daily average (dt = 10 min, so put 6*24=144 points into each average)
veldy = bin_average(vel,144);
veldy.depth = vel.depth;

% monthly averages
velmo = monthly_average(veldy);
velmo.depth = vel.depth;



%% save

disp('saving...')

vel.readme = strvcat(['Velocity from the TAO mooring at ' latstrshort ',' lonstr],...
    ' ',...
    'ADCP is augmented by current meters (if available) at shallow depths',...
    '(>40m) using the method described in Appendix B of Pham, Smyth,',...
    'Sarkar, and Moum (JPO, 2017). Wind stress at the surface,',...
    'current meters at 10m and 25m, and the ADCP are used to find',...
    'shear and velocity to the surface. In cases where no current',...
    'meters were available, only the ADCP data is used.',...
    ' ',...
    'u is zonal velocity [m/s]',...
    'v is meridional velocity [m/s]',...
    'depth in meters',...
    'time is in matlab datenum format',...
    'Ssq is shear squared: (du/dz)^2 + (dv/dz)^2 [s^-2]',...
    'S is shear: sqrt(Ssq) [s^-1]',...
    '',...
    ['created by ' pwd '/load_and_save_tao_adcp_cur.m'],...
    ['by Sally Warner on ' datestr(now,'mmm dd, yyyy')]);

velhr.readme = strvcat(vel.readme,...
    ' ',...
    '10 min data has been binned into hourly averages');

veldy.readme = strvcat(vel.readme,...
    ' ',...
    '10 min data has been binned into daily averages');

velmo.readme = strvcat(vel.readme,...
    ' ',...
    'daily data has been averaged by month. A month nust have data from at',...
    'least 10 days to be averaged');

save([savedir 'vel_' latstrshort '_' lonstr '_10min.mat'],'vel','-v7.3')
save([savedir 'vel_' latstrshort '_' lonstr '_hourly.mat'],'velhr')
save([savedir 'vel_' latstrshort '_' lonstr '_daily.mat'],'veldy')
save([savedir 'vel_' latstrshort '_' lonstr '_monthly.mat'],'velmo')





return

%% plot

dd = dir([savedir 'vel*']);

if ~isempty(dd)
    
for ii = 1:length(dd)
    if strcmp(dd(ii).name(end-8:end),'daily.mat')
        disp(dd(ii).name)
        load([savedir dd(ii).name])
    end
end

vars = {'u';'v';'dudz';'dvdz';'Ssq';'S'};


figure(489)
clf
set(gcf,'position',[376          57        1157        1048])

for ii = 1:length(vars)
    
    ax(ii) = subplot(length(vars),1,ii);
    
    pcolor(veldy.time,veldy.depth,veldy.(vars{ii}))
    shading flat
    axis ij
    colormap(redblue3)
    colorbar
    if ii == 1 | ii == 2
        caxis([-1 1]*3)
    elseif ii == 3 | ii == 4
        caxis([-1 1]*0.04)
    elseif ii == 5
        caxis([0 2]*10^-3)
    else
        caxis([0 0.04])
    end
        
    ylabel(vars{ii})
    datetick
    set(gca,'tickdir','out')
    
    if ii == 1
        title(['Velocity at ' latstrshort ', ' lonstr])
    end
end

linkaxes

export_fig([figdir 'vel.png'],'-r200')

end


%%%%%%%%%%%%%%%%%%%%

dd = dir([savedir 'adcp*']);

if ~isempty(dd)

for ii = 1:length(dd)
    if strcmp(dd(ii).name(end-8:end),'daily.mat')
        disp(dd(ii).name)
        load([savedir dd(ii).name])
    end
end

vars = {'u';'v';'dudz';'dvdz';'Ssq';'S'};


figure(589)
clf
set(gcf,'position',[376          57        1157        1048])

for ii = 1:length(vars)
    
    ax(ii) = subplot(length(vars),1,ii);
    
    pcolor(adcpdy.time,adcpdy.depth,adcpdy.(vars{ii}))
    shading flat
    axis ij
    colormap(redblue3)
    colorbar
    if ii == 1 | ii == 2
        caxis([-1 1]*3)
    elseif ii == 3 | ii == 4
        caxis([-1 1]*0.04)
    elseif ii == 5
        caxis([0 2]*10^-3)
    else
        caxis([0 0.04])
    end
        
    ylabel(vars{ii})
    datetick
    set(gca,'tickdir','out')
    
    if ii == 1
        title(['ADCP Velocity at ' latstrshort ', ' lonstr])
    end
end

linkaxes

export_fig([figdir 'adcp.png'],'-r200')

end



%%%%%%%%%%%%%%%%%%%%


dd = dir([savedir 'cur*']);

if ~isempty(dd)

for ii = 1:length(dd)
    if strcmp(dd(ii).name(end-8:end),'daily.mat')
        disp(dd(ii).name)
        load([savedir dd(ii).name])
    end
end

vars = {'u';'v'};


figure(4689)
clf
set(gcf,'position',[376          57        1157        1048])

for ii = 1:length(curdy.depth)
    
    ax(ii) = subplot(length(curdy.depth),1,ii);
    
    plot(curdy.time,curdy.u(ii,:))
    hold on
    plot(curdy.time,curdy.v(ii,:))
    ylim([-2 2])
        
    ylabel({[num2str(curdy.depth(ii)) ' m'];'u,v [m s^{-1}]'})
    datetick
    set(gca,'tickdir','out')
    
    if ii == 1
        title(['Current Meter Velocity at ' latstrshort ', ' lonstr])
    end
end

linkaxes

export_fig([figdir 'cur.png'],'-r200')


end