% adcp_info_at_0_140W
%
% In this module, specify necessary information about adcps at each
% mooring location.
%
% If all of the ADCP data is already available in:
% ~/ganges/data/TaoTritonPirataRama/high_resolution/hr/adcpXnXXXw_hr.cdf,
% then there is no need to create this m-file for your mooring.
%
% Sally Warner

%% set whether or not you have additional adcp and current meter files.

% set adcpinfo.ynadcp = 1 if you have adcp files that are not in the
% TaoTritonPirataRama database (i.e. high_resolution/hr/adcp_0n140w_hr.cdf
% set adcpinfo.ynadcp = 0 if you do not have any additional ADCP files
adcpinfo.ynadcp = 1;

% set adcpinfo.yncurrent = 1 if you have current meter files that are not in the
% TaoTritonPirataRama database (i.e. high_resolution/hr/cur_0n140w_hr.cdf
% set adcpinfo.yncurrent = 0 if you do not have any additional current 
% meter files
adcpinfo.yncurrent = 0;


%% locations of additional ADCP files

if adcpinfo.ynadcp

% List all of the additional ADCP files that you would like to load. These
% paths should go from the TaoTritonPirataRama_mastercode/mfiles/adcp_info
% folder

% for some reason, when I run this code independently it wants a different
% basedir than if it is called by load_and_save_tao_adcp
% basedir1 = '\\ganges\data\chipod\RAMA17_67E\mooring_data\processed\';
basedir1 = '~/ganges/data/chipod/RAMA17_67E/mooring_data/processed/';

filenames = {'adcp_RAMA_0N67E_Aug2017_to_Sep2019_use_only_for_chipod_processing.mat'};

if exist([basedir1 filenames{1}],'file') == 2
    disp('using basedir1')
    adcpinfo.filenames = [basedir1 filenames{1}];
    adcpinfo.numfiles = 1;
else
    disp('Cannot find auxiliary ADCP file!')
    disp(['looking for: ' basedir1 filenames{1}])
end

%% load each additional adcp file

% load in each file in adcpinfo.filenames and change the variable names so
% that a structure call "ca" can be passed back to the main part of the
% code. ca should contain ALL loaded files and have the following variables:
%     depth (equivalent to taodepth)
%     time (equivalent to taotime)
%     u (in m/s)
%     v (in m/s)
for dd = 1:adcpinfo.numfiles
    
    % load in both .mat and/or netcdf files
    disp(['loading additional ADCP file ' num2str(dd) ' of ' ...
        num2str(adcpinfo.numfiles)])
    if adcpinfo.numfiles > 1
        matorcdf = adcpinfo.filenames{dd}(end-3:end);
        filename = adcpinfo.filenames{dd};
    else
        matorcdf = adcpinfo.filenames(end-3:end);
        filename = adcpinfo.filenames;
    end
        
    if strcmp(matorcdf,'.mat')
        dummy = load(filename);
        raw = dummy.adcp;
    else
        raw = loadnc(filename);
    end
       
    
    %%%%% interpolate to common depth and time grid %%%%%

    % predefine ca structure
    ca(dd).time = time;
    ca(dd).depth = depth;
    ca(dd).u = NaN*ones(length(ca(dd).depth),length(ca(dd).time));
    ca(dd).v = ca(dd).u;

    % first put the raw data on an even depth grid
    disp('interpolating to correct depth grid...')

    % first try a fast way
    [i,j] = size(raw.depth);
    if (i==1 | j==1) & nanmean(diff(raw.depth)) == dz
        cc = 1;
        for ii = 1:length(raw.depth)
            clear inda
            inda = find(ca(dd).depth == raw.depth(ii));
            if ~isempty(inda)
                inddadcp(cc) = inda;
                cc = cc + 1;
            end
        end
        cc = 1;
        for ii = 1:length(ca(dd).depth)
            clear inda
            inda = find(raw.depth == ca(dd).depth(ii));
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
            if size(raw.depth) == size(raw.u)
                raw.uinterpz(:,ii) = interp1(raw.depth(:,ii),raw.u(:,ii),ca(dd).depth);
                raw.vinterpz(:,ii) = interp1(raw.depth(:,ii),raw.v(:,ii),ca(dd).depth);
            else
                raw.uinterpz(:,ii) = interp1(raw.depth,raw.u(:,ii),ca(dd).depth);
                raw.vinterpz(:,ii) = interp1(raw.depth,raw.v(:,ii),ca(dd).depth);
            end
        end
    end

    % loop through all depths to interpolate to common time grid
    disp('interpolating to correct time grid...')
    for ii = 1:length(ca(dd).depth)
        ca(dd).u(ii,:) = interp1(raw.time,raw.uinterpz(ii,:),ca(dd).time);
        ca(dd).v(ii,:) = interp1(raw.time,raw.vinterpz(ii,:),ca(dd).time);
    end

end

clear raw

end



%% current meter info
% load current meter data 

if adcpinfo.yncurrent
    
disp('loading additional current meter files...')

% create structure for current meter data
cc.time = time;
cc.depth = [10 25 45 80 120];
cc.u = NaN*ones(length(cc.depth),length(cc.time));
cc.v = cc.u;

% directory for current meter data
curdir = '../../../tao_array/toApril2017/CURRENTS/';
dd = dir([curdir 'TAO_T0N140W*']);
if isempty(dd)
    curdir = '../../tao_array/toApril2017/CURRENTS/';
    dd = dir([curdir 'TAO_T0N140W*']);
end

% load all files and interpolate to universal time grid (you do not need to
% interpolate to the common depth grid for current meters)
for ii = 1:length(dd)
    disp(dd(ii).name)
    clear raw
    raw = loadnc([curdir dd(ii).name]);
    for jj = 1:length(raw.DEPTH)
        clear indcurd
        indcurd = find(cc.depth == raw.DEPTH(jj));
        if ~isempty(indcurd)
            clear dummyu dummyv
            dummyu = interp1(raw.time,raw.UCUR(jj,:)/100,cc.time);
            dummyv = interp1(raw.time,raw.VCUR(jj,:)/100,cc.time);
            cc.u(indcurd,~isnan(dummyu)) = dummyu(~isnan(dummyu));
            cc.v(indcurd,~isnan(dummyv)) = dummyv(~isnan(dummyv));     
        end
    end
end

clear dd curdir ii raw dummyu dummyv

end