% adcp_info_at_0_100W
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

% List all of the additional ADCP files that you would like to load:

adcpinfo.filenames = {'~/ganges/data/chipod/Pirata17/IRD_data/ADCP-moorings/10W-0N/2012-2013/10W0N_15258_instr_01.mat';...
                      '~/ganges/data/chipod/Pirata17/IRD_data/ADCP-moorings/10W-0N/2014-2015/FR25-10W0N_15258_instr_01.mat';...
                      '~/ganges/data/chipod/Pirata17/IRD_data/ADCP-moorings/10W-0N/2015-2017/10W0N_15258_instr_01.mat';...
                      '~/ganges/data/chipod/Pirata17/IRD_data/ADCP-moorings/10W-0N/2017-2019/0N10W_15258_instr_01.mat'};

adcpinfo.numfiles = length(adcpinfo.filenames);

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
    disp(adcpinfo.filenames{dd})
    
    if adcpinfo.numfiles > 1
        matorcdf = adcpinfo.filenames{dd}(end-3:end);
        filename = adcpinfo.filenames{dd};
    else
        matorcdf = adcpinfo.filenames{1}(end-3:end);
        filename = adcpinfo.filenames{1};
    end
    
%     if strcmp(matorcdf,'.mat')
%         dummy = load(filename);
%         raw = dummy.adcp;
%     else
%         raw = loadnc(filename);
%     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % special instructions for ADCP data from 10W

    % much faster to load in just data structure rather than whole .mat
    % file
    load(filename,'data');
    raw = data;
    raw.timejulian = raw.time;
    raw.time = datenum(datetime(raw.timejulian,'convertfrom','juliandate'));
    raw.instrumentdepth = raw.depth;
    raw.depth = raw.z_bins;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        

       
    
    
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


%%
clear raw

end
end



%% current meter info
% load current meter data (if available)

% if adcpinfo.yncurrent
%     
%     disp('loading additional current meter files...')
% 
%     % create structure for current meter data
%     cc.time = time;
%     cc.depth = [10 25 45 80 120];
%     cc.u = NaN*ones(length(cc.depth),length(cc.time));
%     cc.v = cc.u;
% 
%     % directory for current meter data
%     curdir = '~/ganges/data/tao_array/toApril2017/CURRENTS/';
%     dd = dir([curdir 'TAO_T0N140W*']);
% 
% 
%     % load all files and interpolate to universal time grid (you do not need to
%     % interpolate to the common depth grid for current meters)
%     for ii = 1:length(dd)
%         disp(dd(ii).name)
%         clear raw
%         raw = loadnc([curdir dd(ii).name]);
%         for jj = 1:length(raw.DEPTH)
%             clear indcurd
%             indcurd = find(cc.depth == raw.DEPTH(jj));
%             if ~isempty(indcurd)
%                 clear dummyu dummyv
%                 dummyu = interp1(raw.time,raw.UCUR(jj,:)/100,cc.time);
%                 dummyv = interp1(raw.time,raw.VCUR(jj,:)/100,cc.time);
%                 cc.u(indcurd,~isnan(dummyu)) = dummyu(~isnan(dummyu));
%                 cc.v(indcurd,~isnan(dummyv)) = dummyv(~isnan(dummyv));     
%             end
%         end
%     end
% 
%     clear dd curdir ii raw dummyu dummyv
% 
% end




