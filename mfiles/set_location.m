% set_location 
%
% Use when running the TaoTritonPriataRama gridding suite of codes.
% Set the location where you want to create gridded data.
%
% Format should be a string that includes 'W' or 'E' (eg. '10W' or '140W')
% for longitude, and includes 'N' or 'S' for latitude (e.g. '0N' or '5S').
%
% written by Sally Warner

latstr = '0N';
% latstr = '1.5S';

% lonstr = '10W';
% lonstr = '23W';
% lonstr = '95W';
% lonstr = '110W';
% lonstr = '125W';
% lonstr = '140W';
lonstr = '80.5E';
% lonstr = '67E';


disp(['Processing TaoTritonPirataRama data at ' latstr ', ' lonstr])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% do not change code below %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create necessary strings and numbers for lat and lon

lon = str2num(lonstr(1:end-1));
if strcmp(lonstr(end),'W')
    lon = 360-lon;
end
lonstrnop = lonstr;
for ii = 1:length(lonstr)
    if strcmp(lonstr(ii),'.')
        lonstrnop(ii) = '_';
    end
end


lat = str2num(latstr(1:end-1));
if strcmp(latstr(end),'S')
    lat = -lat;
end
if strcmp(latstr,'0N')
    latstrshort = '0';
else
    latstrshort = latstr;
end
latstrnop = latstr;
for ii = 1:length(latstr)
    if strcmp(latstr(ii),'.')
        latstrnop(ii) = '_';
    end
end
latstrshortnop = latstrshort;
for ii = 1:length(latstrshort)
    if strcmp(latstrshort(ii),'.')
        latstrshortnop(ii) = '_';
    end
end

%% define save and figure directories

time = taotime;

if time(1) < datenum(2005,1,1)
    maindir = ['../../TaoTritonPirataRama_processed_at_' latstrshort '_' lonstr '_from1990'];
else
    maindir = ['../../TaoTritonPirataRama_processed_at_' latstrshort '_' lonstr];
end

if exist(maindir,'dir') == 0
    mkdir(maindir)
end

savedir = [maindir '/processed/'];
if exist(savedir,'dir') == 0
    mkdir(savedir)
end

figdir = [maindir '/figures/'];
if exist(figdir,'dir') == 0
    mkdir(figdir)
end

%% determine the name of the mooring array (for loading daily data)

if lon > 30 & lon <= 120 % Indian
    mooringarray = 'RAMA';
elseif lon > 120 & lon <= 285 % Pacific
    mooringarray = 'TAO_TRITON';
elseif (lon > 285 & lon <= 360) | (lon >= 0 & lon <= 30) % Atlantic
    mooringarray = 'PIRATA';
end



