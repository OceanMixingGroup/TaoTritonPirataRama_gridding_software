% chipod_info_at_0_110W
%
% In this module, specify necessary information about chipods at each
% mooring location.
%
% Sally Warner

%% define depths

% 'depths' defines the depths which you would like to appear in the final
% gridded output. Essentially, this vector should have all the depths at
% which chipods have ever been deployed at this location.
chipodinfo.depths = [9 19 29 40 49];



%% filenames of saved data OLD PROCESSING

% Define parameters for old chipods that were processed with the old 
% processing code (i.e. in the days when Sasha was processing chipods).

% If there are no files from old processing (likely the case unless you're
% processing data from 0, 140W or 0, 110W), set chipodinfo.ynoldchipods = 0;

% set ynoldchipods = 1 if old chipods should be included
% set ynoldchipods = 0 if old chipods should NOT be included
chipodinfo.ynoldchipods = 1;

% paths to old chipod data summary files. If they do not exist (or one
% exists and not the other), set filename for nonexistant files to 'none'.
if chipodinfo.ynoldchipods

    chipodinfo.filenames.oldproc_i = ['../../chipod/tao0N110W_newest_summaries/'...
        'tao0N110W_chipod_summary_10m_12-Mar-2018.mat'];

    chipodinfo.filenames.oldproc_m = 'none';
end


%% depth matrix for OLD PROCESSING
% ignore this part if you don't have old chipods

if chipodinfo.ynoldchipods
% define the index of each depth in cpd structure for each deployment:
%              dpl #: 1
chipodinfo.inddeps = [0; ...    % 9m
                      1; ...    % 19m
                      2; ...    % 29m
                      3; ...    % 40m
                      4]; ...   % 49m

% For this example, during the first deployment (left-most column), chipods
% were deployed at 29m, 49m, and 84m, and they have indices of 1, 2, 3,
% respectively, in the data file. 

end

%% filenames for saved data NEW PROCESSING

% You will almost certainly have new chipods you would like to add to the
% gridded data.

% TAO10_110: lost

% TAO11_110: lost

% TAO13_110: 9m chipod = unit 512

% TAO14_110: lost

% TAO15_110: 9m chipod = unit 734
% TAO15_110: 29m chipod = unit 735
% TAO15_110: 49m chipod = unit 736

% TAO2016: 9m chipod = unit 523
% TAO2016: 29m chipod = unit 525 (batteries exploded)
% TAO2016: 49m chipod = unit 707 (lost to vandalism)

% tao17_110: 9m chipod = unit 734 (bad SD card, no raw data)
% tao17_110: 29m chipod = unit 735 (bad SD card, 1 day raw data)
% tao17_110: 49m chipod = unit 736 (bad SD card, no raw data)

% tao18_110: 9m chipod = unit 811 
% tao18_110: 29m chipod = unit 812
% tao18_110: 49m chipod = unit 813 


% dpl = {'TAO13_110';'TAO15_110';'TAO2016'};
% newdepths = [9    29   49];
% cpds =      [512  0    0 ;        % TAO13_110
%              734  735  736;       % TAO15_110
%              523  0    0;         % TAO16_110
%              0    735  0;         % tao17_110 (no good data)
%              811  812  813];      % tao18_110

chipodinfo.dpl = {'tao13_110';'tao15_110';'tao16_110';'tao18_110'};
chipodinfo.newdepths = [9    29   49];
chipodinfo.cpds =      [512  0    0 ;        % tao13_110
                        734  735  736;       % tao15_110
                        523  0    0;         % tao16_110
                        811  812  813];      % tao18_110


%%%%%%%%%%% directories for each deployment %%%%%%%%%%%
% (This directory should go from the gridded data mfiles folder to the
% folder where all the unit names are saved for that deployment.
for ii = 1:length(chipodinfo.dpl)
    chipodinfo.basedir.(chipodinfo.dpl{ii}) = ['../../chipod/' ...
        chipodinfo.dpl{ii} '/data/'];
end

% if any of the directories are inconsistent with the following structure:
%       ../../chipod/DEPLOYMENT_NAME/data/
% where DEPLOYMENT_NAME = chipodinfo.dpl{ii}, they can be set manually as:
%
% chipodinfo.basedir.tao13_110 = '../../chipod/tao13_110/data/';
% chipodinfo.basedir.tao15_110 = '../../chipod/tao15_110/data/';
% chipodinfo.basedir.tao16_110 = '../../chipod/tao16_110/data/';

