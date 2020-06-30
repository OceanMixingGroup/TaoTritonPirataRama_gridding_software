% chipod_info_at_0_140W
%
% In this module, specify necessary information about chipods at each
% mooring location.
%
% Sally Warner

%% define depths

% 'depths' defines the depths which you would like to appear in the final
% gridded output. Essentially, this vector should have all the depths at
% which chipods have ever been deployed at this location.
chipodinfo.depths = [29 39 49 59 69 89 119];



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

    chipodinfo.filenames.oldproc_i = ['../../chipod/tao0N140W_newest_summaries/'...
        'tao0N140W_chipod_summary_10m_28-Mar-2018.mat'];

    chipodinfo.filenames.oldproc_m = ['../../chipod/tao0N140W_newest_summaries/'...
        'tao0N140W_chipod_summary_10m_N2_28-Mar-2018.mat'];
end


%% depth matrix for OLD PROCESSING
% ignore this part if you don't have old chipods

if chipodinfo.ynoldchipods
% define the index of each depth in cpd structure for each deployment:
% dpl #:   1 2 3 4 5 6 7 8 9
chipodinfo.inddeps = [1 1 1 0 2 1 1 0 1; ...    % 29m
                      0 2 2 2 0 2 2 0 0; ...    % 39m
                      2 3 3 0 3 3 0 1 2; ...    % 49m
                      0 4 4 3 6 4 3 0 3; ...    % 59m
                      0 0 0 4 0 5 0 0 0; ...    % 69m
                      3 0 0 5 7 6 4 0 0; ...    % 84 & 89m
                      0 0 0 6 0 0 0 2 0];       % 124 & 129

% For this example, during the first deployment (left-most column), chipods
% were deployed at 29m, 49m, and 84m, and they have indices of 1, 2, 3,
% respectively, in the data file. 

end

%% filenames for saved data NEW PROCESSING

% You will almost certainly have new chipods you would like to add to the
% gridded data.

% tao14_140: 29m chipod = unit 707
% tao14_140: 49m chipod = unit 713
% tao14_140: 69m chipod = unit 710
% tao14_140: 89m chipod = unit 711
% tao14_140: 119m chipod = unit 712
% (no IC estimates included with tao14_140)

% tao15_140: 29m chipod = unit 726
% tao15_140: 49m chipod = unit 727
% tao15_140: 69m chipod = unit 728
% tao15_140: 89m chipod = unit 729
% tao15_140: 119m chipod = unit 730
% (no IC estimates included with tao15_140)

% tao16_140: 29m chipod = unit 500 - 7 days
% tao16_140: 49m chipod = unit 501 - 3.5 days
% tao16_140: 69m chipod = unit 503 - 1.5 days
% tao16_140: 89m chipod = unit 504 - 1.5 days
% tao16_140: 119m chipod = unit 505 - All data bad
% (IC estimates are included with tao16_140)

% tao17_140: 29m chipod = unit 715 - All data bad
% tao17_140: 49m chipod = unit 716 - All data bad
% tao17_140: 69m chipod = unit 717 - 45 days (with many gaps)
% tao17_140: 89m chipod = unit 718 - No Raw data
% tao17_140: 119m chipod = unit 719 - 4 days of data in a 12-day period, all data bad
% (IC estimates are included with tao17_140 when good)

% tao18_140: 29m chipod = unit 710
% tao18_140: 49m chipod = unit 713
% tao18_140: 69m chipod = unit 714
% tao18_140: 89m chipod = unit 711
% tao18_140: 119m chipod = unit 722
% (IC estimates are included with tao18_140 when good)

% tao19_140: 29m chipod = unit 1111
% tao19_140: 49m chipod = unit 1136
% tao19_140: 69m chipod = unit 1137
% tao19_140: 89m chipod = unit 1138
% tao19_140: 119m chipod = unit 1139 - data recorded but both T1P and T2P are bad, so only have IC estimates
% (IC estimates are included with tao19_140 when good)

chipodinfo.dpl = {'tao14_140';'tao15_140';'tao16_140';...
                  'tao17_140';'tao18_140';'tao19_140'};
chipodinfo.newdepths = [29   49   69   89   119];
chipodinfo.cpds =      [707  713  710  711  712; ...    % tao14_140
                        726  727  728  729  730;        % tao15_140
                        500  501  503  0    0  ;        % tao16_140
                        0    0    717  0    0  ;        % tao17_140
                        710  713  714  711  722;        % tao18_140
                        1111 1136 1137 1138 1139];      % tao19_140


%%%%%%%%%%% directories for each deployment %%%%%%%%%%%
% (This directory should go from the gridded data mfiles folder to the
% folder where all the unit names are saved for that deployment.
% for ii = 1:length(chipodinfo.dpl)
%     chipodinfo.basedir.(chipodinfo.dpl{ii}) = ['../../chipod/' ...
%         chipodinfo.dpl{ii} '/data/'];
% end

for ii = 1:length(chipodinfo.dpl)
    chipodinfo.basedir.(chipodinfo.dpl{ii}) = ['~/ganges/data/chipod/' ...
        chipodinfo.dpl{ii} '/data/'];
end

% if any of the directories are inconsistent with the following structure:
%       ../../chipod/DEPLOYMENT_NAME/data/
% where DEPLOYMENT_NAME = chipodinfo.dpl{ii}, they can be set manually as:
%
% chipodinfo.basedir.tao14_140 = '../../chipod/tao14_140/data/';
% chipodinfo.basedir.tao15_140 = '../../chipod/tao15_140/data/';

