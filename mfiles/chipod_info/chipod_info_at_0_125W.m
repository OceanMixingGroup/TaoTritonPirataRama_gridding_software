% chipod_info_at_0_125W
%
% In this module, specify necessary information about chipods at each
% mooring location.
%
% Sally Warner

%% define depths

% 'depths' defines the depths which you would like to appear in the final
% gridded output. Essentially, this vector should have all the depths at
% which chipods have ever been deployed at this location.
chipodinfo.depths = [30 50 70];



%% filenames of saved data OLD PROCESSING

% Define parameters for old chipods that were processed with the old 
% processing code (i.e. in the days when Sasha was processing chipods).

% If there are no files from old processing (likely the case unless you're
% processing data from 0, 140W or 0, 110W), set chipodinfo.ynoldchipods = 0;

% set ynoldchipods = 1 if old chipods should be included
% set ynoldchipods = 0 if old chipods should NOT be included
chipodinfo.ynoldchipods = 0;

% paths to old chipod data summary files. If they do not exist (or one
% exists and not the other), set filename for nonexistant files to 'none'.
if chipodinfo.ynoldchipods

    chipodinfo.filenames.oldproc_i = 'none';
    chipodinfo.filenames.oldproc_m = 'none';
end


%% depth matrix for OLD PROCESSING
% ignore this part if you don't have old chipods

if chipodinfo.ynoldchipods
% define the index of each depth in cpd structure for each deployment:
%
% see chipod_info_at_0_140W.m for a good example

end

%% filenames for saved data NEW PROCESSING

% You will almost certainly have new chipods you would like to add to the
% gridded data.

% tao14_125: 30m chipod = unit 520
% tao14_125: 50m chipod = unit 521
% tao14_125: 70m chipod = unit 523

% tao15_125: 30m chipod = unit 723 (don't include IC in averages)
% tao15_125: 50m chipod = unit 724 (don't include IC in averages)
% tao15_125: 70m chipod = unit 725 (don't include IC in averages)

% tao16_125: 30m chipod = unit 511 (bad SD cards, only 12 days of data)
% tao16_125: 50m chipod = unit 1112 (bad SD cards, all recorded data was junk)
% tao16_125: 70m chipod = unit 1113 (only IC results are believable)

% tao17_125: 30m chipod = unit 709 (bad SD card, no raw data)
% tao17_125: 50m chipod = unit 720 (chipod lost)
% tao17_125: 70m chipod = unit 726 (chipod lost)

% tao18_125: 30m chipod = unit 723 (3 months good data; IC & VC both good; T1 & T1P good, but T2 & T2P bad)
% tao18_125: 50m chipod = unit 724 (3 months good data; IC good for both T1 & T2; VC for T2P is good, VC for T1P is bad)
% tao18_125: 70m chipod = unit 712 (7 months good data; IC & VC both good; T1 & T1P good, but T2 & T2P bad)

chipodinfo.dpl = {'tao14_125';'tao15_125';'tao16_125';'tao18_125'};
chipodinfo.newdepths = [30    50    70];
chipodinfo.cpds =      [520   521   523 ;      % tao14_125
                        723   724   725;       % tao15_125
                        511   0     1113;      % tao16_125
                        723   724   712];      % tao18_125


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
% chipodinfo.basedir.tao14_140 = '../../chipod/tao14_140/data/';
% chipodinfo.basedir.tao15_140 = '../../chipod/tao15_140/data/';

