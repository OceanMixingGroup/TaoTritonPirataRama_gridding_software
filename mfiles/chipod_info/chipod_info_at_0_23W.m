% chipod_info_at_0_23W
%
% In this module, specify necessary information about chipods at each
% mooring location.
%
% Sally Warner

%% define depths

% 'depths' defines the depths which you would like to appear in the final
% gridded output. Essentially, this vector should have all the depths at
% which chipods have ever been deployed at this location.
chipodinfo.depths = [21 30 35 50 65 81];



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

% Pirata14_23
% 525, 	/data/chipod/Pirata14/,                 -23.0,	30.0

% Pirata15_23
% 718, 	/data/chipod/Pirata15/chipods/23w/718/,     -23.0,	21.0 (Johannes processed with only internal stratification)
% 719, 	/data/chipod/Pirata15/chipods/23w/719/,     -23.0,	35.0 (Johannes processed with only internal stratification)
% 720, 	/data/chipod/Pirata15/chipods/23w/720/,     -23.0,	50.0 (Johannes processed with only internal stratification)
% 721, 	/data/chipod/Pirata15/chipods/23w/721/,     -23.0,	65.0 (Johannes processed with only internal stratification)
% 722, 	/data/chipod/Pirata15/chipods/23w/722/,     -23.0,	81.0 (Johannes processed with only internal stratification)

% Pirata16_23
% 1107, /data/chipod/Pirata16/data/23w/1107/,	-23.0,	21.0
% 1108, /data/chipod/Pirata16/data/23w/1108/,	-23.0,	35.0
% 1109, /data/chipod/Pirata16/data/23w/1108/,	-23.0,	50.0
% 1110, /data/chipod/Pirata16/data/23w/1108/,	-23.0,	65.0
% 1111, /data/chipod/Pirata16/data/23w/1108/,	-23.0,	81.0

% Pirata17_23
% 721   /data/chipod/Pirata17/data/23w/721/,	-23.0,	21.0
% 722   /data/chipod/Pirata17/data/23w/722/,	-23.0,	35.0
% 723   /data/chipod/Pirata17/data/23w/723/,	-23.0,	50.0
% 724   /data/chipod/Pirata17/data/23w/724/,	-23.0,	65.0
% 725   /data/chipod/Pirata17/data/23w/725/,	-23.0,	81.0

% Pirata18_23
% 1101   /data/chipod/Pirata18/data/23w/1101/,	-23.0,	21.0
% 1113   /data/chipod/Pirata18/data/23w/1113/,	-23.0,	35.0
% 1114   /data/chipod/Pirata18/data/23w/1114/,	-23.0,	50.0
% 1115   /data/chipod/Pirata18/data/23w/1115/,	-23.0,	65.0
% 1128   /data/chipod/Pirata18/data/23w/1128/,	-23.0,	81.0

chipodinfo.dpl = {'Pirata14';'Pirata15';'Pirata16';'Pirata17';'Pirata18'};
chipodinfo.newdepths = [21    30    35    50    65    81];
chipodinfo.cpds =      [0     525   0     0     0     0;     % Pirata14
                        718   0   	719   720   721   722;   % Pirata15
                        1107  0     1108  1109  1110  1111;  % Pirata16
                        721   0     722   723   724   725;   % Pirata17
                        1101  0     1113  1114  1115  1128]; % Pirata18


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
chipodinfo.basedir.Pirata14 = '../../chipod/Pirata14/chipods/23w/';
chipodinfo.basedir.Pirata15 = '../../chipod/Pirata15/chipods/23w/';
chipodinfo.basedir.Pirata16 = '../../chipod/Pirata16/data/23w/';
chipodinfo.basedir.Pirata17 = '../../chipod/Pirata17/data/23w/';
chipodinfo.basedir.Pirata18 = '../../chipod/Pirata18/data/23w/';

