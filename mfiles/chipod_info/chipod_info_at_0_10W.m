% chipod_info_at_0_10W
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

% Pirata14_10
% 524,      /data/chipod/Pirata14/,                 30.0

% Pirata15_10
% 709,      /data/chipod/Pirata15/chipods/10w/709/, 21.0 (Johannes only processed with internal stratification)
% 714,      /data/chipod/Pirata15/chipods/10w/714/, 35.0 (ALL DATA BAD)
% 715,      /data/chipod/Pirata15/chipods/10w/715/, 50.0 (Johannes only processed with internal stratification)
% 716,      /data/chipod/Pirata15/chipods/10w/716/, 65.0 (Johannes only processed with internal stratification)
% 717,      /data/chipod/Pirata15/chipods/10w/717/, 81.0 (Johannes only processed with internal stratification)

% Pirata16_10
% 1101, 	/data/chipod/Pirata16/data/10w/1101/,	21.0 (Johannes only processed with internal stratification)
% 1102, 	/data/chipod/Pirata16/data/10w/1102/,	35.0 (Johannes only processed with internal stratification)
% 1104, 	/data/chipod/Pirata16/data/10w/1104/,	50.0 (only 13 days; processed with both Tz_i and Tz_m by Sally in June 2020)
% 1105, 	/data/chipod/Pirata16/data/10w/1105/,	65.0 (only 9 days; processed with both Tz_i and Tz_m by Sally in June 2020)
% 1106, 	/data/chipod/Pirata16/data/10w/1106/,	81.0 (Johannes only processed with internal stratification)

% Pirata17_10
% 710       /data/chipod/Pirata17/data/10w/710/,    21
% 711       /data/chipod/Pirata17/data/10w/711/,    35
% 712       /data/chipod/Pirata17/data/10w/712/,    50 (ALL DATA BAD)
% 713       /data/chipod/Pirata17/data/10w/713/,    65
% 714       /data/chipod/Pirata17/data/10w/714/,    81 (VC BAD; IC GOOD)

% Pirata18_10
% 1123       /data/chipod/Pirata18/data/10w/1123/,    21
% 1124       /data/chipod/Pirata18/data/10w/1124/,    35
% 1125       /data/chipod/Pirata18/data/10w/1125/,    50 (NO GOOD DATA)
% 1126       /data/chipod/Pirata18/data/10w/1126/,    65
% 1127       /data/chipod/Pirata18/data/10w/1127/,    81 

chipodinfo.dpl = {'Pirata14';'Pirata15';'Pirata16';'Pirata17';'Pirata18'};
chipodinfo.newdepths = [21    30    35    50    65    81];
chipodinfo.cpds      = [0     524   0     0     0     0;     % Pirata14
                        709   0     0     715   716   717;   % Pirata15
                        1101  0     1102  1104  1105  1106;  % Pirata16
                        710   0     711   0     713   714;   % Pirata17
                        1123  0     1124  0     1126  1127];  % Pirata18
     % For chipodinfo.cpds, each column corresponds to each depth in
     % chipodinfo.newdepths, and each row corresponds to each deployment in
     % chipodinfo.dpl. Make sure orders stay consistent. Fill with unit #s.

     
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
chipodinfo.basedir.Pirata14 = '../../chipod/Pirata14/chipods/10w/';
chipodinfo.basedir.Pirata15 = '../../chipod/Pirata15/chipods/10w/';
chipodinfo.basedir.Pirata16 = '../../chipod/Pirata16/data/10w/';
chipodinfo.basedir.Pirata17 = '../../chipod/Pirata17/data/10w/';
chipodinfo.basedir.Pirata18 = '../../chipod/Pirata18/data/10w/';


