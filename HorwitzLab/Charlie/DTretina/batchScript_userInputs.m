% This script can be run on any computer (including the shadlen cluster).
% And is designed to probe many color directions, which in this case, are
% specified in this script (as opposed to a dtnt run). I'm also going to
% experiment with running parfor loops on the shadlen cluster.

% clear out the work space
fin

    
%
% define the necessary parameters for the DTcones simulation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.runType = 'dtnt';                         % 'dtnt', or 'absThresh'
params.obsMethod = 'obsMethod_filteredWtFxn';     % 'obsMethod_noClrEqSpace' or 'obsMethod_absThresh' or 'obsMethod_phaseInvariant' or 'obsMethod_filteredWtFxn'
params.Ncones = NaN;                             % set to NaN, except when using the absThresh functionality
params.monCalFile = 'DTcals_equal_bkgnd_1200.mat';                % 'DTcals.mat', 'DTcals_100Hz_framerate.mat', 'DTcals_825Hz_framerate.mat', 'DTcals_apollo_macpig.mat',  or 'DTcals_equal_bkgnd.mat'
params.impulseResponse = 'rieke';                %  'rieke', or 'deltafxn'
params.flatPowerSpect = false;                   % if true, PS is flat w/same integral as the normal PS.
params.enableScones = true;                      % should the S-cones contribute to the pooled response?
params.eyeType = 'monkey';                       % 'monkey' or 'human'
params.coneSampRate = 825;                       % good candidates: [525 600 675 750 825 900 975] These all give rise to nearly an iteger number of 'cone' sampels per monitor refresh
params.colorselection = 'lots';              % could be: 'lots', 'specific', 'guniso'

% define some helpful text files (if necessary), and the paramaters for
% parallel operations
params.DTNT_fname = '';                       % will be created at run time and dynamically on each loop
params.aperatureMosaic = false;               % only for DTV1 experiments
params.saveDir = '';                          % modified below;
params.parallelOperations = true;             % needs to be true for this script (even if parfor is unavailable) so that DTcones knows how to save temporary files

% enable some debugging options if necesary
params.unitTest = false;             % true or false
params.eqMosaic = false;             % for debugging. true or false

% make some notes...
params.notes = 'equal bkgnd for Greg, but normal cone mosaic. 200R*/sec bkgnd';       % notes that should be associated the data file?





%
% In order to probe a bunch of color directions, set up a 'dtnt' struct and
% hijack DTcones ability to import all the relavent gabor parameters.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dtnt.rf_x = -50;     % in tenths of dva
dtnt.rf_y = -35;       % in tenths of dva
dtnt.sigma = 4;    % in tenths of dva
dtnt.nSD = 3;        % number of SDs in the gabor (extends nSD in either direction)
dtnt.theta = 0;
dtnt.gamma = 1;
dtnt.length = .666;  % in seconds
dtnt.speed = 3;
dtnt.sfs = 1;
dtnt.alphas = []; % gets filled in later
dtnt.colorDirs = []; % gets filled in later


if strcmpi(params.colorselection , 'lots')
    
    
    % define the color directions.... just putting points on a sphere for now
    nColors = 100;
    tmp = ceil(sqrt(nColors));
    az = linspace(0, 2*pi, tmp);
    el = linspace(0, pi/2, tmp);
    inds = fullfact([tmp, tmp]);
    dirs_pol = [az(inds(:,1))', el(inds(:,2))'];
    [x,y,z] = sph2cart(dirs_pol(:,1), dirs_pol(:,2), ones(size(dirs_pol,1),1));
    z = z.*3; % scale the s-cone axis
    colorDirs = [x,y,z];
    norms = sqrt(sum(colorDirs.^2,2));
    colorDirs = bsxfun(@rdivide, colorDirs, norms);
    colorDirs(abs(colorDirs)<1e-14) = 0;
    colorDirs = unique(colorDirs, 'rows'); %remove the duplicates
    nColors = size(colorDirs,1);
    alphas = ones(nColors,1) .* .5; % all the same alpha estimate
    
elseif strcmpi(params.colorselection , 'specific')
    
    % In case I want to test a specific set of colors...
    colorDirs = [1 0 0; 0 1 0; 0 0 1];
    colorDirs = bsxfun(@rdivide, colorDirs, sqrt(sum(colorDirs.^2, 2))); % as unit vectors
    alphas = [1, 1, 1];
    nColors = size(colorDirs,1);
    
elseif strcmpi(params.colorselection , 'guniso')
    
    % In case I want to test the gun-iso color directions
    colorDirs = ['ggun'; 'bgun'];
    alphas = [0.98, 0.98];
    nColors = size(colorDirs,1);
end


%
%  CREATE THE NECESSARY TEMPORARY FILES
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create a temporary directory to hold all the data files that will get
% stored as a result of the parfor operation
switch whoami
    case 'hass_mbp'
        params.saveDir = '/Users/charliehass/LabStuff/Huskies/DTcones/Data/tmpBatchData_DTNT/';
    case 'glick_rig1'
        params.saveDir = 'C:\Users\glickfeld_lab\Desktop\Local Data Files\DTcones\tmpBatchData_DTNT\';
    case 'nuke'
        params.saveDir = 'C:\Users\charlie\Desktop\Local Data\DTcones\tmpBatchData_DTNT\';
end
mkdir(params.saveDir);
cd(params.saveDir);

% save a bunch of .mat files and param structs for use during parfor
[pstructs{1:nColors}] = deal(params);
for a = 1:nColors;
    
    % make and save a temp params struct (save in a cell array)
    pstructs{a}.DTNT_fname = [params.saveDir, 'dtnt_batch_', num2str(a)];
    
    % make and save a temp dtnt file
    if ischar(colorDirs(a,:))
        tmpColors = colorDirs(a,:);
    else
        tmpColors = colorDirs(a,:) ./ norm(colorDirs(a,:));
    end
    dtnt.colorDirs = tmpColors;
    dtnt.alphas = alphas(a);
    save(['dtnt_batch_', num2str(a)], 'dtnt');
    
end


%
%
%  PARALLEL FOR LOOP OVER COLOR DIRECTIONS
%
% NOTE: 'parfor' is not guaranteed to iterate in monotonically icreasing
% iteration numbers (e.g., 1:1:N). This means that the order of color
% directions (listed above) might be different than the order of data files
% saved. Since the data files are named based on the time they're saved,
% the order of color directions in the batch data file could be different
% than the order specified by 'colorDirs'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%open a matlabpool
if exist('matlabpool', 'file') == 2;
    poolObj = parpool([2 nColors]);
    %matlabpool open 16
    pause(2)
    fprintf(' *** Using parallel operations *** \n')
end

parfor a = 1:nColors
    % run the simulation
    DTcones(pstructs{a})
end

% close the workers
if exist('matlabpool', 'file') == 2;
    %matlabpool('close')
    delete(poolObj);
end



%
%  REPACKAGE THE DATA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Repackaging the data\n')

% Cycle through the data files and concatenate them into a single set of
% structures.
dataFilePath = [params.saveDir, 'dtnt_batch', '_%d', filesep, 'out_', 'dtnt_batch', '_%d.mat'];

% A hack to make things work on PCs (also, Fileseps on PC are annoying)
idx = regexp(dataFilePath, filesep);
tmpPath = dataFilePath;
tmpPath(idx) = '&'; %changing to a nonsense character that will go through the sprintf statement....
tmpPath = sprintf(tmpPath, 1, 1);
idx = regexp(tmpPath, '&');
tmpPath(idx) = filesep;
% END HACK

load(tmpPath);
nContrasts = numel(gab.contrasts{1});
tmp_colorDirs = [];
tmp_contrasts = repmat({nan(1,nContrasts)}, 1, nColors);

% define the dimensionality of the temporary arrays
warning('need to deal with this!')
switch params.obsMethod
    case 'obsMethod_all'
        [tmp_analyticMean, tmp_analyticVar] = deal(nan(nColors, nContrasts)); % dimensionality will be nColors x nContrasts
        tmp_resp = nan(nColors, nContrasts, gab.nTrials); % dimensionality will be (nColors x nContrasts x nTrials)
    case 'obsMethod_noClrEqSpace'
        [tmp_analyticMean(1:nColors, 1:nContrasts)] = deal({nan(1,3)});
        [tmp_analyticVar(1:nColors, 1:nContrasts)] = deal({nan(1,3)});
        [tmp_resp(1:nColors, 1:nContrasts, 1:gab.nTrials)] = deal({nan(1,3)});
end


for a = 1:nColors
    
    % A hack to make things work on PCs (also, Fileseps on PC are annoying)
    idx = regexp(dataFilePath, filesep);
    tmpPath = dataFilePath;
    tmpPath(idx) = '&'; %changing to a nonsense character that will go through the sprintf statement....
    tmpPath = sprintf(tmpPath, a, a);
    idx = regexp(tmpPath, '&');
    tmpPath(idx) = filesep;
    % END HACK
    
    load(tmpPath);
    tmp_analyticMean(a,:) = idlob.analyticMean;
    tmp_analyticVar(a,:) = idlob.analyticVar;
    tmp_resp(a,:,:) = idlob.resp;
    tmp_colorDirs = cat(1,tmp_colorDirs, gab.colorDirs);
    tmp_contrasts(a) = gab.contrasts;
end

% The data files differ only in color direction, so use the last one as a
% template, but fill in the important things here:
gab.colorDirs = tmp_colorDirs;
gab.contrasts = tmp_contrasts;
idlob.resp = tmp_resp;
idlob.analyticMean = tmp_analyticMean;
idlob.analyticVar = tmp_analyticVar;

% delete all the temporary files and folders
cd([params.saveDir, '/../']) % cd one level up
% rmdir(params.saveDir, 's') % doesn't work on the glick lab PC




%
% save the data
%
%%%%%%%%%%%%%%%%%%%%

fprintf('  Saving the batch data file...\n')

% grab the contents of the .m file used to generate the data (i.e.,
% this file you're reading now)
pathToDTcones = which('DTcones');
fid = fopen(pathToDTcones, 'r');
mfiletext = fscanf(fid, '%c');
fclose(fid);

%cd to where the data should be stored
originalDir = pwd;
switch whoami
    case 'hass_mbp'
        newDir = '/Users/charliehass/LabStuff/Huskies/DTcones/Data';
    case 'glick_rig1'
        newDir = 'C:\Users\glickfeld_lab\Desktop\Local Data Files\DTcones\Data';
    case 'nuke'
        newDir = 'C:\Users\charlie\Desktop\Local Data\DTcones\Data';
end
cd(newDir)

%find a directory, or make one if necessary
d = dir;
suffix = 1;
timestamp = [date, sprintf('_%d', suffix)];
alreadyExists = any(strcmpi(timestamp, {d.name}));
while alreadyExists
    suffix = suffix+1;
    timestamp = [date, sprintf('_%d', suffix)];
    alreadyExists = any(strcmpi(timestamp, {d.name}));
end

%save the data
mkdir(pwd, timestamp)
cd(timestamp)
fname_dat = ['out_', timestamp, '.mat'];
eval(['save ' fname_dat, ' gab cones mon idlob params'])

%save the code of the calling function (i.e., DTcones.m)
fname_mfile = ['DTcones.m_', timestamp];
fid = fopen(fname_mfile, 'w');
fwrite(fid, mfiletext, 'char');
fclose(fid);

%be nice and cd back to where you came from
cd(originalDir)

%%%%%%%%%%%%%%%%%%%%%%%
%
% END SAVE DATA
%

fprintf('   **** All done with the batch process **** \n')

