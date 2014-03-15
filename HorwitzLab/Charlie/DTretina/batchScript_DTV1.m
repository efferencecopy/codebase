% This script can be run on any computer (including the shadlen cluster).
% And is designed to probe many color directions, which in this case, are
% specified in this script (as opposed to a dtnt run). I'm also going to
% experiment with running parfor loops on the shadlen cluster.

% clear out the work space
fin


% define the experiment parameters
params.runType = 'dtv1';                         % 'DTVT', 'DTNT', 'absThresh' or 'default'
params.obsMethod = 'obsMethod_filteredWtFxn';    % obsMethod_noClrEqSpace or obsMethod_absThresh or obsMethod_filteredWtFxn
params.Ncones = nan;                             % set to NaN, except when using the absThresh functionality
params.monCalFile = '';                          % DTcones will use the calibration provided by DTV1 files.
params.impulseResponse = 'rieke';                % params.impulseResponse    => 'rieke', or 'deltafxn'
params.flatPowerSpect = false;
params.enableScones = true;                      % should the S-cones contribute to the pooled response?
params.eyeType = 'monkey';                       % 'monkey' or 'human'
params.coneSampRate = 825;                       % good candidates: [525 600 675 750 825 900 975] These all give rise to nearly an iteger number of 'cone' sampels per monitor refresh
params.aperatureMosaic = true;                   % only for DTV1 experiments




% define some helpful text files (if necessary), and the paramaters for
% parallel operations
params.DTV1_fname = '';             % will be created at run time and dynamically on each loop
params.saveDir = '';                % will be created at run time
params.parallelOperations = true;

% enable some debugging options if necesary
params.unitTest = false;            % true or false
params.eqMosaic = false;            % for debugging. true or false

% make some notes... 
params.notes = '';                  % notes that should be associated the data file?


         
%
%  CREATE THE NECESSARY TEMPORARY FILES
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                  
         
% create a temporary directory to hold all the data files that will get
% stored as a result of the parfor operation
switch whoami
    case 'hass_mbp'
        params.saveDir = '/Users/charliehass/LabStuff/Huskies/DTcones/Data/tmpBatchData_DTV1/';
        DTfnames = fnamesFromTxt2('/Users/charliehass/LabStuff/Huskies/NexFiles/Charlie/Batch Data And Text Files/bothMonkDTPop.txt');
    case 'nuke'
        params.saveDir = 'C:\Users\charlie\Desktop\Local Data\DTcones\tmpBatchData_DTV1\';
        DTfnames = fnamesFromTxt2('\\crash.dhe.duke.edu\charlie\DTcones\Files for server\bothMonkDTPop.txt');
end
mkdir(params.saveDir);
cd(params.saveDir);

% save a bunch param structs for use during parfor
nExpts = numel(DTfnames);
[pstructs{1:nExpts}] = deal(params);
for a = 1:nExpts;
    % make and save a temp params struct (save in a cell array)
    pstructs{a}.DTV1_fname = DTfnames{a}{1};
    pstructs{a}.GTV1_fname = DTfnames{a}{2}; % for aperaturing the cone mosaic according to RF size
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
    %poolObj = parpool('local', 4);
    matlabpool open 6
    pause(2)
    fprintf(' *** Using parallel operations *** \n')
end

parfor a = 1:nExpts
    disp(a)
    
    % run the simulation
    pstructs{a}.DTV1_fname
    DTcones(pstructs{a})
    
    clc
end

% close the workers
if exist('matlabpool', 'file') == 2;
    matlabpool('close')
    %delete(poolObj);
end

% 
% Repackage the data in a way that is similar to the way DTV1 data is
% stored following a call to DTbatchCompile
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load the DTbatch data
switch whoami
    case 'hass_mbp' % charlie's laptop
        load '/Users/charliehass/LabStuff/Huskies/NexFiles/Charlie/Batch Data And Text Files/CardVsInt_06-May-2011.mat'
    case 'nuke'
        load '\\crash.dhe.duke.edu\charlie\DTcones\Files for server\CardVsInt_06-May-2011.mat'
end

ret = [];
for a = 1:numel(out.dat)
    
    
    % Open the appropriate DTcones data file
    expt = out.fnames{a}{1};
    cd([params.saveDir, expt])
    load(['out_', expt, '.mat'])
    ret.fnames{a} = expt;
    
    % Store some useful info
    ret.dat(a).colors = gab.colorDirs;
    if ~all(sign(ret.dat(a).colors) == sign(out.dat(a).expt.standColors)); error('colors are not in the same order'); end
    ret.dat(a).norms = gab.contrasts;
    ret.dat(a).respMean = idlob.analyticMean;
    ret.dat(a).respVar = idlob.analyticVar;
    
    % Compute the ROC values for each color/contrast
    nContrasts = size(idlob.resp,2);
    nColors = size(idlob.resp,1);
    idlob = coneNoiseROC(pstructs{a}, idlob, cones, gab);
    ret.dat(a).roc = idlob.roc;
    ret.dat(a).roc_analytic = idlob.roc_analytic;
    
    % Fit the neurometric functions with a cumulative Weibull. Store the
    % retinometric thresholds and slopes.
    [ret.dat(a).alpha, ret.dat(a).beta] = deal(nan(nColors,1));
    for clr = 1:nColors;
        
        % I don't really know how to fit a Weibull CDF to the analytic
        % data (ML fitting requires number of trials per condition...)
        tmp_norms = gab.contrasts{clr};
        tmp_roc = ret.dat(a).roc_analytic{clr};
        [~, idx] = min(abs(tmp_roc - 0.816));
        alphaGuess = tmp_norms(idx);
        [ret.dat(a).alpha(clr,1), ret.dat(a).beta(clr,1)] = weibullFit(tmp_norms, tmp_roc, 'sse', [alphaGuess 1]);
        
    end
    
end


% make a new directory that will store the two important data structure
% (the one from DTbatchCompile, and the one from the cone mosaic
% simulation). Save a representative .m file in case there are questions
% about what code I ran. Don't discard the original data.
%cd to where the data should be stored
switch whoami
    case 'hass_mbp' % charlie's laptop
        cd '/Users/charliehass/LabStuff/DTcones/Data'
    case 'nuke'
        cd 'C:\Users\charlie\Desktop\Local Data\DTcones\Data';
end

%find a directory, or make one if necessary
d = dir;
suffix = 1;
filename = ['DTV1_batch', sprintf('_%s_%d', date, suffix)];
alreadyExists = any(strcmpi(filename, {d.name}));
while alreadyExists
    suffix = suffix+1;
    filename = ['DTV1_batch', sprintf('_%s_%d', date, suffix)];
    alreadyExists = any(strcmpi(filename, {d.name}));
end

%save the data
mkdir(pwd, filename)
cd(filename)
save  batchData.mat  out ret params % save the batch data

% copy a sample .m file to the batch data directory
d = dir([params.saveDir, expt]); % from the last expt we analyzed in the for loop above
mfileidx = regexpi({d(:).name}, 'DTcones.m_');
mfileidx = cellfun(@(x) ~isempty(x), mfileidx);
copyfile([params.saveDir, expt, filesep, d(mfileidx).name], [pwd, filesep, d(mfileidx).name])


fprintf('   **** All done with the batch process **** \n')
