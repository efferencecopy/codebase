% This script can be run on any computer (including the shadlen cluster)
% and simply runs a simple simulation using the DTcones model.


% clear out the work space
fin


% define the experiment parameters
params.runType = 'default';                      % 'DTVT', 'DTNT', 'absThresh' or 'default'
params.obsMethod = 'obsMethod_filteredWtFxn';     % 'obsMethod_noClrEqSpace' or 'obsMethod_absThresh' or 'obsMethod_phaseInvariant' or 'obsMethod_filteredWtFxn'
params.Ncones = nan;                             % set to NaN, except when using the absThresh functionality
params.monCalFile = 'DTcals.mat';                % 'DTcals.mat', 'DTcals_100Hz_framerate.mat', 'DTcals_825Hz_framerate.mat' or 'DTcals_equal_bkgnd.mat'
params.impulseResponse = 'rieke';                % params.impulseResponse    => 'rieke', or 'deltafxn'
params.flatPowerSpect = false;
params.enableScones = true;                      % should the S-cones contribute to the pooled response?
params.eyeType = 'monkey';                       % 'monkey' or 'human'
params.coneSampRate = 825;                       % good candidates: [525 600 675 750 825 900 975] These all give rise to nearly an iteger number of 'cone' sampels per monitor refresh

% define some helpful text files (if necessary), and the paramaters for
% parallel operations
params.DTV1_fname = '';             % the name of a .nex file for a DT expt
params.DTNT_fname = '';             % The name of a .mat file that has batch data from all the DTNT runs
params.saveDir = '~/Desktop';       % where should the data be saved?
params.parallelOperations = false;

% enable some debugging options if necesary
params.unitTest = true;             % true or false
params.eqMosaic = false;             % for debugging. true or false
params.aperatureMosaic = false;      % only for DTV1 experiments

% make some notes... 
params.notes = '';                  % notes that should be associated the data file?



% run the simulation
DTcones(params)
