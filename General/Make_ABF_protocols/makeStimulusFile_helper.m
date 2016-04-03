%% SIMPLE PULSE TRAINS

fin

params.name = 'Trains_[10_to_60]_SR20kHz_TT2000ms_PW300us';  % the name of the output .atf file
params.type = 'trains';            % 'train', 'pulse'
params.si   = 50e-6;              % the sample INTERVAL (needs to be an iteger)
params.swpDur = 40000;            % The total duration of the sweep IN NUMBERS OF SAMPLES!!!!
params.tStart = 0.500;            % the time of the first pulse
params.pAmp = 1;                % A vector of amplitudes for the pulse height [interleaved variable]
params.pFreq = [10, 20, 40, 60];              % A vector of frequencies for the pulse train [interleaved variable]
params.nPulses = 10;               % The number of pulses in the pulse train. The code will throw a warning if this can not be done in the time aloted

params.pWidth = 350e-6;           % A vector of pulse widths (in seconds)  [interleaved variable]
params.nReps = 25;                % Number of repeates each stimulus should be presented

params.recovery = false;       %  true/false, should there be a recovery pulse?
params.recoveryTime = [0];   %  A vector of numbers corresponding to the recovery time in seconds [interleaved variable]

%
% randomize trials before running the funtion (if desired)
%
% Do the trial randomization outside the makeStimulus function so that the
% randomization can be consistent across multiple files (useful for making
% bi-phasic stimuli on the iso-flex units)
nAmps = numel(params.pAmp);
nFreqs = numel(params.pFreq);
nPulseWidths = numel(params.pWidth);
nRecoveryTimes = numel(params.recoveryTime);
conditions = fullfact([nAmps, nFreqs, nPulseWidths, nRecoveryTimes]);
trlTypes = 1:size(conditions,1);
trlTypes = repmat(trlTypes, 1, params.nReps);
randIdx = randperm(numel(trlTypes)); % randomize the order
trlTypes = trlTypes(randIdx);
params.trlTypes = trlTypes;
params.nSweeps = numel(trlTypes);

%
% generate the .atf file
%

params = makeSweepTemplates_trains(params);
makeAxonTextFile(params, params.templates_trains);



%% 10 SECOND RECOVERY TRAINS + RITs

fin

%
% Define the params that influence the entire data file (trains and RITs)
%
params.name = 'WCSTP_SR20kHz_TT11s_PW300us_';  % the name of the output .atf file
params.si   = 50e-6;              % the sample INTERVAL (needs to be an iteger)
params.swpDur = 20e4;            % The total duration of the sweep IN NUMBERS OF SAMPLES!!!!
params.tStart = 0.500;            % the time of the first pulse
params.pAmp = 1;                  % A vector of amplitudes for the pulse height [interleaved variable]
params.pWidth = 350e-6;           % A vector of pulse widths (in seconds)  [interleaved variable]





%
% make the sweep templates for the pulse trains
%
params.type = 'trains';            % 'train', 'pulse'
params.pFreq = [10, 20, 40];              % A vector of frequencies for the pulse train [interleaved variable]
params.nPulses = 10;
params.recovery = true;       %  true/false, should there be a recovery pulse?
params.recoveryTime = [0.5, 2, 8];   %  A vector of numbers corresponding to the recovery time in seconds [interleaved variable]

params = makeSweepTemplates_trains(params); % templates are stored in params.templates_trains


%
% make the sweep template(s) for the Random impulse trains
%
params.ritFreq = 20;
params.ritHiFreqCut = 58;  % ISIs faster than this will be cutout

params = makeSweepTemplates_poiss(params); % templates are stored in params.templates_poiss



% concatenate templates
allSweepTemplates = cat(2, params.templates_trains, params.templates_poiss);

% present 2 blocks, randomize order within block. Define a block as 1
% repeat of each recovery train, and 3 repeats of the RIT.
recoveryTrain_idx = 1:numel(params.templates_trains);
poissTrain_idx = (1:numel(params.templates_poiss)) + max(recoveryTrain_idx);
blockIdx = [recoveryTrain_idx, repmat(poissTrain_idx, 1, 3)];
trlTypes = [];
for i_block = 1:5;
    randIdx = randperm(numel(blockIdx));
    trlTypes(:,i_block) = blockIdx(randIdx); % each column is a block
end
params.trlTypes = trlTypes(:);
params.nSweeps = numel(trlTypes);

% make the .stf file
makeAxonTextFile(params, allSweepTemplates);





