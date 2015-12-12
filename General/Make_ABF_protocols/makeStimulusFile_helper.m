fin

params.name = 'Trains_[10_to_120]_SR20kHz_TT1250ms_PW100us_three';  % the name of the output .atf file
params.type = 'trains';            % 'train', 'pulse'
params.si   = 50e-6;              % the sample INTERVAL (needs to be an iteger)
params.swpDur = 25000;            % The total duration of the sweep IN NUMBERS OF SAMPLES!!!!
params.tStart = 0.250150;            % the time of the first pulse
params.pAmp = 5.5;                % A vector of amplitudes for the pulse height [interleaved variable]
params.pFreq = [10, 20, 40, 60, 100, 120];              % A vector of frequencies for the pulse train [interleaved variable]
params.nPulses = 7;               % The number of pulses in the pulse train. The code will throw a warning if this can not be done in the time aloted

params.pWidth = 100e-6;           % A vector of pulse widths (in seconds)  [interleaved variable]
params.nReps = 23;                % Number of repeates each stimulus should be presented

params.recovery = false;       %  true/false, should there be a recovery pulse?
params.recoveryTime = [0];   %  A vector of numbers corresponding to the recovery time in seconds [interleaved variable]


%% randomize trials before running the funtion (if desired)

% Do the trial randomization outsid the makeStimulus function so that the
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

%% generate the stimuli

makeStimulusFile(params);

