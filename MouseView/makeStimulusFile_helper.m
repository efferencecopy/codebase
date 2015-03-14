fin

params.name = 'Train_[40_60_100]_SR20kHz_TT12500ms_PW300us.atf';  % the name of the output .atf file
params.type = 'train';            % 'train', 'pulse'
params.si   = 50e-6;              % the sample INTERVAL (needs to be an iteger)
params.swpDur = 25000;            % The total duration of the sweep IN NUMBERS OF SAMPLES!!!!
params.tStart = 0.250;            % the time of the first pulse
params.pAmp = [1];          % A vector of amplitudes for the pulse height [interleaved variable]
params.pFreq = [40, 60, 100];              % A vector of frequencies for the pulse train [interleaved variable]
params.nPulses = 7;               % The number of pulses in the pulse train. The code will throw a warning if this can not be done in the time aloted

params.pWidth = 300e-6;           % A vector of pulse widths (in seconds)  [interleaved variable]
params.nReps = 46;                % Number of repeates each stimulus should be presented



makeStimulusFile(params);