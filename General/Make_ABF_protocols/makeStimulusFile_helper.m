fin

params.name = 'Recovery_Train_[10_20_40Hz]_SR20kHz_TT3725ms_PW1000us.atf';  % the name of the output .atf file
params.type = 'train';            % 'train', 'pulse'
params.si   = 50e-6;              % the sample INTERVAL (needs to be an iteger)
params.swpDur = 74500;            % The total duration of the sweep IN NUMBERS OF SAMPLES!!!!
params.tStart = 0.250;            % the time of the first pulse
params.pAmp = [1];                % A vector of amplitudes for the pulse height [interleaved variable]
params.pFreq = [10, 20, 40];      % A vector of frequencies for the pulse train [interleaved variable]
params.nPulses = 7;               % The number of pulses in the pulse train. The code will throw a warning if this can not be done in the time aloted

params.pWidth = 1e-3;           % A vector of pulse widths (in seconds)  [interleaved variable]
params.nReps = 10;                % Number of repeates each stimulus should be presented

params.recovery = true;       %  true/false, should there be a recovery pulse?
params.recoveryTime = [0.25 0.5 0.75];   %  A vector of numbers corresponding to the recovery time in seconds [interleaved variable]

makeStimulusFile(params);

