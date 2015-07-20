fin

params.name = 'Trains_amps[4_6_8_10]_SR20kHz_TT1250ms_PW300us';  % the name of the output .atf file
params.type = 'train';            % 'train', 'pulse'
params.si   = 50e-6;              % the sample INTERVAL (needs to be an iteger)
params.swpDur = 25000;            % The total duration of the sweep IN NUMBERS OF SAMPLES!!!!
params.tStart = 0.250;            % the time of the first pulse
params.pAmp = [4,6,8,10];                % A vector of amplitudes for the pulse height [interleaved variable]
params.pFreq = [50];              % A vector of frequencies for the pulse train [interleaved variable]
params.nPulses = 7;               % The number of pulses in the pulse train. The code will throw a warning if this can not be done in the time aloted

params.pWidth = 300e-6;           % A vector of pulse widths (in seconds)  [interleaved variable]
params.nReps = 45;                % Number of repeates each stimulus should be presented

params.recovery = false;       %  true/false, should there be a recovery pulse?
params.recoveryTime = [0.25 0.5 0.75];   %  A vector of numbers corresponding to the recovery time in seconds [interleaved variable]

makeStimulusFile(params);

