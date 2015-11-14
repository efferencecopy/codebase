fin

params.name = 'Trains_2_[10_20_40_60_80_100Hz]_SR100kHz_TT1500ms_PW200us.atf';  % the name of the output .atf file
params.type = 'train';            % 'train', 'pulse'
params.si   = 10e-6;              % the sample INTERVAL (needs to be an iteger)
params.swpDur = 150000;            % The total duration of the sweep IN NUMBERS OF SAMPLES!!!!
params.tStart = 0.252;            % the time of the first pulse
params.pAmp = 5;                % A vector of amplitudes for the pulse height [interleaved variable]
params.pFreq = [10,20,40,60,80,100];              % A vector of frequencies for the pulse train [interleaved variable]
params.nPulses = 10;               % The number of pulses in the pulse train. The code will throw a warning if this can not be done in the time aloted

params.pWidth = 200e-6;           % A vector of pulse widths (in seconds)  [interleaved variable]
params.nReps = 5;                % Number of repeates each stimulus should be presented

params.recovery = false;       %  true/false, should there be a recovery pulse?
params.recoveryTime = [0.25 0.5 0.75];   %  A vector of numbers corresponding to the recovery time in seconds [interleaved variable]

makeStimulusFile(params);

