%% MAKE SOME PULSE TRAINS: RECOVERY, RIT, ENVELOPED-RIT

fin

%
% Define the params that influence the entire data file (trains and RITs)
%
params.si   = 50e-6;              % the sample INTERVAL (needs to be an iteger)
params.swpDur = 20e4;             % The total duration of the sweep IN NUMBERS OF SAMPLES!!!!
params.tStart = 0.100;            % the time of the first pulse
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
% make the sweep template(s) for the Random impulse trains with an envelope
%
params.ritFreq = [4,10,20];
params.ritHiFreqCut = 58;  % ISIs faster than this will be cutout
params.rit_Nversions = 1;
params.ritEnvelopeFreq = [1,2,4,8,16];
params = makeSweepTemplates_poiss(params); % templates are stored in params.templates_poiss



% concatenate templates
sweepTemplates = cat(2, params.templates_trains, params.templates_poiss);

