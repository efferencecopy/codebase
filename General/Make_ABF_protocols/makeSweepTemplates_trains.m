function params = makeSweepTemplates_trains(params)

% params should have
%
% params.name           =>  the name of the output .atf file
% params.type           =>  'train', 'pulse'
% params.si             =>  the sample INTERVAL (needs to be an iteger)
% params.swpDur         =>  The total duration of the sweep IN NUMBERS OF SAMPLES!!!!
% params.tStart         =>  the time of the first pulse
% params.pAmp           =>  A vector of amplitudes for the pulse height [interleaved variable]
% params.pFreq          =>  A vector of frequencies for the pulse train [interleaved variable]
% params.nPulses        =>  The number of pulses in a pulse train
% params.pWidth         =>  A vector of pulse widths (in seconds)  [interleaved variable]
% params.recovery       =>  true/false, should there be a recovery pulse?
% params.recoveryTime   =>  A vector of numbers corresponding to the recovery time in seconds [interleaved variable]
% params.nReps          =>  Number of repeates each stimulus should be presented



%
% Generate the stimulus waveforms (one for each unique stimulus type)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%
nAmps = numel(params.pAmp);
nFreqs = numel(params.pFreq);
nPulseWidths = numel(params.pWidth);
nRecoveryTimes = numel(params.recoveryTime);

tt = [0:params.swpDur-1]' .* params.si;
tStartIdx = ceil(params.tStart ./ params.si);
conditions = fullfact([nAmps, nFreqs, nPulseWidths, nRecoveryTimes]);

% initalize the outputs
params.templates_trains = repmat({nan(numel(tt), 1)}, 1, size(conditions, 1));
params.conditions_trains = conditions; % needs to be updated on a cond by cond basis
params.ttype_header_trains = {'pAmp', 'pFreq', 'pWidth', 'tRecov'};

% loop over the conditions and construct the waveform for each sweep
for i_cond = 1:size(conditions, 1)
    
    % zero out everything up to the sample before the first pulse
    params.templates_trains{i_cond}(1:tStartIdx-1) = 0;
    
    % make a pulse "motif" based on the width of the pulse and the
    % amplitude
    tmp_pAmp = params.pAmp(conditions(i_cond, 1));
    tmp_pFreq = params.pFreq(conditions(i_cond, 2));
    tmp_pWidth = params.pWidth(conditions(i_cond, 3));
    samplesPerPulse = ceil(tmp_pWidth ./ params.si); 
    
    % basic error checking
    assert(tmp_pAmp<=10, 'ERROR: pulse amp > 10 volts');
    assert(tmp_pWidth<=1, 'ERROR: pulse amp > 1 second');
    
    if strcmpi(params.type, 'pulse') || tmp_pFreq == 0 % only one pulse
        params.templates_trains{i_cond}(tStartIdx : tStartIdx+samplesPerPulse-1) = tmp_pAmp;
        params.templates_trains{i_cond}(tStartIdx+samplesPerPulse : end) = 0;
        
    else % multiple pulses
        samplesPerPeriod = ceil(1./tmp_pFreq ./ params.si);
        samplesPerIPI = samplesPerPeriod - samplesPerPulse;
        motif = zeros(1,samplesPerPeriod);
        motif(1:samplesPerPulse) = tmp_pAmp;
        
        % quick error checking
        samplesPerTrain = samplesPerPeriod .* params.nPulses;
        idx_trainEnd = samplesPerTrain + tStartIdx;
        assert(idx_trainEnd < numel(tt), 'ERROR: This pulse train can not fit in the sweep-time specified')
        
        % now construct the train pulse by pulse
        idx = tStartIdx;
        for i_pulse = 1:params.nPulses
            params.templates_trains{i_cond}(idx:idx+samplesPerPeriod-1) = motif;
            idx = idx + samplesPerPeriod;
        end
        params.templates_trains{i_cond}(idx:end) = 0; % add the trailing zeros
        
        % add a recovery pulse if desired.
        if params.recovery
            recoveryTime = params.recoveryTime(conditions(i_cond, 4));
            recoveryTimeInSamps = ceil(recoveryTime./params.si) + 1;
            idx = idx - samplesPerIPI + recoveryTimeInSamps;
            
            assert(idx + samplesPerPulse < numel(tt), 'ERROR: Train and recovery pulse can not fit in the sweep-time specified')
            
            params.templates_trains{i_cond}(idx:idx+samplesPerPulse-1) = tmp_pAmp;
        else
            recoveryTime = nan;
        end
        
        
        % define the per-condition ttype array
        params.conditions_trains(i_cond,:) = [tmp_pAmp, tmp_pFreq, tmp_pWidth, recoveryTime];
    
    end

end

