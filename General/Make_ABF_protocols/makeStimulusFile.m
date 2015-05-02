function makeStimulusFile(params)

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


% force default frequency to be zero (only one pulse)
if isempty(params.pFreq)
    params.pFreq = 0;
end

%
% Generate the stimulus waveforms (one for each unique stimulus type)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%
nAmps = numel(params.pAmp);
nFreqs = numel(params.pFreq);
nPulseWidths = numel(params.pWidth);
nRecoveryTimes = numel(params.recoveryTime);
nSweeps = nAmps .* nFreqs .* nPulseWidths .* nRecoveryTimes .* params.nReps;

tt = [0:params.swpDur-1]' .* params.si;
tStartIdx = ceil(params.tStart ./ params.si);
conditions = fullfact([nAmps, nFreqs, nPulseWidths, nRecoveryTimes]);

templates = repmat({nan(numel(tt), 1)}, 1, size(conditions, 1));

% loop over the conditions and construct the waveform for each sweep
for i_cond = 1:size(conditions, 1)
    
    % zero out everything up to the sample before the first pulse
    templates{i_cond}(1:tStartIdx-1) = 0;
    
    % make a pulse "motif" based on the width of the pulse and the
    % amplitude
    tmp_pAmp = params.pAmp(conditions(i_cond, 1));
    tmp_pFreq = params.pFreq(conditions(i_cond, 2));
    tmp_pWidth = params.pWidth(conditions(i_cond, 3));
    samplesPerPulse = ceil(tmp_pWidth ./ params.si) + 1; % need to add 1 b/c this is a 'fence post' problem
    
    % basic error checking
    assert(tmp_pAmp<=10, 'ERROR: pulse amp > 10 volts');
    assert(tmp_pWidth<=1, 'ERROR: pulse amp > 1 second');
    
    if tmp_pFreq == 0 % only one pulse
        templates{i_cond}(tStartIdx : tStartIdx+samplesPerPulse-1) = tmp_pAmp;
        templates{i_cond}(tStartIdx+samplesPerPulse : end) = 0;
        
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
            templates{i_cond}(idx:idx+samplesPerPeriod-1) = motif;
            idx = idx + samplesPerPeriod;
        end
        templates{i_cond}(idx:end) = 0; % add the trailing zeros
        
        % add a recovery pulse if desired.
        if params.recovery
            recoveryTime = params.recoveryTime(conditions(i_cond, 4));
            recoveryTimeInSamps = ceil(recoveryTime./params.si) + 1;
            idx = idx - samplesPerIPI + recoveryTimeInSamps;
            
            assert(idx + samplesPerPulse < numel(tt), 'ERROR: Train and recovery pulse can not fit in the sweep-time specified')
            
            templates{i_cond}(idx:idx+samplesPerPulse-1) = tmp_pAmp;
        end
        
    
    end
    
    
    

end


%
% create a matrix of sweeps to be exported to the atf file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%
trlTypes = 1:size(conditions,1);
trlTypes = repmat(trlTypes, 1, params.nReps);
randIdx = randperm(numel(trlTypes)); % randomize the order
trlTypes = trlTypes(randIdx);
sweeps = nan(numel(tt), numel(trlTypes));
for i_swp = 1:numel(trlTypes)
    sweeps(:,i_swp) = templates{trlTypes(i_swp)}(:);
end

% the atf file needs a time vector, so add that here as the first column
sweeps = cat(2, tt(:), sweeps);



%
% create the header information. Open a new file, and start writing into
% the new file.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%


header{1,:} = {'ATF', '1'};
header{2,:} = {'8', num2str(nSweeps)};
header{3} = {'AcquisitionMode=Episodic Stimulation'};
header{4} = {'Comment='};
header{5} = {'YTop=10'}; % theoretically could be different, but I think this is safe for now
header{6} = {'YBottom=-10'};
header{7} = {'SyncTimeUnits=0.4'};
header{8} = {'SweepStartTimesMS='}; % this gets filled in later...
header{9} = {'SignalsExported=LED_cmd'};

tmp = repmat({'LED_cmd'}, 1, nSweeps);
header{10,:} = {'Signals=', tmp{:}};


tmp = cellfun(@(x,y) sprintf(x,y),...
       repmat({'Trace #%d (V)'}, 1, nSweeps),...
       mat2cell([1:nSweeps]', ones(nSweeps,1))', 'uniformoutput', false);
header{11,:} = {'Time (s)', tmp{:}};



% open a new file
fileID = fopen(params.name,'w');


% iterate over the header, adding line by line. All the entries are
% strings, but I need to append \t and \n characters appropriately
specMotif = '%s \t ';
for row = 1:size(header,1);
    nCols = size(header{row},2);
    formatSpec = repmat(specMotif, 1, nCols-1);
    formatSpec = [formatSpec, '%s \n'];
    fprintf(fileID, formatSpec, header{row}{1:end});
end



% iterate over the sweeps, adding line by line. All the entries are
% doubles, but I need to append \t and \n characters appropriately
specMotif = '%.12f \t ';
nCols = size(sweeps,2);
formatSpec = repmat(specMotif, 1, nCols-1);
formatSpec = [formatSpec, '%.12f \n'];
for row = 1:size(sweeps,1);
    fprintf(fileID, formatSpec, sweeps(row,:));
end

fclose(fileID);



