function makeStimulusFile(params)

% params should have
%
% params.name     =>  the name of the output .atf file
% params.type     =>  'train', 'pulse'
% params.si       =>  the sample INTERVAL (needs to be an iteger)
% params.swpDur   =>  The total duration of the sweep IN NUMBERS OF SAMPLES!!!!
% params.tStart   =>  the time of the first pulse
% params.pAmp     =>  A vector of amplitudes for the pulse height [interleaved variable]
% params.pFreq    =>  A vector of frequencies for the pulse train [interleaved variable]
% params.nPulses  =>  The number of pulses in a pulse train
% params.pWidth   =>  A vector of pulse widths (in seconds)  [interleaved variable]
% params.nReps    =>  Number of repeates each stimulus should be presented


% force default frequency to be zero (only one pulse)
if isempty(params.pFreq)
    params.pFreq = 0;
end


% Generate the stimulus waveforms (one for each unique stimulus type)
nAmps = numel(params.pAmp);
nFreqs = numel(params.pFreq);
nPulseWidths = numel(params.pWidth);
nSweeps = nAmps .* nFreqs .* nPulseWidths .* params.nReps;

tt = [0:params.swpDur-1]' .* params.si;
tStartIdx = ceil(params.tStart ./ params.si);
conditions = fullfact([nAmps, nFreqs, nPulseWidths]);

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
    samplesPerPulse = ceil(tmp_pWidth ./ params.si);
    
    if tmp_pFreq == 0 % only one pulse
        templates{i_cond}(tStartIdx : tStartIdx+samplesPerPulse-1) = tmp_pAmp;
        templates{i_cond}(tStartIdx+samplesPerPulse : end) = 0;
    else
        samplesPerPeriod = ceil(1./tmp_pFreq ./ params.si);
        samplesPerIPI = samplesPerPeriod - samplesPerPulse;
        motif = zeros(1,samplesPerPeriod);
        motif(1:samplesPerPulse-1) = tmp_pAmp;
        idx = tStartIdx;
    end
end




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
fileID = fopen('celldata.atf','w');


% interate over the header, adding line by line. All the entries are
% strings, but I need to append \t and \n characters appropriately
specMotif = '%s \t ';
for row = 1:size(header,1);
    nCols = size(header{row},2);
    formatSpec = repmat(specMotif, 1, nCols-1);
    formatSpec = [formatSpec, '%s \n'];
    fprintf(fileID, formatSpec, header{row}{1:end});
end
fclose(fileID);



