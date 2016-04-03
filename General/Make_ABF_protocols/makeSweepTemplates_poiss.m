function params = makeSweepTemplates_poiss(params, UNITTEST)

% params should have
%
% params.si             =>  the sample INTERVAL (needs to be an iteger)
% params.swpDur         =>  The total duration of the sweep IN NUMBERS OF SAMPLES!!!!
% params.tStart         =>  the time of the first pulse
% params.pAmp           =>  A vector of amplitudes for the pulse height [interleaved variable]
% params.pWidth         =>  A vector of pulse widths (in seconds)  [interleaved variable]
% params.ritFreq        =>  The frequency of the poiss train (approx)
% params.ritHiFreqCut   =>  Cut out the frequencies above this value (in Hz)


if exist('UNITTEST', 'var') && strcmpi(UNITTEST, 'unittest')
    params.pAmp = linspace(1,10, 100); % forces there to be 100 trains
    params.ritFreq = 20;  % forces all the trains to have the same freq
    params.ritFreqCut = 58;
else
    UNITTEST = false;
end


%
% Generate the stimulus waveforms (one for each unique stimulus type)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%
nAmps = numel(params.pAmp);
nFreqs = numel(params.ritFreq);

tStartIdx = ceil(params.tStart ./ params.si);
conditions = fullfact([nAmps, nFreqs]);

params.templates_poiss = repmat({zeros(params.swpDur, 1)}, 1, size(conditions, 1));

% loop over the conditions and construct the waveform for each sweep
for i_cond = 1:size(conditions, 1)
    
    % define some constants for this sweep template
    tmp_pAmp = params.pAmp(conditions(i_cond, 1));
    tmp_pFreq = params.ritFreq(conditions(i_cond, 2));
    samplesPerPulse = ceil(params.pWidth ./ params.si);
    
    % determine the total poiss train duration in samples. assume that it
    % will start at the tStartIdx, and end 500ms before the end of the sweep
    nSampsPerTrain = params.swpDur - tStartIdx - round(0.500 ./ params.si);
    tt = [0:nSampsPerTrain-1] .* params.si;
    
    % now figure out the times of the pulses according to a quasi-poisson
    % process (almost a renewal-process?)
    randNums = unifrnd(0,1, 1, nSampsPerTrain);
    spkThresh = tmp_pFreq .* params.si;
    spikeTrain = randNums < spkThresh;
    spikeTimes = tt(spikeTrain);
    
    minISI = 1 ./ params.ritHiFreqCut;
    tooSoon = [false, diff(spikeTimes)<minISI];
    tooSoonIdx = find(tooSoon);
    spikeIdx = find(spikeTrain);
    spikeIdx(tooSoonIdx) = []; %#ok<FNDSB>
    
    % delay the first (and all subsequent) pulses till the desired time
    spikeIdx = spikeIdx + tStartIdx;
    
    % now make the pulses at the times specified in spikeIdx
    spikeIdx = repmat(spikeIdx(:), 1, samplesPerPulse);
    spikeIdx = bsxfun(@plus, spikeIdx, [0:samplesPerPulse-1]);
    spikeIdx = spikeIdx'; % now each column is a pulse
    params.templates_poiss{i_cond}(spikeIdx) = tmp_pAmp;
    
end

if UNITTEST
    
    for i_cond = 1:numel(params.templates_poiss)
        
        % figure out the ISIs
        aboveThresh = params.templates_poiss{i_cond} > 0.5;
        crossings = diff(aboveThresh)==1;
        isi_samps = diff(find(crossings));
        isi_sec = isi_samps .* params.si;
        all_isi(i_cond,:) = histcounts(isi_sec, [0:0.001:1]);
        
    end
    
    % plot the resulting (average) CDF
    figure
    avgcounts = mean(all_isi, 1);
    stairs([0:0.001:0.999], cumsum(avgcounts./sum(avgcounts)));
    set(gca, 'xscale', 'log')
    axis tight
    
end

