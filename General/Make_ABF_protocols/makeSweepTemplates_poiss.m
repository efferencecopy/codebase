function params = makeSweepTemplates_poiss(params, UNITTEST)

% params should have
%
% params.si              =>  the sample INTERVAL (needs to be an iteger)
% params.swpDur          =>  The total duration of the sweep IN NUMBERS OF SAMPLES!!!!
% params.tStart          =>  the time of the first pulse
% params.pAmp            =>  A vector of amplitudes for the pulse height [interleaved variable]
% params.pWidth          =>  A vector of pulse widths (in seconds)  [interleaved variable]
% params.ritFreq         =>  The frequency of the poiss train (approx)
% params.ritHiFreqCut    =>  Cut out the frequencies above this value (in Hz)
% params.rit_Nversions   =>  Number of different verions of RITs with identical params
% params.ritUseEnvelope  =>  True or false
% params.ritEnvelopeFreq =>  The frequency of a firing rate envelope. Can also be a vector of frequencies


if exist('UNITTEST', 'var') && UNITTEST
    params.rit_Nversions = 100; % forces there to be 100 trains
else
    UNITTEST = false;
end

if ~isfield(params, 'ritUseEnvelope')
    params.ritUseEnvelope = false;
end
if ~isfield(params, 'ritEnvelopeFreq') || isempty(params.ritEnvelopeFreq)
    params.ritEnvelopeFreq = sqrt(-1); % pretty sure this should make things crash if it ends up getting used
end


%
% Generate the stimulus waveforms (one for each unique stimulus type)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%
nAmps = numel(params.pAmp);
nFreqs = numel(params.ritFreq);
nEnvelopes = numel(params.ritEnvelopeFreq);
nVersions = params.rit_Nversions;

tStartIdx = ceil(params.tStart ./ params.si);
conditions = fullfact([nAmps, nFreqs, nEnvelopes, nVersions]);

% initalize the outputs
params.templates_poiss = repmat({zeros(params.swpDur, 1)}, 1, size(conditions, 1));
params.conditions_poiss = conditions; % needs to be updated on a cond by cond basis
params.ttype_header_poiss = {'pAmp', 'pFreq', 'envFreq', 'versionNum'};

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
    
    % make a firing rate envelope (if need be)
    if params.ritUseEnvelope
        tmp_envFreq = params.ritEnvelopeFreq(conditions(i_cond,3));
        envelope = sin(2.*pi.*tt.*tmp_envFreq);
        envelope(envelope<0) = 0;
    else
        envelope = 1;
        tmp_envFreq = nan;
    end
    
    
    % now figure out the times of the pulses according to a quasi-poisson
    % process (almost a renewal-process?)
    randNums = unifrnd(0,1, 1, nSampsPerTrain);
    spkThresh = (envelope .* tmp_pFreq) .* params.si;
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
    
    
    % define the per-condition ttype array
    params.conditions_poiss(i_cond,1:3) = [tmp_pAmp, tmp_pFreq, tmp_envFreq]; % the version number (4th column) is inherited from fullfact as defined above
    
    
end

if UNITTEST
    
    for i_cond = 1:numel(params.templates_poiss)
        
        % figure out the ISIs
        aboveThresh = params.templates_poiss{i_cond} > 0.5;
        crossings = diff(aboveThresh)==1;
        isi_samps = diff(find(crossings));
        isi_sec = isi_samps .* params.si;
        all_isi(i_cond,:) = histc(isi_sec, [0:0.001:2]);
        
        autocorr(:, i_cond) = xcorr(aboveThresh(:));
        
    end
    
    % plot the resulting (average) CDF
    figure
    avgcounts = mean(all_isi, 1);
    stairs([0:0.001:2], cumsum(avgcounts./sum(avgcounts)));
    set(gca, 'xscale', 'log')
    axis tight
    title('average CDF of ISIs')
    
    figure
    hist(sum(all_isi,2))
    title('number of pulses per sweep')
    
    % plot the average auto corr
    figure
    plot(mean(autocorr, 2))
    
end

