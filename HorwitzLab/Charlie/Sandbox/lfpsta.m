% For now, I'm just trying to see if there's anything in the LFP data we
% currently save. As a first pass I'm going to try and create spike
% triggered LFPs

%stealing this code from gratingsAnalysis
GT=nex2stro(findfile('K042710001.nex'));
protocols = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'protocol'));
frameRate = GT.sum.exptParams.framerate;
anlgRate = GT.sum.analog.storeRates{1};
anlgStartIdx = strcmp(GT.sum.rasterCells, 'anlgStartTime');
nFrames = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'nframes'));
fpacq_t = GT.trial(:, strcmp(GT.sum.trialFields(1,:), 'fp_acq'));
stimOn_t = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'stim_on'));
stimOff_t= GT.trial(:,strcmp(GT.sum.trialFields(1,:),'stim_off'));
spikeIdx = strcmp(GT.sum.rasterCells(1,:),'sig001a');
lfpIdx = strcmp(GT.sum.rasterCells(1,:), 'AD01');
Lcc = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'lcont'));
Mcc = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'mcont'));
Scc = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'scont'));


%% trying to apply filtering techniques

sampFreq = anlgRate; %in Hz
time = 0:1/1000:5; %5 seconds long
freq1 = 20;  %in Hz
freq2 = 100; %in Hz
wf1 = sin(2*pi*time*freq1);
wf2 = sin(2*pi*time*freq2);
signal = (wf2+wf1)+normrnd(0,1,1,length(time)); % sum of 2 sinusoids w/additive gaussian noise.
%plot to demonstrate the results:
myfft(signal, sampFreq)


%now design a 60 Hz notch filter using some random code I found on the
%interwebs:
deg = 3; %filter deg
Wn = [59.5*2/sampFreq, 60.5*2/sampFreq];
[B, A] = butter(deg, Wn, 'stop');
newSignal = filtfilt(B, A, signal);
myfft(newSignal, sampFreq); %this does a good job of rejecting just the 60Hz line signal.


%now let's try it on a real set of LFP data
lfp = GT.ras{3,lfpIdx};
newlfp = filtfilt(B,A,lfp);
figure,
hold on,
plot(lfp, 'b')
plot(newlfp, 'r', 'linewidth', 2)

%compare the spectra
figure,
subplot(2,1,1)
myfft(lfp, anlgRate);
subplot(2,1,2)
myfft(newlfp, anlgRate);


%% LFP STA all spikes during stimulus presentation

%analysis specific parameters
preTime = 0.300; %in seconds
postTime = 0.300;
preSamples = preTime .* anlgRate;
postSamples = postTime .* anlgRate;
timeForPlot = (-preSamples:postSamples) .* (1./anlgRate);
denoiseMethod = 'mtm'; %'butter', 'mtm'
detrend = 1;


spLFP = nan(10000, preSamples+postSamples+1);
shuffLFP = nan(60000, preSamples+postSamples+1);
for a = 1:size(GT.ras, 1)
    tOnset = stimOn_t(a);
    tOffset = stimOff_t(a);
    tSpikes = GT.ras{a, spikeIdx};
    validSpikes = (tSpikes > tOnset) & (tSpikes < tOffset);
    tSpikes(~validSpikes) = [];
    tAnlgStart = GT.ras{a, anlgStartIdx};
    lfp = GT.ras{a, lfpIdx};
    timeVec = ((0:length(lfp)-1) .* (1./anlgRate)) + tAnlgStart;
    
    % try to get rid of the 60 cycle noise
    switch denoiseMethod
        case 'butter'
            %design a 60 Hz notch filter using some random code I found on the interwebs:
            deg = 3; %filter deg
            f1 = 55; %low edge of the bandstop
            f2 = 65; %upper edge of the bandstop
            Wn = [f1*2/anlgRate, f2*2/anlgRate];
            [B, A] = butter(deg, Wn, 'stop');
            lfp = filtfilt(B, A, lfp);
        case 'mtm'
            params.tapers = [3, 5];
            params.Fs = anlgRate;
            params.fpass = [0, anlgRate./2];
            params.pad = -1;
            lfp = rmlinesc(lfp, params);
    end
    
    % now try to get rid of slow drifts in the signal
    if detrend
        lfp = locdetrend(lfp', anlgRate, [.500 0.100]);
    end
    
    %now compute the spike triggered LFP
    for j = 1:length(tSpikes);
        spikeTime = tSpikes(j);
        [minErr, minIdx] = min(abs(timeVec - spikeTime));
        if (minIdx-preSamples > 0) && (minIdx+postSamples < length(lfp)) %an ugly hack
            tmp = lfp(minIdx-preSamples : minIdx+postSamples);
            nextRow = sum(~isnan(spLFP(:,1)))+1;
            spLFP(nextRow,:) = tmp';
        end
    end
    
    %now do the same analysis, but with shuffled spike times
    shuffMultiplier = -1;
    nShuffles = length(tSpikes) * shuffMultiplier;
    counter = 0;
    while counter < nShuffles;
        shufSpike = unifrnd(tOnset, tOffset);
        [minErr, minIdx] = min(abs(timeVec - shufSpike));
        if (minIdx-preSamples > 0) && (minIdx+postSamples <= length(lfp)) %an ugly hack
            tmp = lfp(minIdx-preSamples : minIdx+postSamples);
            tmp = tmp - mean(tmp); %subtract off any DC due to non-stationarity.
            nextRow = sum(~isnan(shuffLFP(:,1)))+1;
            shuffLFP(nextRow, :) = tmp';
            counter = counter+1;
        end
    end
end
spLFP(isnan(spLFP(:,1)), :) = [];
shuffLFP(isnan(shuffLFP(:,1)), :) = [];

figure
subplot(1,2,1)
hold on,
plot(timeForPlot, mean(spLFP), 'b');
plot(timeForPlot, mean(shuffLFP), 'r:');
xlabel('time (ms)')
ylabel('????')
hold off,
subplot(1,2,2)
hold on,
sigmaSTA = repmat(std(spLFP, [], 2), 1, size(spLFP,2));
sigmaShuff = repmat(std(shuffLFP, [], 2), 1, size(shuffLFP,2));
plot(timeForPlot, mean(spLFP./sigmaSTA), 'b')
plot(timeForPlot, mean(shuffLFP./sigmaShuff), 'r:')
xlabel('time (ms)')
ylabel('norm somethings???')
title(sprintf('Detrend: %d, Denoise: %s', detrend, denoiseMethod));


%seems like there are trials where the spLFP is zero at all timepoints...
%what's going on?!?!? This also means that the sd for those trials is zero,
%and the normalized spLFP is undefined....  'K042710001.nex' exemplifies
%this well. look for nan entries in [spLFP./sigmaSTA].

%% LFP Peri-Event Time Histogram synched to grating onset
%analysis specific parameters
preTime = 0.3; %in seconds
postTime = 0.7;
preSamples = preTime .* anlgRate;
postSamples = postTime .* anlgRate;
timeForPlot = (-preSamples:postSamples) .* (1./anlgRate);
denoiseMethod = 'butter'; %'butter', 'mtm'
detrend = 1;
onsetType = 'stim_on'; % either fp_acq or stim_on

peth = [];
shuffLFP = [];
shuffLFP = nan(10000, preSamples+postSamples+1);
for a = 1:size(GT.ras, 1)
    switch onsetType
        case 'fp_acq'
            onset_t = fpacq_t(a);
        case 'stim_on'
            onset_t = stimOn_t(a);
    end
    lfp = GT.ras{a, lfpIdx};
    tAnlgStart = GT.ras{a, anlgStartIdx};
    timeVec = ((0:length(lfp)-1) .* (1./anlgRate)) + tAnlgStart;
    [val, idxAtStimOn] = min(abs(timeVec - onset_t));
    
    
    % try to get rid of the 60 cycle noise
    switch denoiseMethod
        case 'butter'
            %design a 60 Hz notch filter using some random code I found on the interwebs:
            deg = 3; %filter deg
            f1 = 59.5; %low edge of the bandstop
            f2 = 60.5; %upper edge of the bandstop
            Wn = [f1*2/anlgRate, f2*2/anlgRate];
            [B, A] = butter(deg, Wn, 'stop');
            lfp = filtfilt(B, A, lfp);
        case 'mtm'
            params.tapers = [3, 5];
            params.Fs = anlgRate;
            params.fpass = [0, anlgRate./2];
            params.pad = -1;
            lfp = rmlinesc(lfp, params);
    end
    
    % now try to get rid of slow drifts in the signal
    if detrend
        lfp = locdetrend(lfp', anlgRate, [.500 0.100]);
    end
    
    %compile the PETH
    peth(a,:) = lfp(idxAtStimOn-preSamples : idxAtStimOn+postSamples);
    
    %now do the same analysis, but with shuffled spike times
    shuffPerTrial = 1;
    counter = 0;
    while counter < shuffPerTrial;
        shuffStart = unifrnd(fpacq_t(a), timeVec(end)); %only take samples when the monk is looking at the monitor
        [minErr, minIdx] = min(abs(timeVec - shuffStart));
        if (minIdx-preSamples > 0) && (minIdx+postSamples <= length(lfp)) %an ugly hack
            tmp = lfp(minIdx-preSamples : minIdx+postSamples);
            nextRow = sum(~isnan(shuffLFP(:,1)))+1;
            shuffLFP(nextRow, :) = tmp';
            counter = counter+1;
        end
    end
end
shuffLFP(isnan(shuffLFP(:,end)), :) = [];


figure
hold on,
plot(timeForPlot, mean(peth,1), 'b')
plot(timeForPlot, mean(shuffLFP, 1), 'r:')
xlim([-preTime, postTime]);
xlabel('Time (seconds)')
ylabel('volts (????)')
title(sprintf('detrend: %d, denoise: %s', detrend, num2str(denoiseMethod)));

%% time frequency spectrogram using the multi-taper method
params.pad = 0;
params.Fs = anlgRate;
params.fpass = [0 150];
params.trialave = 1;
T = 0.1;
TW = 5;
W = TW./T;
params.tapers = [TW, 2*TW-1];
movingwin = [T, 0.050];
[S, t, f] = mtspecgramc(peth', movingwin, params);
figure
subplot(1,2,1)
imagesc(flipud(S'))
ytickInds = get(gca, 'ytick');
set(gca, 'yticklabel', flipud(f(ytickInds)'));
xtickInds = get(gca, 'xtick');
set(gca, 'xticklabel', t(xtickInds)-(preTime+t(1)));
title(sprintf('time/frequency spectrogram\n TW: %d, nTaps: %d, W: %.3f', TW, params.tapers(2), W));


%verifying that the avg specttra has a peak at 75... b/c it's weird that
%the t/f spectrogram does not.
pow = [];
f = [];
for a = 1:size(peth,1)
    [pow(a,:), f(a,:)] = pmtm(peth(a,:), TW, [], anlgRate);
end
subplot(1,2,2)
plot(f(1,:), 10.*log10(pow))
ylabel('Power/freq (dB/Hz)')
xlabel('Frequency');
title('MTM Power Spectral Density')

