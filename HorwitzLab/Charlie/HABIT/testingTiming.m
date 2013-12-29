%% TIMING ANALYSIS FOR HABITUATION PARADIGM
% this code should determine if the expected number of frame refreshes are
% consistent with the number dropped into the rex data file
stro = nex2stro('N:\NexFiles\Charlie\habit_test1.nex'); %oops I deleted this file. Have to make a new one.
habitIndicies

numTrials = size(stro.trial, 1);
for a = 1:numTrials
    startTime = stro.ras{a,3};
    dT = 1/stro.sum.analog.storeRates{1};
    nPts = length(stro.ras{a,2});
    x_tmp = startTime + [0:dT:nPts*dT];
    
    %first the habit stim
    h_ind = (x_tmp > stro.trial(a,repHabitOnInd)) & (x_tmp <= stro.trial(a, repHabitOffInd));
    h_time = x_tmp(h_ind);
    h_sig = stro.ras{a,1}(h_ind);
    h_peaks(a) = sum(diff(sign(h_sig-mean(h_sig))) == 2);
    h_troughs(a) = sum(diff(sign(h_sig-mean(h_sig))) == -2);
    h_conv(a) = sum(conv([-1, -1, -1, -1, 1, 1, 1, 1], sign(h_sig-mean(h_sig))) == 8);    
%     figure(f1)
%     plot(h_time-h_time(1), h_sig, 'b')
    
    
    %now for the gabor
    g_ind = (x_tmp > stro.trial(a,repGaborOnInd)) & (x_tmp <= stro.trial(a, repGaborOffInd));
    g_time = x_tmp(g_ind);
    g_sig = stro.ras{a,1}(g_ind);
    g_peaks(a) = sum(diff(sign(g_sig-mean(g_sig))) == 2);
    g_troughs(a) = sum(diff(sign(g_sig-mean(g_sig))) == -2);
    g_conv(a) = sum(conv([-1, -1, -1, -1, 1, 1, 1, 1], sign(g_sig-mean(g_sig))) == 8);
%     figure(f2)
%     plot(g_time-g_time(1), g_sig, 'k')
    
end
h_dat = [h_peaks(:), h_troughs(:), h_conv(:), stro.trial(:, habitFramesInd)]
g_dat = [g_peaks(:), g_troughs(:), g_conv(:), stro.trial(:, gaborFramesInd)]
    

%% drift rate check
%using the signal from the photodiode to check the drift rate of the gabor.

stro = nex2stro('C:\PlexonData\time test 3.nex');
sampleRate = stro.sum.analog.storeRates{1};
habitIndicies

%find the t2 trials
l_t2Trials = stro.trial(:, gaborXInd) > 0;

t2Trials = find(l_t2Trials);
for a = 1:length(t2Trials)
    gabOnTime = stro.trial(t2Trials(a), repGaborOnInd);
    gabOffTime = stro.trial(t2Trials(a), repGaborOffInd);
    timevec = [0:length(stro.ras{t2Trials(a),1})-1]./sampleRate + stro.ras{t2Trials(a), 3};
    validSamps = (timevec>gabOnTime) & (timevec<gabOffTime);
    tmp = stro.ras{t2Trials(a),1}(validSamps);
    signal2(a,:) = tmp(1:round(sampleRate*9.8));
end

t1Trials = find(~l_t2Trials);
for a = 1:length(t1Trials)
    gabOnTime = stro.trial(t1Trials(a), repGaborOnInd);
    gabOffTime = stro.trial(t1Trials(a), repGaborOffInd);
    timevec = [0:length(stro.ras{t1Trials(a),1})-1]./sampleRate + stro.ras{t1Trials(a), 3};
    validSamps = (timevec>gabOnTime) & (timevec<gabOffTime);
    tmp = stro.ras{t1Trials(a),1}(validSamps);
    signal1(a,:) = tmp(1:round(sampleRate*9.8));
end

diffSig = mean(signal2) - mean(signal1);
[pow, freq] = myfft(diffSig, sampleRate);
[m, ind] = max(pow);
domFreq = freq(ind)

%RESULTS
% when I use a "zero" sf gabor drifting at 1cyc/sec I see a domFreq in the
% power spect at 1Hz. When I use a 1cpd sf gabor drifting at 3cyc/sec the
% power spect is more noisy, but there is a peak at 3hz. If I do this again
% I might consider smoothing the signals before fft to reduce noise in the
% pow spect.
