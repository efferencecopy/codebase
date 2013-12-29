%% trying to tell if low pass filters remove stuff due to spiking.
Fs = 20000; % in samps/sec
time = 0:1/Fs:10;

%make a low pass filter
cutoff = 250;
deg = 4;
[B, A] = butter(deg, cutoff/(Fs/2));
for a = 1:10;
    noise = normrnd(0, 3, 1, 20000);
    [n(a,:), freq] = myfft(filtfilt(B,A,noise), Fs);
end
plot(freq, mean(n,1))


%% make some spikes
REFRACTORY = 0;
FILTER = 0;
spFreq = 80; %firing rate for spikes in Hz
freq1 = 70.2; %sinusoidal lfp
pow = [];
for a = 1:100;
    spTimes = cumsum(exprnd(1./spFreq, 1, 20000));
    spTimes(spTimes>max(time)) = [];
    spTimes = round(spTimes*1000);
    if ~spTimes(1)
        spTimes(1) = spTimes(1)+1; %hack to avoid indexing errors later
    end
    if REFRACTORY
        refPeriod = 1;
        refractoryViolations = find(diff(spTimes)<=refPeriod);
        while ~isempty(refractoryViolations);
            spTimes(refractoryViolations+1) = spTimes(refractoryViolations+1) + refPeriod;
            spTimes(spTimes>max(time*1000)) = [];
            refractoryViolations = find(diff(spTimes)<=refPeriod);
        end
    end
    spwf = zeros(1,length(time));
    sampsPerMs = Fs./1000;
    spwf(spTimes*sampsPerMs) = 100;
    lfp = sin(2*pi*time*freq1+unifrnd(0, 2*pi)); %randomizing the phase
    wf = lfp + spwf;
    if FILTER
        wf = filtfilt(B, A, wf);
    end
    [pow(a,:), freq] = myfft(wf, Fs);
end
figure
plot(freq, 10*log10(mean(pow,1)))
title(sprintf('Refractory <%d>, Filter <%d>', REFRACTORY, FILTER))


% To Do: the while loop bonks b/c one of the spike times is occasionally
% at time zero. Also, write code to compare power spectra with and without
% refractory periods.

%% add the 'lfp' to the 'spikes'
wf = lfp + spwf;
myfft(wf, Fs);

%% Low-pass filter the waveform, then look for power at