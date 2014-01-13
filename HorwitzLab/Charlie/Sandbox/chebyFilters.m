%% trying to demonstrate a knowledge of filter design. Attempting to
% implement a 2nd order, zero-phase, cheby1 filter to produce bandlimited
% LFP data

fin

%generate some white noise data to filter later;
sampFreq = 500; 
nSeconds = 10;
ntrials = 1000;
whiteNoise = normrnd(0, 6, ntrials, sampFreq*nSeconds);



dt = (1/sampFreq);
tt = 0:dt:nSeconds-dt;
tt = repmat(tt, ntrials, 1);
freqs = unifrnd(10, 10, ntrials, 1);
whiteNoise = bsxfun(@times, tt, freqs(:));
whiteNoise = sin(2.*pi.*whiteNoise);



%run the notch filter
% flow = 59.3;
% fhigh = 60.7;
% order = 8;
% r =.1;
% Wn = [flow*2/sampFreq,fhigh*2/sampFreq];
% [B, A] = cheby1(order/2,r, Wn, 'stop');
% Bcell = mat2cell(repmat(B,1000,1), ones(1000,1));
% Acell = mat2cell(repmat(A,1000,1), ones(1000,1));
% tmp = cellfun(@filtfilt, Bcell, Acell, mat2cell(whiteNoise, ones(1000,1)), 'uniformoutput', 0);
% whiteNoise = vertcat(tmp{:});

%now generating a cheby filter
flow = 1;
fhigh = 10;
order = 4;
r = 0.1;
nyquist = sampFreq/2;
Wn = [(flow/nyquist) , fhigh/nyquist]
[B, A] = cheby1(order/2,r, Wn);
%[B, A] = butter(order/2, Wn);

%generating greg's thing
centerfreq = 5; % Hz
cyclespersample = centerfreq./sampFreq;
ncyclesper6sigma = 6;  % Controls the bandwidth
nsamplesper6sigma = ceil(ncyclesper6sigma/cyclespersample);
ncycles = nsamplesper6sigma*cyclespersample;
filtkernel1 = normpdf(linspace(-3,3,nsamplesper6sigma),0,1).*cos(linspace(0,2*pi*ncycles,nsamplesper6sigma));
filtkernel1 = filtkernel1./norm(filtkernel1);

powRaw = [];
powCheby=[];
powGreg=[];
for i = 1:size(whiteNoise,1)
    tmp = conv(whiteNoise(i,:), filtkernel1,'same')./numel(whiteNoise(i,:));
    [powGreg(i,:)] = fft(tmp);
    
    windowed = whiteNoise(i,:) .* hanning(size(whiteNoise,2))';
    [powCheby(i,:)] = fft(filtfilt(B,A,windowed));
    
    powRaw = fft(whiteNoise(i,:));
end


figure
nyquist = sampFreq/2;
freq = linspace(-nyquist, nyquist, sampFreq .* nSeconds);
subplot(3,1,1)
plot(freq, fftshift(mean(abs(powRaw), 1)));
subplot(3,1,2)
plot(freq, fftshift(mean(abs(powCheby),1)));
subplot(3,1,3)
plot(freq, fftshift(mean(abs(powGreg),1)));

figure
posFreq = freq>=0;
subplot(3,1,1)
x = (mean(abs(powRaw).^2, 1));
plot(ifft(x));
xlim([0 80])
subplot(3,1,2)
x = (mean(abs(powCheby).^2, 1));
plot(ifft(x));
xlim([0 80])
subplot(3,1,3)
x = (mean(abs(powGreg).^2, 1));
plot(ifft(x));
xlim([0 80])





