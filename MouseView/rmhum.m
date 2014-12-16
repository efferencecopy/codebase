function out = rmhum(in, sampFreq, winStart, winEnd, lines, showplot)


% check the input arguments
if ~exist('lines', 'var')
    lines = 60;
end

if ~exist('showplot', 'var')
    showplot = false;
end


if ischar(in) && strcmpi(in, 'unittest');
    UNITTEST = 1;
    lines = [60];
    sampFreq = 40e3;
    in = makefakedata(lines, sampFreq, sampFreq .* 2);
    winStart = 1;
    winEnd = numel(in);
else
    UNITTEST = 0;
end


% only anlyzes a subset of the data that is deemed 'quiet'
tmp_in = in(winStart:winEnd);
tmp_in = tmp_in - mean(tmp_in);


% make the frequency axis (this is mostly for making good initial guesses
% for the fminsearch
N = size(tmp_in,1);
if rem(N,2)
    k = -((N-1)./2):((N-1)./2); % pos freqs when N is Odd. For two sided freqs: k = -((N-1)./2):((N-1)./2)
else
    k = -(N./2):((N./2)-1); % 0:(N./2) => abs(neg Freqs) when N is even. For two sided freqs: k = -(N./2):((N./2)-1)
end
ff = (k./N).*sampFreq;



out = tmp_in;
for a = 1:numel(lines)
    
    % look around to see if that's the correct line freq.
    coeffs = fftshift(fft(tmp_in.*hamming(numel(tmp_in))));
    suspect = lines(a);
    valid_ff = ff'>=(suspect-5) & ff'<=(suspect+5);
    peakVals = findpeaks(abs(coeffs(valid_ff)));
    maxVal = max(peakVals);
    idx = (abs(coeffs) == maxVal);
    suspect = ff(idx&valid_ff);
    
    
    
    % find the bins closest idx to the +/-linenoise
    [~, pos_idx] = min(abs(ff - suspect));
    [~, neg_idx] = min(abs(ff + suspect));
    idx = [pos_idx, neg_idx];
    
    
    % pred_mags is easy, just add up the mags for +/- freqs
    pred_mags = sum(abs(coeffs(idx))./N);
    
    % now try to determine the phase of a SINUSOID given by the systhesis
    % direction (ifft), setting all other freq components to zero
    tmp_coeffs = zeros(size(coeffs));
    tmp_coeffs(idx) = coeffs(idx);
    ifft_tmp_in = ifft(ifftshift(tmp_coeffs));
    [~, peakIdx] = findpeaks(ifft_tmp_in);
    tt = [0 : numel(tmp_in)-1]' ./ sampFreq;
    
    timeOfFirstPeak = tt(peakIdx(1));
    period = 1./ff(pos_idx);
    predTimeOfFirstPeak = period/4;
    diffTime = timeOfFirstPeak - predTimeOfFirstPeak;
    
    if diffTime<0
        pred_phase = abs(diffTime./period) .* (2*pi);
    else
        diffTime = predTimeOfFirstPeak - timeOfFirstPeak;
        pred_phase = mod((diffTime./period) .* (2*pi), 2*pi);
    end
    
    % generate some initial guesses (based on the FFT stuff above)
    phi = pred_phase;
    amp = pred_mags;
    omega = suspect;
    initGuess = [amp, phi, omega];
    
    % things that are uesful for fminsearch (are are in the scope of the
    % nested sub-function)
    LB = [0, -inf, omega-1];
    UB =[inf,  inf, omega+1];
    
    % meat and potatoes
    errtype = 'dotcorr'; % does a good job with getting the phase and freq
    initGuess = fminsearch(@minNegCorr, initGuess);
    errtype = 'SSE'; % refines the fit by getting the amplitude
    params = fminsearch(@minNegCorr, initGuess);
    
    % subtract out the line noise, and advance to the next line...
    %  WHEN SUBTRACTING FROM THE FULL WAVEFORM I NEED TO REDEFINE TIME ZERO
    %  TO BE THAT OF THE TEMPLATE...
    tt = [0 : numel(tmp_in)-1]' ./ sampFreq;
    fit = sin(2 .* pi .* (tt+params(2)) .* params(3)) .* params(1);
    out = out - fit;
    tmp_in = tmp_in-fit;



% do some basic plotting if the unittest is turned on
if UNITTEST || showplot
    figure,
    set(gcf, 'position', [745   123   183   636])
    subplot(3,1,1), hold on,
    tt = [0 : numel(tmp_in)-1]' ./ sampFreq;
    plot(tt, tmp_in, 'b')
    plot(tt, fit, 'k', 'linewidth', 2)
    %plot(tt, ifft_tmp_in, 'g', 'linewidth', 2)
    %plot(tt(peakIdx(1)), ifft_tmp_in(peakIdx(1)), 'co', 'markerfacecolor', 'c')
    %xlim([0 0.017]);
    
    subplot(3,1,2), hold on,
    plot(tt, tmp_in, 'b')
    plot(tt, out, 'r')
    %xlim([0 0.017]);
    
    subplot(3,1,3), hold on
    plot(ff, fftshift(abs(fft(out)))./numel(fit), 'r');
    plot(ff, fftshift(abs(fft(tmp_in)))./numel(fit), 'b');
    xlim([omega-5 omega+5])
end

end



% nested numerical solver
    function err = minNegCorr(guess)
        
        A = guess(1);
        P = guess(2);
        F = guess(3);
        
        sinwave = sin(2 .* pi .* (tt+P) .* F) .* A;

        switch errtype
            case 'dotcorr'
                rho = (sinwave' * tmp_in(:)) ./ (norm(sinwave).*norm(tmp_in));
                err = abs(1-rho);
            case 'SSE'
                magicScalar = amp*2+1000;
                err = sum(((tmp_in+magicScalar)-(sinwave+magicScalar)).^2);
                
        end
        
        
        
        
        if any([A,P,F] < LB)
            err = inf;
        end
        if any([A,P,F] > UB)
            err = inf;
        end
        
    end

end % function







% TESTING ROUTINES

function [in, amps] = makefakedata(lines, sampFreq, N)

% generate some white noise
in = normrnd(0,2,N,1);

% randomize the line freq +/- 4
%lines = lines + unidrnd(8,1)-4
lines = lines + unifrnd(-4,4,1)

% generate a sin wave and add to the WN
tt = linspace(0, N./sampFreq, N);
tmp = bsxfun(@times, tt', lines); % inner part (freq * time)
tmp = tmp .* 2*pi; % inner part, now = (2*pi*freq*time)
phase = unifrnd(0, 2*pi, size(lines))
tmp = bsxfun(@plus, tmp, phase);
linenoise = sin(tmp);


% scale the line noise by a random amplitude
amps = unifrnd(1,5,size(lines))
linenoise = bsxfun(@times, linenoise, amps);

% add the line noise to the WN.
in = sum([in, linenoise], 2);


end

