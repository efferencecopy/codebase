function ff = fftfreq(N, samplingRate)

% Calculate the frequncy axis
if rem(N,2)
    k = -((N-1)./2):((N-1)./2); % pos freqs when N is Odd. For two sided freqs: k = -((N-1)./2):((N-1)./2)
else
    k = -(N./2):((N./2)-1); % 0:(N./2) => abs(neg Freqs) when N is even. For two sided freqs: k = -(N./2):((N./2)-1)
end
ff = (k./N).*samplingRate;