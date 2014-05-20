function out = butterFilt(in, sampFreq, cutoff, type)

nyquist = sampFreq./2;
Wn = cutoff./nyquist;

[B,A] = butter(2, Wn, type);

out = filtfilt(B, A, in);
end



