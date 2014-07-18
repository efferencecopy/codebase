function out = butterfilt(in, freqs, sampFreq, type)
%
% out = butterfilt(in, freqs, sampFreq, type)
%
%applies filter using a butterworth filter kernel. Type should be 'high',
%'low', or 'stop'. If type = 'stop', freqs should contain two values.


Wn = freqs./(sampFreq/2);
order = 4;
[B,A] = butter(order, Wn, type);
out = filter(B,A,in);