function out = butterfilt(in, freqs, sampFreq, type, dim)
%
% out = butterfilt(in, freqs, sampFreq, type)
%
%applies filter using a butterworth filter kernel. Type should be 'high',
%'low', or 'stop'. If type = 'stop', freqs should contain two values.


Wn = freqs./(sampFreq/2);
order = 4;
[B,A] = butter(order, Wn, type);
if exist('dim', 'var')
    out = filter(B,A,in,[],dim);
else
    out = filter(B,A,in);
end