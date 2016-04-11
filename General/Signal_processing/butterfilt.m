function out = butterfilt(in, freqs, sampFreq, type, dim)
%
% out = butterfilt(in, freqs, sampFreq, type)
%
%applies filter using a butterworth filter kernel. Type should be 'high',
%'low', or 'stop'. If type = 'stop', freqs should contain two values.


Wn = freqs./(sampFreq/2);
desired_filter_order = 4;
switch lower(type)
    case {'low', 'high'}
        % do nothing, b/c the order argument correctly reflects the
        % resulting filter order
        order = desired_filter_order;
    case {'stop', 'bandpass'}
        order = desired_filter_order ./ 2;
end

[B,A] = butter(order, Wn, type);

if exist('dim', 'var') && (dim==2)
    in = in';
    out = filtfilt(B,A,in);
    out = out';
else
    out = filtfilt(B,A,in);
end