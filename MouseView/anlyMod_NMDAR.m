function params = anlyMod_NMDAR(params)


% The goal is to pull out the peak evoked current, but to ignore any
% transients due to the stimulus. The raw current data comes in synched to
% LED pulse onset (time=0). Ignore the first 1.5 ms worth of data, then
% find the peak. Since I don't know if it's inward or outward, I have to
% trust that transients won't rear their heads...

tvec = params.ivdat.tvec;
l_anlyWindow = tvec > 0.0015;

% iterate over the channels, pull out the peak current. Note the time of
% the peak, so that I can fit an exponential if I want...
