function invgamma = InvertGamma(cal, colourmode)
% function invgamma = InvertGamma(calstruct, colourmode)
%
% Inverts the measured/fitted gamma functions.  These inverted gamma functions
% can be passed to LoadNormalizedGammaTable to linearize the monitor.  This
% does not work with the Bits++ because the digital output of the graphics
% card is < 14 bits (probably 10 or even 8).  If we are using the Bits++
% we'll still want to know these inverted gamma functions to know what
% voltages to pass out to get the desired intensities.
%
% This is kind of ugly, because I hate splining a fit to the data, but
% arguably the fit gamma functions are better indicators of the truth than
% are my noisy measurements.  I've tried it both ways, and the results
% appear very similar.

if colourmode
    nvoltlevels = 2^16;
else
    nvoltlevels = 2^8;
end

invgamma = zeros(nvoltlevels,3);
for gun = 1:3
    if ~isstruct(cal)
        normInt = cal(:,gun);
    else
        %normInt = cal.rawdata.rawGammaTable(:,gun); % Normalized intensity (raw data)
        normInt = cal.gammaTable(:,gun);  % Normalized intensity (using fit)
    end
    normVolt = linspace(0,1,size(normInt,1));  % Normalized voltage
    idx = find(normInt == 0, 1, 'last');  % Getting rid of measurements that were too dim
    if ~isempty(idx)
        normInt = normInt(idx:end);
        normVolt = normVolt(idx:end);
    end
%     interpVolt = spline(normInt, normVolt, linspace(0,1,nvoltlevels));
    interpVolt = interp1(normInt, normVolt, linspace(0,1,nvoltlevels), 'spline'); % 2x speed up using this -- ZLB
    invgamma(:,gun) = interpVolt';
end
