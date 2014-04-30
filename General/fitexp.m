function fit_fxnhand = fitexp(in)

%
% FUNCTION fitexp
%
%  example: fit_fxnHand = fitexp(rawData)
%
% Should accept rawdata inputs. Fits the sum of three exponentials of the
% form:
%
% fit = [A_1 * exp(-tt/tau_1)] + [A_2 * exp(-tt/tau_2)] + [A_3 * exp(-tt/tau_3)] + C
% 
% Returns a function handle that accepts a tt (time) vector, which will
% create the smooth-fitted curve. Times are with refrence to the beginning
% of the raw data set (such that t=0 is the first index to "rawData"
%
%
% C.Hass


% penalize the fit if the asymptotic value is less than the baseline value.
