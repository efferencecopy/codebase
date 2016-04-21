function spd = myPR655parsespdstr(readStr, S_out)
% spd = PR655parsespdstr(readStr,S)
%
% Parse the spectral power distribution string
% returned by the PR655.
% 
% 01/16/09    tbc   Adapted from PR650Toolbox for use with PR655
% 10/2015     Charlie Hass: bug fixes, the original code assumes a
%               consistent length of output from the PR655
%

if nargin < 2 || isempty(S_out)
	S_out = [380 4, 101];
end

idx_returns = regexpi(readStr, '\n');
Nlines = numel(idx_returns);
[S_in, spd] = deal(nan(Nlines-1,1));
for i_line = 1:Nlines-1
    
    % each line is a wavelength reading
    tmpstring = readStr(idx_returns(i_line)+1: idx_returns(i_line+1)-1);  % plus/minus 1 to avoid putting returns in the tmpstring
    
    % stuff before the comma is the wavelength
    idx_comma = regexpi(tmpstring, ',');
    S_in(i_line) = str2double(tmpstring(1:idx_comma-1));
    
    % stuff after the comma is the spd power reading
    spd(i_line) = str2double(tmpstring(idx_comma+1:end));
end

% Convert to our units standard
spd = 4 * spd;

% convert S_in into [start, delta, N] format
N = numel(S_in);
delta = unique(diff(S_in));
assert(numel(delta)==1, 'ERROR: sampling latice is not consistent');
S_in = [S_in(1), delta, N];

% Spline to desired wavelength sampling.
all(S_in == S_out)
if ~all(S_in == S_out)
 spd = SplineSpd(S_in, spd, S_out);
end
