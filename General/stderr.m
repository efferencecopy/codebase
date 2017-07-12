function sem = stderr(data, dim)
%
%  STANDARD ERROR OF THE MEAN 
%
% EXAMPLE:    sem = stderr(data, [dim])
%
% C.Hass 2015

if isempty(data)
    sem = NaN;
    return
end


if ~exist('dim', 'var')
    if isvector(data)
        if isrow(data)
            dim = 2;
        else
            dim = 1;
        end
    else
        dim = 1;
    end
end

sigma = nanstd(data, [], dim);
N = sum(~isnan(data), dim);

sem = sigma ./ sqrt(N);

% fix instances of N=1. Make the SEM for these points == NaN
sem(N==1) = NaN;
