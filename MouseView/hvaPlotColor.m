function [clr_raw, clr_avg] = hvaPlotColor(hva)


switch lower(hva)
    case {'pm', 'py'}
        clr_raw = 'c';
        clr_avg = 'b';
    case {'al', 'som'}
        clr_raw = 'm';
        clr_avg = 'r';
    case {'lm', 'in'}
        clr_raw = [.6 1 .6];
        clr_avg = 'g';
    otherwise
        clr_raw = [.7 .7 .7];
        clr_avg = [.85 .85 .85];
end