function [clr_raw, clr_avg] = hvaPlotColor(hva)


switch lower(hva)
    case 'pm'
        clr_raw = 'c';
        clr_avg = 'b';
    case 'al'
        clr_raw = 'm';
        clr_avg = 'r';
    case 'lm'
        clr_raw = 'k';
        clr_avg = 'k';
    otherwise
        clr_raw = [.8 .8 .8];
        clr_avg = [.8 .8 .8];
end