function [clr_raw, clr_avg] = hvaPlotColor(hva)


switch lower(hva)
    case 'pm'
        clr_raw = 'c';
        clr_avg = 'b';
    case 'al'
        clr_raw = 'm';
        clr_avg = 'r';
    otherwise
        clr_raw = 'k';
        clr_avg = 'k';
end