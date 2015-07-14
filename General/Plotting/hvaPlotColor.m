function [clr_raw, clr_avg] = hvaPlotColor(hva)


switch lower(hva)
    case {'pm', 'py'}
        clr_avg = 'b';
        clr_raw = [0.6 0.6 1];
        
        
    case {'al', 'som'}
        clr_avg = 'r';
        clr_raw = [1 0.6 0.6];
        
        
    case {'lm', 'in'}
        clr_raw = [.6 1 .6];
        clr_avg = 'g';
        
    case 'erc'
        clr_avg = [ 0.5804, 0, 0.8275];
        clr_raw = [0.8549, 0.4392, 0.8392];
        
        
    case 'rl'
        clr_avg = [0.8235, 0.4118, 0.1176];
        clr_raw = [0.8706, 0.7216, 0.5294];
        
        
    otherwise
        clr_raw = [.5 .5 .5];
        clr_avg = [.3 .3 .3];
end