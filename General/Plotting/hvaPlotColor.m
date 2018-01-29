function plt_clr = hvaPlotColor(hva)

switch lower(hva)
    case 'pm'
        plt_clr = [0, 0.4470, 0.7410];
        
    case 'al'
        plt_clr = [0.9290    0.6940    0.1250];
        
    case 'lm'
        plt_clr = [0.8500    0.3250    0.0980];
        
    case 'am'
        plt_clr = [0.4940    0.1840    0.5560];
             
    case 'erc'
        plt_clr = [0 0.5 0];
        
    otherwise
        error('HVA name not recognized')
end