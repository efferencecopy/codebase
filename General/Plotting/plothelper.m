function plothelper(hand, type)

if ~exist('hand', 'var') || isempty(hand)
    hand = gcf;
end

if ~exist('type', 'var') || isempty(type)
    type = 'big';
end


switch lower(type)
    case 'big'
        
        % in case you're about to plot stuff
        set(hand, 'DefaultAxesFontName','Arial')
        set(hand, 'DefaultAxesFontSize', 26)
        set(hand, 'DefaultLineLinewidth', 2)
        
        
        % in case it's already plotted
        xlims = get(gca, 'Xlim');
        ylims = get(gca, 'Ylim');
        fonthands = [gca, cell2mat(get(gca, {'Ylabel', 'Xlabel'}))];
        set(fonthands, 'fontsize', 26)
        set(get(gca, 'title'), 'Fontsize', 28)
        set(gca, 'TickDir', 'out', 'box', 'off')
        set(gca, 'xlim', xlims, 'ylim', ylims)
        
end