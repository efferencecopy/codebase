function plothelper(hand, type)

if ~exist('hand', 'var') || isempty(hand)
    hand = gcf;
end

if ~exist('type', 'var') || isempty(type)
    type = 'big';
end


switch lower(type)
    case 'big'
        set(hand, 'DefaultAxesFontName','Arial')
        set(hand, 'DefaultAxesFontSize', 26)
        set(hand, 'DefaultLineLinewidth', 2)
end