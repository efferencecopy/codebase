function plothelper(hand, type)

if ~exist('hand', 'var') || isempty(hand)
    hand = gcf;
end

if ~exist('type', 'var') || isempty(type)
    type = 'big';
end


switch lower(type)
    case 'big'
        set(hand, 'DefaultAxesFontSize', 20)
        set(hand, 'DefaultLineLinewidth', 2)
end