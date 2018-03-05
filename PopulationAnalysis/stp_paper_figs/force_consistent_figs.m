function hand = force_consistent_figs(hand, type)

switch type
    case 'ax'
        hand.Box = 'off';
        hand.TickDir = 'out';
        hand.FontSize = 16;
    case 'fig'
end