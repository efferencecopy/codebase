% Returns a vector of absolute numbered pins for the bank you request.
function out_pins = pins_rel2abs(bank, in_pins)
switch upper(bank)
    case 'A'
        bankabs = [60 70 29 39 89 99 48 58 7 17 67 77 26 36 86 96 45 55 4 14 64 74 23 33 83 93 42 52 11 21 71 81];
    case 'B'
        bankabs = [50 10 19 9 79 69 38 28 98 88 57 47 16 6 76 66 35 25 95 85 54 44 13 3 73 63 32 22 92 82 61 51];
    case 'C'
        bankabs = [20 30 80 90 49 59 8 18 68 78 27 37 87 97 46 56 5 15 65 75 24 34 84 94 43 53 2 12 62 72 31 41];
    otherwise
        error('Unknown bank requested.');
end
out_pins = bankabs(in_pins);
