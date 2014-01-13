% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pull out the relevant strobed codes (assume 8 uint8's per double)
% and convert to floating point numbers.
function y = codes2num(x)
    if (rem(length(x),8))
        error('Number of elements received in bulk transfer not divisible by 8');
    end
    numElements = length(x)/8;
    uint8s = reshape(x, 8, numElements);
    uint8s = uint8s'; %reshape so that an individual double is represented by each row
    uint8s = fliplr(uint8s);
    y = zeros(numElements, 1);
    for a = 1:numElements
        hex = dec2hex(uint8s(a,:), 2);
        y(a) = hex2num(reshape(hex', 1, 16));
    end
 end