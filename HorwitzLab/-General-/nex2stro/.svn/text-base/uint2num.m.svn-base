function num = uint2num(uint8s, dataType)

% EXAMPLE num = uint2num(uint8s, dataType);
%
%Turns uint8s into numeric values according to their original cast
%(specified by the dataType argument). 'uint8s' should be a vector of one
%or more uint8s. 'dataType' should be a string of one of the following data
%types: int, long, double, or float.


% cah 1/15/08
% cah 3/15/08 error checking for doubles less than 15
% gdlh 8/20/08 adding support for floats

num = [];
if (strcmpi(dataType, 'double')|| strcmpi(dataType, 'float'))
    if (strcmpi(dataType, 'double'))
        nhexcharpernum = 8;
    else % float
        nhexcharpernum = 4;
    end
    numElements = length(uint8s) ./ nhexcharpernum;
    if (numElements ~= round(numElements))
        error(['Excess or missing values (i.e. not an integer multiple of ',num2str(nhexcharpernum),' in uint8 array']);
    end
    uint8s = reshape(uint8s, nhexcharpernum, numElements)'; % each row of 8 represents each double
    uint8s = fliplr(uint8s); %puts the high order bytes back on the left
    uint8s = reshape(uint8s', numElements*nhexcharpernum, 1); %reshape again to retain correct order within each pair of hex numbers later
    hex = dec2hex(uint8s);

    % if the number that's getting converted is less than 15 than there
    % will only be a single column of hex characters (8 chars total)
    % which will make the reshape line bonk. DON'T KNOW WHAT THIS ERROR
    % CHECKING WILL DO WITH NUMBERS THAT ARE NOT SUPPOSED TO BE ZERO,
    % I.E. THE BITE ORDER MAY BE MESSED UP!! (cah 3/15/08)
    if (numel(hex) < 2*numElements*nhexcharpernum)
        hex(:,2) = hex;
        hex(:,1) = 48;
    end
    hex = reshape(hex', 2*nhexcharpernum, numElements)';

    if (strcmpi(dataType, 'double'))
        num = hex2num(hex);
    else % float
        num = singlehex2num(hex);
    end
elseif strcmpi(dataType, 'long');
    numElements = length(uint8s) ./ 4;
    if (numElements ~= round(numElements))
        error('Excess or missing values (i.e. not an integer multiple of 4) in uint8 array');
    end
    uint8s = reshape(uint8s, 4, numElements)'; %each row of 4 represents each long
    num = uint8s * [2^0; 2^8; 2^16; 2^24];
elseif strcmpi(dataType, 'int')
    num = uint8s;
else
    error('Unspecified or inadequate data type specifier');
end
