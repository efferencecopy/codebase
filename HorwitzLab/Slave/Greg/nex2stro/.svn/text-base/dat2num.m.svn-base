function nums = dat2num(vals, codes, offset, dataType, bookend)

%   DAT2NUM
%   
%   EXAMPLE: nums = dat2num(vals, codes, offset, dataType, [bookend]);
%
% Takes the values from a plexon data file and extracts the appropriate
% number of values following each instance of an 8000 level code and
% converts them to a numeric representation for MatLab (i.e. doubles).
% 'vals' should be the vector of event values and ecodes. 'codes' should be
% the vector of ecodes to find in 'vals'. 'offset' should be the offset value
% added during transmission from Rex to the Plexon data stream. 'dataType'
% should be a string defining how Rex represented the data represented by
% 'codes' before they were droped into the Plexon data file. 'dataType' can
% be int, long, or double (case insensitive). 'bookend' should be set to 1
% if the values to be converted are bookended with ecodes. If 'bookend' is
% set to 1 and there are not exactly 2 occurances of the ecode than this
% function skips to the next ecode or returns an empty cell array if called
% with only one ecode. 'nums' is a cell array with the same number of
% elements as 'code' (i.e. one vector for each ecode in 'code')
%

% cah 1/15/08

%******** CONSTANTS **********
lowCutoff = 2000;  %vals below 2000 shouldn't occurr via the low priority code dropping
                   % b/c they'll get confused with real 1000 level ecodes. 
%*****************************

if strcmpi(dataType, 'double')
    uintsPerNumber = 8;
elseif strcmpi(dataType, 'long')
    uintsPerNumber = 4;
elseif strcmpi(dataType, 'float')
    uintsPerNumber = 4;
elseif strcmpi(dataType, 'int')
    uintsPerNumber = 1;
else
    error('Unspecified or inadequate data type specifier.');
end

%find out of range codes (in this case codes that could conflict with 1000
%level codes) and remove them.
vals(vals < lowCutoff) = [];
for j = 1:length(codes)
    nums{j} = [NaN]; %initialize an empty vector to return if there are errors
    codeInd = find(vals == codes(j));
    numCodes = length(codeInd);

    if (numCodes < 1)
        %fprintf('\n *** No <%d> ecodes were found ***\n', codes(j));
        continue %ignore one bad code and continue to the rest
    end

    %Now extract the appropriate number of values following (or between)
    %occurances of the ecode and subtract the offset value.
    if (exist('bookend', 'var') && bookend)
        if (numCodes ~= 2)
            error('\n *** Bookending expected 2 <%d> ecodes but found %d ***\n', codes(j), numCodes);
        end

        %fill up the uint8 vectors and put them into the cell array
        uint8s = vals((codeInd(1)+1) : (codeInd(2)-1))' - offset; %transpose to make column vect
    else
        %make a matrix of indicies. then reshape and extract these values from
        %vals
        ind = nan(numCodes, uintsPerNumber);
        for a = 1:uintsPerNumber
            ind(:,a) = codeInd+a;
        end
        ind = reshape(ind', (numCodes.*uintsPerNumber), 1);
        uint8s = vals(ind)' - offset; %transpose to make column vect
    end

    %now convert to a numeric representation that makes sense for MatLab
    %(i.e., double precision floats)
    nums{j} = uint2num(uint8s, dataType);
end



