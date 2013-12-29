function out = myhex2num(in)

%convert into integers (represented in double precision fp by MatLab).



%determine the number of hex characters per number. Take the character array "in"
%and find the indices to the white spaces, then take the element wise
%difference (diff). The max of this number (minus 1) is the number of hex
%characters per number. Remeber that there are two leading place holders
%(one leading space and one place holder for the sign) so add two to the
%number of hexCharsPerNum to get the number of printing characters per
%number. If some of the numbers are negative than only add 1.
spacesIndex = find(isspace(in) == 1);
hexCharsPerNum = max(diff(spacesIndex))-1;

%determine if the character array represents double or ints, then determine
%the total number of characters per number (including white space and
%signs)
if (hexCharsPerNum == 16)
    doubleData = 1;
    totalCharsPerNum = 17;
elseif (any(find(in == '-'))); %deal with cases where negs are present
    totalCharsPerNum = hexCharsPerNum + 1;
    doubleData = 0;
else
    totalCharsPerNum = hexCharsPerNum + 2;
    doubleData = 0;
end

numElements = length(in) ./ totalCharsPerNum;
bigColumn = reshape(in, totalCharsPerNum, numElements)';

% convert back to numeric representations differently for doubles and ints
if(doubleData)
    bigColumn(:,1) = [];
    out = hex2num(bigColumn);
else
    negatives = bigColumn(:,2) == '-';
    bigColumn(:, 1:2) = []; %take off the leading spaces.
    out = hex2dec(bigColumn);
    out(negatives) = out(negatives) .* -1; %put the negatives back
end