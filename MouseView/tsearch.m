function  [match, names] = tsearch(in, searchString, attribute)
% TSEARCH Search a cell array of text for a specific word(s) or combination
% of phrases.
%
%      match = tsearch(in, phrase)
%
% INPUTS:
%   "in"                a cell array where each cell contains text
%   "searchString"      a string that is the target you're looking for
%   "attribute"         a cell array of other atributres (e.g., filenames)
%                       associated with each entry of "in". For example, each entry
%                       of "in" corresponds to a unique .abf file
%
% OUTPUTS:
%   "match"     the idex to the cells that contain the "searchString"
%   "names"     an optional output, which is the list of file names that
%               correspond to the indicies specified by "match"
%
% C.Hass 02/2014


N = numel(in);
match = cellfun(@regexpi, in, repmat({searchString}, 1, N), repmat({'once'}, 1, N), 'uniformoutput', false);
if iscell(searchString) % complex searches with multiple search targets
    match = cellfun(@(x) ~isempty(cell2mat(x)), match);
else % simple, one word searchString
    match = cellfun(@(x) ~isempty(x), match);
end




% return the optional argument "names" if the optional input "attribute"
% has been specified.
if exist('attribute', 'var')
    names = attribute(match);
end



