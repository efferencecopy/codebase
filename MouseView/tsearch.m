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
bool_AND = ~isempty(regexp(searchString, 'AND', 'once'));
bool_OR = ~isempty(regexp(searchString, 'OR', 'once'));

% figure out if it's a boolean search, parse the string terms, and then do
% the search
if bool_AND || bool_OR
    terms = strsplit(searchString);
    
    % remove the boolean terms
    idx_and = cellfun(@(x) ~isempty(regexp(x, 'AND', 'once')), terms);
    idx_or = cellfun(@(x) ~isempty(regexp(x, 'OR', 'once')), terms);
    terms(idx_and | idx_or) = [];
    
    % run the search
    match = false(N, numel(terms));
    for a = 1:numel(terms)
        tmp = cellfun(@regexpi, in, repmat(terms(a), 1, N), repmat({'once'}, 1, N), 'uniformoutput', false);
        tmp = cellfun(@(x) ~isempty((x)), tmp);
        match(:,a) = tmp;
    end
    
    % perform the boolean operation
    if bool_AND
        match = sum(match,2) == numel(terms);
    elseif bool_OR
        match = sum(match,2)>0;
    end
    
else %simple non-boolean searches
    
    match = cellfun(@regexpi, in, repmat({searchString}, 1, N), repmat({'once'}, 1, N), 'uniformoutput', false);
    match = cellfun(@(x) ~isempty(x), match);
    match = match(:);
end



% return the optional argument "names" if the optional input "attribute"
% has been specified.
if exist('attribute', 'var')
    names = attribute(match);
end



