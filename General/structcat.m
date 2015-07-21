function out = structcat(s, field)

% vertically concatinates the value of the same field across instances of
% structures in a cell array of structures
%
% s should be a cell array of structures
% field should be the field name


out = cellfun(@(x,y) x.(y), s(:), repmat({field}, numel(s), 1), 'uniformoutput', false);
