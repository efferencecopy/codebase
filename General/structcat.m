function out = structcat(s, field)

% vertically concatinates the value of the same field across instances of
% structures in a cell array of structures
%
% s should be a cell array of structures
% field should be the field name or cell array of nested names

if iscell(field)
    N = numel(field);
else
    N = 1;
end

switch N
    case 1
        out = cellfun(@(x,y) x.(y), s(:), repmat({field}, numel(s), 1), 'uniformoutput', false);
    case 2
        out = cellfun(@(x,y) x.(y{1}).(y{2}), s(:), repmat({field}, numel(s), 1), 'uniformoutput', false);
    case 3
        out = cellfun(@(x,y) x.(y{1}).(y{2}).(y{3}), s(:), repmat({field}, numel(s), 1), 'uniformoutput', false);
    case 4
        out = cellfun(@(x,y) x.(y{1}).(y{2}).(y{3}).(y{4}), s(:), repmat({field}, numel(s), 1), 'uniformoutput', false);
    otherwise
        error('indent level is too deep')
end
