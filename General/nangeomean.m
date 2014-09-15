function out = nangeomean(in, dim)

if ~exist('dim', 'var')
    dim = 1;
end


% put the data in cell array
rows = size(in, 1);
cols = size(in, 2);
switch dim
    case 1
        dat = mat2cell(in, rows, ones(1,cols));
    case 2
        dat = mat2cell(in, ones(rows,1), cols);
    otherwise
        error('unknown dimension argument, must be 1 or 2')
end


out = cellfun(@(x) geomean(x(~isnan(x))), dat);