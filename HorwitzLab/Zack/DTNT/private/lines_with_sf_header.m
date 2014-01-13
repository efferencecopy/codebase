function [idx, sfs] = lines_with_sf_header(lines)
idx = find(~cellfun(@isempty, strfind(lines, 'sf:')));
assert(idx(1) == true, ['The text file is malformed: "sf:<number>" must' ...
    ' be the first (uncommented) entry in the text file!']);
sfs = sscanf(char(lines(idx))', 'sf:%f\n');
