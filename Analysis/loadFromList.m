function ax = loadFromList(list)
%
% loads a list of ABF files. The list should be contained in a cell array.
% There are two formats: actual data file names, and:
%   {datePrefix, [numbers]}
% where 'datePrefix' is the string corresonding to the date the data were
% collected, and 'numbers' is the numeric value of the data file (i.e., the
% suffix).
%


% do some error checking. This function tries to be backwards compatiable
% so there will be a few steps.
assert(iscell(list), 'ERROR: input needs to be a cell array');

% determine which form the params.files is in. It could be actual file names, or
% the date prefix, and a list of numbers indicating which file that day.
nchars = cellfun(@numel, list);

fnames = {};
if any(nchars ~= 15) % [prefix, suffix] construction for 'list' input
    for a = 1:numel(list{2})
        suffix = num2str(list{2}(a));
        nZerosNeeded = 4-numel(suffix);
        suffix = [repmat('0',1,nZerosNeeded), suffix];
        fnames{a} = [list{1}, suffix];
    end
else
    fnames = list;
end


ax = {};
fprintf('**** Unpacking %d files:\n', numel(fnames));
for a = 1:numel(fnames)
    fprintf('file %d: %s \n', a, fnames{a})
    ax{a} = abfobj(fnames{a});     
end