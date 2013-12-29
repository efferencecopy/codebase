%ISVALIDNEXFILENAME    Testing NEX filenames for the proper format.
%    ISVALIDNEXFILENAME(nex_filenames) accepts a cell array or char array
%    of nex filenames and returns a logical vector where '1' means the
%    filename conforms to the following standard:
%
%       X or XX, where X is an uppercase letter
%           nnnnnn, where n is a numerical digit (date: mmddyy)
%               nnn, where n is a numerical digit (file #)
%                   .n, (optional) where n is a numerical digit
%                       .nex, (optional, case insensitive) the file extension
%
%    Note: the function treats whitespace as literal. This means if there
%    is stray whitespace anywhere in the filename, the string is considered
%    an invalid NEX filename
%
%    valid_filenames = ISVALIDNEXFILENAME(nex_filenames), outputs a logical
%    vector of the same size and shape as the input where every true
%    element corresponds to a valid filename in the input.
%
%    valid_filenames = ISVALIDNEXFILENAME(..., verbose), when verbose =
%    true (default), the function prints the invalid filenames to the
%    command window.
%
%    [...,tokens] = ISVALIDNEXFILENAME(...), tokens is a cell array that
%    contains the subject ID, month, day, and year strings for every input
%    filename.
%
% 2012/12/26 zlb

function [valid_nexfilenames,nex_tokens] = isvalidnexfilename(nex_file_strs, verbose)
if nargin < 1, return; end

if nargin < 2 || isempty(verbose)
    verbose = true;
end

if ischar(nex_file_strs) % char arrays -> cell array of strings
    nex_file_strs = cellstr(nex_file_strs);
end

%               (subject id)(    month    )(         day          )(year)(###)(  .#  ) ( file extension )
nexfile_pat = '^([A-Z]{1,2})(0[1-9]|1[012])(0[1-9]|[12][0-9]|3[01])(\d\d)\d{3}(?:\.\d)?(?:\.[nN][eE][xX])?$';
nex_tokens = regexp(nex_file_strs, nexfile_pat, 'tokens', 'once');

valid_nexfilenames = ~cellfun(@isempty,nex_tokens);

if verbose && any(~valid_nexfilenames)
    invalid_filenames = nex_file_strs(~valid_nexfilenames);
    fprintf('The following NEX filenames are invalid:\n');
    fprintf('\t%s\n', invalid_filenames{:});
end
