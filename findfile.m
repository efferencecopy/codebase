function outname = findfile(inname, optpath, optsuffix)
%
%    CHARLIE'S LAPTOP VERSION, ONLY GETS CALLED ON MY LAPTOP!!!!!
%
%   This function will take name a nex file name (e.g. "K010309002")
% and find it in the data directory tree (and will append a ".nex"
% if necessary).  The output is the full name (with the path) which can
% be the passed to nex2stro.
%
%   The optional second argument, optpath, can be used to indicate where
% in a directory tree you want to start, but it's mostly for the
% recursion.
%
% GDLH  1/10/09
% CAH   4/2012      Added some functionality for different file types (suffixes)

outname = [];

if (iscell(inname))
    inname = char(inname);
end

if exist('optsuffix', 'var')
    suffix = optsuffix;
else
    suffix = '.NEX';
end

if (length(inname) < length(suffix))
    inname = [inname, suffix];
end
if (~strcmpi(inname(1,end-length(suffix)+1:end),suffix))
    inname = [inname, suffix];
end

if (nargin == 1)
    optpath = '/Users/charliehass/LabStuff/NexFiles';
end

a = dir(optpath);
tmppathlist = {};
for i = 1:size(a,1)
    if (a(i).isdir && ~strcmp(a(i).name,'.') && ~strcmp(a(i).name,'..'))
        tmppathlist{length(tmppathlist)+1} = a(i);
    end
    if (~a(i).isdir)
        if (strcmpi(a(i).name,inname))
            outname = [optpath,'/',a(i).name];
            return;
        end
    end
end

% Recursive part
for i = 1:length(tmppathlist)
    if (isempty(outname))
        outname = findfile(inname, [optpath,'/',tmppathlist{i}.name], suffix);
    end
end
