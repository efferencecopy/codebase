function out = findGTfiles

NEXFILEPATHS = 'nexPaths.mat';
if isunix
    prefix = '/Users/charliehass/LabStuff/NexFiles/Charlie/Batch Data And Text Files/';
elseif ispc
    prefix = 'N:\NexFiles\Charlie\Batch Data And Text Files\';
end
nexPathsPath = [prefix, NEXFILEPATHS];
load(nexPathsPath);

%   FIND VALID GT FILES
out = {};
counter = 1;
for a = 1:length(nexPaths.paths)
    fname = nexPaths.names{a};
    fpath = nexPaths.paths{a};
    
    %make sure that it's my file.
    if isempty(regexpi(fpath, 'charlie'));
        continue
    end
    
    %make sure that the data come from expt using the new LFP board
    exptDate = fname(2:7);
    month = str2num(exptDate(1:2));
    day = str2num(exptDate(3:4));
    year = str2num(exptDate(5:6));
%     try
%         if (year<10) || ((year<=10) && (month<9)) || ((year<=10) && (month<=9) && (day<28))
%             fprintf('Skipping <%s>\n', fname);
%             continue
%         end
%     catch
%         continue
%     end
    
    %make sure it's a graings file
    if getparadigmID(fpath) ~= 150;
        continue
    end
    out{counter} = fname;
    counter = counter+1;    
end


