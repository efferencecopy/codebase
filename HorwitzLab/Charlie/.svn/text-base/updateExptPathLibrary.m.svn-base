function nexPaths = updateExptPathLibrary(startDir, nexPaths)

%extracts all the nex files from the NexFiles directory and stores their
%absolute paths. This should speed things up when doing batch processes.

if ~exist('startDir', 'var')
    startDir = uigetdir;
end
if ~exist('nexPaths', 'var')
    nexPaths.names = {};
    nexPaths.paths = {};
end


% get the contents of the directory.
currentDir = pwd; %cd back to this later.
fprintf('Looking in: <%s>\n', startDir);
cd(startDir);
d = dir;

%append any .nex files into the exptPathLibrary.
%keyboard
names = {d(:).name}';
pat = repmat({'\.nex'}, size(names));
l_nexFiles = cellfun(@(x,y) numel(regexpi(x,y)), names, pat);
nexInDirNames = names(logical(l_nexFiles));
trunc_pwd = pwd;
parent_path_idx = strfind(trunc_pwd, 'NexFiles');
if any(l_nexFiles) && ~isempty(parent_path_idx)
    trunc_pwd(1:parent_path_idx-1) = [];
end
% store relative paths and replace file separator characters with a bogus sequence
% the bogus sequence will be undone when a user polls the nexPaths structure via findNexPath -- Zack
nexInDirPaths = cellfun(@(x) strrep([filesep trunc_pwd filesep x], filesep, '$:$'), nexInDirNames, 'uniformoutput', 0);
nexPaths.names = [nexPaths.names; [nexInDirNames(:)]];
nexPaths.paths = [nexPaths.paths; [nexInDirPaths(:)]];


directories = [d(:).isdir];
for i = find(directories);
    if any(strcmp(d(i).name, {'.', '..', '.DS_Store'}))
        continue % don't recursively search through these directories
    end
    
    %drop into the directory, and recursively search
    nexPaths = updateExptPathLibrary([pwd, filesep, d(i).name], nexPaths);
    cd(startDir)
end

%be nice and cd back to the original directory
if strcmpi(license, '380245') % charlie's laptop
    cd('/Users/charliehass/LabStuff/NexFiles/Charlie/Batch Data And Text Files');
elseif strcmpi(license, '633208') % rig 3 mac
    cd '~/Desktop/Charlie Cone Model Stuff'
elseif strcmp(license, '367516') %lab server
    cd('C:\NO BACKUP\NexFiles\Charlie\Batch Data And Text Files');
elseif ispc
    cd('N:\NexFiles\Charlie\Batch Data And Text Files');
end
save 'nexPaths.mat' nexPaths
cd(currentDir)


