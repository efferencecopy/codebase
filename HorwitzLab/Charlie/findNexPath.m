function filepath = findNexPath(nexPaths, fname)
    %uses a pre-existing mapping of nexfile names to their absolute paths
    %to return the absolute path. Should speed things up when looking for
    %files over the network.
    
    [~, fname] = fileparts(fname); % in case the user passes an absolute path
    fname = [fname '.nex'];
    
    idx = find(strcmpi(nexPaths.names, fname), 2); % find at most 2 true entries
    
    %package the result.
    if ~isempty(idx)
        if length(idx) > 1
            filepath = fname;
            fprintf('File <%s> has none or more than one entry in NexFiles directory\n', fname);
        else
            filepath = strrep(nexPaths.paths{idx}, '$:$', filesep);
%             switch whoami
%                 case 'nuke'
%                     filepath = filepath(2:end); % kinda dumb that the switch/case doesn't do anything...
%             end
            
            %%% OLD VERSION %%%
%             if strcmp(license, '367516') % lab server
%                 filepath = ['C:\NO BACKUP' filepath];
%             elseif ispc % all other lab pc's
%                 filepath = ['N:' filepath];
%             elseif strcmp(license, '380245') % charlie's laptop
%                 filepath = ['/Users/charliehass/LabStuff' filepath];
%             elseif isunix && ~ismac %shadlen cluster
%                 filepath = filepath;                
%             end
        end
    else
        fprintf('    ******  The experiment library needs to be updated!!!  ******\n');
        filepath = findfile(fname);
    end
end