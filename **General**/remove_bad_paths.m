function remove_bad_paths()
bad_path_pattern = ...
    ['\' filesep '([a-zA-Z]\w*\.bundle|\.svn|\.git)(\' filesep '|\' pathsep ')'];

cmdirs = regexp([matlabpath pathsep], ['.[^' pathsep ']*' pathsep], 'match')';
goodpaths = cmdirs(cellfun(@isempty, regexp(cmdirs, bad_path_pattern)));
path([goodpaths{:}]);

warning('Save the path so the changes persist between sessions!');
