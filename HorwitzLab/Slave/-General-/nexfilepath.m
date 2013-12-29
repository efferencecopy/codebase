%NEXFILEPATH   Return the mounted location of our NEX files.

% ZALB 2013/10/07

function filepath = nexfilepath()
filepath = ''; partial_path = '';

host_ip = '128.95.153.12';
win_cmd_str = ['net view ' host_ip];
nix_cmd_str = ['mount | grep ' host_ip ' | sed -n ''s/.*on //;s/ %s.*//p'''];

if strcmp(license, '367516')
    filepath = 'C:\NO BACKUP\NexFiles\';
    return
end

%#ok<*ASGLU>
if ispc % assumes the PC mounted "NO BACKUP" as a drive
    [nil,partial_path] = system(win_cmd_str);
    colon_pos = find(partial_path == ':');
    partial_path = partial_path(colon_pos-1:colon_pos);
elseif ismac % macs have a different `mount` format
    [nil,partial_path] = system(sprintf(nix_cmd_str, '('));
elseif isunix
    [nil,partial_path] = system(sprintf(nix_cmd_str, 'type'));
end

if ~isempty(partial_path)
    filepath = [strtrim(partial_path) filesep 'NexFiles' filesep];
end
