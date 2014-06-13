function markLocalFilesAsCopies(startDir, trashDir)


global GL_DATPATH

% save the top level directory so that recursive operations work
presDir = pwd;


% save the top level directory so that recursive operations work
if ~exist('startDir', 'var');
    startDir = 'C:\Users\glickfeld_lab\Desktop\Local Data Files';
end
cd(startDir)


if ~exist('trashDir', 'var');
    mkdir('filesOnCrash');
    trashDir = [startDir, filesep, 'filesOnCrash'];
end


% load the mdb
mdb = initMouseDB(false, false);

% loop over all the files
d = dir;
for a = 1:numel(d)
    
    fprintf('file %d of %d\n', a, numel(d));
    
    if any(strcmpi(d(a).name, {'.', '..', 'DTcones', 'filesOnCrash'}))
        fprintf('skiping %s \n', d(a).name);
        continue
    end
    
    % the recursive part
    if d(a).isdir
        markLocalFilesAsCopies(d(a).name, trashDir)
    end
    
    
    % look for the file (the quick way);
    fileparts = regexp(d(a).name, '.abf', 'split');
    mouseName = mdb.search(fileparts{1});
    if ~isempty(mouseName);
        tmpPath = findfile(fileparts{1}, [GL_DATPATH, mouseName{1}], '.abf');
    else
        tmpPath = [];
    end
    
    
    % look for it the slow way
    if isempty(tmpPath)
        tmpPath = findfile(fileparts{1}, GL_DATPATH, '.abf');
    end
   
    % if a path has been found, than the file exists on Crash. Move the
    % file to the trashDir
    if ~isempty(tmpPath)
        fprintf('File %s is on Crash\n', d(a).name);
        movefile([startDir, filesep, d(a).name], [trashDir, filesep, d(a).name]);
    end
end

% be nice and CD back to the original dir
cd(presDir)




