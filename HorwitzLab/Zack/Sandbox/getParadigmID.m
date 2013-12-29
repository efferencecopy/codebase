% standalone paradigm identification number getter script
% takes an absolute path to a nex file as an argument
% outputs the paradigm ID as a double
% zlb 2011-04-18
function paradigmID = getParadigmID(filename)
paradigmID = nan;
if nargin == 0 || isempty(filename)
    [fname, pathname] = uigetfile('*.nex', 'Select a NeuroExplorer file');
    filename = strcat(pathname, fname);
end

if isunix
    fid = fopen(filename, 'r', 'l');
elseif ispc
    fid = fopen(filename, 'r');
end

C = onCleanup(@() fclose(fid));

if fid == -1
    error('Unable to open file')
end

magic = fread(fid, 1, 'int32');
if magic ~= 827868494
    error('The file is not a valid .nex file');
end

fread(fid, 260, 'char'); % skip
freq    = fread(fid, 1, 'double');
fread(fid, 8, 'char'); % skip
nvar    = fread(fid, 1, 'int32');
fread(fid, 260, 'char'); % skip

ecodes = [];

for i=1:nvar
    type         = fread(fid, 1, 'int32');
    fread(fid, 68, 'char'); % skip
    offset       = fread(fid, 1, 'int32');
    n            = fread(fid, 1, 'int32');
    fread(fid, 52, 'char'); % skip
    NMarkers     = fread(fid, 1, 'int32'); % how many values are associated with each marker
    MarkerLength = fread(fid, 1, 'int32'); % how many characters are in each marker value
    fread(fid, 8, 'char'); % skip
    filePosition = ftell(fid);
    
    if type == 6 % ecodes
        fseek(fid, offset, 'bof');
        timeStamps = fread(fid, [n 1], 'int32') ./ freq;
        
        if (NMarkers ~= 1)
            error('bad number of markers in file');
        end
        
        mname = fread(fid, [1 64], '*char');
        if ~strncmp('DIO', mname, 3)
            error('unknown marker name');
        end
        
        markers = deblank(fscanf(fid, '%c', [MarkerLength n])'); %transpose, then take off the trailing blank
        markers = markers-48; %convert to numeric
        powersOfTen = 10.^((MarkerLength-2) : -1 : 0)';
        ecodes = [timeStamps markers*powersOfTen];
        %return to where you started in the file!
        fseek(fid, filePosition, 'bof');
    end
    fread(fid, 60, 'char');
end
C.delete;
cellout = dat2num(ecodes, 8999, 4000, 'int', 0);
paradigmID = cellout{1};