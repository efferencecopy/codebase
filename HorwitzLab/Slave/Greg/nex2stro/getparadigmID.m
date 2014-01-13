function paradigmID = getparadigmID(filename)
% paradigmID = getparanum(filename)
% This function will extract the paradigm number from a .nex file.
% GDLH 1/13/10

    paradigmID = nan;
    if (nargin == 0 || isempty(filename))
        [fname, pathname] = uigetfile('*.nex', 'Select a NeuroExplorer file');
        filename = strcat(pathname, fname);
    end

    if isunix
        fid = fopen(filename, 'r', 'l');
    elseif ispc
        fid = fopen(filename, 'r');
    end

    if(fid == -1)
        error('Unable to open file')
    end

    %open the file and make sure that it's a neuroexplorer file
    magic = fread(fid, 1, 'int32');
    if magic ~= 827868494
        error 'The file is not a valid .nex file'
    end

    %read in the header info
    version = fread(fid, 1, 'int32');
    comment = deblank(char(fread(fid, 256, 'char')'));
    freq = fread(fid, 1, 'double');                     %acquisition rate for things other than continuous variables
    tbeg = fread(fid, 1, 'int32') ./ freq;
    tend = fread(fid, 1, 'int32') ./ freq;
    nvar = fread(fid, 1, 'int32');
    fseek(fid, 260, 'cof'); % skip location of next header and padding

    for i=1:nvar
            type         = fread(fid, 1, 'int32');
            varVersion   = fread(fid, 1, 'int32');
            name         = deblank(char(fread(fid, [1 64], 'char')));
            offset       = fread(fid, 1, 'int32');
            n            = fread(fid, 1, 'int32');
            dummy        = fread(fid, 32, 'char'); %skips what's commented out below
            WFrequency   = fread(fid, 1, 'double'); % wf sampling fr.
            ADtoMV       = fread(fid, 1, 'double'); % coeff to convert from AD values to Millivolts.
            NPointsWave  = fread(fid, 1, 'int32'); % number of points in each wave
            NMarkers     = fread(fid, 1, 'int32'); % how many values are associated with each marker
            MarkerLength = fread(fid, 1, 'int32'); % how many characters are in each marker value
            MVOffset    = fread(fid, 1, 'double'); % coeff to shift AD values in Millivolts: mv = raw*ADtoMV+MVOfffset
 
        if (type == 6)
            ecodes = [];
            fseek(fid, offset, 'bof');
            timeStamps = fread(fid, [n 1], 'int32') ./ freq;

            %check to make sure that there are actually markers to
            %retreive before trying
            if (NMarkers ~= 1)
                error('bad number of markers in file');
            end

            %check to make sure that you're about to read in strobed
            %codes
            mname = char(fread(fid, [1 64], 'char'));
            if ~strncmp('DIO', mname, 3)
                error('unknown marker name')
            end

            %now collect the marker (ecode) data. convert them to
            %numeric representations and return them along with timestamps
            markers = deblank(fscanf(fid, '%c', [MarkerLength n])'); %transpose, then take off the trailing blank
            markers = markers-48; %convert to numeric
            powersOfTen = 10.^[(MarkerLength-2) : -1 : 0]';
            ecodes = [timeStamps, markers*powersOfTen];
            paradigmID = cell2mat(dat2num(ecodes, 8999, 4000, 'int', 0));
        end
        dummy = fread(fid, 60, 'char'); %skip some junk to get to the next block of data
    end
    fclose(fid);
end