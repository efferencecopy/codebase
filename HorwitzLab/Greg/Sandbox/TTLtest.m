% Testing TTLtest.d
[fname, pathname] = uigetfile('*.nex', 'Select a NeuroExplorer file');
fid = fopen(strcat(pathname, fname), 'r');

if(fid == -1)
    error('Unable to open file')
    return
end

%open the file and make sure that it's a neuroexplorer file
magic = fread(fid, 1, 'int32');
if magic ~= 827868494
    error 'The file is not a valid .nex file'
end
 
% read in the header info
version = fread(fid, 1, 'int32');
comment = deblank(char(fread(fid, 256, 'char')'));
freq = fread(fid, 1, 'double');                     %acquisition rate for things other than continuous variables
tbeg = fread(fid, 1, 'int32') ./ freq;
tend = fread(fid, 1, 'int32') ./ freq;
nvar = fread(fid, 1, 'int32');
fseek(fid, 260, 'cof'); % skip location of next header and padding


%now read the file putting the appropriate 'type' of Nex info in either
%the ecodes, spikes, or analog vectors. For now ignore the other 'types'

for i=1:nvar
    type         = fread(fid, 1, 'int32');
    varVersion   = fread(fid, 1, 'int32');
    name         = deblank(char(fread(fid, [1 64], 'char')));
    offset       = fread(fid, 1, 'int32');
    n            = fread(fid, 1, 'int32');
    dummy        = fread(fid, 32, 'char'); %sikps what's commented out below
    %WireNumber   = fread(fid, 1, 'int32');
    %UnitNumber   = fread(fid, 1, 'int32');
    %Gain         = fread(fid, 1, 'int32');
    %Filter       = fread(fid, 1, 'int32');
    %XPos         = fread(fid, 1, 'double');
    %YPos         = fread(fid, 1, 'double');
    WFrequency   = fread(fid, 1, 'double'); % wf sampling fr.
    ADtoMV       = fread(fid, 1, 'double'); % coeff to convert from AD values to Millivolts.
    NPointsWave  = fread(fid, 1, 'int32'); % number of points in each wave
    NMarkers     = fread(fid, 1, 'int32'); % how many values are associated with each marker
    MarkerLength = fread(fid, 1, 'int32'); % how many characters are in each marker value
    MVOfffset    = fread(fid, 1, 'double'); % coeff to shift AD values in Millivolts: mv = raw*ADtoMV+MVOfffset
    filePosition = ftell(fid);
    if (type == 5)
        name
        i
        WFrequency;
        fseek(fid, offset, 'bof');
        
        % get the start times of each analog signal snippet
        fragStarts = fread(fid, [n 1], 'int32')./freq;
        
        %get the number of points sammpled during each of the snippets
        nFrag = fread(fid, [n 1], 'int32');
        nFrag(n+1) = NPointsWave;
        nPtsPerFrag = nFrag;
        %now bring in the AtoD data for the entire recording
        data = fread(fid, [NPointsWave 1], 'int16').*ADtoMV + MVOfffset;
        
        fseek(fid, filePosition, 'bof');
    else
        dummy = fread(fid, 60, 'char'); %skip some junk to get to the next block of data
       % Do nothing 
    end
end
fclose(fid);

% Now doing the analysis
idxs = find(diff(data>4))+1; % starts;
tmp = reshape(idxs,2,length(idxs)/2);

% Plotting all the messages on top of each other
firstpulseinblock = tmp(1,find(diff(tmp(1,:)) > 500)+1);
figure; axes; hold on;
for i = 1:length(firstpulseinblock)
    plot([-10:400]/WFrequency*1000,data([-10:400]+firstpulseinblock(i)),'k-')
end
xlabel('Time (ms)'); ylabel('Voltage (V)');
set(gca,'Xlim',[-10, 400]/WFrequency*1000);
set(gca,'Ylim',[-.5 5.5]);
title(['N = ',num2str(length(firstpulseinblock))]);
