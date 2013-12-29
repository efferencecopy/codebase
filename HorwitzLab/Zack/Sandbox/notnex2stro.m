%NOTNEX2STRO   Convert NeuroEXplorer files to a STRO structure.
%   STRO = NOTNEX2STRO() prompts the user for the location of a .nex file and
%   returns a MATLAB structure with fields 'sum', 'trial', 'ras', and 'other'.
%   These fields are populated with information from the .nex file and contain
%   data summaries, trial-by-trial parameters, rasters, and any other data. The
%   STRO structure is the standard MATLAB-based representation of all
%   REX-acquired data (plus spike times) in the Horwitz lab.
%
%   STRO = NOTNEX2STRO(FILENAME) is the same as the above but instead uses the
%   .nex file specified by the absolute path string FILENAME.

% cah       1/21/08 started development...
% gdlh/cah  2/14/08 replaced global variable architecture with nested
%                   functions.
% zalb     10/09/13 removing all non-local variable scopes
    
function stro = notnex2stro(file_name)
if nargin == 0 || isempty(file_name)
    [fname, pathname] = uigetfile('*.nex', 'Select a NeuroEXplorer file');
    if isequal(fname, 0) || isequal(pathname, 0), return; end % user cancelled
    file_name = [pathname fname];
elseif ~exist(file_name, 'file')
    error('Couldn''t find file %s\n', file_name);
end

[stro, ecodes, spikes, anlg, waves, events] = setupStructures();
stro.sum.fileName = file_name;

waitdlg = waitbar(0, 'Opening a .nex file');

% Unpack the .nex file. No arguments are needed b/c all the local variables
% are accessible to the nested functions.
[ecodes, spikes, anlg, waves, events] = unpackNexFile(stro.sum.fileName, ecodes, spikes, anlg, waves, events, waitdlg);
[anlg, events] = cleanupAnalogSignals(anlg, events); % determine if the start times are identical for each channel

% determine which paradigm file to use, then eval the appropriate header
% file. Be sure to return all the appropriate variables (for scoping
% purposes). the feval statement will crash if one of the outputs is
% undefined in paradigmFile
paradigmID = dat2num(ecodes, 8999, 4000, 'int', 0);
stro.sum.paradigmID = paradigmID{1}(1);
paradigmFile = paradigmLibrary(stro.sum.paradigmID);

[C, trialParamsHeader, exptParamsHeader, additionalRasterFields, badLogic, otherInstructions] = feval(paradigmFile);

% Fill up some of the stro.sum contents
stro = getExperimentalParams(stro, ecodes, C.VALOFFSET, exptParamsHeader); % iterates over the exptParamsHeader reconstructing those elements
stro.sum.rasterCells = makeRasterHeader(spikes, anlg, waves, additionalRasterFields); % defined in the 'paradigmFile'.
stro.sum.trialFields = trialParamsHeader; % defined in the 'paradigmFile'
stro.sum.otherFields = otherInstructions;
stro.sum.analog.sigid = anlg.names;
stro.sum.analog.storeRates = anlg.ADFrequency;
stro.sum.analog.ADtoMV = anlg.ADtoMV;
stro.sum.waves.wf_id = waves.names;
stro.sum.waves.storeRates = waves.WFrequency;

% iterate by trials filling up the .raster and .trial fields of stro.
nTrials = length(events.start);
stro.ras = cell(nTrials, length(stro.sum.rasterCells));
stro.trial = nan(nTrials, size(trialParamsHeader, 2));
stro.sum.absTrialNum = nan(nTrials,1);
goodTrial = 0; % the counter of good trials

for b = 1:nTrials % needs descriptive var name
    trialCodes = checkTrial(ecodes, events.start(b), events.stop(b), C, badLogic);
    if ~isempty(trialCodes)
        goodTrial = goodTrial + 1;
        stro = rasterByTrial(stro, spikes, anlg, waves, C.VALOFFSET, events.start(b), events.stop(b), goodTrial, b, additionalRasterFields, trialCodes);
        stro = indexByTrial(stro, trialCodes, trialParamsHeader, goodTrial, C.VALOFFSET);
        stro = otherByTrial(stro, trialCodes, otherInstructions, goodTrial, C);
        stro.sum.absTrialNum(goodTrial) = b; % a mapping b/w actual and good trial nums
    end
    
    if ~rem(b, 3) % update the waitdlg every three trials
        waitbar(b ./ nTrials, waitdlg,  '   ...Processing trials');
    end
end

% now cleanup the bad trials (empty rows) of stro.trial, stro.ras, and
% stro.sum.absTrialNum
waitbar(1, waitdlg,  '   ... cleaning up');
stro.ras((goodTrial+1):end, :) = [];
stro.trial((goodTrial+1):end, :) = [];
stro.sum.absTrialNum(goodTrial+1:end) = [];

% close up the wait dialog
close(waitdlg);

function [stro, ecodes, spikes, anlg, waves, events] = setupStructures()
% this function is necessary because later on we'll test the length of
% spikes.names and anlg.names to iterate over neuronal and anlg data. If
% we decide not to record from an anlg channel, setting the field
% anlg.name to an empty vector (i.e., []) will allow this program to
% function properly. ditto for spike channels. In addition, these
% variables are global, so previous definitions of the variables could
% intrude here.

ecodes = [];

spikes = struct( ...
            'names', {{}}, ...
            'timestamps', {{}} ...
         );

anlg = struct( ...
          'names', {{}}, ...
          'ADFrequency', {{}}, ...
          'ADtoMV', {{}}, ...
          'fragStarts', {{}}, ...
          'nPtsPerFrag', {{}}, ...
          'data', {{}} ...
       );

waves = struct( ...
           'names', {{}}, ...
           'WFrequency', {{}}, ...
           'timestamps', {{}}, ...
           'waveforms', {{}} ...
        );

events = struct();

stro = struct('sum', struct( ...
                        'fileName', [], ...
                        'date', date, ...
                        'paradigmID', [], ...
                        'waves', struct( ...
                                    'wf_id', [], ...
                                    'storeRates', [] ...
                                 ), ...
                        'analog', struct( ...
                                     'sigid', [], ...
                                     'storeRates', [], ...
                                     'ADtoMV', [] ...
                                  ), ...
                        'trialFields', {{}}, ...
                        'rasterCells', {{}}, ...
                        'otherFields', {{}}, ...
                        'exptParams', struct(), ...
                        'absTrialNum', [], ...
                        'concat', [] ...
                     ), ...
              'trial', {{}}, ...
              'ras', {{}}, ...
              'other', {{}} ...
       );

function [ecodes, spikes, anlg, waves, events] = unpackNexFile(fileName, ecodes, spikes, anlg, waves, events, waitdlg)
% open differently for mac and pc:
if isunix
    fid = fopen(fileName, 'r', 'l');
else
    fid = fopen(fileName, 'r');
end

if fid == -1
    error('Unable to open file');
end

% open the file and make sure that it's a neuroexplorer file
magic = fread(fid, 1, 'int32');
if magic ~= 827868494
    error('The file is not a valid .nex file');
end

% read in the header info
fread(fid, 1, 'int32'); % version
fread(fid, 256, 'char'); % comment
freq = fread(fid, 1, 'double'); % acquisition rate for things other than continuous variables
fread(fid, 1, 'int32'); % ./ freq; % tbeg
fread(fid, 1, 'int32'); % ./ freq; % tend
nvar = fread(fid, 1, 'int32');
fseek(fid, 260, 'cof'); % skip location of next header and padding

% now read the file putting the appropriate 'type' of Nex info in either
% the ecodes, spikes, or analog vectors. For now ignore the other
% 'types'

neuronCount = 0;
contSigCount = 0;
waveCount = 0;

for i = 1:nvar
    type = fread(fid, 1, 'int32');
    fread(fid, 1, 'int32'); % varVersion
    name = deblank(fread(fid, [1 64], '*char'));
    offset = fread(fid, 1, 'int32');
    n = fread(fid, 1, 'int32');
    fread(fid, 32, 'char'); % skips what's commented out below
    %WireNumber = fread(fid, 1, 'int32');
    %UnitNumber = fread(fid, 1, 'int32');
    %Gain = fread(fid, 1, 'int32');
    %Filter = fread(fid, 1, 'int32');
    %XPos = fread(fid, 1, 'double');
    %YPos = fread(fid, 1, 'double');
    WFrequency = fread(fid, 1, 'double'); % wf sampling fr.
    ADtoMV = fread(fid, 1, 'double'); % coeff to convert from AD values to Millivolts.
    NPointsWave = fread(fid, 1, 'int32'); % number of points in each wave
    NMarkers = fread(fid, 1, 'int32'); % how many values are associated with each marker
    MarkerLength = fread(fid, 1, 'int32'); % how many characters are in each marker value
    MVOffset = fread(fid, 1, 'double'); % coeff to shift AD values in Millivolts: mv = raw*ADtoMV+MVOffset
    filePosition = ftell(fid);
    
    switch type
        case 0 % neuron
            waitbar(0, waitdlg, '  ...adding spike data');
            neuronCount = neuronCount+1;
            spikes.names(neuronCount) = {name};
            fseek(fid, offset, 'bof');
            spikes.timestamps{neuronCount} = fread(fid, [n 1], 'int32') ./ freq;
            fseek(fid, filePosition, 'bof');
            
        case 1 % event (start/stop)
            waitbar(0, waitdlg, sprintf('  ...adding %s times', name));
            fseek(fid, offset, 'bof');
            ts = fread(fid, [n 1], 'int32') ./ freq;
            events.(lower(name)) = ts;
            fseek(fid, filePosition, 'bof');
            
        case 2 % interval
            % ignore this type for now
            fprintf('Ignoring the ''interval'' type\n');
            
        case 3 % waveform
            waveCount = waveCount+1;
            waves.names(waveCount) = {name};
            waves.WFrequency{waveCount} = WFrequency;
            fseek(fid, offset, 'bof');
            waves.timestamps{waveCount} = fread(fid, [n 1], 'int32') ./ freq;
            waves.waveforms{waveCount} = (fread(fid, [NPointsWave n], 'int16') .* ADtoMV + MVOffset)'; % transpose to make trials go down columns and time accross rows
            fseek(fid, filePosition, 'bof');
            
        case 4 % population vector
            % ignore this type for now
            fprintf('Ignoring the ''population vector'' type\n');
            
        case 5 % continuous variable (i.e. eye signals)
            waitbar(0, waitdlg, sprintf('  ...adding analog data: %s', name));
            contSigCount = contSigCount+1;
            anlg.names(contSigCount) = {name};
            anlg.ADtoMV(contSigCount) = {ADtoMV};
            anlg.ADFrequency{contSigCount} = WFrequency;
            fseek(fid, offset, 'bof');
            
            % get the start times of each analog signal snippet
            anlg.fragStarts{contSigCount} = fread(fid, [n 1], 'int32') ./ freq;
            
            % get the number of points sammpled during each of the snippets
            nFrag = fread(fid, [n 1], 'int32');
            nFrag(n+1) = NPointsWave;
            anlg.nPtsPerFrag{contSigCount} = nFrag;
            % now bring in the AtoD data for the entire recording
            anlg.data{contSigCount} = fread(fid, [NPointsWave 1], 'int16') .* ADtoMV + MVOffset;
            fseek(fid, filePosition, 'bof');
            
        case 6 % marker (i.e. ecodes)
            waitbar(0, waitdlg, '  ...adding ecodes');
            fseek(fid, offset, 'bof');
            timeStamps = fread(fid, [n 1], 'int32') ./ freq;
            
            % check to make sure that there are actually markers to
            % retreive before trying
            if NMarkers ~= 1
                error('bad number of markers in file');
            end
            
            % check to make sure that you're about to read in strobed
            % codes
            mname = fread(fid, [1 64], '*char');
            if ~strncmp('DIO', mname, 3)
                error('unknown marker name');
            end
            
            % now collect the marker (ecode) data. convert them to
            % numeric representations and return them along with timestamps
            markers = deblank(fscanf(fid, '%c', [MarkerLength n])'); % transpose, then take off the trailing blank
            markers = markers - 48; % convert to numeric
            powersOfTen = 10.^((MarkerLength-2):-1:0)';
            ecodes = [ecodes; [timeStamps markers * powersOfTen]];
            
            fseek(fid, filePosition, 'bof');
            
        otherwise
            fprintf('unknown variable type <%s>\n', type);
    end
    fread(fid, 60, 'char'); % skip some junk to get to the next block of data
end
fclose(fid);

function trialCodes = checkTrial(ecodes, trialStart, trialStop, C, badLogic) %#ok<INUSL> % need access to C for `eval`
% pull out the relavent ecodes
ind = ecodes(:,1) >= trialStart & ecodes(:,1) <= trialStop;
trialCodes = ecodes(ind,:);

% test the ecodes to determine if the trial is good or bad. 'badLogic'
% should be defined in the paradigm header file.
bad = eval(badLogic);

% return an empty vector if it's a bad trial
if bad
   trialCodes = [];
end

function stro = rasterByTrial(stro, spikes, anlg, waves, VALOFFSET, start, stop, gtCounter, totTrialCounter, additionalRasterFields, trialCodes)
% fill up the rows of stro.ras trial by trial. treat sike data, anlg
% sigs, and other independently. Format is:
% <neuron 1> <neuron 2> ... <waveForm1>... <anlg 1> <anlg 2> ... <atime> <other>

% each neuron's spikes will go into it's own column:
if ~isempty(spikes.names)
    for ind = 1:length(spikes.names)
        stro.ras{gtCounter, ind} = spikes.timestamps{ind}(spikes.timestamps{ind} >= start & spikes.timestamps{ind} <= stop);
    end
else
    ind = 0;
end

% iterate over the waveform data (if available):
for a = 1:length(waves.names)
    ind = ind + 1;
    whichSpikes = (waves.timestamps{a} >= start & waves.timestamps{a} <= stop);
    stro.ras{gtCounter, ind} = waves.waveforms{a}(whichSpikes, :);
end

% iterate over anlg channels
for a = 1:length(anlg.names)
    ind = ind + 1;
    startInd = anlg.nPtsPerFrag{a}(totTrialCounter) + 1;
    stopInd = anlg.nPtsPerFrag{a}(totTrialCounter+1);
    stro.ras{gtCounter, ind} = anlg.data{a}(startInd:stopInd);
end

% explicitly add the anlg start time. All the analog channels should
% have the same number of start times, so indexing this way is o.k.
if ~isempty(anlg.fragStarts)
    ind = ind + 1;
    stro.ras{gtCounter, ind} = anlg.fragStarts{1}(totTrialCounter);
end

% now look for additional elements as specified by
% additionalRasterFields (which is defined in the 'paradigmFile'
for a = 1:size(additionalRasterFields, 2)
    ind = ind + 1;
    nums = dat2num(trialCodes(:,2), [additionalRasterFields{4, a}], VALOFFSET, additionalRasterFields{2, a}, additionalRasterFields{3, a});
    stro.ras{gtCounter, ind} = nums{1};
end

function stro = indexByTrial(stro, trialCodes, trialParamsHeader, goodTrial, VALOFFSET)
% deal with the ints, longs, and doubles
types = {'int', 'long', 'float', 'double'};
for a = 1:length(types)
    typeInd = strcmp(types{a}, trialParamsHeader(2,:));
    if any(typeInd)
        nums = dat2num(trialCodes(:,2), [trialParamsHeader{4, typeInd}], VALOFFSET, types{a}, 0);
        stro.trial(goodTrial, typeInd) = deal([nums{:}]);
    end
end

% now deal with the 'time' variables
timeInd = find(strcmp('time', trialParamsHeader(2,:)));
timeInd_code = [trialParamsHeader{4, timeInd}];
for a = 1:length(timeInd)
    time = trialCodes(trialCodes(:,2) == timeInd_code(a), 1); % there will be an error if numel(time) > 1
    if ~isempty(time)
        if numel(time) > 1, warning('More than one time code found: using the last!!'); end
        stro.trial(goodTrial, timeInd(a)) = time(end);
    end
end

function stro = otherByTrial(stro, trialCodes, otherInstructions, gtCounter, codeStruct)
ind = 0;
for a = 1:length(otherInstructions)
    ind = ind+1;
    stro.other{gtCounter, a} = feval(otherInstructions{2,a}, trialCodes, codeStruct);
end

function rasterCells = makeRasterHeader(spikes, anlg, waves, additionalRasterFields)
% start with the spike data then add the anlg data. this means that
% stro.raster will be composed of columns:
% <neuron 1> <neuron 2> ... <neuron 3> <anlg 1> <anlg 2> ... <anlg 3> <atime> <other>

rasterCells = {};
if ~isempty(spikes.names)
    for ind = 1:length(spikes.names)
        rasterCells{ind} = spikes.names{ind};
    end
else
    ind = 0; %incase there isn't any spike data
end

% iterate over waveforms (if they exist)
for a = 1:length(waves.names)
    ind = ind + 1;
    rasterCells{ind} = waves.names{a};
end

% now add the analog channels
for a = 1:length(anlg.names)
    ind = ind + 1;
    rasterCells{ind} = anlg.names{a};
end

% explictly add the 'anlgStart' field
if ~isempty(anlg.fragStarts)
    ind = ind + 1;
    rasterCells{ind} = 'anlgStartTime';
end

% now add the additional raster headers (from the paradigm header file)
for a = 1:size(additionalRasterFields, 2)
    ind = ind + 1;
    rasterCells{ind} = additionalRasterFields{1,a};
end

function stro = getExperimentalParams(stro, ecodes, VALOFFSET, exptParamsHeader)
% deal with the entries into the stro.sum.exptParams one by one.
for a = 1:size(exptParamsHeader, 2)
    nums = dat2num(ecodes, [exptParamsHeader{4, a}], VALOFFSET, exptParamsHeader{2, a}, exptParamsHeader{3, a});
    stro.sum.exptParams.(exptParamsHeader{1,a}) = nums{1};
end

function [anlg, events] = cleanupAnalogSignals(anlg, events)
% check to make sure that the number of analog fragments is the same for
% each channel. ditto for the number of AtoD points per fragment

% Ensure an equal number of starts and stops (i.e. check to make sure
% that data acquisition was not arrested part way through a trial)
excessTrial = (length(events.start) - length(events.stop));
if excessTrial == 1
    events.start(end) = [];
    fprintf(' Unequal number of STARTs and STOPs. Deleted the terminal START.\n');
elseif excessTrial == -1 && events.stop(1) < events.start(1)
    events.stop(1) = [];
    fprintf(' Unequal number of STARTs and STOPs. Deleted the initial STOP.\n');
elseif excessTrial ~= 0
    error('  Unequal number of STARTS and STOPS.');
end

% now determine if there are the same number of trial starts/stops as
% there are anlg starts/stops. delete any anlg data that occurs before
% the first trial.
if isempty(anlg.fragStarts)
    return
end

err = abs(events.start(1) - anlg.fragStarts{1});
startIdx = find(err == min(err));
for a = 1:length(anlg.names)
    anlg.fragStarts{a}(1:(startIdx-1)) = [];
    anlg.nPtsPerFrag{a}(1:(startIdx-1)) = [];
end

if numel(anlg.names) > 1
    if ~isequal(anlg.fragStarts{1:end})
        error('frag starts are unequal');
    end
    
    if ~isequal(anlg.nPtsPerFrag{1:end})
        error('points per frag are unequal');
    end
end
