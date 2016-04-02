function makeAxonTextFile(params)

% params should have
%
% params.name           =>  the name of the output .atf file
% params.type           =>  'train', 'pulse'
% params.si             =>  the sample INTERVAL (needs to be an iteger)
% params.swpDur         =>  The total duration of the sweep IN NUMBERS OF SAMPLES!!!!
% params.tStart         =>  the time of the first pulse
% params.pAmp           =>  A vector of amplitudes for the pulse height [interleaved variable]
% params.pFreq          =>  A vector of frequencies for the pulse train [interleaved variable]
% params.nPulses        =>  The number of pulses in a pulse train
% params.pWidth         =>  A vector of pulse widths (in seconds)  [interleaved variable]
% params.recovery       =>  true/false, should there be a recovery pulse?
% params.recoveryTime   =>  A vector of numbers corresponding to the recovery time in seconds [interleaved variable]
% params.nReps          =>  Number of repeates each stimulus should be presented


tt = [0:params.swpDur-1]' .* params.si;


%
% create a matrix of sweeps to be exported to the atf file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(params, 'trlTypes')
    trlTypes = params.trlTypes;
else
    trlTypes = 1:size(conditions,1);
    trlTypes = repmat(trlTypes, 1, params.nReps);
    randIdx = randperm(numel(trlTypes)); % randomize the order
    trlTypes = trlTypes(randIdx);
end

sweeps = nan(numel(tt), numel(trlTypes));
for i_swp = 1:numel(trlTypes)
    sweeps(:,i_swp) = params.templates{trlTypes(i_swp)}(:);
end

% the atf file needs a time vector, so add that here as the first column
sweeps = cat(2, tt(:), sweeps);



%
% create the header information. Open a new file, and start writing into
% the new file.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%


header{1,:} = {'ATF', '1'};
header{2,:} = {'8', num2str(params.nSweeps)};
header{3} = {'AcquisitionMode=Episodic Stimulation'};
header{4} = {'Comment='};
header{5} = {'YTop=10'}; % theoretically could be different, but I think this is safe for now
header{6} = {'YBottom=-10'};
header{7} = {'SyncTimeUnits=0.4'};
header{8} = {'SweepStartTimesMS='}; % this gets filled in later...
header{9} = {'SignalsExported=LED_cmd'};

tmp = repmat({'LED_cmd'}, 1, params.nSweeps);
header{10,:} = {'Signals=', tmp{:}};


tmp = cellfun(@(x,y) sprintf(x,y),...
       repmat({'Trace #%d (V)'}, 1, params.nSweeps),...
       mat2cell([1:params.nSweeps]', ones(params.nSweeps,1))', 'uniformoutput', false);
header{11,:} = {'Time (s)', tmp{:}};


%
% open a new file. Add the waveforms line by line (actually time point by
% time point for all sweeps simultaneously)
%
fileID = fopen(params.name,'w');


% iterate over the header, adding line by line. All the entries are
% strings, but I need to append \t and \n characters appropriately
specMotif = '%s \t ';
for row = 1:size(header,1);
    nCols = size(header{row},2);
    formatSpec = repmat(specMotif, 1, nCols-1);
    formatSpec = [formatSpec, '%s \n'];
    fprintf(fileID, formatSpec, header{row}{1:end});
end



% iterate over the sweeps, adding line by line. All the entries are
% doubles, but I need to append \t and \n characters appropriately
specMotif = '%.12f \t ';
nCols = size(sweeps,2);
formatSpec = repmat(specMotif, 1, nCols-1);
formatSpec = [formatSpec, '%.12f \n'];
for row = 1:size(sweeps,1);
    fprintf(fileID, formatSpec, sweeps(row,:));
end

fclose(fileID);



