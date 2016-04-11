function makeAxonTextFile(params, allTemplates)

% params should have
%
% params.name           =>  the name of the output .atf file
% params.si             =>  the sample INTERVAL (needs to be an iteger)
% params.swpDur         =>  The total duration of the sweep IN NUMBERS OF SAMPLES!!!!
% params.trlTypes       =>  A vector of numbers, which correspond to the
%                           different trial types, if there are multiple instances of the same
%                           number, then that tType will be presented multiple times. This can have a
%                           random sequence (for tType randomization)

tt = [0:params.swpDur-1]' .* params.si;

%
% create a matrix of sweeps to be exported to the atf file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%
sweeps = nan(numel(tt), numel(params.trlTypes));
for i_swp = 1:numel(params.trlTypes)
    sweeps(:,i_swp) = allTemplates{params.trlTypes(i_swp)}(:);
end

% the atf file needs a time vector, so add that here as the first column
sweeps = cat(2, tt(:), sweeps);

%
% create the header information. Open a new file, and start writing into
% the new file.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%


header{1,:} = {'ATF', '1'};
header{2,:} = {'8', num2str(params.nSweeps+1)}; % adding one for the tvec column
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



