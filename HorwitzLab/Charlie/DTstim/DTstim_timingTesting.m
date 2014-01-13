%% Open a DTstim file

fin
DS = dsobj;

%% pulse widths as assessed by the ecodes

widths = [];
for a = 1:size(DS.trial,1)
    if numel(DS.other{a,2}) == numel(DS.other{a,1})
        widths = [widths; (DS.other{a,2}-DS.other{a,1})];
    else
        disp('bad trial')
    end
end

figure
hist(widths, 40)
sum(widths<0.0015) + sum(widths>0.0025)
length(widths)

%% pulse widths as assessed by the continuous records.

thresh = 1; %in volts?
sampFreq = DS.sum.analog.storeRates{1};

widths = [];
laserIdx = strcmp(DS.sum.rasterCells, 'AD13');
for a = 1:size(DS.trial,1);
    tmp = DS.ras{a,laserIdx};
    t = (0:(length(tmp)-1)).*(1./sampFreq); %in seconds
    t = t + DS.ras{a, DS.idx.anlgStart}; %compensate for the offset
    idxUp = false(1,length(t));
    idxDown = false(1,length(t));
    idxUp(2:end) = diff(tmp>thresh) == 1; %notice I'm indexing from 2 on b/c of the diff command.
    idxDown(2:end) = diff(tmp>thresh) == -1;
    
    %assign the times
    tStart = t(idxUp);
    tEnd = t(idxDown);
    
    %calculate the pulse widths
    widths = [widths; (tEnd-tStart)'];
    
    if (1) % looking to see if the ecodes and anlg sig are identical
        pulseTimesCodes = [DS.other{a,1} ; DS.other{a,2}];
        pulseTimesAnlg = [tStart(:); tEnd(:)];
        figure, hold on,
        plot(t, tmp, 'b');
        plot(pulseTimesCodes, ones(length(pulseTimesCodes),1)*(max(tmp)*.5), 'm*')
        plot(pulseTimesAnlg, ones(length(pulseTimesAnlg),1)*(max(tmp)*.5), 'g*')
        ylim([min(tmp)*.98 max(tmp)*1.05])
        keyboard
    end
end

figure
hist(widths, 40)
sum(widths<0.0015) + sum(widths>0.0025)
length(widths)

%% code to open .nex files which lack headers
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

%read in the header info
version = fread(fid, 1, 'int32');
comment = deblank(char(fread(fid, 256, 'char')'));
freq = fread(fid, 1, 'double');                     %acquisition rate for things other than continuous variables
tbeg = fread(fid, 1, 'int32') ./ freq;
tend = fread(fid, 1, 'int32') ./ freq;
nvar = fread(fid, 1, 'int32');
fseek(fid, 260, 'cof'); % skip location of next header and padding


%now read the file putting the appropriate 'type' of Nex info in either
%the ecodes, spikes, or analog vectors. For now ignore the other 'types'
anlg = [];
contSigCount = 0;
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
    
    
    switch name(1:2)
        case 'AD',  % continuous variable (i.e. eye signals)
            disp('here')
            contSigCount = contSigCount+1;
            anlg.names(contSigCount) = {name};
            anlg.ADtoMV(contSigCount) = {ADtoMV};
            anlg.ADFrequency{contSigCount} = WFrequency;
            fseek(fid, offset, 'bof');
            
            % get the start times of each analog signal snippet
            anlg.fragStarts{contSigCount} = fread(fid, [n 1], 'int32')./freq;
            
            %get the number of points sammpled during each of the snippets
            nFrag = fread(fid, [n 1], 'int32');
            nFrag(n+1) = NPointsWave;
            anlg.nPtsPerFrag{contSigCount} = nFrag;
            %now bring in the AtoD data for the entire recording
            anlg.data{contSigCount} = fread(fid, [NPointsWave 1], 'int16').*ADtoMV + MVOfffset;
            
            fseek(fid, filePosition, 'bof');
            
        otherwise
            fprintf('unknown variable type <%s> \n', name);
    end
    dummy = fread(fid, 60, 'char'); %skip some junk to get to the next block of data
end
fclose(fid);

%% SOME CODE TO ASSESS LASER PULSE FIDELITY
fin
DS = dsobj('test082911008.nex');  %test082911001 => 5mW
diodeIdx = strcmp(DS.sum.rasterCells, 'AD13');
anlgRate = DS.sum.analog.storeRates{diodeIdx};
OFFSET = 0.0005;

pulseWidths = DS.trial(:, DS.idx.laserWidth);
intPulseWidths = DS.trial(:, DS.idx.laserIPI);
uniquePulseTypes = unique([pulseWidths, intPulseWidths], 'rows');

%iterate over the data to find the cases with common pulse and int-pulse
%intervals
dat.wfs = {};
dat.pulseTypes = [];
for a = 1:size(uniquePulseTypes,1)
   pulseType = uniquePulseTypes(a,:);
   tList = ismember([pulseWidths, intPulseWidths], pulseType, 'rows');
   tList = find(tList);
   
   dat.pulseTypes(a,:) = pulseType;
   dat.wfs{a} = [];
   dat.avg{a} = [];
   for i = 1:length(tList)
       t = [0:(length(DS.ras{tList(i),diodeIdx})-1)].*(1/anlgRate);
       t = t+DS.ras{tList(i),DS.idx.anlgStart};
       starts = DS.other{tList(i),1};
       stops = DS.other{tList(i),2};
       sig = DS.ras{tList(i),diodeIdx};
       baseline = mean(sig(1:100))
       sig = sig-baseline;
       for p = 1:(length(starts)-1)
           %the entire cycle
           cycleidx = (t>=(starts(p)-OFFSET))&(t<starts(p+1));
           try
           dat.wfs{a}(end+1,:) = sig(cycleidx);
           catch
           end
           
           %just the pulse period
           pulseidx = (t>=(starts(p)))&(t<stops(p));
           pulse = sig(pulseidx);
           dat.avg{a}(end+1) = mean(pulse);
       end
   end
   dat.timeVec{a} = ([0:(size(dat.wfs{a},2)-1)].*(1/anlgRate))-OFFSET;
   
   figure
   subplot(1,2,1)
   plot(dat.timeVec{a}*1000, mean(dat.wfs{a},1))
   xlabel('time (ms)')
   ylabel('Diode Output')
   title(sprintf('Pulse Type: %s', num2str(dat.pulseTypes(a,:))))
   subplot(1,2,2)
   hist(dat.avg{a})
end

figure, hold on,
bar(1:length(dat.avg), cellfun(@mean, dat.avg),'facecolor', 'w')
errorbar(1:length(dat.avg), cellfun(@mean, dat.avg), cellfun(@std, dat.avg), 'k.', 'markersize', 0.1)
for a = 1:size(dat.pulseTypes,1)
    text(a, .1, num2str(dat.pulseTypes(a,1)))
    text(a, .05, num2str(dat.pulseTypes(a,2)))
end
xlabel('pulse type')
ylabel('diode output')
set(gca, 'xticklabel', [])












