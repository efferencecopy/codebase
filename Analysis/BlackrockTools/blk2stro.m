function stro = blk2stro(varargin)

    % syntax for input arguments:
    % stro = blk2stro(pathToNEV, 'NSx', pathToNSx, 'exdef', exptDefStructure)
    
    
    inargs = parseInputs(varargin);
    inargs = inargs.Results;
    assert(~isempty(inargs.exdef), 'ERROR: trial definitions must be defined')
    
    
    %initialize the structure
    stro.sum = [];
    stro.trial = [];
    stro.ras = {};
    stro.other = {};
    
    
    % import the nev.
    nev = openNEV(inargs.nev, 'read', 'nosave', 'nomat');
    nev = forceDouble(nev);
    stro.sum.nev.MetaTags = nev.MetaTags;
    
    % sabotage all current attempts to use the waveforms b/c the units are
    % wrong in the NEV file. there is a different scale factor for chs
    % 1:128 than for the anlg_in channels. find the scale factor in:
    %
    % nev.ElectrodeInfo(ch).DigitalFactor (units in nanoVolts / digi step?)
    % 
    % and see line 746 of openNEV for conversion to 'uV'
    % elecDigiFactors = double(1000./[NEV.ElectrodesInfo(NEV.Data.Spikes.Electrode).DigitalFactor]);
    %
    nev.Data.Spikes.Waveform = [];
    
    % load in the nsx data
    for i_nsx = 1:6

        nsx_type = ['ns', num2str(i_nsx)];
        fpath = inargs.(nsx_type);
        if isempty(fpath)
            continue
        end

        nsx{i_nsx} = openNSx(fpath, 'read', 'precision', 'double');
        nsx{i_nsx} = forceDouble(nsx{i_nsx});
    end
    
    % build the 'sum' field, one for each NSx version
    for i_ns = 1:numel(nsx)
        if isempty(nsx{i_ns})
            continue
        else
            name = sprintf('ns%d', i_ns');
            stro.sum.(name).MetaTags = nsx{i_ns}.MetaTags;
            stro.sum.(name).ElectrodesInfo = nsx{i_ns}.ElectrodesInfo;
        end
    end
    
    
    % build the .Trial and .Ras fields
    switch inargs.exdef.params.method
        case 'nev'
            
            % build the 'trial' field.
            [stro.trial, stro.sum.trialFields] = parseTrialEvents_byNEV(nev, inargs.exdef);
            
            % add the spike times to the .ras array (from the nev file).
            stro = addSpikeTimes(nev, stro);
            
            % add the analog/continuous data to the .ras array (from the nsx files)
            stro = addRasterData(nsx, stro);
            
        case {'ns1', 'ns2', 'ns3', 'ns4', 'ns5', 'ns6'}
            
            % make sure that the designated NSx type exists.
            nsx_num = str2double(inargs.exdef.params.method(end));
            assert(~isempty(nsx{nsx_num}), 'ERROR: could not find the designated NSx data');
            
            % build the 'ras' field with continuous data
            stro = addRasterData(nsx, stro);
            
            % build the .trial field using a thresholding of the analog
            % data
            stro = parseTrialEvents_byNSX(stro, inargs.exdef);
            
            % add the spike times to the .ras array (from the nev file).
            stro = addSpikeTimes(nev, stro);
             
        otherwise
            error('unknown trial definition method')
    end
    
    
    % delete trials that are corrupt. Must have tstart and tstop, and
    % tstart must come before tstop
    idx_start = strcmpi(stro.sum.trialFields, 'trl_start');
    idx_stop =strcmpi(stro.sum.trialFields, 'trl_stop');
    idx = idx_start | idx_stop;
    assert(sum(idx)==2);
    l_nan = isnan(sum(stro.trial(:, idx), 2));
    if any(l_nan)
        stro.trial(l_nan, :) = [];
        stro.ras(l_nan,:) = [];
        assert(size(stro.ras,1) == size(stro.trial,1));
        fprintf('Deleted %d trial(s): undefined start/stop times\n', sum(l_nan));
    end
    l_outoforder = stro.trial(:, idx_start) > stro.trial(:, idx_stop);
    if any(l_outoforder)
        stro.trial(l_outoforder, :) = [];
        stro.ras(l_outoforder,:) = [];
        assert(size(stro.ras,1) == size(stro.trial,1));
        fprintf('Deleted %d trial(s): start time occurs after stop time\n', sum(l_outoforder));
    end
    
    
     % define the indicies to the stro.trial fields
    if isfield(stro.sum, 'trialFields')
        for i_fld = 1:numel(stro.sum.trialFields)
            stro.sum.idx.(stro.sum.trialFields{i_fld}) = i_fld;
        end
    end
    if isfield(stro.sum, 'rasterFields')
        for i_fld = 1:numel(stro.sum.rasterFields)
            stro.sum.idx.(stro.sum.rasterFields{i_fld}) = i_fld;
        end
    end

    
end % main function




function p = parseInputs(optionalInputs)
    p = inputParser;
    p.addParameter('nev', '', @(x) (isempty(x) || (ischar(x) && any(regexpi(x, '.nev')))))
    p.addParameter('ns1', '', @(x) (isempty(x) || (ischar(x) && any(regexpi(x, '.ns1')))))
    p.addParameter('ns2', '', @(x) (isempty(x) || (ischar(x) && any(regexpi(x, '.ns2')))))
    p.addParameter('ns3', '', @(x) (isempty(x) || (ischar(x) && any(regexpi(x, '.ns3')))))
    p.addParameter('ns4', '', @(x) (isempty(x) || (ischar(x) && any(regexpi(x, '.ns4')))))
    p.addParameter('ns5', '', @(x) (isempty(x) || (ischar(x) && any(regexpi(x, '.ns5')))))
    p.addParameter('ns6', '', @(x) (isempty(x) || (ischar(x) && any(regexpi(x, '.ns6')))))
    p.addParameter('rmch', [], @(x) (isempty(x) || isnumeric(x)));
    p.addParameter('exdef', [], @(x) (isempty(x) || isstruct(x)))

    p.parse(optionalInputs{:});
end

function stro = parseTrialEvents_byNSX(stro, exdef)
    
    nsxversion = exdef.params.method;
    trialParams = fieldnames(exdef.tcode);
    
    % iterate over tcode types and determine the event timing. Do this
    % simultaneously for all trials.
    for i_param = 1:numel(trialParams)
        
        % find the analog data that corresponds to the trigger for this
        % trial event
        rasidx = strcmpi(stro.sum.rasterFields, [trialParams{i_param}, '_', nsxversion]);
        assert(sum(rasidx)==1)
        assert(stro.sum.rasterChanIDs(rasidx) == exdef.tcode.(trialParams{i_param})); % cross check index with the exdef field
        
        % pull out the data for this trigger
        trig_anlg = stro.ras(:, rasidx);
        
        % find the threshold crossings. Assume a noiseless signal and make
        % a low threshold. Only store the first threshold crossing.
        threshold = 4000; % in uV.
        trig_crossings = cellfun(@(x) [nan, diff(x>=threshold)]==1, trig_anlg, 'uniformoutput', false);
        trig_crossings = cellfun(@(x) find(x==1, 1, 'first'), trig_crossings, 'uniformoutput', false);
        
        % this line finds instances where there were no threshold crossings
        % (empty vectors) and turns them into nans
        trig_crossings = cellfun(@(x) min([x, nan]), trig_crossings);
        
        % convert to time in seconds
        sampFreq_nsx = stro.sum.(nsxversion).MetaTags.SamplingFreq;
        trig_crossings = trig_crossings ./ sampFreq_nsx;
        
        % add the time from beginning of file
        tfield = [nsxversion, '_start_sec'];
        idx = strcmpi(stro.sum.trialFields, tfield);
        assert(sum(idx)==1);
        offset = stro.trial(:, idx) - (1/sampFreq_nsx); % adding the full offset will be wrong by one sample, so subtract off a sample.
        trig_crossings = trig_crossings + offset;
        
        % add the data to the stro.trial array
        col = size(stro.trial, 2) + 1;
        stro.trial(:, col) = trig_crossings;
        
        % update the stro.sum.trialFields
        stro.sum.trialFields = cat(2, stro.sum.trialFields, trialParams{i_param});
        
        
    end
    
end

function [trial, trialFields] = parseTrialEvents_byNEV(nev, exdef)

    % this function will fill out the stro.trial array with event times,
    % and relies on event times being defined in the NEV file (which is a
    % thresholding of the analog signals for things like 'stimOn' or
    % 'trialStart'. If the NEV file is not reliable or present, then use
    % 'parseTrialEvents_byNSX', but that will require digitizing the analog
    % experimental control signals (stimOn, trialStart etc...)
    
    error('modify this code to use the serial_digital_IO codes instead of the neural ecodes')
    
    % subtracting [(presamps+1) * (1/sampFreq)] should identify the onset
    % of event trigers properly.
    error('fix event times in NEV file b/c nev thresh crossings do not subtract pre-crossing samples')
    
    timeStamps = nev.Data.Spikes.TimeStamp;
    eventsByEtrode = nev.Data.Spikes.Electrode;
    sampRate_ev = nev.MetaTags.SampleRes;
    

    % find the trial starts and stops. make sure there are the same numbers
    ch_start = exdef.tcode.trl_start;
    ch_stop = exdef.tcode.trl_stop;
    
    idx_start = eventsByEtrode == ch_start;
    idx_stop = eventsByEtrode == ch_stop;
    if ~(sum(idx_start) == sum(idx_stop))
        error('I was hoping this day would never come: unequal number of trial start/stops')
    end
    
    % get the times of all the starts and stops.
    tt_samps_start = timeStamps(idx_start);
    tt_samps_stop = timeStamps(idx_stop);
    assert(isfloat(tt_samps_start), 'ERROR: analysis was expecing floats') % this should have been taken care of by 'forceDouble'
    
    tt_sec_start = tt_samps_start ./ sampRate_ev;
    tt_sec_stop = tt_samps_stop ./ sampRate_ev;
    assert(tt_sec_start(1) < tt_sec_stop(1), 'ERROR: starts and stops are out of order')
    
    % define the 'trial' field of the STRO structure
    eventnames = fieldnames(exdef.tcode);
    trialFields = {'trl_start', 'trl_stop'};
    trial = nan(numel(tt_sec_start), numel(eventnames));
    trial(:,1) = tt_sec_start(:);
    trial(:,2) = tt_sec_stop(:);
    
    % now import the remaining trial events and parameters. If there is a
    % 'time' event
    for i_trl = 1:size(trial,1)
        
        % pull out the codes and timestams for each trial
        tstart_samps = tt_samps_start(i_trl);
        tstop_samps = tt_samps_stop(i_trl);
        
        trl_idx = (timeStamps >= tstart_samps) & (timeStamps < tstop_samps);
        trl_eventsByEtrode = eventsByEtrode(trl_idx);
        trl_timeStamps = timeStamps(trl_idx);
                
        for i_ev = 1:numel(eventnames)
            
            if any(strcmpi(eventnames{i_ev}, {'trl_start', 'trl_stop'}))
                % trial start/stop have already been imported
                continue
            end
            
            event_ch = exdef.tcode.(eventnames{i_ev});
            event_idx = trl_eventsByEtrode == event_ch;
            event_tt_samps = trl_timeStamps(event_idx);
            
            if numel(event_tt_samps) > 1
                if strcmpi(eventnames{i_ev}, 'stim_on')
                    event_tt_samps = event_tt_samps(1);
                else
                    error('Expected 1 event but found %d', numel(event_tt_samps))
                end
            end
                    
            % convert to seconds, and add to the trial matrix
            event_tt_sec = event_tt_samps ./ sampRate_ev;
            trial(i_trl, i_ev) = event_tt_sec;
            
            % update the trialFields array
            if i_trl == 1
                trialFields = cat(2, trialFields, eventnames{i_ev});
            end
            
        end
    end
end

function stro = addSpikeTimes(nev, stro)
    
    sampRate_ev = nev.MetaTags.SampleRes;
    timeStamps_samps = nev.Data.Spikes.TimeStamp;
    timeStamps_sec = timeStamps_samps ./ sampRate_ev;
    assert(isfloat(timeStamps_samps), 'ERROR: analysis was expecing floats') % this should have been taken care of by 'forceDouble'
    
    
    Ntrials = size(stro.trial, 1);
    activeChannels = unique(nev.Data.Spikes.Electrode);
    activeChannels = activeChannels(activeChannels<=128); % the neural channels
    if ~isfield(stro.sum, 'rasterFields')
        stro.sum.rasterFields = {}; % nothing here yet
        stro.sum.rasterChanIDs = [];
    end
    
    % stro.ras may or may not already exist. If it doesn't exist, make sure
    % that stro.ras has the correct number of rows by preallocating the
    % first column. This is important b/c if there are no spikes from any
    % channel, then a row of stro.ras could be left out.
    if ~isfield(stro, 'ras') || isempty(stro.ras)
        stro.ras = repmat({[]}, Ntrials, numel(activeChannels));
    end
    
    % pull out the spike times
    for i_trl = 1:Ntrials
        
        trl_start_sec = stro.trial(i_trl, strcmpi(stro.sum.trialFields, 'trl_start'));
        trl_stop_sec = stro.trial(i_trl, strcmpi(stro.sum.trialFields, 'trl_stop'));
        trl_idx = (timeStamps_sec >= trl_start_sec) & (timeStamps_sec < trl_stop_sec);
        trl_etrode = nev.Data.Spikes.Electrode(trl_idx);
        trl_unit = nev.Data.Spikes.Unit(trl_idx);
        trl_timeStamps_sec = timeStamps_sec(trl_idx);
        
        
        for i_ch = 1:numel(activeChannels)
           
            ch = activeChannels(i_ch);
            ch_spikes_idx = trl_etrode == ch;
            ch_units = trl_unit(ch_spikes_idx);
            ch_timeStamps_sec = trl_timeStamps_sec(ch_spikes_idx);

            uniqueUnits = unique(ch_units); % 'zero' unit is MU, ~=0 is SU
            
            for i_unt = 1:numel(uniqueUnits)
                
                unit_idx = ch_units ==  uniqueUnits(i_unt);
                unit_timeStamps_sec = ch_timeStamps_sec(unit_idx);
                
                % look for an existing stro.ras column. If there isn't a
                % label, make one and add the data. SU = single unit, MU =
                % multiunit
                fieldName = nev.ElectrodesInfo(ch).ElectrodeLabel';
                if uniqueUnits(i_unt) == 0
                    fieldName = sprintf('%s_mu', deblank(fieldName)); % MU
                else
                    fieldName = sprintf('%s_su%d', deblank(fieldName), uniqueUnits(i_unt)); % SU
                end
                fieldIdx = strcmpi(fieldName, stro.sum.rasterFields);
                
                
                if sum(fieldIdx) == 0 % the unit has not been added to the array for any trial
                    
                    stro.sum.rasterFields = cat(2, stro.sum.rasterFields, fieldName);
                    stro.sum.rasterChanIDs = cat(2, stro.sum.rasterChanIDs, activeChannels(i_ch));
                    col = numel(stro.sum.rasterFields);
                    
                else
                    assert(sum(fieldIdx) == 1, 'ERROR: duplicate fields in stro.ras');
                    col = fieldIdx;
                end
                
                
                % add the data
                stro.ras{i_trl, col} = unit_timeStamps_sec;
                
            end
        end
    end
    
end

function stro = addRasterData(nsx, stro)

    %initalize the stro.sum.rasterFields header if it doesn't exist
    if ~isfield(stro.sum, 'rasterFields')
        stro.sum.rasterFields = {}; % nothing here yet
        stro.sum.rasterChanIDs = [];
    end

    % iterate over the data and build the 'ras' field.  The dimensions of 'ras'
    % are <Ntrials x Nchannels>, where each channel can have a spike time
    % chanel, a spike WF,  and one or more continuous channels. Deal out each
    % data packet to a cell array.
    for i_ns = 1:numel(nsx);

        % skip nsx version that do not exist
        if isempty(nsx{i_ns})
            continue
        end
        
        
        % currently the code can only deal with nsx files that have
        % already been separated by trials (i.e. external start stop
        % signal for data aquisition).
        assert(iscell(nsx{i_ns}.Data) && iscell(nsx{i_ns}.Data(1)), 'ERROR: nsx data needs to be in cell format to use this module.')
        
        
        % add the nsx start time to the stro.trial array. This is the time at
        % which the first sample of analog data was aquired. Theoretically, the
        % start time could differ b/w nsx versions, but I'd be disappointed if
        % that were true. tbd..
        tstart_nsx_samps = nsx{i_ns}.MetaTags.Timestamp;
        sampRate_nsx = nsx{i_ns}.MetaTags.SamplingFreq;
        tstart_nsx_sec = tstart_nsx_samps ./ sampRate_nsx;
        trl_field_name = sprintf('ns%d_start_sec', i_ns);
        
        % stro.sum.trialFields might not exist yet if the caller specifies
        % 'nsx' as the method to parse .trial field data...
        if ~isfield(stro.sum, 'trialFields')
            stro.sum.trialFields = {};
        end
        stro.sum.trialFields = cat(2, stro.sum.trialFields, trl_field_name);
        
        % check the dimensions of stro.trial and make sure it's the same as
        % should be expected based off the continuous data
        col = numel(stro.sum.trialFields);
        if ~isempty(stro.trial)
            assert(size(stro.trial,2)==(col-1), 'ERROR: stro.trial has unexpected dimensions')
            assert(size(stro.trial,1)==(numel(tstart_nsx_sec)), 'ERROR: wrong number of tstarts for nsx')
        end
        
        % add the start times to the stro.trial array
        stro.trial(:,col) = tstart_nsx_sec(:);

        % add the continuous data
        Nchannels = numel(nsx{i_ns}.ElectrodesInfo);
        for i_ch = 1:Nchannels
            
            % make a rasterFields entry. Look to see if it already exists
            % (it shouldn't). Also error check to make certain that nerual
            % channels have names that correspond precicely to the actual
            % channel number.
            etrode_ID = nsx{i_ns}.ElectrodesInfo(i_ch).ElectrodeID;
            etrode_name = deblank(nsx{i_ns}.ElectrodesInfo(i_ch).Label);
            if etrode_ID <= 128 % neural channels
                assert(regexpi(etrode_name, 'chan')==1, 'ERROR: etrode name must be generic "chan"');
                etrode_num = str2double(etrode_name(5:end));
                assert(etrode_num == etrode_ID, 'ERROR: mismatch b/w etrode name and ID');
            end
                
            etrode_name = sprintf('%s_ns%d', etrode_name, i_ns);
            idx = strcmpi(stro.sum.rasterFields, etrode_name);
            assert(~any(idx), 'ERROR: nsx raster entry alreay exists');
            stro.sum.rasterFields = cat(2, stro.sum.rasterFields, etrode_name);
            stro.sum.rasterChanIDs = cat(2, stro.sum.rasterChanIDs, etrode_ID);
            col = numel(stro.sum.rasterFields);
            
            % allocate the data for all trials for this channel.
            trl_ras = cellfun(@(x) x(i_ch,:), nsx{i_ns}.Data, 'uniformoutput', false);
            
            % the data class should be a double by now, but just check
            dataclass = unique(cellfun(@class, trl_ras, 'uniformoutput', false)');
            assert(numel(dataclass)==1 && strcmpi(dataclass, 'double'), 'ERROR: the data are not doubles!!');
            
            % convert the raw data into micro volts
            digiRange = abs(nsx{i_ns}.ElectrodesInfo(i_ch).MinDigiValue) + nsx{i_ns}.ElectrodesInfo(i_ch).MaxDigiValue;
            analogRange = abs(nsx{i_ns}.ElectrodesInfo(i_ch).MinAnalogValue) + nsx{i_ns}.ElectrodesInfo(i_ch).MaxAnalogValue;
            sampsPerMicrovolt = digiRange ./ analogRange;
            trl_ras = cellfun(@(x,y) x./y, trl_ras, repmat({sampsPerMicrovolt}, size(trl_ras)), 'uniformoutput', false);
            
            % add the data to the stro.ras array
            assert(numel(trl_ras) == size(stro.trial,1), 'ERROR: trial raster data has wrong dimensions')
            stro.ras(:, col) = trl_ras(:);
        end
        
    end

end

function in = forceDouble(in)

    if isstruct(in)
        for i_struct = 1:numel(in);
            fldnames = fieldnames(in(i_struct));
            for i_fld = 1:numel(fldnames)
                in(i_struct).(fldnames{i_fld}) = forceDouble(in(i_struct).(fldnames{i_fld})); % recursive part
            end
        end
    elseif iscell(in)
        for i_idx = 1:numel(in)
            if isnumeric(in{i_idx})
                in{i_idx} = double(in{i_idx});
            end
        end
    elseif ismatrix(in)
        if isnumeric(in)
            in = double(in);
        end
    end


end


