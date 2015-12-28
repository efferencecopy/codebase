function stro = blk2stro(varargin)

    % syntax for input arguments:
    % stro = blk2stro(pathToNEV, 'NSx', pathToNSx, 'trialdef', trialDefStructure)

    % 
    
    
    
    inargs = parseInputs(varargin);
    inargs = inargs.Results;
    assert(~isempty(inargs.trialdef), 'ERROR: trial definitions must be defined')

    
    %initialize the structure
    stro.sum = [];
    stro.trial = [];
    stro.ras = {};
    stro.other = {};
    
    
    % import the nev. Add trial events, and spike times
    nev = openNEV(inargs.nev, 'read', 'nosave', 'nomat');
    nev = forceDouble(nev);
    stro.sum.nev.MetaTags = nev.MetaTags;
    
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
    
    
    % build the 'trial' field. Include the trials start/stop times and any
    % other event timing that was provided during experiment
    [stro.trial, stro.sum.trialFields] = parseTrialEvents_byNEV(nev, inargs.trialdef);
    
    
    % add the spike times to the .ras array (from the nev file).
    stro = addSpikeTimes(nev, stro);
    clear nev % unnecessary now.
    

    % add the analog/continuous data to the .ras array (from the nsx files)
    stro = addRasterData(nsx, stro);
    
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
    p.addParameter('trialdef', [], @(x) (isempty(x) || isstruct(x)))

    p.parse(optionalInputs{:});
end

function [trial, trialFields] = parseTrialEvents_byNEV(nev, trialdef)

    % this function will fill out the stro.trial array with event times,
    % and relies on event times being defined in the NEV file (which is a
    % thresholding of the analog signals for things like 'stimOn' or
    % 'trialStart'. If the NEV file is not reliable or present, then use
    % 'parseTrialEvents_byNSX', but that will require digitizing the analog
    % experimental control signals (stimOn, trialStart etc...)
    
    timeStamps = nev.Data.Spikes.TimeStamp;
    eventsByEtrode = nev.Data.Spikes.Electrode;
    sampRate_ev = nev.MetaTags.SampleRes;
    

    % find the trial starts and stops. make sure there are the same numbers
    ch_start = trialdef.tstart;
    ch_stop = trialdef.tstop;
    
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
    eventnames = fieldnames(trialdef);
    trialFields = {'tstart', 'tstop'};
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
            
            if any(strcmpi(eventnames{i_ev}, {'tstart', 'tstop'}))
                % trial start/stop have already been imported
                continue
            end
            
            event_ch = trialdef.(eventnames{i_ev});
            event_idx = trl_eventsByEtrode == event_ch;
            event_tt_samps = trl_timeStamps(event_idx);
            
            if numel(event_tt_samps) > 1
                if strcmpi(eventnames{i_ev}, 'stimon')
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
    activeChannels = activeChannels(activeChannels<=128);
    stro.sum.rasterFields = {}; % nothing here yet
    
    % make sure that stro.ras has the correct number of rows by
    % preallocating the first column. This is important b/c if there are no
    % spike from any channel, then a row of stro.ras could be left out.
    stro.ras = repmat({[]}, Ntrials, numel(activeChannels));
    
    % pull out the spike times
    for i_trl = 1:Ntrials
        
        trl_start_sec = stro.trial(i_trl, strcmpi(stro.sum.trialFields, 'tstart'));
        trl_stop_sec = stro.trial(i_trl, strcmpi(stro.sum.trialFields, 'tstop'));
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
        % start time could differe b/w nsx versions, but I'd be disappointed if
        % that were true. tbd..
        tstart_nsx_samps = nsx{i_ns}.MetaTags.Timestamp;
        sampRate_nsx = nsx{i_ns}.MetaTags.SamplingFreq;
        tstart_nsx_sec = tstart_nsx_samps ./ sampRate_nsx;
        trl_field_name = sprintf('ns%d_start_sec', i_ns);
        stro.sum.trialFields = cat(2, stro.sum.trialFields, trl_field_name);
        col = numel(stro.sum.trialFields);
        assert(size(stro.trial,2)==(col-1), 'ERROR: stro.trial has unexpected dimensions')
        keyboard
        assert(size(stro.trial,1)==(numel(tstart_nsx_sec)), 'ERROR: wrong number of tstarts for nsx')
        stro.trial(:,col) = tstart_nsx_sec(:);


        % add the data
        Nchannels = numel(nsx{i_ns}.ElectrodesInfo);
        for i_ch = 1:Nchannels
            
            % make a rasterFields entry. Look to see if it already exists
            % (it shouldn't).
            etrode_name = nsx{i_ns}.ElectrodesInfo(i_ch).Label;
            etrode_name = sprintf('%s_ns%d', deblank(etrode_name), i_ns);
            idx = strcmpi(stro.sum.rasterFields, etrode_name);
            assert(~any(idx), 'ERROR: nsx raster entry alreay exists');
            stro.sum.rasterFields = cat(2, stro.sum.rasterFields, etrode_name);
            col = numel(stro.sum.rasterFields);
            
            % allocate the data for all trials for this channel.
            trl_ras = cellfun(@(x,y) x(y,:), nsx{i_ns}.Data, repmat({i_ch}, size(nsx{i_ns}.Data)), 'uniformoutput', false);
            
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


