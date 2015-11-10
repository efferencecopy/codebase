function stro = blk2stro(nev, varargin)

    % syntax for input arguments:
    % stro = blk2stro(pathToNEV, 'NSx', pathToNSx, 'trialdef', trialDefStructure)

    % trial def could be something like: ns4::trl_start

    inargs = parseInputs(nev, varargin);
    inargs = inargs.Results;


    % import the nev
    nev = my_openNEV(nev, 'read', 'nosave', 'nomat');
    stro.sum.nev = nev.MetaTags;


    % build the 'trial' field. Include the trials start/stop times and any
    % other evnet timing that was provided during experiment
    [stro.trial, stro.sum.trialFields] = parseTrialEvents(nev, inargs.trialdef);


    % add the spike times (from the nev file).
    stro = addSpikeTimes(nev, stro);
    clear nev % unnecessary now.


    % load in the nsx data
    for i_nsx = 1:6

        nsx_type = ['ns', num2str(i_nsx)];
        fpath = inargs.(nsx_type);
        if isempty(fpath)
            continue
        end

        nsx{i_nsx} = openNSx(fpath, 'read', 'precision', 'double');

    end


    % add the raster data (from the nsx files)
    stro = addRasterData(nsx, stro);
    
    
    
    % build the 'sum' field. There could be a main exptParams field that comes 
    
    
    
    % define the indicies to the stro.trial fields
    for i_fld = 1:numel(stro.sum.trialFields)
        stro.idx.(stro.sum.trialFields{i_fld}) = i_fld;
    end

    
end % main function




function p = parseInputs(nev, optionalInputs)
    p = inputParser;
    p.addRequired('nev',      @(x) (ischar(x) && any(regexpi(x, '.nev'))))
    p.addParameter('ns1', '', @(x) (ischar(x) && any(regexpi(x, '.ns1'))))
    p.addParameter('ns2', '', @(x) (ischar(x) && any(regexpi(x, '.ns2'))))
    p.addParameter('ns3', '', @(x) (ischar(x) && any(regexpi(x, '.ns3'))))
    p.addParameter('ns4', '', @(x) (ischar(x) && any(regexpi(x, '.ns4'))))
    p.addParameter('ns5', '', @(x) (ischar(x) && any(regexpi(x, '.ns5'))))
    p.addParameter('ns6', '', @(x) (ischar(x) && any(regexpi(x, '.ns6'))))
    p.addParameter('rmch', [], @isnumeric);
    p.addParameter('trialdef', [], @isstruct)

    p.parse(nev, optionalInputs{:});
end

function [trial, trialFields] = parseTrialEvents(nev, trialdef)

    timeStamps = nev.Data.Spikes.TimeStamp;
    eventsByEtrode = nev.Data.Spikes.Electrode;
    sampRate_ev = double(nev.MetaTags.SampleRes);
    

    % find the trial starts and stops. make sure there are the same numbers
    ch_start = trialdef.tstart;
    ch_stop = trialdef.tstop;
    
    idx_start = eventsByEtrode == ch_start;
    idx_stop = eventsByEtrode == ch_stop;
    if ~(sum(idx_start) == sum(idx_stop))
        disp('I was hoping this day would never come: unequal number of trial start/stops')
        keyboard
    end
    
    % get the times of all the starts and stops. import them in 'sample
    % number' units and then convert to time (as a double)
    tt_samps_start = timeStamps(idx_start);
    tt_samps_stop = timeStamps(idx_stop);
    if any(regexpi(class(tt_samps_start), 'uint'))
        tt_samps_start = double(tt_samps_start);
        tt_samps_stop = double(tt_samps_stop);
    end
    
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
            event_tt_sec = double(event_tt_samps) ./ sampRate_ev;
            trial(i_trl, i_ev) = event_tt_sec;
            
            % update the trialFields array
            if i_trl == 1
                trialFields = cat(2, trialFields, eventnames{i_ev});
            end
            
        end
    end
end

function stro = addSpikeTimes(nev, stro)
    
    sampRate_ev = double(nev.MetaTags.SampleRes);
    timeStamps_samps = nev.Data.Spikes.TimeStamp;
    timeStamps_sec = double(timeStamps_samps) ./ sampRate_ev;
    
    Ntrials = size(stro.trial, 1);
    activeChannels = unique(nev.Data.Spikes.Electrode);
    activeChannels = activeChannels(activeChannels<=128);
    stro.sum.rasterFields = {}; % nothing here yet
    stro.ras = {};
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
    % **** need to convert things into uV ******
    % **** use the e: and t: constructs to speed things up?
    % **** if the nsx.Data array is full of cells, do things differently than
    %      if it's not full of cells
    % **** store the start time for the nsx data (on a trial by trial basis).
    %      This is important b/c the start time will depend on the sampling
    %      rate, which could differe between ns1-ns6.
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
        sampRate_nsx = double(nsx{i_ns}.MetaTags.SamplingFreq);
        tstart_nsx_sec = double(tstart_nsx_samps) ./ sampRate_nsx;
        trl_field_name = sprintf('ns%d_start_sec', i_ns);
        stro.sum.trialFields = cat(2, stro.sum.trialFields, trl_field_name);
        col = numel(stro.sum.trialFields);
        assert(size(stro.trial,2)==(col-1), 'ERROR: stro.trial has unexpected dimensions')
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
            digiRange = abs(double(nsx{i_ns}.ElectrodesInfo(i_ch).MinDigiValue)) + double(nsx{i_ns}.ElectrodesInfo(i_ch).MaxDigiValue);
            analogRange = abs(double(nsx{i_ns}.ElectrodesInfo(i_ch).MinAnalogValue)) + double(nsx{i_ns}.ElectrodesInfo(i_ch).MaxAnalogValue);
            sampsPerMicrovolt = digiRange ./ analogRange;
            trl_ras = cellfun(@(x,y) x./y, trl_ras, repmat({sampsPerMicrovolt}, size(trl_ras)), 'uniformoutput', false);
            
            % add the data to the stro.ras array
            assert(numel(trl_ras) == size(stro.trial,1), 'ERROR: trial raster data has wrong dimensions')
            stro.ras(:, col) = trl_ras(:);
            
        end
        
    end

end


