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
clear nev

% define the indicies to the stro.trial fields
for i_fld = 1:numel(stro.sum.trialFields)
    stro.idx.(stro.sum.trialFields{i_fld}) = i_fld;
end





% load in the nsx data
% **** need to convert things into uV ******
% **** use the e: and t: constructs to speed things up.
for i_nsx = 1:6
    
    nsx_type = ['ns', num2str(i_nsx)];
    fpath = inargs.(nsx_type);
    if isempty(fpath)
        continue
    end
    
    nsx{i_nsx} = openNSx(fpath, 'read', 'precision', inargs.precision);

end


% iterate over the data and build the 'ras' field.  The dimensions of 'ras'
% are <Ntrials x Nchannels>, where each channel can have a spike time
% chanel, a spike WF,  and one or more continuous channels. Deal out each
% data packet to a cell array.


% build the 'sum' field. There could be a main exptParams field that comes 



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
    p.addParameter('precision', 'double', @ischar)
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


