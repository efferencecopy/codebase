function out = abfcat(catdim, varargin)

% concatenate ABF files.
%
%  out = abfcat(catdim, axonFiles)
%
% 'axonFiles' should be a cell array of file name (strings) and 'catdim' is
% the dimension along which the data should be concatenated. In episodic
% mode, the 'data' and 'wf' fields are typically 3 dimensional, so 'catdim'
% will be 3 when tall sweeps have the identical command waveform, or 4 when
% sweeps have different command waveforms.


% parse the arguments
if numel(varargin) == 1 && iscell(varargin)
    varargin = varargin{1};
end

for a = 1:numel(varargin)
    ax{a} = abfobj(varargin{a});
end

% use the first file as the standard
out = ax{1};

% check some stuff to make sure that the file in question is idential
% to the standard
for a = 2:numel(ax)
    
    % checking a few (easy) things:
    protocolMatch = strcmpi(out.head.protocolName, ax{a}.head.protocolName); %protocol names should match
    chUnitsMatch = all(strcmpi(out.head.recChUnits, ax{a}.head.recChUnits)); % recorded channel units should match
    chNamesMatch =  all(strcmpi(out.head.recChNames, ax{a}.head.recChNames)); % recorded channel units should match
    %waveformMatch = all(out.wf(:) == ax{a}.wf(:));% the waveforms should be identical
    
    if all([protocolMatch, chUnitsMatch, chNamesMatch])
        
        % find Vclamp channels
        rec_pA = strcmpi(ax{a}.head.recChUnits, 'pa');
        sec_ch = cellfun(@(x) ~isempty(x), regexpi(ax{a}.head.recChNames, 'sec'));
        vclamp_channel = rec_pA & ~sec_ch; % voltage clamp records current (pA) on the primary channel
        
        
        %
        % Checks for voltage clamp experiemnts
        %
        % the file in question passed the initial test, but now consider
        % holding potential
        %
        if any(vclamp_channel)
            rec_mV = strcmpi(ax{a}.head.recChUnits, 'mv');
            vclamp_monitor = rec_mV & sec_ch; % vclamp monitor records voltage (mV) on the secondary channel
            if any(vclamp_monitor)
                monitor_ch = find(vclamp_monitor);
                for i = 1:numel(monitor_ch)
                    idx = monitor_ch(i);
                    delta = mean(out.dat(1:200, idx, 1)) - mean(ax{a}.dat(1:200, idx, 1));
                    holdingMatch(i) = abs(delta) < 2; % allow the Vcmd to be off by 2 mv.
                end
            else
                warning('Holding potentials might not be identical')
                holdingMatch = true;
            end
            
            if all(holdingMatch)
                vclampcheck = true;
            else
                vclampcheck = false;
            end
        else
            vclampcheck = true;
        end
        
        
        %
        % Checks for current clamp experiemnts
        %
        % the file in question passed the initial test, but now consider
        % holding potential
        %
        
        %%% nothing yet %%%
        iclampcheck = true;
        
        
        % if everything checks out, smash the two files together.
        if vclampcheck && iclampcheck
            out.dat = cat(catdim, out.dat, ax{a}.dat);
            out.wf = cat(catdim, out.wf, ax{a}.wf);
        end
    end
    
    if ~all([protocolMatch, chUnitsMatch, chNamesMatch, iclampcheck, vclampcheck])
        fprintf('File <%s> was not included \n', varargin{a})
        keyboard
    end
end
