function out = abfcat(varargin)


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
    waveformMatch = all(out.wf(:) == ax{a}.wf(:));% the waveforms should be identical
    
    if all([protocolMatch, chUnitsMatch, chNamesMatch, waveformMatch])
        
        % the file in question passed the initial test, but no consider
        % holding potential:
        vclamp = strcmpi(ax{a}.head.recChUnits, 'pa');
        if any(vclamp)
            vclamp_monitor = find(strcmpi(ax{a}.head.recChUnits, 'mv'));
            if any(vclamp_monitor)
                for i = vclamp_monitor
                    delta = mean(out.dat(1:200, i, 1)) - mean(ax{a}.dat(1:200, i, 1));
                    holdingMatch(i) = abs(delta) < 2; % allow the Vcmd to be off by 2 mv.
                end
            else
                warning('Holding potentials might not be identical')
                holdingMatch = true;
            end
        end
        
        if all(holdingMatch)
            out.dat = cat(3, out.dat, ax{a}.dat);
            out.wf = cat(3, out.wf, ax{a}.wf);
        end
    end
    
    if ~all([protocolMatch, chUnitsMatch, chNamesMatch, waveformMatch]) || ~all(holdingMatch)
        fprintf('File <%s> was not included \n', varargin{a})
    end
end
