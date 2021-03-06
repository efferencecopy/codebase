function params = anlyMod_EIbalance(params)


% make sure data are present
assert(isfield(params, 'ivdat'), 'ERROR: No optoIV data present')


% define an analysis window based off when the LED pulse turns on. Ignore
% the first 3.5 ms due to LED artifacts, ignore everything after 40 msec
% so that I'm only looking at the time period with the peak
l_anlyWindow = (params.ivdat.tvec > 0.0035) & (params.ivdat.tvec < 0.040);


% iterate over the isolated data field and convert the raw current traces
% to conductances.
for a = 1:size(params.isolatedCurrents, 1) % Num Vholds.
    currentType = params.isolatedCurrents{a, 1};
    experimentalGroup = params.isolatedCurrents{a, 2};
    
    nChanels = size(params.ivdat.(experimentalGroup).vhold, 2);
    for ch = 1:nChanels
        % assign an empty array for now
        params.isolatedData.(currentType).raw_nS{ch} = [];
        
        % look for the correct holding potential
        vholds_exp = params.ivdat.(experimentalGroup).vhold{ch}; % available vHolds
        vhold_req = params.isolatedCurrents{a, 3};
        if numel(vhold_req)>1; % allow the vholds to differ between HS1 and HS2
            vhold_req = vhold_req(ch);
        end
        vHoldAvailable = cellfun(@(x, y) softEq(x, y, 0), vholds_exp, repmat({vhold_req}, size(vholds_exp)), 'uniformoutput', false);
        vHoldAvailable = cellfun(@(x) ~isempty(x) && (x==true), vHoldAvailable);
        
        if ~any(vHoldAvailable); continue; end
        
        % pull out the raw data
        trace_pA = params.ivdat.(experimentalGroup).raw{ch}{vHoldAvailable};
        Erev = params.isolatedCurrents{a, 4};
        if numel(Erev)>1
            Erev = Erev(ch);
        end
        drivingForce_mV = abs(vhold_req - Erev); % in mV
        
        % if isolatedCurrent == 'ampa', then subtract off the NMDA current
        % in the presence of blockers.
        if strcmpi(currentType, 'ampa')
            vholds_nmdaOnly = params.ivdat.nbqxGabazine.vhold{ch};
            idx = cellfun(@(x, y) softEq(x, y, 0), vholds_nmdaOnly, repmat({vhold_req}, size(vholds_nmdaOnly)), 'uniformoutput', false);
            idx = cellfun(@(x) ~isempty(x) && (x==true), idx);
            if any(idx)
                trace_nmdaOnly = params.ivdat.nbqxGabazine.raw{ch}{idx};
                trace_pA = trace_pA - trace_nmdaOnly;
            else
                warning('Could not find control condition for AMPA only current')
                % just use the trace_pA calculated above
            end
        end
        
        % convert to amps and volts. 
        trace_amps = trace_pA ./ 1e12;
        drivingForce_volts = drivingForce_mV ./ 1e3;
        
        % calculate conductance
        trace_siemens = trace_amps ./ drivingForce_volts;
        trace_nS = trace_siemens .* 1e9;
        params.isolatedData.(currentType).raw_nS{ch} = trace_nS;
        
        
        % extract the max (unsigned) peak value
        peakVal_nS = max(abs(trace_nS(l_anlyWindow))); % unsigned value
        peak_idx = (abs(trace_nS) == peakVal_nS) & l_anlyWindow;
        
        %error checking
        assert(sum(peak_idx)==1, sprintf('ERROR: Too  many indicies found <CH %d, currentType: %s>', ch, upper(currentType)))
        if strcmpi(currentType, 'nmda')
            assert(vhold_req>20, 'ERROR: Vhold for NMDA current must be >20mV')
        end
        
        % save the signed value. Some of them will have the wrong sign for
        % the type of current (e.g., negative current for NMDA at
        % Vhold>25mV). Saving the signed nS value will help identify these
        % issues during subsequent data analysis.
        peakVal_nS = trace_nS(peak_idx);

        
        % store the data in the params struct
        params.isolatedData.(currentType).peak_nS{ch} = peakVal_nS;
        params.isolatedData.(currentType).peak_pA{ch} = trace_pA(peak_idx);
        
        % this is redundant with params.ivdat.(exptgroup).raw, but I'm
        % going to store the raw current trace in a form that is easily
        % accessible to down stream analysis that is specifically
        % interested in the isolated currents. Ditto for Vclamp err, and Ra
        params.isolatedData.(currentType).raw_pA{ch} = trace_pA;
        params.isolatedData.(currentType).peakBySweep_pA{ch} = params.ivdat.(experimentalGroup).peakBySweep_pA{ch}{vHoldAvailable};
        params.isolatedData.(currentType).Verr{ch} = params.ivdat.(experimentalGroup).Verr{ch}{vHoldAvailable};
        params.isolatedData.(currentType).Racc{ch} = params.ivdat.(experimentalGroup).Racc{ch}{vHoldAvailable};
        params.isolatedData.(currentType).holdingCurrent{ch} = params.ivdat.(experimentalGroup).holdingCurrent{ch}{vHoldAvailable};
        
        
    end
    
end



% plot the resulting conductances for E/I and NMDA/AMPA
figure % E/I balance
for ch = 1:nChanels
    subplot(nChanels, 1, ch), hold on,
    
    if ~isempty(params.isolatedData.excit.raw_nS{ch});
        plot(params.ivdat.tvec.*1000, params.isolatedData.excit.raw_nS{ch}, 'b', 'linewidth', 2)
    else
        continue
    end
    
    if ~isempty(params.isolatedData.inhib.raw_nS{ch});
        plot(params.ivdat.tvec.*1000, params.isolatedData.inhib.raw_nS{ch}, 'r', 'linewidth', 2)
    else
        continue
    end
    
    xlabel('time (ms)')
    ylabel('conductance (nS)')
    legend('excitation', 'inhibition')
    xlim([-25 300])
    
end


figure % NMDA/AMPA
for ch = 1:nChanels
    subplot(nChanels, 1, ch), hold on,
    if isfield(params.isolatedData, 'ampa') && ~isempty(params.isolatedData.ampa.raw_nS{ch});
        plot(params.ivdat.tvec.*1000, params.isolatedData.ampa.raw_nS{ch}, 'k', 'linewidth', 2)
    else
        continue
    end
    
    if isfield(params.isolatedData, 'nmda') && ~isempty(params.isolatedData.nmda.raw_nS{ch})
        plot(params.ivdat.tvec.*1000, params.isolatedData.nmda.raw_nS{ch}, 'color', [.4 .4 .4], 'linewidth', 2)
    else
        continue
    end
    
    xlabel('time (ms)')
    ylabel('conductance (nS)')
    legend('AMPA', 'NMDA')
    xlim([-25 300])
    
end



