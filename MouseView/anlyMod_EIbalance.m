function params = anlyMod_EIbalance(params)


% make sure data are present
if ~isfield(params, 'ivdat')
    fprintf('No optoIV data present')
end

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
        vHoldAvailable = cellfun(@(x, y) softEq(x, y, 0), vholds_exp, repmat({vhold_req}, size(vholds_exp)), 'uniformoutput', false);
        vHoldAvailable = cellfun(@(x) ~isempty(x) && (x==true), vHoldAvailable);
        
        if ~any(vHoldAvailable); continue; end
        
        % pull out the raw data
        trace_pA = params.ivdat.(experimentalGroup).raw{ch}{vHoldAvailable};
        drivingForce_mV = abs(params.isolatedCurrents{a, 3} - params.isolatedCurrents{a, 4}); % in mV
        
        % convert to amps and volts. 
        trace_amps = trace_pA ./ 1e12;
        drivingForce_volts = drivingForce_mV ./ 1e3;
        
        % calculate conductance
        trace_siemens = trace_amps ./ drivingForce_volts;
        trace_nS = trace_siemens .* 1e9;
        params.isolatedData.(currentType).raw_nS{ch} = trace_nS;
        
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
    if ~isempty(params.isolatedData.ampa.raw_nS{ch});
        plot(params.ivdat.tvec.*1000, params.isolatedData.ampa.raw_nS{ch}, 'k', 'linewidth', 2)
    else
        continue
    end
    
    if ~isempty(params.isolatedData.nmda.raw_nS{ch})
        plot(params.ivdat.tvec.*1000, params.isolatedData.nmda.raw_nS{ch}, 'color', [.4 .4 .4], 'linewidth', 2)
    else
        continue
    end
    
    xlabel('time (ms)')
    ylabel('conductance (nS)')
    legend('AMPA', 'NMDA')
    xlim([-25 300])
end



