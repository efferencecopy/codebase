function params = anlyMod_NMDAR_mac(params)

% need to pull out each file that contains NBQX and GABAZINE and plot peak
% pA as a function of Vhold. 

% what are the files within groups = NMDAR?

% iterate over channels
%   * iterate over Vholds (store in cell array)
%   * include a "corrected" Vhold that incorporates the Vclamp error



% define an analysis window based off when the LED pulse turns on. Ignore
% the first 1.25 ms due to LED artifacts.
l_anlyWindow = params.ivdat.tvec > 0.00125;


% pull out the raw data and determine the peak current
nCh = size(params.ivdat.NMDAR.raw,2);
for ch = 1:nCh;
   nVholds = size(params.ivdat.NMDAR.raw{ch},2);
   snippets = cellfun(@(x,y) x(y), params.ivdat.NMDAR.raw{ch}, repmat({l_anlyWindow}, 1, nVholds), 'uniformoutput', 0);
   max_pA = cellfun(@(x) max(abs(x)), snippets, 'uniformoutput', false); %strictly pos vals
   max_idx = cellfun(@(x,y) find(abs(x)==y), snippets, max_pA, 'uniformoutput', false);
   max_pA = cellfun(@(x,y) x(y), snippets, max_idx);
   NMDAR_ivdat_pA{ch} = max_pA;
   NMDAR_snippets{ch} = snippets; % only for local ploting and debugging.
end


% store the Vholds
NMDAR_ivdat_Vhold = cellfun(@(x) cat(2,x{:}), params.ivdat.NMDAR.vhold, 'uniformoutput', false);


% store the Verr corrected Vholds
for ch = 1:nCh
    Ra_Mohm = params.ivdat.NMDAR.Racc{ch};
    holding_pA = params.ivdat.NMDAR.holdingCurrent{ch};
    holding_nA = cellfun(@(x,y) x./y, holding_pA, repmat({1e3}, size(holding_pA)), 'uniformoutput', false);
    Verr_mV = cellfun(@(x,y) mean(x.*y), Ra_Mohm, holding_nA, 'uniformoutput', true);
    NMDAR_ivdat_Vhold_corrected{ch} = NMDAR_ivdat_Vhold{ch} - Verr_mV;
end


% make sure that the data are in increasing order of Vholds
for ch = 1:nCh
    [~, idx] = sort(NMDAR_ivdat_Vhold{ch});
    NMDAR_ivdat_Vhold{ch} = NMDAR_ivdat_Vhold{ch}(idx);
    NMDAR_ivdat_Vhold_corrected{ch} = NMDAR_ivdat_Vhold_corrected{ch}(idx);
    NMDAR_ivdat_pA{ch} = NMDAR_ivdat_pA{ch}(idx);
end


%
% Plot the results
%
%%%%%%%%%%%%%%%%%%%%%
figure
for a = 1:nCh
    subplot(nCh, 2, ((a-1).*2)+1)
    plot(params.ivdat.tvec(l_anlyWindow), cat(1,NMDAR_snippets{a}{:})')
    
    xlabel('Time from LED onset (sec)')
    ylabel('Current (pA)')
    
    subplot(nCh, 2, a.*2), hold on,
    plot(NMDAR_ivdat_Vhold{a}, NMDAR_ivdat_pA{a}, '-ko')
    plot(NMDAR_ivdat_Vhold_corrected{a}, NMDAR_ivdat_pA{a}, '-bo')
    xlabel('Voltage (mV)')
    ylabel('Current (pA)')
    
    
end




