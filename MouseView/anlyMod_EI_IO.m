function params = anlyMod_EI_IO(params)

% the goal is to have E/I ratios for each pulse (single or train) as a
% function of pulse intensity/freq/width (all the things that can be
% interleaved)

% approach: calculate conductance trace, and peak-by-pulse for each data
% file within a "group". Figure out which ax-files correspond to
% excitation, and to inhibition. Then match the condition types accordingly
% (so that excitation and inhibition are compared for the same stimulus).
% The resulting data could be a list of stimulus conditions and a list of
% E/I ratios associated with each condition.

% remember that each channel can have a unique holding potential for E or
% I.


% Each pulse should have its own analysis window, from 3 ms after the
% pulse, to about 15 ms after the pulse. the upper limit should not encroch
% into the next pulse epoch (if present).



% initialize the output arguments. this is important b/c I'll concatenate
% into these arrays, and they must exist prior to the first concatenation.
nGroups = size(params.isolatedCurrents, 1);
for i_group = 1:nGroups;
    group = params.isolatedCurrents{i_group,1};
    params.avg.peak_nS.(group).cond = {[] []};
    params.avg.peak_nS.(group).vals = {{} {}};
end


% loop through each file, condition, and channel.
nFiles = numel(params.files{2});
for i_fid = 1:nFiles;
    
    nConds = size(params.tdict{i_fid}.conds, 1);
    for i_cond = 1:nConds;
        
        % find the pulse on indicies, and identify an analysis window that
        % is appropriate for this condition
        t_anlyWindowStart = 0.002; % in sec
        t_anlyWindowEnd = 0.010;
        t_anlyStartIdx = floor(t_anlyWindowStart .* params.ax{i_fid}.head.sampRate);
        t_anlyEndIdx = ceil(t_anlyWindowEnd .* params.ax{i_fid}.head.sampRate);
        
        nCh = size(params.avg.trace_pA{i_fid},2);
        for i_ch = 1:nCh
            
            % figure out what the driving force is, and make sure that the
            % Vhold is close enough to the desired value specified by the
            % "isolatedCurrents" field in the physology_notes.m
            vhold = params.avg.vhold{i_fid}{i_cond, i_ch};
            drivingForce = [];
            fieldName = [];
            match = false;
            i_rev = 1;
            while ~match && i_rev<=size(params.isolatedCurrents,1);
                Erev = params.isolatedCurrents{i_rev, 3};
                if numel(Erev)>1
                    Erev = Erev(i_ch);
                end
                
                match = abs(vhold-Erev)<1;
                if match
                    drivingForce = Erev - params.isolatedCurrents{i_rev, 4};
                    drivingForce = abs(drivingForce); % in mV
                    fieldName = params.isolatedCurrents{i_rev, 1};
                end
                
                i_rev = i_rev+1;
            end
            
            if ~match
                continue
            end
            
            % pull out the max conductance for each pulse
            pOnIdx = params.tdict{i_fid}.pOnIdx{i_cond, i_ch};
            nPulses = numel(pOnIdx);
            maxval_nS = nan(1, nPulses);
            for i_pulse = 1:nPulses
                
                % define an analysis window for each pulse
                if nPulses>1
                    assert(t_anlyEndIdx < (pOnIdx(2)-pOnIdx(1)), 'ERROR: analysis window includes following pulse interval');
                end
                anlyWindow = pOnIdx(i_pulse)+t_anlyStartIdx : pOnIdx(i_pulse)+t_anlyEndIdx;
                snippet = params.avg.trace_pA{i_fid}{i_cond, i_ch}(anlyWindow);
                
                % now determine the maximal deviation from zero (positive
                % or negative). 
                [~, maxidx] = max(abs(snippet));
                maxval_pA = snippet(maxidx);
                
                % calculate the bkgnd pA
                bkgnd_pA = mean(params.avg.trace_pA{i_fid}{i_cond, i_ch}(anlyWindow(1)-21:anlyWindow(1)));
                
                % make sure the actual current takes into account the
                % decaying current it may be riding on
                maxval_pA = maxval_pA - bkgnd_pA;
                
                % convert to conductance
                maxval_amps = maxval_pA ./ 1e12;
                drivingForce_volts = drivingForce ./ 1e3;
                maxval_siemans = maxval_amps ./ drivingForce_volts;
                maxval_nS(i_pulse) = maxval_siemans .* 1e9;
                
                
            end
            
            % store the maxval conductance once all the pulses have been
            % accounted for
            idx = size(params.avg.peak_nS.(fieldName).cond{i_ch},1) + 1;
            params.avg.peak_nS.(fieldName).cond{i_ch}(idx,:) = params.tdict{i_fid}.conds(i_cond,:);
            params.avg.peak_nS.(fieldName).vals{i_ch}{idx,1} = maxval_nS;
            
            
        end % i_ch
    end % i_cond
end % i_fid







% 
% OPTIONAL SUMMARY FIGURE 1: PULSE AMPLITUDE IS THE 
% ONLY VARIABLE THAT'S INTERLEAVED
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure out if pulse amp is the only variable interleaved
bigcondmtx = [];
for i_group = 1:nGroups;
    for i_ch = 1:2;
        group = params.isolatedCurrents{i_group,1};
        bigcondmtx = cat(1, bigcondmtx, params.avg.peak_nS.(group).cond{i_ch});
    end
end
nAmps = numel(unique(bigcondmtx(:,1)));
nWidths = numel(unique(bigcondmtx(:,2)));
nFreqs = numel(unique(bigcondmtx(:,3)));

onlyPulseAmp = (nAmps>1) && (nWidths==1) && (nFreqs==1);
if onlyPulseAmp

    for i_group = 1:nGroups
        
        group = params.isolatedCurrents{i_group,1};
        
        figure
        set(gcf, 'name', group, 'position', [373    17   452   789]);
        plothelper(gcf);
        
        
        for i_ch = 1:2
            subplot(2,1,i_ch)
            if ~isempty(params.avg.peak_nS.(group).cond{i_ch})
                
                % sort the data and then plot
                x = params.avg.peak_nS.(group).cond{i_ch}(:,1);
                [x, idx] = sort(x);
                y = cell2mat(params.avg.peak_nS.(group).vals{i_ch});
                y = y(idx);
                plot(abs(x), abs(y), '-ko', 'markerfacecolor', 'k');
                
            end
            
            % putting this outside the 'if' statement so that axes get titles
            % even when there's no data...
            xlabel('pulse amplitude (V)')
            ylabel('conductance (nS)')
            title(sprintf('Channel %d', i_ch))
            
            
        end
    end
    
    
    % plot the E/I ratio
    figure
    set(gcf, 'name', 'EI ratio', 'position', [373    17   452   789]);
    plothelper(gcf);
    
    for i_ch = 1:2
        subplot(2,1,i_ch)
        if ~isempty(params.avg.peak_nS.excit.vals{i_ch})
            
            excit_nS = cell2mat(params.avg.peak_nS.excit.vals{i_ch});
            excit_ledV = params.avg.peak_nS.excit.cond{i_ch}(:,1);
            inhib_nS = cell2mat(params.avg.peak_nS.inhib.vals{i_ch});
            inhib_ledV = params.avg.peak_nS.inhib.cond{i_ch}(:,1);
            
            % make sure that the nS values are sorted identically for the excit and
            % inhib conditions
            [excit_ledV, idx] = sort(excit_ledV, 'ascend');
            inhib_ledV = inhib_ledV(idx);
            assert(all(inhib_ledV == excit_ledV), 'ERROR: no correspondence b/w excit and inhib LED pow');
            
            ei_ratio = abs(excit_nS(idx))./abs(inhib_nS(idx));
            plot(excit_ledV, ei_ratio, '-ko', 'markerfacecolor', 'k');
            set(gca, 'yscale', 'log')
        end
        xlabel('Pulse Voltage (V)')
        ylabel('EI ratio')
        title(sprintf('Channel %d', i_ch))
    end

end



% 
% OPTIONAL SUMMARY FIGURE 1: MULTIPLE FREQUENCIES AND POSSIBLY MULTIPLE
% PULSE AMPLITUDES. SAME PULSE WIDTH.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if (nWidths==1) && (nFreqs>1);
    
    
    % simple plot of Pn:P1 conductance as a function of pulse
    % number
    for i_group = 1:nGroups
        
        group = params.isolatedCurrents{i_group,1};
        
        figure
        set(gcf, 'name', group, 'position', [373    17   452   789]);
        plothelper(gcf);
        
        
        for i_ch = 1:2
            subplot(2,1,i_ch), hold on,
            if ~isempty(params.avg.peak_nS.(group).cond{i_ch})
                
                % pull out the data 
                tmp_cond = params.avg.peak_nS.(group).cond{i_ch};
                tmp_vals = params.avg.peak_nS.(group).vals{i_ch}; % a cell array
                
                
                nconds = size(tmp_cond,1);
                legtext = {};
                cmap = colormap(jet);
                cidx = round(linspace(1, size(cmap,1),nconds));
                clrs = cmap(cidx,:);
                for i_cond = 1:nconds
                    y = tmp_vals{i_cond};
                    plot(1:numel(y), y./y(1), '-o', 'color', clrs(i_cond,:), 'markerfacecolor', clrs(i_cond,:));
                    legtext = cat(2, legtext, sprintf('%.1f V, %.3f ms, %.0f Hz',...
                                  tmp_cond(i_cond,1), tmp_cond(i_cond, 2).*1000, tmp_cond(i_cond,3)));
                end
                set(gca, 'yscale', 'log')
                legend(legtext)
                legend boxoff
            end
            
            xlabel('pulse number')
            ylabel('conductance (nS)')
            title(sprintf('Channel %d', i_ch))
            
        end
        
    end
    
    
    
    % plot the E/I ratio as a function of pulse number
    figure
    set(gcf, 'name', 'DI ratio', 'position', [373    17   452   789]);
    plothelper(gcf);
    for i_ch = 1:2
        subplot(2,1,i_ch), hold on,
        if ~isempty(params.avg.peak_nS.excit.vals{i_ch})
            
            % pull out the data
            excit_nS = params.avg.peak_nS.excit.vals{i_ch};
            excit_conds = params.avg.peak_nS.excit.cond{i_ch};
            inhib_nS = params.avg.peak_nS.inhib.vals{i_ch};
            inhib_conds = params.avg.peak_nS.inhib.cond{i_ch};
            
            % make sure the experimental conditions were identical between
            % the excitation and inhibition files, and then sort the raw
            % data accordingly
            assert(size(inhib_conds,1) == size(excit_conds,1), 'ERROR: condition mismatch between excit and inhib')
            sortidx = nan(size(inhib_conds,1),1);
            for a = 1:size(excit_conds,1)
               idx = ismember(excit_conds,inhib_conds(a,:), 'rows');
               assert(sum(idx)==1, 'ERROR: too many matches found')
               sortidx(a) = find(idx);
            end
            
            % since I'm using the excit_conds as a template, reorder just
            % the inhib_conds and _vals
            inhib_nS = inhib_nS(sortidx);
            inhib_conds = inhib_conds(sortidx, :);
            
            
            % now calculate the E/I ratio and plot
            nconds = size(excit_conds,1);
            legtext = {};
            cmap = colormap(jet);
            cidx = round(linspace(1, size(cmap,1),nconds));
            clrs = cmap(cidx,:);
            for i_cond = 1:nconds
                
                
                ei_ratio = abs(excit_nS{i_cond})./abs(inhib_nS{i_cond});
                plot(1:numel(ei_ratio), ei_ratio, '-o', 'color', clrs(i_cond,:), 'markerfacecolor', clrs(i_cond,:));
                
                
                legtext = cat(2, legtext, sprintf('%.1f V, %.3f ms, %.0f Hz',...
                    tmp_cond(i_cond,1), tmp_cond(i_cond, 2).*1000, tmp_cond(i_cond,3)));
            end
            set(gca, 'yscale', 'log')
            legend(legtext, 'location', 'northwest')
            legend boxoff
        end
        
        xlabel('pulse number')
        ylabel('E/I ratio')
        title(sprintf('Channel %d', i_ch))
        
    end
    
    
end

