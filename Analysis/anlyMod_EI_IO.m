function params = anlyMod_EI_IO(params)

% the goal is to have E/I ratios for each pulse (single or train) as a
% function of pulse intensity/freq/width (all the things that can be
% interleaved)

nExptGroups = size(params.groups,1);
for i_group = 1:nExptGroups;
    
    % assign a field name for the structure
    group = params.groups{i_group, 1};
    
    % initialize the output arguments. this is important b/c I'll
    % concatenate into these arrays, and they must exist prior to the first
    % concatenation. The reason I concatenate is b/c I don't know the
    % dimensionality a priori
    n_pA_types = size(params.isolatedCurrents, 1);
    for i_pAType = 1:n_pA_types;
        pA_type = params.isolatedCurrents{i_pAType,1};
        params.(group).avg.peak_nS.(pA_type).cond = {[] []};
        params.(group).avg.peak_nS.(pA_type).vals = {{} {}};
    end
    
    
    % loop through each file, condition, and channel.
    nFiles = numel(params.(group).tdict);
    for i_fid = 1:nFiles;
        
        nConds = size(params.(group).tdict{i_fid}.conds, 1);
        for i_cond = 1:nConds;
            
            % find the pulse on indicies, and identify an analysis window that
            % is appropriate for this condition
            t_anlyWindowStart = 0.002; % in sec
            t_anlyWindowEnd = 0.010;
            t_anlyStartIdx = floor(t_anlyWindowStart .* params.(group).head{i_fid}.sampRate);
            t_anlyEndIdx = ceil(t_anlyWindowEnd .* params.(group).head{i_fid}.sampRate);
            
            nCh = size(params.(group).avg.trace_pA{i_fid},2);
            for i_ch = 1:nCh
                
                % figure out what the driving force is, and make sure that the
                % Vhold is close enough to the desired value specified by the
                % "isolatedCurrents" field in the physology_notes.m
                vhold = params.(group).avg.vhold{i_fid}{i_cond, i_ch};
                if isempty(vhold)
                    continue % no data for this channel/condition
                end
                drivingForce = [];
                pA_type = [];
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
                        pA_type = params.isolatedCurrents{i_rev, 1};
                    end
                    
                    i_rev = i_rev+1;
                end
                
                
                if ~match
                    continue % no data for this channel/condition were collected at the correct vhold...
                end
                
                % pull out the max conductance for each pulse
                pOnIdx = params.(group).tdict{i_fid}.pOnIdx{i_cond, i_ch};
                nPulses = numel(pOnIdx);
                maxval_nS = nan(1, nPulses);
                for i_pulse = 1:nPulses
                    
                    % define an analysis window for each pulse
                    if nPulses>1
                        assert(t_anlyEndIdx < (pOnIdx(2)-pOnIdx(1)), 'ERROR: analysis window includes following pulse interval');
                    end
                    anlyWindow = pOnIdx(i_pulse)+t_anlyStartIdx : pOnIdx(i_pulse)+t_anlyEndIdx;
                    snippet = params.(group).avg.trace_pA{i_fid}{i_cond, i_ch}(anlyWindow);
                    
                    % now determine the maximal deviation from zero (positive
                    % or negative).
                    [~, maxidx] = max(abs(snippet));
                    maxval_pA = snippet(maxidx);
                    
                    % calculate the bkgnd pA
                    switch group
                        case 'pulses'
                            Nbkgndpts = ceil(5e-3 .* params.(group).head{i_fid}.sampRate);
                        case 'trains'
                            Nbkgndpts = ceil(500e-6 .* params.(group).head{i_fid}.sampRate);
                    end
                    bkgnd_pA = mean(params.(group).avg.trace_pA{i_fid}{i_cond, i_ch}(pOnIdx-Nbkgndpts : pOnIdx));
                    
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
                idx = size(params.(group).avg.peak_nS.(pA_type).cond{i_ch},1) + 1;
                params.(group).avg.peak_nS.(pA_type).cond{i_ch}(idx,:) = params.(group).tdict{i_fid}.conds(i_cond,:);
                params.(group).avg.peak_nS.(pA_type).vals{i_ch}{idx,1} = maxval_nS;
                
                
            end % i_ch
        end % i_cond
    end % i_fid
end % i_group


%
% Stop here and make sure that there are no duplicate conditions. This
% could happen if I collect two datasets with overlapping conditions. When
% there are duplicates, consolidate the datasets, but plot a figure of each
% condition separately just as a sanity check.
%
%%%%%%%%%%%%%%%%%%%%%%
for i_group = 1:nExptGroups;
    group = params.groups{i_group, 1};
    
    for i_pAType = 1:n_pA_types
        pA_type = params.isolatedCurrents{i_pAType,1};
        
        for i_ch = 1:2
            tmp_conds = params.(group).avg.peak_nS.(pA_type).cond{i_ch};
            tmp_raw = params.(group).avg.peak_nS.(pA_type).vals{i_ch};
            
            %look for duplicates
            unique_conds = unique(tmp_conds, 'rows');
            out_raw = {};
            for i_cond = 1:size(unique_conds,1)
                list = ismember(tmp_conds, unique_conds(i_cond,:), 'rows');
                if sum(list)>1
                    warning('Averaging across duplicate conditions')
                end
                
                valsToAverage = cell2mat(tmp_raw(list));
                out_raw = cat(1, out_raw, {mean(valsToAverage, 1)});
            end
            
            % repackage the 'conds' and 'vals' fields
            params.(group).avg.peak_nS.(pA_type).cond{i_ch} = unique_conds;
            params.(group).avg.peak_nS.(pA_type).vals{i_ch} = out_raw;
        end
    end
end


%
%
%  OPTIONAL FIGURES
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i_group = 1:nExptGroups;
    
    group = params.groups{i_group, 1};
    
    switch group
        case 'pulses'
            summaryFigure_1(params, group)
        case 'trains'
            summaryFigure_2(params, group)
    end
    
end

end %  MAIN FUNCTION



function [idx_excit, idx_inhib] = sortconditions(excit_conds, inhib_conds)
    
    % make a 1to1 correspondence b/w conditions in the excitation and
    % inhibition data arrays. This might be difficult becuase not all
    % trial types will exist and they may be out of order
    
    % by default, use the excit_conds as the template, but switch to
    % inhib_conds if there were fewer inhib_conds
    if size(excit_conds,1) <= size(inhib_conds,1)
        template_one = excit_conds;
        template_two = inhib_conds;
        normaloutput = true;
    else
        template_one = inhib_conds;
        template_two = excit_conds;
        normaloutput = false;
    end
    
    
    % iterate through "template_one" looking for the corresponding
    % condition in "template_two"
    nConds = size(template_one, 1);
    newidx = nan(nConds,1);
    for i_cond = 1:nConds
        list = ismember(template_two, template_one(i_cond,:), 'rows');
        assert(sum(list)==1, 'ERROR: more than 1 match was found');
        newidx(i_cond) = find(list==1);
    end

    
    % asemble the outputs
    if normaloutput
        idx_excit = [1:size(excit_conds,1)]';
        idx_inhib = newidx;
    else
        idx_inhib = [1:size(inhib_conds,1)]';
        idx_excit = newidx;
    end

end

function summaryFigure_1(params, group)

    n_pA_types = size(params.isolatedCurrents,1);
    for i_pAType = 1:n_pA_types
        
        pA_type = params.isolatedCurrents{i_pAType,1};
        
        figure
        set(gcf, 'name', [group, ': ', pA_type], 'position', [373    17   452   789]);
        
        
        for i_ch = 1:2
            subplot(2,1,i_ch)
            if ~isempty(params.(group).avg.peak_nS.(pA_type).cond{i_ch})
                
                % sort the data and then plot
                x = params.(group).avg.peak_nS.(pA_type).cond{i_ch}(:,1);
                [x, idx] = sort(x);
                y = cell2mat(params.(group).avg.peak_nS.(pA_type).vals{i_ch});
                y = y(idx);
                plot(x, abs(y), '-ko', 'markerfacecolor', 'k');
                
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
    set(gcf, 'name', [group, ': E/I ratio'], 'position', [373    17   452   789]);
    for i_ch = 1:2
        subplot(2,1,i_ch)
        excitPresent = ~isempty(params.(group).avg.peak_nS.excit.vals{i_ch});
        inhibPresent = ~isempty(params.(group).avg.peak_nS.inhib.vals{i_ch});
        if excitPresent && inhibPresent
            
            excit_nS = cell2mat(params.(group).avg.peak_nS.excit.vals{i_ch});
            excit_ledV = params.(group).avg.peak_nS.excit.cond{i_ch}(:,1);
            inhib_nS = cell2mat(params.(group).avg.peak_nS.inhib.vals{i_ch});
            inhib_ledV = params.(group).avg.peak_nS.inhib.cond{i_ch}(:,1);
            
            % make sure that the nS values are sorted identically for the excit and
            % inhib conditions
            [idx_excit, idx_inhib] = sortconditions(excit_ledV, inhib_ledV);
            excit_nS = excit_nS(idx_excit);
            excit_ledV = excit_ledV(idx_excit);
            inhib_nS = inhib_nS(idx_inhib);
            inhib_ledV = inhib_ledV(idx_inhib);
            
            assert(all(excit_ledV==inhib_ledV), 'ERROR: conditions are mismatched')
            
            ei_ratio = (abs(excit_nS)-abs(inhib_nS)) ./ (abs(excit_nS)+abs(inhib_nS)) ;
            plot(excit_ledV, ei_ratio, '-ko', 'markerfacecolor', 'k');
        end
        xlabel('Pulse Voltage (V)')
        ylabel('(E-I)/(E+I)')
        title(sprintf('Channel %d', i_ch))
    end


end

function summaryFigure_2(params, group)
 
    % simple plot of Pn:P1 conductance as a function of pulse
    % number
    n_pA_types = size(params.isolatedCurrents,1);
    for i_pAType = 1:n_pA_types
        
        pA_type = params.isolatedCurrents{i_pAType,1};
        
        figure
        set(gcf, 'name', [group, ': ', pA_type], 'position', [373    17   452   789]);
        
        
        for i_ch = 1:2
            subplot(2,1,i_ch), hold on,
            if ~isempty(params.(group).avg.peak_nS.(pA_type).cond{i_ch})
                
                % pull out the data 
                tmp_cond = params.(group).avg.peak_nS.(pA_type).cond{i_ch};
                tmp_vals = params.(group).avg.peak_nS.(pA_type).vals{i_ch}; % a cell array
                
                
                nconds = size(tmp_cond,1);
                legtext = {};
                cmap = colormap(jet);
                cidx = round(linspace(1, size(cmap,1),nconds));
                clrs = cmap(cidx,:);
                for i_cond = 1:nconds
                    y = tmp_vals{i_cond};
                    TF = tmp_cond(i_cond, 3);
                    dt = 1/TF;
                    tt = 0 : dt : (dt*(numel(y)-1));
                    tt = tt.* 1000;
                    plot(tt, y./y(1), '-o', 'color', clrs(i_cond,:), 'markerfacecolor', clrs(i_cond,:));
                    legtext = cat(2, legtext, sprintf('%.1f V, %.3f ms, %.0f Hz',...
                                  tmp_cond(i_cond,1), tmp_cond(i_cond, 2).*1000, tmp_cond(i_cond,3)));
                end
                set(gca, 'yscale', 'log')
                legend(legtext)
                legend boxoff
            end
            
            xlabel('pulse time (ms)')
            ylabel('conductance (nS)')
            title(sprintf('Channel %d', i_ch))
            
        end
        
    end
    
    
    
    % plot the E/I ratio as a function of pulse number
    figure
    set(gcf, 'name', [group, ': E/I ratio'], 'position', [373    17   452   789]);
    for i_ch = 1:2
        subplot(2,1,i_ch), hold on,
        excitPresent = ~isempty(params.(group).avg.peak_nS.excit.vals{i_ch});
        inhibPresent = ~isempty(params.(group).avg.peak_nS.inhib.vals{i_ch});
        if excitPresent && inhibPresent
            
            % pull out the data
            excit_nS = params.(group).avg.peak_nS.excit.vals{i_ch};
            excit_conds = params.(group).avg.peak_nS.excit.cond{i_ch};
            inhib_nS = params.(group).avg.peak_nS.inhib.vals{i_ch};
            inhib_conds = params.(group).avg.peak_nS.inhib.cond{i_ch};
            
            % make sure the experimental conditions were identical between
            % the excitation and inhibition files, and then sort the raw
            % data accordingly
            [idx_excit, idx_inhib] = sortconditions(excit_conds, inhib_conds);
            excit_nS = excit_nS(idx_excit);
            excit_conds = excit_conds(idx_excit,:);
            inhib_nS = inhib_nS(idx_inhib);
            inhib_conds = inhib_conds(idx_inhib,:);
            
            assert(all(excit_conds(:)==inhib_conds(:)), 'ERROR: conditions mismatch');
            
            % now calculate the E/I ratio and plot
            nconds = size(excit_conds,1);
            legtext = {};
            cmap = colormap(jet);
            cidx = round(linspace(1, size(cmap,1),nconds));
            clrs = cmap(cidx,:);
            for i_cond = 1:nconds
                
                tmp_e = abs(excit_nS{i_cond});
                tmp_i = abs(inhib_nS{i_cond});
                ei_ratio = (tmp_e - tmp_i) ./ (tmp_e + tmp_i);
                
                TF = excit_conds(i_cond,3);
                dt = 1/TF;
                tt = 0 : dt : (dt*(numel(y)-1));
                tt = tt.* 1000;
                
                plot(tt, ei_ratio, '-o', 'color', clrs(i_cond,:), 'markerfacecolor', clrs(i_cond,:));
                
                
                legtext = cat(2, legtext, sprintf('%.1f V, %.3f ms, %.0f Hz',...
                    excit_conds(i_cond,1), excit_conds(i_cond, 2).*1000, excit_conds(i_cond,3)));
            end
            legend(legtext, 'location', 'northwest')
            legend boxoff
        end
        
        xlabel('pulse number')
        ylabel('(E-I)/(E+I)')
        title(sprintf('Channel %d', i_ch))
        
    end
    
end



