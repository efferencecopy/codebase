function fig_1_tau(dat, plotgroups)



groupdata.tau_neg = repmat({[]}, 1, size(plotgroups, 1));
groupdata.tau_pos = repmat({[]}, 1, size(plotgroups, 1));

for i_ex = 1:numel(dat)
    
    for i_ch = 1:2
        
        % check to make sure this neuron was defined
        hs_name = sprintf('dcsteps_hs%d', i_ch);
        isvalid = ~strncmp(dat{i_ex}.info.fid.(hs_name)(end-8:end), [filesep, 'none.abf'], 9);
        if ~isvalid
            continue
        end
        
        % check to make sure the cell was healthy (as based on presence of
        % holding current)
        max_holding_pa = abs(dat{i_ex}.dcsteps.pA_holding{i_ch}(2));
        if max_holding_pa > 20; continue; end
        
        % check the attributes
        ch_attribs = {dat{i_ex}.info.cellType{i_ch}, dat{i_ex}.info.brainArea, dat{i_ex}.info.opsin};
        group_idx = groupMatcher(plotgroups, ch_attribs);
        if sum(group_idx) == 0; continue; end
        
        % add the statistics to the appropriate cell array
        if isfield(dat{i_ex}.dcsteps, 'tau') && ~isempty(dat{i_ex}.dcsteps.tau{i_ch})
            tmp_all_tau = dat{i_ex}.dcsteps.tau{i_ch}(:,3);
            tmp_Icmd = dat{i_ex}.dcsteps.tau{i_ch}(:,1);
            
            if any(tmp_Icmd==30) && any(tmp_Icmd == -30)
                tau_pos_30pA = tmp_all_tau(tmp_Icmd==30);
                tau_neg_30pA = tmp_all_tau(tmp_Icmd == -30);
                
                groupdata.tau_neg{group_idx} = cat(1, groupdata.tau_neg{group_idx}, tau_neg_30pA);
                groupdata.tau_pos{group_idx} = cat(1, groupdata.tau_pos{group_idx}, tau_pos_30pA);
            end
        end
    end
end


f = figure; hold on,
f.Position = [440   156   561   642];
allTaus = cat(1, groupdata.tau_pos{:}, groupdata.tau_neg{:});
edges_taus = linspace(min(allTaus)-0.002, 0.050, 50).*1000;
legtext = {};
for i_group = 1:numel(groupdata.tau_pos)
    hva_name = plotgroups{i_group, 3};
    grp_clr = hvaPlotColor(hva_name);
    
    subplot(2,1,1), hold on,
    N = histcounts(groupdata.tau_pos{i_group}.*1000, edges_taus);
    N = [0, N];
    cdf_vals = cumsum(N)./sum(N);
    stairs(edges_taus, cdf_vals, '-', 'color', grp_clr)
    num_cells = numel(groupdata.tau_pos{i_group});
    legtext{i_group} = sprintf('%s, %s, %s, n=%d', plotgroups{i_group,1},plotgroups{i_group, 2},plotgroups{i_group,3}, num_cells);
    force_consistent_figs(gca, 'ax');
    xlabel('Time constant (ms)')
    ylabel('Proportion')
    title('Tau from +30 pA cmd')
    
    
    subplot(2,1,2), hold on
    N = histcounts(groupdata.tau_neg{i_group}.*1000, edges_taus);
    N = [0, N];
    cdf_vals = cumsum(N)./sum(N);
    stairs(edges_taus, cdf_vals, '-', 'color', grp_clr)
    force_consistent_figs(gca, 'ax');
    xlabel('Time constant (ms)')
    ylabel('Proportion')
    title('Tau from -30 pA cmd')
end
hl = legend(legtext, 'Location', 'best');
hl.Interpreter = 'none';
legend boxoff






