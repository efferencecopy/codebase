function fig_tau(dat, plotgroups)



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
            tmp_all_tau = dat{i_ex}.dcsteps.tau{i_ch}(:,4);
            tmp_Icmd = dat{i_ex}.dcsteps.tau{i_ch}(:,1);
            
            l_pos_cmd = (tmp_Icmd > 25) & (tmp_Icmd < 35);
            if any(l_pos_cmd) 
                tau_pos_30pA = tmp_all_tau(l_pos_cmd);
                groupdata.tau_pos{group_idx} = cat(1, groupdata.tau_pos{group_idx}, tau_pos_30pA);
            end 
            l_neg_cmd = (tmp_Icmd < -25) & (tmp_Icmd > -35);
            if any(l_neg_cmd)
                tau_neg_30pA = tmp_all_tau(l_neg_cmd);
                groupdata.tau_neg{group_idx} = cat(1, groupdata.tau_neg{group_idx}, tau_neg_30pA);
            end
        end
    end
end


f = figure; hold on,
f.Position = [440   156   561   642];
allTaus = cat(1, groupdata.tau_pos{:}, groupdata.tau_neg{:});
edges_taus = linspace(min(allTaus)-0.002, 0.050, 50).*1000;
legtext_neg = {};
legtext_pos = {};
for i_group = 1:numel(groupdata.tau_pos)
    hva_name = plotgroups{i_group, 3};
    grp_clr = hvaPlotColor(hva_name);
    
    subplot(2,1,1), hold on,
    N = histcounts(groupdata.tau_pos{i_group}.*1000, edges_taus);
    N = [0, N];
    cdf_vals = cumsum(N)./sum(N);
    stairs(edges_taus, cdf_vals, '-', 'color', grp_clr)
    num_cells = numel(groupdata.tau_pos{i_group});
    legtext_pos{i_group} = sprintf('%s, %s, %s, n=%d', plotgroups{i_group,1},plotgroups{i_group, 2},plotgroups{i_group,3}, num_cells);
    force_consistent_figs(gca, 'ax');
    xlabel('Time constant (ms)')
    ylabel('Proportion')
    title('Tau from +30 pA cmd')
    
    
    subplot(2,1,2), hold on
    N = histcounts(groupdata.tau_neg{i_group}.*1000, edges_taus);
    N = [0, N];
    cdf_vals = cumsum(N)./sum(N);
    stairs(edges_taus, cdf_vals, '-', 'color', grp_clr)
    num_cells = numel(groupdata.tau_neg{i_group});
    legtext_neg{i_group} = sprintf('%s, %s, %s, n=%d', plotgroups{i_group,1},plotgroups{i_group, 2},plotgroups{i_group,3}, num_cells);
    force_consistent_figs(gca, 'ax');
    xlabel('Time constant (ms)')
    ylabel('Proportion')
    title('Tau from -30 pA cmd')
end
subplot(2,1,1)
hl = legend(legtext_pos, 'Location', 'best');
hl.Interpreter = 'none';
legend boxoff
subplot(2,1,2)
h2 = legend(legtext_neg, 'Location', 'best');
h2.Interpreter = 'none';
legend boxoff

%
%  box plot of Rin
%
%%%%%%%%%%%%%%%%%%%%%
box_groups = [];
box_colors = [];
box_data = [];
for i_group = 1:numel(groupdata.tau_neg)
    
    hva_name = plotgroups{i_group, 3};
    box_colors = cat(1, box_colors, hvaPlotColor(hva_name));
    
    N = numel(groupdata.tau_neg{i_group});
    box_data = cat(1, box_data, groupdata.tau_neg{i_group} .* 1000);
    box_groups = cat(1, box_groups, repmat(hva_name, N, 1));
end
hf = figure;
hb = boxplot(box_data, box_groups,...
             'grouporder', plotgroups(:,3),...
             'colors', box_colors,...
             'notch', 'on',...
             'symbol', '',...
             'whisker', 0);
set(gca, 'tickdir', 'out', 'fontsize', 24, 'ylim', [0, 18])
set(hb, 'linewidth', 3)
ylabel(sprintf('Membrane Time Constant\n(ms)'))
legend(gca, plotgroups(:,3)')



