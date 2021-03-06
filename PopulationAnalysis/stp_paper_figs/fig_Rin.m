function fig_Rin(dat, plotgroups)

groupdata.Rin_peak = repmat({[]}, 1, size(plotgroups, 1)); % should only have N cells, where N = size(plotgroups, 1). Each cell has a matrix with a cononicalGrid:
groupdata.Rin_asym = repmat({[]}, 1, size(plotgroups, 1)); % should only have N cells, where N = size(plotgroups, 1). Each cell has a matrix with a cononicalGrid:
groupdata.R_pos_30pA = repmat({[]}, 1, size(plotgroups, 1));
groupdata.R_neg_30pA = repmat({[]}, 1, size(plotgroups, 1));


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
        groupdata.Rin_peak{group_idx} = cat(1, groupdata.Rin_peak{group_idx}, dat{i_ex}.dcsteps.IVpeak.Rin{i_ch});
        groupdata.Rin_asym{group_idx} = cat(1, groupdata.Rin_asym{group_idx}, dat{i_ex}.dcsteps.IVasym.Rin{i_ch});
        
        % now calculate the +/- 30pA Rins
        v_rest = dat{i_ex}.dcsteps.Vrest{i_ch};
        Icmd = dat{i_ex}.dcsteps.IVasym.raw{i_ch}(:,1);
        delta_mV = dat{i_ex}.dcsteps.IVasym.raw{i_ch}(:,2) - v_rest;
        MOhm = abs(delta_mV  ./ (Icmd./1000));
        
        l_pos_cmd = (Icmd > 25) & (Icmd < 35);
        if any(l_pos_cmd)
            groupdata.R_pos_30pA{group_idx} = cat(1, groupdata.R_pos_30pA{group_idx}, MOhm(l_pos_cmd));
        end
        
        l_neg_cmd = (Icmd < -25) & (Icmd > -35);
        if any(l_neg_cmd)
            groupdata.R_neg_30pA{group_idx} = cat(1, groupdata.R_neg_30pA{group_idx}, MOhm(l_neg_cmd));
        end
    end
end


%
% plot histograms of Input resistance measured two different ways
%
%%%%%%%%%%%%%%%%%%%%%%%%%
f=figure;
Ngroups = (numel(groupdata.Rin_peak)); % adding one for "summary" cdf fig
N_plt_rows = Ngroups+1;
f.Position = [506   158   838   749];
allNums = cat(1, groupdata.Rin_peak{:}, groupdata.Rin_asym{:});
edges = linspace(min(allNums)-1, max(allNums)+1, 50);
for i_group = 1:numel(groupdata.Rin_asym)
    
    hva_name = plotgroups{i_group, 3};
    grp_clr = hvaPlotColor(hva_name);
    
    %
    % current clamp, peak vals
    %%%%%%%%%%%%%%%%%%
    pltidx = sub2ind([2, N_plt_rows], 1, i_group);
    subplot(N_plt_rows, 2, pltidx)
    
    h = histogram(groupdata.Rin_peak{i_group}, edges);
    h.FaceColor = grp_clr;
    xbar = nanmean(groupdata.Rin_peak{i_group});
    hold on,
    plot(xbar, 0.5, 'kv', 'markerfacecolor', 'k', 'markersize', 5)
    N = numel(groupdata.Rin_peak{i_group});
    legtext = sprintf('%s, %s, %s, n=%d', plotgroups{i_group,1},plotgroups{i_group, 2},plotgroups{i_group,3}, N);
    hl = legend(legtext, 'Location', 'northeast');
    legend boxoff
    hl.Interpreter = 'none';
    xlabel('Input R (MOhm)')
    ylabel('counts')
    if i_group == 1
        title('Rin (iClamp, peak)')
    end
    
    % CDF summary
    pltidx = sub2ind([2, N_plt_rows], 1, Ngroups+1);
    subplot(N_plt_rows, 2, pltidx), hold on,
    N = histcounts(groupdata.Rin_peak{i_group}, edges);
    N(end+1) = 0;
    cdf_vals = cumsum(N)./sum(N);
    stairs(edges, cdf_vals, '-', 'color', grp_clr)
    xlabel('Input R (MOhm)')
    ylabel('Proportion')
    
    %
    % current clamp asym vals
    %%%%%%%%%%%%%%%%%%%%%%%%
    pltidx = sub2ind([2, N_plt_rows], 2, i_group);
    subplot(N_plt_rows, 2, pltidx)
    
    h = histogram(groupdata.Rin_asym{i_group}, edges);
    h.FaceColor = grp_clr;
    xbar = nanmean(groupdata.Rin_asym{i_group});
    hold on,
    plot(xbar, 0.5, 'kv', 'markerfacecolor', 'k', 'markersize', 5)
    xlabel('Input R (MOhm)')
    ylabel('counts')
    if i_group == 1
        title('Rin (iClamp, asym)')
    end
    
    % CDF summary
    pltidx = sub2ind([2, N_plt_rows], 2, Ngroups+1);
    subplot(N_plt_rows, 2, pltidx), hold on,
    N = histcounts(groupdata.Rin_asym{i_group}, edges);
    N(end+1) = 0;
    cdf_vals = cumsum(N)./sum(N);
    stairs(edges, cdf_vals, '-', 'color', grp_clr)
    xlabel('Input R (MOhm)')
    ylabel('Proportion')
end



%
% CDF plots for +/- 30pA conditions
%
%%%%%%%%%%%%%%%%%%%%%%%%
f=figure;
f.Units = 'Pixels';
f.Position = [182   110   575   648];
Ngroups = (numel(groupdata.Rin_peak)); % adding one for "summary" cdf fig

allNums = cat(1, groupdata.R_pos_30pA{:}, groupdata.R_neg_30pA{:});
edges = linspace(min(allNums)-1, max(allNums)+1, 50);
legtext_pos = {};
legtext_neg = {};
for i_group = 1:numel(groupdata.Rin_asym)
    
    hva_name = plotgroups{i_group, 3};
    grp_clr = hvaPlotColor(hva_name);
    
    % +30 pA
    subplot(2,1,1), hold on,
    N = histcounts(groupdata.R_pos_30pA{i_group}, edges);
    N = [0, N];
    cdf_vals = cumsum(N)./sum(N);
    stairs(edges, cdf_vals, '-', 'color', grp_clr)
    N = numel(groupdata.R_pos_30pA{i_group});
    legtext_pos{i_group} = sprintf('%s, n=%d', plotgroups{i_group,3}, N);
    force_consistent_figs(gca, 'ax');
    title('Rin from +30pA step');
    xlabel('Input R (MOhm)');
    ylabel('Proportion');
    
    % -30 pA
    subplot(2,1,2), hold on,
    N = histcounts(groupdata.R_neg_30pA{i_group}, edges);
    N = [0, N];
    cdf_vals = cumsum(N)./sum(N);
    stairs(edges, cdf_vals, '-', 'color', grp_clr);
    N = numel(groupdata.R_neg_30pA{i_group});
    legtext_neg{i_group} = sprintf('%s, n=%d', plotgroups{i_group,3}, N);
    force_consistent_figs(gca, 'ax');
    title('Rin from -30pA step');
    xlabel('Input R (MOhm)');
    ylabel('Proportion');
end
subplot(2,1,1)
h1 = legend(legtext_pos);
subplot(2,1,2)
h2 = legend(legtext_neg);


%
% bar plot of Rin
%
%%%%%%%%%%%%%%%%%%%%%%%%
f=figure;
hold on,
f.Units = 'Pixels';
f.Position = [182   110   575   648];
Ngroups = (numel(groupdata.Rin_peak));

legtext = {};
xtick_labels = {};
hb = [];
for i_group = 1:numel(groupdata.Rin_asym)
    
    hva_name = plotgroups{i_group, 3};
    grp_clr = hvaPlotColor(hva_name);
    
    % -30 pA
    xmedian = median(groupdata.R_neg_30pA{i_group});
    N = length(groupdata.R_neg_30pA{i_group});
    btstrp_inds = unidrnd(N, [N, 10000]);
    btstrp_data = groupdata.R_neg_30pA{i_group}(btstrp_inds);
    se = std(median(btstrp_data, 1));
    hb(i_group) = bar(i_group, xmedian, 'facecolor', grp_clr, 'edgecolor', grp_clr);
    plot([i_group, i_group], [xmedian, xmedian+se], 'color', grp_clr, 'linewidth', 3)
    
    N = numel(groupdata.R_neg_30pA{i_group});
    legtext{i_group} = sprintf('%s, n=%d', hva_name, N);
    xtick_labels{i_group} = hva_name;
end
set(gca, 'xtick', [1:4], 'xticklabel', xtick_labels)
title('Rin from -30pA step');
ylabel(sprintf('Input R (MOhm)\n(Median +/- SE)'));
legend(hb, legtext)
legend boxoff

% %
% %  box plot of Rin
% %
% %%%%%%%%%%%%%%%%%%%%%
% box_groups = [];
% box_colors = [];
% box_data = [];
% for i_group = 1:numel(groupdata.R_neg_30pA)
%     
%     hva_name = plotgroups{i_group, 3};
%     box_colors = cat(1, box_colors, hvaPlotColor(hva_name));
%     
%     N = numel(groupdata.R_neg_30pA{i_group});
%     box_data = cat(1, box_data, groupdata.R_neg_30pA{i_group});
%     box_groups = cat(1, box_groups, repmat(hva_name, N, 1));
% end
% hf = figure;
% hb = boxplot(box_data, box_groups,...
%              'grouporder', plotgroups(:,3),...
%              'colors', box_colors,...
%              'notch', 'on',...
%              'symbol', '',...
%              'whisker', 0);
% set(gca, 'tickdir', 'out', 'fontsize', 24, 'ylim', [30, 100])
% set(hb, 'linewidth', 3)
% ylabel('Input Resistance')
% legend(gca, plotgroups(:,3)')
% 
% 

% Do an analysis of variance (kruskal wallis)
X = [];
group = [];
for i_group = 1:numel(groupdata.Rin_asym)
    
    hva_name = plotgroups{i_group, 3};
    
    % -30 pA
    X = cat(1, X, groupdata.R_neg_30pA{i_group});
    N = length(groupdata.R_neg_30pA{i_group});
    group = cat(1, group, repmat({hva_name}, N, 1));
end
[~,~,STATS] = kruskalwallis(X, group);
comparison = multcompare(STATS);
p_sig = comparison(:,end) < 0.05;
sig_comps = comparison(p_sig, 1:2);
h_pmtx = figure;
hold on,
for i_sig = 1:size(sig_comps, 1)
   x=sig_comps(i_sig, 1);
   y=sig_comps(i_sig, 2);
   text(x, y, 'X') 
end
ylim([0, 5])
xlim([0, 5])
labels = {plotgroups{:,3}};
set(gca, 'XTickLabel', labels, 'XTick', [1,2,3,4])
set(gca, 'YTickLabel', labels, 'YTick', [1,2,3,4])



