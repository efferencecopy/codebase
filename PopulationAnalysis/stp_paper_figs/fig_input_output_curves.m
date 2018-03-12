function fig_input_output_curves(dat, plotgroups)

Ngroups = size(plotgroups, 1);
groupdata = [];
groupdata.fi = repmat({[]}, 1, Ngroups);

% % latice one
% bin_edge_l = [20, 140, 290, 590];
% bin_edge_r = [40, 160, 310, 610];
% bin_center = [30, 150, 300, 600];
% 
% % latice two
% bin_edge_l = [90,  220, 335, 465, 710];
% bin_edge_r = [110, 240, 355, 485, 740];
% bin_center = [100, 230, 345, 470, 725];
 
% both latices
bin_edge_l = [20, 90,  140, 220, 290, 335, 465, 590, 710];
bin_edge_r = [40, 110, 160, 240, 310, 355, 485, 610, 740];
bin_center = [30, 100, 150, 230, 300, 345, 470, 600, 725];

for i_ex = 1:numel(dat)
    
    for i_ch = 1:2
        
        % check to make sure this neuron was defined
        hs_name = sprintf('dcsteps_hs%d', i_ch);
        isvalid = ~strncmp(dat{i_ex}.info.fid.(hs_name)(end-8:end), [filesep, 'none.abf'], 9);
        if ~isvalid
            continue
        end
        
        % check the attributes
        ch_attribs = {dat{i_ex}.info.cellType{i_ch}, dat{i_ex}.info.brainArea, dat{i_ex}.info.opsin};
        group_idx = groupMatcher(plotgroups, ch_attribs);
        if sum(group_idx) == 0; continue; end
        
        % pull out the values that correspond to the pre-defined Icmd bins
        pa = dat{i_ex}.dcsteps.fi_curve{i_ch}(:,1);
        fr = dat{i_ex}.dcsteps.fi_curve{i_ch}(:,2);
        l_icmd = nan(size(bin_center));
        bined_fr = nan(size(bin_center));
        for i_cmd = 1:numel(bin_center)
            tmp_idx = (pa > bin_edge_l(i_cmd)) & (pa <= bin_edge_r(i_cmd));
            assert(sum(tmp_idx) <=1, 'ERROR: found too many matches')
            if any(tmp_idx)
                bined_fr(i_cmd) = fr(find(tmp_idx==1));
            end
        end
        
        if sum(~isnan(bined_fr)) >= 3
            % add the statistics to the appropriate cell array
            groupdata.fi{group_idx} = cat(1, groupdata.fi{group_idx}, {bined_fr});
        else
            continue
        end

    end
end




% -------  plot the spike rates as scatter plot  -----------
figure, hold on,
legtext = {};
leghands = [];
group_all_frs = {};
group_all_pa = {};
for i_group = 1:Ngroups
    tmp_dat = groupdata.fi{i_group};
    group_all_frs{i_group} = cat(1, tmp_dat{:});
end

for i_group = 1:Ngroups
    hva_name = plotgroups{i_group, 3};
    grp_clr = hvaPlotColor(hva_name);
    
    X = bin_center;
    Y = group_all_frs{i_group};
    
    hp = my_errorbar(X, nanmean(Y,1), stderr(Y, 1), 's',...
                                                    'color', grp_clr,...
                                                    'linewidth', 2,...
                                                    'markersize', 8,...
                                                    'markeredgecolor', grp_clr,...
                                                    'markerfacecolor', grp_clr);
    
    designMtx = bin_center(:);
    response = nanmean(Y, 1)';
    b0 = [0, 200, .001];
    [params, liklihood] = halfSquareFit(designMtx, response, b0, 'simple');
    fit = params(1)+params(3)*(max(bin_center-params(2),0).^2);
    plot(X, fit, '-', 'color', grp_clr, 'linewidth', 2)
    
    
    leghands(i_group) = hp;
    N = size(Y, 1);
    legtext{i_group} =  sprintf('%s, %s, %s, N=%d', plotgroups{i_group,1},plotgroups{i_group, 2},plotgroups{i_group,3}, N );
    
    
end
force_consistent_figs(gca, 'ax');
xlabel('Current injection (pA)')
ylabel('Avg Spike Rate (Hz)')
legend(leghands, legtext, 'location', 'northwest')
legend boxoff


