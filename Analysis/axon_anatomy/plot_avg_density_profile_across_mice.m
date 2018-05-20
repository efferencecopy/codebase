function plot_avg_density_profile_across_mice(hva_aggs, all_hvas)

% make one figure for each mouse_layer_line
hf = figure();
mouse_lines = fieldnames(hva_aggs);
for i_mline = 1:numel(mouse_lines)
    subplot(1, numel(mouse_lines), i_mline)
    hold on,
    legend_text = {};
    max_x = 0;
    max_y = 0;
    for i_hva = 1:numel(all_hvas)
        % check for data
        if ~isfield(hva_aggs.(mouse_lines{i_mline}), all_hvas{i_hva})
            continue
        end
        
        plt_clr = hvaPlotColor(all_hvas{i_hva});
        N = hva_aggs.(mouse_lines{i_mline}).avg_profiles.(all_hvas{i_hva}).N;
        legend_text{end+1} = sprintf('%s, N=%d', upper(all_hvas{i_hva}), N);
        xx = hva_aggs.(mouse_lines{i_mline}).avg_profiles.(all_hvas{i_hva}).xx_common;
        yy = hva_aggs.(mouse_lines{i_mline}).avg_profiles.(all_hvas{i_hva}).yy_avg;
        plot(xx, yy, '-', 'linewidth', 2, 'color', plt_clr)
        max_x = max(max_x, max(xx));
        max_y = max(max_y, max(yy));
    end
    legend(legend_text, 'location', 'best')
    xlim([-125, max_x])
    ylim([0, max_y])
    xlabel('depth from L1/L23 boundary')
    ylabel('fluorescence')
    title(sprintf('Average for mouse type: %s', mouse_lines{i_mline}))
    drawnow
end

end







