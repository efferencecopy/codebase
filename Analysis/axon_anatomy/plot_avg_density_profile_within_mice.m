function plot_avg_density_profile_within_mice(dat, all_hvas, NORMALIZE)

mouse_names = fieldnames(dat);
for mouse = mouse_names'
    mouse = mouse{1};
    hf = figure; hold on,
    legend_text = {};
    for i_hva = 1:numel(all_hvas)
        if isempty(dat.(mouse).(all_hvas{i_hva}).images)
            continue % no data to plot
        end
        legend_text{end+1} = all_hvas{i_hva};
        
        % compute yy_profile average across different lengthed profiles
        [xx_avg, yy_avg] = compute_avg_profile(dat, mouse, all_hvas{i_hva});
        
        if NORMALIZE
            yy_avg = yy_avg ./ max(yy_avg(xx_avg>0));
        end
        plt_clr = hvaPlotColor(all_hvas{i_hva});
        plot(xx_avg, yy_avg, '-', 'color', plt_clr, 'linewidth', 2)
    end
    hf.Name = mouse;
    ylabel(all_hvas{i_hva});
    legend(legend_text, 'location', 'best')
    axis tight
    
    % add the layer boundary information
    y_val_delta = diff(get(gca, 'ylim')) ./ 10;
    y_val_for_layers = max(get(gca, 'ylim')) + y_val_delta;
    bar_width = y_val_delta .* 0.75;
    for i_hva = 1:numel(all_hvas)
        if isempty(dat.(mouse).(all_hvas{i_hva}).images)
            continue % no data to plot
        end
        layers = dat.(mouse).(all_hvas{i_hva}).layers;
        plt_clr = hvaPlotColor(all_hvas{i_hva});
        plot_layer_boundaries(layers, y_val_for_layers, bar_width, plt_clr, hf);
        
        % increment the y_val
        y_val_for_layers = y_val_for_layers + y_val_delta;
    end
    drawnow
    
end











