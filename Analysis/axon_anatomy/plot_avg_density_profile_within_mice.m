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
        xx_profile = {};
        yy_profile = {};
        common_right_edge = inf;
        common_left_edge = -inf;
        n_imgs = numel(dat.(mouse).(all_hvas{i_hva}).images);
        for i_img = 1:n_imgs
            % raw profile
            yy_profile{i_img} = dat.(mouse).(all_hvas{i_hva}).profiles{i_img};
            
            % layer boundaries
            l1_yy = dat.(mouse).(all_hvas{i_hva}).layers.L1_L23{i_img};
            l1_idx = round(mean(l1_yy));

            % plot, but define x=0 as the L1/L23 boundary
            xx_profile{i_img} = [1:numel(yy_profile{i_img})] - l1_idx;
            
            common_left_edge = max(common_left_edge, xx_profile{i_img}(1));
            common_right_edge = min(common_right_edge, xx_profile{i_img}(end));
            
        end
        
        
        new_yy = cellfun(@(x,y) x(find(y==common_left_edge) : find(y==common_right_edge)), yy_profile, xx_profile, 'uniformoutput', false);
        new_yy = cat(2, new_yy{:});
        avg_yy = mean(new_yy, 2);
        new_xx = [common_left_edge : common_right_edge];
        if NORMALIZE
            avg_yy = avg_yy ./ max(avg_yy(new_xx>0));
        end
        plt_clr = hvaPlotColor(all_hvas{i_hva});
        plot(new_xx, avg_yy, '-', 'color', plt_clr, 'linewidth', 2)
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


