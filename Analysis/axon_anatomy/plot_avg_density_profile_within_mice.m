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
        n_imgs = numel(dat.(mouse).(all_hvas{i_hva}).images);
        xx_raw = {};
        yy_raw = {};
        for i_img = 1:n_imgs
            % raw profile
            yy_raw{i_img} = dat.(mouse).(all_hvas{i_hva}).profiles{i_img};
            
            % layer boundaries
            l1_yy = dat.(mouse).(all_hvas{i_hva}).layers.L1_L23{i_img};
            l1_idx = round(mean(l1_yy));
            
            % now, pixel values are relative to the L1 boundary
            xx_raw{i_img} = [1:numel(yy_raw{i_img})] - l1_idx;
        end
        xx_min = min(cellfun(@(x) x(1), xx_raw));
        xx_max = max(cellfun(@(x) x(end), xx_raw));
        
        % aggregate the profiles in a subfunction
        xx_common = xx_min : xx_max;
        yy_common = pad_for_agg_sum(xx_common, xx_raw, yy_raw);
        yy_avg = nanmean(yy_common);
        
        if NORMALIZE
            yy_avg = yy_avg ./ max(yy_avg(xx_common>0));
        end
        plt_clr = hvaPlotColor(all_hvas{i_hva});
        plot(xx_common, yy_avg, '-', 'color', plt_clr, 'linewidth', 2)
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
end

function yy_common = pad_for_agg_sum(xx_common, xx_raw, yy_raw)

% since the outputs will have the same number of samples (with padding), I
% can output a matrix and not a cell array
n_imgs = numel(xx_raw);
n_xx = numel(xx_common);

% define the output variable. fill it up with padding (nans)
yy_common = nan(n_imgs, n_xx);

% iterate over the yy_raw profiles and put them into the yy_common output,
% keeping in mind how many pads they should have on the front and back.
for i_img = 1:n_imgs
    idx_start = xx_raw{i_img}(1) - xx_common(1) + 1;
    idx_stop = idx_start + numel(yy_raw{i_img}) - 1;
    yy_common(i_img, idx_start : idx_stop) = yy_raw{i_img};
end

end








