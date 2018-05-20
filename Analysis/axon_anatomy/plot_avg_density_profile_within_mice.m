function plot_avg_density_profile_within_mice(dat, all_hvas, NORMALIZE, PROFILE_TYPE)

switch PROFILE_TYPE
    case 'baselined'
        ptype = 'baselined_profiles';
        yy_ptype = 'yy_baselined';
    case 'average'
        ptype = 'avg_profiles';
        yy_ptype = 'yy_avg';
end


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
        
        % grab the pre-computed y-profile
        if isfield(dat.(mouse), ptype)
            xx_avg = dat.(mouse).(ptype).(all_hvas{i_hva}).xx_common;
            yy_avg = dat.(mouse).(ptype).(all_hvas{i_hva}).(yy_ptype);
        else
            [xx_avg, yy_avg] = compute_avg_profile(dat, mouse, all_hvas{i_hva});
        end
        
        switch NORMALIZE
            case 'max'
                l_L23 = xx_avg>0 & xx_avg<400;
                yy_avg = yy_avg ./ max(yy_avg(l_L23));
            case 'integral'
                l_ctx = xx_avg>-100 & xx_avg<850;
                yy_avg = yy_avg ./ nansum(yy_avg(l_ctx));
        end
        plt_clr = hvaPlotColor(all_hvas{i_hva});
        plot(xx_avg, yy_avg, '-', 'color', plt_clr, 'linewidth', 2)
    end
    hf.Name = sprintf('%s, %s', mouse, dat.(mouse).layer_line);
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











