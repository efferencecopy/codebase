function plot_density_profiles_from_all_slices(dat, all_hvas)

mouse_names = fieldnames(dat);
for i_mouse = 1:numel(mouse_names)
    mouse = mouse_names{i_mouse};
    
    hf = figure;
    hf.Name = sprintf('%s, %s', mouse, dat.(mouse).layer_line);
    hf.Units = 'normalized';
    hf.Position = [0.3049   -0.1722    0.3618    1.0578];
    max_y = 0;
    max_x = 0;
    min_x = inf;
    ha = [];
    for i_hva = 1:numel(all_hvas)
        hva = all_hvas{i_hva};
        n_imgs = numel(dat.(mouse).(hva).profiles);
        ha(i_hva) = subplot(numel(all_hvas),1,i_hva); hold on,
        ylabel(hva);
        for i_img = 1:n_imgs
            
            % raw profile
            profile_raw = dat.(mouse).(hva).profiles{i_img};
            
            % layer boundaries
            l1_yy = dat.(mouse).(hva).layers.L1_L23{i_img};
            l1_idx = round(mean(l1_yy));
            
            % plot, but define x=0 as the L1/L23 boundary
            xx = [1:numel(profile_raw)] - l1_idx;
            hp = plot(xx, profile_raw);
            
%             % add a dashed line denoting the bottom of the cortex
%             L5_L6_yy = dat.(mouse).(hva).layers.L5_L6{i_img};
%             L5_L6_idx = mean(L5_L6_yy) - l1_idx;
%             plot([L5_L6_idx, L5_L6_idx], [min(profile_raw), max(profile_raw)], '--k')
            
            % store the max axis vals
            max_y = max(max_y, max(profile_raw));
            max_x = max(max_x, xx(end));
            min_x = min(min_x, xx(1));
            
        end
    end
    set(ha, 'YLim', [0, max_y], 'XLim', [min_x, max_x])
    
end


