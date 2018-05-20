function plot_raw_images_with_layer_boundaries(dat, all_hvas)

mouse_names = fieldnames(dat);
for mouse = mouse_names'
    mouse = mouse{1};
    
    nrows = 0;
    for i_hva = 1:numel(all_hvas)
        nrows = max(nrows, numel(dat.(mouse).(all_hvas{i_hva}).images));
    end
    
    hf = figure;
    hf.Name = sprintf('%s, %s', mouse, dat.(mouse).layer_line);
    hf.Units = 'normalized';
    hf.Position = [0.1542    0.0089    0.3042    0.9033];
    for i_hva = 1:numel(all_hvas)
        n_imgs = numel(dat.(mouse).(all_hvas{i_hva}).images);
        for i_img = 1:n_imgs
            plt_idx = sub2ind([nrows, numel(all_hvas)], i_img, i_hva);
            ha = subplot(numel(all_hvas), nrows,plt_idx);
            img_raw = dat.(mouse).(all_hvas{i_hva}).images{i_img};
            hi = image(max(img_raw, [], 3));
            colormap((gray(256)))
            ha.XTick = [];
            ha.YTick = [];
            hold on,
            l1_yy = dat.(mouse).(all_hvas{i_hva}).layers.L1_L23{i_img};
            l23_yy = dat.(mouse).(all_hvas{i_hva}).layers.L23_L4{i_img};
            l4_yy = dat.(mouse).(all_hvas{i_hva}).layers.L4_L5{i_img};
            l5_yy = dat.(mouse).(all_hvas{i_hva}).layers.L5_L6{i_img};
            lbot_yy = dat.(mouse).(all_hvas{i_hva}).layers.bottom{i_img};
            plot([1:size(img_raw,2)], l1_yy, 'w--')
            plot([1:size(img_raw,2)], l23_yy, 'c--')
            plot([1:size(img_raw,2)], l4_yy, 'b--')
            plot([1:size(img_raw,2)], l5_yy, 'g--')
            plot([1:size(img_raw,2)], lbot_yy, 'r--')
            if i_img == 1
                ha.YLabel.String = all_hvas{i_hva};
            end
        end
    end    
end

