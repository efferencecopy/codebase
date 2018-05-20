function dat = get_average_layer_boundaries(dat, all_hvas)

all_mice = fieldnames(dat);
for i_mouse = 1:numel(all_mice)
    for i_hva = 1:numel(all_hvas)
        if isempty(dat.(all_mice{i_mouse}).(all_hvas{i_hva}).images)
            continue % no data to plot
        end
        
        % get layers for the 4 representative images for a specific HVA
        layers = dat.(all_mice{i_mouse}).(all_hvas{i_hva}).layers;
        layer_names = fieldnames(layers);
        
        % need to loop over images and measure distance to L1 in that img
        n_img = numel(dat.(all_mice{i_mouse}).(all_hvas{i_hva}).images);
        n_layers = numel(layer_names);
        depth_relative_to_L1 = nan(n_layers, n_img);
        for i_img = 1:n_img
            for i_layer = 1:n_layers
                val = mean(layers.(layer_names{i_layer}){i_img});  % avg Y-pos for a single layer/img
                L1 = mean(layers.L1_L23{i_img});
                depth_relative_to_L1(i_layer, i_img) = val - L1;
            end
        end
        
        % store the avg layer boundaries (from L1)
        dat.(all_mice{i_mouse}).(all_hvas{i_hva}).layer_avg = mean(depth_relative_to_L1, 2);
        
    end
    
end


end
