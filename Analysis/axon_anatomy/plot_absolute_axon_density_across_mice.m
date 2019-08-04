function all_norm_dat = plot_absolute_axon_density_across_mice(hva_aggs, all_hvas, norm_method)
% make a figure that has N_lines x N_layers cells. Each element
% represents the avg_unnormed fluorescence intensity in that layer.
% plot the results as imagesc (norm and un-norm to L23)
%
% will need to keep better track of layer boundaries:
%  * keep track of fluorescence/area for each mouse, present the avg of
%     these metrics. This will account for differences in area across mice
%     and lines


mouse_lines = fieldnames(hva_aggs);
n_lines = numel(mouse_lines);
all_norm_dat = struct();
abs_f_dat = struct(); % init a container to fill up
for i_mline = 1:n_lines
   
    abs_f_dat.(mouse_lines{i_mline}) = []; % hard coding layer count
    n_mice = -1; % a counter updated later
    for i_hva = 1:numel(all_hvas)
        
        % pull out the values for all mice and average across mice
        all_f_per_pix = hva_aggs.(mouse_lines{i_mline}).(all_hvas{i_hva}).f_per_pix;
        abs_f_dat.(mouse_lines{i_mline})(:, i_hva, :) = all_f_per_pix;
        if i_hva == 1
            n_mice = size(all_f_per_pix, 3);
        else
            assert(size(all_f_per_pix, 3) == n_mice, 'ERROR: n-mice differs across HVAs')
        end
        
    end
end

%
% a heat map for each layer_line
%%%%%%%%%%%%%%%%%%%
hf = figure;
hf.Position = [36         513        1361         236];
max_abs_val = -Inf;
min_abs_val = Inf;
for i_mline = 1:n_lines
    subplot(1, n_lines, i_mline)
    raw_dat = abs_f_dat.(mouse_lines{i_mline});
    assert(size(raw_dat, 2) == 4, 'ERROR: some HVAs are missing')
    
    switch norm_method
        case {'lm', 'al', 'pm', 'am'}
            % norm within a mouse to L2/3 of LM
            norm_idx_hva = strcmpi(all_hvas, norm_method);
            norm_dat = nan(size(raw_dat));
            for i_mouse = 1:size(raw_dat, 3)
                norm_dat(:,:,i_mouse) = raw_dat(:,:,i_mouse) ./ raw_dat(2, norm_idx_hva, i_mouse);
            end
        case 'self'
            norm_dat = nan(size(raw_dat));
            for i_mouse = 1:size(raw_dat, 3)
                norm_dat(:,:,i_mouse) = bsxfun(@rdivide, raw_dat(:,:,i_mouse), raw_dat(2, :, i_mouse));
            end
    end
    
    % remove mice that aren't sampled in all areas
    l_good_mice = false(1,1,size(raw_dat, 3));
    for i_mouse = 1:size(raw_dat, 3)
        l_good_mice(1,1,i_mouse) = ~any(isnan(reshape(norm_dat(:,:,i_mouse), 1, [])));
    end
    norm_dat = norm_dat(:,:,l_good_mice);
    all_norm_dat.(mouse_lines{i_mline}) = norm_dat; % raw_dat(:,:,l_good_mice);
    avg_img = mean(norm_dat, 3);
    
    imagesc(avg_img);
    colorbar
    title(sprintf('F/pix relative to L2/3 in %s\n %s, n=%d', upper(norm_method), mouse_lines{i_mline}, size(norm_dat, 3)))
    layer_labels = {'L1', 'L2/3', 'L4', 'L5', 'L6'};
    set(gca, 'ytick', [1,2,3,4,5], 'yticklabels', layer_labels);
    set(gca, 'xtick', [1,2,3,4], 'xticklabels', all_hvas, 'tickdir', 'out');
    
    % keep track of cmap limits
    max_abs_val = max(max_abs_val, max(avg_img(:)));
    min_abs_val = min(min_abs_val, min(avg_img(:)));
end

for i_plt = 1:n_lines
    subplot(1, n_lines, i_plt)
    caxis([min_abs_val, max_abs_val])
end











