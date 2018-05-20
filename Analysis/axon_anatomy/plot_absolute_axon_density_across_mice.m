function plot_absolute_axon_density_across_mice(hva_aggs, all_hvas)
% make a figure that has N_lines x N_layers cells. Each element
% represents the avg_unnormed fluorescence intensity in that layer.
% plot the results as imagesc (norm and un-norm to L23)
%
% will need to keep better track of layer boundaries:
%  * keep track of fluorescence/area for each mouse, present the avg of
%     these metrics. This will account for differences in area across mice
%     and lines


all_hvas = {'lm', 'pm', 'al', 'am'}; % hard coding so that plot order is specified
mouse_lines = fieldnames(hva_aggs);
n_lines = numel(mouse_lines);
abs_f_dat = struct(); % init a container to fill up
for i_mline = 1:n_lines
   
    n_hvas = numel(all_hvas);
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
    
    % norm within a mouse to L2/3 of LM
    norm_idx_hva = strcmpi(all_hvas, 'lm');
    norm_dat = nan(size(raw_dat));
    for i_mouse = 1:size(raw_dat, 3)
        norm_dat(:,:,i_mouse) = raw_dat(:,:,i_mouse) ./ raw_dat(2, norm_idx_hva, i_mouse);
    end
    
    % remove mice that aren't sampled in all areas
    l_good_mice = false(1,1,size(raw_dat, 3));
    for i_mouse = 1:size(raw_dat, 3)
        l_good_mice(1,1,i_mouse) = ~any(isnan(reshape(norm_dat(:,:,i_mouse), 1, [])));
    end
    norm_dat = norm_dat(:,:,l_good_mice);
    avg_img = nanmean(norm_dat, 3);
    
    imagesc(avg_img);
    colorbar
    title(sprintf('F/pix relative to L2/3 in LM\n %s, n=%d', mouse_lines{i_mline}, size(norm_dat, 3)))
    layer_labels = {'L1', 'L2/3', 'L4', 'L5', 'L6'};
    set(gca, 'ytick', [1,2,3,4,5], 'yticklabels', layer_labels);
    set(gca, 'xtick', [1,2,3,4], 'xticklabels', all_hvas);
    
    % keep track of cmap limits
    max_abs_val = max(max_abs_val, max(avg_img(:)));
    min_abs_val = min(min_abs_val, min(avg_img(:)));
end

for i_plt = 1:n_lines
    subplot(1, n_lines, i_plt)
    caxis([min_abs_val, max_abs_val])
end



% %
% % line plots of absolute F for L2/3 only, across lines
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% error('need to collect nan values so that I can normalize within mouse
% before average')
% hf = figure;
% hold on,
% legend_text = {};
% for i_mline = 1:n_lines
%     layer_23 = abs_f_dat.(mouse_lines{i_mline})(2,:);
%     assert(numel(layer_23) == 4, 'ERROR: some HVAs are missing')
%     plot(1:numel(all_hvas), layer_23, '-', 'linewidth', 2)
%     legend_text{end+1} = mouse_lines{i_mline};
% end
% set(gca, 'xtick', 1:numel(all_hvas), 'xticklabels', all_hvas);
% xlim([0.5, numel(all_hvas)+0.5])
% title('F/pix for L2/3 across areas and genotypes')
% ylabel('F/pix for L2/3')
% legend(legend_text, 'location', 'best')








