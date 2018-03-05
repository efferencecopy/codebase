function  [xx_avg, yy_avg] = compute_avg_profile(dat, mouse_name, hva_name)

n_imgs = numel(dat.(mouse_name).(hva_name).images);
xx_raw = {};
yy_raw = {};
for i_img = 1:n_imgs
    % raw profile
    yy_raw{i_img} = dat.(mouse_name).(hva_name).profiles{i_img};
    
    % layer boundaries
    l1_yy = dat.(mouse_name).(hva_name).layers.L1_L23{i_img};
    l1_idx = round(mean(l1_yy));
    
    % now, pixel values are relative to the L1 boundary
    xx_raw{i_img} = [1:numel(yy_raw{i_img})] - l1_idx;
end
xx_min = min(cellfun(@(x) x(1), xx_raw));
xx_max = max(cellfun(@(x) x(end), xx_raw));

% aggregate the profiles in a subfunction
xx_avg = xx_min : xx_max;
yy_common = pad_for_agg_sum(xx_avg, xx_raw, yy_raw);
yy_avg = nanmean(yy_common);

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