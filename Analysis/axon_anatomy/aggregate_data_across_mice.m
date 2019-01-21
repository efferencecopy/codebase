function hva_aggs = aggregate_data_across_mice(dat, all_hvas, NORMALIZE, PROFILE_TYPE)

switch PROFILE_TYPE
    case 'baselined'
        ptype = 'baselined_profiles';
        yy_ptype = 'yy_baselined';
    case 'average'
        ptype = 'avg_profiles';
        yy_ptype = 'yy_avg';
end

% make an empty container for the aggregate data
hva_aggs = struct();

% loop over mice and hvas, save data to the agg container
mouse_names = fieldnames(dat);
for i_mouse = 1:numel(mouse_names)
    mouse = mouse_names{i_mouse};
    layer_line = dat.(mouse).layer_line;
    for i_hva = 1:numel(all_hvas)
        
        % init the output containers
        if ~isfield(hva_aggs, layer_line)
            hva_aggs.(layer_line) = struct();
        end
        if ~isfield(hva_aggs.(layer_line), all_hvas{i_hva})
            hva_aggs.(layer_line).(all_hvas{i_hva}).xx_vals = {};
            hva_aggs.(layer_line).(all_hvas{i_hva}).yy_vals = {};
        end
        if ~isfield(hva_aggs.(layer_line).(all_hvas{i_hva}), 'f_per_pix')
            hva_aggs.(layer_line).(all_hvas{i_hva}).f_per_pix = [];
        end
        
        % add NaN values to the dataset if real data are not present
        if isempty(dat.(mouse).(all_hvas{i_hva}).images)
            hva_aggs.(layer_line).(all_hvas{i_hva}).f_per_pix = cat(...
                3,...
                hva_aggs.(layer_line).(all_hvas{i_hva}).f_per_pix,...
                nan(5,1));
            continue
        end
        
        % grab the pre-computed y-profile, which is the average profile
        % across 4 images (within a mouse)
        if isfield(dat.(mouse), ptype)
            xx_avg = dat.(mouse).(ptype).(all_hvas{i_hva}).xx_common;
            yy_avg = dat.(mouse).(ptype).(all_hvas{i_hva}).(yy_ptype);
        else
            % if PROFILE_TYPE is 'average', and it doesn't currently exist,
            % we can just call compute_avg_profile and use the return args.
            % But if we want 'baseline' profiles we should throw an error
            % b/c the return of compute_avg_profile is 'average', not
            % 'baseline'
            if strcmpi(PROFILE_TYPE, 'average')
                [xx_avg, yy_avg] = compute_avg_profile(dat, mouse, all_hvas{i_hva});
            else
                error('baseline metrics have not been computed')
            end
        end
        
        switch NORMALIZE
            case 'max'
                L23_xx = dat.(mouse).(all_hvas{i_hva}).layer_avg(2);
                l_L23 = xx_avg>0 & xx_avg<L23_xx;
                yy_avg = yy_avg ./ max(yy_avg(l_L23));
            case 'integral'
                l_ctx = xx_avg>-75 & xx_avg<850; % roughly all of cortex
                yy_avg = yy_avg ./ nansum(yy_avg(l_ctx));
            case 'none'
                % do nothing to the yy_avg
            otherwise
                error('unknown normalize method')
        end
        
        % add these data to the aggregated dataset
        hva_aggs.(layer_line).(all_hvas{i_hva}).xx_vals = cat(1, hva_aggs.(layer_line).(all_hvas{i_hva}).xx_vals, xx_avg);
        hva_aggs.(layer_line).(all_hvas{i_hva}).yy_vals = cat(1, hva_aggs.(layer_line).(all_hvas{i_hva}).yy_vals, yy_avg);
        
        
        % compute the absolute fluroescence in each layer
        layer_bounds = dat.(mouse).(all_hvas{i_hva}).layer_avg;
        layer_bounds = cat(1, -75, layer_bounds); % so that L1 can get calculated too

        n_layers = numel(layer_bounds);
        f_per_pix = nan(n_layers-1, 1);
        for i_layer = 1:n_layers-1
            l_in_layer = (xx_avg > layer_bounds(i_layer)) & (xx_avg < layer_bounds(i_layer+1));
            f = sum(yy_avg(l_in_layer));
            img_width = size(dat.(mouse).(all_hvas{i_hva}).images{1}, 2);
            area_pix = sum(l_in_layer) * img_width;
            f_per_pix(i_layer) = f ./ area_pix;
        end
        % add these data to the output array
        hva_aggs.(layer_line).(all_hvas{i_hva}).f_per_pix = cat(...
                    3,...
                    hva_aggs.(layer_line).(all_hvas{i_hva}).f_per_pix,...
                    f_per_pix);
        
    end
end


% adjust the xx_avg and yy_avg arrays so that they are the same dimensions
mouse_lines = fieldnames(hva_aggs);
for i_mline = 1:numel(mouse_lines)
    for i_hva = 1:numel(all_hvas)
        % check for data
        if ~isfield(hva_aggs.(mouse_lines{i_mline}), all_hvas{i_hva})
            continue
        end
        
        xx_raw = hva_aggs.(mouse_lines{i_mline}).(all_hvas{i_hva}).xx_vals;
        yy_raw = hva_aggs.(mouse_lines{i_mline}).(all_hvas{i_hva}).yy_vals;
        min_x = min(cellfun(@(x) min(x), xx_raw));
        max_x = max(cellfun(@(x) max(x), xx_raw));
        xx_common = min_x : max_x;
        yy_common = pad_for_agg_sum(xx_common, xx_raw, yy_raw);
        yy_avg_across_mice = nanmean(yy_common, 1);
        
        % add the agg_avg to the array
        if ~isfield(hva_aggs.(mouse_lines{i_mline}), 'avg_profiles')
            hva_aggs.(mouse_lines{i_mline}).avg_profiles = struct();
        end
        hva_aggs.(mouse_lines{i_mline}).avg_profiles.(all_hvas{i_hva})= struct(...
            'xx_common', xx_common, 'yy_avg', yy_avg_across_mice, 'N', size(yy_raw, 1));
    end
end

end




function yy_common = pad_for_agg_sum(xx_common, xx_raw, yy_raw)

% since the outputs will have the same number of samples (with padding), I
% can output a matrix and not a cell array
n_mice = numel(xx_raw);
n_xx = numel(xx_common);

% define the output variable. fill it up with padding (nans)
yy_common = nan(n_mice, n_xx);

% iterate over the yy_raw profiles and put them into the yy_common output,
% keeping in mind how many pads they should have on the front and back.
for i_img = 1:n_mice
    idx_start = xx_raw{i_img}(1) - xx_common(1) + 1;
    idx_stop = idx_start + numel(yy_raw{i_img}) - 1;
    yy_common(i_img, idx_start : idx_stop) = yy_raw{i_img};
end

end



