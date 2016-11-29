function l_match = groupMatcher(group_params, cell_params)

    % group_params: {cell_type, layer, brain_area, opsin}
    %               the characteristics of each group. One row for EACH group.
    %
    % cell_params: {cell_type, brain_area, opsin}
    %               the characteristics of the neuron . Only one row
    %
    % l_match: a logical vector where 1 indicates a group that matches the
    % neuron's characteristics. There can be no more than one match.


    % prep the cell_params array. Turn NaN into string('nan') so that regex works
    l_nan = cellfun(@(x) all(isnan(x)), cell_params);
    cell_params(l_nan) = cellfun(@num2str, cell_params(l_nan), 'uniformoutput', false);

    % does the cell_params cell type match any of the group_params options?
    [cell_type, cell_layer] = parse_type_layer(cell_params{1});
    l_celltype_match = subfxn_celltype_match(group_params(:,1), cell_type);
    l_layer_match = subfxn_layer_match(group_params(:,2), cell_layer);
    l_brainarea_match = subfxn_brainarea_match(group_params(:,3), cell_params{2});
    l_opsin_match = subfxn_opsin_match(group_params(:,4), cell_params{3});

    l_match = l_celltype_match & l_layer_match & l_brainarea_match & l_opsin_match;
    assert(sum(l_match)<=1, 'ERROR: found too many group indicies')

end


function l_celltype_match = subfxn_celltype_match(group_celltype, cell_celltype)
    % modify the celltype to include special flags such as 'all_pv', or 'all_som'
    switch lower(cell_celltype)
        case {'pvtom', 'fs'}
            cell_celltype = [cell_celltype, ' all_pv'];
        case {'somtom', 'npvin'}
            cell_celltype = [cell_celltype, ' all_som'];
    end
    l_group_is_any = strcmpi(group_celltype, 'any');
    l_group_matches_input = cellfun(@(x) ~isempty(regexp(cell_celltype, x, 'once')), group_celltype);
    l_celltype_match = l_group_is_any | l_group_matches_input;
    
end

function l_layer_match = subfxn_layer_match(group_layers, cell_layer)
    l_layer_is_any = strcmpi(group_layers, 'any');
    l_layer_matches_input = cellfun(@(x) ~isempty(regexp(cell_layer, x, 'once')), group_layers);
    l_layer_match = l_layer_is_any | l_layer_matches_input;
end

function l_brainarea_match = subfxn_brainarea_match(group_brainarea, cell_brainarea)
    % modify the cell_brainarea to include special flags, such as 'med' and  'lat'
    switch lower(cell_brainarea)
        case {'am', 'pm', 'am/pm', 'pm/am'}
            cell_brainarea = [cell_brainarea, ' med'];
        case {'al', 'lm', 'al/lm', 'lm/al'}
            cell_brainarea = [cell_brainarea, ' lat'];
    end
    l_brainarea_is_any = strcmpi(group_brainarea, 'any');
    l_brainarea_matches_input = cellfun(@(x) ~isempty(regexp(cell_brainarea, x, 'once')), group_brainarea);
    l_brainarea_match = l_brainarea_is_any | l_brainarea_matches_input;
end

function l_opsin_match = subfxn_opsin_match(group_opsins, cell_opsin)
    l_opsin_is_any = strcmpi(group_opsins, 'any');
    l_opsin_matches_input = cellfun(@(x) ~isempty(regexpi(cell_opsin, x)), group_opsins);
    l_opsin_match = l_opsin_matches_input | l_opsin_is_any;
end

function [cell_type, cell_layer] = parse_type_layer(type_layer)
	idx_flag = regexp(type_layer, '_L');
    cell_type = type_layer(1:idx_flag-1);
    cell_layer = type_layer(idx_flag+1:end);
end


