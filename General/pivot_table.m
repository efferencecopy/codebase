function [pivot_table] = pivot_table(t, aggregate_by, ...
        data, function_handle, varargin)
% PIVOT_TABLE   create a pivot table (similar to that created by Excel), using
% an input table, one or more columns to aggregate by, one or more data columns,
% and one or more functions to apply to the aggregated data.
%
%   EXAMPLE
% p = pivot_table(t, aggregate_by, data, function_handle)
% p = pivot_table(t, aggregate_by, data, function_handle, weight_by)
% p = pivot_table(t, aggregate_by, data, function_handle, weight_by, 'StableOrder', true)
%
%   INPUTS
%   t: a table to create the pivot table from
%   aggregate_by: either a string containing the table column name in t to aggregate over, or a cellstr of multiple table columns.
%   data: string containing name of table columns in t to use for data, or a cellstr of multiple table variables.
%   function_handle: a function handle to apply to the aggregated data
%
%   weight_by (optional): a string containing the table column name to weight the aveage

%   OPTIONS
%   StableOrder (default false)
%   should the output aggregated values be sorted or not

    p = inputParser;
    p.addRequired('t', @(x) istable(x));
    p.addRequired('aggregate_by', @(x) ischar(x) || iscellstr(x));
    p.addRequired('data', @(x) ischar(x) || iscellstr(x));
    p.addRequired('function_handle', @(x) isa(function_handle, 'function_handle'));
    p.addOptional('weight_by', [], @(x) ischar(x) || iscellstr(x) || isempty(x));
    p.addParameter('StableOrder', false, @(x) isa(x, 'logical'));
    p.parse(t, aggregate_by, data, function_handle, varargin{:});

    weight_by = p.Results.weight_by;
    use_stable_order = p.Results.StableOrder;

    pivot_table = table;

    if use_stable_order
        [unique_values, ~, ic] = unique(t(:, aggregate_by), 'stable');
    else
        [unique_values, ~, ic] = unique(t(:, aggregate_by));
    end
    pivot_table(:, aggregate_by) = unique_values;

    % handle single or multiple output columns
    has_multiple_data_columns = iscell(data);
    if has_multiple_data_columns
        columns = data;
    else
        columns = {data};
    end

    % aggregate the data and apply the function_handle
    for col_i = 1:length(columns)
        if isempty(weight_by)
            output_variable{col_i} = sprintf('%s_of_%s', func2str(function_handle), columns{col_i});
            x = accumarray(ic, t{:, columns{col_i}}, [], function_handle);
        else
            output_variable{col_i} = sprintf('w_mean_of_%s', columns{col_i});
            if isequal(@mean, function_handle)
                y = accumarray(ic, t{:, columns{col_i}} .* t{:, weight_by}, [], @sum);
                z = accumarray(ic, t{:, weight_by}, [], @sum);
                x = y ./ z;
            else
                error('weight_by only valid for @mean aggregation function');
            end
        end

        % Some function names can not be turned into strings, especially if
        % an anonymous function was passed as the function_handle.
        % In this case, use f_1, f_2 etc. for the output columns names
        try
            pivot_table{:, output_variable{col_i}} = x;
        catch
            pivot_table{:, sprintf('f_%d', col_i)} = x;
        end
    end
end

