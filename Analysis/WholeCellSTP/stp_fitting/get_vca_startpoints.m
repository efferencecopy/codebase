function custpts = get_vca_startpoints(lower_bounds, upper_bounds, pts_per_dim)
%get_vca_startpoints 
%
%   EXAMPLE custpts = get_vca_startpoints(bounds, pts_per_dim)
%
%   return custom start points for MultiStart solver. The custom start
%   points will be uniformly distributed across the parameter space. None
%   of the start points will be positioned exactly at the UB/LB of the
%   parameter. 

% linearly space the start points along each dimension
N_dims = numel(upper_bounds);
start_points_1d = {};
for i_dim = 1:N_dims
    tmp_pts = linspace(lower_bounds(i_dim), upper_bounds(i_dim), pts_per_dim+2); % adding 2 so that the start points are not on the bounds
    start_points_1d{i_dim} = tmp_pts(2:end-1); % hacking off the first/last points b/c they are on the bounds.
end

% make an ndgrid
[grid{1:N_dims}] = ndgrid(start_points_1d{:});

% convert each grid matrix to a column vector, concatenate into a
% CustomStartPointSet
grid = cellfun(@(x) x(:), grid, 'uniformoutput', false);
grid = cat(2, grid{:});
custpts = CustomStartPointSet(grid);

end

