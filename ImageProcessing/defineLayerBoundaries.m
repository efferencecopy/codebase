function defineLayerBoundaries

% game plan
%
% Make gui with some buttons
%   * define L1 L2/3 line
%   * define L2/3 L4 line
%   * define L4 L5/6 line
%   * define bottom of L6
%   * change color chanel
%   * export data
%   * scroll through focal planes

[raw, mask] = io_importData;


end



function [raw, mask] = io_importData

[fnames, fpaths] = uigetfile({['*'], ['all files']}, 'Select the raw image and mask', 'MultiSelect', 'on');

% make sure there were only two selections
assert(numel(fnames)==2, 'ERROR: Please select exactly 2 files')


% the raw image will be a .tiff and the mask will be a .mat
idx_mat = cellfun(@regexpi, fnames, repmat({'.mat'}, size(fnames)), 'uniformoutput', false);
idx_mat = cellfun(@(x) ~isempty(x), idx_mat);

idx_tiff = cellfun(@regexpi, fnames, repmat({'.tif'}, size(fnames)), 'uniformoutput', false);
idx_tiff = cellfun(@(x) ~isempty(x), idx_tiff);


% load the mask
load([fpaths, fnames{idx_mat}]); % loads a structure called "mask"


end