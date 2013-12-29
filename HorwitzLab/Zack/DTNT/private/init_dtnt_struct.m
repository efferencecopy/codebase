% This function creates an empty, `nelements`-long trial specs structure. By
% default, this function resets the globally-defined counters. Pass `false` as
% the second argument to prevent this.

function ts = init_dtnt_struct(nelements, reset_counters)
global gl

if nargin < 2 || isempty(reset_counters)
    reset_counters = true;
end

ts_fields = {'colordir' 'predictedthreshold' 'parentvertices' 'parentOOGs' 'measuredthreshold' 'OOG'};
ts = cell2struct(cell(numel(ts_fields),nelements), ts_fields);

if reset_counters
    gl.dtnt.nspecs = 0;
    gl.dtnt.nfiles = 0;
    gl.dtnt.filenames = {};
end
