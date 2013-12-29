function tempname = mwpath(fname)
%
% tempname=meshtemppath(fname)
%
% get full temp-file name by prepend working-directory and current session name
%
% author: Qianqian Fang (fangq <at> nmr.mgh.harvard.edu)
%
% input:
%    fname: input, a file name string
%
% output:
%    tempname: output, full file name located in the working directory
%
%    if global variable ISO2MESH_TEMP is set in 'base', it will use it
%    as the working directory; otherwise, will use matlab function tempdir
%    to return a working directory.
%
%    if global variable ISO2MESH_SESSION is set in 'base', it will be
%    prepended for each file name, otherwise, use supplied file name.
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%
    username = ['iso2mesh-' getenv('username')]; % for windows
    tdir = tempdir();
    if tdir(end) ~= filesep, tdir = [tdir filesep]; end
    tdir = [tdir username filesep];
    if ~exist(tdir, 'dir'), mkdir(tdir); end
    tempname = [tdir fname];
end
