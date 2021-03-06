function saveoff(v, f, fname)
%
% saveoff(v,f,fname)
%
% save a surface mesh to Geomview Object File Format (OFF)
%
% author: Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
% date: 2007/03/28
%
% input:
%      v: input, surface node list, dimension (nn,3)
%      f: input, surface face element list, dimension (be,3)
%      fname: output file name
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%
    fid = fopen(fname, 'w');
    if fid == -1
        error('You do not have permission to save mesh files.');
    end
    fprintf(fid, 'OFF\n');
    fprintf(fid, '%d %d %d\n', length(v), length(f), 0);
    fprintf(fid, '%f %f %f\n', v');
    bigdude = [repmat(size(f, 2), size(f, 1), 1) f-1];
    fprintf(fid, '%d %d %d %d\n', bigdude');
    fclose(fid);
end

