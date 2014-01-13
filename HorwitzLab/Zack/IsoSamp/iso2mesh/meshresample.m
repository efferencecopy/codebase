function [node,elem] = meshresample(v, f, keepratio)
%
% [node,elem]=meshresample(v,f,keepratio)
%
% resample mesh using CGAL mesh simplification utility
%
% author: Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
% date: 2007/11/12
%
% input:
%    v: list of nodes
%    f: list of surface elements (each row for each triangle)
%    keepratio: decimation rate, a number less than 1, as the percentage
%               of the elements after the sampling
%
% output:
%    node: the node coordinates of the sampled surface mesh
%    elem: the element list of the sampled surface mesh
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%
    [node,elem] = domeshsimplify(v, f, keepratio);

    if isempty(node)
        warning('meshre:defects', 'Mesh defects detected... attempting to repair');
        [vnew,fnew] = meshcheckrepair(v, f);
        [node,elem] = domeshsimplify(vnew, fnew, keepratio);
    end

    [node,~,J] = unique(node, 'rows');
    elem = J(elem);
    saveoff(node, elem, mwpath('post_remesh.off'));
end


% function to perform the actual resampling
function [node,elem] = domeshsimplify(v, f, keepratio)
    exesuff = '.exe';
    saveoff(v, f, mwpath('pre_remesh.off'));
    deletemeshfile(mwpath('post_remesh.off'));
    system([' "' mcpath('cgalsimp2') exesuff '" "' mwpath('pre_remesh.off') '" ' num2str(keepratio) ' "' mwpath('post_remesh.off') '" > nul']);
    [node,elem] = readoff(mwpath('post_remesh.off'));
end
