function p = getvarfrom(ws, name)
%
% p=getvarfrom(ws,name)
%
% get variable value by name from specified work-space
%
% author: Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
%
% input:
%    ws: name of the work-space, for example, 'base'
%    name: name string of the variable
%
% output:
%    p: the value of the specified variable, if the variable does not
%       exist, return empty array
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%
    p = [];
    if evalin(ws, ['exist(''' name ''',''var'')']);
        p = evalin(ws, name);
    end
end

