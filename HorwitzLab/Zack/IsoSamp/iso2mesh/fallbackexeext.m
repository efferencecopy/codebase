function exesuff=fallbackexeext(exesuffix, exename)
%
% exesuff=fallbackexeext(exesuffix, exename)
%
% get the fallback external tool extension names for the current platform
%
% author: Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
%
% input:
%     exesuffix: the output executable suffix from getexeext
%     exename: the executable name
%
% output:
%     exesuff: file extension for iso2mesh tool binaries
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

exesuff=exesuffix;
existresult = exist([mcpath(exename) exesuff],'file');
if(strcmp(exesuff,'.mexa64') && ~existresult) % fall back to i386 linux
    exesuff='.mexglx';
    return;
end
if(strcmp(exesuff,'.mexmaci64') && ~existresult) % fall back to i386 mac
    exesuff='.mexmaci';
end
if(strcmp(exesuff,'.mexmaci') && ~existresult) % fall back to ppc mac
    exesuff='.mexmac';
end

