function [oog]=isItOog(currentLMS,stimIntensity)
global M
global bkgndLMS

%currentLMS=currentLMS*stimIntensity;
%currentLMS=bkgndLMS.*(currentLMS+1);
%rgb=currentLMS*inv(M);
rgb=inv(M)*(bkgndLMS.*(currentLMS'*stimIntensity+1));

if min(rgb)>0 && max(rgb)<1
    oog=0;
else
    oog=1;
end
%if (oog)
%    keyboard
%end