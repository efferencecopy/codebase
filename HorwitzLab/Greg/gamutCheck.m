function [h,scalar] = gamutCheck(cc, bkgndrgb, M, tail)
% function [h,scalar] = gamutCheck(cc, bkgndrgb, M, tail)
%
% Will determine whether a requested cone contrast is physically obtainable
% on the monitor.
%
%   OUTPUT
%       h: 1 = in gamut, 0 = not in gamut
%       scalar: the largest scalar by which cc can be multiplied by and still
%   fit inside the gamut.
%
%   INPUT
%       bkgndrgb: RGBs of the background in normalized intensity units
%       cc: Cone contrast desired (direction and amplitude)
%       M: The 'M' matrix obtained by taking the dotproducts of the
%       phosphor emission spectra and the cone fundamentals
%       tail: 'both' means we're checking for a biphasic stimulus (like a
%       grating).  'single' means we're checking for a monophasic stimulus.

if (~strcmp(tail,'both') & ~strcmp(tail,'single'))
    error('Third argument must be either ''both'' or ''single''');
end

if (size(bkgndrgb,1) < size(bkgndrgb,2))
   bkgndrgb = bkgndrgb';  % forcing a column vector
end
if (size(cc,1) < size(cc,2))
   cc = cc';  % forcing a column vector
end
% cc
% bkgndrgb
% M
% tail

bkgndlms = M*bkgndrgb;
rgb = inv(M)*(bkgndlms.*(cc));  % rgb is is delta units from bkgnd

scalefactors = [(1-bkgndrgb(1))/rgb(1);...
    (1-bkgndrgb(2))/rgb(2);...
    (1-bkgndrgb(3))/rgb(3);...
    (-bkgndrgb(1))/rgb(1);...
    (-bkgndrgb(2))/rgb(2);...
    (-bkgndrgb(3))/rgb(3)];

% No flipping about the origin if we're just considering monopolar stimuli
if (strcmp(tail,'single'))
    scalefactors(scalefactors < 0) = [];
end

signs = sign(scalefactors);
scalefactors = abs(scalefactors*(1-2*eps));  % A hack to avoid round off errors
% abs to avoid flipping the polarity
collisionpt = [];
for i = 1:size(scalefactors,1)
    collisionpt(i,:) = rgb*signs(i)*scalefactors(i)+bkgndrgb;
end
Limpossible = any(collisionpt<0,2) | any(collisionpt>1,2);
scalar = min(abs(scalefactors(~Limpossible)));
h = scalar > 1;
end

% For testing
% M =[0.0715    0.1314    0.0187;...
%     0.0259    0.1366    0.0275;...
%     0.0022    0.0104    0.1042];
% 
% bkgndrgb = [.57 .48 .41]'
% cc = [.25 .25 .25]';
% 

% Sanity checks.  
% rgb intensity at peak
%inv(M)*(bkgndlms+stimlms)
% rgb intensity at trough
%inv(M)*(bkgndlms-stimlms) 

% % Should be out
% inv(M)*(bkgndlms.*(1+stimcc))
% inv(M)*(bkgndlms.*(1-stimcc))
%
%
% stimcc./cc and scalar should be equal
% stimrgb = rgb*scalar;
% stimlms = M*stimrgb;
% stimcc = stimlms./bkgndlms;
% stimcc./cc
% 
% scalar


