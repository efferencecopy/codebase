function out = ConvertConeContrastBasis(Min, Mout, bkgndrgb, cc)
% lmsout = ConvertConeContrastBasis(Min, Mout, bkgndrgb, cc)
%
% Converts contrasts in one basis to contrasts in another basis.
% Useful for looking at NeuroThresh data in cone-contrast spaces 
% defined by different cone fundamentals.
%
% cc should be nx3.

if (size(cc,2) ~= 3)
    cc = cc';
end
if (size(cc,2) ~= 3)
    error('cc should be n x 3');
end
if (size(bkgndrgb,1) ~=3 & size(bkgndrgb,2) ~=3)
    error('bkgndrgb should have 3 elements');
end
if (size(bkgndrgb,1) ~=3)
    bkgndrgb = bkgndrgb';
end

Mconv = Mout/Min;  % matrix that converts cone excitations in stro file
n = size(cc,1);
% irrespective of which fundamentals were used) to 10 deg cone excitations
bkgndlms = Min*bkgndrgb;
lms = cc.*repmat(bkgndlms',n,1)+repmat(bkgndlms',n,1);
lms_out = lms*Mconv';
bkgndlms_out = Mconv*bkgndlms;

out = (lms_out-repmat(bkgndlms_out',n,1))./repmat(bkgndlms_out',n,1);
end