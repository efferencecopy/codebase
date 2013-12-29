function out = rgb2RGB(in, invGammaTable)
% function out = rgb2RGB(in, invGammaTable)
%
%      Takes an image in normalized rgb intensity units (0,1) and returns
% the normalized RGB voltages (0,1) needed to achieve those intensities.
% Doesn't include the ambient right now.
%
%      "in" is an image which can be Nx3 or NxMx3 or what have you (one of
% the dimensions has to have size , other than that, anything goes).
%      "invGammaTable" is the inverse gamma table (converts intensities to 
% voltages) which can be obtained via InvertGamma.m.
%
%      This function treats the fist dimension of "in" with a size '3' as   
% the rgb dimension.  This could get tricky if your image has three pixels
% in one of the two spatial dimensions.
%
% GDLH 6/21/07

sizein = size(in);
rgbdim = min(find(sizein == 3));
if (rgbdim == 0)
    error('No dimension of input has size 3');
end
nonrgbdims = 1:ndims(in);
nonrgbdims(rgbdim) = [];
reshapedin = permute(in, [nonrgbdims rgbdim]);  % After this, the rgb dimension comes last
dims = size(reshapedin);
reshapedin = reshape(reshapedin, length(reshapedin(:))/3, 3);
[rgb, idx1, idx2] = unique(reshapedin, 'rows');
maxgamidx = size(invGammaTable,1)-1;
RGBim = zeros(size(reshapedin,1),3);
for i = 1:size(rgb,1)
    RGB = [invGammaTable(round(rgb(i,1)*maxgamidx)+1,1),...
           invGammaTable(round(rgb(i,2)*maxgamidx)+1,2),...
           invGammaTable(round(rgb(i,3)*maxgamidx)+1,3)];
    L = logical(idx2 == i);
    RGBim(L,:) = repmat(RGB,sum(L),1);
end
reshapedout = reshape(RGBim, [sizein(nonrgbdims), 3]);
out = reshape(RGBim,sizein);
