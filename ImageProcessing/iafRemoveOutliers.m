function img_out = iafRemoveOutliers(img_in, stdMultiplier)


stdAllData = std(img_in(:));
stdByPix = std(img_in, [] ,3);
stdTooBig = stdByPix > (stdMultiplier .* stdAllData);
normFact = (stdByPix ./ stdAllData);

normFact = normFact .* stdTooBig;
normFact = max(cat(3,normFact, ones(size(normFact))), [], 3);
img_out = bsxfun(@rdivide, img_in, normFact);