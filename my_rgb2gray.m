function I = my_rgb2gray(img, info, coef)
%RGB2GRAY Convert RGB image or colormap to grayscale.


if (ndims(img)~=3);
    error('Image must be 3D')
end

% convert to double if need be
wasfloat = true;
if ~(isa(img, 'double') || isa(img, 'single'))
    wasfloat = false;
    img = double(img)./double(info.MaxSampleValue(1));
end

% Shape input matrix so that it is a n x 3 array and initialize output matrix
origSize = size(img);
img = reshape(img(:),origSize(1)*origSize(2),3);
sizeOutput = [origSize(1), origSize(2)];
I = img * coef';

% normalize so that the values are back to 0 to 1 range
I = I./max(I);

%Make sure that the output matrix has the right size
I = reshape(I,sizeOutput);

% bring things back
if ~wasfloat
    I = round(I.*double(info.MaxSampleValue(1))); %bring it back to int world
    I = uint8(I);
end
