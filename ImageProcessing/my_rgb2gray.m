function I = my_rgb2gray(img, info, coef)
%RGB2GRAY Convert RGB image or colormap to grayscale.


% do some error checking
if (ndims(img)~=3);
    error('Image must be 3D')
end
if ~strcmpi('truecolor', info.ColorType)
    error('Image must be truecolor')
end
if ~isinteger(img)
    error('Image must be a truecolor (UINT8/16) to process')
end


% define the coefficients for conversion to gray scale. The default values
% in matlab's ind2rgb are 
if ~exist('coef', 'var')
    if any(regexp(info.Filename, '_red'))
        coef = [0.5850, 0.3250, 0];
    elseif any(regexp(info.Filename, '_green'))
        coef = [0.2989, 0.5870, 0.1140];
    end
end

% this line was stolen from Matlab's rgb2gray, it's identical to a matrix
% mulitply, but takes into consideration the class of 'img' (in this case
% it's an INTEGER) and does some rounding under the hood to force the
% scaled values back into integers
I = imapplymatrix(coef, img, class(img));
