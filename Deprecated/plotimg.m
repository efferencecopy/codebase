function hand = plotimg(img, tags, clrmap)


% generate the plot, assuming that the image is truecolor, forcing the
% scaling to be between the min and max dac values
bits = tags.BitsPerSample(1);
maxdac = 2^bits - 1;
hand = imshow(img, [0, maxdac]);

if exist('clrmap', 'var')
    colormap(clrmap)
else 
    colormap('jet')
end