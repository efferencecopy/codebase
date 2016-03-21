function info = loadTIFF_info(fpath)

warnID = 'MATLAB:imagesci:tiffmexutils:libtiffWarning';
warning('off', warnID);

tifflink = Tiff(fpath, 'r');

info.height_pix = tifflink.getTag('ImageLength');
info.width_pix = tifflink.getTag('ImageWidth');
info.bitdepth = tifflink.getTag('BitsPerSample');

tifflink.close();