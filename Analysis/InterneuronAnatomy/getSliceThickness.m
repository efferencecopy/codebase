function thickness = getSliceThickness(fpath, info)

% get the z resolution
ztext = info(1).ImageDescription;
idx = regexpi(ztext, 'spacing=', 'end');
zresolution = str2double(ztext(idx+1 : idx+8));



% unpack each frame of the image
img = nan(info(1).Height, info(1).Width, numel(info), 3);
for i_fr = 1:numel(info)
    tmp = imread(fpath, 'tif', 'index', i_fr);
    for i_clr = 1:3
        img(:,:,i_fr,i_clr) = tmp(:,:,i_clr);
    end
end

% x-project, and then y-project
xproj = mean(img, 1);
yproj = mean(xproj, 2);

curves = squeeze(yproj);
curves = bsxfun(@minus, curves, min(curves, [], 1));
curves = bsxfun(@rdivide, curves, max(curves, [], 1));
curves = curves'; % now each row is a color channel

% set a threshold at 0.4 (40% of max height)
thresh = 0.5;
crossing_up = cat(2, false(3,1), diff(curves > thresh, 1, 2) == 1);
crossing_down = cat(2, false(3,1), diff(curves > thresh, 1, 2) == -1);

X = [1:numel(info)] .* zresolution;
thickness = nan(1,3);
for i_clr = 1:3
    low = X(crossing_up(i_clr,:));
    hi = X(crossing_down(i_clr,:));
    
    % check for mirroring
    if numel(low)>1 || numel(hi)>1
        warning('Mirroring found for ch %d file: %s', i_clr, fpath)
        continue
    end
    
    % check for non-threshold crossings
    if isempty(low) || isempty(hi)
        warning('No thresh crossing found for ch %d file: %s', i_clr, fpath)
        continue
    end
    
    thickness(i_clr) = hi-low;
end


