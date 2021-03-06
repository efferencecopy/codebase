function detectionContours(colors, sfs, data, plane)
% PLOT DETECTION CONTOURS
%
% EXAMPLE: detectionContours(colors, sfs, data, [plane])
% Plots detection ellipses for data generated by quest experiments. The
% optional "plane" argument specifies which plane in color space to plot.
% The default is the isoluminant plane.
    
    if ~exist('plane', 'var')
        plane = 'isolum';
    end
    
    % strip out the achromatic color directions
    ach = [1 1 1]./norm([1 1 1]);
    achIdx = softEq(ach, colors, [], 'rows');
    if any(achIdx)
        colors(achIdx, :) = [];
        data(achIdx, :, :) = [];
    end
    
    %
    % now plot the contours
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    figure;
    SEM = nanstd(data, 0, 3) ./ sqrt(sum(~isnan(data), 3));
    lowErr = nanmean(data,3)-SEM;
    highErr = nanmean(data,3)+SEM;
    whichSfs = sort(sfs);
    whichSfs = (sfs == whichSfs(1));
    
    switch lower(plane)
        case 'isolum'
            basisOne = [1 -1 0] ./ norm([1 -1 0]); %angles wrt l-m
            basisTwo = [0 0 1];

        case 'lm'
            basisOne = [1 0 0];
            basisTwo = [0 1 0];
    end
    
    x = diag(colors * repmat(basisOne', 1, size(colors, 1)));
    y = diag(colors * repmat(basisTwo', 1, size(colors, 1)));
    thetas = atan(y./x);
    thetas = [thetas; thetas+pi];
    r = nanmean(data(:, whichSfs, :), 3);
    r = [r;r];
    h.fig = polar(thetas, r, 'b.');
    
    %plot 1sd error bars
    hold on,
    lowTmp = [lowErr(:,whichSfs)', lowErr(:,whichSfs)'];
    highTmp = [highErr(:,whichSfs)', highErr(:,whichSfs)'];
    errs = [lowTmp ; highTmp];
    thetas = [thetas' ; thetas'];
    polar(thetas, errs, 'b-')
    set(get(gca, 'children'), 'MarkerSize', 15, 'LineWidth', 2);
    hold off,
    
    %update the uicontrols
    sfsAsText = num2str(sort(sfs)', 3);
    h.uimenu = uicontrol('style', 'popupmenu', 'string', sfsAsText, 'callback', @updatePlot, 'position', [20, 30, 70, 35]);
    uicontrol('style', 'text', 'string', 'Spt. Freq', 'position', [22, 70, 68, 15])
    
    %package relevant data into the user field
    d.colors = colors;
    d.data = data;
    d.sfs = sfs;
    d.basisOne = basisOne;
    d.basisTwo = basisTwo;
    d.hands = h;
    set(gcf, 'userdata', d)
end


%
%   PLOT THE DETECTION CONTOUR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function updatePlot(varargin);

    %determine which sfs to plot
    d = get(gcf, 'UserData');
    uiVal = get(d.hands.uimenu, 'value');
    sortSfs = sort(d.sfs);
    sfsToPlot = sortSfs(uiVal);
    sfsIdx = [d.sfs == sfsToPlot];
      
    %update the plot
    x = diag(d.colors * repmat(d.basisOne', 1, size(d.colors, 1)));
    y = diag(d.colors * repmat(d.basisTwo', 1, size(d.colors, 1)));
    thetas = atan(y./x);
    thetas = [thetas; thetas+pi];
    r = nanmean(d.data(:, sfsIdx, :), 3);
    r = [r;r];
    h.fig = polar(thetas, r, 'b.');
    
    %update the 1 SEM error bars
    SEM = nanstd(d.data, 0, 3) ./ sqrt(sum(~isnan(d.data), 3));
    lowErr = nanmean(d.data,3)-SEM;
    highErr = nanmean(d.data,3)+SEM;
    lowTmp = [lowErr(:,sfsIdx)', lowErr(:,sfsIdx)'];
    highTmp = [highErr(:,sfsIdx)', highErr(:,sfsIdx)'];
    errs = [lowTmp ; highTmp];
    thetas = [thetas' ; thetas'];
    hold on,
    polar(thetas, errs, 'b-')
    set(get(gca, 'children'), 'MarkerSize', 15, 'LineWidth', 2);
    hold off,
end



