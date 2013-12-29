function threshSurfPlot(colors, alphas, viewangle, plottype, fpar)

% unpack the color directions used. turn them into unit vectors.
if iscell(colors)
    colors = cat(1, colors{:});
end
colordirs = bsxfun(@rdivide, colors, sqrt(sum(colors.^2, 2))); % make sure they're unit vecs
colordirs = [colordirs; -colordirs];

% unpack the threshold estimates.
if iscell(alphas)
    alphas = cat(1, alphas{:}); %in CC b/w 0 and 100
end
alphas = [alphas; alphas];



%  Plotting
switch lower(plottype)
    case '3d'
        % raw data first
        cordinates = bsxfun(@times, colordirs, alphas);
        
        hold on,
        set(gca, 'linewidth', 1, 'fontsize', 22)
        plot3(cordinates(:,1), cordinates(:,2), cordinates(:,3),'ko', 'markersize', 7, 'markerfacecolor', 'k')
        
        % Plotting the fit
        plotlims = max(abs(cordinates)).*1.1;
        [x,y,z] = meshgrid(linspace(-plotlims(1),plotlims(1),150),linspace(-plotlims(2),plotlims(2),150),linspace(-plotlims(3),plotlims(3),50));
        tmp = [x(:) y(:) z(:)];
        v = sum(abs(tmp * reshape(fpar(2:end),3,3)).^fpar(1),2);
        fv = isosurface(x,y,z,reshape(v,size(x,1),size(x,2),size(x,3)),1);
        patch(fv,'EdgeAlpha',0,'FaceAlpha',.2);
        
        % annotations and such
        axis equal
        xlabel('\DeltaL/L')
        ylabel('\DeltaM/M')
        zlabel('\DeltaS/S')
        
        switch viewangle
            case 'lm plane'
                newview = [0 0 1];
            case 'isolum'
                newview = [1 1 0];
            case 's vs l+m'
                newview = [1 -1 0];
        end
        view(newview)
        
    case '2d'
        
        % determine what the basis vectors are
        switch lower(viewangle)
            case 'isolum'
                xvec = [1/sqrt(2), -1/sqrt(2), 0];
                yvec = [0 0 1];
            case 'lm plane'
                xvec = [1, 0, 0];
                yvec = [0, 1, 0];
            case 's vs l+m'
                xvec = [1/sqrt(2), 1/sqrt(2), 0];
                yvec = [0 0 1];
        end
        
        % project the raw data onto the appropriate plane
        cordinates = bsxfun(@times, colordirs, alphas);
        monkThreshCords = cordinates * [xvec(:), yvec(:)];
        
        % determine the modeled threshold in that plane
        angles = linspace(0,2*pi, 1000);
        tmp = (cos(angles(:)) * xvec) + (sin(angles(:)) * yvec);
        v = sum(abs(tmp * reshape(fpar(2:end),3,3)).^fpar(1),2);
        modThresh = (1./v).^(1./fpar(1));
        
        % project the model threshold onto the plane
        modThreshCords = bsxfun(@times, tmp, modThresh); %scale the color dirs by the vector norm
        modThreshCords = modThreshCords * [xvec(:), yvec(:)];
        
        % plot the data
        hold on,
        set(gca, 'linewidth', 1 , 'fontsize', 22)
        plot(monkThreshCords(:,1), monkThreshCords(:,2), 'k.', 'markersize', 8)
        plot(modThreshCords(:,1), modThreshCords(:,2), '-', 'linewidth', 3, 'color', [.7 .7 .7]);
        axis equal
        
end








