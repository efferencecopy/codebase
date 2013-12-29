function im = DrawGaborEdge(bkgndrgb, gaborrgb, edgergb, theta, lambda, sigma, gamma, phi, xoff, yoff, etheta, edisp, gausslim, pixperdeg); 
    %keyboard   
    stimsizeindeg = norminv([1-gausslim gausslim],0,1)*abs(sigma)/min([1 gamma]);
    %stimsizeindeg = [min(stimsizeindeg(1),-.65) max(stimsizeindeg(2),.65)];
    stimsizeinpix = round(2*stimsizeindeg(2)*pixperdeg);
    interval = linspace(stimsizeindeg(1), stimsizeindeg(2), stimsizeinpix);
  %  [X, Y] = meshgrid(interval-xoff/2,interval+yoff/2);
    size(interval)
  [X, Y] = meshgrid(interval-xoff,interval+yoff); 
    xprime = X.*cos(-theta)+Y.*sin(-theta);
    yprime = -X.*sin(-theta)+Y.*cos(-theta);
    fittedgabor = exp(-(xprime.^2+gamma.^2.*yprime.^2)./(2.*sigma.^2)).*cos(2.*pi.*yprime./lambda-phi);
    orientation = theta+etheta;
    if (mod(orientation, pi/2) < 1e-5 || mod (orientation, pi/2) > pi/2-1e-5)
        orientation = pi/2*round(orientation/(pi/2));
    end
    edgeim = cosd(orientation*180/pi)*(Y+edisp*sigma*cosd(orientation*180/pi))+...
            sind(orientation*180/pi)*(X+edisp*sigma*sind(orientation*180/pi));
    edgeim(edgeim == 0) = eps;
    edgeim = sign(edgeim).*exp(-(xprime.^2+gamma.^2.*yprime.^2)./(2.*sigma.^2));
    im = zeros(stimsizeinpix,stimsizeinpix,3);
    for plane = 1:3
        im(:,:,plane) = (fittedgabor.*gaborrgb(plane))+(edgeim*edgergb(plane))+bkgndrgb(plane);
    end
end
