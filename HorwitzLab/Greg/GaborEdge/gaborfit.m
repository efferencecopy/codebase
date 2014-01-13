% gaborfit.m
%
% Fits a gabor function to a square 2-D grayscale image (like an STA).
% Requires the error function gaborfiterr.m which it calls from
% fminsearch.
%
% input argument should be a square matrix.
% output is a structure with several fields:
%   out.theta: orientation in radians (0 horizonal, pi/2 vertical) 
%   out.lambda: period in pixels/cycle 
%   out.phi: phase in radians
%   out.sigma: sigma of Gaussian component
%   out.gamma: distention factor orthogonal to theta 
%   out.xoffset: x offset in pixels
%   out.yoffset: y offset in pixels
%   out.exitflag: 1 = converged, 0 = did not converge
%
% See comment at the end of the function to see how to sample the fitted gabor.
%
% GDLH 8/5/08

function out = gaborfit(im)
    if (size(im,1) ~= size(im,2))
        error('input matrix must be square');
    end
    nstixperside = size(im,1);
    im = im./max(abs(im(:)));

    % Getting initguess parameters
    imsq = im.^2;
    interval = [1:nstixperside];
    centerinterval = ceil(median(interval));
    interval = interval-centerinterval;
    initguess.xoffset = (interval*sum(imsq)')./sum(imsq(:));
    initguess.yoffset = -(interval*sum(imsq')')./sum(imsq(:));

    % Looking at the power spectrum
    powerspectra = abs(fftshift(fft2(im))).^2;
    powerspectra = powerspectra .* (powerspectra > .1*max(powerspectra(:)));
    powerspectra = powerspectra./sum(powerspectra(:));
    phases = angle(fftshift(fft2(im)));

    [i,j] = meshgrid(interval,interval);
    angles = atan2(i,j)-pi/2;
    amp = sqrt(i.^2+j.^2);
    [peaksfy,peaksfx] = ind2sub(size(powerspectra),find(powerspectra == max(powerspectra(:)),1));
    initguess.lambda = nstixperside./sum(sum(amp.*powerspectra));
    cosines = powerspectra.*cos(angles);
    sines = powerspectra.*sin(angles);
    if (peaksfy<=ceil(median(interval)))
        meansin = mean(mean(sines(1:centerinterval,:)));
    else
        meansin = mean(mean(sines(centerinterval:nstixperside,:)));
    end
    if (peaksfx<=ceil(median(interval)))
        meancos = mean(mean(cosines(:,1:centerinterval)));
    else
        meancos = mean(mean(cosines(:,centerinterval:nstixperside)));
    end
    initguess.theta = atan2(meansin,meancos);
    initguess.theta = mod(initguess.theta + pi/2,2*pi);

    if (initguess.lambda>2*nstixperside)
        initguess.lambda = 5*nstixperside;
        initguess.sigma = sqrt(sum(abs(im(:)) > 0.5));
        initguess.gamma = 1;
        initguess.phi = 0;
    else
        initguess.sigma = initguess.lambda/3.5;
        initguess.gamma = 1.4;
        proj = [initguess.xoffset initguess.yoffset]*[cos(pi/2-initguess.theta); -sin(pi/2-initguess.theta)];
        [i,j] = ind2sub(size(powerspectra),find(max(powerspectra(:))==powerspectra,1));
        initguess.phi = phases(i,j);
        err1 = gaborfiterr([initguess.theta,initguess.lambda,initguess.phi,initguess.sigma,initguess.gamma,initguess.xoffset,initguess.yoffset],im);
        err2 = gaborfiterr([initguess.theta,initguess.lambda,initguess.phi+pi,initguess.sigma,initguess.gamma,initguess.xoffset,initguess.yoffset],im);
        if (err2<err1)
            initguess.phi = initguess.phi+pi;
        end
    end

    % Preparing for fitting
    theta = initguess.theta;
    lambda = initguess.lambda;  % pixels per cycle
    phi = initguess.phi;
    sigma = initguess.sigma;
    gamma = initguess.gamma;
    xoffset = initguess.xoffset;
    yoffset = initguess.yoffset;

    % Doing the fitting
    options = optimset('MaxIter',2000,'MaxFunEvals',3000,'Display','off');
    [x,fval,exitflag] = fminsearch(@(x)gaborfiterr(x,im),[theta,lambda,phi,sigma,gamma,xoffset,yoffset],options);
    out.theta = x(1);
    out.lambda = x(2);
    out.phi = x(3);
    out.sigma = x(4);
    out.gamma = x(5);
    out.xoffset = x(6);
    out.yoffset = x(7);
    out.exitflag = exitflag;
    
    % Calculating the fitted gabor
    % [X, Y] = meshgrid(interval(1:2:end),interval);
    % X = X-xoffset; Y = Y+yoffset; % So negative numbers mean down
    % xprime = X.*cos(-theta)+Y.*sin(-theta);
    % yprime = -X.*sin(-theta)+Y.*cos(-theta);
    % fittedgabor = exp(-(xprime.^2+gamma.^2.*yprime.^2)./(2.*sigma.^2)).*cos(2.*pi.*yprime./lambda-phi);

end
