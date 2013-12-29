function error = gaborfiterr(params,im)
    theta = params(1);
    lambda = params(2);
    phi = params(3);
    sigma = params(4);
    gamma = params(5);
    xoffset = params(6);
    yoffset = params(7);
       
    interval = [1:size(im,1)];
    interval = interval-ceil(median(interval));
	[X, Y] = meshgrid(interval,interval);
    X = X-xoffset; Y = Y+yoffset; % So negative numbers mean down
    xprime = X.*cos(-theta)+Y.*sin(-theta);
    yprime = -X.*sin(-theta)+Y.*cos(-theta);
    gabor = exp(-(xprime.^2+gamma.^2.*yprime.^2)./(2.*sigma.^2)).* ...
        cos(2.*pi.*yprime./lambda-phi);
    error = sum((gabor(:)-im(:)).^2);
end