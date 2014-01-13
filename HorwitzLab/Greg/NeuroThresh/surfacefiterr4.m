function out = surfacefiterr4(data, params, Loog)

    % An attempt at making a unifying framework for fitting pairs of
    % parallel planes, ellipsoids, and bent planes.  Trying to take
    % advantage of the fact that a pair of parallel planes can be written
    % as a quadratic form.  We can use this as an initial guess and tweak
    % the coefficients around until the sum of the (log, radial) squared
    % errors are minimized.  This worked well.
    %
    % GDLH 3/14/11
    
    nparams = length(params);
    if (nparams ~=6)
        error('Wrong number of parameters passed to surfacefiterr4');
    end

    if (nargin < 3)
        Loog = false(size(data,1),1);
    else
        Loog = logical(Loog);
    end
    
    a = params(1);
    b = params(2);
    c = params(3);
    d = params(4);
    e = params(5);
    f = params(6);
    
    x = data(:,1);
    y = data(:,2);
    z = data(:,3);

    % First converting data into polar coordinates
    [th,ph,r] = cart2sph(x,y,z);
    
    predr2 = 1./(a.*(cos(ph).*cos(th)).^2 +...
                b.*(cos(ph).*sin(th)).^2 +...
                c.*sin(ph).^2 + ...
                2*d.*cos(ph).*cos(th).*cos(ph).*sin(th) +...
                2*e.*cos(ph).*cos(th).*sin(ph) +...
                2*f.*cos(ph).*sin(th).*sin(ph));
    if (any(predr2<0))
       % disp('negative predicted r2');
    	predr2(predr2<0) = inf;
    end
    predr = sqrt(predr2);
    resid = log(predr)-log(r);
    resid(Loog&(predr > r)) = 0;
    
    out = sum(resid.^2);
end