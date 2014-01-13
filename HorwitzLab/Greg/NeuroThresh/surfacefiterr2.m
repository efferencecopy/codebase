function out = surfacefiterr2(data, params, Loog)

    % Function for computing the error of a parametric function fit to 
    % 3-D data assuming radial error.  We fit symmetric pairs of surfaces 
    % to the points (no need to preassign the points to two clusters).
    % Two classes of fit are supported: 
    % 1) Planar: ax+by+cz+1 = 0
    % 2) Quadratic: 
    %       ax+by+cz+dx^2+ey^2+fxy+1 = 0;
    % This quadratic form assumes that the data have been rotated such that
    % they are spread out mostly in x and y (z is basically a function on x
    % and y).

    % Which fit is used depends on how many parameters are provided in the
    % PARAMS parameter vector.  If PARAMS has three elements, we use a
    % plane, if it has six elements, we use a quadratic surface.
    
    nparams = length(params);
    if (nparams ~=3 && nparams ~=6)
        error('Wrong number of parameters passed to surfacefiterr');
    end

    if (nargin < 3)
        Loog = false(size(data,1),1);
    else
        Loog = logical(Loog);
    end
    
    fittype = 0;  % 0 = plane; ax+by+cz = 1
    a = params(1);
    b = params(2);
    c = params(3);
    if (nparams > 3)
        fittype = 1;  % 1 = quad; ax+by+cz+dx^2+ey^2+fz^2 = 1
        d = params(4);
        e = params(5);
        f = params(6);
    end
    
    x = data(:,1);
    y = data(:,2);
    z = data(:,3);

    % First converting data into polar coordinates
    [th,ph,r] = cart2sph(x,y,z);
    
    if fittype == 0
        predr = -1./(params(1).*cos(ph).*cos(th)+params(2).*cos(ph).*sin(th)+params(3).*sin(ph));  
        resid = log(abs(predr))-log(r);
    else
        % Getting coefficients for quadratic equation in r.
        % C is all ones.
        A = (d*x.^2+e*y.^2+f*x.*y)./(r.^2);
        B = (a*x+b*y+c*z)./r;
        
        predr1 = 2./(-B+sqrt(B.^2-4.*A));
        predr2 = 2./(-B-sqrt(B.^2-4.*A));
        
        % These are obviously bad solutions - return a big error
        Lbad1 = ~(real(predr1)==predr1) | isinf(predr1);
        Lbad2 = ~(real(predr2)==predr2) | isinf(predr2);
        
        if any(~Loog&Lbad1 & Lbad2)
            out = 10^10;
            return;
        end
        predr1(Lbad1) = Inf;
        predr2(Lbad2) = Inf;
%keyboard        
        predr = nan*ones(size(predr1));
        predr(Lbad1) = predr2(Lbad1);
        predr(Lbad2) = predr1(Lbad2);
        L = ~Lbad1 & ~Lbad2;
        predr(L) = min([abs(predr1(L)), abs(predr2(L))],[],2);
        resid = log(abs(predr))-log(r);
    end
    resid(Loog&(abs(predr) > r)) = 0;
    if ~isreal(resid)
        keyboard
    end
    out = sum(resid.^2);
end