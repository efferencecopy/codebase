function out = surfacefiterr(data, params, Loog)

    % Function for computing the error of a parametric function fit to 
    % 3-D data assuming radial error.  
    % Three classes of fit are supported: 
    % 1) Planar: ax+by+cz+1 = 0
    % 2) Quadratic without cross-terms: 
    %       ax+by+cz+dx^2+ey^2+fz^2+1 = 0;
    % 3) Quadratic with cross-terms:
    %       ax+by+cz+dx^2+ey^2+fz^2+gxy+hyz+ixz+1 = 0;

    % Which fit is used depends on how many parameters are provided in the
    % PARAMS parameter vector.  If PARAMS has three elements, we use a
    % plane, if it has six elements, we use a quadratic surface without
    % cross-terms, if it has nine elements, we use a full quadratic fit.
    
    % The optional "Loog" argument is for a binary vector where 0 means
    % that the corresponding point in the "data" argument was in the gamut
    % of the monitor and a 1 means that it was not.
    
    nparams = length(params);
    if (nparams ~=3 && nparams ~=6 & nparams ~=9)
        error('Wrong number of parameters passed to surfacefiterr');
    end
    if (nargin < 3)
        Loog = false(size(data,1),1);
    else
        Loog = logical(Loog);
    end
    fittype = 0;
    a = params(1);
    b = params(2);
    c = params(3);
    if (nparams > 3)
        fittype = 1;
        d = params(4);
        e = params(5);
        f = params(6);
    end
    if (nparams > 6)
        fittype = 2;
        g = params(7);
        h = params(8);
        i = params(9);
    end
    
    x = data(:,1);
    y = data(:,2);
    z = data(:,3);
    keyboard
    % First converting data into polar coordinates
    % At the moment, using Will Kleiber's convention, but this should be
    % changed to the native matlab convention.
    r = sqrt(x.^2+y.^2+z.^2);
    th = acos(y./r);
    ph = atan2(z,x);
    L = logical(ph<0);
    ph(L) = -ph(L);
    ph(~L) = 2*pi-ph(~L);
    
    if fittype == 0
        predr = -1./(a.*sin(th).*cos(ph)+b.*cos(th)-c.*sin(th).*sin(ph));
       % if (any(predr < 0))
       %     keyboard
       %     out = 10^10;
       %     return
       % else
        predr = abs(predr);
        resid = log(predr)-log(r);
        %end
    else
        % Getting coefficients for quadratic equation in r
        % C is all ones.
        if (fittype == 1)
            A = (d*x.^2+e*y.^2+f*z.^2)./(r.^2);
        else
            A = (d*x.^2+e*y.^2+f*z.^2+g.*x.*y+h.*y.*z+i.*x.*z)./(r.^2);
        end
        B = (a*x+b*y+c*z)./r;
        
        predr1 = 2./(-B+sqrt(B.^2-4.*A));
        predr2 = 2./(-B-sqrt(B.^2-4.*A));
        
        % These are obviously bad solutions - return a big error
        %Lbad1 = ~Loog&(predr1<0 | ~(real(predr1)==predr1) | isinf(predr1));
        %Lbad2 = ~Loog&(predr2<0 | ~(real(predr2)==predr2) | isinf(predr2));
        Lbad1 = predr1<0 | ~Loog&(~(real(predr1)==predr1) | isinf(predr1));
        Lbad2 = predr2<0 | ~Loog&(~(real(predr2)==predr2) | isinf(predr2));

        if any(Lbad1 & Lbad2)
            out = 10^10;
            return;
        end

        predr = nan*ones(size(predr1));
        predr(~Lbad1) = predr1(~Lbad1);
        predr(~Lbad2) = predr2(~Lbad2);
        L = ~Lbad1 & ~Lbad2;
        predr(L) = min([predr1(L), predr2(L)],[],2);
        resid = log(predr)-log(r);
      %  if (any(any(real(resid)~=resid & Loog)))
      %      keyboard
      %  end
        resid(real(resid)~=resid & Loog) = 0;  % imaginary solutions are OK for OOGs
    end
    resid(Loog&(predr > r)) = 0;
  %  resid(Loog&(predr < r)) = 10*resid(Loog&(predr < r));
    
    if ~isreal(resid)
        keyboard
    end
    out = sum(resid.^2);
end