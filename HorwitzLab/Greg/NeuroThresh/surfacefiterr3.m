function out = surfacefiterr3(data, params, Loog)

    % Function for to fit paired (or multiple?) NT data sets from the same
    % neuron under the constraint that the isoresponse surfaces have to be
    % scaled versions of each other. 
    % IN PROGRESS
    % GDLH 6/2/11
    
    nparams = length(params);
    if (nparams ~= 7) % This will need to change when were' doing >2 col dirs
        error('Wrong number of parameters passed to surfacefiterr3');
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
    scalefactors = [1 params(7:end)];  % scalefactor for first group is always '1'
    x = data(:,1);
    y = data(:,2);
    z = data(:,3);
    whichgroup = data(:,4);
    
    if (length(unique(whichgroup)) ~= length(scalefactors))
        keyboard
        error('Wrong number of scalefactors passed in.');
    end
        
    % First converting data into polar coordinates
    resid = [];
    for i = 1:2  % hardcoding two groups
        L = whichgroup == i;
        [th,ph,r] = cart2sph(x(L),y(L),z(L));
        
        predr2 = 1./(a.*(cos(ph).*cos(th)).^2 +...
            b.*(cos(ph).*sin(th)).^2 +...
            c.*sin(ph).^2 + ...
            2*d.*cos(ph).*cos(th).*cos(ph).*sin(th) +...
            2*e.*cos(ph).*cos(th).*sin(ph) +...
            2*f.*cos(ph).*sin(th).*sin(ph))/scalefactors(i);
        predr = sqrt(abs(predr2));
        tmpresid = log(abs(predr))-log(r);
        tmpresid(Loog(L)&(predr > r)) = 0;
        resid = [resid; tmpresid];
    end
    out = sum(resid.^2);
end