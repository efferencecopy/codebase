function [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog)
    
% [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog)
%
% Function for fitting surfaces to NeuroThresh data.
% 
% INPUTS
% Scaled is an n x 3 matrix of staircase termination points.
% Loog is a binary vector indicating which color directions went out of gamut.
%   (1 = OOG, 0 = not OOG) 
%
% OUTPUTS
% planeparams = [a b c] where |ax+by+cz| = 1;
% quadparams = [a b c d e f] where 
%       ax^2+by^2+cz^2+dxy+exz+fyz = 1
%
% planeSSE and quadSSE are the sum of squared errors associated with the
% fits.
% planeparams and quadparams are passed back in a coordinate frame that
% differs from the original coordinate frame by a transformation captured 
% by the xformmat matrix (a whitening matrix).
%
% Updated to use the new general framework for fitting both parallel pairs
% of planes and quadric surfaces.
%
% GDLH 3/20/11

    [v,d] = eig(cov([scaled(~Loog,:); -scaled(~Loog,:)]));
    d = diag(d);
    if (min(d) < 2*eps)
        disp('Too few data points for whitening');
        whtmat = eye(3);
    else
        whtmat = v*diag(sqrt(1./d));
        scaled = scaled*whtmat;
    end
    [th,ph,r] = cart2sph(scaled(:,1),scaled(:,2),scaled(:,3));
 
    errs = zeros(50,50);
    tmp = linspace(0,pi,size(errs,1));
    for i = 1:length(tmp)
        for j = 1:length(tmp)
            [tmpa, tmpb, tmpc] = sph2cart(tmp(i),tmp(j),1);
            predr = 1./(tmpa.*cos(ph).*cos(th)+tmpb.*cos(ph).*sin(th)+tmpc.*sin(ph));
            predr = predr.*geomean(r./abs(predr));
            resid = log(abs(predr))-log(r);
            resid(Loog&(abs(predr) > r)) = 0;
            errs(i,j) = sum(resid.^2);
        end
    end
    [cand_i, cand_j] = ind2sub(size(errs),find(errs(:) < prctile(errs(:),5)));
    planeparams = [];
    planeSSE = sum(r.^2);
    for i = 1:size(cand_i,1)
        [a, b, c] = sph2cart(tmp(cand_i(i)), tmp(cand_j(i)),1);
        predr = 1./(a.*cos(ph).*cos(th)+b.*cos(ph).*sin(th)+c.*sin(ph));
        initguess = [a;b;c].*geomean(abs(predr)./r); % more accurate than previous shiftfactor way
        
        % Doing gradient descent for planar fit.
        options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-6,'TolX',10^-6,'Display','off');
        warning('off'); % with fminunc, adding "'LargeScale','off'" to the options will squlech the warnings
        [tmpplaneparams, tmpplaneSSE, exitflag, out4] = fminsearch(@(x) surfacefiterr2(scaled, x, Loog),initguess,options);
        warning('on');
        if (~exitflag)
            keyboard;
        end
        
        if (tmpplaneSSE < planeSSE)
            planeparams = tmpplaneparams;
            planeSSE = tmpplaneSSE;
        end
    end
    initguess = [planeparams(1)^2 planeparams(2)^2 planeparams(3)^2 planeparams(1)*planeparams(2) planeparams(1)*planeparams(3) planeparams(2)*planeparams(3)]';
    [quadparams, quadSSE, exitflag, out4] = fminsearch(@(x) surfacefiterr4(scaled, x, Loog),initguess,options);
    if (~exitflag)
        keyboard;
    end
    xformmat = whtmat;
end
