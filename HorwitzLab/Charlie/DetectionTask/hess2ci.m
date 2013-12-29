function ci = hess2ci(monk, cell, confLevel);

% This function accepts DTmocs data and creates confidence intervals for
% threshold ratios. I'm using code that takes the hessian of weibull
% fitting as an estimate of SE for thresholds. Then I convert these to a
% couchy. The CI's are based off the couchy... I'm also robbing code from
% greg
%
% CAH and GDLH 11/09

TRs = cell.alpha ./ monk.alpha;  %one TR for each clr/sfs condition
alpha = 1-confLevel;
for clr = 1:size(TRs,1)
    for sfs = 1:size(TRs,2)
        ci.up(clr, sfs) = nan; %these get replaced with legit values
        ci.lo(clr, sfs) = nan; %if TRs and Hessian estimates are o.k.
        goodTRs = ~any(isnan(TRs(clr,sfs)));
        goodHessErrs = ~any(isnan([cell.err.mle.alpha(clr, sfs), monk.err.mle.alpha(clr, sfs)]));
        if (goodTRs && goodHessErrs)
            muCell = cell.alpha(clr, sfs);
            muMonk = monk.alpha(clr, sfs);
            sigmaCell = cell.err.mle.alpha(clr, sfs);
            sigmaMonk = monk.err.mle.alpha(clr, sfs);
            
            
            % First, a Monte Carlo simulation to determine the domain to
            % evaluate
            niter = 10000;
            X = normrnd(muCell, sigmaCell, niter, 1);
            Y = normrnd(muMonk, sigmaMonk, niter, 1);
            Z = X./Y;
            
            % Now using the formulas from Kamerund 1978
            z = linspace(min(Z).*0.8, max(Z).*1.2, 2000);
            w = (sigmaMonk/sigmaCell)*z;  % w is a scaled version of z
            s = 1./sqrt(w.^2+1);
            k = (muCell/sigmaCell.*w+muMonk/sigmaMonk).*s.^2;
            M = -.5*(muMonk/sigmaMonk.*w-muCell/sigmaCell).^2.*s.^2;
            Q = k.*s.*sqrt(2*pi).*(1-2.*normcdf(-k./s))+(2*s.^2.*exp(-k.^2./(2*s.^2)));
            
            fg = 1/(2*pi).*Q.*exp(M);
            fz = (sigmaMonk/sigmaCell)*fg;
            fz = fz./sum(fz);
            
            % CDF
            Fz = cumsum(fz);
            ci.up(clr, sfs) = interp1(Fz, z, [1-alpha]);
            ci.lo(clr, sfs) = interp1(Fz, z, [alpha]);
        end
    end
end
