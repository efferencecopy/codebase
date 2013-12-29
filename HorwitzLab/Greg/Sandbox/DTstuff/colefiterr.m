function err = colefiterr(params, lms, Loog, forceellipsoid)
% fminsearch only searches over a single vector of parameters
% so even though we want to minimize the error over 'beta' and
% 'mechanisms', which are logically distinct, we have to bundle
% them together.
if (forceellipsoid)
    beta = 2;
else
    beta = params(1);    
end

[th,phi,r] = cart2sph(lms(:,1),lms(:,2),lms(:,3));
predr = zeros(length(r),1);
options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-10);

if (length(params) == 3) % we're fitting a 1-mechanism (linear) model
    wholesum = abs([cos(phi(:)).*cos(th(:)) cos(phi(:)).*sin(th(:)) sin(phi(:))]*params(:));
    predr = (1./wholesum);
end
if (length(params) == 7) % we're fitting a 2-mechanism model
    wholesum = sum(abs([cos(phi(:)).*cos(th(:)) cos(phi(:)).*sin(th(:)) sin(phi(:))]*[reshape(params(2:end),3,2), [0; 0; 0]]).^params(1),2);
    predr = (1./wholesum).^(1/params(1));
end
if (length(params) == 10) % we're fiting the full-blown three parameter model
    wholesum = sum(abs([cos(phi(:)).*cos(th(:)) cos(phi(:)).*sin(th(:)) sin(phi(:))]*reshape(params(2:end),3,3)).^params(1),2);
    predr = (1./wholesum).^(1/params(1));
end

resid = log(predr)-log(r);
resid(Loog&(predr > r)) = 0;
err = sum(abs(resid));
%err = sum(resid.^2);