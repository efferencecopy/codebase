% Script to test out ideas for making sure that linear combinations of
% truncated standard normals stay in the gamut.  Why isn't it sufficient to
% make sure that the variances of all of the marginals are one?  Somehow,
% even though the marginal is a linear combination of truncated standard
% normals, and it's constrained to have a variance of '1', it is not
% constrained to be truncated. e.g. it might fall outside of the gamut.

ndims = 300;
nbasisvects = 20;
x = normrnd(0,1,ndims, nbasisvects);
M = x*x';
sds = sqrt(diag(M));
M = M.*((1./sds)*(1./sds)');
sqM = real(sqrtm(M));
niter = 1000;

for i = 1:niter
    i
    p = zeros(ndims,1);
    while any(p < 0.005 | p > 0.995)
        disp('rejecting');
        x = normrnd(0,1,300,1);
        p = normcdf(x,0,1);
    end
    y = sqM*x;
    p = normcdf(y,0,1);
    if any(p < 0.005 | p > 0.995)
        keyboard
    end
end

%%% Hmmm... not sure what to do.  I have to find a way of manipulating M 
% (and the distribution of x) so that the marginal distributions of the
% elements of M*x are truncated Gaussians.  

%%
% OK, I can do it by starting with normals, going through sqM, going
% through the normCDF and then the truncated normal inverseCDF.  That
% sounds like a lot of operations.

% Creating the truncated normal inverse cdf
% This may be really slow.  There's probably a better way.
p = linspace(0,1,10000);
ninv = norminv(p,0,1);
smallidx = find(p<0.005,1,'last')
largeidx = find(p>0.995,1,'first')
truncnorm.x = linspace(0,1,largeidx-smallidx-1);
truncnorm.inv = ninv(smallidx+1:largeidx-1);

ndims = 300;
nbasisvects = 20;
x = normrnd(0,1,ndims, nbasisvects);
M = x*x';
sds = sqrt(diag(M));
M = M.*((1./sds)*(1./sds)');
sqM = real(sqrtm(M));
niter = 3000;

data = [];
for i = 1:niter
    i
    x = normrnd(0,1,300,1);
    y = sqM*x;
    p = normcdf(y,0,1);
    y = interp1(truncnorm.x, truncnorm.inv, p);
    p = normcdf(y,0,1);
    if any(p < 0.005 | p > 0.995)
        keyboard
    end
    data = [data; var(x), var(y)];
end

%%
% Doing it again, but this time using an algorithm that I might want to use
% in an actual experiment: 1) pick uniform random numbers using EJ's routine.
% 2) Convert them to normal using the normal inverse cdf. 3) Matrix
% multiply. 4) Convert to gun voltages by going through a function that is
% the concatenation of (a) normcdf (to get uniform rvs), (b) inverse
% truncated normal cdf (to get truncated normals), (c) inverse gamma
% function (to get gun volatges).  Whew.
% Just doing it one step at a time here.

% Loading a datafile just to get the gamma functions
stro = nex2stro(findfile('K051309001.nex'));
gammaTable = stro.sum.exptParams.gamma_table;
gammaTable = reshape(gammaTable, length(gammaTable)/3, 3);
gammaTable1 = interp1(linspace(0,255,256),gammaTable,linspace(0,255,65536), 'spline');
invgamma = InvertGamma(gammaTable, 1);

gausslocut = 0.005;
gausshicut = 0.995;
minval = norminv(gausslocut);
maxval = norminv(gausshicut);

% Gun intensity units relative to background
sigma = [.15 .15 .15];
bkgndRGB = stro.trial(1,[end-3, end-2, end-1]);
mu = [gammaTable(bkgndRGB(1),1) gammaTable(bkgndRGB(2),2) gammaTable(bkgndRGB(3),3)];

% Making a truncated normal inverse CDF by using the standard normal
% inverse CDF and then stretching it slightly in x so that
% f(0) = minval and f(1) = maxval.
p = linspace(gausslocut,gausshicut,2^16);
truncnorm.inv = norminv(p);
truncnorm.x = linspace(0,1,2^16);

% Finding a function that converts standard normals into truncated normals.
% f(y) = truncnorminv(normcdf(y)).  This is mostly for debugging.
%normconvert.x = linspace(-6,6,2^16);  % Standard normals essentially never leave the range [-6:6]
%normconvert.y = interp1(truncnorm.x, truncnorm.inv, normcdf(normconvert.x));
%plot(normconvert.x,normconvert.y);

% Finding a function that converts standard normals into truncated normals.
% with the appropriate mean and variance.
normconvert.x = linspace(-6,6,2^16)';  % Standard normals essentially never leave the range [-6:6]
y = interp1(truncnorm.x, truncnorm.inv, normcdf(normconvert.x))';
for gun = 1:3
    normconvert.y(:,gun) = y*sigma(gun)+mu(gun);
end
plot(normconvert.x,normconvert.y);

% Now integrating the inverse gamma functions with the normconvert
% function.
for gun = 1:3
    gaussgamma(:,gun) = round(invgamma(round(normconvert.y(:,gun)*2^16),gun)*2^16);
end
    
% Setting up a matrix for introducing correlation among the stimulus elements
nstix = 200;
nbasisvects = 20;
niter = 100;
x = normrnd(0,1,3*nstix, nbasisvects);
M = x*x';
sds = sqrt(diag(M));
M = M.*((1./sds)*(1./sds)');
M = eye(3*nstix);   % for debugging
sqM = real(sqrtm(M));
for i = 1:niter
    % Creating the stimuli
    x = getEJrandnums(3*nstix,unidrnd(2^16));  % uniform [0:2^16-1]
    xx = norminv((x+1)./(2^16+1)); % quantized, truncated normals [-4.1696:4.1696]
    y = sqM*xx; % y's are approximately standard normal
    
    % Converting y's into indicies into the "gaussgamma" lookup table
    slope = length(normconvert.x)/(normconvert.x(end)-normconvert.x(1));
    intercept = 2^15;   % 2^16/2
    yidxs = round(slope*y+intercept);
    z = [];
    for gun = 1:3
        idxs = [1:nstix]+nstix*(gun-1);
        z(idxs) = gaussgamma(yidxs(idxs),gun);
    end
end
figure;
subplot(3,1,1);
hist(xx,20); % should be quantized, truncated normal
subplot(3,1,2);
hist(normcdf(y),20)   % Should be approximately uniform
% Undoing the monitor gamma
intensities = [];
for gun = 1:3
    idxs = [1:nstix]+nstix*(gun-1);
    intensities(:,gun) = interp1(linspace(0,1,size(gammaTable,1)), gammaTable(:,gun), z(idxs)/2^16);
end
subplot(3,1,3);
hist(intensities,20)   % Should be truncated normal
mean(intensities)
std(intensities)
