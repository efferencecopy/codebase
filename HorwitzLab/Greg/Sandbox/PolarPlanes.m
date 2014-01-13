%%
% Contents
% Section 1: An initial stab at fitting planes to data, minimizing error in 'r'.
% Includes code to test the dependence of the fit on linear transformations of the
% data.
%
% Section 2: Testing the dependence of a quadratic fit on linear
% transformations of the data.
%
% Section 2.1: Testing the dependence of linear and quadratic fit (of real
% data) with different data rotations.  Using the new NTsurfacefit.m which 
% treats the pair of planes as a special case of the general quadric
% surface.
%
% Section 3: Making a few plots for Paul Sampson (statistical consulting
% service) showing fits to an example dataset.
%
% Section 4: Fitting a plane and quadratic surfaces to simulated data and
% plotting these surfaces.
%
% Section 5: Trying a bootstrap test on an actual dataset (need to run some
% code in NeuroThresh.m first)
%
% Section 6: Testing the rejection probability of the bootstrap test by
% simulating data consistent with the null hypothesis.  Ideally, the 
% probability of rejecting the null hypothesis would be alpha.
%
% Section 7: Testing the rejection probability of the bootstrap test on a
% simple univariate linear regression.
%
% Section 8: Trying an idea to try to estimate SNR as a function of color
% direction from NeuroThresh data in an attempt to compare neural responses
% to DTspot data.  The basic idea is to fit the isoresponse surface 
% and the baseline rate, assume Poisson spiking and identical
% contrast response functions across color directions (normalized by
% distance to isoresponse surface).  Now we have a full distribution of
% spike counts at each location in color space.  Find the scale factor for
% which ROC area = 0.82.
%
% Section 8.1: How good is the approximation of a single contrast-response
% function in all color directions (using NeuroThresh)?
%
% Section 8.2: Comparing contrast-response functions measured in DTspot and
% NeuroThresh.  IN PROGRESS
%
% Section 9: What does the white noise could look like in LMS space?  How
% do the preferred color directions as assessed by WN, GT, and NT compare?
%
% Section 10: Playing around with a few models.  Planar isoresponse surface
% + nonlinear contrast response functions.  Also experimenting with
% non-planar isoresponse surfaces parameterized in spherical coordinates.
%%
% Section 1
% Playing around with planes in x,y,z and theta,phi,r

a = .1;
b = .1;
c = .1;
n = 10;

[x,y] = meshgrid(linspace(-3,3,100),linspace(-3,3,100));
% Using equation for plane: ax+by+cz+1 =0;
z= -(1+a*x+b*y)/c;

figure; axes; hold on;
surf(x,y,z);
plot3(0,0,0,'k*');

% Using Will Kleiber's spherical coordinates conventions
% can always switch back later.
% theta is elevation, but 0 means straight up
% phi is counterclockwise angle in the x, z plane
% Surface of plane in th, ph, log(r)
r = sqrt(x.^2+y.^2+z.^2);
th = acos(y./r);
ph = atan2(z,x);
L = logical(ph<0);
ph(L) = -ph(L);
ph(~L) = 2*pi-ph(~L);

figure; axes; hold on;
surf(th,ph,log(r));

% Creating random data and plotting them on same axes as surface
x = normrnd(0,1,n,1);
y = normrnd(0,1,n,1);
z = -(1+a*x+b*y)/c;
r = sqrt(x.^2+y.^2+z.^2);
r = r+normrnd(0,r/10);  % adding some noise.

th = acos(y./r);
ph = atan2(z,x);
L = logical(ph<0);
ph(L) = -ph(L);
ph(~L) = 2*pi-ph(~L);
plot3(th,ph,log(r),'m.');

% Bringing the noisy points back into xyz space
x = r.*sin(th).*cos(ph);
y = r.*cos(th);
z = -r.*sin(th).*sin(ph);

% Trying to fit the plane
% Creating an initial guess

pcs = princomp([x y z]);
mn = mean([x y z]);
initxyz = -sign(mn*pcs(:,3))* pcs(:,3)./norm(mn*pcs(:,3));
lights = [x y z];
[pointertoplane, out2, out3, out4] = fminsearch(@(x) surfacefiterr(lights,x),initxyz);

% Checking to see whether a linear transformation of the data leads to an
% estimate of the plane that is a linear transformation away from the
% original estimate of the plane (it is).
% 
% Here's the logic: 
% The equation for a plane is: [a b c]*[x y z]' = -1;
% Where a,b,c are parameters of the plane and x,y,z are spatial dimensions.
% if I transform [x y z] by a matrix A, I can rewrite this as:
% [a b c]*inv(A)*A*[x y z]'.  So if I feed A*[x y z]' into planefiterr
% instead of [x y z]' I get plane coefficients a',b',c' which are equal to
% [a b c]*inv(A).  To recover [a b c] from [a' b' c'] I have to
% postmultiply by A, or alternatively premultiply by A'.

pointerstoplane = [];
for i = 1:10
    A = normrnd(0,1,3,3);
    transformedxyz = [x y z]*A';
    pcs = princomp(transformedxyz);
    mn = mean(transformedxyz);
    initxyz = -sign(mn*pcs(:,3))* pcs(:,3)./norm(mn*pcs(:,3));
    [newpointertoplane, out2, out3, out4] = fminsearch(@(x) surfacefiterr(transformedxyz,x),initxyz);
    pointerstoplane = [pointerstoplane, A'*newpointertoplane];
end

pointerstoplane;


%%
% Section 2
% Testing the dependence of a quadratic fit on linear transformations of
% the data.  The SSE from the quadratic fit *does* depend somewhat on the
% coordinates in which the data are represented.  It's not as
% bad as the thin plate spline though.

a = normrnd(0,1)/10;
b = normrnd(0,1)/10;
c = normrnd(1,1)/10;

% simulating data points
n = 40;
x = normrnd(0,1,n,1);
y = normrnd(0,1,n,1);
z = -(1+a*x+b*y)/c;
[th,ph,r] = cart2sph(x,y,z);
noise = normrnd(1,.02,n,1);
r = exp(log(r)+log(noise));

% Bringing the noisy points back into xyz space
[x, y, z] = sph2cart(th,ph,r);

% Checking to see whether a linear transformation of the data changes the
% ratio of the SSE between the linear and the quadratic fits.

data = [];
for i = 1:40
    A = normrnd(0,40,3,3);
    % A = eye(3);
    transformedxyz = [x y z]*A';
    [v,d] = eig(cov([transformedxyz; -transformedxyz]));
    d = diag(d);
    whtmat = diag(sqrt(1./d))*v';
    transformedxyz = transformedxyz*whtmat';
    [th,ph,r] = cart2sph(transformedxyz(:,1),transformedxyz(:,2),transformedxyz(:,3));

    errs = zeros(500,500);
    tmp = tan(linspace(-pi/2.001,pi/2.001,size(errs,1)));
    for i = 1:length(tmp)
        for j = 1:length(tmp)
            % Parametrically varying a and b and leaving c=1
            predr = -1./(tmp(i).*cos(ph).*cos(th)+tmp(j).*cos(ph).*sin(th)+sin(ph));
            resid = log(abs(predr))-log(r);
            resid = resid-mean(resid);  % That extra degree of freedom
            errs(i,j) = sum(resid.^2);
        end
    end
    [i,j] = ind2sub(size(errs),find(errs(:) == min(errs(:))));
    predr = -1./(tmp(i).*cos(ph).*cos(th)+tmp(j).*cos(ph).*sin(th)+sin(ph));
    shiftfactor = mean(log(abs(predr))-log(r));
    atmp = tmp(i)*exp(shiftfactor);
    btmp = tmp(j)*exp(shiftfactor);
    ctmp = exp(shiftfactor);
    
    % Changing coordinate frames for gradient descent fits of planar and
    % quadratic surfaces.
    rotmat = MakeOrtho([[atmp; btmp; ctmp;] normrnd(0,1,3,2)]);
    rotmat = [rotmat(:,[2 3]) rotmat(:,1)];
    transformedxyz = transformedxyz*rotmat;
    
    rotabc = inv(rotmat)*[atmp;btmp;ctmp];
    options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-6,'TolX',10^-6);
    
    [planeparams, planeSSE, out3, out4] = fminsearch(@(x) surfacefiterr2(transformedxyz,x),rotabc, options);
    [quadparams, quadSSE, out3, out4] = fminsearch(@(x) surfacefiterr2(transformedxyz,x),[planeparams',0,0,0],options);
    data = [data; planeSSE quadSSE]
end


std(data)./mean(data)  % coefficient of variation
% The quality of the quadratic fit appears to depend very little on the
% initial linear transformation of the data.

%%
% Section 2.1
% Looking at the influence of data rotations on the plane and quadratic
% fits of real data (which has OOGs).  


%NTfilename = 'Sedna/S022210002.nex';  % pan color
%NTfilename = 'Sedna/S050410006.nex';  % planar
NTfilename = 'Kali/K091010004.nex'; % Bowl
NTfilename = 'Sedna/S033010003.nex';  % Cylinder?
NT = nex2stro(['/Users/greghorwitz/NexFiles/Greg/',NTfilename]);
out = NTpreprocess(NT,.4,2);
scaled = out(:,[2:4]).*repmat(out(:,5), 1,3);
Loog = logical(out(:,7));
data = [];
for i = 1:10
    A = normrnd(0,1,3,3);
    warning('off');
    [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled*A, Loog, 'newway');
    warning('on');
    data = [data; planeSSE quadSSE planeparams' quadparams]
end

%%
% Section 3
% Making figures for Paul Sampson et al.
% Need to have "data" from NeuroThreshStuff.m and Laccept
% This is a good file: 'K082509007.nex'

% 1) Points and fitted plane in xyz
Loog = logical(data(Laccept,6));
scaled = data(Laccept,[1:3]) .*repmat(data(Laccept,4), 1,3);


% Not sure how to handle OOG points.  For now, ignoring them.
[pcs, score, d] = princomp(scaled(~Loog,:));
mn = mean(scaled(~Loog,:))';
initxyz = -sign(mn'*pcs(:,3))* pcs(:,3)./norm(mn'*pcs(:,3));
options = optimset('MaxFunEvals',10000,'MaxIter',8000);
[pointertoplane, planeSSE, exitflag, out4] = fminsearch(@(x) surfacefiterr(scaled(~Loog,:),x),initxyz, options);
if (exitflag == 0)
    disp('Trying again with new initial guess');
    [pointertoplane, planeSSE, exitflag, out4] = fminsearch(@(x) surfacefiterr(scaled(~Loog,:),x),mn, options);
end

% Plotting
lims = [min(scaled(~Loog,:)); max(scaled(~Loog,:))];
figure; axes; hold on;
[x,y] = meshgrid(linspace(lims(1,1),lims(2,1),2),linspace(lims(1,2),lims(2,2),2));
% Using equation for plane: ax+by+cz+1 =0;
z=-(1+pointertoplane(1)*x+pointertoplane(2)*y)/pointertoplane(3);
h=surf(x,y,z);
set(h,'FaceAlpha',.75,'FaceColor',[.5 .5 .5]);
plot3(scaled(~Loog,1),scaled(~Loog,2), scaled(~Loog,3),'r.','MarkerSize',16);


%%
% Section 3 cont.
% Thin-plate smoothing spline
% Not sure this rotating thing is doing what I want it to do.
% Might need to work on this.
scaled = data(Laccept,[1:3]) .*repmat(data(Laccept,4), 1,3);
rotdat = scaled(~Loog,:);

r = sqrt(rotdat(:,1).^2+rotdat(:,2).^2+rotdat(:,3).^2);
th = acos(rotdat(:,2)./r);
ph = atan2(rotdat(:,3),rotdat(:,1));
L = logical(ph<0);
ph(L) = -ph(L);
ph(~L) = 2*pi-ph(~L);

[F,p] = tpaps([th ph]',log(r'));
[a,b] = meshgrid(linspace(min(th),max(th),100),linspace(min(ph),max(ph),100));
tmp = reshape(permute(cat(3,a,b),[3,1,2]),[2,size(a,1)*size(b,1)]);
splinepoints = fnval(F,tmp);
figure; axes; hold on;
h = surf(a,b,reshape(splinepoints,100,100));
plot3(th,ph,log(r),'r.');
resid = fnval(F,[th ph]')-log(r)';
splineSSE = sum(resid.^2)
title(num2str(splineSSE));

% Getting planeSSE
[pcs, score, d] = princomp(scaled(~Loog,:));
mn = mean(scaled(~Loog,:))';
initxyz = -sign(mn'*pcs(:,3))* pcs(:,3)./norm(mn'*pcs(:,3));
options = optimset('MaxFunEvals',10000,'MaxIter',8000);
[tmp, planeSSE, exitflag] = fminsearch(@(x) surfacefiterr(scaled(~Loog,:),x),initxyz, options);

% Now bootstrapping the residuals and fitting planes
niter = 2;
nulldist = [];
for i = 1:niter
    tmpresid = exp(resid(unidrnd(length(resid),[1 length(resid)])));
    l = r.*tmpresid'.*sin(th).*cos(ph);
    m = r.*tmpresid'.*cos(th);
    s = -r.*tmpresid'.*sin(th).*sin(ph);
    [pcs, score, d] = princomp([l m s]);
    mn = mean([l m s])';
    initxyz = -sign(mn'*pcs(:,3))* pcs(:,3)./norm(mn'*pcs(:,3));

    [tmp, SSE1, exitflag] = fminsearch(@(x) surfacefiterr([l m s],x),initxyz, options);
    Fboot = tpaps([th ph]',log(r'+tmpresid));
    SSE2 = sum((fnval(F,[th ph]')-log(r')-log(tmpresid)).^2);

    nulldist = [nulldist; SSE1./SSE2];
end
p = sum(nulldist<planeSSE/splineSSE)./niter;
%%
% Section 3 cont.
% Getting an empirical distribution of splineSSE and smoothing parameters
% for random rotations of the data.  The bottom line from this code is that
% the splineSSE depends alot on the coordinate frame in which the data are
% represented.  This is not good and steers me against using the thin plate
% spline.
niter = 2000;
scaled = data(Laccept,[1:3]) .*repmat(data(Laccept,4), 1,3);
out = nan*ones(niter,2);
for i = 1:niter
    rotdat = scaled(~Loog,:)*orth(normrnd(0,1,3,3))';

    r = sqrt(rotdat(:,1).^2+rotdat(:,2).^2+rotdat(:,3).^2);
    th = acos(rotdat(:,2)./r);
    ph = atan2(rotdat(:,3),rotdat(:,1));
    L = logical(ph<0);
    ph(L) = -ph(L);
    ph(~L) = 2*pi-ph(~L);

    [F,p] = tpaps([th ph]',log(r'),.99);
    resid = fnval(F,[th ph]')-log(r)';
    splineSSE = sum(resid.^2);
    out(i,:) = [splineSSE p];
end
figure;
subplot(2,2,1);
hist(out(:,1));
title('SplineSSE');
subplot(2,2,2);
hist(out(:,2));
title('Smoothing parameter');
subplot(2,2,3); hold on;
plot(out(:,1),out(:,2),'k.');
plot([planeSSE planeSSE],[min(out(:,2)) max(out(:,2))],'r-','linewidth',2);
xlabel('SplineSSE');
ylabel('Smoothing parameter');
%%
% Section 3 cont.
% 2) replotting in PCA space

M = pcs*diag(1./sqrt(d));
a = scaled*M;
cov(a(~Loog,:))  % Sanity check
lims = [min(a(~Loog,:)); max(a(~Loog,:))];
[x,y] = meshgrid(linspace(lims(1,1),lims(2,1),2),linspace(lims(1,2),lims(2,2),2));
b = pointertoplane'*pcs*diag(sqrt(d));
z=-(1+b(1)*x+b(2)*y)/b(3);

figure; axes; hold on;
plot3(a(~Loog,1),a(~Loog,2), a(~Loog,3),'r.','MarkerSize',16);
h=surf(x,y,z);
set(h,'FaceAlpha',.75,'EdgeAlpha',0,'FaceColor',[.5 .5 .5]);
%%
% Section 3 cont.
% Plotting in spherical coordinates
r = data(:,4);
th = acos(data(:,2));
ph = atan2(data(:,3),data(:,1));
L = logical(ph<0);
ph(L) = -ph(L);
ph(~L) = 2*pi-ph(~L);

figure; axes; hold on;
plot3(th(~Loog),ph(~Loog),log(r(~Loog)),'r.','MarkerSize',16);

% Now the surface
lims = [min(scaled(~Loog,[1 2])); max(scaled(~Loog,[1 2]))];
[x,y] = meshgrid(linspace(lims(1,1),lims(2,1),50),linspace(lims(1,2),lims(2,2),50));

a = pointertoplane(1);
b = pointertoplane(2);
c = pointertoplane(3);
z = -(1+a*x+b*y)/c;
r = sqrt(x.^2+y.^2+z.^2);
th = acos(y./r);
ph = atan2(z,x);
L = logical(ph<0);
ph(L) = -ph(L);
ph(~L) = 2*pi-ph(~L);
h = surf(th,ph,log(r));
set(h,'FaceAlpha',.5,'EdgeAlpha',0);
xlabel('theta');
ylabel('phi');
zlabel('r');


%%
% Section 3 cont.
% Residuals
r = data(:,4);
th = acos(data(:,2));
ph = atan2(data(:,3),data(:,1));
L = logical(ph<0);
ph(L) = -ph(L);
ph(~L) = 2*pi-ph(~L);

predr = -1./(pointertoplane(1).*sin(th).*cos(ph)+pointertoplane(2).*cos(th)-pointertoplane(3).*sin(th).*sin(ph));
figure;
subplot(2,1,1);
plot(predr(~Loog),r(~Loog)-predr(~Loog),'k.');
xlabel('fitted');
ylabel('data-fitted');

subplot(2,1,2);
plot(log(predr(~Loog)),log(r(~Loog))-log(predr(~Loog)),'k.');
xlabel('log(fitted)');
ylabel('log(data-fitted)');

%%
% Section 4
% Fitting a plane and a quadratic surface to simulated data. 
a = .2;
b = -.2;
c = -.1;

n = 40;
x = normrnd(0,1,n,1);
y = normrnd(0,1,n,1);
z = -(1+a*x+b*y)/c;
r = sqrt(x.^2+y.^2+z.^2);
r = r+normrnd(0,r/4);  % adding some noise.

th = acos(y./r);
ph = atan2(z,x);
L = logical(ph<0);
ph(L) = -ph(L);
ph(~L) = 2*pi-ph(~L);

% Bringing the noisy points back into xyz space
x = r.*sin(th).*cos(ph);
y = r.*cos(th);
z = -r.*sin(th).*sin(ph);

% Trying to fit the plane
% Creating an initial guess
mn = mean([x y z]);
initxyz = -pinv(mn);
initxyz = [a b c]';
options = optimset('MaxFunEvals',50000,'MaxIter',50000);
[planeparams, planeerr, out3, out4] = fminsearch(@(fn) surfacefiterr([x y z],fn),initxyz,options);
[quadparams1, quaderr1, out3, out4] = fminsearch(@(fn) surfacefiterr([x y z],fn),[planeparams;0;0;0],options);
[quadparams2, quaderr2, out3, out4] = fminsearch(@(fn) surfacefiterr([x y z],fn),[planeparams;0;0;0;0;0;0],options);

% Plotting the plane
figure; axes; hold on;
plot3(0,0,0,'y*','MarkerSize',10);
plot3(x,y,z,'r.','MarkerSize',10);
[xx yy] = meshgrid(linspace(-2,2,2),linspace(-2,2,2));
planez = -(1+planeparams(1)*xx+planeparams(2)*yy)/planeparams(3);
h = surf(xx,yy,planez);
set(h,'FaceAlpha',.2);


% Plotting the quadratic fit without the crossterms
% The quadratic surface: 
% ax+by+cz+dx^2+ey^2+fz^2+gxy+hyz+ixz+1 = 0

% Rearranging terms 
% fz^2+(c+hy+ix)z+(1+ax+by+dx^2+ey^2+gxy) = 0
% General solution to quadratic equation: az^2+bz+c = 0
% z = 2c/(-b+-4ac)
cprime = 1+quadparams1(1).*mn(1)+quadparams1(2).*mn(2)+quadparams1(4).*mn(1).^2+quadparams1(5).*mn(2).^2;
z1 = 2*cprime./(-quadparams1(3)+sqrt(quadparams1(3).^2-4*quadparams1(6)*cprime));
z2 = 2*cprime./(-quadparams1(3)-sqrt(quadparams1(3).^2-4*quadparams1(6)*cprime));
if (abs(mn(3)-z1) < abs(mn(3)-z2))
    signz = 1;
else
    signz = -1;
end

[xx yy] = meshgrid(linspace(-2,2,100),linspace(-2,2,100));
cprime = 1+quadparams(1).*xx+quadparams(2).*yy+quadparams(4).*xx.^2+quadparams(5).*yy.^2;
quadz = 2*cprime./(-quadparams(3)+signz.*sqrt(quadparams(3).^2-4*quadparams(6)*cprime));
Limag = real(quadz) ~= quadz;
xx(Limag) = nan;
yy(Limag) = nan;
quadz(Limag) = nan;
h = surf(xx,yy,quadz);
set(h,'FaceAlpha',.2,'EdgeAlpha',.2,'edgecolor','red');

% Figuring out which of the two solutions (surfaces) to use for the
% quadratic fit with the crossterms
cprime = 1+quadparams2(1).*mn(1)+quadparams2(2).*mn(2)+quadparams2(4).*mn(1).^2+quadparams2(5).*mn(2).^2+quadparams2(7).*mn(1).*mn(2);
bprime = quadparams2(3)+quadparams2(8).*mn(2)+quadparams2(9).*mn(1);
z1 = 2*cprime./(-bprime+sqrt(bprime.^2-4*quadparams2(6)*cprime));
z2 = 2*cprime./(-bprime-sqrt(bprime.^2-4*quadparams2(6)*cprime));
if (abs(mn(3)-z1) < abs(mn(3)-z2))
    signz = 1;
else
    signz = -1;
end
cprime = 1+quadparams2(1).*xx+quadparams2(2).*yy+quadparams2(4).*xx.^2+quadparams2(5).*yy.^2+quadparams2(7).*xx.*yy;
bprime = quadparams2(3)+quadparams2(8).*yy+quadparams2(9).*xx;

quadz = 2*cprime./(-bprime+signz.*sqrt(bprime.^2-4*quadparams2(6)*cprime));
Limag = real(quadz) ~= quadz;
xx(Limag) = nan;
yy(Limag) = nan;
quadz(Limag) = nan;
h = surf(xx,yy,quadz);
set(h,'FaceAlpha',.2,'EdgeAlpha',.2,'edgecolor','blue');

title(['SSEs: ',num2str(planeerr),' ',num2str(quaderr1),' ',num2str(quaderr2)]);
%%
% Section 5
% Now trying a bootstrapping test comparing a planar model to a quadratic
% model.
% Need to have "data" from NeuroThreshStuff.m and Laccept
% and you should have gone through the part in which a single "group"
% of points is identified.
% This is a good file: 'K082509007.nex'
% Still need a good way of making the initial guess.

% Fitting the plane and the quadratic surface
Loog = logical(data(Laccept,6));
scaled = data(Laccept,[1:3]) .*repmat(data(Laccept,4), 1,3);
mn = mean(scaled(~Loog,:))';
[pcs, score, d] = princomp(scaled(~Loog,:));
mn = mean(scaled(~Loog,:))';
initxyz = -sign(mn'*pcs(:,3))* pcs(:,3)./norm(mn'*pcs(:,3));
options = optimset('MaxFunEvals',10000,'MaxIter',8000);
[planeparams, planeSSE, out3, out4] = fminsearch(@(fn) surfacefiterr(scaled(~Loog,:),fn),initxyz,options);
[quadparams, quadSSE, out3, out4] = fminsearch(@(fn) surfacefiterr(scaled(~Loog,:),fn),[planeparams;0;0;0;0;0;0],options);

% Plotting data
figure; axes; hold on;
plot3(scaled(~Loog,1),scaled(~Loog,2),scaled(~Loog,3),'r.','MarkerSize',10);
plot3(scaled(Loog,1),scaled(Loog,2),scaled(Loog,3),'y.','MarkerSize',10);
plot3(0,0,0,'y*','MarkerSize',10);
title(num2str([planeSSE quadSSE]));

% Plotting plane
[xx yy] = meshgrid(linspace(-.1,.1,2),linspace(-.1,.1,2));
planez = -(1+planeparams(1)*xx+planeparams(2)*yy)/planeparams(3);
h = surf(xx,yy,planez);
set(h,'FaceAlpha',.2);

% Plotting quadratic surface
[xx yy] = meshgrid(linspace(-.2,.2,100),linspace(-.2,.2,100));
cprime = 1+quadparams(1).*xx+quadparams(2).*yy+quadparams(4).*xx.^2+quadparams(5).*yy.^2+quadparams(7).*xx.*yy;
bprime = quadparams(3)+quadparams(8).*yy+quadparams(9).*xx;
for signz = [-1 1];
    quadz = 2*cprime./(-bprime+signz.*sqrt(bprime.^2-4*quadparams(6)*cprime));
    Limag = real(quadz) ~= quadz;
    xx(Limag) = nan;
    yy(Limag) = nan;
    quadz(Limag) = nan;
    h = surf(xx,yy,quadz);
    set(h,'FaceAlpha',.2,'EdgeAlpha',.2,'edgecolor','blue');
end
% Getting the residuals from the quadratic fit
r = sqrt(sum(scaled(:,[1 2 3]).^2,2));

A = ([scaled(~Loog,[1 2 3]).^2,...
    scaled(~Loog,1).*scaled(~Loog,2),...
    scaled(~Loog,2).*scaled(~Loog,3),...
    scaled(~Loog,1).*scaled(~Loog,3)]...
        *quadparams([4 5 6 7 8 9]))./(r(~Loog,:).^2);

B = (scaled(~Loog,[1 2 3])*quadparams([1 2 3]))./r(~Loog,:);

predr1 = 2./(-B+sqrt(B.^2-4.*A));
predr2 = 2./(-B-sqrt(B.^2-4.*A));

% Clunky code below
predr = nan*ones(size(predr1));
predr(predr1 < 0) = predr2(predr1 < 0);
predr(predr2 < 0) = predr1(predr2 < 0);
L = isnan(predr); % color directions with two positive solutions
predr(L) = min([predr1(L), predr2(L)],[],2);
resid = log(r)-log(predr);

exp(resid)

% Finding r's that lie on the best fit plane.
th = acos(scaled(~Loog,2)./r(~Loog));
ph = atan2(scaled(~Loog,3),scaled(~Loog,1));
L = logical(ph<0);
ph(L) = -ph(L);
ph(~L) = 2*pi-ph(~L);
planer = -1./(planeparams(1).*sin(th).*cos(ph)+planeparams(2).*cos(th)-planeparams(3).*sin(th).*sin(ph))

niter = 10;
nulldist = nan*ones(niter,1);
drawnow;
for i = 1:niter
    i
    tmpresid = exp(resid(unidrnd(length(resid),[1 length(resid)])));
    l = planer.*tmpresid.*sin(th).*cos(ph);
    m = planer.*tmpresid.*cos(th);
    s = -planer.*tmpresid.*sin(th).*sin(ph);

    [tmp, SSE1, exitflag] = fminsearch(@(x) surfacefiterr([l m s],x),planeparams, options);
    [tmp, SSE2, exitflag] = fminsearch(@(x) surfacefiterr([l m s],x),[tmp;0;0;0;0;0;0], options);
    nulldist(i)= SSE1./SSE2;
end
p = sum(nulldist>planeSSE/quadSSE)./niter;
figure; axes; hold on;
hist(nulldist,40);
title(['p=',num2str(p)]);
plot(planeSSE/quadSSE,0,'m*');
%%
% Section 6
% Testing the alpha level of the bootstrap test
% Generating points on a plane plus some radial noise.
% Then doing the bootstrap test, comparing the qualities of plane to 
% quadratic fits to see how often we get a significant result.

% Simulated plane
niter = 1000;
nbootiter = 2000;
ps = [];
allplaneparams = [];
planeSSEs = [];
quadSSEs = [];
allquadparams = [];
types = [];
for i  = 1:niter
    % Inititial rotation of the data is probably important.
    % 'C' can't be too small. *(Is this really a problem?  I don't know.)
    a = normrnd(0,1)/10;
    b = normrnd(0,1)/10;
    c = normrnd(1,1)/10;

    % simulating data points
    n = 40;
    x = normrnd(0,1,n,1);
    y = normrnd(0,1,n,1);
    z = -(1+a*x+b*y)/c;
    r = sqrt(x.^2+y.^2+z.^2);
    noise = normrnd(1,.02,n,1);
    r = exp(log(r)+log(noise));

    th = acos(y./r);
    ph = atan2(z,x);
    L = logical(ph<0);
    ph(L) = -ph(L);
    ph(~L) = 2*pi-ph(~L);

    % Bringing the noisy points back into xyz space
    x = r.*sin(th).*cos(ph);
    y = r.*cos(th);
    z = -r.*sin(th).*sin(ph);

    % Trying to fit the plane
    % Creating an initial guess
    %mn = mean([x y z])';
    %[pcs, score, d] = princomp([x y z]);
    %initxyz = -sign(mn'*pcs(:,3))* pcs(:,3)./norm(mn'*pcs(:,3));
    initxyz = [a b c]';  % Using the true values as an initial guess.
    options = optimset('MaxFunEvals',50000,'MaxIter',50000);
    disp('Fitting plane');
    [planeparams, planeSSE, out3, out4] = fminsearch(@(fn) surfacefiterr2([x y z],fn),initxyz,options);
    disp('Fitting quadratic surface');
    %[quadparams, quadSSE, out3, out4] = fminsearch(@(fn) surfacefiterr([x y z],fn),[planeparams;0;0;0],options);
    [quadparams, quadSSE, out3, out4] = fminsearch(@(fn) surfacefiterr2([x y z],fn),[planeparams;0;0;0],options);
     
%     A = ([x.^2 y.^2 x.*y]*quadparams([4 5 6]))./(r.^2);
%     B = ([x y z]*quadparams([1 2 3]))./r;
% 
%     quadpredr1 = 2./(-B+sqrt(B.^2-4.*A));
%     quadpredr2 = 2./(-B-sqrt(B.^2-4.*A));
% 
%     resid = nan*ones(size(quadpredr1));
%     resid(quadpredr1 < 0) = log(r(quadpredr1 < 0))-log(quadpredr2(quadpredr1 < 0));
%     resid(quadpredr2 < 0) = log(r(quadpredr2 < 0))-log(quadpredr1(quadpredr2 < 0));
%     L = isnan(resid); % color directions with two positive solutions
%     resid(L) = min([abs(log(r(L))-log(quadpredr1(L))) abs(log(r(L))-log(quadpredr2(L)))],[],2);

    % Finding r's that lie on the best fit plane.
    th = acos(y./r);
    ph = atan2(z,x);
    L = logical(ph<0);
    ph(L) = -ph(L);
    ph(~L) = 2*pi-ph(~L);
    planer = -1./(planeparams(1).*sin(th).*cos(ph)+planeparams(2).*cos(th)-planeparams(3).*sin(th).*sin(ph));

    % Debugging
     resid = log(r)-log(planer);
    
    nulldist = nan*ones(nbootiter,1);
    wait_h = waitbar(0,'Bootstrapping...');
    disp(['Bootstrapping.  PlaneSSE = ',num2str(planeSSE),' QuadSSE = ',num2str(quadSSE)]);
    SSEs = [];
    for j = 1:nbootiter
        waitbar(j/nbootiter, wait_h);
        tmpresid = exp(resid(unidrnd(length(resid),[1 length(resid)])));
        
        % Debugging
        %tmpresid = normrnd(1,.02,n,1);
     
        l = planer.*tmpresid.*sin(th).*cos(ph);
        m = planer.*tmpresid.*cos(th);
        s = -planer.*tmpresid.*sin(th).*sin(ph);

        [tmp, SSE1, exitflag1] = fminsearch(@(x) surfacefiterr2([l m s],x),planeparams, options);
        % [tmp, SSE2, exitflag2] = fminsearch(@(x) surfacefiterr([l m s],x),[planeparams;0;0;0], options);
        [tmp, SSE2, exitflag2] = fminsearch(@(x) surfacefiterr2([l m s],x),[planeparams;0;0;0], options);
        if (~exitflag1)
            disp('Bonked on a plane fit');
        end
        if (~exitflag2)
            disp('Bonked on a quadratic fit');
        end
        if (~exitflag1 | ~exitflag2)
            j = j-1;
        end
        nulldist(j)= (SSE1-SSE2)./SSE2;
        SSEs(j,:) = [SSE1 SSE2];
    end
    close(wait_h);
    p = sum(nulldist>(planeSSE-quadSSE)/quadSSE)./nbootiter;
    planeSSEs = [planeSSEs; planeSSE mean(SSEs(:,1)) std(SSEs(:,1))];
    quadSSEs = [quadSSEs; quadSSE mean(SSEs(:,2)) std(SSEs(:,2))];
    ps = [ps; p]
    types = [types; std(resid)];
    allplaneparams = [allplaneparams; planeparams'];
    allquadparams = [allquadparams; quadparams'];
end

figure;
hist(ps);
figure;
fittdiff = sum((allplaneparams-allquadparams(:,[1 2 3])).^2,2);
hist(ps(fittdiff < 100));

%%
% Section 7
% OK, now something even simpler
% Taking data that are Gaussian deviations away from a line.  Fitting a
% line and a quadratic by OLS.  Resampling from the quadratic residuals,
% sticking them on the best fit line and calculating the SSE of the line
% and quadratic fits.

% The distribution of p-values generated by this code is impressively
% uniform!

% Simulated plane
niter = 1000;
nbootiter = 2000;
ps = [];
for i  = 1:niter
    b = normrnd(0,1,2,1)/10;

    % simulating data points
    n = 40;
    x = normrnd(0,1,n,1);
    y = b(1)+b(2)*x+normrnd(0,.1,n,1);
    b_linear = [ones(n,1) x]\y;
    yhat_linear = [ones(n,1) x]*b_linear;
    lineSSE = sum((y-yhat_linear).^2);
    b_quad = [ones(n,1) x x.^2]\y;
    yhat_quad = [ones(n,1) x x.^2]*b_quad;
    quadSSE = sum((y-yhat_quad).^2);
    resid = y-yhat_quad;
    
    nulldist = nan*ones(nbootiter,1);
    wait_h = waitbar(0,'Bootstrapping...');
    disp(['Bootstrapping.  PlaneSSE = ',num2str(lineSSE),' QuadSSE = ',num2str(quadSSE)]);
    for j = 1:nbootiter
        waitbar(j/nbootiter, wait_h);
        tmpresid = resid(unidrnd(length(resid),[1 length(resid)]));
        y = yhat_linear+tmpresid;
        
        b_linear = [ones(n,1) x]\y;
        yhat_linear = [ones(n,1) x]*b_linear;
        SSE1 = sum((y-yhat_linear).^2);
        b_quad = [ones(n,1) x x.^2]\y;
        yhat_quad = [ones(n,1) x x.^2]*b_quad;
        SSE2 = sum((y-yhat_quad).^2);
        
        nulldist(j)= SSE1./SSE2;
    end
    close(wait_h);
    p = sum(nulldist>lineSSE/quadSSE)./nbootiter;
    ps = [ps; p]
 end

figure;
hist(ps);

%%
% Section 8
% trying to get a meaningful (if model-bound) estimate of neurometric
% threshold from NeuroThresh data. (Run the begining part of
% NeuroThreshStuff first so that we have "planeparams" and "quadparams").
% Not assuming a linear contrast response function.

% Getting raw spike counts.
nspikes = zeros(size(stro.trial,1),1);
for i = 1:size(stro.trial,1)
    spiketimes = stro.ras{i,spikeidx};
    nspikes(i) = sum(spiketimes > stimon_t(i)+stro.sum.exptParams.latency/1000 & spiketimes < stimoff_t(i));
    dur(i)= stimoff_t(i)-stimon_t(i)-stro.sum.exptParams.latency/1000;
end

tmp1 = [];
tmp2 = [];
for i = find(Laccept)'
    if (data(i,end) == 1)
        continue
    end
    % Getting predicted distances to plane (in transformed space)
    coloridx = uniquecoloridxs(i);
    colordir = data(i,[1 2 3]);
    rotcol =  colordir*xformmat; % To get into the same space as planeparams/quadparams  
    
    r = sqrt(sum(rotcol.^2,2));
    A = (quadparams(4).*rotcol(:,1).^2+quadparams(5).*rotcol(:,2).^2+quadparams(6).*rotcol(:,1).*rotcol(:,2))./(r.^2);
    B = (quadparams(1).*rotcol(:,1)+quadparams(2).*rotcol(:,2)+quadparams(3).*rotcol(:,3))./r;

    quadpredr1 = 2./(-B+sqrt(B.^2-4.*A));
    quadpredr2 = 2./(-B-sqrt(B.^2-4.*A));
    quadpredr1(real(quadpredr1) ~= quadpredr1) = Inf;
    quadpredr2(real(quadpredr2) ~= quadpredr2) = Inf;
    quadpredr = min([abs(quadpredr1), abs(quadpredr2)]);
    quadpredr = quadpredr./r;  
    % Quadpredr is now multiplier to get a unit vector in color direction to surface
    
    linpredr = (planeparams'*planeparams)/abs(rotcol*planeparams); % distance to plane (multiplier)
    % (predr*rotcol)*planeparams == norm(planeparams).^2
    % predr is the scaling factor that puts rotcol on the plane
    tmp1 = [tmp1;quadpredr data(i,4)];

    L = coloridxs == coloridx;
    lms = stro.trial(L,lmsidxs);
    [u,s,v] = svd(lms);
    unitvector = v(:,1);
    contrasts = abs(lms*unitvector);
    resp = nspikes(L);
    tmp2 = [tmp2; repmat(coloridx,sum(L),1) contrasts/quadpredr, resp];
   
  % color_id, normalized contrast, response
end
figure; axes; hold on;
plot(log(tmp1(:,1)),log(tmp1(:,2)),'k.');
plot([min(log(tmp1(:,1))),max(log(tmp1(:,1)))],[min(log(tmp1(:,2))),max(log(tmp1(:,2)))],'k:');
xlabel('log predr');
ylabel('log r');

figure; axes; hold on;
for i = unique(tmp2(:,1))'
    L = tmp2(:,1) == i;
    h = plot(tmp2(L,2), tmp2(L,3),'k.');
    set(h,'Color',unifrnd(0,1,1,3));
   % drawnow;
   % pause;
   % cla
end
beta = glmfit([tmp2(:,2) tmp2(:,2).^2],tmp2(:,3),'poisson','link','identity');
%plot([0 2],beta(1)+beta(2)*[0 2],'k-');
tmp = logspace(-.5, .3);
plot(tmp,beta(1)+beta(2)*tmp+beta(3)*tmp.^2,'k-');
thr = stro.sum.exptParams.threshold*mean(dur);
plot([min(tmp2(:,2)) 1 1],[thr thr 0],'k-');
set(gca,'xscale','log')

% Now have to find the "normalized contrast" at which the ROC area is 0.82
% Getting baseline response
 L = logical(~isnan(fpacq_t));
 baseline = [];
 for i = find(L)'
     spiketimes = stro.ras{i,spikeidx};
     baseline = [baseline; sum(spiketimes > stimon_t(i)-stro.sum.exptParams.latency/1000 & spiketimes < stimon_t(i))];
 end
 [baselinepdf,n] = hist(baseline,[0:5*max(nspikes)]);
 baselinepdf = baselinepdf./sum(baselinepdf);
 %baselinepdf = poisspdf(0:5*max(nspikes), mean(baseline));
 %baselinepdf = poisspdf(0:5*max(nspikes), beta(1));
 tmp = [];
 for testcontrast = 0:.01:2
     testlambda = beta(1)+beta(2)*testcontrast+beta(3)*testcontrast.^2;
     testpdf = poisspdf(0:5*max(nspikes),testlambda);
     testcdf = poisscdf(0:5*max(nspikes),testlambda);
     rocarea = baselinepdf*(1-testcdf') + .5* testpdf*baselinepdf';
    tmp  = [tmp; testcontrast, rocarea];
    % Logic for the above line:
    % P(X>Y) = P(X>0|y=1)P(Y=1) + ...
	% P(X=Y) = P(X=0|y=0)P(Y=0) + ...
end
err = (tmp(:,2) - .82).^2;
% Normalized contrast (0 = blank, 1 = on isoresponse surface) 
% at which modeled neurometric threshold is 0.82
threshnormcontrast = tmp(find(err == min(err)),1);
% Now estimating threshold for a particular color direction in cone
% contrast units
colordir = [0 0 1];
colordir = colordir./norm(colordir);
rotcol =  colordir*xformmat; % To get into the same space as planeparams
%predr = (planeparams'*planeparams)/abs(rotcol*planeparams); % linear fit
r = sqrt(sum(rotcol.^2,2));
A = (quadparams(4).*rotcol(:,1).^2+quadparams(5).*rotcol(:,2).^2+quadparams(6).*rotcol(:,1).*rotcol(:,2))./(r.^2);
B = (quadparams(1).*rotcol(:,1)+quadparams(2).*rotcol(:,2)+quadparams(3).*rotcol(:,3))./r;

quadpredr1 = 2./(-B+sqrt(B.^2-4.*A));
quadpredr2 = 2./(-B-sqrt(B.^2-4.*A));
quadpredr1(real(quadpredr1) ~= quadpredr1) = Inf;
quadpredr2(real(quadpredr2) ~= quadpredr2) = Inf;
quadpredr = min([abs(quadpredr1), abs(quadpredr2)]);
quadpredr = quadpredr./r;  % Quadpredr is now multiplier to get a unit vector in color dierciton to surface

predthresh = quadpredr*threshnormcontrast;
title(['Color Dir: [',num2str(colordir,2),'] Estimated threshold: ',num2str(predthresh)])

%%
% Testing out parametric ROC calculation
rocs = [];
for i = 1:1000
    x = poissrnd(mean(baseline), 100,1);
    y = poissrnd(testnormcontrast, 100,1);
    rocs = [rocs; roc(x,y)];
end
hist(rocs)
title(num2str(mean(rocs)))

%%
% Section 8.1
% Comparing contrast-response functions in multiple directions in color
% space in NeuroThresh.  How good is the approximation of a single
% underlying contrast-response function that works in all color directions
% but the domain of which is scaled?
NTfilename = 'S030510003';
NT = nex2stro(findfile(NTfilename));
out = NTpreprocess(NT,.4, 1);
out(out(:,end) == 1,:) = [];  % Getting rid of the OOGs
spikeidx = strcmp(NT.sum.rasterCells(1,:),getSpikenum(NT));
coloridxs = NT.trial(:,strcmp(NT.sum.trialFields(1,:),'coloridx'));
stimon_t = NT.trial(:,strcmp(NT.sum.trialFields(1,:),'stim_on'));
stimoff_t = NT.trial(:,strcmp(NT.sum.trialFields(1,:),'stim_off'));
lms = NT.trial(:,[find(strcmp(NT.sum.trialFields(1,:),'lcont'))...
            find(strcmp(NT.sum.trialFields(1,:),'mcont'))...
            find(strcmp(NT.sum.trialFields(1,:),'scont'))]);

for i = 1:size(NT.trial,1)
    spiketimes = NT.ras{i,spikeidx};
    nspikes(i) = sum(spiketimes > stimon_t(i)+NT.sum.exptParams.latency/1000 & spiketimes < stimoff_t(i));
end

uniquecoloridxs = out(:,1);
uniquecolordirs = out(:,[2:4]);

% First, getting the color directions, contrasts, and responses into an
% easy to use format.  Also plotting contrast response functions in vector
% norms;
figure; subplot(2,1,1); hold on;
data = [];
for i = uniquecoloridxs'
    L = coloridxs == i;
    contrast = abs(uniquecolordirs(uniquecoloridxs == i,:)*lms(L,:)');
    data(length(data)+1).contrast = contrast';
    data(length(data)).response = nspikes(L);
    data(length(data)).lms = uniquecolordirs(uniquecoloridxs == i,:);
    
    h = plot(contrast, nspikes(L),'ko');
    set(h,'MarkerFaceColor',unifrnd(0,1,3,1));
end
set(gca,'XScale','log');
[beta, dev, stats] = glmfit(cat(1,data.contrast),cat(1,data.response),'poisson');
subplot(2,1,2);
plot(cat(1,data.contrast), stats.resid,'k.');
set(gca,'XScale','log');

% Now, fitting contrast-response curves (lines for now) separately for each
% color direction.  No constant because it has to be the same for all color
% directions.
beta = [];
for i = 1:length(data)
    beta(i) = glmfit(data(i).contrast,data(i).response,'poisson','constant','off');
end
figure; axes; hold on;
for i = 1:length(data)
   h = plot(beta(i)*data(i).contrast,data(i).response,'ko');
   set(h,'MarkerFaceColor',unifrnd(0,1,3,1));
end
set(gca,'XScale','log','YScale','log');


% Now fitting a gigantic model to the whole dataset
% with different slopes for each color direction but a common y-intercept.
r = cat(1,data.response);
designmatrix = zeros(length(r),length(data));
for i = 1:length(data)
    if i == 1
        tmp = [1 length(data(i).response)];
    else
        tmp = length(cat(1,data(1:i-1).response))+[1 length(data(i).response)];
    end
    designmatrix(tmp(1):tmp(2),i) = data(i).contrast;
end
[beta,dev,stats] = glmfit(designmatrix,r,'poisson');
figure; subplot(2,1,1); hold on;
yhat = [];
for i = 1:length(data)
   h = plot(beta(i+1)*data(i).contrast+beta(1),data(i).response,'ko');
   set(h,'MarkerFaceColor',unifrnd(0,1,3,1));
   yhat= [yhat; repmat(i,size(data(i).contrast,1),1), beta(i+1)*data(i).contrast+beta(1)];
end
set(gca,'XScale','log');

subplot(2,1,2); hold on;
for i = 1:length(data)
    L = yhat(:,1) == i;
    h = plot(yhat(L,2),stats.resid(L),'ko');
    set(h,'MarkerFaceColor',unifrnd(0,1,3,1));
end
set(gca,'XScale','log');

% Plotting beta as a function of color direction
scaled = uniquecolordirs./repmat(beta(2:end,:),1,3);
figure; axes; hold on;
plot3(scaled(:,1),scaled(:,2),scaled(:,3),'k.')
plot3(-scaled(:,1),-scaled(:,2),-scaled(:,3),'k.')
plot3(0,0,0,'y.')
%%
% Section 8.2
% Comparing contrast-response functions measured in DTspot and NeuroThresh.
DTfilename = 'S012710006';
NTfilename = 'S012710004';
DT = nex2stro(findfile(DTfilename));
NT = nex2stro(findfile(NTfilename));

% First getting contrast-response function from DTspot
ntrials = size(DT.trial,1);
DTcolordirs = reshape(DT.sum.exptParams.RF_colors,[3 3]);
stimon_t = DT.trial(:,strcmp(DT.sum.trialFields(1,:),'flash_on'))
stimoff_t = DT.trial(:,strcmp(DT.sum.trialFields(1,:),'flash_off'))
dt = stimoff_t-stimon_t;
DTspikerate = zeros(ntrials,1);
for i = 1:ntrials
    spiketimes = DT.ras{i,strcmp(DT.sum.rasterCells,'sig001a')};
    nspikes = sum(spiketimes > stimon_t(i) & spiketimes < stimoff_t(i)) 
    DTspikerate(i) = nspikes./dt(i)
end




%%
% Section 9
% How does the STA compare with the isoresponse surface fit by NeuroThresh?

WNfile = 'K072009001'; GTfile = 'K072009002'; NTfile = 'K072009003'; % Lum simple cell
WNfile = 'K090309001'; GTfile = 'K090309003'; NTfile = 'K090309006';%Pancolor complex
%WNfile = 'K102209008'; GTfile = 'K102209006'; NTfile = 'K102209007';% lum
%WNfile = 'K010110001'; GTfile = 'K010110003'; NTfile = 'K010110005'; %Nice lum simple cell
%WNfile = 'S012710001'; GTfile = 'S012710002'; NTfile = 'S012710003';%color
%WNfile = 'S033010001'; GTfile = 'S033010002'; NTfile = 'S033010003';%blue
WNfile = 'S071410006'; GTfile = 'S071410004'; NTfile = 'S071410005';% DO - discrep.
WNfile = 'S021810011'; GTfile = 'S021810007'; NTfile = 'S021810008';% 

%WNfile = 'S030410001'; GTfile = 'S030410002'; NTfile = 'S030410003';%magenta pan


WN=nex2stro(findfile(WNfile));
GT=nex2stro(findfile(GTfile));
NT=nex2stro(findfile(NTfile));
GTout = getGratingTuning(GT, 1);
GTconeweights = GTout.color.prefcolor/sum(abs(GTout.color.prefcolor));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Getting cone weights from white noise
framerate = WN.sum.exptParams.framerate;
nstixperside = WN.sum.exptParams.nstixperside;
ntrials = length(WN.sum.absTrialNum);
stimonidx = find(strcmp(WN.sum.trialFields(1,:),'stim_on'));
stimoffidx = find(strcmp(WN.sum.trialFields(1,:),'all_off'));
nframesidx = find(strcmp(WN.sum.trialFields(1,:),'num_frames'));
noisetypeidx = find(strcmp(WN.sum.trialFields(1,:),'noise_type'));
sigmaidxs = strmatch('sigma',WN.sum.trialFields(1,:));

hepidx = find(strcmp(WN.sum.rasterCells(1,:),'AD11'));
vepidx = find(strcmp(WN.sum.rasterCells(1,:),'AD12'));
anlgStartTimeidx = find(strcmp(WN.sum.rasterCells(1,:),'anlgStartTime'));
eyestart_t = [WN.ras{:,anlgStartTimeidx}]';
eyesampperiod = 1/WN.sum.analog.storeRates{1};

gammaTable = WN.sum.exptParams.gamma_table;
gammaTable = reshape(gammaTable, length(gammaTable)/3, 3);
gammaTable1 = interp1(linspace(0,255,256),gammaTable,linspace(0,255,65536), 'spline');
invgamma = InvertGamma(gammaTable, 0);

% Reconstructing the M matrix
fundamentals = WN.sum.exptParams.fundamentals;
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = WN.sum.exptParams.mon_spd;
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;

% Getting the background rgb/lms
ridx = find(strcmp(WN.sum.trialFields(1,:),'bkgnd_r'));
gidx = find(strcmp(WN.sum.trialFields(1,:),'bkgnd_g'));
bidx = find(strcmp(WN.sum.trialFields(1,:),'bkgnd_b'));
bkgndRGB = [mode(WN.trial(:,ridx)), mode(WN.trial(:,gidx)), mode(WN.trial(:,bidx))];
bkgndrgb = [gammaTable(bkgndRGB(1)+1,1); gammaTable(bkgndRGB(2)+1,2); gammaTable(bkgndRGB(3)+1,3)];
bkgndlms = M*bkgndrgb;

spikename = getSpikenum(WN);
spikeidx = find(strcmp(WN.sum.rasterCells(1,:),spikename));
maxT = 9;
Lgunnoise = WN.trial(:,noisetypeidx) == 1;
Lconenoise = WN.trial(:,noisetypeidx) == 2;
WN.ras(Lconenoise,:) = [];
WN.trial(Lconenoise,:) = [];
out = getWhtnsStats(WN,maxT,'STCOVmex',{nstixperside^2, 3, maxT},spikename);
STAs = out{1};
nspikes = out{3};
normfactor = 1./(2.5*max(abs(STAs(:))));
muvect = reshape(repmat([.5 .5 .5],nstixperside^2,1),nstixperside^2*3,1);

% Plotting
figure;
for i = 1:size(STAs,2)
    STA = normfactor*(STAs(:,i)-muvect)+muvect;
    STA = reshape(STA,[nstixperside nstixperside 3]);
    subplot(1,size(STAs,2),i);
    image(STA);
    set(gca,'XTick',[],'YTick',[]); axis square;
end

STAgunmat = reshape(STAs,[nstixperside^2 3 maxT]);
STAgunmat = permute(STAgunmat,[2 1 3]);
STAgunmat = reshape(STAgunmat,[3, nstixperside^2*maxT]);

energy = sum(STAgunmat.^2);
whichpix = logical(energy>3*std(energy));
[u,s,v] = svd(STAgunmat(:,whichpix));
u = mean(STAgunmat(:,whichpix),2)
%u(:,1) = STAgunmat(:,energy == max(energy)); % Single hottest pixel
if (sum(v(:,1)) < 0)
    u = -u;
end
WNgunweights = u(:,1)./sum(abs(u(:,1)));
WNconeweights = inv((diag(1./bkgndlms)*M)')*u(:,1);
% Do I have to do something about bkgndlms here, above?
WNconeweights = WNconeweights'./sum(abs(WNconeweights));
if (WNconeweights*GTconeweights' < 0)
    WNconeweights = -WNconeweights;
    WNgunweights = -WNgunweights;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Getting cone weights from NeuroThresh
NTout = NTpreprocess(NT,.4,2);
NTout(:,1) = [];
Laccept = true(size(NTout,1),1);

scaled = NTout(:,[1:3]) .*repmat(NTout(:,4), 1,3);
Loog = logical(NTout(:,6));

[v,d] = eig(cov([scaled(~Loog,:); -scaled(~Loog,:)]));
d = diag(d);
whtmat = diag(sqrt(1./d))*v';
scaled = scaled*whtmat';
[th,ph,r] = cart2sph(scaled(~Loog,1),scaled(~Loog,2),scaled(~Loog,3));

errs = zeros(500,500);
tmp = tan(linspace(-pi/2.001,pi/2.001,size(errs,1)));
for i = 1:length(tmp)
    for j = 1:length(tmp)
        % Parametrically varying a and b and leaving c=1
        predr = -1./(tmp(i).*cos(ph).*cos(th)+tmp(j).*cos(ph).*sin(th)+sin(ph));
        resid = log(abs(predr))-log(r);
        resid = resid-mean(resid);  % That extra degree of freedom
        errs(i,j) = sum(resid.^2);
    end
end
[i,j] = ind2sub(size(errs),find(errs(:) == min(errs(:))));
a = tmp(i); b = tmp(j); c = 1;
predr = -1./(tmp(i).*cos(ph).*cos(th)+tmp(j).*cos(ph).*sin(th)+sin(ph));
shiftfactor = mean(log(abs(predr))-log(r));
a = a*exp(shiftfactor);
b = b*exp(shiftfactor);
c = c*exp(shiftfactor);

% Changing coordinate frames for gradient descent fits of planar and
% quadratic surfaces.
rotmat = MakeOrtho([[a; b; c;] normrnd(0,1,3,2)]);
rotmat = [rotmat(:,[2 3]) rotmat(:,1)];
xyz = scaled(:,[1 2 3])*rotmat;
rotabc = inv(rotmat)*[a;b;c];

% Doing gradient descent for planar fit.
options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-6,'TolX',10^-6);
[planeparams, planeSSE, exitflag, out4] = fminsearch(@(x) surfacefiterr2(xyz, x, Loog),rotabc,options);
[quadparams, quadSSE, exitflag, out4] = fminsearch(@(x) surfacefiterr2(xyz, x, Loog),[rotabc;0;0;0],options);
NTconeweights = planeparams'*rotmat'*whtmat;
NTconeweights = NTconeweights./sum(abs(NTconeweights));
if (NTconeweights*GTconeweights' < 0)
    NTconeweights = -NTconeweights;
end

GTconeweights
WNconeweights
NTconeweights

% Plotting isoresponse plane and stimulus cloud from white noise.
scaled = NTout(:,[1:3]) .*repmat(NTout(:,4), 1,3);
lim = 6;
[x y] = meshgrid(linspace(-lim,lim,10),linspace(-lim,lim,10));
z = -(1+planeparams(1)*x+planeparams(2)*y)/planeparams(3);
%z = -(1+quadparams(1)*x+quadparams(2)*y+quadparams(4)*x.^2+quadparams(5)*y.^2+quadparams(6).*x.*y)/quadparams(3);

tmp = permute(cat(3,x,y,z),[3 1 2]);
tmp = reshape(tmp,3,size(x,2)*size(y,2));
xformM = inv(rotmat'*whtmat);

xformed = xformM*tmp;
xformed = reshape(xformed, 3, size(x,2), size(y,2));
xformed = permute(xformed,[2 3 1]);

figure; hold on;
color = [.5 .5 .5];
h1 = surf(xformed(:,:,1),xformed(:,:,2),xformed(:,:,3));
set(h1,'CDataMapping','direct','FaceColor',color,'EdgeColor',color,'EdgeAlpha',.2,'FaceAlpha',.2);
h2 = surf(-xformed(:,:,1),-xformed(:,:,2),-xformed(:,:,3));
set(h2,'CDataMapping','direct','FaceColor',color,'EdgeColor',color,'EdgeAlpha',.2,'FaceAlpha',.2);
xlabel('L'); ylabel('M'); zlabel('S');
plot3(scaled(~Loog,1),scaled(~Loog,2),scaled(~Loog,3),'k.')
plot3(-scaled(~Loog,1),-scaled(~Loog,2),-scaled(~Loog,3),'k.')
plot3(scaled(Loog,1),scaled(Loog,2),scaled(Loog,3),'y.')
plot3(-scaled(Loog,1),-scaled(Loog,2),-scaled(Loog,3),'y.')

% Now the stimulus cloud
sigma_rgb = [unique(WN.trial(:,strcmp(WN.sum.trialFields(1,:), 'sigma1'))),...
        unique(WN.trial(:,strcmp(WN.sum.trialFields(1,:), 'sigma2'))),...
        unique(WN.trial(:,strcmp(WN.sum.trialFields(1,:), 'sigma3')))]/1000;
sigma_rgb(all(sigma_rgb == 0,2),:) = [];
cov_rgb = diag(sigma_rgb.^2);
cov_lms = M*cov_rgb*M';
[v,d] = eig(cov_lms);
[x,y,z] = ellipsoid(0,0,0,sqrt(d(1,1))./bkgndlms(1),sqrt(d(2,2))./bkgndlms(2),sqrt(d(3,3))./bkgndlms(2));
tmpmat = v*[x(:) y(:) z(:)]';
newxyz = reshape(tmpmat',[size(x,1),size(x,2),3]);
h = surf(newxyz(:,:,1),newxyz(:,:,2),newxyz(:,:,3));
set(h,'EdgeAlpha',0,'FaceAlpha',.3);

h(1) = plot3([0 WNconeweights(1)],[0 WNconeweights(2)],[0 WNconeweights(3)],'r-');
h(2) = plot3([0 GTconeweights(1)],[0 GTconeweights(2)],[0 GTconeweights(3)],'k-');
h(3) = plot3([0 NTconeweights(1)],[0 NTconeweights(2)],[0 NTconeweights(3)],'m-');
legend(h,{'WN','GT','NT'})

% Now we do the same thing, but plot everything in RGB space
figure; axes; hold on;
[x,y,z] = ellipsoid(0,0,0,sigma_rgb(1),sigma_rgb(2),sigma_rgb(3));
h = surf(x,y,z);
set(h,'EdgeAlpha',0,'FaceAlpha',.3);
% NT.  Need to figure out how to deal with bkgndlms
rgb = (inv(M)*(scaled.*repmat(bkgndlms',size(scaled,1),1))')';
plot3(rgb(~Loog,1),rgb(~Loog,2),rgb(~Loog,3),'k.')
plot3(-rgb(~Loog,1),-rgb(~Loog,2),-rgb(~Loog,3),'k.')
plot3(rgb(Loog,1),rgb(Loog,2),rgb(Loog,3),'y.')
plot3(-rgb(Loog,1),-rgb(Loog,2),-rgb(Loog,3),'y.')
axis equal;
plot3([0 WNgunweights(1)],[0 WNgunweights(2)],[0 WNgunweights(3)],'r-');
NTgunweights = M'*NTconeweights';
NTgunweights = NTgunweights./sum(abs(NTgunweights));
plot3([0 NTgunweights(1)],[0 NTgunweights(2)],[0 NTgunweights(3)],'m-');

GTgunweights = M'*GTconeweights';
GTgunweights = GTgunweights./sum(abs(GTgunweights));
plot3([0 GTgunweights(1)],[0 GTgunweights(2)],[0 GTgunweights(3)],'k-');

set(gca,'Xlim',[-.4 .4],'YLim',[-.4 .4],'Zlim',[-.4 .4])

lim = 1;
[x y] = meshgrid(linspace(-lim,lim,10),linspace(-lim,lim,10));
z = -(1+planeparams(1)*x+planeparams(2)*y)/planeparams(3);
%z = -(1+quadparams(1)*x+quadparams(2)*y+quadparams(4)*x.^2+quadparams(5)*y.^2+quadparams(6).*x.*y)/quadparams(3);

tmp = permute(cat(3,x,y,z),[3 1 2]);
tmp = reshape(tmp,3,size(x,2)*size(y,2));
xformM = inv(diag(1./bkgndlms)*M)*inv(rotmat'*whtmat);

xformed = xformM*tmp;
xformed = reshape(xformed, 3, size(x,2), size(y,2));
xformed = permute(xformed,[2 3 1]);

color = [.5 .5 .5];
h1 = surf(xformed(:,:,1),xformed(:,:,2),xformed(:,:,3));
set(h1,'CDataMapping','direct','FaceColor',color,'EdgeColor',color,'EdgeAlpha',.2,'FaceAlpha',.2);
h2 = surf(-xformed(:,:,1),-xformed(:,:,2),-xformed(:,:,3));
set(h2,'CDataMapping','direct','FaceColor',color,'EdgeColor',color,'EdgeAlpha',.2,'FaceAlpha',.2);

%%
% Section 10
% Playing around with a few models of color tuning.  One hope is to fit
% models to the NeuroThresh data that use *all* the data - not just the
% terminations of the staircases.  Can we fit models that are separable in
% contrast and color direction?  Is this misguided?
figure;
weights = [1 1 1]';
[x,y,z] = meshgrid([-1:.1:1],[-1:.1:1],[-1:.1:1]);
%stim = [x(:) y(:) z(:)];
[th,ph,r] = cart2sph(x(:),y(:),z(:));
dotprods = [cos(ph).*cos(th) cos(ph).*sin(th) sin(ph)]*weights;
%dotprods = [x(:) y(:) z(:)]*weights;

%plot3(x(:),y(:),z(:),'.');
pow = .5
colors = {'red','green'};
vals = [.5 1];
for i = 1:length(vals)
    val = vals(i);
    fv = isosurface(x,y,z,reshape(r.^pow.*dotprods,size(x)),val)
    h = patch(fv);
    set(h,'FaceColor',colors{i},'EdgeColor',colors{i});
    fv = isosurface(x,y,z,reshape(r.^pow.*dotprods,size(x)),-val)
    h = patch(fv);
    set(h,'FaceColor',colors{i},'EdgeColor',colors{i});
end

% The model:
% r.^n*[cos(ph)cos(th) cos(ph)sin(th) sin(ph)]*[a b c]' 
% includes a plane and some bent surfaces.  It appears to be good for
% modeling a single isoresponse surface, but the contrast response
% functions get weird as we change 'n'.  

% GLM.  How do I express a GLM in spherical coordinates?
%%
% Doing everything in 2-D for now
% Hoping to stumble onto a good parametrization for a full model by
% flailing around in polar coordinates.  Basically unsuccessful.

weights = [.1 .2]';
[x,y] = meshgrid([-1:.1:1],[-1:.1:1]);
[th,r] = cart2pol(x(:),y(:));

pow1 = 1;
pow2 = 1;
dotprods = [cos(th) sin(th)].^pow1*weights;
dotprods = [cos(th).^2 sin(th).^2 2*cos(th).*sin(th)]*[.1 .1 0]';
fr = reshape(r.^pow2.*dotprods,size(x));

%fr = exp(fr);
%fr = 1-.3.*fr+2.*max(fr,0).^2
figure; subplot(2,1,1);
contour(x,y,fr); axis square;
%set(gca,'XScale','log','YScale','log')
subplot(2,1,2);
surf(x,y,fr); axis square;
%set(gca,'XScale','log','YScale','log')


%%
% Can we even *make* a model neuron that has identical contrast response
% functions in all directions (up to a horizontal scale factor) an smooth,
% nonplanar isoresponse contours?  Are they nested?


% Generic quadric surface 
% ax2+by2+cz2+dxy+eyz+fxz+gx+hy+iz = 1;
%coefficients = [1 0 0 0 0  0 0 0 0]; Gives parallel planes.


%coefficients = [1/2 1/2 0  1 0 0  0 0 0];
%[x,y,z] = meshgrid([-1:.2:1],[-1:.2:1],[-1:.2:1]);
%variables = [x(:).^2 y(:).^2 z(:).^2 x(:).*y(:) y(:).*z(:) x(:).*z(:) x(:) y(:) z(:)];
%fr = variables*coefficients';
%cla; hold on;
%isosurface(x,y,z,reshape(fr,size(x)), 1)
%isosurface(-x,-y,-z,reshape(fr,size(x)), 1)
%isosurface(x,y,z,reshape(fr,size(x)), .1)
%isosurface(-x,-y,-z,reshape(fr,size(x)), .5)
%isosurface(x,y,z,reshape(fr,size(x)), 0)
%set(gca,'XLim',[-1 1],'YLim',[-1 1],'ZLim',[-1 1]);

% 2-D case, pure rotation (norm v = norm u = 1);
% alpha can be from -1 to 1;
alpha = 0;
[x,y] = meshgrid([-1:.1:1],[-1:.1:1]);
variables = [x(:).^2 y(:).^2 x(:).*y(:)];
coefficients = [alpha.^2 1-alpha.^2 2*alpha*sqrt(1-alpha.^2)];
fr = variables*coefficients';
cla; hold on;
contour(x,y,reshape(fr,size(x)));
axis square

% 2-D case with the freedom to change the scale of x and y
%alpha = 1; beta = -3;
%coefficients = [alpha.^2 beta.^2 2*alpha*beta];
%fr = variables*coefficients';
%figure; cla; hold on;
%contour(x,y,reshape(fr,size(x)));
%axis square


%%
% 3-D case
alpha = 0;
beta = 0;
gamma = 1-alpha^2-beta^2;
[x,y,z] = meshgrid([-1:.1:1],[-1:.1:1],[-1:.1:1]);
coefficients = [alpha^2 beta^2 gamma^2 2*alpha*beta 2*alpha*gamma 2*beta*gamma];
variables = [x(:).^2 y(:).^2 z(:).^2 x(:).*y(:) x(:).*z(:) y(:).*z(:)];
fr = variables*coefficients';
cla; hold on;
isosurface(x,y,z,reshape(fr,size(x)), .5)
set(gca,'XLim',[-1 1],'YLim',[-1 1],'ZLim',[-1 1]);

%%
% Converting a pair of planes represented as |ax+by+cz| = 1 to the 
% ax^2+by^2+cz^2+dxy+exz+fyz = 1 representation

abc = [-1 1 -2];
[x,y,z] = meshgrid([-2:.2:2],[-2:.2:2],[-2:.2:2]);
variables = [x(:) y(:) z(:)];
fr = variables*abc';

% The first way
figure; subplot(2,1,1); hold on;
p = patch(isosurface(x,y,z,reshape(fr,size(x)), 1));
set(p, 'FaceColor', 'red', 'EdgeColor', 'none');
p = patch(isosurface(-x,-y,-z,reshape(fr,size(x)), 1));
set(gca,'XLim',[-2 2],'YLim',[-2 2],'ZLim',[-2 2]);
set(p, 'FaceColor', 'red', 'EdgeColor', 'none');
axis square;

% The second way
alpha = abc(1); beta = abc(2); gamma = abc(3);
coefficients = [alpha^2 beta^2 gamma^2 2*alpha*beta 2*alpha*gamma 2*beta*gamma];
variables = [x(:).^2 y(:).^2 z(:).^2 x(:).*y(:) x(:).*z(:) y(:).*z(:)];
fr = variables*coefficients';
%subplot(2,1,2);
p = patch(isosurface(x,y,z,reshape(fr,size(x)), 1))
set(p, 'FaceColor', 'blue', 'EdgeColor', 'none');
set(gca,'XLim',[-2 2],'YLim',[-2 2],'ZLim',[-2 2]);
axis square;

%%
% Seeing how tweaking the coefficients/eigenvalues change the shape of the
% surface

abc = [-1 1 -2];
[x,y,z] = meshgrid([-2:.2:2],[-2:.2:2],[-2:.2:2]);
variables = [x(:) y(:) z(:)];
fr = variables*abc';

% plotting the parallel planes
alpha = abc(1); beta = abc(2); gamma = abc(3);
coefficients = [alpha^2 beta^2 gamma^2 2*alpha*beta 2*alpha*gamma 2*beta*gamma];
variables = [x(:).^2 y(:).^2 z(:).^2 x(:).*y(:) x(:).*z(:) y(:).*z(:)];
fr = variables*coefficients';
figure;
p = patch(isosurface(x,y,z,reshape(fr,size(x)), 1))
set(p, 'FaceColor', 'blue', 'EdgeColor', 'none');
set(gca,'XLim',[-2 2],'YLim',[-2 2],'ZLim',[-2 2]);
axis square;
A = [alpha^2 alpha*beta alpha*gamma;...
    alpha*beta beta^2 beta*gamma;...
    alpha*gamma beta*gamma gamma^2];
[v,d] = eig(A)

% Tweaking the eigenvalues
d(1,1) = -d(3,3)/10;
d(2,2) = d(3,3)/10;
newcoeffs = v*d*v';
coefficients = [newcoeffs(1,1) newcoeffs(2,2) newcoeffs(3,3) 2*newcoeffs(1,2) 2*newcoeffs(1,3) 2*newcoeffs(2,3)];
fr = variables*coefficients';

p = patch(isosurface(x,y,z,reshape(fr,size(x)), 1))
set(p, 'FaceColor', 'blue', 'EdgeColor', 'none');
set(gca,'XLim',[-2 2],'YLim',[-2 2],'ZLim',[-2 2]);
axis square;

%%
% Trying to fit real data using this new way of representing a pair or
% parallel planes as a quadratic.
NTfilename = 'S080910003.nex';  % Beautiful funnel
%NTfilename = 'Sedna/S050410006.nex';  % planar
%NTfilename = 'Sedna/S033010003.nex';  % Cylinder?
%NTfilename = 'Kali/K070110005.nex';
%NTfilename = 'Kali/K091010004.nex'; % Beautiful funnel
% NTfilename = 'Kali/K102010006.nex'; % Pan color
%NTfilename = 'S0.nex';  


NT = nex2stro(findfile(NTfilename));
out = NTpreprocess(NT,.4, 0);
scaled = out(:,[2:4]).*repmat(out(:,5),[1 3]);
Loog = logical(out(:,end));
[planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);

rotateddata = scaled*xformmat;
figure; axes; hold on;
plot3(rotateddata(~Loog,1),rotateddata(~Loog,2),rotateddata(~Loog,3),'k.');
plot3(-rotateddata(~Loog,1),-rotateddata(~Loog,2),-rotateddata(~Loog,3),'k.');
plot3([zeros(sum(Loog),1) rotateddata(Loog,1)]',[zeros(sum(Loog),1) rotateddata(Loog,2)]',[zeros(sum(Loog),1) rotateddata(Loog,3)]','y-');
plot3([zeros(sum(Loog),1) -rotateddata(Loog,1)]',[zeros(sum(Loog),1) -rotateddata(Loog,2)]',[zeros(sum(Loog),1) -rotateddata(Loog,3)]','y-');

% Plotting the planar fits
tmp1 = linspace(-mean(abs(rotateddata(:,1))),mean(abs(rotateddata(:,1))),2);
tmp2 = linspace(-mean(abs(rotateddata(:,2))),mean(abs(rotateddata(:,2))),2);
[xx yy] = meshgrid(tmp1,tmp2);
planez = -(1+planeparams(1)*xx+planeparams(2)*yy)/planeparams(3); % plane 1
h = surf(xx,yy,planez);
set(h,'FaceAlpha',.2);
planez = (1-planeparams(1)*xx-planeparams(2)*yy)/planeparams(3); % plane 2
h = surf(xx,yy,planez);
set(h,'FaceAlpha',.2);

% Plotting quadratic fit
m = max(max(abs(rotateddata(~Loog,:))));
domain = 3*linspace(-m,m,30);
[x,y,z] = meshgrid(domain,domain,domain);
variables = [x(:).^2 y(:).^2 z(:).^2 2*x(:).*y(:) 2*x(:).*z(:) 2*y(:).*z(:)];
fr = variables*quadparams;
p = patch(isosurface(x,y,z,reshape(fr,size(x)), 1));
set(p,'FaceAlpha',.2,'FaceColor','green','Edgealpha',0);
axis vis3d;
camlight;
lighting phong;
axis square;


% % can we get a better fit with a different initial guess?
% % Not obviously
% 
% A = [planeparams(1).^2 planeparams(1)*planeparams(2) planeparams(1)*planeparams(3);...
%      planeparams(1)*planeparams(2) planeparams(2).^2 planeparams(2)*planeparams(3);...
%      planeparams(1)*planeparams(3) planeparams(2)*planeparams(3) planeparams(3).^2];
%  
% % Doing a grid search in a 3-D space spanned by the eigenvectors
% [v,d] = eig(A);
% [th,ph,r] = cart2sph(rotateddata(:,1),rotateddata(:,2),rotateddata(:,3));
% out = [];
% eigvals = [-.5:.025:1];
% for i = 1:length(eigvals)
%     for j = 1:length(eigvals)
%         for k = 1:length(eigvals)
%             dd = d;
%             dd(1,1) = eigvals(i);
%             dd(2,2) = eigvals(j);
%             dd(3,3) = d(3,3)+eigvals(k);
%             newcoeffs = v*dd*v';
%             
%             predr2 = 1./(newcoeffs(1,1).*(cos(ph).*cos(th)).^2 +...
%                 newcoeffs(2,2).*(cos(ph).*sin(th)).^2 +...
%                 newcoeffs(3,3).*sin(ph).^2 + ...
%                 2*newcoeffs(1,2).*cos(ph).*cos(th).*cos(ph).*sin(th) +...
%                 2*newcoeffs(1,3).*cos(ph).*cos(th).*sin(ph) +...
%                 2*newcoeffs(2,3).*cos(ph).*sin(th).*sin(ph));
%             if (any(predr2<0))
%                 predr2(predr2<0) = inf;
%             end
%             predr = sqrt(predr2);
%             
%             resid = log(abs(predr))-log(r);
%             resid(Loog&(abs(predr) > r)) = 0;
%             out(i,j,k) = sum(resid.^2);
%         end
%     end
% end
% 
% %figure;
% %imagesc(out);
% min(out(:))
% [i,j,k] = ind2sub([length(eigvals), length(eigvals), length(eigvals)],find(out == min(out(:))))
% dd = d;
% dd(1,1) = eigvals(i);
% dd(2,2) = eigvals(j);
% dd(3,3) = d(3,3)+eigvals(k);
% newcoeffs = v*dd*v';
% initguess = [newcoeffs(1,1) newcoeffs(2,2) newcoeffs(3,3) newcoeffs(1,2) newcoeffs(1,3) newcoeffs(2,3)];
% m = max(rotateddata(:));
% domain = linspace(-m,m,20);
% [x,y,z] = meshgrid(domain,domain,domain);
% variables = [x(:).^2 y(:).^2 z(:).^2 2*x(:).*y(:) 2*x(:).*z(:) 2*y(:).*z(:)];
% fr = variables*initguess';
% p = patch(isosurface(x,y,z,reshape(fr,size(x)), 1));
% set(p,'FaceAlpha',.2,'FaceColor','red','Edgealpha',0);
% axis vis3d;
% camlight;
% lighting phong;
% axis square;
% 
% [quadparams, quadSSEtobeat, exitflag, out4] = fminsearch(@(x) surfacefiterr4(rotateddata, x, Loog),initguess,options);
% quadSSEtobeat
% fr = variables*quadparams';
% p = patch(isosurface(x,y,z,reshape(fr,size(x)), 1));
% set(p,'FaceAlpha',.2,'FaceColor','blue','Edgealpha',0);
% 
%%
% What do the eigenvalues/vectors of 1-sheet hyperboloid mean?
% I assume the negative eigenvector is pointing down the hollow part of the
% surface, but it would be good to know for sure.  Confirmed.

NT = nex2stro(findfile('K082409005'));
out = NTpreprocess(NT,.4, 2);
scaled = out(:,[2:4]).*repmat(out(:,5),[1 3]);
Loog = logical(out(:,end));
[planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);

% Plotting quadratic fit
figure; axes; hold on;
rotateddata = scaled*xformmat;
m = max(max(abs(rotateddata(~Loog,:))));
domain = 3*linspace(-m,m,30);
[x,y,z] = meshgrid(domain,domain,domain);
variables = [x(:).^2 y(:).^2 z(:).^2 2*x(:).*y(:) 2*x(:).*z(:) 2*y(:).*z(:)];
fr = variables*quadparams;
p = patch(isosurface(x,y,z,reshape(fr,size(x)), 1));
set(p,'FaceAlpha',.2,'FaceColor','green','Edgealpha',0);
axis vis3d;
camlight;
lighting phong;
axis square;

A = [quadparams(1) quadparams(4) quadparams(5);...
    quadparams(4) quadparams(2) quadparams(6);...
    quadparams(5) quadparams(6) quadparams(3)];
[v,d] = eig(A);

for i = 1:3
    h = plot3([0 v(1,i)],[0 v(2,i)],[0 v(3,i)],'y-','Linewidth',2);
    if d(i,i)<0
        set(h,'Color','red')
    end
end
