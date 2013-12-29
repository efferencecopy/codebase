%%
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
    nspikes(i,:) = sum(spiketimes > stimon_t(i)+NT.sum.exptParams.latency/1000 & spiketimes < stimoff_t(i));
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
plot(cat(1,data.contrast),exp(beta(1)*cat(1,data.contrast)+beta(2)))
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
% Fitting a plane and plotting contrast-response functions normalized to
% distance from plane.  

load ('T_cones_smj10');

%NTfilename = 'S020310002' % Lum
%NTfilename = 'S021010011'; % lum
%NTfilename = 'S021810003'; % lum
%NTfilename = 'S021810006'; % lum
NTfilename = 'S021810008'; % lum
%NTfilename = 'S033010003'; % DO complicated

%NTfilename = 'S032710003'; % Lum
%NTfilename = 'S040210003'; % L-M
%NTfilename = 'S050410006'; % Lum
%NTfilename = 'S080610002' % Complicated funnel
%NTfilename = 'S080910003' % L-M+S funnel (insensitive)
%NTfilename = 'S081010005' % Complicated cylinder
%NTfilename = 'S090210005' % L-M
%NTfilename = 'S090910005' % L-M 
%NTfilename = 'S091310003' % L-M (funnel)


NT = nex2stro(findfile(NTfilename));
sf = NT.trial(:,strcmp(NT.sum.trialFields(1,:),'sf')); % sf wasn't dropped correctly in early files!
if (sf(1) == 0)
    error('No sf');
end
out = NTpreprocess(NT,.4, 1);
NTscaled = out(:,[2:4]).*repmat(out(:,5), 1,3);  % Not the same as DTspot convention (100 x smaller)
fundamentals = NT.sum.exptParams.fundamentals;
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = NT.sum.exptParams.mon_spd;
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;
M10 = T_cones_smj10*mon_spd;
NTscaled = NTscaled*(M10/M)';
[planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(NTscaled, zeros(size(NTscaled,1),1));

% Doing some plotting (for sanity checks)
figure; axes; hold on;
plot3(NTscaled(:,1),NTscaled(:,2),NTscaled(:,3),'y.')
plot3(-NTscaled(:,1),-NTscaled(:,2),-NTscaled(:,3),'y.')
coneweights = xformmat*planeparams;
%plot3([0 coneweights(1)]/100,[0 coneweights(2)]/100,[0 coneweights(3)]/100)
lim = 2;
[x y] = meshgrid(linspace(-lim,lim,10),linspace(-lim,lim,10));
z = -(1+planeparams(1)*x+planeparams(2)*y)/planeparams(3);
tmp = permute(cat(3,x,y,z),[3 1 2]);
tmp = reshape(tmp,3,size(x,2)*size(y,2));
xformed = xformmat'\tmp;
xformed = reshape(xformed, 3, size(x,2), size(y,2));
xformed = permute(xformed,[2 3 1]);
color = [.5 .5 .5];
h1 = surf(xformed(:,:,1),xformed(:,:,2),xformed(:,:,3));
set(h1,'CDataMapping','direct','FaceColor',color,'EdgeColor',color,'EdgeAlpha',.2,'FaceAlpha',.2);
h2 = surf(-xformed(:,:,1),-xformed(:,:,2),-xformed(:,:,3));
set(h2,'CDataMapping','direct','FaceColor',color,'EdgeColor',color,'EdgeAlpha',.2,'FaceAlpha',.2);

% Computing projection of each staircase termination onto conweights
threshproj = NTscaled*coneweights;
figure;
hist(abs(threshproj),logspace(-1,1,20))
xlabel('Projection onto coneweights');
set(gca,'XScale','log')

% First plotting normalized contrast-response functions.
spikeidx = strcmp(NT.sum.rasterCells(1,:),getSpikenum(NT));
coloridxs = NT.trial(:,strcmp(NT.sum.trialFields(1,:),'coloridx'));
stimon_t = NT.trial(:,strcmp(NT.sum.trialFields(1,:),'stim_on'));
stimoff_t = NT.trial(:,strcmp(NT.sum.trialFields(1,:),'stim_off'));
lms = NT.trial(:,[find(strcmp(NT.sum.trialFields(1,:),'lcont'))...
            find(strcmp(NT.sum.trialFields(1,:),'mcont'))...
            find(strcmp(NT.sum.trialFields(1,:),'scont'))]);
lms = lms*(M10/M)';

for i = 1:size(NT.trial,1)
    spiketimes = NT.ras{i,spikeidx};
    nspikes(i,:) = sum(spiketimes > stimon_t(i)+NT.sum.exptParams.latency/1000 & spiketimes < stimoff_t(i));
end

uniquecoloridxs = out(:,1);
uniquecolordirs = out(:,[2:4]);  % unit vectors
data = [];
for i = uniquecoloridxs'
    L = coloridxs == i;
    cdir = uniquecolordirs(uniquecoloridxs == i,:);
    contrast = abs(cdir*lms(L,:)');  % vector norms
    contrast = abs(lms(L,:)*xformmat*planeparams)';  % normalized so 1 = on plane
    data(length(data)+1).contrast = contrast;
    data(length(data)).lms = lms(L,:);
    data(length(data)).cdir = cdir;
    data(length(data)).response = nspikes(L);
    data(length(data)).Loog = out(uniquecoloridxs == i,end);
end

figure; axes; hold on;
for i = 1:length(data)
    if (~data(i).Loog)
        h = plot(data(i).contrast, data(i).response,'o');
    else
        h = plot(data(i).contrast, data(i).response,'+'); 
    end
    c = unifrnd(0,1,1,3);
    set(h,'MarkerFaceColor',c,'MarkerEdgeColor',c);
  %  i
  %  data(i).cdir
  %  pause
end
set(gca,'Xscale','log');

% Fitting a contrast response function
cr = [];
for i = 1:length(data)
    cr = [cr; data(i).contrast', data(i).response];
end
cr = sortrows(cr,1);
% Getting a coarse (nonparametric) view of the contrast response fn.
winlen = 20;
tmp = [];
for i = 1:winlen:size(cr,1)-winlen
    tmp = [tmp; median(cr([i:i+winlen-1],1)),median(cr([i:i+winlen-1],2))];
end
figure; axes; hold on;
plot(tmp(:,1),tmp(:,2),'b-')
set(gca,'Xscale','log');

% As a first pass, just fitting a linear contrast response function by LS
b = regress(cr(:,2),[log10(cr(:,1)) ones(size(cr,1),1)]);
plot(cr(:,1),b(1)*log10(cr(:,1))+b(2),'k:');
%b = glmfit(log10(cr(:,1)),cr(:,2),'poisson','link','identity');
%plot(cr(:,1),b(1)*log10(cr(:,1))+b(2),'k-');

% The moral of the story appears to be that contrast-response functions
% are not the same even when distance to the plane is taken into
% account.  Bleh.  Maybe the difference between L+M and L-M cells is big
% enough that this won't matter.
%
% Continuing from above
% Estimating response at psychophysical threshold
% Just asserting fitted iso-detection parameters for now
DTNTsfs = [0.5008 0.9919 1.9839 3.9677];
err = abs(DTNTsfs-sf(1));
whichsf = err == min(err);
params =[ 1.1817    2.8475    2.8556      3.0997;
        0.1203    0.0147   -0.0035      0.0012;
       -0.0052   -0.0033    0.0072      0.0126;
        0.1445    0.2288    0.1570      0.0491;
       -0.0367   -0.5428    0.2218     -0.2422;
       -0.3620   -0.2225    0.4693      0.2491;
       0.1134    0.0706    0.0446       0.0476;
       0.7796    1.2509   -0.7186       -0.3305;
       -0.6624   -1.2128    0.5492     -0.2676;
      -0.0508   -0.0515    0.0766      -0.0105];
params = params(:,whichsf);
mechs = reshape(params(2:end),3,3);

% Trying one color direction first
% cdir = mkbasis([1 1 0]);
% [th,phi,r] = cart2sph(cdir(1),cdir(2),cdir(3));
% options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-10);
% dist2PsychoSurface = fminsearch(@(r) colefitpredr(th,phi,r,params),r,options);
% sum(abs(dist2PsychoSurface*cdir*mechs).^params(1))  % sanity check should equal 1
% dist2PsychoSurface*cdir  % In 100*cone contrast units, this is on the isodetection surface
% dist2NeuroSurface = 1/(cdir*coneweights);
% dist2NeuroSurface*cdir % in cone contrast units, this is on the isoresponse plane
% dist2NeuroSurface*cdir*coneweights % Sanity check: should equal 1
% reldist = abs((dist2PsychoSurface/100)/dist2NeuroSurface)  % larger distance to psycho surface?!
% estcellresp = b(1)*log10(reldist)+b(2)

% Trying it in a few color directions now
options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-10);
[x,y,z] = sphere(60);
xyz = [x(:) y(:) z(:)];
estresp = zeros(size(x));
for i = 1:size(xyz,1)
    cdir = mkbasis(xyz(i,:));
    [th,phi,r] = cart2sph(cdir(1),cdir(2),cdir(3));
    dist2PsychoSurface = fminsearch(@(r) colefitpredr(th,phi,r,params),r,options);
    dist2NeuroSurface = abs(1/(cdir*coneweights));
    reldist = abs((dist2PsychoSurface/100)/dist2NeuroSurface);
    estresp(i) = max([eps b(1)*log10(reldist)+b(2)]);
    [x(i),y(i),z(i)] = sph2cart(th,phi,dist2PsychoSurface);
end
figure; axes; hold on;
h = surf(x,y,z,estresp./max(estresp(:)));
set(h,'EdgeColor','none')

% Plotting the planes
lim = .5;
[x y] = meshgrid(linspace(-lim,lim,10),linspace(-lim,lim,10));
z = -(1+planeparams(1)*x+planeparams(2)*y)/planeparams(3);
tmp = permute(cat(3,x,y,z),[3 1 2]);
tmp = reshape(tmp,3,size(x,2)*size(y,2));
xformed = xformmat'\tmp;
xformed = reshape(xformed, 3, size(x,2), size(y,2));
xformed = permute(xformed,[2 3 1])*100;
color = [.5 .5 .5];
h1 = surf(xformed(:,:,1),xformed(:,:,2),xformed(:,:,3),zeros(size(xformed(:,:,3))));
set(h1,'EdgeAlpha',.2,'FaceAlpha',.2);
h2 = surf(-xformed(:,:,1),-xformed(:,:,2),-xformed(:,:,3),zeros(size(xformed(:,:,1))));
set(h2,'FaceColor',color,'EdgeColor',color,'EdgeAlpha',.2,'FaceAlpha',.2);
colorbar('YTick',[0 .25 .5 .75 1],'YTickLabel',num2str(max(estresp(:))*[0 .25 .5 .75 1]',2))
title([num2str(DTNTsfs(whichsf),2),' cpd']);

% How far is the orthogonal vector that just touches the plane from the
% isodection surface?  This might provide a number we can use to decide
% whether the contrasts used in NeuroThresh data are sufficiently low that
% it makes sense to compare them to the DTNT data.
planevect =  coneweights'/norm(coneweights).^2*100;  % x100 to adhere to DTNT contrast convention
planevect*coneweights  % Sanity check: should equal 100
[th,phi,r1] = cart2sph(planevect(1),planevect(2),planevect(3));
options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-10);
r2 = abs(fminsearch(@(r) colefitpredr(th,phi,r,params),r1,options));
[x,y,z] = sph2cart(th,phi,r2);
% v = sum(abs([x y z] *reshape(params(2:end),3,3)).^params(1),2);
% "normfact" is the scale factor that you have to multiply planeparams by
% to get to the isodection surface
normfact = r2./r1
v = sum(abs(normfact*planevect*reshape(params(2:end),3,3)).^params(1),2)  % Sanity check should equal 1



%%
% Plotting the contrasts tested in NeuroThresh superimposed on the 3-D
% isodetection surface.  This should tell us if we're testing the same
% range of contasts in the two paradigms or how much they differ.
% What this shows is that we basically have *no* data inside of the
% isodetection contour.  All firing rate estimates at detection threshold 
% are necessarily extrapolations.

load ('T_cones_smj10');

NTfilename = 'K070809005'; % Lum
%NTfilename = 'K071709006'; % Lum
%NTfilename = 'K082109009'; % Lum 
%NTfilename = 'S032710003'; % Lum
%NTfilename = 'S050410006'; % Lum
%NTfilename = 'S040210003'; % L-M
NTfilename = 'S020310002' % Lum
%NTfilename = 'S090210005' % L-M
%NTfilename = 'S091310003' % L-M (funnel)
%NTfilename = 'S080610002' % Complicated funnel
%NTfilename = 'S081010005' % Complicated cylinder
%NTfilename = 'S090910005' % L-M 

NT = nex2stro(findfile(NTfilename));
sf = NT.trial(:,strcmp(NT.sum.trialFields(1,:),'sf')); % sf wasn't dropped correctly in early files!
if (sf(1) == 0)
    error('No sf');
end
out = NTpreprocess(NT,.4, 1);
NTscaled = out(:,[2:4]).*repmat(out(:,5), 1,3);  % Not the same as DTspot convention (100 x smaller)
fundamentals = NT.sum.exptParams.fundamentals;
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = NT.sum.exptParams.mon_spd;
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;
M10 = T_cones_smj10*mon_spd;
NTscaled = NTscaled*(M10/M)'*100;  % x 100 to match DTspot contrast convention

% Getting individual trial contrasts
coloridxs = NT.trial(:,strcmp(NT.sum.trialFields(1,:),'coloridx'));
lms = NT.trial(coloridxs>2,[find(strcmp(NT.sum.trialFields(1,:),'lcont'))...
            find(strcmp(NT.sum.trialFields(1,:),'mcont'))...
            find(strcmp(NT.sum.trialFields(1,:),'scont'))]);
lms = lms*(M10/M)'*100;  % x 100 to match DTspot contrast convention

DTNTsfs = [0.5008 0.9919 1.9839 3.2238 3.9677];
err = abs(DTNTsfs-sf(1));
whichsf = err == min(err);
params =[ 1.1817    2.8475    2.8556         0    3.0997;
        0.1203    0.0147   -0.0035         0    0.0012;
       -0.0052   -0.0033    0.0072         0    0.0126;
        0.1445    0.2288    0.1570         0    0.0491;
       -0.0367   -0.5428    0.2218         0   -0.2422;
       -0.3620   -0.2225    0.4693         0    0.2491;
       0.1134    0.0706    0.0446         0    0.0476;
       0.7796    1.2509   -0.7186         0   -0.3305;
       -0.6624   -1.2128    0.5492         0   -0.2676;
      -0.0508   -0.0515    0.0766         0   -0.0105];
params = params(:,whichsf);
% Plotting the surface
figure; axes; hold on;
[x,y,z] = meshgrid(linspace(-3,3,50),linspace(-3,3,50),linspace(-9,9,50));
tmp = [x(:) y(:) z(:)];
v = sum(abs(tmp *reshape(params(2:end),3,3)).^params(1),2);
fv = isosurface(x,y,z,reshape(v,size(x,1),size(x,2),size(x,3)),1);
h = patch(fv);
set(h, 'FaceColor','green','EdgeColor','green');
set(h,'FaceAlpha',.5);
plot3(lms(:,1),lms(:,2),lms(:,3),'k.');
plot3(-lms(:,1),-lms(:,2),-lms(:,3),'k.');
plot3(NTscaled(:,1),NTscaled(:,2),NTscaled(:,3),'m*');
plot3(-NTscaled(:,1),-NTscaled(:,2),-NTscaled(:,3),'m*');
set(gca,'Xlim',[-3 3],'Ylim',[-3 3],'Zlim',[-9 9]);
title(['sf = ',num2str(DTNTsfs(whichsf),2)]);

% How many points are inside the volume?
v = sum(abs(lms *reshape(params(2:end),3,3)).^params(1),2);
disp('# points inside surface  Total # pts.');
disp([sum(v < 1) length(v)])
disp('# points inside 5x surface  Total # pts.');
disp([sum(v < 5) length(v)])

%%
%  Comparing contrast sensitivity in NeuroThresh and DTspot for those few
%  cells for which we have data in both paradigms.

%DTfilename = 'S030310011'; NTfilename = 'S030310010';  % pan color
DTfilename = 'K091010003'; NTfilename = 'K091010004';  % tightly tuned
%DTfilename = 'K102010002'; NTfilename = 'K102010006';  % pan color
%DTfilename = 'K020311002'; NTfilename = 'K020311003'; % L-M plane
%DTfilename = 'K020311005'; NTfilename = 'K020311006'; % pan color very weak resp
%DTfilename = 'S022210003'; NTfilename = 'S022210002'; % Little DTspot data - high baseline cell

% Loading NT file ----------------------
NT = nex2stro(findfile(NTfilename));
NTout = NTpreprocess(NT,.4, 1);
Loog = logical(NTout(:,7));
NTscaled = NTout(:,[2:4]).*repmat(NTout(:,5), 1,3);
fundamentals = NT.sum.exptParams.fundamentals;
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = NT.sum.exptParams.mon_spd;
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
mon_spd = SplineSpd([380:4:780]', mon_spd, [380:5:780]');
M_NT = fundamentals'*mon_spd;
[planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(NTscaled, Loog,'ellipsoid');
NTconeweights = xformmat*planeparams;
% Done loading NT file ----------------------

% Now loading DTspot file ----------------------
DT = nex2stro(findfile(DTfilename));
flashside = sign(DT.trial(:,strcmp(DT.sum.trialFields(1,:),'flash_x')))*sign(DT.sum.exptParams.rf_x);  % + means in RF
gtrials = logical(DT.trial(:,strcmp(DT.sum.trialFields(1,:),'trial_type')));
% Getting rid of T2 and grating trials
L = flashside == -1 | gtrials;
DT.trial(L,:) = [];
DT.ras(L,:) = [];
stimon_t = DT.trial(:,strcmp(DT.sum.trialFields(1,:),'flash_on'));
stimoff_t = DT.trial(:,strcmp(DT.sum.trialFields(1,:),'flash_off'));
bkgndrgb = [DT.sum.exptParams.bkgnd_r, DT.sum.exptParams.bkgnd_g, DT.sum.exptParams.bkgnd_b];
flashR = DT.trial(:,strcmp(DT.sum.trialFields(1,:),'flash_R'));
flashG = DT.trial(:,strcmp(DT.sum.trialFields(1,:),'flash_G'));
flashB = DT.trial(:,strcmp(DT.sum.trialFields(1,:),'flash_B'));
coloridx = DT.trial(:,strcmp(DT.sum.trialFields(1,:),'color_dir'));
M_DT = reshape(DT.sum.exptParams.m_mtx,3,3);
monspect = reshape(DT.sum.exptParams.mon_spect,81,3);
% monitor spectra differ by a factor of 5/4 if comparing SplineSpd to
% SplineRaw.  This doesn't matter if we're looking at cone contrasts.
DTcoldirs = mkbasis(reshape(DT.sum.exptParams.RF_colors,3,3));
DTcoldirs(:,3) = [];
bkgndlms = M_DT * bkgndrgb';
x = 0:255; %the normal range of the gamma look up table
xx = linspace(0, 255, 2^16); %the desired quantization of the gammaTable
g1 = reshape(DT.sum.exptParams.gamma_table, 256, 3);
gammaTable = [spline(x, g1(:,1), xx)', spline(x, g1(:,2), xx)', spline(x, g1(:,3), xx)'];
RGB = [flashR flashG flashB]+1;
rgb = [gammaTable(RGB(:,1), 1), gammaTable(RGB(:,2), 2), gammaTable(RGB(:,3), 3)];
LMS = (M_DT*rgb')-repmat(bkgndlms,1,size(rgb,1));
DTccs = LMS' ./ repmat(bkgndlms',size(rgb,1),1);
spikeidx = 1; % lame
spikerates = nan(size(DT.trial,1),1);
for i = 1:size(DT.trial,1)
    spiketimes = DT.ras{i,spikeidx};
    nspikes = sum(spiketimes > stimon_t(i)+NT.sum.exptParams.latency/1000 & spiketimes < stimoff_t(i));
    spikerates(i) = nspikes./(stimoff_t(i)-stimon_t(i)-NT.sum.exptParams.latency/1000);
end
% Done loading and preprocessing DTspot file ----------------------

% Plotting data from NeuroThresh (points on isoresponse surface)
figure; axes; hold on;
plot3(NTscaled(~Loog,1),NTscaled(~Loog,2),NTscaled(~Loog,3),'k.');
plot3(-NTscaled(~Loog,1),-NTscaled(~Loog,2),-NTscaled(~Loog,3),'k.');
plot3([zeros(sum(Loog),1) NTscaled(Loog,1)]'/10,[zeros(sum(Loog),1) NTscaled(Loog,2)]'/10,[zeros(sum(Loog),1) NTscaled(Loog,3)]'/10,'y-');
plot3([zeros(sum(Loog),1) -NTscaled(Loog,1)]'/10,[zeros(sum(Loog),1) -NTscaled(Loog,2)]'/10,[zeros(sum(Loog),1) -NTscaled(Loog,3)]'/10,'y-');
plot3(DTccs(coloridx == 1,1),DTccs(coloridx == 1,2),DTccs(coloridx == 1,3),'m.')
plot3(DTccs(coloridx == 2,1),DTccs(coloridx == 2,2),DTccs(coloridx == 2,3),'g.')

% Plotting data in whitened space
%[v,d] = eig(cov([NTscaled(~Loog,:); -NTscaled(~Loog,:)]));
%tmp  = NTscaled*v*1/sqrt(d);
%figure; axes; hold on;
%plot3(tmp(~Loog,1),tmp(~Loog,2),tmp(~Loog,3),'k.');
%plot3(-tmp(~Loog,1),-tmp(~Loog,2),-tmp(~Loog,3),'k.');
%plot3([zeros(sum(Loog),1) tmp(Loog,1)]'/10,[zeros(sum(Loog),1) tmp(Loog,2)]'/10,[zeros(sum(Loog),1) tmp(Loog,3)]'/10,'y-');
%plot3([zeros(sum(Loog),1) -tmp(Loog,1)]'/10,[zeros(sum(Loog),1) -tmp(Loog,2)]'/10,[zeros(sum(Loog),1) -tmp(Loog,3)]'/10,'y-');
%tmp  = DTccs*v*1/sqrt(d);
%plot3(tmp(coloridx == 1,1),tmp(coloridx == 1,2),tmp(coloridx == 1,3),'m.')
%plot3(tmp(coloridx == 2,1),tmp(coloridx == 2,2),tmp(coloridx == 2,3),'g.')

% Plotting contrast-response functions from DT data 
% and showing a predicted point from NeuroThresh (using a plane fit)
figure; axes; hold on; plotcols = ['m','g'];
vectornormcc = nan(size(coloridx));
for i = 1:2
    L_DT = coloridx == i;
    vectornormcc(L_DT) = DTccs(L_DT,:)*DTcoldirs(:,i);
    plot(vectornormcc(L_DT),spikerates(L_DT),[plotcols(i),'.']);
    % Finding the cone contrast vector norm on the plane
    cc = 1./abs(DTcoldirs(:,i)'*NTconeweights);
    plot(cc,NT.sum.exptParams.threshold,[plotcols(i),'*'])
end

% Plotting contrast-response functions from DT data
% and showing a predicted point from NeuroThresh (using a ellipsoid fit)
for i = 1:2
    xfDTcoldir = DTcoldirs(:,i)'*xformmat;
    tmp = xfDTcoldir'*xfDTcoldir;
    v = [tmp(1) tmp(5) tmp(9) 2*tmp(2) 2*tmp(3) 2*tmp(6)];
    cc = 1./abs(v*quadparams);
    plot(sqrt(cc),NT.sum.exptParams.threshold,[plotcols(i),'s'])
    % Testing  (should equal 1)
    tmp = xfDTcoldir'*cc*xfDTcoldir;
    v = [tmp(1) tmp(5) tmp(9) 2*tmp(2) 2*tmp(3) 2*tmp(6)];
    abs(v*quadparams)
end