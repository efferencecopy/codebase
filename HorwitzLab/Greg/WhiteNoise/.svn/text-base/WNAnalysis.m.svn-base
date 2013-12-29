% Scripts to analyze white noise data
%
% Loading the data

WN=nex2stro;
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

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stuff for looking at spike-triggered averages, etc.
spikename = getSpikenum(WN);
spikeidx = find(strcmp(WN.sum.rasterCells(1,:),spikename));
maxT = 9;
Lgunnoise = WN.trial(:,noisetypeidx) == 1;
Lconenoise = WN.trial(:,noisetypeidx) == 2;
if (sum(Lconenoise) == 0)
    whichnoisetype = 1;
    disp('Gun noise only');
elseif (sum(Lgunnoise) == 0)
    whichnoisetype = 2;
    disp('Cone noise only');
else
    querystr = ['Which noisetype? 1=gun (n=',num2str(sum(Lgunnoise)),') 2=cone (n=',num2str(sum(Lconenoise)),'): '];
    whichnoisetype = input(querystr);
end
tmpstro = WN;
if (whichnoisetype == 2) && (any(Lgunnoise))  % Removing gun noise
    disp(['Getting rid of ',num2str(sum(Lgunnoise)),' gun noise trials!']);
    tmpstro.ras(Lgunnoise,:) = [];
    tmpstro.trial(Lgunnoise,:) = [];
elseif (whichnoisetype == 1) && (any(Lconenoise))  % Removing cone noise
    disp(['Getting rid of ',num2str(sum(Lconenoise)),' cone noise trials!']);
    tmpstro.ras(Lconenoise,:) = [];
    tmpstro.trial(Lconenoise,:) = [];
end
out = getWhtnsStats(tmpstro,maxT,'STCOVmex', {nstixperside^2, 3, maxT}, spikename);
tmpstro = [];
STAs = out{1};
STCs = out{2};
nspikes = out{3};
%%
% Plotting the STA and PCs
if (~exist('mask'))
    mask = zeros(nstixperside, nstixperside);   % Pixel mask.  Someday make a tool to make this non-zero.
end

Lmask = logical(repmat(~mask(:),[3 1]));
PCs = [];
for i = 1:size(STCs,2)
    STC = reshape(STCs(:,i), 3*nstixperside^2, 3*nstixperside^2);
    subSTC = STC(Lmask, Lmask);
    subSTA = STAs(Lmask,i);
    P = eye(size(subSTC))-subSTA*inv(subSTA'*subSTA)*subSTA';
    subSTC = P*subSTC*P';
    [tmp,d] = eig(subSTC);
    v = repmat(double(Lmask),[1 size(tmp,2)]);
    v(Lmask,:) = tmp;
    [newd, idxs] = sort(diag(d));
    v = v(:,idxs);
    v = v(:,[end end-1 end-2 2]);  % Not taking last (trivially zero) vector
    PCs = cat(3, PCs, v);
end
% Normalizing images
nstixperside = WN.sum.exptParams.nstixperside;
template = reshape([1:nstixperside^2],nstixperside,nstixperside);
edgepixels = [template(:,1); template(1,[2:end-1])'; template(end,[2:end-1])'; template(:,end)];
if (all(mask(edgepixels)))
    PCs = PCs.*max(STAs(:));
else
    edgepixelidxs = [edgepixels; edgepixels+nstixperside^2; edgepixels+2*(nstixperside^2)];
    PCelements = PCs(edgepixelidxs,:,:);
    PCsds = std(PCelements);    % One std calculated per PC
    PCs = PCs.*repmat(std(STAs(:,1))./PCsds,[size(PCs,1),1,1]);
end

rowidxs = reshape([1:3*nstixperside^2],[nstixperside^2 3]);
maxes = []; mins = [];
imagevectors = [STAs, reshape(PCs,[300 size(PCs,2)*size(PCs,3)])];
for i = 1:3
    maxes = [maxes; max(max(imagevectors(rowidxs(:,i),:)))];
    mins = [mins; min(min(imagevectors(rowidxs(:,i),:)))];
end
potentialnormfactors = [(1-[.5; .5; .5]-eps)./(maxes-[.5; .5; .5]); (-[.5; .5; .5]+eps)./(mins-[.5; .5; .5])];
% 'eps' in above line is a kludge that is required for avoiding
% out of bounds errors.
potentialnormfactors(potentialnormfactors < 0) = []; % if min > mu or max < mu
normfactor = min(potentialnormfactors);

%muvect = reshape(repmat(bkgndrgb',nstixperside^2,1),nstixperside^2*3,1);
muvect = reshape(repmat([.5 .5 .5],nstixperside^2,1),nstixperside^2*3,1);

% Plotting
figure;
for i = 1:size(STAs,2)
    STA = normfactor*(STAs(:,i)-muvect)+muvect;
    STA = reshape(STA,[nstixperside nstixperside 3]);
    subplot(6,size(STAs,2),i);
    image(STA);
    set(gca,'XTick',[],'YTick',[]); axis square;
    if (i == 1)
        ylabel('STA');
    end
    for j = 1:size(v,2)
        PC = normfactor*(PCs(:,j,i)-muvect)+muvect;
        PC = reshape(PC,[nstixperside nstixperside 3]);
        subplot(6,size(STAs,2), j*size(STAs,2)+i);
        image(PC);
        set(gca,'XTick',[],'YTick',[]); axis square;
    end
    STC = reshape(STCs(:,i),[sqrt(length(STCs(:,i))),sqrt(length(STCs(:,i)))]);
    STV = reshape(diag(STC),[nstixperside nstixperside 3]);
    STV = (STV-mean(STV(:)))./(2*range(STV(:)))+.5;
    subplot(6,size(STAs,2),5*size(STAs,2)+i);
    image(STV);
    set(gca,'XTick',[],'YTick',[]); axis square;
    if (i == 1)
        ylabel('STV');
    end
end

%%
% Significance testing on the STAs, STVs, broken down by color
L = WN.trial(:,noisetypeidx)==1;
mu1idx = find(strcmp(WN.sum.trialFields(1,:),'mu1'));
mu2idx = find(strcmp(WN.sum.trialFields(1,:),'mu2'));
mu3idx = find(strcmp(WN.sum.trialFields(1,:),'mu3'));
sigma1idx = find(strcmp(WN.sum.trialFields(1,:),'sigma1'));
sigma2idx = find(strcmp(WN.sum.trialFields(1,:),'sigma2'));
sigma3idx = find(strcmp(WN.sum.trialFields(1,:),'sigma3'));
muvect = unique(WN.trial(L,[mu1idx mu2idx mu3idx]),'rows')/1000;
sigmavect = unique(WN.trial(L,[sigma1idx sigma2idx sigma3idx]),'rows')/1000;
sigmavect(all(any(sigmavect == 0),2),:) = [];
gausslims = [WN.sum.exptParams.gauss_locut WN.sum.exptParams.gauss_hicut]/1000;
mumat = repmat(reshape(repmat(muvect,nstixperside^2,1),[nstixperside^2*3, 1]),[1,size(STAs,2)]);
sigmamat = repmat(reshape(repmat(sigmavect,nstixperside^2,1),[nstixperside^2* 3, 1]),[1,size(STAs,2)]);
zscoremeans = (STAs-mumat)./(sigmamat*sqrt(nspikes));

% Doing the calculations for the variances
% Only considering one correction factor per dimension
% (assuming variances on green and blue guns are same as red gun)

NPOINTS = 65536;
x = linspace(gausslims(1),gausslims(2),NPOINTS);
Fx = norminv(x)*sigmavect(1);
sigmacorrectionfactor = std(Fx)./sigmavect(1);
for i = 1:size(STCs,2)
    STC = reshape(STCs(:,i),[sqrt(length(STCs(:,i))),sqrt(length(STCs(:,i)))]);
    STVs(:,i) = diag(STC)./nspikes;
end
muvar = (sigmavect(1)*sigmacorrectionfactor)^2;
sigmavar = muvar*sqrt(2/nspikes);
zscorevars = (STVs-muvar)./sigmavar;
maxzscore = max(abs([zscoremeans(:);zscorevars(:)]));
alpha = 0.01;
crit = norminv(1-alpha,0,1);
% Plotting
figure;
for i = 1:size(STAs,2)
    % First STAs
    zmat = reshape(zscoremeans(:,i),[nstixperside nstixperside 3]);
    for j = 1:3
        pmat = logical(abs(zmat(:,:,j))>crit);
        im = zmat(:,:,j)./(2*maxzscore)+.5;
        im = repmat(im,[1 1 3]);
        sigidxs = find(pmat);
        im(sigidxs) = .5;  % red to .5 where sig.  Looks red on dark and cyan on bright.
        subplot(6,size(STAs,2),(j-1)*size(STAs,2)+i);
        image(im);
        axis image;
        set(gca,'XTick',[],'YTick',[]);
    end
    % Then STVs
    zmat = reshape(zscorevars(:,i),[nstixperside nstixperside 3]);
    for j = 4:6
        pmat = logical(abs(zmat(:,:,j-3))>crit);
        im = zmat(:,:,j-3)./(2*maxzscore)+.5;
        im = repmat(im,[1 1 3]);
        sigidxs = find(pmat);
        im(sigidxs) = .5;
        subplot(6,size(STAs,2),(j-1)*size(STAs,2)+i);
        image(im);
        axis image;
        set(gca,'XTick',[],'YTick',[]);
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Looking at the responses to synthetic images
% Finding which trials were synthetic image trials
imageidx = find(strcmp(WN.sum.rasterCells(1,:),'synth_image'));

allimages = WN.ras(:,imageidx);
L = [];
for i = 1:length(allimages)
    L = [L; ~isnan(allimages{i}(1))];
end
figure;
plot(L);
[x,y] = ginput(2);
close;
imageidxs = find(L);
imageidxs = imageidxs(imageidxs<max(x) & imageidxs>min(x));
% Pulling out the images and putting them through the gamma functions
% to get them in normalized intensity units.
nelements = length(allimages{imageidxs(1)});
npix = nelements/3;


npixperside = sqrt(npix);
images = nan(nelements,length(imageidxs));
for i = 1:length(imageidxs)
    im = allimages{imageidxs(i)};
    for j = 1:3
        idxs = [1:npix]+npix*(j-1);
        images(idxs,i) = gammaTable1(im(idxs),j);
    end
end

% Subtracting the background
bkgndrgbmat = repmat(bkgndrgb', npix, 1);
bkgndrgbmat = repmat(bkgndrgbmat(:), 1, length(imageidxs));
[uniqueimages, i, whichimage] = unique((images-bkgndrgbmat)','rows');
nuniqueims = max(whichimage);
spikerates = {};
meanfr = [];
stdfr = [];
ns = [];
for i = 1:max(whichimage)
    Lcont = whichimage == i;
    stimon_t = WN.trial(imageidxs(Lcont),stimonidx);
    numframes = WN.trial(imageidxs(Lcont),nframesidx)
    stimoff_t = stimon_t+numframes/framerate;
    spikes = WN.ras(imageidxs(Lcont),spikeidx);
    stimdurations = stimoff_t-stimon_t;
    for j = 1:sum(Lcont)
        spikerates{i}(j) = sum((spikes{j}>stimon_t(j)) & (spikes{j} < stimoff_t(j)))/stimdurations(j);
    end
    meanfr = [meanfr; mean(spikerates{i})];
    stdfr = [stdfr; std(spikerates{i})];
    ns = [ns; sum(Lcont)];
end
sem = stdfr./sqrt(ns);
sem(ns < 2) = nan;

[u,s,v] = svd(uniqueimages');
if ((size(s,2) == 1) || s(1,1)./s(2,2) > 1000)    % 1-D images
    [y,contrastsort] = sort(v(:,1));
    figure;    
    subplot(2,1,1); hold on;
    plot(v(contrastsort,1)*s(1),meanfr(contrastsort),'ko-', 'Linewidth', 1,'MarkerSize',4)
    plot([v(contrastsort,1) v(contrastsort,1)]*s(1),[meanfr(contrastsort)-sem(contrastsort) meanfr(contrastsort)+sem(contrastsort)],'k-.')
    ylabel('sp/sec');
    xlabel('Projection');
    subplot(2,2,3);
    avgim = v(contrastsort(1))*u(:,1)+bkgndrgbmat(:,1);
    image(reshape(avgim,npixperside, npixperside, 3));
    axis image; set(gca,'Xtick',[],'Ytick',[]);
    subplot(2,2,4);
    avgim = v(contrastsort(end))*u(:,1)+bkgndrgbmat(:,1);
    image(reshape(avgim,npixperside, npixperside, 3));
    axis image; set(gca,'Xtick',[],'Ytick',[]);
    
    figure;  % Subplots may not be arranged in a meaningful order
    for i = 1:nuniqueims
        subplot(1,nuniqueims,i);
        scalefactor = .5/max(abs(uniqueimages(:)));
        image(reshape(scalefactor*uniqueimages(i,:)',[nstixperside, nstixperside, 3])+.5);
        title(num2str(meanfr(i),2));
        set(gca,'Xtick',[],'YTick',[]);
        axis image;
    end
    
else  % 2-D images
    % Rotating coordinate axes
    dists = sum(v(:,1).^2+v(:,2).^2,2);
    corneridx = find(dists == max(dists),1);
    thetas = [atan(v(corneridx,2)./v(corneridx,1))+pi/4;...
                 atan(v(corneridx,2)./v(corneridx,1))-pi/4];
    theta = thetas(abs(thetas) == min(abs(thetas)));
    rotmat = [cos(theta) sin(theta); -sin(theta) cos(theta)];
    basis = u(:,[1 2])*rotmat' * rotmat*s([1 2], [1 2])*rotmat';
    basis = pinv(basis)';  % im = b'*w so w = im*pinv(b)
    norms = sqrt(sum(basis.^2));
    basis = basis./repmat(norms, size(basis,1),1);
    % This is flawed because the noise in the STA contributes to 
    % its norm but we're zeroing out nuisance elements in the PC1 so the 
    % when the are constrained to have the same norm, the "signal" in the
    % PC1 vector is huge, and so the weights onto the PC1 vector are
    % artificially small (relative to projections onto the STA).
    %weights = (rotmat*s([1 2], [1 2])*rotmat') * rotmat*v(:,[1 2])';

    weights = uniqueimages*basis;
    normweights = [(weights(:,1)-min(weights(:,1)))/(max(weights(:,1))-min(weights(:,1))),...
                   (weights(:,2)-min(weights(:,2)))/(max(weights(:,2))-min(weights(:,2)))];
    rankweights = round(normweights*(sqrt(nuniqueims)-1))+1;
    frim = zeros(sqrt(nuniqueims));
    semim = zeros(sqrt(nuniqueims));
    figure;
    for i = 1:nuniqueims;
        frim(rankweights(i,1),rankweights(i,2)) = meanfr(i);
        semim(rankweights(i,1),rankweights(i,2)) = sem(i);
        ind = sub2ind(max(rankweights), rankweights(i,1), rankweights(i,2))
        subplot(max(rankweights(:,1)), max(rankweights(:,2)), ind);
        scalefactor = .5/max(abs(uniqueimages(:)));
        image(reshape(scalefactor*uniqueimages(i,:)',[nstixperside, nstixperside, 3])+.5);
        title(num2str(meanfr(i),2));
        set(gca,'Xtick',[],'YTick',[]);
        axis image;
    end
    frim = frim';
    semim = semim';
    figure; imagesc(frim); axis xy; axis image;
    colormap(gray);
end

%%
% For L vs M replays
% L is on the x axis, M is on the y axis
LMScontrast = [];
for i = 1:nuniqueims
    im = reshape(uniqueimages(i,:),[nstixperside^2 3]);
    [u,s,v] = svd(im);
    if (max(abs(u(:,1))) ~= max(u(:,1)))
       s(1) = -s(1);
    end
    LMScontrast(:,i) = (M*im(find(abs(u(:,1))==max(abs(u(:,1)))),:)')./bkgndlms;
end
rmscontrast = unique(sqrt(mean(LMScontrast.^2)));  % in case we're interested in this
totalcontrast = sqrt(sum(LMScontrast.^2));
uniquecontrasts = unique(round(totalcontrast*100))/100;
data = [];
for i = 1:length(uniquecontrasts)
    L = abs(totalcontrast-uniquecontrasts(i)) < 1/100;
    if (sum(L) == 8)
       frs = meanfr(L);
       theta = atan2(LMScontrast(2,L), LMScontrast(1,L))';
       [theta,j] = sort(theta);
       data = cat(3,data,[theta, frs(j)]);
       disp('theta    contrast');
       disp([theta, totalcontrast(L)']); 
    end
end
figure; subplot(3,3,5);
posbColors = {'k', 'b', 'm', 'c', 'g', 'r'};
for i = size(data,3):-1:1
    theta = data(:,1,i);
    resp = data(:,2,i);
    polar([theta; theta(1)],[resp; resp(1)],[posbColors{i}, '-']);
    hold on;
end
ims = uniqueimages(L,:);
theta = squeeze(data(:,1,end));
resp = squeeze(data(:,2,end));
for i = 1:8
    binnedtheta = round(mod(theta(i),2*pi)./(pi/4));
    if (binnedtheta == 8 | binnedtheta == 0)
        plotnum = 6;
    elseif (binnedtheta > 4)
        plotnum = binnedtheta +2;
    elseif (binnedtheta == 4)
        plotnum = 4;
    else
        plotnum = 4-binnedtheta;
    end
    subplot(3,3,plotnum);
    image(reshape(ims(j(i),:)',[nstixperside, nstixperside, 3])+.5);
    title(num2str(resp(i),2));
    set(gca,'Xtick',[],'YTick',[]);
    axis image;
end
%%
% For L vs M replays at multiple contrasts
if (size(data,3) > 1)
    posbColors = {'k', 'b', 'm', 'c', 'g', 'r'};
    figure; axes; hold on; 
    for i = 1:size(data,3);
        plot(data(:,1,i),data(:,2,i),[posbColors{i}, 'o-'])
    end
end
set(gca,'YScale','log');
%%
% Rasters for the replay images

figure; axes; hold on;
counter = 0;
for i = 1:size(uniqueimages,2)
    Lcont = whichimage == i;
    stimon_t = WN.trial(imageidxs(Lcont),stimonidx);
    numframes = WN.trial(imageidxs(Lcont),nframesidx);
    stimoff_t = stimon_t+numframes/framerate;
    spikes = WN.ras(imageidxs(Lcont),spikeidx);
 
    for j = 1:sum(Lcont)
        nsp = length(spikes{j})
        if (nsp > 0)
            plot(0,counter,'g*');
            plot(stimoff_t(j)-stimon_t(j),counter,'r*');
        else
            plot(0,counter,'k*');
            plot(stimoff_t(j)-stimon_t(j),counter,'k*');
        end
        plot([spikes{j} spikes{j}]'-stimon_t(j),[zeros(nsp,1) .5*ones(nsp,1)]'+counter,'k-')
        counter = counter + 1;
    end
    plot([-.2 1],[counter counter],'k-');
end
set(gca,'Xlim',[-.2 1]); xlabel('time(s)');
title(WN.sum.fileName);


%%
% Getting the projections onto the PC1 and STA
% Just experimenting with the code

whichframe = 4;  % in the middle
nframestoconsider = 2;
msperframe = 1000/WN.sum.exptParams.framerate;
%%%%%%%%%%%%
% Skip this part if using a basis that came from somewhere else
%%%%%%%%%%
if (~exist('basis'))
    STC = reshape(STCs(:,whichframe), 3*nstixperside^2, 3*nstixperside^2);
    [v,d] = eig(STC);
    basis = [STAs(:,whichframe) v(:,end)];
    basis = basis./repmat(sqrt(sum(basis.^2)), size(basis,1),1);
end

nframestotal = sum(WN.trial(:,nframesidx));
initargs = {mkbasis(repmat(basis,nframestoconsider,1)) 3, nframestotal, [nstixperside^2 3 nframestoconsider]};
out = getWhtnsStats(WN,whichframe-ceil(nframestoconsider/2),'STPROJmod',initargs);
projs = out{1};
Lspike = out{2};

% 1-D firing rate plot
whichvect = 2;
nbins = 5;
x = linspace(prctile(projs(:,whichvect), 5), prctile(projs(:,whichvect), 95),nbins+2);
Na = hist(projs(:,whichvect), x);
Ns = zeros(size(Na));
for i = 1:max(Lspike)
    Ns = Ns + hist(projs(Lspike >=i,whichvect), x);
end
figure;
plot(x(2:end-1),Ns(2:end-1)./Na(2:end-1)/msperframe*1000,'ko');

% 2-D firing rate plot
if (round(sqrt(nuniqueims)) == sqrt(nuniqueims));
    nbins = sqrt(nuniqueims);
else
    nbins = 4;
end
x = [nbins nbins; prctile(projs, 10); prctile(projs, 90)];
%x = [nbins nbins; min(weights); max(weights)];

Na = hist2(projs, x);
Ns = zeros(size(Na));
for i = 1:max(Lspike)
    Ns = Ns + hist2(projs(Lspike >=i,:), x);
end

figure;
imageprops = mkbasis([-1 1; 0 1; 1 1; -1 0; 0 0; 1 0; -1 -1; 0 -1; 1 -1]')';
for i = [1 2 3 4 6 7 8 9] 
    subplot(3,3,i);
    tmpim = basis/2*imageprops(i,:)'+.5;  % just for a quick rendering
    image(reshape(tmpim,[nstixperside nstixperside 3]));
    axis image; set(gca,'XTick',[]); set(gca,'Ytick',[]);
end
subplot(3,3,5)
imagesc(Ns./Na);
axis image; set(gca,'XTick',[]); set(gca,'Ytick',[]);
colormap gray;
figure;
plot(frim,(Ns./Na)/msperframe*1000,'k.');

%%
% Rotating the basis vectors onto which we project the spike-triggered 
% stimuli.  Ideally, we'll see the best correlation between the projected
% firing rate plot and the empirical firing rate plot when we don't rotate
% at all.
nbins = sqrt(length(frim(:)));
x = [nbins nbins; min(weights); max(weights)];
data = [];
for theta = linspace(-pi/4,pi/4,5)
    theta
    rotmat = [cos(theta) sin(theta); -sin(theta) cos(theta)];
    initargs = {basis*rotmat', 3, nframestotal, [nstixperside^2 3 1]};
    out = getWhtnsStats(WN,whichframe,'STPROJmod',initargs);
    projs = out{1};
    Lspike = out{2};
    Na = hist2(projs, x);
    Ns = zeros(size(Na));
    for i = 1:max(Lspike)
        Ns = Ns + hist2(projs(Lspike >=i,:), x);
    end
    [rho,p] = corr(frim(:),Ns(:)./Na(:),'type','Kendall')
    data = [data; theta rho p];
end
figure;
plot(data(:,1)*180/pi,data(:,2),'k-.');
xlabel('Deg axes rotation')
ylabel('r')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Everything below this point is more special-purpose code...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% Looking for missed frames
% (Late peak is from the synth image trials)
stimon_t = WN.trial(:,stimonidx);
stimoff_t = WN.trial(:,stimoffidx);
framerate = WN.sum.exptParams.framerate;
predstimoff_t = stimon_t+WN.trial(:,nframesidx)/framerate;
x = stimoff_t-predstimoff_t;

hist(x(x<0),100);
plot(WN.trial(:,nframesidx), x,'k.')


%%
% Analyzing fixational saccades
sacstats = getSacData(WN);
L = logical(WN.trial(:,noisetypeidx)~= 0);
amplitudes = [];
directions = [];
peakv = [];
pathlengths = [];
durations = [];
BINWIDTH = .01;
timebins = [-.25:BINWIDTH:.5];
PSTHmat = zeros(sum(L),length(timebins));
for i = find(L')
    stimon_t = WN.trial(i,stimonidx);
    numframes = WN.trial(i,nframesidx);
    stimoff_t = stimon_t+numframes/WN.sum.exptParams.framerate;
    st = sacstats.starttimes{i};
    Lsac = (st > stimon_t) & (st < stimoff_t-.1); %.1 to avoid a saccade that leaves the window
    if any(Lsac)
        if any(sacstats.amplitudes{i}(Lsac) > 2)
            keyboard
        end
        amplitudes = [amplitudes; sacstats.amplitudes{i}(Lsac)];
        directions = [directions; sacstats.directions{i}(Lsac)];
        peakv = [peakv; sacstats.peakv{i}(Lsac)];
        pathlengths = [pathlengths; sacstats.pathlengths{i}(Lsac)];
        durations = [durations; sacstats.durations{i}(Lsac)];
        spiketimes = [];
        for j = find(Lsac')
            tmp = WN.ras{i,1}-sacstats.starttimes{i}(j);
            spiketimes = [spiketimes; tmp((tmp > timebins(1)-BINWIDTH/2) & (tmp < timebins(end)+BINWIDTH/2))];
        end
        [n,x] = hist(spiketimes, timebins);
        PSTHmat(i,:) = n./(BINWIDTH*sum(Lsac));
    end
end
figure;
subplot(3,2,1);
hist(amplitudes,100);
xlabel('amplitude'); ylabel('count');
subplot(3,2,3);
plot(amplitudes,peakv,'k.'); 
xlabel('amplitude'); ylabel('peak vel.');
subplot(3,2,4);
plot(amplitudes,pathlengths,'k.');
xlabel('amplitude'); ylabel('traj. length');
subplot(3,2,5);
[rho, theta]= hist(directions,20);
polar([theta theta(1)],[rho rho(1)],'k-')
subplot(3,2,6); hold on;
plot(timebins, mean(PSTHmat),'k-');
plot([0 0],[0 max(mean(PSTHmat))],'b:');
set(gca,'Xlim',[min(timebins) max(timebins)]);
xlabel('time wrt saccade onset (s)'); ylabel('sp/sec');
%%
% Cross validation to find static nonlinearities
% 1-D version

niters = 2;
CVprop = .7;
whichframe = 4;
nbins = 5;
Lgaussnoise = WN.trial(:,noisetypeidx) == 1;
binlim = 2*mode(WN.trial(Lgaussnoise,9))/1000
ntrials = size(WN.trial,1);
ntrialsSTC = round(ntrials*CVprop);
frs = [];
vs = [];
for iter = 1:niters
    trialidxs = randperm(ntrials);
    tmpstro = WN;
    tmpstro.trial = WN.trial(trialidxs(1:ntrialsSTC),:);
    tmpstro.ras = WN.ras(trialidxs(1:ntrialsSTC),:);
    out = getWhtnsStats(tmpstro, maxT, 'STCOVmex', {nstixperside^2, 3, whichframe});
    STAs = out{1};
    STCs = out{2};
    STC = reshape(STCs(:,whichframe), 3*nstixperside^2, 3*nstixperside^2);
    STC = STC.*~Lmask;
    [v,d] = eig(STC);
    d = diag(d);
    
    basis = [STAs(:,whichframe) v(:,end)];
    basis = basis./sqrt(sum(basis.^2));
    if (basis'*synthimvect < 0)% synthimvect is created by the chunk of code above
        basis = -basis;
    end
    tmpstro = WN;
    tmpstro.trial = WN.trial(trialidxs(ntrialsSTC+1:end),:);
    tmpstro.ras = WN.ras(trialidxs(ntrialsSTC+1:end),:);
    nframestotal = sum(WN.trial(:,nframesidx));
    initargs = {basis, whichframe-1, nframestotal, [nstixperside^2 3 1]};
    out = getWhtnsStats(tmpstro,whichframe,'STPROJmod',initargs);
    projs = out{1};
    Lspike = out{2};

    x = linspace(-binlim, binlim, nbins+2);
    Na = hist(projs(:,1), x);
    Ns = zeros(size(Na));
    for i = 1:max(Lspike)
        Ns = Ns + hist(projs(Lspike >=i,1), x);
    end
    frs = [frs; Ns./Na]
    vs = [vs,basis];
end

figure; subplot(2,1,1);
plot(x,mean(frs));
subplot(2,2,4);
image(reshape(mean(vs,2),[nstixperside, nstixperside, 3])+.5);  % Not rendered accurately
axis image; set(gca,'XTick',[],'YTick',[]);
subplot(2,2,3);
image(reshape(-mean(vs,2),[nstixperside, nstixperside, 3])+.5);  % Not rendered accurately
axis image; set(gca,'XTick',[],'YTick',[]);


%%
% Cross validation to find static nonlinearities
% 2-D version

niters = 10;
CVprop = .7;
whichframe = 4;
nbins = 10;
Lgaussnoise = WN.trial(:,noisetypeidx) == 1;
binlim = 1*mode(WN.trial(Lgaussnoise,9))/1000;
ntrials = size(WN.trial,1);
ntrialsSTC = round(ntrials*CVprop);
Nas = [];
Nss = [];
for iter = 1:niters
    trialidxs = randperm(ntrials);
    tmpstro = WN;
    tmpstro.trial = WN.trial(trialidxs(1:ntrialsSTC),:);
    tmpstro.ras = WN.ras(trialidxs(1:ntrialsSTC),:);
    out = getWhtnsStats(tmpstro, maxT, 'STCOVmex', {nstixperside^2, 3, whichframe});
    STAs = out{1};
    STCs = out{2};
    STC = reshape(STCs(:,whichframe), 3*nstixperside^2, 3*nstixperside^2);
    [v,d] = eig(STC);
    
    basis = mkbasis([STAs(:,whichframe) v(:,end)]);
    basis(:,2) = sign(basis(:,1)'*basis(:,2))*basis(:,2);   % PC constrained to have a + proj onto STA
    tmpstro = WN;
    tmpstro.trial = WN.trial(trialidxs(ntrialsSTC+1:end),:);
    tmpstro.ras = WN.ras(trialidxs(ntrialsSTC+1:end),:);
    nframestotal = sum(WN.trial(:,nframesidx));
    initargs = {basis, whichframe-1, nframestotal, [nstixperside^2 3 1]};
    out = getWhtnsStats(tmpstro,whichframe,'STPROJmod',initargs);
    projs = out{1};
    Lspike = out{2};

    x = [nbins nbins; -binlim -binlim; binlim binlim];
    Na = hist2(projs, x);
    Ns = zeros(size(Na));
    for i = 1:max(Lspike)
        Ns = Ns + hist2(projs(Lspike >=i,:), x);
    end
    Nas = cat(3,Nas,Na);
    Nss = cat(3,Nss,Ns);
end
fr = Nss./Nas;
figure;
imagesc(mean(fr,3));

%% Looking at STA from the multielectrode array (no STC analysis because stimulus is too big).
% 

WN=nex2stro;
LFPsamprate = WN.sum.analog.storeRates{1};  % Assuming all channels are sampled at the same rate
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
maxT = 9;

potentialspikeidxs = find(strncmp(WN.sum.rasterCells,'sig',3))
for spikeidx = potentialspikeidxs
    out = getWhtnsStats(WN,maxT,'STAmex', {nstixperside^2, 3, maxT}, WN.sum.rasterCells{spikeidx});
    nspikes = out{3};
    STAs = out{1};
    %STVs = (out{2}-out{1}.^2)./nspikes;
    STVs = out{2}./nspikes;
    
    % Significance testing on the STAs, STVs, broken down by color
    L = WN.trial(:,noisetypeidx)==1; % 1 = gun noise
    mu1idx = find(strcmp(WN.sum.trialFields(1,:),'mu1'));
    mu2idx = find(strcmp(WN.sum.trialFields(1,:),'mu2'));
    mu3idx = find(strcmp(WN.sum.trialFields(1,:),'mu3'));
    sigma1idx = find(strcmp(WN.sum.trialFields(1,:),'sigma1'));
    sigma2idx = find(strcmp(WN.sum.trialFields(1,:),'sigma2'));
    sigma3idx = find(strcmp(WN.sum.trialFields(1,:),'sigma3'));
    muvect = unique(WN.trial(L,[mu1idx mu2idx mu3idx]),'rows')/1000;
    sigmavect = unique(WN.trial(L,[sigma1idx sigma2idx sigma3idx]),'rows')/1000;
    sigmavect(all(any(sigmavect == 0),2),:) = [];
    gausslims = [WN.sum.exptParams.gauss_locut WN.sum.exptParams.gauss_hicut]/1000;
    if (sum(L) > sum(~L))  % If more gun noise than cone noise
        mumat = repmat(reshape(repmat(muvect,nstixperside^2,1),[nstixperside^2*3, 1]),[1,size(STAs,2)]);
        sigmamat = repmat(reshape(repmat(sigmavect,nstixperside^2,1),[nstixperside^2* 3, 1]),[1,size(STAs,2)]);
        zscoremeans = (STAs-mumat)./(sigmamat*sqrt(nspikes));
        
        % Doing the calculations for the variances
        % Only considering one correction factor per dimension
        % (assuming variances on green and blue guns are same as red gun)
        
        NPOINTS = 65536;
        x = linspace(gausslims(1),gausslims(2),NPOINTS);
        Fx = norminv(x)*sigmavect(1);
        sigmacorrectionfactor = std(Fx)./sigmavect(1);
        muvar = (sigmavect(1)*sigmacorrectionfactor)^2;
        sigmavar = muvar*sqrt(2/nspikes);
        zscorevars = (STVs-muvar)./sigmavar;
    else
        zscoremeans = (STAs-mumat)./sqrt(nspikes);
        zscorevars = zeros(size(zscoremeans));
    end
    maxzscore = max(abs([zscoremeans(:);zscorevars(:)]));
    alpha = 0.01;
    crit = norminv(1-alpha,0,1);
    
    % Plotting
    figure;
    for i = 1:size(STAs,2)
        % First STAs
        zmat = reshape(zscoremeans(:,i),[nstixperside nstixperside 3]);
        for j = 1:3
            pmat = logical(abs(zmat(:,:,j))>crit);
            im = zmat(:,:,j)./(2*maxzscore)+.5;
            im = repmat(im,[1 1 3]);
            sigidxs = find(pmat);
            im(sigidxs) = .5;  % red to .5 where sig.  Looks red on dark and cyan on bright.
            subplot(6,size(STAs,2),(j-1)*size(STAs,2)+i);
            image(im);
            axis image;
            set(gca,'XTick',[],'YTick',[]);
        end
        % Then STVs
        zmat = reshape(zscorevars(:,i),[nstixperside nstixperside 3]);
        for j = 4:6
            pmat = logical(abs(zmat(:,:,j-3))>crit);
            im = zmat(:,:,j-3)./(2*maxzscore)+.5;
            im = repmat(im,[1 1 3]);
            sigidxs = find(pmat);
            im(sigidxs) = .5;
            subplot(6,size(STAs,2),(j-1)*size(STAs,2)+i);
            image(im);
            axis image;
            set(gca,'XTick',[],'YTick',[]);
        end
    end
    set(gcf,'Name',WN.sum.rasterCells{spikeidx});
    drawnow;
end

%%
% Looking at LFP-triggered averages (filtered, interpolated LFPs)
centerfreq = 100; % Hz
cyclespersample = centerfreq./LFPsamprate;
ncyclesper6sigma = 3;  % Controls the bandwidth
nsamplesper6sigma = ceil(ncyclesper6sigma/cyclespersample);
ncycles = nsamplesper6sigma*cyclespersample;
filtkernel1 = normpdf(linspace(-3,3,nsamplesper6sigma),0,1).*cos(linspace(0,2*pi*ncycles,nsamplesper6sigma));
filtkernel2 = normpdf(linspace(-3,3,nsamplesper6sigma),0,1).*sin(linspace(0,2*pi*ncycles,nsamplesper6sigma));
filtkernel1 = filtkernel1./norm(filtkernel1);
filtkernel2 = filtkernel2./norm(filtkernel2);
for whichADchan = 1:32;
    ntrials = size(WN.trial,1);
    nstixperside = WN.sum.exptParams.nstixperside; %get number of stixels per side
    msperframe = 1000/WN.sum.exptParams.framerate; %calculate msec per frame
    secsperframe = 1/WN.sum.exptParams.framerate;
    gammaTable = WN.sum.exptParams.gamma_table; %get gamma_table
    gammaTable = reshape(gammaTable,length(gammaTable)/3,3); %reshapse gamma_table into three columns
    invgammaTable = InvertGamma(gammaTable,1); %invert gamma_table
    ngammasteps = size(invgammaTable,1); %get number of rows of gamma_table (65536)
    seedidx = strcmp(WN.sum.trialFields(1,:),'seed'); %get seed index from trialFields
    nframesidx = strcmp(WN.sum.trialFields(1,:),'num_frames'); %get nframes index from trialFields
    stimonidx = strcmp(WN.sum.trialFields(1,:),'stim_on'); %get stimon index from trialFields
    muidxs = [find(strcmp(WN.sum.trialFields(1,:),'mu1')),... %get mu indices into vector from trialFields
        find(strcmp(WN.sum.trialFields(1,:),'mu2')),...
        find(strcmp(WN.sum.trialFields(1,:),'mu3'))];
    sigmaidxs = [find(strcmp(WN.sum.trialFields(1,:),'sigma1')),... %get sigma indices into vector from trialFields
        find(strcmp(WN.sum.trialFields(1,:),'sigma2')),...
        find(strcmp(WN.sum.trialFields(1,:),'sigma3'))];
    LFPstarttimes = [WN.ras{:,strcmp(WN.sum.rasterCells,'anlgStartTime')}];
    
    feval('STAmex','init', {nstixperside^2, 3, maxT});
    % This code was largely lifted from getWhtnsStats
    hwait = waitbar(0,'Finding stimuli...'); %display progress bar
    L = WN.trial(:,noisetypeidx)==1; % 1 = gun noise
    for i = 1:ntrials %get values and insert into given column
        if (WN.trial(i,noisetypeidx)==1)
            if (sum(L) > sum(~L))
                noisetype = 1;
            else
                continue
            end
        end
        
        if (WN.trial(i,noisetypeidx)==2)
            if (sum(~L) > sum(L))
                noisetype = 2;
            else
                continue
            end
        end

        seed = WN.trial(i,seedidx);
        nframes = WN.trial(i,nframesidx);
        mu = WN.trial(i,muidxs)/1000;
        sigma = WN.trial(i,sigmaidxs)/1000;
        
        if (noisetype == 1)  % Gaussian gun
            x = linspace(WN.sum.exptParams.gauss_locut/1000, WN.sum.exptParams.gauss_hicut/1000, ngammasteps);%dividing gauss by 1000 and making equal intervals so that there are 65536 values
            for gun = 1:3
                invnormcdf(:,gun) = norminv(x)*sigma(gun)+mu(gun);
            end
            randnums = getEJrandnums(3*nstixperside^2*nframes, seed);
            randnums = reshape(randnums, [nstixperside^2*3, nframes]);
            for j = 1:3
                idxs = [1:nstixperside^2]+nstixperside^2*(j-1);
                randnums(idxs,:) = reshape(invnormcdf(randnums(idxs,:)+1,j),[length(idxs),nframes]);
            end
        elseif (noisetype == 2)  % Binary cone
            colordirlms = sign(fullfact([2,2,2])-1.5);
            randnums = getEJrandnums(nstixperside^2*nframes, seed);
            randnums = reshape(randnums, [nstixperside^2, nframes]);
            randnums = mod(randnums, size(colordirlms,1))+1;
            tmp = colordirlms(randnums,:);
            tmp = reshape(tmp, [nstixperside^2 nframes 3]);
            tmp = permute(tmp,[1 3 2]);
            randnums = reshape(tmp,[nstixperside^2*3 nframes]);
        end
        
        t_stimon = WN.trial(i,stimonidx);
        LFP = WN.ras{i,strcmp(WN.sum.rasterCells,['AD',num2str(16+whichADchan)])};
        filteredLFP = conv(LFP, filtkernel1,'same').^2+conv(LFP, filtkernel2,'same').^2;
        LFPtimes = LFPstarttimes(i)+[0:length(filteredLFP)-1]/LFPsamprate; % sec
        frametimes = t_stimon+(linspace(0, nframes*secsperframe, nframes)+(secsperframe/2)'); % sec
        rectLFP = interp1(LFPtimes, filteredLFP, frametimes); 
        %   plot(LFPtimes,abs(LFP)); hold on; plot(frametimes,rectLFP,'m-');
        feval('STAmex',randnums(:),rectLFP);
        if (ishandle(hwait))
            waitbar(i/ntrials,hwait);
        else
            break;
        end
    end
    close(hwait);
    out = feval('STAmex','return');
    eval('clear STAmex');
    
    % Now doing some plotting
    STAs = out{1};
    STVs = out{2};
    STAlims = [min(STAs(:)) max(STAs(:))];
    STVlims = [min(STVs(:)) max(STVs(:))];
    STAs = (STAs-STAlims(2))/(STAlims(1)-STAlims(2));
    STVs = (STVs-STVlims(2))/(STVlims(1)-STVlims(2));
    
    figure; colormap(gray(255)); 
    set(gcf,'Name',['LFP ',num2str(whichADchan)]);
    for i = 1:size(STAs,2)
        % First STAs
        mat = reshape(STAs(:,i),[nstixperside nstixperside 3]);
        for j = 1:3
            im = mat(:,:,j);
            subplot(6,size(STAs,2),(j-1)*size(STAs,2)+i);
            image(im*255);
            axis image;
            set(gca,'XTick',[],'YTick',[]);
        end
        % Then STVs
        mat = reshape(STVs(:,i),[nstixperside nstixperside 3]);
        for j = 1:3
            im = mat(:,:,j);
            subplot(6,size(STAs,2),(j-1+3)*size(STAs,2)+i);
            image(im*255);
            axis image;
            set(gca,'XTick',[],'YTick',[]);
        end
    end
end

%%
% Looking at rectified LFPs at the time of stimus onset 
ADchans = find(strncmp('AD',WN.sum.rasterCells,2));
LFPchans = ADchans(3:end);  % Getting rid of the eye position channels
hwait = waitbar(0,'Computing stimulus-locked LFPs...'); %display progress bar
poststim_t = 0.5; % secs
prestim_t = 0.3; % secs
nprestimsamples = prestim_t*LFPsamprate;
npoststimsamples = poststim_t*LFPsamprate;
nsamples = nprestimsamples+npoststimsamples+1;
PSTHs = zeros(length(LFPchans),nsamples);
for i = 1:ntrials %get values and insert into given column
    t_stimon = WN.trial(i, stimonidx);    
    for ADchan1 = LFPchans
        LFP = WN.ras{i,ADchan1};
        LFPtimes = LFPstarttimes(i)+[0:length(LFP)-1]/LFPsamprate; % sec
        err = (LFPtimes-t_stimon).^2;
        startidx = find(err == min(err),1);
        if (startidx > nprestimsamples & length(LFP)-startidx > npoststimsamples)
            LFPclip = LFP(startidx+[-nprestimsamples:npoststimsamples]);
            PSTHs(ADchan1 == LFPchans,:) = PSTHs(ADchan1 == LFPchans,:) + abs(LFPclip)';
        end
    end
    waitbar(i/ntrials,hwait);
end
close(hwait);
figure; subplot(2,1,1);
imagesc(PSTHs);
subplot(2,1,2);
tvect = linspace(-prestim_t, poststim_t, nsamples);
plot(tvect,mean(PSTHs));
% Cool.  I think we're seeing the stimulus refresh in the LFP.
L = tvect>.3;
freqs = linspace(0, LFPsamprate/2, sum(L)/2);
[Pxx,F] = periodogram(mean(PSTHs(:,L)),[],[],LFPsamprate) 
figure; plot(F,Pxx); set(gca,'YScale','log');

%%
% Looking at the auto- and cross-correlation functions of the LFPs

hwait = waitbar(0,'Finding stimuli...'); %display progress bar
maxlag = 0.2; % sec
delayfromstimon = 0.25;
maxlagsamples = maxlag*LFPsamprate;
nsamples = 2*maxlagsamples+1;
data = zeros(length(LFPchans),length(LFPchans), nsamples);
for i = 1:ntrials % get values and insert into given column
    t_stimon = WN.trial(i, stimonidx);
    t_stimoff = WN.trial(i, stimoffidx);
    for idx1 = 1:length(LFPchans)
        for idx2 = idx1:length(LFPchans)
            LFP1 = WN.ras{i,LFPchans(idx1)};
            LFP2 = WN.ras{i,LFPchans(idx2)};
            LFPtimes = LFPstarttimes(i)+[0:length(LFP1)-1]/LFPsamprate; % sec
            err = (LFPtimes-t_stimon-delayfromstimon).^2;
            startidx = find(err == min(err),1);
            err = (LFPtimes-t_stimoff).^2;
            endidx = find(err == min(err),1);
            LFP1 = LFP1(startidx:endidx);
            LFP2 = LFP2(startidx:endidx);
            if (length(LFP1) >= nsamples)
                tmp = xcorr(LFP1,LFP2,maxlagsamples,'coeff');
                data(idx1,idx2,:) = tmp;
            end
        end
    end
    waitbar(i/ntrials,hwait);
end
close(hwait);

% First looking at all the autocorrelation functions
tvect = [-maxlagsamples:maxlagsamples]/LFPsamprate;
tmp = zeros(length(LFPchans),nsamples);
for i = 1:length(LFPchans)
    tmp(i,:) = data(i,i,:);
end
figure;
subplot(2,1,1);
title('Autocorrelations');
imagesc(tmp); set(gca,'Xtick',[]);
subplot(2,1,2);
plot(tvect,tmp');
xlabel('lag secs');

% Now looking at cross-correlations
tmp = [];
for sig1 = 1:length(LFPchans)
    for sig2 = sig1+1:length(LFPchans)
        tmp = [tmp;squeeze(data(sig1,sig2,:))'];
    end
end
figure;
imagesc(tmp);

% % Looking at cross correlograms of adjacent and non-adjacent pairs of electrodes
% % Only appropriate for most medial bank of 32 (port C).
tmp1 = []; tmp2 = [];
for i = 1:length(LFPchans)-1
    if (mod(i,2))  % if i is odd
        tmp1 = [tmp1; squeeze(data(i,i+1,:))'];
    else  % if i is even
        tmp2 = [tmp2; squeeze(data(i,i+1,:))'];
    end
end
figure;
subplot(1,2,1); % Nearby pairs
plot(tvect, tmp1);
set(gca,'Ylim',[min([tmp1(:);tmp2(:)])*.9 max([tmp1(:);tmp2(:)])*1.1]); 
subplot(1,2,2); % Distant pairs
plot(tvect, tmp2);
set(gca,'Ylim',[min([tmp1(:);tmp2(:)])*.9 max([tmp1(:);tmp2(:)])*1.1]); 

%%
% Cranking through a bunch of frequencies to find where we get the best
% signal in the LFP-TA and LFP-TV
LFPsamprate = WN.sum.analog.storeRates{1};  % Assuming all channels are sampled at the same rate
whichADchan = 4;  % Pick something good
nfreqs = 10;
lowfreq = 30;
highfreq = 300;
%highfreq = LFPsamprate/2;  % Nyquist
ncyclesper6sigma = 3;  % Controls the bandwidth
ntrials = size(WN.trial,1);
nstixperside = WN.sum.exptParams.nstixperside; %get number of stixels per side
msperframe = 1000/WN.sum.exptParams.framerate; %calculate msec per frame
secsperframe = 1/WN.sum.exptParams.framerate;
gammaTable = WN.sum.exptParams.gamma_table; %get gamma_table
gammaTable = reshape(gammaTable,length(gammaTable)/3,3); %reshapse gamma_table into three columns
invgammaTable = InvertGamma(gammaTable,1); %invert gamma_table
ngammasteps = size(invgammaTable,1); %get number of rows of gamma_table (65536)
seedidx = strcmp(WN.sum.trialFields(1,:),'seed'); %get seed index from trialFields
nframesidx = strcmp(WN.sum.trialFields(1,:),'num_frames'); %get nframes index from trialFields
stimonidx = strcmp(WN.sum.trialFields(1,:),'stim_on'); %get stimon index from trialFields
muidxs = [find(strcmp(WN.sum.trialFields(1,:),'mu1')),... %get mu indices into vector from trialFields
    find(strcmp(WN.sum.trialFields(1,:),'mu2')),...
    find(strcmp(WN.sum.trialFields(1,:),'mu3'))];
sigmaidxs = [find(strcmp(WN.sum.trialFields(1,:),'sigma1')),... %get sigma indices into vector from trialFields
    find(strcmp(WN.sum.trialFields(1,:),'sigma2')),...
    find(strcmp(WN.sum.trialFields(1,:),'sigma3'))];
LFPstarttimes = [WN.ras{:,strcmp(WN.sum.rasterCells,'anlgStartTime')}];

for centerfreq = logspace(log10(lowfreq),log10(highfreq),nfreqs); % Hz
    cyclespersample = centerfreq./LFPsamprate;
    nsamplesper6sigma = ceil(ncyclesper6sigma/cyclespersample);
    ncycles = nsamplesper6sigma*cyclespersample;
    filtkernel1 = normpdf(linspace(-3,3,nsamplesper6sigma),0,1).*cos(linspace(0,2*pi*ncycles,nsamplesper6sigma));
    filtkernel2 = normpdf(linspace(-3,3,nsamplesper6sigma),0,1).*sin(linspace(0,2*pi*ncycles,nsamplesper6sigma));
    filtkernel1 = filtkernel1./norm(filtkernel1);
    filtkernel2 = filtkernel2./norm(filtkernel2);
    
    feval('STAmex','init', {nstixperside^2, 3, maxT});
    % This code was largely lifted from getWhtnsStats
    hwait = waitbar(0,'Finding stimuli...'); %display progress bar
    L = WN.trial(:,noisetypeidx)==1; % 1 = gun noise
    for i = 1:ntrials %get values and insert into given column
         if (WN.trial(i,noisetypeidx)==1)
            if (sum(L) > sum(~L))
                noisetype = 1;
            else
                continue
            end
        end
        
        if (WN.trial(i,noisetypeidx)==2)
            if (sum(~L) > sum(L))
                noisetype = 2;
            else
                continue
            end
        end
        
        seed = WN.trial(i,seedidx);
        nframes = WN.trial(i,nframesidx);
        mu = WN.trial(i,muidxs)/1000;
        sigma = WN.trial(i,sigmaidxs)/1000;
        
        if (noisetype == 1)  % Gaussian gun
            x = linspace(WN.sum.exptParams.gauss_locut/1000, WN.sum.exptParams.gauss_hicut/1000, ngammasteps);%dividing gauss by 1000 and making equal intervals so that there are 65536 values
            for gun = 1:3
                invnormcdf(:,gun) = norminv(x)*sigma(gun)+mu(gun);
            end
            randnums = getEJrandnums(3*nstixperside^2*nframes, seed);
            randnums = reshape(randnums, [nstixperside^2*3, nframes]);
            for j = 1:3
                idxs = [1:nstixperside^2]+nstixperside^2*(j-1);
                randnums(idxs,:) = reshape(invnormcdf(randnums(idxs,:)+1,j),[length(idxs),nframes]);
            end
        elseif (noisetype == 2)  % Binary cone
            colordirlms = sign(fullfact([2,2,2])-1.5);
            randnums = getEJrandnums(nstixperside^2*nframes, seed);
            randnums = reshape(randnums, [nstixperside^2, nframes]);
            randnums = mod(randnums, size(colordirlms,1))+1;
            tmp = colordirlms(randnums,:);
            tmp = reshape(tmp, [nstixperside^2 nframes 3]);
            tmp = permute(tmp,[1 3 2]);
            randnums = reshape(tmp,[nstixperside^2*3 nframes]);
        end
        
        t_stimon = WN.trial(i, stimonidx);
        
        LFP = WN.ras{i,strcmp(WN.sum.rasterCells,['AD',num2str(16+whichADchan)])};
        filteredLFP = conv(LFP, filtkernel1,'same').^2+conv(LFP, filtkernel2,'same').^2;
        LFPtimes = LFPstarttimes(i)+[0:length(filteredLFP)-1]/LFPsamprate; % sec
        frametimes = t_stimon+(linspace(0, nframes*secsperframe, nframes)+(secsperframe/2)'); % sec
        rectLFP = interp1(LFPtimes, abs(filteredLFP), frametimes);  % This is probably not the best
        %   plot(LFPtimes,abs(LFP)); hold on; plot(frametimes,rectLFP,'m-');
        feval('STAmex',randnums(:),rectLFP);
        if (ishandle(hwait))
            waitbar(i/ntrials,hwait);
        else
            break;
        end
    end
    close(hwait);
    out = feval('STAmex','return');
    eval('clear STAmex');
    
    % Now doing some plotting
    STAs = out{1};
    STVs = out{2};
    STAlims = [-max(abs(STAs(:))) max(abs(STAs(:)))];
    STVlims = [min(STVs(:)) max(STVs(:))];
    STAs = (STAs-STAlims(2))/(STAlims(1)-STAlims(2));
    STVs = (STVs-STVlims(2))/(STVlims(1)-STVlims(2));
    
    figure; colormap(gray(255));
    set(gcf,'Name',['LFP ',num2str(whichADchan),' ',num2str(centerfreq),' Hz']);
    for i = 1:size(STAs,2)
        % First STAs
        mat = reshape(STAs(:,i),[nstixperside nstixperside 3]);
        for j = 1:3
            im = mat(:,:,j);
            subplot(6,size(STAs,2),(j-1)*size(STAs,2)+i);
            image(im*255);
            axis image;
            set(gca,'XTick',[],'YTick',[]);
        end
        % Then STVs
        mat = reshape(STVs(:,i),[nstixperside nstixperside 3]);
        for j = 1:3
            im = mat(:,:,j);
            subplot(6,size(STAs,2),(j-1+3)*size(STAs,2)+i);
            image(im*255);
            axis image;
            set(gca,'XTick',[],'YTick',[]);
        end
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The stuff below is basically obsolete.  It's stuff I was working on just
% to make sure that could do offline analysis at all.
%
% First, looking at eye movements
figure;
for i = 1:ntrials
    hep = WN.ras{i,hepidx}*4096/400;
    vep = WN.ras{i,vepidx}*4096/400;
    nsamps = length(hep);
    t = [0:1:nsamps-1]*eyesampperiod+eyestart_t(i);
    subplot(3,1,1); hold on;
    plot(hep,vep);
   % subplot(3,1,2); hold on;
   % plot(t-stimon_t(i), hep,'m-');
   % plot(t-stimon_t(i), vep,'g-');
end

subplot(3,1,1)
axis square
% Theory on the eye movements:
% in nex2stro.m we multiply the raw AD values by 10/2^12 to get to
% "millivolts".  But we know that on REX 40 AD values is 1 degree and we
% know that REX and Plexon have similar A/D cards.  Both Rex and Plexon
% appear to want A/D signals from -5V to +5V (they both saturate outside of
% this range).  
%
% 10 V/102.4 deg * 4096 steps/10 V = 40 steps/deg

%%
% Select a pixel for analysis
[x,y] = ginput;
whichpix = [round(x) round(y)];
idx = (whichpix(1)-1)*nstixperside+whichpix(2);
STA_t = STAs([idx, idx+nstixperside^2, idx+2*nstixperside^2],:);
figure; axes; hold on;
plot(STA_t(1,:),'r-','LineWidth',2);
plot(STA_t(2,:),'g-','LineWidth',2);
plot(STA_t(3,:),'b-','LineWidth',2);

%%
% Looking at the gun intensities of STA and PCs
% at a single frame
colors = {'red','green','blue'};
rowidxs = reshape([1:3*nstixperside^2],[nstixperside^2, 3]);
whichframe = input('Which frame would you like to see?');
figure;
subplot(6,1,1); hold on;
maxval = max(abs([STAs(:,whichframe);PCs(:,1,whichframe)]));
for gun = 1:3
    plot(STAs(rowidxs(:,gun),whichframe),'Color',colors{gun});
end
set(gca,'Ylim', [-1.5 1.5]*maxval);
for i = 1:4
    subplot(6,1,1+i); hold on;
    for gun = 1:3
        plot(PCs(rowidxs(:,gun),i,whichframe),'Color',colors{gun});
    end
    set(gca,'Ylim', [-1.5 1.5]*maxval);
end
STC = reshape(STCs(:,whichframe),[sqrt(length(STCs(:,whichframe))),sqrt(length(STCs(:,whichframe)))]);
STV = reshape(diag(STC),[nstixperside nstixperside 3]);
STV = (STV-mean(STV(:)));
subplot(6,1,6); hold on;
for gun = 1:3
    plot(reshape(STV(:,:,gun),[nstixperside^2 1]),'Color',colors{gun});
end

%% Looking trial by trial at eye movements, spikes, events, etc.
figure; axes; hold on;
framerate = WN.sum.exptParams.framerate;
predstimoff_t = stimon_t+WN.trial(:,nframesidx)/framerate;
for i = 1:ntrials
    timeanchor = stimon_t(i);
    hep = WN.ras{i,hepidx}*4096/400;
    vep = WN.ras{i,vepidx}*4096/400;
    spiketimes = WN.ras{i,spikeidx}-timeanchor;
    nspikes = length(spiketimes);
    nsamps = length(hep);
    t = [0:1:nsamps-1]*eyesampperiod+eyestart_t(i)-timeanchor;
    minep = min([hep; vep]);
    
    plot(t, hep,'m-');
    plot(t, vep,'g-');
    plot(stimon_t(i)-timeanchor,1,'g*');
    plot(stimoff_t(i)-timeanchor,1,'r*');
    plot(predstimoff_t(i)-timeanchor,1,'m*');
    plot([spiketimes spiketimes]',[zeros(nspikes,1) .1*ones(nspikes,1)]'+minep,'b-')
    
    pause;
    cla;
end

%% Looking at a PSTH of gun noise trials synced on stimulus onset
binwidth = .025;
bins =[-.5:binwidth:1.5];
psthmat = nan(ntrials, length(bins));
L = WN.trial(:,noisetypeidx) == 1;  % gun noise
stimon_t = WN.trial(:,stimonidx);
numframes = WN.trial(:,nframesidx)
stimoff_t = stimon_t+numframes/framerate;
for i = find(L)'
    if (stimoff_t(i)-stimon_t(i) > bins(end))
        spiketimes = WN.ras{i,spikeidx}-stimon_t(i);
        psthmat(i,:) = histc(spiketimes, bins);
    end
end

plot((bins+binwidth/2),nanmean(psthmat)/binwidth,'k-','Linewidth',2)
xlabel('Time(sec)');
ylabel('sp/sec');
set(gca,'Xlim',[bins(1) bins(end)]);
%% Code for making 1-D firing rate functions from 2-D data
figure;
subplot(2,1,1); hold on;
idx = round(size(frim,2)/2);
plot(frim(idx,:),'k-');
plot(frim(idx,:)+semim(idx,:),'k-.');
plot(frim(idx,:)-semim(idx,:),'k-.');
ylim = get(gca,'YLim');
set(gca,'XTickLabel',[],'XTick',[1:size(frim,1)]);

subplot(2,1,2); hold on;
idx = round(size(frim,1)/2);
plot(frim(:,idx),'k.-');
plot(frim(:,idx)+semim(:,idx),'k-.');
plot(frim(:,idx)-semim(:,idx),'k-.');
set(gca,'Ylim',ylim);
set(gca,'XTickLabel',[],'XTick',[1:size(frim,2)]);

figure;
imagesc(frim');
axis image;
set(gca,'XTick',[],'YTick',[]);
colormap(gray);
colorbar;