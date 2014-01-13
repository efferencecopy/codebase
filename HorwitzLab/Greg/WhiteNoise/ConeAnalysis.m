% Script for comparing STAs measured under gun noise and cone noise

stro=nex2stro;
framerate = stro.sum.exptParams.framerate;
ntrials = length(stro.sum.absTrialNum);
stimonidx = find(strcmp(stro.sum.trialFields(1,:),'stim_on'));
stimoffidx = find(strcmp(stro.sum.trialFields(1,:),'all_off'));
nframesidx = find(strcmp(stro.sum.trialFields(1,:),'num_frames'));
noisetypeidx = find(strcmp(stro.sum.trialFields(1,:),'noise_type'));
sigmaidxs = strmatch('sigma',stro.sum.trialFields(1,:));

hepidx = find(strcmp(stro.sum.rasterCells(1,:),'AD11'));
vepidx = find(strcmp(stro.sum.rasterCells(1,:),'AD12'));
anlgStartTimeidx = find(strcmp(stro.sum.rasterCells(1,:),'anlgStartTime'));
eyestart_t = [stro.ras{:,anlgStartTimeidx}]';
eyesampperiod = 1/stro.sum.analog.storeRates{1};

gammaTable = stro.sum.exptParams.gamma_table;
gammaTable = reshape(gammaTable, length(gammaTable)/3, 3);
gammaTable1 = interp1(linspace(0,255,256),gammaTable,linspace(0,255,65536), 'spline');
invgamma = InvertGamma(gammaTable, 0);

% Reconstructing the M matrix
fundamentals = stro.sum.exptParams.fundamentals;
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = stro.sum.exptParams.mon_spd;
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;

% Getting the background rgb/lms
ridx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_r'));
gidx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_g'));
bidx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_b'));
bkgndRGB = [mode(stro.trial(:,ridx)), mode(stro.trial(:,gidx)), mode(stro.trial(:,bidx))];
bkgndrgb = [gammaTable(bkgndRGB(1)+1,1); gammaTable(bkgndRGB(2)+1,2); gammaTable(bkgndRGB(3)+1,3)];
bkgndlms = M*bkgndrgb;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Getting STAs
spikename = getSpikenum(stro);
spikeidx = find(strcmp(stro.sum.rasterCells(1,:),spikename));
maxT = 8;
nstixperside = stro.sum.exptParams.nstixperside;
for noisetype = 1:2
    tmpstro = stro;
    L = stro.trial(:,noisetypeidx) == noisetype;
    tmpstro.ras(~L,:) = [];
    tmpstro.trial(~L,:) = [];
    out = getWhtnsStats(tmpstro,maxT,'STCOVmex', {nstixperside^2, 3, maxT}, spikename);
    if (noisetype == 1)
        STAs_gun = out{1};
        STCs_gun = out{2};
        nspikes_gun = out{3};
    elseif (noisetype == 2)
        STAs_cone = out{1};
        STCs_cone = out{2};
        nspikes_cone = out{3};        
    end
end
tmpstro = [];


%%
% Cone weights maps from cone noise
figure;
rng = [min(STAs_cone(:)) max(STAs_cone(:))];
STAs_cone_z = STAs_cone./sqrt(nspikes_cone);
alpha = [0.001];  % max of 3 critical values
crit = norminv(1-alpha,0,1);
for i = 1:size(STAs_cone,2)
    for j = 1:3
        rowidxs = [(j-1)*nstixperside^2+1, j*nstixperside^2];
        subplot(3,size(STAs_cone,2),i+size(STAs_cone,2)*(j-1));
        im = reshape(STAs_cone(rowidxs(1):rowidxs(2),i),[nstixperside nstixperside]);    
        im = (im-rng(1))/(rng(2)-rng(1));  % rescaling
        im = repmat(im, [1 1 3]);
        
        zim = reshape(STAs_cone_z(rowidxs(1):rowidxs(2),i),[nstixperside nstixperside]);
        for k = 1:length(alpha)
            pmat = logical(abs(zim)>crit(k));
            sigidxs = find(pmat);
            im(sigidxs+nstixperside^2*(k-1)) = .5;
        end
        image(im);
        
        axis image;
        set(gca,'XTick',[],'YTick',[]);
    end
end

%%
% Comparing estimates of cone weights from gun noise and cone noise.
% At the moment, not normalizing either one by the background lms since
% that will have the same effect on both estimates and so isn't useful
% for the comparison between them.
% Finding whichpixels
sdthresh = 3;
tmpSTA = reshape(STAs_gun, [nstixperside^2 3 maxT]);
tmpSTA = permute(tmpSTA, [2 1 3]);
STAgunmat = reshape(tmpSTA,[3 nstixperside^2*maxT]);
tmpSTA = reshape(STAs_cone, [nstixperside^2 3 maxT]);
tmpSTA = permute(tmpSTA, [2 1 3]);
STAconemat = reshape(tmpSTA,[3 nstixperside^2*maxT]);
conesigmas = stro.trial(stro.trial(:,noisetypeidx) == 2, sigmaidxs)/1000;  % cone excitation units

whichpixgun = logical(sum(STAgunmat.^2)>sdthresh*std(sum(STAgunmat.^2)));
whichpixcone = logical(sum(STAconemat.^2)>sdthresh*std(sum(STAconemat.^2)));
whichpix = find(whichpixcone | whichpixgun);

% First, the gun noise STA
[u,s,v] = svd(STAgunmat(:,whichpix));
if (sum(v(:,1)) < 0)
    u = -u;
end
oldv = v(:,1);
coneweights = inv(M')*u(:,1);
%coneweights = coneweights.*bkgndlms;   % Cone specific adaptation (on the weights?)
%coneweights = coneweights.*(M*[1 1 1]'); 
coneweights_gun = coneweights./sum(abs(coneweights));

% Cone weights from cone noise expressed in contrast units
%[u,s,v] = svd(STAconemat(:,whichpix)./repmat(lmscontrasts',[1 length(whichpix)]));%
[u,s,v] = svd(STAconemat(:,whichpix));
if (v(:,1)'*oldv< 0 | sum(oldv) < 0)
    u = -u;
end
u(:,1) = u(:,1)./mode(abs(conesigmas))';
coneweights_cone = u(:,1)./sum(abs(u(:,1)));
figure; axes; hold on;
bar([coneweights_cone coneweights_gun]);
legend({'cone','gun'});
set(gca,'XTick',[1 2 3],'XTickLabel',{'L','M','S'});

%%
% Rasters for gun noise and cone noise

figure; axes; hold on;
counter = 0;
maxstimduration = 0;
for i = 1:2
    Lnoise = stro.trial(:,noisetypeidx) == i;
    stimon_t = stro.trial(Lnoise,stimonidx);
    numframes = stro.trial(Lnoise,nframesidx);
    stimoff_t = stimon_t+numframes/framerate;
    spikes = stro.ras(Lnoise,spikeidx);
    stimdurations = stimoff_t-stimon_t;
    maxstimduration = max([maxstimduration; stimdurations]);
    [junk, sortedidxs] = sort(stimdurations);
    
    for j = sortedidxs'
        plot(0,counter,'g*');
        plot(stimoff_t(j)-stimon_t(j),counter,'r*');
        nspikes = length(spikes{j});
        plot([spikes{j} spikes{j}]'-stimon_t(j),[zeros(nspikes,1) .5*ones(nspikes,1)]'+counter,'k-')
        counter = counter + 1;
    end
    plot([-.2 maxstimduration+.2],[counter counter],'k-');
end
set(gca,'Xlim',[-.2 maxstimduration+.2]); xlabel('time(s)');
title(stro.sum.fileName);

%%
% Looking at the PCs from the cone noise data.
% Top three rows are largest PC, bottom three rows are smallest
% Projecting orthogonal to STA - GDLH 9/28/09
figure;
ndims = sqrt(size(STCs_cone,1));
for i = 1:size(STCs_cone,2)
   STC = reshape(STCs_cone(:,i),[ndims ndims]);
   P = eye(size(STAs_cone,1),size(STAs_cone,1))-STAs_cone(:,i)*inv(STAs_cone(:,i)'*STAs_cone(:,i))*STAs_cone(:,i)';
   [v,d] = eig(P*STC*P');
   [sorted, sortidxs] = sort(diag(d));
   v = v(:,sortidxs);
   tmp = reshape(v(:,end),[sqrt(ndims/3) sqrt(ndims/3) 3])
   tmp = tmp./(2*max(abs(tmp(:))));
   tmp = tmp+.5;
   for j = 1:3
       subplot(7,size(STCs_cone,2),(j-1)*size(STCs_cone,2)+i)
       image(255*tmp(:,:,j));
       colormap(gray(255));
       axis image;
       set(gca,'Xtick',[],'Ytick',[]);
   end
   tmp = reshape(v(:,2),[sqrt(ndims/3) sqrt(ndims/3) 3])
   tmp = tmp./(2*max(abs(tmp(:))));
   tmp = tmp+.5;
   for j = 1:3
       subplot(7,size(STCs_cone,2),(j-1)*size(STCs_cone,2)+i+(size(STCs_cone,2)*4))
       image(255*tmp(:,:,j));
       colormap(gray(255));
       axis image;
       set(gca,'Xtick',[],'Ytick',[]);
   end
end
%%
% Rendering an RGB image from the LMS STA.  Scaling cone STA by contrast.
% ignoring the possiblity that sigma or mu might have changed across
% trials.  Or that mu might not be zero.
gunsigmas = stro.trial(stro.trial(:,noisetypeidx) == 1, sigmaidxs)/1000;  % normalized gun intensity units
conesigmas = stro.trial(stro.trial(:,noisetypeidx) == 2, sigmaidxs)/1000;  % cone excitation units
lmscontrasts = mode(conesigmas)./bkgndlms';

tmpSTA = reshape(STAs_cone, [nstixperside^2 3 maxT]);
tmpSTA = permute(tmpSTA, [2 1 3]);
STAconemat = reshape(tmpSTA,[3 nstixperside^2*maxT]);

% Plotting a few transformations of the cone STAs
figure; set(gcf,'DefaultAxesColorOrder',[1 0 0; 0 1 0; 0 0 1]);
subplot(3,3,1); plot(STAconemat');  % raw lms STA
subplot(3,3,2); plot((STAconemat./repmat(mode(conesigmas)',[1 size(STAconemat,2)]))');  % lms STA normalized by sigma
subplot(3,3,3); plot((STAconemat./repmat(lmscontrasts',[1 size(STAconemat,2)]))');  % lms STA normalized by contrast

subplot(3,3,4); plot((M'*STAconemat)');  % raw lms->rgb STA
subplot(3,3,5); plot((M'*(STAconemat./repmat(mode(conesigmas)',[1 size(STAconemat,2)])))');  % lms STA normalized by sigma converted to rgb
subplot(3,3,6); plot((M'*(STAconemat./repmat(lmscontrasts',[1 size(STAconemat,2)])))');  % lms STA normalized by contrast converted to rgb

subplot(3,3,8);
tmpSTA = reshape(STAs_gun, [nstixperside^2 3 maxT]);
tmpSTA = permute(tmpSTA, [2 1 3]);
STAgunmat = reshape(tmpSTA,[3 nstixperside^2*maxT]);
plot(STAgunmat')

% Plotting the cone STAs converted to rgb
figure;
STAlmstorgb = M'*STAconemat./repmat(lmscontrasts',[1 size(STAconemat,2)]);
STAlmstorgb = STAlmstorgb./(2*range(STAlmstorgb(:)));
STAlmstorgb = STAlmstorgb + repmat([.5 .5 .5]',[1 size(STAlmstorgb,2)]);
STAlmstorgb = permute(reshape(STAlmstorgb,[3 nstixperside^2 maxT]),[2 1 3]);
for i = 1:maxT
    subplot(2,maxT,i);
    image(reshape(STAlmstorgb(:,:,i),[nstixperside, nstixperside, 3]))
    axis image;
    set(gca,'XTick',[],'Ytick',[]);
end

% Plotting the measured rgb STA alongside the computed one
STArgb = STAgunmat./(2*range(STAgunmat(:)));
STArgb = STArgb + repmat([.5 .5 .5]',[1 size(STArgb,2)]);
STArgb = permute(reshape(STArgb,[3 nstixperside^2 maxT]),[2 1 3]);
for i = 1:maxT
    subplot(2,maxT,i+maxT);
    image(reshape(STArgb(:,:,i),[nstixperside, nstixperside, 3]))
    axis image;
    set(gca,'XTick',[],'Ytick',[]);
end
