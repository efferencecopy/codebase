% Script for looking at data collected with the GaborEdge paradigm
stro = nex2stro;
ntrials = size(stro.trial,1);
ntrlspercond = stro.sum.exptParams.trialspercond;
fpacqidx = find(strcmp(stro.sum.trialFields(1,:),'fp_acq'));
stimonidx = find(strcmp(stro.sum.trialFields(1,:),'stim_on'));
stimoffnidx = find(strcmp(stro.sum.trialFields(1,:),'stim_off'));
fpacq_t = stro.trial(:,fpacqidx);
stimon_t = stro.trial(:,stimonidx);
numframes = stro.sum.exptParams.nframes;
framerate = stro.sum.exptParams.framerate;
stimoff_t = stimon_t+numframes/framerate;
thetaidx = find(strcmp(stro.sum.trialFields(1,:),'g_theta'));
lambdaidx = find(strcmp(stro.sum.trialFields(1,:),'g_lambda'));
phiidx = find(strcmp(stro.sum.trialFields(1,:),'g_phi'));
sigmaidx = find(strcmp(stro.sum.trialFields(1,:),'g_sigma'));
gammaidx = find(strcmp(stro.sum.trialFields(1,:),'g_gamma'));
xoffsetidx = find(strcmp(stro.sum.trialFields(1,:),'g_xoff'));
yoffsetidx = find(strcmp(stro.sum.trialFields(1,:),'g_yoff'));
gaborcontrastidx = find(strcmp(stro.sum.trialFields(1,:),'g_cont'));
gabortypeidx = find(strcmp(stro.sum.trialFields(1,:),'g_colortype'));
gaborrgbidxs = [find(strcmp(stro.sum.trialFields(1,:),'g_r'));...
                find(strcmp(stro.sum.trialFields(1,:),'g_g'));...
                find(strcmp(stro.sum.trialFields(1,:),'g_b'))];
gaborlmsidxs = [find(strcmp(stro.sum.trialFields(1,:),'g_l'));...
                find(strcmp(stro.sum.trialFields(1,:),'g_m'));...
                find(strcmp(stro.sum.trialFields(1,:),'g_s'))];
edgeangleidx = find(strcmp(stro.sum.trialFields(1,:),'e_theta'));
edgedispidx = find(strcmp(stro.sum.trialFields(1,:),'e_disp'));
edgecontrastidx = find(strcmp(stro.sum.trialFields(1,:),'e_cont'));
edgetypeidx = find(strcmp(stro.sum.trialFields(1,:),'e_colortype'));
edgergbidxs = [find(strcmp(stro.sum.trialFields(1,:),'e_r'));...
                find(strcmp(stro.sum.trialFields(1,:),'e_g'));...
                find(strcmp(stro.sum.trialFields(1,:),'e_b'))];
edgelmsidxs = [find(strcmp(stro.sum.trialFields(1,:),'e_l'));...
                find(strcmp(stro.sum.trialFields(1,:),'e_m'));...
                find(strcmp(stro.sum.trialFields(1,:),'e_s'))];

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
bkgndrgb = stro.sum.exptParams.bkgndrgb;
bkgndlms = M*bkgndrgb;

gaborcolortypes = stro.trial(:,gabortypeidx);
gaborcontrasts = stro.trial(:,gaborcontrastidx);
edgecontrasts = stro.trial(:,edgecontrastidx);
edgeangles = stro.trial(:,edgeangleidx);
edgedisplacements = stro.trial(:,edgedispidx);

spikename = getSpikenum(stro);
spikeidx = find(strcmp(stro.sum.rasterCells(1,:),spikename));
spikes = stro.ras(:,spikeidx);

%%
% Looking at rasters, condition by condition
nconds = size(unique([gaborcolortypes gaborcontrasts, edgecontrasts, edgeangles, edgedisplacements],'rows'),1);
means = [];
sds = [];
ns = [];
figure;
axescounter = 1;
for h = unique(gaborcolortypes)' 
    for i = unique(gaborcontrasts)' 
        for j = unique(edgecontrasts)'
            for k = unique(edgeangles)'
                for l = unique(edgedisplacements)'
                    L = (gaborcolortypes == h) & (gaborcontrasts == i) & (edgecontrasts == j) &...
                        (edgeangles == k) & (edgedisplacements == l);
                    if (sum(L) > 0)
                        tmp = [];
                        subplot(ceil(sqrt(nconds)),ceil(sqrt(nconds)), axescounter);
                        hold on;
                        trlidxs = find(L);
                        for counter = 1:sum(L)
                            trlidx = trlidxs(counter);
                            plot(0,counter,'g*');
                            plot(stimoff_t(trlidx)-stimon_t(trlidx),counter,'r*');
                            plot(fpacq_t(trlidx)-stimon_t(trlidx),counter,'m*');
                            nspikestot = length(spikes{trlidx});
                            plot([spikes{trlidx} spikes{trlidx}]'-stimon_t(trlidx),[zeros(nspikestot,1) .5*ones(nspikestot,1)]'+counter,'k-');
                            title(sprintf('%1.1f %1.2f %1.1f ',i, j, k));
                            nspikes = sum((spikes{trlidx} < stimoff_t(trlidx)) &...
                                          (spikes{trlidx} > stimon_t(trlidx)));
                            tmp = [tmp; nspikes];
                        end
                        means = [means; mean(tmp)];
                        sds = [sds; sqrt(var(tmp))];
                        ns = [ns; length(tmp)];
                        axescounter = axescounter+1;
                        set(gca,'XTick',[]','YTick',[]);
                        set(gca,'XLim',[fpacq_t(trlidx)-stimon_t(trlidx)-.1, stimoff_t(trlidx)-stimon_t(trlidx)+.1]);
                    end
                end
            end
        end
    end
end
figure;
errorbar(means, sds./sqrt(ns));

% gray scale plot of firing rates
figure;
for i = 1:axescounter-1
    subplot(ceil(sqrt(nconds)),ceil(sqrt(nconds)), i);
    image(255*(means(i)-min(means))/(max(means)-min(means)));
    colormap(gray(255));
    axis image;
    set(gca,'Visible','off');
end
set(gcf,'Color',[0 0 0]);
%%
% A visual representation of each gabor/edge stimulus.
% Assumes a single theta, lambda, phi, sigma, gamma,
% gabor colortype, and edge color type.
theta = unique(stro.trial(:,thetaidx));
lambda = unique(stro.trial(:,lambdaidx));
phi = unique(stro.trial(:,phiidx));
sigma = unique(stro.trial(:,sigmaidx));
gamma = unique(stro.trial(:,gammaidx));
gausslim = stro.sum.exptParams.gausslim/1000;
pixperdeg = stro.sum.exptParams.pixperdeg;
nconds = size(unique([gaborcolortypes gaborcontrasts, edgecontrasts, edgeangles, edgedisplacements],'rows'),1);

figure;
axescounter = 1;
for h = unique(gaborcolortypes)'
    for i = unique(gaborcontrasts)'
        for j = unique(edgecontrasts)'
            for k = unique(edgeangles)'
                for l = unique(edgedisplacements)'
                    L = (gaborcolortypes == h) & (gaborcontrasts == i) & (edgecontrasts == j) &...
                        (edgeangles == k) & (edgedisplacements == l);
                    if (sum(L) > 0)
                        subplot(ceil(sqrt(nconds)),ceil(sqrt(nconds)), axescounter);
                        gaborrgb = unique(stro.trial(L,gaborrgbidxs),'rows');
                        edgergb = unique(stro.trial(L,edgergbidxs),'rows');
                        im = DrawGaborEdge(bkgndrgb, gaborrgb, edgergb, theta, lambda, sigma, gamma, phi, 0, 0, k, l, gausslim, pixperdeg);
                        image(im);
                        set(gca,'Visible','off');
                        axis image; axis tight;
                        axescounter = axescounter +1;
                    end
                end
            end
        end
    end
end
%%
% Looking at 1-D marginal averages
spikerates = [];
for i = 1:ntrials
    spiketimes = stro.ras{i,spikeidx};
    nspikes = sum(spiketimes > stimon_t(i) & spiketimes < stimoff_t(i));
    spikerates = [spikerates; nspikes./(stimoff_t(i)-stimon_t(i))];
end

% Gabor contrast
labels = {'STA','orth1','orth2','orth3','Achrom','L-M','mech'};
cols = {'m','b','g','c','k','r','y'};
figure; axes; hold on; 
xlabel('Gabor rms contrast'); ylabel('spikes/sec');
for j = unique(gaborcolortypes)'
    data = [];
    for i = unique(gaborcontrasts)'
        L = (gaborcolortypes == j) & (gaborcontrasts == i);
        tmp = spikerates(L);
        data = [data; mean(tmp), sqrt(var(tmp)./length(tmp))];
    end
    if (size(data,1) > 1)
        h = [h; errorbar(unique(gaborcontrasts),data(:,1),data(:,2),cols{j+1})];
    end
end
if (length(unique(gaborcolortypes)) > 1)
    legend (labels{unique(gaborcolortypes)+1});
end

% Edge contrast
data = [];
for i = unique(edgecontrasts)'
    L = (edgecontrasts == i);
    tmp = spikerates(L);
    data = [data; mean(tmp), sqrt(var(tmp)./length(tmp))];
end
if (size(data,1) > 1)
    figure; axes; errorbar(unique(edgecontrasts),data(:,1),data(:,2)); 
    xlabel('Edge contrast');
    ylabel('spikes/sec');
end

% Edge orientation
data = [];
for i = unique(edgeangles)'
    for j = [1,-1]
        L = (edgeangles == i) & (sign(edgecontrasts) == j);
        tmp = spikerates(L);
        data = [data; mean(tmp), sqrt(var(tmp)./length(tmp))];
    end
end
angs = unique(edgeangles);
angs = [angs, angs+pi]';
angs = angs(:);
[angs,i] = sort(angs);
data= data(i,:);
if (size(data,1) > 1)
    figure; axes; 
    errorbar(angs,data(:,1),data(:,2)); 
    xlabel('Edge angles');
    ylabel('spikes/sec');
end

% Edge displacements
data = [];
for i = unique(edgedisplacements)'
    L = (edgedisplacements == i);
    tmp = spikerates(L);
    data = [data; mean(tmp), sqrt(var(tmp)./length(tmp))];
end
if (size(data,1) > 1)
    figure; axes; errorbar(unique(edgedisplacements),data(:,1),data(:,2));
    xlabel('Edge displacements');
    ylabel('spikes/sec');
end

%%
% Firing rate surface collapsing over edge orientation and displacement
uniquegaborcolortypes = unique(gaborcolortypes);
if (length(uniquegaborcolortypes) > 1)
    a = listdlg('PromptString','Which gabor colortype?','SelectionMode','Single','ListSize',[150 50],'ListString',num2str(uniquegaborcolortypes));
    whichgaborcolortype = uniquegaborcolortypes(a);
else
    whichgaborcolortype = uniquegaborcolortypes;
end

spikecounts = [];
for i = 1:ntrials
    spiketimes = stro.ras{i,spikeidx};
    nspikes = sum(spiketimes > stimon_t(i) & spiketimes < stimoff_t(i));
    spikecounts = [spikecounts; nspikes];
end

uniquegcont= unique(gaborcontrasts)';
uniqueecont = unique(edgecontrasts)';
frim = nan*ones(length(uniqueecont), length(uniquegcont));
for i = 1:length(uniquegcont)
    for j = 1:length(uniqueecont)
        L = (gaborcontrasts == uniquegcont(i)) & (edgecontrasts == uniqueecont(j));
        L = L & (gaborcolortypes == whichgaborcolortype)
        frim(j,i) = mean(spikecounts(L));
    end
end
figure
imagesc(frim);
colormap(gray);
set(gca,'YDir','normal');
xlabel('gabor contrast');
ylabel('edge contrast');
set(gca,'XTick',[1:length(uniquegcont)],'XTickLabel',uniquegcont);
set(gca,'YTick',[1:length(uniqueecont)],'YTickLabel',uniqueecont);

%%
% 2-D contrast response functions conditioned on "edge style"
uniquegaborcolortypes = unique(gaborcolortypes);
if (length(uniquegaborcolortypes) > 1)
    a = listdlg('PromptString','Which gabor colortype?','SelectionMode','Single','ListSize',[150 50],'ListString',num2str(uniquegaborcolortypes));
    whichgaborcolortype = uniquegaborcolortypes(a);
else
    whichgaborcolortype = uniquegaborcolortypes;
end


spikecounts = [];
for i = 1:ntrials
    spiketimes = stro.ras{i,spikeidx};
    nspikes = sum(spiketimes > stimon_t(i) & spiketimes < stimoff_t(i));
    spikecounts = [spikecounts; nspikes];
end

uniquegcont= unique(gaborcontrasts)';
uniqueecont = unique(edgecontrasts)';
designmatrix = [edgeangles edgedisplacements];
trialtypes = unique(designmatrix,'rows');
figure;
frim = zeros(length(uniqueecont), length(uniquegcont), size(trialtypes,1));
tmp = [];
for k = 1:size(trialtypes,1)
    for i = 1:length(uniquegcont)
        for j = 1:length(uniqueecont)
            L = (gaborcontrasts == uniquegcont(i)) & (edgecontrasts == uniqueecont(j));
            if (uniqueecont(j) ~= 0)
                L = L & (edgeangles == trialtypes(k,1)) & (edgedisplacements == trialtypes(k,2));
            end
            if (uniquegcont(i) ~= 0)
               L = L & (gaborcolortypes == whichgaborcolortype);
            end
            frim(j,i,k) = mean(spikecounts(L));
            tmp = [tmp; find(L)];
        end
    end
end

gaborrgb = unique(stro.trial(L,gaborrgbidxs),'rows');
maxfr = max(frim(:));
for k = 1:size(trialtypes,1)
    subplot(length(unique(edgeangles)), length(unique(edgedisplacements)),k)
    image(255*frim(:,:,k)./maxfr);
    colormap(gray(255));
    set(gca,'YDir','normal');
    set(gca,'XTick',[1:length(uniquegcont)],'XTickLabel',uniquegcont);
    set(gca,'YTick',[1:length(uniqueecont)],'YTickLabel',uniqueecont);

    if (k == 1)
        xlabel('gabor contrast');
        ylabel('edge contrast');
    end
    %title(['A: ',num2str(trialtypes(k,1)),' D: ',num2str(trialtypes(k,2))]);
    axis square;
end
h = colorbar('EastOutside')
set(h,'YTick',[1 size(colormap,1)],'YTickLabel',[0 maxfr])

% Now plotting images
figure;
for k = 1:size(trialtypes,1)
    im = DrawGaborEdge(bkgndrgb, gaborrgb, edgergb, theta, lambda, sigma, gamma, phi, 0, 0, trialtypes(k,1), trialtypes(k,2), gausslim, pixperdeg);
    subplot(length(unique(edgeangles)), length(unique(edgedisplacements)),k)
    image(im);
    set(gca,'Visible','off');
    axis image; axis tight;
end


%%
% Randomization tests
% Before publishing anything from here, better make sure it's working
% properly.
spikecounts = [];
for i = 1:ntrials
    spiketimes = stro.ras{i,1};
    nspikes = sum(spiketimes > stimon_t(i) & spiketimes < stimoff_t(i));
    spikecounts = [spikecounts; nspikes];
end
alltrials = [gaborcontrasts edgecontrasts edgeangles edgedisplacements];
strings = {'Gabor Contrast','Edge Contrast','Edge Angle','Edge Displacement'};

niter = 1000;
for whichtest = 1:size(alltrials,2)
    trialtypes = unique(alltrials,'rows');
    subidxs = [1:size(alltrials,2)];
    subidxs(whichtest) = [];
    subtrialtypes = unique(alltrials(:,subidxs),'rows');

    data = [];
    order = [1:length(spikecounts)]';
    for iter = 1:niter
        tmp = [];
        for i = unique(trialtypes(:,whichtest))'
            L = (alltrials(:,whichtest) == i);
            tmp = [tmp; mean(spikecounts(order(L))) var(spikecounts(order(L)))];
        end
        F = var(tmp(:,1))./mean(tmp(:,2));
        data = [data; F];
        
        % Shuffling data
        order = [1:length(spikecounts)]';
        for i = 1:size(subtrialtypes,1)
            L = ones(size(spikecounts,1),1);
            for j = 1:size(subtrialtypes,2)
                L = L & alltrials(:,subidxs(j)) == subtrialtypes(i,j);
            end
            whichtrials = find(L);
            order(whichtrials) = whichtrials(randperm(length(whichtrials)));
        end
    end
    p = sum(data(2:end) >= data(1))/(size(data,1)-1);
    disp([strings{whichtest},' p:',num2str(p)]);
end

%%
% PSTH of activity leading up to stimulus on - I'm pretty sure that at
% least some of these cells act like they're "anticipating" the stimulus
% presentation.

binwidth = .01;
bins = [-1:binwidth:.2];
psth = zeros(1,length(bins));
t_fpacq = zeros(length(ntrials),1);

for trlidx = 1:ntrials
    spiketimes = spikes{trlidx}-stimon_t(trlidx);
    spiketimes(spiketimes<bins(1)-binwidth/2) = [];
    spiketimes(spiketimes>bins(end)+binwidth/2) = [];
    [n,x] = hist(spiketimes,bins);
    psth = psth+n;
    t_fpac(trlidx) = fpacq_t(trlidx)-stimon_t(trlidx);
end
 
figure;
plot(bins,psth);
hold on;
plot(median(t_fpac),max(psth),'m*');
plot([min(t_fpac) max(t_fpac)],[max(psth) max(psth)],'m-');

%%
% Verifying that the RGB projection stuff is all working out
% First reconstructing the cone contrasts
gaborrgb = stro.trial(:,gaborrgbidxs);
gaborlms = stro.trial(:,gaborlmsidxs);
stim = gaborrgb+repmat(bkgndrgb',ntrials,1);
for i = 1:ntrials
    conexc = M*stim(i,:)';
    stim(i,:) = (conexc-bkgndlms)./bkgndlms;
end
figure;
plot(stim-gaborlms) % residuals should be small (like 10^-16)

% Then looking at the projections
rawrgb = stro.sum.exptParams.gaborrawrgb;
disp(['gabor colortype is ',num2str(unique(stro.trial(:,gabortypeidx)))]);
disp(['edge colortype is ',num2str(unique(stro.trial(:,edgetypeidx)))]);
if (unique(stro.trial(:,gabortypeidx)) == 1)
    plot(gaborlms*[1 1 0]');
end

% Checking the gabor contrast
residual = sqrt(mean(stim.^2,2))-stro.trial(:,gaborcontrastidx);
figure;
plot(residual); % residuals should be small (like 10^-9)

% Checking the relationship between gabor RGBs and raw rgb
tmp = gaborrgb./repmat(sqrt(sum(gaborrgb.^2,2)),1,3);
residual = tmp - repmat(rawrgb',ntrials,1);
figure;
plot(residual) % residuals should be small (like 10^-16)
%%
% Checking to make sure the cone contrasts are being reported correctly
% First, getting the M matrix
spd = stro.sum.exptParams.mon_spd;
fund = stro.sum.exptParams.fundamentals;
spd = reshape(spd,length(spd)/3,3);
fund = reshape(fund,length(fund)/3,3);
spd = SplineSpd(SToWls([380 4 101]),spd, SToWls([380 5 81]));
M = fund'*spd;
bkgndlms = M*stro.sum.exptParams.bkgndrgb;

% Calculating cone contrasts from the rgbs reported
deltacone = (M*stro.trial(:,[15 16 17 ])')';
a = deltacone./repmat(bkgndlms',size(deltacone,1),1);
b = stro.trial(:,[18 19 20]);
plot(a(:)-b(:))  % residuals should be small, like 10^-17

% cone contrasts as reported
% compared to rms cone contrast reported
cc = stro.trial(:,[18 19 20]);
a = sqrt(mean(cc.^2,2));
b = stro.trial(:,13);
[a b] 
plot(a-b); % residuals should be small, like 10^-9
%%
% Getting color directions in cone contrast space
labels = {'STA','orth1','orth2','orth3','Achrom','L-M','mech'};

for i = unique(gaborcolortypes)'
    L = gaborcolortypes == i;
    ccs = stro.trial(L,[18 19 20]);
    colordirs = ccs./repmat(sum(sqrt(ccs.^2),2),1,3);
    if(any(std(abs(colordirs)) > 10^-10))
        error(['More than one color direction detected for gabor color condition ',num2str(i)]);
    end
    colordir = nanmean(colordirs);
    disp([labels{i+1},' ',num2str(colordir)]);
end
