% modelling stuff for specific aim 1
% The key here is going to be to find two cells: one that is very sensitive
% to L-M and another that responds to L-M, but weakly.  Then, fit
% neurometric models to each and see how many of the less sensitive cells
% we have to pool (with realistic correlations) to rival the sensitivity of
% a single very sensitive neuron.

filenames = {'K111308003.nex','K111908003.nex'};
spikenames = {'sig001b', 'sig001a'};
plotcol = [0 0 0; .4 0 1];
figure; axes; hold on;
ylabel('Area under ROC curve');
xlabel('L-M rms cone contrast');
poststimtime = 0.05; % seconds after stimoff to include
thresholds = [];
for fileidx = 1:2
    stro = nex2stro(findfile(filenames{fileidx}));
    spikename = spikenames{fileidx};

    stimonidx = find(strcmp(stro.sum.trialFields(1,:),'stim_on'));
    numframes = stro.sum.exptParams.nframes;
    framerate = stro.sum.exptParams.framerate;
    gaborcontrastidx = find(strcmp(stro.sum.trialFields(1,:),'g_cont'));
    gabortypeidx = find(strcmp(stro.sum.trialFields(1,:),'g_colortype'));
    gaborlmsidxs = [find(strcmp(stro.sum.trialFields(1,:),'g_l'));...
                    find(strcmp(stro.sum.trialFields(1,:),'g_m'));...
                    find(strcmp(stro.sum.trialFields(1,:),'g_s'))];
    spikeidx = find(strcmp(stro.sum.rasterCells(1,:),spikename));
    gaborcolortypes = stro.trial(:,gabortypeidx);
    Lcolortype = gaborcolortypes == 5;  % L-M
    Lcontrast = stro.trial(:,gaborcontrastidx) <= 0.05;
    % Cutting out all trials except the ones we are interested in 
    stro.trial = stro.trial(Lcolortype & Lcontrast,:);
    stro.ras = stro.ras(Lcolortype & Lcontrast,:);

    gaborcontrasts = stro.trial(:,gaborcontrastidx);
    Lnull = all(stro.trial(:,[18 19 20]) == 0,2);  % blanks
    gaborcontrasts(Lnull) = 0;

    prestimspikes = nan*ones(sum(Lcolortype & Lcontrast),1);
    spikes = nan*ones(sum(Lcolortype & Lcontrast),1);
    for i = 1:length(spikes)
        spikes(i) = sum(stro.ras{i,spikeidx} > stro.trial(i,stimonidx) & ...
            stro.ras{i,spikeidx} < stro.trial(i,stimonidx)+numframes./framerate+poststimtime);
        prestimspikes(i) = sum(stro.ras{i,spikeidx} < stro.trial(i,stimonidx) & ...
            stro.ras{i,spikeidx} > stro.trial(i,stimonidx)-numframes./framerate-poststimtime);
    end

    spikedistns = {};
    uniquegaborcontrasts = unique(gaborcontrasts);
    ns = [];
    for i = 1:length(uniquegaborcontrasts)
        L = gaborcontrasts == uniquegaborcontrasts(i);
        ns(i) = sum(L);
        spikedistns{i} = spikes(L);
    %    idxs = fullfact(repmat(sum(L),ncells,1));
    %    tmpspikes = spikes(L);
    %    spikedistns{i} = mean(tmpspikes(idxs),2);
    end

    rocs = [];
    for i = 1:length(spikedistns)
        rocs(i) = roc(prestimspikes, spikedistns{i});
    end

    [a,b,g,s] = weibullFit(uniquegaborcontrasts, rocs, 'sse',[mean(uniquegaborcontrasts), 1 1]);
    [a,b,g,s] = weibullFit(uniquegaborcontrasts, [ns'.*rocs',ns'.*(1-rocs)'], 'mle', [a b]);

    thresholds(fileidx) = a;
    plot(uniquegaborcontrasts,rocs,'o','MarkerSize',4,...
        'MarkerFaceColor',plotcol(fileidx,:),'MarkerEdgeColor',plotcol(fileidx,:))
    x = linspace(.001,.05,1000);
    plot(x, g +(0.5-g).*exp(-((x./a).^b)),'-','Linewidth',2,'Color',plotcol(fileidx,:));
end

%%
% Plotting a rasters for a few selected conditions for plot insets.
% Showing that some neurons don't fire at contrasts where other neurons  
% give a reliable signal.
datapath = 'C:\NexFiles\Kali';
filenames = {'K111308003.nex','K111908003.nex'};
spikenames = {'sig001b', 'sig001a'};
whichcontrasts = {[2  6],[1 5]}; % from lowest = 1
axescounter = 1;
figure;
poststimtime = 0;
for fileidx = 1:length(filenames)
    stro = nex2stro([datapath,'/',filenames{fileidx}]);
    spikename = spikenames{fileidx};

    fpacqidx = find(strcmp(stro.sum.trialFields(1,:),'fp_acq'));
    fpacq_t = stro.trial(:,fpacqidx);
    stimonidx = find(strcmp(stro.sum.trialFields(1,:),'stim_on'));
    stimoffnidx = find(strcmp(stro.sum.trialFields(1,:),'stim_off'));
    stimon_t = stro.trial(:,stimonidx);
    numframes = stro.sum.exptParams.nframes;
    framerate = stro.sum.exptParams.framerate;
    stimoff_t = stimon_t+numframes/framerate;   gaborcontrastidx = find(strcmp(stro.sum.trialFields(1,:),'g_cont'));
    gabortypeidx = find(strcmp(stro.sum.trialFields(1,:),'g_colortype'));
    gaborlmsidxs = [find(strcmp(stro.sum.trialFields(1,:),'g_l'));...
                    find(strcmp(stro.sum.trialFields(1,:),'g_m'));...
                    find(strcmp(stro.sum.trialFields(1,:),'g_s'))];
    spikeidx = find(strcmp(stro.sum.rasterCells(1,:),spikename));
    gaborcolortypes = stro.trial(:,gabortypeidx);
    Lcolortype = gaborcolortypes == 5;  % L-M
    gaborcontrasts = stro.trial(:,gaborcontrastidx);
    uniquegaborcontrasts = unique(gaborcontrasts);
    for contrastidx = whichcontrasts{fileidx}
        Lcontrast = stro.trial(:,gaborcontrastidx) == uniquegaborcontrasts(contrastidx);
        subplot(2,2,axescounter); hold on;
        trlidxs = find(Lcontrast & Lcolortype);
        for counter = 1:sum(Lcontrast & Lcolortype)
            trlidx = trlidxs(counter);
            spikes = stro.ras{trlidx,spikeidx};
            spikes(spikes < fpacq_t(trlidx)) = [];
            spikes(spikes > stimoff_t(trlidx)+.2) = [];
            nspikestot = length(spikes);
            plot([spikes spikes]'-stimon_t(trlidx),[zeros(nspikestot,1) .5*ones(nspikestot,1)]'+counter-1,'k-');
        end
        set(gca,'Xlim',[-.2 numframes./framerate+.2]);
        set(gca,'Ylim',[0 counter],'Ytick',[]);
        plot([0 0],[0 counter],'g-');
        plot([numframes./framerate numframes./framerate],[0 counter],'r-');
        plot([numframes./framerate+poststimtime numframes./framerate+poststimtime],[0 counter],'r:');
        
        axescounter = axescounter+1;
    end
end
%%
% A figure showing the STA, fitted gabor (ortho1 and achrom) and contrast
% response functions.  Assumes there's no cone noise in the WN file.

WNstro = nex2stro('C:\NexFiles\Kali\K111308001.nex');
GBstro = nex2stro('C:\NexFiles\Kali\K111308002.nex');
spikename = 'sig001b';

%WNstro = nex2stro('C:\NexFiles\Kali\K111908002.nex');
%GBstro = nex2stro('C:\NexFiles\Kali\K111908003.nex');
%spikename = 'sig001a';

% Getting STA
nstixperside = WNstro.sum.exptParams.nstixperside;
out = getWhtnsStats(WNstro,9,'STCOVmex', {nstixperside^2, 3, 9}, spikename);
STAs = out{1};
energy = sum(STAs.^2);
peakframe = find(energy == max(energy));
STA = reshape(STAs(:,peakframe),[nstixperside, nstixperside, 3]);
maxes = squeeze(max(max(STA)));
mins = squeeze(min(min(STA)));
potentialnormfactors = [(1-[.5; .5; .5]-eps)./(maxes-[.5; .5; .5]); (-[.5; .5; .5]+eps)./(mins-[.5; .5; .5])];
potentialnormfactors(potentialnormfactors < 0) = []; % if min > mu or max < mu
normfactor = min(potentialnormfactors);
mumat = .5*ones(nstixperside, nstixperside, 3);
STAim = normfactor*(STA-mumat)+mumat;

% Getting gabors and contrast response functions
thetaidx = find(strcmp(GBstro.sum.trialFields(1,:),'g_theta'));
lambdaidx = find(strcmp(GBstro.sum.trialFields(1,:),'g_lambda'));
phiidx = find(strcmp(GBstro.sum.trialFields(1,:),'g_phi'));
sigmaidx = find(strcmp(GBstro.sum.trialFields(1,:),'g_sigma'));
gammaidx = find(strcmp(GBstro.sum.trialFields(1,:),'g_gamma'));
xoffsetidx = find(strcmp(GBstro.sum.trialFields(1,:),'g_xoff'));
yoffsetidx = find(strcmp(GBstro.sum.trialFields(1,:),'g_yoff'));
gaborcontrastidx = find(strcmp(GBstro.sum.trialFields(1,:),'g_cont'));
gabortypeidx = find(strcmp(GBstro.sum.trialFields(1,:),'g_colortype'));
gaborrgbidxs = [find(strcmp(GBstro.sum.trialFields(1,:),'g_r'));...
                find(strcmp(GBstro.sum.trialFields(1,:),'g_g'));...
                find(strcmp(GBstro.sum.trialFields(1,:),'g_b'))];

theta = unique(GBstro.trial(:,thetaidx));
lambda = unique(GBstro.trial(:,lambdaidx));
phi = unique(GBstro.trial(:,phiidx));
sigma = unique(GBstro.trial(:,sigmaidx));
gamma = unique(GBstro.trial(:,gammaidx));
xoffset = unique(GBstro.trial(:,xoffsetidx));
yoffset = unique(GBstro.trial(:,yoffsetidx));
colortype = GBstro.trial(:,gabortypeidx);
contrast = GBstro.trial(:,gaborcontrastidx);
gausslim = GBstro.sum.exptParams.gausslim/1000;
pixperdeg = GBstro.sum.exptParams.pixperdeg;

L = colortype == 4 & contrast == max(contrast);  % Achromatic
gaborrgb1 = unique(GBstro.trial(L,gaborrgbidxs),'rows');
L = colortype == 1 & contrast == max(contrast);  % Isoluminant
gaborrgb2 = unique(GBstro.trial(L,gaborrgbidxs),'rows');

%figure;
%subplot(2,1,1);
Achim = DrawGaborEdge([.5 .5 .5], gaborrgb1*2, [0 0 0], theta, lambda, sigma, gamma, phi, xoffset, yoffset, 0, 0, gausslim, pixperdeg);
%image(Achim)
%subplot(2,1,2);
Isoim = DrawGaborEdge([.5 .5 .5], gaborrgb2*2, [0 0 0], theta, lambda, sigma, gamma, phi, xoffset, yoffset, 0, 0, gausslim, pixperdeg);
%image(Isoim)
spikerates = [];
% Contrast response function
stimon_t = GBstro.trial(:,find(strcmp(GBstro.sum.trialFields(1,:),'stim_on')));
numframes = GBstro.sum.exptParams.nframes;
framerate = GBstro.sum.exptParams.framerate;
stimoff_t = stimon_t+numframes/framerate;
spikerates = [];
for i = 1:size(GBstro.trial,1)
    spiketimes = GBstro.ras{i,find(strcmp(GBstro.sum.rasterCells(1,:),spikename))};
    nspikes = sum(spiketimes > stimon_t(i) & spiketimes < stimoff_t(i));
    spikerates = [spikerates; nspikes./(stimoff_t(i)-stimon_t(i))];
end

figure; axeshandle = axes; hold on;
% Gabor contrast
cols = {'c',[],[],'k'};
xlabel('Gabor rms contrast'); ylabel('spikes/sec');
for j = [1,4]
    data = [];
    for i = unique(contrast)'
        L = (colortype == j) & (contrast == i);
        tmp = spikerates(L);
        data = [data; mean(tmp), sqrt(var(tmp)./length(tmp))];
    end
    if (size(data,1) > 1)
        errorbar(unique(contrast),data(:,1),data(:,2),cols{j},'linewidth',2);
    end
end
axes('units','inches','position',[1 3.5 .5 .5]);
image(STAim);
set(gca,'XTick',[],'YTick',[]);
axes('units','inches','position',[1.75 3.5 .5 .5]);
image(Achim);
set(gca,'XTick',[],'YTick',[]);
axes('units','inches','position',[2.5 3.5 .5 .5]);
image(Isoim);
set(gca,'XTick',[],'YTick',[]);


%%
% Adding a single extra point for the dark blob at 10% contrast
% for cell K11190800X.
GBstro = nex2stro('C:\NexFiles\Kali\K111908007.nex');
spikename = 'sig001a';

colortype = GBstro.trial(:,gabortypeidx);
contrast = GBstro.trial(:,gaborcontrastidx);
L = colortype == 4 & contrast == min(contrast);  % Achromatic
stimon_t = GBstro.trial(:,find(strcmp(GBstro.sum.trialFields(1,:),'stim_on')));
stimoff_t = stimon_t+numframes/framerate;
spikerates = [];
for i = 1:size(GBstro.trial,1)
    spiketimes = GBstro.ras{i,find(strcmp(GBstro.sum.rasterCells(1,:),spikename))};
    nspikes = sum(spiketimes > stimon_t(i) & spiketimes < stimoff_t(i));
    spikerates = [spikerates; nspikes./(stimoff_t(i)-stimon_t(i))];
end
mn = mean(spikerates(L));
sem = sqrt(var(spikerates(L))/sum(L));
axes(axeshandle);
errorbar(-min(contrast),mn,sem,'Color',[.5 .5 .5],'linewidth',2);
a = get(gcf,'Children');
set(gcf,'Children',flipud(a))
%%
% OK, here's where I'm really getting into the modeling stuff.
% Taking a single cell, fitting a parametric model to the firing rates
% and simulating many cells from this example one.
stro = nex2stro('C:\NexFiles\Kali\2008\K111308003.nex');
spikename = 'sig001b';

%stro = nex2stro('C:\NexFiles\Kali\K110708004.nex');
%spikename = getSpikenum(stro);

stimonidx = find(strcmp(stro.sum.trialFields(1,:),'stim_on'));
poststimtime = 0.05; % seconds after stimoff to include
gaborcontrastidx = find(strcmp(stro.sum.trialFields(1,:),'g_cont'));
gabortypeidx = find(strcmp(stro.sum.trialFields(1,:),'g_colortype'));
spikeidx = find(strcmp(stro.sum.rasterCells(1,:),spikename));
gaborcolortypes = stro.trial(:,gabortypeidx);
Lcolortype = gaborcolortypes == 5;  % 5 is L-M
Lcontrast = stro.trial(:,gaborcontrastidx) <= 1; %0.05;
% Cutting out all trials except the ones we are interested in 
stro.trial = stro.trial(Lcolortype & Lcontrast,:);
stro.ras = stro.ras(Lcolortype & Lcontrast,:);

prestimspikes = nan*ones(sum(Lcolortype & Lcontrast),1);
spikes = nan*ones(sum(Lcolortype & Lcontrast),1);
for i = 1:length(spikes)
    spikes(i) = sum(stro.ras{i,spikeidx} > stro.trial(i,stimonidx) & ...
        stro.ras{i,spikeidx} < stro.trial(i,stimonidx)+numframes./framerate+poststimtime);
    prestimspikes(i) = sum(stro.ras{i,spikeidx} < stro.trial(i,stimonidx) & ...
        stro.ras{i,spikeidx} > stro.trial(i,stimonidx)-numframes./framerate-poststimtime);
end
gaborcontrasts = stro.trial(:,gaborcontrastidx);
uniquegaborcontrasts = unique(gaborcontrasts);
means = [];
nonparametricrocs = [];
for i = 1:length(uniquegaborcontrasts)
    L = gaborcontrasts == uniquegaborcontrasts(i);
    ns(i) = sum(L);
    means(i) = mean(spikes(L));
    nonparametricrocs(i) = roc(prestimspikes,spikes(L));
end
meanprestim = mean(prestimspikes);
% Assuming that all the prestimspikerates are '0' which is almost true
% for the cell I'm considering.
% Let's see what kind of neurometric curve I get fron approximating the
% firing rate dist as normal...
%vars = 1.5*means;
%rocs = 1-normcdf(0,means,sqrt(vars));  % only works for cells with zero baseline
%rocs = normroc(meanprestim, sqrt(1.5*meanprestim), means, sqrt(vars));
%rocs(isnan(rocs)) = 0.5;
%plot(uniquegaborcontrasts, nonparametricrocs,'o-')

% % Testing out generating correlated random variables
% ncells = 50;
% nreps = 1000;
% r = .2;
% S2 = r*ones(ncells);
% S2 = S2+(1-r)*eye(ncells);
% Q = sqrtm(S2);
% Z = normrnd(0, 1, ncells, nreps);
% A = Q*Z;
% figure; 
% subplot(2,2,1); hold on;
% plot(mean(A'),var(A'),'k.');
% plot(mean(Z'),var(Z'),'b.');
% subplot(2,2,2);
% imagesc(corrcoef(A'));
% subplot(2,2,3); hold on;
% plot(A(1,:), A(2,:),'k.');
% subplot(2,2,4); hold on;
% plot(Z(1,:),Z(2,:),'b.');
% corrcoef(A(1,:),A(4,:))  % Should be close to r
% corrcoef(Z(1,:),Z(4,:))  % Should be close to 0
% % Seems to work OK.

%%
% Evaluate this to fit a parametric model to the contrast-response
% function.
% R = rmax*i^n/(i^n+c^n);
[modelparams, resid] = nlinfit(uniquegaborcontrasts', means, @nakarushton, [10*means(end), 1, max(uniquegaborcontrasts)])

figure; axes; hold on;
plot(uniquegaborcontrasts, means,'o');
x = linspace(min(uniquegaborcontrasts),max(uniquegaborcontrasts),100)
plot(x, nakarushton(modelparams,x),'-');
xlabel('rms cone contrast');
ylabel('spike count');
SStoaccountfor = sum((means-mean(means)).^2);
[means] = nakarushton(modelparams, uniquegaborcontrasts);
meanprestim = 0;
R2 = 1-sum(resid.^2)./SStoaccountfor
%%
% Create random responses from a simulated ensemble of neurons with the
% option of adding pooling noise
% n = 1 parametric rocs
% n = 1000, r = 0.2
% n = 1000, r = 0.1, pooling noise = 0

figure; axes; hold on
plot(uniquegaborcontrasts, nonparametricrocs,'+-')

ncells = 100;
nreps = 5000;
vartomean = 1.5;   % 1.5
poolingVMR = .6;  % .6
gainreductionfactor = 1.6; %1.6 % 
r = 0.2;
S2 = r*ones(ncells);
S2 = S2+(1-r)*eye(ncells);
Q = sqrtm(S2);
data = [];
for i = 0:length(uniquegaborcontrasts)
    Z = normrnd(0, 1, ncells, nreps);
    A = Q*Z;
    if (i == 0) % Doing the prestimspikes;
        A = A.*sqrt(vartomean*meanprestim)+meanprestim;        
    else
        A = A.*sqrt(vartomean*means(i)/gainreductionfactor)+means(i)/gainreductionfactor;
    end
    data = [data; mean(A,1)];
end
data = data+normrnd(0,sqrt(max(data,1)*poolingVMR));
%data(data<0) = 0;
correct = data(2:size(data,1),:)>repmat(data(1,:),size(data,1)-1,1);
%correct = correct+.5*(data(2:size(data,1),:)==repmat(data(1,:),size(data,1)-1,1));

rocs = sum(correct,2)./size(data,2);
%plot(uniquegaborcontrasts, rocs, 'ko');

[a,b,g,s] = weibullFit(uniquegaborcontrasts, rocs, 'sse',[mean(uniquegaborcontrasts), 1 1]);
[a,b,g,s] = weibullFit(uniquegaborcontrasts, [ns'.*rocs',ns'.*(1-rocs)'], 'mle', [a b]);
%plot(uniquegaborcontrasts,rocs,'o','MarkerSize',2,'MarkerFaceColor',[.4 .4 .4]);
%x = linspace(min(uniquegaborcontrasts),max(uniquegaborcontrasts),100);
%plot(x, 1-0.5.*exp(-((x./a).^b)),'-','Linewidth',1,'color','black');

% For plotting on the axes with the real data
plot(uniquegaborcontrasts,rocs,'o','MarkerSize',4,'MarkerFaceColor',[.7 .7 .7],'MarkerEdgeColor',[.7 .7 .7]);
plot(x, 1-0.5.*exp(-((x./a).^b)),'-','Linewidth',2,'color',[.7 .7 .7]);

%%
% NM2S figure
% Two panels showing performance with two polarities of gradient.

datapath = 'C:\PlexonData\Kali Psychophysics';
% S+ gradients
Sincfilenames = {'KP122408002.nex','KP122408005.nex'};  % Good
% S- gradients
Sdecfilenames = {'KP122208004.nex','KP122108003.nex'};  % Good

contrastlim = 0.05;
colors = [.8 .8 0; .5 .5 .5; 0 0 1];
figure;
set(gcf,'position',[463 170 308 543]);
for panel = 1:2
    subplot(2,1,panel); hold on;
    if (panel == 1)
        filenames = Sincfilenames;
        whichgradpolarities = [0 1];
    else
        filenames = Sdecfilenames;
        whichgradpolarities = [-1 0];
    end
    stro ={};
    for fileidx = 1:length(filenames)
        tmpstro = nex2stro([datapath,'/',filenames{fileidx}]);
        if (isempty(stro))
            stro = tmpstro;
        else
            stro = strocat(stro,tmpstro)
        end
    end
    
    ntrials = size(stro.trial,1);
    correctidx = find(strcmp(stro.sum.trialFields(1,:),'correct'));
    nonmatchsideidx = find(strcmp(stro.sum.trialFields(1,:),'nonmatchside'));
    correct = stro.trial(:,correctidx);
    nonmatchside = stro.trial(:,nonmatchsideidx);

    T1ccidxs = [find(strcmp(stro.sum.trialFields(1,:),'T1_Lcc'));...
        find(strcmp(stro.sum.trialFields(1,:),'T1_Mcc'));...
        find(strcmp(stro.sum.trialFields(1,:),'T1_Scc'))];
    T2ccidxs = [find(strcmp(stro.sum.trialFields(1,:),'T2_Lcc'));...
        find(strcmp(stro.sum.trialFields(1,:),'T2_Mcc'));...
        find(strcmp(stro.sum.trialFields(1,:),'T2_Scc'))];
    Grad1ccidxs = [find(strcmp(stro.sum.trialFields(1,:),'Grad1_Lcc'));...
        find(strcmp(stro.sum.trialFields(1,:),'Grad1_Mcc'));...
        find(strcmp(stro.sum.trialFields(1,:),'Grad1_Scc'))];
    Grad2ccidxs = [find(strcmp(stro.sum.trialFields(1,:),'Grad2_Lcc'));...
        find(strcmp(stro.sum.trialFields(1,:),'Grad2_Mcc'));...
        find(strcmp(stro.sum.trialFields(1,:),'Grad2_Scc'))];

    T1ccs = stro.trial(:,T1ccidxs);
    T2ccs = stro.trial(:,T2ccidxs);
    nonmatchccs = T2ccs;
    nonmatchccs(find(nonmatchside),:) = T1ccs(find(nonmatchside),:);
    matchccs = T1ccs;
    matchccs(find(nonmatchside),:) = T2ccs(find(nonmatchside),:);
    samplecc = unique(matchccs,'rows');
    deltaccs = nonmatchccs-repmat(samplecc,ntrials,1);

    gradccs = stro.trial(:,Grad2ccidxs);
    gradccs(find(nonmatchside),:) = stro.trial(find(nonmatchside),Grad1ccidxs);
    gradpolarity = sign(gradccs(:,3));

    contrasts = deltaccs(:,3);  % This is a difference between cone contrasts (relative to standard background)
    uniquecontrasts = unique(contrasts);
    x = [];
    n = [];
    for j = 1:length(whichgradpolarities)
        for i = 1:length(uniquecontrasts)
            L = gradpolarity == whichgradpolarities(j);
            L = L & contrasts == uniquecontrasts(i);
            n(:,i,j) = sum(L);
            x(:,i,j) = sum(correct(L));
        end
    end
    for j = 1:length(whichgradpolarities)
        col = colors(find(whichgradpolarities(j) == [-1 0 1]),:);
        p = squeeze(x(:,:,j)./n(:,:,j));
        se = sqrt((p.*(1-p))./n(:,:,j));
        errorbar(uniquecontrasts, p, se,'Color',col,'Linewidth',2);
    end
    set(gca,'YTick',[.4:.1:1],'XLim',[min(uniquecontrasts) max(uniquecontrasts)],'Ylim',[.4 1.01])
    plot([min(uniquecontrasts) max(uniquecontrasts)],[.5 .5],'k:')
    set(gca,'XTick',[-.05:.025:.05]);
    xtick = get(gca,'XTick');
    set(gca,'XTickLabel',xtick*100);
    ytick = get(gca,'YTick');
    set(gca,'YTickLabel',ytick*100);
    ylabel('% Correct');
    if (panel == 2)
        xlabel('S-cone contrast difference');
    end
end

%%
% Neurometric/psychometric comparison
% Using data from DTspot

whichcolordir = 2;
filename = 'C:\NexFiles\Sedna\2008\S122808011.nex';
stro=nex2stro(filename);
correct = stro.trial(:,find(strcmp(stro.sum.trialFields(1,:),'correct')));
colordir = stro.trial(:,find(strcmp(stro.sum.trialFields(1,:),'color_dir')));
contrast = stro.trial(:,find(strcmp(stro.sum.trialFields(1,:),'cntrst_lev')));
stimon_t = stro.trial(:,find(strcmp(stro.sum.trialFields(1,:),'flash_on')));
nframes = stro.trial(:,find(strcmp(stro.sum.trialFields(1,:),'numframes')));
stimoff_t = stimon_t+nframes/75;
flashside = sign(stro.trial(:,8))== -1;  % 1 = T1, 0 = T2
colordirs = reshape(stro.sum.exptParams.RF_colors,[3,3])/100;
framerate = stro.sum.exptParams.frame_rate; 
M = reshape(stro.sum.exptParams.m_mtx,[3 3]);
xx = linspace(0, 255, 2^16); %the desired quantization of the gammaTable
g1 = reshape(stro.sum.exptParams.gamma_table, 256, 3);
gammaTable = [spline([0:255], g1(:,1), xx)', spline([0:255], g1(:,2), xx)', spline([0:255], g1(:,3), xx)'];
bkgndrgb = [stro.sum.exptParams.bkgnd_r stro.sum.exptParams.bkgnd_g stro.sum.exptParams.bkgnd_b];
flash_R = stro.trial(:,find(strcmp(stro.sum.trialFields(1,:),'flash_R')));
flash_G = stro.trial(:,find(strcmp(stro.sum.trialFields(1,:),'flash_G')));
flash_B = stro.trial(:,find(strcmp(stro.sum.trialFields(1,:),'flash_B')));
pred_rgb = [gammaTable(flash_R+1, 1), gammaTable(flash_G+1, 2), gammaTable(flash_B+1, 3)];
pred_rgb = pred_rgb - repmat(bkgndrgb, size(pred_rgb,1), 1);
bkgndLMS = M*bkgndrgb';
pred_lms = [(M * pred_rgb') ./ repmat(bkgndLMS, 1, size(pred_rgb,1))]';

uniquecontrasts = unique(contrast)';
psychdata = [];
for i = 2:max(uniquecontrasts)
    L = colordir == whichcolordir;
    L = L & contrast  == i;
    cc = unique(pred_lms(L,:),'rows');
    psychdata = [psychdata; sqrt(mean(cc.^2)) sum(correct(L)), sum(L)];
end

[aSSE, bSSE, gSSE] = weibullFit(psychdata(:,1),psychdata(:,2)./psychdata(:,3) , 'sse' ,[median(psychdata(:,1)) 1]);
[aMLE, bMLE, gMLE, success] = weibullFit(psychdata(:,1),[psychdata(:,2),psychdata(:,3)-psychdata(:,2)], 'mle', [aSSE, bSSE]);
x = logspace(log10(min(psychdata(:,1))), log10(max(psychdata(:,1))),100);
ypsych = 1 + (0.5 - 1).*exp(-(x./aMLE).^bMLE);

% Now doing neurometric fn.
% First just looking at spike counts
% later doing something more exotic
spikecounts = nan*ones(length(stimon_t),1);
for i = 1:length(stimon_t)
   spikecounts(i) = sum(stro.ras{i} > stimon_t(i)+.3 & stro.ras{i} < stimoff_t(i));
end
%nulldistn = spikecounts(flashside == 0);
nulldistn = spikecounts(contrast == 1);
neurodata = [];
spikerates = [];
for i = 2:max(uniquecontrasts)
    L = colordir == whichcolordir;
    L = L & contrast  == i;
    cc = unique(pred_lms(L,:),'rows');
    T1distn = spikecounts(L & flashside == 1);
    rocA = roc(nulldistn, T1distn);
    neurodata = [neurodata; sqrt(mean(cc.^2)), rocA*sum(L), sum(L)];
    spikerates = [spikerates; mean(T1distn)];
end
[aSSE, bSSE, gSSE] = weibullFit(neurodata(:,1),neurodata(:,2)./neurodata(:,3) , 'sse' ,[median(neurodata(:,1)) 1]);
[aMLE, bMLE, gMLE, success] = weibullFit(neurodata(:,1),[neurodata(:,2),neurodata(:,3)-neurodata(:,2)], 'mle', [aSSE, bSSE]);
yneuro = 1 + (0.5 - 1).*exp(-(x./aMLE).^bMLE);

figure;
axes; hold on;
plot(psychdata(:,1),psychdata(:,2)./psychdata(:,3),'k.','MarkerSize',10);
plot(x,ypsych,'k-','LineWidth',2);
plot(neurodata(:,1),neurodata(:,2)./neurodata(:,3),'b.','MarkerSize',10);
plot(x,yneuro,'b-','LineWidth',2);

set(gca,'XScale','log');
set(gca,'XLim',[psychdata(1,1)/1.1 psychdata(end,1)*1.1]);
xlabel('RMS Cone Contrast');
ylabel('% Correct');

% Rasters

for i = 0:1
    figure; axes; hold on;
    counter = 0;
    for j = 1:max(uniquecontrasts)
        L = colordir == whichcolordir;
        L = L & contrast  == j;
        L = L & flashside == i;
        whichtrials = find(L);
        spikes = stro.ras(L);
        for k = 1:length(spikes)
            nspikestot = length(spikes{k});
            plot([spikes{k} spikes{k}]'-stimon_t(whichtrials(k)),[zeros(nspikestot,1) .5*ones(nspikestot,1)]'+counter,'k-');
            counter = counter+1;
        end
        plot([0 1],[counter counter],'k:');
    end
    plot([stimoff_t(1)-stimon_t(1) stimoff_t(1)-stimon_t(1)],[0 counter],'k-');
    set(gca,'Xlim',[-.5 1]);
end
%%
% Now some lame modeling simulate a pool of neurons.
% In progress...
[modelparams, resid] = nlinfit(neurodata(:,1), spikerates, @nakarushton, [10*spikerates(end), 1, max(uniquecontrasts)])

figure; axes; hold on;
plot(neurodata(:,1), spikerates,'o');
x = linspace(min(neurodata(:,1)),max(neurodata(:,1)),100)
plot(x, nakarushton(modelparams,x),'-');
xlabel('rms cone contrast');
ylabel('spike count');
SStoaccountfor = sum((spikerates-mean(spikerates)).^2);
[means] = nakarushton(modelparams, neurodata(:,1));
meanprestim = 0;
R2 = 1-sum(resid.^2)./SStoaccountfor
means = means+mean(mean(nulldistn));  % Only do this once!

% 20 independent cells or 1000 cells with r = 0.1!
ncells = 1000;  % 20
nreps = 500;
vartomean = 1.5;   % 1.5
poolingVMR = .6;  % .6
gainreductionfactor = 1.6; %1.6 % 
r = 0.1; % 0
S2 = r*ones(ncells);
S2 = S2+(1-r)*eye(ncells);
Q = sqrtm(S2);
data = [];
for i = 0:size(neurodata,1)
    Z = normrnd(0, 1, ncells, nreps);
    A = Q*Z;
    if (i == 0) % Doing the prestimspikes;
        A = A.*sqrt(vartomean*meanprestim)+meanprestim;        
    else
        A = A.*sqrt(vartomean*means(i)/gainreductionfactor)+means(i)/gainreductionfactor;
    end
    data = [data; mean(A,1)];
end
data = data+normrnd(0,sqrt(max(data,1)*poolingVMR));
%data(data<0) = 0;
correct = data(2:size(data,1),:)>repmat(data(1,:),size(data,1)-1,1);
%correct = correct+.5*(data(2:size(data,1),:)==repmat(data(1,:),size(data,1)-1,1));

rocs = sum(correct,2)./size(data,2);
plot(neurodata(:,1), rocs, 'ko');

[a,b,g,s] = weibullFit(neurodata(:,1), rocs, 'sse',[mean(neurodata(:,1)), 1 1]);
[a,b,g,s] = weibullFit(neurodata(:,1), [neurodata(:,3).*rocs,neurodata(:,3).*(1-rocs)], 'mle', [a b]);
plot(neurodata(:,1),rocs,'o','MarkerSize',2,'MarkerFaceColor',[.4 .4 .4]);
x = linspace(min(neurodata(:,1)),max(neurodata(:,1)),100);
plot(x, 1-0.5.*exp(-((x./a).^b)),'-','Linewidth',1,'color','black');

% For plotting on the axes with the real data
plot(neurodata(:,1,rocs),'o','MarkerSize',4,'MarkerFaceColor',[.7 .7 .7],'MarkerEdgeColor',[.7 .7 .7]);
plot(x, 1-0.5.*exp(-((x./a).^b)),'-','Linewidth',2,'color',[.7 .7 .7]);
%%
% Gabor fit for STA.
% This stuff also derives estimates of preferred orientation, spatial
% frequency, and spatial phase from the fitted gabor that can be compared
% to measurements made by the grating paradigm.

stro=nex2stro(findfile('K122908001.nex'));
whichframe = 4;
framerate = stro.sum.exptParams.framerate;
nstixperside = stro.sum.exptParams.nstixperside;
ntrials = length(stro.sum.absTrialNum);
stimonidx = find(strcmp(stro.sum.trialFields(1,:),'stim_on'));
stimoffidx = find(strcmp(stro.sum.trialFields(1,:),'all_off'));
nframesidx = find(strcmp(stro.sum.trialFields(1,:),'num_frames'));
noisetypeidx = find(strcmp(stro.sum.trialFields(1,:),'noise_type'));
sigmaidxs = strmatch('sigma',stro.sum.trialFields(1,:));

out = getWhtnsStats(stro,whichframe,'STCOVmex', {nstixperside^2, 3, whichframe}, getSpikenum(stro));
tmpstro = [];
STAs = out{1};
STCs = out{2};
nspikes = out{3};

% fitting the gabor
im = reshape(STAs(:,end),nstixperside*nstixperside,3);
[u,s,v] = svd(im');
if (max(abs(v(:,1))) ~= max(v(:,1)))
    v = -v;
    u = -u;
end
im = reshape(v(:,1),[nstixperside nstixperside]);
struct = gaborfit(im);

pixperstix = 10;
npix = pixperstix*nstixperside;
margin = 1-(1/pixperstix);
interval = linspace(1-margin,nstixperside+margin,npix)-ceil(median(1:nstixperside));
[X, Y] = meshgrid(interval,interval);
X = X-struct.xoffset;
Y = Y+struct.yoffset;
xprime = X.*cos(-struct.theta)+Y.*sin(-struct.theta);
yprime = -X.*sin(-struct.theta)+Y.*cos(-struct.theta);
fittedgabor = exp(-(xprime.^2+struct.gamma.^2.*yprime.^2)./(2.*struct.sigma.^2)).*cos(2.*pi.*yprime./struct.lambda-struct.phi);

colorgabor = u(:,1)*reshape(fittedgabor,1,npix^2)/3;
colorgabor = colorgabor+repmat([.5 .5 .5]',1,npix^2);
colorgabor = reshape(colorgabor,[3, npix, npix]);
colorgabor = permute(colorgabor,[2 3 1]);
image(colorgabor);
axis image;
set(gca,'XTick',[],'YTick',[]);

% Here are the fitted orientation and spatial frequency from the gabor
disp('orientation estimate is (deg): ');
mod([struct.theta*180/pi struct.theta*180/pi+180],360)
% lambda is wavelength in stixels.
mondist = 100; % cm
screenwidth = 35; %cm
screenwidthinpix = 1024; % singlewide pixels
pixpercm = screenwidthinpix/screenwidth;
cmperdeg = screenwidth/(2*atan2(screenwidth/2, mondist)*180/pi);
pixperdeg = pixpercm*cmperdeg;
stixperdeg = pixperdeg/(2*stro.sum.exptParams.npixperstix); % "2*" because npixperstix is in doublewides.

%1/struct.lambda is in cycles/stixel
disp('sf estimate is (cpd): ');
(1/struct.lambda)*stixperdeg

% Trying to get a quantitive estimate of phase from grating data (protocol 3) 
% that we can compare with the phi parameter from the gabor fit.
% The phase calculation in the grating paradigm does not take the x and y
% offset into account (the gabor fit does).  Basically, the grating phase
% is off from the gabor phase by an amount due to xoffset and yoffset.
gratingstro=nex2stro('C:\NexFiles\Kali\K122908005.nex');

diam = gratingstro.trial(end,5); % degrees
orient = gratingstro.trial(end,11); % rad

pixpercycle = pixperdeg/((1/struct.lambda)*stixperdeg); % singlewide pixels (y)
sizeinpix = round(pixperdeg*diam/2)*2;
xinc = (2*pi/pixpercycle)*cos((pi/2)-orient);
yinc = (2*pi/pixpercycle)*sin((pi/2)-orient);
[xramp, yramp] = meshgrid(xinc*([0:2:sizeinpix-1]), yinc*([0:sizeinpix-1]));
a = cos(xramp+yramp+phase);
imagesc(a); colormap(gray); 
% the way grating works, the origin is at the upper left corner

% Finding phase in the center of the display
centerx = (sizeinpix/2)-6*struct.xoffset; % 6 = pixels per stixel
centery = (sizeinpix/2)+6*struct.yoffset;
hold on; plot(centerx/2, centery,'y*');
angletorf = atan2(-centery, centerx);
distancetorf = sqrt(centerx^2+centery^2);

ncycles = (1/pixpercycle)*sin(mod(orient,pi)-angletorf)*distancetorf % how many cycles are between the corner and the center
phasedelta = rem(ncycles,1)*2*pi
disp('preferred phase estimate is (rad): ');
struct.phi-phasedelta;

%%
% Figure showing STAs, fitted Gabors, and positioned squares for NM2S
% paradigm

filenames = {'K032608004', 'K032008005'}
configurations = [2 1];
spikenames = {'sig001a','sig001a'};
%stro = nex2stro(findfile('K090108006'));
%stro = nex2stro(findfile('K110608001'));
%stro = nex2stro(findfile('K032608004'));  % Good for config 2
%stro = nex2stro(findfile('K050508003'));
%stro = nex2stro(findfile('K110408001'));
%stro = nex2stro(findfile('K032008005'));
%stro = nex2stro(findfile('K032608003'));
%stro = nex2stro(findfile('K111908002'));  % Lowpass cell, lambda is huge
%stro = nex2stro(findfile('K090508003'));
%stro = nex2stro(findfile('K090408002'));
%stro = nex2stro(findfile('K091008001'));
%stro = nex2stro(findfile('K091608005'));
%stro = nex2stro(findfile('K010809002'));
%stro = nex2stro(findfile('K052308001'));
%stro = nex2stro(findfile('K053008003'));
%stro = nex2stro(findfile('K041008001'));

% constants
maxT = 8;
mondist = 100; % cm
screenwidth = 36; %cm
screenwidthinpix = 1024; % singlewide pixels
pixpercm = screenwidthinpix/screenwidth;
cmperdeg = screenwidth/(2*atan2(screenwidth/2, mondist)*180/pi);
pixperdeg = pixpercm*cmperdeg;
gl_stimsize = 1.2;  % 1 degree on a side square

for i = 1:length(filenames)
    stro = nex2stro(findfile(filenames{i}));

    framerate = stro.sum.exptParams.framerate;
    nstixperside = stro.sum.exptParams.nstixperside;
    ntrials = length(stro.sum.absTrialNum);
    stimonidx = find(strcmp(stro.sum.trialFields(1,:),'stim_on'));
    stimoffidx = find(strcmp(stro.sum.trialFields(1,:),'all_off'));
    nframesidx = find(strcmp(stro.sum.trialFields(1,:),'num_frames'));
    noisetypeidx = find(strcmp(stro.sum.trialFields(1,:),'noise_type'));
    sigmaidxs = strmatch('sigma',stro.sum.trialFields(1,:));
    stixperdeg = pixperdeg/(2*stro.sum.exptParams.npixperstix);

    L = stro.trial(:,noisetypeidx) == 1;  % Getting rid of cone noise
    stro.ras(~L,:) = [];
    stro.trial(~L,:) = [];
    spikename = spikenames{i};

    out = getWhtnsStats(stro,maxT,'STCOVmex', {nstixperside^2, 3, maxT}, spikename);
    STAs = out{1};
    STCs = out{2};
    nspikes = out{3};

    energy = sum(STAs.^2);
    whichframe = find(energy == max(energy));

    % fitting the gabor
    STA = reshape(STAs(:,whichframe),nstixperside*nstixperside,3);
    [u,s,v] = svd(STA');
    im = reshape(v(:,1),[nstixperside nstixperside]);
    struct = gaborfit(im);  % GH: having trouble with the way PTB and Matlab set up the y-axis (increasing or decreasing)
    % Converting struct fields to dva
    struct.lambda = struct.lambda/stixperdeg;
    struct.sigma = struct.sigma/stixperdeg;
    struct.xoffset = struct.xoffset/stixperdeg;
    struct.yoffset = struct.yoffset/stixperdeg;
 
    STAim = reshape(STA./(3*max(abs(STA(:))))+.5,[nstixperside nstixperside 3])
    im = DrawGaborEdge([.5 .5 .5], u(:,1)/4, [0 0 0], struct.theta, struct.lambda, struct.sigma, struct.gamma, struct.phi, struct.xoffset, struct.yoffset, 0, 0, .995, pixperdeg);

    figure;
    subplot(2,1,1);
    image(STAim);
    set(gca,'XTick',[],'YTick',[]); axis square;
    subplot(2,1,2);
    image(im);
    hold on; set(gca,'XTick',[],'YTick',[]);
    axis square; %axis xy;

    phasealignshift = ((-2/pi)*mod(struct.phi,pi)+1)*struct.lambda/4;
    xlims = get(gca,'Xlim');
    ylims = get(gca,'Ylim');
    stim.x = median(xlims)/pixperdeg; % in dva
    stim.y = median(ylims)/pixperdeg; % in dva

    if (configurations(i) == 1)
        x = stim.x+struct.xoffset;
        y = stim.y-struct.yoffset;
    elseif (configurations(i) == 2)
        x1 = sin(struct.theta)*(gl_stimsize/2-phasealignshift);
        y1 = -cos(struct.theta)*(gl_stimsize/2-phasealignshift);

        x2 = -sin(struct.theta)*(gl_stimsize/2+phasealignshift);
        y2 = cos(struct.theta)*(gl_stimsize/2+phasealignshift);

        if (x1^2+y1^2 < x2^2+y2^2)  % Move center by the short distance
            x = stim.x+struct.xoffset+x1;
            y = stim.y-(struct.yoffset+y1);  % '-' to account for difference in plotting convention
        else
            x = stim.x+struct.xoffset+x2;
            y = stim.y-(struct.yoffset+y2);
        end
    end
    stimcenterinpix = [x y]*pixperdeg;
    % Picking points on the square
    rotmat = [cos(-struct.theta) -sin(-struct.theta); sin(-struct.theta) cos(-struct.theta) ];
    % points = [linspace(-.5,.5,100)' -.5*ones(100,1);
    %     .5*ones(100,1) linspace(-.5,.5,100)';
    %     linspace(-.5,.5,100)' .5*ones(100,1);
    %     -.5*ones(100,1) linspace(-.5,.5,100)']*gl_stimsize;
    points = [-.5 -.5; .5 -.5; .5 .5; -.5 .5; -.5 -.5]*gl_stimsize;
    if (configurations(i) == 1)
        transformedpoints = points'.*pixperdeg+repmat(stimcenterinpix',1,size(points,1));
    elseif (configurations(i) == 2)
        transformedpoints = (rotmat*points').*pixperdeg+repmat(stimcenterinpix',1,size(points,1));
    end
    plot(transformedpoints(1,:),transformedpoints(2,:),'k-','Linewidth',2)
end

