stro = nex2stro;
ntrials = size(stro.trial,1);
ntrlspercond = stro.sum.exptParams.trialspercond;
fpacqidx = find(strcmp(stro.sum.trialFields(1,:),'fp_acq'));
stimonidx = find(strcmp(stro.sum.trialFields(1,:),'stim_on'));
fpoffidx = find(strcmp(stro.sum.trialFields(1,:),'fp_off'));
saccadeidx = find(strcmp(stro.sum.trialFields(1,:),'saccade'));

fpacq_t = stro.trial(:,fpacqidx);
fpoff_t = stro.trial(:,fpoffidx);
stimon_t = stro.trial(:,stimonidx);
framerate = stro.sum.exptParams.framerate;
stimoff_t = stimon_t+stro.sum.exptParams.nframes_discriminanda/framerate;

saccade_t = stro.trial(:,saccadeidx);

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

spikeidx = find(strcmp(stro.sum.rasterCells(1,:),getSpikenum(stro)));
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
%%
% Basic performance
colors = [.8 .8 0; .5 .5 .5; 0 0 1];
gradpolarity = sign(gradccs(:,3));
contrasts = deltaccs(:,3);  % This is a difference between cone contrasts (relative to standard background)
uniquecontrasts = unique(contrasts);
uniquegrads = unique(gradpolarity);
figure; axes; hold on;
x = [];
n = [];
for j = 1:length(uniquegrads)
    for i = 1:length(uniquecontrasts)
        L = gradpolarity == uniquegrads(j);
        L = L & contrasts == uniquecontrasts(i);
        n(:,i,j) = sum(L);
        x(:,i,j) = sum(correct(L));
    end
end
for j = 1:length(uniquegrads)
    col = colors(find(uniquegrads(j) == [-1 0 1]),:);
    p = squeeze(x(:,:,j)./n(:,:,j));
    se = sqrt((p.*(1-p))./n(:,:,j));
    errorbar(uniquecontrasts, p, se,'Color',col,'Linewidth',2);
end
ps =[];
for i = 1:length(uniquecontrasts)
    tmpx = squeeze(x(:,i,:));
    tmpn = squeeze(n(:,i,:));
    [h,p] = equalproptest(tmpx,tmpn,0.05);
    ps = [ps; p];
    if (h)
        plot(uniquecontrasts(i),1,'k*');
    end
end
% Overall performance
[h,p] = equalproptest(sum(x),sum(n),0.05);
if (h)
    disp('Alert!  Gradient influences overall percent correct!');
end
disp(['Test on overall performance: ',num2str(p)]);
disp('Contrast by contrast % correct comparisons');
disp(ps')
title(stro.sum.fileName);
%%
% performance as a function of contrast and side
contrasts = deltaccs(:,3);
uniquecontrasts = unique(contrasts);
gradpolarity = sign(gradccs(:,3));
uniquegrads = unique(gradpolarity);

x = [];
n = [];
for k = 1:length(uniquegrads)
    for i = 1:length(uniquecontrasts)
        for j = 0:1 % 0 means T2 is correct, 1 means T1 is correct
            L = contrasts == uniquecontrasts(i);
            L = L & nonmatchside == j;
            L = L & gradpolarity == uniquegrads(k);
            n(i,j+1,k) = sum(L);
            x(i,j+1,k) = sum(correct(L));
        end
    end
end

figure; 
for k = 1:length(uniquegrads)
    subplot(length(uniquegrads),2,2*(k-1)+1); hold on;
    for i = 1:size(x,1)
        h = plot([1 0], x(i,:,k)./n(i,:,k),'k.-');
        col = [.5 .5 .5] + uniquecontrasts(i)*[-1 -1 1];
        set(h,'Linewidth',3,'MarkerSize',10,'Color',col);
    end
    set(gca,'Xtick',[0 1],'Xticklabel',{'T1','T2'});
    ylabel('Prop. correct');

    subplot(length(uniquegrads),2,2*(k-1)+2); hold on;
    plot(uniquecontrasts, x(:,1,k),'k-'); 
    plot(uniquecontrasts, x(:,2,k),'b-');
    if (k == length(uniquegrads))
        legend('T2 correct','T1 correct');
    end
end
%%
% Looking for any unexpected effects of the gradient

gradpolarity = sign(gradccs(:,3));
contrasts = deltaccs(:,3);
uniquecontrasts = unique(contrasts);
uniquegrads = unique(gradpolarity);
lat = saccade_t-fpoff_t;
if (length(uniquegrads) > 1)
    figure; axes; hold on; xlabel('Contrast'); ylabel('Latency');
    for i = 1:length(uniquegrads)
        L = gradpolarity == uniquegrads(i);
        meanlat = mean(lat(L));
        disp(['Grad ',num2str(uniquegrads(i)),' % correct: ',num2str(100*sum(correct(L)/sum(L)))]);
        disp(['Grad ',num2str(uniquegrads(i)),' mean latency: ',num2str(meanlat)]);
        disp(['-------------']);
        data = [];
        for j = 1:length(uniquecontrasts)
            L = gradpolarity == uniquegrads(i);
            L = L & contrasts == uniquecontrasts(j);
            data = [data, [mean(lat(L)); sqrt(var(lat(L))/sum(L))]];
        end
        col = colors(find(uniquegrads(i) == [-1 0 1]),:);
        errorbar(uniquecontrasts, data(1,:),data(2,:),'Color',col,'Linewidth',2);
    end
end
       
%%
% Running mean of performance

lat = saccade_t-fpoff_t;

figure;
k = 20;
kernel = ones(1,k)./k;
x = conv(kernel, correct);
x(end-k:end) = [];
x(1:k-1) = [];
subplot(2,1,1);
plot(x);
x = conv(kernel, lat);
x(end-k:end) = [];
x(1:k-1) = [];
subplot(2,1,2);
plot(x);
r = corrcoef([correct, lat])
title(['r = ',num2str(r(1,2))]);

%%
% Looking at eye movements
EPx = stro.ras(:,find(strcmp(stro.sum.rasterCells, 'AD11')));
EPy = stro.ras(:,find(strcmp(stro.sum.rasterCells, 'AD12')));
EPstart_t = [stro.ras{:,find(strcmp(stro.sum.rasterCells, 'anlgStartTime'))}]';
rates = [stro.sum.analog.storeRates{:}];
POSTSACMARGIN = 0.05;
EPnormfactor = ((10/4096)*40)^-1;
%(10 V/2^12 AD values)*40 AD values per degree;

figure; axes; hold on; axis equal;
for i = 1:size(EPx,1)
    t = [0:length(EPx{i})-1]*1/rates(1)+EPstart_t(i);
    x = EPx{i};
    y = EPy{i};
    L = t > fpacq_t(i) & t < saccade_t(i)+POSTSACMARGIN;
    if (nonmatchside(i) & correct(i));
        col = [1 0 0];
    elseif (nonmatchside(i) & ~correct(i));
        col = [0 .5 0];
    elseif (~nonmatchside(i) & correct(i));
        col = [0 1 0];
    elseif (~nonmatchside(i) & ~correct(i));
        col = [.5 0 0];
    end
    plot(EPnormfactor*x(L),EPnormfactor*y(L),'Color',col);
   % pause;
end
t1x = stro.sum.exptParams.targdistprop*(stro.sum.exptParams.rf_x/10);
t1y = stro.sum.exptParams.targdistprop*(stro.sum.exptParams.rf_y/10);
t2x = -t1x;
t2y = -t1y;
plot(t1x, t1y,'ko','MarkerSize',10);
plot(t2x, t2y,'ko','MarkerSize',10);


%%
% Latency analysis
lat = saccade_t-fpoff_t;
figure; subplot(2,2,1);
[n,x] = hist(lat, 20);
bar(x,n);
set(gca,'XLim',[min(x), max(x)]);

subplot(2,2,2); hold on;
[n1,x] = hist(lat(correct == 1), x);
[n2,x] = hist(lat(correct == 0), x);
bar(x,n1,'g');
bar(x,n2,'r');
set(gca,'XLim',[min(x), max(x)]);

data = []; sem = [];
for i = 1:length(uniquecontrasts)
    for j = 0:1 % 0 means T2 is correct, 1 means T1 is correct
        for k = 0:1% incorrect/correct
            L = contrasts == uniquecontrasts(i);
            L = L & nonmatchside == j;
            Lcorrect = logical(correct == k);
            L = L & Lcorrect;
            data(i,j+1,k+1) = mean(lat(L));
            sem(i,j+1,k+1) = sqrt(var(lat(L))/sum(L));
        end
    end
end
subplot(2,2,3); hold on; % Correct trials
errorbar(uniquecontrasts, data(:,2,1), sem(:,2,1),'k.-'); % looking to the left, incorrectly
errorbar(uniquecontrasts, data(:,2,2), sem(:,2,2),'m.-'); % looking to the left, correctly
set(gca,'Ylim',[min(data(:)), max(data(:))]);

subplot(2,2,4); hold on; % Incorrect trials
errorbar(uniquecontrasts, data(:,1,1), sem(:,1,1),'k.-'); % looking to the right, incorrectly
errorbar(uniquecontrasts, data(:,1,2), sem(:,1,2),'m.-'); % looking to the right, correctly
set(gca,'Ylim',[min(data(:)), max(data(:))]);
legend('inc','cor');

%% 
% Fixation saccades
sacstats = getSacData(stro);
saccadeoccurred = nan*ones(ntrials,1);
for i = 1:ntrials
    st = sacstats.starttimes{i};
    Lsac = (st > stimon_t(i)) & (st < stimoff_t(i));
   % Lsac = (st < stimon_t(i)) & (st > stimon_t(i)-(stimoff_t(i)-stimon_t(i))); 
    if (any(Lsac))
        saccadeoccurred(i) = 1;
    else
        saccadeoccurred(i) = 0;
    end
end
corrcoef([saccadeoccurred, correct])

contrasts = deltaccs(:,3);
uniquecontrasts = unique(contrasts);
data = [];
for i = uniquecontrasts'
    L = contrasts == i;
    r = corrcoef([saccadeoccurred(L), correct(L)]);
    data = [data; r(1,2)];
end
figure; axes; hold on;
plot(uniquecontrasts,data,'ko-');
%%
% More fixational sacade stuff
%sacstats = getSacData(stro);
amplitudes = [];
directions = [];
peakv = [];
pathlengths = [];
durations = [];
BINWIDTH = .01;
timebins = [-.25:BINWIDTH:.9];
PSTHmat = [];
for i = 1:ntrials
    st = sacstats.starttimes{i};
    Lsac = (st > fpacq_t(i)) & (st < fpoff_t(i)-.1); 
    %Lsac = ones(length(st),1);
    if any(Lsac)
        if any(sacstats.amplitudes{i}(Lsac) > 50)
            keyboard
        end
        amplitudes = [amplitudes; sacstats.amplitudes{i}(Lsac)];
        directions = [directions; sacstats.directions{i}(Lsac)];
        peakv = [peakv; sacstats.peakv{i}(Lsac)];
        pathlengths = [pathlengths; sacstats.pathlengths{i}(Lsac)];
        durations = [durations; sacstats.durations{i}(Lsac)];
        
        sactimes = [];
        for j = find(Lsac')
            tmp = sacstats.starttimes{i}(j)-stimon_t(i);
            sactimes = [sactimes; tmp((tmp > timebins(1)-BINWIDTH/2) & (tmp < timebins(end)+BINWIDTH/2))];
        end
        [n,x] = hist(sactimes, timebins);
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
xlabel('time wrt stimulus onset (s)'); ylabel('saccades/sec');
stimdur = stro.sum.exptParams.nframes_discriminanda/stro.sum.exptParams.framerate;
plot([stimdur stimdur],[0 max(mean(PSTHmat))],'g:');

%%
% Looking at effects of microstimulation
correctidx = find(strcmp(stro.sum.trialFields(1,:),'correct'));
nonmatchsideidx = find(strcmp(stro.sum.trialFields(1,:),'nonmatchside'));
mstimidx = find(strcmp(stro.sum.trialFields(1,:),'mstim'));
correct = stro.trial(:,correctidx);
nonmatchside = stro.trial(:,nonmatchsideidx);
mstim = stro.trial(:,mstimidx);

contrasts = deltaccs(:,3);
uniquecontrasts = unique(contrasts);

x = [];
n = [];
for k = 0:1
    for i = 1:length(uniquecontrasts)
        for j = 0:1 % 0 means T2 is correct, 1 means T1 is correct
            L = contrasts == uniquecontrasts(i);
            L = L & nonmatchside == j;
            L = L & mstim == k;
            n(i,j+1,k+1) = sum(L);
            x(i,j+1,k+1) = sum(correct(L));
        end
    end
end

figure; axes; hold on;
plot(uniquecontrasts, x(:,1,1)./n(:,1,1),'ko-');
plot(uniquecontrasts, x(:,2,1)./n(:,2,1),'ro-');
plot(uniquecontrasts, x(:,1,2)./n(:,1,2),'ko:');
plot(uniquecontrasts, x(:,2,2)./n(:,2,2),'ro:');
legend('Color on right nostim','Color on left nostim','Color on right stim','Color on left stim')
xlabel('Delta contrast');
ylabel('% chosen');

%%
% Looking for a tuning curve to the stimulus in the RF
psthendbins = [-.1 .2];
countingendbins = [.04 .2];

nbins = 20;
bins = linspace(psthendbins(1),psthendbins(2),nbins);
correctidx = find(strcmp(stro.sum.trialFields(1,:),'correct'));
nonmatchsideidx = find(strcmp(stro.sum.trialFields(1,:),'nonmatchside'));
mstimidx = find(strcmp(stro.sum.trialFields(1,:),'mstim'));
correct = stro.trial(:,correctidx);
nonmatchside = stro.trial(:,nonmatchsideidx);
uniquecontrasts = unique(deltaccs,'rows');
grads = stro.trial(:,find(strcmp(stro.sum.trialFields(1,:),'Grad1_Lcc'))+[0 1 2]);
uniquegrads = unique(grads,'rows');

spikename = getSpikenum(stro);
spikeidx = find(strcmp(stro.sum.rasterCells(1,:),spikename));
spikes = stro.ras(:,spikeidx);
figure;
data = {};
titles = {'T1 is correct','T2 is correct'};
for k = 1:size(uniquegrads,1)
    for i = 1:size(uniquecontrasts,1)
        for j = 0:1 % 0 means T2 is correct, 1 means T1 is correct
            L = all(deltaccs == repmat(uniquecontrasts(i,:),[size(deltaccs,1) 1]),2);
            L = L & all(grads == repmat(uniquegrads(k,:),[size(grads,1) 1]),2);
            L = L & nonmatchside == j;
            trlidxs = find(L);
            if (isempty(trlidxs))
                continue;
            end
        
            subplot(3,2,(-j+2)+2*(k-1)); hold on; set(gca,'Color',[.7 0 .6] )
            title(titles(-j+2));
            psth = zeros(1,length(bins));
            spikecounts = [];
            for counter = 1:sum(L)
                trlidx = trlidxs(counter);
                tmpspikes = spikes{trlidx};
                Lpsth = tmpspikes > stimon_t(trlidx)+psthendbins(1) &...
                    tmpspikes < stimon_t(trlidx)+psthendbins(2);
                Lanalysis = tmpspikes > stimon_t(trlidx)+countingendbins(1) &...
                    tmpspikes < stimon_t(trlidx)+countingendbins(2);
            
                spikecounts = [spikecounts; sum(Lanalysis)];
                [n,x] = hist(tmpspikes(Lpsth)-stimon_t(trlidx),bins);
                psth = psth+n;
            end
            data(i,-j+2, k) = {spikecounts./(countingendbins(2)-countingendbins(1))};
            psth = psth/sum(L)/(x(2)-x(1));
            h = plot(x,psth);
            set(h,'Color',[i i i]/length(uniquecontrasts),'LineWidth',2);
            set(gca,'Xlim',psthendbins);
        end
    end
end
ymax = 0;
naxes = length(get(gcf,'Children'));
for i = 1:naxes
    subplot(3,2,i);
    tmp = get(gca,'Ylim');
    ymax = max([ymax tmp(2)]);
end
for i = 1:naxes
    subplot(3,2,i)
    set(gca,'Ylim',[0 ymax]);
end

% First column of data is colored square at T1
% Second column of data is colored square at T2

means = [];
ses = [];
for k = 1:size(uniquegrads,1)
    for i = 1:length(uniquecontrasts)
        for j = 1:size(data,2)
            means(i,j,k) = mean(data{i,j,k});
            ses(i,j,k) = sqrt(var(data{i,j,k})/length(data{i,j,k}));
        end
    end
end
idx = find(var(uniquecontrasts) == max(var(uniquecontrasts)),1);
cols = [.5 .5 0; .5 .5 .5; .5 .5 1];
linestyles = ['-',':'];
figure; hold on;
for k = 1:size(data,3)
    for j = 1:size(data,2)
        errorbar(uniquecontrasts(:,idx), means(:,j,k), ses(:,j,k),...
            'linewidth',2,'Color',cols(k,:),'LineStyle',linestyles(j));
    end
end

% ANOVAs
ys = get(gca,'Ylim');
for k = 1:size(data,3)
for i = 1:size(data,2)
    tmp = [];
    for j= 1:length(uniquecontrasts)
        tmp = [tmp; data{j,i,k}, j*ones(size(data{j,i,k},1),1)];
    end
    p = anova1(tmp(:,1),tmp(:,2),'off');
    text(0,ys(2)*(1-(k-1)/(2*size(data,3))),['p = ',num2str(p)]);
end
end

% Regression
tmp = [];
for k = [2 3]  % which gradients to look at
for i = 1:size(data,2)
    for j= 1:length(uniquecontrasts)
        tmp = [tmp; data{j,i,k},...
            j*ones(size(data{j,i,k},1),1), k*ones(size(data{j,i,k},1),1), j*k*ones(size(data{j,i,k},1),1)];
    end
end
end
tmp = [tmp, ones(size(tmp,1),1)];
[b, bint, r, rint, stats] = regress(tmp(:,1),tmp(:,2:end));
%%
% Looking for an effect of choice when the same gray square is at the RF
% In all of these trials the gray square is at the RF.
figure;
data = [];
plotcols = [.7 .5 .7; .5 .7 .5];
for i = 1:length(uniquecontrasts)
    for j = 0:1
        L = all(deltaccs == repmat(uniquecontrasts(i,:),[size(deltaccs,1) 1]),2);
        L = L & nonmatchside == 0; % 0 means T2 is correct (i.e. RF sees gray square)
        L = L & correct == j;
        trlidxs = find(L);
        sum(L)
        subplot(2,2,j+1); hold on; set(gca,'Color',plotcols(j+1,:))
        psth = zeros(1,length(bins));
        for counter = 1:sum(L)
            trlidx = trlidxs(counter);
            tmpspikes = spikes{trlidx};
            tmpspikes(tmpspikes < stimon_t(trlidx)+psthendbins(1)) = [];
            tmpspikes(tmpspikes > stimon_t(trlidx)+psthendbins(2)) = [];
            [n,x] = hist(tmpspikes-stimon_t(trlidx),bins);
            psth = psth+n;
        end
        psth = psth/sum(L)/(x(2)-x(1));
        h = plot(x,psth);
        set(h,'Color',[i i i]/length(uniquecontrasts),'LineWidth',2);
        data(i,j+1) = mean(psth);
    end
end
ymax = 0;
for i = 1:2
    subplot(2,2,i);
    tmp = get(gca,'Ylim');
    ymax = max([ymax tmp(2)]);
end
for i = 1:2
    subplot(2,2,i)
    set(gca,'Ylim',[0 ymax]);
end

subplot(2,2,3); hold on;
plot(uniquecontrasts,data(:,1),'Color',plotcols(1,:),'Linewidth',2);
plot(uniquecontrasts, data(:,2),'Color',plotcols(2,:),'Linewidth',2);

subplot(2,2,1);
title('Went to T1 (incorrectly)');
subplot(2,2,2);
title('Went to T2 (correctly)');


%%
% Signal detection theory-style simulations.
% How should I measure the size of the behavioral effect in the face of
% possible criterion shifts?
%
% x is the noise distribution, y is the signal distribution
% Model of decision making is: pick y if (y-crit)^2 > (x-crit)^2
%
% 1) Size of the effect can be estimated by the displacement between the two
% curves.  
% 2) Criterion shifts can force p(choose y) < 0.5.
% 3) Effect-specific noise can manifest as a trough with nadir > 0.5.  (If I
% were to see this in the data I would assume that she is able to see the
% gradient and can use it to her advantage.)

nsamples = 2000;
nreps = 10;
crit = 0;
effects = [-1 0 1]; % effects of gradient
effectaddnoise = .2; % effect specific noise
signalstrengths = [-4:.1:4];
figure; axes; hold on;
for j = 1:length(effects)
    data =[];
    for i = 1:length(signalstrengths)
        x = normrnd(0,1,nsamples,nreps);
        y = normrnd(signalstrengths(i)+effects(j),1+effectaddnoise*abs(effects(j)),nsamples,nreps);
        L = sum((y-crit).^2 >(x-crit).^2);
        data = [data; L/nsamples];
    end
    plot(signalstrengths,data);
end

%%
% Tuning curve
S = unique(stro.trial(:,13));
grad = unique(stro.trial(:,13));
data = [];

for i = 1:length(S)
    L = stro.trial(:,13) == S(i);
    spikes = stro.ras(L,spikeidx);
    

end
