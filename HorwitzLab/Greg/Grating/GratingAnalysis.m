% Script for playing around with gratings data
GT=nex2stro;
orients = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'orient'));
sfs = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'sf'));
diams = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'diam'));
protocols = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'protocol'));
framerate = GT.sum.exptParams.framerate;
nframes = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'nframes'));
stimon_t = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'stim_on'));
stimoff_t= GT.trial(:,strcmp(GT.sum.trialFields(1,:),'stim_off'));
spikeidx = strcmp(GT.sum.rasterCells(1,:),getSpikenum(GT));
% GH: Note this stimoff_t is the time that the n+1st frame would have
% appeared had one more frame been included.  Because the communication is
% fast, STIMFOFFCD gets dropped before this.

Lcc = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'lcont'));
Mcc = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'mcont'));
Scc = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'scont'));
colordirections = unique([Lcc Mcc Scc],'rows');

spikerates = [];
baselines = [];  baseline_t = 0.5;
for i = 1:size(GT.trial,1)
    spiketimes = GT.ras{i,spikeidx};
    nspikes = sum(spiketimes > stimon_t(i) & spiketimes < stimoff_t(i));
    spikerates = [spikerates; nspikes./(stimoff_t(i)-stimon_t(i))];
    nspikes = sum(spiketimes > stimon_t(i)-baseline_t & spiketimes < stimon_t(i));
    baselines = [baselines; nspikes./baseline_t];
end

% Looking at orientation tuning curve(s)
orienttrials = find(protocols == 1);
protocolswitches = orienttrials(find(diff(orienttrials) > 1)+1);
starts = [find(orienttrials == 1,1,'first'); protocolswitches];
stops = [];
for i = 1:length(starts)
    if (i == length(starts))
        stops = [stops; max(orienttrials(orienttrials > starts(i)))];
    else
        stops = [stops; max(orienttrials(orienttrials > starts(i) & orienttrials < starts(i+1)))];
    end
end

% Skipping earliest protocol 1 trials (if user clicks reset)
expectedntrials = length(unique(orients(orienttrials)))*GT.sum.exptParams.ntrlspercond;
if (stops(1)-starts(1)+1 > expectedntrials)
    starts(1) = stops(1)-expectedntrials+1;
end

figure; axes; hold on;
for i = 1:length(starts)
    trlidxs = [starts(i):stops(i)];
    x = orients(trlidxs);
    y = spikerates(trlidxs);
    Ltmp = x == min(x);
    y = [y; y(Ltmp)];
    x = [x; x(Ltmp)+2*pi];
        
    pp = csape(x,y,'periodic');
    xx = linspace(0,2*pi,100);
    fit = ppval(pp,xx);
    subplot(1,length(starts),i);
    hold on;
    plot(180/pi*x,y,'k.');
    plot(180/pi*xx,fit,'b-');
    sf = unique(sfs(trlidxs));  % For axis label
    if (length(sf > 1))
        disp('Error: Parameters changed during protocol 1 but "reset" not clicked');
    end
    title(['SF: ',num2str(sf),' cyc/deg']);
    xlabel('orientation (deg)');
    ylabel('response (sp/sec)');
    set(gca,'YLim',[0 max(spikerates)]);
    mu = []; sem = [];
    for j = unique(x)'
        mu(j == unique(x))= mean(y(x==j));
        sem(j == unique(x)) = std(y(x==j))/sqrt(sum(x==j));
    end
end
errorbar(unique(x)*180/pi,mu,sem);
plot([0 max(x)*180/pi], repmat(mean(baselines),1,2),'k:');
cd = unique([Lcc(stops(end)-starts(1)) Mcc(stops(end)-starts(1)) Scc(stops(end)-starts(1))],'rows')
text(0,1.1*max(mu),['color direction: ',num2str(cd*100,'%2.1d, %2.1d, %2.1d')])
axis tight

figure; axes; 
polar(unique(x)',mu+sem);
hold on;
h = polar(unique(x)',mu);
set(h,'LineWidth',2);
polar(unique(x)',mu-sem);
polar(linspace(0,2*pi,100),repmat(mean(baselines),1,100));

% Looking at spatial frequency tuning curve(s)
sftrials = find(protocols == 2);
protocolswitches = sftrials(find(diff(sftrials) > 1)+1);
starts = [min(sftrials); protocolswitches];
stops = [];
for i = 1:length(starts)
    if (i == length(starts))
        stops = [stops; max(sftrials(sftrials > starts(i)))];
    else
        stops = [stops; max(sftrials(sftrials > starts(i) & sftrials < starts(i+1)))];
    end
end
figure; axes; hold on;
for i = 1:length(starts)
    trlidxs = [starts(i):stops(i)];
    x = sfs(trlidxs);
    y = spikerates(trlidxs);
    pp = csape(x,y,'variational');
    
    xx = linspace(min(sfs),max(sfs),100);
    fit = ppval(pp,xx);
    subplot(1,length(starts),i);
    hold on;
    plot(x,y,'k.');
    plot(xx,fit,'b-');
    orient = unique(orients(trlidxs)); % for axis labeling
    title(['orientation: ',num2str(orient*180/pi),' deg']);
    set(gca,'XScale','log');
    xlabel('spatial frequency (cyc/deg)');
    ylabel('response (sp/sec)');
    set(gca,'YLim',[0 max(spikerates)]);
    mu = []; sem = [];
    for j = unique(x)'
        mu(j == unique(x))= mean(y(x==j));
        sem(j == unique(x)) = std(y(x==j))/sqrt(sum(x==j));
    end
end
errorbar(unique(x),mu,sem);
plot([min(x) max(x)], repmat(mean(baselines),1,2),'k:');
axis tight

%%
% Looking at area summation curve
if (any(protocols == 5))
    Lprotocol = protocols == 5;
    uniquesizes = unique(diams(Lprotocol));
end
mu = []; sem = [];
for i = 1:length(uniquesizes)
    L = Lprotocol & diams == uniquesizes(i);
    mu = [mu; mean(spikerates(L))];
    sem = [sem; sqrt(var(spikerates(L))/sum(L))];
end
figure;
errorbar(uniquesizes, mu, sem);
xlabel('Aperture diameter (deg)');
ylabel('sp/sec');

%%
% Looking at protocol 4 stuff (color tuning)
Lprotocol = protocols == 4;
if (any(Lprotocol))
    Lcc = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'lcont'));
    Mcc = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'mcont'));
    Scc = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'scont'));
    
    colordirections = [Lcc Mcc Scc];
    uniquecolordirs = unique(colordirections(Lprotocol,:),'rows');
    data = [];
    for i = 1:size(uniquecolordirs,1)
        L = Lprotocol &...
            colordirections(:,1) == uniquecolordirs(i,1) & ...
            colordirections(:,2) == uniquecolordirs(i,2) & ...
            colordirections(:,3) == uniquecolordirs(i,3);
        data = [data; mean(spikerates(L)) sqrt(var(spikerates(L))/sum(L))];
    end
    % Plotting
    basis = [1 1 0; 1 -1 0; 0 0 1];
    DKLs = (basis*uniquecolordirs')';
    coordinates = sign(DKLs);
    figure;
    peakfr = max(data(:,1));
    labels = {'L+M','L-M';'S','L-M';'S','L+M'};
    for i = 1:3
        subplot(1,3,i); polar(0,peakfr); hold on;
        polar(linspace(0,2*pi,100), repmat(mean(baselines),1,100));
        text(0,2*peakfr,labels{i,1},'HorizontalAlignment','center');
        text(.9*peakfr,.1*peakfr,labels{i,2});
    end
    
    for i = 1:size(uniquecolordirs,1)
        subplot(1,3,1);
        if (coordinates(i,3) == 0)
            theta = atan2(coordinates(i,1),coordinates(i,2));
            polar(theta, data(i,1),'k.');
            polar(theta+pi, data(i,1),'k.');
        end
        subplot(1,3,2);
        if (coordinates(i,1) == 0)
            theta = atan2(coordinates(i,3),coordinates(i,2));
            polar(theta, data(i,1),'k.');
            polar(theta+pi, data(i,1),'k.');
        end
        subplot(1,3,3);
        if (coordinates(i,2) == 0)
            theta = atan2(coordinates(i,3),coordinates(i,1));
            polar(theta, data(i,1),'k.');
            polar(theta+pi, data(i,1),'k.');
        end
    end
    for i = 1:3
        subplot(1,3,i);
        polar(linspace(0,2*pi,100),repmat(mean(baselines),1,100),'k:');
    end
end
%%
% Rasters for the color tuning measurements
Lprotocol = protocols == 4;
if (any(Lprotocol))
    figure;
    set(gcf,'Position',[11 56 516 887]);
    Lcc = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'lcont'));
    Mcc = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'mcont'));
    Scc = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'scont')); 
    colordirections = [Lcc Mcc Scc];
    uniquecolordirs = unique(colordirections(Lprotocol,:),'rows');
    for i = 1:size(uniquecolordirs,1)
        L = Lprotocol &...
            colordirections(:,1) == uniquecolordirs(i,1) & ...
            colordirections(:,2) == uniquecolordirs(i,2) & ...
            colordirections(:,3) == uniquecolordirs(i,3);
        trlidxs = find(L);
        subplot(size(uniquecolordirs,1),1,i); hold on;
        for j = trlidxs'
            plotidx = find(j==trlidxs);
            sp = GT.ras{j,1}-stimon_t(j);
            plot([sp sp]',[plotidx*ones(length(sp),1) plotidx*ones(length(sp),1)+.5]','k-');
            plot([stimoff_t(j)-stimon_t(j) stimoff_t(j)-stimon_t(j)],[plotidx plotidx+.5],'m-','Linewidth',2);
        end
        set(gca,'YLim',[0 plotidx+1],'YTick',[],'XTick',[],'XLim',[-.1 1]);
        title(num2str(uniquecolordirs(i,:)));
    end
    set(gca,'XTick',[0 .5 1]);
end

%%
% Experimenting with fitting a linear model to protocol 4 data and picking
% a "preferred color direction".
Lprotocol = protocols == 4;
if (any(Lprotocol))
    Lcc = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'lcont'));
    Mcc = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'mcont'));
    Scc = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'scont'));
    
    colordirections = [Lcc(Lprotocol) Mcc(Lprotocol) Scc(Lprotocol)];
    uniquecolordirs = unique(colordirections,'rows');
    responses = spikerates(Lprotocol);
    maxrespidx = find(responses == max(responses),1);
    %stimnorms = sqrt(sum(colordirections.^2,2));
    initguess = colordirections(maxrespidx,:);
    initguess = initguess.*(responses(maxrespidx)./norm(initguess))
    
    [coneweights, fval, exitflag] = fminsearch(@(x) linmodfiterr(colordirections, responses, x), initguess);

    % Let's see how well we did.
    if (exitflag)
        [tmp1,tmp2, idxs] = unique(colordirections,'rows');
        predictions = abs(colordirections*coneweights');
        figure
        subplot(2,2,1);
        plot(predictions,responses,'k.');
        axis equal;
        axis square;
        xlabel('prediction'); ylabel('response');
        title(coneweights);
        subplot(2,2,3); hold on;
        for i = 1:max(idxs)
            plot(i,responses(idxs==i)-predictions(idxs==i),'k.')
        end
        axis equal
        axis square
        xlabel('color direction'); ylabel('residual');
        subplot(2,2,4);
        for i = 1:max(idxs)
            text(1,i,num2str(tmp1(i,:)));
        end
        set(gca,'YLim',[0 max(idxs)+1], 'Xlim',[0 10],'XTick',[],'YTick',[]);
    end
end

%%
% Estimating F1/F0 on a condition by condition basis
% Using protocol 4: All colors that drive the cell decently well.
% There are a few files that have protocol 6 trials that are 
% especially designed for estimating f1/f0
binwidth = .0025;
% This smoothing filter may be unnecessary
%filterwidth = .000012; % s  
%smoothingfilter = normpdf(linspace(-3,3,6*filterwidth./binwidth),0,1);
%smoothingfilter = smoothingfilter./sum(smoothingfilter);

Lprotocol = protocols == 4;
start_t = stimon_t(Lprotocol);
end_t = stimoff_t(Lprotocol);
minstimdur = min(end_t-start_t);
if (any(protocols == 4))
    Lcc = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'lcont'));
    Mcc = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'mcont'));
    Scc = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'scont'));
    
    colordirections = [Lcc(Lprotocol) Mcc(Lprotocol) Scc(Lprotocol)];
    uniquecolordirs = unique(colordirections,'rows');
    bins = [0:binwidth:minstimdur];
    data = [];
    for i = 1:size(uniquecolordirs,1)
        L = Lprotocol &...
            Lcc == uniquecolordirs(i,1) & ...
            Mcc == uniquecolordirs(i,2) &... 
            Scc == uniquecolordirs(i,3);
        start_t = stimon_t(L);
        end_t = stimoff_t(L);
        spiketimes = GT.ras(L,1);
        totalspikes = [];
        for j = 1:sum(L)
            totalspikes = [totalspikes; spiketimes{j}-start_t(j)];
        end
        totalspikes(totalspikes > minstimdur) = [];
        psth = histc(totalspikes',bins)/binwidth/sum(L);

     %   psth = conv(psth,smoothingfilter,'full');
     %   psth(1:length(smoothingfilter)/2-1) = [];
     %   psth(end:-1:end-length(smoothingfilter)/2) = [];

        k = unique(GT.trial(L,strcmp(GT.sum.trialFields(1,:), 'tf')));
        if (length(k) > 1)
            error('Too many temporal frequencies!');
        end
        basis0 = ones(length(psth),1);
        basis1 = exp(-2*pi*sqrt(-1)*k*[0:length(psth)-1]/length(psth))';
        F0 = abs((psth-mean(baselines))*basis0);
        F1 = 2*abs((psth-mean(baselines))*basis1);
        % Not sure what the logic for the the "2*" in the line above
        % but it makes a half-wave rectified sinewave have a modulation
        % ratio of pi/2, which is what it's supposed to have.
     %   F0 = abs(sum(psth));
     %   F1 = abs(sum((psth)*basis1));

        data = [data; F1 F0];
    end
    figure; axes; hold on;
    plotcols = uniquecolordirs
    for i = 1:size(uniquecolordirs,1)
        col = uniquecolordirs(i,:).*[1 1 .2];
        col = col./(2*max(abs(col))) + .5;        
        h = plot(data(i,2)/length(psth),data(i,1)./data(i,2),'ko','Markersize',10,'MarkerFaceColor',col);
    end
    xlabel('mean spike rate');
    ylabel('modulation ratio');
    legend(num2str(uniquecolordirs));
  
    L = data(:,2)/length(psth) > 10;
    if (any(L))
        unweighted = mean(data(L,1)./data(L,2));
        weighted = sum(data(L,2).*(data(L,1)./data(L,2)))/sum(data(L,2));
        title(['unweighted: ',num2str(unweighted,3),' weighted: ',num2str(weighted,3)]);
    else
        title('No conditions drove cell > 10 sp/sec');
    end
end

%%
% Looking at Gabor contrast-response functions
if (any(protocols == 7))
    figure; axes; 
    Lprotocol = protocols == 7;
    responses = spikerates(Lprotocol);
    LMSs = [Lcc(Lprotocol) Mcc(Lprotocol) Scc(Lprotocol)];
    Lblank = all(LMSs == 0,2);
    LvsM = dot(LMSs, repmat([1./sqrt(2) -1./sqrt(2) 0],size(LMSs,1),1),2);
    S = dot(LMSs, repmat([0 0 1],size(LMSs,1),1),2)
    uniqueS = unique(S(S>0));
    uniqueLvsM = unique(LvsM(LvsM>0));
    
    % First S
    mn = [];
    se = [];
    for i = 1:size(uniqueS)
        L = S == uniqueS(i);
        mn(i)= mean(responses(L));
        se(i) = std(responses(L)./sqrt(sum(L)));
    end
    subplot(2,1,1); hold on;
    plot(uniqueS,mn);
    errorbar(uniqueS, mn, se);
    title('S');


    % Now L vs M 
    mn = [];
    se = [];
    for i = 1:size(uniqueLvsM)
        L = LvsM == uniqueLvsM(i);
        mn(i)= mean(responses(L));
        se(i) = std(responses(L)./sqrt(sum(L)));
    end
    subplot(2,1,2); hold on;
    plot(uniqueLvsM,mn);
    errorbar(uniqueLvsM, mn, se);
    title('L-M');
end
%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%

%%
% Estimating F1/F0 from drifting gratings
% Doing it on a trial-by0trial basis.
% This stuff didn't work that well
sacstats = getSacData(GT);
bindur = .0001;  % s
filterwidth = .01; % s  
smoothingfilter = normpdf(linspace(-3,3,6*filterwidth./bindur),0,1);
smoothingfilter = smoothingfilter./sum(smoothingfilter);

% Setting up an L vector of trials to analyze
% Just temporary for now.
tfs = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'tf'));
%L = spikerates-mean(baselines) > 20;

L = protocols == 6 & spikerates-mean(baselines) > 10;
if (sum(L)==0)
    L = protocols == 4 & spikerates-mean(baselines) > 10;
end

trlidxs = find(L);

% % Visualizing the rasters and eye movements
% % for diagnostics
% figure;
clf;
set(gcf,'Position',[44 111 367 846]);
for i = trlidxs'
    subplot(length(trlidxs),1,find(i == trlidxs)); hold on;
    starteye = GT.ras{i,strcmp(GT.sum.rasterCells,'anlgStartTime')};
    eyeh = GT.ras{i,strcmp(GT.sum.rasterCells,'AD11')};
    eyev = GT.ras{i,strcmp(GT.sum.rasterCells,'AD12')};
    t = linspace(starteye, starteye+length(eyeh)/500,length(eyeh));
    sp = GT.ras{i,1}-stimon_t(i);
    plot([sp sp]',[zeros(length(sp),1) .2*ones(length(sp),1)]','k-');
    plot([0 0],[0 .5],'m-','Linewidth',2);
    plot([stimoff_t(i)-stimon_t(i) stimoff_t(i)-stimon_t(i)],[0 .5],'m-','Linewidth',2);
    
    Leye = t>stimon_t(i) & t<stimoff_t(i);
    plot(t(Leye)-stimon_t(i),eyeh(Leye)-mean(eyeh(Leye))*5,'r-');
    plot(t(Leye)-stimon_t(i),eyev(Leye)-mean(eyev(Leye))*5,'g-');
    set(gca,'YLim',[-.1 .5],'YTick',[],'XTick',[],'XLim',[-.1 1]);
end

% Now doing the analysis on a trial by trial basis
F0 = [];
F1 = [];
snippetlen = [];
for i = trlidxs'
    startsac_t = sacstats.starttimes{i};
    endsac_t = sacstats.endtimes{i}+.1; % cutting out from sac start to 100 ms after sac end
    L = startsac_t < stimon_t(i) | endsac_t > stimoff_t(i); % getting rid of saccades that span stimon or stimoff
    startsac_t(L) = [];
    endsac_t(L) = [];
    % Counterintuitively, we want to start epochs at the *ends* 
    % of saccades and end them at the *beginings* of saccades.
    start_t = [stimon_t(i); endsac_t]; 
    end_t = [startsac_t;stimoff_t(i)]; 
    
    % getting intersaccade times
    epochdur = end_t-start_t;
 
    for j = 1:length(epochdur)
        if (epochdur(j) < .25/tfs(i))  % Not using very short epochs
            continue;
        end
        psth = histc(GT.ras{i},[start_t(j):bindur:end_t(j)]);
        psth = conv(psth,smoothingfilter,'full');
        psth = psth/bindur;  % To get back to sp/sec (which is the units of "baseline"). 
        k = (end_t(j)-start_t(j))/(1/tfs(i)); % k is the number of cycles in the snippet
        basis1 = exp(-2*pi*sqrt(-1)*k*[0:length(psth)-1]/length(psth))';
        F0(length(F0)+1) = abs(sum(psth-mean(baselines)));
        F1(length(F1)+1) = abs(sum((psth-mean(baselines))'*basis1));
        snippetlen(length(snippetlen)+1) = length(psth);        
    end
end
modRatio = nanmean(F1./F0);
subplot(length(trlidxs),1,1);
title(['Modulation ratio = ',num2str(modRatio)]);


%%
% Looking at spatial phase tuning curve a la Williams and Shapley.
% Protocol 3.  This stuff never worked that well.

if (any(protocols == 3))
    spikerates = [];
    offsets = [0.02 .1];  % Both relative to stimon
    for i = 1:size(GT.trial,1)
        spiketimes = GT.ras{i,1};
        nspikes = sum(spiketimes > stimon_t(i)+offsets(1) & spiketimes < stimon_t(i)+offsets(2));
        spikerates = [spikerates; nspikes./(offsets(2)-offsets(1))];
    end

    % Getting eye positions
    % Fixation saccades
    sacstats = getSacData(GT);
    saccadeoccurred = nan*ones(size(GT.trial,1),1);
    for i = 1:size(GT.trial,1)
        st = sacstats.starttimes{i};
        Lsac = (st > stimon_t(i)-.02) & (st < stimoff_t(i));
        if (any(Lsac))
            saccadeoccurred(i) = 1;
        else
            saccadeoccurred(i) = 0;
        end
    end

    phasetrials = find(protocols == 3);
    phases = GT.trial(:,find(strcmp(GT.sum.trialFields(1,:),'phase')));
    uniquephases = unique(phases(phasetrials));
    means = [];
    sems = [];
    for i = 1:length(uniquephases)
        L = protocols == 3 & phases == uniquephases(i) & ~saccadeoccurred;
        means = [means; mean(spikerates(L))];
        sems = [sems; sqrt(var(spikerates(L))/sum(L))];
    end
    %means = [means; means(1)]
    %sems = [sems; sems(1)]
    %uniquephases = [uniquephases; uniquephases(1)+2*pi];

    bl = mean(baselines(protocols == 3));
    figure;
    errorbar(uniquephases, means-bl, sems);
    ylabel('response (sp/sec)');
    xlabel('phase (rad)');
    basis1 = exp(-2*pi*sqrt(-1)*[0:length(uniquephases)-1]/length(uniquephases))';
    F1 = abs(sum((means-bl)'*basis1));
    F0 = abs(sum(means-bl));

    %abs(fft(means-bl))
    title(['F1/F0 = ',num2str(F1/F0)]);
end

b = exp(-2*pi*i*[0:6]/7)
abs(sum(a))

% For some reason I have to multiply F0 by 0.5 so that a halfwave
% rectified sinewave give f1/f0 = pi/2

% This stuff is all consistent, but doesn't give pi/2 for halfwave
% rectified sinwave
%abs(fft(means))
%F0 = sum(means);
%F1 = sqrt(sum(((means)'*basis).^2));

% Looking at the rasters.
% This should give an impression about the quality of the eye position
% compensation

spikename = getSpikenum(GT);
spikeidx = find(strcmp(GT.sum.rasterCells(1,:),spikename));
spikes = GT.ras(:,spikeidx);
stimon_t = GT.trial(:,2);
%stimoff_t = GT.trial(:,3);

phasetrials = find(protocols == 3);
if (any(phasetrials))
phases = GT.trial(:,find(strcmp(GT.sum.trialFields(1,:),'phase')));
uniquephases = unique(phases(phasetrials));
figure;
for i = 1:length(uniquephases)
    subplot(length(uniquephases),1,i);
    hold on;
    L = protocols == 3 & phases == uniquephases(i) & ~saccadeoccurred;
    trlidxs = find(L);
    for counter = 1:sum(L)
        trlidx = trlidxs(counter);
        %plot(0,counter,'g*');
        plot(stimoff_t(trlidx)-stimon_t(trlidx),counter,'r*');
        nspikestot = length(spikes{trlidx});
        plot([spikes{trlidx} spikes{trlidx}]'-stimon_t(trlidx),[zeros(nspikestot,1) .5*ones(nspikestot,1)]'+counter,'k-');
    end
    plot([0 0],[0 counter+1],'g-');
    plot([offsets(1) offsets(1)],[0 counter+1],'k:');
    plot([offsets(2) offsets(2)],[0 counter+1],'k:');  % Spike counting windows
    
    set(gca,'XLim',[-.1 .4],'YLim',[0 counter+1],'YTick',[]);
    if (i<length(uniquephases))
       set(gca,'XTick',[]) 
    end
end
end
