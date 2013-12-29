%  Code for analyzing FixStim data.
%
% Section 1: Rasters and PSTHs for FixStim data
% Section 2: Eye position plots for FixStim data
% Section 3: Endpoints of saccades made within a particular window
% Section 4: Power series 
% Section 5: Analysis of FixStim_interfreq data
% Section 6: Fixational eyemovements around the time of optostim
% Section 7: Looking at optostim-driven spikes vs microsaccades times
% Section 8: FixStim_contrastresponse analyses

%%
% Section 1: Rasters and PSTHs

stro = nex2stro;
synctimestr = 'FPOFF';
stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
fpoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fpoff_t'));
stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimoff_t'));
targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targon_t'));
targoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targoff_t'));
opt_stimfreq = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'opt_stimfreq'));
targ_shown = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_shown'));
stim_type = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_type'));
optfreq = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'opt_stimfreq'));

if (strcmp(synctimestr,'FPOFF'))
    sync_t = fpoff_t; disp('All trials aligned to FP off');
else
    sync_t = stimon_t;  % Ecode for alignment
    L = isnan(sync_t);
    sync_t(L) = targon_t(L); % sync of targon if no optical stim
    L = isnan(sync_t);
    sync_t(L) = fpoff_t(L); % sync on fpoff if no optical nor visual stim
end
offset = [-.3 .7];  % pre and post time
Lspikechans = strncmp(stro.sum.rasterCells,'sig0',4);

if (any(Lspikechans))
    for whichspike = find(Lspikechans)
        whichspike
        spikes = stro.ras(:,whichspike);
        binwidth = .005;
        bins = offset(1):binwidth:offset(2);
        PSTH = zeros(4,length(bins));
        figure; subplotcntr = 0;
        for i = 0:1  % Target appearing or not
            for j = unique(stim_type)' % optical stim or not
                subplotcntr=subplotcntr+1; subplot(2,2,subplotcntr); hold on;
                L = targ_shown == i & stim_type == j;
                trlidxs = find(L);
                for counter = 1:sum(L)
                    trlidx = trlidxs(counter);
                    tmpspikes = spikes{trlidx}-sync_t(trlidx);
                    tmpspikes(tmpspikes < offset(1) | tmpspikes > offset(2)) = [];
                    nspikestot = length(tmpspikes);
                    plot([tmpspikes tmpspikes]',[zeros(nspikestot,1) .5*ones(nspikestot,1)]'+counter,'k-');
                    PSTH(subplotcntr,:) = PSTH(subplotcntr,:) + hist(tmpspikes, bins);
                end
                if (j > 0) % Plotting the time course of optical stimulation
                    dur = mode(stimoff_t-stimon_t);
                    secspercycle = 1/unique(optfreq(optfreq > 0));
                    transitions = 0:secspercycle/2:dur;
                    x = [transitions; transitions];
                    x = [offset(1); x(:); offset(2)];
                    y = [repmat([0 1],1,length(transitions)/2) 0]*2-2;
                    y = [y;y];
                    if (length(x(:)) == length(y(:)))
                        plot(x,y(:)','k-','linewidth',2);
                        %  else
                        %      keyboard;
                    end
                end
                PSTH(subplotcntr,:) = PSTH(subplotcntr,:)/binwidth/sum(L);
                set(gca,'XLim', offset,'Ytick',[],'YLim',[-4 sum(L)+5]);
                titlestr = [];
                if (i == 1)
                    titlestr = cat(2, titlestr,'Targ');
                    xlabel('Time (s)');
                end
                if (i == 1 & j ~= 0)
                    titlestr = cat(2, titlestr,' + ');
                end
                if (j == 2)
                    titlestr = cat(2, titlestr,' optical stim (',[num2str(unique(optfreq(optfreq > 0))),' Hz)']);
                end
                title(titlestr);
            end
        end
        
        % PSTHs
        figure;
        for i = 1:4
            subplot(2,2,i);
            plot(bins,PSTH(i,:),'k-','LineWidth',2);
            set(gca,'YLim',[0 10*ceil(max(PSTH(:)/10))]);
            set(gca,'Xlim',offset);
            if (i == 1) 
                title(['All trials aligned to ',synctimestr]);
            end
        end
    end
end
%%
% Section 2)
% Eye position plots for the various trial types

ALIGNEPS = 0; % Set eye_h and eye_v = 0 at t = sync_t 
H = stro.ras(:,strcmp(stro.sum.rasterCells,'AD11'));
V = stro.ras(:,strcmp(stro.sum.rasterCells,'AD12'));
ep_t = stro.ras(:,strcmp(stro.sum.rasterCells,'anlgStartTime'));
ADfreq = stro.sum.analog.storeRates{1};

sync_t = stimon_t;  % Ecode for alignment
L = isnan(sync_t);
sync_t(L) = targon_t(L); % sync of targon if no optical stim
L = isnan(sync_t);
sync_t(L) = fpoff_t(L); % sync on fpoff if no optical nor visual stim

% Aligning on the appearance of the target
%sync_t(logical(targ_shown)) = targon_t(logical(targ_shown)); 
% Aligning on the electrical stimulation
%sync_t(~isnan(stimon_t)) = stimon_t(~isnan(stimon_t)); 

% Xt and Yt plots
subplotcntr = 0;
figure;  % ({x,y},t)
for i = 0:1  % Target appearing or not
    for j = [0 1 2]  % No stim, electrical stim, optical stim
        L = stim_type == j & targ_shown == i;
        trlidxs = find(L);
        subplotcntr=subplotcntr+1;
        subplot(2,3,subplotcntr); hold on;
        for counter = 1:sum(L)
            trlidx = trlidxs(counter);
            h = H{trlidx}*4096/400;
            v = V{trlidx}*4096/400;
            t = [0:length(h)-1]./ADfreq+ep_t{trlidx};
            t = t-sync_t(trlidx);
            if (ALIGNEPS)
                idx = find(t.^2 == min(t.^2));
                h = h-h(idx);
                v = v-v(idx);
            end
            plot(t,h,'g-');
            plot(t,v,'r-');
        end
%        set(gca,'Xlim',offset,'Ylim',[-4 3]);
        set(gca,'Xlim',offset);
    end
end
%equatesubplotaxeslims;


% Looking at the X,Y components of saccade. Taking EP from 100 to 0 ms before
% FPOFF and subtracting EP at 200 to 300 ms.
data = {};
for i = 0:1  % Target appearing or not
    for j = [0 1 2]  % No stim, electrical stim, optical stim
        L = targ_shown == i & stim_type == j;
        trlidxs = find(L);
        tmp = [];
        for counter = 1:sum(L)
            trlidx = trlidxs(counter);
            h = H{trlidx}*4096/400;
            v = V{trlidx}*4096/400;
            t = [0:length(h)-1]./ADfreq+ep_t{trlidx};
            t = t-sync_t(trlidx);
            L1 = t > -.1 & t < 0;
            L2 = t > .2 & t < .3;
            hcomp= mean(h(L2))-mean(h(L1));
            vcomp= mean(v(L2))-mean(v(L1));
            tmp = [tmp; hcomp, vcomp];
        end
        data{i+1, j+1} = tmp;
    end
end
figure; subplotcntr = 0;
for i = 0:1  % Target appearing or not
    for j = [0 1 2]  % No stim, electrical stim, optical stim
        subplotcntr=subplotcntr+1;
        subplot(2,3,subplotcntr); hold on;
        tmp = data{i+1,j+1};
        if (~isempty(tmp))
            plot(tmp(:,1),tmp(:,2),'k.');
        end
    end
end
%equatesubplotaxeslims;

% XY plot
subplotcntr = 0;
figure;  % (x,y)
for i = 0:1  % Target appearing or not
    for j = [0 1 2]  % No stim, electrical stim, optical stim
        L = targ_shown == i & stim_type == j;
        trlidxs = find(L);
        subplotcntr=subplotcntr+1;
        subplot(2,3,subplotcntr); hold on;
        for counter = 1:sum(L)
            trlidx = trlidxs(counter);
            h = H{trlidx}*4096/400;
            v = V{trlidx}*4096/400;
            t = [0:length(h)-1]./ADfreq+ep_t{trlidx};
            t = t-sync_t(trlidx);
            h(t<0) = [];
            v(t<0) = [];
            
            if (ALIGNEPS)
                idx = find(t.^2 == min(t.^2));
                h = h-h(idx);
                v = v-v(idx);
            end
            plot(h,v,'k-');
        end
        if (i == 0)
            set(gca,'Xlim',[-30 30],'Ylim',[-30 30]);
%                        set(gca,'Xlim',[-4 4],'Ylim',[-4 4]);   
        else
            set(gca,'Xlim',[-4 4],'Ylim',[-4 4]);            
        end
    end
end
%equatesubplotaxeslims;



%%
% Section 3
% Metrics of saccades made within 300 ms of fixation point offset.
% Noise in the eye traces will be a pain.
sacstats = getSacData(stro); close;
data = cell(3,2);
SACTWIN = [.1 .45];
MINAMP = .5;
tmp = [];
for i = 1:size(stro.trial,1)
    latencies = sacstats.starttimes{i}-fpoff_t(i);
    amps = sacstats.amplitudes{i};
    idx = find(latencies > SACTWIN(1) & latencies < SACTWIN(2) & amps > MINAMP,1,'first');
    tmp = data{stim_type(i)+1,targ_shown(i)+1};
    if (isempty(idx))
        tmp = [tmp; 0 0];
    else
        [x,y] = pol2cart(sacstats.directions{i}(idx),sacstats.amplitudes{i}(idx));
        tmp = [tmp; x y];
    end
    data{stim_type(i)+1,targ_shown(i)+1} = tmp;
end

targx = stro.sum.exptParams.rf_x/10;
targy = stro.sum.exptParams.rf_y/10;
unitvecttotarg = [targx; targy]./norm([targx;targy]);
figure;
tmplims = [];
for i = 1:numel(data)
    subplot(2,3,i); hold on;
    if (~isempty(data{i}))
        compass(data{i}(:,1),data{i}(:,2));
    end
    plot(targx, targy,'ms','MarkerSize',5);
    axis equal;
    tmplims = [tmplims; get(gca,'XLim') get(gca,'YLim')];
end
PROJTHRESH = 1;
for i = 1:numel(data)
    subplot(2,3,i); hold on;
    set(gca,'XLim',[min(tmplims(:,1)) max(tmplims(:,2))]);
    set(gca,'YLim',[min(tmplims(:,3)) max(tmplims(:,4))]);
    plot([min(tmplims(:,1)) max(tmplims(:,2))],...
        [(PROJTHRESH-unitvecttotarg(1)*min(tmplims(:,1)))/unitvecttotarg(2),...
        (PROJTHRESH-unitvecttotarg(1)*max(tmplims(:,2)))/unitvecttotarg(2)],'k:');
end

% getting numbers of saccades into the relevant quadrant per condition
nsacstotarg = zeros(3,2);

for i = 1:numel(data)
    if (any(data{i}))
        projections = data{i}*unitvecttotarg;
        nsacstotarg(i) = sum(projections > PROJTHRESH);
    end
end

% Bar plot showing numbrs of saccades made to relevant area
figure; axes; hold on;
bar(nsacstotarg')
set(gca,'XTick',[1 2],'XTicklabel',{'targ absent','targ present'});
ylabel('# saccades into relevant quadrant');
legend({'no stim','electrical','optical'},'location','NorthWest');
set(gcf,'Name',stro.sum.fileName);

% Comparing optical to no stim
for i = 1:2
    ns = [size(data{1,i},1) size(data{3,i},1)]
    [h,p] = equalproptest([nsacstotarg(1,i) nsacstotarg(3,i)],ns,0.01);
    if (h)
       plot(i,max(nsacstotarg(:,i)),'k*','MarkerSize',5);
    end
end


%%
% % Section 4
% % Power series on behavior
% stro = nex2stro;
% fpoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fpoff_t'));
% targ_shown = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_shown'));
% SACTWIN = [.1 .5]; % Saccade latencies must be between first and last
% MINAMP = .5; % Minimum saccade amplitude
% 
% sacstats = getSacData(stro); close;
% laserpowers = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'laser_power'));
% uniquepowers = unique(laserpowers);
% 
% data = cell(length(uniquepowers),2);
% tmp = [];
% for i = 1:size(stro.trial,1)
%     latencies = sacstats.starttimes{i}-fpoff_t(i);
%     amps = sacstats.amplitudes{i};
%     idx = find(latencies > SACTWIN(1) & latencies < SACTWIN(2) & amps > MINAMP,1,'first');
%     tmp = data{uniquepowers==laserpowers(i),targ_shown(i)+1};
%     if (isempty(idx))
%         tmp = [tmp; nan nan];
%  %       tmp = [tmp; 0 0]; % Count a "no-saccade" trial as a saccade to (0,0)?
%     else
%         [x,y] = pol2cart(sacstats.directions{i}(idx),sacstats.amplitudes{i}(idx));
%         tmp = [tmp; x y];
%     end
%     data{uniquepowers==laserpowers(i),targ_shown(i)+1} = tmp;
% end
% 
% figure; axes; hold on;
% cmap = hot(size(data,1));
% for i = 1:size(data,1)
%     h = plot(data{i,1}(:,1),data{i,1}(:,2),'ko');
%     set(h,'MarkerFaceColor',cmap(i,:),'MarkerEdgeColor','none');
% end
% set(gca,'Color',[.2 .2 .2])
% 
% 
% % Mehrdad's mean to variance plot
% rf = [stro.sum.exptParams.rf_x stro.sum.exptParams.rf_y]/10;
% %if (rf(1) < 0) rf = -rf; end  % RF is always on the right (for Sedna)
% rf =  [1.2000   -2.6000]
% dists = [];
% for i = 1:size(data,1)
%     xy = [data{i,1}(:,1),data{i,1}(:,2)];
%     xydiff = xy-repmat(rf,size(xy,1),1);
%     dists = [dists; repmat(uniquepowers(i),size(xy,1),1) sqrt(sum(xydiff.^2,2))];
% end
% % Figure 1
% figure; axes; hold on;
% boxplot(dists(:,2),dists(:,1));
% %set(gca,'Color',[.2 .2 .2])
% set(h,'MarkerFaceColor',cmap(i,:),'MarkerEdgeColor','none');
% 
% % Figure 2
% figure; axes; hold on;
% for i = 1:size(data,1)
%     xy = [data{i,1}(:,1),data{i,1}(:,2)];
%     xydiff = xy-repmat(rf,size(xy,1),1);
%     dist = sqrt(sum(xydiff.^2,2));
%     xydiff2 = xy-repmat(nanmean(xy),size(xy,1),1);
%     dist2 = sqrt(sum(xydiff2.^2,2));
%     h = plot(nanmean(dist),nanmean(dist2),'ko');
%     set(h,'MarkerFaceColor',cmap(i,:),'MarkerEdgeColor','none');
% end
% set(gca,'Color',[.2 .2 .2])
% %xlabel('median distance'); ylabel('IQR of distances');
% 
% % Figure 3
% figure;
% for i = 1:size(data,1)
%     xy = [data{i,1}(:,1),data{i,1}(:,2)];
%    % xydiff = xy-repmat(rf,size(xy,1),1);
%    % dist = sqrt(sum(xydiff.^2,2));
%      dist = norm(nanmean(xy)-rf);
%     
%     xydiff2 = xy-repmat(nanmean(xy),size(xy,1),1);
%     dist2 = sqrt(sum(xydiff2.^2,2));
%     biasval(i) = nanmean(dist);
%     stdval(i) = nanmean(dist2);    
% end    
% subplot(2,1,1);hold on;
% plot(uniquepowers/255*50,biasval,'k');
% subplot(2,1,2);hold on;
% plot(uniquepowers/255*50,stdval,'k');
% 
% for i = 1:size(data,1)
%     xy = [data{i,2}(:,1),data{i,2}(:,2)];
%    % xydiff = xy-repmat(rf,size(xy,1),1);
%    % dist = sqrt(sum(xydiff.^2,2));
%      dist = norm(nanmean(xy)-rf);
%     
%     xydiff2 = xy-repmat(nanmean(xy),size(xy,1),1);
%     dist2 = sqrt(sum(xydiff2.^2,2));
%     biasval(i) = nanmean(dist);
%     stdval(i) = nanmean(dist2);    
% end    
% subplot(2,1,1);hold on;
% plot(uniquepowers/255*50,biasval,'r');
% subplot(2,1,2);hold on;
% plot(uniquepowers/255*50,stdval,'r');

% subplot(2,1,1);
% ylabel('Bias (deg)'); %xlabel('Power (mW)');
% subplot(2,1,2);
% ylabel('Standard deviation (deg)'); xlabel('Power (mW)');
% 
% Code for analyzing data from Zack's Fixstim_freq data

%%
%{
% Section 5
% Hacked together analysis of Fixstim_freq data
stro = nex2stro;
offset = [-.1 .1];  % pre and post time wrt stimon/stimoff
stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
fpoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fpoff_t'));
stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimoff_t'));
targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targon_t'));
targoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targoff_t'));
opt_stimfreq = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'opt_stimfreq'));
targ_shown = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_shown'));
optfreq = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'elec_stimfreq'));
uniqfreqs = unique(optfreq);

% Hack below 
if (isnan(targ_shown(1)))
    targ_shown(:,1) = 0;
end
dur = mode(stimoff_t-stimon_t);

sync_t = stimon_t;  % Ecode for alignment
Lspikechans = strncmp(stro.sum.rasterCells,'sig0',4);

if (any(Lspikechans))
    for whichspike = find(Lspikechans)
        spikes = stro.ras(:,whichspike);
        binwidth = .005;
        bins = offset(1):binwidth:dur+offset(2);
        PSTH = zeros(1,length(bins));
        for j = uniqfreqs'
            PSTH = zeros(1,length(bins));
            figure; subplot(2,1,1); hold on;
            L = optfreq == j;
            trlidxs = find(L);
            for counter = 1:sum(L)
                trlidx = trlidxs(counter);
                tmpspikes = spikes{trlidx}-sync_t(trlidx);
                tmpspikes(tmpspikes < offset(1) | tmpspikes > dur+offset(2)) = [];
                nspikestot = length(tmpspikes);
                plot([tmpspikes tmpspikes]',[zeros(nspikestot,1) .5*ones(nspikestot,1)]'+counter,'k-');
                PSTH = PSTH + hist(tmpspikes, bins);
            end
            if (j > 0) % Plotting the time course of optical stimulation
                secspercycle = 1/unique(j(j > 0));
                transitions = 0:secspercycle/2:dur;
                if (ceil(length(transitions)/2) ~= floor(length(transitions)/2))
                    transitions(end+1) = dur; %automatic shutoff
                end
                x = [transitions; transitions];
                x = [offset(1); x(:); offset(2)];
                y = [repmat([0 1],1,length(transitions)/2) 0]*2-2;
                y = [y;y];
                if (length(x(:)) == length(y(:)))
                    plot(x,y(:)','k-','linewidth',2);
                end
            end
            set(gca,'XLim', [offset(1) dur+offset(2)],'Ytick',[],'YLim',[-4 sum(L)+5]);
            title(['Frequency: ',num2str(j)]);
            
            % PSTH
            subplot(2,1,2); hold on;
            plot(bins,PSTH,'k-','LineWidth',2);
            set(gca,'YLim',[0 10*ceil(max(PSTH(:)/10))]);
            set(gca,'Xlim',[offset(1) dur+offset(2)]);
        end
    end
end

%%
% Section 6
% Looking at fixational eye movements synced to optostim pulses
% Assuming Fixstim_freq paradigm for now

H = stro.ras(:,strcmp(stro.sum.rasterCells,'AD11'));
V = stro.ras(:,strcmp(stro.sum.rasterCells,'AD12'));
ep_t = stro.ras(:,strcmp(stro.sum.rasterCells,'anlgStartTime'));
ADfreq = stro.sum.analog.storeRates{1};
sync_t = stimon_t;  % Someday align on every pulse, not just the first

% Xt and Yt plots

for j = uniqfreqs'
    figure; axes; hold on;
    hs = [];
    vs = [];
    L = optfreq == j;
    trlidxs = find(L);
    % Plotting the time course of optical stimulation
    % doing it first because we want it in the background
    if (j > 0)
        dur = mode(stimoff_t-stimon_t);
        secspercycle = 1/unique(j(j > 0));
        transitions = 0:secspercycle/2:dur;
        if (ceil(length(transitions)/2) ~= floor(length(transitions)/2))
            transitions(end+1) = dur; %automatic shutoff
        end
        x = [transitions; transitions];
        x = [offset(1); x(:); offset(2)];
        y = [repmat([0 1],1,length(transitions)/2) 0]/4 - .5;
        y = [y;y];
        if (length(x(:)) == length(y(:)))
            plot(x,y(:)','k-','linewidth',1);
        end
    end
    for counter = 1:sum(L)
        trlidx = trlidxs(counter);
        hep = H{trlidx}*4096/400;
        vep = V{trlidx}*4096/400;
        t = [0:length(hep)-1]./ADfreq+ep_t{trlidx};
        t = t-sync_t(trlidx);
        h = plot(t,hep,'g-'); set(h,'Color',[.75 1 .75])
        h = plot(t,vep,'r-'); set(h,'Color',[1 .75 .75])
        hs(counter,:) = hep(t > offset(1) & t < offset(2));
        vs(counter,:) = vep(t > offset(1) & t < offset(2));
    end
    plot_t = linspace(offset(1),offset(2),size(hs,2));
    plot(plot_t,mean(hs),'g-','Linewidth',2);
    plot(plot_t,mean(vs),'r-','Linewidth',2);
    set(gca,'Xlim',offset,'Ylim',[-.6 .6]);
    title(['Frequency: ',num2str(j)]);
end

%%
% Section 7
% Comparing counts of optostim produced spikes after microsaccades 
stro = nex2stro;
sacstats = getSacData(stro); close;
stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targon_t'));
L = ~isnan(stimon_t); % & isnan(targon_t); including target trials for now 
stimtrials = find(L)';
data = [];
for i = stimtrials
    spiketimes = stro.ras{i,1}-stimon_t(i);
    nspikes = sum(spiketimes > 0 & spiketimes < .01);
    sactimes = sacstats.starttimes{i}-stimon_t(i);
    lastsactime = max(sactimes(sactimes < 0));
    if (isempty(lastsactime))
        lastsactime = nan;
    end
    data = [data; lastsactime nspikes];
end
data(isnan(data(:,1)),:) = [];
plot(data(:,1),data(:,2),'k.'); lsline;
[r,p] = corrcoef(data);
title(['r = ',num2str(r(1,2)),' p = ',num2str(p(1,2))]);


%%
% Section 8
% FixStim_contrastresponse analyses
stro = nex2stro;
offset = [-.1 .75];
stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimoff_t'));
targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targon_t'));
targoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targoff_t'));
targ_shown = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_shown'));
laser_power = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'laser_power'));
contrasts = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'contrast'));
uniquecontrasts = unique(contrasts);
uniquepowers = unique(laser_power);

sync_t = targon_t;  % Ecode for alignment
Lspikechans = strncmp(stro.sum.rasterCells,'sig0',4);

if (any(Lspikechans))
    for whichspike = find(Lspikechans)
        spikes = stro.ras(:,whichspike);
        binwidth = .005;
        bins = offset(1):binwidth:offset(2);
        PSTH = zeros(1,length(bins));
        for j = uniquecontrasts'
            figure;
            peakPSTH = 0;
            for i = uniquepowers'
                PSTH = zeros(1,length(bins));
                subplot(2,length(uniquepowers),find(i == uniquepowers)); hold on;
                L = contrasts == j & laser_power == i;
                trlidxs = find(L);
                for counter = 1:sum(L)
                    trlidx = trlidxs(counter);
                    tmpspikes = spikes{trlidx}-sync_t(trlidx);
                    tmpspikes(tmpspikes < offset(1) | tmpspikes > offset(2)) = [];
                    nspikestot = length(tmpspikes);
                    plot([tmpspikes tmpspikes]',[zeros(nspikestot,1) .5*ones(nspikestot,1)]'+counter,'k-');
                    PSTH = PSTH + hist(tmpspikes, bins);
                end
                set(gca,'XLim', offset,'Ytick',[],'YLim',[-4 sum(L)+5]);
                title(['Contrast: ',num2str(j),' Power: ',num2str(i)]);
                
                % PSTH
                subplot(2,length(uniquepowers),length(uniquepowers)+find(i == uniquepowers)); hold on;
                plot(bins,PSTH,'k-','LineWidth',2);
                set(gca,'YLim',[0 10*ceil(max(PSTH(:)/10))]);
                set(gca,'Xlim',offset);
                peakPSTH = max([peakPSTH, PSTH]);
            end
            for i = uniquepowers'
                subplot(2,length(uniquepowers),length(uniquepowers)+find(i == uniquepowers)); 
                set(gca,'Ylim',[0 peakPSTH*1.1])
            end
        end
    end
end

%%
% Caculating luminance values for targets
R = stro.sum.exptParams.targ_r;
G = stro.sum.exptParams.targ_g;
B = stro.sum.exptParams.targ_b;
load Dell4BitsCal;
cal = cals{end};
cal.bgColor
p = SplineRaw([380:4:780]', cal.P_device, [380:5:780]');

load T_xyzJuddVos
vl = T_xyzJuddVos(2,:);
bkgndlum = vl*(p*cal.bgColor);  % this is 90 cd/m^2
vl = vl*90/bkgndlum; %scling vlambda
vl*SplineRaw([380:4:780]', cal.P_ambient, [380:5:780]')  % The black fixation point

% Computing the luminance of the dim target
r =cal.gammaTable(R+1,1);
g =cal.gammaTable(G+1,2);
b =cal.gammaTable(B+1,3);
dimlum = vl*(p*[r;g;b])

%%
% U-probe analyses LFP
% Single trial browser
whichtrial = 60

LADchan = strncmp(stro.sum.rasterCells,'AD',2);
fpoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fpoff_t'));
stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimoff_t'));
targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targon_t'));
targoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targoff_t'));

sig = [];
for whichchan = find(LADchan)
    sig(find(LADchan) == whichchan,:) = stro.ras{whichtrial,whichchan}';
end
AD_t = stro.ras{whichtrial,end}+[0:2:2*size(sig,2)-1]/5000;
meanLFP = mean(sig);
colors = hot(sum(LADchan));

figure; axes; hold on; set(gca,'Color',[.5 .5 .5])
for i = 1:size(sig,1)
   h = plot(AD_t,sig(i,:)-meanLFP,'k-');
   set(h,'Color',colors(i,:))
   plot([stimon_t(whichtrial) stimoff_t(whichtrial)],[0 0],'b-','LineWidth',3)
   plot([targon_t(whichtrial) targoff_t(whichtrial)],[0 0]+.01,'k-','LineWidth',3)
   plot(fpoff_t(whichtrial),0,'k*','LineWidth',3)
end


%%
% Filtered LFP for one trial

if (~rem(size(sig,2),2)) % if numel sig is even
    sig(:,end) = [];
    AD_t(end) = [];
end
sigdur_s = size(sig,2)/5000;
freqs = linspace(0,2500,size(sig,2)/2);
L = freqs<60 | freqs>120; % frequencies to omit
%if (rem(size(sig,2),2)) % if numel sig is odd
    L = logical([1, L, fliplr(L)]); % '1' assumes that we want to omit DC
%else
%    L = logical([L, fliplr(L)]);    
%end
fftout = fft(sig');
fftout(L,:) = 0;

figure; subplot(2,1,1);
plot(AD_t,real(ifft(fftout)))
subplot(2,1,2);
plot(AD_t,imag(ifft(fftout)))

%%
% Filtered LFP for all trials
% (Mean across electrodes subtracted for plotting)
LADchan = strncmp(stro.sum.rasterCells,'AD',2);
fpoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fpoff_t'));
stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimoff_t'));
targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targon_t'));
targoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targoff_t'));
chancols = find(strncmp(stro.sum.rasterCells,'AD',2));
ntrials = size(stro.trial,1);
types = [];
nsamps = 999;
data =nan*ones(ntrials,nsamps,16);
sampfreqHz = stro.sum.analog.storeRates{whichcontact};
for whichcontact = 1:16
    for whichtrial = 1:ntrials
        whichtrial
        sig = [];
        if isnan(stimon_t(whichtrial)) & isnan(targon_t(whichtrial))
            type = 1;
            sync_t = fpoff_t(whichtrial);
        elseif ~isnan(stimon_t(whichtrial)) & isnan(targon_t(whichtrial))
            type = 2;
            sync_t = stimon_t(whichtrial);
        elseif isnan(stimon_t(whichtrial)) & ~isnan(targon_t(whichtrial))
            type = 3;
            sync_t = targon_t(whichtrial);
        else
            type = 4;
            sync_t = stimon_t(whichtrial);
        end
        types(whichtrial) = type;
        
        sig = stro.ras{whichtrial,chancols(whichcontact)}';
        AD_t = stro.ras{whichtrial,end}+[0:2:2*size(sig,2)-1]/sampfreqHz;
        startidx = find(AD_t-sync_t<0,1,'last');
        endidx = startidx+nsamps-1;
        
        sig = sig(startidx:endidx);
        AD_t = AD_t(startidx:endidx);
        
        if (~rem(size(sig,2),2)) % if numel sig is even
            sig(:,end) = [];
            AD_t(end) = [];
        end
        
        sigdur_s = size(sig,2)/sampfreqHz;
        freqs = linspace(0,sampfreqHz/2,size(sig,2)/2);
        L = freqs<0; % frequencies to omit
        L = logical([0, L, fliplr(L)]); % assuming we want to omit the DC
        
        fftout = fft(sig');
        fftout(L) = 0;
        filteredsig = real(ifft(fftout));
        
        data(whichtrial,:,whichcontact) = filteredsig;
    end
end

mn = mean(data,3);

for i = 1:16
    figure;
    for j = 1:4
        subplot(2,2,j);
        clow = min(min(data(types == j,:,i)-mn(types == j,:)));
        chigh = max(max(data(types == j,:,i)-mn(types == j,:)));
        imagesc(data(types == j,:,i)-mn(types == j,:),[clow chigh]);
    end
end

%%
% Looking one electrode at a time
% trial types: (1)blank, (2)stim, (3)targ, (4)stim+targ
% Syncing to (1)fpoff, (2)stimon, (3)targon, (4)stimon 
% This didn't show much

whichcontact = 8;
chancols = find(strncmp(stro.sum.rasterCells,'AD',2));

types = [];
nsamps = 1000;
data = zeros(length(stimon_t),nsamps)
for i = 1:length(stimon_t)
    if isnan(stimon_t(i)) & isnan(targon_t(i))
        type = 1;
        sync_t = fpoff_t(i);
    elseif ~isnan(stimon_t(i)) & isnan(targon_t(i))
        type = 2;
        sync_t = stimon_t(i);
    elseif isnan(stimon_t(i)) & ~isnan(targon_t(i))
        type = 3;
        sync_t = targon_t(i);
    else
        type = 4;
        sync_t = stimon_t(i);
    end
    types(i) = type;
    sig = stro.ras{i,chancols(whichcontact)}';
    t = stro.ras{i,end}+[0:2:2*size(sig,2)-1]/stro.sum.analog.storeRates{whichcontact};
    startidx = find(t-sync_t <0,1,'last');
    endidx = startidx+nsamps-1;
    data(i,:) = sig(startidx:endidx);
end
figure;
axes;
L = types == 2;
plot(data(L,:)');


%%
% Analyses of the cells in Merhdad's figures 
% Kali: K102711014 & K102711015
% K102711016 is opposite side target control
% Sedna: S081111007 & S081111009
% [Sedna: S072712001 Sort of a negative result as is S072911003]

if (isempty(stro))
    stro = strocat(nex2stro,nex2stro);
end
fpoff_t = stro.trial(:,3);
EP_t = [stro.ras{:,end}]';
ADfreq = stro.sum.analog.storeRates{1};
sacstats = getSacData(stro);
for i = 0:1
    for j = [0 2]
        L = stro.trial(:,8) == i & stro.trial(:,9) == j;
        figure; axes; hold on;
        title(['i = ',num2str(i),' j = ',num2str(j)]);
        sum(L)
        for k = find(L)'
            H = stro.ras{k,3};
            V = stro.ras{k,4};
            t = [0:length(H)-1]*(1/ADfreq)+EP_t(k);
            
            sacidx = find(sacstats.starttimes{k} > fpoff_t(k),1,'first')
            amp = sacstats.amplitudes{k}(sacidx);
            dir = sacstats.directions{k}(sacidx);
            [sacx,sacy] = pol2cart(dir,amp);
            if (sacstats.starttimes{k}(sacidx)<fpoff_t(k)+.35)  % latency criterion
                h = plot([0 sacx],[0 sacy],'k-');
                set(h,'Color',[.5 .5 .5]);
                h = plot(sacx,sacy,'ko');
                set(h,'MarkerFaceColor','black','MarkerSize',5,'MarkerEdgeColor',[.5 .5 .5]);
            end
        end
        set(gca,'Xlim',[-10 10],'Ylim',[-10 10]);
        axis square;
    end
end

%}