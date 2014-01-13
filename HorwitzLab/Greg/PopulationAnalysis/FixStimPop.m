% FixStim population analyses

% Section 1: Looking at fixation around the time of optical stimulation
%
% Section 2: Tabulating the magnitude of the optical stimulation effect on neural activity and RF location. 
%
% Section 3: Looking for relatiohships between microsaccades and opto-responses in SC

%%
% Section 1: Looking at fixation around the time of optical stimulation

[fnames, spikeIdx] = fnamesFromTxt2();
offset = [-.2 .2];  % pre and post time
H_data = [];
V_data = [];
infomat = [];

for filecounter = 1:size(fnames,1)
    filename = findfile(char(fnames{filecounter}));
    paradigmID = getparadigmID(filename);
    if paradigmID ~= 105
        error([fnames{counter},' is not a FixStim file']);
    else
        stro = nex2stro(filename);
        ADfreq = stro.sum.analog.storeRates{1};
        sampleperiod = 1/ADfreq;
        H = stro.ras(:,strcmp(stro.sum.rasterCells,'AD11'));
        V = stro.ras(:,strcmp(stro.sum.rasterCells,'AD12'));
        ep_t = stro.ras(:,strcmp(stro.sum.rasterCells,'anlgStartTime'));
        stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
        fpoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fpoff_t'));
        % Finding time that stimulation *would have* come on in a no-stim
        % trial. Fix point goes off first. 
        targ_shown = logical(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_shown')));
        avgtargstimon_t = nanmean(stimon_t(targ_shown)-fpoff_t(targ_shown));
        avgnotargstimon_t = nanmean(stimon_t(~targ_shown)-fpoff_t(~targ_shown));
        sync_t = stimon_t;
        sync_t(targ_shown & isnan(stimon_t)) = avgtargstimon_t+fpoff_t(targ_shown & isnan(stimon_t));
        sync_t(~targ_shown & isnan(stimon_t)) = avgnotargstimon_t+fpoff_t(~targ_shown & isnan(stimon_t));

        for cond = 1:4
            if (cond == 1)
                L = ~targ_shown & isnan(stimon_t); % No targ, no stim 
            elseif (cond == 2)
                L = ~targ_shown & ~isnan(stimon_t); % No targ, stim
            elseif (cond == 3)
                L = targ_shown & isnan(stimon_t); % Targ, no stim 
            elseif (cond == 4)
                L = targ_shown & ~isnan(stimon_t); % Targ, stim 
            end
            trlidxs = find(L);
            hs = []; vs = [];
            for counter = 1:sum(L)
                trlidx = trlidxs(counter);
                hep = H{trlidx}*4096/400;
                vep = V{trlidx}*4096/400;
                t = [0:length(hep)-1]*sampleperiod+ep_t{trlidx};
                t = t-sync_t(trlidx);
                L_t = t > offset(1)-sampleperiod/2 & t < offset(2)+sampleperiod/2;
                % hack
                if (~isempty(H_data))
                    while (sum(L_t) ~= size(H_data,2))
                        if sum(L_t) > size(H_data,2)
                            L_t(find(L_t,1,'last')) = 0;
                        elseif sum(L_t) < size(H_data,2)
                            L_t(find(L_t,1,'last')+1) = 1;
                        end
                    end
                end
                hs(counter,:) = hep(L_t)';
                vs(counter,:) = vep(L_t)';
            end
            H_data = [H_data; hs];
            V_data = [V_data; vs];
            infomat = [infomat; repmat([cond filecounter],sum(L),1)];
        end
    end
end

figure;
plot_t = linspace(offset(1),offset(2),size(H_data,2));
for cond = 1:4
    subplot(2,2,cond); hold on;
    L = infomat(:,1) == cond;
    plot(plot_t,mean(H_data(L,:)),'r-','Linewidth',2)
    plot(plot_t,std(H_data(L,:))'*[-1 1]+repmat(mean(H_data(L,:))',1,2),'r-')
    
    plot(plot_t,nanmean(V_data(L,:)),'b-','Linewidth',2)
    plot(plot_t,std(V_data(L,:))'*[-1 1]+repmat(mean(V_data(L,:))',1,2),'b-')
end

figure;
for cond = 1:4
    subplot(2,2,cond); hold on;
    L = infomat(:,1) == cond;
    imagesc(H_data(L,:))
    axis tight;
end
figure;
for cond = 1:4
    subplot(2,2,cond); hold on;
    L = infomat(:,1) == cond;
    imagesc(V_data(L,:));
    axis tight;
end


%%
% Section 2
% Looking at optogenetic neural activation as a function of RF position
% FixStimPop is a good list.

[fnames, spikeIdx] = fnamesFromTxt2();
data = [];
for filecounter = 1:size(fnames,1)
    filename = findfile(char(fnames{filecounter}));
    paradigmID = getparadigmID(filename);
    if paradigmID ~= 105
        error([fnames{counter},' is not a FixStim file']);
    else
        stro = nex2stro(filename);
        rfx = stro.sum.exptParams.rf_x;
        rfy = stro.sum.exptParams.rf_y;
        ras = stro.ras(:,strcmp(stro.sum.rasterCells, 'sig001a'));
        stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
        stimtrials = find(~isnan(stimon_t))';
        tmp = zeros(length(stimtrials),2);
        for i = stimtrials
            spiketimes = ras{i}-stimon_t(i);
            tmp(i == stimtrials,:) = [sum(spiketimes > -.01 & spiketimes < 0) sum(spiketimes > 0 & spiketimes < .01)];
        end
        data = [data; rfx rfy mean(tmp,1) filecounter];
    end
end

L = data(:,1) > 0 & data(:,2) < 0; % Only targets in the appropriate quadrant
data(~L,:) = [];
L = data(:,3) > 0 & data(:,4) > 0;
data(~L,:) = [];

spikecount_ratio = data(:,4)./data(:,3);
if (isnan(spikecount_ratio) | spikecount_ratio == 0)
    spikecount_ratio = eps;
elseif (isinf(spikecount_ratio))
    spikecount_ratio = 1;
end

% Plotting a symbol for every experiment
figure; axes; hold on;
unique_rfxys = unique(data(:,[1 2]),'rows');
for i = 1:size(data,1)
 %   h = plot(data(i,1),data(i,2),'ko','MarkerSize',10*spikecount_ratio(i),'MarkerFaceColor',[spikecount_ratio(i)/nanmax(spikecount_ratio) 0 0]);
 %   h = plot(data(i,1),data(i,2),'ko','MarkerSize',10*spikecount_ratio(i),'MarkerEdgeColor',[spikecount_ratio(i)/nanmax(spikecount_ratio) 0 0],'LineWidth',2);
    h = plot(data(i,1)/10,data(i,2)/10,'k.','MarkerSize',10*spikecount_ratio(i),'MarkerEdgeColor',[spikecount_ratio(i)/max(spikecount_ratio) 0 0]);    
end
set(gca,'Xlim',[0 8],'Ylim',[-8 0])
axis equal;
disp(['Max Spike count ratio is ',num2str(max(spikecount_ratio))]);

% Plotting one symbol for each RF location
figure; axes; hold on;
unique_rfxys = unique(data(:,[1 2]),'rows');
for i = 1:size(unique_rfxys,1)
    L = data(:,1) == unique_rfxys(i,1) & data(:,2) == unique_rfxys(i,2);
    scr = max(spikecount_ratio(L));
    h = plot(unique_rfxys(i,1)/10,unique_rfxys(i,2)/10,'k.','MarkerSize',10*scr,'MarkerEdgeColor',[scr/max(spikecount_ratio) 0 0]);    
end
set(gca,'Xlim',[0 8],'Ylim',[-8 0])
axis equal;

% A key
figure; axes; hold on;
set(gca,'Xlim',[0 8],'Ylim',[-8 0])
axis equal;
scr_steps = [1 2 5 max(spikecount_ratio)];
for i = 1:length(scr_steps)
    i
    h = plot(i,-4,'k.','MarkerSize',10*scr_steps(i),'MarkerEdgeColor',[scr_steps(i)/max(spikecount_ratio) 0 0]);
    text(i,-3,num2str(scr_steps(i)));
end





%%
% Section 3
% Looking for relatiohships between microsaccades and opto-responses in SC
INCLUDETARGTRIALS = 1;
[fnames, spikeIdx] = fnamesFromTxt2();
data = [];
for filecounter = 1:size(fnames,1)
    filename = findfile(char(fnames{filecounter}));
    paradigmID = getparadigmID(filename);
    if paradigmID ~= 105
        error([fnames{counter},' is not a FixStim file']);
    else
        stro = nex2stro(filename);
        sacstats = getSacData(stro); close;
        stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
        targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targon_t'));
        L = ~isnan(stimon_t);
        if (INCLUDETARGTRIALS)
            L = L & isnan(targon_t);
        end
        stimtrials = find(L)';
        tmp = [];
        for i = stimtrials
            spiketimes = stro.ras{i,1}-stimon_t(i);
            nspikes1 = sum(spiketimes > -.01 & spiketimes < 0);
            nspikes2 = sum(spiketimes > 0 & spiketimes < .02);
            sactimes = sacstats.starttimes{i}-stimon_t(i);
            lastsactime = max(sactimes(sactimes < 0));
            if (isempty(lastsactime))
                lastsactime = nan;
            end
            tmp = [tmp; lastsactime nspikes1 nspikes2];
        end
        tmp(isnan(tmp(:,1)),:) = [];
        p = signrank(tmp(:,2), tmp(:,3));
        if p < 0.001
            %tmp(tmp(:,1) < -0.5,:) = []; % getting rid of saccades that are too long ago
            plot(tmp(:,1),tmp(:,2),'k.'); lsline;
            [r,p] = corr(tmp,'type','spearman');
            title(['r = ',num2str(r(1,2)),' p = ',num2str(p(1,2))]);
            data = [data; r(1,2) p(1,2) filecounter];
            stro.sum.fileName
            keyboard
        end
    end
end
hist(data(:,1));

