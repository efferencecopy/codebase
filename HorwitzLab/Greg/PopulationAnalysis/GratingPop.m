% GratingPop.m
%
% Contents:
% ---------
% Section 1)
% Correlation analyses on the ranked responses across color directions.
% This analysis will tell us whether cells that respond well to L+M (or L-M)
% also respond well to the intermediate color direction that drives the
% greater response (L-M-S).  It will also be useful for seeing what "types"
% of color tuning we find with the grating paradigm.
%
% Section 2) 
% Getting F1/F0 modulation ratios for a bunch of cells.
%
% Section 3)
% Taking a look at F1/F0 modulation in different color directions.
%
% Section 3.1)
% Looking at F1/F0 modulation in different color directions like section 3,
% above, but now calculating a "bootstrap F-statistic" to quantify now much
% the F1/F0 ratio varies with color direction.
%
% Section 4) 
% Checking to see whether the algorithm I'm using to find preferred color
% direction from gratings data depends on initial guess. (it does).
%
% Section 5) 
% Finding double-opponent cells based on the relative phase of response to
% L- and M-cone isolating stimuli.
%
% Section 5.1)
% Finding double-opponent cells based on criteria similar to those of
% Johnson and Shapley.
%
% Section 6)
% Collecting orientation tuning curves and receptive field locations for
% Jochen Triech
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% SECTION 1
% Correlation analyses on the ranked responses across color directions.
% This analysis will tell us whether cell that respond well to L+M (or L-M)
% also respond well to the intermediate color direction that drives the
% greater response (L-M-S).  It will also be useful for seeing what "types"
% of color tuning we find with the grating paradigm.
%%%%%%%%%%%%%%%%%%%%%%%%%%

data = {};
[fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\Gratings.txt');
for a = 1:size(fnames,1)
    stro = {};
    for i = 1:size(fnames{a},2)
        tmpstro = nex2stro(findfile(char(fnames{a}(i))));
        if (isempty(stro))
            stro = tmpstro;
        else
            stro = strocat(stro, tmpstro);
        end
    end

    out = getGratingTuning(stro,spikeIdx(a));
    data{a} = out;
end

% % Lennie plot
% tmp = [];
% for i = 1:length(data)
%     tmp = [tmp; data{i}.prefSF data{i}.preforient data{i}.prefcolor];
% end
% 
% tmp(tmp(:,1) > 2,:) = [];
% 
% Lscone = tmp(:,5) > 0;
% figure; axes; hold on;
% plot(tmp(Lscone,3),tmp(Lscone,4),'k.');
% plot(tmp(~Lscone,3),tmp(~Lscone,4),'ko');
% plot(-tmp(~Lscone,3),-tmp(~Lscone,4),'k.');
% plot(-tmp(Lscone,3),-tmp(Lscone,4),'ko');
% plot([0 1 0 -1 0],[1 0 -1 0 1],'k-');
% plot([0 0],[-1 1],'k:');
% plot([-1 1],[0 0],'k:');
% axis equal;

% Correlation analysis of color responses
rawresponses = [];
rankedresponses = [];
prefsfs = [];
for i = 1:length(data)
    if all(all(data{i}.color.colors == data{end}.color.colors))
        prefsfs = [prefsfs; data{i}.sf.prefSF];
        rawresponses = [rawresponses; data{i}.color.colresp(:,1)'];
        rankedresponses = [rankedresponses; tiedrank(data{i}.color.colresp(:,1))'];
    end
end
% 
% L = prefsfs < prctile(prefsfs,33);
% % L = prefsfs < prctile(prefsfs,66) & prefsfs >= prctile(prefsfs,33);
% % L = prefsfs >= prctile(prefsfs,66);
%  rankedresponses = rankedresponses(L,:);

figure;
subplot(3,1,1);
for i = 1:size(data{end}.color.colors,1)
    text(.1,1-(i/(size(data{end}.color.colors,1)+1)), [num2str(i),') ',num2str(data{end}.color.colors(i,:))]);
end
set(gca,'XTick',[],'YTick',[]);

subplot(3,1,2); 
[rho, p] = corr(rankedresponses,'type','spearman');
sigs = p < 0.01.*(1-eye(9));
imagesc(rho); hold on;
[i,j] = ind2sub(size(sigs),find(sigs));
plot(i,j,'y*');
axis square; colormap(gray);
subplot(3,1,3);
plot(mean(rankedresponses),'k.-')
ylabel('mean rank'); xlabel('color');

% Do cells that respond well to L+M exhibit a large *differential*
% response to the two intermediate color directions?
figure; axes; hold on;
[n,x] = hist2([rankedresponses(:,8),rankedresponses(:,3)-rankedresponses(:,4)],9)
imagesc(n); colormap(gray);
xlabel('response to L+M');
ylabel('intermediate response asymmetry');
corr(rankedresponses(:,8),rankedresponses(:,3)-rankedresponses(:,4));
exgiveny = linspace(1,9,9)*(n./repmat(sum(n),9,1))
plot(exgiveny,'y-'); axis image;
% Asymmetry in responses to intermediated color directions does not appear 
% to be strongly related to luminance tuning.


[u,s,v] = svd(rankedresponses-repmat(mean(rankedresponses),size(rankedresponses,1),1));
%[u,s,v] = svd(rankedresponses);
figure;
subplot(3,1,1);
plot(diag(s)./sum(diag(s)),'k.');
subplot(3,1,2); hold on;
plot(-v(:,1),'k.-');
plot(-v(:,2),'m.-');
subplot(3,1,3); hold on;
plot(-u(:,1),-u(:,2),'k.');
% The 1st PC is the mean (sort of indicating pancolor?) and the 2nd is a
% color vs/luminance axis.  Scatterplot on 1st and 2nd is oriented.  
% What does this mean?

% NOTICED:
% Nice anticorrelation between L-M and L+M responses, as hoped.
% The stronger L-M is, the weaker L+M-S
% Cells that like L tend not to like S
% Cells that like L-M tend not to like S-(L+M)
% Cells that like L+M tend not to like S

% Pairwise scatterplots.
for i =1:9
    for j = 1:i-1
        figure; axes; hold on;
        [n,x] = hist2([rankedresponses(:,i),rankedresponses(:,j)],9);
        imagesc(n);
        axis ij;
        exgiveny = linspace(1,9,9)*(n./repmat(sum(n),9,1));
        plot(exgiveny,'y-');
        colormap(gray);
        axis tight;
        xlabel(num2str(data{end}.color.colors(i,:)));
        ylabel(num2str(data{end}.color.colors(j,:)));
        set(gca,'YTickLabel', flipud(get(gca,'YTickLabel')));
        r = corrcoef([rankedresponses(:,i),rankedresponses(:,j)]);
        title(['r = ',num2str(r(1,2))]);
    end
end

% Does the bias in the preference for the L-M-S intermediate color
% direction over the L-M+S intermediate color direction depend on spatial
% frequency (like psychophysical detection does?)
MwithSidx = 3;
MoppSidx = 4;
%(A-B)/(A+B)
intermedbias = (rawresponses(:,MwithSidx)-rawresponses(:,MoppSidx))./(rawresponses(:,MwithSidx)+rawresponses(:,MoppSidx));
figure; subplot(1,3,1);
plot(prefsfs, intermedbias,'k.');
lsline;
xlabel('pref sf'); ylabel('bias for L-M-S');
% Marginal hist
subplot(1,3,2);
hist(intermedbias,30);
mn = nanmean(intermedbias);
[h,p] = ttest(intermedbias);
title(['grand mean: ',num2str(mn),' p: ',num2str(p)]);
% Just looking at the 0.5 cpd neurons
L = prefsfs == 0.5;
subplot(1,3,3);
hist(intermedbias(L),30);
mn = nanmean(intermedbias(L));
[h,p] = ttest(intermedbias(L));
title(['0.5 cpd mean: ',num2str(mn),' p: ',num2str(p)]);





%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% SECTION 2
% Looking at a histogram of modulation ratios
%%%%%%%%%%%%%%%%%%%%%%%%%%

data = {};
[fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\Gratings.txt');
for a = 1:size(fnames,1)
    stro = {};
    for i = 1:size(fnames{a},2)
        tmpstro = nex2stro(findfile(char(fnames{a}(i))));
        if (isempty(stro))
            stro = tmpstro;
        else
            stro = strocat(stro, tmpstro);
        end
    end

    out = getGratingTuning(stro,spikeIdx(a));
    data{a} = out;
end

x = [];
for i = 1:length(data)
    x = [x; data{i}.modulationratio];
end
figure;
hist(x,20)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% SECTION 3
% Looking at modulation ratios in a different color directions.
%%%%%%%%%%%%%%%%%%%%%%%%%%
data = {};
[fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\Gratings.txt');
for a = 1:size(fnames,1)
    stro = {};
    for i = 1:size(fnames{a},2)
        tmpstro = nex2stro(findfile(char(fnames{a}(i))));
        if (isempty(stro))
            stro = tmpstro;
        else
            stro = strocat(stro, tmpstro);
        end
    end
    
    framerate = stro.sum.exptParams.framerate;
    nframes = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'nframes'));
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_on'));
    stimoff_t= stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_off'));
    tfs = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'tf'));
    protocols = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'protocol'));
    Lcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'lcont'));
    Mcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'mcont'));
    Scc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'scont'));
    
    L4 = protocols == 4;
    baselines = [];  baseline_t = 0.5;
    for i = find(L4)'
        spiketimes = stro.ras{i,spikeIdx(a)};
        nspikes = sum(spiketimes > stimon_t(i)-baseline_t & spiketimes < stimon_t(i));
        baselines = [baselines; nspikes./baseline_t]; % spikes/sec
    end
    
    if (sum(L4)<3*9)  % Too few (or no) protocol 4 trials.
        continue
    end
    minstimdur = min(stimoff_t(L4)-stimon_t(L4))
    if (minstimdur < .9)  % Stimulus presentation too short. Must be at least 900 ms
        continue
    end
    
    binwidth = .001; % Does this matter?  Principled way of choosing this?
    bins = [0:binwidth:minstimdur];  % Bad assuming stimulus dur
    colordirections = [Lcc(L4) Mcc(L4) Scc(L4)];
    uniquecolordirections = unique(colordirections,'rows');
    tmp = [];
    for i = 1:size(uniquecolordirections,1)
        L = L4 &...
            Lcc == uniquecolordirections(i,1) & ...
            Mcc == uniquecolordirections(i,2) &...
            Scc == uniquecolordirections(i,3);
        start_t = stimon_t(L);
        spiketimes = stro.ras(L,spikeIdx(a));
        totalspikes = [];
        for j = 1:sum(L)
            totalspikes = [totalspikes; spiketimes{j}-start_t(j)];
        end
        psth = histc(totalspikes',bins)/binwidth/sum(L);
        k = unique(stro.trial(L,strcmp(stro.sum.trialFields(1,:), 'tf')));
        if (length(k) > 1)
            error('Too many temporal frequencies!');
        end
        if (k ~= 3)
            error('TF not 3 Hz');
        end
        basis0 = ones(length(psth),1);
        basis1 = exp(-2*pi*sqrt(-1)*k*[0:length(psth)-1]/length(psth))';
        F0 = abs((psth-mean(baselines))*basis0);
        F1 = 2*abs((psth-mean(baselines))*basis1);
        % Not sure what the logic for the the "2*" in the line above
        % but it makes a half-wave rectified sinewave have a modulation
        % ratio of pi/2, which is what it's supposed to have.
        tmp = [tmp; F1 F0];
    end
    idx = length(data)+1;
    data(idx).colors = uniquecolordirections;
    data(idx).F1 = tmp(:,1);
    data(idx).F0 = tmp(:,2);
    data(idx).filename = fnames{a};
    data(idx).meanfr = tmp(:,2)/length(psth);
    data(idx).baseline = mean(baselines); 
end
%%
% Now unpacking this data stucture into something useable

COLORDIRS = data(end).colors;
out = nan*ones(9,length(data));
for a = 1:length(data)
    if(isempty(data(a)))
        continue;
    end
    L = [ismember(data(a).colors,-COLORDIRS,'rows'),...
        ismember(data(a).colors,COLORDIRS,'rows')];
    if all(any(L,2))
        out(:,a) =  data(a).F1./data(a).F0;
        out(data(a).meanfr < 10, a) = nan;  % 10 sp/sec or more
    else
        disp('Wrong colors');
    end
end

figure;
subplot(2,1,1);
plot(nanmean(out,2))
ylabel('F1/F0'); xlabel('Color direction');
subplot(2,1,2);
plot(nanstd(out),'.')
ylabel('std(F1/F0)'); xlabel('Cell index');

% This worked better than standard deviation or resampling statistic
figure;
plot(nanmax(out)./nanmin(out),'k.');
ylabel('fold change in F1/F0'); xlabel('Cell index');
set(gca,'Yscale','log');

%%
% Section 3.1
% Quantifying the ratio of between color variance in F1/F0 to the within
% color variance of F1/F0.  Using resampling.  This didn't work that well.
data = [];
[fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\Gratings.txt');
for a = 1:size(fnames,1)
    stro = {};
    for i = 1:size(fnames{a},2)
        tmpstro = nex2stro(findfile(char(fnames{a}(i))));
        if (isempty(stro))
            stro = tmpstro;
        else
            stro = strocat(stro, tmpstro);
        end
    end
    
    framerate = stro.sum.exptParams.framerate;
    nframes = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'nframes'));
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_on'));
    stimoff_t= stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_off'));
    tfs = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'tf'));
    protocols = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'protocol'));
    Lcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'lcont'));
    Mcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'mcont'));
    Scc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'scont'));
    
    L4 = protocols == 4;
    baselines = [];  baseline_t = 0.5;
    for i = find(L4)'
        spiketimes = stro.ras{i,spikeIdx(a)};
        nspikes = sum(spiketimes > stimon_t(i)-baseline_t & spiketimes < stimon_t(i));
        baselines = [baselines; nspikes./baseline_t]; % spikes/sec
    end
    
    if (sum(L4)<3*9)  % Too few (or no) protocol 4 trials.
        continue
    end
    minstimdur = min(stimoff_t(L4)-stimon_t(L4));
    if (minstimdur < .9)  % Stimulus presentation too short. Must be at least 900 ms
        continue
    end
    binwidth = .001; % Does this matter?  Principled way of choosing this?
    bins = [0:binwidth:minstimdur];  % Bad assuming stimulus dur
    colordirections = [Lcc(L4) Mcc(L4) Scc(L4)];
    uniquecolordirections = unique(colordirections,'rows');

    meanfiringrate = []; fftbycolor = {};
    for i = 1:size(uniquecolordirections,1)
        L = L4 &...
            Lcc == uniquecolordirections(i,1) & ...
            Mcc == uniquecolordirections(i,2) &...
            Scc == uniquecolordirections(i,3);
        start_t = stimon_t(L);
        spiketimes = stro.ras(L,spikeIdx(a));
        k = unique(stro.trial(L,strcmp(stro.sum.trialFields(1,:), 'tf')));
        if (length(k) > 1)
            error('Too many temporal frequencies!');
        end
        if (k ~= 3)
            error('TF not 3 Hz');
        end
        basis0 = ones(length(bins),1);
        basis1 = exp(-2*pi*sqrt(-1)*k*[0:length(bins)-1]/length(bins))';
        % Yes, I can get the same F1 and F0 doing the calculation
        % on a trial-to-trial basis.  Hooray for linearity.
        psth = [];
        for j = 1:sum(L)
            psth(j,:) = histc(spiketimes{j}-start_t(j),bins)/binwidth;
            if (size(psth,1) > size(psth,2))  % Annoying transpose issue
                psth = psth';
            end
            F0(j) = (psth(j,:)-mean(baselines))*basis0;
            F1(j) = 2*(psth(j,:)-mean(baselines))*basis1;
        end
        fftbycolor{i} = [F1' F0'];
        meanfiringrate(i) = sum(psth(:))./length(psth)./sum(L);
    end
    
    
    % Bootstrapping to estimate between-color-direction and
    % within-color-direction variance in F1/F0
    % Generating resamples of data
    niter = 200;
    tmp = zeros(niter,length(fftbycolor),2);
    for i = 1:niter
        for j = 1:length(fftbycolor)
            n = size(fftbycolor{j},1);
            tmp(i,j,:) = abs(sum(fftbycolor{j}(unidrnd(n,n,1),:)));
        end
    end
    % Computing statistics from resampled data
    bootdata = [];
    for i = 1:length(fftbycolor)
        modratios = tmp(:,i,1)./tmp(:,i,2);
        Q1 = mean(tmp(:,i,1)./tmp(:,i,2));
        Q2 = var(tmp(:,i,1)./tmp(:,i,2));
        bootdata(i,:) = [Q1 Q2];
    end
    % Bootdata: first column is mean, second is 
    % Figuring out which color directions have mean firing rates > 10 sp/s.
    
    L =  meanfiringrate > 10;
    if (sum(L) > 2)
        data(a,:) = [var(bootdata(L,1)) mean(bootdata(L,2)) sum(L)]
    end
    % First column is between var, second column is within var
end

figure
plot(data(:,1)./data(:,2),'k.')
% High numbers are potentially interesting cells.

%%
% Section 4
% Checking to see whether algorthim for finding preferred color
% directions from gratings data is robust to the initial guess.
% Answer: not very robust.
niter = 100;  % initial guesses per cell
data = [];
[fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\Gratings.txt');
for a = 1:size(fnames,1)
    stro = nex2stro(findfile(char(fnames{a})));
    protocols = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'protocol'));
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_on'));
    stimoff_t= stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_off'));
    spikeidx = strcmp(stro.sum.rasterCells(1,:),getSpikenum(stro));

    spikerates = [];
    for i = 1:size(stro.trial,1)
        spiketimes = stro.ras{i,spikeidx};
        nspikes = sum(spiketimes > stimon_t(i) & spiketimes < stimoff_t(i));
        spikerates = [spikerates; nspikes./(stimoff_t(i)-stimon_t(i))];
    end
    Lprotocol = protocols == 4;
    if (any(Lprotocol))
        tmp = nan(niter,3);
        for i = 1:niter
            Lcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'lcont'));
            Mcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'mcont'));
            Scc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'scont'));
            
            colordirections = [Lcc(Lprotocol) Mcc(Lprotocol) Scc(Lprotocol)];
            uniquecolordirs = unique(colordirections,'rows');
         %   responses = spikerates(Lprotocol);
         %   maxrespidx = find(responses == max(responses),1);
            
            initguess = normrnd(0,1,1,3);
          %  initguess = initguess*(responses(maxrespidx)./norm(initguess));
             initguess = initguess*(mean(responses)./norm(initguess));
            
            [coneweights, fval, exitflag] = fminsearch(@(x) linmodfiterr(colordirections, responses, x), initguess);
            if (~exitflag)
                continue
            end
            tmp(i,:) = coneweights;
        end
        b = unique(round(tmp),'rows')
        n = zeros(size(b,1),1);
        for i = 1:size(b,1)
            n(i) = sum(all(round(tmp) == repmat(b(i,:),niter,1),2))
        end
        coneweights = b(find(n== max(n),1),:);
        Lconsistent = all(b == repmat(coneweights,size(b,1),1),2) | all(b == -repmat(coneweights,size(b,1),1),2);
        nerrs = sum(n(~Lconsistent));
    end
    data = [data; a nerrs coneweights]
end

%%
% Section 5
% Finding double-opponent cells from grating data.
% Based on the relative phase of responses to L- and M-cone isolating stimuli 

[fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\Gratings.txt');
data = zeros(size(fnames,1),6);
for a = 1:size(fnames,1)
    stro = nex2stro(findfile(char(fnames{a})));
    protocols = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'protocol'));
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_on'));
    stimoff_t= stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_off'));
    spikeidx = strcmp(stro.sum.rasterCells(1,:),getSpikenum(stro,'first'));
    sfs = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'sf'));

    spikerates = [];
    for i = 1:size(stro.trial,1)
        spiketimes = stro.ras{i,spikeidx};
        nspikes = sum(spiketimes > stimon_t(i) & spiketimes < stimoff_t(i));
        spikerates = [spikerates; nspikes./(stimoff_t(i)-stimon_t(i))];
    end
    L4 = protocols == 4;
    if sum(L4) >= 9*3
        Lcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'lcont'));
        Mcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'mcont'));
        Scc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'scont'));
        colordirections = [Lcc(L4) Mcc(L4) Scc(L4)];
        uniquecolordirs = unique(colordirections,'rows');
       
        baselines = [];  baseline_t = 0.5;
        for i = find(L4)'
            spiketimes = stro.ras{i,spikeIdx(a)};
            nspikes = sum(spiketimes > stimon_t(i)-baseline_t & spiketimes < stimon_t(i));
            baselines = [baselines; nspikes./baseline_t]; % spikes/sec
        end
        
        minstimdur = min(stimoff_t(L4)-stimon_t(L4));
        if (minstimdur < .9)  % Stimulus presentation too short. Must be at least 900 ms
            continue
        end
        
        binwidth = .001; % Does this matter?  Principled way of choosing this?
        bins = [0:binwidth:minstimdur];  % Bad assuming stimulus dur
        colordirections = [Lcc(L4) Mcc(L4) Scc(L4)];
        uniquecolordirections = unique(colordirections,'rows');
        tmp = [];
        for i = 1:size(uniquecolordirections,1)
            L = L4 &...
                Lcc == uniquecolordirections(i,1) & ...
                Mcc == uniquecolordirections(i,2) &...
                Scc == uniquecolordirections(i,3);
            start_t = stimon_t(L);
            spiketimes = stro.ras(L,spikeIdx(a));
            totalspikes = [];
            for j = 1:sum(L)
                totalspikes = [totalspikes; spiketimes{j}-start_t(j)];
            end
            psth = histc(totalspikes',bins)/binwidth/sum(L);
            k = unique(stro.trial(L,strcmp(stro.sum.trialFields(1,:), 'tf')));
            if (length(k) > 1)
                error('Too many temporal frequencies!');
            end
            if (k ~= 3)
                error('TF not 3 Hz');
            end
            basis0 = ones(length(psth),1);
            basis1 = exp(-2*pi*sqrt(-1)*k*[0:length(psth)-1]/length(psth))';
            F0 = abs((psth-mean(baselines))*basis0);
            F1 = 2*abs((psth-mean(baselines))*basis1);
            theta = angle((psth-mean(baselines))*basis1);
            tmp = [tmp; F1 F0 theta];
        end        
        % Pull out the L and M-cone isolating stimuli
        Ll = uniquecolordirections(:,1) ~= 0 & uniquecolordirections(:,2) == 0 & uniquecolordirections(:,3) == 0;
        Lm = uniquecolordirections(:,1) == 0 & uniquecolordirections(:,2) ~= 0 & uniquecolordirections(:,3) == 0;
        frl = tmp(Ll,2)/length(psth);
        frm = tmp(Lm,2)/length(psth);
        if (frl > 5) & (frm > 5)
            data(a,1:6) = [tmp(Ll,1)./tmp(Ll,2) frl tmp(Ll,3) tmp(Lm,1)./tmp(Lm,2) frm tmp(Lm,3)]
        end
        
        % ------------------------------------
        % Pulling out info about bandpass-ness
        % ------------------------------------
        if (~any(protocols == 2))
             data(a,:) = 0;
             continue;
        end
        lastorienttrial = find(protocols == 1,1,'last');
        if (all(protocols(1:lastorienttrial) == 1))
            firstorienttrial = 1;
        else
            firstorienttrial = find(protocols(1:lastorienttrial) ~= 1,1,'last')+1;
        end
        lastSFtrial = find(protocols == 2,1,'last');
        if (all(protocols(1:lastSFtrial) == 1))
            firstSFtrial = 1;
        else
            firstSFtrial = find(protocols(1:lastSFtrial) ~= 2,1,'last')+1;
        end
        trlidxs = [firstSFtrial:lastSFtrial];
        x = sfs(trlidxs);
        y = spikerates(trlidxs);
        tmp = [];
        for i = sort(unique(x)')
            tmp = [tmp; mean(y(x == i))];
        end
        if (max(tmp) == tmp(1))  % if lowest SF produces the largest response
            data(a,:) = 0;
        end
    end
end

figure; axes; hold on;
plot(data(:,1),data(:,4),'k.');
Lmodratio = data(:,1) >= 1 & data(:,4) >= 1;  % Pretty strict.  Might want to examine L-Ms instead?
angulardifference = data(:,3)-data(:,6)
Lphase = abs(angulardifference) > pi/2 & abs(angulardifference) < 3*pi/2;
plot(data(Lphase,1),data(Lphase,4),'m.');

[fnames{Lmodratio & Lphase}]'


%%
% Section 5.1
% Finding double-opponent cells from grating data.
% Based on the modulation ratio across all conditions, color sensitivity,
% and bandpassness

[fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\Gratings.txt');
data = zeros(size(fnames,1), 4);
for a = 1:size(fnames,1)
    stro = nex2stro(findfile(char(fnames{a})));
    GTstruct = getGratingTuning(stro, 1);
    if(~isfield(GTstruct.sf,'stim'))
        Lsf = 0;
    elseif any(GTstruct.sf.resp(2:end,1) >  GTstruct.sf.resp(1,1))
        Lsf = 1;
    else
        Lsf = 0;
    end
    if (~isfield(GTstruct,'color'))
        continue;
    end
    if (isnan(GTstruct.color.colors))
        continue;
    end
    dotprods = abs(GTstruct.color.colors*[1 -1 0]');
    rgidx = dotprods == max(dotprods) & GTstruct.color.colors(:,3) == 0;
    rg = GTstruct.color.colresp(rgidx,1);
    lumidx = all(GTstruct.color.colors == repmat(abs(GTstruct.color.colors(rgidx,:)),size(GTstruct.color.colors,1),1),2);
    lum = GTstruct.color.colresp(lumidx,1);
    data(a,:) = [rg lum Lsf GTstruct.modulationratio];
end

figure; axes; hold on;
L = data(:,1) > data(:,2) & data(:,3) & data(:,4) >= 1;
[fnames{L}]'


%%
% Section 6
% Collecting orientation tuning curves and receptive field locations for
% Jochen Triech
data = {};
[fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\Gratings.txt');
for a = 1:size(fnames,1)
    stro = {};
    for i = 1:size(fnames{a},2)
        tmpstro = nex2stro(findfile(char(fnames{a}(i))));
        if (isempty(stro))
            stro = tmpstro;
        else
            stro = strocat(stro, tmpstro);
        end
    end

    out = getGratingTuning(stro,spikeIdx(a));
    data{a} = out;
    data{a}.rf = [stro.sum.exptParams.rf_x stro.sum.exptParams.rf_y]/10;
end

ForJochen = []
for i = 1:length(data)
    ForJochen(i,:) = [data{i}.rf data{1}.orient.preforient data{i}.orient.stim(1:8)' data{i}.orient.resp((1:8),1)' data{i}.orient.resp((1:8),2)'];
end

plot(ForJochen(:,1),ForJochen(:,2),'k.');

% Is the effect he's looking for there?
angles = [];
for i = 1:length(data)
    angles(i,:) = [atan2(data{i}.rf(2),data{i}.rf(1)) data{i}.orient.preforient];
end
stat = mod(-angles(:,1)-angles(:,2),pi); 
hist(stat)
1-binocdf(sum(stat<pi/2),length(stat),.5)

