% CONTENTS
% --------
%
% Section 1:
% Summary plots for the statistics of fixational saccades.
%
% Section 1.1:
% Microsaccade counts/frequencies during various trial epochs.
%
% Section 1.2:
% Comparison of errors and corrects as a function of target choice.
%
% Section 2:
% Probability of detection (from a list of DTspot files) as a function of
% saccade occurrence (binary), spatial frequency, and color direction.
%
% Section 2.5:
% Correlation between saccade occurrence and correct as a function of
% choice (in a single counting window).
%
% Section 3:
% Sliding window analysis looking at correlation between saccade occurrence
% and correct (or target choice).
%
% Section 4:
% Distribution of saccade directions (relative to correct target)
% as a function of choice or correct/incorrect at different points
% during the trial.
%
% Section 5:
% Distribution of saccade amplitudes as a function of choice or
% correct/incorrect at different points during the trial.
%
% Section 5.1:
% Distribution of saccade endpoints (amplitude and direction) as a function
% of time in the trial.  Maybe this will allow us to identify
% "anticipatory saccades".
%
% Section 5.2:
% Distribution of saccade endpoints (amplitude and direction) on *error*
% trials.  Are saccades directed toward the flash (as predicted by
% microsaccade inhibition) or toward the chosen target (more consistent
% with microsaccadic suppression of visual detection).
%
% Section 5.3:
% Microsaccade frequency, as a function of time in the trial, conditional
% on color.
%
% Section 6:
% Looking at detection thresholds for a bunch of experiments.  Trying to
% weed out the bad ones (bad behavior days).
%
% Section 7:
% Comparison of psychometric thresholds on trials with an without an eye
% movement.  Pooling data across multiple sessions prior to calculating
% threshold.
%
% Section 8:
% Probability correct as a function of when a microsaccade was made.
% Option to trim saccades based on direction (toward one or the other
% target).
%
% Section 8.1
% Threshold as a function of when a microsaccade was made.
%
% Section 9:
% Probablity correct on trials in which a saccade was made during a
% specific window.  As a function of spatial frequency and color.  Also
% chisquared tests for equality of proportions (% correct as a function of
% color, etc.)
%
% Section 9.1
% Probability correct on trials in which a saccade was made during a
% specific time window, conditional on target choice.
%
% Section 9.2
% Probability correct on trials in which a saccade was made during a
% specific time window, conditional on color and only using sessions in
% which more errors were made on Achromatic trials
%
% Section 10:
% Comparing parameters (amplitude, direction, frequency) from DTspot and
% white noise.
%
% Section 10.1:
% Getting RF eccentricities from white noise experiments.
%
% Section 11:
% Probability correct and target choice as a function of when a microsaccade was made.
% Analyzing upwards and downward saccades separately (since the stimulus
% drifts upward, an upward eye movement would serve to decrease its retinal
% speed).
%
% Section 12:
% Weibull-fitted psychometric functions for individual color/sf conditions
% computed from trials with a fixational saccade and from trials without a
% fixational saccade (during the stimulus period).  Also permutation test
% for significance in difference between thresholds.  Option to eliminate
% trials on a saccade amplitude criterion.
%
% Section 12.1:
% Looking for a relationship between contrast and microsaccade occurrence.
%
% Section 12.2
% Looking at the Weibull-fitted psychmetric thresholds: comparing no
% saccade conditions to saccade conditions but subject to an amplitude
% criterion.  (So trials with saccades that don't meet the criterion don't
% get used in this analysis at all.)  Based on Section 12.
%
% Section 13:
% Looking at eye position traces synced to frames on.  Testing whether the
% delay between fixation achieved and frames on was always at least 300 ms.
%
% Section 14:
% Looking at saccade-triggered PSTHs from cells studied in DTspot recording
% experiments.
%
% Section 14.1:
% Saccade-triggered PSTHs from cells studied in DTspot recording
% experiments.  L-M only. (DTspotL-M.txt)
%
% Section 15:
% Comparing performance (and microsaccadic suppression) across potential
% S-cone isolation color directions (using 2 and 10 deg fundamentals)
% at high and low temporal frequency.
%
% Section 16
% Is there a correlation between the time of the last microsaccade and the
% latency to the operant saccade (the reaction time?)
%
% Section 17
% How does the magnitude of microsaccadic suppression vary with amplitude?
% (During the stimulus plateau only.)
%
% Section 17.1
% How does the magnitude of microsaccadic suppression vary with amplitude?
% As a function of time in the trial.
%
% Section 18
% How do behavioral/microsaccadic measures vary from the earliest data
% collection sessions to the most recent?
%%
% Section 1
% Generates a few summary plots of saccadic eye movements made in the
% fixation window.

filenames = fnamesFromTxt;
amplitudes = [];
directions = [];
peakv = [];
pathlengths = [];
durations = [];
BINWIDTH = .01;
timebins = [-.5:BINWIDTH:.767];
PSTH = zeros(2,length(timebins));  % histogram of saccade initiations, not spikes

totalnumtrials = [0 0];
for a = 1:size(filenames,1)
    stro = nex2stro(findfile(filenames(a,:)));
    %    stro = DTfilterquesttrials(stro,'PaperDefaults');
    if isempty(stro.sum.analog.sigid) || isempty(stro.trial)
        continue;
    end
    ntrials = size(stro.trial,1);
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_on'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_off'));
    frameon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'frame_on'));
    targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_on'));
    correct = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
    sacstats = getSacData(stro);
    close;
    
    for i = 1:ntrials
        st = sacstats.starttimes{i};
        Lsac = (sacstats.starttimes{i} > frameon_t(i)) & (sacstats.endtimes{i} < targon_t(i));
        if any(Lsac)
            if any(sacstats.amplitudes{i}(Lsac) > 2)
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
            PSTH(correct(i)+1,:) = PSTH(correct(i)+1,:) + n;
        end
        totalnumtrials(correct(i)+1) = totalnumtrials(correct(i)+1) + 1;
    end
end

figure;
subplot(3,2,1);
hist(amplitudes,50);
xlabel('amplitude'); ylabel('count');

subplot(3,2,2);
[rho, theta]= hist(directions,20);
polar([theta theta(1)],[rho rho(1)],'k-');

subplot(3,2,3);
[n,x] = hist2([amplitudes,peakv],[50 40]);
imagesc(flipud(n)); colormap(gray); axis xy;
m = (x{1}(end)-x{1}(1))./diff(get(gca,'XLim')+[.5 -.5]);
b = x{1}(1)-m;  % assuming first bin is at '1'
set(gca,'XTick',([0 .25 .5 .75 1]-b)./m,'XTickLabel',[0 .25 .5 .75 1]);
m = (x{2}(end)-x{2}(1))./diff(get(gca,'YLim')+[.5 -.5]);
b = x{2}(1)-m;  % assuming first bin is at '1'
bins = 25:50:150;
set(gca,'YTick',(bins-b)./m,'YTickLabel',bins);
xlabel('amplitude (deg)'); ylabel('peak vel. (deg/sec)');

subplot(3,2,4);
[n,x] = hist2([amplitudes,pathlengths],50);
imagesc(flipud(n)); colormap(gray); axis xy;
m = (x{1}(end)-x{1}(1))./diff(get(gca,'XLim')+[.5 -.5]);
b = x{1}(1)-m;  % assuming first bin is at '1'
set(gca,'XTick',([0 .25 .5 .75 1]-b)./m,'XTickLabel',[0 .25 .5 .75 1]);

m = (x{2}(end)-x{2}(1))./diff(get(gca,'YLim')+[.5 -.5]);
b = x{2}(1)-m;  % assuming first bin is at '1'
set(gca,'YTick',([0:.5:1.5]-b)./m,'YTickLabel',[0:.5:1.5]);
xlabel('amplitude (deg)'); ylabel('traj. length (deg)');


subplot(3,2,5); hold on;
plot(timebins, sum(PSTH)./(BINWIDTH*sum(totalnumtrials)),'k-');
plot([0 0],[0 3],'b:');
set(gca,'Xlim',[min(timebins) max(timebins)]);
xlabel('time wrt stimulus onset (s)'); ylabel('saccades/sec');

subplot(3,2,6); hold on;
plot(timebins, PSTH(1,:)./(BINWIDTH*totalnumtrials(1)),'k-');
plot(timebins, PSTH(2,:)./(BINWIDTH*totalnumtrials(2)),'b-');
legend('inc','cor');
plot([0 0],[0 3],'b:');
set(gca,'Xlim',[min(timebins) max(timebins)]);
xlabel('time wrt stimulus onset (s)'); ylabel('saccades/sec');

figure;
plot(timebins,PSTH'./repmat(totalnumtrials,length(timebins),1));

% What is the duration of a 0.16 deg sacccade?
plot(amplitudes,durations,'.');
lsline
%%
% Section 1.1:
% Microsaccade counts/frequencies during various trial epochs.
% Separated by correct/incorrect

filenames = fnamesFromTxt;
% epoch 1: frames on to stim on
% epoch 2: stim
% epoch 3: delay
% epoch 4: whole trial
data = [];
for a = 1:size(filenames,1)
    stro = nex2stro(findfile(filenames(a,:)));
    stro = DTfilterquesttrials(stro, 'PaperDefaults');
    if isempty(stro.sum.analog.sigid) || isempty(stro.trial)
        continue;
    end
    ntrials = size(stro.trial,1);
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_on'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_off'));
    frameon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'frame_on'));
    targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_on'));
    correct = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
    sacstats = getSacData(stro);
    close;
    
    for i = 1:ntrials
        st = sacstats.starttimes{i};
        sacrate1 = sum(st > frameon_t(i) & st < stimon_t(i))./(stimon_t(i)-frameon_t(i));
        sacrate2 = sum(st > stimon_t(i) & st < stimoff_t(i))./(stimoff_t(i)-stimon_t(i));
        sacrate3 = sum(st > stimoff_t(i) & st < targon_t(i))./(targon_t(i)-stimoff_t(i));
        sacrate4 = sum(st > frameon_t(i) & st < targon_t(i))./(targon_t(i)-frameon_t(i));
        saccount = sum(st > frameon_t(i) & st < targon_t(i));
        data = [data; correct(i) saccount sacrate1 sacrate2 sacrate3 sacrate4];
    end
end

Lcorrect = logical(data(:,1) == 1);
for i = 2:size(data,2)
    [h,p] = ttest2(data(Lcorrect,i), data(~Lcorrect,i))
end
mean(data(Lcorrect,:))-mean(data(~Lcorrect,:))

%%
% Section 1.2:
% Corrects and errors as a function of target choice

filenames = fnamesFromTxt;
data = [];
offsets = [.167 .5];   % relative to stim on, when to look for saccades
for a = 1:size(filenames,1)
    stro = nex2stro(findfile(filenames(a,:)));
    %    stro = DTfilterquesttrials(stro,'PaperDefaults');
    if isempty(stro.sum.analog.sigid) || isempty(stro.trial)
        continue;
    end
    ntrials = size(stro.trial,1);
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_on'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_off'));
    correct = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
    flashside = -1*sign(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_x')));
    flashside = max(flashside, 0); % 1 is to the left
    choice = (correct & flashside) | (~correct & ~flashside);
    sacstats = getSacData(stro);
    close;
    sacrate = zeros(ntrials,1);
    for i = 1:ntrials
        st = sacstats.starttimes{i};
        sacrate(i) = sum(st > stimon_t(i)+offsets(1) & st < stimon_t(i)+offsets(2))./...
            (stimon_t(i)+offsets(2)-stimon_t(i)-offsets(1));
    end
    
    data = [data; sum(correct & choice)./sum(choice) sum(correct & ~choice)./sum(~choice) mean(sacrate(logical(choice))) mean(sacrate(logical(~choice)))];
end
[h,p] = ttest(data(:,1)-data(:,2))
mean(data)
[h,p] = ttest(data(:,3)-data(:,4))

%%
% Section 2
% Reading in a list of DTspot data files and looking at detection
% performance as a function of spatial frequency, color direction, and the
% occurrence of fixational eye movements.
sigmalims = [0 9];  % minimum and maximum sigma (in DVA x10)

filenames = fnamesFromTxt;
data = [];
offset = [.167 .5];
for i = 1:size(filenames,1)
    stro = nex2stro(findfile(filenames(i,:)));
    stro = DTfilterquesttrials(stro,'PaperDefaults');
    
    sigmas = stro.trial(:,strcmp(stro.sum.trialFields(1,:), 'gabor_sigma'));
    if (any(sigmas < min(sigmalims)) | any(sigmas > max(sigmalims)))
        continue;
    end
    if isempty(stro.sum.analog.sigid) || isempty(stro.trial)
        continue;
    end
    ntrials = size(stro.trial,1);
    lms = stro.sum.exptParams.RF_colors;
    lms = reshape(lms,[3,size(lms,1)/3]);
    lms(:,all(lms == 0)) = [];
    whichcol = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'color_dir'));
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_on'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_off'));
    spatialPeriods = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'gabor_lambda'));
    sfs = stro.sum.exptParams.pixperdeg ./ spatialPeriods;
    correct = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
    sacstats = getSacData(stro);
    close;
    % Looping over trials
    saccadeoccurred = nan*ones(ntrials,1);
    for j = 1:ntrials
        st = sacstats.starttimes{j};
        Lsac = (st > stimon_t(j)+offset(1)) & (st < stimon_t(j)+offset(2));
        if (any(Lsac))
            saccadeoccurred(j) = 1;
        else
            saccadeoccurred(j) = 0;
        end
    end
    trialtypes = unique([sfs whichcol],'rows');
    tmpdata = zeros(size(trialtypes,1),1);
    for j = 1:size(trialtypes,1)
        L = sfs == trialtypes(j,1) & whichcol == trialtypes(j,2);
        r = corrcoef(correct(L), saccadeoccurred(L));
        tmpdata(j) = r(1,2);  % Could be Nan if one of the two vectors has var = 0
    end
    data = [data; lms(:,trialtypes(:,2))' trialtypes(:,1) tmpdata];
end

% Preprocessing
cols = data(:,[1:3]);
norms = sqrt(sum(cols'.^2))';
cols = cols./repmat(norms,1,3);
cols(cols(:,1) < 0,:) = -cols(cols(:,1) < 0,:); % convention: L-cone always positive
uniquecols = unique(cols,'rows')
for i = 1:size(uniquecols,1)
    L = softEq(uniquecols(i,:),uniquecols,10,'rows');
    L = ismember(cols,uniquecols(L,:),'rows');
    data(L,[1 2 3]) = repmat(uniquecols(i,:),sum(L),1);
end
uniquecols = unique(data(:,[1 2 3]),'rows')
uniquesfs = unique(data(:,4));

% Plotting
figure; subplot(1,2,1); hold on;
hist(data(:,5),20);
plot(nanmean(data(:,5)),0,'w*','MarkerSize',10);
[h,p] = ttest(data(:,5));
title(['p = ',num2str(p)]);
xlabel('Correlation btn sac. & correct');
ylabel('count');
subplot(1,2,2); hold on;
set(gcf,'InvertHardCopy','off');
set(gca,'Color',[.5 .5 .5]); set(gca,'XScale','log','XLim',[.2 4])
for i = 1:size(uniquecols,1)
    color = uniquecols(i,:);
    Lcolor = all(softEq(cols, repmat(color,size(data,1),1)),2);
    if (color(1) == color(2) && color(2) == color(3))
        % It's achromatic so do nothing
    else
        color(:,[1, 2]) = color(:,[1, 2])*3;
    end
    color = color./(2*max(abs(color)));
    tmp =[];
    for j = 1:length(uniquesfs)
        Lsf = data(:,4) == uniquesfs(j);
        tmp = [tmp; nanmean(data(Lcolor&Lsf,5)) nanstd(data(Lcolor&Lsf,5))...
            sum(Lcolor&Lsf)-sum(isnan(data(Lcolor&Lsf,5)))];
    end
    errorbar(uniquesfs,tmp(:,1),tmp(:,2)./sqrt(tmp(:,3)),'k-','color',.5+color,'linewidth',2);
end


% ANOVA
groups = nan*ones(size(data,1),2);
groups(:,2) = data(:,4);  % sf
for i = 1:size(uniquecols,1)
    color = uniquecols(i,:);
    L = all(softEq(data(:,[1 2 3]), repmat(color,size(data,1),1)),2);
    groups(L,1) = i;
end

[p,t,stats,terms]=anovan(data(:,end),{groups(:,1) data(:,4)},...
    'model',2, 'sstype',2, 'varnames',{'Color', 'SF'});

% Another plot
figure; axes; hold on;
set(gcf,'InvertHardCopy','off');
set(gca,'Color',[.5 .5 .5]); set(gca,'XScale','log','XLim',[.2 4])
for i = 1:size(data,1)
    color = data(i,[1 2 3]);
    if (color(1) == color(2) && color(2) == color(3))
        % It's achromatic so do nothing
    else
        color(:,[1, 2]) = color(:,[1, 2])*3;
    end
    color = color./(2*max(abs(color)));
    plot(data(i,4),data(i,5),'k.','color',.5+color,'MarkerSize',15);
end
plot([min(data(:,4)),max(data(:,4))],[0 0],'k:');
xlabel('SF (cpd)');
ylabel('Corr btn sac & correct');
text(min(data(:,4)), .32, {['Col.: p=',num2str(p(1),2)],...
    ['SF: p=',num2str(p(2),2)],...
    ['ColxSF: p=',num2str(p(3),2)]})

%%
% Section 2.5
% Looking at the mean correlation between saccade occurrence and correct as
% a function of choice.

sigmalims = [0 9];  % minimum and maximum sigma (in DVA x10)
filenames = fnamesFromTxt;
data = [];
offset = [.2 -.2];
for i = 1:size(filenames,1)
    stro = nex2stro(findfile(filenames(i,:)));
    stro = DTfilterquesttrials(stro,0);
    
    sigmas = stro.trial(:,strcmp(stro.sum.trialFields(1,:), 'gabor_sigma'));
    if (any(sigmas < min(sigmalims)) | any(sigmas > max(sigmalims)))
        continue;
    end
    
    if isempty(stro.sum.analog.sigid) || isempty(stro.trial)
        continue;
    end
    ntrials = size(stro.trial,1);
    lms = stro.sum.exptParams.RF_colors;
    lms = reshape(lms,[3,size(lms,1)/3]);
    lms(:,all(lms == 0)) = [];
    whichcol = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'color_dir'));
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_on'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_off'));
    spatialPeriods = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'gabor_lambda'));
    sfs = stro.sum.exptParams.pixperdeg ./ spatialPeriods;
    correct = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
    flashside = -1*sign(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_x')));
    flashside = max(flashside, 0); % 1 is to the left
    choice = (correct & flashside) | (~correct & ~flashside);
    sacstats = getSacData(stro);
    close;
    % Looping over trials
    saccadeoccurred = nan*ones(ntrials,1);
    for j = 1:ntrials
        st = sacstats.starttimes{j};
        Lsac = (st > stimon_t(j)+offset(1)) & (st < stimoff_t(j)+offset(2));
        if (any(Lsac))
            saccadeoccurred(j) = 1;
        else
            saccadeoccurred(j) = 0;
        end
    end
    trialtypes = unique([sfs whichcol],'rows');
    tmpdata = zeros(size(trialtypes,1),2);
    for j = 1:size(trialtypes,1)
        L = sfs == trialtypes(j,1) & whichcol == trialtypes(j,2);
        r1 = corrcoef(correct(L&choice), saccadeoccurred(L&choice));
        r2 = corrcoef(correct(L&~choice), saccadeoccurred(L&~choice));
        tmpdata(j,:) = [r1(1,2), r2(1,2)];
    end
    data = [data; lms(:,trialtypes(:,2))' trialtypes(:,1) tmpdata];
end
figure;
for i = 1:2
    subplot(2,1,i);
    hist(data(:,4+i),30); hold on;
    plot(nanmean(data(:,4+i)),0,'y*','Markersize',10);
    [h,p] = ttest(data(:,4+i));
    title(['p = ',num2str(p)]);
end
xlabel('Correlation btn sac. & correct');
%%
% Section 3
%
% Sliding window for calculating correlation between saccade occurrence and
% correct.  This code can also do the identical analysis on choice
% direction.  (Or correct as a function of choice direction).

whichanalysis = 'both';  % 'choice', 'correct', 'both'
filenames = fnamesFromTxt;
binwidth = .1; %sec
offsets = [-.2:binwidth:.6];
data = [];
for i = 1:size(filenames,1)
    stro = nex2stro(findfile(filenames(i,:)));
    %  stro = DTfilterquesttrials(stro,'PaperDefaults');
    if (isempty(stro.sum.analog.sigid) || isempty(stro.trial))
        continue;
    end
    ntrials = size(stro.trial,1);
    lms = stro.sum.exptParams.RF_colors;
    lms = reshape(lms,[3,size(lms,1)/3]);
    whichcol = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'color_dir'));
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_on'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_off'));
    spatialPeriods = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'gabor_lambda'));
    sfs = stro.sum.exptParams.pixperdeg ./ spatialPeriods;
    correct = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
    flashside = -1*sign(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_x')));
    flashside = max(flashside, 0); % 1 is to the left
    choice = (correct & flashside) | (~correct & ~flashside);
    sacstats = getSacData(stro);
    close;
    % Looping over trials
    saccadeoccurred = nan*ones(ntrials,1);
    tmpdata = [];
    for k = 1:length(offsets)
        for j = 1:ntrials
            st = sacstats.starttimes{j};
            Lsac = (st > stimon_t(j)+offsets(k)-binwidth/2) & (st < stimon_t(j)+offsets(k)+binwidth/2);
            if (any(Lsac))
                saccadeoccurred(j) = 1;
            else
                saccadeoccurred(j) = 0;
            end
        end
        if strcmp(whichanalysis,'choice')
            r = corrcoef(choice, saccadeoccurred);
            tmpdata(k) = r(1,2);
        elseif strcmp(whichanalysis,'correct')
            r = corrcoef(correct, saccadeoccurred);
            tmpdata(k) = r(1,2);
        elseif strcmp(whichanalysis,'both')
            r1 = corrcoef(correct(choice), saccadeoccurred(choice));
            r2 = corrcoef(correct(~choice), saccadeoccurred(~choice));
            tmpdata(k,:) = [r1(1,2), r2(1,2)];
        end
    end
    if (strcmp(whichanalysis,'both'))
        tmpdata = tmpdata(:)';  % concatenating T1 and T2 choices
    end
    data = [data; tmpdata]
end

figure; axes; hold on;
errorbar(offsets,nanmean(data(:,[1:length(offsets)])),...
    1.96*nanstd(data(:,[1:length(offsets)]))./sqrt(sum(~isnan(data(:,[1:length(offsets)])))));
if (size(data,2) > length(offsets))
    errorbar(offsets,nanmean(data(:,[length(offsets)+1:end])),...
        1.96*nanstd(data(:,[length(offsets)+1:end]))./sqrt(sum(~isnan(data(:,[length(offsets)+1:end])))),'k-');
end
plot([offsets(1) offsets(end)],[0 0],'k:');
if strcmp(whichanalysis,'choice')
    ylabel('Correlation saccade and choice');
elseif strcmp(whichanalysis,'correct')
    ylabel('Correlation saccade and correct');
end
xlabel('time');

%%
% Section 4
%
% Distribution of saccade directions at different points during the trial
% as a function of choice direction (shows if microsaccades tend to be
% directed toward the to-be-chosen target) or correct/incorrect (shows if
% microsaccades in a particular direction tend to be associated with target
% visibility/invisibility).

whichanalysis = 'choice';  % 'choice', 'correct', 'choice|correct'
filenames = fnamesFromTxt;
binwidth = .1; %sec
offsets = [-.2:binwidth:.8];
data1 = cell(size(offsets)); % T1 or correct
data2 = cell(size(offsets)); % T2 or incorrect
n1 = 0;
n2 = 0;

for i = 1:size(filenames,1)
    stro = nex2stro(findfile(filenames(i,:)));
    stro = DTfilterquesttrials(stro);
    
    if (isempty(stro.sum.analog.sigid) || isempty(stro.trial))
        continue;
    end
    ntrials = size(stro.trial,1);
    whichcol = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'color_dir'));
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_on'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_off'));
    targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_on'));
    correct = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
    flashside = -1*sign(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_x')));
    flashside = max(flashside, 0); % 1 is to the left
    sacstats = getSacData(stro);
    close;
    
    % Looping over trials
    tmpdata1 = [];  tmpdata2 = [];
    for k = 1:length(offsets)
        for j = 1:ntrials
            dirs = sacstats.directions{j};
            st = sacstats.starttimes{j};
            st(sacstats.endtimes{j} >= targon_t(j)) = [];  % Omitting saccades that terminate after targets appear
            Lsac = (st > stimon_t(j)+offsets(k)-binwidth/2) & (st < stimon_t(j)+offsets(k)+binwidth/2);
            if (any(Lsac))
                if (strcmp(whichanalysis,'choice'))
                    if ((correct(j) && flashside(j)) || (~correct(j) && ~flashside(j)))
                        tmpdata1 = [tmpdata1; dirs(Lsac)];
                    else
                        tmpdata2 = [tmpdata2; dirs(Lsac)];
                    end
                elseif (strcmp(whichanalysis,'correct'))
                    if (correct(j))
                        tmpdata1 = [tmpdata1; dirs(Lsac)];
                    else
                        tmpdata2 = [tmpdata2; dirs(Lsac)];
                    end
                elseif (strcmp(whichanalysis,'choice|correct'))
                    if (correct(j) && flashside(j))
                        tmpdata1 = [tmpdata1; dirs(Lsac)];
                    elseif (correct(j) && ~flashside(j))
                        tmpdata2 = [tmpdata2; dirs(Lsac)];
                    end
                end
            end
        end
        data1{k} = [data1{k}; tmpdata1];
        data2{k} = [data2{k}; tmpdata2];
        if (strcmp(whichanalysis,'choice'))
            nt1ch = sum((correct & flashside) | (~correct & ~flashside));
            n1 = n1 + nt1ch;
            n2 = n2 + ntrials - nt1ch;
        elseif (strcmp(whichanalysis,'correct'))
            n1 = n1 + sum(correct);
            n2 = n2 + sum(~correct);
        elseif (strcmp(whichanalysis,'choice|correct'))
            n1 = n1 + sum(correct & flashside);
            n2 = n2 + sum(correct & ~flashside);
        end
    end
end

% Plotting
maxextent = 0;
bins = linspace(-pi,pi,20);
figure('Position', [-103 745 1284 99]);
for i = 1:length(offsets)
    subplot(1,length(offsets),i); hold on;
    [count1,x] = hist(data1{i},bins);
    [count2,x] = hist(data2{i},bins);
    plot([count1 count1(1)]./n1.*cos([bins bins(1)]), [count1 count1(1)]./n1.*sin([bins bins(1)]),'b-');
    plot([count2 count2(1)]./n2.*cos([bins bins(1)]), [count2 count2(1)]./n2.*sin([bins bins(1)]),'k-');
    maxextent = max([maxextent; get(gca,'XLim')';get(gca,'YLim')']);
    xlabel(num2str(offsets(i)));
    [h,p] = WatsonU2Test(data1{i}, data2{i});
    title(['p = ',num2str(p)]);
end
if (strcmp(whichanalysis,'choice'))
    legend('T1','T2','Location','EastOutside');
elseif (strcmp(whichanalysis,'correct'))
    legend('cor','inc','Location','EastOutside');
end
for i = 1:length(offsets)
    subplot(1,length(offsets),i);
    set(gca,'Xlim',[-maxextent maxextent],'Ylim',[-maxextent maxextent]);
    set(gca,'XTick',[],'YTick',[]);
end

%%
% Section 5
%
% Distribution of saccade amplitudes at different points during the trial
% as a function of choice direction or correct/incorrect (shows if
% microsaccades of particular amplitude tend to be associated with target
% visibility/invisibility).

whichanalysis = 'correct';  % 'choice', 'correct'
filenames = fnamesFromTxt;

binwidth = .1; %sec
offsets = [-.2:binwidth:.8];
data1 = cell(size(offsets)); % T1 or correct
data2 = cell(size(offsets)); % T2 or incorrect
n1 = 0;
n2 = 0;

for i = 1:size(filenames,1)
    stro = nex2stro(findfile(filenames(i,:)));
    %   stro = DTfilterquesttrials(stro);
    
    if (isempty(stro.sum.analog.sigid) || isempty(stro.trial))
        continue;
    end
    ntrials = size(stro.trial,1);
    whichcol = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'color_dir'));
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_on'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_off'));
    targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_on'));
    correct = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
    flashside = -1*sign(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_x')));
    flashside = max(flashside, 0); % 1 is to the left
    sacstats = getSacData(stro);
    close;
    
    % Looping over trials
    tmpdata1 = [];  tmpdata2 = [];
    for k = 1:length(offsets)
        for j = 1:ntrials
            amps = sacstats.amplitudes{j};
            st = sacstats.starttimes{j};
            st(sacstats.endtimes{j} >= targon_t(j)) = [];  % Omitting saccades that terminate after targets appear
            Lsac = (st > stimon_t(j)+offsets(k)-binwidth/2) & (st < stimon_t(j)+offsets(k)+binwidth/2);
            if (any(Lsac))
                if (strcmp(whichanalysis,'choice'))
                    if ((correct(j) && flashside(j)) || (~correct(j) && ~flashside(j)))
                        tmpdata1 = [tmpdata1; amps(Lsac)];
                    else
                        tmpdata2 = [tmpdata2; amps(Lsac)];
                    end
                elseif (strcmp(whichanalysis,'correct'))
                    if (correct(j))
                        tmpdata1 = [tmpdata1; amps(Lsac)];
                    else
                        tmpdata2 = [tmpdata2; amps(Lsac)];
                    end
                end
            end
        end
        data1{k} = [data1{k}; tmpdata1];
        data2{k} = [data2{k}; tmpdata2];
        if (strcmp(whichanalysis,'choice'))
            nt1ch = sum((correct & flashside) | (~correct & ~flashside));
            n1 = n1 + nt1ch;
            n2 = n2 + ntrials - nt1ch;
        elseif (strcmp(whichanalysis,'correct'))
            n1 = n1 + sum(correct);
            n2 = n2 + sum(~correct);
        end
    end
end

% Plotting
bins = linspace(0,2,20);
figure('Position', [-103 745 1284 99]);
for i = 1:length(offsets)
    subplot(1,length(offsets),i); hold on;
    [count1,x] = hist(data1{i},bins);
    [count2,x] = hist(data2{i},bins);
    plot(bins,count1/n1,'b-');
    plot(bins,count2/n2,'k-');
    xlabel(num2str(offsets(i)));
    [h,p] = ttest2(data1{i}, data2{i});
    title(['p = ',num2str(p)]);
    set(gca,'XTick',[],'YTick',[]);
end
if (strcmp(whichanalysis,'choice'))
    legend('T1','T2','Location','EastOutside');
elseif (strcmp(whichanalysis,'correct'))
    legend('cor','inc','Location','EastOutside');
end


%%
% Section 5.1
%
% Distribution of saccade endpoints (amplitudes and direction)
% as a function of time during the trial.

filenames = fnamesFromTxt;

t_binwidth = .1; %sec
s_binwidth = .1; %deg
t_bins = [-.2:t_binwidth:.8];
s_bins = [-.5:s_binwidth:.5];
cordata = zeros(length(s_bins),length(s_bins),length(t_bins));
incdata = zeros(length(s_bins),length(s_bins),length(t_bins));
T1data = zeros(length(s_bins),length(s_bins),length(t_bins));
T2data = zeros(length(s_bins),length(s_bins),length(t_bins));

for i = 1:size(filenames,1)
    stro = nex2stro(findfile(filenames(i,:)));
    %   stro = DTfilterquesttrials(stro);
    
    if (isempty(stro.sum.analog.sigid) || isempty(stro.trial))
        continue;
    end
    ntrials = size(stro.trial,1);
    whichcol = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'color_dir'));
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_on'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_off'));
    targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_on'));
    correct = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
    flashside = -1*sign(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_x')));
    flashside = max(flashside, 0); % 1 is to the left
    choice = (correct & flashside) | (~correct & ~flashside);
    sacstats = getSacData(stro);
    close;
    
    % Looping over trials
    for k = 1:length(t_bins)
        for j = 1:ntrials
            dirs = sacstats.directions{j};
            amps = sacstats.amplitudes{j};
            st = sacstats.starttimes{j};
            L = logical(sacstats.endtimes{j} >= targon_t(j));
            st(L) = [];    % Omitting saccades that terminate after targets appear
            amps(L) = [];  % Omitting saccades that terminate after targets appear
            dirs(L) = [];  % Omitting saccades that terminate after targets appear
            Lsac = (st > stimon_t(j)+t_bins(k)-t_binwidth/2) & (st < stimon_t(j)+t_bins(k)+t_binwidth/2);
            try
                if (any(Lsac))
                    [x,y] = pol2cart(dirs(Lsac),amps(Lsac));
                    xidxs = find(histc(x,s_bins));
                    yidxs = find(fliplr(histc(y,s_bins)));
                    if (correct(j))
                        cordata(yidxs,xidxs,repmat(k,sum(Lsac),1)) =  cordata(yidxs,xidxs,repmat(k,sum(Lsac),1)) + 1;
                    else
                        incdata(yidxs,xidxs,repmat(k,sum(Lsac),1)) =  incdata(yidxs,xidxs,repmat(k,sum(Lsac),1)) + 1;
                    end
                    if (choice(j))
                        T1data(yidxs,xidxs,repmat(k,sum(Lsac),1)) =  T1data(yidxs,xidxs,repmat(k,sum(Lsac),1)) + 1;
                    else
                        T2data(yidxs,xidxs,repmat(k,sum(Lsac),1)) =  T2data(yidxs,xidxs,repmat(k,sum(Lsac),1)) + 1;
                    end
                    
                end
            catch
                keyboard
            end
        end
    end
end

% Differential pics
% T1-T2 (normalized for number of choices)
ntimepoints = length(t_bins);
figure;
data = T1data./sum(T1data(:))-T2data./sum(T2data(:));
lim = max(abs(data(:)));
for i = 1:ntimepoints
    subplot(ceil(sqrt(ntimepoints)),ceil(sqrt(ntimepoints)),i)
    imagesc(data(:,:,i),[-lim lim]);
    axis square;
    colormap(jet);
    set(gca,'XTick',[],'YTick',[]);
    title(num2str(t_bins(i)));
end
set(gcf,'Name','T1-T2');


% correct-incorrect (normalized for number of correct/incorrect)
figure;
data = cordata./sum(cordata(:))-incdata./sum(incdata(:));
lim = max(abs(data(:)));
for i = 1:ntimepoints
    subplot(ceil(sqrt(ntimepoints)),ceil(sqrt(ntimepoints)),i)
    imagesc(data(:,:,i),[-lim lim]);
    axis square;
    colormap(jet);
    set(gca,'XTick',[],'YTick',[]);
    title(num2str(t_bins(i)));
end
set(gcf,'Name','Correct-incorrect');


%%
% Section 5.2
%
% Distribution of saccade endpoints (amplitudes and direction)
% as a function of time during the trial on error trials specifically.

filenames = fnamesFromTxt;

t_binwidth = .1; %sec
s_binwidth = .1; %deg
t_bins = [-.2:t_binwidth:.8];
s_bins = [-.5:s_binwidth:.5];
T1cordata = zeros(length(s_bins),length(s_bins),length(t_bins));
T2cordata = zeros(length(s_bins),length(s_bins),length(t_bins));
T1incdata = zeros(length(s_bins),length(s_bins),length(t_bins));
T2incdata = zeros(length(s_bins),length(s_bins),length(t_bins));

for i = 1:size(filenames,1)
    stro = nex2stro(findfile(filenames(i,:)));
    stro = DTfilterquesttrials(stro,'PaperDefaults');
    
    if (isempty(stro.sum.analog.sigid) || isempty(stro.trial))
        continue;
    end
    ntrials = size(stro.trial,1);
    whichcol = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'color_dir'));
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_on'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_off'));
    targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_on'));
    correct = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
    flashside = -1*sign(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_x')));
    flashside = max(flashside, 0); % 1 is to the left
    choice = (correct & flashside) | (~correct & ~flashside);
    sacstats = getSacData(stro);
    close;
    
    % Looping over trials
    for k = 1:length(t_bins)
        for j = 1:ntrials
            dirs = sacstats.directions{j};
            amps = sacstats.amplitudes{j};
            st = sacstats.starttimes{j};
            L = logical(sacstats.endtimes{j} >= targon_t(j));
            st(L) = [];    % Omitting saccades that terminate after targets appear
            amps(L) = [];  % Omitting saccades that terminate after targets appear
            dirs(L) = [];  % Omitting saccades that terminate after targets appear
            Lsac = (st > stimon_t(j)+t_bins(k)-t_binwidth/2) & (st < stimon_t(j)+t_bins(k)+t_binwidth/2);
            try
                if (any(Lsac))
                    [x,y] = pol2cart(dirs(Lsac),amps(Lsac));
                    xidxs = find(histc(x,s_bins));
                    yidxs = find(fliplr(histc(y,s_bins)));
                    if (correct(j))
                        if (choice(j))
                            T1cordata(yidxs,xidxs,repmat(k,sum(Lsac),1)) =  T1cordata(yidxs,xidxs,repmat(k,sum(Lsac),1)) + 1;
                        else
                            T2cordata(yidxs,xidxs,repmat(k,sum(Lsac),1)) =  T2cordata(yidxs,xidxs,repmat(k,sum(Lsac),1)) + 1;
                        end
                    else
                        if (choice(j))
                            T1incdata(yidxs,xidxs,repmat(k,sum(Lsac),1)) =  T1incdata(yidxs,xidxs,repmat(k,sum(Lsac),1)) + 1;
                        else
                            T2incdata(yidxs,xidxs,repmat(k,sum(Lsac),1)) =  T2incdata(yidxs,xidxs,repmat(k,sum(Lsac),1)) + 1;
                        end
                    end
                end
            catch
                keyboard
            end
        end
    end
end

cordata = log10(T1cordata./T2cordata);
cordata(T1cordata < 2 | T2cordata < 2) = 0;
incdata = log10(T1incdata./T2incdata);
incdata(T1incdata < 2 | T2incdata < 2) = 0;
lim = max(abs([cordata(:); incdata(:)]));

figure
for i = 1:length(s_bins)
    subplot(2,length(s_bins),i);
    imagesc(cordata(:,:,i),[-lim lim]);
    set(gca,'Xtick',[],'Ytick',[]);
    title(num2str(t_bins(i)));
    subplot(2,length(s_bins),i+length(s_bins));
    imagesc(incdata(:,:,i),[-lim lim]);
    set(gca,'Xtick',[],'Ytick',[]);
end
colormap(jet)

% Calculating the difference in endpoint mean (T1 vs T2 incorrects)
for i = 1:length(s_bins)
    pT1 = T1incdata(:,:,i)./sum(sum(T1incdata(:,:,i)));
    pT2 = T2incdata(:,:,i)./sum(sum(T2incdata(:,:,i)));
    mnT1 = [sum([length(s_bins):-1:1]*pT1) sum([1:length(s_bins)]*pT1')];
    mnT2 = [sum([length(s_bins):-1:1]*pT2) sum([1:length(s_bins)]*pT2')];
    differencevector = mnT1-mnT2;
    [t_bins(i) differencevector]
end

% Resampling test on errors (assuming an empirical multinomial distribution for endpoints).
niter = 2000;
tmp = zeros(niter,1);
p = []; differencevector = [];
for frame = 1:length(s_bins)
    tmp = zeros(niter,1);
    nT1 = sum(sum(T1incdata(:,:,frame)));
    nT2 = sum(sum(T2incdata(:,:,frame)));
    empdistn = (T1incdata(:,:,frame)+T2incdata(:,:,frame))./(nT1+nT2);
    
    % Bug workaround from Mathworks site
    % (http://www.mathworks.com/support/bugreports/644205)
    e = min([0 cumsum(empdistn(:))'],1);
    e(end) = 1;
    empdistn = reshape(diff(e),size(empdistn));
    
    for iter = 1:niter
        if (iter == 1)
            resampT1 = T1incdata(:,:,frame);
            resampT2 = T2incdata(:,:,frame);
        else
            resampT1 = reshape(mnrnd(nT1,empdistn(:)),size(empdistn));
            resampT2 = reshape(mnrnd(nT2,empdistn(:)),size(empdistn));
        end
        resampT1 = resampT1./sum(sum(resampT1));
        resampT2 = resampT2./sum(sum(resampT2));
        mnT1 = [sum([length(s_bins):-1:1]*resampT1) sum([1:length(s_bins)]*resampT1')];
        mnT2 = [sum([length(s_bins):-1:1]*resampT2) sum([1:length(s_bins)]*resampT2')];
        tmp(iter) = norm(mnT1-mnT2);
    end
    differencevector(frame,:) = mnT1-mnT2;
    p(frame) = sum(tmp > tmp(1))/niter;
end
[t_bins; p]


% Calculating the difference in endpoint mean (corrects vs incorrect -
% optional: conditional on target)
cordata = T1cordata+T2cordata;
incdata = T1incdata+T2incdata;

% Resampling test on errors (assuming an empirical multinomial distribution for endpoints).
niter = 2000;
tmp = zeros(niter,1);
p = []; differencevector = [];
for frame = 1:length(s_bins)
    tmp = zeros(niter,1);
    ncor = sum(sum(cordata(:,:,frame)));
    ninc = sum(sum(incdata(:,:,frame)));
    empdistn = (cordata(:,:,frame)+incdata(:,:,frame))./(ninc+ncor);
    
    % Bug workaround from Mathworks site
    % (http://www.mathworks.com/support/bugreports/644205)
    e = min([0 cumsum(empdistn(:))'],1);
    e(end) = 1;
    empdistn = reshape(diff(e),size(empdistn));
    
    for iter = 1:niter
        if (iter == 1)
            resampcor = cordata(:,:,frame);
            resampinc = incdata(:,:,frame);
        else
            resampcor = reshape(mnrnd(ncor,empdistn(:)),size(empdistn));
            resampinc = reshape(mnrnd(ninc,empdistn(:)),size(empdistn));
        end
        resampcor = resampcor./sum(sum(resampcor));
        resampinc = resampinc./sum(sum(resampinc));
        mncor = [sum([length(s_bins):-1:1]*resampcor) sum([1:length(s_bins)]*resampcor')];
        mninc = [sum([length(s_bins):-1:1]*resampinc) sum([1:length(s_bins)]*resampinc')];
        tmp(iter) = norm(mncor-mninc);
    end
    differencevector(frame,:) = mncor-mninc;
    p(frame) = sum(tmp > tmp(1))/niter;
end
[t_bins; p]

%%
% Section 5.3
% Microsaccade frequency as a function of color

filenames = fnamesFromTxt;
colordirs = mkbasis([1 1 1; 1 -1 0]');
BINWIDTH = 0.05;
timebins = [-.5:BINWIDTH:.767];
PSTH = zeros(2,length(timebins));  % histogram of saccade initiations, conditional on color

totalnumtrials = [0 0];
for a = 1:size(filenames,1)
    stro = nex2stro(findfile(filenames(a,:)));
    stro = DTfilterquesttrials(stro,'PaperDefaults');
    if isempty(stro.sum.analog.sigid) || isempty(stro.trial)
        continue;
    end
    ntrials = size(stro.trial,1);
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_on'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_off'));
    frameon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'frame_on'));
    targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_on'));
    correct = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
    lms = stro.sum.exptParams.RF_colors;
    lms = mkbasis(reshape(lms,[3,size(lms,1)/3]));
    lms(:,all(lms == 0)) = [];
    whichcol = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'color_dir'));
    sacstats = getSacData(stro);
    close;
    
    for i = 1:ntrials
        st = sacstats.starttimes{i};
        colidx = find(lms(:,whichcol(i))'*colordirs > .99);
        sactimes = [];
        tmp = st-stimon_t(i);
        tmp((tmp < timebins(1)-BINWIDTH/2) | (tmp > timebins(end)+BINWIDTH/2)) = [];
        [n,x] = hist(tmp, timebins);
        PSTH(colidx,:) = PSTH(colidx,:) + n;
        totalnumtrials(colidx) = totalnumtrials(colidx) + 1;
    end
end

figure; axes; hold on;
plot(timebins, PSTH(1,:)./(BINWIDTH*totalnumtrials(1)),'k-');
plot(timebins, PSTH(2,:)./(BINWIDTH*totalnumtrials(2)),'r-');

L = timebins > .167 & timebins < .5
PSTH*L'  %  microsaccade counts during stimulus plateau

%%
% Section 6
% Looking at thresholds and percent correct across a bunch of experiments.  Are there outliers?
colordirs = mkbasis([1 1 1; 0 0 1; 1 -1 0]');
periods = [206 58 16];
filenames = fnamesFromTxt;
thresholds = nan*ones(length(colordirs), length(periods), size(filenames,1));
ncor = nan*ones(length(colordirs), length(periods), size(filenames,1));
ntot = nan*ones(length(colordirs), length(periods), size(filenames,1));

for a = 1:size(filenames,1)
    stro = nex2stro(findfile(filenames(a,:)));
    if isempty(stro.sum.analog.sigid) || isempty(stro.trial)
        continue;
    end
    
    [tmpthresholds, tmpcolorDirs, tmpsfs] = DTquestUnpackGH(stro, 'mode');
    %Color is on the rows, sf is on the columns of tmpthresholds
    stro = DTfilterquesttrials(stro,'PaperDefaults');
    correct = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
    spatialPeriods = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'gabor_lambda'));
    whichcol = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'color_dir'));
    
    tmpcolorDirs = mkbasis(tmpcolorDirs')';
    tmpperiods = stro.sum.exptParams.pixperdeg./tmpsfs;
    for i = 1:size(colordirs,1)
        for j = 1:length(periods)
            coloridx = find(abs(colordirs(:,i)'*mkbasis(tmpcolorDirs')) > .99);
            periodidx = find(periods(j) == round(tmpperiods));
            if (~isempty(coloridx) & ~isempty(periodidx))
                thresholds(i,j,a) = tmpthresholds(coloridx,periodidx);
                L = spatialPeriods == periods(j) & whichcol == coloridx;
                ntot(i,j,a) = sum(L);
                ncor(i,j,a) = sum(L & logical(correct));
            end
        end
    end
end

figure;
stdevs = nanstd(log10(thresholds),[],3);
for i = 1:size(colordirs,1)
    for j = 1:length(periods)
        whichaxes = sub2ind([size(colordirs,1) length(periods)],i,j);
        subplot(size(colordirs,1),length(periods),whichaxes);
        hist(log10(squeeze(thresholds(i,j,:))))
        set(gca,'Xlim',[-2 2]);
        title(num2str(stdevs(i,j)));
    end
end

% For Sedna
figure; axes; hold on;
plot(squeeze(thresholds(3,3,:)))
whichfiles = squeeze(thresholds(3,3,:)) > 8;
plot(find(whichfiles),squeeze(thresholds(1,1,whichfiles)),'m*')
filenames(whichfiles,:)  % Omit these

% Comparing percent correct for L-M and achromatic
pcor = ncor./ntot;
figure;
for i = 1:3
    subplot(3,1,i); hold on;
    plot(squeeze(pcor(1,i,:)),'k-');
    plot(squeeze(pcor(3,i,:)),'r-');
    ach = squeeze(pcor(1,i,:));
    LvM = squeeze(pcor(3,i,:));
    [h,p] = ttest(ach-LvM);
    [h,p] = equalproptest(sum(ncor([1 3],i,:),3),sum(ntot([1 3],i,:),3),.05)
    title(num2str(p));
end

% Overall proportion correct within sf
nansum(ncor,3)./nansum(ntot,3)

% Proportion correct integrating across sfs
nansum(nansum(ncor,3),2)./nansum(nansum(ntot,3),2)


%%
% Section 7
% Calculating thresholds with and without fixational saccades.  Pooling all
% the data prior to estimating thresholds (for each color/spatial frequency
% condition).

filenames = fnamesFromTxt;
data = {};
offset = [.2 .4];
for i = 1:size(filenames,1)
    stro = nex2stro(findfile(filenames(i,:)));
    %  stro = DTfilterquesttrials(stro,[],nan);
    
    if isempty(stro.sum.analog.sigid) || isempty(stro.trial)
        continue;
    end
    DTindicies;
    ntrials = size(stro.trial,1);
    lms = stro.sum.exptParams.RF_colors;
    lms = reshape(lms,[3,size(lms,1)/3]);
    lms(:,all(lms == 0)) = [];
    whichcol = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'color_dir'));
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_on'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_off'));
    spatialPeriods = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'gabor_lambda'));
    sfs = stro.sum.exptParams.pixperdeg ./ spatialPeriods;
    uniquesfs = unique(sfs);
    correct = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
    exptColors = reshape(stro.sum.exptParams.RF_colors, 3, 3)';
    exptColors = exptColors.*repmat(sign(sum(sign(exptColors+eps),2)),1,size(exptColors,2)); %standardizing color directions
    noColor = sum(abs(exptColors), 2) == 0;
    exptColors(noColor,:) = [];
    norms = sqrt(sum(exptColors.^2, 2));
    exptColors = exptColors./repmat(norms, 1, 3);
    bkgndrgb = [stro.sum.exptParams.bkgnd_r, stro.sum.exptParams.bkgnd_g, stro.sum.exptParams.bkgnd_b];
    M = reshape(stro.sum.exptParams.m_mtx, 3, 3);
    bkgndlms = M * bkgndrgb';
    x = 0:255; %the normal range of the gamma look up table
    xx = linspace(0, 255, 2^16); %the desired quantization of the gammaTable
    g1 = reshape(stro.sum.exptParams.gamma_table, 256, 3);
    gammaTable = [spline(x, g1(:,1), xx)', spline(x, g1(:,2), xx)', spline(x, g1(:,3), xx)'];
    
    RGB = stro.trial(:,flashRInd:flashBInd)+1;
    rgb = [gammaTable(RGB(:,1), 1), gammaTable(RGB(:,2), 2), gammaTable(RGB(:,3), 3)];
    LMS = (M*rgb')-repmat(bkgndlms,1,size(rgb,1));
    ccs = LMS' ./ repmat(bkgndlms',size(rgb,1),1);
    
    
    sacstats = getSacData(stro);
    close;
    % Looping over trials
    saccadeoccurred = nan*ones(ntrials,1);
    for j = 1:ntrials
        st = sacstats.starttimes{j};
        Lsac = (st > stimon_t(j)+offset(1)) & (st < stimon_t(j)+offset(2));
        if (any(Lsac))
            saccadeoccurred(j) = 1;
        else
            saccadeoccurred(j) = 0;
        end
    end
    for j = 1:size(exptColors,1)
        for k = 1:length(uniquesfs)
            L = whichcol == j & sfs == uniquesfs(k);
            if (sum(L) == 0)
                continue;
            end
            contrast = sqrt(sum(ccs(L,:).^2,2));
            tmpdata = [contrast, correct(L), saccadeoccurred(L)];
            foundit = 0;
            for m = 1:length(data)
                if (all(data{m}.color == exptColors(j,:)) & data{m}.sf == uniquesfs(k))
                    data{m}.perf = [data{m}.perf; tmpdata];
                    foundit = 1;
                end
            end
            if (~foundit)
                m = length(data)+1;
                data{m}.color = exptColors(j,:);
                data{m}.sf = uniquesfs(k);
                data{m}.perf = tmpdata;
            end
        end
    end
end
% Now, computing psychometric thresholds
for i = 1:length(data)
    tmp = data{i}.perf;
    contrast = tmp(:,1)*100;
    correct = tmp(:,2);
    saccadeoccurred = logical(tmp(:,3));
    thresholds = [];
    for j = 1:2
        if (j == 1);
            L = ~saccadeoccurred;
        else
            L = saccadeoccurred;
        end
        if (sum(correct(L)) < 8 | sum(~correct(L)) < 8)
            aML = nan;
        else
            [aSSE, bSSE, gSSE] = weibullFit(contrast(L),correct(L), 'sse', [mean(contrast),5]);
            [aML, bML, gML, success] = weibullFit(contrast(L), [correct(L), ~correct(L)], 'mle', [aSSE, bSSE]);
            if (~success | abs(aML - aSSE)>10)
                aML = nan;
            end
        end
        thresholds = [thresholds; aML];
        data{i}.thresholds = thresholds;
    end
end
% Only looking at the "standard" sfs
standardperiods = [206, 58, 16];
omitlist = [];
for i = 1:length(data)
    if (~any(abs(stro.sum.exptParams.pixperdeg./data{i}.sf-standardperiods) < .0001))
        data{i}.sf
        omitlist = [omitlist; i];
    end
end
data(omitlist) = [];
% data{i}.perf = contrast, correct, saccade occurred
% Looking for an effect of "saccade".
% This GLM fit is flawed!  Psychometric function goes from 0 to 1 instead
% of 0.5 to 1!
% link = @(x)(-log(1./(2.*x-1)-1));
% invlink = @(x)(.5*(1./(1+exp(-x))+1));
% derivlink = @(x)(-1./(2.*x.^2-3.*x+1));
%
% % From stattestlink.m below
% link = @(mu) log(mu ./ (1-mu));
% derivlink = @(mu) 1 ./ (mu .* (1-mu));
% invlink = @(eta) 1 ./ (1 + exp(-eta));
%
% statsmat = [];
% for i = 1:length(data)
%     [b,dev,stats] = glmfit(data{i}.perf(:,[1 3]),data{i}.perf(:,2),'binomial','link',{link, derivlink, invlink});
%
%    % b(2) is threshold w/o saccade, b(2)+b(3) is threshold w/ saccade
%     statsmat(i,1) = b(3)./(b(2)+b(3));
%     statsmat(i,2) = stats.p(3);
% end
% statsmat = reshape(statsmat,[3,3,2]);
% % Negative coefficients mean that a saccade (saccade occurred = 1) tends to contribute to
% % incorrect responses (correct = 0).

% Now extracting the threshold ratios
tmpdata = [];
for i = 1:length(data)
    if all(~isnan(data{i}.thresholds))
        tmpdata = [tmpdata; data{i}.color, data{i}.sf, data{i}.thresholds(1), data{i}.thresholds(2)]
    end
end
figure;
hist(log(tmpdata(:,5)./tmpdata(:,6)),20);  % without a saccade/with a saccade
geomean(tmpdata(:,5)./tmpdata(:,6))  % a 10% loss of sensitivity?
exp(median(log((tmpdata(:,5)./tmpdata(:,6)))))
[p, h] = signrank(tmpdata(:,5)-tmpdata(:,6));
xlabel('Log threshold ratio');

%% Section 8:
% Probability correct as a function of time between saccade time and go
% signal.  Option to omit saccades based on direction (towards a target)

filenames = fnamesFromTxt;
NOSACTOWARDTARG = 0;
WEDGEWIDTH = pi/2;
% CONDITIONEDONSAC = 1;  % Thsis should always be '1'

binwidth = .05; %sec
bins = [-.2:binwidth:.8];
x = zeros(length(bins),3);
n = zeros(length(bins),3);
colordirs = mkbasis([1 1 1; 1 -1 0; 0 0 1]');

allsacdirections = [];
for i = 1:size(filenames,1)
    stro = nex2stro(findfile(filenames(i,:)));
    stro = DTfilterquesttrials(stro,10,nan);
    
    if (isempty(stro.sum.analog.sigid) || isempty(stro.trial))
        continue;
    end
    ntrials = size(stro.trial,1);
    T1angle = atan2(stro.sum.exptParams.rf_x, stro.sum.exptParams.rf_y);
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_on'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_off'));
    targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_on'));
    correct = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
    whichcol = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'color_dir'));
    lms = stro.sum.exptParams.RF_colors;
    lms = reshape(lms,[3,size(lms,1)/3]);
    lms(:,all(lms == 0)) = [];
    lms = mkbasis(lms);
    sacstats = getSacData(stro);
    close;
    
    
    % Looping over trials
    for k = 1:length(bins)
        for j = 1:ntrials
            coloridx = find(abs(lms(:,whichcol(j))'*colordirs) > .99);  % Avoiding roundoff errors
            if (isempty(coloridx))
                continue;
            end
            dirs = sacstats.directions{j};
            amps = sacstats.amplitudes{j};
            st = sacstats.starttimes{j};
            L = logical(sacstats.endtimes{j} >= targon_t(j));   % Omitting saccades that terminate after targets appear
            if (NOSACTOWARDTARG) % Omitting saccades towards a target
                rotsacdir = mod(sacstats.directions{j}-T1angle, 2*pi);
                Lt1 = logical(rotsacdir < WEDGEWIDTH/2 | rotsacdir >  mod(-WEDGEWIDTH/2,2*pi));
                Lt2 = logical(rotsacdir > pi-WEDGEWIDTH/2 & rotsacdir <  pi+WEDGEWIDTH/2);
                %  disp(['Eliminating ',num2str(sum(Lt1|Lt2)),' saccades'])
                L = L | Lt1 | Lt2;
            end
            st(L) = [];
            amps(L) = [];
            dirs(L) = [];
            
            Lsac = (st > stimon_t(j)+bins(k)-binwidth/2) & (st < stimon_t(j)+bins(k)+binwidth/2);
            %             if (~CONDITIONEDONSAC)
            %                 Lsac = ~Lsac;   % condition on "no saccade occurred"
            %             end
            if (any(Lsac))
                if (correct(j))
                    x(k,coloridx) = x(k,coloridx)+1;
                end
                n(k,coloridx) = n(k,coloridx)+1;
            end
        end
    end
end

% All together
colors = {'k-','r-','b-'};
figure;
subplot(2,1,1); hold on;
p = sum(x,2)./sum(n,2);
se = sqrt((p.*(1-p))./sum(n,2))
h = patch([bins fliplr(bins)],[p+se; flipud(p-se)]','b')
set(h,'FaceColor',[.5 .5 .5])
plot(bins, p);

plot([0 .167 .5 .666],[0.6 .75 .75 0.6],'k:');
ylabel('Probability correct');
set(gca,'Xlim',[bins(1) bins(end)]);
subplot(2,1,2); hold on;
for i = 1:3
    plot(bins, sum(n,2),char(colors(i)))
end
ylabel('# saccades');
xlabel('Time wrt stim on');
set(gca,'Xlim',[bins(1) bins(end)]);


% Broken down by color direction
figure;
subplot(2,1,1); hold on;
for i = 1:3
    plot(bins, x(:,i)./n(:,i),char(colors(i)));
end

ylabel('Probability correct');
set(gca,'Xlim',[bins(1) bins(end)]);
subplot(2,1,2); hold on;
for i = 1:3
    plot(bins, n(:,i),char(colors(i)))
end
ylabel('# saccades');
xlabel('Time wrt stim on');
set(gca,'Xlim',[bins(1) bins(end)]);

%%
% Section 8.1
% Threshold as a function of when a microsaccade was made.
% Based largely on section 11

filenames = fnamesFromTxt;
binwidth = .05; %sec
bins = [-.2:binwidth:.8];
data = [];
for i = 1:size(filenames,1)
    stro = nex2stro(findfile(filenames(i,:)));
    stro = DTfilterquesttrials(stro,'PaperDefaults');
    
    if isempty(stro.sum.analog.sigid) || isempty(stro.trial)
        continue;
    end
    ntrials = size(stro.trial,1);
    lms = stro.sum.exptParams.RF_colors;
    lms = reshape(lms,[3,size(lms,1)/3]);
    lms(:,all(lms == 0)) = [];
    lms = mkbasis(lms);
    whichcol = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'color_dir'));
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_on'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_off'));
    targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_on'));
    spatialPeriods = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'gabor_lambda'));
    correct = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
    flashR = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_R'));
    flashG = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_G'));
    flashB = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_B'));
    bkgndrgb = [stro.sum.exptParams.bkgnd_r, stro.sum.exptParams.bkgnd_g, stro.sum.exptParams.bkgnd_b];
    M = reshape(stro.sum.exptParams.m_mtx, 3, 3);
    bkgndlms = M * bkgndrgb';
    x = 0:255; %the normal range of the gamma look up table
    xx = linspace(0, 255, 2^16); %the desired quantization of the gammaTable
    g1 = reshape(stro.sum.exptParams.gamma_table, 256, 3);
    gammaTable = [spline(x, g1(:,1), xx)', spline(x, g1(:,2), xx)', spline(x, g1(:,3), xx)'];
    RGB = [flashR, flashG, flashB];
    rgb = [gammaTable(RGB(:,1)+1, 1), gammaTable(RGB(:,2)+1, 2), gammaTable(RGB(:,3)+1, 3)];
    cc = ((M*rgb')'-repmat(bkgndlms',size(RGB,1),1)) ./ repmat(bkgndlms',size(RGB,1),1);
    ccnorms = sqrt(sum(cc.^2,2));
    sacstats = getSacData(stro);
    close;
    
    for j = 1:ntrials
        tmp = zeros(1,length(bins));
        for k = 1:length(bins)
            if (stro.trial(j,strcmp(stro.sum.trialFields(1,:),'numframes'))  == 0)
                continue;  % Why are there zero frame trials?
            end
            st = sacstats.starttimes{j};
            tmp(k) = any((st > stimon_t(j)+bins(k)-binwidth/2) & (st < stimon_t(j)+bins(k)+binwidth/2));
        end
        data = [data; tmp correct(j) ccnorms(j)];
    end
end

cd 'C:\Matlab Code\Analysis\Sandbox\Greg\DTemstuff'
thresholds = [];
for k = 1:length(bins)
    Lsac = logical(data(:,k));
    [fittedparams, success(1)] = weibullFitGH(data(:,end), data(:,end-1), 'sse', [mean(data(:,end)), 1]);
    [fittedparams_full, success(2)] = weibullFitGH([data(:,end), Lsac], data(:,end-1), 'mle', [fittedparams(1),fittedparams(2)]);    thresholds(k) = log10((fittedparams_full(1)+fittedparams_full(4))/fittedparams_full(1))
end
plot(thresholds,'m-')
%% Section 9
% Probablity correct on trials in which a saccade was made during a
% specific time window.  As a function of spatial frequency and color.

filenames = fnamesFromTxt;
offsets = [.167 .5];   % relative to stim on, when to look for saccades
colordirs = mkbasis([1 1 1; 1 -1 0]');
%periods = [16 23 33 48 58 69 100 143 206];
periods = [16 58 206];
NOSACTOWARDTARG = 0;
WEDGEWIDTH = pi/2;
OMITCURVEDSACS = 0;

x_sac = zeros(size(colordirs,2), length(periods));
n_sac = zeros(size(colordirs,2), length(periods));
x_nosac = zeros(size(colordirs,2), length(periods));
n_nosac = zeros(size(colordirs,2), length(periods));
curvedomitcount = [0 0];
for i = 1:size(filenames,1)
    stro = nex2stro(findfile(filenames(i,:)));
    stro = DTfilterquesttrials(stro,'PaperDefaults');
    
    if isempty(stro.sum.analog.sigid) || isempty(stro.trial)
        continue;
    end
    ntrials = size(stro.trial,1);
    T1angle = atan2(stro.sum.exptParams.rf_x, stro.sum.exptParams.rf_y);
    lms = stro.sum.exptParams.RF_colors;
    lms = reshape(lms,[3,size(lms,1)/3]);
    lms(:,all(lms == 0)) = [];
    lms = mkbasis(lms);
    whichcol = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'color_dir'));
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_on'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_off'));
    targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_on'));
    spatialPeriods = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'gabor_lambda'));
    correct = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
    sacstats = getSacData(stro);
    close;
    
    % Looping over trials
    for j = 1:ntrials
        coloridx = find(abs(lms(:,whichcol(j))'*colordirs) > .99);  % Avoiding roundoff errors
        if (isempty(coloridx))
            disp('Got here 1');
            continue;
        end
        periodidx = find(spatialPeriods(j) == periods);
        if (isempty(periodidx))
            disp('Got here 2');
            %  keyboard
            continue;
        end
        st = sacstats.starttimes{j};
        L = logical(sacstats.endtimes{j} >= targon_t(j));
        if (NOSACTOWARDTARG) % Omitting saccades towards a target
            rotsacdir = mod(sacstats.directions{j}-T1angle, 2*pi);
            Lt1 = logical(rotsacdir < WEDGEWIDTH/2 | rotsacdir >  mod(-WEDGEWIDTH/2,2*pi));
            Lt2 = logical(rotsacdir > pi-WEDGEWIDTH/2 & rotsacdir <  pi+WEDGEWIDTH/2);
            %  disp(['Eliminating ',num2str(sum(Lt1|Lt2)),' saccades'])
            L = L | Lt1 | Lt2;
        end
        % From DTFixSac, Section 1: amplitude ranges from [0 0.8367], peakv ranges from [20 125]
        % equation of line: peakv = (125-20)/0.8376 * amp + 20
        if (OMITCURVEDSACS)
            a = sacstats.amplitudes{j};
            v = sacstats.peakv{j};
            v_hat = (125-20)/.8367*a+20;
            Lcurve = v > v_hat;
            %             if (any(Lcurve))
            %                 disp('got a curved saccade');
            %             end
            n_curved = sum(Lcurve & ~L);
            n_potentialcurve = sum(ones(size(Lcurve)) & ~L);
            curvedomitcount = curvedomitcount + [n_curved n_potentialcurve]; % curved and not already chucked out
            L = L | Lcurve;
        end
        
        st(L) = [];    % Omitting saccades that terminate after targets appear (+ other criteria)
        Lsac = (st > stimon_t(j)+offsets(1)) & (st < stimon_t(j)+offsets(2));
        if (any(Lsac))
            if (correct(j))
                x_sac(coloridx, periodidx) = x_sac(coloridx,periodidx)+1;
            end
            n_sac(coloridx,periodidx) = n_sac(coloridx,periodidx)+1;
        else
            if (correct(j))
                x_nosac(coloridx, periodidx) = x_nosac(coloridx,periodidx)+1;
            end
            n_nosac(coloridx,periodidx) = n_nosac(coloridx,periodidx)+1;
        end
    end
end

% Getting rid of non-existent conditions
sfs = stro.sum.exptParams.pixperdeg ./ periods;
L = sum(n_sac) == 0;
% Only looking at the "standard" periods
L = L | ~ismember(periods,[16 58 206]);
x_sac(:,L) = []; n_sac(:,L) = []; x_nosac(:,L) = []; n_nosac(:,L) = [];
sfs(L) = [];

% x_sac and n_sac: Color on rows, SF on columns
% For color dummy variable Achr = [0 0], L-M = [0 1], S = [1 0];
% design matrix:
%   0) (Implicit intercept)
%   1) 1 = sac, 0 = no sac
%   2) 1 = ~L-M, 0 = ~L-M
%   3) 1 = S, 0 = ~S
%   4) SF

sfmat = repmat(sfs,size(colordirs,2),1);
designmat = [repmat([0 0; 1 0; 0 1],length(sfs),1) sfmat(:)];
designmat = [[ones(size(designmat,1),1); zeros(size(designmat,1),1)], [designmat; designmat]];
interactions = [designmat(:,1).*designmat(:,2) designmat(:,1).*designmat(:,3) designmat(:,1).*designmat(:,4)];
data = [[x_sac(:),n_sac(:)];[x_nosac(:),n_nosac(:)]];
[b,dev,stats] = glmfit([designmat interactions],data,'binomial');

parameters = {'Intercept','Saccade','L-M dummy','S dummy','SF','Saccade*L-M','Saccade*S','Saccade*SF'};
for i = 1:length(parameters)
    str = sprintf('%s: \t%f\t p=%f ',char(parameters{i}),stats.beta(i), stats.p(i));
    if (stats.p(i) < 0.05)
        str = [str,'*'];
    end
    if (stats.p(i) < 0.01)
        str = [str,'*'];
    end
    disp(str);
end

% Looking at the effect of color across SFs
figure;
subplot(1,3,1); hold on;
p1 = sum(x_sac,2)./sum(n_sac,2);
p2 = sum(x_nosac,2)./sum(n_nosac,2);
p =  (sum(x_sac,2)+sum(x_nosac,2))./(sum(x_sac,2)+sum(x_nosac,2));
se = sqrt(p1.*(1-p1)./sum(n_sac,2) + p2.*(1-p2)./sum(n_nosac,2));
bar((p2-p1)./p2,'black');
for i = 1:3
    plot([i i],[-se(i) se(i)]+(p2(i)-p1(i))./p2(i),'k-','Linewidth',2);
end
set(gca,'XTickLabel',{'Ach','L-M','S'});

% Joint data
subplot(1,3,2)
%surfl((x_nosac./n_nosac)-(x_sac./n_sac));
bar3((x_nosac./n_nosac)-(x_sac./n_sac));

set(gca,'XTick',[1 2 3],'XTickLabel',num2str(sfs',2));
set(gca,'YTick',[1 2 3],'YTickLabel',{'Ach','L-M','S'});
colormap(gray);

% Looking at the effect of color across SFs
subplot(1,3,3); hold on;
p1 = sum(x_sac,1)./sum(n_sac,1);
p2 = sum(x_nosac,1)./sum(n_nosac,1);
se = sqrt(p1.*(1-p1)./sum(n_sac,1) + p2.*(1-p2)./sum(n_nosac,1));
bar((p2-p1)./p2,'black');
for i = 1:3
    plot([i i],[-se(i) se(i)]+(p2(i)-p1(i))./p2(i),'k-','Linewidth',2);
end
set(gca,'XTick',[1 2 3],'XTickLabel',num2str(sfs',2));

disp(['Proportion no sac trials ',num2str(1-sum(n_sac(:))/(sum(n_nosac(:))+sum(n_sac(:))))]);
% Does performance vary as a function of color?
[h,p] = equalproptest(sum(x_sac,2)+sum(x_nosac,2),sum(n_sac,2)+sum(n_nosac,2),0.05);
disp(['% correct for different colors: ',num2str([(sum(x_sac,2)+sum(x_nosac,2))./(sum(n_sac,2)+sum(n_nosac,2))]')])
disp(['% correct is the same for different colors: ',num2str(p)]);
% Does saccade probability vary as a function of color?
[h,p] = equalproptest(sum(n_sac,2),sum(n_sac,2)+sum(n_nosac,2),0.05);
disp(['probability of saccade is the same for different colors: ',num2str(p)]);
% Does performance differ between sac and nosac (conditional on color)?
for i = 1:size(n_sac,1)
    [h,p(i)] = equalproptest([sum(x_sac(i,:)) sum(x_nosac(i,:))],[sum(n_sac(i,:)) sum(n_nosac(i,:))],0.05);
    [sum(x_sac(i,:)) sum(x_nosac(i,:))]./[sum(n_sac(i,:)) sum(n_nosac(i,:))]
end
disp(['performance is the same for sac and nosac (conditional on color): ',num2str(p)]);
disp(['performance on saccade trials: ',num2str(sum(x_sac(:))./sum(n_sac(:)))]);
disp(['performance on no saccade trials: ',num2str(sum(x_nosac(:))./sum(n_nosac(:)))]);
[h,p] = equalproptest([sum(x_sac(:)) sum(x_nosac(:))],[sum(n_sac(:)) sum(n_nosac(:))],.05);
disp(['Comparision between these two proportions: ',num2str(p)]);

% Does performance vary as a function of spatial frequency?
[h,p] = equalproptest(sum(x_sac,1)+sum(x_nosac,1),sum(n_sac,1)+sum(n_nosac,1),0.05);
disp(['% correct is the same for different spatial frequencies: ',num2str(p)]);
% Does saccade probability vary as a function of spatial frequency?
[h,p] = equalproptest(sum(n_sac,1),sum(n_sac,1)+sum(n_nosac,1),0.05);
disp(['saccade probability is the same for different spatial frequencies: ',num2str(p)]);
% Does performance differ between sac and nosac (conditional on spatial
% frequency)?  (Two sample test of proportions)
for i = 1:size(n_sac,1)
    [h,p(i)] = equalproptest([sum(x_sac(:,i)) sum(x_nosac(:,i))],[sum(n_sac(:,i)) sum(n_nosac(:,i))],0.05);
end
disp(['performance is the same for sac and nosac (conditional on spatial freq): ',num2str(p)]);

% Looking at the effect of color for individual SFs
figure;
tmp = ((x_nosac./n_nosac)-(x_sac./n_sac));
lims = [min(tmp(:))*.9 max(tmp(:))*1.1];
for i = 1:length(sfs)
    subplot(length(sfs),1,i); % Looking at the effect of color at various SFs
    p1 = x_sac(:,i)./n_sac(:,i);
    p2 = x_nosac(:,i)./n_nosac(:,i);
    bar((p2-p1)./p2);
    title(num2str(sfs(i)));
    set(gca,'XTick',[],'YLim',lims);
    if (i == 2)
        ylabel('Difference in prob. correct (nosac - sac)');
    end
end
set(gca,'XTick',[1 2 3],'XTickLabel',{'Ach','L-M','S'});
disp(['% correct overall: ',num2str(sum(x_sac(:)+x_nosac(:))/sum(n_sac(:)+n_nosac(:)))]);


% Sedna get slightly more corrects on achromatic.  Yikes.
%%
% Section 9.1
% Probablity correct (and sac/cor correlation) on trials in which a saccade
% was made during a specific time window conditional on target choice.

filenames = fnamesFromTxt;
offsets = [.167 .5];   % relative to stim on, when to look for saccades
x_sac = [0 0];
n_sac = [0 0];
x_nosac = [0 0];
n_nosac = [0 0];
ACHROMONLY = 1;
for i = 1:size(filenames,1)
    stro = nex2stro(findfile(filenames(i,:)));
    stro = DTfilterquesttrials(stro,'PaperDefaults');
    
    if isempty(stro.sum.analog.sigid) || isempty(stro.trial)
        continue;
    end
    ntrials = size(stro.trial,1);
    whichcol = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'color_dir'));
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_on'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_off'));
    correct = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
    flashside = -1*sign(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_x')));
    flashside = max(flashside, 0); % 1 is to the left
    choice = (correct & flashside) | (~correct & ~flashside);
    choice = abs(choice-1)+1;  % 1 = T1, 2 = T2
    
    sacstats = getSacData(stro);
    close;
    
    % Looping over trials
    for j = 1:ntrials
        if (ACHROMONLY & whichcol(j) ~= 1)
            continue;
        end
        st = sacstats.starttimes{j};
        Lsac = (st > stimon_t(j)+offsets(1)) & (st < stimon_t(j)+offsets(2));
        if (any(Lsac))
            if (correct(j))
                x_sac(choice(j)) = x_sac(choice(j))+1;
            end
            n_sac(choice(j)) = n_sac(choice(j))+1;
        end
        if (~any(Lsac))
            if (correct(j))
                x_nosac(choice(j)) = x_nosac(choice(j))+1;
            end
            n_nosac(choice(j)) = n_nosac(choice(j))+1;
        end
    end
end


for choice = 1:2
    disp(['performance on T',num2str(choice),' choice saccade trials: ',num2str(x_sac(choice)./n_sac(choice))]);
    disp(['performance on T',num2str(choice),' choice no saccade trials: ',num2str(x_nosac(choice)./n_nosac(choice))]);
    [h,p] = equalproptest([x_sac(choice) x_nosac(choice)],[n_sac(choice) n_nosac(choice)],.05);
    disp(['Comparision between these two proportions: ',num2str(p)]);
end
disp([num2str(sum(n_sac)),' saccade trials out of ',num2str(sum(n_sac+n_nosac)),' trials total']);

disp(['prop. trials with a saccade ', num2str(n_sac./n_nosac)]);
disp(['prop error trials ',num2str(1-(x_sac+x_nosac)./(n_sac+n_nosac))])
disp(['percent correct ',num2str((x_sac+x_nosac)./(n_sac+n_nosac))])

%        Cor Inc
% Sac
% no Sac

for choice = 1:2
    conttab = [x_sac(choice) n_sac(choice)-x_sac(choice); x_nosac(choice) n_nosac(choice)-x_nosac(choice)]
    % [p,x2] = chisquarecont(conttab)
    
    data = [repmat([1 1],conttab(1,1),1); repmat([1 0],conttab(1,2),1); repmat([0 1],conttab(2,1),1); repmat([0 0],conttab(2,2),1)];
    corrcoef(data)
end

%%
% Section 9.2
% Probability correct on trials in which a saccade was made during a specific
% time window conditional on color - only looking at experiments in which
% percent correct was greater on achromatic trials.  This is an issue for
% Sedna because she tends to have greater percent correct on achromatics.
% (?!)

filenames = fnamesFromTxt;
offsets = [.167 .5];   % relative to stim on, when to look for saccades 
x_sac = [0 0];
n_sac = [0 0];
x_nosac = [0 0];
n_nosac = [0 0];
colordirs = mkbasis([1 1 1; 1 -1 0]');
for i = 1:size(filenames,1)
    stro = nex2stro(findfile(filenames(i,:)));
    stro = DTfilterquesttrials(stro,'PaperDefaults');
  
    if isempty(stro.sum.analog.sigid) || isempty(stro.trial)
        continue;
    end
    ntrials = size(stro.trial,1);
    whichcol = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'color_dir'));
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_on'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_off'));
    correct = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
    lms = stro.sum.exptParams.RF_colors;
    lms = reshape(lms,[3,size(lms,1)/3]);
    lms(:,all(lms == 0)) = [];
    lms = mkbasis(lms);

    tmp = [];
    for j = 1:size(colordirs,2)
        coloridx = find(abs(lms(:,j)'*colordirs) > .99);  % Avoiding roundoff errors
        L = logical(whichcol == coloridx);
        sum(correct(L))
        sum(L)
        
        tmp = [tmp; sum(correct(L))./sum(L)];
    end
    if (tmp(1) > tmp(2))
        disp('Skipping a file');
        continue;
    end
    sacstats = getSacData(stro);
    close;
    % Looping over trials
    for j = 1:ntrials
        coloridx = find(abs(lms(:,whichcol(j))'*colordirs) > .99);  % Avoiding roundoff errors
        if (isempty(coloridx))
            continue;
        end
        st = sacstats.starttimes{j};
        Lsac = (st > stimon_t(j)+offsets(1)) & (st < stimon_t(j)+offsets(2));
        if (any(Lsac))
            if (correct(j))
                x_sac(coloridx) = x_sac(coloridx)+1;
            end
            n_sac(coloridx) = n_sac(coloridx)+1;
        end
        if (~any(Lsac))
            if (correct(j))
                x_nosac(coloridx) = x_nosac(coloridx)+1;
            end
            n_nosac(coloridx) = n_nosac(coloridx)+1;
        end
    end
end
% Verifying higher proportion correct in L-M
(x_sac+x_nosac)./(n_sac+n_nosac)

% Sign of the effect is correct.
(x_nosac./n_nosac)-(x_sac./n_sac)

% effect of microsaccades on achromatic trials not significant
[h,p] = equalproptest([x_sac(1) x_nosac(1)],[n_sac(1) n_nosac(1)],0.05);
[h,p] = equalproptest([x_sac(2) x_nosac(2)],[n_sac(2) n_nosac(2)],0.05)



%%
% Section 10
% Comparing fixational eyemovements made during DTspot and white noise.

WNfilelists = {'N:\NexFiles\nexfilelists\Greg\ColorOpponent.txt','N:\NexFiles\nexfilelists\Greg\Lum.txt'};
WNfilenames = {};
WNspikenums = [];
for i = 1:length(WNfilelists)
    [tmpnames, tmpnums] = fnamesFromTxt2(WNfilelists{i});
    WNfilenames = cat(1,tmpnames, WNfilenames);
    WNspikenums = cat(1, tmpnums, WNspikenums);
end

DTfilelists = {'N:\NexFiles\nexfilelists\Greg\KaliQuestLrg.txt'};
DTfilenames = {};
DTspikenums = [];
for i = 1:length(DTfilelists)
    [tmpnames, tmpnums] = fnamesFromTxt2(DTfilelists{i});
    DTfilenames = cat(1,tmpnames, DTfilenames);
    DTspikenums = cat(1,tmpnums, DTspikenums);
end

amplitudes = [];
directions = [];
peakv = [];
pathlengths = [];
durations = [];
whichlist = [];

data = [];
for a = 1:2
    if (a == 1)
        filenames = WNfilenames;
        spikenums = WNspikenums;
    else
        filenames = DTfilenames;
        spikenums = DTspikenums;
    end
    
    for i = 1:size(filenames,1)
        stro = {};
        for f = 1:size(filenames{i},1)
            tmpstro = nex2stro(findfile(char(filenames{i}(f))));
            if (isempty(stro))
                stro = tmpstro;
            else
                stro = strocat(stro, tmpstro);
            end
        end
        if (a == 1)
            stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_on'));
            stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'all_off')); 
            % Need to strip out replay trials
        elseif (a == 2)
            stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_on'));
            stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_off'));
        end
        sacstats = getSacData(stro);
        close;
        
        nsac = 0;
        nsec = 0;
        for j = 1:size(stro.trial,1)
            st = sacstats.starttimes{j};
            Lsac = (sacstats.starttimes{j} > stimon_t(j)) & (sacstats.endtimes{j} < stimoff_t(j)-0.1);
            
            if any(Lsac)
                lastsacidx = find(Lsac,1,'last'); % Omitting final sac if long
                if (sacstats.amplitudes{j}(lastsacidx) > 2)
                    disp('Omitting last saccade');
                    sacstats.amplitudes{j}(lastsacidx) = [];
                    sacstats.directions{j}(lastsacidx) = [];
                    sacstats.peakv{j}(lastsacidx) = [];
                    sacstats.pathlengths{j}(lastsacidx) = [];
                    sacstats.durations{j}(lastsacidx) = [];
                    Lsac(lastsacidx) = 0;
                end                    
                if any(sacstats.amplitudes{j}(Lsac) > 2)
                    keyboard  % This only happens once on a pathological trial
                end
                amplitudes = [amplitudes; sacstats.amplitudes{j}(Lsac)];
                directions = [directions; sacstats.directions{j}(Lsac)];
                peakv = [peakv; sacstats.peakv{j}(Lsac)];
                pathlengths = [pathlengths; sacstats.pathlengths{j}(Lsac)];
                durations = [durations; sacstats.durations{j}(Lsac)];
                whichlist = [whichlist; repmat(a,sum(Lsac),1)];
            end
            nsac = nsac+sum(Lsac);
            nsec = nsec+sum(stimoff_t(j)-stimon_t(j));
        end
        data = [data; nsac, nsec, a];
    end
end

L = whichlist == 1; % 1 = white noise
figure;
% Amplitude comparison
bins = linspace(0,2,20);
[n1,x] = hist(amplitudes(L),bins);
[n2,x] = hist(amplitudes(~L),bins);
subplot(3,1,1); hold on;
plot(bins, n1./sum(n1),'r-');  % white noise in red
plot(bins, n2./sum(n2),'b-');  % DTspot in blue
[h,p] = ttest2(amplitudes(L), amplitudes(~L))

% Direction comparison
subplot(3,1,2); hold on;
[rho, theta]= hist(directions(L),20);
polar([theta theta(1)],[rho rho(1)]./sum(rho),'r-');
[rho, theta]= hist(directions(~L),20);
polar([theta theta(1)],[rho rho(1)]./sum(rho),'b-');
axis square
[h,p] = WatsonU2Test(directions(L),directions(~L))

% Frequency comparison
WNsacrate = data(data(:,3) == 1,1)./data(data(:,3) == 1,2)
DTsacrate = data(data(:,3) == 2,1)./data(data(:,3) == 2,2)
bins = linspace(0,2,20);
[n1,x] = hist(WNsacrate,bins);
[n2,x] = hist(DTsacrate,bins);
subplot(3,1,3); hold on;
plot(bins, n1./sum(n1),'r-');  % white noise in red
plot(bins, n2./sum(n2),'b-');  % DTspot in blue
[h,p] = ttest2(WNsacrate, DTsacrate)

%% 
% Section 10
% Receptive field eccentricities from white noise experiments.

filelists = {'N:\NexFiles\nexfilelists\Greg\ColorOpponent.txt','N:\NexFiles\nexfilelists\Greg\Lum.txt'};
filenames = {};
spikenums = [];
L = [];
for i = 1:length(WNfilelists)
    [tmpnames, tmpnums] = fnamesFromTxt2(filelists{i});
    filenames = cat(1,tmpnames, filenames);
    spikenums = cat(1,tmpnums, spikenums);
    L = cat(1,repmat(i,length(tmpnames),1),L);
end

data = [];
for i = 1:size(filenames,1)
    stro = {};
    for f = 1:size(filenames{i},1)
        tmpstro = nex2stro(findfile(char(filenames{i}(f))));
        if (isempty(stro))
            stro = tmpstro;
        else
            stro = strocat(stro, tmpstro);
        end
    end
    data = [data; stro.sum.exptParams.rf_x  stro.sum.exptParams.rf_y];
end
figure; axes; hold on;
cols = [1 0 0; 0 0 0];
for i = unique(L)'
    h = plot(data(L == i,1),data(L == i,2),'*');
    set(h,'color',cols(i,:));
end
set(gca,'XLim',[-100 100],'YLim',[-100 100])
plot(-50,-35,'y*')
mean(data)

%%
% Section 11
% Probablity correct on trials in which a saccade was made as a function of 
% time, color, spatial frequency, and whether the saccade had an up or downward
% component.  This will help answer the question, is there something about
% the saccade-induced change in stimulus speed that makes it easier or
% harder to detect.

filenames = fnamesFromTxt;
colordirs = mkbasis([1 1 1; 1 -1 0; 0 0 1]');
periods = [16 58 206];
%binwidth = .1; %sec
%bins = [-.2:binwidth:.8];
binwidth = .1; %sec
bins = 0.7170;  % For asking how often an upward microsaccade precedes an upward choice

x = zeros(length(bins),3,3,2);  % bins, color, sf, EM (up, down)
n = zeros(length(bins),3,3,2);
m = zeros(length(bins),3,3,2);
% m = 1:T1 choices. x = 1:correct.

for i = 1:size(filenames,1)
    stro = nex2stro(findfile(filenames(i,:)));
    stro = DTfilterquesttrials(stro,'PaperDefaults');
  
    if isempty(stro.sum.analog.sigid) || isempty(stro.trial)
        continue;
    end
    ntrials = size(stro.trial,1);
    lms = stro.sum.exptParams.RF_colors;
    lms = reshape(lms,[3,size(lms,1)/3]);
    lms(:,all(lms == 0)) = [];
    lms = mkbasis(lms);
    whichcol = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'color_dir'));
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_on'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_off'));
    targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_on'));
    spatialPeriods = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'gabor_lambda'));
    correct = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
    flashside = -1*sign(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_x')));
    flashside = max(flashside, 0); % 1 is to the left
    choice = (correct & flashside) | (~correct & ~flashside);
    sacstats = getSacData(stro);
    close;
  
    
     % Looping over trials
    for k = 1:length(bins)
        for j = 1:ntrials
            coloridx = find(abs(lms(:,whichcol(j))'*colordirs) > .99);  % Avoiding roundoff errors
            if (isempty(coloridx))
                continue;
            end
            periodidx = find(spatialPeriods(j) == periods);
            if (isempty(periodidx))
                continue;
            end      
            dirs = sacstats.directions{j};
            amps = sacstats.amplitudes{j};
            st = sacstats.starttimes{j};
            Lsac = logical((st > stimon_t(j)+bins(k)-binwidth/2) & (st < stimon_t(j)+bins(k)+binwidth/2));
            diridx = 1+(mod(dirs,2*pi) > pi); % upward = 1, downward = 2
            if (any(Lsac))
                if (correct(j))
                    x(k,coloridx,periodidx,diridx(Lsac)) = x(k,coloridx,periodidx,diridx(Lsac))+1;
                end
                if (choice(j))
                    m(k,coloridx,periodidx,diridx(Lsac)) = m(k,coloridx,periodidx,diridx(Lsac))+1;
                end
                n(k,coloridx,periodidx,diridx(Lsac)) = n(k,coloridx,periodidx,diridx(Lsac))+1;
            end
        end
    end
end

% first collapsing across stimuli and just looking at eye movement
% direction

p1 = sum(squeeze(sum(x(:,:,:,1),2)),2)./sum(squeeze(sum(n(:,:,:,1),2)),2);
p2 = sum(squeeze(sum(x(:,:,:,2),2)),2)./sum(squeeze(sum(n(:,:,:,2),2)),2);

% What's the probability that an upward microsaccade precedes an upward
% choice?
nn = squeeze(sum(sum(n,2),3))
mm = nn-squeeze(sum(sum(m,2),3));
p = mm./nn
p(1)./p(2)
figure;
axes; hold on;
plot(bins, p1,'g-');
plot(bins, p2,'r-');
legend({'Upward saccade','Downward saccade'});

% Now breaking it down by color (only)
p = squeeze(sum(sum(x,3),4))./squeeze(sum(sum(n,3),4));
figure; set(gcf,'DefaultAxesColorOrder',[0 0 0; 1 0 0; 0 0 1])
axes; hold on;
plot(bins, p,'-');

% Now breaking it down by spatial frequency and eye movement direction
p1 = squeeze(sum(x(:,:,:,1),2))./squeeze(sum(n(:,:,:,1),2));
p2 = squeeze(sum(x(:,:,:,2),2))./squeeze(sum(n(:,:,:,2),2));

figure;
axes; hold on;
plot(bins, p2,'-');
plot(bins, p1,':');

legend(num2str(stro.sum.exptParams.pixperdeg ./ periods'));

% Now breaking it down by color (and eye movement direction)
p1 = squeeze(sum(x(:,:,:,1),3))./squeeze(sum(n(:,:,:,1),3));
p2 = squeeze(sum(x(:,:,:,2),3))./squeeze(sum(n(:,:,:,2),3));

figure; set(gcf,'DefaultAxesColorOrder',[0 0 0; 1 0 0; 0 0 1])
axes; hold on;
plot(bins, p2,'-');
plot(bins, p1,':');



% How much do these microsaccades change the speed of the stimulus?
sfs = stro.sum.exptParams.pixperdeg ./ periods';
tf = mode(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'gabor_speed')));
stimspeed = tf./sfs; % deg/sec
sacspeed = cat(1,sacstats.amplitudes{:})./cat(1,sacstats.durations{:}); % deg/sec
figure; axes; hold on;
plot(sacspeed,cat(1,sacstats.peakv{:}),'k.')
plot([0 1000],[0 1000]);
xlabel('average speed (deg/sec)'); ylabel('peak speed (deg/sec)');
% The speed of the stimulus is lower than the speed of even a small
% saccade. 
sum(sacspeed > max(stimspeed))./length(sacspeed)
% But then again, even a fast saccade nearly orthogonal to the drift direction
% will have a slow component in the direction of the drift.

%%
% Section 12
% Psychometric thresholds on trials with and without a fixational saccade.
% Permutation test for significant effect of saccade.

filenames = fnamesFromTxt;
colordirs = mkbasis([1 1 1; 1 -1 0; 0 0 1]');
periods = [206 58 16];  % Low to high spatial frequency
offset = [.167 .5];
AMPLITUDERANGE = [0 .4];  % Only considering microsaccades within this range
    
data = cell(length(periods), size(c,1));
for i = 1:size(filenames,1)
    stro = nex2stro(findfile(filenames(i,:)));
    stro = DTfilterquesttrials(stro,'PaperDefaults');
  
    if isempty(stro.sum.analog.sigid) || isempty(stro.trial)
        continue;
    end
    ntrials = size(stro.trial,1);
    lms = stro.sum.exptParams.RF_colors;
    lms = reshape(lms,[3,size(lms,1)/3]);
    lms(:,all(lms == 0)) = [];
    lms = mkbasis(lms);
    whichcol = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'color_dir'));
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_on'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_off'));
    targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_on'));
    spatialPeriods = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'gabor_lambda'));
    correct = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
    flashR = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_R'));
    flashG = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_G'));
    flashB = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_B'));
    bkgndrgb = [stro.sum.exptParams.bkgnd_r, stro.sum.exptParams.bkgnd_g, stro.sum.exptParams.bkgnd_b];
    M = reshape(stro.sum.exptParams.m_mtx, 3, 3);
    bkgndlms = M * bkgndrgb';
    x = 0:255; %the normal range of the gamma look up table
    xx = linspace(0, 255, 2^16); %the desired quantization of the gammaTable
    g1 = reshape(stro.sum.exptParams.gamma_table, 256, 3);
    gammaTable = [spline(x, g1(:,1), xx)', spline(x, g1(:,2), xx)', spline(x, g1(:,3), xx)'];

    sacstats = getSacData(stro);
    close;
  
    saccadeoccurred = nan*ones(ntrials,1);
    for j = 1:ntrials
        st = sacstats.starttimes{j};
        Lsac = (st > stimon_t(j)+offset(1)) & (st < stimon_t(j)+offset(2));
        if (~isempty(AMPLITUDERANGE))
           Lsac = Lsac & (sacstats.amplitudes{j} > AMPLITUDERANGE(1)) & (sacstats.amplitudes{j} < AMPLITUDERANGE(2));
        end
        if (any(Lsac))
            saccadeoccurred(j) = 1;
        else
            saccadeoccurred(j) = 0;
        end
    end
    for a = 1:size(colordirs,1)
        for b = 1:length(periods)
            coloridx = find(abs(lms'*colordirs(:,a)) > .99);  % Avoiding roundoff errors
            if (isempty(coloridx))
                continue;
            end
            L = whichcol == coloridx;
            L = L & spatialPeriods == periods(b);
            L = L & stro.trial(:,strcmp(stro.sum.trialFields(1,:),'numframes')) > 0;
            % Above, need to eliminate 0 frame trials.  What are these?
            RGB =[flashR(L), flashG(L), flashB(L)];
            rgb = [gammaTable(RGB(:,1)+1, 1), gammaTable(RGB(:,2)+1, 2), gammaTable(RGB(:,3)+1, 3)];
            cc = ((M*rgb')'-repmat(bkgndlms',sum(L),1)) ./ repmat(bkgndlms',sum(L),1);
            ccnorms = sqrt(sum(cc.^2,2));
            if (any(ccnorms > 1.5))
                keyboard
            end
            tmp = [saccadeoccurred(L), correct(L), ccnorms];
            data{b,a} = cat(1,data{b,a},tmp);
        end
    end
end

colordirs = mkbasis([1 1 1; 1 -1 0]');

% How many trials of each type?
% Ten times fewer microsaccade trials across the board
figure
for a = 1:size(colordirs,2)
    for b = 1:length(periods)
        subplot(size(data,1),size(data,2),sub2ind(size(data),a, b)); hold on;
        n_sac = sum(data{b,a}(:,1) == 1);
        n_nosac = sum(data{b,a}(:,1) == 0);
        bar([0 1],[n_nosac n_sac],'k');
    end
end


% binning
h = figure
plotcols = ['k','m'];
for a = 1:size(colordirs,2)
    for b = 1:length(periods)
        subplot(size(data,1),size(data,2),sub2ind(size(data),a, b)); hold on;
        for c = 1:-1:0
            Lcor = logical(data{b,a}(:,2));
            Lsac = logical(data{b,a}(:,1) == c);
            bins = logspace(log10(min(data{b,a}(:,3))),log10(max(data{b,a}(:,3))),5);
            [n1,x1] = hist(data{b,a}(Lcor&Lsac,3),bins);
            [n2,x2] = hist(data{b,a}(~Lcor&Lsac,3),bins);
            h = plot(bins, n1./(n1+n2),'o');
            set(h,'markerfacecolor',plotcols(c+1));
        end
        set(gca,'XScale','log');
    end
end

% Doing a significance test using the asymptotic -2*log likelihood thing
p = [];
for a = 1:size(colordirs,2)
    for b = 1:length(periods)
        tmp = data{b,a};
        Lsac = logical(tmp(:,1));
        [fittedparams, success(1)] = weibullFitGH(tmp(:,3), tmp(:,2), 'sse', [mean(tmp(:,3)), 1]);
        [fittedparams_full, success(2)] = weibullFitGH([tmp(:,3), Lsac], tmp(:,2), 'mle', [fittedparams(1),fittedparams(2)]);
        [fittedparams_nested, success(2)] = weibullFitGH(tmp(:,3), tmp(:,2), 'mle', [fittedparams(1),fittedparams(2)]);
        pred_nested = 1-.5*exp(-(tmp(:,3)./fittedparams_nested(1)).^fittedparams_nested(2));
        pred_full = 1-.5*exp(-(tmp(:,3)./(fittedparams_full(1)+tmp(:,1)*fittedparams_full(4))).^fittedparams_full(2));
        dev_nested = sum(log(sum([tmp(:,2).*pred_nested ~tmp(:,2).*(1-pred_nested)],2)));
        dev_full = sum(log(sum([tmp(:,2).*pred_full ~tmp(:,2).*(1-pred_full)],2)));
        teststat = -2*(dev_nested-dev_full);
        p(a,b) = 1-chi2cdf(teststat,1);
        out(b,a,1) = fittedparams_full(1)+fittedparams_full(4);  % Threshold with saccade would be high
        out(b,a,2) = fittedparams_full(1);  % Threshold with saccade would be high
        out(b,a,3) = fittedparams_nested(1); % Threshold without saccade would be low
    end
end
% Now in "log units"
effectsize = log10(out(:,:,1)./out(:,:,2));
bins = [0:.05:.8];
figure; axes; hold on;
[x,n] = hist(effectsize(:),bins);
bar(n,x,'k')
[x,n] = hist(reshape(effectsize(p < 0.01),sum(p(:)<0.01),1),bins);
h = bar(n,x,'r');

% Doing a significance test using a permutation test
% Very similar results to chi-squared test
niter = 200;
p = zeros(size(colordirs,1),length(periods));
se = zeros(size(colordirs,1),length(periods));
for a = 1:size(colordirs,1)
    for b = 1:length(periods)
        teststat = zeros(niter,1);
        statforse = zeros(niter,1);
        for c = 1:niter
            tmp = data{b,a};
            Lsac = logical(tmp(:,1));
            if (c > 1)
               Lsac = Lsac(randperm(length(Lsac)));
            end
            [fittedparams, success(1)] = weibullFitGH(tmp(:,3), tmp(:,2), 'sse', [mean(tmp(:,3)), 1]);
            [fittedparams_full, success(2)] = weibullFitGH([tmp(:,3), Lsac], tmp(:,2), 'mle', [fittedparams(1),fittedparams(2)]);
            [fittedparams_nested, success(2)] = weibullFitGH(tmp(:,3), tmp(:,2), 'mle', [fittedparams(1),fittedparams(2)]);
            pred_nested = 1-.5*exp(-(tmp(:,3)./fittedparams_nested(1)).^fittedparams_nested(2));
            pred_full = 1-.5*exp(-(tmp(:,3)./(fittedparams_full(1)+tmp(:,1)*fittedparams_full(4))).^fittedparams_full(2));
            dev_nested = sum(log(sum([tmp(:,2).*pred_nested ~tmp(:,2).*(1-pred_nested)],2)));
            dev_full = sum(log(sum([tmp(:,2).*pred_full ~tmp(:,2).*(1-pred_full)],2)));
            teststat(c) = fittedparams_full(4);
            statforse(c) = log10((fittedparams_full(1)+fittedparams_full(4))./fittedparams_full(1));
            if (c == 1)
                out(b,a,1) = fittedparams_full(1)+fittedparams_full(4);  % Threshold with saccade would be high
                out(b,a,2) = fittedparams_nested(1); % Threshold without saccade would be low
            end
        end
        p(a,b) = 1-sum(teststat(1) > teststat(2:end))./size(teststat,1);
        se(a,b) = std(statforse)
    end
end

% 
% % Comparing thresholds with fitted cumulative weibulls 
% % and doing a permutation test to assess significance.
% niter = 20;
% outsac = nan*ones(length(periods), size(colordirs,1), niter);
% outnosac = nan*ones(length(periods), size(colordirs,1), niter);
% 
% for a = 1:size(colordirs,1)
%     for b = 1:length(periods)
%         for i = 1:niter
%             tmp = data{b,a};
%             Lsac = logical(tmp(:,1));
%             if (i > 1)
%                 Lsac = Lsac(randperm(length(Lsac)));
%             end
%             [aSSEnosac, bSSE, gSSE] = weibullFit(tmp(~Lsac,3), tmp(~Lsac,2), 'sse', [mean(tmp(~Lsac,3)), 1]);
%             [aMLnosac, bML, gML, success] = weibullFit(tmp(~Lsac,3), [tmp(~Lsac,2) ~tmp(~Lsac,2)], 'mle', [aSSEnosac, bSSE]);
%             
%             [aSSEsac, bSSE, gSSE] = weibullFit(tmp(Lsac,3), tmp(Lsac,2), 'sse', [mean(tmp(Lsac,3)), 1]);
%             [aMLsac, bML, gML, success] = weibullFit(tmp(Lsac,3), [tmp(Lsac,2) ~tmp(Lsac,2)], 'mle', [aSSEsac, bSSE]);
%             if (~success)
%                 keyboard
%             end
%             outsac(b,a,i) = aMLsac;
%             outnosac(b,a,i) = aMLnosac;
%           i
%         end
%     end
%     b
% end
% 
% out = outsac-outnosac;
% threshchange = 1-outnosac(:,:,1)./outsac(:,:,1)
% p = [];
% for a = 1:size(colordirs,1)
%     for b = 1:length(periods)
%         p(a,b) = sum(out(b,a,:) > out(b,a,1))./size(out,3);
%     end
% end
% % 
% Comparing thresholds with fitted lines
niter = 100;
outsac = nan*ones(length(periods), size(colordirs,1), niter);
outnosac = nan*ones(length(periods), size(colordirs,1), niter);
for a = 1:size(colordirs,1)
    for b = 1:length(periods)
        for i = 1:niter
            tmp = data{b,a};
            Lsac = logical(tmp(:,1));
            Lcor = logical(tmp(:,2));
            if (i > 1)
                Lsac = Lsac(randperm(length(Lsac)));
            end
            [bs,bint,r,rint,stats]=regress(tmp(:,2),[ones(size(tmp,1),1) Lsac log10(tmp(:,3))]);
            outsac(b,a,i) = 10^(.75-bs(1)-bs(2))/bs(3);
            outnosac(b,a,i) = 10^(.75-bs(1))/bs(3);
            if (i == 1)
                subplot(size(data,1),size(data,2),sub2ind(size(data),a, b)); hold on;
                plot(10.^linspace(-4,0,2),bs(1)+bs(3)*linspace(-4,0,2))
            end
        end
    end
    bs
end
out = outsac-outnosac;
threshchange = 1-outnosac(:,:,1)./outsac(:,:,1);
p = [];
for a = 1:size(colordirs,1)
    for b = 1:length(periods)
        p(a,b) = sum(out(b,a,:) > out(b,a,1))./size(out,3);
    end
end
% Now in "log units"
log10(outsac(:,:,1)./outnosac(:,:,1))
% 
%
% threshchange = [];
% p = [];
% % Comparing thresholds just looking at medians (means)
% for a = 1:size(colordirs,1)
%     for b = 1:length(periods)
%         tmp = data{b,a};
%         Lsac = logical(tmp(:,1));
%         Lcorrect = logical(tmp(:,2));
%         % Threshold with a saccade should be low, threshold without a
%         % saccade should be high, so ratio should be smaller than 1
%         % Positive threshchanges agree with prediction
%         threshchange(a,b) = 1-(mean(tmp(Lsac,3))./mean(tmp(~Lsac,3)));
%         p(a,b) = ranksum(tmp(~Lsac,3), tmp(Lsac,3));
%     end
% end
% p, threshchange: each row is a color, each column a SF


%%
% Section 12.1
% Looking for a relationship between contrast and microsaccade occurrence

filenames = fnamesFromTxt;
colordirs = mkbasis([1 1 1; 1 -1 0]');
periods = [206 58 16];  % Low to high spatial frequency
offset = [.167 .5];
data = cell(length(periods), size(colordirs,2));

for fileidx = 1:size(filenames,1)
    stro = nex2stro(findfile(filenames(fileidx,:)));
    stro = DTfilterquesttrials(stro,'PaperDefaults');
    if isempty(stro.sum.analog.sigid) || isempty(stro.trial)
        continue;
    end
    ntrials = size(stro.trial,1);
    lms = stro.sum.exptParams.RF_colors;
    lms = reshape(lms,[3,size(lms,1)/3]);
    lms(:,all(lms == 0)) = [];
    lms = mkbasis(lms);
    whichcol = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'color_dir'));
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_on'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_off'));
    targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_on'));
    spatialPeriods = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'gabor_lambda'));
    correct = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
    flashR = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_R'));
    flashG = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_G'));
    flashB = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_B'));
    bkgndrgb = [stro.sum.exptParams.bkgnd_r, stro.sum.exptParams.bkgnd_g, stro.sum.exptParams.bkgnd_b];
    M = reshape(stro.sum.exptParams.m_mtx, 3, 3);
    bkgndlms = M * bkgndrgb';
    x = 0:255; %the normal range of the gamma look up table
    xx = linspace(0, 255, 2^16); %the desired quantization of the gammaTable
    g1 = reshape(stro.sum.exptParams.gamma_table, 256, 3);
    gammaTable = [spline(x, g1(:,1), xx)', spline(x, g1(:,2), xx)', spline(x, g1(:,3), xx)'];

    sacstats = getSacData(stro);
    close;
  
    saccadeoccurred = nan*ones(ntrials,1);
    for j = 1:ntrials
        st = sacstats.starttimes{j};
        Lsac = (st > stimon_t(j)+offset(1)) & (st < stimon_t(j)+offset(2));
        if (any(Lsac))
            saccadeoccurred(j) = 1;
        else
            saccadeoccurred(j) = 0;
        end
    end
    for a = 1:size(colordirs,2)
        for b = 1:length(periods)
            coloridx = find(abs(lms'*colordirs(:,a)) > .99);  % Avoiding roundoff errors
            if (isempty(coloridx))
                continue
            end
            L = whichcol == coloridx;
            L = L & spatialPeriods == periods(b);
            
            RGB =[flashR(L), flashG(L), flashB(L)];
            rgb = [gammaTable(RGB(:,1)+1, 1), gammaTable(RGB(:,2)+1, 2), gammaTable(RGB(:,3)+1, 3)];
            cc = ((M*rgb')'-repmat(bkgndlms',sum(L),1)) ./ repmat(bkgndlms',sum(L),1);
            ccnorms = sqrt(sum(cc.^2,2));
            tmp = [saccadeoccurred(L), correct(L), ccnorms];
            data{b,a} = cat(1,data{b,a},tmp);
        end
    end
end

% Do microsaccades occur more frequently on low contrast trials?
for i = 1:size(data,1)
    for j = 1:size(data,2)
        Lsac = logical(data{i,j}(:,1));
        Lcor = logical(data{i,j}(:,2));
        % sanity check
        [h,p(i,j)] = ttest2(data{i,j}(Lcor,3), data{i,j}(~Lcor,3));
        effect(i,j) = mean(data{i,j}(Lcor,3))-mean(data{i,j}(~Lcor,3));
        
        % The real analysis
        [h,p(i,j)] = ttest2(data{i,j}(Lsac,3), data{i,j}(~Lsac,3));
        effect(i,j) = mean(data{i,j}(Lsac,3))-mean(data{i,j}(~Lsac,3));
        % Large effects mean higher contrasts with saccades
    end
end
effect
p

% Does the relationship between microsaccades and correct depend on the
% contrast of the stimulus (microsaccadic inhibition suggests that high
% contrast stimuli should be more effective at suppressing microsaccades).

for i = 1:size(data,1)
    for j = 1:size(data,2)
        Lsac = logical(data{i,j}(:,1));
        Lcor = logical(data{i,j}(:,2));
        contrast = data{i,j}(:,3);
        L = contrast >= prctile(contrast,75);
        p_sac = sum(Lcor(L&Lsac))./sum(L&Lsac);
        p_nosac = sum(Lcor(L&~Lsac))./sum(L&~Lsac);
        relsup_hi(i,j) = p_nosac - p_sac;

        p_sac = sum(Lcor(~L&Lsac))./sum(~L&Lsac);
        p_nosac = sum(Lcor(~L&~Lsac))./sum(~L&~Lsac);
        relsup_lo(i,j) = p_nosac - p_sac;
    end
end
relsup_hi
relsup_lo

%%
% Section 12.2
% Psychometric thresholds on trials with and without a fixational saccade.
% Omitting from the analysis trials with a saccade outside of the tolerance
% range.

filenames = fnamesFromTxt;
colordirs = mkbasis([1 1 1; 1 -1 0; 0 0 1]');
periods = [206];  % Low to high spatial frequency
offset = [.167 .5];
AMPLITUDERANGE = [0 .4];  % Only considering microsaccades within this range

data = cell(length(periods), size(colordirs,1));
for i = 1:size(filenames,1)
    stro = nex2stro(findfile(filenames(i,:)));
    stro = DTfilterquesttrials(stro,'PaperDefaults');
  
    if isempty(stro.sum.analog.sigid) || isempty(stro.trial)
        continue;
    end
    ntrials = size(stro.trial,1);
    lms = stro.sum.exptParams.RF_colors;
    lms = reshape(lms,[3,size(lms,1)/3]);
    lms(:,all(lms == 0)) = [];
    lms = mkbasis(lms);
    whichcol = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'color_dir'));
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_on'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_off'));
    targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_on'));
    spatialPeriods = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'gabor_lambda'));
    correct = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
    flashR = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_R'));
    flashG = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_G'));
    flashB = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_B'));
    bkgndrgb = [stro.sum.exptParams.bkgnd_r, stro.sum.exptParams.bkgnd_g, stro.sum.exptParams.bkgnd_b];
    M = reshape(stro.sum.exptParams.m_mtx, 3, 3);
    bkgndlms = M * bkgndrgb';
    x = 0:255; %the normal range of the gamma look up table
    xx = linspace(0, 255, 2^16); %the desired quantization of the gammaTable
    g1 = reshape(stro.sum.exptParams.gamma_table, 256, 3);
    gammaTable = [spline(x, g1(:,1), xx)', spline(x, g1(:,2), xx)', spline(x, g1(:,3), xx)'];

    sacstats = getSacData(stro);
    close;
  
    saccadeoccurred = nan*ones(ntrials,1);
    withinamprange = nan*ones(ntrials,1);
    for j = 1:ntrials
        st = sacstats.starttimes{j};
        Lsac = (st > stimon_t(j)+offset(1)) & (st < stimon_t(j)+offset(2));
        Lamp = (sacstats.amplitudes{j}(Lsac) > AMPLITUDERANGE(1)) & (sacstats.amplitudes{j}(Lsac) < AMPLITUDERANGE(2));
        if (any(Lsac))
            saccadeoccurred(j) = 1;
        else
            saccadeoccurred(j) = 0;
        end
        if (any(Lamp))
            withinamprange(j) = 1;
        else
            withinamprange(j) = 0;
        end
    end
    for a = 1:size(colordirs,1)
        for b = 1:length(periods)
            coloridx = find(abs(lms'*colordirs(:,a)) > .99);  % Avoiding roundoff errors
            if (isempty(coloridx))
                continue;
            end
            L = whichcol == coloridx;
            L = L & spatialPeriods == periods(b);
            L = L & stro.trial(:,strcmp(stro.sum.trialFields(1,:),'numframes')) > 0;
            % Above, need to eliminate 0 frame trials.  What are these?
            RGB =[flashR(L), flashG(L), flashB(L)];
            rgb = [gammaTable(RGB(:,1)+1, 1), gammaTable(RGB(:,2)+1, 2), gammaTable(RGB(:,3)+1, 3)];
            cc = ((M*rgb')'-repmat(bkgndlms',sum(L),1)) ./ repmat(bkgndlms',sum(L),1);
            ccnorms = sqrt(sum(cc.^2,2));
            if (any(ccnorms > 1.5))
                keyboard
            end
            tmp = [saccadeoccurred(L), correct(L), ccnorms, withinamprange(L)];
            data{b,a} = cat(1,data{b,a},tmp);
        end
    end
end

% Doing a significance test using the asymptotic -2*log likelihood thing
p = [];
for a = 1:size(colordirs,2)
    for b = 1:length(periods)
        tmp = data{b,a};
        tmp(logical(tmp(:,1)) & ~logical(tmp(:,4)),:) = [];  % omitting saccade trials with amplitudes outside of the range
        Lsac = logical(tmp(:,1)); 
        [fittedparams, success(1)] = weibullFitGH(tmp(:,3), tmp(:,2), 'sse', [mean(tmp(:,3)), 1]);
        [fittedparams_full, success(2)] = weibullFitGH([tmp(:,3), Lsac], tmp(:,2), 'mle', [fittedparams(1),fittedparams(2)]);
        [fittedparams_nested, success(2)] = weibullFitGH(tmp(:,3), tmp(:,2), 'mle', [fittedparams(1),fittedparams(2)]);
        pred_nested = 1-.5*exp(-(tmp(:,3)./fittedparams_nested(1)).^fittedparams_nested(2));
        pred_full = 1-.5*exp(-(tmp(:,3)./(fittedparams_full(1)+tmp(:,1)*fittedparams_full(4))).^fittedparams_full(2));
        dev_nested = sum(log(sum([tmp(:,2).*pred_nested ~tmp(:,2).*(1-pred_nested)],2)));
        dev_full = sum(log(sum([tmp(:,2).*pred_full ~tmp(:,2).*(1-pred_full)],2)));
        teststat = -2*(dev_nested-dev_full);
        p(a,b) = 1-chi2cdf(teststat,1);
        out(b,a,1) = fittedparams_full(1)+fittedparams_full(4);  % Threshold with saccade would be high
        out(b,a,2) = fittedparams_full(1);  % Threshold with saccade would be high
        out(b,a,3) = fittedparams_nested(1); % Threshold without saccade would be low
    end
end
% Now in "log units"
effectsize = log10(out(:,:,1)./out(:,:,2))
p
out(1,1,[1 2])

%%
% Section 13
% Looking at eye movements synced to frames on.
% This code verifies that the monkey fixated consistently 300 ms before frames on.

filenames = fnamesFromTxt;
data = [];
delta_t = 1/stro.sum.analog.storeRates{1};
for i = 1:size(filenames,1)
    stro = nex2stro(findfile(filenames(i,:)));
    if isempty(stro.sum.analog.sigid) || isempty(stro.trial)
        continue;
    end
    figure; axes; hold on;
    EPstart_t = [stro.ras{:,strcmp(stro.sum.rasterCells,'anlgStartTime')}]';
    frameon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'frame_on'));
    data = [data; min(frameon_t-EPstart_t) mean(frameon_t-EPstart_t)];
    for j = 1:size(stro.trial,1)
        x = stro.ras{j,strcmp(stro.sum.rasterCells,'AD11')};
        y = stro.ras{j,strcmp(stro.sum.rasterCells,'AD12')};
        t = [EPstart_t(j):delta_t:EPstart_t(j)+(length(x)-1)*delta_t];
        t= t-frameon_t(j);
        L = t > .6 | t < -1;
        x(L) = [];
        y(L) = [];
        t(L) = [];
        plot(t,x,'g-');
        plot(t,y,'r-');
    end
    plot([-.3 -.3], [-1 1],'k-');
end

%%
% Section 14
% Looking at saccade triggered PSTHs for a bunch of DTspot files.
% Is there a difference between cells that are sensitive to S and cells
% that are sensitive to L-M with regard to their response to microsaccades?

preTime = .100;
postTime = .200;
binSize = .010;
edges = -preTime:binSize:postTime;
filenames = fnamesFromTxt2;
whichcardinal = []; sfs = [];
data = [];
whichfiles = {};
counter = 1;  % Can't use expt as counter because some files may not have EP data.
for expt = 1:length(filenames)
    disp(filenames{expt}{1})
    stro = nex2stro(findfile(filenames{expt}{1}));
    if isempty(stro.sum.analog.sigid);
        continue
    end
    stro = stripOutGratingTrials(stro);
    sacstats = getSacData(stro);
    close; %getSacData spits out a figure
    contrastLevels = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'cntrst_lev'));
    whichColor = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'color_dir'));
    spatialPeriods = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'gabor_lambda'));
   
    colors = mkbasis(reshape(stro.sum.exptParams.RF_colors,[3 3]));
    if (any(abs([0 0 1]*colors)>.99))
        whichcardinal(counter) = 0;
    elseif (any(abs([1/sqrt(2) -1/sqrt(2) 0]*colors)>.99))
        whichcardinal(counter) =1;
    else
        whichcardinal(counter) =2;
    end
    sfs(counter) = unique(stro.sum.exptParams.pixperdeg ./ spatialPeriods);
    
    correct = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
    flashside = -1*sign(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_x')));
    flashside = max(flashside, 0); % 1 is to the left
    choice = (correct & flashside) | (~correct & ~flashside);
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_on'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_off'));
    
    psth = [];
    for clr = 1:length(unique(whichColor));
        for cntrst = 2:length(unique(contrastLevels));
            tList = (whichColor == clr) & (contrastLevels == cntrst) & (flashside); %only consider in RF trials
            tList = find(tList);
            counts = [];
            for trl = 1:length(tList);
                sacStarts = sacstats.starttimes{tList(trl)};
                sacStarts(sacStarts < stimon_t(tList(trl))) = [];  %only count saccads during the flash
                sacStarts(sacStarts > stimoff_t(tList(trl))) = [];
                for a = 1:length(sacStarts);
                    spikeTimes = stro.ras{tList(trl), 1} - sacStarts(a);
                    spikeTimes(spikeTimes < -preTime) = [];
                    spikeTimes(spikeTimes > postTime) = [];
                    tmp_counts = histc(spikeTimes, edges);
                    tmp_counts(end) = [];
                    counts(end+1, :) = tmp_counts;
                end
            end
            if (~isempty(counts))
                tmp_psth = mean(counts,1) ./ binSize;
                psth(end+1,:) = tmp_psth;
            end
        end
    end
    data(counter, :) = mean(psth);
    whichfiles{counter} = filenames{expt}(1);
    counter = counter+1;
%     figure
%     hold on,
%     bar(edges(1:end-1), mean(psth,1), 1)
%     plot([0 0], [0, max(mean(psth,1))+.1], 'k', 'linewidth', 2)
%     xlim([min(edges)-0.025, max(edges)+0.025])
%     title(sprintf('Cell: %s', filenames{expt}{1}));
end
L = sfs < 5;
meanmat = repmat(mean(data,2),1,size(data,2));
figure;
subplot(2,3,1)
imagesc(data(L,:)-meanmat(L,:));
subplot(2,3,2);
imagesc(data(L & whichcardinal == 0,:)-meanmat(L & whichcardinal == 0,:));  % S
title('S');
subplot(2,3,3);
imagesc(data(L & whichcardinal == 1,:)-meanmat(L & whichcardinal == 1,:));  % L-M
title('L-M');

subplot(2,3,4)
plot(mean(data(L,:)-meanmat(L,:)));
subplot(2,3,5)
plot(mean(data(L & whichcardinal == 0,:)-meanmat(L & whichcardinal == 0,:)));
subplot(2,3,6)
plot(mean(data(L & whichcardinal == 1,:)-meanmat(L & whichcardinal == 1,:)));
% The following neurons tested in DTspot are the "best" by eye (looking at the
% associated gratings data).
% S040809004 (S)
% S042109006 (Pan)
% S051409003 (S)
% S052009003 (L-M)
% S060409004 (L-M)
% S081709002 (Pan)
% S081809001 (S)

% Best L-M cells: S042109007, S052009004, S060409005.
% Clear signs of microsaccade-related suppression in these cells.
% Arguably more than S-cone dominated/tested cells.

%%
% Section 14.1
% Looking more closely specifically at the L-M DTspot files.  Is there an effect of
% microsaccades on spike rate?

BINWIDTH = .01;
timebins = [-.25:BINWIDTH:.5];
filenamelist = 'N:\NexFiles\nexfilelists\Greg\DTEM\DTspotLvsM.txt';
[filenames, spikenums] = fnamesFromTxt2(filenamelist);
data = nan*ones(length(filenames), length(timebins));
for i = 1:size(filenames,1)
    stro = {};
    for f = 1:size(filenames{i},1)
        tmpstro = nex2stro(findfile(char(filenames{i}(f))));
        if (isempty(stro))
            stro = tmpstro;
        else
            stro = strocat(stro, tmpstro);
        end
    end
    
    if (isempty(stro.sum.analog.sigid))
        continue;
    end
    sacstats = getSacData(stro);
    close;
    
    whichColor = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'color_dir'));
    colors = mkbasis(reshape(stro.sum.exptParams.RF_colors,[3 3]));
    Lcol = whichColor == find(abs([1/sqrt(2) -1/sqrt(2) 0]*colors)>.99);
    
    flashside = sign(stro.sum.exptParams.rf_x) == sign(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_x')));
    % flashside == 1 means "in RF"
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'frame_on'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_off'));
    PSTH = zeros(1,length(timebins));
    L = Lcol;
    for j = find(L')
        st = sacstats.starttimes{j};
        Lsac = (st > stimon_t(j)) & (st < stimoff_t(j));
        if any(Lsac)
            spiketimes = [];
            for k = find(Lsac')
                tmp = stro.ras{j,1}-sacstats.starttimes{j}(k);
                spiketimes = [spiketimes; tmp((tmp > timebins(1)-BINWIDTH/2) & (tmp < timebins(end)+BINWIDTH/2))];
            end
            [n,x] = hist(spiketimes, timebins);
            PSTH = PSTH + n./(BINWIDTH*sum(Lsac));
        end
    end
    data(i,:) = PSTH./sum(L);
end
PSTHs = data./repmat(max(data,[],2),1,size(data,2));

figure;
subplot(2,1,1);
imagesc(PSTHs);
subplot(2,1,2);
plot(timebins,mean(PSTHs));
set(gca,'Ylim',[0.1 .55]);
set(gca,'Xlim',[timebins(1) timebins(end)])

%%
% Section 15
% Comparing performance (and microsaccadic suppression) across potential
% S-cone isolation color directions at high and low temporal frequency.

filenames = fnamesFromTxt();
offset = [.167 -.167];
allcolors = [];
x_sac = zeros(size(filenames,1),3);
n_sac = zeros(size(filenames,1),3);
x_nosac = zeros(size(filenames,1),3);
n_nosac = zeros(size(filenames,1),3);
data = [];

for a = 1:size(filenames,1)
    stro = nex2stro(findfile(filenames(a,:)));
    sacstats = getSacData(stro);
    close;
    ntrials = size(stro.trial,1);
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_on'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_off'));
    correct = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
    
    [thresholds, colorDirs, sfs] = DTquestUnpack(stro,'mode');
    % close;
    speed = stro.trial(:,strcmp(stro.sum.trialFields(1,:), 'gabor_speed'));
    whichcol = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'color_dir'));
    
    colors = reshape(stro.sum.exptParams.RF_colors,[3 3]);
    if (isempty(allcolors))
        allcolors = colors;
    else
        if (~all(colors == allcolors))
            error('Disagreement in colors across experiments')
        end
    end
    
    % Looping over trials
    saccadeoccurred = nan*ones(ntrials,1);
    for j = 1:ntrials
        st = sacstats.starttimes{j};
        Lsac = (st > stimon_t(j)+offset(1)) & (st < stimoff_t(j)+offset(2));
        if (any(Lsac))
            saccadeoccurred(j) = 1;
        else
            saccadeoccurred(j) = 0;
        end
    end
    trialtypes = unique(whichcol);
    tmpdata = zeros(size(trialtypes,1),1);
    for j = 1:size(trialtypes,1)
        L = whichcol == trialtypes(j);
        r = corrcoef([correct(L), saccadeoccurred(L)]);
        tmpdata(j) = r(1,2);
        
        if (all(speed == 3))
            x_sac(a,j) = x_sac(a,j)+sum(correct & saccadeoccurred & L);
            n_sac(a,j) = n_sac(a,j)+sum(saccadeoccurred & L);
            x_nosac(a,j) = x_nosac(a,j)+sum(correct & ~saccadeoccurred & L);
            n_nosac(a,j) = n_nosac(a,j)+sum(~saccadeoccurred & L);
        end
    end
    
    if (length(thresholds) == 2)  % only 2 color directions tested
        thresholds = [thresholds; nan];
        tmpdata = [tmpdata; nan];
    end
    if (any(tmpdata == 1))
        keyboard
    end
    data = [data; thresholds', tmpdata', speed(1)];
end
speeds = unique(data(:,7));

% All possible threshold ratios (3 Hz/not 3 Hz, within color directions)
% Too bad these things aren't independent
figure; axes; hold on;
for i = 1:size(allcolors,2)
    L = data(:,7) == 3;
    a = find(L);
    b = find(~L);
    designmat = fullfact([sum(L), sum(~L)]);
    for j = 1:size(designmat,1)
        plot(i,data(a(designmat(j,1)),i)./data(b(designmat(j,2)),i),'k.')
    end
end

% Now looking at paired threshold ratios.  Lame, but whatever.
L = data(:,7) == 3;
a = find(L);
b = find(~L);
a = a(1:min(length(a),length(b)))
b = b(1:min(length(a),length(b)))

ratios = data(a,[1 2 3])./data(b,[1 2 3]);
figure; subplot(2,1,1); hold on;
plot(ratios(:,1),'b.-')
plot(ratios(:,2),'y.-')
legend({'2°','10°'});
xlabel('experiment'); ylabel('threshold ratio 3/15 Hz');
[h,p] = ttest(ratios(:,1)-ratios(:,2))
subplot(2,1,2);
hist(ratios(:,1)-ratios(:,2));
xlabel('threshold ratio 3/15 Hz');
title(['p = ',num2str(p)]);


% Getting summary statistics on the behavioral performance
out = [];
for i = 1:size(allcolors,2)
    for j = 1:length(speeds)
        spd = speeds(j);
        L = data(:,7) == spd;
        out(j,:,i) = [mean(data(L,i)) std(data(L,i)) sum(L)];
    end
end

% Making a plot of the threshold ratios
% Using the formulas from Kamerund 1978
figure; axes; hold on;
for i = 1:size(out,3)
    muX = out(1,1,i);
    sigmaX = out(1,2,i)./sqrt(out(1,3,i));
    muY = out(2,1,i);
    sigmaY = out(2,2,i)./sqrt(out(2,3,i));
    
    z = linspace(0.2,.6,1000);
    w = (sigmaY/sigmaX)*z;  % w is a scaled version of z
    s = 1./sqrt(w.^2+1);
    k = (muX/sigmaX.*w+muY/sigmaY).*s.^2;
    M = -.5*(muY/sigmaY.*w-muX/sigmaX).^2.*s.^2;
    Q = k.*s.*sqrt(2*pi).*(1-2.*normcdf(-k./s))+(2*s.^2.*exp(-k.^2./(2*s.^2)));
    
    fg = 1/(2*pi).*Q.*exp(M);
    fz = (sigmaY/sigmaX)*fg;
    % plot(z,fz./max(fz),'m-');
    
    % CDF
    Fz = cumsum(fz);
    lastidx = find(diff(Fz) == 0,1,'first');
    if (~isempty(lastidx))
        Fz = Fz(1:lastidx);
        z = z(1:lastidx);
    end
    Fz = Fz./Fz(end);
    CI = interp1(Fz,z,[.05 .95]);
    plot(i,muX/muY,'k*')
    plot([i,i],CI);
end
set(gca,'Xlim',[0 4]);
xlabel('Color direction'); ylabel('Threshold ratio 3/15 Hz');


% Looking at the microsaccadic suppression effect
figure;
for j = 1:length(speeds)
    spd = speeds(j);
    L = data(:,7) == spd;
    subplot(length(speeds),1,j)
    boxplot(data(L,[4 5 6]));
    hold on;
    title(['speed: ',num2str(spd)]);
    for i = 1:3
        [h,p] = ttest(data(L,i+3));
        if (h)
            plot(i,0,'m*');
        end
    end
    ylabel('Microsaccade effect');
end

xlabel('Color direction');
% Comparing microsaccade effects between S-cone isolating
% directions at 3 Hz.
L = data(:,7) == 3;
[h,p] = ttest(data(L,4)-data(L,5))
nansum(x_sac)./nansum(n_sac)
nansum(x_nosac)./nansum(n_nosac)

% Paired comparisons
for i = 1:3
    effect = (sum(x_nosac(:,i))./sum(n_nosac(:,i)))-(sum(x_sac(:,i))./sum(n_sac(:,i)));
    [h,p] = equalproptest([sum(x_sac(:,i)) sum(x_nosac(:,i))],[sum(n_sac(:,i)) sum(n_nosac(:,i))],0.05);
    disp(['Color dir ',num2str(i),' ',num2str(effect),' p = ',num2str(p)]);
    L = any(n_nosac > 0,2);
    [h,p] = ttest(x_sac(L,i)./n_sac(L,i)-x_nosac(L,i)./n_nosac(L,i));
    disp(['One sample t-test on relative suppression: ',num2str(p)]);
    disp(' ');
end

%%
% Section 16
% Looking for a correlation between the time of the last microsaccade
% (relative to fixation point off) and the latency to the operant saccade.

filenames = fnamesFromTxt;
data = [];

for a = 1:size(filenames,1)
    stro = nex2stro(findfile(filenames(a,:)));
    stro = DTfilterquesttrials(stro,'PaperDefaults');
    if isempty(stro.sum.analog.sigid) || isempty(stro.trial)
        continue;
    end
    ntrials = size(stro.trial,1);
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_on'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_off'));
    frameon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'frame_on'));
    targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_on'));
    correct = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
    sacstats = getSacData(stro);
    close;
    
    for i = 1:ntrials
        lastmsac = find(sacstats.endtimes{i} < targon_t(i),1,'last');
        if any(lastmsac)
            if (length(sacstats.starttimes{i}) >= lastmsac+1)
                msaclat = targon_t(i)-sacstats.starttimes{i}(lastmsac);
                saclat = sacstats.starttimes{i}(lastmsac+1)-targon_t(i);
                data = [data; msaclat saclat];
            end
        end
    end
end
figure;
subplot(2,2,1);
plot(data(:,1),data(:,2),'k.')
xlabel('\musac latency');
ylabel('saccade latency');
subplot(2,2,2);
hist(data(:,1),20);
xlabel('\musac latency');
subplot(2,2,3);
hist(data(:,2),20);
xlabel('saccade latency');

[h,p]= corrcoef(data)

%%
% Section 17
% Does the magnitude of saccadic suppression vary with microsaccade
% amplitude?  Looking at achromatic only.  During starting plateau only.
filenames = fnamesFromTxt;
offsets = [.167 .5];   % relative to stim on, when to look for saccades
data = [];
for b = 1:size(filenames,1)
    stro = nex2stro(findfile(filenames(b,:)));
    stro = DTfilterquesttrials(stro,'PaperDefaults');

    if isempty(stro.sum.analog.sigid) || isempty(stro.trial)
        continue;
    end
        
    ntrials = size(stro.trial,1);
    lms = stro.sum.exptParams.RF_colors;
    lms = reshape(lms,[3,size(lms,1)/3]);
    lms(:,all(lms == 0)) = [];
    lms = mkbasis(lms);
    whichcol = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'color_dir'));
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_on'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_off'));
    targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_on'));
    correct = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
    flashR = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_R'));
    flashG = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_G'));
    flashB = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_B'));
    bkgndrgb = [stro.sum.exptParams.bkgnd_r, stro.sum.exptParams.bkgnd_g, stro.sum.exptParams.bkgnd_b];
    M = reshape(stro.sum.exptParams.m_mtx, 3, 3);
    bkgndlms = M * bkgndrgb';
    x = 0:255; %the normal range of the gamma look up table
    xx = linspace(0, 255, 2^16); %the desired quantization of the gammaTable
    g1 = reshape(stro.sum.exptParams.gamma_table, 256, 3);
    gammaTable = [spline(x, g1(:,1), xx)', spline(x, g1(:,2), xx)', spline(x, g1(:,3), xx)'];
    
    sacstats = getSacData(stro);
    close;
    
    % Looping over trials
    % And getting data for precent correct analysis
    for j = 1:ntrials
        coloridx = find(abs(lms(:,whichcol(j))'*[.5774 .5774 .5774]') > .99);  % Avoiding roundoff errors
        if (coloridx ~= 1)
            continue;
        end
        st = sacstats.starttimes{j};
        L = logical(sacstats.endtimes{j} >= targon_t(j));
        st(L) = [];    % Omitting saccades that terminate after targets appear
        Lsac = (st > stimon_t(j)+offsets(1)) & (st < stimon_t(j)+offsets(2));
        if (any(Lsac))
            k = find(Lsac);
            data = [data; correct(j) max(sacstats.amplitudes{j}(k))];
        else
            data = [data; correct(j) nan];            
        end
    end
end
binedges = prctile(data(:,2),[0:10:100]);
tmpdata = [];
for j = 2:length(binedges)
    L = data(:,2)>binedges(j-1) & data(:,2)<binedges(j);
    binmiddle = (binedges(j)+binedges(j-1))/2;
    tmpdata = [tmpdata; binmiddle sum(data(L,1)) sum(L)];
end
figure; axes; hold on;
pcor = tmpdata(:,2)./tmpdata(:,3);
plot(tmpdata(:,1),pcor,'.-')
plot(tmpdata(:,1),pcor + sqrt((pcor.*(1-pcor))./tmpdata(:,2)) ,'-');
plot(tmpdata(:,1),pcor - sqrt((pcor.*(1-pcor))./tmpdata(:,2)) ,'-');
Lnosac = isnan(data(:,2));
pcornosac = sum(data(Lnosac,1))/sum(Lnosac);
plot([0 1],[pcornosac pcornosac],'k-');
[h,p] = equalproptest(tmpdata(:,2),tmpdata(:,3),0.05)

% Testing large (>0.5 deg) saccades vs the others
L_bigsac = data(:,2) >= 0.5;
sum(L_bigsac)/length(L_bigsac)  % What fraction of trials have big microsaccades?
sum(data(L_bigsac,1))./sum(L_bigsac)
[h,p] = equalproptest([sum(data(L_bigsac,1)) sum(data(~L_bigsac,1))],[sum(L_bigsac) sum(~L_bigsac)],0.05)


%%
% Section 17.1
% Does the magnitude of saccadic suppression vary with microsaccade
% amplitude?  Achromatic only.  As a function of time in the trial.
filenames = fnamesFromTxt;
binwidth = .05; %sec
bins = [-.2:binwidth:.767];
x_sac = zeros(2,length(bins));
n_sac = zeros(2,length(bins));
x_nosac = 0;
n_nosac = 0;
for b = 1:size(filenames,1)
    stro = nex2stro(findfile(filenames(b,:)));
    stro = DTfilterquesttrials(stro,'PaperDefaults');

    if isempty(stro.sum.analog.sigid) || isempty(stro.trial)
        continue;
    end
        
    ntrials = size(stro.trial,1);
    lms = stro.sum.exptParams.RF_colors;
    lms = reshape(lms,[3,size(lms,1)/3]);
    lms(:,all(lms == 0)) = [];
    lms = mkbasis(lms);
    whichcol = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'color_dir'));
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_on'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_off'));
    targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_on'));
    correct = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
    flashR = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_R'));
    flashG = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_G'));
    flashB = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_B'));
    bkgndrgb = [stro.sum.exptParams.bkgnd_r, stro.sum.exptParams.bkgnd_g, stro.sum.exptParams.bkgnd_b];
    M = reshape(stro.sum.exptParams.m_mtx, 3, 3);
    bkgndlms = M * bkgndrgb';
    x = 0:255; %the normal range of the gamma look up table
    xx = linspace(0, 255, 2^16); %the desired quantization of the gammaTable
    g1 = reshape(stro.sum.exptParams.gamma_table, 256, 3);
    gammaTable = [spline(x, g1(:,1), xx)', spline(x, g1(:,2), xx)', spline(x, g1(:,3), xx)'];
    
    sacstats = getSacData(stro);
    close;
    
    amplitudethreshold = 0.5;  % Reviewer 1's suggested threshold
  %  amplitudethreshold = 0.3;  % near the median

    
    % Looping over trials
    % And getting data for precent correct analysis
    for j = 1:ntrials
        coloridx = find(abs(lms(:,whichcol(j))'*[.5774 .5774 .5774]') > .99);  % Avoiding roundoff errors
        if (coloridx ~= 1)
            continue;
        end
        st = sacstats.starttimes{j};
        L = logical(sacstats.endtimes{j} >= targon_t(j));
        st(L) = [];    % Omitting saccades that terminate after targets appear
        Lsac = (st > stimon_t(j)+bins(1)-binwidth/2) & (st < stimon_t(j)+bins(end)+binwidth/2);
        if (any(Lsac))
            for k = find(Lsac)'
                whichbin = logical(hist(sacstats.starttimes{j}(k)-stimon_t(j),bins));
                amp = sacstats.amplitudes{j}(k);
                if (amp > amplitudethreshold)
                    ampidx = 1;
                else
                    ampidx = 2;
                end
                
                if (correct(j))
                    x_sac(ampidx,whichbin) = x_sac(ampidx,whichbin)+1;
                end
                n_sac(ampidx,whichbin) = n_sac(ampidx,whichbin)+1;
            end
        end
        if (~any(Lsac))
            if (correct(j))
                x_nosac = x_nosac+1;
            end
            n_nosac = n_nosac+1;
        end
    end
end
figure
plot(bins,x_sac'./n_sac')
%%
% Section 18
% Comparing microsaccade parameters as a function of time during data
% collection.

filenames = fnamesFromTxt;
amplitudes = [];
directions = [];
peakv = [];
starttimes = [];
pertrialdata = [];
totaltime = [0 0 0];
nsacs = 2000;
tmptime = 0;

for a = 1:size(filenames,1)
    stro = nex2stro(findfile(filenames(a,:)));
    stro = DTfilterquesttrials(stro,'PaperDefaults');
    if isempty(stro.sum.analog.sigid) || isempty(stro.trial)
        continue;
    end
    ntrials = size(stro.trial,1);
    lms = stro.sum.exptParams.RF_colors;
    lms = reshape(lms,[3,size(lms,1)/3]);
    lms(:,all(lms == 0)) = [];
    lms = mkbasis(lms);
    whichcol = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'color_dir'));
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_on'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_off'));
    frameon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'frame_on'));
    targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_on'));
    correct = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
    sacstats = getSacData(stro);
    close;
    
    for i = 1:ntrials
        st = sacstats.starttimes{i};
        Lsac = (st > frameon_t(i)) & (st < targon_t(i)-.2);  % ad hoc margin
        tmptime = tmptime + (targon_t(i)-.2-frameon_t(i));        
        totaltime(2) = totaltime(2)+ (targon_t(i)-.2-frameon_t(i));
        if (length(amplitudes) < nsacs)
            totaltime(1) = totaltime(1)+ (targon_t(i)-.2-frameon_t(i));
            tmptime = 0;
        end
        if any(Lsac)
            if any(sacstats.amplitudes{i}(Lsac) > 2)
                keyboard
            end
            amplitudes = [amplitudes; sacstats.amplitudes{i}(Lsac)];
            directions = [directions; sacstats.directions{i}(Lsac)];
            peakv = [peakv; sacstats.peakv{i}(Lsac)];
            starttimes = [starttimes; st(Lsac)-stimon_t(i)];
            if (length(amplitudes) > nsacs)
                totaltime(3) = totaltime(3)+tmptime;
                tmptime = 0;
            end
        end
        coloridx = find(abs(lms(:,whichcol(i))'*[.5774 .5774 .5774]') > .99);  % Avoiding roundoff errors
        if (coloridx)
            Lplateau = (st > stimon_t(i) + .167) & (st < stimon_t(i)+.5);
            pertrialdata = [pertrialdata; correct(i) any(Lplateau)];
        end
    end
end

figure;
early = 1:nsacs;
late = length(amplitudes)-nsacs+1:length(amplitudes);
[n1,x1] = hist(amplitudes(early),[0:.05:1])
[n2,x2] = hist(amplitudes(late),[0:.05:1])
subplot(2,2,1); hold on;
plot(x1,n1,'k-','LineWidth',2);
plot(x2,n2,'b-','LineWidth',2);
set(gca,'XTick',[0 .25 .5 .75 1]);
mean(amplitudes(early))
mean(amplitudes(late))
[h,p] = ttest2(amplitudes(early),amplitudes(late))
xlabel('Amplitude (deg)');
ylabel('# saccades');

subplot(2,2,2); 
[r1, t1]= hist(directions(early),[-pi:pi/8:pi]);
[r2, t2]= hist(directions(late),[-pi:pi/8:pi]);
h(1) = polar([t1 t1(1)],[r1 r1(1)],'k-');hold on;
h(2) = polar([t2 t2(1)],[r2 r2(1)],'b-');
set(h,'Linewidth',2);
[h,p] = WatsonU2Test (directions(early),directions(late))

subplot(2,2,3); hold on;
binwidth = 0.05;
[n1,x1] = hist(starttimes(early),[-.2:binwidth:1])
[n2,x2] = hist(starttimes(late),[-.2:binwidth:1])
plot(x1,n1/binwidth/totaltime(1),'k-','LineWidth',2);
plot(x2,n2/binwidth/(totaltime(2)-totaltime(3)),'b-','LineWidth',2);
xlabel('Time from stim on (s)');
ylabel('Saccades per sec');
set(gca,'XLim',[-.2 1]);

subplot(2,2,4); hold on;
ntrials = 2000;
tmpearly = pertrialdata(1:ntrials,[1 2]);
Lsac = logical(tmpearly(:,2));
p(1) = sum(tmpearly(Lsac))./sum(Lsac);
se(1) = p(1)*(1-p(1))/sqrt(sum(Lsac));
p(2) = sum(tmpearly(~Lsac))./sum(~Lsac);
se(2) = p(2)*(1-p(2))/sqrt(sum(~Lsac));
tmplate = pertrialdata(length(pertrialdata)-ntrials+1:length(pertrialdata),[1 2]);
Lsac = logical(tmplate(:,2));
p(3) = sum(tmplate(Lsac))./sum(Lsac);
se(3) = p(3)*(1-p(3))/sqrt(sum(Lsac));
p(4) = sum(tmplate(~Lsac))./sum(~Lsac);
se(4) = p(4)*(1-p(4))/sqrt(sum(~Lsac));
subplot(2,2,4); hold on;
bar([1 2 3 4], p);
for i =1:4
    plot([i i],p(i)+[se(i) -se(i)],'k-','LineWidth',2);
end
set(gca,'Ylim',[.6 .8],'XLim',[.1 4.9])
  