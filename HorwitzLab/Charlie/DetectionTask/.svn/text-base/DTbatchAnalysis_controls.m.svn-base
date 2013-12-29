%% (22) POPULATION WIDE SAC TRIGGERED PSTH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
global cardVsIntBatchPath blpBatchPath
preprocessDTbatchData

rawWeigths = nan(length(out.dat),3);
for a = 1:length(out.dat)
    tmp = out.dat(a).grating.color.prefcolor;
    tmp = tmp ./ sum(abs(tmp)); %normalize the weights
    %make sure that the points all fall in the upper part of the diamond
    if sum(sign(tmp(1:2))) == 0 %i.e., L-M or -L+M
        if tmp(1) > 0;
            tmp = tmp.*-1;
        end
    elseif abs(sum(sign(tmp(1:2)))) == 2 %i.e., L+M or -L-M
        if tmp(1) < 0;
            tmp = tmp.*-1;
        end
    end
    rawWeights(a,:) = tmp;
end

%iterate over the DT files to compute the sac-triggered psth
l_validExpts = ~(commonExclusions | out.errors(:,neuroThresh1Ind) | out.errors(:,neuroThresh2Ind));
fileList = find(l_validExpts);
filtWeights = rawWeights(l_validExpts,:);
NSHUFFS = 100;
preTime = .250;
postTime = .500;
BINWIDTH = .020;
timebins = -preTime:BINWIDTH:postTime;
data = nan(length(fileList), length(timebins));
shuffData = nan(length(fileList), length(timebins));
shuffNormToMaxData = nan(length(fileList), length(timebins));
lm = nan(length(fileList), 2);
fileName = {};
for a = 1:length(fileList)
    expt = fileList(a);
    out.fnames{expt}
    
    %open a file
    DT = dtobj(out.fnames{expt}{1});
    DT = stripOutGratingTrials(DT);
    l_inRF = DT.trial(:, DT.idx.flashX) == DT.sum.exptParams.rf_x;
    frameon_t = DT.trial(:,DT.idx.frameOn);
    targon_t = DT.trial(:, DT.idx.targOn);
    choice_t = DT.trial(:,DT.idx.choiceTime);
    
    %get the analog data. if it's not there than skip this file.
    if ~isempty(DT.sum.analog.sigid)
        sacstats = getSacData(DT);
        close(gcf);
    else
        fprintf('   ****** File <%d> has no analog data *****', a)
        continue
    end
    
    %setup the shuffled data correspondences
    shuffInd = nan(NSHUFFS+1,size(DT.trial,1));
    shuffInd(1,:) = [1:size(DT.trial,1)];
    for i = 2:NSHUFFS+1;
        shuffInd(i,:) = randperm(size(DT.trial,1));
    end
    
    
    trl_psth = nan(size(DT.trial,1), length(timebins));
    expt_psth = nan(NSHUFFS+1, length(timebins));
    for shuff = 1:NSHUFFS+1;
        for j = 1:size(DT.trial,1);
            anlyStart = frameon_t(j);
            anlyEnd = frameon_t(j)+0.200+0.666;
            
            sacStarts = sacstats.starttimes{shuffInd(shuff,j)};
            sacStarts = sacStarts - DT.ras{shuffInd(shuff,j), DT.idx.anlgStart}; %sacStarts from time zero
            sacStarts = sacStarts + DT.ras{j, DT.idx.anlgStart}; %sacStarts from trial(trl) anlgStart
            sacStarts(sacStarts < anlyStart) = [];
            sacStarts(sacStarts > anlyEnd) = [];
            if any(sacStarts)
                sac_psth = [];
                for k = 1:length(sacStarts)
                    spiketimes = DT.ras{j,1}-sacStarts(k);
                    spiketimes(spiketimes < -preTime-BINWIDTH/2) = [];
                    spiketimes(spiketimes > postTime+BINWIDTH/2) = [];
                    sac_psth(end+1,:) = hist(spiketimes, timebins);
                    sac_psth(end,timebins+sacStarts(k) > choice_t(j)) = nan;
                end
                trl_psth(j,:) = nanmean(sac_psth,1);
            end
        end
        expt_psth(shuff,:) = nanmean(trl_psth,1)./BINWIDTH;
    end
    
    %complie the psth's, coneweights, and filenames.
    data(a,:) = expt_psth(1,:);
    tmp = expt_psth(2:end,:);
    shuffData(a,:) = nanmean(tmp,1);
    shuffNormToMaxData(a,:) = nanmean(tmp./repmat(max(tmp,[],2),1,size(tmp,2)),1);
    lm(a,:) = sign(filtWeights(a, [1,2]));
    fileName{a} = out.fnames{a}{1};
end

NORMMETH = 'baseline';
switch lower(NORMMETH)
    case 'max'
        normPSTH = data./repmat(max(data,[],2), 1, size(data,2));
        normShuffPSTH = shuffNormToMaxData;
        titleString = 'Normalized to Max';
    case 'baseline'
        baselineBins = timebins <= 0;
        normPSTH = data./repmat(nanmean(data(:,baselineBins),2), 1, size(data,2));
        normShuffPSTH = shuffData./repmat(nanmean(shuffData(:,baselineBins),2), 1, size(shuffData,2));
        titleString = 'Normalized to Baseline';
end


%start by ploting the data in non-normalized units
figure
hold on,
plot(timebins, nanmean(data,1), 'k')
plot(timebins, nanmean(shuffData,1), 'r:')
title('Non-Normalized Average PSTH')
xlabel('Time From Microsaccade')
ylabel('Sp/Sec')
hold off


%now plot normalized
figure, hold on,
normSEM = nanstd(normPSTH,[],1)./sqrt(size(normPSTH,1));
shuffNormSEM = nanstd(normShuffPSTH,[],1)./sqrt(size(normShuffPSTH,1));
patch([timebins, fliplr(timebins)], [nanmean(normShuffPSTH,1)+shuffNormSEM, fliplr(nanmean(normShuffPSTH,1)-shuffNormSEM)], 'r', 'facecolor', [1 0.5 0.5])
plot(timebins, nanmean(normShuffPSTH,1), 'r')
patch([timebins, fliplr(timebins)], [nanmean(normPSTH,1)+normSEM, fliplr(nanmean(normPSTH,1)-normSEM)], 'k', 'facecolor', [0.6 0.6 0.6])
plot(timebins, nanmean(normPSTH,1), 'k')
xlabel('Time From Microsaccade')
ylabel('Percent Change')
title(titleString)
hold off

%seperate the psth's by cone weights
l_RG = sum(lm, 2) == 0;
l_LUM = sum(lm, 2) == 2;
figure
subplot(1,2,1), hold on, %L-M
plot(timebins, nanmean(normPSTH(l_RG,:),1), 'k')
plot(timebins, nanmean(normShuffPSTH(l_RG,:),1), 'r:')
xlabel('Time From Microsaccade')
ylabel(titleString)
title('L-M')
hold off
subplot(1,2,2), hold on, %L+M
plot(timebins, nanmean(normPSTH(l_LUM,:),1), 'k')
plot(timebins, nanmean(normShuffPSTH(l_LUM,:),1), 'r:')
xlabel('Time From Microsaccade')
ylabel(titleString)
title('L+M')
hold off

%seperate by monkey
for a = 1:length(fileName)
    if ~isempty(fileName{a})
        initials(a,1) = fileName{a}(1);
    else
        initials(a,1) = ' ';
    end
end
l_sedna = initials=='S';
l_kali = initials=='K';
figure
subplot(1,2,1), hold on, %Sedna
sem = nanstd(normPSTH(l_sedna,:),[],1) ./ sqrt(sum(l_sedna));
avg = mean(normPSTH(l_sedna,:),1);
patch([timebins, fliplr(timebins)], [avg+sem, fliplr(avg-sem)], 'k', 'facecolor', [0.6, 0.6, 0.6])
plot(timebins, avg, 'k')
semShuff = nanstd(normShuffPSTH(l_sedna,:),[],1) ./ sqrt(sum(l_sedna));
avgShuff = mean(normShuffPSTH(l_sedna,:),1);
patch([timebins, fliplr(timebins)], [avgShuff+semShuff, fliplr(avgShuff-semShuff)], 'r', 'facecolor', [1, 0.5, 0.5])
plot(timebins, avgShuff, 'r')
xlabel('Time From Microsaccade')
ylabel(titleString)
title('Sedna')
xlim([min(timebins), max(timebins)])
hold off
subplot(1,2,2), hold on, %kali
sem = nanstd(normPSTH(l_kali,:),[],1) ./ sqrt(sum(l_kali));
avg = mean(normPSTH(l_kali,:),1);
patch([timebins, fliplr(timebins)], [avg+sem, fliplr(avg-sem)], 'k', 'facecolor', [0.6, 0.6, 0.6])
plot(timebins, avg, 'k')
semShuff = nanstd(normShuffPSTH(l_kali,:),[],1) ./ sqrt(sum(l_kali));
avgShuff = mean(normShuffPSTH(l_kali,:),1);
patch([timebins, fliplr(timebins)], [avgShuff+semShuff, fliplr(avgShuff-semShuff)], 'r', 'facecolor', [1, 0.5, 0.5])
plot(timebins, avgShuff, 'r')
xlabel('Time From Microsaccade')
ylabel(titleString)
title('Kali')
xlim([min(timebins), max(timebins)])
hold off

%% (22.2) FIRING RATE BETWEEN SAC AND NOSAC TRIALS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
global cardVsIntBatchPath blpBatchPath
preprocessDTbatchData

sacZscores = [];
nosacZscores = [];
for a = 1:length(out.dat);
    out.fnames{a}
    if commonExclusions(a)
        disp('Does not meet inclusion criteria')
        continue
    end
    
    %open a file
    DT = dtobj(out.fnames{a}{1});
    DT = stripOutGratingTrials(DT);
    l_inRF = DT.trial(:, DT.idx.flashX) == DT.sum.exptParams.rf_x;
    contrast = DT.trial(:, DT.idx.cntrstLev);
    colorDir = DT.trial(:, DT.idx.colorDir);
    gaborOn_t = DT.trial(:, DT.idx.flashOn);
    gaborOff_t = DT.trial(:,DT.idx.flashOff);
    spikes = DT.ras(:, DT.idx.spikes);
    tStart = mat2cell(gaborOn_t, ones(length(gaborOn_t),1));
    tEnd = mat2cell(gaborOff_t, ones(length(gaborOff_t),1));
    nSpikes = cellfun(@(x,y,z)(sum((x>y)&(x<=z))), spikes, tStart, tEnd);
    
    %get the analog data. if it's not there than skip this file.
    if ~isempty(DT.sum.analog.sigid)
        sacstats = getSacData(DT);
        nSacDuringGabor = cellfun(@(x,y,z)(sum((x>y)&(x<=z))), sacstats.starttimes', tStart, tEnd);
        sacTrials = nSacDuringGabor>0;
        close(gcf);
    else
        fprintf('   ****** File <%d> has no analog data *****', a)
        continue
    end
    
    for clr = 1:max(colorDir);
        %only look at the data if the neuron fired spikes to a near
        %threshold stimulus
        if out.dat(a).c.alpha(clr,1)
            disp('nan neurothresh')
            continue
        end
        
        %don't look at non-zero contrasts or out or RF trials.
        for cntrst = 2:max(contrast);
            if (clr==1)&&(cntrst==1) %i.e., just do this once
                sacTList = (~l_inRF | colorDir==1) & sacTrials;
                noSacTList = (~l_inRF | colorDir==1) & ~sacTrials;
            else
                sacTList = l_inRF & (colorDir==clr) & (contrast==cntrst) & sacTrials;
                noSacTList = l_inRF & (colorDir==clr) & (contrast==cntrst) & ~sacTrials;
            end
            sacCounts = nSpikes(sacTList);
            noSacCounts = nSpikes(noSacTList);
            sacCounts = nSpikes(sacTList);
            noSacCounts = nSpikes(noSacTList);
            
            %create a pooled estimate of sigma and mu.
            n1 = numel(sacCounts);
            n2 = numel(noSacCounts);
            if~n1
                if ~n1 && ~n2; keyboard; end
                sigma = std(noSacCounts);
                mu = mean(noSacCounts);
            elseif ~n2
                sigma = std(sacCounts);
                mu = mean(sacCounts);
            else
                sigma = sqrt(((n1-1)*var(sacCounts) + (n2-1)*var(noSacCounts))/(n1+n2-2));
                mu = (n1.*mean(sacCounts)+n2.*mean(noSacCounts))./ (n1+n2);
            end
            tmpSacZ = (sacCounts-mu) ./ (sigma+eps);
            tmpNoSacZ = (noSacCounts-mu) ./ (sigma+eps);
            sacZscores = [sacZscores(:);tmpSacZ(:)];
            nosacZscores = [nosacZscores(:);tmpNoSacZ(:)];
        end
    end
end

%a simple t test on the z scores
[h,p] = ttest2(sacZscores, nosacZscores)

%a permutation test
NPERMS = 10000;
permScores = nan(1,NPERMS);
combDat = [sacZscores(:); nosacZscores(:)];
for a = 1:NPERMS
    permInd = randperm(length(combDat));
    tmpSac = combDat(permInd(1:length(sacZscores)));
    tmpNoSac = combDat(permInd(length(sacZscores)+1:end));
    sigma = sqrt(((length(tmpSac)-1)*var(tmpSac) + (length(tmpNoSac)-1)*var(tmpNoSac))/(length(tmpSac)+length(tmpNoSac)-2));
    permScores(a) = (mean(tmpNoSac)-mean(tmpSac))./sigma;
end
exptPoolSigma = sqrt(((length(sacZscores)-1)*var(sacZscores) + (length(nosacZscores)-1)*var(nosacZscores))/(length(sacZscores)+length(nosacZscores)-2));
exptScore = (mean(nosacZscores)-mean(sacZscores))./exptPoolSigma
critVals = prctile(permScores, [0.025 0.975])
p = sum((permScores<-abs(exptScore)) | (permScores>abs(exptScore)))./length(permScores)

%% (25) 2º vs 10º fundamentals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, clc
global cardVsIntBatchPath blpBatchPath
preprocessDTbatchData

l_tenDegFund = out.errors(:, tenDegFundInd);
l_valid = l_tenDegFund & ~(commonExclusions | out.errors(:, cardTieInd) | out.errors(:, intTieInd));


figure %all the TRs
minTRs = min(rawTRs,[],2);
hist(minTRs(l_valid))
xlabel('TR for 10deg Fund')
ylabel('Counts')


% card vs int
card = rawCardTRs(l_valid);
int = rawIntTRs(l_valid);
valfornan = max([card(:);int(:)]) .* 1.2;
card(isnan(card)) = valfornan;
int(isnan(int)) = valfornan;
sedna = l_sedna(l_valid);
kali = l_kali(l_valid);
figure, hold on,
plot(card(sedna), int(sedna), 'kv', 'markerfacecolor', 'g', 'markersize', 8)
plot(card(kali), int(kali), 'ko', 'markerfacecolor', 'b', 'markersize', 8)
plot([.2, 1.5*valfornan], [.2, 1.5*valfornan], 'k')
set(gca, 'yscale', 'log', 'xscale', 'log')
legend('Sedna', 'Kali')
xlim([0.5 1.5*valfornan])
ylim([0.5 1.5*valfornan])
xlabel('Card TR')
ylabel('Int TR')

%% (26) FRAMES VS. NO-FRAMES (MOCS Spikes)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, clc
global cardVsIntBatchPath blpBatchPath
preprocessDTbatchData

PRETIME = -0.200;
POSTTIME = 0.200;
BINSIZE = 0.010;
bins = -0.200:BINSIZE:0.200;
psth = nan(length(out.dat), length(bins));
nosac_psth = nan(length(out.dat), length(bins));
sacsPerSec = nan(length(out.dat), length(bins));
framePresent = ones(length(out.dat),1); %initialize to ones, than change to zeros as needed
for a = 1:length(out.dat);
    DT = dtobj(out.fnames{a}{1});
    
    % flag expts that didn't use frames
    if all(DT.trial(:, DT.idx.framesPresent) == 0)
        framePresent(a) = 0;
    end
    
    %get the counts pre and post frame presentation (but only when the
    %frames were present)
    frameOnTime = DT.trial(:, DT.idx.frameOn);
    frameOnTime = mat2cell(frameOnTime, ones(1, length(frameOnTime)));
    spikeTimesFromFrames = cellfun(@(x,y) x-y, DT.ras(:, DT.idx.spikes), frameOnTime, 'uniformoutput', 0);
    zeroAsCells = mat2cell(zeros(length(frameOnTime),1), ones(length(frameOnTime),1));
    preTimeAsCells = mat2cell(repmat(PRETIME, length(frameOnTime),1), ones(length(frameOnTime),1));
    postTimeAsCells = mat2cell(repmat(POSTTIME, length(frameOnTime),1), ones(length(frameOnTime),1));
    preFrameCounts = cellfun(@(x,y,z) sum((x>y) & (x<=z)), spikeTimesFromFrames, preTimeAsCells, zeroAsCells);
    postFrameCounts = cellfun(@(x,y,z) sum((x>y) & (x<=z)), spikeTimesFromFrames, zeroAsCells, postTimeAsCells);
    
    meanOfDiffs(a) = mean(postFrameCounts-preFrameCounts);
    [h(a), p(a)] = ttest(postFrameCounts, preFrameCounts);
    
    %compile the psth synced to frame on (or the time when the frames would
    %have come on had there been frames.
    edges = mat2cell(repmat(bins, length(spikeTimesFromFrames),1), ones(length(spikeTimesFromFrames),1));
    trl_psth = cellfun(@histc, spikeTimesFromFrames, edges, 'uniformoutput', 0);
    trl_psth = cellfun(@(x) x(:)', trl_psth, 'uniformoutput', 0);
    trl_psth = vertcat(trl_psth{:});
    psth(a,:) = mean(trl_psth, 1) ./ BINSIZE;
    
    %now looking at eyemovements
    if ~isempty(DT.sum.analog.sigid)
        sacstats = getSacData(DT, 0.3);
        close(gcf);
        
        %calculating the frame triggered saccade average
        sacTimesFromFrames = cellfun(@(x,y) x-y, sacstats.starttimes', frameOnTime, 'uniformoutput', 0);
        trl_sacsPerSec = cellfun(@histc, sacTimesFromFrames, edges, 'uniformoutput', 0);
        trl_sacsPerSec = cellfun(@(x) x(:)', trl_sacsPerSec, 'uniformoutput', 0);
        trl_sacsPerSec = vertcat(trl_sacsPerSec{:});
        sacsPerSec(a,:) = mean(trl_sacsPerSec,1) ./ BINSIZE;
        
        %calculating a spike psth for trials w/o saccades during the
        %analysis window
        sacTrials = sum(trl_sacsPerSec,2) > 0;
        nosac_psth(a,:) = mean(trl_psth(~sacTrials,:),1) ./ BINSIZE;
        
        % ttest on counts (pre vs. post frames) for no-sac trials
        [h_nosac(a), p_nosac(a)] = ttest(postFrameCounts(~sacTrials), preFrameCounts(~sacTrials));
    else
        fprintf('   ****** File <%d> has no analog data *****', a)
        continue
    end
    
    
end

%normalize the spike psth's by their background (pre frame) rate
NORMMETH = 'max';
switch NORMMETH
    case 'baseline'
        prebins = bins<0;
        norm_psth = psth ./ repmat(mean(psth(:,prebins),2),1,size(psth,2));
        norm_nosac_psth = nosac_psth ./ repmat(mean(nosac_psth(:,prebins),2),1,size(nosac_psth,2));
    case 'max'
        norm_psth = psth ./ repmat(max(psth,[],2),1,size(psth,2));
        norm_nosac_psth = nosac_psth ./ repmat(max(nosac_psth,[],2),1,size(nosac_psth,2));
end

%plot psths for instances of sig increase w/out frames.
l_sigIncreseNoFrames = (h_nosac'>0) & (meanOfDiffs' > 0) & ~framePresent;
figure, hold on,
stairs(bins, nanmean(norm_psth((l_sigIncreseNoFrames&l_kali),:)), 'k')
stairs(bins, nanmean(norm_nosac_psth((l_sigIncreseNoFrames&l_kali),:)), 'k--')
stairs(bins, nanmean(norm_psth((l_sigIncreseNoFrames&l_sedna),:)), 'b')
stairs(bins, nanmean(norm_nosac_psth((l_sigIncreseNoFrames&l_sedna),:)), 'b--')
title('Sig Increase, No Frames')
xlabel('time from frames on')
ylabel('rate')
legend('Sacs Kali', 'No Sacs Kali', 'Sacs Sedna', 'No Sacs Sedna', 'location', 'northwest')
xlim([bins(1), bins(end)])

%plot psths for instances of sig decrease w/out frames.
l_sigDecreaseNoFrames = (h_nosac'>0) & (meanOfDiffs' < 0) & ~framePresent;
figure, hold on,
stairs(bins, nanmean(norm_psth((l_sigDecreaseNoFrames&l_kali),:)), 'k')
stairs(bins, nanmean(norm_nosac_psth((l_sigDecreaseNoFrames&l_kali),:)), 'k--')
stairs(bins, nanmean(norm_psth((l_sigDecreaseNoFrames&l_sedna),:)), 'b')
stairs(bins, nanmean(norm_nosac_psth((l_sigDecreaseNoFrames&l_sedna),:)), 'b--')
title('Sig Decrease, No Frames')
xlabel('time from frames on')
ylabel('rate')
legend('Sacs Kali', 'No Sacs Kali', 'Sacs Sedna', 'No Sacs Sedna', 'location', 'northeast')
xlim([bins(1), bins(end)])

%plot sedna's data side by side w/w/out frames.
figure,
subplot(2,2,1),hold on, %no frames
stairs(bins, nanmean(norm_psth(~framePresent&l_sedna,:),1), 'k')
stairs(bins, nanmean(norm_nosac_psth(~framePresent&l_sedna,:),1), 'k--')
title('Sedna: no frames')
ylabel('Rate')
subplot(2,2,2)%no frames
hist(meanOfDiffs(~framePresent&l_sedna),20)
[hPop, pPop] = ttest(meanOfDiffs(~framePresent&l_sedna));
title(sprintf('(post-pre) counts\n h=%d, p=%.3f', hPop, pPop))
xlabel('diff')
subplot(2,2,3),hold on, %with frames
stairs(bins, nanmean(norm_psth(framePresent&l_sedna,:),1), 'k')
stairs(bins, nanmean(norm_nosac_psth(framePresent&l_sedna,:),1), 'k--')
title('Sedna: with frames')
ylabel('Rate')
subplot(2,2,4) %with frames
hist(meanOfDiffs(framePresent&l_sedna),20)
[hPop, pPop] = ttest(meanOfDiffs(framePresent&l_sedna));
title(sprintf('(post-pre) counts\n h=%d, p=%.3f', hPop, pPop))
Xlabel('diff')

%plot kali's data side by side w/w/out frames.
figure,
subplot(2,2,1),hold on, %no frames
stairs(bins, nanmean(norm_psth(~framePresent&l_kali,:),1), 'k')
stairs(bins, nanmean(norm_nosac_psth(~framePresent&l_kali,:),1), 'k--')
title('Kali: no frames')
ylabel('Rate')
subplot(2,2,2)%no frames
hist(meanOfDiffs(~framePresent&l_kali),20)
[hPop, pPop] = ttest(meanOfDiffs(~framePresent&l_kali));
title(sprintf('(post-pre) counts\n h=%d, p=%.3f', hPop, pPop))
xlabel('diff')
subplot(2,2,3),hold on, %with frames
stairs(bins, nanmean(norm_psth(framePresent&l_kali,:),1), 'k')
stairs(bins, nanmean(norm_nosac_psth(framePresent&l_kali,:),1), 'k--')
title('Kali: with frames')
ylabel('Rate')
subplot(2,2,4) %with frames
hist(meanOfDiffs(framePresent&l_kali),20)
[hPop, pPop] = ttest(meanOfDiffs(framePresent&l_kali));
title(sprintf('(post-pre) counts\n h=%d, p=%.3f', hPop, pPop))
Xlabel('diff')


%% (27) FRAMES VS. NO-FRAMES (MoCS Behavior)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, clc
global cardVsIntBatchPath blpBatchPath
preprocessDTbatchData

% use sedna's MoCS data to construct sCSFs w/w/out frames
[color, sptFreq, thresh, frames, sedna] = deal([]);
for a = 1:length(out.dat)
    %monkey?
    if lower(out.fnames{a}{1}(1)) == 's'
        tmp_sedna = 1;
    else
        continue;
    end
    
    DT = dtobj(out.fnames{a}{1});
    [m, c, expt] = DTunpack(DT,1);
    close all;
    
    %frames?
    if all(DT.trial(:, DT.idx.framesPresent) == 0)
        tmp_frames = 0;
    else
        tmp_frames = 1;
    end
    
    %color, sf, thresh
    for sf =1:size(m.alpha,2);
        for clr = 1:size(m.alpha,1);
            color = [color; expt.standColors(clr,:)];
            sptFreq = [sptFreq; expt.sfs];
            thresh = [thresh; m.alpha(clr, sf)];
            frames = [frames; tmp_frames];
            sedna = [sedna; tmp_sedna];
        end
    end
    
end

%plot the sCSFs by color (only for sedna)
figure,
subplot(2,2,1), hold on, %lvm
l_lvm = ismember(sign(color), [1 -1 0], 'rows');
l_sedna_LvM_Frames = l_lvm & sedna & frames;
l_sedna_LvM_NoFrames = l_lvm & sedna & ~frames;
plot(sptFreq(l_sedna_LvM_Frames), 1./thresh(l_sedna_LvM_Frames), 'rs')
plot(sptFreq(l_sedna_LvM_NoFrames), 1./thresh(l_sedna_LvM_NoFrames), 'rs', 'markerfacecolor', 'r')
legend('Frames', 'No Frames', 'location', 'southwest')
set(gca, 'xscale', 'log', 'yscale', 'log')
subplot(2,2,2), hold on, %s-iso
l_s = ismember(sign(color), [0 0 1], 'rows');
l_sedna_S_Frames = l_s & sedna & frames;
l_sedna_S_NoFrames = l_s & sedna & ~frames;
plot(sptFreq(l_sedna_S_Frames), 1./thresh(l_sedna_S_Frames), 'bs')
plot(sptFreq(l_sedna_S_NoFrames), 1./thresh(l_sedna_S_NoFrames), 'bs', 'markerfacecolor', 'b')
legend('Frames', 'No Frames', 'location', 'southwest')
set(gca, 'xscale', 'log', 'yscale', 'log')
subplot(2,2,3), hold on, %swm
l_swm = ismember(sign(color), [1 -1 -1], 'rows');
l_sedna_SwM_Frames = l_swm & sedna & frames;
l_sedna_SwM_NoFrames = l_lvm & sedna & ~frames;
plot(sptFreq(l_sedna_SwM_Frames), 1./thresh(l_sedna_SwM_Frames), 'gs')
plot(sptFreq(l_sedna_SwM_NoFrames), 1./thresh(l_sedna_SwM_NoFrames), 'gs', 'markerfacecolor', 'g')
legend('Frames', 'No Frames', 'location', 'southwest')
set(gca, 'xscale', 'log', 'yscale', 'log')
subplot(2,2,4), hold on, %swl
l_swl = ismember(sign(color), [1 -1 1], 'rows');
l_sedna_SwL_Frames = l_swl & sedna & frames;
l_sedna_SwL_NoFrames = l_swl & sedna & ~frames;
plot(sptFreq(l_sedna_SwL_Frames), 1./thresh(l_sedna_SwL_Frames), 'ms')
plot(sptFreq(l_sedna_SwL_NoFrames), 1./thresh(l_sedna_SwL_NoFrames), 'ms', 'markerfacecolor', 'm')
legend('Frames', 'No Frames', 'location', 'southwest')
set(gca, 'xscale', 'log', 'yscale', 'log')


%% (34) LOOKING FOR NON-STATIONARITIES: FIRING RATE VS. TIME

clear, clc
global cardVsIntBatchPath blpBatchPath
preprocessDTbatchData

MAKEPLOT = 1;
t = {};
spk = {};
rho = nan(length(out.dat),1);
p = nan(length(out.dat),1);
slope = nan(length(out.dat),1);
for a = 1:length(out.dat);
    if any(isnan(rawTRs(a,:)))
        continue
    end
    disp(a)
    DT = dtobj(out.fnames{a}{1});
    frameRate = DT.sum.exptParams.frame_rate;
    stimOn = DT.trial(:, DT.idx.flashOn); %in seconds
    stimOff = stimOn + (DT.trial(:, DT.idx.nFrames)./frameRate);
    t_on = mat2cell(stimOn, ones(length(stimOn),1));
    t_off = mat2cell(stimOff, ones(length(stimOff),1));
    counts = cellfun(@(on,off,ras) sum((ras>on)&(ras<=off)), t_on, t_off, DT.ras(:,DT.idx.spikes));
    [rho(a), p(a)] = corr(stimOn, counts, 'type', 'spearman');
    b = [stimOn, ones(length(stimOn),1)] \ counts;
    slope(a) = b(1);
    
    %package away for group comparisons.
    t{a} = stimOn;
    spk{a} = counts;
    
    %plot if need be
    if MAKEPLOT && p(a)<0.05
        figure
        plot((stimOn-stimOn(1))/60, counts, 'k.');
        title(sprintf('rho = %.3f, p = %.3f, slope = %.3f', rho(a), p(a), slope(a)));
        ylabel('counts')
        xlabel('time (min)')
    end
end


%% (12) SURROUND SUPPRESSION AND THRESHOLD RATIOS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear
global cardVsIntBatchPath blpBatchPath
preprocessDTbatchData

MAKEPLOT = 0;
suppIdx = nan(length(out.dat),1);
for a = 1:length(out.dat);
    diams = out.dat(a).grating.areasummation.stim;
    resp = out.dat(a).grating.areasummation.resp(:,1);
    [respAtPeak, idxToPrefDiam] = max(resp);
    [~, bigIdx] = max(diams);
    respToBiggest = resp(bigIdx);
    suppIdx(a) = (respAtPeak-respToBiggest) ./ respAtPeak;
    
    if MAKEPLOT
        figure, hold on,
        plot(diams, resp, 'ko-')
        plot(diams(idxToPrefDiam), respAtPeak, 'ko', 'markerfacecolor', 'r')
        plot(diams(bigIdx), respToBiggest, 'ko', 'markerfacecolor', 'y')
        title(sprintf('suppression idx = %.3f', suppIdx(a)))
    end
end

% (1) histogram of supp indicies
figure
hist(suppIdx, 25)
xlabel('Suppression Index')
ylabel('counts')



% (2) surround suppression vs. threshold ratio
l_validConds = ~commonExclusions;
filtTRs = rawTRs(l_validConds,:);
minFiltTRs = min(filtTRs, [],2);
l_nans = isnan(minFiltTRs);
valForNans = max(minFiltTRs) .* 1.2;
minFiltTRs(l_nans) = valForNans;
filtSuppIdx = suppIdx(l_validConds);
figure, hold on,
plot(filtSuppIdx, minFiltTRs, 'ko', 'markerfacecolor', 'k')
set(gca,'yscale', 'log')
ylim([min(minFiltTRs).*.98, valForNans.*1.2])
xlim([-0.05 1.05])
xlabel('suppression index')
ylabel('minimum TR')
[rho, p] = corr(filtSuppIdx(:), minFiltTRs(:), 'type', 'spearman');
title(sprintf('spearmans rho=%.3f, p=%.3f', rho, p))

% (3) surround suppression vs. CSI
filtCSIs = rawCSIs(l_validConds);
figure
plot(filtCSIs, filtSuppIdx, 'ko', 'markerfacecolor', 'k')
xlabel('CSI')
ylabel('Supression Index')
