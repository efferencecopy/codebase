%% (1)  BEGIN: DESGINATE A BATCH DATA FILE

% clear out the currrent workspace
fin

% designate a text file for the batch analysis:
SVSLM = 'SvsLM_31-Aug-2011.mat';
CARDVSINT = 'CardVsInt_06-May-2011.mat'; % the typical one: 'CardVsInt_06-May-2011.mat' or 'CardVsInt_NewCPCriteria.mat' for different CP low cutoff
BLPLFP = 'blp_DTbatch.mat';


global cardVsIntBatchPath blpBatchPath cardVsCardBatchPath
if ismac %charlie's laptop
    prefix = '~/LabStuff/Huskies/NexFiles/Charlie/Batch Data And Text Files/';
elseif strcmpi(license,  '367516') % when logged into cahass at wanprc on the server
    prefix = 'C:\NO BACKUP\NexFiles\Charlie\Batch Data And Text Files\';
else %assume all other pc's
    prefix = 'N:\NexFiles\Charlie\Batch Data And Text Files\';
end
cardVsIntBatchPath = [prefix, CARDVSINT];
blpBatchPath = [prefix, BLPLFP];
cardVsCardBatchPath = [prefix, SVSLM];


%% (2) SUMMARY OF PREFERED COLORS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
global cardVsIntBatchPath blpBatchPath %#ok<REDEF>
preprocessDTbatchData


lvmAxis = [1 -1 0] ./ norm([1 -1 0]);
SAxis = [0 0 .2];
l_validConds = ~(commonExclusions | out.errors(:,prefTieInd));
x = prefIsolum(l_validConds,:) * lvmAxis';
y = prefIsolum(l_validConds,:) * SAxis';
thetas = atan(y./x);
bins = linspace(-pi/4,pi/2,4);
figure;
subplot(1,2,1)
rose([thetas;thetas+pi], [bins(:); bins(:)+pi]);
set(get(gca, 'children'), 'linewidth', 2)
title('Prefered Color')
nTies = sum(out.errors(:, prefTieInd));
xlabel(sprintf('N = %d expts \n %d tie(s)', size(prefIsolum,1), nTies));

%bar chart of pref / int colors
cardinals = [1 -1 0; 0 0 1];
l_validConds = ~(commonExclusions | out.errors(:,cardTieInd));
cardText = sprintf('cardinals \n nValid = %d', sum(l_validConds));
for a = 1:size(cardinals,1);
    conditionList = ismember(sign(prefCards(l_validConds,:)), cardinals(a,:), 'rows');
    nCardinals(a) = sum(conditionList);
end

intermediates = [1 -1 1; 1 -1 -1];
l_validConds = ~(commonExclusions | out.errors(:,intTieInd));
intText = sprintf('intermediate \n nValid = %d', sum(l_validConds));
for a = 1:size(intermediates,1);
    conditionList = ismember(sign(prefInts(l_validConds,:)), intermediates(a,:), 'rows');
    nIntermediates(a) = sum(conditionList);
end
subplot(1,2,2)
bar([nCardinals, nIntermediates])
set(gca, 'xticklabel', num2str([cardinals; intermediates]));
maxy = max(get(gca, 'ylim'));
hold on,
plot([2.5 2.5], [0 maxy+3], 'k', 'linewidth', 2)
axis tight
ylim([0 maxy+3])
text(1, maxy+2, cardText);
text(2.8, maxy+2, intText)
ylabel('counts')


%% (3) ECCENTRICITY
%%%%%%%%%%%%%%%%%%%%%%

clear
global cardVsIntBatchPath blpBatchPath
preprocessDTbatchData

%Cartesian plot of RF centers
kCounter = 0;
sCounter = 0;
figure, hold on,
validConds = ~commonExclusions;
for a = find(validConds)';
    monkey = out.fnames{a}{1}(1); %the first letter of the file name
    switch lower(monkey)
        case 'k'
            plot(out.dat(a).expt.rfpos(1)/10, out.dat(a).expt.rfpos(2)/10, 'ko', 'markersize', 10*out.dat(a).grating.areasummation.prefsize)
            kCounter = kCounter+1;
        case 's'
            plot(out.dat(a).expt.rfpos(1)/10, out.dat(a).expt.rfpos(2)/10, 'co', 'markersize', 10*out.dat(a).grating.areasummation.prefsize)
            sCounter = sCounter+1;
    end
end
plot([-0.5, 0; 0.5, 0], [0, -0.5; 0, 0.5], 'r', 'linewidth', 2)
maxEcc = max(abs([get(gca, 'xlim'), get(gca, 'ylim')]));
ylim([-maxEcc, maxEcc])
xlim([-maxEcc, maxEcc])
axis square
title('RF Locations')
xlabel('Degrees of visual angle')
ylabel('Degrees of visual angle')
text(-5, 6, sprintf('K cells = %d \nS cells = %d', kCounter, sCounter))


% histogram of eccentricities
l_validConds = ~commonExclusions;
filtEcc = rawEccentricity(l_validConds);
rawMinTR = min(rawTRs, [], 2);
filtMinTR = rawMinTR(l_validConds);
meanTR = nanmean(filtMinTR);
valfornan = max(filtMinTR).*1.5;
l_nans = isnan(filtMinTR);
filtMinTR(l_nans) = valfornan;

figure,
set(gcf, 'position', [-1107 203 997 426])
subplot(1,2,1)
hold on,
hist(filtEcc)
plot(mean(filtEcc), 1, 'wv', 'markersize', 6, 'markerfacecolor', 'm')
title('eccentricity')
xlabel('Degrees of visual angle')
ylabel('counts')
hold off

subplot(1,2,2); %Ecc vs. TR
hold on,
plot(filtEcc, filtMinTR, 'ko')
plot(filtEcc(l_nans), filtMinTR(l_nans), 'ko', 'markerfacecolor', 'k')
[rho, p] = corr(filtEcc(:), filtMinTR(:), 'type', 'spearman');
[rhoNoNan, pNoNan] = corr(filtEcc(~l_nans), filtMinTR(~l_nans), 'type', 'spearman');
xlabel('Eccentricity (DVA)')
ylabel('Min TR')
title(sprintf('Rho = %.3f, p = %.3f \n rhoNoNan = %.3f, pNoNan = %.3f', rho, p, rhoNoNan, pNoNan));


%do the Eccentricities differ b/w monkey?
if ~any([sum(l_kali), sum(l_sedna)]==0)
    [h, p] = ttest2(rawEccentricity(l_sedna&validConds), rawEccentricity(l_kali&validConds))
    avgSednaEcc = mean(rawEccentricity(l_sedna&validConds))
    avgKaliEcc = mean(rawEccentricity(l_kali&validConds))
end


%% (4) POPULATION HISTOGRAM OF THRESHOLD RATIOS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, clc, close all
global cardVsIntBatchPath blpBatchPath
preprocessDTbatchData


l_validConds = ~commonExclusions;
filtTRs = rawTRs(l_validConds,:); %pull out the valid experiments
filtTRs = min(filtTRs, [], 2); %just pile all the data into a single vector
valForNan = 1.5 * max(filtTRs);
nNans = sum(isnan(filtTRs));
filtTRs(isnan(filtTRs)) = valForNan;
figure
edges = logspace(log10(0.8), log10(10), 20);
counts = histc(filtTRs, edges);
bar(edges, counts, 'histc')
set(gca, 'xscale', 'log','xlim', [0.6, 6])
set(gca, 'tickdir', 'out')
title(sprintf('n = %d, nNans = %d', length(filtTRs), nNans))
ylim([0, 20])


% Threshold ratios by monkey?
if ~any([sum(l_kali), sum(l_sedna)]==0)
    minTRs = min(rawTRs,[],2);
    kTRs = minTRs(l_kali & l_validConds);
    sTRs = minTRs(l_sedna & l_validConds);
    n_kTRs = length(kTRs);
    n_sTRs = length(sTRs);
    [h,p] = ttest2(sTRs, kTRs);
    avgKaliTR = nanmean(kTRs)
    semKaliTR = nanstd(kTRs) ./ sqrt(sum(~isnan(kTRs)))
    semSednaTR = nanstd(sTRs) ./ sqrt(sum(~isnan(sTRs)))
    avgSednaTR = nanmean(sTRs)
    kTRs(isnan(kTRs)) = max(minTRs) .* 1.5;
    sTRs(isnan(sTRs)) = max(minTRs) .* 1.5;
    edges = logspace(log10(0.8), log10(10), 20);
    k_counts = histc(kTRs, edges);
    s_counts = histc(sTRs, edges);
    
    figure, hold on,
    k = bar(edges, k_counts,'histc');
    set(k, 'facecolor', [0.4 0.4 0.4], 'edgeColor', 'k', 'lineWidth', 3)
    s = bar(edges, s_counts,'histc');
    set(s, 'facecolor', [0.2 0.2 1], 'facealpha', 0.2, 'edgeColor', 'b', 'lineWidth', 3);
    set(gca, 'xscale', 'log', 'xlim', [edges(1).*.9, edges(end).*1.1])
    delete(findobj('marker', '*'))
    legend('Kali', 'Sedna')
    text(2, 10, sprintf('Kali: %d \nSedna: %d', n_kTRs, n_sTRs));
    xlabel('Minimum TR')
    ylabel('Counts')
    title(sprintf('TRs by Monkey. One TR per neuron. p = %.3f', p))
end



%% (5) THRESHOLD RATIOS: (1)Card vs. Int (2) pref1 vs pref2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
global cardVsIntBatchPath blpBatchPath
preprocessDTbatchData

% compile the data
l_validConds = ~(commonExclusions);
cardTRs = rawCardTRs(l_validConds);
intTRs = rawIntTRs(l_validConds);
nBothNans = sum(isnan(cardTRs) & isnan(intTRs));
valForNan = max([cardTRs ; intTRs]) .* 1.2;
cardTRs(isnan(cardTRs)) = valForNan;
intTRs(isnan(intTRs)) = valForNan;
filtCSIs = rawCSIs(l_validConds);
filtColLum = (filtCSIs >= 0.5) & (filtCSIs < 2);
filtColOnly = filtCSIs >= 2;
filtLumOnly = filtCSIs < 0.5;

%compute CI's if desired
plotCI = 0;
CItype = 'hess'; % "boot" or "hess"
nStraps = 5000;
confInterval = .95;
if plotCI
    for a = 1:length(out.dat)
        fprintf('computing expt <%d> out of <%d>\n', a, length(out.dat));
        switch lower(CItype)
            case 'boot'
                CI = boot2ci(out.dat(a).m, out.dat(a).c, out.dat(a).expt, nStraps, confInterval);
            case 'hess'
                CI = hess2ci(out.dat(a).m, out.dat(a).c, confInterval);
        end
        %determine which is for card color
        sidx = ismember(sign(out.dat(a).expt.standColors), [0 0 1], 'rows');
        lvmidx = ismember(sign(out.dat(a).expt.standColors), [1 -1 0], 'rows');
        cardidx = sidx|lvmidx;
        cardCI(a,:) = [CI.lo(cardidx), CI.up(cardidx)];
        intCI(a,:) = [CI.lo(~cardidx), CI.up(~cardidx)];
    end
    filtCardCI = cardCI(l_validConds,:);
    filtIntCI = intCI(l_validConds,:);
end


%FIG 1, Card vs. Int TRs... grouped by CSI
figure,
hold on,
if plotCI
    plot(filtCardCI', [intTRs'; intTRs'], 'b')
    plot([cardTRs'; cardTRs'], filtIntCI', 'b')
end
plot(cardTRs(filtColLum), intTRs(filtColLum), 'ko', 'markerfacecolor', 'y', 'markersize', 8)
plot(cardTRs(filtColOnly), intTRs(filtColOnly), 'ko', 'markerfacecolor', 'r', 'markersize', 8)
plot(cardTRs(filtLumOnly), intTRs(filtLumOnly), 'ko', 'markerfacecolor', 'k', 'markersize', 8)
plot([0.2 1.5*valForNan], [0.2 1.5*valForNan], 'k')
xlim([0.5 1.5*valForNan])
ylim([0.5 1.5*valForNan])
set(gca, 'xscale', 'log', 'yscale', 'log')
xlabel('Cardinal TR')
ylabel('Intermediate TR')
text(2, 0.8, sprintf('n = %d \n n both-nans = %d', length(cardTRs), nBothNans))
title('Card Vs. Int By CSI')
axis square
hold off,


if ~any([sum(l_kali), sum(l_sedna)]==0)
    tmpCardTR = rawCardTRs;
    tmpCardTR(isnan(tmpCardTR)) = max(rawTRs(:).*1.2);
    tmpIntTR = rawIntTRs;
    tmpIntTR(isnan(tmpIntTR)) = max(rawTRs(:).*1.2);
    kCardTR = tmpCardTR(l_validConds & l_kali);
    kIntTR = tmpIntTR(l_validConds & l_kali);
    sCardTR = tmpCardTR(l_validConds & l_sedna);
    sIntTR = tmpIntTR(l_validConds & l_sedna);
    
    figure, hold on,
    plot(kCardTR, kIntTR, 'ko', 'markerfacecolor', 'k', 'markersize', 8)
    plot(sCardTR, sIntTR, 'go', 'markerfacecolor', 'g', 'markersize', 8)
    plot([0.2, max(rawTRs(:)).*1.2*1.5], [0.2, max(rawTRs(:)).*1.2*1.5], 'k:')
    xlim([0.5 max(rawTRs(:)).*1.2*1.1])
    ylim([0.5 max(rawTRs(:)).*1.2*1.1])
    set(gca, 'xscale', 'log', 'yscale', 'log')
    axis square
    legend('Kali', 'Sedna')
    xlabel('Card TR')
    ylabel('Int TR')
    title('Card Vs. Int By Monkey')
end


%% (6) THRESHOLD RATIOS AND NTs BY COLOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, clc, close all,
global cardVsIntBatchPath blpBatchPath
preprocessDTbatchData

%DEfine if the analysis should be restricted to a single color per neuron,
%or if both colors should be analyzed
BOTHTRS = 'two'; % 'one', or 'two'
EXCLUDENANS = 1;
COLORFILTER = 'mintr'; %'prefcard', or 'mintr'

%cycle through the expts and pull out the indicies to the TRs for each
%color direction
[sTR,lvmTR,swmTR,swlTR] = deal([]);
[sNT,lvmNT,swmNT,swlNT] = deal([]);
l_validConds = ~(commonExclusions | out.errors(:, cardMismatchInd) | out.errors(:,intMismatchInd));
for a = find(l_validConds)';
    NTs = out.dat(a).c.alpha;
    TRs = NTs ./ out.dat(a).m.alpha;
    switch BOTHTRS
        case 'one'
            if strcmpi(COLORFILTER, 'mintr')
                [~, idx] = min(TRs);
            elseif strcmpi(COLORFILTER, 'prefisolum')
                [~, idx] = ismember(sign(prefIsolum(a,:)), sign(out.dat(a).expt.standColors), 'rows');
                if ~idx; continue; end
            end
            tmpColor = sign(out.dat(a).expt.standColors(idx,:));
            NTs = NTs(idx);
            TRs = TRs(idx);
        case 'two'
            tmpColor = sign(out.dat(a).expt.standColors);
    end
    
    %deal with the cardinal TR
    if ismember([0 0 1], tmpColor, 'rows')
        [~,idx] = ismember([0 0 1], tmpColor, 'rows');
        if (EXCLUDENANS && (~isnan(TRs(idx)))) || ~EXCLUDENANS
            sTR(end+1,1) = TRs(idx);
            sNT(end+1,1) = NTs(idx);
        end
    elseif ismember([1 -1 0], tmpColor, 'rows')
        [~,idx] = ismember([1 -1 0], tmpColor, 'rows');
        if (EXCLUDENANS && (~isnan(TRs(idx)))) || ~EXCLUDENANS
            lvmTR(end+1,1) = TRs(idx);
            lvmNT(end+1,1) = NTs(idx);
        end
    end
    
    %deal with the intermediate TR
    if ismember([1 -1 1], tmpColor, 'rows')
        [~,idx] = ismember([1 -1 1], tmpColor, 'rows');
        if (EXCLUDENANS && (~isnan(TRs(idx)))) || ~EXCLUDENANS
            swlTR(end+1,1) = TRs(idx);
            swlNT(end+1,1) = NTs(idx);
        end
    elseif ismember([1 -1 -1], tmpColor, 'rows')
        [~,idx] = ismember([1 -1 -1], tmpColor, 'rows');
        if (EXCLUDENANS && (~isnan(TRs(idx)))) || ~EXCLUDENANS
            swmTR(end+1,1) = TRs(idx);
            swmNT(end+1,1) = NTs(idx);
        end
    end
end

% Box plot of threshold ratios
figure
hold on,
popMedian = nanmedian([sTR;lvmTR;swmTR;swlTR]);
plot([0.5 4.5], [popMedian, popMedian], 'k:')
tmp_TRs = [sTR; lvmTR; swmTR; swlTR];
group = char(repmat('Siso', length(sTR),1), repmat('LvM', length(lvmTR),1), repmat('SwM', length(swmTR),1), repmat('SwL', length(swlTR),1));
boxplot(tmp_TRs, group)
set(gca, 'yscale', 'log')
ylabel('Threshold Ratio')
%ylim([0.5, 4])
hold off
if strcmpi(BOTHTRS, 'one')
    [p_TR, atab_TR, stats_TR] = anova1(log(tmp_TRs), group, 'off');
    title(sprintf('Main effect of color, p = %g', p_TR))
end
    

% Box plot of neurometric thresholds
figure
hold on,
popMedian = nanmedian([sNT;lvmNT;swmNT;swlNT]);
plot([0.5 4.5], [popMedian, popMedian], 'k:')
tmp_NTs = [sNT; lvmNT; swmNT; swlNT];
group = char(repmat('Siso', length(sNT),1), repmat('LvM', length(lvmNT),1), repmat('SwM', length(swmNT),1), repmat('SwL', length(swlNT),1));
boxplot(tmp_NTs, group)
ylabel('Neurometric Threshold')
%ylim([0, 0.2])
hold off
if strcmpi(BOTHTRS, 'one')
    [p_NT, atab_NT, stats_NT] = anova1(tmp_NTs, group, 'off');
    title(sprintf('Main effect of color, p = %g', p_NT))
    figure
    [c,m,~,nms] = multcompare(stats_NT);
end


% % run a permutation test on the neurons that were sensitive in both color
% % directions tested
% [p, fstats] = TRpermtest();




%% (7) LINEAR MECHANISMS MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, clc
global cardVsIntBatchPath blpBatchPath
preprocessDTbatchData

%cycle through the valid neurons and compare the actual and predicted TRs
%of the intermediate color direction based on a linear model
l_validConds = ~(commonExclusions | out.errors(:, cardMismatchInd) | out.errors(:, intMismatchInd));
validCells = find(l_validConds);
for a = 1:length(validCells)
    cellInd = validCells(a);
    exptColors = out.dat(cellInd).expt.standColors;
    cardIdx = ismember(sign(exptColors), [1 -1 0], 'rows') | ismember(sign(exptColors), [0 0 1], 'rows');
    cardNT = out.dat(cellInd).c.alpha(cardIdx);
    cardColor = exptColors(cardIdx,:);
    cardUnit = cardColor ./ norm(cardColor);
    intColor = exptColors(~cardIdx,:);
    intUnit = intColor ./ norm(intColor);
    
    %make the prediction based on a linear model
    predIntNT(a) = cardNT / abs(intUnit * cardUnit');
    
    %determine the measured neurometric threshold
    measIntNT(a) = out.dat(cellInd).c.alpha(~cardIdx);
end

%plot the data
figure
hold on,
valForNan = max([predIntNT(:) ; measIntNT(:)]) .* 1.2;
tmpPredNT = predIntNT;
tmpPredNT(isnan(tmpPredNT)) = valForNan;
tmpMeasNT = measIntNT;
tmpMeasNT(isnan(tmpMeasNT)) = valForNan;
plot(tmpMeasNT, tmpPredNT, 'bo', 'markersize', 9, 'linewidth', 2)
xlabel('measured Int NT')
ylabel('predicted Int NT')
plot([0 valForNan.*1.2], [0, valForNan.*1.2], 'k-')
nNans = sum(isnan(measIntNT(:)+predIntNT(:))); %either the pred or meas TR is a nan
[p,h] = signtest(tmpPredNT, tmpMeasNT);
text(valForNan*.6, valForNan.*.2, sprintf('p=%.3f', p));
title(sprintf('n = %d, nNan = %d', sum(l_validConds), nNans));
xlim([0 valForNan*1.1])
ylim([0 valForNan*1.1])
axis square
hold off

%histogram of ratios (ratios are invariant to linear transformations of the
%color space). THIS ANALYSIS WILL ONLY ILLUSTRATE RATIOS FOR NEUORONS IN
%WHICH BOTH NEUROMETRIC THRESHOLDS WERE NON-NAN
figure
ratio = predIntNT./measIntNT;
ratio(isnan(ratio)) = [];
edges = logspace(log10(0.5), log10(3), 12);
counts = histc(ratio, edges);
bar(edges, counts, 'style', 'histc')
set(gca, 'xscale', 'log')
hold on,
plot(mean(ratio), 1, 'r*')
xlabel(sprintf('n = %d', length(ratio)))
[~,p] = ttest(ratio-1);
title(sprintf('p = %.5f', p));

%% (8) CRFs FOR INT DATA PROJECTED ONTO CARD AXIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



clear, clc
global cardVsIntBatchPath blpBatchPath
preprocessDTbatchData


MAKEPLOT = 0;
TRIALLENGTH = 0.666;  % in seconds
CARD = 0; %grouping variable for the GLM
INT = 1;
l_validConds = ~(commonExclusions | (out.errors(:, neuroThresh1Ind) & out.errors(:,neuroThresh2Ind)));
exptList = find(l_validConds);
[boffset, bcontrast, bcolor, binteract, pcontrast, pinteract, pcolor] = deal(nan(1, length(l_validConds)));
[popCounts, popYhat, popContrast, popCards] = deal([]);
for a = 1:length(exptList)
    ex = exptList(a);
    
    %deal with the cardinal light data (counts not rates)
    s_idx = ismember(sign(out.dat(ex).expt.standColors), [0 0 1], 'rows');
    lvm_idx = ismember(sign(out.dat(ex).expt.standColors), [1 -1 0], 'rows');
    card_Idx = s_idx | lvm_idx;
    if sum(card_Idx) ~= 1; keyboard; end
    counts = [];
    contrasts = [];
    colors = [];
    for cntrst = 1:length(out.dat(ex).c.crfIn{card_Idx})
        tmp = out.dat(ex).c.crfIn{card_Idx}{cntrst}.*TRIALLENGTH;
        counts = [counts; tmp(:)];
        contrasts = [contrasts; repmat(out.dat(ex).expt.norms{card_Idx}(cntrst), length(tmp), 1)];
        colors = [colors; repmat(CARD, length(tmp), 1)];
    end
    
    %deal with the intermediate light data
    int_Idx = ~card_Idx;
    for cntrst = 1:length(out.dat(ex).c.crfIn{int_Idx})
        tmp = out.dat(ex).c.crfIn{int_Idx}{cntrst}.*TRIALLENGTH;
        counts = [counts; tmp(:)];
        cardColor = out.dat(ex).expt.standColors(card_Idx,:);
        cardUnit = cardColor ./ norm(cardColor);
        intColor = out.dat(ex).expt.standColors(int_Idx,:);
        intUnit = intColor./norm(intColor);
        intLMS = intUnit .* out.dat(ex).expt.norms{int_Idx}(cntrst);
        projOntoCard = abs(intLMS * cardUnit(:));
        contrasts = [contrasts; repmat(projOntoCard, length(tmp), 1)];
        colors = [colors; repmat(INT, length(tmp), 1)];
    end
    
    
    %run the possion regression with three terms: Color, Contrast,
    %Interaction...
    counts = round(counts); %multiplying by .666 is a hack b/c I don't compute the trial length exactly. As such, some of the counts aren't integers.
    interaction = contrasts.*colors;
    predictors = [contrasts, colors, interaction];
    [beta, ~, stats] = glmfit(predictors, counts, 'poisson', 'link', 'log');
    boffset(ex) = beta(1);
    bcontrast(ex) = beta(2);
    bcolor(ex) = beta(3);
    binteract(ex) = beta(4);
    pcontrast(ex) = stats.p(2);
    pcolor(ex) = stats.p(3);
    pinteract(ex) = stats.p(4);
    
    
    %compute some things for a population figure. normalize by card
    %threshold so that the CRFs are computed in threshold units.
    popCounts = [popCounts; counts];
    popYhat = [popYhat; glmval(beta, predictors, 'log')];
    popContrast = [popContrast; contrasts(:) ./ out.dat(ex).m.alpha(card_Idx)]; %scaled to card thresh
    popCards = [popCards; colors==CARD]; %so that I can plot card and int seperately
    
    maxContrast = 3.5;
    nInterps = 400;
    interpContrast = linspace(0, maxContrast, nInterps)' .* out.dat(ex).m.alpha(card_Idx);
    cardInterp(a,:) = glmval(beta, [interpContrast, ones(nInterps,1).*CARD, (interpContrast.*(ones(nInterps,1).*CARD))], 'log');
    intInterp(a,:) = glmval(beta, [interpContrast, ones(nInterps,1).*INT, (interpContrast.*(ones(nInterps,1).*INT))], 'log');
    
    %plot the raw data along with the glm fit data. Poisson regresson
    %assumes a linear relationship b/w contrast and the log(counts) so plot
    %on a semilogy axis
    if MAKEPLOT
        figure,hold on,
        l_card = colors==CARD;
        plot(contrasts(l_card), counts(l_card), 'bo')
        plot(contrasts(~l_card), counts(~l_card), 'ko')
        yhat = glmval(beta, predictors, 'log');
        plot(contrasts(l_card), yhat(l_card), 'b-')
        plot(contrasts(~l_card), yhat(~l_card), 'k-')
        legend('Cardinal', 'Intermediate')
        title(sprintf('Cntrst = %.3f, Color = %.3f, Interact = %.3f', pcontrast(ex), pcolor(ex), pinteract(ex)))
        xlabel('Projection Onto Cardinal Direction')
        ylabel('Counts')
        set(gca, 'xscale', 'log')
        set(gcf, 'name', sprintf('%s, exp = %d', out.fnames{ex}{1}, a))
        ylim([0, max(counts).*1.1])
        
        
        %         %trying to see if the interpolation thing works...
        %         cardPT = out.dat(ex).m.alpha(card_Idx);
        %         plot(linspace(0, maxContrast, nInterps).*cardPT, cardInterp(a,:), 'cv')
        %         plot(linspace(0, maxContrast, nInterps).*cardPT, intInterp(a,:), 'gv')
        hold off
    end
end


% a lame way to summarize the general linear model fits. Basically plotting
% the beta values and demonstrating that as a population, the betas
% associated with contrast and interaction (the 2 slope determinates) are
% significantly different than zero
figure
subplot(2,2,1), hold on, % beta0... offset
hist(boffset)
plot(mean(boffset), 5, 'vm', 'markerfacecolor', 'm')
[~, pOffset] = ttest(boffset);
title(sprintf('pOffset = %.3f', pOffset));
axis tight
hold off
subplot(2,2,2), hold on, %beta1... contrast
hist(bcontrast)
plot(mean(bcontrast), 5, 'vm', 'markerfacecolor', 'm')
[~, pContrast] = ttest(bcontrast);
title(sprintf('pContrast = %.3f', pContrast));
axis tight
hold off
subplot(2,2,3), hold on, %beta2... color
hist(bcolor)
plot(mean(bcolor), 5, 'vm', 'markerfacecolor', 'm')
[~, pColor] = ttest(bcolor);
title(sprintf('pColor = %.3f', pColor));
axis tight
hold off
subplot(2,2,4), hold on,  %beta3... interaction
hist(binteract)
plot(mean(binteract), 5, 'vm', 'markerfacecolor', 'm');
[~, pInteract] = ttest(binteract);
title(sprintf('pInteract = %.3f', pInteract))
axis tight
hold off


%now show a population CRF in threshold units.
figure
subplot(1,2,1), hold on, %yhat vs contrast
popCards = logical(popCards);
plot(popContrast(popCards), popYhat(popCards), 'b.')
plot(linspace(0, maxContrast, nInterps), mean(cardInterp,1), 'b', 'linewidth', 3)
plot(popContrast(~popCards), popYhat(~popCards), 'k.')
plot(linspace(0,maxContrast, nInterps), mean(intInterp,1), 'k', 'linewidth', 3)
ylabel('yHat')
xlabel('Contrast (x card thresh)')
ylim([0 120])
subplot(1,2,2), hold on %counts vs contrast
contrastBins = linspace(0,maxContrast, 15);
cardCountsPerBin = nan(1,length(contrastBins)-1);
intCountsPerBin = nan(1,length(contrastBins)-1);
for a = 1:length(contrastBins)-1;
    l_inBin = (popContrast>contrastBins(a)) & (popContrast<=contrastBins(a+1));
    cardCountsPerBin(a) = median(popCounts(l_inBin & popCards));
    intCountsPerBin(a) = median(popCounts(l_inBin & ~popCards));
end
plot(contrastBins(1:end-1), cardCountsPerBin, 'b', 'linewidth', 3)
plot(contrastBins(1:end-1), intCountsPerBin, 'k', 'linewidth', 3)
plot(linspace(0, maxContrast, nInterps), median(cardInterp,1), 'b')
plot(linspace(0,maxContrast, nInterps), median(intInterp,1), 'k')
ylim([0, 100])
ylabel('counts')
xlabel('contrast (x card thres)')

% histograms of ratios: (B2 vs. B2+B4). B2 corresponds to CARD alone
figure,
subplot(1,2,1), hold on,
b_card = bcontrast;
b_int = bcontrast+binteract;
alpha = 0.05;
l_sig = pcontrast<alpha & pinteract<alpha;
maxXY = max(b_int);
plot(b_card, b_int, 'bo')
plot(b_card(l_sig), b_int(l_sig), 'bo', 'markerfacecolor', 'b')
plot([0, maxXY], [0, maxXY],'k')
axis tight
subplot(1,2,2)
hist(b_card./b_int, 30)

%% (9) SLOPES
%
%   (1) neurometric vs. psychometric (slope ratio)
%   (2) slope ratio by color (with attention to S-iso vs. L-M)
%   (3) correlation b/w TR and SR?
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, clc
global cardVsIntBatchPath blpBatchPath
preprocessDTbatchData

l_validExpts = ~(commonExclusions);% | (out.errors(:, neuroThresh1Ind) & out.errors(:, neuroThresh2Ind)));
exptList = find(l_validExpts);
[SR_S,SR_LvM,SR_SwM,SR_SwL] = deal([]); %initialize vectors for SRs
[neuroSlope_S, neuroSlope_LvM, neuroSlope_SwM, neuroSlope_SwL] = deal([]); %initialze vectors for neuro slopes
[psySlope_S, psySlope_LvM, psySlope_SwM, psySlope_SwL] = deal([]); %initialze vectors for psy slopes
rawSRs = nan(length(exptList), 2);
filtCSIs = nan(length(exptList),1);
for a = 1:length(exptList)
    ex = exptList(a);
    for clr = 1:size(out.dat(a).expt.standColors, 1)
        
        % exclude color conditions where the neuron was insensitive
        if isnan(out.dat(ex).c.alpha(clr));
            continue
        end
        
        %determine some basic info
        color = sign(out.dat(ex).expt.standColors(clr,:));
        neuroSlope = out.dat(ex).c.beta(clr);
        psySlope = out.dat(ex).m.beta(clr);
        ratio = neuroSlope ./ psySlope;
        
        %a hack to exclude rediculously steep slopes.
        if any([neuroSlope, psySlope]>150)
            disp('************ EXCLUDING DATA *****************')
            continue
        end
        
        
        %compile the SRs by color direction (lumping the 2 SRs for each
        %expt in the same analysis)
        if all(color == [0 0 1]) %S-iso
            SR_S(end+1) = ratio;
            neuroSlope_S(end+1) = neuroSlope;
            psySlope_S(end+1) = psySlope;
        elseif all(color == [1 -1 0]) %LvM
            SR_LvM(end+1) = ratio;
            neuroSlope_LvM(end+1) = neuroSlope;
            psySlope_LvM(end+1) = psySlope;
        elseif all(color == [1 -1 -1]) %SwM
            SR_SwM(end+1) = ratio;
            neuroSlope_SwM(end+1) = neuroSlope;
            psySlope_SwM(end+1) = psySlope;
        elseif all(color == [1 -1 1]) %SwL
            SR_SwL(end+1) = ratio;
            neuroSlope_SwL(end+1) = neuroSlope;
            psySlope_SwL(end+1) = psySlope;
        else
            error('Unidentified color')
        end
        
        %now compile a <nExpts x 2> big matrix to house all the raw SRs
        rawSRs(a, clr) = ratio;
        filtCSIs(a,1) = rawCSIs(ex);
    end
end


% (1) simple histogram of Slope Ratios
allSRs = rawSRs(:);
l_nans = isnan(allSRs);
valForNan = max(allSRs) .* 3;
allSRs(l_nans) = valForNan;
edges = logspace(log10(min(allSRs).*.98), log10(max(allSRs).*1.15), 20);
n = histc(allSRs, edges);
figure, hold on,
bar(edges, n, 'type', 'histc')
avg = mean(allSRs(~l_nans));
plot(avg, 7, 'vw', 'markerfacecolor', 'w')
set(gca, 'xscale', 'log', 'fontsize', 14)
xlabel('slope ratio')
ylabel('count')
title(sprintf('nanval = %.3f', valForNan));
hold off


% (2) All types of slopes broken down by color
figure
subplot(1,3,1),hold on, %slope Ratios
plot(ones(length(SR_S)), SR_S, 'bo')
plot(ones(length(SR_LvM)).*2, SR_LvM, 'ro')
plot(ones(length(SR_SwL)).*3, SR_SwL, 'mo')
plot(ones(length(SR_SwM)).*4, SR_SwM, 'go')
plot([0.6 1.4], repmat(median(SR_S),1,2), 'b', 'linewidth', 3)
plot([1.6 2.4], repmat(median(SR_LvM),1,2), 'r', 'linewidth', 3)
plot([2.6 3.4], repmat(median(SR_SwL),1,2), 'm', 'linewidth', 3)
plot([3.6 4.4], repmat(median(SR_SwM),1,2), 'g', 'linewidth', 3)
xlim([0,5])
ylabel('Slope Ratio')
xlabel('Color')
pCard = ranksum(SR_S, SR_LvM);
pInt = ranksum(SR_SwL, SR_SwM);
title(sprintf('P_{card} = %.3f, P_{int} = %.3f', pCard, pInt))
hold off
subplot(1,3,2), hold on, %neurometric slopes
plot(ones(length(neuroSlope_S)), neuroSlope_S, 'bo')
plot(ones(length(neuroSlope_LvM)).*2, neuroSlope_LvM, 'ro')
plot(ones(length(neuroSlope_SwL)).*3, neuroSlope_SwL, 'mo')
plot(ones(length(neuroSlope_SwM)).*4, neuroSlope_SwM, 'go')
plot([0.6 1.4], repmat(median(neuroSlope_S),1,2), 'b', 'linewidth', 3)
plot([1.6 2.4], repmat(median(neuroSlope_LvM),1,2), 'r', 'linewidth', 3)
plot([2.6 3.4], repmat(median(neuroSlope_SwL),1,2), 'm', 'linewidth', 3)
plot([3.6 4.4], repmat(median(neuroSlope_SwM),1,2), 'g', 'linewidth', 3)
xlim([0,5])
ylabel('Neurometric Slope')
xlabel('Color')
pCard = ranksum(neuroSlope_S, neuroSlope_LvM);
pInt = ranksum(neuroSlope_SwL, neuroSlope_SwM);
title(sprintf('P_{card}=%.3f, P_{int}=%.3f', pCard, pInt))
hold off
subplot(1,3,3), hold on, %psychometric slopes
plot(ones(length(psySlope_S)), psySlope_S, 'bo')
plot(ones(length(psySlope_LvM)).*2, psySlope_LvM, 'ro')
plot(ones(length(psySlope_SwL)).*3, psySlope_SwL, 'mo')
plot(ones(length(psySlope_SwM)).*4, psySlope_SwM, 'go')
plot([0.6 1.4], repmat(median(psySlope_S),1,2), 'b', 'linewidth', 3)
plot([1.6 2.4], repmat(median(psySlope_LvM),1,2), 'r', 'linewidth', 3)
plot([2.6 3.4], repmat(median(psySlope_SwL),1,2), 'm', 'linewidth', 3)
plot([3.6 4.4], repmat(median(psySlope_SwM),1,2), 'g', 'linewidth', 3)
xlim([0,5])
ylabel('Psychometricmetric Slope')
xlabel('Color')
pCard = ranksum(psySlope_S, psySlope_LvM);
pInt = ranksum(psySlope_SwL, psySlope_SwM);
title(sprintf('P_{card}=%.3f, P_{int}=%.3f', pCard, pInt))
hold off


% (3) Diffs in neurometric and psychometric slopes for L-M and S
sDiff = psySlope_S - neuroSlope_S;
lvmDiff = psySlope_LvM - neuroSlope_LvM;
p = ranksum(sDiff(:), lvmDiff(:));
figure, hold on,
plot([0.7, 1.3], repmat(median(sDiff), 2,1), 'b', 'linewidth', 3)
plot([1.7, 2.3], repmat(median(lvmDiff), 2,1), 'r', 'linewidth', 3)
plot(ones(length(sDiff)), sDiff, 'bo', 'markerfacecolor', 'b')
plot(ones(length(lvmDiff)).*2, lvmDiff, 'ro', 'markerfacecolor', 'r')
legend('S Diff', 'L-M Diff')
title(sprintf('(Psy - Neuro) slopes \n Wilcoxon p=%.3f', p));
ylabel('Diff in slope')
hold off

% (4) SRs for CSIs
allSRs = max(rawSRs,[], 2);
l_nans = isnan(allSRs);
l_CO = filtCSIs >= 2;
l_CL = (filtCSIs >= 0.5) & (filtCSIs < 2);
l_LO = filtCSIs < 0.5;
[rho, p] = corr(allSRs(~l_nans), filtCSIs(~l_nans), 'type', 'spearman');
bigSR = max(allSRs).*1.2;
smallSR = min(allSRs).*0.90;
figure, hold on,
plot(filtCSIs(l_CO&~l_nans), allSRs(l_CO&~l_nans), 'ko', 'markerfacecolor', 'r')
plot(filtCSIs(l_CL&~l_nans), allSRs(l_CL&~l_nans), 'ko', 'markerfacecolor', 'y')
plot(filtCSIs(l_LO&~l_nans), allSRs(l_LO&~l_nans), 'ko', 'markerfacecolor', 'k')
plot([0.5, 2; 0.5, 2], [smallSR, smallSR; bigSR, bigSR], 'k--')
set(gca, 'xscale', 'log', 'yscale', 'log')
ylim([smallSR, bigSR])
title(sprintf('speaman''s rho = %.3f, p=%.3f', rho, p))
xlabel('CSI')
ylabel('Max Slope Ratio')
hold off



%% (10) COLOR TUNING WIDTH VS. THRESHOLD RATIOS (ISOLUM ONLY)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, clc, close all
global cardVsIntBatchPath blpBatchPath
preprocessDTbatchData

%find the ratio of rate to the prefered and orthogonal gratings.
sScaleFactor = 5; %necessary for finding the orth. Int color dir
l_validConds = ~(commonExclusions | (out.errors(:, neuroThresh1Ind) & out.errors(:,neuroThresh2Ind)));
exptList = find(l_validConds);
[raisedExp, opRatio, planeStatDT, circVarDT] = deal(nan(length(exptList),1));
plotFigs = 0;
PLOTGTPOP = 1;
for a = 1:length(exptList)
    %compute the o/p ratio
    expt = exptList(a);
    colors = out.dat(expt).grating.color.colors;
    colors(:,3) = colors(:,3)./sScaleFactor;
    l_isolum = ~(colors * [1;1;0]);
    prefRate = max(out.dat(expt).grating.color.colresp(l_isolum, 1));
    if prefRate<3;
        disp('DT no resp')
        continue
    end
    
    if sScaleFactor == 5 %this is the only way there will be orthogonals for the intermediates
        if prefRate == 0;
            opRatio(a) = nan;
        else
            prefIdx = (out.dat(expt).grating.color.colresp(:,1) == prefRate) & l_isolum;
            [~, orthIdx] = min(abs(colors(l_isolum,:) * colors(prefIdx,:)'));
            isolumColors = colors(l_isolum,:);
            orthIdx = ismember(colors, isolumColors(orthIdx,:), 'rows');
            orthRate = out.dat(expt).grating.color.colresp(orthIdx,1);
            opRatio(a) = orthRate./prefRate;
        end
    end
    
    
    %compute the plane stat
    if prefRate == 0
        planeStatDT(a) = nan;
    else
        planeStatDT(a) = sum(out.dat(expt).grating.color.colresp(l_isolum,1)) ./ (sum(l_isolum).*prefRate);
    end
    
    %compute the circular variance
    isolumColors = colors(l_isolum,:); %colors has been scaled....
    x = [isolumColors(:,1:2)] * ([1 -1]./norm([1 -1]))';
    y = [isolumColors(:,3)];
    theta = atan(y./x);
    isolumRates = out.dat(expt).grating.color.colresp(l_isolum,1);
    if sum(isolumRates) == 0;
        circVarDT(a) = nan;
    else
        circVarDT(a) = circleVar(isolumRates, theta);
    end
    
    
    %compute the exponent of a raised sinusoid
    resp = out.dat(expt).grating.color.colresp(l_isolum,1);
    x = [isolumColors(:,1:2)] * ([1 -1]./norm([1 -1]))';
    y = [isolumColors(:,3)];
    theta = atan(y./x);
    [prefTheta, raisedExp(a), gain] = raisedCos(theta, resp);
    
    if plotFigs && any(strcmp(out.fnames{expt}{1}, {'K070510002', 'S060809002', 'S022111004', 'S012910004'}));
        
        [theta, ind] = sort(theta);
        isolumRates = isolumRates(ind);
        polar([theta; theta+pi; theta(1)], [isolumRates; isolumRates; isolumRates(1)], '-bo')
        set(get(gca, 'children'), 'markerfacecolor', 'b')
        title(sprintf('expnt: %.3f, CV: %.3f plane:%.3f', raisedExp(a), circVarDT(a), planeStatDT(a)));
        xlabel(sprintf(out.fnames{expt}{1}))
        hold on,
        
        % for the rasied exponents
        xx = 0:0.001:(2*pi);
        model = abs(cos(xx-prefTheta).^raisedExp(a));
        model = gain .* model;
        polar(xx, model, 'k-')
        polar(prefTheta, max(resp), 'r+')
        hold off
        
        keyboard
        cla
    end
end

%load in a mass of GT data and find the sphere stat for them as well
gt = load('/Users/charliehass/LabStuff/NexFiles/Charlie/Batch Data And Text Files/allGTFiles.mat');
gt = gt.out;
if PLOTGTPOP
    [planeStatGT, circVarGT, raisedExpGT, rawModRatioGT, rawRaisedExpGT, rawPlaneStatGT] = deal(nan(numel(gt.dat), 1));
    for a = 1:numel(gt.dat)
        try
            nTrials = gt.dat(a).gratings.color.colresp(:,3);
            if min(nTrials)<6
                continue
            end
            
            colors = gt.dat(a).gratings.color.colors;
            colors(:,3) = colors(:,3)./sScaleFactor;
            l_isolum = ~(colors * [1;1;0]);
            prefRate = max(gt.dat(a).gratings.color.colresp(l_isolum, 1));
            
            if prefRate<3;
                disp('GT no resp')
                continue
            end
            
            %compute the plane stat
            if prefRate == 0
                planeStatGT(a) = nan;
            else
                planeStatGT(a) = sum(gt.dat(a).gratings.color.colresp(l_isolum,1)) ./ (sum(l_isolum).*prefRate);
            end
            
            %compute the circular variance
            isolumColors = colors(l_isolum,:); %colors has been scaled....
            x = [isolumColors(:,1:2)] * ([1 -1]./norm([1 -1]))';
            y = [isolumColors(:,3)];
            theta = atan(y./x);
            isolumRates = gt.dat(a).gratings.color.colresp(l_isolum,1);
            if all(isolumRates == 0);
                circVarGT(a) = nan;
            else
                circVarGT(a) = circleVar(isolumRates, theta);
            end
            
            %find the best fitting raised sinusoid
            resp = gt.dat(a).gratings.color.colresp(l_isolum,1);
            x = [isolumColors(:,1:2)] * ([1 -1]./norm([1 -1]))';
            y = [isolumColors(:,3)];
            theta = atan(y./x);
            if all(isolumRates == 0)
                raisedExpGT(a) = nan;
            else
                [prefTheta, raisedExpGT(a), gain] = raisedCos(theta, resp);
            end
            
            %store the modulation ratio
            rawModRatioGT(a) = gt.dat(a).gratings.modulationratio;
            rawRaisedExpGT(a) = raisedExpGT(a);
            rawPlaneStatGT(a) = planeStatGT(a);
            
        catch
            continue
        end
    end
    planeStatGT(isnan(planeStatGT)) = [];
    circVarGT(isnan(circVarGT)) = [];
    raisedExpGT(isnan(raisedExpGT)) = [];
end


% ***************************************** %
%           O/P Ratio Figures
% ***************************************** %
if sScaleFactor == 5;
    figure
    set(gcf, 'position', [9         404        1385         380])
    subplot(1,2,1), hold on,
    %find the linear predictions across a range of phase shifts
    phi = linspace(0, pi, 10000);
    phi(end) = [];
    colors = gt.dat(7).gratings.color.colors;
    colors(:,3) = colors(:,3)./sScaleFactor;
    l_isolum = ~(colors * [1;1;0]);
    colors = colors(l_isolum, :);
    pred_opRatio = nan(size(phi));
    for a = 1:numel(phi)
        mecVec = (cos(phi(a)).*[1/sqrt(2); -1/sqrt(2);0])+(sin(phi(a)).*[0;0;1]);
        R = abs(colors * mecVec);
        [prefResp, prefIdx] = max(R);
        prefColor = colors(prefIdx,:);
        [~, orthIdx] = min(abs(prefColor * colors'));
        orthResp = R(orthIdx);
        pred_opRatio(a) = orthResp ./ prefResp;
    end
    edges1 = linspace(0,1,20);
    opRatioCounts = histc(opRatio, edges1);
    edges2 = linspace(0,1,200);
    pred_opRatioPDF = histc(pred_opRatio, edges2)./numel(pred_opRatio);
    [ax h1, h2] = plotyy(edges1, opRatioCounts, edges2, pred_opRatioPDF, 'bar', 'plot');
    set(ax, 'xLim', [0,1]);
    set(h2, 'color', 'r', 'linewidth', 2)
    set(ax(2), 'ycolor', 'r')
    set(get(ax(1), 'ylabel'), 'string', 'Raw Data', 'fontsize', 14, 'fontname', 'times')
    set(get(ax(2), 'ylabel'), 'string', 'Linear Prediction', 'fontsize', 14, 'fontname', 'times')
    xlabel('Orth/Pref Ratio')
    hold off
    
    subplot(1,2,2), hold on,
    filtTRs = rawTRs(l_validConds,:);
    filtMinTRs = min(filtTRs, [], 2);
    l_nans = isnan(filtMinTRs);
    valForNan = 1.5*max(filtMinTRs);
    filtMinTRs(l_nans) = valForNan;
    l_zeros = opRatio == 0;
    valForZeros = min(opRatio(opRatio>0)) .* 0.98;
    tmpOPRatio = sum([opRatio(:), l_zeros.*valForZeros], 2);
    plot(tmpOPRatio, filtMinTRs, 'o', 'markerfacecolor', 'b')
    plot(tmpOPRatio(l_nans), filtMinTRs(l_nans), 'ko', 'markerfacecolor', 'k')
    plot(tmpOPRatio(l_zeros), filtMinTRs(l_zeros), 'go', 'markerfacecolor', 'g')
    set(gca, 'yscale', 'log')
    ylim([.98*min(filtMinTRs), 1.1*valForNan])
    xlim([min(tmpOPRatio)*0.9, max(tmpOPRatio)*1.1])
    [rho, p] = corr(opRatio(~l_nans), filtMinTRs(~l_nans), 'type', 'spearman');
    xlabel('Orth/Pref rate for color')
    ylabel('Min TR')
    title(sprintf('Rho = %.3f, p = %.3f', rho, p));
    hold off
end


% ***************************************** %
%           Circular Variance
% ***************************************** %
figure
set(gcf, 'position', [9         404        1385         380])
subplot(1,2,1), hold on,
%find the prediction based on a linear neuron w/phase shifts
phi = linspace(0, pi, 5000);
phi(end) = [];
colors = gt.dat(7).gratings.color.colors;
colors(:,3) = colors(:,3)./sScaleFactor;
l_isolum = ~(colors * [1;1;0]);
colors = colors(l_isolum, :);
x = [colors(:,1:2)] * ([1 -1]./norm([1 -1]))';
y = [colors(:,3)];
theta = atan(y./x);
pred_CV = nan(size(phi));
for a = 1:numel(phi)
    mecVec = (cos(phi(a)).*[1/sqrt(2); -1/sqrt(2);0])+(sin(phi(a)).*[0;0;1]);
    R = abs(colors * mecVec);
    pred_CV(a) = circleVar(R, theta);
end
%do the plotting
edges1 = linspace(0,1,20);
CVCounts = histc(circVarDT, edges1);
edges2 = linspace(0,1,200);
pred_CVPDF = histc(pred_CV, edges2)./numel(pred_CV);
if PLOTGTPOP
    CVCounts_GT = histc(circVarGT, edges1);
    bar(edges1, CVCounts_GT, 'facecolor', 'k')
end
[ax, h1, h2] = plotyy(edges1, CVCounts, edges2, pred_CVPDF, @bar, @plot);
set(ax, 'xLim', [0 ,1]);
%set(ax(1), 'ylim', [0 60], 'ytick', [0:5:15])
set(h2, 'color', 'r', 'linewidth', 2)
set(get(h1, 'children'), 'faceColor', [.40, .40, .40])
set(ax(2), 'ycolor', 'r')
set(get(ax(1), 'ylabel'), 'string', 'Raw Data', 'fontsize', 14, 'fontname', 'times')
set(get(ax(2), 'ylabel'), 'string', 'Linear Prediction', 'fontsize', 14, 'fontname', 'times')
xlabel('Circular Variance')
hold off

subplot(1,2,2), hold on,
filtTRs = rawTRs(l_validConds,:);
filtTRs = min(filtTRs,[],2);
l_nans = isnan(filtTRs)  | isnan(circVarDT);
valfornan = max(filtTRs)*1.1;
filtTRs(l_nans) = valfornan;
plot(circVarDT, filtTRs, 'bo')
plot(circVarDT(l_nans), filtTRs(l_nans), 'ko');
set(gca, 'yscale', 'log')
[rho, p] = corr(circVarDT(~l_nans), filtTRs(~l_nans), 'type', 'spearman');
title(sprintf('rho = %.3f, p = %.3f', rho, p))
hold off





% ***************************************** %
%           Iso-lum plane stat
% ***************************************** %
figure,
set(gcf, 'position', [9         404        1385         380])
subplot(1,2,1), hold on
%figuring out the prediction based on sinusoids for plane stat
phi = linspace(0, pi, 5000);
phi(end) = [];
colors = gt.dat(7).gratings.color.colors;
colors(:,3) = colors(:,3)./sScaleFactor;
l_isolum = ~(colors * [1;1;0]);
colors = colors(l_isolum, :);
predStat = nan(size(phi));
for a = 1:length(phi)
    mecVec = (cos(phi(a)).*[1/sqrt(2); -1/sqrt(2);0])+(sin(phi(a)).*[0;0;1]);
    R = abs(colors * mecVec);
    maxR = max(R);
    predStat(a) = sum(R) ./ (maxR .* numel(R));
end
edges1 = linspace(0,1,20);
planeStatDTCounts = histc(planeStatDT, edges1);
edges2 = linspace(0,1,200);
predStatPDF = histc(predStat, edges2)./numel(predStat);
if PLOTGTPOP
    planeStatGTCounts = histc(planeStatGT, edges1);
    bar(edges1, planeStatGTCounts, 'facecolor', 'k')
end
[ax, h1, h2] = plotyy(edges1, planeStatDTCounts, edges2, predStatPDF, @bar, @plot);
set(ax, 'xLim', [0.15,1]);
%set(ax(1), 'ylim', [0, 70])
set(h2, 'color', 'r', 'linewidth', 2)
set(ax(2), 'ycolor', 'r')
set(get(h1, 'children'), 'faceColor', [.40, .40, .40])
set(get(ax(1), 'ylabel'), 'string', 'Raw Data', 'fontsize', 14, 'fontname', 'times')
set(get(ax(2), 'ylabel'), 'string', 'Linear Prediction', 'fontsize', 14, 'fontname', 'times')
xlabel('Isolumplane Statistic')
hold off

subplot(1,2,2), hold on,
filtTRs = rawTRs(l_validConds,:);
filtMinTRs = min(filtTRs, [], 2);
l_nanTRs = isnan(filtMinTRs);
valForNan = 1.5*max(filtMinTRs);
filtMinTRs(l_nans) = valForNan;
plot(planeStatDT, filtMinTRs, 'o', 'markerfacecolor', 'b')
plot(planeStatDT(l_nans), filtMinTRs(l_nans), 'ko', 'markerfacecolor', 'k')
set(gca, 'yscale', 'log')
ylim([.98*min(filtMinTRs), 1.1*valForNan])
xlim([min(planeStatDT)*0.9, max(planeStatDT)*1.1])
l_nan_twiAndTR = l_nanTRs | isnan(raisedExp);
[rho, p] = corr(planeStatDT(~l_nan_twiAndTR), filtMinTRs(~l_nan_twiAndTR), 'type', 'spearman');
xlabel('Isolumplane Statistic')
ylabel('Min TR')
title(sprintf('Rho = %.3f, p = %.3f', rho, p));
hold off


% ***************************************** %
%           Raised Exponent
% ***************************************** %
%find the prediction based on a linear neuron w/phase shifts
phi = linspace(0, pi, 50);
phi(end) = [];
colors = gt.dat(7).gratings.color.colors;
colors(:,3) = colors(:,3)./sScaleFactor;
l_isolum = ~(colors * [1;1;0]);
colors = colors(l_isolum, :);
x = [colors(:,1:2)] * ([1 -1]./norm([1 -1]))';
y = [colors(:,3)];
theta = atan(y./x);
pred_expnt = nan(size(phi));
for a = 1:numel(phi)
    mecVec = (cos(phi(a)).*[1/sqrt(2); -1/sqrt(2);0])+(sin(phi(a)).*[0;0;1]);
    R = abs(colors * mecVec);
    [predPrefTheta, pred_expnt(a), predGain] = raisedCos(theta, R);
end

figure,
set(gcf, 'position', [9         404        1385         380])
subplot(1,2,1), hold on
if PLOTGTPOP
    tmp = [raisedExp(:); raisedExpGT(:)];
    edges = logspace(log10(min(tmp)*.9), log10(max(tmp)*1.1), 35);
    %edges = logspace(log10(10e-3), log10(10e2), 35); %for sfn poster
    GTcounts = histc(raisedExpGT, edges);
    DTcounts = histc(raisedExp, edges);
    maxY = max([DTcounts(:);GTcounts(:)]);
else
    tmp = raisedExp(:);
    edges = logspace(log10(min(tmp)*.9), log10(max(tmp)*1.1), 35);
    DTcounts = histc(raisedExp, edges);
    maxY = max(DTcounts);
end
if PLOTGTPOP
    h_GT = bar(edges, GTcounts, 'type', 'histc');
    set(h_GT, 'facecolor', 'k')
end
edges2 = logspace(log10(min(pred_expnt)*.9), log10(max(pred_expnt)*1.1),200);
predExpntPDF = histc(pred_expnt, edges2)./numel(pred_expnt);
h_DT = bar(edges, DTcounts, 'type', 'histc');
set(h_DT, 'facecolor', [0.4 0.4 0.4])
[h_ax, ~, h_Pred] = plotyy([nan], [nan], edges2, predExpntPDF, @plot, @plot);
set(h_Pred, 'color', 'r', 'linewidth', 2)
set(h_ax(2), 'ycolor', 'r')
set([gca, h_ax], 'xscale', 'log', 'xlim', [edges(1)./10, edges(end)*10])
%set([gca, h_ax], 'tickdir', 'out', 'xscale', 'log', 'xlim', [10^(-2.5),10^(2)]) %for sfn poster
set(gca, 'ytick', [0:10:maxY], 'ylim', [0, maxY*1.02])
xlabel('exponent')
hold off

subplot(1,2,2), hold on,
filtTRs = rawTRs(l_validConds,:);
filtMinTRs = min(filtTRs, [], 2);
l_nanTRs = isnan(filtMinTRs);
valForNan = 1.5*max(filtMinTRs);
filtMinTRs(l_nanTRs) = valForNan;
plot(raisedExp, filtMinTRs, 'bo')
plot(raisedExp(l_nanTRs), filtMinTRs(l_nanTRs), 'ko')
set(gca, 'yscale', 'log', 'xscale', 'log')
%xlim([0.03, 15])
ylim([.7 valForNan*1.05])
l_nan_expAndTR = l_nanTRs | isnan(raisedExp);
[rho, p] = corr(raisedExp(~l_nan_expAndTR), filtMinTRs(~l_nan_expAndTR), 'type', 'spearman');
title(sprintf('Rho = %.3f, p = %.3f', rho, p));
xlabel('Exponent')
ylabel('Threshold Ratio')
hold off

% ***************************************** %
%           compare the three values
% ***************************************** %
figure, hold on,
plot3(opRatio, circVarDT, planeStatDT, 'k.')
plot3([0, 1], [0 1], [0 1], 'k-')
xlabel('opRatio')
ylabel('circVarDT')
zlabel('planeStatDT')
hold off

% ***************************************** %
%    Modulation ratio vs. tuning width
% ***************************************** %
filtModRatios = rawModRatios(l_validConds);
figure
subplot(1,2,1), hold on,
if PLOTGTPOP
    plot(rawModRatioGT, rawRaisedExpGT, 'go')
end
plot(filtModRatios, raisedExp, 'bo')
l_nans = isnan(filtModRatios + raisedExp);
[rho, p] = corr(filtModRatios(~l_nans), raisedExp(~l_nans), 'type', 'spearman')
set(gca, 'xscale', 'log', 'yscale', 'log')
title(sprintf('Exponent, r = %.3f, p = %.3f', rho, p))
xlabel('modulation ratio')
ylabel('exponent')
hold off

subplot(1,2,2), hold on,
if PLOTGTPOP
    loglog(rawModRatioGT, rawPlaneStatGT, 'go');
end
loglog(filtModRatios, planeStatDT, 'bo')
l_nans = isnan(filtModRatios + planeStatDT);
[rho, p] = corr(filtModRatios(~l_nans), planeStatDT(~l_nans), 'type', 'spearman')
set(gca, 'xscale', 'log', 'yscale', 'log')
title(sprintf('Plane Stat, r = %.3f, p = %.3f', rho, p))
xlabel('modulation ratio')
ylabel('planestat')
hold off

% ***************************************** %
%    CSI vs. tuning width
% ***************************************** %
filtCSI = rawCSIs(l_validConds);
figure
subplot(1,2,1)
semilogx(filtCSI, planeStatDT, 'bo')
l_nocorr = isnan(filtCSI + planeStatDT) | isinf(filtCSI + planeStatDT);
[rho, p] = corr(filtCSI(~l_nocorr), planeStatDT(~l_nocorr), 'type', 'spearman');
title(sprintf('rho = %.3f p = %.3f', rho, p));
subplot(1,2,2)
loglog(filtCSI, raisedExp, 'bo')
l_nocorr = isnan(filtCSI + raisedExp) | isinf(filtCSI + raisedExp);
[rho, p] = corr(filtCSI(~l_nocorr), planeStatDT(~l_nocorr), 'type', 'spearman');
title(sprintf('rho = %.3f p = %.3f', rho, p));


% ***************************************** %
%    CSI vs. pref color
% ***************************************** %
l_s = ismember(sign(prefIsolum(l_validConds,:)), [0 0 1], 'rows');
l_lvm = ismember(sign(prefIsolum(l_validConds,:)), [1 -1 0], 'rows');
l_swm = ismember(sign(prefIsolum(l_validConds,:)), [1 -1 -1], 'rows');
l_swl = ismember(sign(prefIsolum(l_validConds,:)), [1 -1 1], 'rows');

figure,
subplot(1,2,1), hold on,
plot(ones(sum(l_s),1), planeStatDT(l_s), 'bo')
plot(ones(sum(l_lvm),1)*2, planeStatDT(l_lvm), 'ro')
plot(ones(sum(l_swm),1)*3, planeStatDT(l_swm), 'go')
plot(ones(sum(l_swl),1)*4, planeStatDT(l_swl), 'mo')
xlim([0 5])
subplot(1,2,2), hold on,
semilogy(ones(sum(l_s),1), raisedExp(l_s), 'bo')
semilogy(ones(sum(l_lvm),1)*2, raisedExp(l_lvm), 'ro')
semilogy(ones(sum(l_swm),1)*3, raisedExp(l_swm), 'go')
semilogy(ones(sum(l_swl),1)*4, raisedExp(l_swl), 'mo')
set(gca, 'yscale', 'log')
xlim([0 5])

%doing a kruskal walis
group = (l_s) + (l_lvm*2) + (l_swm*3) + (l_swl*4);
p = kruskalwallis(planeStatDT, group);

%% (10.2) COLOR TUNING WIDTH VS THRESHOLD RATIO (ALL GT COLORS)
clear, clc, close all
global cardVsIntBatchPath blpBatchPath
preprocessDTbatchData

% a new analysis sum(resp) ./ (max(resp).*9)
sScaleFactor = 5; %necessary for finding the orth. Int color dir
l_validConds = ~commonExclusions;
exptList = find(l_validConds);
sphereStat_DT = nan(1, length(exptList));
plotSurf = 0;
for a = 1:length(exptList);
    ex = exptList(a);
    colors = out.dat(ex).grating.color.colors;
    colors(:,3) = colors(:,3)./sScaleFactor;
    norms = sqrt(sum(colors.^2,2));
    unitVecs = bsxfun(@ldivide, colors, norms);
    colResp = out.dat(ex).grating.color.colresp(:,1);
    if max(colResp) == 0;
        continue
    else
        sphereStat_DT(a) = sum(colResp) ./ (numel(colResp) .* max(colResp));
    end
    
    if plotSurf
        figure, hold on,
        cordinates = unitVecs .* repmat(colResp, 1, 3);
        cordinates = [cordinates; cordinates.*-1];
        if any(isnan(cordinates))|any(isinf(cordinates)), continue, end
        sphere = unitVecs .* repmat(max(colResp), length(colResp), 3);
        sphere = [sphere; sphere.*-1];
        tuningTri = DelaunayTri(cordinates);
        tuningTri = convexHull(tuningTri);
        trisurf(tuningTri, cordinates(:,1), cordinates(:,2), cordinates(:,3), 'facecolor', 'c', 'facealpha', 0.5, 'linewidth', 0.5)
        sphereTri = DelaunayTri(sphere);
        sphereTri = convexHull(sphereTri);
        trisurf(sphereTri, sphere(:,1), sphere(:,2), sphere(:,3), 'facecolor', 'k', 'facealpha', 0.2, 'linestyle', ':')
        plot3([3,0,0;-3,0,0], [0,3,0;0,-3,0], [0,0,3;0,0,-3], 'r-', 'linewidth', 3)
        xlabel('L-cone contrast')
        ylabel('M-cone contrast')
        zlabel('S-cone contrast')
        axis equal square
        title(sprintf('sphereStat = %.3f', sphereStat_DT(a)))
    end
end


%load in a mass of GT data and find the sphere stat for them as well
gt = load('/Users/charliehass/LabStuff/NexFiles/Charlie/Batch Data And Text Files/allGTFiles.mat');
gt = gt.out;
sphereStat_GT = nan(numel(gt.dat), 1);
for a = 1:numel(gt.dat)
    try
        nTrials = gt.dat(a).gratings.color.colresp(:,3);
        if min(nTrials)<6
            continue
        end
        R = gt.dat(a).gratings.color.colresp(:,1);
        sphereStat_GT(a) = sum(R)./(numel(R).*max(R));
    catch
        continue
    end
end
sphereStat_GT(isnan(sphereStat_GT)) = [];



%determine what the distribution of sphereStats would be assuming a linear
%model
colors = bsxfun(@times, out.dat(50).grating.color.colors, [1, 1, 1./sScaleFactor]);
nThetas = 1000;
azimuth = linspace(0, pi, nThetas);
elevation = linspace(0, pi, nThetas);
azimuth(end) = [];
elevation(end) = [];
comb = fullfact([nThetas-1, nThetas-1]);
pred_sphereStat = nan(size(comb,1), 1);
for a = 1:size(comb,1)
    theta = azimuth(comb(a,1));
    phi = elevation(comb(a,2));
    unitVec = [sin(theta)*cos(phi); sin(theta)*sin(phi); cos(theta)];
    R = abs(colors * unitVec);
    pred_sphereStat(a) = sum(R) ./ (max(R).*numel(R));
end


%plot the results
figure,
set(gcf, 'position', [9         404        1385         380])
subplot(1,2,1), hold on,
edges1 = linspace(0,1, 30);
sphereStatCounts_DT = histc(sphereStat_DT, edges1);
sphereStatCounts_GT = histc(sphereStat_GT, edges1);
edges2 = linspace(0,1, 200);
predStatPDF = histc(pred_sphereStat, edges2)./numel(pred_sphereStat);
bar(edges1, sphereStatCounts_GT, 'FaceColor', 'k');
[ax, h1, h2] = plotyy(edges1, sphereStatCounts_DT, edges2, predStatPDF, @bar, @plot);
set(ax, 'xLim', [0,1]);
%set(ax(1), 'ylim', [0,13], 'ytick', [0:2:12])
set(h2, 'color', 'r', 'linewidth', 2)
set(get(h1, 'children'), 'faceColor', [.40, .40, .40])
set(ax(2), 'ycolor', 'r')
set(get(ax(1), 'ylabel'), 'string', 'Raw Data', 'fontsize', 14, 'fontname', 'times')
set(get(ax(2), 'ylabel'), 'string', 'Linear Prediction', 'fontsize', 14, 'fontname', 'times')
xlabel('Sphere Statistic')
title(sprintf('Ngt = %d, Ndt = %d', numel(sphereStat_GT), numel(sphereStat_DT)));
hold off

subplot(1,2,2), hold on,
filtTRs = rawTRs(l_validConds,:);
filtMinTRs = min(filtTRs, [], 2);
l_nans = isnan(filtMinTRs);
valForNan = 1.5*max(filtMinTRs);
filtMinTRs(l_nans) = valForNan;
plot(sphereStat_DT, filtMinTRs, 'o', 'markerfacecolor', 'b')
plot(sphereStat_DT(l_nans), filtMinTRs(l_nans), 'ko', 'markerfacecolor', 'k')
set(gca, 'yscale', 'log')
ylim([.98*min(filtMinTRs), 1.1*valForNan])
xlim([min(sphereStat_DT)*0.9, max(sphereStat_DT)*1.1])
[rho, p] = corr(sphereStat_DT(~l_nans)', filtMinTRs(~l_nans), 'type', 'spearman');
xlabel('Spherical Statistic')
ylabel('Min TR')
title(sprintf('Rho = %.3f, p = %.3f', rho, p));


%% (11.1) VARIABLILITY AS A FUNCTION OF COLOR DIRECTION
%
% HERE I'LL ASK THE QUESTION: DOES VARIABILITY GROW WITH FIRING RATE
% SIMILARLY FOR L-M AND S?
%
%****************************************************************%



%cells with strange mean/var relationships:
%K110510002
%K102010002
%K101210004
%K090210005
%K071110004
%K070110002
%S031011005
%S081709004
%S051409008b
%S042109007



clear, clc
global cardVsIntBatchPath blpBatchPath
preprocessDTbatchData

l_valid = ~commonExclusions;
validExpts = find(l_valid);
makePlot = 1;
[l_s, l_lvm] = deal(false(length(validExpts),1));
[bwls_raw, bols_raw, btls_raw] = deal(nan(length(validExpts),1));
[bwls_log, bols_log] = deal(nan(length(validExpts),2));
sigmaAtThresh = nan(length(validExpts),5);
weightedLeastSquares = @(A, B, WT) ((A'*WT*A)\A'*WT*B); %function definition for WLS
for a = 1:length(validExpts);
    ex = validExpts(a);
    exptColors = sign(out.dat(ex).expt.standColors);
    l_s(a) = ismember([0 0 1], exptColors, 'rows');
    l_lvm(a) = ismember([1 -1 0], exptColors, 'rows');
    if l_s(a) && l_lvm(a)
        error('something is wrong')
    elseif l_s(a)
        datIdx = ismember(exptColors, [0 0 1], 'rows');
    elseif l_lvm(a)
        datIdx = ismember(exptColors, [1 -1 0], 'rows');
    else
        fprintf('cell <%d> has neither card color\n', a)
    end
    
    %bail on this data set if the neurometric threshold is a nan (and thus,
    %the neuron's CRF was flat.
    if isnan(out.dat(ex).c.alpha(datIdx))
        l_s(a) = false;
        l_lvm(a) = false;
        continue
    end
    
    
    %plot variance as a function of mean counts, and fit the points with a
    %line (assuming a proportional relationship between variance and rate.
    TRIALLENGTH = 0.666;
    crf = out.dat(ex).c.crfIn{datIdx}';
    nTrials = cellfun(@length, crf);
    tlengths = mat2cell(repmat(TRIALLENGTH,size(crf)), ones(size(crf)));
    countsPerTrial = cellfun(@(x,y) round(x.*y), crf, tlengths, 'uniformoutput', 0); %round so that you get integer number of spikes
    meanCountsByContrast = cellfun(@mean, countsPerTrial);
    varsByContrast = cellfun(@var, countsPerTrial);
    
    
    %First: find the weights for the weighted regression
    wt_raw = [];
    wt_log = [];
    for w = 1:length(meanCountsByContrast); %i.e., the number of contrasts
        inds = unidrnd(nTrials(w), nTrials(w), 1e3);
        btSamps = countsPerTrial{w}(inds);
        wt_raw(w) = var(var(btSamps));
        logRaw = log10(var(btSamps));
        wt_log(w) = var(logRaw(~isinf(logRaw))); %this line is a hack, and prevents zero variance cases from occuring...
    end
    
    %Second: compute the WLS on the raw data
    X = meanCountsByContrast(:);
    Y = varsByContrast(:);
    l_zeroMean = meanCountsByContrast(:)==0;
    X(l_zeroMean) = [];
    Y(l_zeroMean) = [];
    wt_raw(l_zeroMean) = [];
    W = diag(1./wt_raw,0);
    bwls_raw(a) = weightedLeastSquares(X,Y,W);
    bols_raw(a) = X\Y;
    
    %Third: compute the WLS on the log data
    X = [log10(X), ones(length(X),1)];
    Y = log10(Y);
    wt_log(l_zeroMean) = [];
    W = diag(1./wt_log,0);
    bwls_log(a,:) = weightedLeastSquares(X,Y,W);
    bols_log(a,:) = X\Y;
    
    %Fourth: TLS via PCA
    X = meanCountsByContrast(:);
    Y = varsByContrast(:);
    l_zeroMean = meanCountsByContrast(:)==0;
    X(l_zeroMean) = [];
    Y(l_zeroMean) = [];
    yx = [(Y(:)), (X(:))];
    covMtx = (yx' * yx) ./ (size(yx,1)-1);
    [vec, val] = eig(covMtx);
    PC1 = vec(:,2);
    btls_raw(a) = PC1(1)./PC1(2);
    
    
    if makePlot && any(sigmaAtThresh(a,:)>10)
        if l_s(a)
            colors = 'bo';
        elseif l_lvm(a)
            colors = 'ro';
        end
        
        
        figure
        set(gcf, 'position', [157   102   923   608], 'name', out.fnames{ex}{1})
        subplot(2,2,1), hold on, %the raw data
        x = linspace(0, max(meanCountsByContrast),50);
        plot(meanCountsByContrast,varsByContrast,colors)
        plot(x, bols_raw(a).*x, 'k')
        plot(x, bwls_raw(a).*x, 'c')
        plot(x, btls_raw(a).*x, 'g')
        legend('raw data', 'OLS', 'WLS', 'TLS', 'location', 'northwest')
        xylim = [min([meanCountsByContrast(:);varsByContrast(:)])*.9, max([meanCountsByContrast(:);varsByContrast(:)])*1.2];
        title(sprintf('slope = %.3f', bwls_raw(a)))
        
        subplot(2,2,2), hold on, %log data
        plot(meanCountsByContrast(meanCountsByContrast>0),varsByContrast(meanCountsByContrast>0),colors)
        yols_log = 10^bols_log(a,2).*x.^bols_log(a,1);
        plot(x, yols_log, 'k')
        ywls_log = 10^bwls_log(a,2).*x.^bwls_log(a,1);
        plot(x, ywls_log, 'c')
        title(sprintf('exponent = %.3f, slope = %.3f', bwls_log(a,1), 10^bwls_log(a,2)))
        
        subplot(2,2,3), hold on, %the neurometric function
        roc = out.dat(ex).c.roc{datIdx};
        norms = out.dat(ex).expt.norms{datIdx};
        cc = [norms(2)*0.98 : 0.001 : norms(end)*1.04];
        alpha = out.dat(ex).c.alpha(datIdx);
        beta = out.dat(ex).c.beta(datIdx);
        mod = 1 - 0.5.*exp(-(cc./alpha).^beta);
        plot(norms, roc, colors)
        plot(cc, mod, 'k')
        set(gca, 'xscale', 'log')
        ylim([0.4 1.03]);
        xlim([cc(1), cc(end)])
        title(sprintf('alpha = %.3f', alpha));
        
        subplot(2,2,4), hold on, %contrast response function
        sem = sqrt(varsByContrast)./sqrt(nTrials);
        errorbar(norms, meanCountsByContrast, sem, [colors,'-'])
        xlim([0, norms(end)*1.02])
    end
end

%For now, just look at the WLS from log transformed data (i.e., assuming a
%power law relationship)
figure
exponent = bwls_log(:,1);
slope = 10.^bwls_log(:,2);
l_nans = isnan(exponent);
l_negExp = exponent<0;
s_slope = slope(l_s & ~l_nans & ~l_negExp);
lvm_slope = slope(l_lvm & ~l_nans & ~l_negExp);
s_exp = exponent(l_s & ~l_nans & ~l_negExp);
lvm_exp = exponent(l_lvm & ~l_nans & ~l_negExp);
subplot(1,2,1), hold on, %exponent
plot(ones(length(lvm_exp),1), lvm_exp, 'ro')
plot(ones(length(s_exp),1)*3, s_exp, 'bo')
xlim([0,4])
[p,h] = ranksum(s_exp, lvm_exp);
title(sprintf('wilcoxon p = %.3f', p));
ylabel('Exponent')
subplot(1,2,2), hold on, %slope
plot(ones(length(lvm_slope),1), lvm_slope, 'ro')
plot(ones(length(s_slope),1)*3, s_slope, 'bo')
xlim([0,4])
[p,h] = ranksum(s_slope, lvm_slope);
title(sprintf('wilcoxon p = %.3f', p));
ylabel('Slope')

figure, hold on,
plot(s_slope, s_exp, 'bo')
plot(lvm_slope, lvm_exp, 'ro')
xlabel('Slope')
ylabel('Exponent')

%now look at the WLS data from non-log transformed data (i.e., assuming a
%linear relationship)
figure, hold on,
slope = bwls_raw;
l_nans = isnan(slope);
s_slope = slope(l_s & ~l_nans);
lvm_slope = slope(l_lvm & ~l_nans);
plot(ones(length(lvm_slope),1), lvm_slope, 'ro')
plot(ones(length(s_slope),1)*3, s_slope, 'bo')
xlim([0,4])
[p,h] = ranksum(s_slope, lvm_slope);
title(sprintf('wilcoxon p = %.3f', p));
ylabel('Slope')

%now look at sigmaAtThresh for S and L-M neurons
nonan = ~isnan(sigmaAtThresh(:,1));
figure, hold on,
%ols
plot([1], sigmaAtThresh(l_lvm,1), 'r.')
plot([2], sigmaAtThresh(l_s,1), 'b.')
ranksum(sigmaAtThresh(l_lvm&nonan,1), sigmaAtThresh(l_s&nonan,1))

%ols_log
plot([5], sigmaAtThresh(l_lvm,2), 'r.')
plot([6], sigmaAtThresh(l_s,2), 'b.')
ranksum(sigmaAtThresh(l_lvm&nonan,2), sigmaAtThresh(l_s&nonan,2))

%wls
plot([9], sigmaAtThresh(l_lvm,3), 'r.')
plot([10], sigmaAtThresh(l_s,3), 'b.')
ranksum(sigmaAtThresh(l_lvm&nonan,3), sigmaAtThresh(l_s&nonan,3))

%wls_log
plot([13], sigmaAtThresh(l_lvm,4), 'r.')
plot([14], sigmaAtThresh(l_s,4), 'b.')
ranksum(sigmaAtThresh(l_lvm&nonan,4), sigmaAtThresh(l_s&nonan,4))

%tls
plot([17], sigmaAtThresh(l_lvm,5), 'r.')
plot([18], sigmaAtThresh(l_s,5), 'b.')
ranksum(sigmaAtThresh(l_lvm&nonan,5), sigmaAtThresh(l_s&nonan,5))
%% (11.2) VARAINCE VS. MEAN FOR L-M NEURONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, clc
global cardVsIntBatchPath blpBatchPath
preprocessDTbatchData

%function definition:
weightedLeastSquares = @(A, B, WT) ((A'*WT*A)\A'*WT*B); %function definition for WLS

fitMeth = 'wls'; % wls or tls
nShuffles = 1e4;
l_valid = ~commonExclusions;
actualSlopeDiffs = nan(length(out.dat),1);
p_diff = nan(length(out.dat),1);
for a = find(l_valid)'
    lvmIdx = ismember(sign(out.dat(a).expt.standColors), [1 -1 0], 'rows');
    if ~any(lvmIdx); continue; end %only analyze the expts with L-M
    if all(lvmIdx); keyboard; end  % this shouldn't happen.
    %if all(isnan(out.dat(a).c.alpha)); continue; end
    
    %iterate over each color direction and determine the relationship b/w
    %mean and variance
    slope = [];
    combinedCounts = [];
    combinedVars = [];
    combinedWts = [];
    colorIdentity = [];
    TRIALLENGTH = 0.666;
    figure
    set(gcf, 'position', [ -1076 255 1054 281], 'name', out.fnames{a}{2})
    for clr = 1:2;
        crf = out.dat(a).c.crfIn{clr}';
        crf(1) = []; %eliminate the zero contrast condition
        nTrials = cellfun(@length, crf);
        tlengths = repmat({TRIALLENGTH},size(crf));
        countsPerTrial = cellfun(@(x,y) round(x.*y), crf, tlengths, 'uniformoutput', 0); %round so that you get integer number of spikes
        meanCountsByContrast = cellfun(@mean, countsPerTrial);
        varsByContrast = cellfun(@var, countsPerTrial);
        
        
        %calculate the slope via weighted least squares
        wt_raw = [];
        for w = 1:length(meanCountsByContrast); %i.e., the number of contrasts
            inds = unidrnd(nTrials(w), nTrials(w), 1e4);
            btSamps = countsPerTrial{w}(inds);
            wt_raw(w) = var(var(btSamps));
        end
        
        X = meanCountsByContrast(:);
        Y = varsByContrast(:);
        l_zeroMean = meanCountsByContrast(:)==0;
        X(l_zeroMean) = [];
        Y(l_zeroMean) = [];
        switch fitMeth
            case 'wls'
                W = wt_raw(~l_zeroMean);
                W = diag(1./W,0);
                slope(clr) = weightedLeastSquares(X,Y,W);
            case 'tls'
                slope(clr) = tls(X,Y);
        end
        
        %aggregate the data across colors:
        combinedCounts = [combinedCounts ; meanCountsByContrast(:)];
        combinedVars = [combinedVars ; varsByContrast(:)];
        combinedWts = [combinedWts ; wt_raw(:)];
        colorIdentity = [colorIdentity; ones(length(meanCountsByContrast),1)*clr];
        
        %plot the resultant mean to variance relationships
        color = {'k', 'r'};
        subplot(1,3,1), hold on,
        plot(X,Y, 'o', 'color', color{(find(lvmIdx)==clr)+1})
        plot(X, X*slope(clr), '-', 'color', color{(find(lvmIdx)==clr)+1})
        xlabel('counts')
        ylabel('variance')
        hold off
        
        %plot the contrast response functions
        subplot(1,3,2), hold on,
        errorbar(out.dat(a).expt.norms{clr}(2:end), meanCountsByContrast, sqrt(varsByContrast)./sqrt(nTrials), '.-', 'color',color{(find(lvmIdx)==clr)+1})
        xlabel('contrast')
        ylabel('counts')
        set(gca, 'xscale', 'log')
        axis tight
        hold off
    end
    
    %store the empirical difference in slopes
    actualSlopeDiffs(a) = slope(lvmIdx)-slope(~lvmIdx);
    
    %now randomly permute the data sets and build up a distribution of
    %slope differences
    shuffSlopeDiffs = nan(nShuffles, 1);
    nColorOne = sum(colorIdentity==1);
    if nColorOne~=7; keyboard; end %shouldn't happen
    for shuff = 1:nShuffles;
        mapping = randperm(length(combinedCounts));
        l = mapping<=nColorOne;
        
        shuffSlope = [];
        for clr = 1:2
            if clr == 1;
                l_shuffClr = l;
            elseif clr == 2;
                l_shuffClr = ~l;
            end
            
            X = combinedCounts(l_shuffClr);
            Y = combinedVars(l_shuffClr);
            l_zeroMean = X==0;
            X(l_zeroMean) = [];
            Y(l_zeroMean) = [];
            switch fitMeth
                case 'wls'
                    W = combinedWts(l_shuffClr);
                    W(l_zeroMean) = [];
                    W = diag(1./W,0);
                    shuffSlope(clr) = weightedLeastSquares(X,Y,W);
                case 'tls'
                    shuffSlope(clr) = tls(X,Y);
            end
        end
        %store the difference b/w shuffled slopes
        shuffSlopeDiffs(shuff) = shuffSlope(1) - shuffSlope(2);
    end
    
    %plot the resulting distribution of shuffled slope values
    if actualSlopeDiffs(a)>0
        prcnt = sum(shuffSlopeDiffs>actualSlopeDiffs(a)) ./ sum(~isnan(shuffSlopeDiffs));
    elseif actualSlopeDiffs(a)<0
        prcnt = sum(shuffSlopeDiffs<actualSlopeDiffs(a)) ./ sum(~isnan(shuffSlopeDiffs));
    end
    subplot(1,3,3), hold on,
    hist(shuffSlopeDiffs, 30)
    plot(actualSlopeDiffs(a), 10, 'wv', 'markerfacecolor', 'r')
    title(sprintf('p = %.3f', prcnt));
    
end


%% (13) CSI VS. THRESHOLD RATIO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, clc, close all
global cardVsIntBatchPath blpBatchPath
preprocessDTbatchData


l_validConds = ~(commonExclusions | out.errors(:, prefTieInd));
[minTRs, ~] = min(rawTRs, [], 2);
filtMinTRs = minTRs(l_validConds);
l_nanTRs = isnan(filtMinTRs);
filtMinTRs(l_nanTRs) = max(filtMinTRs) .* 1.5;

filtCSIs = rawCSIs(l_validConds);
l_zeroCSIs = filtCSIs==0;
l_infCSIs = isinf(filtCSIs);
filtCSIs(l_infCSIs) = max(filtCSIs(~l_infCSIs)) .* 1.5;
filtCSIs(l_zeroCSIs) = min(filtCSIs(~l_zeroCSIs)) .* 0.95;

filtJSCSIs = rawJSCSIs(l_validConds);
l_zeroJSCSIs = filtJSCSIs==0;
l_infJSCSIs = isinf(filtJSCSIs);
filtJSCSIs(l_infCSIs) = max(filtJSCSIs(~l_infCSIs)) .* 1.5;
filtJSCSIs(l_zeroCSIs) = min(filtJSCSIs(~l_zeroCSIs)) .* 0.95;


figure;
subplot(2,2,1), hold on, %regular deff of CSI
plot(filtCSIs, filtMinTRs, 'bo')
plot(filtCSIs(l_nanTRs), filtMinTRs(l_nanTRs), 'ro')
plot(filtCSIs(l_infCSIs), filtMinTRs(l_infCSIs), 'go')
ylim([min(filtMinTRs)-0.1, max(filtMinTRs)+0.3])
set(gca, 'xscale', 'log', 'yscale', 'log')
ypts = get(gca, 'ylim');
plot([0.5 0.5], ypts, 'k:')
plot([2 2], ypts, 'k:')
xlabel('Color Sensitivity Index')
ylabel('Min Threshold Ratio')
[rho, p] = corr(filtCSIs(:), filtMinTRs(:), 'type', 'spearman');
title(sprintf('rho = %.3f; p = %.3f', rho, p))
text(0.04, 0.9, sprintf('n = %d', length(filtMinTRs)))
hold off

subplot(2,2,3); %regular deff by color
hold on,
filtPrefIsolum = prefIsolum(l_validConds, :);
filtFnames = out.fnames(l_validConds);
colorTypes = [1 -1 0; 0 0 1; 1 -1 1; 1 -1 -1];
markerColors = {'r', 'b', 'm', 'g'};
label1 = xlabel('hello');
for a = 1:length(filtMinTRs)
    idx = ismember(colorTypes, sign(filtPrefIsolum(a,:)), 'rows');
    l = loglog(filtCSIs(a), filtMinTRs(a), 'o');
    set(l, 'MarkerEdgeColor', markerColors{idx}, 'buttondownfcn', ['printFname(', num2str(a), ',label1)'])
end
set(gca, 'xscale', 'log', 'yscale', 'log')
ylim([min(filtMinTRs)-0.2, max(filtMinTRs)+0.6])
ypts = get(gca, 'ylim');
plot([0.5 0.5], ypts, 'k:')
plot([2 2], ypts, 'k:')
hold off,
printFname = @(ind, hand) set(hand, 'string', filtFnames{ind}(1));


subplot(2,2,2) %johnson and shapley CSI
hold on,
plot(filtJSCSIs, filtMinTRs, 'bo');
plot(filtJSCSIs(l_nanTRs), filtMinTRs(l_nanTRs), 'ro')
plot(filtJSCSIs(l_infJSCSIs), filtMinTRs(l_infJSCSIs), 'go')
plot(filtJSCSIs(l_zeroJSCSIs), filtMinTRs(l_zeroJSCSIs), 'mo')
ylim([min(filtMinTRs)-0.2, max(filtMinTRs)+0.6])
set(gca, 'xscale', 'log', 'yscale', 'log')
ypts = get(gca, 'ylim');
plot([0.5 0.5], ypts, 'k:')
plot([2 2], ypts, 'k:')
xlabel('Color Sensitivity Index (Johnson & Shapley)')
ylabel('Min Threshold Ratio')
[rho, p] = corr(filtJSCSIs(:), filtMinTRs(:), 'type', 'spearman');
title(sprintf('rho = %.3f; p = %.3f', rho, p))
hold off


subplot(2,2,4); %j/s deff by color
hold on,
label2 = xlabel('hello2');
for a = 1:length(filtMinTRs)
    idx = ismember(colorTypes, sign(filtPrefIsolum(a,:)), 'rows');
    l = loglog(filtJSCSIs(a), filtMinTRs(a), 'o');
    set(l, 'MarkerEdgeColor', markerColors{idx}, 'buttondownfcn', ['printFname(', num2str(a), ',label2)'])
end
set(gca, 'xscale', 'log', 'yscale', 'log')
ylim([min(filtMinTRs)-0.2, max(filtMinTRs)+0.6])
ypts = get(gca, 'ylim');
plot([0.5 0.5], ypts, 'k:')
plot([2 2], ypts, 'k:')
hold off,


bins = logspace(log10(0.04), log10(35), 25);
N = histc(filtCSIs, bins);
figure, 
bar(bins, N, 'histc')
set(gca, 'xscale', 'log', 'tickDir', 'out')
box off
xlim([0.04 40])



%% (23) NT, PT, and TR VS. SPATIAL FREQUENCY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, clc, close all
global cardVsIntBatchPath blpBatchPath
preprocessDTbatchData

x = log10(linspace(0.45,4.3,100));
yticloc = {log10([0.0025, .005, 0.01, 0.02, 0.04, 0.08, .16, .32]);...
    log10([0.0025, .005, 0.01, 0.02, 0.04, 0.08, .16, .32]);...
    log10([0.5, 1, 2, 4, 8])};
xticloc = log10([0.5, 1, 2, 4]);
pred = [log10([rawSF;rawSF]), ones(2*numel(rawSF),1)];
sfs = log10(rawSF');

figure
for a = 1:3
    switch a
        case 1
            tmp = rawPTs;
        case 2
            tmp = rawNTs;
        case 3
            tmp = rawTRs;
    end
    
    %run the regression
    resp = log10(tmp(:)); %tmp is <nx2>, convert to a col vec
    l_nan = isnan(resp) | isnan(pred(:,1));
    b = pred(~l_nan,:) \ resp(~l_nan);
    
    
    %plot the results
    subplot(3,1,a), hold on,
    plot(sfs, log10(tmp'), '.');
    plot(x, b(1).*x+b(2), 'k-');
    xlim([x(1), x(end)]);
    set(gca, 'xtick', xticloc, 'xticklabel', cellfun(@(x) num2str(10.^(x)), mat2cell(xticloc(:), ones(size(xticloc(:)))), 'uniformoutput', 0))
    set(gca, 'ytick', yticloc{a}, 'ytickLabel', cellfun(@(x) num2str(10.^(x)), mat2cell(yticloc{a}', ones(size(yticloc{a}'))), 'uniformoutput', 0))
    [rho, p] = corr(pred(~l_nan), resp(~l_nan), 'type', 'spearman')
    title(sprintf('rho = %.3f, p = %.3f', rho, p))
end


%% (38) COLOR CODING CSI VS. TR BY SPATIAL FREQ
clear, clc, close all
global cardVsIntBatchPath
preprocessDTbatchData

%a one line function
dogPred = @(R0, SF, k_c, mu_c, sigma_c, k_s, mu_s, sigma_s) (R0 + (k_c .* exp(-((SF-mu_c)./(2.*sigma_c)).^2)) - (k_s .* exp(-((SF-mu_s)./(2.*sigma_s)).^2)));

PLOT = 1;
l_validConds = ~(commonExclusions);
lowpass = false(length(out.dat),1);
l_isolum = false(length(out.dat),1);
[allSfs, allResp] = deal(nan(numel(l_validConds), 4));
[prefSF_dog, bandwidth] = deal(nan(numel(out.dat),1));
for a = 1:length(out.dat)
    if ~l_validConds(a)
        continue
    end
    
    GT = gtobj(out.fnames{a}{2});
    t_on = GT.trial(:, GT.idx.stimon);
    nFrames = GT.trial(:, GT.idx.nFrames);
    t_off = t_on + (nFrames ./ GT.sum.exptParams.framerate);
    t_fpon = GT.trial(:, GT.idx.fpon);
    
    %use cell arrays to get the firing rates
    t_fpon = mat2cell(t_fpon, ones(size(t_fpon)));
    t_on = mat2cell(t_on, ones(size(t_on)));
    t_off = mat2cell(t_off, ones(size(t_off)));
    rates = cellfun(@(b, e, sp) sum((sp>=b)&(sp<e))./(e-b), t_on, t_off, GT.ras(:, GT.idx.spikes));
    bkgndRate = cellfun(@(b,e,sp) sum((sp>=b)&(sp<=e))./(e-b), t_fpon, t_on, GT.ras(:, GT.idx.spikes));
    
    
    %now find the pref sf
    l_sf = GT.trial(:, GT.idx.protocol) == 2;
    if sum(l_sf) < 3;
        disp('Not enough SF trials!!!')
        continue
    end
    
    sptFreqs = unique(GT.trial(l_sf,GT.idx.sf));
    meanRates = [];
    SEM = [];
    for i = 1:length(sptFreqs)
        tList = l_sf & (GT.trial(:,GT.idx.sf)==sptFreqs(i));
        if sum(tList)==0; keyboard; end
        meanRates(i) = mean(rates(tList));
        SEM(i) = std(rates(tList))./sqrt(sum(tList));
    end
    
    if all(meanRates < 3)
        disp('Firing rates too low!!!')
        continue
    end
    
    
    %compare the lowest and pref sf
    [~, maxRateInd] =  max(meanRates);
    [~, smallStimInd] = min(sptFreqs);
    if meanRates(smallStimInd) >= ((meanRates(maxRateInd))./2);
        lowpass(a) = true;
    end
    
    
    %determine bandwidth using the fitted SF tuning functions
    R0 = mean(bkgndRate);
    [Kc, SIGMAc, Ks, SIGMAs] = DoGfit(R0, GT.trial(l_sf, GT.idx.sf), rates(l_sf));
    SFpred = 0.5:0.001:4;
    Rpred = dogPred(R0, SFpred, Kc, 0, SIGMAc, Ks, 0, SIGMAs);
    [maxRate, maxIdx] = max(Rpred);
    prefSF_dog(a) = SFpred(maxIdx);
    err = abs(Rpred-(maxRate./2));
    [~,lowSideIdx] = min(err(1:maxIdx));
    [~,highSideIdx] = min(err(maxIdx:end));
    highSideIdx = highSideIdx+(maxIdx-1);
    
    if lowSideIdx == 1; %low pass
        bandwidth(a) = inf;
    elseif highSideIdx == numel(Rpred); %highpass
        bandwidth(a) = -1;
    elseif (lowSideIdx == 1) && (highSideIdx == numel(Rpred)) %very broad
        bandwith(a) = nan;
    else
        bandwidth(a) = (log2(SFpred(highSideIdx)) - log2(SFpred(lowSideIdx)))./2;
    end
    
    
    %flag the instances where SF was measured with isoluminant stimuli
    LMS = GT.trial(:, [GT.idx.lcc | GT.idx.mcc | GT.idx.scc]);
    if all((LMS(l_sf,:)*[1;1;0])==0)
        l_isolum(a) = true;
        sf_color = unique(LMS(l_sf,:), 'rows');
    end
    
    
    plotList = {'S060409005', 'S032911002', 'S022111004'};
    if PLOT && any(strcmpi(out.fnames{a}{1}, plotList));
        figure
        hold on,
        plot(GT.trial(l_sf, GT.idx.sf), rates(l_sf), 'b.');
        errorbar(sptFreqs, meanRates, SEM, 'ko', 'markerfacecolor', 'k');
        plot(SFpred, Rpred, 'r-')
        plot(rawSF(a), dogPred(R0, rawSF(a), Kc, 0, SIGMAc, Ks, 0, SIGMAs), 'g*');
        plot(prefSF_dog(a), dogPred(R0, prefSF_dog(a), Kc, 0, SIGMAc, Ks, 0, SIGMAs), 'k*')
        if ~(isempty(lowSideIdx) && isempty(highSideIdx))
        plot([SFpred(lowSideIdx), SFpred(highSideIdx)], [(maxRate./2),(maxRate./2)], 'm--')
        end
        [maxPred, ind] = max(Rpred);
        maxSF = SFpred(ind);
        plot(maxSF, maxPred, 'rx')
        set(gca, 'xscale', 'log')
        title(sprintf('%s TR: %.3f CSI: %.3f', out.fnames{a}{1}, min(rawTRs(a,:), [],2), rawCSIs(a)));
        text(0.5, 2, num2str([Kc, SIGMAc, Ks, SIGMAs]))
        hold off,
    end
    allSfs(a,:) = sptFreqs(:)';
    allResp(a,:) = meanRates(:)';
end


%
% Comparing CSI for lowpass and non-lowpass cells
%%%%%%%%%%%%%%%%%%%%%%%
l_validConds = ~commonExclusions;
tmpCSI = rawCSIs;
valForInf = max(tmpCSI(~isinf(tmpCSI))) .* 1.1;
valForZero = min(tmpCSI(tmpCSI>0)) .* .9;
tmpCSI(isinf(tmpCSI)) = valForInf;
tmpCSI(tmpCSI == 0) = valForZero;
csi_lowpass = rawCSIs(l_validConds & lowpass);
csi_others = rawCSIs(l_validConds & ~lowpass);
edges = logspace(log10(min(tmpCSI).*.9), log10(max(tmpCSI).*1.05), 20);
N_lowpass = histc(csi_lowpass, edges);
N_others = histc(csi_others, edges);
figure, hold on,
h_low = bar(edges, N_lowpass, 'type', 'histc');
h_oth = bar(edges, N_others, 'type', 'histc');
set(h_low, 'facecolor', 'k')
set(h_oth, 'facecolor', 'b', 'facealpha', 0.4)
set(gca, 'xscale', 'log')
legend('Lowpass', '','Others')
[p,h] = ranksum(csi_lowpass, csi_others);
xlabel('CSI')
ylabel('Count')



%
% Looking at the population SF tuning
% for color sensitive vs. insensitive neurons
%%%%%%%%%%%%%%%%%%%%%%
l_insensitive = l_validConds & out.errors(:, neuroThresh1Ind) & out.errors(:, neuroThresh2Ind);
l_sensitive = l_validConds & ~(out.errors(:, neuroThresh1Ind) & out.errors(:, neuroThresh2Ind));

figure
subplot(1,3,1) %first the insensitive neurons
tmpSfs = allSfs(l_insensitive,:);
tmpResp = allResp(l_insensitive,:);
l_nan = isnan(tmpSfs(:,1));
tmpResp(l_nan,:) = [];
tmpSfs(l_nan,:) = [];
[Kc, SIGMAc, Ks, SIGMAs] = DoGfit(0, tmpSfs(:), tmpResp(:));
SFpred = logspace(log10(0.4), log10(5), 100);
Rpred = dogPred(0, SFpred, Kc, 0, SIGMAc, Ks, 0, SIGMAs);
hold on,
%plot(allSfs(l_insensitive,:)', allResp(l_insensitive,:)', '-.', 'color', [1,1,1].*0.5)
errorbar(allSfs(3,:), nanmean(allResp(l_insensitive,:)), nanstd(allResp(l_insensitive,:))./sqrt(sum(l_insensitive)), 'bo', 'markerfacecolor', 'b');
plot(SFpred, Rpred, 'b', 'linewidth', 2)
set(gca, 'xscale', 'log')
ylabel('Firing Rate (sp/sec')
xlabel('Spatial Frequency')
title('Insensitive Neurons')

subplot(1,3,2) %now for the sensitive neurons
hold on,
tmpSfs = allSfs(l_sensitive,:);
tmpResp = allResp(l_sensitive,:);
l_nan = isnan(tmpSfs(:,1));
tmpResp(l_nan,:) = [];
tmpSfs(l_nan,:) = [];
[Kc, SIGMAc, Ks, SIGMAs] = DoGfit(0, tmpSfs(:), tmpResp(:));
SFpred = logspace(log10(0.4), log10(5), 100);
Rpred = dogPred(0, SFpred, Kc, 0, SIGMAc, Ks, 0, SIGMAs);
%plot(allSfs(l_sensitive,:)', allResp(l_sensitive,:)', '-.', 'color', [1,1,1].*0.5)
errorbar(allSfs(3,:), nanmean(allResp(l_sensitive,:)), nanstd(allResp(l_sensitive,:))./sqrt(sum(l_sensitive)), 'ko', 'markerfacecolor', 'k');
plot(SFpred, Rpred, 'k', 'linewidth', 2)
xlim([0.35 4.2])
set(gca, 'xscale', 'log')
title('Sensitive')

subplot(1,3,3), hold on, %normalized to the max, sensitive and insensitive plotted on the same axes
%sensitive neurons
l_nan = isnan(allSfs(:,1));
tmpSfs = allSfs(l_sensitive&~l_nan,:);
tmpResp = allResp(l_sensitive&~l_nan,:);
tmpResp = bsxfun(@rdivide, tmpResp, max(tmpResp, [], 2));
[Kc, SIGMAc, Ks, SIGMAs] = DoGfit(0, tmpSfs(:), tmpResp(:));
Rpred = dogPred(0, SFpred, Kc, 0, SIGMAc, Ks, 0, SIGMAs);
meanNorm = mean(tmpResp, 1);
errorbar(tmpSfs(1,:), meanNorm, std(tmpResp,1)./sqrt(size(tmpResp,1)), 'ko', 'markersize', 10)
plot(SFpred, Rpred, 'k', 'linewidth', 2)
% insensitive neurons
tmpSfs = allSfs(l_insensitive&~l_nan,:);
tmpResp = allResp(l_insensitive&~l_nan,:);
tmpResp = bsxfun(@rdivide, tmpResp, max(tmpResp, [], 2));
[Kc, SIGMAc, Ks, SIGMAs] = DoGfit(0, tmpSfs(:), tmpResp(:));
Rpred = dogPred(0, SFpred, Kc, 0, SIGMAc, Ks, 0, SIGMAs);
meanNorm = mean(tmpResp, 1);
errorbar(tmpSfs(1,:), meanNorm, std(tmpResp,1)./sqrt(size(tmpResp,1)), 'bo', 'markersize', 10)
plot(SFpred, Rpred, 'b', 'linewidth', 2)
xlim([0.35 5.2])
set(gca, 'xscale', 'log')
ylim([0.1 0.9])
legend('', 'Sensitive', '', 'Insensitive')

%
% Threshold ratios vs. bandwidth
%%%%%%%%%%%%%%%%%%%%%%
minTRs = min(rawTRs,[],2);
minTRs(isnan(minTRs)) = max(minTRs) .* 1.05;
tmp = bandwidth;
valForInf = max(tmp(~isinf(tmp))) .* 1.1;
tmp(isinf(tmp)) = valForInf;
figure
hold on
plot(minTRs(l_validConds), tmp(l_validConds), 'ko', 'markerfacecolor', 'k');
plot(minTRs(l_validConds&l_isolum), tmp(l_validConds & l_isolum), 'ro', 'markerfacecolor', 'r')
plot([min(minTRs(l_validConds))*.95, max(minTRs(l_validConds)).*1.05], [0,0], ':r')
plot([min(minTRs(l_validConds))*.95, max(minTRs(l_validConds)).*1.05], [valForInf, valForInf], 'r:')
xlim([min(minTRs(l_validConds))*.95, max(minTRs(l_validConds)).*1.05])
ylim([min(tmp)*1.05, max(tmp).*1.05])
set(gca, 'xscale', 'log')
[rho, p] = corr(tmp(tmp>0), minTRs(tmp>0), 'type', 'spearman');
title(sprintf('rho: %.3f p: %.3f', rho, p))
xlabel('Threshold Ratio')
ylabel('Bandwidth')


%% (40) MULTI-DEMINSIONAL CORRELATION AND PARTIAL CORRELATION OF TUNING PARAMETERS

clear, clc, close all
global cardVsIntBatchPath
preprocessDTbatchData

%plan:
%make a scatter plot of exponents, SF, CSI, TR, DS, and orientation bandwidth
%calculate a covariance matrix (do my data meet the assumptions for covariance?)
%partial correlations?
%PCA, then project the raw data onto the first two PCs, or 3d scatter (PC1, PC2, TR)

% GLOBALS
l_valid = ~commonExclusions;
sScaleFactor = 5;
includeInfAndNan = 1; %1 => changes the values to numbers, 0 => leaves them as is (corr probably won't work)
PLOT = 0;

% calculate the exponents and orientation tuning bandwidth
exponents = nan(size(l_valid));
bandwidths_orient = nan(size(l_valid));
bandwidths_sf = []; %going to be a structure b/c I'll estimate the bandwidth two ways...
RFsize = nan(size(l_valid));
for a = find(l_valid)';
    fprintf('Now analyzing cell <%d>\n', a)
    
    
    %
    %RF size
    %%%%
    RFsize(a) = out.dat(a).grating.areasummation.prefsize;
    
    %
    %exponents
    %%%%
    colors = out.dat(a).grating.color.colors;
    colors(:,3) = colors(:,3)./sScaleFactor;
    l_isolum = ~(colors * [1;1;0]);
    isolumColors = colors(l_isolum, :);
    clrResp = out.dat(a).grating.color.colresp(l_isolum,1);
    x = [isolumColors(:,1:2)] * ([1 -1]./norm([1 -1]))';
    y = [isolumColors(:,3)];
    theta = atan(y./x);
    [~, exponents(a), ~] = raisedCos(theta, clrResp);
    
    %
    %orientation tuning bandwidth (half width at half height)
    %%%%
    orients = out.dat(a).grating.orient.stim;
    oriResp = out.dat(a).grating.orient.resp(:,1);
    twoPiIdx = orients == (2*pi); %getting rid of the redundant data point
    if sum(twoPiIdx) == 0; error('could not find the 2pi value'); end
    orients(twoPiIdx) = [];
    oriResp(twoPiIdx) = [];
    
    %center the tuning fxn at pi
    [~, maxIdx] = max(oriResp);
    maxOrient = orients(maxIdx);
    newOrients = orients - (maxOrient-pi);
    newOrients = mod(newOrients, 2*pi);
    [~, sortInd] = sort(newOrients);
    newOrients = newOrients(sortInd);
    newOriResp = oriResp(sortInd);
    
    %find the max, and then interpolate the response at half height (where
    %half height is resp/sqrt(2)
    [maxResp, maxIdx] = max(newOriResp); %the index will have shifted
    if maxResp <= 0; error('max orient rate <= zero'); end
    halfHeight = maxResp./sqrt(2);
    highSideIdx = find(newOriResp(maxIdx+1:end)<halfHeight, 1, 'first') + maxIdx; %the index to the first point to the right that's below half max
    lowSideIdx = find(newOriResp(1:maxIdx)<halfHeight, 1, 'last'); %the index to the first point to the left that's below the halfheight
    if isempty(highSideIdx) || isempty(lowSideIdx)
        bandwidths_orient(a) = inf;
    else
        highBeta = [newOrients(highSideIdx-1), 1; newOrients(highSideIdx), 1] \ newOriResp(highSideIdx-1:highSideIdx); %slope and intercept for linear interpolation
        lowBeta = [newOrients(lowSideIdx), 1; newOrients(lowSideIdx+1), 1] \ newOriResp(lowSideIdx:lowSideIdx+1);
        highOrient = (halfHeight-highBeta(2))./highBeta(1);
        lowOrient = (halfHeight-lowBeta(2))./lowBeta(1);
        bandwidths_orient(a) = (highOrient-lowOrient)./2;
    end
    
    
    %
    %spatial frequency tuning bandwidth
    % (half width at half height)
    % For the linear interpolation, a value of -1 implies low-pass, and a
    % value of inf implies unknown but consistent with high-pass
    %%%%
    if isfield(out.dat(a).grating.sf, 'stim')
        sf = out.dat(a).grating.sf.stim;
        sfResp = out.dat(a).grating.sf.resp(:,1);
        
        %linear interpolation first:
        [maxResp_sf, maxInd_sf] = max(sfResp);
        if maxResp_sf<=0; error('max SF rate <= zero'); end
        halfHeight_sf = maxResp_sf ./ sqrt(2);
        highSideIdx_sf = find(sfResp(maxInd_sf+1:end)<halfHeight_sf, 1, 'first') + maxInd_sf; %the index to the first point to the right that's below half max
        lowSideIdx_sf = find(sfResp(1:maxInd_sf)<halfHeight_sf, 1, 'last'); %the index to the first point to the left that's below the halfheight
        if isempty(highSideIdx_sf) && isempty(lowSideIdx_sf)
            bandwidths_sf(a).lin = nan;
            fprintf('#### %d\n', a);
        elseif isempty(lowSideIdx_sf)
            bandwidths_sf(a).lin = inf;
        elseif isempty(highSideIdx_sf)
            bandwidths_sf(a).lin = -1;
        else
            highBeta_sf = [sf(highSideIdx_sf-1), 1; sf(highSideIdx_sf), 1] \ sfResp(highSideIdx_sf-1:highSideIdx_sf); %slope and intercept for linear interpolation
            lowBeta_sf = [sf(lowSideIdx_sf), 1; sf(lowSideIdx_sf+1), 1] \ sfResp(lowSideIdx_sf:lowSideIdx_sf+1);
            highSF = (halfHeight_sf-highBeta_sf(2))./highBeta_sf(1);
            lowSF = (halfHeight_sf-lowBeta_sf(2))./lowBeta_sf(1);
            bandwidths_sf(a).lin = (log2(highSF)-log2(lowSF))/2;
        end
        
    else %there's no SF data
        disp('No SF data')
        bandwidths_sf(a).lin = nan;
    end
    
    %trying to catch bugs....
    fprintf('bandwidth for cell %d. (orient) = %.3f, (sf) = %.3f\n', a, bandwidths_orient(a), bandwidths_sf(a).lin)
    
    if PLOT
        %orientation bandwidth
        subplot(1,2,1)
        cla,
        hold on,
        plot(newOrients, newOriResp, '-bo')
        plot(orients, ones(size(newOrients)).*out.dat(a).grating.baselines(1), 'k:')
        plot(newOrients, ones(size(newOrients)).*halfHeight, 'm--')
        if ~isempty(highSideIdx) && ~isempty(lowSideIdx)
            plot([highOrient, highOrient], [0, halfHeight], 'm--')
            plot([lowOrient, lowOrient], [0, halfHeight], 'm--')
        else
            disp('No Linear Orient Bandwidth Estimate')
        end
        ylim([0, maxResp*1.1])
        
        %spatial frequency bandwidth
        subplot(1,2,2)
        cla,
        hold on,
        plot(sf, sfResp, '-bo')
        plot(sf, ones(size(sf)).*out.dat(a).grating.baselines(1), 'k:')
        plot(sf, ones(size(sf)).*halfHeight_sf, 'm--')
        if ~isempty(highSideIdx_sf) && ~isempty(lowSideIdx_sf)
            plot([highSF, highSF], [0, halfHeight_sf], 'm--')
            plot([lowSF, lowSF], [0, halfHeight_sf], 'm--')
        else
            disp('No Linear Orient Bandwidth Estimate')
        end
        ylim([0, maxResp_sf*1.1])
        drawnow
        pause(0.5)
    end
end

%run the analysis
params = {'TR', 'Orient Bandwidth', 'Exponent', 'CSI', 'pref SF', 'SF bandwidth', 'Mod Ratio', 'Eccentricity', 'RF Size'};
corMtx = nan(numel(l_valid), numel(params));
for a = 1:numel(params)
    tmp = [];
    switch params{a}
        case 'TR'
            tmp = min(rawTRs,[],2);
            if includeInfAndNan
                l_nans = isnan(tmp);
                tmp(l_nans) = max(tmp).*1.1;
            end
        case 'Orient Bandwidth'
            tmp = bandwidths_orient(:);
            if includeInfAndNan
                l_inf = isinf(tmp);
                tmp(l_inf) = max(tmp(~l_inf)) .* 1.1;
            end
        case 'Exponent'
            tmp = exponents(:);
            l_little = exponents<0.001;
            tmp(l_little) = 0.001;
        case 'CSI'
            tmp = rawCSIs(:);
            if includeInfAndNan
                l_inf = isinf(tmp);
                tmp(l_inf) = max(tmp(~l_inf)).*1.1;
            end
        case 'pref SF'
            tmp = rawSF(:);
        case 'SF bandwidth'
            for i = 1:numel(bandwidths_sf);
                if ~isempty(bandwidths_sf(i).lin)
                    tmp(i) = bandwidths_sf(i).lin;
                else
                    tmp(i) = nan; %happens when the sf tuning fxn doesn't go below half height
                end
            end
            valForLowpass = max(tmp(~isinf(tmp))) .* 1.1;
            tmp(isinf(tmp)) = valForLowpass;
            tmp(tmp<0) = nan; % bnadwidths that fall off the high side are degenerate
        case 'Mod Ratio'
            tmp = rawModRatios(:); %has nan values, but I can't do anything with them b/c it just didn't fire hard enough
            l_nan = isnan(tmp);
            tmp(l_nan) = 1;
        case 'Eccentricity'
            tmp = rawEccentricity(:);
        case 'RF Size'
            tmp = RFsize(:);
        otherwise
            error('do not know case: %s', params{a});
    end
    corMtx(:,a) = tmp;
end

nParams = numel(params);
[rho, p] = deal(nan(nParams, nParams));
offset = 1;
counter = 0;
figure
for r = 1:nParams
    offset = offset+1;
    for c = offset:nParams
        counter = counter+1;
        subplot(4, 9, counter)
        l_nans = sum(isnan(corMtx(:, [r,c])),2);
        plot(corMtx(l_valid,r), corMtx(l_valid,c), 'bo');
        [rho(r,c), p(r,c)] = corr(corMtx(l_valid&~l_nans,r), corMtx(l_valid&~l_nans,c), 'type', 'spearman');
        if any(strcmpi(params{r}, {'TR', 'Exponent', 'CSI', 'SF', 'Mod Ratio'}))
            set(gca, 'xscale', 'log')
        end
        if any(strcmpi(params{c}, {'TR', 'Exponent', 'CSI', 'SF', 'Mod Ratio'}))
            set(gca, 'yscale', 'log')
        end
        xlabel(params{r})
        ylabel(params{c})
        title(sprintf('%.2f, %.2f', rho(r,c), p(r,c)))
        
        %color the background according to the pvalue
        if p(r,c)<0.05
            set(gca, 'color', [1 .9 .9])
        elseif (p(r,c)<0.2) && (p(r,c)>0.05)
            set(gca, 'color', [.9 1 .9])
        elseif isnan(p(r,c))
            set(gca, 'color', [.85 .85 1])
        end
        axis tight
        
    end
end

%% (41) ISOLUM RESULTS CONDITIONED ON CP AND TR


clear, clc, close all
global cardVsIntBatchPath
preprocessDTbatchData

MAKEPLOT = 1;
ANALYZESF = 0;


%recover the grand CP estimates:
grandCP = nan(size(out.dat));
for a = 1:numel(out.dat)
    grandCP(a) = out.dat(a).c.cp.grandCP.val;
end


%pull up a figure of CP vs. TR. Allow the user to select the neurons to
%test further.
l_valid = ~commonExclusions;
TRs = min(rawTRs,[],2);
TRs(isnan(TRs)) = max(TRs(:))*1.1;
figure, hold on,
plot(TRs(l_valid), grandCP(l_valid), 'ko', 'markerfacecolor', 'k', 'markersize', 10);
set(gca, 'xscale', 'log')
xlim([min(TRs(:))*.95, max(TRs(:))*1.05])
ylim([min(grandCP)*.95, max(grandCP)*1.05])
set(gcf, 'position', [142 93 1072 691])
xlabel('Threshold Ratio')
ylabel('Choice Probability')
title('Select a group of cells')
[x,y] = ginput(1);
plot([x,x], [y, max(get(gca, 'ylim'))], 'r-', 'linewidth', 2)
plot([min(get(gca, 'xlim')),x], [y,y], 'r-', 'linewidth', 2)

%find the points that corresponed to the inside of the contour and plot
%them on the same axes as a sanity check
critVal_cp = y
critVal_tr = x
l_trs = sum(TRs < critVal_tr, 2) > 0; %taking the sum b/c one of the two TRs could fit the criteria, doing the == to turn it into a logical
l_cp = grandCP(:) > critVal_cp;
l_meetsCriteria = l_trs & l_cp & l_valid;

%verify that the cells chosen for the rest of the analysis are actually the
%ones you wanted to select with the ginput feature
plot(TRs(l_meetsCriteria), grandCP(l_meetsCriteria), 'bo', 'markerfacecolor', 'b', 'markersize', 8)


%
% now do compare TRs for the different colors
% cycle through the expts and pull out the indicies to the TRs for each
% color direction
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[sTR,lvmTR,swmTR,swlTR] = deal([]);
for a = find(l_meetsCriteria)';
    
    %only analyze one TR per neuron so that a kruskal wallis will work
    [~, idxToMin] = min(rawTRs(a,:));
    tmpColor = out.dat(a).expt.standColors(idxToMin,:);
    
    
    %deal with the cardinal TR
    if all([0 0 1] == sign(tmpColor))
        sTR(end+1,1) = rawCardTRs(a);
    elseif all([1 -1 0] == sign(tmpColor))
        lvmTR(end+1,1) = rawCardTRs(a);
    elseif all([1 -1 1] == sign(tmpColor));
        swlTR(end+1,1) = rawIntTRs(a);
    elseif all([1 -1 -1] == sign(tmpColor));
        swmTR(end+1,1) = rawIntTRs(a);
    end
end

%trying for a box plot
figure
set(gcf, 'position', [128         365        1218         419])
subplot(1,2,1)
hold on,
filtTRs = rawTRs(l_meetsCriteria,:);
popMedian = nanmedian(filtTRs(:));
plot([0.5 4.5], [popMedian, popMedian], 'k:')
tmp_TRs = [sTR; lvmTR; swmTR; swlTR];
group = char(repmat('Siso', length(sTR),1), repmat('LvM', length(lvmTR),1), repmat('SwM', length(swmTR),1), repmat('SwL', length(swlTR),1));
boxplot(tmp_TRs, group)
set(gca, 'yscale', 'log')
ylim([0.5, 4])
[p, atab, stats] = anova1(log(tmp_TRs), group, 'off');
title(sprintf('Main effect of color: p = %g', p));
ylabel('Threshold Ratio')
hold off

subplot(1,2,2)
hold on,
plot(1, sTR, 'bo', 'markerfacecolor', 'b', 'markersize', 8)
plot(1, geomean(sTR(~isnan(sTR))), 'ks', 'markerfacecolor', 'k', 'markersize', 10)
plot(2, lvmTR, 'ro', 'markerfacecolor', 'r', 'markersize', 8)
plot(2, geomean(lvmTR(~isnan(lvmTR))), 'ks', 'markerfacecolor', 'k', 'markersize', 10)
plot(3, swmTR, 'go', 'markerfacecolor', 'g', 'markersize', 8)
plot(3, geomean(swmTR(~isnan(swmTR))), 'ks', 'markerfacecolor', 'k', 'markersize', 10)
plot(4, swlTR, 'mo', 'markerfacecolor', 'm', 'markersize', 8)
plot(4, geomean(swlTR(~isnan(swlTR))), 'ks', 'markerfacecolor', 'k', 'markersize', 10)




%
% looking at the spatial frequency tuning of the most
% sensitive cells
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ANALYZESF
    
    dogPred = @(R0, SF, k_c, mu_c, sigma_c, k_s, mu_s, sigma_s) (R0 + (k_c .* exp(-((SF-mu_c)./(2.*sigma_c)).^2)) - (k_s .* exp(-((SF-mu_s)./(2.*sigma_s)).^2)));
    [bandwidths, prefSF] = deal(nan(size(out.dat)));
    [allSF, allResp] = deal([]);
    for a = find(l_meetsCriteria')
        GT = gtobj(out.fnames{a}{2});
        rates = GT.spikeResp('rate');
        R0 = out.dat(a).grating.baselines(1);
        l_sf = GT.trial(:, GT.idx.protocol) == 2;
        [Kc, SIGMAc, Ks, SIGMAs] = DoGfit(R0, GT.trial(l_sf, GT.idx.sf), rates(l_sf));
        SFpred = 0.3:0.001:5.2;
        Rpred = dogPred(R0, SFpred, Kc, 0, SIGMAc, Ks, 0, SIGMAs);
        
        
        %find the bandwidth
        [maxRate, maxIdx] = max(Rpred);
        err = abs(Rpred-(maxRate./2));
        [~,lowSideIdx] = min(err(1:maxIdx));
        [~,highSideIdx] = min(err(maxIdx:end));
        highSideIdx = highSideIdx+(maxIdx-1);
        if lowSideIdx == 1;
            bandwidths(a) = inf;
        elseif highSideIdx == numel(Rpred);
            bandwidths(a) = -1;
        elseif (lowSideIdx == 1) && (highSideIdx == numel(Rpred))
            bandwith(a) = nan;
        else
            bandwidths(a) = (log2(SFpred(highSideIdx)) - log2(SFpred(lowSideIdx)));
        end
        
        %store the mean responses away for future plotting
        allSF = [allSF; out.dat(a).grating.sf.stim(:)'];
        allResp = [allResp; out.dat(a).grating.sf.resp(:,1)'];
        
        %plot the tuning of each neuron
        if MAKEPLOT
            figure, hold on,
            plot(GT.trial(l_sf, GT.idx.sf), rates(l_sf), 'ko')
            plot(SFpred, Rpred, 'k', 'linewidth', 2)
            if ~(isempty(lowSideIdx) && isempty(highSideIdx))
                plot([SFpred(lowSideIdx), SFpred(highSideIdx)], [(maxRate./2),(maxRate./2)], 'm--')
            end
            set(gca, 'xscale', 'log')
            xlabel('Spatial Frequency')
            ylabel('Firing Rate')
            title(sprintf('CC: %s,  BW: %.3f', num2str(unique(GT.trial(l_sf, [GT.idx.lcc | GT.idx.mcc | GT.idx.scc]), 'rows')), bandwidths(a)))
        end
    end
    
    
    %plot the TRs vs Bandwidth, along side the average SF tuning.
    figure
    set(gcf, 'position', [203 338 1125 446])
    subplot(1,3,1), hold on,
    plot(allSF', allResp', '--', 'color', [1 1 1].*0.4)
    plot(allSF(1,:), mean(allResp,1), 'ko', 'markerfacecolor', 'k', 'markersize', 10)
    [Kc, SIGMAc, Ks, SIGMAs] = DoGfit(R0, allSF(:), allResp(:));
    SFpred = 0.3:0.001:5.2;
    Rpred = dogPred(R0, SFpred, Kc, 0, SIGMAc, Ks, 0, SIGMAs);
    plot(SFpred, Rpred, 'k', 'linewidth', 2)
    set(gca, 'xscale', 'log')
    title('(sub)population average')
    xlabel('Spatial Frequency')
    ylabel('Firing Rate')
    axis tight
    hold off
    
    subplot(1,3,2)
    normResp = bsxfun(@rdivide, allResp, max(allResp, [], 2));
    plot(allSF(1,:), mean(normResp,1), '-ko', 'markersize', 10, 'markerfacecolor', 'k', 'linewidth', 2)
    set(gca, 'xscale', 'log')
    title('Normalized SF tuning')
    xlabel('Spatial Frequency')
    ylabel('Norm Resp')
    
    
    subplot(1,3,3), hold on,
    tmp_bw = bandwidths;
    valForInf = max(tmp_bw(~isinf(tmp_bw))) * 1.1;
    tmp_bw(isinf(tmp_bw)) = valForInf;
    plot(TRs(l_meetsCriteria), tmp_bw(l_meetsCriteria), 'bo', 'markersize', 10, 'markerfacecolor', 'b')
    plot([min(TRs)*.95, max(TRs).*1.05], [0,0], ':r')
    plot([min(TRs)*.95, max(TRs).*1.05], [valForInf*.95, valForInf*.95], 'r:')
    xlim([min(TRs)*.95, max(TRs).*1.05])
    ylim([min(tmp_bw)*1.05, max(tmp_bw).*1.05])
    set(gca, 'xscale', 'log')
    xlabel('Threshold Ratio')
    ylabel('SF Bandwidth')
    [rho, p] = corr(TRs(l_meetsCriteria&(tmp_bw(:)>0)), bandwidths(l_meetsCriteria&(tmp_bw(:)>0))', 'type', 'spearman');
    title(sprintf('rho = %.3f, p = %.3f', rho, p))
end %if ANALYZESF

%% (42) CONDITION THE BOX AND WHISKER PLOTS ON DIFFERENT TR LEVELS

clear, clc, close all
global cardVsIntBatchPath
preprocessDTbatchData

%recover the grand CP estimates:
grandCP = nan(size(out.dat));
for a = 1:numel(out.dat)
    grandCP(a) = out.dat(a).c.cp.grandCP.val;
end

l_valid = ~commonExclusions;
TRs = min(rawTRs,[],2);
TRs(isnan(TRs)) = max(TRs)*1.1;
critTR = max(TRs);
deltaTR = 0.2;
timeout = 3;

figure
set(gcf, 'position', [128         365        1218         419])
while(1)
    if sum(TRs(l_valid) < critTR) == 0
        disp('done')
        break
    end
    
    exptList = l_valid & (TRs<critTR);
    
    [sTR,lvmTR,swmTR,swlTR] = deal([]);
    for a = find(exptList)';
        %deal with the cardinal TR
        if ismember([0 0 1], sign(out.dat(a).expt.standColors), 'rows')
            sTR(end+1,1) = rawCardTRs(a);
        elseif ismember([1 -1 0], sign(out.dat(a).expt.standColors), 'rows')
            lvmTR(end+1,1) = rawCardTRs(a);
        end
        
        %deal with the intermediate TR
        if ismember([1 -1 1], sign(out.dat(a).expt.standColors), 'rows')
            swlTR(end+1,1) = rawIntTRs(a);
        elseif ismember([1 -1 -1], sign(out.dat(a).expt.standColors), 'rows')
            swmTR(end+1,1) = rawIntTRs(a);
        end
    end
    
    try
        subplot(1,3,1)
        cla
        hold on,
        plot(TRs(l_valid), grandCP(l_valid), 'bo', 'markerfacecolor', 'b', 'markersize', 10)
        plot([critTR, critTR], [max(grandCP(l_valid)).*0.93, min(grandCP(l_valid))*1.07], '-r', 'linewidth', 2)
        set(gca, 'xscale', 'log')
        xlim([min(TRs(l_valid))*.95, max(TRs(l_valid))*1.1])
        hold off
        
        %plot the results for each round
        subplot(1,3,2)
        cla
        hold on,
        filtTRs = rawTRs(exptList,:);
        popAvg = nanmedian(filtTRs(:));
        plot([0.5 4.5], [popAvg, popAvg], 'k:')
        tmp_TRs = [sTR; lvmTR; swmTR; swlTR];
        group = char(repmat('Siso', length(sTR),1), repmat('LvM', length(lvmTR),1), repmat('SwM', length(swmTR),1), repmat('SwL', length(swlTR),1));
        boxplot(tmp_TRs, group)
        set(gca, 'yscale', 'log')
        ylim([0.5, 4])
        try
            pCard = ranksum(sTR(~isnan(sTR)), lvmTR(~isnan(lvmTR)));
            pInt = ranksum(swmTR(~isnan(swmTR)), swlTR(~isnan(swlTR)));
            title(sprintf('P_{card}=%.3f, P_{int}=%.3f', pCard, pInt));
        catch
        end
        hold off
        
        subplot(1,3,3)
        cla
        hold on,
        plot(1, sTR, 'bo', 'markerfacecolor', 'b', 'markersize', 8)
        plot(1, geomean(sTR(~isnan(sTR))), 'ks', 'markerfacecolor', 'k', 'markersize', 10)
        plot(2, lvmTR, 'ro', 'markerfacecolor', 'r', 'markersize', 8)
        plot(2, geomean(lvmTR(~isnan(lvmTR))), 'ks', 'markerfacecolor', 'k', 'markersize', 10)
        plot(3, swmTR, 'go', 'markerfacecolor', 'g', 'markersize', 8)
        plot(3, geomean(swmTR(~isnan(swmTR))), 'ks', 'markerfacecolor', 'k', 'markersize', 10)
        plot(4, swlTR, 'mo', 'markerfacecolor', 'm', 'markersize', 8)
        plot(4, geomean(swlTR(~isnan(swlTR))), 'ks', 'markerfacecolor', 'k', 'markersize', 10)
        title(sprintf('TR < %.3f', critTR))
        hold off
    catch
        disp('tried but failed...')
        break
    end
    
    %pause and update the critval
    critTR = critTR - deltaTR;
    keyboard
end

%% (43) CRFs FOR INT DATA PROJECTED ONTO CARD AXIS (HALF SQUARE FIT)
%
% This analysis cell does a bunch of stuff:
%
% (1) Fits CRFs with a model that assums a simple scaling between the card
% and intermediate directions. It starts by fitting them independently
% (full model), and then compares the liklihood to a fit of the CRFs jointly
% (reduced model).
%
% (2) Calculates the (interplated) firing rate at psychometric threshold,
% and the fitted firing rate at zero contrast. Plots these data as box and
% whiskers
%
% (3) Calculates 95% confidence intervals for the scaleFactor when using
% yoked fits.
%
% (4) Compares Beta3s to CSI, TR, and CP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, clc, close all
global cardVsIntBatchPath blpBatchPath
preprocessDTbatchData


MAKEPLOT = 0;
CONTRAST_DEF = 'proj onto card'; % could be: 'proj onto card' or 'psy thresh'
TRIALLENGTH = 0.666;  % in seconds
CARD = 0; %grouping variable
INT = 1;
l_validConds = ~(commonExclusions | (out.errors(:, neuroThresh1Ind) & out.errors(:,neuroThresh2Ind)));
exptList = find(l_validConds);
[devstat, devstat_shared, p, p_shared, scaleFactor] = deal(size(exptList));
CI = nan(numel(exptList),2);
[driven.S, driven.LvM, driven.SwL, driven.SwM] = deal([]); %for a comparison of baseline vs. driven rates
[Ro.S, Ro.LvM, Ro.SwL, Ro.SwM] = deal([]);
[TR.S, TR.LvM, TR.SwL, TR.SwM] = deal([]);
for a = 1:numel(exptList)
    ex = exptList(a);
    
    %deal with the cardinal light data (counts not rates)
    s_idx = ismember(sign(out.dat(ex).expt.standColors), [0 0 1], 'rows');
    lvm_idx = ismember(sign(out.dat(ex).expt.standColors), [1 -1 0], 'rows');
    card_Idx = s_idx | lvm_idx;
    if sum(card_Idx) ~= 1; error('didn''t get the correct number of cardinal colors'); end
    [counts, contrasts, colors] = deal([]); %for the halfSquareFit
    [cardContrasts, cardCounts, cardSEM] = deal([]); %for estimating spike threshold
    for cntrst = 1:length(out.dat(ex).c.crfIn{card_Idx})
        tmp_counts = out.dat(ex).c.crfIn{card_Idx}{cntrst}.*TRIALLENGTH;
        tmp_counts = round(tmp_counts);
        counts = [counts; tmp_counts(:)];
        switch lower(CONTRAST_DEF)
            case 'proj onto card'
                tmp_contrast = out.dat(ex).expt.norms{card_Idx}(cntrst);
            case 'psy thresh'
                tmp_contrast = out.dat(ex).expt.norms{card_Idx}(cntrst);
                tmp_contrast = tmp_contrast ./ out.dat(ex).m.alpha(card_Idx);
        end
        contrasts = [contrasts; repmat(tmp_contrast, length(tmp_counts), 1)];
        colors = [colors; repmat(CARD, length(tmp_counts), 1)];
        cardContrasts = [cardContrasts, tmp_contrast];
        cardCounts = [cardCounts, mean(tmp_counts)];
        cardSEM = [cardSEM, std(tmp_counts)./sqrt(numel(tmp_counts))];
    end
    
    %deal with the intermediate light data
    int_Idx = ~card_Idx;
    [intContrasts, intCounts, intSEM] = deal([]); %for estimating spike threshold
    for cntrst = 1:length(out.dat(ex).c.crfIn{int_Idx})
        tmp_counts = out.dat(ex).c.crfIn{int_Idx}{cntrst}.*TRIALLENGTH;
        tmp_counts = round(tmp_counts);
        counts = [counts; tmp_counts(:)];
        cardColor = out.dat(ex).expt.standColors(card_Idx,:);
        switch lower(CONTRAST_DEF)
            case 'proj onto card'
                cardUnit = cardColor ./ norm(cardColor);
                intColor = out.dat(ex).expt.standColors(int_Idx,:);
                intUnit = intColor./norm(intColor);
                intLMS = intUnit .* out.dat(ex).expt.norms{int_Idx}(cntrst);
                tmp_contrast = abs(intLMS * cardUnit(:));
            case 'psy thresh'
                tmp_contrast = out.dat(ex).expt.norms{int_Idx}(cntrst);
                tmp_contrast = tmp_contrast ./ out.dat(ex).m.alpha(int_Idx);
        end
        contrasts = [contrasts; repmat(tmp_contrast, length(tmp_counts), 1)];
        colors = [colors; repmat(INT, length(tmp_counts), 1)];
        intContrasts = [intContrasts, tmp_contrast];
        intCounts = [intCounts, mean(tmp_counts)];
        intSEM = [intSEM, std(tmp_counts)./sqrt(numel(tmp_counts))];
    end
    
    
    %estimate the spike threshold
    bkgndCounts = (out.dat(ex).c.crfIn{1}{1}.*TRIALLENGTH); %crfIn and crfOut have both in/out trials for zero contrast
    spkThresh = mean(bkgndCounts) + (std(bkgndCounts)./sqrt(numel(bkgndCounts)));
    cardThresh = cardContrasts(max(find(cardCounts > spkThresh, 1, 'first')-1, 1));
    intThresh = intContrasts(max(find(intCounts > spkThresh, 1, 'first')-1, 1));

    
    %start by fitting the CRFs separately
    l_card = colors == CARD;
    l_aboveThresh = contrasts>=cardThresh;
    tList = l_card&l_aboveThresh;
    b0_card = [ones(sum(tList),1), contrasts(tList).^2] \ counts(tList);
    params0 = [mean(bkgndCounts), cardThresh, b0_card(2)];
    [fit_card, lik_card] = halfSquareFit(contrasts(l_card), counts(l_card), params0, 'simple');
    
    l_int = colors == INT;
    l_aboveThresh = contrasts>=intThresh;
    tList = l_int&l_aboveThresh;
    b0_int = [ones(sum(tList),1), contrasts(tList).^2] \ counts(tList);
    params0 = [mean(bkgndCounts), intThresh, b0_int(2)];
    [fit_int, lik_int] = halfSquareFit(contrasts(l_int), counts(l_int), params0, 'simple');
    
    %fit the CRFs jointly (4 parameter fit)
    yoke_b0 = mean([fit_card(1), fit_int(1)]); %guess for baseline counts
    yoke_b1 = fit_card(2); %guess for spike threshold
    yoke_b2 = fit_card(3); %guess for gain
    yoke_b3 = sqrt(fit_int(3)/fit_card(3)); %guess for contrast scaling
    params0 = [yoke_b0, yoke_b1, yoke_b2, yoke_b3];
    [fit_yoke, lik_yoke, CI(a,:)] = halfSquareFit([contrasts, colors], counts, params0, 'yoked');
    scaleFactor(a) = fit_yoke(4);
    
    %fit the CRFs separately (again) but using the yoked fit as the initial
    %guess.
    params0 = fit_yoke(1:3);
    [fit_card2, lik_card2] = halfSquareFit(contrasts(l_card), counts(l_card), params0, 'simple');
    params0 = [fit_yoke(1), fit_yoke(2)./fit_yoke(3), fit_yoke(3).*fit_yoke(4)^2];
    [fit_int2, lik_int2] = halfSquareFit(contrasts(l_int), counts(l_int), params0, 'simple');
    
    
    %select the best separate fits based off their liklihoods
    if lik_int2 < lik_int
        lik_int = lik_int2;
        fit_int = fit_int2;
    end
    if lik_card2< lik_card
        lik_card = lik_card2;
        fit_card = fit_card2;
    end
        
    % deviance for yoked vs separate fits
    devstat(a) = 2*(lik_yoke - (lik_card+lik_int));
    p(a) = 1- chi2cdf(devstat(a),2);
    
    %fit both sets with a model that has a shared baseline
    params0 = [mean([fit_card(1), fit_int(1)]), fit_card(2:3), fit_int(2:3)]; 
    [fit_shared, lik_shared] = halfSquareFit([contrasts, colors], counts, params0, 'sharedBaseline');
    devstat_shared(a) = 2*(lik_yoke - lik_shared);
    p_shared(a) = 1- chi2cdf(devstat_shared(a),1);
    
    plotlist_bigScale = {'S030411002', 'S082709002', 'S052009004'};
    plotlist_ScaleEqOne = {'S042809003', 'K020311002'};
    if MAKEPLOT && any(strcmpi(plotlist_ScaleEqOne, out.fnames{ex}{1}))
        %set up the figure
        if strcmpi(CONTRAST_DEF, 'proj onto card')
            f(a) = figure;
        else
            figure(f(a))
        end
        
        %plot the data on normal contrast axes
        set(gcf, 'position', [213 573 1311 435])
        subplot(1,3,1), cla
        hold on
        pltClr = {'r', 'b'};
        for j = 1:2;
            colorDir{j} = num2str(out.dat(ex).expt.standColors(j,:));
            ccs = out.dat(ex).expt.norms{j};
            counts = cellfun(@(x,y) x.*y, out.dat(ex).c.crfIn{j}, repmat({0.666}, size(out.dat(ex).c.crfIn{j})), 'uniformoutput', 0);
            avg = cellfun(@mean, counts);
            sem = cellfun(@(x) std(x)./sqrt(numel(x)), counts);
            errorbar(ccs, avg, sem, 'color', pltClr{j})
        end
        set(gca, 'xscale', 'log')
        xlabel('vector norm CC')
        ylabel('counts')
        jointContrasts = [out.dat(ex).expt.norms{1}(:) ; out.dat(ex).expt.norms{2}(:)];
        minx = min(jointContrasts(jointContrasts~=0));
        maxx = max(jointContrasts);
        set(gca, 'xscale', 'log', 'xlim', [minx, maxx])
        legend(colorDir{1}, colorDir{2})
        
        %define contrast as projOnCard, or psyThresh. Then plot the data
        if strcmpi(CONTRAST_DEF, 'proj onto card')
            subplot(1,3,2), cla
            xlabel('proj onto card')
        else
            subplot(1,3,3), cla
            xlabel('psy thresh')
        end
        hold on
        
        
        %plot the raw data first
        errorbar(cardContrasts, cardCounts, cardSEM, 'ks', 'markerfacecolor', 'k')
        errorbar(intContrasts, intCounts, intSEM, 'rs', 'markerfacecolor','r')
        
        %CRFs fit separately
        xx_card = linspace(0, max(contrasts(l_card)), 100);
        xx_int = linspace(0, max(contrasts(l_int)), 100);
        pred_card = fit_card(1)+fit_card(3)*(max(xx_card-fit_card(2),0).^2);
        plot(xx_card, pred_card, 'k--')
        pred_int = fit_int(1)+fit_int(3)*(max(xx_int-fit_int(2),0).^2);
        plot(xx_int, pred_int, 'r--')
        
        % CRFs fit jointly
        pred_yoked_card = fit_yoke(1)+fit_yoke(3)*(max(xx_card-fit_yoke(2),0).^2);
        pred_yoked_int = fit_yoke(1)+fit_yoke(3)*(max(fit_yoke(4).*xx_int-fit_yoke(2),0).^2);
        plot(xx_card, pred_yoked_card, '-m')
        plot(xx_int, pred_yoked_int, '-m')
        
        % CRFs with sharedBaseline fit
        pred_shared_card = fit_shared(1)+fit_shared(3)*(max(xx_card-fit_shared(2),0).^2);
        pred_shared_int = fit_shared(1)+fit_shared(5)*(max(xx_int-fit_shared(4),0).^2);
        plot(xx_card, pred_shared_card, 'c')
        plot(xx_int, pred_shared_int, 'c')
        ylabel('Spike Count')
        %title(sprintf('Dev = %.3f p-val %.3f', devstat(a), p(a)))
        title(sprintf('Scale Factor: %.3f', scaleFactor(a)));
        set(gcf, 'name', out.fnames{ex}{1})
        
        %set the axes
        jointContrasts = [cardContrasts(:) ; intContrasts(:)];
        minx = min(jointContrasts(jointContrasts~=0));
        maxx = max([cardContrasts(:) ; intContrasts(:)]);
        set(gca, 'xscale', 'log', 'xlim', [minx, maxx])
        
        if strcmpi(CONTRAST_DEF, 'psy thresh')
            subplot(1,3,2)
            yylims = get(gca, 'ylim');
            for j = 1:3;
                subplot(1,3,j);
                set(gca, 'ylim', yylims)
            end
        end
            
    end
    
    if strcmpi(CONTRAST_DEF, 'psy thresh')
        %keep track of baseline rates and driven responses across colors
        [minTR, minTR_idx] = min(out.dat(ex).c.alpha ./ out.dat(ex).m.alpha);
        tmpColor = sign(out.dat(ex).expt.standColors(minTR_idx,:));
        if all(tmpColor == [0 0 1])
            Ro.S(end+1,1) = fit_shared(1);
            driven.S(end+1,1) = (fit_shared(1)+fit_shared(3)*(max(1-fit_shared(2),0).^2)) - fit_shared(1);
            TR.S(end+1,1) = minTR;
        elseif all(tmpColor == [1 -1 0])
            Ro.LvM(end+1,1) = fit_shared(1);
            driven.LvM(end+1,1) = (fit_shared(1)+fit_shared(3)*(max(1-fit_shared(2),0).^2)) - fit_shared(1);
            TR.LvM(end+1,1) = minTR;
        elseif all(tmpColor == [1 -1 1])
            Ro.SwL(end+1,1) = fit_shared(1);
            driven.SwL(end+1,1) = (fit_shared(1)+fit_shared(5)*(max(1-fit_shared(4),0).^2)) - fit_shared(1);
            TR.SwL(end+1,1) = minTR;
        elseif all(tmpColor == [1 -1 -1])
            Ro.SwM(end+1,1) = fit_shared(1);
            driven.SwM(end+1,1) = (fit_shared(1)+fit_shared(5)*(max(1-fit_shared(4),0).^2)) - fit_shared(1);
            TR.SwM(end+1,1) = minTR;
        end
    end
end

% plot the p-values and a histogram of the deviances
figure
subplot(2,2,1)
hist(devstat)
title('devstat yoked')
subplot(2,2,2)
hist(p)
title('p yoked')
subplot(2,2,3)
hist(devstat_shared)
title('devstat shared')
subplot(2,2,4)
hist(p_shared)
title('p shared')

%plot the scaleFactors on a log axis
figure, hold on,
bins = logspace(log10(min(scaleFactor)*.97), log10(max(scaleFactor)*1.04), 15);
nAll = histc(scaleFactor, bins);
nNonSig = histc(scaleFactor(p_shared>0.05), bins);
h_all = bar(bins, nAll, 'style', 'histc');
h_nonSig = bar(bins, nNonSig, 'style', 'histc');
set(h_all, 'edgecolor', 'w');
set(h_nonSig, 'edgecolor', 'w', 'facecolor', 'y')
set(gca, 'xscale', 'log', 'tickDir', 'out')
child = get(gca, 'children');
set(child([1,3]), 'visible', 'off')
legend('All neurons', 'Consistent with model')
xlabel('Scale Factor')
ylabel('Counts')
[h_popHist,p_popHist] = ttest(log10(scaleFactor(p_shared>0.05)));
title(sprintf('Different than one? p = %g', p_popHist))

%plot the scaleFactors with their confidence intervals
figure, hold on,
plot(CI', repmat([1:size(CI,1)], 2, 1), 'b')
plot(scaleFactor, [1:size(CI,1)], 'b.');
plot([1, 1], [0, size(CI,1)], 'm:')
xlim([-2, 4])
ylim([-2, size(CI,1)+2])
ylabel('Cell Number')
xlabel('95% CI for Scale Factor')
insideCI = sum((CI(:,1)<=1)&(CI(:,2)>=1));
biggerThanOne = sum(sum(CI>1,2) == 2);
smallerThanOne = sum(sum(CI<1,2) == 2);
title(sprintf('CIs that include 1: %g', insideCI))
hold off

%print out some usefull info
fprintf('\n\n\n\t **********  SUMMARY  ***********\n')
fprintf('Total number of neurons: %d\n', numel(out.dat))
fprintf('Number of neurons subjected to halfSquare analysis: %d\n', numel(devstat));
fprintf('Number of neurons consistent with model %d\n', sum(p_shared>0.05));
fprintf('Geometric mean scale factor for neurons that are consistent with model %.3f\n', geomean(scaleFactor(p>0.05)));


%make a boxplot of driven rates by color if need be.
if strcmpi(CONTRAST_DEF, 'psy thresh')
    tmp_driven = [driven.S ; driven.LvM ; driven.SwL ; driven.SwM];
    tmp_Ro =  [Ro.S ; Ro.LvM ; Ro.SwL ; Ro.SwM];
    tmp_TRs = [TR.S ; TR.LvM ; TR.SwL ; TR.SwM];
    groups = [repmat('Siso', numel(Ro.S), 1) ; repmat(' LvM', numel(Ro.LvM), 1); repmat(' SwL', numel(Ro.SwL), 1); repmat(' SwM', numel(Ro.SwM), 1)];
    
    figure
    
    subplot(1,3,1)
    boxplot(tmp_Ro, groups)
    [p_Ro, atab_Ro, stats_Ro] = anova1(tmp_Ro, groups, 'off');
    title(sprintf('Baseline Rate, p = %g', p_Ro))
    ylabel('Spike Count')
    
    subplot(1,3,2)
    boxplot(tmp_driven, groups)
    [p_driven, atab_driven, stats_driven] = anova1(tmp_driven, groups, 'off');
    title(sprintf('Driven Rate p = %g', p_driven))
    ylabel('Spike Count')
    
    subplot(1,3,3)
    boxplot(tmp_TRs, groups)
    [p_TRs, atab_TRs, stats_TRs] = anova1(tmp_TRs, groups, 'off');
    title(sprintf('Threshold Ratios p = %g', p_TRs))
    ylabel('Threshold Ratio')
    set(gca, 'yscale', 'log')
end

% make a plot of the Beta3 (i.e., the scale factor) as a function of TR and
% CSI
TRs = min(rawTRs(l_validConds,:),[],2);
CSI = rawCSIs(l_validConds);
figure
subplot(1,2,1)
plot(scaleFactor, TRs, 'k.')
set(gca, 'yscale', 'log', 'xscale', 'log')
xlabel('Beta 3')
ylabel('TR')
[rho_TR,p_TR] = corr(TRs, scaleFactor')
subplot(1,2,2)
plot(scaleFactor, CSI, 'k.')
set(gca, 'yscale', 'log', 'xscale', 'log')
xlabel('Beta 3')
ylabel('CSI')
[rho_CSI,p_CSI] = corr(CSI, scaleFactor', 'type', 'spearman')


%%

logB3 = log10(scaleFactor);
SD = std(logB3);
xbar = mean(logB3);

l_consistent = (logB3 > -SD) & (logB3 <= SD)
l_bigger = logB3 > SD

[h_CSI, p_CSI] = ranksum(CSI(l_consistent), CSI(l_bigger))
[h_TR, p_TR] = ttest2(TRs(l_consistent), TRs(l_bigger))


l_CO = CSI>=2;
l_other = CSI<2;
[h_B3, p_B3] = ttest2(logB3(l_CO), logB3(l_other))
    

%% CRF for example cell
fin
DT = dtobj('S042809003');
[m, c, expt] = DTunpack(DT);

figure, hold on,
pltClr = {'r', 'b'}
for a = 1:2;
    colorDir{a} = num2str(expt.standColors(a,:));
    ccs = expt.norms{a};
    counts = cellfun(@(x,y) x.*y, c.crfIn{a}, repmat({0.666}, size(c.crfIn{a})), 'uniformoutput', 0);
    avg = cellfun(@mean, counts);
    sem = cellfun(@(x) std(x)./sqrt(numel(x)), counts);
    errorbar(ccs, avg, sem, 'color', pltClr{a})
end

set(gca, 'xscale', 'log')
legend(colorDir{1}, colorDir{2})
axis tight


%% CELLS THAT ARE INSENSITIVE TO ONE COLOR...

clear, clc, close all
global cardVsIntBatchPath
preprocessDTbatchData


l_oneColor = sum(isnan(rawNTs),2) == 1;

for a = find(l_oneColor)'
    figure, hold on,
    for clr = 1:2
        if all(sign(out.dat(a).prefCard) == sign(out.dat(a).expt.standColors(clr,:)))
            plot(out.dat(a).expt.norms{clr}(2:end), out.dat(a).c.roc{clr}(2:end), '-ro')
        else
            plot(out.dat(a).expt.norms{clr}(2:end), out.dat(a).c.roc{clr}(2:end), '-ko')
        end
    end
    set(gca,'xscale', 'log')
end



