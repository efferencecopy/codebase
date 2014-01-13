%% (1)  BEGIN: DESGINATE A BATCH DATA FILE

% clear out the currrent workspace
fin

% designate a text file for the batch analysis:
CARDVSINT = 'CardVsInt_06-May-2011.mat'; % the typical one: 'CardVsInt_06-May-2011.mat'
global cardVsIntBatchPath

if isunix %charlie's laptop
    prefix = '/Users/charliehass/LabStuff/NexFiles/Charlie/Batch Data And Text Files/';
elseif strcmpi(license,  '367516') % when logged into cahass at wanprc on the server
    prefix = 'C:\NO BACKUP\NexFiles\Charlie\Batch Data And Text Files\';
else %assume all other pc's
    prefix = 'N:\NexFiles\Charlie\Batch Data And Text Files\';
end
cardVsIntBatchPath = [prefix, CARDVSINT];



%% HALF SQUARE FIT
%
% One of the reviewers asked us to split the data set from the cardinal
% color direction into two groups and to fit our model. Under these
% conditions, the distribution of Beta3 should be centered at one.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, clc, close all
global cardVsIntBatchPath %#ok<*REDEF>
preprocessDTbatchData

NITERS = 1;
WHICHCOLOR = 'card';
CLR1 = 0;  % just a grouping variable.
CLR2 = 1;
TRIALLENGTH = 0.666;  % in seconds
l_validConds = ~(commonExclusions | (out.errors(:, neuroThresh1Ind) & out.errors(:,neuroThresh2Ind)));
exptList = find(l_validConds);
CI = nan(numel(exptList),2);
meanScaleFactors = nan(numel(exptList),1);

for a = 1:numel(exptList)
    ex = exptList(a);
    
    % determine which color to use
    switch WHICHCOLOR
        case 'card'
            s_idx = ismember(sign(out.dat(ex).expt.standColors), [0 0 1], 'rows');
            lvm_idx = ismember(sign(out.dat(ex).expt.standColors), [1 -1 0], 'rows');
            clr_idx = s_idx | lvm_idx;
            if sum(clr_idx) ~= 1; error('didn''t get the correct number of cardinal colors'); end
        case 'int'
            swm_idx = ismember(sign(out.dat(ex).expt.standColors), [1 -1 -1], 'rows');
            swl_idx = ismember(sign(out.dat(ex).expt.standColors), [1 -1 1], 'rows');
            clr_idx = swm_idx | swl_idx;
            if sum(clr_idx) ~= 1; error('didn''t get the correct number of cardinal colors'); end
    end
    
    [devstat, p, scaleFactors] = deal( nan(NITERS,1));
    for iter = 1:NITERS
        [counts, contrasts, colors] = deal([]); %for the halfSquareFit
        [contrasts_clr1, contrasts_clr2, counts_clr1, counts_clr2] = deal([]); %for estimating spike threshold
        for cntrst = 1:length(out.dat(ex).c.crfIn{clr_idx})
            
            % extract the data and split into two roughly equal groups
            tmp_contrast = out.dat(ex).expt.norms{clr_idx}(cntrst);
            tmp_counts = out.dat(ex).c.crfIn{clr_idx}{cntrst}.*TRIALLENGTH;
            tmp_counts = round(tmp_counts);
            nCounts = numel(tmp_counts);
            tmp_counts = tmp_counts(randperm(nCounts)); % so that multipule iterations don't give rise to the same grouping b/w clr1 and clr2
            
            
            % allocate some data to "color 1"
            N_clr1 = round((nCounts+unifrnd(-0.5, 0.5))/2); % adding a randome number b/w +/- 0.5 so that clr1 doesn't always get more trials (when there are an odd number of trials)
            contrasts = [contrasts; repmat(tmp_contrast, N_clr1, 1)];
            colors = [colors; repmat(CLR1, N_clr1, 1)];
            counts = [counts; tmp_counts(1:N_clr1)'];
            contrasts_clr1 = [contrasts_clr1, tmp_contrast];
            counts_clr1 = [counts_clr1, mean(tmp_counts(1:N_clr1))];
            
            
            % allocate some data to "color 2"
            N_clr2 = nCounts - N_clr1;
            contrasts = [contrasts; repmat(tmp_contrast, N_clr2, 1)];
            colors = [colors; repmat(CLR2, N_clr2, 1)];
            counts = [counts; tmp_counts(N_clr1+1:end)'];
            contrasts_clr2 = [contrasts_clr2, tmp_contrast];
            counts_clr2 = [counts_clr2, mean(tmp_counts(N_clr1+1:end))];
            
        end
        
        
        %estimate the spike threshold
        bkgndCounts = (out.dat(ex).c.crfIn{1}{1}.*TRIALLENGTH); %crfIn and crfOut have both in/out trials for zero contrast
        spkThresh = mean(bkgndCounts) + (std(bkgndCounts)./sqrt(numel(bkgndCounts)));
        thresh_clr1 = contrasts_clr1(max(find(counts_clr1 > spkThresh, 1, 'first')-1, 1));
        if isempty(thresh_clr1); thresh_clr1 = contrasts_clr1(4); end %arbitrary threshold when the neruon doesn't respond over the range tested
        thresh_clr2 = contrasts_clr2(max(find(counts_clr2 > spkThresh, 1, 'first')-1, 1));
        if isempty(thresh_clr2); thresh_clr2 = contrasts_clr2(4); end %arbitrary threshold when the neruon doesn't respond over the range tested
        
        %start by fitting the CRFs separately
        l_clr1 = colors == CLR1;
        l_aboveThresh = contrasts>=thresh_clr1;
        tList = l_clr1 & l_aboveThresh;
        b0_clr1 = [ones(sum(tList),1), contrasts(tList).^2] \ counts(tList);
        params0 = [mean(bkgndCounts), thresh_clr1, b0_clr1(2)];
        [fit_clr1, lik_clr1] = halfSquareFit(contrasts(l_clr1), counts(l_clr1), params0, 'simple');
        
        l_clr2 = colors == CLR2;
        l_aboveThresh = contrasts >= thresh_clr2;
        tList = l_clr2 & l_aboveThresh;
        b0_clr2 = [ones(sum(tList),1), contrasts(tList).^2] \ counts(tList);
        params0 = [mean(bkgndCounts), thresh_clr2, b0_clr2(2)];
        [fit_clr2, lik_clr2] = halfSquareFit(contrasts(l_clr2), counts(l_clr2), params0, 'simple');
        
        %fit the CRFs jointly (4 parameter fit)
        yoke_b0 = mean([fit_clr1(1), fit_clr2(1)]); %guess for baseline counts
        yoke_b1 = fit_clr1(2); %guess for spike threshold
        yoke_b2 = fit_clr1(3); %guess for gain
        yoke_b3 = sqrt(fit_clr2(3)/fit_clr1(3)); %guess for contrast scaling
        params0 = [yoke_b0, yoke_b1, yoke_b2, yoke_b3];
        [fit_yoke, lik_yoke, CI(a,:)] = halfSquareFit([contrasts, colors], counts, params0, 'yoked');
        scaleFactors(iter) = fit_yoke(4);
    end
    meanScaleFactors(a) = mean(scaleFactors);
end

%plot the scaleFactors on a log axis
figure
edges = logspace(log10(min(meanScaleFactors).*0.9), log10(max(meanScaleFactors).*1.05), 20);
countsPerBin = histc(meanScaleFactors, edges);
b = bar(edges, countsPerBin, 'type', 'histc');
set(gca, 'xscale', 'log', 'tickDir', 'out')
children = get(gca, 'children');
set(children(1), 'visible', 'off')
xlabel('mean \beta_{3} across iterations')
ylabel('count')
title(sprintf('N iters = %d, color dir = %s', NITERS, WHICHCOLOR))

%% CARDINAL MECHANISMS MODEL ASSESSED BEHAVIORALLY
%
% One of the reviewers wanted to see a test of the cardinal model just
% using our behavioral data. Would the thresholds for intermediate color
% directions be predicted on the basis of cardinal mechanisms at all SFs
% for the contrast sensitivity functions we measured?
%
%%%%%%%%%%%%%%%%%%%%


%clear out all the junk that's in the workspace.
clear, close all, clc

% PICK AN OBSERVER:
observer{1} = '/Users/charliehass/LabStuff/NexFiles/Charlie/Kali/text files/questCSFdata.txt'; %kali
observer{2} = '/Users/charliehass/LabStuff/NexFiles/Charlie/Sedna/text files/quest.txt'; %sedna
observer{3} = '/Users/charliehass/LabStuff/NexFiles/Charlie/CharliePsychophysics/Text Files/charlieTrainingSet.txt'; %charlie

for a = 1:length(observer)
    
    % STEP TWO: batch process the data:
    nTrials = 15;
    perfRange = []; %don't filter on the baisis of performance
    [colors{a}, sfs{a}, data{a}] = questBatchProcess(observer{a}, 'mode', nTrials, perfRange);
    
    % STEP THREE: remove spatial frequencies that don't fall on the traditional latice
    ind = softEq(0.8893, sfs{a}, 4) | softEq(0.9919, sfs{a}, 4);
    sfs{a};
    sfs{a}(ind) = [];
    data{a}(:,ind,:) = [];
end

%
% STEP FIVE: plot the CSF and make sure that there are enough experiments
% per condition
%
%%%%%%%%%%%%%%%%%%%%
for a = 1:length(observer)
    order = 2;
    f(a) = fitCSF(colors{a}, sfs{a}, data{a}, order, 'all');
    nExpts = sum(~isnan(data{a}), 3);
    title(observer{a});
    xlabel(sprintf('Min Num Expts: %d', min(nExpts(:))))
end


%adjust the axes so that all of them are identical
SEM = cellfun(@(x) nanstd((1./(x./100)),0,3)./sqrt(sum(~isnan(x), 3)), data, 'uniformoutput', 0);
upperBound = cellfun(@(x,y) nanmean((1./(x./100)),3)+y, data, SEM, 'uniformoutput', 0);
lowerBound = cellfun(@(x,y) nanmean((1./(x./100)),3)-y, data, SEM, 'uniformoutput', 0);
maxY = max(cellfun(@(x) max(max((x))), upperBound));
minY = min(cellfun(@(x) min(min((x))), lowerBound));
for a = 1:length(f)
    figure(f(a).hand)
    set(gca, 'ylim', [minY.*0.98, maxY.*1.1], 'linewidth', 2, 'fontSize', 24)
    box off
    children = get(gca, 'children');
    set(children, 'linewidth', 4, 'markersize', 38)
end

%
% STEP SIX: cycle through the observers and compare the mean QUEST
% sensitivity for the INT colors to what we predict under the cardinal
% model. Do this separately for each spatial frequency
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
names = {'Monkey K', 'Monkey S', 'Human'};
figure
set(gcf, 'position', [-26           5        1467         824])
for a = 1:numel(observer);
   % find the mean thresholds
   meanThresh = nanmean( data{a}./100 , 3 );
   
   % find indicies to the color directions
   s_idx = ismember(sign(colors{a}), [0 0 1], 'rows');
   lvm_idx = ismember(sign(colors{a}), [1 -1 0], 'rows');
   swm_idx = ismember(sign(colors{a}), [1 -1 -1], 'rows');
   swl_idx = ismember(sign(colors{a}), [1 -1 1], 'rows');
   
   % convert the intermediate thresholds into vectors with the appropriate
   % length, and then project them onto the cardinal colors. Divide the
   % projection by the threshold in the cardinal color direction
   swm_vecs = [colors{a}(swm_idx,:)' * meanThresh(swm_idx,:)]';
   swmProjOnS = abs(swm_vecs * [0;0;1]) ./ meanThresh(s_idx,:)';
   swmProjOnLvM = abs(swm_vecs * colors{a}(lvm_idx,:)') ./ meanThresh(lvm_idx,:)';
   
   swl_vecs = [colors{a}(swl_idx,:)' * meanThresh(swl_idx,:)]';
   swlProjOnS = abs(swl_vecs * [0;0;1]) ./ meanThresh(s_idx,:)';
   swlProjOnLvM = abs(swl_vecs * colors{a}(lvm_idx,:)') ./ meanThresh(lvm_idx,:)';   
   
   % order the data by sf.
   [~,inds] = sort(sfs{a}, 2, 'ascend');
   [~, maxSfIdx] = max(sfs{a});
    
   subplot(1,3,a)
   set(gca, 'fontsize', 18)
   hold on,
   h = area([0,1], [1,1]);
   set(h, 'faceColor', [.9 .9 .9])
   plot(swmProjOnLvM(inds), swmProjOnS(inds), '-kv', 'markersize', 10, 'linewidth', 1, 'markerfacecolor', 'none')
   plot(swlProjOnLvM(inds), swlProjOnS(inds), '--ks', 'markersize', 10, 'linewidth', 1, 'markerfacecolor', 'k')
   maxVal = max([swmProjOnLvM; swlProjOnLvM; swmProjOnS; swlProjOnS]) .*1.07;
   set(gca, 'xlim', [0, maxVal], 'ylim', [0, maxVal])
   axis square
   if a == 2; xlabel('Projection onto the L-M axis (threshold units)'); end
   if a == 1; ylabel('Projection onto the S-axis (threshold units)'); end
   title(names{a})
   legend('', 'S with M', 'S with L', 'location', 'northwest')
   legend boxoff
   text(swlProjOnLvM(maxSfIdx).*1.05, swlProjOnS(maxSfIdx), [num2str(sfs{a}(maxSfIdx), 2), ' cpd'])
   
   hold off
end


%% CHOICE PROBABILITY
%
% Trying to justify why we didn't see significant CP in more cells
%
%%%%%%%%%%%%%%%%%%%%

clear, clc, close all
global cardVsIntBatchPath %#ok<*REDEF>
preprocessDTbatchData

% cycle through the data and determine the mean firing rate for the
% conditions in which CP was calculated.
meanRatesForCP = nan(numel(out.dat), 16);
choiceRatios = nan(numel(out.dat), 16);
for a = 1:numel(out.dat)
    meanRates = [cellfun(@mean, out.dat(a).c.crfIn{1}), cellfun(@mean, out.dat(a).c.crfIn{2})];
    validConds = ~[isnan(out.dat(a).c.cp.in{1}), isnan(out.dat(a).c.cp.in{2})];
    meanRates(~validConds) = nan;
    meanRatesForCP(a,:) = meanRates;
    
    T1T2ratio = [out.dat(1).c.cp.N_T1T2ratio{1}, out.dat(1).c.cp.N_T1T2ratio{2}];
    T1T2ratio(~validConds) = nan;
    choiceRatios(a,:) = T1T2ratio;
end

rateAmongCPConditions = nanmean(nanmean(meanRatesForCP,2))
ratioAmongCPConditions = nanmean(nanmean(choiceRatios,2))




%% PSTH's JUSTIFY THE USE OF F0 AMPLITUDE

clear, clc, close all
global cardVsIntBatchPath %#ok<*REDEF>
preprocessDTbatchData

CELLTYPES = 'simple';


l_valid = ~(commonExclusions | out.errors(:, neuroThresh1Ind) | out.errors(:, neuroThresh2Ind) | out.errors(:, lt5GTtrialsInd) | out.errors(:, prefTieInd) | out.errors(:, cardMismatchInd) | out.errors(:, intMismatchInd));
switch CELLTYPES
    case 'simple'
        l_cellType = rawModRatios > 1;
    case 'complex'
        l_cellType = rawModRatios < 1;
end
l_valid = l_valid & l_cellType;

for a = find(l_valid)'
    
    % open the GT file
    GT = gtobj(out.fnames{a}{2});
    LMS = GT.trial(:,[GT.idx.lcc | GT.idx.mcc | GT.idx.scc]);
    l_color = ismember(sign(LMS), sign(out.dat(a).prefIsolum), 'rows');
    l_p4 = getTrialList(GT, 4);
    tList = l_p4 & l_color;
    zeroTimes = mat2cell(GT.trial(tList, GT.idx.stimon), ones(sum(tList),1),1);
    spikeTimes_gt = cellfun(@(x,y)(x-y), GT.ras(tList, 1), zeroTimes, 'UniformOutput', 0);
    
    figure
    set(gcf, 'position', [29   585   834   244])
    subplot(1,2,1)
    hold on,
    counter = mat2cell([0:sum(tList)-1]', ones(sum(tList), 1), 1);
    cellfun(@(x, y)plot([x, x]', [zeros(1,length(x))+y; [ones(1, length(x)).*0.8 + y]], 'k'), spikeTimes_gt, counter);
    xlim([-0.1, 1.1])
    hold off
    
    
    
    % open the DT file
    DT = dtobj(out.fnames{a}{1});
    prefColorIdx = find(ismember(sign(out.dat(a).expt.standColors), sign(out.dat(a).prefIsolum), 'rows'));
    l_color = DT.trial(:, DT.idx.colorDir) == prefColorIdx;
    l_maxContrast = DT.trial(:, DT.idx.cntrstLev) == max(DT.trial(:, DT.idx.cntrstLev));
    rfx = DT.sum.exptParams.rf_x;
    rfy = DT.sum.exptParams.rf_y;
    l_inRF = (DT.trial(:,DT.idx.flashX) == rfx) & (DT.trial(:,DT.idx.flashY) == rfy);
    tList = l_color & l_maxContrast & l_inRF;
    
    % get the DT spike times, make the psth
    % pull out the spike times
    zeroTimes = mat2cell(DT.trial(tList, DT.idx.repFlashOn), ones(sum(tList),1),1);
    spikeTimes = cellfun(@(x,y)(x-y), DT.ras(tList, 1), zeroTimes, 'UniformOutput', 0);
    edges = linspace(-0.200, 0.800, 70);
    binwidth = edges(2)-edges(1);
    countsPerBin = cellfun(@(x,y) histc(x,y), spikeTimes, repmat({edges}, sum(tList), 1), 'uniformoutput', 0);
    countsPerBin = cellfun(@(x)(x(:)'), countsPerBin, 'uniformoutput', 0); %allign for vertcat
    countsPerBin = vertcat(countsPerBin{:});
    avgCountsPerBin = mean(countsPerBin) ./ binwidth;
    
    
    subplot(1,2,2)
    bar(edges, avgCountsPerBin, 'type', 'histc')
    xlim([-0.2, 0.7])
    
end









