function [thresholds, colorDirs, sfs] = DTquestUnpack(stro, fitMeth, WINSIZE, PERFRANGE)
%
%   EXAMPLE:  [thresholds, colorDirs, sfs] = DTquestUnpack(stro, fitMeth, [WINSIZE],[PERFRANGE]);
%
% unpacks QUEST experiments and returns the threshold estimate as specified
% by 'fitMeth' (a string argument either 'mode' or 'weibull'). WINSIZE is
% an optional argument specifing the number of trials over which a running
% average of behavioral performance should be computed. The default size is
% 10 trials. PERFRANGE is an optional argument that allows the user to
% filter out QUEST runs based on the subject's asympotic performance. For
% example, if PERFRANGE is set to 0.1 than QUEST runs with asympotic
% performace outside the range of 82 +/- 10 percent correct will be
% rejected.
%
% CAH 03/09


%define the optional argumentS if they're not given
if ~exist('WINSIZE', 'var')
    WINSIZE = 10;
end
if ~exist('PERFRANGE', 'var')
    PERFRANGE = []; %by default, don't filter
end

%unpack the necessary experimental params.
DTindicies;
bkgndrgb = [stro.sum.exptParams.bkgnd_r, stro.sum.exptParams.bkgnd_g, stro.sum.exptParams.bkgnd_b];
M = reshape(stro.sum.exptParams.m_mtx, 3, 3);
bkgndlms = M * bkgndrgb';
x = 0:255; %the normal range of the gamma look up table
xx = linspace(0, 255, 2^16); %the desired quantization of the gammaTable
g1 = reshape(stro.sum.exptParams.gamma_table, 256, 3);
gammaTable = [spline(x, g1(:,1), xx)', spline(x, g1(:,2), xx)', spline(x, g1(:,3), xx)'];

colorDirs = reshape(stro.sum.exptParams.RF_colors, 3,3)';
colorDirs(softEq([0 0 0], colorDirs, [], 'rows'), :) = [];
colorRanges = reshape(stro.sum.exptParams.quest_ranges, 3, 3)';
colorRanges(softEq([0 0 0], colorRanges, [], 'rows'), :) = [];
spatialPeriods = sort(unique(stro.trial(:,gaborLambdaInd)));
sfs = stro.sum.exptParams.pixperdeg ./ spatialPeriods;
figure;
for clr = 1:size(colorDirs, 1);
    yRange = [0 0]; %for plotting
    for spfr = 1:length(spatialPeriods);
        %reconstruct the prior
        sigma = stro.sum.exptParams.quest_sigma;
        beta = stro.sum.exptParams.quest_beta;
        range = [colorRanges(clr, 1):colorRanges(clr, 3):colorRanges(clr, 2)];
        guess = stro.sum.exptParams.low_scalars(clr);
        prior = (1./sigma*sqrt(2*pi)) .* exp(-(range-guess).^2 ./ (2*sigma^2));
        Q = log(prior);
        
        %pull out the appropriate trials
        l_sfs = stro.trial(:,gaborLambdaInd) == spatialPeriods(spfr);
        l_clr = stro.trial(:,colorDirInd) == clr;
        tList = find(l_sfs & l_clr);
        nTrials = length(tList);
        threshByTrial = nan(1, nTrials);
        
        %reconstruct the trials
        unitVec = colorDirs(clr,:) ./ norm(colorDirs(clr,:));
        CCtrl = nan(1,nTrials);
        for trl = 1:nTrials;
            [m, idx] = max(Q);
            CCsim = range(idx); %the expected CC based on my reconstruction
            RGB = stro.trial(tList(trl),flashRInd:flashBInd)+1;
            rgb = [gammaTable(RGB(1), 1), gammaTable(RGB(2), 2), gammaTable(RGB(3), 3)];
            LMS = ((M*rgb')-bkgndlms) ./ bkgndlms;
            CCtrl(trl) = norm(LMS).*100; %the actual CC presented
            
            %update the quest function
            liklihood = 1 - 0.5.*exp(-(CCtrl(trl)./range).^beta);
            if ~stro.trial(tList(trl), correctInd)
                liklihood = 1-liklihood;
            end
            Q = Q + log(liklihood);
            threshByTrial(trl) = CCsim;
        end
        
        %assign the final threshold estimate based on the modal value of
        %the quest function
        [m, idx] = max(Q);
        questAlpha(clr, spfr) = range(idx);
        
        %estimate threshold a second time based on a fit to a cumulative weibull.
        CCtmp = CCtrl;%./100; %for compatibility wth mocs weibull fitting
        [aSSE, bSSE, gSSE] = weibullFit(CCtmp, stro.trial(tList, correctInd), 'sse' ,[questAlpha(clr, spfr), 85]);
        correctByContrast = stro.trial(tList, correctInd);
        wrongByContrast = ~correctByContrast;
        [weibullAlpha(clr, spfr), beta(clr, spfr), gamma(clr, spfr), success] = weibullFit(CCtmp, [correctByContrast(:), wrongByContrast(:)], 'mle', [aSSE, bSSE]);
        if success <= 0;
            weibullAlpha(clr, spfr) = NaN;
        end
        
        %plot the raw data
        filter = ones(1,WINSIZE)./WINSIZE;
        perf = conv(stro.trial(tList, correctInd), filter);
        perf(end-(WINSIZE-2):end) = []; %convolution is zero padded on the tails
        perf(1:(WINSIZE-1)) = nan; %these values don't reflect a true running avg
        subplot(size(colorDirs, 1), length(sfs), (clr-1).*length(sfs)+spfr)
        hands{spfr} = plotyy([1:nTrials], threshByTrial, [1:nTrials], perf);
        set(get(min(hands{spfr}), 'children'), 'color', 'k', 'marker', '.', 'linewidth', 0.25);
        set(get(max(hands{spfr}), 'children'), 'color', 'b', 'marker', '.', 'linewidth', 0.25);
        set(hands{spfr}, 'ycolor', 'k')
        yLimTmp = get(min(hands{spfr}), 'ylim');
        if yLimTmp(1) < yRange(1)
            yRange(1) = yLimTmp(1);
        end
        if yLimTmp(2) > yRange(2)
            yRange(2) = yLimTmp(2);
        end
        
        %filter out runs where the subjects performace was outside some
        %predetermined range
        if ~isempty(PERFRANGE)
            if (perf(end) < (0.8-PERFRANGE)) || (perf(end) > (0.8+PERFRANGE))
                questAlpha(clr, spfr) = nan;
                weibullAlpha(clr, spfr) = nan;
            end
        end
        
        %now assign the output appropriately
        switch lower(fitMeth)
            case 'mode'
                thresholds(clr, spfr) = questAlpha(clr, spfr);
            case 'weibull'
                thresholds(clr, spfr) = weibullAlpha(clr, spfr);
        end

    end
    %set all the axes to have the same range
    for ax = 1:length(spatialPeriods);
        set(min(hands{ax}), 'box', 'off', 'ylim', yRange, 'ytickmode', 'auto');
        set(max(hands{ax}), 'box', 'off', 'ylim', [0.4 1], 'ytickmode', 'auto');
        title = sprintf('[%.2f, %.2f, %.2f] %.3f cyc/deg', unitVec(1), unitVec(2), unitVec(3), sfs(ax));
        set(get(max(hands{ax}), 'title'), 'string', title);
        if ax < length(spatialPeriods)
            set(max(hands{ax}), 'yticklabel', [])
        elseif ax > 1
            set(min(hands{ax}), 'yticklabel', [])
        end
        if ax == 1;
            set(get(min(hands{ax}), 'ylabel'), 'string', 'Thresh. Estimate (%CC)');
        elseif ax == length(spatialPeriods)
            set(get(max(hands{ax}), 'ylabel'), 'string', 'Performance');
        end
    end
end
    

