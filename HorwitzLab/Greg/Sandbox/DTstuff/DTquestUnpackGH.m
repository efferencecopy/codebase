function [thresholds, colorDirs, sfs, QuestTrajectories] = DTquestUnpackGH(stro, fitMeth, WINSIZE, PERFRANGE);
%
%   EXAMPLE:  [thresholds, colorDirs, sfs, QuestTrajectories] = DTquestUnpackGH(stro, fitMeth, [WINSIZE],[PERFRANGE]);
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
colorDirs(sum(abs(colorDirs), 2)==0, :) = [];
colorRanges = reshape(stro.sum.exptParams.quest_ranges, 3, 3)';
colorRanges(sum(colorRanges, 2)==0, :) = [];
spatialPeriods = sort(unique(stro.trial(:,gaborLambdaInd)));
sfs = stro.sum.exptParams.pixperdeg ./ spatialPeriods;
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
        if (nTrials < 6)
            disp('got here')
            thresholds(clr, spfr) = nan;
            continue;
        end
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
        [aSSE, bSSE, gSSE] = weibullFit(CCtmp, stro.trial(tList, correctInd), 'sse' ,[questAlpha(clr, spfr), 1]);
        correctByContrast = stro.trial(tList, correctInd);
        wrongByContrast = ~correctByContrast;
        if (sum(correctByContrast) < 3 | sum(wrongByContrast) < 3)
            disp('Got here 2');
            thresholds(clr, spfr) = nan;
            continue;
        end
        [weibullAlpha(clr, spfr), beta(clr, spfr), gamma(clr, spfr), success] = weibullFit(CCtmp, [correctByContrast(:), wrongByContrast(:)], 'mle', [aSSE, bSSE]);
        if (weibullAlpha(clr, spfr) > 100)
            disp('Yuck - crazy threshold');
           % keyboard
        end
        if success <= 0;
            weibullAlpha(clr, spfr) = NaN;
        end
        
        
        %now assign the output appropriately
        switch lower(fitMeth)
            case 'mode'
                thresholds(clr, spfr) = questAlpha(clr, spfr);
            case 'weibull'
                thresholds(clr, spfr) = weibullAlpha(clr, spfr);
        end
        QuestTrajectories{clr,spfr} = threshByTrial;
    end
end
    

