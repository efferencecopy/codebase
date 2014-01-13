function [thresholds, colorDirs, sfs] = habitUnpack(stro, fitMeth, winsize, perfrange)
%  UNPACKS HABITUATION EXPERIMENTS
%
%   EXAMPLE: [thresh, colors, sfs] = habitUnpack(stro)
%

if ~exist('winsize', 'var')
    winsize = 20;
end
if ~exist('perfrange', 'var')
    perfrange = [];
end


%first some preliminaries
habitIndicies
nColors = length(unique(stro.trial(:, colorDirInd)));
spatialPeriods = sort(unique(stro.trial(:,gaborLambdaInd)));
nSfs = length(unique(stro.trial(:, gaborLambdaInd)));
sfs = stro.sum.exptParams.pixperdeg ./ spatialPeriods;
colorDirs = reshape(stro.sum.exptParams.gaborColor, 3,3)';
emptyColors = sum(abs(colorDirs), 2)==0;
colorDirs(emptyColors, :) = [];
questDomain = reshape(stro.sum.exptParams.questDomain, 3, 3)';
questDomain(emptyColors, :) = [];
stro.sum.exptParams.alphaGuess(emptyColors) = [];
x = 0:255; %the normal domain of the gamma look up table
xx = linspace(0, 255, 2^16); %the desired quantization of the gammaTable
g1 = reshape(stro.sum.exptParams.gammaTable, 256, 3);
gammaTable = [spline(x, g1(:,1), xx)', spline(x, g1(:,2), xx)', spline(x, g1(:,3), xx)'];
bkgndrgb = stro.sum.exptParams.bkgndrgb;
M = reshape(stro.sum.exptParams.Mmtx, 3, 3);
bkgndlms = M * bkgndrgb;

figure
for clr = 1:nColors
    yRange = [0 0]; %for plotting
    for sf = 1:nSfs
        %reconstruct the prior
        sigma = stro.sum.exptParams.questSigma;
        beta = stro.sum.exptParams.questBeta;
        domain = [questDomain(clr, 1):questDomain(clr, 3):questDomain(clr, 2)];
        guess = stro.sum.exptParams.alphaGuess(clr);
        prior = (1./sigma*sqrt(2*pi)) .* exp(-(domain-guess).^2 ./ (2*sigma^2));
        Q = log(prior);
        
        %pull out the appropriate trials
        l_sfs = stro.trial(:,gaborLambdaInd) == spatialPeriods(sf);
        l_clr = stro.trial(:,colorDirInd) == clr;
        tList = find(l_sfs & l_clr);
        nTrials = length(tList);
        
        %reconstruct the trials
        CCtrl = nan(1,nTrials);
        threshByTrial = nan(1, nTrials);
        for trl = 1:nTrials;
            [m, idx] = max(Q);
            threshByTrial(trl) = domain(idx); %the expected CC based on my reconstruction
            RGB = stro.trial(tList(trl),actGaborRInd:actGaborBInd)+1;
            rgb = [gammaTable(RGB(1), 1), gammaTable(RGB(2), 2), gammaTable(RGB(3), 3)];
            LMS = ((M*rgb')-bkgndlms) ./ bkgndlms;
            CCtrl(trl) = norm(LMS).*100; %the actual CC presented
            
            %update the quest function
            liklihood = 1 - 0.5.*exp(-(CCtrl(trl)./domain).^beta);
            if ~stro.trial(tList(trl), correctInd)
                liklihood = 1-liklihood;
            end
            Q = Q + log(liklihood);
        end
        
        %check to make sure that the modal value of the quest function was
        %presented on each trial, then assign a threshold estimate based on
        %the modal value of the quest
        prcntDiff = abs(threshByTrial - CCtrl)./threshByTrial;
        if any(prcntDiff>0.1)
            [CCtrl(:), threshByTrial(:)]
            %error('the modal value of the quest function at least 0.1%CC different than that presented') %#ok<CTPCT>
        end
        [dum, idx] = max(Q);
        questAlpha(clr, sf) = domain(idx);
        
        %estimate threshold a second time based on a fit to a cumulative weibull.
        [aSSE, bSSE, gSSE] = weibullFit(CCtrl, stro.trial(tList, correctInd), 'sse' ,[questAlpha(clr, sf), 85]);
        correctByContrast = stro.trial(tList, correctInd);
        wrongByContrast = ~correctByContrast;
        [weibullAlpha(clr, sf), beta, gamma, success] = weibullFit(CCtrl, [correctByContrast(:), wrongByContrast(:)], 'mle', [aSSE, bSSE]);
        if success <= 0;
            weibullAlpha(clr, sf) = NaN;
        end
        
        %plot the raw data
        filter = ones(1,winsize)./winsize;
        perf = conv(stro.trial(tList, correctInd), filter);
        perf(end-(winsize-2):end) = []; %convolution is zero padded on the tails
        perf(1:(winsize-1)) = nan; %these values don't reflect a true running avg
        subplot(size(colorDirs, 1), length(sfs), (clr-1).*length(sfs)+sf)
        hands = plotyy([1:nTrials], threshByTrial, [1:nTrials], perf);
        unit = LMS./norm(LMS);
        title(sprintf('[%.3f, %.3f, %.3f]', unit(1), unit(2), unit(3)));
        set(get(min(hands), 'children'), 'color', 'k', 'marker', 'none', 'linewidth', 3);
        set(get(max(hands), 'children'), 'color', 'b', 'marker', '.', 'linewidth', 0.25);
        nYlabels = 5;
        set(max(hands), 'yLim', [0,1], 'yTick', linspace(0,1,nYlabels), 'ytickLabel', linspace(0,1,nYlabels));
        Y1lims = get(min(hands), 'ylim');
        set(min(hands),'yTick', linspace(Y1lims(1),Y1lims(2),nYlabels), 'ytickLabel', linspace(Y1lims(1),Y1lims(2),nYlabels));
        set(hands, 'ycolor', 'k')
               
        %filter out runs where the subjects performace was outside some
        %predetermined range
        if ~isempty(perfrange)
            if (perf(end) < (0.82-perfrange)) || (perf(end) > (0.82+perfrange))
                questAlpha(clr, sf) = nan;
                weibullAlpha(clr, sf) = nan;
            end
        end
        
        %now assign the output appropriately
        switch lower(fitMeth)
            case 'mode'
                thresholds(clr, sf) = questAlpha(clr, sf);
            case 'weibull'
                thresholds(clr, sf) = weibullAlpha(clr, sf);
        end
    end %spt freq
end %colors