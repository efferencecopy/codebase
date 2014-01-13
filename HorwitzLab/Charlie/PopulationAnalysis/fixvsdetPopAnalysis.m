%% Fixation vs. Detection Population analysis
%

%first just compile the data
analParams.cellNum = 1;
analParams.start = 0; %gabor onset
analParams.end = 0.666; %gabor offset
analParams.lowCutoff = 5; %min number of trials to compute cp
analParams.meth = 2; %rate
[fnames, spikenums] = fnamesFromTxt2('N:\NexFiles\Charlie\Kali\text files\DTvsGTpopList.txt');
out.fnames = fnames;
for ex = 1:length(fnames)
    disp(fnames{ex})
    GT = gtobj(fnames{ex}{1});
    out.dat{ex}.gt{1} = getGratingTuning(GT);
    DT = nex2stro(findfile(fnames{ex}{2}));
    [DT. fixTrials] = stripOutGratingTrials(DT);
    [out.dat{ex}.dt.m, out.dat{ex}.dt.c, out.dat{ex}.dt.expt] = DTunpack(DT,analParams);
    if length(fnames{ex})==3
        GT = gtobj(fnames{ex}{3});
        out.dat{ex}.gt{2} = getGratingTuning(GT);
    end
end

%save the compiled data
saveDir = 'N:\NexFiles\Charlie\Kali\batchData';
presDir = pwd;
str = ['fixvdet_batch_' date];
cd(saveDir)
eval(['save ' str ' out']);
fprintf('**** Batch data saved to directory: %s ****\n', saveDir);
cd(presDir)

%% Neurmetric functions from CRFs in GT and DT


%compile NTs
load('N:\NexFiles\Charlie\Kali\batchData\fixvdet_batch_06-Jun-2010.mat')
nExpts = length(out.dat);
for ex = 1:nExpts
    nGTfiles = length(out.dat{ex}.gt);
    for a = 1:nGTfiles;
        nColors = size(out.dat{ex}.gt{a}.crf.colors,1);
        for clr = 1:nColors
            nContrasts = length(out.dat{ex}.gt{a}.crf.rates{clr});
            for cntrst = 1:nContrasts;
                noise = out.dat{ex}.gt{a}.crf.rates{clr}{1};
                signal = out.dat{ex}.gt{a}.crf.rates{clr}{cntrst};
                out.dat{ex}.gt{a}.crf.roc{clr,1}(cntrst) = roc(noise, signal);
            end
            
            %fit a neurometric function to the roc areas
            fprintf('ex: %d, clr: %d \n', ex, clr);
            errs = abs(0.82-out.dat{ex}.gt{a}.crf.roc{clr,1});
            aGuess = out.dat{ex}.gt{a}.crf.norms{clr,1}(find(errs == min(errs), 1, 'last'));
            nTrialsByCntrst = cellfun(@length, out.dat{ex}.gt{a}.crf.rates{clr,1});
            correctByContrast = out.dat{ex}.gt{a}.crf.roc{clr,1}.*nTrialsByCntrst;
            wrongByContrast = (nTrialsByCntrst - correctByContrast);
            [out.dat{ex}.gt{a}.crf.alpha(clr,1), out.dat{ex}.gt{a}.crf.beta(clr,1), out.dat{ex}.gt{a}.crf.gamma(clr,1), success(2), modErrs] = weibullFit(out.dat{ex}.gt{a}.crf.norms{clr,1}, [correctByContrast(:), wrongByContrast(:)], 'mle', [aGuess, 2]);
            out.dat{ex}.gt{a}.crf.errs.alpha(clr,1) = modErrs(1);
        end
    end
end


%iterate over the expts and plot the data
plotSE = 0;
nExpts = length(out.dat);
for ex = 1:nExpts
    figure;
   % set(gcf, 'position', [-1100, 22, 884, 695])
    %L-M neurofuns from GT
    gNorms = {};
    gX = {};
    gAlpha = [];
    gsemAlpha = [];
    gMod = {};
    gRoc = {};
    gRates = {};
    gSEM = {};
    for a = 1:length(out.dat{ex}.gt)
        lvmIdx = ismember(out.dat{ex}.gt{a}.crf.colors, [1 -1 0], 'rows');
        gNorms{a} = out.dat{ex}.gt{a}.crf.norms{lvmIdx};
        gX{a} = linspace(gNorms{a}(2).*.98, gNorms{a}(end).*1.2, 500);
        gAlpha(a) = out.dat{ex}.gt{a}.crf.alpha(lvmIdx,1);
        gsemAlpha(a) = out.dat{ex}.gt{a}.crf.errs.alpha(lvmIdx,1);
        gBeta = out.dat{ex}.gt{a}.crf.beta(lvmIdx,1);
        gGamma = out.dat{ex}.gt{a}.crf.gamma(lvmIdx,1);
        gMod{a} = gGamma - (gGamma - 0.5).*exp(-((gX{a}./gAlpha(a)).^gBeta));
        gRoc{a} = out.dat{ex}.gt{a}.crf.roc{lvmIdx,1};
        gRates{a} = cellfun(@mean, out.dat{ex}.gt{a}.crf.rates{lvmIdx});
        gSigma = cellfun(@std, out.dat{ex}.gt{a}.crf.rates{lvmIdx});
        nTrials = cellfun(@length, out.dat{ex}.gt{a}.crf.rates{lvmIdx});
        gSEM{a} = gSigma ./ sqrt(nTrials);
    end
    
    %L-M neurofuns from DT
    lvmIdx = ismember(sign(out.dat{ex}.dt.expt.standColors), [1 -1 0], 'rows');
    dNorms = out.dat{ex}.dt.expt.norms{lvmIdx,1};
    dX = linspace(dNorms(2).*.98, dNorms(end), 500);
    dAlpha = out.dat{ex}.dt.c.alpha(lvmIdx,1);
    dsemAlpha = out.dat{ex}.dt.c.err.mle.alpha(lvmIdx,1);
    dBeta = out.dat{ex}.dt.c.beta(lvmIdx,1);
    dGamma = out.dat{ex}.dt.c.gamma(lvmIdx,1);
    dMod = dGamma - (dGamma - 0.5).*exp(-((dX./dAlpha).^dBeta));
    dRoc = out.dat{ex}.dt.c.roc{lvmIdx};
    dRates = cellfun(@mean, out.dat{ex}.dt.c.crfIn{lvmIdx});
    dSigma = cellfun(@std, out.dat{ex}.dt.c.crfIn{lvmIdx});
    nTrials = cellfun(@length, out.dat{ex}.dt.c.crfIn{lvmIdx});
    dSEM = dSigma ./ sqrt(nTrials);
    
    subplot(2,2,1) %neurofuns for l-m
    hold on,
    pltClr = {'r', 'm'};
    for a = 1:length(out.dat{ex}.gt)
        plot(gX{a}, gMod{a}, pltClr{a});
        plot(gNorms{a}, gRoc{a}, [pltClr{a}, 'o']);
    end
    plot(dX, dMod, 'k');
    plot(dNorms, dRoc, 'ko')
    if plotSE
        for a = 1:length(out.dat{ex}.gt);
            plot([gAlpha(a)-gsemAlpha(a), gAlpha(a)+gsemAlpha(a)], [0.82, 0.82], [pltClr{a}, '.-']);
        end
        plot([dAlpha-dsemAlpha, dAlpha+dsemAlpha], [0.82, 0.82], 'k.-');
    end
    set(gca, 'xscale', 'log')
    ylabel('p(Correct)')
    xlabel('cone contrast')
    title('neurmetric functions')
    axis tight
    ylim([0.4 1.1])
    hold off
    
    subplot(2,2,3) %mean rates for l-m
    hold on,
    for a = 1:length(out.dat{ex}.gt)
        errorbar(gNorms{a}, gRates{a}, gSEM{a}, pltClr{a})
    end
    errorbar(dNorms, dRates, dSEM, 'k')
    set(gca, 'xscale', 'log')
    xlabel('cone contrast')
    ylabel('mean rate')
    title('contrast response fxn')
    axis tight
    hold off
    
    %S-iso neurofuns from GT
    gNorms = {};
    gX = {};
    gAlpha = [];
    gsemAlpha = [];
    gMod = {};
    gRoc = {};
    gRates = {};
    gSEM = {};
    for a = 1:length(out.dat{ex}.gt)
        sIdx = ismember(out.dat{ex}.gt{a}.crf.colors, [0 0 1], 'rows');
        gNorms{a} = out.dat{ex}.gt{a}.crf.norms{sIdx};
        gX{a} = linspace(gNorms{a}(2).*.98, gNorms{a}(end).*1.2, 500);
        gAlpha(a) = out.dat{ex}.gt{a}.crf.alpha(sIdx,1);
        gsemAlpha(a) = out.dat{ex}.gt{a}.crf.errs.alpha(sIdx,1);
        gBeta = out.dat{ex}.gt{a}.crf.beta(sIdx,1);
        gGamma = out.dat{ex}.gt{a}.crf.gamma(sIdx,1);
        gMod{a} = gGamma - (gGamma - 0.5).*exp(-((gX{a}./gAlpha(a)).^gBeta));
        gRoc{a} = out.dat{ex}.gt{a}.crf.roc{sIdx,1};
        gRates{a} = cellfun(@mean, out.dat{ex}.gt{a}.crf.rates{sIdx});
        gSigma = cellfun(@std, out.dat{ex}.gt{a}.crf.rates{sIdx});
        nTrials = cellfun(@length, out.dat{ex}.gt{a}.crf.rates{sIdx});
        gSEM{a} = gSigma ./ sqrt(nTrials);
    end
    
    %S-iso neurofuns from DT
    sIdx = ismember(sign(out.dat{ex}.dt.expt.standColors), [0 0 1], 'rows');
    dNorms = out.dat{ex}.dt.expt.norms{sIdx,1};
    dX = linspace(dNorms(2).*.98, dNorms(end), 500);
    dAlpha = out.dat{ex}.dt.c.alpha(sIdx,1);
    dsemAlpha = out.dat{ex}.dt.c.err.mle.alpha(sIdx,1);
    dBeta = out.dat{ex}.dt.c.beta(sIdx,1);
    dGamma = out.dat{ex}.dt.c.gamma(sIdx,1);
    dMod = dGamma - (dGamma - 0.5).*exp(-((dX./dAlpha).^dBeta));
    dRoc = out.dat{ex}.dt.c.roc{sIdx};
    dRates = cellfun(@mean, out.dat{ex}.dt.c.crfIn{sIdx});
    dSigma = cellfun(@std, out.dat{ex}.dt.c.crfIn{sIdx});
    nTrials = cellfun(@length, out.dat{ex}.dt.c.crfIn{sIdx});
    dSEM = dSigma ./ sqrt(nTrials);
    
    subplot(2,2,2) %neurofuns for S-iso
    hold on,
    pltClr = {'b', 'c'};
    for a = 1:length(out.dat{ex}.gt)
        plot(gX{a}, gMod{a}, pltClr{a});
        plot(gNorms{a}, gRoc{a}, [pltClr{a}, 'o']);
    end
    plot(dX, dMod, 'k');
    plot(dNorms, dRoc, 'ko')
    if plotSE
        for a = 1:length(out.dat{ex}.gt)
            plot([gAlpha(a)-gsemAlpha(a), gAlpha(a)+gsemAlpha(a)], [0.82, 0.82], [pltClr{a}, '.-']);
        end
        plot([dAlpha(a)-dsemAlpha(a), dAlpha(a)+dsemAlpha(a)], [0.82, 0.82], 'k.-');
    end
    set(gca, 'xscale', 'log')
    ylabel('p(Correct)')
    xlabel('cone contrast')
    title('neurmetric functions')
    axis tight
    ylim([0.4 1.1])
    hold off
    
    subplot(2,2,4) %mean rates for S-iso
    hold on,
    for a = 1:length(out.dat{ex}.gt)
        errorbar(gNorms{a}, gRates{a}, gSEM{a}, pltClr{a})
    end
    errorbar(dNorms, dRates, dSEM, 'k')
    set(gca, 'xscale', 'log')
    xlabel('cone contrast')
    ylabel('mean rate')
    title('contrast response fxn')
    axis tight
    hold off
    
    
    %labels some specifics about each cell.
    gColors = out.dat{ex}.gt{1}.color.colors;
    projOnLum = gColors * [1;1;0];
    isolumIdx = projOnLum == 0;
    lumIdx = ismember(sign(out.dat{ex}.gt{1}.color.colors), [1 1 0], 'rows');
    isolumRates = out.dat{ex}.gt{1}.color.colresp(isolumIdx, 1);
    lumRate = out.dat{ex}.gt{1}.color.colresp(lumIdx, 1);
    csi = max(isolumRates) ./ lumRate;
    sf = out.dat{ex}.gt{1}.sf.prefSF;
    figName = sprintf('%s, csi: %.3f, sf: %.3f', out.fnames{ex}{1}, csi, sf);
    set(gcf, 'numbertitle', 'off', 'name', figName)
end

%% trying to implement poisson regression (for now on the S-iso cases)


[fnames, spikenums] = fnamesFromTxt2('N:\NexFiles\Charlie\Kali\text files\DTvsGTpopList.txt');
pvals = nan(2,4,length(fnames));
for a = 1:length(fnames)
    figure
    set(gcf, 'position', [142   246   915   476]);
    %open up the GT data
    gt = gtobj(fnames{a}{1});
    stimon = mat2cell(gt.trial(:,gt.idx.stimon), ones(size(gt.trial,1),1));
    stimoff = mat2cell(gt.trial(:, gt.idx.stimoff), ones(size(gt.trial,1),1));
    gtcounts = cellfun(@(x,y,z) length(x((x>y)&(x<z))), gt.ras(:, gt.idx.spikes), stimon, stimoff);
    LMS = [gt.trial(:, gt.idx.lcc), gt.trial(:, gt.idx.mcc), gt.trial(:, gt.idx.scc)];
    
    %open up the DT data
    dt = dtobj(fnames{a}{2});
    stimon = mat2cell(dt.trial(:,dt.idx.flashOn), ones(size(dt.trial,1),1));
    stimoff = mat2cell(dt.trial(:, dt.idx.flashOff), ones(size(dt.trial,1),1));
    dtcounts = cellfun(@(x,y,z) length(x((x>y)&(x<z))), dt.ras(:, dt.idx.spikes), stimon, stimoff);
    dtcolors = reshape(dt.sum.exptParams.RF_colors, 3, 3)';
    dtcolors((dtcolors(:,1)<0),:) = dtcolors((dtcolors(:,1)<0),:)*-1;
    nolm = sum(abs(dtcolors(:,1:2)),2)==0;
    negS = dtcolors(:,3)<0;
    dtcolors((nolm&negS),:) = dtcolors((nolm&negS),:) * -1;
    l_zerocntrst = sum(dt.LMS,2) == 0;
    rfx = dt.sum.exptParams.rf_x;
    rfy = dt.sum.exptParams.rf_y;
    l_inRF = (dt.trial(:,dt.idx.flashX) == rfx) & (dt.trial(:,dt.idx.flashY) == rfy);
    
    colors = [0 0 1; 1 -1 0];
    for clr = 1:2
        % initialize the stuff for glmfit
        x = []; %the predictors: <task, contrast, task*contrast>
        y = []; %the dependent variable: <counts>
        
        %deal with GT
        l_p7s = gt.getTrialList(7, colors(clr,:));
        l_bkgnd = gt.getTrialList(7, [0 0 0]);
        l_gtTrials = l_p7s | l_bkgnd;
        gT = zeros(sum(l_gtTrials), 1);
        gC = LMS(l_gtTrials, :);
        gC = sqrt(sum(gC.^2, 2));
        gTC = gT.*gC;
        x = [x; [gT, gC, gTC]];
        y = [y; gtcounts(l_gtTrials)];
        
        %now for DT.... a bigger pain.
        clrIdx = find(ismember(sign(dtcolors), colors(clr,:), 'rows'));
        l_color = dt.trial(:,dt.idx.colorDir) == clrIdx;
        l_dtTrials = (l_color | l_zerocntrst) & l_inRF;
        avgNTrials = sum(l_dtTrials) ./ 8;
        if avgNTrials < 8;
            disp('too few trials')
            continue
        end
        dT = ones(sum(l_dtTrials), 1);
        dC = sqrt(sum(dt.LMS(l_dtTrials, :).^2, 2));
        dTC = dT.*dC;
        x = [x;[dT,dC,dTC]];
        y = [y; dtcounts(l_dtTrials)];
        
        %now for the regression... cross your fingers.
        [betas, dev, stats] = glmfit(x, y, 'poisson', 'link', 'log');
        yhat = glmval(betas, x, 'log');
        pvals(clr,:,a) = stats.p;
        
        
        %plot some stuff
        try
            plotColors = {'b', 'r'};
            subplot(1,2,clr)
            hold on,
            scatter(gC, gtcounts(l_gtTrials),[plotColors{clr},'+'])
            scatter(dC, dtcounts(l_dtTrials), 'k.')
            l_gt = x(:,1) == 0;
            gtYhat = unique(yhat(l_gt));
            gtnorms = unique(x(l_gt,2));
            plot(gtnorms, gtYhat, plotColors{clr})
            l_dt = x(:,1) == 1;
            dtYhat = unique(yhat(l_dt));
            dtnorms = unique(x(l_dt,2));
            plot(dtnorms, dtYhat, 'k--');
            ylabel('counts')
            xlabel('contrast')
            legend('Gratings', 'Detection', 'location', 'northwest')
            title(sprintf('Poisson Regression \n pvals = [offset: %.3f, Task: %.3f, Contrast: %.3f, Interact: %.3f]', stats.p(1), stats.p(2),stats.p(3),stats.p(4)))
            set(gcf, 'name', gt.name, 'numbertitle', 'off')
        catch
            wtf
            gt.name
        end
    end
end





