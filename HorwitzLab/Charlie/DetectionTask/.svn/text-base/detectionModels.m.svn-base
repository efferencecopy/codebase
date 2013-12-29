
fin
% designate a text file for the batch analysis:
CARDVSINT = 'CardVsInt_06-May-2011.mat';


global cardVsIntBatchPath
if isunix
    prefix = '/Users/charliehass/LabStuff/NexFiles/Charlie/Batch Data And Text Files/';
elseif ispc
    prefix = 'N:\NexFiles\Charlie\Batch Data And Text Files\';
end
cardVsIntBatchPath = [prefix, CARDVSINT];
preprocessDTbatchData

% constants
MAKEPLOT = 0;
MEANTOCCSMETH = 'exponential';
MEANTOVARMETH = 'ols';
COUNTSDISTRIBUTION = 'normal'; %could be: normal, poisson
POOLINGMETHOD = 'sum'; % could be: average, sum
SPTFREQCUTOFF = 2;
MAXPOOLSIZE = 602;
possiblePoolSizes = [[1,3,5,7,10], [11:10:101], [102:100:MAXPOOLSIZE]];
NUMPOOLS = 3000; %only consider this many pools per color/poolsize condition
EXCLUDENANS = 0; %should the analysis consider insensitive neurons?
RHO = 0.10;
NUMTRIALS = 500;

% Goal: determine, on average, how many neurons would be needed to
% recapitulate the behavioral detection thresholds using the neurons tested
colorDirs = [1 -1 -1;...
    1 -1 1;...
    1 -1 0;...
    0 0 1];
countsAtThresh = nan(length(out.dat), 4);
sigmaAtThresh = nan(length(out.dat),4);
bkgndCounts = nan(length(out.dat),4);
tic
for clr = 1:4
    fprintf('Color: <%d>\n', clr)
    %initialize the pool
    validUnits = [];

    % STEP ONE: select the neurons that will make up the pool
    for ex = 1:length(out.dat)
        clridx = ismember(sign(out.dat(ex).expt.standColors), colorDirs(clr,:), 'rows');
        if sum(clridx)~=1
            continue %bail if this color isn't represented, or if for some dumb reason it's represented twice (which is likely a problem)
        elseif out.dat(ex).grating.sf.prefSF > SPTFREQCUTOFF;
            continue
        elseif commonExclusions(ex)
            continue
        end
        
        if EXCLUDENANS && isnan(out.dat(ex).c.alpha(clridx))
            continue
        end
            

        %now extract the contrast response function
        n = numel(validUnits)+1;
        validUnits(n).crf = out.dat(ex).c.crfIn{clridx};
        validUnits(n).norms = out.dat(ex).expt.norms{clridx};
        validUnits(n).psychThresh = out.dat(ex).m.alpha(clridx);
        validUnits(n).neuroThresh = out.dat(ex).c.alpha(clridx);
        validUnits(n).roc = out.dat(ex).c.roc{clridx};
        validUnits(n).prefSF = out.dat(ex).expt.sfs;
        validUnits(n).unitNum = ex;
    end
    fprintf('Finished compiling the data\n')
    fprintf('Population has %d neurons. %d (%.1f%%) are insensitive\n', numel(validUnits), sum(isnan([validUnits.neuroThresh])), (sum(isnan([validUnits.neuroThresh]))./numel(validUnits))*100)


    %STEP TWO: Fit the mean/contrast and mean/var relationships
    for n = 1:length(validUnits) %redefining "n"

        %first the mean to intensity relationship via poisson regression
        counts = [];
        ccs = [];
        for i = 1:length(validUnits(n).norms)
            tmp = round(validUnits(n).crf{i}' .* 0.666); %converting from rates to counts
            counts = [counts; tmp];
            ccs = [ccs; ones(numel(validUnits(n).crf{i}),1).*validUnits(n).norms(i)];
        end

        switch MEANTOCCSMETH
            case 'exponential'
                %not sure this method adequately captures an "additive"
                %baseline...
                validUnits(n).crfBeta = glmfit(ccs, counts, 'poisson', 'link', 'log');
        end

        %now the mean to variance relationship
        switch MEANTOVARMETH
            case 'ols'
                xbar = cellfun(@(x,y) mean(round(x.*y)), validUnits(n).crf, repmat({0.666}, size(validUnits(n).crf)));
                svar = cellfun(@(x,y) var(round(x.*y)), validUnits(n).crf, repmat({0.666}, size(validUnits(n).crf)));
                zeroMean = xbar == 0;
                xbar(zeroMean) = [];
                svar(zeroMean) = [];
                validUnits(n).meanVarSlope = xbar(:) \ svar(:);
        end
        
        %estimate the average number of spikes and the SD at behaviroal
        %detection threshold.
        countsAtThresh(n, clr) = exp(validUnits(n).crfBeta(1) + validUnits(n).crfBeta(2).*validUnits(n).psychThresh);
        sigmaAtThresh(n, clr) = sqrt(validUnits(n).meanVarSlope .* countsAtThresh(n,clr));
        bkgndCounts(n,clr) = exp(validUnits(n).crfBeta(1));

        if MAKEPLOT
            figure
            subplot(1,2,1)
            cla, hold on,
            plot(ccs, counts, 'k.')
            plot(ccs, glmval(validUnits(n).crfBeta, ccs, 'log'))
            subplot(1,2,2)
            cla, hold on,
            plot(xbar, svar, 'k.')
            plot([0;xbar(:)], validUnits(n).meanVarSlope.*[0;xbar(:)], 'k')
            xlabel('Mean Counts')
            ylabel('Var Counts')
            hold off
        end
    end
    fprintf('Done fitting CRFs and Mean vs Var\n')

    %STEP THREE: run the population model
    avgBehavioralThresh(clr) = mean([validUnits.psychThresh]);
    allNorms = [validUnits.norms];
    allNorms(allNorms==0) = [];
    synthNorms = [0, logspace(log10(avgBehavioralThresh(clr)*.3), log10(avgBehavioralThresh(clr)*1.75), 7)];
    for p = 1:numel(possiblePoolSizes)
        poolSize = possiblePoolSizes(p);
        fprintf('Pool size %d of %d\n', p, numel(possiblePoolSizes))

        neuronPools= unidrnd(numel(validUnits), NUMPOOLS, poolSize); %selecting from all possible pool types with replacement
        pooledThreshold = nan(size(neuronPools,1),1);
        for a = 1:size(neuronPools,1) %iterating over pools

            %cycle through each neuron in the pool and draw from their
            %(parameterized) CRFs to create new distributions
            switch COUNTSDISTRIBUTION
                case 'poisson'
                    synthUnitCounts = nan(NUMTRIALS, 8, poolSize); %initialize the synthetic pool
                case 'normal'
                    synthUnitMeans = nan(poolSize, 8);
                    synthUnitVars = nan(poolSize, 8);
            end

            for i = 1:size(neuronPools,2) %iterating over neurons in pool

                unitIdx = neuronPools(a,i);

                switch MEANTOCCSMETH
                    case 'exponential'
                        crfBetas = validUnits(unitIdx).crfBeta;
                        unitMeans = exp(crfBetas(1) + (crfBetas(2).*synthNorms));
                        if any(unitMeans>500); disp('high rates'); keyboard; end %doing a sanity check to make sure improper extrapolation of the CRF doesn't occur
                end

                switch MEANTOVARMETH
                    case 'ols'
                        slope = validUnits(unitIdx).meanVarSlope;
                        unitVars = slope.*unitMeans;
                end

                %compile the population counts, or means and vars....
                switch COUNTSDISTRIBUTION
                    case 'normal'
                        synthUnitMeans(i,:) = unitMeans; %store them here so that I can analytically derive the pooled mean and var later
                        synthUnitVars(i,:) = unitVars;
                    case 'poisson'
                        tmpMeans = repmat(unitMeans, NUMTRIALS, 1); %store the counts for a monte-carlo style simulation later
                        synthPoolCounts(:,:,i) = poissrnd(tmpMeans);
                end
            end

            %pool the responses acording to a specified rule.
            %"pooledResponse" should have deminsions [nTrials x nContrasts]
            switch POOLINGMETHOD
                case 'average'
                    switch COUNTSDISTRIBUTION
                        case 'normal'
                            pooledMean = mean(synthUnitMeans,1);
                            pooledVar = pooledVarsByContrast(synthUnitVars, RHO);
                            pooledVar = pooledVar ./ poolSize^2; %the output of pooledVarsByContast assumes 'sum', so compinsate for pooling via 'mean'
                            pooledStd = sqrt(pooledVar);
                        case 'poisson'
                            pooledResponse = mean(synthPoolCounts,3);
                    end
                case 'sum'
                    switch COUNTSDISTRIBUTION
                        case 'normal'
                            pooledMean = sum(synthUnitMeans,1); %could be wrong (should be sum(...))
                            pooledVar = pooledVarsByContrast(synthUnitVars, RHO);
                            pooledStd = sqrt(pooledVar);
                        case 'poisson'
                            pooledResponse = sum(synthPoolCounts,3);
                    end
            end

            %run the ROC analysis
            switch COUNTSDISTRIBUTION
                case 'normal'
                    %find the domian over which to calculate the noise dist
                    if pooledStd(1) == 0
                        noiseVals = [pooledMean(1), pooledMean(1)]; %just consider the region around the mean
                    else
                        noiseVals = norminv([0.001, 0.999], pooledMean(1), pooledStd(1));
                    end
                    
                    %iterate over contrasts
                    for j = 1:size(pooledMean,2);
                        if pooledStd(j) == 0
                            sigVals = [pooledMean(j), pooledMean(j)];
                        else
                            sigVals = norminv([0.001,0.999], pooledMean(j), pooledStd(j));
                        end
                        
                        lowVal = min([noiseVals(1), sigVals(1)]);
                        highVal = max([noiseVals(2), sigVals(2)]);
                        x = linspace(lowVal, highVal, 500);
                        pFA = [1, 1-normcdf(x, pooledMean(1), pooledStd(1)), 0]; %manually add the p=0 and p=1 condition
                        pHit = [1, 1 - normcdf(x, pooledMean(j), pooledStd(j)), 0];
                        pooledROC(j) = -trapz(pFA, pHit);
                    end
                case 'poisson'
                    for j = 1:size(pooledResponse,2);
                        pooledROC(j) = roc(pooledResponse(:,1), pooledResponse(:,j));
                    end
            end
            
            %some basic error checking
            if any(isnan(pooledROC)) || (numel(pooledROC)<8)
                keyboard
            end

            %fit the neurometric function with a cumulative weibull
            l_zero = synthNorms==0;
            fitNorms = synthNorms(~l_zero);
            fitROC = pooledROC(~l_zero);
            errs = abs(0.82-fitROC);
            aGuess = fitNorms(find(errs == min(errs), 1, 'last'));
            correctByContrast = (fitROC .* NUMTRIALS);
            wrongByContrast = (NUMTRIALS - correctByContrast);
            [aSSE, bSSE] = weibullFit(fitNorms, fitROC, 'sse', [aGuess 1]);
            [alpha, ~, success] = weibullFit(fitNorms, [correctByContrast(:), wrongByContrast(:)], 'mle', [aSSE, bSSE]);

            if success
                pooledThreshold(a) = alpha;
            end
        end

        avgThreshByPoolSize(clr, p) = mean(pooledThreshold);
        semThreshByPoolSize(clr, p) = std(pooledThreshold) ./ sqrt(length(pooledThreshold));
        normalizedThreshold = pooledThreshold./avgBehavioralThresh(clr);
        avgNormThreshByPoolSize(clr, p) = mean(normalizedThreshold);
        semNormThreshByPoolSize(clr, p) = std(normalizedThreshold) ./ sqrt(length(normalizedThreshold));
    end
end
toc

figure, hold on,
pltClr = {'g', 'm', 'r', 'b'};
PLOTTYPE = 'normalized';
nPlots = size(avgThreshByPoolSize,1);
for clr = 1:nPlots
    switch PLOTTYPE
        case 'normalized'
            errorbar(possiblePoolSizes, avgNormThreshByPoolSize(clr,:), semNormThreshByPoolSize(clr,:), ['-', pltClr{clr},'.'], 'linewidth', 2);
            plot([0.5, (max(possiblePoolSizes)+0.5)], [1;1], 'k:', 'linewidth', 2);
        case 'regular'
            errorbar(possiblePoolSizes, avgThreshByPoolSize(clr,:), semThreshByPoolSize(clr,:), ['-', pltClr{clr},'.'], 'linewidth', 2);
            plot([0.5, (max(possiblePoolSizes+0.5))], [avgBehavioralThresh(clr); avgBehavioralThresh(clr)], 'k:', 'linewidth', 2);
    end
end
set(gca, 'xscale', 'log', 'yscale', 'log')
xlabel('Pool Size')
ylabel('Pooled Threshold')
title(sprintf('Rho = %.3f, Including Insensitive: %d', RHO, ~EXCLUDENANS))
axis tight

% I've added new functionality for insensitive neurons. Some of them have
% std = 0. This is a problem for norminv and normpdf, but not for normcdf.
% I use norminv to figure out the correct domain of analysis during the
% computation of area under roc. Instead of letting norminv bonking, I flag
% these cases and manually adjust the domain.

%Also, some of the CRFs shoot off to huge values if you extrapolate too
%much on the contrast axis. I'm now looking for cases where the fitted
%firing rate exceedes 500sp/sec (just as a sanity check) and I'm only
%considering contrasts that are 0.3 to 1.75x threshold.
