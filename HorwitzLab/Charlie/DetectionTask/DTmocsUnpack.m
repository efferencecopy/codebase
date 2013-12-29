function [monkFits, cellFits, expt] = DTmocsUnpack(stro, defaultParams)

%allow a calling function to suppress the ui prompt window
if nargin<2
    defaultParams = 0;
end

% This code isn't backwards compatible, so make sure things are o.k.
datecheck(stro);

%remove any grating catch trials that may be present
stro = stripOutGratingTrials(stro);

% create the relavant indicies and recompute the monitorcalibration data
bkgndrgb = [stro.sum.exptParams.bkgnd_r, stro.sum.exptParams.bkgnd_g, stro.sum.exptParams.bkgnd_b];
M = reshape(stro.sum.exptParams.m_mtx, 3, 3);
bkgndLMS = M * bkgndrgb';
x = 0:255; %the normal range of the gamma look up table
xx = linspace(0, 255, 2^16); %the desired quantization of the gammaTable
g1 = reshape(stro.sum.exptParams.gamma_table, 256, 3);
gammaTable = [spline(x, g1(:,1), xx)', spline(x, g1(:,2), xx)', spline(x, g1(:,3), xx)'];


% determine some relavent experimental characteristics
numTrials = size(stro.trial, 1);
sptPeriods = unique(stro.trial(:,stro.idx.gaborLambda)); %in pixels
expt.sfs = 1./(sptPeriods./stro.sum.exptParams.pixperdeg);
nSptFreqs = size(sptPeriods, 1);
nColors = max(stro.trial(:,stro.idx.colorDir));
contrasts = unique(stro.trial(:, stro.idx.cntrstLev));
nContrasts = max(contrasts);
nCellChannels = length(strmatch('sig00', strvcat(stro.sum.rasterCells{:})));
if numel(stro.sum.rasterCells) > 0
    for a = 1:length(stro.sum.rasterCells)
        wf_channels(a) = strcmpi(stro.sum.rasterCells{a}(end-2:end), '_wf');
    end
    nCellChannels = nCellChannels-sum(wf_channels); %don't want to count the wf column as a spike channel
    nFigRows = 1 + (nCellChannels>0); %add a figure row for neural data
else
    nCellChannels = 0;
    nFigRows = 1;
end

if nCellChannels>0;
    [cellNum, anlyStart, anlyEnd, lowCutoff, analysisMeth] = getAnalysisParams(nCellChannels, defaultParams);
    colIndsToPossibleCells = strmatch('sig001', stro.sum.rasterCells);
    expt.cellNum = colIndsToPossibleCells(cellNum);
else
    expt.cellNum = 0;
end

% standardize the color directions. This is an important step, and allows
% population analyses to index data in the same way from potentially
% heterogeneous experimental sessions. The convention will be to force the
% L-cone sign to be positive. Also, the S-iso condition will be constrained
% to be positive.
expt.standColors = reshape(stro.sum.exptParams.RF_colors, 3, 3)';
expt.standColors(sum(abs(expt.standColors),2) == 0, :) = []; %remove nonexistent colors
norms = sqrt(sum(expt.standColors.^2, 2));
expt.standColors = expt.standColors ./ repmat(norms, 1, 3);
expt.standColors(expt.standColors(:,1)<0, :) = expt.standColors(expt.standColors(:,1)<0, :) .* -1; %flip the signs on color dirs with -Lcones
l_siso = softEq([0,0,1], abs(expt.standColors), [], 'rows');
if sum(expt.standColors(l_siso,:)) < 0;
    expt.standColors(l_siso,:) = expt.standColors(l_siso,:) .* -1;
end

%convert from DAC gun space to cone contrasts. I have to add one to all the
%RGB values dropped in the data file b/c those values are DAC voltages are
%b/w 0 and maxdac-1... but I need them to be an index b/w 1 and maxdac.
actGunVals = [stro.trial(:,stro.idx.actflashR), stro.trial(:,stro.idx.actflashG), stro.trial(:,stro.idx.actflashB)] + 1;
act_rgb = [gammaTable(actGunVals(:,1), 1), gammaTable(actGunVals(:,2), 2), gammaTable(actGunVals(:,3), 3)];
act_rgb = act_rgb - repmat(bkgndrgb, numTrials, 1);
act_lms = [(M * act_rgb') ./ repmat(bkgndLMS, 1, numTrials)]';

%roll through the different experimental conditions and fit psycho- and
%neurometric functions
monkFits = [];
cellFits = [];
grandZscores{1,2} = [];
for a = 1:nColors
    figure; %open a new figure for each color direction
    colorTList = [stro.trial(:,stro.idx.colorDir) == a];

    for b = 1:nSptFreqs
        sfsTList = [stro.trial(:, stro.idx.gaborLambda) == sptPeriods(b)];

        for c = 1:nContrasts;
            cntrstList = [stro.trial(:, stro.idx.cntrstLev) == c];
            %Currently, the first cntrstLevel is the 'no stimulus'
            %condition. This condition doesn't really have a sfs or color,
            %so ignore these aspects of this condition.
            if (c == 1)
                tList = cntrstList;
            else
                tList = (colorTList & sfsTList & cntrstList);
            end
            t_lms = act_lms(tList,:);
            expt.norms{a,b}(c) = norm(unique(t_lms, 'rows'));
            monkFits.performance{a,b}(c) = sum(stro.trial(tList, stro.idx.correct)) ./ sum(tList);
            monkFits.nTrials{a,b}(c) = sum(tList);
        end %cyc through contrasts

        % now fit psycho- and neurometric functions. Plot as need be.
        monkFits = psychFun(expt.norms, monkFits.nTrials{a,b}, a, b, nFigRows, nSptFreqs, monkFits);
        if expt.cellNum > 0
            conditionTList = (colorTList & sfsTList);
            [cellFits, grandZscores] = neuroFun(stro, expt.norms, expt.cellNum, conditionTList, a, b, nFigRows, nSptFreqs, analysisMeth, anlyStart, anlyEnd, lowCutoff, cellFits, grandZscores);
        end
    end %cyc through sptFreqs
    
    %add some titles to the figures once everything has been calculated for
    %a single color direction
    colorLabels = unique(act_lms(tList,:), 'rows') ./ norm(unique(act_lms(tList,:), 'rows'));
    if numel(colorLabels) ~= 3; error('RGB values may differ within a single contrast'); end
    
    dots = strfind(stro.sum.fileName, '.');
    if ispc
        slash = strfind(stro.sum.fileName, '\');
    elseif ismac
        slash = strfind(stro.sum.fileName, '/');
    end
    set(gcf, 'Name', sprintf('%s ==> [%.2f, %.2f, %.2f]', stro.sum.fileName(slash(end)+1:dots(end)-1), colorLabels(1), colorLabels(2), colorLabels(3)));
    set(gcf,'NumberTitle', 'off');
    for j = 1:nSptFreqs;
        subplot(nFigRows, nSptFreqs, j);
        if ~isempty(cellFits)
            title(sprintf('sptFreq: %.2f cyc/deg \n \\alpha_{p}=%.3f, \\alpha_{n}=%.3f \n \\beta_{p}=%.3f, \\beta_{n}=%.3f', expt.sfs(j), monkFits.alpha(a,j), cellFits.alpha(a,j), monkFits.beta(a,j), cellFits.beta(a,j)));
        else
            title(sprintf('sptFreq: %.2f cyc/deg \n psy: %.3f', expt.sfs(j), monkFits.alpha(a,j)));
        end
    end
end % cyc through colors

%add the CP for the entire cell (pooled across color/contrast/sf)
if (numel(grandZscores{1}) > 0) && (numel(grandZscores{2}) > 0)
    cellFits.cp.grandCP.val = choiceProb(grandZscores, [], lowCutoff);
    [cellFits.cp.grandCP.p,~] = ranksum(grandZscores{1}, grandZscores{2});
else
    cellFits.cp.grandCP.p = nan;
    cellFits.cp.grandCP.val = nan;
end

end %function

                %**************************************%
                %**          SUBFUNCTIONS            **%
                %**************************************%

function monkFits = psychFun(norms, nTrialsByCntrst, a, b, nFigRows, nSptFreqs, monkFits)
    %fits a psychometric function to the choice data and coughs up
    %parameter estimates. Works on a case by case basis where each
    %"case" is a unique sfs/color combination
    
    %fit the psychometric function
    zeroInd = norms{a,b} == 0; %don't consider the zero contrast condition
    tmpNorms = norms{a,b}(~zeroInd);
    errs = abs(0.82-monkFits.performance{a,b}(~zeroInd));
    aGuess = tmpNorms(find(errs == min(errs), 1, 'last'));
    [aSSE, bSSE, ~, success(1)] = weibullFit(norms{a,b}, monkFits.performance{a,b}, 'sse', [aGuess 1]);
    correctByContrast = (monkFits.performance{a,b}.*nTrialsByCntrst);
    wrongByContrast = (nTrialsByCntrst - correctByContrast);
    [monkFits.alpha(a,b), monkFits.beta(a,b), monkFits.gamma(a,b), success(2), modErrs] = weibullFit(norms{a,b}, [correctByContrast(:), wrongByContrast(:)], 'mle', [aSSE, bSSE]);
    dNorm = min(diff(norms{a,b}))./100;
    modelMLE = monkFits.gamma(a,b) + (0.5 - monkFits.gamma(a,b)).*exp(-(([0:dNorm:max(norms{a,b})]./monkFits.alpha(a,b)).^monkFits.beta(a,b)));

    if all(success);
        %determine errors for the parameter estimates if the SSE and MLE
        %searches have been sucessful
        monkFits.err.mle.alpha(a,b) = modErrs(1);
        monkFits.err.mle.beta(a,b) = modErrs(2);
    else
        %nans are used to denote unsucessful fits during subsequent
        %analyses
        monkFits.alpha(a,b) = NaN;
        monkFits.beta(a,b) = NaN;
        monkFits.err.boot.alpha(a,b) = nan;
        monkFits.err.boot.beta(a,b) = nan;
    end
    
    subplot(nFigRows, nSptFreqs, b)
    h1 = semilogx(norms{a,b}, monkFits.performance{a,b}, 'k.');
    hold on
    h2 = semilogx([0:dNorm:max(norms{a,b})], modelMLE, 'k');
    hold off
    hg = hggroup;  % So that legend works
    set(get(get(hg,'Annotation'),'LegendInformation'),'IconDisplayStyle','on')
    set(h1,'Parent',hg);
    set(h2,'Parent',hg);
    xlabel('Cone Contrast');
    ylabel('p(Correct)');
    xlim([(norms{a,b}(2)./3), max(norms{a,b})]); %don't display the zero contrast
    ylim([min(monkFits.performance{a,b})-.02, 1.02]);
    hold off
end %psychFun

function [cellFits, grandZscores] = neuroFun(stro, norms, cellNum, conditionTList, a, b, nFigRows, nSptFreqs, analysisMeth, anlyStart, anlyEnd, lowCutoff, cellFits, grandZscores)
    %gets called for each spatial frequency/color type.
    %some preliminaries and experimental params:
    contrasts = unique(stro.trial(:, stro.idx.cntrstLev));
    nContrasts = max(contrasts);
    l_rfX = stro.trial(:,stro.idx.flashX) == stro.sum.exptParams.rf_x;
    l_rfY = stro.trial(:,stro.idx.flashY) == stro.sum.exptParams.rf_y;
    inRFlist = (l_rfX & l_rfY);
    
    %figure out what target the monkey picked. T1 => in RF. This is
    %necessary to compute CP (which occurs below)
    corrects = stro.trial(:, stro.idx.correct);
    choices = nan(size(stro.trial, 1), 1);
    choices(corrects & inRFlist) = 1; %in the RF and correct
    choices(~corrects & ~inRFlist) = 1; %incorrect but not in the RF
    choices(isnan(choices)) = 2; %everything else is T2 choice
    
    %****************************************************%
    %*    RASTER PLOTS AND COMPILING SPIKE COUNTS       *%
    %****************************************************%
    subplot(nFigRows, nSptFreqs, b+(nSptFreqs));
    hold on, %the raster subplot
    xlim([-0.3, 0.8]);
    set(gca, 'XTick', [-0.3, 0, .5, 0.8], 'XTickLabel', [-0.3, 0, .5, 0.8]);
    
    %first compile the noise distributions, adding lines to the raster as you go.
    counter = 0;
    noiseTlist = find(stro.trial(:, stro.idx.cntrstLev) == 1); %cntrstLev = 1 implies no stimulus was shown
    noiseChoices = choices(noiseTlist); %for CP analysis
    cellFits.nTrialsIn{a,b}(1) = length(noiseTlist);
    for k = 1:size(noiseTlist,1);
        [spikes, tFlashOn, tFlashOff, noiseStats(k)] = getTrialSpikes(stro, noiseTlist(k), cellNum, analysisMeth, anlyStart, anlyEnd);
        frameOnTime = stro.trial(noiseTlist(k), stro.idx.frameOn);
        rasterLine(tFlashOn, tFlashOff, frameOnTime, spikes, counter);
        counter = counter+1;
    end
    
    %now compile the "signal" distributions on a contrast by contrast basis. Do this seperately for trials with the stimulus inside and outside the RF. Computing the
    %ROC area and choice probability as we go.
    cellFits.crfIn{a,b}{1} = noiseStats;
    cellFits.crfOut{a,b}{1} = noiseStats;
    for c = 2:nContrasts; %2 => be careful not to include the 'zero' contrast conditions
        %add a bold line to separate the noise raster from the rest of the
        %stuff to come.
        plot([-0.3, 0.8], [counter, counter], 'k');
        counter = counter+1;
        
        %now cycle through the contrasts for trials with stimuli in the RF.
        %Plot these trials to the rasters.
        cntrstList = [stro.trial(:, stro.idx.cntrstLev) == c];
        tList = find(conditionTList & cntrstList & inRFlist); %make sure it's in the RF, and has the correct color/sf/contrast!!
        inRFchoices{c} = choices(tList); %for CP analysis
        cellFits.nTrialsIn{a,b}(c) = length(tList);
        for k = 1:size(tList, 1); %cyc through trials
            [spikes, tFlashOn, tFlashOff, inRFstats{c}(k)] = getTrialSpikes(stro, tList(k), cellNum, analysisMeth, anlyStart, anlyEnd);
            frameOnTime = stro.trial(tList(k), stro.idx.frameOn);
            rasterLine(tFlashOn, tFlashOff, frameOnTime, spikes, counter)
            counter = counter+1; 
        end %cyc through 'in RF' trials
        
        %now for trials with stimuli outside the RF.
        tList = find(conditionTList & cntrstList & ~inRFlist);
        outRFchoices{c} = choices(tList); %for CP analysis
        for k = 1:size(tList, 1); %cyc through trials
            [~, ~, ~, outRFstats{c}(k)] = getTrialSpikes(stro, tList(k), cellNum, analysisMeth, anlyStart, anlyEnd);
        end %cyc through 'out RF' trials
        
        % compile the average neural response to each stimulus. this will
        % either be the mean counts, rate, or f1 amplitude. "inRFstats" is
        % a vector of responses on each trial (either rate, counts, of f1
        % amps.
        cellFits.crfIn{a,b}{c} = inRFstats{c}; % 'c' goes from 2-nContrasts
        cellFits.crfOut{a,b}{c} = outRFstats{c};
    end %cyc through contrasts
    hold off, %the raster subplot
    ylim([0, counter])
    
    %**************************************************%
    %*       ROC ANALYSIS AND NEUROFUN PLOTING        *%
    %**************************************************%
    %calculate area under the ROC
    for c = 1:nContrasts
        if(c==1) %the 'zero contrast' condition
            cellFits.roc{a,b}(c) = roc(noiseStats, noiseStats);
        else
            cellFits.roc{a,b}(c) = roc(noiseStats, inRFstats{c});
        end
    end
    
    %start by fitting the neurometric:
    zeroList = norms{a,b} == 0;
    tmpNorms = norms{a,b}(~zeroList);
    errs = abs(0.82-cellFits.roc{a,b}(~zeroList));
    aGuess = tmpNorms(find(errs == min(errs), 1, 'last'));
    [aSSE, bSSE, ~, success(1)] = weibullFit(norms{a,b}, cellFits.roc{a,b}, 'sse', [aGuess 1]);
    correctByContrast = (cellFits.roc{a,b}.*cellFits.nTrialsIn{a,b});
    wrongByContrast = (cellFits.nTrialsIn{a,b} - correctByContrast);
    [cellFits.alpha(a,b), cellFits.beta(a,b), cellFits.gamma(a,b), success(2), modErrs] = weibullFit(norms{a,b}, [correctByContrast(:), wrongByContrast(:)], 'mle', [aSSE, bSSE]);
    neurometric = cellFits.gamma(a,b) + (0.5 - cellFits.gamma(a,b)).*exp(-(([0:0.001:max(norms{a,b})]./cellFits.alpha(a,b)).^cellFits.beta(a,b)));
    
    if ~all(success) || (sum(cellFits.roc{a,b} > 0.80) == 0);
        cellFits.alpha(a,b) = NaN;
        cellFits.beta(a,b) = NaN;
        cellFits.err.mle.alpha(a,b) = NaN;
        cellFits.err.mle.beta(a,b) = NaN;
    else
        cellFits.err.mle.alpha(a,b) = modErrs(1);
        cellFits.err.mle.beta(a,b) = modErrs(2);
    end
    
    subplot(nFigRows, nSptFreqs, b)
    hold on, %don't erase what's already there
    h1 = semilogx(norms{a,b}, cellFits.roc{a,b}, 'bo');
    h2 = semilogx([0:0.001:max(norms{a,b})], neurometric, 'b');
    hg = hggroup;  % So that legend works
    set(get(get(hg,'Annotation'),'LegendInformation'),'IconDisplayStyle','on')
    set(h1,'Parent',hg);
    set(h2,'Parent',hg);
    ylim([min(cellFits.roc{a,b})-0.02, 1.01])
    hold off,
    legend({'Behavior','Neuron'},'location','northwest');
    
    %**********************************************%
    %*        CHOICE PROBABILITY ANALYSIS         *%
    %**********************************************%
    
    %measure CP at each contrast
    zT1in = [];
    zT2in = [];
    zT1out = [];
    zT2out = [];
    for c = 1:nContrasts
        choiceRatio = nan; %initialize this.
        if(c==1)
            [cellFits.cp.in{a,b}(c), T1noise, T2noise] = choiceProb(noiseStats, noiseChoices, lowCutoff);
            cellFits.cp.out{a,b}(c) = cellFits.cp.in{a,b}(c); %just for book keeping
            choiceRatio = numel(T1noise)./numel(T2noise);
        else
            %first for trials with stim in RF
            [cellFits.cp.in{a,b}(c), T1tmp, T2tmp] = choiceProb(inRFstats{c}, inRFchoices{c}, lowCutoff);
            if all([length(T1tmp), length(T2tmp)]>=lowCutoff)
                [zT1_tmp, zT2_tmp] = normalizeDists(T1tmp, T2tmp);
                zT1in = [zT1in; zT1_tmp(:)];
                zT2in = [zT2in; zT2_tmp(:)];
                choiceRatio = numel(T1tmp) ./ numel(T2tmp);
            end
            
            
            %now for trials with stim outside RF
            [cellFits.cp.out{a,b}(c), T1tmp, T2tmp] = choiceProb(outRFstats{c}, outRFchoices{c}, lowCutoff);
            if all([length(T1tmp), length(T2tmp)]>=lowCutoff)
                [zT1_tmp, zT2_tmp] = normalizeDists(T1tmp, T2tmp);
                zT1out = [zT1out; zT1_tmp(:)];
                zT2out = [zT2out; zT2_tmp(:)];
            end
        end
        
        % trying to calculate the ratio of T1:T2 choices for each CP condition
        cellFits.cp.N_T1T2ratio{a,b}(c) = choiceRatio;
    end
    
    %accumulate the zscores across all conditions. Add the zero contrast
    %trials only once.
    grandZscores = {[grandZscores{1}; zT1in], [grandZscores{2}; zT2in]};
    if a == 1;
        if all([length(T1noise), length(T2noise)]>=lowCutoff)
            [zT1_noise, zT2_noise] = normalizeDists(T1noise, T2noise);
            grandZscores = {[grandZscores{1}; zT1_noise(:)], [grandZscores{2}; zT2_noise(:)]};
        end
    end
    
    %measure the CP pooled across contrasts (within a specific sfs/color combo)
    cellFits.cp.poolConIn.val(a,b) = choiceProb({zT1in, zT2in}, [], lowCutoff);
    if ~(isempty(zT1in) || isempty(zT2in))
        [cellFits.cp.poolConIn.p(a,b), cellFits.cp.poolConIn.h(a,b)] = ranksum(zT1in, zT2in);
    else
        cellFits.cp.poolConIn.p(a,b) = nan;
        cellFits.cp.poolConIn.h(a,b) = 0; %ranksum returns this as a logical so nan's aren't allowed
    end
    cellFits.cp.poolConOut.val(a,b) = choiceProb({zT1out, zT2out}, [], lowCutoff);
    if ~(isempty(zT1out) || isempty(zT2out))
        [cellFits.cp.poolConOut.p(a,b), cellFits.cp.poolConOut.h(a,b)] = ranksum(zT1out, zT2out);
    else
        cellFits.cp.poolConOut.p(a,b) = nan;
        cellFits.cp.poolConOut.h(a,b) = 0;
    end
    
    %measure the OUT of RF CP without normalizing. Include the zero
    %contrast trials too.
    outStats = horzcat(outRFstats{:})';
    outStats = [outStats; noiseStats(:)];
    outChoices = vertcat(outRFchoices{:});
    outChoices = [outChoices; noiseChoices];
    cellFits.cp.poolConOut_nonNorm(a,b) = choiceProb(outStats, outChoices, lowCutoff);
end %neuroFun

function [CP, T1dist, T2dist] = choiceProb(counts, choices, lowCutoff)
    %
    %   EXAMPLE : [CP, T1dist, T2dist] = choiceProb(dist, [choices], lowCutoff)
    %
    %This function takes neural data and computes CP. If dist is a vector
    %than this function expects a vector of choices. In this case, the
    %function will parse dist into a T1 and a T2 distribution. If dist is
    %an (2x1) cell array and choices is unspecified, than the funtion assumes
    %that dist is already parsed b/w T1 and T2.
    
    %have the T1 and T2 dists already been parsed?
    if ~iscell(counts)
        T1dist = counts(choices == 1);
        T2dist = counts(choices == 2);
    elseif iscell(counts)
        T1dist = counts{1};
        T2dist = counts{2};
    end
    
    %make sure that there are an adequate number of T1 and T2 decisions
    if ((length(T1dist) < lowCutoff) || (length(T2dist) < lowCutoff))
        CP = NaN;
    else
        CP = roc(T2dist, T1dist);%is this correct?
    end
end

function [zT1out, zT2out] = normalizeDists(T1dist, T2dist)
    %normalizes the counts distributions according to a pooled estimate of
    %sigma and mu.
    n1 = length(T1dist);
    n2 = length(T2dist);
    if ~n1
        sigma = std(T2dist);
        mu = mean(T2dist);
    elseif ~n2
        sigma = std(T1dist);
        mu = mean(T1dist);
    else
        %sigma = sqrt(((n1-1)*var(T1dist) + (n2-1)*var(T2dist))/(n1+n2-2)); % normal version
        sigma = sqrt( ((n1*var(T1dist) + n2*var(T2dist))./(n1+n2)) + ((n1*n2*(mean(T1dist)-mean(T2dist))^2)./(n1+n2)^2) ); % the Maunsell version
        mu = (n1.*mean(T1dist)+n2.*mean(T2dist))./ (n1+n2);
    end
    zT1out = (T1dist-mu) ./ (sigma+eps);
    zT2out = (T2dist-mu) ./ (sigma+eps);
end

function [tSpikes, tFlashOn, tFlashOff, trlStat] = getTrialSpikes(stro, trialNum, cellNum, analysisMeth, anlyStart, anlyEnd)
    %all this function does is pull out the spike times for an individual
    %trial and represent them as times from stimulus onset. Then it
    %computes the appropriate trial spiking statistic
    
    tFlashOn = stro.trial(trialNum, stro.idx.repFlashOn);
    tFlashOff = stro.trial(trialNum, stro.idx.nFrames)./stro.sum.exptParams.frame_rate+tFlashOn;
    tSpikes = stro.ras{trialNum, cellNum};
    
    %anlyStart and anlyEnd could be strings or numbers. Convert to numbers
    %as needed:
    if ischar(anlyStart)
        switch lower(anlyStart)
            case 'gabor on'
                anlyStart = tFlashOn;
            case 'frames on'
                anlyStart = stro.trial(trialNum, stro.idx.frameOn);
            otherwise
                error(sprintf('start event <%s> not yet supported'));
        end
    elseif isnumeric(anlyStart)
        anlyStart = anlyStart + tFlashOn;
    end
    if ischar(anlyEnd)
        switch lower(anlyEnd)
            case 'gabor off'
                anlyEnd = tFlashOff; %time from stim onset (simply the duration of the stimulus as calculated from fFrames)
            case 'targ on'
                anlyEnd = stro.trial(trialNum, stro.idx.targOn);
            otherwise
                error(sprintf('stop event <%s> not yet supported'));
        end
    elseif isnumeric(anlyEnd)
        anlyEnd = anlyEnd + tFlashOn;
    end
    
    %now return the relavant trial spike statistic
    switch lower(analysisMeth);
        case 'counts' %counts
            trlStat = sum((tSpikes>anlyStart) & (tSpikes<=anlyEnd));
        case 'rate' %rate in spikes/sec
            trlStat = sum((tSpikes > anlyStart) & (tSpikes<=anlyEnd)) ./ (anlyEnd-anlyStart);
        case 'f1'; %f1 response amplitude.
            if ~isempty(tSpikes)
                binwidth = 0.025;
                driftRate = stro.trial(trialNum, stro.idx.driftRate);
                edges = [anlyStart:binwidth:anlyEnd];
                psth = histc(tSpikes, edges) ./ binwidth; %psth in sp/sec
                basis1 = exp(-2*pi*sqrt(-1)*driftRate*edges);
                trlStat = 2*abs(psth(:)'*basis1(:));
            else
                trlStat = 0;
            end
    end      
end

function rasterLine(tFlashOn, tFlashOff, frameOnTime, spikes, counter)
    plot(0, counter+0.4, 'g*', 'MarkerSize', 3);
    plot(tFlashOff, counter+0.4, 'r*', 'MarkerSize', 3);
    plot((frameOnTime-tFlashOn), counter+0.4, 'c*', 'MarkerSize', 3);
    spikes = spikes-tFlashOn; %allign all the plots (and times) relative to stim on.
    spikes(spikes<-0.3) = []; %show 300ms before the stim
    spikes(spikes>0.8) = []; %only include the 134ms after the stimulus
    plot([spikes(:), spikes(:)]', [zeros(1,length(spikes))+counter; [ones(1, length(spikes)).*0.8 + counter]], 'k');
end

function datecheck(stro)
    %starting 12/18/08, I introduced zero contrast trials. This analysis
    %code is incompatible with expts prior to this update.
    monthDay = str2num(stro.sum.fileName(end-12:end-9));
    year = str2num(stro.sum.fileName(end-8:end-7));
    err = (monthDay < 1218) && (year <= 08);
    if err
        error(sprintf('\n DTmocsUnpack isn''t backwards compatabile with experiments prior to 12/18/08'));
    end
end

function [cellNum, anlyStart, anlyEnd, lowCutoff, analysisMeth] = getAnalysisParams(nCellChannels, defaultParams)
    
    %DTmocsUnpack can be called from the command window, some other
    %function, or by DTgui. Depending on which calling method is used,
    %default params may be a scalar(outside function) or a structure(gui
    %update window). If default params is unspecified, than popup the
    %following userinput menu.
    if isnumeric(defaultParams) && defaultParams > 0;
        cellNum = 1;
        anlyStart = 'gabor on';
        anlyEnd = 'gabor off';
        lowCutoff = 5;
        analysisMeth = 'rate';
        return
    elseif isstruct(defaultParams)
        cellNum = defaultParams.cellNum;
        anlyStart = defaultParams.start;
        anlyEnd = defaultParams.end;
        lowCutoff = defaultParams.lowCutoff;
        analysisMeth = defaultParams.meth;
        return
    end
    
    name = 'DT analysis params';
    prompt = {'Unit Number', 'Start Time From Gabor Onset(sec)', 'End Time From Gabor Onset (sec)', 'min trials for CP', 'analysis method (''counts''; ''rate''; ''F1'')'};
    numlines = 1;
    if nCellChannels > 1
        defaultAnswers = {sprintf('%d units available', nCellChannels), 'gabor on', 'gabor off', '5', 'rate'};
    else
        defaultAnswers = {'1', 'gabor on', 'gabor off', '5', 'rate'};
    end
    options.Resize='on';
    options.WindowStyle='normal';
    options.Interpreter='tex';
    params = inputdlg(prompt, name, numlines, defaultAnswers, options);
    
    %determine which neuron to analyze (if more than one is present)
    cellNum = str2double(params{1});
    if isempty(cellNum) || (cellNum>nCellChannels);
        error('Unit number is unspecified or unavailable');
    end
    
    
    %deal with the start time. Set to one of the trial events for a precise
    %(trial by trial) time. Otherwise, specify a time from onset.
    acceptedOnsetEvents = {'gabor on', 'frames on'};
    if any(strcmp(params{2}, acceptedOnsetEvents))
        anlyStart = params{2};
    elseif ~isempty(str2double(params{2}))
        anlyStart = str2double(params{2});
    else
        error('unknown start time designation')
    end
    
    
    %deal with the analysis offset time
    acceptedOffsetEvents = {'gabor off', 'targ on', 'sac on'};
    if any(strcmp(params{3}, acceptedOffsetEvents))
        anlyEnd = params{3};
    elseif ~isnan(str2double(params{3}))
        anlyEnd = str2double(params{3});
    else
        error('unknown end time designation')
    end
       
    lowCutoff = str2double(params{4});
    analysisMeth = params{5};
    disp('done with defaults')
end



%****************************************************************%
%
%                       TESTING CODE
%
%****************************************************************%


function testROC()
% CODE FOR TESTING DTMOCSUNPACK
% MAKES A TEST DATA FILE WITH A KNOWN ROC AREA

%START BY OPENING A DATA FILE AND GUTTING ALL THE REAL DATA.
    stro = dtobj('S061609005');
    [stro.ras{:,1}] = deal([]);
    l_color2 = stro.trial(:, stro.idx.colorDir) == 2;
    stro.trial(l_color2, :) = [];
    stro.ras(l_color2, :) = [];
    stro.trial = [stro.trial; stro.trial; stro.trial; stro.trial; stro.trial; stro.trial; stro.trial; stro.trial];
    stro.ras = [stro.ras; stro.ras; stro.ras; stro.ras; stro.ras; stro.ras; stro.ras; stro.ras];



    %setting up the psychometric results
    color1 = abs(stro.sum.exptParams.RF_colors(1:3))+eps;
    triplets(1,2:8) = logspace(log10(0.25.*color1(1)), log10(color1(1)), 7);
    triplets(2,2:8) = logspace(log10(0.25.*color1(2)), log10(color1(2)), 7);
    triplets(3,2:8) = logspace(log10(0.25.*color1(3)), log10(color1(3)), 7);
    contrasts = sqrt(sum(triplets.^2, 1));
    alpha = contrasts(5);
    beta = 7;
    binoProbs = 1 - 0.5.*exp(-(contrasts./alpha).^beta);
    semilogx(contrasts, binoProbs, 'b');
    stro.trial(:, stro.idx.correct) = nan;
    for a = 1:length(unique(stro.trial(:, stro.idx.cntrstLev)));
        l_trial = stro.trial(:, stro.idx.cntrstLev) == a;
        nCorrect = binornd(sum(l_trial), binoProbs(a));
        tNums = find(l_trial);
        correctTrials = tNums(1:nCorrect);
        incorrectTrials = tNums(nCorrect+1:end);
        stro.trial(correctTrials, stro.idx.correct) = 1;
        stro.trial(incorrectTrials, stro.idx.correct) = 0;
    end


    %setting up the neurometric results. Modeling spike counts on each trial as
    %a draw from a poisson distribution. Determining the outcome by creating a
    %lookup table.
    lambdas = logspace(log10(1), log10(30), 30);
    lambdas(end+1) = 50;
    distributions = poissrnd(repmat(lambdas, 10000, 1));
    for a = 1:size(distributions,2);
        rocVal(a) = roc(distributions(:,1), distributions(:,a));
    end
    l_redundant = find(rocVal == 1);
    lambdas(l_redundant(2:end)) = [];
    rocVal(l_redundant(2:end)) = [];
    meanCounts = interp1(rocVal, lambdas, binoProbs);

    %populate the trials with spike counts.
    for a = 1:size(stro.trial,1)
        startTime = stro.trial(a, stro.idx.flashOn);
        stopTime = stro.trial(a, stro.idx.flashOff);
        cntrst = stro.trial(a, stro.idx.cntrstLev);
        nSpikes = poissrnd(meanCounts(cntrst));
        if stro.trial(a,stro.idx.flashX) == stro.sum.exptParams.rf_x;
            stro.ras{a,1} = linspace(startTime+0.010, stopTime-0.010, nSpikes);
        else
            stro.ras{a,1} = linspace(startTime+0.010, stopTime-0.010, poissrnd(meanCounts(1)));
        end
    end


end

function testCP()
% this script will create a data file with no dependence between firing
% rates and choices (i.e., cp should = 0.5). Also, firing rates are drawn
% from two poissions. One for the bottom four contrasts and another for the
% highest four contrasts. The consequence should be a step function in the
% neurometric function. Lastly, corrects (in the stro.trial field) are
% assigned randomly (and hence T1 choices). Thus the psychometric function
% should be flat and = 0.5.

    % How many times should I iterate over this simulation (to get a feel
    % for the mean CP as a result of the sim)
    nIters = 200;

    %open a data file and gut the data from one of the color directions.
    stro = dtobj('S061609005');
    [stro.ras{:,1}] = deal([]);
    l_color2 = stro.trial(:, stro.idx.colorDir) == 2;
    stro.trial(l_color2, :) = [];
    stro.ras(l_color2, :) = [];
    
    %make the data file four times as long.
    stro.trial = [stro.trial; stro.trial; stro.trial; stro.trial; stro.trial; stro.trial; stro.trial; stro.trial];
    stro.ras = [stro.ras; stro.ras; stro.ras; stro.ras; stro.ras; stro.ras; stro.ras; stro.ras];
    
    for i = 1:nIters;
        disp(i)
        
        %allocate choices randomly.
        stro.trial(:, stro.idx.correct) = unidrnd(2, size(stro.trial,1), 1) -1;
        
        %loop through trials, assign spike times. I'll draw spike counts from
        %two different distributions: one for the lowest four contrasts and
        %another for the highest four contrasts.
        meanCounts = [3, 30];
        for a = 1:size(stro.trial,1);
            startTime = stro.trial(a, stro.idx.flashOn);
            stopTime = stro.trial(a, stro.idx.flashOff);
            normCntrst = (stro.trial(a, stro.idx.cntrstLev) > 4) + 1; %cntrsts are spilt between the 4 highest and lowest
            nSpikes = poissrnd(meanCounts(normCntrst));
            if stro.trial(a,stro.idx.flashX) == stro.sum.exptParams.rf_x;
                stro.ras{a,1} = sort(unifrnd(startTime+0.010, stopTime-0.020, nSpikes, 1));
            else
                stro.ras{a,1} = sort(unifrnd(startTime+0.010, stopTime-0.020, poissrnd(meanCounts(1)), 1));
            end
        end
        
% %         %     plot histograms of counts for T1 and T2 choices by contrast
% %         l_corrects = stro.trial(:, stro.idx.correct);
% %         l_inRF = stro.trial(:,flashX) == stro.sum.exptParams.rf_x;
% %         l_t1 = (l_corrects & l_inRF) | (~l_corrects & ~l_inRF);
% %         figure
% %         for a = 1:length(unique(stro.trial(:, stro.idx.cntrstLev)));
% %             l_tlist = (stro.trial(:,stro.idx.cntrstLev) == a) & l_inRF;
% %             
% %             subplot(length(unique(stro.trial(:, stro.idx.cntrstLev))), 1, a)
% %             hold on,
% %             tmp1 = cellfun(@length, stro.ras(l_t1 & l_tlist, 1));
% %             [cnts1, pos1] = hist(tmp1);
% %             bar(pos1, cnts1, 'b')
% %             tmp2 = cellfun(@length, stro.ras(~l_t1 & l_tlist, 1));
% %             [cnts2, pos2] = hist(tmp2);
% %             bar(pos2, cnts2, 'r')
% %             plot(mean(tmp1), max([cnts1(:);cnts2(:)])+3, 'bv')
% %             plot(mean(tmp2), max([cnts1(:);cnts2(:)])+3, 'rv')
% %             hold off
% %         end
        
        %unpack the data and hold onto the cp.
        [m, c, expt] = DTunpack(stro,1);
        close all
        cp(i,:) = c.cp.in{1};
        poolIn_p(i) = c.cp.poolConIn.p;
        poolIn_val(i) = c.cp.poolConIn.val;
    end
    
    %CONCLUSIONS: If I iterate 100 times, the resulting mean CP is between
    %~0.5 with and SD = 0.05. Distributions of CPs appear symetric about
    %the mean. For pooled CP, the mean val was ~0.5 and SD = 0.02.
    %Importantly, about 5% of the simulated T1 dists were significantly
    %differnt than the T2 dist (i.e. poolConIn.p < 0.05). I think that this
    %is what we expect for a random process. The distribution of p values
    %for pooled CP looked like it could/would be uniform with more
    %iterations.

end
