function out = dtcp(DT, permuteChoices, alternateStat)
%takes a DT object and computes the CP for various conditions in the
%experiment.

if ~exist('permuteChoices', 'var')
    permuteChoices = 0;
end

if ~exist('alternateStat', 'var');
    alternateStat = [];
end

%first, strip out the grating repeat trials
DT = stripOutGratingTrials(DT);

%now determine the number of spikes elicited during the stimulus period for
%each trial
if isempty(alternateStat)
    tStart = DT.trial(:, DT.idx.repFlashOn);
    nFrames = DT.trial(:, DT.idx.nFrames);
    tEnd = (nFrames ./ DT.sum.exptParams.frame_rate) + tStart;
    tEnd_cell = mat2cell(tEnd, ones(length(tEnd),1));
    tStart_cell = mat2cell(tStart, ones(length(tStart),1));
    counts = cellfun(@(r,s,e) sum((r>s)&(r<=e)), DT.ras(:, DT.idx.spikes), tStart_cell, tEnd_cell);
    rates = counts ./ (tEnd-tStart);
else
    rates = alternateStat;
end

% unpack the colors used during the experiment, and the LMS values for each
% trial type
exptColors = reshape(DT.sum.exptParams.RF_colors, 3, 3)';
exptColors(sum(abs(exptColors),2)==0,:) = [];
exptColors = exptColors ./ repmat(sqrt(sum(exptColors.^2,2)), 1, 3);
trialColors = unique(DT.LMS ,'rows');



% determine which target the monkey chose
l_inRF = (DT.trial(:, DT.idx.flashX) == DT.sum.exptParams.rf_x) & (DT.trial(:, DT.idx.flashY) == DT.sum.exptParams.rf_y);
l_correct = DT.trial(:, DT.idx.correct);
l_T1 = (l_inRF & l_correct) | (~l_inRF & ~l_correct);
l_T2 = (l_inRF & ~l_correct) | (~l_inRF & l_correct);

%cycle through the trial types and determine the CP. Also develop a pooled
%estimate of CP.
val{1,2} = [];
norm{1,2} = [];
nTrials{1,2} = [];
rocval{1,2} = [];
zT1 = [];
zT2 = [];
for clr = 1:size(trialColors,1)
    l_trialLMS = ismember(DT.LMS, trialColors(clr,:), 'rows');
    if trialColors(clr,:) == [0 0 0];
        tList = l_trialLMS & ~isnan(rates); %the alternate stat can be a nan. remove those here
    else
        tList = l_inRF & l_trialLMS & ~isnan(rates);
    end
    
    %compute the ROC val
    noiseTrials = ismember(DT.LMS, [0 0 0], 'rows');
    noiseRates = rates(noiseTrials & ~isnan(rates));
    trialRates = rates(tList);
    if all([numel(noiseRates), numel(trialRates)]>=8)
        rocest = roc(noiseRates, trialRates);
    else
        rocest = nan;
    end
    
    %determine the T1 and T2 distributions of rates. Permute if need be.
    T1dist = rates(tList&l_T1);
    T2dist = rates(tList&l_T2);
    if permuteChoices
        nT1 = numel(T1dist);
        nT2 = numel(T2dist);
        tmp = [T1dist(:); T2dist(:)];
        ind = randperm(numel(tmp));
        shuff = tmp(ind);
        T1dist = shuff(1:nT1);
        T2dist = shuff(nT1+1:end);
        
        %error checking...
        if (numel(T1dist) ~= nT1) || (numel(T2dist) ~= nT2)
            keyboard
        end
        if all(tmp == [T1dist;T2dist]) && all([nT1, nT2]>0) && std(tmp)
            %keyboard
        end
    end
    
    
    %compute the CP estimate for a given color/contrast condition
    if any([numel(T1dist), numel(T2dist)]<5)
        cpest = NaN;
    else
        cpest = roc(T2dist, T1dist);%is this correct?
    end
    
    %compute the zscores for this color/cont condition
    if all([numel(T1dist), numel(T2dist)]>=5)
        n1 = length(T1dist);
        n2 = length(T2dist);
        if ~n1
            sigma = std(T2dist);
            mu = mean(T2dist);
        elseif ~n2
            sigma = std(T1dist);
            mu = mean(T1dist);
        else
            sigma = sqrt(((n1-1)*var(T1dist) + (n2-1)*var(T2dist))/(n1+n2-2));
            mu = (n1.*mean(T1dist)+n2.*mean(T2dist))./ (n1+n2);
        end
        zT1tmp = (T1dist-mu) ./ (sigma+eps);
        zT2tmp = (T2dist-mu) ./ (sigma+eps);
        
        zT1 = [zT1; zT1tmp];
        zT2 = [zT2; zT2tmp];
    end
    
    %store the data temporarly    
    [~, ind] = max(exptColors * trialColors(clr,:)');
    if sum(abs(trialColors(clr,:)))>0
        val{ind}(numel(val{ind})+1) = cpest;
        norm{ind}(numel(norm{ind})+1) = exptColors(ind,:) * trialColors(clr,:)';
        nTrials{ind}(numel(nTrials{ind})+1) = sum(tList);
        rocval{ind}(numel(rocval{ind})+1) = rocest;
    else
        zeroContCP = cpest;
        zeroContNtrials = sum(tList);
        zeroContROC = rocest;
    end
end %cycling through trial color types

%rearrange the cp, roc, and ntrials in ascending order of contrast
for a = 1:2
    [tmpNorms, inds] = sort(norm{a});
    out.norms(a,:) = [0,tmpNorms];
    out.cp(a,:) = [zeroContCP,val{a}(inds)];
    out.nTrials(a,:) = [zeroContNtrials, nTrials{a}(inds)];
    out.roc(a,:) = [zeroContROC, rocval{a}(inds)];
end

%add some other miscelaneous CP values:
out.cpout = computeCP((~l_inRF&l_T1&~isnan(rates)), (~l_inRF&l_T2&~isnan(rates)), rates);
if all([numel(zT1), numel(zT2)]>0)
    out.zScoreIn = roc(zT2, zT1);
    [out.zScoreIn_P,~] = ranksum(zT2, zT1);
else
    out.zScoreIn = nan;
    out.zScoreIn_P = nan;
end




% NESTED SUBFUNCTION:
    function out = computeCP(T1, T2, FR)
        T1rates = FR(T1);
        T2rates = FR(T2);
        if any([length(T1rates), length(T2rates)]<5)
            out = NaN;
        else
            out = roc(T2rates, T1rates);%is this correct?
        end
        
    end
    


end
