%% (24) PREDICTING TRs TUNING ON BASIS OF GT TUNING
%
% Here I'll test whether neurons tuning during GT agrees with the TR's
% estimated during DT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, clc
global cardVsIntBatchPath blpBatchPath
preprocessDTbatchData


l_validConds = ~(commonExclusions | out.errors(:, prefTieInd));
expList = find(l_validConds);
[DT_cardNT, DT_cardTR, GT_cardResp, DT_intNT, DT_intTR, GT_intResp] = deal(nan(length(expList),1));
for a = 1:length(expList)
    ex = expList(a);
    
    gt = gtobj(out.fnames{a}{2});
    stimon = mat2cell(gt.trial(:, gt.idx.stimon), ones(size(gt.trial,1),1));
    stimoff = mat2cell(gt.trial(:, gt.idx.stimoff), ones(size(gt.trial,1),1));
    nSpikes = cellfun(@(x,y,z) sum((x>y)&(x<z)), gt.ras(:,gt.idx.spikes), stimon, stimoff);
    l_Card = gt.getTrialList(4, sign(out.dat(ex).prefCard));
    l_Int = gt.getTrialList(4, sign(out.dat(ex).prefInt));
    [p, h] = ranksum(nSpikes(l_Card), nSpikes(l_Int));
    if (p<1)
        for clr = 1:2
            DTclr = sign(out.dat(ex).expt.standColors(clr,:));
            tmpNT = out.dat(ex).c.alpha(clr,1);
            tmpTR = out.dat(ex).c.alpha(clr,1) ./ out.dat(ex).m.alpha(clr,1);
            tmpTR = min(tmpTR, inf); %trying to deal with the nan TRs. I think this will make the lt/gt statements below work
            indForGT = ismember(sign(out.dat(ex).grating.color.colors), DTclr, 'rows');
            tmpGTresp = out.dat(ex).grating.color.colresp(indForGT,1);
            iscard = all(DTclr == [1 -1 0]) || all(DTclr == [0 0 1]);
            isint = all(DTclr == [1 -1 1]) || all(DTclr == [1 -1 -1]);
            if iscard
                DT_cardTR(a) = tmpTR;
                DT_cardNT(a) = tmpNT;
                GT_cardResp(a) = tmpGTresp;
            elseif isint
                DT_intTR(a) = tmpTR;
                DT_intNT(a) = tmpNT;
                GT_intResp(a) = tmpGTresp;
            else
                error('color not recognized')
            end
        end
    end
end

%setting up a contingency table to compute the Fisher Exact Test, whose
%null hypothesis is that the two outcomes are independent (i.e., the
%proportion of DTcards is independent of a neurons' behavior in GT).
%Rejecting the null hypothesis (i.e., that the two outcomes are dependent)
%would suggest that GT data predicts DT data.
%   H0 => the two marginal distributions are independent, and a neuron's
% DTcardinal status is independent of its GT cardinal status.
%   Halt => a neuron's DT status is contingent on its GT status.
f11 = sum((DT_cardTR < DT_intTR) & (GT_cardResp > GT_intResp));
f21 = sum((DT_cardTR > DT_intTR) & (GT_cardResp > GT_intResp));
f12 = sum((DT_cardTR < DT_intTR) & (GT_cardResp < GT_intResp));
f22 = sum((DT_cardTR > DT_intTR) & (GT_cardResp < GT_intResp));

alpha = 0.05;
[h, p] = fisherExact(f11,f21,f12,f22, alpha)


figure %the TRs vs GT
tmpDT_cardTR = DT_cardTR;
tmpDT_intTR = DT_intTR;
tmpDT_cardTR(isinf(tmpDT_cardTR)) = nan;
tmpDT_intTR(isinf(tmpDT_intTR)) = nan;
l_bothNan = isnan(min([tmpDT_cardTR, tmpDT_intTR], [],2));
valfornan = max([tmpDT_cardTR(:);tmpDT_intTR(:)])*1.5;
tmpDT_cardTR(isnan(tmpDT_cardTR)) = valfornan;
tmpDT_intTR(isnan(tmpDT_intTR)) = valfornan;
DTratio = tmpDT_cardTR./tmpDT_intTR;
DTratio(l_bothNan) = nan; %convert the ones to nans.
GTratio = GT_cardResp./GT_intResp;
loglog(DTratio, GTratio, 'o')
[rho, p] = corr(DTratio(~l_bothNan), GTratio(~l_bothNan), 'type', 'spearman');
xlabel('DT ratio of TRs');
ylabel('GT ratio of resps');
title(sprintf('rho=%.3f, p=%.3f', rho, p))

figure %the NTs vs GT
tmpDT_cardNT = DT_cardNT;
tmpDT_intNT = DT_intNT;
l_bothNan = isnan(min([tmpDT_cardNT, tmpDT_intNT], [],2));
valfornan = max([tmpDT_cardNT(:);tmpDT_intNT(:)])*1.5;
tmpDT_cardNT(isnan(tmpDT_cardNT)) = valfornan;
tmpDT_intNT(isnan(tmpDT_intNT)) = valfornan;
DTratio = tmpDT_cardNT./tmpDT_intNT;
DTratio(l_bothNan) = nan; %convert the ones to nans.
GTratio = GT_cardResp./GT_intResp;
loglog(DTratio, GTratio, 'o')
[rho, p] = corr(DTratio(~l_bothNan), GTratio(~l_bothNan), 'type', 'spearman');
xlabel('DT ratio of NTs');
ylabel('GT ratio of resps');
title(sprintf('rho=%.3f, p=%.3f', rho, p))

%% (35) LOOKING FOR THE MECHANISM DIRECTION THAT BEST PREDICTS THE DT CRFs
% The results of this analysis are odd. For S-iso, the model for the joint
% CRF produces an effect where the card and int CRF are identical. for L-M,
% the joint CRF model produces an effect where the card and int CRFs are
% different, and the ordinal relationship changes (int on top vs. card on
% top)...


clear, clc
global cardVsIntBatchPath
preprocessDTbatchData

MAKEPLOT = 1;
SCALE = 0; %a binary flag to scale the S axis?
TRIALLENGTH = 0.666; %in sec
CARD = 0;
INT = 1;
l_valid = ~(commonExclusions | (out.errors(:, neuroThresh1Ind) & out.errors(:,neuroThresh2Ind)));
exptList = find(l_valid);
[mechVect, betaMech1, betaMech2, pMech1, pMech2] = deal(nan(length(exptList), 3));
exitFlag = false(numel(exptList),1);
plotList = {'S030411005','S061209002', 'K120310002'};
for a = 1:length(exptList);
    ex = exptList(a);
    fprintf('Expermient <%s>, ex:%d \n', out.fnames{ex}{1}, ex);
    
    %find the rates in response to Card and Int
    sIdx = ismember(sign(out.dat(ex).expt.standColors), [0 0 1], 'rows');
    lvmIdx = ismember(sign(out.dat(ex).expt.standColors), [1 -1 0], 'rows');
    cardIdx = sIdx|lvmIdx;
    [cardRates, intRates, cardNorms, intNorms] = deal([]); %initialize these variables
    for cntrst = 1:length(out.dat(ex).expt.norms{1})
        %card
        cardRates = [cardRates , out.dat(ex).c.crfIn{cardIdx}{cntrst}];
        nCard = length(out.dat(ex).c.crfIn{cardIdx}{cntrst});
        cardNorms = [cardNorms, repmat(out.dat(ex).expt.norms{cardIdx}(cntrst), 1, nCard)];
        %int
        intRates = [intRates , out.dat(ex).c.crfIn{~cardIdx}{cntrst}];
        nInt = length(out.dat(ex).c.crfIn{~cardIdx}{cntrst});
        intNorms = [intNorms, repmat(out.dat(ex).expt.norms{~cardIdx}(cntrst), 1, nInt)];
    end
    
    %prepare for the gradient descent.
    counts = ([cardRates(:); intRates(:)].*TRIALLENGTH);
    intUnit = out.dat(ex).expt.standColors(~cardIdx,:);
    intUnit = intUnit./norm(intUnit);
    cardUnit = out.dat(ex).expt.standColors(cardIdx,:);
    cardUnit = cardUnit./norm(cardUnit);
    intLMS = (intUnit(:) * intNorms(:)')'; %converting from norms to LMS triplets
    cardLMS = (cardUnit(:) * cardNorms(:)')';
    LMS = [cardLMS ; intLMS];
    if SCALE
        LMS = bsxfun(@times, LMS, [1,1,.2]);
    end
    [mechGuesses, exitVals, flags] = DTfindMechVect(LMS, counts, 'deviance');
    
    %DTfindMechVect yields 2 guesses. Pick the one with the lowest
    %deviance.
    idx = false(size(mechGuesses,1),1);
    [~, tmp] = min(exitVals);
    idx(tmp) = true;
    mechVect(a,:) = mechGuesses(idx,:);
    exitFlag(a) = flags(idx);
    
    % now do the poisson regression in two ways, as a projection onto the
    % cardinal axis, and as a projection onto the mechanism vector....
    % Start with the cardinal axis.
    if exitFlag(a)
        mech1Proj = abs(LMS*mechGuesses(idx,:)');
        mech2Proj = abs(LMS*mechGuesses(~idx,:)');
        color = [repmat(CARD, length(cardRates), 1); repmat(INT, length(intRates), 1)];
        
        %mechinism direction case 1
        mech1Interact = mech1Proj .* color;
        mech1Predictors = [mech1Proj, mech1Interact];
        [betaMech1(a,:), ~, statsMech1] = glmfit(mech1Predictors, counts, 'poisson', 'link', 'log');
        pMech1(a,:) = statsMech1.p';
        
        %mechinism direction case 2
        mech2Interact = mech2Proj .* color;
        mech2Predictors = [mech2Proj, mech2Interact];
        [betaMech2(a,:), ~, statsMech2] = glmfit(mech2Predictors, counts, 'poisson', 'link', 'log');
        pMech2(a,:) = statsMech2.p';
    end
    
    
    %do some ploting if desired
    if MAKEPLOT;% && (rem(a,3) == 0); %&& any(strcmp(out.fnames{ex}{1}, plotList))
        l_card = color == CARD;
        l_int = color == INT;
        cardAlphaLMS = out.dat(ex).m.alpha(cardIdx) * cardUnit;
        intAlphaLMS = out.dat(ex).m.alpha(~cardIdx) * intUnit;
        
        
        mech1Interp = logspace(log10(min(mech1Proj(mech1Proj~=0))*.95), log10(max(mech1Proj).*1.07), 150);
        mech1InterpCard = glmval(betaMech1(a,:)', [mech1Interp(:), mech1Interp(:).*CARD], 'log');
        mech1CardAlphaNorm = abs(cardAlphaLMS * mechGuesses(idx,:)');
        mech1CardAlphaCount = glmval(betaMech1(a,:)', [mech1CardAlphaNorm, mech1CardAlphaNorm*CARD], 'log');
        mech1InterpInt = glmval(betaMech1(a,:)', [mech1Interp(:), mech1Interp(:).*INT], 'log');
        mech1IntAlphaNorm = abs(intAlphaLMS * mechGuesses(idx,:)');
        mech1IntAlphaCount = glmval(betaMech1(a,:)', [mech1IntAlphaNorm, mech1IntAlphaNorm*INT], 'log');
        
        mech2Interp = logspace(log10(min(mech2Proj(mech2Proj~=0))*.95), log10(max(mech2Proj).*1.07), 150);
        mech2InterpCard = glmval(betaMech2(a,:)', [mech2Interp(:), mech2Interp(:).*CARD], 'log');
        mech2CardAlphaNorm = abs(cardAlphaLMS * mechGuesses(~idx,:)');
        mech2CardAlphaCount = glmval(betaMech2(a,:)', [mech2CardAlphaNorm, mech2CardAlphaNorm*CARD], 'log');
        mech2InterpInt = glmval(betaMech2(a,:)', [mech2Interp(:), mech2Interp(:).*INT], 'log');
        mech2IntAlphaNorm = abs(intAlphaLMS * mechGuesses(~idx,:)');
        mech2IntAlphaCount = glmval(betaMech2(a,:)', [mech2IntAlphaNorm, mech2IntAlphaNorm*INT], 'log');
        
        
        figure
        set(gcf, 'position', [83 357 1031 427]);
        subplot(1,3,1), hold on; %proj onto mechVect1
        plot(mech1Proj(l_card), counts(l_card), 'k.')
        plot(mech1Proj(l_int), counts(l_int), 'b.')
        plot(mech1Interp, mech1InterpCard, 'k')
        plot(mech1Interp, mech1InterpInt, 'b')
        plot(mech1CardAlphaNorm, mech1CardAlphaCount, 'kv', 'markersize', 8, 'markerfacecolor', 'k');
        plot(mech1IntAlphaNorm, mech1IntAlphaCount, 'b^', 'markersize', 8, 'markerfacecolor', 'b');
        ylim([0, max(counts).*1.1])
        xlim([min(mech1Proj(mech1Proj~=0))*.95, max(mech1Proj).*1.07])
        set(gca, 'xscale', 'log')
        title(sprintf('Mech1: [%s]', num2str(mechGuesses(idx,:))))
        xlabel('Proj onto mech vect1')
        ylabel('Counts')
        legend('Card', 'Int', 'location', 'northwest');
        legend boxoff
        hold off,
        
        subplot(1,3,2), hold on; %proj onto mechVect2
        plot(mech2Proj(l_card), counts(l_card), 'k.')
        plot(mech2Proj(l_int), counts(l_int), 'b.')
        plot(mech2Interp, mech2InterpCard, 'k')
        plot(mech2Interp, mech2InterpInt, 'b')
        plot(mech2CardAlphaNorm, mech2CardAlphaCount, 'kv', 'markersize', 8, 'markerfacecolor', 'k');
        plot(mech2IntAlphaNorm, mech2IntAlphaCount, 'b^', 'markersize', 8, 'markerfacecolor', 'b');
        ylim([0, max(counts).*1.1])
        xlim([min(mech2Proj(mech2Proj~=0))*.95, max(mech2Proj).*1.07])
        set(gca, 'xscale', 'log')
        title(sprintf('Mech2: [%s]', num2str(mechGuesses(~idx,:))))
        xlabel('Proj onto mech vect2')
        ylabel('Counts')
        legend('Card', 'Int', 'location', 'northwest');
        legend boxoff
        hold off,
        
        %plot the deviance function and the color directions tested, along
        %with the supposed mechanism direction
        lvmUnit = [1 -1 0]./norm([1 -1 0]);
        sUnit = [0 0 1];
        theta = linspace(0, pi, 300);
        dev = nan(size(theta));
        for t = 1:length(theta)
            clrVec = (cos(theta(t))*lvmUnit(:)) + (sin(theta(t))*sUnit(:));
            clrVec = clrVec./norm(clrVec); %make it a unit vector
            proj = abs(LMS*clrVec(:));
            
            %mechinism direction case
            [~, dev(t)] = glmfit(proj, counts, 'poisson', 'link', 'log');
        end
        cardTheta = mod(atan(cardUnit(3)./(cardUnit*lvmUnit(:))), pi);
        intTheta = mod(atan(intUnit(3)./(intUnit*lvmUnit(:))), pi);
        mechTheta = mod(atan(mechGuesses(:,3)./(mechGuesses*lvmUnit(:))), pi);
        subplot(1,3,3), hold on,
        plot(theta*(180/pi), dev);
        plot(cardTheta*(180/pi), 0, 'k*')
        plot(intTheta*(180/pi), 0, 'b*')
        plot(mechTheta(idx)*(180/pi), exitVals(idx), 'g^')
        plot(mechTheta(~idx)*(180/pi), exitVals(~idx), 'mv')
        xlim([0,180])
    end
    
end


%summary of beta values
figure, hold on,
scatter(betaMech1(:,2), betaMech1(:,2)+betaMech1(:,3), 'bo');
x = get(gca, 'xlim');
y = get(gca, 'ylim');
maxxy =max([x,y]);
plot([0,maxxy], [0, maxxy], 'k')
axis tight
hold off

%summary of mechanism directions
figure
x = mechVect * ([1/sqrt(2), -1/sqrt(2),0])';
y = mechVect * [0,0,1]';
theta =atan(y./x);
r1 = rose(mod(theta, pi));
% hold on
% xx = prefIsolum(l_valid,:) * ([1/sqrt(2), -1/sqrt(2),0])';
% yy = prefIsolum(l_valid,:) * [0,0,1]';
% thetaPrefIsolum = atan(yy./xx);
% r2 = rose(mod(thetaPrefIsolum, pi));
% set(r2, 'color', [1 0 0])

%% (36) COLOR TUNING ON THE BASIS OF GT
% taking the responses to the nine stimuli used during GT protocol 4 and
% fitting a linear model. Using this fit as an estimate of the neurons
% preferred color direction. Various scalings of the S-axis are possible.



clear, clc
global cardVsIntBatchPath
preprocessDTbatchData

SCALEFACTOR = 0.2;
l_valid = ~(commonExclusions | (out.errors(:, neuroThresh1Ind) & out.errors(:,neuroThresh2Ind)));
exptList = find(l_valid);
prefColor = nan(length(exptList), 3);
for a = 1:length(exptList)
    ex = exptList(a);
    fprintf('Analyzing <%s>, file %d out of %d\n', out.fnames{ex}{2}, a, length(exptList));
    
    %unpack the GT file
    GT = gtobj(out.fnames{ex}{2});
    L4 = GT.trial(:, GT.idx.protocol) == 4;
    Lcc = GT.trial(:, GT.idx.lcc);
    Mcc = GT.trial(:, GT.idx.mcc);
    Scc = GT.trial(:, GT.idx.scc) .* SCALEFACTOR;
    stimon_t = GT.trial(:, GT.idx.stimon);
    stimoff_t = stimon_t + (GT.trial(:, GT.idx.nFrames) ./ GT.sum.exptParams.framerate);
    stimon_t = mat2cell(stimon_t, ones(size(stimon_t)));
    stimoff_t = mat2cell(stimoff_t, ones(size(stimoff_t)));
    spikerates = cellfun(@(b,e,r) sum((r>b)&(r<=e))./(e-b), stimon_t, stimoff_t, GT.ras(:, GT.idx.spikes));
    
    %run the fitting routine
    errorflag = 0;
    if sum(L4) > 8
        colordirections = [Lcc(L4) Mcc(L4) Scc(L4)];
        uniquecolordirections = unique(colordirections,'rows');
        [normresps,rawresps,sds,ns] = deal([]); %initialize the vectors
        for i = 1:size(uniquecolordirections,1)
            cdir = uniquecolordirections(i,:);
            L = Lcc == cdir(1) &  Mcc == cdir(2) & Scc == cdir(3);
            normresps = [normresps; mean(spikerates(L&L4))];   % NOT normalizing by cone contrast
            rawresps = [rawresps;mean(spikerates(L&L4))];
            sds = [sds; std(spikerates(L&L4))];
            ns = [ns; sum(L&L4)];
        end
        responses = spikerates(L4);
        maxrespidx = find(responses == max(responses),1);
        initguess = colordirections(maxrespidx,:);
        initguess = initguess.*(responses(maxrespidx)./norm(initguess));
        options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-6,'TolX',10^-6);
        [mechVect, fval, exitflag] = fminsearch(@(x) linmodfiterr(colordirections, responses, x), initguess, options);
        
        if exitflag
            prefColor(a,:) = mechVect./norm(mechVect);
        end
    end
end




%plot a polar histogram of preferred color directions
lvmUnit = [1./sqrt(2), -1./sqrt(2), 0];
x = prefColor * lvmUnit(:);
y = prefColor(:,3);
theta = atan(y./x);
rose(mod(theta, pi))

%% COLOR SENSITIVITY VS. DIRECTION SELECTIVITY


clear, clc, close all
global cardVsIntBatchPath
preprocessDTbatchData



% loop through the GT files and compute the direction selectivity index.
% Try using the raw data, and some type of parametric fit to the data.

figure
DS_raw = nan(numel(out.dat),1);
for a = 1:numel(out.dat)
    if any(out.dat(a).grating.orient.resp(:,3) < 3); disp('got one'); continue; end
    
    orients = out.dat(a).grating.orient.stim;
    resp = out.dat(a).grating.orient.resp(:,1);
    
    % DS estimated from the raw data
    [pref_resp, idx] = max(resp);
    pref_orient = orients(idx);
    null_orient = mod(pref_orient + pi, 2*pi);
    null_resp = resp(softEq(orients,null_orient));
    DS_raw(a) = 1-(null_resp ./ pref_resp);
    
    % DS estimated by spline fits
    Ltmp = orients == min(orients);
    y = [resp; resp(Ltmp)];
    x = [orients; orients(Ltmp)+2*pi];
    pp = csape(x,y,'periodic');
    xx = linspace(0,2*pi,100);
    fit = ppval(pp,xx);
    
    
    cla,
    hold on,
    set(gca, 'fontsize', 14, 'box', 'off')
    plot(orients, resp, '-o')
    %plot(xx, fit, 'b--')
    plot(pref_orient, pref_resp, 'rs', 'markersize', 13, 'markerfacecolor', 'r')
    plot(null_orient, null_resp, 'bs', 'markersize', 13, 'markerfacecolor', 'b')
    xlabel('Orientation')
    ylabel('Mean Firing Rate')
    xlim([0, 2*pi])
    
end



% Plot the histogram of DS
l_valid = ~commonExclusions & ~isnan(DS_raw);
figure
hist(DS_raw(l_valid))
set(gca, 'fontsize', 14, 'box', 'off', 'tickDir', 'out')
xlabel('Direction Selectivity Index')
ylabel('Count')


% Plot CSIs vs. DS
figure
hold on,
plot(rawCSIs(l_valid), DS_raw(l_valid), 'bo')
set(gca, 'xscale', 'log', 'fontsize', 14)
ylabel('Direction Selectivity Index')
xlabel('Color Sensitivity Index')
[r,p] = corr(rawCSIs(l_valid), DS_raw(l_valid), 'type', 'spearman');
t = text(0.02, .8, sprintf('r = %.3f\np = %.3f', r,p));
set(t, 'fontSize', 14)


% now for TRs
minTRs = min(rawTRs,[],2);
l_nan = isnan(minTRs);
minTRs(l_nan) = max(minTRs).*1.2;
figure,
plot(minTRs(l_valid), DS_raw(l_valid), 'bo')
set(gca, 'xscale', 'log', 'fontsize', 14)
xlabel('Threshold Ratio')
ylabel('Direction Selectivity Index')
[r,p] = corr(minTRs(l_valid), DS_raw(l_valid), 'type', 'spearman');
[r_nonan, p_nonan] = corr(minTRs(l_valid & ~l_nan), DS_raw(l_valid & ~l_nan), 'type', 'spearman')
t = text(0.8, .8, sprintf('r = %.5f\np = %.3f', r,p));
set(t, 'fontSize', 14)



