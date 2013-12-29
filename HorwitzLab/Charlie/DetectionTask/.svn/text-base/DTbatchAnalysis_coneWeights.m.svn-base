%% (16) CONE WEIGHTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
global cardVsIntBatchPath blpBatchPath
preprocessDTbatchData

[l_SwM, l_SwL, l_S, l_LvM, l_other] = deal(logical(zeros(length(out.dat),1)));
for a = 1:length(out.dat)
    tmp = out.dat(a).grating.color.prefcolor;
    tmp = tmp ./ sum(abs(tmp)); %normalize the weights
    %make sure that the points all fall in the upper part of the diamond
    if sum(sign(tmp(1:2))) == 0 %i.e., L-M or -L+M
        if tmp(1) > 0;
            tmp = tmp.*-1;
        end
    elseif abs(sum(sign(tmp(1:2)))) == 2 %i.e., L+M or -L-M
        if tmp(1) < 0;
            tmp = tmp.*-1;
        end
    end
    rawWeights(a,:) = tmp;
    
    %determine which pref color goes with which cone weights
    prefIsolum = out.dat(a).prefIsolum;
    if all(sign(prefIsolum) == [0 0 1])
        l_S(a) = 1;
    elseif all(sign(prefIsolum) == [1 -1 0])
        l_LvM(a) = 1;
    elseif all(sign(prefIsolum) == [1 -1 1])
        l_SwL(a) = 1;
    elseif all(sign(prefIsolum) == [1 -1 -1])
        l_SwM(a) = 1;
    else
        l_other(a) = 1;
    end
end

%plot all the data
negS = rawWeights(:,3) < 0;
figure, hold on,
plot(rawWeights(~negS,1), rawWeights(~negS,2), 'bo')
plot(rawWeights(negS,1), rawWeights(negS,2), 'bo', 'markerfacecolor', 'b')
plot([0, 1], [1, 0], 'k')
plot([0, -1], [1, 0], 'k')
plot([0 0], [1 0], 'k')
plot([1 -1], [0 0], 'k')
plot([-.5 .5], [.5 .5], 'k+', 'markersize', 15)
axis equal
title('Normalized Cone Weights')
hold off

%now sized by TR, the biger the marker size, the lower the TR (i.e., more
%sensitive at threshold
MINSIZE = 8;
minTRs = min(rawTRs, [], 2);
markSize = (1./(minTRs./max(minTRs))).* 5;
figure, hold on,
title('Sized by 1/TR')
counter = 0;
for a = 1:length(minTRs);
    
    %filter out expts that don't meet certain criteria
    if commonExclusions(a)
        continue
    end
    
    if ~isnan(minTRs(a))
        if rawWeights(a,3) < 0;
            plot(rawWeights(a,1), rawWeights(a,2), 'ro', 'markersize', markSize(a))
        else
            plot(rawWeights(a,1), rawWeights(a,2), 'bo', 'markersize', markSize(a))
        end
    else
        plot(rawWeights(a,1), rawWeights(a,2), 'ko', 'markerfacecolor', 'k', 'markersize', max([1,MINSIZE-2]))
    end
end
plot([0, 1], [1, 0], 'k')
plot([0, -1], [1, 0], 'k')
plot([0 0], [1 0], 'k')
plot([1 -1], [0 0], 'k')
plot([-.5 .5], [.5 .5], 'k+', 'markersize', 15)
axis equal
hold off

%now according to pref isolum color
figure, hold on,
plot(rawWeights(l_S,1), rawWeights(l_S,2), 'ko', 'markerfacecolor', 'b');
plot(rawWeights(l_LvM,1), rawWeights(l_LvM,2), 'ko', 'markerfacecolor', 'r');
plot(rawWeights(l_SwL,1), rawWeights(l_SwL,2), 'ko', 'markerfacecolor', 'm');
plot(rawWeights(l_SwM,1), rawWeights(l_SwM,2), 'ko', 'markerfacecolor', 'g');
plot([0, 1], [1, 0], 'k')
plot([0, -1], [1, 0], 'k')
plot([0 0], [1 0], 'k')
plot([1 -1], [0 0], 'k')
plot([-.5 .5], [.5 .5], 'k+', 'markersize', 15)
axis equal
hold off
title('Pref Isolum Color')

%now according to CSI
l_CO = rawCSIs(:) >= 2;%color only
l_CL = (rawCSIs(:) >=0.5) & (rawCSIs(:) < 2); %color luminance
l_LO = rawCSIs(:) < 0.5;
figure, hold on,
plot(rawWeights(l_CO,1), rawWeights(l_CO,2), 'ko', 'markerfacecolor', 'r');
plot(rawWeights(l_CL,1), rawWeights(l_CL,2), 'ko', 'markerfacecolor', 'y');
plot(rawWeights(l_LO,1), rawWeights(l_LO,2), 'ko', 'markerfacecolor', 'k');
plot([0, 1], [1, 0], 'k')
plot([0, -1], [1, 0], 'k')
plot([0 0], [1 0], 'k')
plot([1 -1], [0 0], 'k')
plot([-.5 .5], [.5 .5], 'k+', 'markersize', 15)
axis equal
hold off
title('Pref Isolum Color')

% bar chart depicting the percentage of CSI types for each pref isolum color
figure, hold on,
bar(1,sum([l_CO & l_S])./sum(l_S),0.9, 'r')
bar(2,sum([l_CL & l_S])./sum(l_S),0.9, 'y')
bar(3,sum([l_LO & l_S])./sum(l_S),0.9, 'k')
bar(5,sum([l_CO & l_SwL])./sum(l_SwL),0.9, 'r')
bar(6,sum([l_CL & l_SwL])./sum(l_SwL),0.9, 'y')
bar(7,sum([l_LO & l_SwL])./sum(l_SwL),0.9, 'k')
bar(9,sum([l_CO & l_SwM])./sum(l_SwM),0.9, 'r')
bar(10,sum([l_CL & l_SwM])./sum(l_SwM),0.9, 'y')
bar(11,sum([l_LO & l_SwM])./sum(l_SwM),0.9, 'k')
bar(13,sum([l_CO & l_LvM])./sum(l_LvM),0.9, 'r')
bar(14,sum([l_CL & l_LvM])./sum(l_LvM),0.9, 'y')
bar(15,sum([l_LO & l_LvM])./sum(l_LvM),0.9, 'k')
legend('Color Only', 'Color-Lum', 'Lum-Only')
set(gca, 'xtick', [2 6 10 14], 'XTickLabel', {'S-iso', 'SwL', 'SwM', 'L-M'})
ylim([0,1])
title('Cell Type by Prefered Color')



% bar chart depicting the percentage of color types for each CSI type
figure, hold on,
bar(1,sum([l_CO & l_S])./sum(l_CO),0.9, 'b')
bar(2,sum([l_CO & l_SwL])./sum(l_CO),0.9, 'm')
bar(3,sum([l_CO & l_SwM])./sum(l_CO),0.9, 'g')
bar(4,sum([l_CO & l_LvM])./sum(l_CO),0.9, 'r')

bar(6,sum([l_CL & l_S])./sum(l_CL),0.9, 'b')
bar(7,sum([l_CL & l_SwL])./sum(l_CL),0.9, 'm')
bar(8,sum([l_CL & l_SwM])./sum(l_CL),0.9, 'g')
bar(9,sum([l_CL & l_LvM])./sum(l_CL),0.9, 'r')

bar(11,sum([l_LO & l_S])./sum(l_LO),0.9, 'b')
bar(12,sum([l_LO & l_SwL])./sum(l_LO),0.9, 'm')
bar(13,sum([l_LO & l_SwM])./sum(l_LO),0.9, 'g')
bar(14,sum([l_LO & l_LvM])./sum(l_LO),0.9, 'r')
legend('S-iso', 'SwL', 'SwM', 'L-M')
set(gca, 'xtick', [2.5 7.5 12.5], 'XTickLabel', {'Color Only', 'Color Lum', 'Lum Only'})
ylim([0,1])
title('Pref Color by CSI Type')

%% (16.3) CONE WEIGHTS USING SOME NEW CODE
%implementing a gradient descent with different initial conditions and
%representations of contrasts.... The fact that cone weights change with
%initial conditions is unsettling, but there's no obvious way to solve the
%problem using the normal equations...

clear, clc
global cardVsIntBatchPath
preprocessDTbatchData

nGT = length(out.fnames);
oldWts = nan(nGT, 3);
newWts = nan(nGT, 3);
newWts_CC = nan(nGT, 3);
for a = 1:nGT;
    GT = gtobj(out.fnames{a}{2});
    LMS = GT.trial(:,[GT.idx.lcc | GT.idx.mcc | GT.idx.scc]);
    l_p4 = GT.trial(:, GT.idx.protocol)==4;
    stimOn = GT.trial(:, GT.idx.stimon);
    frameRate = GT.sum.exptParams.framerate;
    nFrames = GT.trial(:, GT.idx.nFrames);
    stimOff = stimOn + (nFrames/frameRate);
    
    %convert to cell arrays and compute the cone weights
    t_on = mat2cell(stimOn, ones(length(stimOn),1));
    t_off = mat2cell(stimOff, ones(length(stimOff),1));
    counts = cellfun(@(on, off, ras) sum((ras>on)&(ras<=off)), t_on, t_off, GT.ras(:,GT.idx.spikes));
    
    %run the regression to determine cone weights (use the counts as the
    %response, and the actual color directions as the predictors. This
    %gives you the color direction which when dot producted against the
    %actual colors best predicts the actual responses)
    response = counts(l_p4);
    predictors = LMS(l_p4,:);
    %newWts(a,:) = predictors \ response;
    [~, maxidx] = min(response);
    initguess = predictors(maxidx,:);
    options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-6,'TolX',10^-6);
    [coneweights, fval, exitflag] = fminsearch(@(x) linmodfiterr(predictors, response, x), initguess, options);
    if (exitflag)
        newWts(a,:) = coneweights;
    else
        newWts(a,:) = [nan nan nan];
    end
    oldWts(a,:) = out.dat(a).grating.color.prefcolor;
    
    % re-run the regression using unit vectors for color directions and
    % counts/contrast unit as the response.
    norms = sqrt(sum(LMS.^2,2));
    response = counts(l_p4)./norms(l_p4,:);
    predictors = LMS(l_p4,:)./repmat(norms(l_p4,:), 1,3);
    [~, maxidx] = min(response);
    initguess = predictors(maxidx,:);
    options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-6,'TolX',10^-6);
    [coneweights, fval, exitflag] = fminsearch(@(x) linmodfiterr(predictors, response, x), initguess, options);
    if (exitflag)
        newWts_CC(a,:) = coneweights;
    else
        newWts_CC(a,:) = [nan nan nan];
    end
    %newWts_CC(a,:) = predictors \ response;
end

%convert to normalized rates, then make all the M-weights positive
normNew = newWts ./ repmat(sum(abs(newWts),2), 1, 3);
l_negM = normNew(:,2)<0;
normNew(l_negM,:) = normNew(l_negM,:)*-1;
normOld = oldWts ./ repmat(sum(abs(oldWts),2), 1, 3);
l_negM = normOld(:,2)<0;
normOld(l_negM,:) = normOld(l_negM,:)*-1;
normNewCC = newWts_CC ./ repmat(sum(abs(newWts_CC),2), 1, 3);
l_negM = normNewCC(:,2)<0;
normNewCC(l_negM,:) = normNewCC(l_negM,:)*-1;


figure
subplot(1,3,1), hold on,
plot(normNew(:,1), normNew(:,2), 'ko');
plot([0, 1], [1, 0], 'k')
plot([0, -1], [1, 0], 'k')
plot([0 0], [1 0], 'k')
plot([1 -1], [0 0], 'k')
plot([-.5 .5], [.5 .5], 'k+', 'markersize', 15)
ylim([0,1])
hold off
subplot(1,3,2), hold on,
plot(normOld(:,1), normOld(:,2), 'ko');
plot([0, 1], [1, 0], 'k')
plot([0, -1], [1, 0], 'k')
plot([0 0], [1 0], 'k')
plot([1 -1], [0 0], 'k')
plot([-.5 .5], [.5 .5], 'k+', 'markersize', 15)
ylim([0,1])
hold off
subplot(1,3,3), hold on,
plot(normNewCC(:,1), normNewCC(:,2), 'ko');
plot([0, 1], [1, 0], 'k')
plot([0, -1], [1, 0], 'k')
plot([0 0], [1 0], 'k')
plot([1 -1], [0 0], 'k')
plot([-.5 .5], [.5 .5], 'k+', 'markersize', 15)
ylim([0,1])
hold off

%% (16.4) CONE WEIGHTS IN THE L/M PLANE VS TRs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, clc, close all
global cardVsIntBatchPath blpBatchPath
preprocessDTbatchData

PLOT = 1;
rawWeights = nan(numel(out.dat), 3);
thetas_old = nan(numel(out.dat), 1); %from getGRatingsTuning
thetas_exp = nan(numel(out.dat), 1); %from the raisedExp fitting
raisedExp = nan(numel(out.dat),1);
for a = 1:numel(out.dat)
    %find the pref color using the fitting routine from getGratingTuning
    tmp = out.dat(a).grating.color.prefcolor;
    tmp = tmp./norm(tmp); %convert to unit vector.
    if sum(sign(tmp(1:2))) == 0 %i.e., L-M or -L+M
        if tmp(1) > 0;
            tmp = tmp.*-1;
        end
    elseif abs(sum(sign(tmp(1:2)))) == 2 %i.e., L+M or -L-M
        if tmp(1) < 0;
            tmp = tmp.*-1;
        end
    end
    rawWeights(a,:) = tmp;
    tmp = atan(tmp(2)./tmp(1));
    thetas_old(a) = mod(tmp, pi);
    
    %now estimate the prefColor from the raisedExp routine (in the L/M
    %plane), this will also help identify the cases where color tuning is
    %quite broad (and thus not particularly useful for this analysis)
    l_LMplane = out.dat(a).grating.color.colors(:,3) == 0;
    LMresp = out.dat(a).grating.color.colresp(l_LMplane,1);
    LMsem = out.dat(a).grating.color.colresp(l_LMplane,2)./sqrt(out.dat(a).grating.color.colresp(l_LMplane,3));
    LMcolors = out.dat(a).grating.color.colors(l_LMplane,:);
    Lcc = LMcolors(:,1);
    Mcc = LMcolors(:,2);
    LMthetas = atan(Mcc./Lcc);
    [tmp_theta, raisedExp(a), gain] = raisedCos(LMthetas, LMresp);
    thetas_exp(a) = mod(tmp_theta, pi);
    
    if PLOT && thetas_exp(a)>(pi/2) && ~commonExclusions(a) && all(sign(prefCards(a,:)) == [1 -1 0])
        figure
        plotTheta = [LMthetas; LMthetas+pi]';
        plotLMResp = [LMresp; LMresp]';
        perr = polar([plotTheta; plotTheta], [plotLMResp-[LMsem', LMsem']; plotLMResp+[LMsem', LMsem']], 'k');
        hold on
        p1 = polar(plotTheta, plotLMResp, 'ko');
        set(p1, 'markerfacecolor', 'k', 'markersize', 10)
        set(perr, 'linewidth', 3)
        
        
        %add the fits and prefColor estimates
        t = 0:0.001:(2*pi);
        pred = gain .* abs(cos(t - thetas_exp(a))).^raisedExp(a);
        p2 = polar(t, pred, 'b');
        set(p2, 'linewidth', 2)
        p3 = polar(thetas_exp(a), gain, 'gx');
        p4 = polar(thetas_old(a), gain, 'rx');
        set([p3, p4], 'markersize', 13, 'linewidth', 3)
        hold off
        
        %annotate the figure.
        title(sprintf('%s exp: %.3f TR: %.3f', out.fnames{a}{1}, raisedExp(a), rawCardTRs(a)));
    end
end

%compile the sensitivities to the L-M color direction as assesed by DTspot
l_lvm = ismember(sign(prefCards), [1 -1 0], 'rows');
cardSensitivities = 1./ rawCardTRs;
valForNans = min(cardSensitivities) .* .5;
filtSensitivities = cardSensitivities;
filtSensitivities(isnan(filtSensitivities)) = valForNans;

%now plot the TR's as a function of preferred color in the L/M plane
figure
p1 = polar(thetas_exp(l_lvm), filtSensitivities(l_lvm), 'ko');
hold on,
polar([thetas_exp(l_lvm), thetas_old(l_lvm)]', [filtSensitivities(l_lvm), filtSensitivities(l_lvm)]', 'k')
p2 = polar(thetas_old(l_lvm), filtSensitivities(l_lvm), 'ro');
set([p1, p2], 'markersize', 6, 'linewidth', 2)
x = 0:0.001:(2.*pi);
p3 = polar(x, abs(cos(x-(3*pi/4))), 'b');
set(p3, 'linewidth', 2)
title('Old (red) and Exp (black) prefColors')


%now bin the data to get a better estimate of the average TR as a function
%of pref color (first for the l-m color direction)
expCutoff = 0; %i.e less than this won't get counted...
edges = linspace(0,pi,7);
[~, bin] = histc(thetas_exp, edges);
avgSen = [];
SEM = [];
N = [];
for a = 1:max(bin);
    l_bin = (bin == a) & ~commonExclusions;
    l_expCrit = raisedExp > expCutoff;
    
    %first without excluding on the basis of exponent
    tmp = cardSensitivities(l_lvm & l_bin);
    avgSen(a,1) = nanmean(tmp);
    N(a,1) = sum(~isnan(tmp));
    SEM(a,1) = nanstd(tmp)./sqrt(N(a,1));
    
    
    %now culling data on the basis of exponent
    tmp = cardSensitivities(l_lvm & l_bin & l_expCrit);
    avgSen(a,2) = nanmean(tmp);
    N(a,2) = sum(~isnan(tmp));
    SEM(a,2) = nanstd(tmp)./sqrt(N(a,2));
end
N
SEM
figure
plotEdges = edges(2:end)-(diff(edges)/2);
p2 = polar([plotEdges; plotEdges], [(avgSen(:,2)+SEM(:,2))'; (avgSen(:,2)-SEM(:,2))'], 'r');
hold on,
p1 = polar([plotEdges; plotEdges], [(avgSen(:,1)+SEM(:,1))'; (avgSen(:,1)-SEM(:,1))'], 'k');
set([p1, p2], 'linewidth', 2)
p4 = polar(plotEdges, avgSen(:,2)', 'ro');
p3 = polar(plotEdges, avgSen(:,1)', 'ko');
set(p4, 'linewidth', 2, 'markersize', 8, 'markerfacecolor', 'r')
set(p3, 'linewidth', 2, 'markersize', 8, 'markerfacecolor', 'k')
x = 0:0.001:(2*pi);
modScaleFactor = avgSen(5,1);
p5 = polar(x, abs(cos(x-(3*pi/4)))*modScaleFactor, 'b');
set(p5, 'linewidth', 2)
hold off
title('TRs to L-M stimuli as a fxn of prefColor')


%now bin the data to get a better estimate of the average TR as a function
%of pref color (now for the SwM intermediate)
l_SwM = ismember(sign(prefInts), [1 -1 -1], 'rows');
intSensitivities = 1./rawIntTRs;
avgSen = [];
SEM = [];
N = [];
for a = 1:max(bin);
    l_bin = bin == a;
    l_expCrit = raisedExp > expCutoff;
    
    %first without excluding on the basis of exponent
    tmp = intSensitivities(l_SwM & l_bin);
    avgSen(a,1) = nanmean(tmp);
    N(a,1) = sum(~isnan(tmp));
    SEM(a,1) = nanstd(tmp)./sqrt(N(a,1));
    
    %now culling data on the basis of exponent
    tmp = intSensitivities(l_SwM & l_bin & l_expCrit);
    avgSen(a,2) = nanmean(tmp);
    N(a,2) = sum(~isnan(tmp));
    SEM(a,2) = nanstd(tmp)./sqrt(N(a,2));
end
N
SEM
%Generating a linear prediction on the assumption that the neuron is
%"tunend" to a direction in the L/M plane, but always has the same S-cone
%component as the SwM intermediate
SwMunit = [-0.1386   0.1386   0.9806]; %notice that it's adjusted for the 2nd quadrant of the L/M plane
C = sqrt(2*0.1386^2);
L = linspace(C, -C, 100)'; %space from C to -C, which corresponds to all L to all -L
L(end) = [];
M = sqrt(C^2 - L.^2);
synthVecs = [L(:), M(:), repmat(SwMunit(3), numel(L), 1)];
%sqrt(sum(synthVecs.^2,2)) %should be a bunch of ones
synthPred = synthVecs * SwMunit(:);
synthTheta = mod(atan(M./L),pi);

figure
plotEdges = edges(2:end)-(diff(edges)/2);
p1 = polar([plotEdges; plotEdges], [(avgSen(:,1)+SEM(:,1))'; (avgSen(:,1)-SEM(:,1))'], 'k');
hold on,
p2 = polar([plotEdges; plotEdges], [(avgSen(:,2)+SEM(:,2))'; (avgSen(:,2)-SEM(:,2))'], 'r');
set([p1, p2], 'linewidth', 2)
p3 = polar(plotEdges, avgSen(:,1)', 'ko');
p4 = polar(plotEdges, avgSen(:,2)', 'ro');
set(p3, 'linewidth', 2, 'markersize', 8, 'markerfacecolor', 'k')
set(p4, 'linewidth', 2, 'markersize', 8, 'markerfacecolor', 'r')
x = 0:0.001:(2*pi);
p5 = polar(x, abs(cos(x-(3*pi/4)))*.8, 'b'); %assuming an L-M mech
p6 = polar(synthTheta, synthPred.*.8, 'b:'); %assuming S-cone input
set([p5, p6], 'linewidth', 2)
hold off
title('TRs to SwM stimuli as a fxn of prefColor')

