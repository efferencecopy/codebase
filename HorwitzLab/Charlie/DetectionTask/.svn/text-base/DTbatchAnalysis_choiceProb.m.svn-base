%% (18) CHOICE PROBABILITY: ZERO CONTRAST TRIALS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
global cardVsIntBatchPath blpBatchPath
preprocessDTbatchData

for a = 1:length(out.dat);
    rawBkgndCP(a) = out.dat(a).c.cp.in{1}(1);
end
tmp_CSIs = rawCSIs; %change this to switch between CSI definitions

l_validConds = ~(commonExclusions);% | (out.errors(:, neuroThresh1Ind) & out.errors(:, neuroThresh2Ind)));
filt_bkgndCP = rawBkgndCP(l_validConds);
avgCP = nanmean(filt_bkgndCP);
figure,
subplot(1,4,1)% all cells
hold on,
hist(filt_bkgndCP,10)
plot(avgCP, 2, 'm*')
xlabel('All Cells')
title(sprintf('avg CP = %.3f', avgCP))
axis tight
hold off

subplot(1,4,2) %color-luminance cells
hold on,
l_validConds = ~commonExclusions(:) & (tmp_CSIs>0.5) & (tmp_CSIs<2);
filtCP = rawBkgndCP(l_validConds);
hist(filtCP, 7)
avgCP = nanmean(filtCP);
plot(avgCP, 0.5, 'm*')
axis tight
xlabel('Color Luminance Cells')
title(sprintf('avg CP = %.3f', avgCP))
hold off

subplot(1,4,3) %color only cells
hold on,
l_validConds = ~commonExclusions & (tmp_CSIs>2);
filtCP = rawBkgndCP(l_validConds);
hist(filtCP, 7)
avgCP = nanmean(filtCP);
plot(avgCP, 0.5, 'm*')
axis tight
xlabel('Color-Only Cells')
title(sprintf('avg CP = %.3f', avgCP))
hold off

subplot(1,4,4) %lum only cells
hold on,
l_validConds = ~commonExclusions & (tmp_CSIs<0.5);
filtCP = rawBkgndCP(l_validConds);
hist(filtCP, 7)
avgCP = nanmean(filtCP);
plot(avgCP, 0.5, 'm*')
axis tight
xlabel('Luminance-Only Cells')
title(sprintf('avg CP = %.3f', avgCP))
hold off

%% (19) CHOICE PROBABILITY AS A FUNCTION OF CONTRAST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
global cardVsIntBatchPath blpBatchPath
preprocessDTbatchData


for a = 1:length(out.dat)
    inCP_clr1(a,:) = out.dat(a).c.cp.in{1};
    inCP_clr2(a,:) = out.dat(a).c.cp.in{2};
    outCP_clr1(a,:) = out.dat(a).c.cp.out{1};
    outCP_clr2(a,:) = out.dat(a).c.cp.out{2};
end
inCP = vertcat(inCP_clr1, inCP_clr2);
outCP = vertcat(outCP_clr1, outCP_clr2);

figure
hold on,
plot(repmat((1:8),size(outCP,1), 1)', outCP', 'ob')
plot(repmat((1:8),size(inCP,1), 1)', inCP', 'or')
plot([1, 8], [0.5, 0.5], 'k:')
plot((1:8), nanmean(outCP), 'sc', 'markersize', 8,'markerfacecolor', 'c')
plot((1:8), nanmean(inCP), 'sm', 'markersize', 8,'markerfacecolor', 'm')
hold off,
xlim([0.5 8.5])
l = legend('out of RF', 'in RF');
c = get(l, 'children');
set(c(1), 'color', 'r', 'markerfacecolor', 'r')
set(c(4), 'color', 'b', 'markerfacecolor', 'b')
xlabel('Contrast Level')
ylabel('Choice Probability')
title('CP for all colors/contrasts')



%% (20) CHOICE PROBABILITY USING SOME NEW CODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, clc
global cardVsIntBatchPath blpBatchPath
preprocessDTbatchData

batchProcess =0;
permuteChoices = 1;
if batchProcess
    for a = 1:length(out.dat)
        out.fnames{a}{1}
        cp(a) = dtcp(dtobj(out.fnames{a}{1}), permuteChoices);
    end
    cp_dtcp = [cp(:).zScoreIn];
    p_dtcp = [cp(:).zScoreIn_P];
else
    if isunix
        load '/Users/charliehass/LabStuff/NexFiles/Charlie/Batch Data And Text Files/dtcp_batch_May-02.mat'
    elseif ispc
        error('add a path for pc')
    end
    cp_dtcp = [cp(:).zScoreIn];
    p_dtcp = [cp(:).zScoreIn_P];
end


figure, hold on,
l_valid = ~(commonExclusions | (out.errors(:, neuroThresh1Ind) & out.errors(:, neuroThresh2Ind)));
plot([0.5 0.5], [0, 42], '--k', 'linewidth', 2)
hist(cp_dtcp(l_valid))
h = findobj(gca,'Type','patch');
set(h,'FaceColor','k','EdgeColor','w')
plot(nanmean(cp_dtcp(l_valid)),42, 'kv', 'markerfacecolor', 'k', 'markersize', 15)
[p,h] = signtest(cp_dtcp(l_valid)-0.5);
title(sprintf('DTCP p=%.3f', p));
set(gca, 'fontsize', 20, 'linewidth', 2, 'tickDir', 'out')
ylim([0, 42])
ylabel('Count')
xlabel('Pooled Choice Probability')
xlim([.3 .8])

for a = 1:length(out.dat)
    if isfield(out.dat(a).c.cp, 'grandCP')
        cp_DTunpack(a) = out.dat(a).c.cp.grandCP.val;
        p_DTunpack(a) = min(1, out.dat(a).c.cp.grandCP.p); %coverts the nans to ones
    else
        cp_DTunpack(a) = nan;
        p_DTunpack(a) = 1;
    end
end
figure, hold on,
hist(cp_DTunpack(l_valid))
plot(nanmean(cp_DTunpack(l_valid)),4, 'mv')
plot([0.5 0.5], [0 18], 'w')
[p, h] = signtest(cp_DTunpack(l_valid)-0.5);
title(sprintf('DTunpack p=%.3f', p))


%as a funtion of contrast (new code vs old code)
cntrst_dtcp = [];
cntrst_DTunpack = [];
for a = 1:length(out.dat)
    for i = 1:2
        cntrst_dtcp = [cntrst_dtcp; cp(a).cp(i,:)'];
        cntrst_DTunpack = [cntrst_DTunpack; out.dat(a).c.cp.in{i}'];
    end
end
figure, hold on,
plot(cntrst_dtcp, cntrst_DTunpack, '.')
plot([0 1], [0 1], 'k-')
title('DTCP vs. DTmocsUnpack')
axis square

%TR vs CP
trs = rawTRs(l_valid,:);
trs(isnan(trs)) = max(trs(:))*1.1;
grandCPs = cp_dtcp(l_valid);
l_sig = p_dtcp(l_valid)<0.05;
figure, hold on,
plot(trs, grandCPs, 'o') %plotting both TR estimates
plot(trs(l_sig), grandCPs(l_sig), 'ro', 'markerfacecolor', 'r')
set(gca, 'xscale', 'log')
ylabel('choice probability')
xlabel('Threshold Ratio')
[~, p] = corr(trs(~isnan(grandCPs'),:), grandCPs(~isnan(grandCPs))', 'type', 'spearman');
title(sprintf('Spearman''s p=%.3f', mean(p)))
hold off

%% (20) PLOTTING CP BY COLOR

clear, clc, close all
global cardVsIntBatchPath blpBatchPath
preprocessDTbatchData

%assign the CP measurements by color direction tested
[s, lvm, swm, swl] = deal([]);
for a = 1:numel(out.dat)
    colors = out.dat(a).expt.standColors;
    idx_s = ismember(sign(colors), [0 0 1], 'rows');
    idx_lvm = ismember(sign(colors), [1 -1 0], 'rows');
    idx_swm = ismember(sign(colors), [1 -1 -1], 'rows');
    idx_swl = ismember(sign(colors), [1 -1 1], 'rows');
    if ~any(idx_s | idx_lvm | idx_swm | idx_swl)
        error('Something is messed up')
    end
    
    if any(idx_s)
        s = [s, out.dat(a).c.cp.poolConIn.val(idx_s)];
    end
    if any(idx_lvm)
        lvm = [lvm, out.dat(a).c.cp.poolConIn.val(idx_lvm)];
    end
    if any(idx_swm)
        swm = [swm, out.dat(a).c.cp.poolConIn.val(idx_swm)];
    end
    if any(idx_swl)
        swl = [swl, out.dat(a).c.cp.poolConIn.val(idx_swl)];
    end
end


figure, 
hold on,
plot(1, s, 'bo')
plot(1, nanmean(s), 'bs', 'markersize', 10, 'markerfacecolor', 'b')
plot(2, lvm, 'ro')
plot(2, nanmean(lvm), 'rs', 'markersize', 10, 'markerfacecolor', 'r')
plot(3, swm, 'go')
plot(3, nanmean(swm), 'gs', 'markersize', 10, 'markerfacecolor', 'g')
plot(4, swl, 'mo')
plot(4, nanmean(swl), 'ms', 'markersize', 10, 'markerfacecolor', 'm')




%% (21) GRAND CP SEPERATED BY CSI

clear, clc, close all
global cardVsIntBatchPath blpBatchPath
preprocessDTbatchData

[cp, p] = deal(nan(length(out.dat),1));
for a = 1:length(out.dat)
    if isfield(out.dat(a).c.cp, 'grandCP')
        cp(a) = out.dat(a).c.cp.grandCP.val;
        p(a) = min(1, out.dat(a).c.cp.grandCP.p); %coverts the nans to ones
    end
end

l_valid = sum(out.errors(:, [orientMismatchInd, posMismatchInd]),2) == 0;
l_sig = p<0.05;
figure, hold on,
set(gca, 'fontsize', 18)
plot(rawCSIs(l_valid), cp(l_valid), 'ko', 'markerfacecolor', 'k')
%plot(rawCSIs(l_valid&l_sig), cp(l_valid&l_sig), 'ro', 'markerfacecolor', 'r')
set(gca, 'xscale', 'log')
xlabel('CSI')
ylabel('Grand CP')
corrConds = l_valid & ~isnan(cp);
[rho, p_rho] = corr(rawCSIs(corrConds), cp(corrConds), 'type', 'spearman');
title(sprintf('rho: %.3f, p:%.3f', rho, p_rho))


