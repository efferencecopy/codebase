%% (29) BANDLIMITED POWER LOCAL FIELD POTENTIALS

clear, clc, close all
global cardVsIntBatchPath blpBatchPath
preprocessDTbatchData

%cycle through the DT files and compute the CP on the bases of BLP LFP
makePlot = 1;
for a = 1:length(out.dat);
    
    %make sure that the data come from expt using the new LFP board
    exptDate = out.fnames{a}{1}(2:7);
    month = str2num(exptDate(1:2));
    day = str2num(exptDate(3:4));
    year = str2num(exptDate(5:6));
    if (year<10) || ((year<=10) && (month<9)) || ((year<=10) && (month<=9) && (day<28))
        fprintf('Skipping <%s>\n', out.fnames{a}{1});
        continue
    end
    
    
    DT = dtobj(out.fnames{a}{1});
    lfpRasIdx = DT.idx.lfp;
    anlgSumIdx = strcmpi(DT.sum.analog.sigid, DT.sum.rasterCells{lfpRasIdx});
    
    
    %open a figure if need be
    if makePlot
        figure
        figName = sprintf('%s: avg TR: %.3f', out.fnames{a}{1}, nanmean(rawTRs(a,:)));
        set(gcf, 'name', figName, 'position', [-1262 39 1226 672]);
    end
    
    band = {[8,15], [20 45], [50 70], [80 95]}; %for cheby
    %band = {10, 30, 60, 85}; %for greg's thing
    for bnd = 1:length(band);
        
        %determine the band limited power
        blp = blplfp(DT, lfpRasIdx, band{bnd});
        l_oob = cellfun(@(x) any(isnan(x)), blp);
        if sum(l_oob) == length(l_oob); close, continue; end
        
        % create and avereage across the gabor presentation
        gaborOn = mat2cell(DT.trial(:,DT.idx.repFlashOn), ones(size(DT.trial,1),1));
        nFrames = mat2cell(DT.trial(:, DT.idx.nFrames), ones(size(DT.trial,1),1));
        frameRate = repmat({DT.sum.exptParams.frame_rate}, size(DT.trial,1),1);
        anlgRate = repmat({DT.sum.analog.storeRates{anlgSumIdx}}, size(blp));
        gaborOnSamp = cellfun(@(gbOn, anSt, anRt)  round((gbOn-anSt)*anRt), gaborOn, DT.ras(:, DT.idx.anlgStart), anlgRate, 'uniformoutput', 0);
        nSampsPerStim = repmat({round(0.666*DT.sum.analog.storeRates{anlgSumIdx})}, size(blp));
        gaborOffSamp = cellfun(@(gbOnSamp, n) gbOnSamp+n, gaborOnSamp, nSampsPerStim, 'uniformoutput', 0);
        
        %compute the LFP synched to gabor onset
        preSamp = repmat({200}, size(blp));
        gaborLFP = cellfun(@(x, b, pre, e) x(b-pre:e), blp, gaborOnSamp, preSamp, gaborOffSamp, 'uniformoutput', 0);
        gaborLFP(l_oob) = [];
        
        %compute CRFs for each color contrast condition
        trlAvg = cellfun(@(x, b, e) nanmean(x(b:e)), blp, gaborOnSamp, gaborOffSamp);
        l_inRF = (DT.trial(:, DT.idx.flashX) == DT.sum.exptParams.rf_x) & (DT.trial(:, DT.idx.flashY) == DT.sum.exptParams.rf_y);
        nCont = length(unique(DT.trial(:, DT.idx.cntrstLev)));
        xbar = nan(1,nCont);
        SEM = nan(1, nCont);
        for cntrst = 1:nCont
            l_cntrst = DT.trial(:, DT.idx.cntrstLev) == cntrst;
            tList = l_cntrst & l_inRF;
            xbar(cntrst) = nanmean(trlAvg(tList));
            SEM(cntrst) = nanstd(trlAvg(tList)) ./ sqrt(sum(~isnan(trlAvg(tList))));
        end
        
        %compute the CP for the BLP
        blpStruct{a}(bnd) = dtcp(DT, 0, trlAvg);
        if sum(isnan(trlAvg)) ~= sum(l_oob); keyboard; end
        
        if makePlot
            %plot the raw data as an image
            subplot(3, length(band), bnd)
            tmp = [gaborLFP{:}]';
            [~, mapping] = sort(DT.trial(~l_oob, DT.idx.cntrstLev));
            imagesc(tmp(mapping,:))%order the raw data by contrast (zero on top)
            set(gca, 'xticklabel', num2str(str2num(get(gca, 'xtickLabel'))-preSamp{1}))
            title(sprintf('band: %s', num2str(band{bnd})));
            colormap jet
            %plot the average blp synched to gabor onset
            subplot(3, length(band), bnd+length(band))
            errorbar(nanmean([gaborLFP{:}]'), nanstd([gaborLFP{:}]')./sqrt(length(gaborLFP)))
            xlim([0 866])
            set(gca, 'xticklabel', num2str(str2num(get(gca, 'xtickLabel'))-preSamp{1}))
            %plot the roc values for BLP and single units
            subplot(3,length(band), bnd+(length(band)*2))
            hold on,
            plot(blpStruct{a}(bnd).norms(:,2:end)', blpStruct{a}(bnd).roc(:,2:end)', 'o-')
            spike_roc = vertcat(out.dat(a).c.roc{:});
            spike_norms = vertcat(out.dat(a).expt.norms{:});
            plot(spike_norms(:,2:end)', spike_roc(:,2:end)', 'v:')
            set(gca, 'xscale', 'log')
            axis tight
            ylim([.4 1])
            xlabel(sprintf('CP: %.3f, p=%.3f', blpStruct{a}(bnd).zScoreIn, blpStruct{a}(bnd).zScoreIn_P));
            title(sprintf('nans: %d', sum(isnan(trlAvg))))
        end
    end
end

%define a structure to save:
lfp.bands = band;
lfp.dat = blpStruct;
str = ['BLP_LFP_' date];
eval(['save ' str ' lfp']);
presDir = pwd;
fprintf('**** Batch data saved to directory: %s ****\n', presDir);


%% (30) COMPARING SPIKES AND BLP-LFP
%
% Trying to determine if the area under ROC changes b/w spike and lfp.
% ** looks like the insensitive ones (in spikes) are more sensitive in lfp


clear, clc
global cardVsIntBatchPath blpBatchPath
preprocessDTbatchData

figure
diffOrRatio = 1; %1=>diff 0=>ratio
for bnd = 1:4
    diffs = [];
    for ex = 1:length(lfp.dat)
        if isempty(lfp.dat{ex}); continue; end %bail if there's no lfp data
        if commonExclusions(ex); continue; end %bail when there wasn't much data
        if ~(out.errors(ex, neuroThresh1Ind) && out.errors(ex, neuroThresh2Ind)); continue; end
        
        spike_roc = vertcat(out.dat(ex).c.roc{:});
        spike_norm = vertcat(out.dat(ex).expt.norms{:});
        lfp_roc = lfp.dat{ex}(bnd).roc;
        lfp_norm = lfp.dat{ex}(bnd).norms;
        
        %errorchecking:
        if any(abs(spike_norm-lfp_norm) > 0.01); keyboard; end %make sure I'm comparing the correct conditions
        
        %do the comparison
        if diffOrRatio
            tmp_diff = spike_roc - lfp_roc;
        else
            tmp_diff = spike_roc ./ lfp_roc;
        end
        diffs = [diffs; tmp_diff(:)];
    end
    
    %plot the results.
    subplot(1,4,bnd), hold on,
    if diffOrRatio
        hist(diffs, 15)
        plot(nanmean(diffs), 6, 'rv', 'markerfacecolor', 'r')
        plot([0 0], [0 max(hist(diffs, 15))], 'w-')
        [h,p] = ttest(diffs);
        xlabel('spike-lfp')
    else
        edges = logspace(log10(0.2), log10(3), 20);
        counts = histc(diffs, edges);
        bar(edges, counts, 'histc')
        set(gca, 'xscale', 'log','xlim', [0.2, 3])
        plot(nanmean(diffs), 6, 'rv', 'markerfacecolor', 'r')
        plot([1 1], [0 max(counts)], 'w-')
        [h,p] = ttest(diffs-1);
        xlabel('spike/lfp')
    end
    title(sprintf('band: %s p=%.3f', num2str(lfp.bands{bnd}), p))
end

%% (31) NTs FOR SPIKES VS. BLP LFP

clear, clc
global cardVsIntBatchPath blpBatchPath
preprocessDTbatchData

makePlot = 0;
nBands = length(lfp.bands);
[lfpAlpha{1:nBands}, spikeAlpha{1:nBands}] = deal(nan(length(lfp.dat),2));
[lfpBeta{1:nBands}, spikeBeta{1:nBands}] = deal(nan(length(lfp.dat),2));
for ex = 1:length(lfp.dat)
    if isempty(lfp.dat{ex}); continue; end %bail if there's no lfp data
    if commonExclusions(ex); continue; end %bail when there wasn't much data
    %if (out.errors(ex, neuroThresh1Ind) && out.errors(ex, neuroThresh2Ind)); continue; end
    
    if makePlot
        figure
        set(gcf, 'position', [-1142 305 1130 386], 'name', num2str(ex));
    end
    for bnd = 1:nBands;
        for clr = 1:2;
            %errorchecking:
            spike_norm = out.dat(ex).expt.norms{clr};
            lfp_norm = lfp.dat{ex}(bnd).norms(clr,:);
            if any(abs(spike_norm-lfp_norm) > 0.01); keyboard; end %make sure I'm comparing the correct conditions
            
            %compute the NT based on the LFP data (using the same process
            %as in DTmocsUnpack
            lfp_roc = lfp.dat{ex}(bnd).roc(clr,:);
            if any(isnan(lfp_roc)); continue; end
            nTrialsIn = out.dat(ex).c.nTrialsIn{clr};
            l_nans = isnan(lfp_roc);
            lfp_roc(l_nans) = [];
            nTrialsIn(l_nans) = [];
            lfp_norm(l_nans) = [];
            correctByContrast = (lfp_roc.*nTrialsIn);
            wrongByContrast = (nTrialsIn - correctByContrast);
            zeroList = lfp_norm == 0;
            tmpNorms = lfp_norm(~zeroList);
            errs = abs(0.82-lfp_roc(~zeroList));
            aGuess = tmpNorms(find(errs == min(errs), 1, 'last'));
            [aSSE, bSSE, ~, success(1)] = weibullFit(lfp_norm, lfp_roc, 'sse', [aGuess 1]);
            [tmpAlpha, tmpBeta, tmpGamma, success(2)] = weibullFit(lfp_norm, [correctByContrast(:), wrongByContrast(:)], 'mle', [aSSE, bSSE]);
            if ~all(success) || ~any(lfp_roc > 0.80);
                tmpAlpha = NaN;
                tmpBeta = NaN;
            end
            
            %keep track of the results
            lfpAlpha{bnd}(ex,clr) = tmpAlpha;
            lfpBeta{bnd}(ex, clr) = tmpBeta;
            spikeAlpha{bnd}(ex, clr) = out.dat(ex).c.alpha(clr);
            spikeBeta{bnd}(ex, clr) = out.dat(ex).c.beta(clr);
            
            if makePlot
                subplot(2, nBands, bnd+(nBands*(clr-1))), hold on,
                ccs = logspace(log10(0.001), log10(max(lfp_norm)*1.1), 100);
                sG = out.dat(ex).c.gamma(clr);
                sB = out.dat(ex).c.beta(clr);
                sA = out.dat(ex).c.alpha(clr);
                spikeMod = sG + (0.5 - sG).*exp(-((ccs./sA).^sB));
                lfpMod = tmpGamma + (0.5 - tmpGamma).*exp(-((ccs./tmpAlpha).^tmpBeta));
                plot(lfp_norm, lfp_roc, 'ro');
                plot(spike_norm, out.dat(ex).c.roc{clr}, 'ko');
                plot(ccs, lfpMod, 'r')
                plot(ccs, spikeMod, 'k')
                set(gca, 'xscale', 'log')
                ylim([0.3 1.03])
                xlim([lfp_norm(2)*0.98 lfp_norm(end)*1.07])
                title(sprintf('Spike = %.3f, LFP = %.3f', sA, tmpAlpha))
                xlabel('Cone Contrast')
                ylabel('Prob(correct)')
                if bnd+(nBands*(clr-1)) == 1
                    legend('LFP', 'Spikes')
                end
            end
        end
    end
end

figure
for bnd = 1:nBands
    subplot(2,nBands,bnd), hold on, %thresholds
    plot(lfpAlpha{bnd}', spikeAlpha{bnd}', 'bo-')
    maxAlpha = max([lfpAlpha{bnd}(:); spikeAlpha{bnd}(:)]) .* 1.1;
    plot([0 maxAlpha], [0 maxAlpha], 'k')
    axis tight
    xlabel('LFP')
    ylabel('SPIKES')
    title(sprintf('NT in band: [%s]', num2str(lfp.bands{bnd})))
    
    subplot(2,nBands,bnd+nBands), hold on, %slopes
    plot(lfpBeta{bnd}', spikeBeta{bnd}', 'bo-')
    maxBeta = 15;%max([lfpBeta{bnd}(:); spikeBeta{bnd}(:)]) .* 1.1
    plot([0 maxBeta], [0 maxBeta], 'k')
    xlim([0 maxBeta])
    ylim([0 maxBeta])
    xlabel('LFP')
    ylabel('SPIKES')
    title(sprintf('Slope in band: [%s]', num2str(lfp.bands{bnd})))
end

figure
for bnd = 1:nBands
    plt = bnd;
    subplot(2,nBands,plt), hold on, %thresholds
    edges = logspace(log10(0.01), log10(8), 20);
    counts = histc([lfpAlpha{bnd}(:)./spikeAlpha{bnd}(:)], edges);
    bar(edges, counts, 'type', 'histc')
    plot([1 1], [0 max(counts)], 'g')
    plot(nanmedian([lfpAlpha{bnd}(:)./spikeAlpha{bnd}(:)]), 1, 'rv', 'markerfacecolor', 'r')
    set(gca, 'xscale', 'log')
    xlim([min(edges(counts>0))*.8, max(edges(counts>0))*2])
    title(sprintf('Band: [%s]', num2str(lfp.bands{bnd})))
    if plt == 1;
        ylabel('thresholds')
    end
    
    plt = bnd+nBands;
    subplot(2,nBands,plt), hold on, %slopes
    counts = hist([lfpBeta{bnd}(:)./spikeBeta{bnd}(:)], edges);
    bar(edges, counts, 'type', 'histc')
    plot([1 1], [0 max(counts)], 'g')
    plot(nanmedian([lfpBeta{bnd}(:)./spikeBeta{bnd}(:)]), 1, 'rv', 'markerfacecolor', 'r')
    set(gca, 'xscale', 'log')
    xlim([min(edges(counts>0))*.8, max(edges(counts>0))*2])
    xlabel('LFP / SPIKES')
    if plt == nBands+1;
        ylabel('slope')
    end
end

%NTs vs. CSI
figure
for bnd = 1:nBands
    minNTs = min(lfpAlpha{bnd}, [],2);
    [~, p_rho] = corr(minNTs(~isnan(minNTs)), rawCSIs(~isnan(minNTs)), 'type', 'spearman');
    subplot(1,nBands, bnd)
    plot(rawCSIs, minNTs, 'bo')
    set(gca, 'xscale', 'log')
    title(sprintf('Band: [%s], p=%.3f', num2str(lfp.bands{bnd}), p_rho))
    xlabel('CSI')
    ylabel('LFP NT')
end

%% (32) BAND LIMITED POWER CP
%
% this analysis is pretty busted. I'm not sure why there's so little data
% in the BLP CP
%
clear, clc
global cardVsIntBatchPath blpBatchPath
preprocessDTbatchData

figure
spike_cp = nan(length(out.dat),1);
cp = nan(length(lfp.dat),1);
p = nan(length(lfp.dat),1);
for bnd = 1:4
    for ex = 1:length(lfp.dat)
        if isempty(lfp.dat{ex}); continue; end %bail if there's no lfp data
        if commonExclusions(ex); continue; end %bail when there wasn't much data
        if (out.errors(ex, neuroThresh1Ind) && out.errors(ex, neuroThresh2Ind)); continue; end
        
        lfp_cp = lfp.dat{ex}(bnd).zScoreIn;
        lfp_p = lfp.dat{ex}(bnd).zScoreIn_P;
        
        cp(ex) = lfp_cp;
        p(ex) = lfp_p;
        spike_cp(ex) = out.dat(ex).c.cp.grandCP.val;
    end
    
    %plot the results.
    set(gcf, 'name', 'spike/lfp')
    subplot(1,4,bnd), hold on,
    hist(cp,7)
    plot(nanmean(cp), 6, 'rv', 'markerfacecolor', 'r')
    plot([0.5 0.5], [0 max(hist(cp))], 'w-')
    [h,pttest] = ttest(cp-0.5);
    title(sprintf('band: %s p=%.3f', num2str(lfp.bands{bnd}), pttest))
end


%% (33) BLP in gratings S vs. L-M

clear, clc
global cardVsIntBatchPath blpBatchPath
preprocessDTbatchData

%cycle through the GT files
makePlot = 1;
[sSpike lvmSpike sLFP{1:4} lvmLFP{1:4}] = deal(nan(length(out.fnames),1));
for a = 1:length(out.fnames);
    
    %make sure that the data come from expt using the new LFP board
    exptDate = out.fnames{a}{2}(2:7);
    month = str2num(exptDate(1:2));
    day = str2num(exptDate(3:4));
    year = str2num(exptDate(5:6));
    if (year<10) || ((year<=10) && (month<9)) || ((year<=10) && (month<=9) && (day<28))
        fprintf('Skipping <%s>\n', out.fnames{a}{2});
        continue
    end
    
    
    GT = gtobj(out.fnames{a}{2});
    nTrials = size(GT.trial,1);
    lfpRasIdx = GT.idx.lfp;
    anlgSumIdx = strcmpi(GT.sum.analog.sigid, GT.sum.rasterCells{lfpRasIdx});
    if ~any(anlgSumIdx); continue; end
    stimOn = mat2cell(GT.trial(:, GT.idx.stimon), ones(nTrials,1));
    nFrames = mat2cell(GT.trial(:, GT.idx.nFrames), ones(nTrials,1));
    frameRate = repmat({GT.sum.exptParams.framerate}, nTrials,1);
    stimOff = cellfun(@(x,y,z) (x/y)+z, nFrames, frameRate, stimOn, 'uniformoutput', 0);
    counts = cellfun(@(x,y,z) sum((x>y)&(x<=z)), GT.ras(:,GT.idx.spikes), stimOn, stimOff);
    l_p4 = GT.trial(:, GT.idx.protocol) == 4;
    ccs = GT.trial(:, (GT.idx.lcc | GT.idx.mcc | GT.idx.scc));
    l_s = ismember(sign(ccs), [0 0 1], 'rows');
    l_lvm = ismember(sign(ccs), [1 -1 0], 'rows');
    l_swl = ismember(sign(ccs), [1 -1 1], 'rows');
    l_swm = ismember(sign(ccs), [1 -1 -1], 'rows');
    l_lum = ismember(sign(ccs), [1 1 0], 'rows');
    
    %what's the csi?
    lvmCount = mean(counts(l_p4&l_lvm));
    sCount= mean(counts(l_p4&l_s));
    swlCount= mean(counts(l_p4&l_swl));
    swmCount= mean(counts(l_p4&l_swm));
    lumCount= mean(counts(l_p4&l_lum));
    if any(([lvmCount, sCount, swlCount, swmCount]./lumCount) < 2)
        max([lvmCount, sCount, swlCount, swmCount]./lumCount)
        continue
    end
    
    
    %calculate the spike response
    sSpike(a) = mean(counts(l_p4&l_s));
    lvmSpike(a) = mean(counts(l_p4&l_lvm));
    nStrials = sum(l_p4&l_s);
    nLVMtrials = sum(l_p4&l_lvm);
    
    band = {[8,15], [20 45], [50 70], [80 95]}; %for cheby
    %band = {10, 30, 60, 85}; %for greg's thing
    for bnd = 1:length(band)
        %determine the band limited power
        blp = blplfp(GT, lfpRasIdx, band{bnd});
        l_oob = cellfun(@(x) any(isnan(x)), blp);
        if sum(l_oob) == length(l_oob); close, continue; end
        
        % create and avereage across the gratings presentation
        anlgRate = repmat({GT.sum.analog.storeRates{anlgSumIdx}}, size(blp));
        stimOnSamp = cellfun(@(on, anSt, anRt)  round((on-anSt)*anRt), stimOn, GT.ras(:, GT.idx.anlgStart), anlgRate, 'uniformoutput', 0);
        nSampsPerStim = cellfun(@(x,y,z) round(x-y)*z, stimOff, stimOn, anlgRate, 'uniformoutput', 0);
        stimOffSamp = cellfun(@(onSamp, n) onSamp+n, stimOnSamp, nSampsPerStim, 'uniformoutput', 0);
        trlAvg = cellfun(@(x, b, e) nanmean(x(b:e)), blp, stimOnSamp, stimOffSamp);
        
        sblp = trlAvg(l_p4&l_s);
        lvmblp = trlAvg(l_p4&l_lvm);
        if numel(sblp) >=nStrials;
            sLFP{bnd}(a) = nanmean(sblp);
        end
        if numel(lvmblp) >= nLVMtrials;
            lvmLFP{bnd}(a) = nanmean(lvmblp);
        end
    end
end


figure
nbands = length(band);
for bnd = 1:nbands;
    subplot(1,nbands,bnd); hold on,
    plot(sSpike, sLFP{bnd}, 'bo')
    plot(lvmSpike, lvmLFP{bnd}, 'ro')
    c = [sLFP{bnd} ./ lvmLFP{bnd}];
    [h,p] = ttest(c(~isnan(c))-1)
    nanmean(c)
end
