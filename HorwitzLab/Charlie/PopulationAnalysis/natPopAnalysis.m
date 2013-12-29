% Noise At Thresh (NAT) compile batch data
%
%This script should take GT files and analyze noise at threshold
%between L-M and S-iso colors. I'll start by writing code for compiling
%batch data.

[fnames, spikenums] = fnamesFromTxt2('N:\NexFiles\Charlie\Kali\text files\noiseAtThreshPopList.txt');
for ex = 1:length(fnames)
    disp(fnames{ex}{1})
    GT = gtobj(fnames{ex}{1});
    out.GT{ex} = getGratingTuning(GT);
end

%save the compiled data
saveDir = 'N:\NexFiles\Charlie\Kali\batchData';
presDir = pwd;
str = ['nat_batch_' date];
cd(saveDir)
eval(['save ' str ' out']);
fprintf('**** Batch data saved to directory: %s ****\n', saveDir);
cd(presDir)

%% Load a pre-saved data structure
global currentNatBatch currentGTvsDTbatch
currentNatBatch = 'N:\NexFiles\Charlie\Kali\batchData\nat_batch_06-Jun-2010.mat';
currentGTvsDTbatch = 'N:\NexFiles\Charlie\Kali\batchData\fixvdet_batch_05-Jun-2010.mat';


%% Plot the SD as a function of contrast, and as a function of mean rate
clear
load(currentNatBatch)
figure;
for ex = 1:length(out.GT)
    sIdx = ismember(out.GT{ex}.crf.colors, [0 0 1], 'rows');
    sNorms = out.GT{ex}.crf.norms{sIdx};
    sMu = cellfun(@mean, out.GT{ex}.crf.rates{sIdx});
    sSigma = cellfun(@std, out.GT{ex}.crf.rates{sIdx});
    sVar = cellfun(@var, out.GT{ex}.crf.rates{sIdx});
    
    lvmIdx = ismember(out.GT{ex}.crf.colors, [1 -1 0], 'rows');
    lvmNorms = out.GT{ex}.crf.norms{lvmIdx};
    lvmMu = cellfun(@mean, out.GT{ex}.crf.rates{lvmIdx});
    lvmSigma = cellfun(@std, out.GT{ex}.crf.rates{lvmIdx});
    lvmVar = cellfun(@var, out.GT{ex}.crf.rates{lvmIdx});
    
    %sigma vs. contrast
    subplot(2,1,1)
    hold on,
    plot(sNorms, sSigma, 'b')
    plot(lvmNorms, lvmSigma, 'r')
    hold off
    xlabel('cone contrast')
    ylabel('sigma')
    
    %sigma vs mean rate
    subplot(2,1,2)
    hold on,
    plot(sMu, sSigma, 'b')
    plot(lvmMu, lvmSigma, 'r')
    hold off
    xlabel('mean rate')
    ylabel('sigma')
end


%% number of spikes "in v1" due to S-iso and L-M at detection threshold

%first, based on GT data
clear
global currentNatBatch %#ok<REDEF>
load(currentNatBatch)
load('C:\Documents and Settings\horwitzlab\Desktop\MatlabCode\Slave\REX slave programs\sCSF_Fits_Kali.mat')
sRates = nan(1,length(out.GT));
lvmRates = nan(1,length(out.GT));
for ex = 1:length(out.GT);
    colors = [0, 0, 1;1, -1, 0];
    for clr = 1:2;
        clrIdx = ismember(out.GT{ex}.crf.colors, colors(clr,:), 'rows');
        clrNorms = out.GT{ex}.crf.norms{clrIdx};
        clrMu = cellfun(@mean, out.GT{ex}.crf.rates{clrIdx});
        prefSF = out.GT{ex}.sf.prefSF;
        [~, idx] = ismember(sign(colors(clr,:)), sign(fp.colors), 'rows');
        thresh = 1./(polyval(fp.polyCoeff(idx,:), prefSF));
        if (thresh>min(clrNorms)) && (thresh < max(clrNorms))
            tmp = interp1(clrNorms, clrMu, thresh);
            if clr ==1
                sRates(ex) = tmp;
            elseif clr == 2;
                lvmRates(ex) = tmp;
            end
        else
            disp('out range')
        end
    end
end
figure
hold on,
plot([lvmRates; sRates], 'k-')
plot(ones(length(lvmRates)), lvmRates, 'r.', 'markersize', 5)
plot(1, nanmean(lvmRates), 'm*')
plot(ones(length(sRates)).*2, sRates, 'b.', 'markersize', 5)
plot(2, nanmean(sRates), 'm*')
xlim([0.5 2.5])
xlabel('Based on GT P7 data')
ylabel('Mean Rate')
invalid = isnan(sum([lvmRates(:) + sRates(:)],2));
valsForTtest = [lvmRates(:), sRates(:)];
valsForTtest(invalid,:) = [];
[h,p] = ttest(valsForTtest(:,1), valsForTtest(:,2));
title(sprintf('Paired Ttest, p: %.3f', p));

%now based on DT data.
clear
global currentGTvsDTbatch %#ok<REDEF>
filterOutFlatCRF = 0;
filterOutNoTrials = 1;
load(currentGTvsDTbatch)
sRates = nan(1,length(out.dat));
lvmRates = nan(1,length(out.dat));
for ex = 1:length(out.dat)
    colors = [0, 0, 1;1, -1, 0];
    for clr = 1:2;
        clrIdx = ismember(sign(out.dat{ex}.dt.expt.standColors), colors(clr,:), 'rows');
        clrNorms = out.dat{ex}.dt.expt.norms{clrIdx};
        clrMu = cellfun(@mean, out.dat{ex}.dt.c.crfIn{clrIdx});
        thresh = out.dat{ex}.dt.m.alpha(clrIdx);
        if (thresh>min(clrNorms)) && (thresh < max(clrNorms))
            tmp = interp1(clrNorms, clrMu, thresh);
            if clr ==1
                sRates(ex) = tmp;
            elseif clr == 2;
                lvmRates(ex) = tmp;
            end
        end
        if filterOutFlatCRF
            [Y, C] = deal([]);
            for cntrst = 1:length(clrNorms)
                Y = [Y, out.dat{ex}.dt.c.crfIn{clrIdx}{cntrst}];
                C = [C, repmat(clrNorms(cntrst), 1, length(out.dat{ex}.dt.c.crfIn{clrIdx}{cntrst}))];
            end
            [b, ~, stats] = glmfit(C, round(Y.*.666), 'poisson', 'link', 'log');
            stats.p;
            if stats.p(2) > 0.05
                if clr == 1
                    sRates(ex) = nan;
                elseif clr == 2;
                    lvmRates(ex) = nan;
                end
            end
        end
        if filterOutNoTrials
            nTrials = cellfun(@length, out.dat{ex}.dt.c.crfIn{clrIdx});
            if min(nTrials) < 8
                if clr == 1
                    sRates(ex) = nan;
                elseif clr == 2;
                    lvmRates(ex) = nan;
                end
            end
        end
    end
end
figure
hold on,
plot([lvmRates; sRates], 'k-')
plot(ones(length(lvmRates)), lvmRates, 'r.', 'markersize', 5)
plot(ones(length(sRates)).*2, sRates, 'b.', 'markersize', 5)
plot(1, nanmean(lvmRates), 'm*')
plot(2, nanmean(sRates), 'm*')
xlim([0.5 2.5])
xlabel('Based on DT data')
ylabel('Mean Rate')
[h,pt] = ttest(lvmRates(:), sRates(:));
pw = ranksum(lvmRates(~isnan(lvmRates)), sRates(~isnan(sRates)));
title(sprintf('Paired Ttest, p: %.3f \n Wilcoxon, p: %.3f', pt, pw));
set(gcf, 'name', sprintf('filtFlatCRF: %d, filtNoTrials: %d', filterOutFlatCRF, filterOutNoTrials))


%% S-iso vs L-M "spikes in v1" for protocol 4 GT trials

%the logic here is that the P4 stimuli are roughly equated for
%detectability and may show consistent results to the Mullen paper (i.e.,
%more spikes for S-iso than for L-M). 

clear
global currentNatBatch %#ok<REDEF>
load(currentNatBatch)
sRates = nan(1,length(out.GT));
lvmRates = nan(1,length(out.GT));
colors = [1, -1, 0; 0, 0, 1];
rates = nan(length(out.GT), 2);
for ex = 1:length(out.GT)
    for clr = 1:2
        clrIdx = ismember(sign(out.GT{ex}.color.colors), colors(clr,:), 'rows');
        rates(ex,clr) = out.GT{ex}.color.colresp(clrIdx, 1);
    end
end
figure
hold on,
plot(rates', 'k')
plot(ones(size(rates,1)), rates(:,1), 'r.')
plot(ones(size(rates,1))*2, rates(:,2), 'b.')
plot(1, mean(rates(:,1)), 'm*')
plot(2, mean(rates(:,2)), 'm*')
xlim([0.5, 2.5])
[h, p] = ttest(rates(:,1), rates(:,2));
title(sprintf('KALI\nTtest p: %.3f', p));
xlabel('data from P4 GT trials')


% now do the same analysis, but for Sedna.
clear
fnames = fnamesFromTxt2('N:\NexFiles\Charlie\Sedna\text files\allGrating.txt');
colors = [1, -1, 0; 0, 0, 1];
rates = nan(length(fnames), 2);
for ex = 1:length(fnames)
     out = getGratingTuning(gtobj(fnames{ex}));
     for clr = 1:2
         clrIdx = ismember(sign(out.color.colors), colors(clr,:), 'rows');
         rates(ex,clr) = out.color.colresp(clrIdx, 1);
     end
end
figure
hold on,
plot(rates', 'k')
plot(ones(size(rates,1)), rates(:,1), 'r.')
plot(ones(size(rates,1))*2, rates(:,2), 'b.')
plot(1, mean(rates(:,1)), 'm*')
plot(2, mean(rates(:,2)), 'm*')
xlim([0.5, 2.5])
[h, p] = ttest(rates(:,1), rates(:,2));
title(sprintf('SEDNA\nTtest p: %.3f', p));
xlabel('data from P4 GT trials')


%% L-M vs S TRs for DT
clear
global currentGTvsDTbatch
load(currentGTvsDTbatch)
[sTR, lvmTR, sNT, lvmNT] = deal(nan(length(out.dat),1));
colors = [0 0 1;1 -1 0];
for ex = 1:length(out.dat)
    for clr = 1:2;
        colorInd = ismember(sign(out.dat{ex}.dt.expt.standColors), colors(clr,:), 'rows');
        nt = out.dat{ex}.dt.c.alpha(colorInd);
        pt = out.dat{ex}.dt.m.alpha(colorInd);
        
        %filter out the expts where alpha wasn't in the range tested
        norms = out.dat{ex}.dt.expt.norms{colorInd};
        if ~((pt>min(norms)) && (pt < max(norms)))
            continue
        end
        
        %filter out expts where there are insufficient trials
        nTrials = cellfun(@length, out.dat{ex}.dt.c.crfIn{colorInd});
        if min(nTrials) < 8
            continue
        end
        
        %if you've made it here that means that the data are good.
        if clr == 1;
            sTR(ex) = nt/pt;
            sNT(ex) = nt;
        elseif clr == 2
            lvmTR(ex) = nt/pt;
            lvmNT(ex) = nt;
        end
    end
end

valForNan = 2;
tmp_sTR = sTR;
tmp_sTR(isnan(tmp_sTR)) = valForNan;
tmp_lvmTR = lvmTR;
tmp_lvmTR(isnan(tmp_lvmTR))= valForNan;
figure
hold on,
plot(tmp_lvmTR, tmp_sTR, 'k.')
plot([0.35 valForNan*1.02], [0.35 valForNan*1.02], 'k')
title('TR for L-M and S')
xlabel('L-M TR')
ylabel('S-iso TR')
set(gca, 'xscale', 'log', 'yscale', 'log')
xlim([0.35 2.1])
ylim([0.35 2.1])





