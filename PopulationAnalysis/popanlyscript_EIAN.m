%% IMPORT THE NECESSARY DATA INTO A datDB STRUCTURE

fin

mdb = initMouseDB('new'); % clear out the junk from previous datDB runs.

% buid a structure of params from physiology notes. The globals are used
% during execution of the physiology_notes script
global GL_ADD_TO_MDB GL_SUPPRESS_ANALYSIS
GL_ADD_TO_MDB = true; %#ok<NASGU>
GL_SUPPRESS_ANALYSIS = true; %#ok<NASGU>
physiology_notes % run to store prams in MDB.

% cd to the directory of workbooks and import the spreadsheet of cells
workbookpath = [GL_DOCUPATH, filesep, 'Other_workbooks', filesep, 'Cell_Library'];
workbooksheet = 'E_I and A_N ratios';
[~, ~, raw] = xlsread(workbookpath, workbooksheet);
header = raw(1,:);
mouseNames = raw(2:end,1);
siteNumber = raw(2:end,2);
goodNeurons(:,1) = cat(1,raw{2:end, 3});
goodNeurons(:,2) = cat(1,raw{2:end, 4});
goodNeurons = logical(goodNeurons);
neuronType = raw(2:end, 6:7);
HVA = raw(2:end, 5);
layer = [raw(2:end, 8); raw(2:end, 9)];


%
% build a population data structure by performing the appropriate analysis.
% I'll use the params.name field to build each new entry into the DB.
% Finish by saving the structure (possibly with a date string)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initalize things that I care about
[dat.ampa.peak_nS, dat.nmda.peak_nS, dat.excit.peak_nS, dat.inhib.peak_nS] = deal(nan(numel(mouseNames), 2));
[dat.ampa.peak_pA, dat.nmda.peak_pA, dat.excit.peak_pA, dat.inhib.peak_pA] = deal(nan(numel(mouseNames), 2));
[dat.NMDAR.ivcurve.mV, dat.NMDAR.ivcurve.mv_corrected, dat.NMDAR.ivcurve.pA] = deal(repmat({[] []}, numel(mouseNames), 1));
dat.cellDepth = deal(nan(numel(mouseNames), 2));


% initalize things for control analyses.
[dat.ampa.Verr, dat.nmda.Verr, dat.excit.Verr, dat.inhib.Verr] = deal(repmat({[] []}, numel(mouseNames), 1));
[dat.ampa.holding, dat.nmda.holding, dat.excit.holding, dat.inhib.holding] = deal(repmat({[] []}, numel(mouseNames), 1));
[dat.ampa.peakBySweep, dat.nmda.peakBySweep, dat.excit.peakBySweep, dat.inhib.peakBySweep] = deal(repmat({[] []}, numel(mouseNames), 1));
[dat.ampa.raw_pA, dat.nmda.raw_pA, dat.excit.raw_pA, dat.inhib.raw_pA] = deal(repmat({[] []}, numel(mouseNames), 1));

% iterate over the mice in the cell library (some mice get analyzed
% multiple times if there were multiple experiments per mouse)
for ex = 1:numel(mouseNames)
    
    ex_mouseName = mouseNames{ex};
    ex_siteNum = siteNumber{ex};
    mdb = initMouseDB(false, true);
    [~, idx] = mdb.search(ex_mouseName);
    
    params = mdb.mice{idx}.popAnly{ex_siteNum};
    params.fxns = {@anlyMod_optoIV, @anlyMod_EIbalance, @anlyMod_NMDAR};
    params = invitroAnalysisOverview(params);
    close all; drawnow
    
    % store the cell depth measurements
    dat.cellDepth(ex,:) = params.celldepth;
    
    %
    % add the data to an array for the AMPA/NMDA
    % Excitation/Inhibition ratio stuff. Adding 'NMDAR' to the 'group'
    % array allows the script to pull out NMDAR IV curve data
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    group = {'ampa', 'nmda', 'excit', 'inhib'};
    for g = 1:numel(group)
        
        switch group{g}
            case {'ampa', 'nmda', 'excit', 'inhib'}
                if isfield(params.isolatedData, group{g})
                    for ch = 1:2
                        
                        if ~isempty(params.isolatedData.(group{g}).raw_nS{ch});
                            dat.(group{g}).peak_nS(ex, ch) = params.isolatedData.(group{g}).peak_nS{ch};
                            dat.(group{g}).peak_pA(ex, ch) = params.isolatedData.(group{g}).peak_pA{ch};
                            dat.(group{g}).raw_pA{ex, ch} = params.isolatedData.(group{g}).raw_pA{ch};
                            dat.tvec{ex} = params.ivdat.tvec;
                            dat.(group{g}).Verr{ex, ch} = params.isolatedData.(group{g}).Verr{ch};
                            dat.(group{g}).peakBySweep{ex, ch} = params.isolatedData.(group{g}).peakBySweep_pA{ch};
                            dat.(group{g}).holding{ex,ch} = params.isolatedData.(group{g}).holdingCurrent{ch};
                        end
                    end
                end
                
            case 'NMDAR'
                if isfield(params.ivdat, 'NMDAR')
                    for ch = 1:2
                        
                        if ~isempty(params.ivdat.NMDAR.ivcurve.mV{ch});
                            dat.NMDAR.ivcurve.mV{ex, ch} = params.ivdat.NMDAR.ivcurve.mV{ch};
                            dat.NMDAR.ivcurve.mv_corrected{ex, ch} = params.ivdat.NMDAR.ivcurve.mV_corrected{ch};
                            dat.NMDAR.ivcurve.pA{ex, ch} = params.ivdat.NMDAR.ivcurve.pA{ch};
                        end
                        
                    end
                end
                
        end
        
        
    end % iterate over groups
    
end % iterate over recording sites

% package all the useful things into a single structure
dat.hva = HVA;
dat.goodNeurons = goodNeurons;
dat.neuronType = neuronType;
dat.mice = mouseNames;
dat.siteNum = cat(1, siteNumber{:});
dat.layer = layer;


% create grouping lists for HVAs
hvas = repmat(dat.hva, 2,1);
hvaList.('pm') = cellfun(@(x) ~isempty(x), regexpi(hvas, 'pm'));
hvaList.('lm') = cellfun(@(x) ~isempty(x), regexpi(hvas, 'lm'));
hvaList.('al') = cellfun(@(x) ~isempty(x), regexpi(hvas, 'al'));
hvaList.('und') = cellfun(@(x) ~isempty(x), regexpi(hvas, 'und'));

%
% create grouping lists for neuron type
neuronType = dat.neuronType(:);
neuronType = cellfun(@num2str, neuronType, 'uniformoutput', false);
l_IN = cellfun(@(x) ~isempty(x), regexpi(neuronType, 'in'));
l_SOM = cellfun(@(x) ~isempty(x), regexpi(neuronType, 'som'));
typeList.IN = l_IN | l_SOM;
typeList.SOM = l_SOM;
typeList.PY = cellfun(@(x) ~isempty(x), regexpi(neuronType, 'py'));
typeList.und = ~typeList.IN & ~typeList.PY;


%
% create a grouping list for layer
layer = dat.layer(:);
layer = cellfun(@num2str, layer, 'uniformoutput', false);
layerList.L_23 = cellfun(@(x) ~isempty(x), regexpi(layer, '2/3'));
layerList.L_4 = cellfun(@(x) ~isempty(x), regexpi(layer, '4'));
layerList.L_5 = cellfun(@(x) ~isempty(x), regexpi(layer, '5'));


% save the data and the grouping lists
originalDir = pwd;
cd(GL_POPDATPATH);
save('popAnly_EIAN.mat', 'dat', 'hvaList', 'typeList', 'layerList')
cd(originalDir);

% be nice and return these variables to their default values
GL_ADD_TO_MDB = false;
GL_SUPPRESS_ANALYSIS = false;

%% EXCITATION VS. INHIBITION

fin

% load in the pre-saved population data
load([GL_POPDATPATH, 'popAnly_EIAN.mat']);
l_23 = layerList.L_23;
l_valid = dat.goodNeurons(:) & l_23;


%
% plot the raw currents for isolated excitation and inhibition
%
%%%%%%%%%%%%%%%%%%%%%%%%%
nPlots = sum(l_valid);
nRows = ceil(sqrt(nPlots));
tmp_inhib_raw = cat(1, dat.inhib.raw_pA(:));
tmp_excit_raw = cat(1, dat.excit.raw_pA(:));
tmp_names = repmat(dat.mice, 2,1);
tmp_tvec = repmat(dat.tvec', 2,1);
listToPlot = find(l_valid);
figure
set(gcf, 'position', [7 18 1396 795])
for i = 1:nPlots;
    
    a = listToPlot(i);
    if any([isempty(tmp_inhib_raw{a}), isempty(tmp_excit_raw{a})])
        continue
    end
    
    subplot(nRows, nRows, i)
    hold on,
    plot(tmp_tvec{a}, tmp_inhib_raw{a}, 'r')
    plot(tmp_tvec{a}, tmp_excit_raw{a}, 'b')
    t = title(tmp_names{a});
    set(t, 'interpreter', 'none')
    axis tight
    hold off
    
end

% pull out peak conductances. do some error checking, and make sure
% all the nS values are positive.
excit_nS_signed = dat.excit.peak_nS(:);
excit_nS_signed = excit_nS_signed(l_valid);
assert(all(excit_nS_signed<0), 'ERROR: some excitatory currents are positive...')
excit_nS_unsigned = abs(excit_nS_signed);

inhib_nS_signed = dat.inhib.peak_nS(:);
inhib_nS_signed = inhib_nS_signed(l_valid);
assert(all(inhib_nS_signed>0), 'ERROR: some inhibitory currents are negative...')
inhib_nS_unsigned = abs(inhib_nS_signed);



%
% plot peak conductances for all cells
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp_mice = [dat.mice(:); dat.mice(:)];
tmp_mice = tmp_mice(l_valid);
tmp_siteNum = [dat.siteNum(:); dat.siteNum(:)];
tmp_siteNum = tmp_siteNum(l_valid);
figure, hold on,
set(gca, 'fontsize', 25)
for a = 1:numel(excit_nS_unsigned)
    p = plot(excit_nS_unsigned(a), inhib_nS_unsigned(a), 'ko', 'markerfacecolor', 'k', 'markersize', 10);
    printTitle = @(a,b,c) title(sprintf('%s, cell %d',c{1},c{2}));
    set(p, 'buttonDownFcn', {printTitle, {tmp_mice{a}, tmp_siteNum(a)}})
    t = get(get(p, 'parent'), 'title');
end
l_nan = isnan(excit_nS_unsigned) | isnan(inhib_nS_unsigned);
betas = [excit_nS_unsigned(~l_nan), ones(size(excit_nS_unsigned(~l_nan)))] \ inhib_nS_unsigned(~l_nan);
xvals_all = get(gca, 'xlim');
yvals_all = get(gca, 'ylim');
plot(xvals_all, [xvals_all(:), ones(2,1)]*betas(:), '--k', 'linewidth', 3)
xlabel('Excit conductance (nS)')
ylabel('Inhib conductance (nS)')
title(sprintf('ALL DATA: E/I ratio = %.2f', betas(1)))
set(t, 'interpreter', 'none');

%
% plot the data, but color code each HVA.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure, hold on,
title('By HVAs, only PY cells')
set(gca, 'fontsize', 25)
groups = {'und', 'pm', 'lm', 'al'};
h_fit = [];
for a = 1:numel(groups);
    l_PY = typeList.PY;
    l_HVA = hvaList.(groups{a});
    l_toPlot = l_PY & l_HVA;
    l_toPlot(~l_valid) = [];
    tmp_excit = excit_nS_unsigned(l_toPlot);
    tmp_inhib = inhib_nS_unsigned(l_toPlot);
    [clr_fit, clr_raw] = hvaPlotColor(groups{a});
    plot(tmp_excit, tmp_inhib, 'o', 'markerfacecolor', clr_raw, 'markeredgecolor', clr_raw, 'markersize', 10)
    l_nan = isnan(tmp_excit) | isnan(tmp_inhib);
    betas = [tmp_excit(~l_nan), ones(size(tmp_excit(~l_nan)))] \ tmp_inhib(~l_nan);
    xvals = get(gca, 'xlim');
    h = plot(xvals, [xvals(:), ones(2,1)]*betas(:), '--', 'color', clr_fit, 'linewidth', 3);
    h_fit = [h_fit,h];
    xlim(xvals_all)
    ylim(yvals_all)
    xlabel('Excit conductance (nS)')
    ylabel('Inhib conductance (nS)')
    
end
legend(h_fit, groups, 'location', 'southeast')



%
% Plot all E/I ratios, color code by neuron type
%
%%%%%%%%%%%%%%%%%%
groupTypes = {'IN', 'PY', 'SOM'};
figure, hold on,
title('By cell type')
set(gca, 'fontsize', 25)
h_fit = [];
for a = 1:numel(groupTypes);
    l_toPlot = typeList.(groupTypes{a});
    l_toPlot(~l_valid) = [];
    tmp_excit = excit_nS_unsigned(l_toPlot);
    tmp_inhib = inhib_nS_unsigned(l_toPlot);
    [clr_fit, clr_raw] = hvaPlotColor(groupTypes{a});
    plot(tmp_excit, tmp_inhib, 'o', 'markerfacecolor', clr_raw, 'markeredgecolor', clr_raw, 'markersize', 10)
    l_nan = isnan(tmp_excit) | isnan(tmp_inhib);
    betas = [tmp_excit(~l_nan), ones(size(tmp_excit(~l_nan)))] \ tmp_inhib(~l_nan);
    xvals = get(gca, 'xlim');
    h = plot(xvals, [xvals(:), ones(2,1)]*betas(:), '--', 'color', clr_fit, 'linewidth', 3);
    h_fit = [h_fit, h];
    xlim(xvals_all)
    ylim(yvals_all)
    xlabel('Excit conductance (nS)')
    ylabel('Inhib conductance (nS)')
end
legend(h_fit, groupTypes, 'location', 'southeast')



%
% A general figure that has one subplot for each cell type and shows data
% across HVAs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nHVAs = numel(fieldnames(hvaList));
nCellTypes = numel(fieldnames(typeList));
HVATypes = fieldnames(hvaList);
cellTypes = fieldnames(typeList);
figure
for iArea = 1:nHVAs
    
    % initialize the plot
    subplot(1,nHVAs, iArea), hold on,
    title(sprintf('HVA: %s', HVATypes{iArea}));
    
    h_fit = [];
    l_leg = true(nCellTypes,1);
    for iType = 1:nCellTypes
        % make a list of the appropriate cells
        l_type = typeList.(cellTypes{iType});
        l_hva = hvaList.(HVATypes{iArea});
        l_toPlot = l_type & l_hva;
        l_toPlot(~l_valid) = [];
        if ~any(l_toPlot);
            l_leg(iType) = false;
            continue
        end
        
        tmp_excit = excit_nS_unsigned(l_toPlot);
        tmp_inhib = inhib_nS_unsigned(l_toPlot);
        [clr_fit, clr_raw] = hvaPlotColor(cellTypes{iType});
        plot(tmp_excit, tmp_inhib, 'o', 'markerfacecolor', clr_raw, 'markeredgecolor', clr_raw, 'markersize', 10)
        l_nan = isnan(tmp_excit) | isnan(tmp_inhib);
        betas = [tmp_excit(~l_nan), ones(size(tmp_excit(~l_nan)))] \ tmp_inhib(~l_nan);
        xvals = get(gca, 'xlim');
        h = plot(xvals, [xvals(:), ones(2,1)]*betas(:), '--', 'color', clr_fit, 'linewidth', 3);
        h_fit = [h_fit, h];
        xlim(xvals_all)
        ylim(yvals_all)
    end
    legend(h_fit, cellTypes{l_leg}, 'location', 'northwest')
end




%% E/I RATIOS BY CELL DEPTH

fin

% load in the pre-saved population data
load([GL_POPDATPATH, 'popAnly_EIAN.mat']);
l_valid = dat.goodNeurons(:) & typeList.PY & layerList.L_23;


% pull out peak conductances. do some error checking, and make sure
% all the nS values are positive.
excit_nS_signed = dat.excit.peak_nS(:);
excit_nS_signed = excit_nS_signed(l_valid);
assert(all(excit_nS_signed<0), 'ERROR: some excitatory currents are positive...')
excit_nS_unsigned = abs(excit_nS_signed);

inhib_nS_signed = dat.inhib.peak_nS(:);
inhib_nS_signed = inhib_nS_signed(l_valid);
assert(all(inhib_nS_signed>0), 'ERROR: some inhibitory currents are negative...')
inhib_nS_unsigned = abs(inhib_nS_signed);

% define the E/I ratio
ei_ratio = excit_nS_unsigned ./ inhib_nS_unsigned;

% define the cell depth measurements
cellDepth = dat.cellDepth(:);
cellDepth = cellDepth(l_valid);




%
% plot all the data and make points clickable
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp_mice = [dat.mice(:); dat.mice(:)];
tmp_mice = tmp_mice(l_valid);
tmp_siteNum = [dat.siteNum(:); dat.siteNum(:)];
tmp_siteNum = tmp_siteNum(l_valid);
figure, hold on,
set(gca, 'fontsize', 25)
for a = 1:numel(ei_ratio)
    p = plot(cellDepth(a), ei_ratio(a), 'ko', 'markerfacecolor', 'k', 'markersize', 10);
    printTitle = @(a,b,c) title(sprintf('%s, cell %d',c{1},c{2}));
    set(p, 'buttonDownFcn', {printTitle, {tmp_mice{a}, tmp_siteNum(a)}})
    t = get(get(p, 'parent'), 'title');
end
xlabel('cell depth')
ylabel('ei_ratio')
set(t, 'interpreter', 'none');
set(gca, 'yscale', 'log')



%
% color code points by HVA or by laminar location
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% grouping list for layer
l_L23 = layerList.L_23;
l_L23 = l_L23(l_valid);
l_L4 = layerList.L_4;
l_L4 = l_L4(l_valid);
l_L5 = layerList.L_5;
l_L5 = l_L5(l_valid);


figure, hold on,
plot(cellDepth(l_L23), ei_ratio(l_L23), 'bo', 'markerfacecolor', 'b')
plot(cellDepth(l_L4), ei_ratio(l_L4), 'ko', 'markerfacecolor', 'k')
plot(cellDepth(l_L5), ei_ratio(l_L5), 'ro', 'markerfacecolor', 'r')
set(gca, 'yscale', 'log')
legend('Layer 2/3', 'Layer 4', 'Layer 5')
ylabel('Ex / Inhib ratio')
xlabel('Cell Depth')
title('Grouped by Laminar Location')


% now by area
l_AL = hvaList.al;
l_AL = l_AL(l_valid);
l_PM = hvaList.pm;
l_PM = l_PM(l_valid);

figure, hold on,
plot(cellDepth(l_AL), ei_ratio(l_AL), 'ro', 'markerfacecolor', 'r')
plot(cellDepth(l_PM), ei_ratio(l_PM), 'bo', 'markerfacecolor', 'b')
set(gca, 'yscale', 'log')
legend('Area AL', 'Area PM')
ylabel('Ex / Inhib ratio')
xlabel('Cell Depth')
title('Grouped by brain area')


%% AMPA to NMDA RATIOS

fin

% load in the pre-saved population data
load([GL_POPDATPATH, 'popAnly_EIAN.mat'])
l_valid = dat.goodNeurons(:);

% pull out peak conductances. do some error checking, and then convert the
% nS values into positive numbers.
ampa_nS_signed = dat.ampa.peak_nS(:);
nmda_nS_signed = dat.nmda.peak_nS(:);

l_valid = l_valid & ~isnan(ampa_nS_signed) & ~isnan(nmda_nS_signed);
ampa_nS_signed = ampa_nS_signed(l_valid);
nmda_nS_signed = nmda_nS_signed(l_valid);

assert(all(ampa_nS_signed<0), 'ERROR: some AMPA currents are positive or nan')
assert(all(nmda_nS_signed>0), 'ERROR: some NMDA currents are negative or nan')
ampa_nS_unsigned = abs(ampa_nS_signed);
nmda_nS_unsigned = abs(nmda_nS_signed);




%
% plot peak conductances for all cells
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure, hold on,
set(gca, 'fontsize', 25)
plot(ampa_nS_unsigned, nmda_nS_unsigned, 'ko', 'markerfacecolor', 'k', 'markersize', 10)
l_nan = isnan(ampa_nS_unsigned) | isnan(nmda_nS_unsigned);
betas = [ampa_nS_unsigned(~l_nan), ones(size(ampa_nS_unsigned(~l_nan)))] \ nmda_nS_unsigned(~l_nan);
xvals_all = get(gca, 'xlim');
yvals_all = get(gca, 'ylim');
plot(xvals_all, [xvals_all(:), ones(2,1)]*betas(:), '--k', 'linewidth', 3)
xlabel('AMPA conductance (nS)')
ylabel('NMDA conductance (nS)')
title(sprintf('ALL DATA: A/N ratio = %.2f', betas(1)))



%
% plot the data, but color code each HVA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure, hold on,
set(gca, 'fontsize', 25)
groups = {'pm', 'lm', 'al', 'und'};
title('A/N ratio by HVA, PY cells only')
h_fit =[];
for a = 1:numel(groups);
    l_PY = typeList.PY;
    l_HVA = hvaList.(groups{a});
    l_toPlot = l_PY & l_HVA;
    l_toPlot(~l_valid) = [];
    tmp_ampa = ampa_nS_unsigned(l_toPlot);
    tmp_nmda = nmda_nS_unsigned(l_toPlot);
    [clr_fit, clr_raw] = hvaPlotColor(groups{a});
    plot(tmp_ampa, tmp_nmda, 'o', 'markerfacecolor', clr_raw, 'markeredgecolor', clr_raw, 'markersize', 10)
    l_nan = isnan(tmp_ampa) | isnan(tmp_nmda);
    betas = [tmp_ampa(~l_nan), ones(size(tmp_ampa(~l_nan)))] \ tmp_nmda(~l_nan);
    xvals = get(gca, 'xlim');
    h = plot(xvals, [xvals(:), ones(2,1)]*betas(:), '--', 'color', clr_fit, 'linewidth', 3);
    h_fit = [h_fit, h];
    xlim(xvals_all)
    ylim(yvals_all)
    xlabel('AMPA conductance (nS)')
    ylabel('NMDA conductance (nS)')
end
legend(h_fit, groups, 'location', 'northwest')


groupTypes = {'IN', 'PY', 'SOM'};
figure, hold on,
set(gca, 'fontsize', 25)
title('By cell type')
h_fit = [];
for a = 1:numel(groupTypes);
    l_toPlot = typeList.(groupTypes{a});
    l_toPlot(~l_valid) = [];
    tmp_ampa = ampa_nS_unsigned(l_toPlot);
    tmp_nmda = nmda_nS_unsigned(l_toPlot);
    [clr_fit, clr_raw] = hvaPlotColor(groupTypes{a});
    plot(tmp_ampa, tmp_nmda, 'o', 'markerfacecolor', clr_raw, 'markeredgecolor', clr_raw, 'markersize', 10)
    l_nan = isnan(tmp_ampa) | isnan(tmp_nmda);
    betas = [tmp_ampa(~l_nan), ones(size(tmp_ampa(~l_nan)))] \ tmp_nmda(~l_nan);
    xvals = get(gca, 'xlim');
    h = plot(xvals, [xvals(:), ones(2,1)]*betas(:), '--', 'color', clr_fit, 'linewidth', 3);
    h_fit = [h_fit,h];
    xlim(xvals_all)
    ylim(yvals_all)
    xlabel('AMPA conductance (nS)')
    ylabel('NMDA conductance (nS)')
end
legend(h_fit, groupTypes, 'location', 'northwest')


%% A/N RATIOS BY CELL DEPTH

fin

% load in the pre-saved population data
load([GL_POPDATPATH, 'popAnly_EIAN.mat']);
l_valid = dat.goodNeurons(:) & typeList.IN;


% pull out peak conductances. do some error checking, and make sure
% all the nS values are positive.
ampa_nS_signed = dat.ampa.peak_nS(:);
ampa_nS_signed = ampa_nS_signed(l_valid);
assert(all(ampa_nS_signed<0 | isnan(ampa_nS_signed)), 'ERROR: some AMPA currents are positive...')
ampa_nS_unsigned = abs(ampa_nS_signed);

nmda_nS_signed = dat.nmda.peak_nS(:);
nmda_nS_signed = nmda_nS_signed(l_valid);
assert(all(nmda_nS_signed>0 | isnan(nmda_nS_signed)), 'ERROR: some NMDA currents are negative...')
nmda_nS_unsigned = abs(nmda_nS_signed);

% define the E/I ratio
an_ratio = ampa_nS_unsigned ./ nmda_nS_unsigned;

% define the cell depth measurements
cellDepth = dat.cellDepth(:);
cellDepth = cellDepth(l_valid);




%
% plot all the data and make points clickable
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp_mice = [dat.mice(:); dat.mice(:)];
tmp_mice = tmp_mice(l_valid);
tmp_siteNum = [dat.siteNum(:); dat.siteNum(:)];
tmp_siteNum = tmp_siteNum(l_valid);
figure, hold on,
set(gca, 'fontsize', 25)
for a = 1:numel(an_ratio)
    p = plot(cellDepth(a), an_ratio(a), 'ko', 'markerfacecolor', 'k', 'markersize', 10);
    printTitle = @(a,b,c) title(sprintf('%s, cell %d',c{1},c{2}));
    set(p, 'buttonDownFcn', {printTitle, {tmp_mice{a}, tmp_siteNum(a)}})
    t = get(get(p, 'parent'), 'title');
end
xlabel('cell depth')
ylabel('A/N ratio')
set(t, 'interpreter', 'none');
set(gca, 'yscale', 'log')



%
% color code points by HVA or by laminar location
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% grouping list for layer
l_L23 = layerList.L_23;
l_L23 = l_L23(l_valid);
l_L4 = layerList.L_4;
l_L4 = l_L4(l_valid);
l_L5 = layerList.L_5;
l_L5 = l_L5(l_valid);


figure, hold on,
plot(cellDepth(l_L23), an_ratio(l_L23), 'bo', 'markerfacecolor', 'b')
plot(cellDepth(l_L4), an_ratio(l_L4), 'ko', 'markerfacecolor', 'k')
plot(cellDepth(l_L5), an_ratio(l_L5), 'ro', 'markerfacecolor', 'r')
set(gca, 'yscale', 'log')
legend('Layer 2/3', 'Layer 4', 'Layer 5')
ylabel('AMPA / NMDA ratio')
xlabel('Cell Depth')
title('Grouped by Laminar Location')


% now by area
group = {'al', 'pm', 'lm', 'und'};
figure, hold on,
for a = 1:numel(group)
    l_hva = hvaList.(group{a});
    l_hva(~l_valid) = [];
    [clr_fit, clr_raw] = hvaPlotColor(group{a});
    plot(cellDepth(l_hva), an_ratio(l_hva), 'o', 'markeredgecolor', clr_raw, 'markerfacecolor', clr_raw, 'markersize', 8)
end
set(gca, 'yscale', 'log')
legend(group)
ylabel('AMPA / NMDA ratio')
xlabel('Cell Depth')
title('Grouped by brain area')

%% A/N vs E/I ratios

fin

% load in the pre-saved population data
load([GL_POPDATPATH, 'popAnly_EIAN.mat']);
l_valid = dat.goodNeurons(:) & typeList.PY & layerList.L_23;


% pull out peak conductances. do some error checking, and make sure
% all the nS values are positive.
ampa_nS_signed = dat.ampa.peak_nS(:);
ampa_nS_signed = ampa_nS_signed(l_valid);
assert(all(ampa_nS_signed<0 | isnan(ampa_nS_signed)), 'ERROR: some AMPA currents are positive...')
ampa_nS_unsigned = abs(ampa_nS_signed);

nmda_nS_signed = dat.nmda.peak_nS(:);
nmda_nS_signed = nmda_nS_signed(l_valid);
assert(all(nmda_nS_signed>0 | isnan(nmda_nS_signed)), 'ERROR: some NMDA currents are negative...')
nmda_nS_unsigned = abs(nmda_nS_signed);

% define the A/N ratio
an_ratio = ampa_nS_unsigned ./ nmda_nS_unsigned;



% pull out peak conductances. do some error checking, and make sure
% all the nS values are positive.
excit_nS_signed = dat.excit.peak_nS(:);
excit_nS_signed = excit_nS_signed(l_valid);
assert(all(excit_nS_signed<0), 'ERROR: some excitatory currents are positive...')
excit_nS_unsigned = abs(excit_nS_signed);

inhib_nS_signed = dat.inhib.peak_nS(:);
inhib_nS_signed = inhib_nS_signed(l_valid);
assert(all(inhib_nS_signed>0), 'ERROR: some inhibitory currents are negative...')
inhib_nS_unsigned = abs(inhib_nS_signed);

% define the E/I ratio
ei_ratio = excit_nS_unsigned ./ inhib_nS_unsigned;


%
% plot all the data and make points clickable
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp_mice = [dat.mice(:); dat.mice(:)];
tmp_mice = tmp_mice(l_valid);
tmp_siteNum = [dat.siteNum(:); dat.siteNum(:)];
tmp_siteNum = tmp_siteNum(l_valid);
figure, hold on,
set(gca, 'fontsize', 25)
for a = 1:numel(an_ratio)
    p = plot(ei_ratio(a), an_ratio(a), 'ko', 'markerfacecolor', 'k', 'markersize', 10);
    printTitle = @(a,b,c) title(sprintf('%s, cell %d',c{1},c{2}));
    set(p, 'buttonDownFcn', {printTitle, {tmp_mice{a}, tmp_siteNum(a)}})
    t = get(get(p, 'parent'), 'title');
end
x = xlabel('E/I ratio');
y = ylabel('A/N ratio');
set(gca, 'yscale', 'log', 'xscale', 'log')
set([t,x,y], 'interpreter', 'none');

%
% plot by HVA
%
groups = {'al', 'pm', 'lm', 'und'};
figure, hold on,
for a = 1:numel(groups)
    l_hva = hvaList.(groups{a});
    l_hva(~l_valid) = [];
    [clr_fit, clr_raw] = hvaPlotColor(groups{a});
    
    plot(ei_ratio(l_hva), an_ratio(l_hva), 'o', 'markeredgecolor', clr_raw, 'markerfacecolor', clr_raw, 'markersize', 10);
end
x = xlabel('E/I ratio');
y = ylabel('A/N ratio');
set(gca, 'yscale', 'log', 'xscale', 'log')




%% CONTROL ANALYSES

fin

% load in the pre-saved population data
load([GL_POPDATPATH, 'popAnly_EIAN.mat'])


% create grouping lists
l_valid = dat.goodNeurons(:);

% grab the data (vclamp errors)
tmp_ampa_verr = dat.ampa.Verr(:);
tmp_ampa_verr = tmp_ampa_verr(l_valid);
tmp_nmda_verr = dat.nmda.Verr(:);
tmp_nmda_verr = tmp_nmda_verr(l_valid);
tmp_inhib_verr = dat.inhib.Verr(:);
tmp_inhib_verr = tmp_inhib_verr(l_valid);
Nplts = ceil(sqrt(numel(tmp_nmda_verr)));

figure
for a = 1:numel(tmp_ampa_verr)
    
    if isempty(tmp_ampa_verr{a}); continue; end
    
    subplot(Nplts, Nplts, a)
    hold on,
    plot(tmp_ampa_verr{a}, 'k', 'linewidth', 2)
    plot(tmp_nmda_verr{a}, 'b', 'linewidth', 2)
    plot(tmp_inhib_verr{a}, 'r', 'linewidth', 2)
    set(gcf, 'name', 'Voltage Clamp Errors')
end





% grab the data (holding current)
tmp_ampa_hold = dat.ampa.holding(:);
tmp_ampa_hold = tmp_ampa_hold(l_valid);
tmp_nmda_hold = dat.nmda.holding(:);
tmp_nmda_hold = tmp_nmda_hold(l_valid);
tmp_inhib_hold = dat.inhib.holding(:);
tmp_inhib_hold = tmp_inhib_hold(l_valid);
Nplts = ceil(sqrt(numel(tmp_nmda_verr)));

figure
for a = 1:numel(tmp_ampa_hold)
    
    if isempty(tmp_ampa_hold{a}); continue; end
    
    subplot(Nplts, Nplts, a)
    hold on,
    plot((tmp_ampa_hold{a}-tmp_ampa_hold{a}(1))./tmp_ampa_hold{a}(1), 'k', 'linewidth', 2)
    plot((tmp_nmda_hold{a}-tmp_nmda_hold{a}(1))./tmp_nmda_hold{a}(1), 'b', 'linewidth', 2)
    plot((tmp_inhib_hold{a}-tmp_inhib_hold{a}(1))./tmp_inhib_hold{a}(1), 'r', 'linewidth', 2)
    set(gcf, 'name', 'Holding Current')
end




% grab the data (peak current by sweep)
tmp_ampa_peak = dat.ampa.peakBySweep(:);
tmp_ampa_peak = tmp_ampa_peak(l_valid);
tmp_nmda_peak = dat.nmda.peakBySweep(:);
tmp_nmda_peak = tmp_nmda_peak(l_valid);
tmp_inhib_peak = dat.inhib.peakBySweep(:);
tmp_inhib_peak = tmp_inhib_peak(l_valid);
Nplts = ceil(sqrt(numel(tmp_nmda_verr)));

figure
for a = 1:numel(tmp_ampa_peak)
    
    if isempty(tmp_ampa_peak{a}); continue; end
    
    subplot(Nplts, Nplts, a)
    hold on,
    plot((tmp_ampa_peak{a}-tmp_ampa_peak{a}(1))./tmp_ampa_peak{a}(1), 'k', 'linewidth', 2)
    plot((tmp_nmda_peak{a}-tmp_nmda_peak{a}(1))./tmp_nmda_peak{a}(1), 'b', 'linewidth', 2)
    plot((tmp_inhib_peak{a}-tmp_inhib_peak{a}(1))./tmp_inhib_peak{a}(1), 'r', 'linewidth', 2)
    set(gcf, 'name', 'Peak Currents')
end














