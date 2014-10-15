%% IMPORT THE NECESSARY DATA INTO A datDB STRUCTURE

fin

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


%
% build a population data structure by performing the appropriate analysis.
% I'll use the params.name field to build each new entry into the DB.
% Finish by saving the structure (possibly with a date string)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% iterate over the mice in the cell library (some mice get analyzed
% multiple times if there were multiple experiments per mouse)
[dat.ampa.peak_nS, dat.nmda.peak_nS, dat.excit.peak_nS, dat.inhib.peak_nS] = deal(nan(numel(mouseNames), 2));
[dat.ampa.peak_pA, dat.nmda.peak_pA, dat.excit.peak_pA, dat.inhib.peak_pA] = deal(nan(numel(mouseNames), 2));
[dat.ampa.Verr, dat.nmda.Verr, dat.excit.Verr, dat.inhib.Verr] = deal(repmat({[] []}, numel(mouseNames), 1));
[dat.ampa.holding, dat.nmda.holding, dat.excit.holding, dat.inhib.holding] = deal(repmat({[] []}, numel(mouseNames), 1));
[dat.ampa.peakBySweep, dat.nmda.peakBySweep, dat.excit.peakBySweep, dat.inhib.peakBySweep] = deal(repmat({[] []}, numel(mouseNames), 1));
[dat.ampa.raw_pA, dat.nmda.raw_pA, dat.excit.raw_pA, dat.inhib.raw_pA] = deal(repmat({[] []}, numel(mouseNames), 1));
for ex = 1:numel(mouseNames)
    
    ex_mouseName = mouseNames{ex};
    ex_siteNum = siteNumber{ex};
    mdb = initMouseDB(false, true);
    [~, idx] = mdb.search(ex_mouseName);
    
    params = mdb.mice{idx}.popAnly{ex_siteNum};
    params.fxns = {@anlyMod_optoIV, @anlyMod_EIbalance};
    params = invitroAnalysisOverview(params);
    close all; drawnow
    
    %
    % add the data to an array for the AMPA/NMDA
    % Excitation/Inhibition ratio stuff
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    idx = goodNeurons(ex,:);
    group = {'ampa', 'nmda', 'excit', 'inhib'};
    for g = 1:numel(group)
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
    end
    
end

% package all the useful things into a single structure (and then save the
% structure)
dat.hva = HVA;
dat.goodNeurons = goodNeurons;
dat.neuronType = neuronType;
dat.mice = mouseNames;
dat.siteNum = cat(1, siteNumber{:});

originalDir = pwd;
cd(GL_POPDATPATH);
save popAnly_EIAN.mat dat
cd(originalDir);

% be nice and return these variables to their default values
GL_ADD_TO_MDB = false;
GL_SUPPRESS_ANALYSIS = false;

%% EXCITATION VS. INHIBITION

fin

% load in the pre-saved population data
load([GL_POPDATPATH, 'popAnly_EIAN.mat'])


% create grouping lists
l_valid = dat.goodNeurons(:);
hvas = repmat(dat.hva, 2,1);
hvas = hvas(l_valid);
hvalist.('pm') = cellfun(@(x) ~isempty(x), regexpi(hvas, 'pm'));
hvalist.('lm') = cellfun(@(x) ~isempty(x), regexpi(hvas, 'lm'));
hvalist.('al') = cellfun(@(x) ~isempty(x), regexpi(hvas, 'al'));
hvalist.('und') = cellfun(@(x) ~isempty(x), regexpi(hvas, 'und'));

% now group by IN vs PY
neuronType = dat.neuronType(:);
neuronType = cellfun(@num2str, neuronType, 'uniformoutput', false);
neuronType = neuronType(l_valid);
typeList.IN = cellfun(@(x) ~isempty(x), regexpi(neuronType, 'in'));
typeList.PY = cellfun(@(x) ~isempty(x), regexpi(neuronType, 'py'));
typeList.und = ~typeList.IN & ~typeList.PY;

% plot the raw currents for isolated excitation and inhibition
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

% pull out peak conductances
raw_excit = dat.excit.peak_pA(:);
raw_excit = abs(raw_excit(l_valid));
raw_inhib = dat.inhib.peak_pA(:);
raw_inhib = raw_inhib(l_valid);

% plot peak conductances for all cells
figure, hold on,
set(gca, 'fontsize', 25)
plot(raw_excit, raw_inhib, 'ko', 'markerfacecolor', 'k', 'markersize', 10)
l_nan = isnan(raw_excit);
betas = [raw_excit(~l_nan), ones(size(raw_excit(~l_nan)))] \ raw_inhib(~l_nan);
xvals_all = get(gca, 'xlim');
yvals_all = get(gca, 'ylim');
plot(xvals_all, [xvals_all(:), ones(2,1)]*betas(:), '--k', 'linewidth', 3)
xlabel('Excit conductance (nS)')
ylabel('Inhib conductance (nS)')
title(sprintf('ALL DATA: E/I ratio = %.2f', betas(1)))

% plot the data, but color code each HVA
figure, hold on,
title('By HVAs')
set(gca, 'fontsize', 25)
groups = {'und', 'pm', 'lm', 'al'};
for a = 1:numel(groups);
    tmp_excit = raw_excit(hvalist.(groups{a}));
    tmp_inhib = raw_inhib(hvalist.(groups{a}));
    [clr_fit, clr_raw] = hvaPlotColor(groups{a});
    plot(tmp_excit, tmp_inhib, 'o', 'markerfacecolor', clr_raw, 'markeredgecolor', clr_raw, 'markersize', 10)
    l_nan = isnan(tmp_excit) | isnan(tmp_inhib);
    betas = [tmp_excit(~l_nan), ones(size(tmp_excit(~l_nan)))] \ tmp_inhib(~l_nan);
    xvals = get(gca, 'xlim');
    plot(xvals, [xvals(:), ones(2,1)]*betas(:), '--', 'color', clr_raw, 'linewidth', 3)
    xlim(xvals_all)
    ylim(yvals_all)
    xlabel('Excit conductance (nS)')
    ylabel('Inhib conductance (nS)')
end


% color code by neuron type
groupTypes = {'IN', 'PY', 'und'};
pltClr = {'r', 'b', [.85 .85 .85]};
figure, hold on,
title('By cell type')
set(gca, 'fontsize', 25)
for a = 1:numel(groupTypes);
    tmp_excit = raw_excit(typeList.(groupTypes{a}));
    tmp_inhib = raw_inhib(typeList.(groupTypes{a}));
    clr_raw = pltClr{a};
    plot(tmp_excit, tmp_inhib, 'o', 'markerfacecolor', clr_raw, 'markeredgecolor', clr_raw, 'markersize', 10)
    l_nan = isnan(tmp_excit);
    betas = [tmp_excit(~l_nan), ones(size(tmp_excit(~l_nan)))] \ tmp_inhib(~l_nan);
    xvals = get(gca, 'xlim');
    plot(xvals, [xvals(:), ones(2,1)]*betas(:), '--', 'color', clr_raw, 'linewidth', 3)
    xlim(xvals_all)
    ylim(yvals_all)
    xlabel('Excit conductance (nS)')
    ylabel('Inhib conductance (nS)')
end





%% AMPA to NMDA RATIOS

fin

% load in the pre-saved population data
load([GL_POPDATPATH, 'popAnly_EIAN.mat'])


% create grouping lists
l_valid = dat.goodNeurons(:);
hvas = repmat(dat.hva, 2,1);
hvas = hvas(l_valid);
hvalist.('pm') = cellfun(@(x) ~isempty(x), regexpi(hvas, 'pm'));
hvalist.('lm') = cellfun(@(x) ~isempty(x), regexpi(hvas, 'lm'));
hvalist.('und') = cellfun(@(x) ~isempty(x), regexpi(hvas, 'und'));


% pull out peak conductances
raw_ampa = dat.ampa.peak_nS(:);
raw_ampa = raw_ampa(l_valid);
raw_nmda = dat.nmda.peak_nS(:);
raw_nmda = raw_nmda(l_valid);

% plot peak conductances for all cells
figure, hold on,
set(gca, 'fontsize', 25)
plot(raw_ampa, raw_nmda, 'ko', 'markerfacecolor', 'k', 'markersize', 10)
l_nan = isnan(raw_ampa);
betas = [raw_ampa(~l_nan), ones(size(raw_ampa(~l_nan)))] \ raw_nmda(~l_nan);
xvals_all = get(gca, 'xlim');
yvals_all = get(gca, 'ylim');
plot(xvals_all, [xvals_all(:), ones(2,1)]*betas(:), '--k', 'linewidth', 3)
xlabel('AMPA conductance (nS)')
ylabel('NMDA conductance (nS)')
title(sprintf('ALL DATA: A/N ratio = %.2f', betas(1)))

% plot the data, but color code each HVA
figure, hold on,
set(gca, 'fontsize', 25)
groups = {'und', 'pm', 'lm'};
for a = 1:numel(groups);
    tmp_ampa = raw_ampa(hvalist.(groups{a}));
    tmp_nmda = raw_nmda(hvalist.(groups{a}));
    [clr_fit, clr_raw] = hvaPlotColor(groups{a});
    plot(tmp_ampa, tmp_nmda, 'o', 'markerfacecolor', clr_raw, 'markeredgecolor', clr_raw, 'markersize', 10)
    l_nan = isnan(tmp_ampa);
    betas = [tmp_ampa(~l_nan), ones(size(tmp_ampa(~l_nan)))] \ tmp_nmda(~l_nan);
    xvals = get(gca, 'xlim');
    plot(xvals, [xvals(:), ones(2,1)]*betas(:), '--', 'color', clr_raw, 'linewidth', 3)
    xlim(xvals_all)
    ylim(yvals_all)
    xlabel('AMPA conductance (nS)')
    ylabel('NMDA conductance (nS)')
end



% now group by IN vs PY
neuronType = dat.neuronType(:);
neuronType = cellfun(@num2str, neuronType, 'uniformoutput', false);
neuronType = neuronType(l_valid);
typeList.IN = cellfun(@(x) ~isempty(x), regexpi(neuronType, 'in'));
typeList.PY = cellfun(@(x) ~isempty(x), regexpi(neuronType, 'py'));
typeList.und = ~typeList.IN & ~typeList.PY;

groupTypes = {'IN', 'PY', 'und'};
pltClr = {'r', 'b', [.85 .85 .85]};
figure, hold on,
set(gca, 'fontsize', 25)
title('By cell type')
for a = 1:numel(groupTypes);
    tmp_ampa = raw_ampa(typeList.(groupTypes{a}));
    tmp_nmda = raw_nmda(typeList.(groupTypes{a}));
    clr_raw = pltClr{a};
    plot(tmp_ampa, tmp_nmda, 'o', 'markerfacecolor', clr_raw, 'markeredgecolor', clr_raw, 'markersize', 10)
    l_nan = isnan(tmp_ampa);
    betas = [tmp_ampa(~l_nan), ones(size(tmp_ampa(~l_nan)))] \ tmp_nmda(~l_nan);
    xvals = get(gca, 'xlim');
    plot(xvals, [xvals(:), ones(2,1)]*betas(:), '--', 'color', clr_raw, 'linewidth', 3)
    xlim(xvals_all)
    ylim(yvals_all)
    xlabel('AMPA conductance (nS)')
    ylabel('NMDA conductance (nS)')
end



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














