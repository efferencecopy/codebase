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
                end
            end
        end
    end
    
end

% package all the useful things into a single structure (and then save the
% structure)
dat.hva = HVA;
dat.goodNeurons = goodNeurons;
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
hvalist.('und') = cellfun(@(x) ~isempty(x), regexpi(hvas, 'und'));

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
plot(raw_excit, raw_inhib, 'ko', 'markerfacecolor', 'k')
l_nan = isnan(raw_excit);
betas = [raw_excit(~l_nan), ones(size(raw_excit(~l_nan)))] \ raw_inhib(~l_nan);
xvals = get(gca, 'xlim');
plot(xvals, [xvals(:), ones(2,1)]*betas(:), '--k')
xlabel('Excit conductance (nS)')
ylabel('Inhib conductance (nS)')
title(sprintf('ALL DATA: E/I ratio = %.2f', betas(1)))

% plot the data, but color code each HVA
figure, hold on,
groups = {'und', 'pm', 'lm'};
for a = 1:numel(groups);
    tmp_excit = raw_excit(hvalist.(groups{a}));
    tmp_inhib = raw_inhib(hvalist.(groups{a}));
    [clr_fit, clr_raw] = hvaPlotColor(groups{a});
    plot(tmp_excit, tmp_inhib, 'o', 'markerfacecolor', clr_raw, 'markeredgecolor', clr_raw)
    l_nan = isnan(tmp_excit);
    betas = [tmp_excit(~l_nan), ones(size(tmp_excit(~l_nan)))] \ tmp_inhib(~l_nan);
    xvals = get(gca, 'xlim');
    plot(xvals, [xvals(:), ones(2,1)]*betas(:), '--', 'color', clr_raw)
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
plot(raw_ampa, raw_nmda, 'ko', 'markerfacecolor', 'k')
l_nan = isnan(raw_ampa);
betas = [raw_ampa(~l_nan), ones(size(raw_ampa(~l_nan)))] \ raw_nmda(~l_nan);
xvals = get(gca, 'xlim');
plot(xvals, [xvals(:), ones(2,1)]*betas(:), '--k')
xlabel('AMPA conductance (nS)')
ylabel('NMDA conductance (nS)')
title(sprintf('ALL DATA: A/N ratio = %.2f', betas(1)))

% plot the data, but color code each HVA
figure, hold on,
groups = {'und', 'pm', 'lm'};
for a = 1:numel(groups);
    tmp_ampa = raw_ampa(hvalist.(groups{a}));
    tmp_nmda = raw_nmda(hvalist.(groups{a}));
    [clr_fit, clr_raw] = hvaPlotColor(groups{a});
    plot(tmp_ampa, tmp_nmda, 'o', 'markerfacecolor', clr_raw, 'markeredgecolor', clr_raw)
    l_nan = isnan(tmp_ampa);
    betas = [tmp_ampa(~l_nan), ones(size(tmp_ampa(~l_nan)))] \ tmp_nmda(~l_nan);
    xvals = get(gca, 'xlim');
    plot(xvals, [xvals(:), ones(2,1)]*betas(:), '--', 'color', clr_raw)
    xlabel('AMPA conductance (nS)')
    ylabel('NMDA conductance (nS)')
end


%% CONTROL ANALYSES

fin

% load in the pre-saved population data
load([GL_POPDATPATH, 'popAnly_EIAN.mat'])


% create grouping lists
l_valid = dat.goodNeurons(:);

% grab the data
tmp_ampa_verr = dat.ampa.Verr(:);
tmp_ampa_verr = tmp_ampa_verr(l_valid);



