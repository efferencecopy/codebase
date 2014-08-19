%% IMPORT THE NECESSARY DATA INTO A datDB STRUCTURE

fin

% buid a structure of params from physiology notes
global GL_ADD_TO_MDB GL_SUPPRESS_ANALYSIS
GL_ADD_TO_MDB = true;
GL_SUPPRESS_ANALYSIS = true;
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


% build a datDB_EIAN by performing the appropriate analysis. I'll use the
% params.name field to build each new entry into the DB. Finish by saving
% the DB (possibly with a date string)
% Things to save
%  * raw conductances
%  * raw currents?
%  * peak conductances
%  * list of errors
%      * Vclamp error too big

% iterate over the mice in the cell library (some mice get analyzed
% multiple times if there were multiple experiments per mouse)
[dat.ampa.peak.nS, dat.nmda.peak.nS, dat.excit.peak.nS, dat.inhib.peak.nS] = deal(nan(numel(mouseNames), 2));
[dat.ampa.peak.pA, dat.nmda.peak.pA, dat.excit.peak.pA, dat.inhib.peak.pA] = deal(nan(numel(mouseNames), 2));
for ex = 1:numel(mouseNames)
    
    ex_mouseName = mouseNames{ex};
    ex_siteNum = siteNumber{ex};
    mdb = initMouseDB(false, true);
    [~, idx] = mdb.search(ex_mouseName);
    
    params = mdb.mice{idx}.popAnly{ex_siteNum};
    params.fxns = {@anlyMod_optoIV, @anlyMod_EIbalance};
    params = invitroAnalysisOverview(params);
    close all
    
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
                    dat.(group{g}).peak.nS(ex, ch) = params.isolatedData.(group{g}).peak_nS{ch};                    
                    dat.(group{g}).peak.nA(ex, ch) = params.isolatedData.(group{g}).peak_pA{ch};
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

%% EXCITATION VS. INHIBITION

fin

% load in the pre-saved population data
load([GL_POPDATPATH, 'popAnly_EIAN.mat'])

% pull out raw data
l_valid = dat.goodNeurons(:);
raw_excit = dat.excit.peak.nS(:);
raw_excit = raw_excit(l_valid);
raw_inhib = dat.inhib.peak.nS(:);
raw_inhib = raw_inhib(l_valid);

% create grouping lists
hvas = repmat(dat.hva, 2,1);
hvas = hvas(l_valid);
hvalist.('pm') = cellfun(@(x) ~isempty(x), regexpi(hvas, 'pm'));
hvalist.('lm') = cellfun(@(x) ~isempty(x), regexpi(hvas, 'lm'));
hvalist.('und') = cellfun(@(x) ~isempty(x), regexpi(hvas, 'und'));


% plot the most raw form of the data
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

% pull out raw data
l_valid = dat.goodNeurons(:);
raw_ampa = dat.ampa.peak(:);
raw_ampa = raw_ampa(l_valid);
raw_nmda = dat.nmda.peak(:);
raw_nmda = raw_nmda(l_valid);

% create grouping lists
hvas = repmat(dat.hva, 2,1);
hvas = hvas(l_valid);
hvalist.('pm') = cellfun(@(x) ~isempty(x), regexpi(hvas, 'pm'));
hvalist.('lm') = cellfun(@(x) ~isempty(x), regexpi(hvas, 'lm'));
hvalist.('und') = cellfun(@(x) ~isempty(x), regexpi(hvas, 'und'));


% plot the most raw form of the data
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











