%% SPECIFY WHICH EXPERIMENTS SHOULD CONTRIBUTE, LOAD THE DATA

fin


% decide what experiment to run
EXPTTYPE = 2;
switch EXPTTYPE
    case 1
        EXPTTYPE = 'all';
    case 2
        EXPTTYPE = 'Manifold';
    case 3
        EXPTTYPE = 'RIT_test';
    case 4
        EXPTTYPE = 'Passive_Props';
    case 5
        EXPTTYPE = 'LGN';
    case 6
        EXPTTYPE = 'IN_strength';
    case 7
        EXPTTYPE = 'Zucker_Stim';
    case 8
        EXPTTYPE = 'IN_powers';
end



%%%% DEFINE THE ANALYSIS PARAMS %%%%

params.pretime.vclamp = 0.002;     % seconds before pulse onset
params.posttime.vclamp = 0.015;    % seconds after pulse onset
params.pretime.dcsteps = 0.100;    % seconds before current onset
params.posttime.dcsteps = 0.300;   % seconds after current offset 
params.pretime.iclamp = 0.003;
params.posttime.iclamp = 0.019;    % actually a minimum value, real value stored in the dat structure, and depends on ISI
params.expttype = EXPTTYPE;
switch EXPTTYPE
    case {'IN_strength', 'IN_powers'}
        params.force_epsc = true; % allows small amplitude SOM EPSCs to be analyzed correctly
    otherwise
        params.force_epsc = false;
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% figure out which experiments should be analyzed
wb_path = [GL_DOCUPATH, 'Other_workbooks', filesep, 'wholeCellSTPCellList.xlsx'];
[~, ~, wb_expt] = xlsread(wb_path, 2);
if ~strcmpi(EXPTTYPE, 'all')
    header_idx = strcmpi(EXPTTYPE, wb_expt(1,:));
    assert(sum(header_idx) == 1)
    l_to_analyze = cellfun(@(x) numel(x)==1 && (x==1 || x=='1'), wb_expt(:, header_idx));
else
    l_to_analyze = true(size(wb_expt,1), 1);
    l_to_analyze(1) = false; % don't analyze the header row
end
expt_idx = find(l_to_analyze);

% load in the workbook that contains all the experimental information
[~,~,wb_info] = xlsread(wb_path, 1);

% generate the header index info
for i_atrib = 1:size(wb_info,2)
    fldname = wb_info{1,i_atrib};
    fldname(isspace(fldname)) = [];
    hidx.(fldname) = i_atrib;
end

% throw an error if the wb_expt and wb_info have different numbers of
% neurons
assert(size(wb_expt,1)==size(wb_info,1), 'ERROR: unequal numbers of neurons in wb_info and wb_expt')
    
% now that the header is formed, delete the first row from the wb_info and
% subtract one from the expt_idx to "delete" the first row from the wb_expt
wb_info(1,:) = [];
expt_idx = expt_idx - 1;

% convert the file names into fully qualified paths so that par-for can run
% without calling a global
dcsteps_HS1_fpath = cellfun(@(x,y) strcat(GL_DATPATH, x, filesep, 'Physiology', filesep, y, '.abf'), wb_info(:,hidx.MouseName), wb_info(:, hidx.ABFDCstepsHS1), 'uniformoutput', false);
dcsteps_HS2_fpath = cellfun(@(x,y) strcat(GL_DATPATH, x, filesep, 'Physiology', filesep, y, '.abf'), wb_info(:,hidx.MouseName), wb_info(:, hidx.ABFDCstepsHS2), 'uniformoutput', false);
vclamp_fpath = cellfun(@(x,y) strcat(GL_DATPATH, x, filesep, 'Physiology', filesep, y, '.abf'), wb_info(:,hidx.MouseName), wb_info(:, hidx.ABFOptostimVclamp), 'uniformoutput', false);
iclamp_fpath = cellfun(@(x,y) strcat(GL_DATPATH, x, filesep, 'Physiology', filesep, y, '.abf'), wb_info(:,hidx.MouseName), wb_info(:, hidx.ABFOptostimIclamp), 'uniformoutput', false);
wb_info(:,hidx.ABFDCstepsHS1) = dcsteps_HS1_fpath;
wb_info(:,hidx.ABFDCstepsHS2) = dcsteps_HS2_fpath;
wb_info(:,hidx.ABFOptostimVclamp) = vclamp_fpath;
wb_info(:,hidx.ABFOptostimIclamp) = iclamp_fpath;
                   
% make each row it's own cell array so that it can be passed as a single
% argument to the function that does the major unpacking
attributes = {};
for i_ex = 1:numel(expt_idx)
    attributes{i_ex,1} = wb_info(expt_idx(i_ex),:);
end


%
% LOAD THE DATA FILES
%
dat = {};
Nexpts = numel(attributes);

pool = gcp('nocreate');
if isempty(pool)
    pool = parpool();
end
for i_ex_par = 200:Nexpts
    dat{i_ex_par} = wcstp_compile_data(attributes{i_ex_par}, hidx, params);
    close all
end
fprintf('All done importing data\n')

%% FORCE THE DATASET TO ONLY INCLUDE "COMPLETE" SETS

% define a complete set as:
% having 12, 25, 50 Hz
% having 3 recovery times [0.33, 1, 5.5]
% a minimum of 3 trials per condition (33 total)
%
% enforce this by changing the "isvalid" flag
% if both channels are ~ valid, then delete the cell array in "dat"

l_good_expts = false(numel(dat), 1);
any_is_valid_vclamp = false(numel(dat), 1);
for i_ex = 1:numel(dat)
    % reassemble the tdict array, and pull out trial counts for each channel
    conds = fieldnames(dat{i_ex}.expt);
    tdict = nan(numel(conds), 5);
    has_3trls = nan(numel(conds), 2);
    for i_cond = 1:numel(conds)
        tdict(i_cond,:) = dat{i_ex}.expt.(conds{i_cond}).tdict;
        
        has_data = cellfun(@(x) ~isempty(x), dat{i_ex}.expt.(conds{i_cond}).stats.EPSCamp);
        gt_3trls = cellfun(@(x) size(x, 3), dat{i_ex}.expt.(conds{i_cond}).stats.EPSCamp) >= 3;
        has_3trls(i_cond,:) = has_data & gt_3trls;
    end
    has_3trls = all(has_3trls, 1);
    
    % make sure the t-types are correct
    unique_tfs = unique(tdict(:,3));
    unique_tfs(unique_tfs==0) = [];
    if numel(unique_tfs) == 3
        correct_tfs = all(unique_tfs(:) == [12; 25; 50]);
    else
        correct_tfs = false;
    end
    correct_tfs = [correct_tfs, correct_tfs];
    
    unique_recovs =unique(tdict(:,4));
    unique_recovs(unique_recovs == 0) = [];
    if numel(unique_recovs) == 3
        correct_recovs = all(unique_recovs(:) == [333; 1000; 5500]);
    else
        correct_recovs = false;
    end
    correct_recovs = [correct_recovs, correct_recovs];
    
    ch_good = has_3trls & correct_tfs & correct_recovs;
    l_good_expts(i_ex) = any(ch_good);
    
    % change the identity of the is_valid flag, but only for the cases
    % where the "validitiy" was initially true
    tmp_is_valid_vclamp = dat{i_ex}.info.HS_is_valid_Vclamp;
    for i_ch = 1:2
        if tmp_is_valid_vclamp(i_ch)
            if ~ch_good(i_ch)
                tmp_is_valid_vclamp(i_ch) = false;
            end
        end
    end
    dat{i_ex}.info.HS_is_valid_Vclamp = tmp_is_valid_vclamp;
    any_is_valid_vclamp(i_ex) = any(tmp_is_valid_vclamp);
end

dat = dat(any_is_valid_vclamp);
fprintf('deleted %d experiments\n', sum(~any_is_valid_vclamp))

%% FILTER OUT MICE THAT DO NOT HAVE MULTIPLE HVAs REPRESENTED

all_mouse_names = {};
for i_ex = 1:numel(dat)
    all_mouse_names = cat(1, all_mouse_names, dat{i_ex}.info.mouseName);
end

unique_mouse_names = unique(all_mouse_names);
l_mice_with_2HVAs = false(size(attributes));
for i_mouse = 1:numel(unique_mouse_names)
    l_mouse = strcmpi(all_mouse_names, unique_mouse_names{i_mouse});
    
    % by definition, only 1 HVA is  no good
    if sum(l_mouse) <= 1; continue; end
    
    % make a matrix that shows hva & celltype, for each valid Vclamp
    % recording
    [ex_hva, ex_celltype] = deal({});
    for i_ex = find(l_mouse)'
        for i_ch = 1:2
            isvalid = dat{i_ex}.info.HS_is_valid_Vclamp(i_ch);
            if isvalid
                ex_hva = cat(1, ex_hva, dat{i_ex}.info.brainArea);
                ex_celltype = cat(1, ex_celltype, dat{i_ex}.info.cellType{i_ch});
            end
        end
    end
    
    % eliminate the non-PY cells
    l_py_L23 = strcmpi(ex_celltype, 'py_l23');
    ex_hva_from_py_cells = ex_hva(l_py_L23);
    l_lat = cellfun(@(x) ~isempty(x), regexpi(ex_hva_from_py_cells, 'lm|al'));
    l_med = cellfun(@(x) ~isempty(x), regexpi(ex_hva_from_py_cells, 'pm|am'));
    
    % count this as a good mouse if there is at least one of each med/lat
    if sum(l_lat)>0 && sum(l_med)>0
        l_mice_with_2HVAs(l_mouse) = true;
    end
end

dat = dat(l_mice_with_2HVAs);
fprintf('deleted %d experiments\n', sum(~l_mice_with_2HVAs))


%% QULAITY CONTROL PLOTS

close all

for i_ex = 180:numel(dat)

    f = figure;
    f.Name = sprintf('Mouse %s, site %s, HVA: %s', dat{i_ex}.info.mouseName, dat{i_ex}.info.siteNum, dat{i_ex}.info.brainArea);
    f.Position = [332  96 1259 665];
    
    for i_ch = 1:2
        
        % plot the current step data set to help identify cell types
        if ~isempty(dat{i_ex}.dcsteps.Vm_raw{i_ch})
            Vm = dat{i_ex}.dcsteps.Vm_raw{i_ch};
            Icmd = dat{i_ex}.dcsteps.Icmd{i_ch};
            N = size(Vm,2);
            tt = ([0:N-1] ./ dat{i_ex}.info.sampRate.dcsteps) - dat{i_ex}.info.pretime.dcsteps;
            if i_ch == 1; col = 1; else col = 3; end
            
            if any(Icmd<0)
                pltidx = sub2ind([4,3], col, 3);
                ha = subplot(3,4,pltidx); hold on,
                neg_cmds = unique(Icmd(Icmd < 0));
                for i_cmd = 1:numel(neg_cmds)
                    l_cmd = Icmd == neg_cmds(i_cmd);
                    plot(tt, mean(Vm(l_cmd,:), 1))
                end
                if isfield(dat{i_ex}.dcsteps, 'Ih_sag') && ~isempty(dat{i_ex}.dcsteps.Ih_sag.sag{i_ch})
                    peak_vals = repmat(dat{i_ex}.dcsteps.Ih_sag.peak_Vm{i_ch}, 2, 1);
                    sag_amps = repmat(dat{i_ex}.dcsteps.Ih_sag.sag{i_ch}, 2, 1);
                    asym_vals = repmat(dat{i_ex}.dcsteps.Ih_sag.Vm_asym{i_ch}, 2, 1);
                    plot([tt(1), tt(end)], peak_vals, 'm:')
                    plot([tt(1), tt(end)], peak_vals-sag_amps, 'r:')
                    plot([tt(1), tt(end)], asym_vals, 'g:')
                    plot(repmat([12e-3, 25e-3], 2, 1), [min(Vm(:)), dat{i_ex}.dcsteps.Vrest{i_ch}], 'k--')
                end
                axis tight
                ha.Box = 'off';
                ha.TickDir = 'out';
                xlabel('Seconds')
                ylabel('mV')
            end
            
            vmidx = (Icmd>0) & (Icmd<350);
            if any(vmidx);
                pltidx = sub2ind([4,3], col, 2);
                ha = subplot(3,4,pltidx);
                plot(tt, Vm(vmidx,:)');
                axis tight
                ha.Box = 'off';
                ha.TickDir = 'out';
                ylabel('mV')
            end
            
            vmidx = (Icmd >= 350);
            if any(vmidx);
                pltidx = sub2ind([4,3], col, 1);
                ha = subplot(3,4,pltidx);
                plot(tt, Vm(vmidx,:)');
                axis tight
                ha.Box = 'off';
                ha.TickDir = 'out';
                ylabel('mV')
                title(dat{i_ex}.info.cellType{i_ch})
            end
        end
        
        % figure out the subplot column number for the following axes
        if i_ch == 1; col = 2; else col = 4; end
        
        % series resistance
        if isfield(dat{i_ex}.qc, 'Rs') && dat{i_ex}.info.HS_is_valid_Vclamp(i_ch)
            pltidx = sub2ind([4,3], col, 1);
            subplot(3,4,pltidx)
            tmp = squeeze(dat{i_ex}.qc.Rs{i_ch});
            plot(tmp)
            ylim([min([0,min(tmp(:))]), max(tmp(:))*1.5])
            ylabel('R_{s} (MOhm)')
            xlabel('trial number')
            title(sprintf('Channel %d', i_ch))
        end
        
        % verr
        if isfield(dat{i_ex}.qc, 'Verr') && dat{i_ex}.info.HS_is_valid_Vclamp(i_ch)
            pltidx = sub2ind([4,3],col,2);
            subplot(3,4,pltidx)
            tmp = squeeze(dat{i_ex}.qc.verr{i_ch});
            plot(tmp)
            ylim([min([0,min(tmp(:))]), max(tmp(:))*1.5])
            ylabel('SS Verr (mV)')
            xlabel('trial number')
            
        end
        
        % p1amps
        if isfield(dat{i_ex}.qc, 'p1amp') && dat{i_ex}.info.HS_is_valid_Vclamp(i_ch)
            pltidx = sub2ind([4,3], col, 3);
            subplot(3,4,pltidx), hold on,
            tmp = squeeze(dat{i_ex}.qc.p1amp{i_ch});
            tmp_norm = squeeze(dat{i_ex}.qc.p1amp_norm{i_ch});
            plot(tmp)
            plot(tmp_norm, 'r')
            ylim([min([0,min(tmp(:))]), max(tmp(:))*1.5])
            ylabel('P1 Amp')
            xlabel('trial number')
        end
            
        
    end
    drawnow
    

    
end


%% SUMMARY OF PASSIVE PROPERTIES

close all; clc

% define a set of attributes for each analysis group
% {CellType, Layer,  BrainArea,  OpsinType}
% Brain Area can be: 'AL', 'PM', 'AM', 'LM', 'any', 'med', 'lat'. CASE SENSITIVE
plotgroups = {
    'PY', 'L23', 'AM', 'any';...
    'PY', 'L23', 'PM', 'any';...
    'PY', 'L23', 'AL', 'any';...
    'PY', 'L23', 'LM', 'any';...
    };

groupdata.Rin_peak = repmat({[]}, 1, size(plotgroups, 1)); % should only have N cells, where N = size(plotgroups, 1). Each cell has a matrix with a cononicalGrid:
groupdata.Rin_asym = repmat({[]}, 1, size(plotgroups, 1));
groupdata.tau = repmat({[]}, 1, size(plotgroups, 1));
groupdata.Vrest = repmat({[]}, 1, size(plotgroups, 1));
groupdata.Depth = repmat({[]}, 1, size(plotgroups, 1));
groupdata.IVcurve_peak = repmat({{}}, 1, size(plotgroups,1));
groupdata.IVcurve_asym = repmat({{}}, 1, size(plotgroups,1));
groupdata.Ih_sag = repmat({{}}, 1, size(plotgroups,1));
groupdata.Ih_Vm = repmat({{}}, 1, size(plotgroups,1));
groupdata.Ih_asym = repmat({{}}, 1, size(plotgroups,1));
groupdata.starttime = repmat({[]}, 1, size(plotgroups,1));


for i_ex = 1:numel(dat)
    
    for i_ch = 1:2
        
        % check to make sure this neuron was defined
        hs_name = sprintf('dcsteps_hs%d', i_ch);
        isvalid = ~strncmp(dat{i_ex}.info.fid.(hs_name)(end-8:end), [filesep, 'none.abf'], 9);
        if ~isvalid
            continue
        end
        
        % check to make sure the cell was healthy (as based on presence of
        % holding current)
        max_holding_pa = abs(dat{i_ex}.dcsteps.pA_holding{i_ch}(2));
        if max_holding_pa > 20; disp(max_holding_pa); continue; end
        
        % check the attributes
        ch_attribs = {dat{i_ex}.info.cellType{i_ch}, dat{i_ex}.info.brainArea, dat{i_ex}.info.opsin};
        group_idx = groupMatcher(plotgroups, ch_attribs);
        if sum(group_idx) == 0; continue; end
        
        % add the statistics to the appropriate cell array
        groupdata.Rin_peak{group_idx} = cat(1, groupdata.Rin_peak{group_idx}, dat{i_ex}.dcsteps.IVpeak.Rin{i_ch});
        groupdata.Rin_asym{group_idx} = cat(1, groupdata.Rin_asym{group_idx}, dat{i_ex}.dcsteps.IVasym.Rin{i_ch});
        groupdata.Vrest{group_idx} = cat(1, groupdata.Vrest{group_idx}, dat{i_ex}.dcsteps.Vrest{i_ch});
        groupdata.Depth{group_idx} = cat(1, groupdata.Depth{group_idx}, dat{i_ex}.info.cellDepth_um(i_ch));
        groupdata.IVcurve_peak{group_idx} = cat(1, groupdata.IVcurve_peak{group_idx}, {dat{i_ex}.dcsteps.IVpeak.raw{i_ch}});
        groupdata.IVcurve_asym{group_idx} = cat(1, groupdata.IVcurve_asym{group_idx}, {dat{i_ex}.dcsteps.IVasym.raw{i_ch}});
        groupdata.starttime{group_idx} = cat(1, groupdata.starttime{group_idx}, dat{i_ex}.info.fileStartTime_24hrs);
        if isfield(dat{i_ex}.dcsteps, 'Ih_sag')
            groupdata.Ih_sag{group_idx} = cat(1, groupdata.Ih_sag{group_idx}, [dat{i_ex}.dcsteps.Ih_sag.sag{i_ch}]);
            groupdata.Ih_Vm{group_idx} = cat(1, groupdata.Ih_Vm{group_idx}, [dat{i_ex}.dcsteps.Ih_sag.peak_Vm{i_ch}]);
            groupdata.Ih_asym{group_idx} = cat(1, groupdata.Ih_asym{group_idx}, [dat{i_ex}.dcsteps.Ih_sag.Vm_asym{i_ch}]);
        end
        if isfield(dat{i_ex}.dcsteps, 'tau') && ~isempty(dat{i_ex}.dcsteps.tau{i_ch})
            tmp_all_tau = dat{i_ex}.dcsteps.tau{i_ch}(:,3);
            l_neg = dat{i_ex}.dcsteps.tau{i_ch}(:,1) < 0;
            groupdata.tau{group_idx} = cat(1, groupdata.tau{group_idx}, mean(tmp_all_tau(l_neg)));
        end
    end
end


%
% plot histograms of Input resistance measured two different ways
%
%%%%%%%%%%%%%%%%%%%%%%%%%
f=figure;
Ngroups = (numel(groupdata.Rin_peak)); % adding one for "summary" cdf fig
N_plt_rows = Ngroups+1;
f.Position = [506   158   838   749];
allNums = cat(1, groupdata.Rin_peak{:}, groupdata.Rin_asym{:});
edges = linspace(min(allNums)-1, max(allNums)+1, 50);
for i_group = 1:numel(groupdata.Rin_asym)
    
    hva_name = plotgroups{i_group, 3};
    grp_clr = hvaPlotColor(hva_name);
    
    %
    % current clamp, peak vals
    %%%%%%%%%%%%%%%%%%
    pltidx = sub2ind([2, N_plt_rows], 1, i_group);
    subplot(N_plt_rows, 2, pltidx)
    
    h = histogram(groupdata.Rin_peak{i_group}, edges);
    h.FaceColor = grp_clr;
    xbar = nanmean(groupdata.Rin_peak{i_group});
    hold on,
    plot(xbar, 0.5, 'kv', 'markerfacecolor', 'k', 'markersize', 5)
    N = numel(groupdata.Rin_peak{i_group});
    legtext = sprintf('%s, %s, %s, n=%d', plotgroups{i_group,1},plotgroups{i_group, 2},plotgroups{i_group,3}, N);
    hl = legend(legtext, 'Location', 'northeast');
    legend boxoff
    hl.Interpreter = 'none';
    xlabel('Input R (MOhm)')
    ylabel('counts')
    if i_group == 1
        title('Rin (iClamp, peak)')
    end
    
    % CDF summary
    pltidx = sub2ind([2, N_plt_rows], 1, Ngroups+1);
    subplot(N_plt_rows, 2, pltidx), hold on,
    N = histcounts(groupdata.Rin_peak{i_group}, edges);
    N(end+1) = 0;
    cdf_vals = cumsum(N)./sum(N);
    stairs(edges, cdf_vals, '-', 'color', grp_clr)
    xlabel('Input R (MOhm)')
    ylabel('Proportion')
    
    %
    % current clamp asym vals
    %%%%%%%%%%%%%%%%%%%%%%%%
    pltidx = sub2ind([2, N_plt_rows], 2, i_group);
    subplot(N_plt_rows, 2, pltidx)
    
    h = histogram(groupdata.Rin_asym{i_group}, edges);
    h.FaceColor = grp_clr;
    xbar = nanmean(groupdata.Rin_asym{i_group});
    hold on,
    plot(xbar, 0.5, 'kv', 'markerfacecolor', 'k', 'markersize', 5)
    xlabel('Input R (MOhm)')
    ylabel('counts')
    if i_group == 1
        title('Rin (iClamp, asym)')
    end
    
    % CDF summary
    pltidx = sub2ind([2, N_plt_rows], 2, Ngroups+1);
    subplot(N_plt_rows, 2, pltidx), hold on,
    N = histcounts(groupdata.Rin_asym{i_group}, edges);
    N(end+1) = 0;
    cdf_vals = cumsum(N)./sum(N);
    stairs(edges, cdf_vals, '-', 'color', grp_clr)
    xlabel('Input R (MOhm)')
    ylabel('Proportion')
end

%
% membrane time constants
%
%%%%%%%%%%%%%%%%%%%%%%%%%%
figure, hold on,
allTaus = cat(1, groupdata.tau{:});
edges_taus = linspace(min(allTaus)-0.002, 0.050, 50).*1000;
legtext = {};
for i_group = 1:numel(groupdata.tau)
    hva_name = plotgroups{i_group, 3};
    grp_clr = hvaPlotColor(hva_name);
    
    N = histcounts(groupdata.tau{i_group}.*1000, edges_taus);
    N(end+1) = 0;
    cdf_vals = cumsum(N)./sum(N);
    stairs(edges_taus, cdf_vals, '-', 'color', grp_clr)
    num_cells = numel(groupdata.tau{i_group});
    legtext{i_group} = sprintf('%s, %s, %s, n=%d', plotgroups{i_group,1},plotgroups{i_group, 2},plotgroups{i_group,3}, num_cells);
end
hl = legend(legtext, 'Location', 'northeast');
legend boxoff
hl.Interpreter = 'none';
xlabel('Time constant (ms)')
ylabel('Proportion')

%
% plot a histograms of Vrest and depth for each group
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=figure;
Ngroups = (numel(groupdata.Rin_peak));
f.Position = [466   464   835   440];
for i_group = 1:numel(groupdata.Vrest)
    
    hva_name = plotgroups{i_group, 3};
    grp_clr = hvaPlotColor(hva_name);
    
    % Vrest
    allNums = cat(1, groupdata.Vrest{:});
    edges = linspace(min(allNums)-10, max(allNums)+10, 30);
    
    pltidx = sub2ind([2, Ngroups], 1, i_group);
    subplot(Ngroups, 2, pltidx)
    
    h = histogram(groupdata.Vrest{i_group}, edges);
    h.FaceColor = grp_clr;
    xbar = nanmean(groupdata.Vrest{i_group});
    hold on,
    plot(xbar, 0.5, 'kv', 'markerfacecolor', 'k', 'markersize', 5)
    legtext = sprintf('%s, %s, %s', plotgroups{i_group,1},plotgroups{i_group, 2},plotgroups{i_group,3} );
    legend(legtext, 'Location', 'northeast')
    legend boxoff
    xlabel('Vrest(mV)')
    ylabel('counts')
    if i_group == 1
        title('Vrest')
    end
    
    % Cell depth
    allNums = cat(1, groupdata.Depth{:});
    edges = linspace(min(allNums)-10, max(allNums)+10, 30);
    
    pltidx = sub2ind([2, Ngroups], 2, i_group);
    subplot(Ngroups, 2, pltidx)
    
    h = histogram(groupdata.Depth{i_group}, edges);
    h.FaceColor = grp_clr;
    xbar = nanmean(groupdata.Depth{i_group});
    hold on,
    plot(xbar, 0.5, 'kv', 'markerfacecolor', 'k', 'markersize', 5)
    xlabel('Cell Depth (um))')
    ylabel('counts')
    if i_group == 1
        title('Cell Depth')
    end
    
    
    
end



%
% scatter plots of cell depth vs. Rin, Vrest
%
%%%%%%%%%%%%%%%%%%%%%%%%%
f=figure;
Ngroups = (numel(groupdata.Rin_asym));
f.Position = [49         305        1333         602];
for i_group = 1:numel(groupdata.Rin_asym)
    
    hva_name = plotgroups{i_group, 3};
    grp_clr = hvaPlotColor(hva_name);
    
    % for axis standardization
    allRin = cat(1, groupdata.Rin_asym{:});
    allDepth = cat(1, groupdata.Depth{:});
    allTimes = cat(1, groupdata.starttime{:});
    allVrest = cat(1, groupdata.Vrest{:});
    
    
    % cell depth vs. Rin
    pltidx = sub2ind([3, Ngroups], 1, i_group);
    subplot(Ngroups, 3, pltidx), hold on,
    X = groupdata.Depth{i_group};
    Y = groupdata.Rin_asym{i_group};
    [B,BINT] = regress(Y(:), [ones(size(X(:))), X(:)]);
    [top_int, bot_int, x_mod, y_mod] = regression_line_ci(0.05, B, X, Y);
    CI_up = top_int - y_mod;
    CI_down = y_mod - bot_int;
    plot(X, Y, 'o', 'color', grp_clr, 'markerfacecolor', grp_clr, 'markersize', 5);
    shadedErrorBar(x_mod, y_mod, [CI_up ; CI_down], {'color', grp_clr});
    legtext = sprintf('%s, %s, %s', plotgroups{i_group,1},plotgroups{i_group, 2},plotgroups{i_group,3} );
    legend(legtext, 'Location', 'best')
    legend boxoff
    xlabel('Cell Depth (um)')
    ylabel('Input Resistance')
    if i_group == 1
        title('Rin vs. Depth')
    end
    xlim([min(allDepth).*.95, max(allDepth).*1.05])
    ylim([min(allRin).*.95, max(allRin).*1.05])
    
    % incubation period vs. Rin
    pltidx = sub2ind([3, Ngroups], 2, i_group);
    subplot(Ngroups, 3, pltidx), hold on,
    
    X = groupdata.starttime{i_group};
    Y = groupdata.Rin_asym{i_group};
    [B,BINT] = regress(Y(:), [ones(size(X(:))), X(:)]);
    [top_int, bot_int, x_mod, y_mod] = regression_line_ci(0.05, B, X, Y);
    CI_up = top_int - y_mod;
    CI_down = y_mod - bot_int;
    plot(X, Y, 'o', 'color', grp_clr, 'markerfacecolor', grp_clr, 'markersize', 5)
    shadedErrorBar(x_mod, y_mod, [CI_up ; CI_down], {'color', grp_clr});
    xlabel('start time (24hrs)')
    ylabel('Input Resistance')
    if i_group == 1
        title('Rin vs. Incubation Time')
    end
    xlim([min(allTimes).*.95, max(allTimes).*1.05])
    ylim([min(allRin).*.95, max(allRin).*1.05])
    
    % scatter plot of Vrest vs depth
    pltidx = sub2ind([3, Ngroups], 3, i_group);
    subplot(Ngroups, 3, pltidx), hold on,
    

    Y = groupdata.Rin_asym{i_group};
    X = groupdata.Vrest{i_group};
    [B,BINT] = regress(Y(:), [ones(size(X(:))), X(:)]);
    [top_int, bot_int, x_mod, y_mod] = regression_line_ci(0.05, B, X, Y);
    CI_up = top_int - y_mod;
    CI_down = y_mod - bot_int;
    plot(X, Y, 'o', 'color', grp_clr, 'markerfacecolor', grp_clr, 'markersize', 5)
    shadedErrorBar(x_mod, y_mod, [CI_up ; CI_down], {'color', grp_clr});
    xlabel('V rest (mV)')
    ylabel('Rin asym (Mohm)')
    if i_group == 1
        title('Rin vs. Vrest')
    end
    ylim([min(allRin).*.95, max(allRin).*1.05])
    xlim([min(allVrest).*.95, max(allVrest).*1.05])
    
end



%
% Line plots of I-V curves, and peak vs. steady state
%
%%%%%%%%%%%%%%%%%%%%%%%%%
f=figure;
Ngroups = (numel(groupdata.IVcurve_peak));
f.Position = [114 426 1667 456];
legtext = {};
leghand = [];

% I-V curves for peak values
subplot(1,3,1), hold on,
for i_group = 1:numel(groupdata.IVcurve_peak)
    hva_name = plotgroups{i_group, 3};
    grp_clr = hvaPlotColor(hva_name);
    
    for i_cell = 1:numel(groupdata.IVcurve_peak{i_group})
        X = groupdata.IVcurve_peak{i_group}{i_cell}(:,1);
        Y = groupdata.IVcurve_peak{i_group}{i_cell}(:,2);
        hp = plot(X, Y, '-', 'color', grp_clr);
    end
    legtext{i_group} = sprintf('%s, %s, %s', plotgroups{i_group,1},plotgroups{i_group, 2},plotgroups{i_group,3} );
    leghand(i_group) = hp;
end
title('IV curve, peak vals')
xlabel('current injection (pA)')
ylabel('voltage response (mV)')
legend(leghand, legtext, 'location', 'best')
legend boxoff

% I-V curves for asym values
subplot(1,3,2), hold on,
for i_group = 1:numel(groupdata.IVcurve_asym)
    hva_name = plotgroups{i_group, 3};
    grp_clr = hvaPlotColor(hva_name);
    
    for i_cell = 1:numel(groupdata.IVcurve_asym{i_group})
        X = groupdata.IVcurve_asym{i_group}{i_cell}(:,1);
        Y = groupdata.IVcurve_asym{i_group}{i_cell}(:,2);
        plot(X, Y, '-', 'color', grp_clr)
    end
end
title('IV curve, asym vals')
xlabel('current injection (pA)')
ylabel('voltage response (mV)')

% directly compare peak vs. steady state DC response 
subplot(1,3,3), hold on,
edges = [-825:25:625];
group_iv_vals = {};
for i_group = 1:numel(groupdata.IVcurve_peak)
    hva_name = plotgroups{i_group, 3};
    grp_clr = hvaPlotColor(hva_name);
    
    vals = repmat({[]}, 1, numel(edges));
    for i_cell = 1:numel(groupdata.IVcurve_peak{i_group})
        X = groupdata.IVcurve_peak{i_group}{i_cell}(:,1);
        Y_peak = groupdata.IVcurve_peak{i_group}{i_cell}(:,2);
        Y_asym = groupdata.IVcurve_asym{i_group}{i_cell}(:,2);
        diff_val = Y_asym-Y_peak;
        plot(X, diff_val, '-', 'color', grp_clr)
        for i_val = 1:numel(diff_val)
            idx = histc(X(i_val), edges);
            vals{logical(idx)}(end+1) = diff_val(i_val);
        end
    end
    group_iv_vals{i_group} = vals;
end
axis tight
title('Steady State - Peak')
xlabel('DC injection (pA)')
ylabel('Difference (mV)')


% ----- regression of IV curve --------------
hf = figure;
hold on,
for i_group = 1:numel(groupdata.IVcurve_peak)
    hva_name = plotgroups{i_group, 3};
    grp_clr = hvaPlotColor(hva_name);
    
    [X,Y] = deal([]);
    for i_val = 1:numel(edges);
        if edges(i_val)<=0; continue; end
        tmpvals = group_iv_vals{i_group}{i_val};
        Y = cat(1,Y, tmpvals(:));
        X = cat(1,X, ones(numel(tmpvals), 1).*edges(i_val));
    end
    [B,BINT] = regress(Y(:), [ones(size(X(:))), X(:)]);
    [top_int, bot_int, x_mod, y_mod] = regression_line_ci(0.05, B, X, Y);
    CI_up = top_int - y_mod;
    CI_down = y_mod - bot_int;
    plot(X, Y, '.', 'color', grp_clr)
    shadedErrorBar(x_mod, y_mod, [CI_up ; CI_down], {'color', grp_clr});
end
xlabel('Current Injection (pA)')
ylabel('Difference (mV)')
title('Steady State - Peak')


% -----  compare Ih sag --------- 
legtext = {};
leghand = [];
f = figure;
f.Position = [384 353 1071 415];
subplot(1,2,1), hold on,
edges = [-150:10:-80];
vals = repmat({[]}, 1, numel(edges));
for i_group = 1:numel(groupdata.Ih_sag)
    hva_name = plotgroups{i_group, 3};
    grp_clr = hvaPlotColor(hva_name);
    N = numel(groupdata.Ih_sag{i_group});
    for i_ex = 1:N
        hp = plot(groupdata.Ih_Vm{i_group}{i_ex}, groupdata.Ih_sag{i_group}{i_ex}, '.-', 'color', grp_clr);
    end
    legtext{i_group} = sprintf('%s, %s, %s', plotgroups{i_group,1},plotgroups{i_group, 2},plotgroups{i_group,3} );
    leghand(i_group) = hp;
end
for i_group = 1:numel(groupdata.Ih_sag)
    hva_name = plotgroups{i_group, 3};
    grp_clr = hvaPlotColor(hva_name);
    N = numel(groupdata.Ih_sag{i_group});
    for i_ex = 1:N
        tmpval = groupdata.Ih_sag{i_group}{i_ex};
        tmpmv = groupdata.Ih_Vm{i_group}{i_ex};
        for i_val= 1:numel(tmpval)
            idx = histc(tmpmv(i_val), edges);
            vals{logical(idx)}(end+1) = tmpval(i_val);
        end
    end
    l_enough = cellfun(@numel, vals) > 5;
    avg_ih = cellfun(@nanmean, vals);
    plot(edges(l_enough), avg_ih(l_enough), '-', 'color', grp_clr, 'linewidth', 4)
end
xlabel('Membrane Potential At Sag Peak (mV)')
ylabel('Sag Amplitude (mV)')
legend(leghand, legtext, 'location', 'best')
legend boxoff

subplot(1,2,2), hold on,
for i_group = 1:numel(groupdata.Ih_asym)
    hva_name = plotgroups{i_group, 3};
    grp_clr = hvaPlotColor(hva_name);
    N = numel(groupdata.Ih_asym{i_group});
    for i_ex = 1:N
        sag_Vm = groupdata.Ih_Vm{i_group}{i_ex} - groupdata.Ih_sag{i_group}{i_ex};
        asym_Vm =  groupdata.Ih_asym{i_group}{i_ex};
        rebound = sag_Vm - asym_Vm;
        plot(groupdata.Ih_Vm{i_group}{i_ex}, rebound, 'o-', 'color', grp_clr)
    end
end

edges = [-150:10:-80];
group_iv_vals = {};
for i_group = 1:numel(groupdata.Ih_asym)
    hva_name = plotgroups{i_group, 3};
    grp_clr = hvaPlotColor(hva_name);
    
    vals = repmat({[]}, 1, numel(edges));
    N = numel(groupdata.Ih_asym{i_group});
    for i_ex = 1:N
        sag_Vm = groupdata.Ih_Vm{i_group}{i_ex} - groupdata.Ih_sag{i_group}{i_ex};
        asym_Vm =  groupdata.Ih_asym{i_group}{i_ex};
        rebound = sag_Vm - asym_Vm;
        tmpmv = groupdata.Ih_Vm{i_group}{i_ex};
        for i_val= 1:numel(rebound)
            idx = histc(tmpmv(i_val), edges);
            vals{logical(idx)}(end+1) = rebound(i_val);
        end
    end
    l_enough = cellfun(@numel, vals) > 5;
    avg_rebound = cellfun(@nanmean, vals);
    plot(edges(l_enough), avg_rebound(l_enough), '-', 'color', grp_clr, 'linewidth', 4)
    group_iv_vals{i_group} = vals;
end
xlabel('Membrane Potential At Sag Peak (mV)')
ylabel('Sag Rebound (mV)')

% ----- Regression for sag rebound --------
hf = figure;
hold on,
for i_group = 1:numel(groupdata.Ih_asym)
    hva_name = plotgroups{i_group, 3};
    grp_clr = hvaPlotColor(hva_name);
    [X,Y] = deal([]);
    for i_val = 1:numel(edges);
        tmpvals = group_iv_vals{i_group}{i_val};
        Y = cat(1,Y, tmpvals(:));
        X = cat(1,X, ones(numel(tmpvals), 1).*edges(i_val));
    end
    [B,BINT] = regress(Y(:), [ones(size(X(:))), X(:)]);
    [top_int, bot_int, x_mod, y_mod] = regression_line_ci(0.05, B, X, Y);
    CI_up = top_int - y_mod;
    CI_down = y_mod - bot_int;
    plot(X, Y, '.', 'color', grp_clr)
    shadedErrorBar(x_mod, y_mod, [CI_up ; CI_down], {'color', grp_clr});
end
xlabel('Membrane Potential At Sag Peak (mV)')
ylabel('Sag Rebound (mV)')


%% SUMMARY OF SPIKES: THRESHOLDS, RESET, F-I CURVES


% trial cutoff seems busted for lateral areas, would have thought there
% would be roughly equal trials for med and lat
%
% WHY IS THERE SUCH A LARGE DIFFERENCE IN NEURON COUNTS B/W THE SPIKE
% FREQUENCY STUFF AND THE THRESHOLD STUFF. I WOULD THINK THAT N FOR HZ AND 
% THRESHOLD WOULD BE THE SAME (NOT FOR RESET, WHICH REQUIRES 2 SPIKES...)
%
% MAYBE DEFINE N's FOR FIRST SPIKE ONLY (INSTEAD OF FOR FIRST AND LAST 
% SEPARATELY)
%
% Need to compute spike threshold on a spike by spike basis, or compute AHP
% between spike threshold crossings and not dvdt crossings..


close all; clc

trl_cutoff_num = 5;

% define a set of attributes for each analysis group
% {CellType, Layer,  BrainArea,  OpsinType}
% Brain Area can be: 'AL', 'PM', 'AM', 'LM', 'any', 'med', 'lat'. CASE SENSITIVE
plotgroups = {
    'PY', 'L23', 'med', 'any';...
    'PY', 'L23', 'lat', 'any';...
    };

groupdata = [];
groupdata.Vthresh = repmat({[]}, 1, size(plotgroups, 1)); % should only have N cells, where N = size(plotgroups, 1).
groupdata.Vreset = repmat({[]}, 1, size(plotgroups, 1));
groupdata.fi = repmat({[]}, 1, size(plotgroups, 1));

Ngroups = (numel(groupdata.Vthresh));
for i_ex = 1:numel(dat)
    
    for i_ch = 1:2
        
        % check to make sure this neuron was defined
        hs_name = sprintf('dcsteps_hs%d', i_ch);
        isvalid = ~strncmp(dat{i_ex}.info.fid.(hs_name)(end-8:end), [filesep, 'none.abf'], 9);
        if ~isvalid
            continue
        end
        
        % check the attributes
        ch_attribs = {dat{i_ex}.info.cellType{i_ch}, dat{i_ex}.info.brainArea, dat{i_ex}.info.opsin};
        group_idx = groupMatcher(plotgroups, ch_attribs);
        if sum(group_idx) == 0; continue; end
        
        % add the statistics to the appropriate cell array
        groupdata.Vthresh{group_idx} = cat(1, groupdata.Vthresh{group_idx}, dat{i_ex}.dcsteps.Vthresh(i_ch));
        groupdata.Vreset{group_idx} = cat(1, groupdata.Vreset{group_idx}, dat{i_ex}.dcsteps.Vreset(i_ch));
        groupdata.fi{group_idx} = cat(1, groupdata.fi{group_idx}, dat{i_ex}.dcsteps.fi_curve(i_ch));
    end
end


% fix the data so that the pA is on a consistent latice
%pA_cmd_bins = [-10;0;25;50;100;200;250;300;400;500;600;700;860];
pA_cmd_bins = [-20:25:875];
for i_group = 1:Ngroups
    
    tmp_vthresh_data = groupdata.Vthresh{i_group};
    tmp_vreset_data = groupdata.Vreset{i_group};
    tmp_fi_data = groupdata.fi{i_group};
    
    
    % reorganize the data [Nneurons, NpA]. One for first, one for last
    N_neurons = size(groupdata.Vthresh{i_group},1);
    [out_vthresh_first, out_vthresh_last] = deal(nan(N_neurons, numel(pA_cmd_bins)));
    [out_vreset_first, out_vreset_last] = deal(nan(N_neurons, numel(pA_cmd_bins)));
    out_fi_data = nan(N_neurons, numel(pA_cmd_bins));
    for i_neuron = 1:N_neurons
        
        % make sure each used the same pA_cmd for all measurements
        pA_thresh = tmp_vthresh_data{i_neuron}(:,1);
        pA_reset = tmp_vreset_data{i_neuron}(:,1);
        pA_fi = tmp_fi_data{i_neuron}(:,1);
        assert(all(pA_thresh==pA_reset))
        assert(all(pA_thresh==pA_fi))
        
        % bin the experimental commands into the ones specified. Only
        % consider positive current injections.
        [~,~,bin_idx]=histcounts(pA_thresh, pA_cmd_bins);
        for i_pA = 1:numel(pA_cmd_bins);
            
            
            l_in_bin = bin_idx == i_pA;
            
            out_vthresh_first(i_neuron, i_pA) = nanmean(tmp_vthresh_data{i_neuron}(l_in_bin, 2));
            out_vthresh_last(i_neuron, i_pA) = nanmean(tmp_vthresh_data{i_neuron}(l_in_bin, 3));
            out_vreset_first(i_neuron, i_pA) = nanmean(tmp_vreset_data{i_neuron}(l_in_bin, 2));
            out_vreset_last(i_neuron, i_pA) = nanmean(tmp_vreset_data{i_neuron}(l_in_bin, 3));
            out_fi_data(i_neuron, i_pA) = nanmean(tmp_fi_data{i_neuron}(l_in_bin, 2));
        end
        
        % here is where I could interpolate the fi_data, and set the values
        % above the max pA to nans. could also cull experiments that lacked
        % a high pA condition (e.g, >600pA), but this might unnecessarily
        % cull the IN recordings)
        old_dat = out_fi_data(i_neuron,:);
        l_nan = isnan(old_dat);
        old_dat = old_dat(~l_nan);
        old_pa = pA_cmd_bins(~l_nan);
        interp_dat = interp1(old_pa, old_dat, pA_cmd_bins, 'linear');
        out_fi_data(i_neuron, :) = interp_dat;
        
    end
    
    groupdata.Vthresh_ordered_first{i_group} = out_vthresh_first;
    groupdata.Vthresh_ordered_last{i_group} = out_vthresh_last;
    groupdata.Vreset_ordered_first{i_group} = out_vreset_first;
    groupdata.Vreset_ordered_last{i_group} = out_vreset_last;
    groupdata.fi_ordered{i_group} = out_fi_data;
    
end



% first the spike thresholds
figure, hold on,
legtext = {};
leghands = [];
for i_group = 1:Ngroups
    
    hva_name = plotgroups{i_group, 3};
    grp_clr = hvaPlotColor(hva_name);
    
    % first spike
    tmp_dat = groupdata.Vthresh_ordered_first{i_group};
    N = sum(~isnan(tmp_dat),1);
    
    l_N_enough = N>= trl_cutoff_num;
    xbar = nanmean(tmp_dat(:,l_N_enough), 1);
    sem = stderr(tmp_dat(:,l_N_enough), 1);
    plt_pA = pA_cmd_bins(l_N_enough);
    
    hp = shadedErrorBar(plt_pA, xbar, sem, {'o-', 'color', grp_clr});
    leghands(i_group) = hp.mainLine;
    legtext{i_group} =  sprintf('first sipke: %s, %s, %s', plotgroups{i_group,1},plotgroups{i_group, 2},plotgroups{i_group,3} );
    
    % last spike
    tmp_dat = groupdata.Vthresh_ordered_last{i_group};
    N = sum(~isnan(tmp_dat),1);
    
    l_N_enough = N>=trl_cutoff_num;
    xbar = nanmean(tmp_dat(:,l_N_enough), 1);
    sem = stderr(tmp_dat(:,l_N_enough), 1);
    plt_pA = pA_cmd_bins(l_N_enough);
    
    shadedErrorBar(plt_pA, xbar, sem, {'s--', 'color', grp_clr});
    
end
xlabel('Current injection (pA)')
ylabel('Avg Spike Threshold (mV)')
legend(leghands, legtext, 'location', 'northwest')
legend boxoff

% second, AHP reset Vm
figure, hold on,
legtext = {};
leghands = [];
for i_group = 1:Ngroups
    
    hva_name = plotgroups{i_group, 3};
    grp_clr = hvaPlotColor(hva_name);
    
    % first spike
    tmp_dat = groupdata.Vreset_ordered_first{i_group};
    N = sum(~isnan(tmp_dat),1);
    
    l_N_enough = N>= trl_cutoff_num;
    xbar = nanmean(tmp_dat(:,l_N_enough), 1);
    sem = stderr(tmp_dat(:,l_N_enough), 1);
    plt_pA = pA_cmd_bins(l_N_enough);
    
    hp = shadedErrorBar(plt_pA, xbar, sem, {'o-', 'color', grp_clr});
    leghands(i_group) = hp.mainLine;
    legtext{i_group} =  sprintf('first sipke: %s, %s, %s', plotgroups{i_group,1},plotgroups{i_group, 2},plotgroups{i_group,3} );
   
    
    % last spike
    tmp_dat = groupdata.Vreset_ordered_last{i_group};
    N = sum(~isnan(tmp_dat),1);
    
    l_N_enough = N>=trl_cutoff_num;
    xbar = nanmean(tmp_dat(:,l_N_enough), 1);
    sem = stderr(tmp_dat(:,l_N_enough), 1);
    plt_pA = pA_cmd_bins(l_N_enough);
    
    shadedErrorBar(plt_pA, xbar, sem, {'s--', 'color', grp_clr});
    
end
xlabel('Current injection (pA)')
ylabel('Avg AHP voltage (mV)')
legend(leghands, legtext, 'location', 'northwest')
legend boxoff

%% F-I CURVES FOR CURRENT INJECTIONS
close all; clc

% define a set of attributes for each analysis group
% {CellType, Layer,  BrainArea,  OpsinType}
% Brain Area can be: 'AL', 'PM', 'AM', 'LM', 'any', 'med', 'lat'. CASE SENSITIVE
plotgroups = {
    'PY', 'L23', 'med', 'any';...
    'PY', 'L23', 'lat', 'any';...
    };

groupdata = [];
groupdata.Vthresh = repmat({[]}, 1, size(plotgroups, 1)); % should only have N cells, where N = size(plotgroups, 1).
groupdata.Vreset = repmat({[]}, 1, size(plotgroups, 1));
groupdata.fi = repmat({[]}, 1, size(plotgroups, 1));

Ngroups = (numel(groupdata.Vthresh));
for i_ex = 1:numel(dat)
    
    for i_ch = 1:2
        
        % check to make sure this neuron was defined
        hs_name = sprintf('dcsteps_hs%d', i_ch);
        isvalid = ~strncmp(dat{i_ex}.info.fid.(hs_name)(end-8:end), [filesep, 'none.abf'], 9);
        if ~isvalid
            continue
        end
        
        % check the attributes
        ch_attribs = {dat{i_ex}.info.cellType{i_ch}, dat{i_ex}.info.brainArea, dat{i_ex}.info.opsin};
        group_idx = groupMatcher(plotgroups, ch_attribs);
        if sum(group_idx) == 0; continue; end
        
        % add the statistics to the appropriate cell array
        groupdata.Vthresh{group_idx} = cat(1, groupdata.Vthresh{group_idx}, dat{i_ex}.dcsteps.Vthresh(i_ch));
        groupdata.Vreset{group_idx} = cat(1, groupdata.Vreset{group_idx}, dat{i_ex}.dcsteps.Vreset(i_ch));
        groupdata.fi{group_idx} = cat(1, groupdata.fi{group_idx}, dat{i_ex}.dcsteps.fi_curve(i_ch));
    end
end

% -------  plot the spike rates as scatter plot  -----------
figure, hold on,
legtext = {};
leghands = [];
group_all_frs = {};
group_all_pa = {};
for i_group = 1:Ngroups
    tmp_dat = groupdata.fi{i_group};
    frs = [];
    pa = [];
    for i_ex = 1:numel(tmp_dat)
        l_spk = tmp_dat{i_ex}(:,1)>=0;
        pa = cat(1, pa, tmp_dat{i_ex}(l_spk,1));
        frs = cat(1, frs, tmp_dat{i_ex}(l_spk,2));
    end
    group_all_pa{i_group} = pa;
    group_all_frs{i_group} = frs;
end

for i_group = 1:Ngroups
    hva_name = plotgroups{i_group, 3};
    grp_clr = hvaPlotColor(hva_name);
    
    X = group_all_pa{i_group};
    Y = group_all_frs{i_group};
    
    [B,BINT] = regress(Y(:), [ones(size(X(:))), X(:)]);
    [top_int, bot_int, x_mod, y_mod] = regression_line_ci(0.05, B, X, Y);
    CI_up = top_int - y_mod;
    CI_down = y_mod - bot_int;
    
    hp = plot(X, Y, '.', 'color', grp_clr);
    shadedErrorBar(x_mod, y_mod, [CI_up ; CI_down], {'color', grp_clr});
    
    leghands(i_group) = hp;
    legtext{i_group} =  sprintf('%s, %s, %s', plotgroups{i_group,1},plotgroups{i_group, 2},plotgroups{i_group,3} );
    
    
end
xlabel('Current injection (pA)')
ylabel('Avg Spike Rate (Hz)')
legend(leghands, legtext, 'location', 'northwest')
legend boxoff


% -------  bin the data and plot the average  -----------
figure
edges = 0:150:900;
for i_group = 1:Ngroups
    
    hva_name = plotgroups{i_group, 3};
    grp_clr = hvaPlotColor(hva_name);
    
    X = group_all_pa{i_group};
    Y = group_all_frs{i_group};
    [counts, bin_idx] = histc(X, edges);
    bined_fi_curve_avg = nan(1, numel(edges));
    bined_fi_curve_sem = nan(1, numel(edges));
    for i_bin = 1:numel(edges)
        bined_fi_curve_avg(i_bin) = mean(Y(bin_idx==i_bin));
        bined_fi_curve_sem(i_bin) = stderr(Y(bin_idx==i_bin));
    end
    l_nans = isnan(bined_fi_curve_avg);
    l_enough_data = counts'>1;
    l_good = l_enough_data & ~l_nans;
    bined_fi_curve_avg = bined_fi_curve_avg(l_good);
    bined_fi_curve_sem = bined_fi_curve_sem(l_good);
    
    hse = shadedErrorBar(edges(l_good), bined_fi_curve_avg, [bined_fi_curve_sem ; bined_fi_curve_sem], {'color', grp_clr});
    leghands(i_group) = hse.mainLine;
    legtext{i_group} =  sprintf('%s, %s, %s', plotgroups{i_group,1},plotgroups{i_group, 2},plotgroups{i_group,3} );
end
xlabel('Current injection (pA)')
ylabel('Avg Spike Rate (Hz)')
hl = legend(leghands, legtext, 'location', 'northwest');
legend boxoff
hl.Interpreter = 'none';
set(gca, 'fontsize', 14)

%% SCATTER PLOT OF DIFFERENT TAU ESTIMATES
close all;
figure, hold on,
for i_ex = 1:numel(dat)
    if ~isfield(dat{i_ex}.dcsteps, 'tau')
        dat{i_ex}.info.fid
        continue
    end
    tau = dat{i_ex}.dcsteps.tau;
    if numel(tau)>=1 && ~isempty(tau{1})
        plot(tau{1}(:,3), tau{1}(:,4), 'ow', 'markerfacecolor', 'b')
    end
    if numel(tau)>=2 && ~isempty(tau{2})
        plot(tau{2}(:,3), tau{2}(:,4), 'ow', 'markerfacecolor', 'b')
    end
end
xlim([0, 0.07])
ylim([0, 0.07])
plot([0, 0.07], [0, 0.07], 'k')
set(gca, 'fontsize', 16)
xlabel('tau from exponential fit')
ylabel('tau from quick estimate')
axis square


%% PLOT THE RAW VCLAMP WAVEFORMS FOLLOWING EACH PULSE

close all


PLOT_ALL_TRIALS = false;
DEBUG_MEAN = true;
DEBUG_ALL = false;
NORM_TO_SMOOTH_P1 = false;
PLOT_RIT = false;


% goal: to plot the snipets, one after another, on a separate subplot for
% each stimulus type. All data will be plotted (if all trials are present),
% along with the mean

for i_ex = 1:numel(dat)
    
    f = figure;
    f.Units = 'Normalized';
    f.Position = [0.1 0.01 0.8 0.9];
    f.Name = sprintf('%s, site: %s, opsin: %s', dat{i_ex}.info.mouseName, dat{i_ex}.info.siteNum, dat{i_ex}.info.opsin);
    
    for i_ch = 1:2
        
        % skip past recording channels that have no data
        isvalid = dat{i_ex}.info.HS_is_valid_Vclamp(i_ch);
        if ~isvalid
            continue
        end
        
        % skip non-py cells
        if ~strcmpi(dat{i_ex}.info.cellType{i_ch}, 'PY_L23');
            continue
        end
        
        % skip non-ochief
        if ~regexpi(dat{i_ex}.info.opsin, 'chief'); continue; end
        
        
        conds = fieldnames(dat{i_ex}.expt);
        if PLOT_RIT
           plotList = 1:numel(conds); 
        else
           plotList = find(~strncmpi(conds, 'ritv', 4));
        end
        
        N_subplots = numel(plotList);
        for i_cond = 1:N_subplots
            idx = plotList(i_cond);
            
            snips = dat{i_ex}.expt.(conds{idx}).raw.snips{i_ch};
            if isempty(snips)
                continue
            end
            
            % transform the array for easy plotting. Add some nans to
            % separate the pulses.
            samps_per_pulse = size(snips, 2);
            Ntrials = size(snips, 3);
            Npulses = size(snips, 1);
            Nnans = round(samps_per_pulse./12);
            snips = cat(2, snips, nan(Npulses, Nnans, Ntrials));
            new_snip_length = size(snips, 2);
            snips = permute(snips, [2,1,3]);
            snips = reshape(snips, [], Ntrials);
            
            if NORM_TO_SMOOTH_P1
                trl_nums = dat{i_ex}.expt.(conds{idx}).realTrialNum{i_ch};
                p1amps = dat{i_ex}.qc.p1amp_norm{i_ch}(trl_nums);
                snips = bsxfun(@rdivide, snips, p1amps);
            end
            
            
            % plot the raw wave forms and the average waveform
            pltNum = sub2ind([2, N_subplots], i_ch, i_cond);
            hs = subplot(N_subplots, 2, pltNum); hold on,
            if PLOT_ALL_TRIALS
                plot(snips, '-', 'linewidth', 0.25)
            end
            plot(nanmean(snips,2), 'k-', 'linewidth', 2)
            hs.XTick = [];
            axis tight
            
            % add a point for the peak value
            peakVals = -1 .* dat{i_ex}.expt.(conds{idx}).stats.EPSCamp{i_ch};
            assert(~isempty(peakVals), 'ERROR: could not locate EPSCs')
            peakVals = permute(peakVals, [1,3,2]); % Npulses x Nsweeps
            
            if NORM_TO_SMOOTH_P1
                peakVals = bsxfun(@rdivide, peakVals, p1amps);
            end
            
            peakTimes = dat{i_ex}.expt.(conds{idx}).stats.latency{i_ch};
            peakTimes = permute(peakTimes, [1,3,2]); % Npulses x Nsweeps
            peakTimes_idx = round(peakTimes .* dat{i_ex}.info.sampRate.vclamp);
            presamps = round(dat{i_ex}.info.pretime.vclamp .* dat{i_ex}.info.sampRate.vclamp);
            offset = presamps + (((1:Npulses)-1) .* new_snip_length);
            peakTimes_idx = bsxfun(@plus, peakTimes_idx, offset(:));
            if DEBUG_ALL
                plot(peakTimes_idx, peakVals, 'co', 'markerfacecolor', 'c')
            end
            if DEBUG_MEAN
                plot(mean(peakTimes_idx, 2), mean(peakVals, 2), 'ro', 'markerfacecolor', 'r')
            end
            
            
        end
        
        drawnow
        
    end
    
    % close figures that are blank
    if isempty(f.Children); close(f); end
    
end
    

%% PLOT THE RAW ICLAMP WAVEFORMS FOLLOWING EACH PULSE

close all


PLOT_ALL_TRIALS = false;
DEBUG_MEAN = true;
DEBUG_ALL = false;
PLOT_RIT = false;


% goal: to plot the snipets, one after another, on a separate subplot for
% each stimulus type. All data will be plotted (if all trials are present),
% along with the mean

for i_ex = 1:numel(dat)
    
    % skip experiments with no WCSTP data
    if ~isfield(dat{i_ex}, 'iclamp')
        continue
    end
    
    f = figure;
    f.Units = 'Normalized';
    f.Position = [0.1 0.01 0.8 0.9];
    f.Name = sprintf('%s, site: %s, opsin: %s', dat{i_ex}.info.mouseName, dat{i_ex}.info.siteNum, dat{i_ex}.info.opsin);
    
    for i_ch = 1:2
        
        % skip past recording channels that have no data
        isvalid = dat{i_ex}.info.HS_is_valid_Iclamp(i_ch);
        if ~isvalid
            continue
        end
        
        
        conds = fieldnames(dat{i_ex}.iclamp);
        if PLOT_RIT
           plotList = 1:numel(conds); 
        else
           plotList = find(~strncmpi(conds, 'ritv', 4));
        end
        
        N_subplots = numel(plotList);
        for i_cond = 1:N_subplots
            idx = plotList(i_cond);
            
            % pull out the snips and add back the baseline
            snips = dat{i_ex}.iclamp.(conds{idx}).raw.snips{i_ch};
            %baselines = dat{i_ex}.iclamp.(conds{idx}).raw.bkgndVm{i_ch};
            %snips = bsxfun(@plus, snips, baselines);
            
            if isempty(snips)
                continue
            end
            
            % transform the array for easy plotting. Add some nans to
            % separate the pulses.
            samps_per_pulse = size(snips, 2);
            Ntrials = size(snips, 3);
            Npulses = size(snips, 1);
            Nnans = round(samps_per_pulse./12);
            snips = cat(2, snips, nan(Npulses, Nnans, Ntrials));
            new_snip_length = size(snips, 2);
            snips = permute(snips, [2,1,3]);
            snips = reshape(snips, [], Ntrials);
            

            % plot the raw wave forms and the average waveform
            pltNum = sub2ind([2, N_subplots], i_ch, i_cond);
            hs = subplot(N_subplots, 2, pltNum); hold on,
            if PLOT_ALL_TRIALS
                plot(snips, '-', 'linewidth', 0.25)
            end
            plot(nanmean(snips,2), 'k-', 'linewidth', 2)
            hs.XTick = [];
            axis tight
            
            % add a point for the peak value
            peakVals = dat{i_ex}.iclamp.(conds{idx}).stats.EPSPamp{i_ch};
            peakVals = permute(peakVals, [1,3,2]); % Npulses x Nsweeps
            
            
            peakTimes = dat{i_ex}.iclamp.(conds{idx}).stats.latency{i_ch};
            peakTimes = permute(peakTimes, [1,3,2]); % Npulses x Nsweeps
            peakTimes_idx = round(peakTimes .* dat{i_ex}.info.sampRate.iclamp);
            presamps = round(dat{i_ex}.info.pretime.iclamp .* dat{i_ex}.info.sampRate.iclamp);
            offset = presamps + (((1:Npulses)-1) .* new_snip_length);
            peakTimes_idx = bsxfun(@plus, peakTimes_idx, offset(:));
            if DEBUG_ALL
                plot(peakTimes_idx, peakVals, 'co', 'markerfacecolor', 'c')
            end
            if DEBUG_MEAN
                plot(mean(peakTimes_idx, 2), mean(peakVals, 2), 'ro', 'markerfacecolor', 'r')
            end
            
            
        end
        
        drawnow
        
    end
    
    % close figures that are blank
    if isempty(f.Children); close(f); end
    
end
    

%% COMPARISON OF FIRST PULSE RESPONSES ACROSS OPSINS

NORM_VALS = true;
PLOT_RAW = false;
CELL_TYPE = 'PY_L23';

p1pop.opsin = {};
p1pop.avgwf_vclamp = [];
p1pop.power = [];
p1pop.ND_status = {};
p1pop.p1amp = [];
for i_ex = 1:numel(dat)
    for i_ch = 1:2
        
        % check the cell type
        correct_cell_type = strcmpi(CELL_TYPE, dat{i_ex}.info.cellType{i_ch});
        if ~correct_cell_type; continue; end
        
        isvalid = dat{i_ex}.info.HS_is_valid_Vclamp(i_ch);
        if isvalid
           
            condnames = fieldnames(dat{i_ex}.expt);
            tmp_wfs = [];
            for i_cond = 1:numel(condnames)
                first_pulse = dat{i_ex}.expt.(condnames{i_cond}).raw.snips{i_ch}(1,:,:);
                if isempty(first_pulse); continue; end
                first_pulse = permute(first_pulse, [3,2,1]); % Npulses x Ntime
                
                if NORM_VALS
                    normvals = dat{i_ex}.expt.(condnames{i_cond}).stats.EPSCamp{i_ch};
                    if isempty(normvals); continue; end
                    normvals = permute(normvals(1,1,:), [3,2,1]);
                    first_pulse = bsxfun(@rdivide, first_pulse, normvals);
                end
                
                tmp_wfs = cat(1, tmp_wfs, first_pulse);
            end
            
            p1pop.opsin{end+1} = dat{i_ex}.info.opsin;
            p1pop.avgwf_vclamp = cat(1, p1pop.avgwf_vclamp, mean(tmp_wfs, 1));
            p1pop.tt = ([0:size(tmp_wfs,2)-1] ./ dat{i_ex}.info.sampRate.vclamp) - dat{i_ex}.info.pretime.vclamp;
            p1pop.power(end+1) = dat{i_ex}.expt.(condnames{1}).tdict(1);
            p1pop.ND_status{end+1} = dat{i_ex}.info.laser_power_supply_ND_filt;
            p1pop.p1amp(end+1) = nanmean(dat{i_ex}.qc.p1amp{i_ch}, 3);
            
        end
        
    end    
end

% -------- P1 amp waveforms ---------% 

f = figure; hold on,
clrs.chronos = 'b';
clrs.chief = 'r';
for opsin = {'chronos', 'chief'}
    l_opsin = cellfun(@(x) ~isempty(x), regexpi(p1pop.opsin, opsin));
    wf_avg = mean(p1pop.avgwf_vclamp(l_opsin,:), 1);
    wf_sem = stderr(p1pop.avgwf_vclamp(l_opsin, :), 1);
    if PLOT_RAW
        plot(p1pop.tt.*1000, p1pop.avgwf_vclamp(l_opsin,:)', '-', 'color',  clrs.(opsin{1}))
    end
    shadedErrorBar(p1pop.tt.*1000, wf_avg, wf_sem, {'linewidth', 2, 'color', clrs.(opsin{1})});
end
xlabel('time (ms)')
legend('chronos', 'chief')
if NORM_VALS
    ylabel('norm amplitude')
else
    ylabel('average waveform (pA)')
end


% -------- P1 amp vs. Laser power ---------% 

% adjust the nominal laser power to account for a lack of an ND filter
l_10 = cellfun(@(x) ~isempty(regexpi(x, '10')), p1pop.ND_status);
l_ND_filter = cellfun(@(x) isempty(regexpi(x, 'none')), p1pop.ND_status);
pow_scalar = ones(numel(p1pop.ND_status), 1);
pow_scalar(~l_ND_filter) = 5;
l_to_plot = l_ND_filter;


f = figure; hold on,
clrs.chronos = 'b';
clrs.chief = 'r';
for opsin = {'chronos', 'chief'}
    l_opsin = cellfun(@(x) ~isempty(x), regexpi(p1pop.opsin, opsin));
    l_tmp = l_opsin;
    tmp_pow = p1pop.power(l_tmp) .* pow_scalar(l_tmp)';
    tmp_amps = p1pop.p1amp(l_tmp);
    plot(tmp_pow, tmp_amps, 'o', 'markerfacecolor', clrs.(opsin{1}), 'markeredgecolor', clrs.(opsin{1}))
end
xlabel('Laser cmd voltage')
ylabel('EPSC amplitude')


%% PPRs VS. DEPTH FOR EPSCs ONTO PY CELLS

clc, close all

ppr_vs_depth.pnp1 = [];
ppr_vs_depth.depth = [];
ppr_vs_depth.brainarea = {};
ppr_vs_depth.opsin = {};
for i_ex = 1:numel(dat)
    for i_ch = 1:2
        % is there Vclamp data,
        isvalid = dat{i_ex}.info.HS_is_valid_Vclamp(i_ch);
        if ~isvalid; continue; end
        
        % is this a PY cell?
        is_py_cell = strcmpi(dat{i_ex}.info.cellType{i_ch}, 'py_l23');
        if ~is_py_cell; continue; end
        
        % for 25Hz only: get the average P1 through P10 amplitudes
        dat_25hz = [];
        condnames = fieldnames(dat{i_ex}.expt);
        for i_cond = 1:numel(condnames)
            ex_tf = dat{i_ex}.expt.(condnames{i_cond}).tdict(3);
            tf_is_valid =  ex_tf == 25;
            is_not_rit = dat{i_ex}.expt.(condnames{i_cond}).tdict(5) == 0;
            if tf_is_valid && is_not_rit
                tmp_dat = permute(dat{i_ex}.expt.(condnames{i_cond}).stats.EPSCamp{i_ch}, [3,1,2]);
                real_trl_nums = dat{i_ex}.expt.(condnames{i_cond}).realTrialNum{i_ch};
                smooth_p1_amps = dat{i_ex}.qc.p1amp_norm{i_ch}(real_trl_nums);
                tmp_dat(:,1) = smooth_p1_amps(:);
                dat_25hz = cat(1, dat_25hz, tmp_dat(:,1:10));
            end
        end
        
        if ~isempty(dat_25hz)
            pprs = bsxfun(@rdivide, dat_25hz, dat_25hz(:,1));
            avg_ppr = mean(pprs, 1);
            
            ppr_vs_depth.pnp1 = cat(1, ppr_vs_depth.pnp1, avg_ppr);
            ppr_vs_depth.depth = cat(1,  ppr_vs_depth.depth, dat{i_ex}.info.cellDepth_um(i_ch));
            ppr_vs_depth.brainarea = cat(1, ppr_vs_depth.brainarea, dat{i_ex}.info.brainArea);
            ppr_vs_depth.opsin = cat(1, ppr_vs_depth.opsin, dat{i_ex}.info.opsin);
        end
        
    end
end

l_opsin = cellfun(@(x) ~isempty(x), regexpi(ppr_vs_depth.opsin, 'chief'));
l_has_depth = ~isnan(ppr_vs_depth.depth);
l_med = cellfun(@(x) ~isempty(x), regexpi(ppr_vs_depth.brainarea, 'am|pm'));
l_lat = cellfun(@(x) ~isempty(x), regexpi(ppr_vs_depth.brainarea, 'al|lm'));

p2p1 = ppr_vs_depth.pnp1(:,2);
p4p1 = ppr_vs_depth.pnp1(:,4);
p10p1 = ppr_vs_depth.pnp1(:,10);

ymax = 410;
figure, hold on
set(gcf, 'position', [440         459        1084         416])
subplot(2,3,1), hold on,
l_plt = l_med & l_opsin & l_has_depth;
X = ppr_vs_depth.depth(l_plt);
Y = p2p1(l_plt);
plot(X, Y, 'bo')
[B,~] = regress(Y(:), [ones(size(X(:))), X(:)]);
[top_int, bot_int, x_mod, y_mod] = regression_line_ci(0.05, B, X, Y);
CI_up = top_int - y_mod;
CI_down = y_mod - bot_int;
shadedErrorBar(x_mod, y_mod, [CI_up ; CI_down], {'color', 'b'});
xlim([150, ymax])
[rho, p] = corr(X, Y, 'type', 'pearson')
text(155, max(Y), sprintf('r = %.2f, p = %.2f', rho, p))
title('P2:P1', 'fontsize', 14)
ylabel('Med Areas', 'fontsize', 14)

subplot(2,3,2), hold on,
Y = p4p1(l_plt);
plot(X, Y, 'bo')
[B,~] = regress(Y(:), [ones(size(X(:))), X(:)]);
[top_int, bot_int, x_mod, y_mod] = regression_line_ci(0.05, B, X, Y);
CI_up = top_int - y_mod;
CI_down = y_mod - bot_int;
shadedErrorBar(x_mod, y_mod, [CI_up ; CI_down], {'color', 'b'});
xlim([150, ymax])
[rho, p] = corr(X, Y, 'type', 'pearson');
text(155, max(Y), sprintf('r = %.2f, p = %.2f', rho, p))
title('P4:P1', 'fontsize', 14)

subplot(2,3,3), hold on,
Y = p10p1(l_plt);
plot(X, Y, 'bo')
[B,~] = regress(Y(:), [ones(size(X(:))), X(:)]);
[top_int, bot_int, x_mod, y_mod] = regression_line_ci(0.05, B, X, Y);
CI_up = top_int - y_mod;
CI_down = y_mod - bot_int;
shadedErrorBar(x_mod, y_mod, [CI_up ; CI_down], {'color', 'b'});
xlim([150, ymax])
[rho, p] = corr(X, Y, 'type', 'pearson');
text(155, max(Y), sprintf('r = %.2f, p = %.2f', rho, p))
title('P10:P1', 'fontsize', 14)

subplot(2,3,4), hold on,
l_plt = l_lat & l_opsin & l_has_depth;
X = ppr_vs_depth.depth(l_plt);
Y = p2p1(l_plt);
plot(X, Y, 'ro')
[B,~] = regress(Y(:), [ones(size(X(:))), X(:)]);
[top_int, bot_int, x_mod, y_mod] = regression_line_ci(0.05, B, X, Y);
CI_up = top_int - y_mod;
CI_down = y_mod - bot_int;
shadedErrorBar(x_mod, y_mod, [CI_up ; CI_down], {'color', 'r'});
xlim([150, ymax])
[rho, p] = corr(X, Y, 'type', 'pearson');
text(155, max(Y), sprintf('r = %.2f, p = %.2f', rho, p))
xlabel('depth (um)', 'fontsize', 14)
ylabel('Lat Areas', 'fontsize', 14)

subplot(2,3,5), hold on,
Y = p4p1(l_plt);
plot(X, Y, 'ro')
[B,~] = regress(Y(:), [ones(size(X(:))), X(:)]);
[top_int, bot_int, x_mod, y_mod] = regression_line_ci(0.05, B, X, Y);
CI_up = top_int - y_mod;
CI_down = y_mod - bot_int;
shadedErrorBar(x_mod, y_mod, [CI_up ; CI_down], {'color', 'r'});
xlim([150, ymax])
[rho, p] = corr(X, Y, 'type', 'pearson');
text(155, max(Y), sprintf('r = %.2f, p = %.2f', rho, p))
xlabel('depth (um)', 'fontsize', 14)

subplot(2,3,6), hold on,
Y = p10p1(l_plt);
plot(X, Y, 'ro')
[B,~] = regress(Y(:), [ones(size(X(:))), X(:)]);
[top_int, bot_int, x_mod, y_mod] = regression_line_ci(0.05, B, X, Y);
CI_up = top_int - y_mod;
CI_down = y_mod - bot_int;
shadedErrorBar(x_mod, y_mod, [CI_up ; CI_down], {'color', 'r'});
xlim([150, ymax])
[rho, p] = corr(X, Y, 'type', 'pearson');
text(155, max(Y), sprintf('r = %.2f, p = %.2f', rho, p))
xlabel('depth (um)', 'fontsize', 14)

% plot all data aggregated together
figure, hold on,
l_plt = l_opsin & l_has_depth;
X = ppr_vs_depth.depth(l_plt);
Y = p2p1(l_plt);
plot(X, Y, 'ko')
[B,~] = regress(Y(:), [ones(size(X(:))), X(:)]);
[top_int, bot_int, x_mod, y_mod] = regression_line_ci(0.05, B, X, Y);
CI_up = top_int - y_mod;
CI_down = y_mod - bot_int;
shadedErrorBar(x_mod, y_mod, [CI_up ; CI_down], {'color', 'k'});
xlim([150, ymax])
[rho, p] = corr(X, Y, 'type', 'pearson');
text(155, max(Y), sprintf('r = %.2f, p = %.2f', rho, p))
xlabel('depth (um)', 'fontsize', 14)
ylabel('All Cells', 'fontsize', 14)
%% STP SUMMARY FOR EACH RECORDING


NORMAMPS = true;
PLOTRAW = false;

for i_ex = 1:numel(dat)
    
    if ~isfield(dat{i_ex}, 'expt')
        continue
    end
    
    f = figure;
    f.Name = sprintf('Mouse %s, site %s, opsin %s', dat{i_ex}.info.mouseName, dat{i_ex}.info.siteNum, dat{i_ex}.info.opsin);
    f.Units = 'normalized';
    f.Position = [0 0    1    1];
    
    for i_ch = 1:2
        
        isvalid = dat{i_ex}.info.HS_is_valid_Vclamp(i_ch);
        if ~isvalid
            continue
        end
        
        conds = fieldnames(dat{i_ex}.expt);
        Nconds = numel(conds);
        for i_cond = 1:Nconds
            if isempty(dat{i_ex}.expt.(conds{i_cond}).stats.EPSCamp{i_ch})
                continue
            end
            
            subplot(Nconds,2, 2.*(i_cond-1) + i_ch);
            
            tmp = dat{i_ex}.expt.(conds{i_cond}).stats.EPSCamp{i_ch};
            
            if NORMAMPS
                p1 = squeeze(dat{i_ex}.qc.p1amp{i_ch});
                normfact = squeeze(dat{i_ex}.qc.p1amp_norm{i_ch});
                realTrlNums = dat{i_ex}.expt.(conds{i_cond}).realTrialNum{i_ch};
                normfact = normfact(realTrlNums);
                normfact = permute(normfact, [3,1,2]);
                tmp = bsxfun(@rdivide, tmp, normfact);
            end
            
            tt = dat{i_ex}.expt.(conds{i_cond}).pOnTimes;
            if PLOTRAW
                plot(tt, permute(tmp, [1,3,2]), '-')
            else
                xbar = mean(tmp,3);
                sem = stderr(tmp,3);
                my_errorbar(tt(1:end-1), xbar(1:end-1), sem(1:end-1), '-ok', 'markersize', 2, 'linewidth', 2);
            end
            
            xlim([0, dat{i_ex}.info.sweepLength.vclamp]);
            if i_cond < Nconds
                set(gca, 'xticklabel', []);
            end
             
            if ~NORMAMPS
                set(gca, 'yticklabel', []);
            end
            
        end
    end
    drawnow
    
end

%% STP: LOOKING FOR CLASS 1b AND 2 INPUTS

NORM_TO_RUNNING_AVG = true;
PULSE_NUM = 4; % PnP1

% initialize the figure
f = figure; hold on,

for i_ex = 1:numel(dat)
    
    if ~isfield(dat{i_ex}, 'expt')
        continue
    end
    
    
    for i_ch = 1:2
        
        isvalid = dat{i_ex}.info.HS_is_valid_Vclamp(i_ch);
        if ~isvalid
            continue
        end
        
        conds = fieldnames(dat{i_ex}.expt);
        Nconds = numel(conds);
        pnp1 = [];
        isi = [];
        for i_cond = 1:Nconds
            
            is_rit = strcmp(conds{i_cond}(1:4), 'RITv');
            no_data = isempty(dat{i_ex}.expt.(conds{i_cond}).stats.EPSCamp{i_ch});
            is_not_py_cell = ~strcmpi(dat{i_ex}.info.cellType{i_ch}(1:3), 'py_');
            if is_rit || no_data || is_not_py_cell
                continue
            end
            
            % grab the data
            tmp = dat{i_ex}.expt.(conds{i_cond}).stats.EPSCamp{i_ch};
            
            if NORM_TO_RUNNING_AVG
                % deal with non-stationarities
                normfact = squeeze(dat{i_ex}.qc.p1amp_norm{i_ch});
                realTrlNums = dat{i_ex}.expt.(conds{i_cond}).realTrialNum{i_ch};
                normfact = normfact(realTrlNums);
                normfact = permute(normfact, [3,1,2]);
                tmp = bsxfun(@rdivide, tmp, normfact);
            end
            
            % always re-normalize to the first pulse
            normfact = tmp(1,:,:);
            tmp = bsxfun(@rdivide, tmp, normfact);
            
            
            % store the PNP1 data
            tmp = mean(tmp, 3);
            if numel(tmp)<PULSE_NUM; continue; end
            pnp1 = cat(1, pnp1, tmp(PULSE_NUM));
            
            % store the isi
            isi = cat(1, isi, round(1000 ./ dat{i_ex}.expt.(conds{i_cond}).tdict(3)));
        end
        
        % sometimes there will be multiple versions of the same isi.
        % average across these
        unique_isi = unique(isi);
        plot_dat = [];
        plot_isi = [];
        for i_cond = 1:numel(unique_isi)
            plot_dat = cat(1, plot_dat, mean(pnp1(isi == unique_isi(i_cond))));
            plot_isi = cat(1, plot_isi, unique_isi(i_cond));
        end
        
        % plot the data
        switch dat{i_ex}.info.opsin
            case {'chief_flx', 'chief_cit', 'chief'}
                pltclr = 'r';
            case {'chronos', 'chronos_flx', 'chronos_gfp', 'chronos_tom', 'chronos_9_gfp'}
                pltclr = 'g';
            otherwise
                error('unknown opsin type')
        end
        
        figure(f)
        plot(plot_isi, plot_dat, '-', 'color', pltclr)
        xlabel('ISI (ms)')
        ylabel(sprintf('P%d:Pn ratio', PULSE_NUM))
        set(gca, 'xscale', 'log', 'yscale', 'log')
        axis tight
        
        
    end
    
end


%% I-CLAMP POPULATION ANALYSIS (DATA COLLECTION)

% loop through the experiments. Pull out the trains data. Ignore the
% recovery train (if present) and aggregate across recovery conditions.
iclamp_pop = [];
iclamp_pop.TFsAllExpts = [];
iclamp_pop.MaxNPulses = 0;
for i_ex = 1:numel(dat)
    
    if ~isfield(dat{i_ex}, 'iclamp');
        continue
    end
    
    % find the normal trains. Assume the field name is NOT 'ritv'
    condnames = fieldnames(dat{i_ex}.iclamp);
    l_trains = ~strncmp(condnames, 'RITv', 4);
    if sum(l_trains)==0; continue; end % no trains data
    trainParams = cellfun(@(x) dat{i_ex}.iclamp.(condnames{x}).tdict, mat2cell(find(l_trains), ones(sum(l_trains),1), 1), 'uniformoutput', false);
    trainParams = cat(1, trainParams{:});
    
    % make sure the pulse amplitude and width were identical across ttypes
    assert(numel(unique(trainParams(:,1)))==1, 'ERROR: more than one pulse amplitude')
    assert(numel(unique(trainParams(:,2)))==1, 'ERROR: more than one pulse width');
    
    % identify the unique TF conditions for this experiment, and update the
    % running log of TFs used across all experiments
    uniqueTFs = unique(trainParams(:,3));
    tmp = cat(1, iclamp_pop.TFsAllExpts, uniqueTFs);
    iclamp_pop.TFsAllExpts = unique(tmp);
    
    % store the stim params for all stim types, which will be useful for
    % indexing later.
    allStimParams = cellfun(@(x) dat{i_ex}.iclamp.(condnames{x}).tdict, mat2cell((1:numel(l_trains))', ones(numel(l_trains),1), 1), 'uniformoutput', false);
    allStimParams = cat(1, allStimParams{:});
    
    % store some metadata
    iclamp_pop.info{i_ex} = dat{i_ex}.info;
    
    
    % aggregate data within TF conditons (Separately for each
    % recording channel
    for i_tf = 1:numel(uniqueTFs);
        
        tfidx = find(trainParams(:,3) == uniqueTFs(i_tf)); % condnames that contain the train with a particular TF
        
        for i_ch = 1:2
            
            % check to make sure this neuron was defined
            isvalid = dat{i_ex}.info.HS_is_valid_Iclamp(i_ch);
            if ~isvalid
                iclamp_pop.dat{i_ex}.xbar{i_ch}{i_tf} = [];
                iclamp_pop.tfs{i_ex}{i_ch} = [];
                continue
            end
            
            % iterate over the trains with the same freq. 
            catdat = [];
            for i_cond = 1:numel(tfidx)
                
                % pull out the data
                condIdx = ismember(allStimParams, trainParams(tfidx(i_cond),:), 'rows');
                assert(sum(condIdx)==1, 'ERROR, found zero or more than 1 instance of the trial type')
                tmpdat = dat{i_ex}.iclamp.(condnames{condIdx}).raw.snips{i_ch}; %[Npulses x Ntime x Nsweeps]
                tmpdat = mean(tmpdat,3); % mean across sweeps
                if ~isempty(tmpdat)
                    % normalize by the p1Amp
                    p1_amps = squeeze(dat{i_ex}.iclamp.(condnames{condIdx}).stats.EPSPamp{i_ch}(1,:,:));
                    normfact = mean(p1_amps);
                    tmpdat = tmpdat ./ normfact;
                    
                    % store in a matrix for averaging later
                    catdat = cat(3, catdat, tmpdat);
                end
            end
            
            % check to make sure there are data for these conditions. Even
            % though this recording channel should be defined (See above),
            % it's possible there are no data due to deletion of single
            % sweeps
            if isempty(catdat)
                iclamp_pop.dat{i_ex}.xbar{i_tf}{i_ch} = [];
                iclamp_pop.tfs{i_ex}{i_ch} = [];
            else
                
                % now average across trails, and re-normalize to the first
                % pulse. Store in the population data structure.
                avg = mean(catdat,3); % now [Npulses, Ntime]
                
                % delete the last pulse (if it's a recovery pulse, but
                % store the recovery pulse in a different field of the
                % population structure
                isrecovery = trainParams(tfidx(i_cond),4) > 0;
                if isrecovery
                    avg(end,:) = [];
                end
                
                iclamp_pop.dat{i_ex}.xbar{i_tf}{i_ch} = avg;
                iclamp_pop.tfs{i_ex}{i_ch} = uniqueTFs;
                iclamp_pop.MaxNPulses = max([iclamp_pop.MaxNPulses, size(avg,1)]);

            end
        end
    end
end


%% I-CLAMP POPULATION ANALYSIS (PLOTTING)

clc; close all

PLOTERRBAR = true;

% define a set of attributes for each analysis group
% {CellType, Layer,  BrainArea,  OpsinType}
% Brain Area can be: 'AL', 'PM', 'AM', 'LM', 'any', 'med', 'lat'. CASE SENSITIVE
plotgroups = {
    'PY', 'L23', 'LM', 'chief';...
    'PY', 'L23', 'AL', 'chief';...
    'PY', 'L23', 'PM', 'chief';...
    'PY', 'L23', 'AM', 'chief';...
    };

% initalize the population structure
allTFs = iclamp_pop.TFsAllExpts;
for i_grp = 1:size(plotgroups, 1)
   groupdata_raw{i_grp} =  repmat({[]}, numel(allTFs), 1);
end

% iterate over the experiments. For each recording channel, determine what
% the attributes are, and place the data in the correct ploting group.
for i_ex = 1:numel(dat)
    
    for i_ch = 1:2
        
        % check to make sure this neuron was defined
        isvalid = dat{i_ex}.info.HS_is_valid_Iclamp(i_ch);
        if ~isvalid
            continue
        end
        
        % check to make sure this neuron has iClamp population data, which
        % it might not if it was only RITs or something like that.
        if isempty(iclamp_pop.dat{i_ex})
            continue
        end
        
        % check the attributes
        ch_attribs = {dat{i_ex}.info.cellType{i_ch}, upper(dat{i_ex}.info.brainArea), dat{i_ex}.info.opsin}; % force the brain area to be uppercase
        group_idx = groupMatcher(plotgroups, ch_attribs);
        if sum(group_idx) == 0; continue; end
        
        
        % add data to the appropriate group data array. This will probably
        % bonk if the sampling rate or number of pulses is different across experiments
        ch_tfs = iclamp_pop.tfs{i_ex}{i_ch};
        for i_tf = 1:numel(ch_tfs)
            tf_idx = allTFs == ch_tfs(i_tf);
            tmpdat = iclamp_pop.dat{i_ex}.xbar{i_tf}{i_ch}; % [Npulses x Ntime]
            
            % add some nans for separating the pulses during plotting.
            whitesamps = round(size(tmpdat,2) ./ 17);
            tmpdat = cat(2, tmpdat, nan(size(tmpdat,1), whitesamps));
            tmpdat = reshape(tmpdat', [], 1)'; % notice the extra transpose to make this a row vector
            
            
            % add to the population structure
            if ~isempty(tmpdat)
                groupdata_raw{group_idx}{tf_idx} = cat(1, groupdata_raw{group_idx}{tf_idx}, tmpdat);
            end
        end

    end
end


% plot the data, all line series
f = figure;
f.Units = 'normalized';
f.Position = [0.1, 0.01, 0.3, 0.9];
Ntfs = numel(allTFs);
for i_tf = 1:Ntfs
    subplot(Ntfs, 1, i_tf), hold on,
    for i_grp = 1:numel(groupdata_raw)
        tmp = groupdata_raw{i_grp}{i_tf};
        xbar = nanmean(tmp, 1);
        N = sum(~isnan(tmp),1);
        unique(N)
        sem = nanstd(tmp, [], 1) ./ sqrt(N);
        plt_clr = hvaPlotColor(plotgroups{i_grp, 3});
        if PLOTERRBAR
            shadedErrorBar(1:size(tmp,2), xbar, sem, {'color', plt_clr, 'linewidth', 2});
        else
            plot(tmp', '-', 'color', plt_clr)
            %plot(mean(tmp,1), '-', 'color', groupcolors{i_grp}, 'linewidth', 3)
        end
    end
    ylim([0,2])
    ylabel('EPSP amplitude')
end




%% V-CLAMP POPULATION ANALYSIS (DATA COLLECTION)

% loop through the experiments. Pull out the trains data. Ignore the
% recovery train (if present) and aggregate across recovery conditions.
MIN_TRL_COUNT = 0; % set to 0 to turn off
MIN_TFS_COUNT = 0; % set to 0 to turn off
MIN_RECOV_COUNT = 0; % set to 0 to turn off
recovpop = [];
recovpop.TFsAllExpts = [];
recovpop.recoveryTimesAllExpts = [];
for i_ex = 1:numel(dat)
    
    % find the normal trains. Assume the field name is NOT 'ritv'
    condnames = fieldnames(dat{i_ex}.expt);
    l_trains = ~strncmp(condnames, 'RITv', 4);
    if sum(l_trains)==0; continue; end % no trains data
    
    condnames_trains = condnames(l_trains);
    trainParams = cellfun(@(x) dat{i_ex}.expt.(x).tdict, condnames_trains, 'uniformoutput', false);
    trainParams = cat(1, trainParams{:});
    recovpop.trainParams{i_ex} = trainParams;
    unique_tfs = unique(trainParams(:,3));
    unique_recov = unique(trainParams(:,4));
    unique_recov(unique_recov == 0) = [];
    if isempty(unique_recov); continue; end
    
    % make sure the pulse amplitude and width were identical across ttypes
    assert(numel(unique(trainParams(:,1)))==1, 'ERROR: more than one pulse amplitude')
    assert(numel(unique(trainParams(:,2)))==1, 'ERROR: more than one pulse width');
    
    % aggregate data across TF conditons and recovery times (Separately for each
    % recording channel)
    recovpop.dat{i_ex}.psc_amps = {};
    recovpop.dat{i_ex}.psc_wfs = {};
    recovpop.dat{i_ex}.ignore={[],[]};
    Ntfs = numel(unique_tfs);
    Nrecov = numel(unique_recov);
    for i_ch = 1:2
        
        % remove instances where there is insufficient data to constrain
        % the fit
        insufficient_data = size(dat{i_ex}.qc.p1amp{i_ch}, 3) < MIN_TRL_COUNT;
        insufficient_tfs = Ntfs < MIN_TFS_COUNT;
        insufficient_recovs = Nrecov < MIN_RECOV_COUNT;
        if insufficient_data || insufficient_tfs || insufficient_recovs
            recovpop.dat{i_ex}.ignore{i_ch} = true;
        end
        
        
        %flag instances of potential bad datafiles that should be excluded
        %later
        if strcmpi(dat{i_ex}.info.cellType{i_ch}, 'PY_L23')
            p1amp_avg = nanmean(dat{i_ex}.qc.p1amp{i_ch});
%             if p1amp_avg < 200
%                 recovpop.dat{i_ex}.ignore{i_ch} = true;
%                 warning('culling data on basis of P1 amp')
%             end
        end
        
        % initialize the outputs. I need to concatenate similar TFs
        % together, but keep each of their recovery conditions separate
        %
        % recovpop.dat.psc_amps = {{TF}, {psc_amps}, {{recov_time}, {recov_psc_amp}}}
        % recovpop.dat.psc_wfs = {{TF}, {psc_wfs}, {{recov_time}, {recov_psc_wfs}}}
        recovpop.dat{i_ex}.psc_amps{i_ch} = {num2cell(unique_tfs), repmat({[]}, Ntfs,1), repmat({num2cell(unique_recov), repmat({[]}, Nrecov,1)}, Ntfs,1)};
        recovpop.dat{i_ex}.psc_wfs{i_ch} =  {num2cell(unique_tfs), repmat({[]}, Ntfs,1), repmat({num2cell(unique_recov), repmat({[]}, Nrecov,1)}, Ntfs,1)};

        for i_tf = 1:Ntfs
            
            for i_recov = 1:Nrecov
                
                cond_idx = (trainParams(:,3)==unique_tfs(i_tf)) & (trainParams(:,4)==unique_recov(i_recov));
                
                % check to make sure this neuron was defined, and that the
                % stimulus train is a is a recovery train. also check to make
                % sure there are data for these conditions. Even though this
                % recording channel should be defined (See above), it's
                % possible there are no data due to deletion of single sweeps
                isvalid = dat{i_ex}.info.HS_is_valid_Vclamp(i_ch);
                isrecovery = trainParams(cond_idx,4) > 0;
                no_data = isempty(dat{i_ex}.expt.(condnames_trains{cond_idx}).stats.EPSCamp{i_ch});
                if ~isvalid || ~isrecovery || no_data
                    continue
                end
                
                % pull out the EPSC amplitudes (mean across sweeps), norm to
                % first pulse
                tmp_psc = mean(dat{i_ex}.expt.(condnames_trains{cond_idx}).stats.EPSCamp{i_ch}, 3);
                norm_fact = tmp_psc(1);
                tmp_psc = tmp_psc ./ norm_fact;
                
                % concatenate 1:10 into 'trains', and the 11th into 'recov'
                recovpop.dat{i_ex}.psc_amps{i_ch}{2}{i_tf} = cat(1, recovpop.dat{i_ex}.psc_amps{i_ch}{2}{i_tf}, tmp_psc(1:end-1)');
                recovpop.dat{i_ex}.psc_amps{i_ch}{3}{i_tf,2}{i_recov} = cat(1, recovpop.dat{i_ex}.psc_amps{i_ch}{3}{i_tf,2}{i_recov}, tmp_psc(end));

                % pull out the EPSC waveforms (mean across sweeps)
                tmp_wfs = mean(dat{i_ex}.expt.(condnames_trains{cond_idx}).raw.snips{i_ch}, 3);% Npulses x Ntime
                norm_fact = min(tmp_wfs(1,:), [], 2);
                tmp_wfs = tmp_wfs ./ -norm_fact; 
                
                % store the waveforms in the population structures
                recovpop.dat{i_ex}.psc_wfs{i_ch}{2}{i_tf} = cat(3, recovpop.dat{i_ex}.psc_wfs{i_ch}{2}{i_tf}, tmp_wfs(1:end-1, :));
                recovpop.dat{i_ex}.psc_wfs{i_ch}{3}{i_tf,2}{i_recov} = cat(3, recovpop.dat{i_ex}.psc_wfs{i_ch}{3}{i_tf,2}{i_recov}, tmp_wfs(end,:));
                
            end
        end
    end
    
    
    
    % identify the unique TF conditions for this experiment, and update the
    % running log of TFs used across all experiments. This will be used by
    % the plotting routines to figure out the "canonical grid"
    tmp = cat(1, recovpop.TFsAllExpts, unique_tfs);
    recovpop.TFsAllExpts = unique(tmp);
    
    % identify the unique recovery times for this experiment, and update the
    % running log of recov times used across all experiments. This will be used by
    % the plotting routines to figure out the "canonical grid"
    tmp = cat(1, recovpop.recoveryTimesAllExpts, unique_recov);
    recovpop.recoveryTimesAllExpts = unique(tmp);

    
end

%% V-CLAMP POPULATION ANALYSIS (PLOTS)
% clc; close all
% plot recovery pulse amplitude (normalized) vs. train frequency. Do this
% separately for cell types, brain areas, opsins, etc...

PLOT_RAW_DATA = false;
PLOT_AVG_DATA = true;
PLOT_MANIFOLD_OF_AVG_RAW_DATA = false;
PLOT_AVG_MANIFOLD = true;

% define a set of attributes for each analysis group
% {CellType, Layer,  BrainArea,  OpsinType}
% Brain Area can be: 'AL', 'PM', 'AM', 'LM', 'any', 'med', 'lat'. CASE SENSITIVE
man_plotgrps = {
    'PY', 'L23', 'med', 'chief';...
    'PY', 'L23', 'lat', 'chief';...
    };

allTFs = recovpop.TFsAllExpts;
Ntfs = numel(allTFs);
allRecoveryTimes = recovpop.recoveryTimesAllExpts';
Nrecovs = numel(allRecoveryTimes);

% set up the output arrays
%
%  {{trains_(tf=X)}{recov=A, tf=X}{recov=B, tf=X}...} one cell for each recovery pulse, one row for each TF
%
template = repmat({[]}, Ntfs, Nrecovs+1);
groupdata = [];
groupdata.wfs = repmat({template}, size(man_plotgrps,1), 1);
groupdata.amps = repmat({template}, size(man_plotgrps,1), 1);

empty_array = repmat({[]}, 1, size(man_plotgrps, 1)); % should only have N cells, where N = size(man_plotgrps, 1). Each cell has a matrix with a cononicalGrid:
[groupdata_smooth, group_params] = deal(empty_array);

% iterate over the experiments. For each recording channel, determine what
% the attributes are, and place the data in the correct ploting group.
for i_ex = 1:numel(dat)
    
    % skip experiments that have no data (i.e., trains but no recovery)
    if isempty(recovpop.dat{i_ex}); continue; end
    
    for i_ch = 1:2
        % check to make sure this neuron was defined
        isvalid = dat{i_ex}.info.HS_is_valid_Vclamp(i_ch);
        if ~isvalid
            continue
        end
        
        % Ignore the data (could be due to small EPSCs or insufficient
        % data)
        if recovpop.dat{i_ex}.ignore{i_ch}
            continue
        end
        
        % check the attributes
        ch_attribs = {dat{i_ex}.info.cellType{i_ch}, dat{i_ex}.info.brainArea, dat{i_ex}.info.opsin};
        group_idx = groupMatcher(man_plotgrps, ch_attribs);
        if sum(group_idx) == 0; continue; end
        
        % pull out the waveform data
        %
        % Need to order things by ALL TFS and RECOVs, and not w/r/t an
        % individual data file's lattice
        %
        wf_top_level = recovpop.dat{i_ex}.psc_wfs{i_ch};
        amps_top_level = recovpop.dat{i_ex}.psc_amps{i_ch};
        ex_tfs = cell2mat(wf_top_level{1});
        
        wfs_train = cellfun(@(x) mean(x, 3), wf_top_level{2}, 'uniformoutput', false); % average across replicates
        amps_train = cellfun(@(x) mean(x, 1), amps_top_level{2}, 'uniformoutput', false); % average across replicates
        for i_tf = 1:numel(ex_tfs)
            
            % grab the data from cell array containing average across sweeps
            tf_wf = wfs_train{i_tf};
            tf_wf = cat(2, tf_wf, nan(size(tf_wf, 1), 25)); % add some nans for separation
            tf_wf = reshape(tf_wf', 1, []); % a row vector
            tf_amps = amps_train{i_tf};
            
            
            % push the wfs data into the correct slot in the groupdata
            row_idx = allTFs == ex_tfs(i_tf);
            tmp_group = groupdata.wfs{group_idx}{row_idx, 1};
            tmp_group = cat(1, tmp_group, tf_wf); % Ntraces x Ntime
            groupdata.wfs{group_idx}{row_idx, 1} = tmp_group;
            
            % push the amps data into the correct slot in the groupdata
            tmp_group = groupdata.amps{group_idx}{row_idx, 1};
            tmp_group = cat(1, tmp_group, tf_amps); % Nneurons x Npulses
            groupdata.amps{group_idx}{row_idx, 1} = tmp_group;
            
            % loop over the recovery conditions and deal with those
            ex_recovs = cell2mat(wf_top_level{3}{i_tf, 1});
            wfs_recovs = cellfun(@(x) mean(x, 3), wf_top_level{3}{i_tf,2}, 'uniformoutput', false); % average across replicates
            amps_recovs = cellfun(@(x) mean(x, 3), amps_top_level{3}{i_tf, 2}, 'uniformoutput', false); % average across replicates. Each recovTime gets a scalar
            for i_recov = 1:numel(ex_recovs)
                col_idx = [false, allRecoveryTimes == ex_recovs(i_recov)]; % leading False to offset the "trains" position
                tmp_group = groupdata.wfs{group_idx}{row_idx, col_idx};
                tmp_group = cat(1, tmp_group, wfs_recovs{i_recov});
                groupdata.wfs{group_idx}{row_idx, col_idx} = tmp_group;
                
                % deal with amps
                tmp_group = groupdata.amps{group_idx}{row_idx, col_idx};
                tmp_group = cat(1, tmp_group, amps_recovs{i_recov});
                groupdata.amps{group_idx}{row_idx, col_idx} = tmp_group;
            end
        end
        if PLOT_AVG_MANIFOLD
            groupdata_smooth{group_idx} = cat(3, groupdata_smooth{group_idx}, pprpop.smoothManifold{i_ex}{i_ch});
            group_params{group_idx} = cat(1, group_params{group_idx}, pprpop.params{i_ex}{i_ch});
        end
    end
end



%
% Plots of the waveforms
%
%%%%%%%%%%%%
f = figure;
f.Units = 'Normalized';
f.Position = [0.0865, 0, 0.3161, 1];
f.Color = 'w';
Ngroups = size(man_plotgrps,1);
legtext = repmat({{}}, 1, Ntfs);
leghand = repmat({[]}, 1, Ntfs);
for i_group = 1:size(man_plotgrps, 1)
    plt_clr = hvaPlotColor(man_plotgrps{i_group, 3});
    
    for i_tf = 1:Ntfs
        
        subplot(Ntfs, 1, i_tf), hold on,
        
        % check to make sure there are data for this TF condition
        data_exist = ~isempty(groupdata.wfs{i_group}{i_tf,1});
        if data_exist
            
            
            % concatenate the time series, padding with nans
            N_expts = cellfun(@(x) size(x,1), groupdata.wfs{i_group}(i_tf,:));
            N_expts = max(N_expts);
            all_wfs = cellfun(@(x) cat(1, x, nan(N_expts-size(x,1), size(x,2))), groupdata.wfs{i_group}(i_tf,:), 'uniformoutput', false);
            all_wfs = cellfun(@(x) cat(2, x, nan(N_expts, 200)), all_wfs, 'uniformoutput', false);
            all_wfs = cat(2, all_wfs{:});
            
            % determine the mean across experiments
            xbar = nanmean(all_wfs, 1);
            sem = stderr(all_wfs, 1);
            
            % mean +/- sem
            if PLOT_AVG_DATA
                hp = shadedErrorBar(1:numel(xbar), xbar, sem, {'color', plt_clr});
                leghand{i_tf}(end+1) = hp.mainLine;
                legtext{i_tf}{end+1} = sprintf('%s, %s, %s, N=%d', man_plotgrps{i_group,1}, man_plotgrps{i_group,2}, man_plotgrps{i_group,3}, N_expts);
            end
            
            
            % all the data
            if PLOT_RAW_DATA
                plot(1:numel(xbar), all_wfs', '-', 'color', plt_clr)
            end
            
        end
        
    end
end
% indicate the TF, and recovery times, add a legend
for i_tf = 1:Ntfs
    subplot(Ntfs, 1, i_tf)
    %ylim([-1.4, 0.2])
    plt_train_tf = allTFs(i_tf);
    l_recov_exist = cellfun(@(x) ~isempty(x), groupdata.wfs{i_group}(i_tf,:));
    l_recov_exist(1) = []; % ignore the first array b/c it's not recov pulse
    plt_recov_ms = allRecoveryTimes(l_recov_exist);
    plt_string = sprintf('Train = %d Hz, Recovery (in ms): %s', plt_train_tf, num2str(plt_recov_ms));
    txt_hand = text(750, 0.2, plt_string);
    txt_hand.FontSize = 9;
    legend(leghand{i_tf}, legtext{i_tf}, 'location', 'best', 'interpreter', 'none')
    legend boxoff
end



%
% Plots of the normalized amplitudes
%
%%%%%%%%%%%%
f = figure;
f.Units = 'Normalized';
f.Position = [0.0865, 0, 0.3161, 1];
f.Color = 'w';
Ngroups = size(man_plotgrps,1);
legtext = repmat({{}}, 1, Ntfs);
leghand = repmat({[]}, 1, Ntfs);
for i_group = 1:size(man_plotgrps, 1)
    plt_clr = hvaPlotColor(man_plotgrps{i_group, 3});
    
    for i_tf = 1:Ntfs
        
        hs = subplot(Ntfs, 1, i_tf); hold on,
        hs.TickDir = 'out';
        
        % check to make sure there are data for this TF condition
        data_exist = ~isempty(groupdata.amps{i_group}{i_tf,1});
        if data_exist
            
            
            % grab the data for this tf condition
            N_expts = cellfun(@(x) size(x,1), groupdata.amps{i_group}(i_tf,:));
            N_expts = max(N_expts);
            train_amps = groupdata.amps{i_group}{i_tf,1};
            recov_amps = cellfun(@(x) cat(1, x, nan(N_expts-size(x,1), size(x,2))), groupdata.amps{i_group}(i_tf,2:end), 'uniformoutput', false);
            recov_amps = cat(2, recov_amps{:});
            
            
            % determine the mean across experiments
            xbar_trains = nanmean(train_amps, 1);
            sem_trains = stderr(train_amps, 1);
            xx_trains = 1:numel(xbar_trains);
            xbar_recov = nanmean(recov_amps, 1);
            sem_recov = stderr(recov_amps, 1);
            xx_recov = numel(xbar_trains)+1 : numel(xbar_trains)+numel(xbar_recov);
            
            % plot
            hp = errorbar(xx_trains, xbar_trains, sem_trains, 'color', plt_clr);
            errorbar(xx_recov, xbar_recov, sem_recov, 'o', 'markerfacecolor', plt_clr, 'markeredgecolor', plt_clr)
            leghand{i_tf}(end+1) = hp(1);
            legtext{i_tf}{end+1} = sprintf('%s, %s, %s, N=%d', man_plotgrps{i_group,1}, man_plotgrps{i_group,2}, man_plotgrps{i_group,3}, N_expts);
            
            % add predictions from the fits to the averaged raw data
            if PLOT_MANIFOLD_OF_AVG_RAW_DATA
                % check that the plot_group defined here is the same as the
                % plot-group defined by the group_fit_to_grand_avg dataset
                data_matches = all(strcmpi(man_plotgrps(i_group,:), group_fit_to_grand_avg{i_group}.grouptypes));
                assert(data_matches, 'ERROR: plot groups do not match')
                
                % get the fit parameters and predict the amplitudes for
                % each TF/recov condition
                stp_params = group_fit_to_grand_avg{i_group}.params;
                model = group_fit_to_grand_avg{i_group}.model;
                [d, tau_d, f, tau_f] = parse_vca_model_coeffs(stp_params, model);
                A0 = 1;
                
                N_train_pulses = numel(xbar_trains);
                train_tf = allTFs(i_tf);
                l_recovs = cellfun(@(x) ~isempty(x), groupdata.amps{i_group}(i_tf,2:end));
                recov_times = allRecoveryTimes(l_recovs);
                
                % iterate over TF/recov conds and generate predictions
                % based on the model
                train_amps = [];
                recov_amps = [];
                for i_recov = 1:numel(recov_times)
                    ipi_sec = 1/train_tf;
                    ptimes = cumsum(ones(1, N_train_pulses)*ipi_sec) - ipi_sec;
                    recov_time_sec = recov_times(i_recov) ./ 1000;
                    ptimes(end+1) = ptimes(end)+recov_time_sec;
                    
                    pred_amps = predict_vca_psc(ptimes, d, tau_d, f, tau_f, A0);
                    train_amps(1,:) = pred_amps(1:N_train_pulses);
                    recov_amps(i_recov) = pred_amps(end);
                end
                plot([xx_trains, xx_recov], [train_amps, recov_amps], '--', 'color', plt_clr)
            end
            if PLOT_AVG_MANIFOLD
                train_isi_ms = 1000 ./ allTFs(i_tf);
                [~, smooth_tf_idx] = min(abs(pprpop.smoothManifold_isi_ms - train_isi_ms));
                N_pulses_smooth = pprpop.smoothManifold_numPulses;
                
                grid_average = nanmean(groupdata_smooth{i_group}, 3);
                h_smooth = plot(1:N_pulses_smooth, grid_average(:, smooth_tf_idx), '--', 'color', plt_clr);
            end
        end
        
    end
end
% indicate the TF, and recovery times, add a legend
%figure(f)
for i_tf = 1:Ntfs
    subplot(Ntfs, 1, i_tf)
    plt_train_tf = allTFs(i_tf);
    l_recov_exist = cellfun(@(x) ~isempty(x), groupdata.wfs{i_group}(i_tf,:));
    l_recov_exist(1) = []; % ignore the first array b/c it's not recov pulse
    plt_recov_ms = allRecoveryTimes(l_recov_exist);
    plt_string = sprintf('Train = %d Hz, Recovery (in ms): %s', plt_train_tf, num2str(plt_recov_ms));
    txt_hand = text(2, 1.65, plt_string);
    txt_hand.FontSize = 9;
    legend(leghand{i_tf}, legtext{i_tf}, 'location', 'best', 'interpreter', 'none')
    legend boxoff
end


%
% Plots expanded predictions of model
%
%%%%%%%%%%%%
if PLOT_MANIFOLD_OF_AVG_RAW_DATA
    f = figure;
    f.Units = 'Normalized';
    f.Position = [0.0865, 0, 0.3161, 1];
    f.Color = 'w';
    for i_group = 1:size(man_plotgrps, 1)
        plt_clr = hvaPlotColor(man_plotgrps{i_group, 3});
        for i_tf = 1:Ntfs
            
            hs = subplot(Ntfs, 1, i_tf); hold on,
            hs.TickDir = 'out';
            
            
            % check that the plot_group defined here is the same as the
            % plot-group defined by the group_fit_to_grand_avg dataset
            data_matches = all(strcmpi(man_plotgrps(i_group,:), group_fit_to_grand_avg{i_group}.grouptypes));
            assert(data_matches, 'ERROR: plot groups do not match')
            
            % get the fit parameters and predict the amplitudes for
            % each TF/recov condition
            stp_params = group_fit_to_grand_avg{i_group}.params;
            model = group_fit_to_grand_avg{i_group}.model;
            [d, tau_d, f, tau_f] = parse_vca_model_coeffs(stp_params, model);
            
            N_train_pulses = 10;
            train_tf = allTFs(i_tf);
            recov_times_ms = 100:100:15000;
            
            % iterate over TF/recov conds and generate predictions
            % based on the model
            train_amps = [];
            recov_amps = [];
            for i_recov = 1:numel(recov_times_ms)
                ipi_sec = 1/train_tf;
                ptimes = cumsum(ones(1, N_train_pulses)*ipi_sec) - ipi_sec;
                recov_time_sec = recov_times_ms(i_recov) ./ 1000;
                ptimes(end+1) = ptimes(end)+recov_time_sec;
                
                pred_amps = predict_vca_psc(ptimes, d, tau_d, f, tau_f, A0);
                train_amps(1,:) = pred_amps(1:N_train_pulses);
                recov_amps(i_recov) = pred_amps(end);
            end
            xx_trains = 1:numel(train_amps);
            plot(xx_trains, train_amps, '-', 'color', plt_clr)
            xx_recov = xx_trains(end) + 1 : xx_trains(end) + numel(recov_amps);
            plot(xx_recov, recov_amps, 'o',...
                                       'markeredgecolor', plt_clr,...
                                       'markerfacecolor', plt_clr,...
                                       'markersize', 4)
            t = title(sprintf('%d Hz: d=[%.1f,%.1f], dtau=[%.1f,%.1f], f=[%.1f,%.1f], ftau=[%.1f,%.1f]', train_tf, d, tau_d, f, tau_f));
            t.FontSize=10;
        end
    end
end


%% BUILD AN MULTI-FACTOR ANOVA ANALYSIS



anovaTFs = [12, 25, 50];
[all_pnum_labels, all_log_pprs, all_tf_labels, all_hva_labels] = deal([]);
for i_group = 1:size(man_plotgrps, 1)
    for i_tf = 1:numel(anovaTFs)
        tf_idx = recovpop.TFsAllExpts == anovaTFs(i_tf);
        
        % get the PPRs for the trains only (no recov conds)
        tmp_pprs = log(groupdata.amps{i_group}{tf_idx,1});
        
        % make a pulse number index vector
        tmp_pnum_labels = repmat(1:size(tmp_pprs,2), size(tmp_pprs, 1), 1);
        tmp_pnum_labels = tmp_pnum_labels(:);
        
        % make the other grouping variables
        tmp_tf_labels = ones(size(tmp_pnum_labels)) .* anovaTFs(i_tf);
        tmp_hva_labels = repmat(man_plotgrps{i_group, 3}, numel(tmp_pnum_labels), 1);
        
        % concatenate the tmp labels into the population vectors
        all_log_pprs = cat(1, all_log_pprs, tmp_pprs(:));
        all_pnum_labels = cat(1, all_pnum_labels, tmp_pnum_labels);
        all_tf_labels = cat(1, all_tf_labels, tmp_tf_labels);
        all_hva_labels = cat(1, all_hva_labels, tmp_hva_labels);
    end
end

[p, tbl, stats] = anovan(all_log_pprs,...
                        {all_hva_labels, all_tf_labels, all_pnum_labels},...
                        'model', 'interaction',...
                        'varnames', {'all_hva_labels', 'all_tf_labels', 'all_pnum_labels'});
multcompare(stats, 'Dimension', [3])



%%  WITHIN-MOUSE STP COMPARISON

clc; close all;


% iterate over mice and store the PPRs for each HVA tested
%
% 1) I can ignore mice that only have one line in the sheet
freq_tags = {'_f12_', '_f25_', '_f50_'};
wincomp_pop = struct();
for i_ex = 1:numel(dat)
    
    ex_mouse_name = dat{i_ex}.info.mouseName;
    ex_hva = dat{i_ex}.info.brainArea;
    

    % is this mouse in the pop array? If not, add it to the array
    is_in_array = isfield(wincomp_pop, ex_mouse_name);
    if ~is_in_array
        wincomp_pop.(ex_mouse_name) = struct('LM', {[]}, 'PM', {[]}, 'AL', {[]}, 'AM', {[]}, 'RL', {[]}, 'med', {[]}, 'lat', {[]});
    end
      
    % loop over cells and TF conditions. Keep the PPRs in an array that is 
    % Npulses x NTFs x Ncells
    for i_ch = 1:2
        
        % ignore some of the experiments
        is_valid = dat{i_ex}.info.HS_is_valid_Vclamp(i_ch);
        if ~is_valid; continue; end
        is_chief = any(regexpi(dat{i_ex}.info.opsin, 'chief'));
        is_pyl23 = any(regexpi(dat{i_ex}.info.cellType{i_ch}, 'PY_L23'));
        if ~is_chief || ~is_pyl23; continue; end
        
        % find the fields that have the correct frequencies and iterate
        % through them in order
        all_conds = fieldnames(dat{i_ex}.expt);
        cell_data = nan(10,3); % pulses x TF
        for i_tf = 1:numel(freq_tags)
            tf_recov_cond_inds = cellfun(@(x) ~isempty(x), regexpi(all_conds, freq_tags{i_tf}));
            tf_recov_conds = all_conds(tf_recov_cond_inds);
            tf_recov_data_sum = zeros(10, 1);
            n_sweeps = 0;
            for i_recov = 1:numel(tf_recov_conds)
                tmp_dat = dat{i_ex}.expt.(tf_recov_conds{i_recov}).stats.EPSCamp{i_ch};
                if isempty(tmp_dat); continue; end
                tmp_dat = tmp_dat(1:10,:,:); % remove the recovery pulse
                tf_recov_data_sum = tf_recov_data_sum + sum(tmp_dat, 3);
                n_sweeps = n_sweeps + size(tmp_dat, 3);
            end
            
            % compute the mean for this TF condition (across recov times).
            % Normalize to the first pulse
            tf_data_mean = tf_recov_data_sum ./ n_sweeps;
            pprs = tf_data_mean ./ tf_data_mean(1);
            
            % add to the cell_data matrix
            cell_data(:, i_tf) = pprs;
        end
        
        % add the cell_data matrix to the appropriate HVA field in the
        % population array
        tmp_array = wincomp_pop.(ex_mouse_name).(ex_hva);
        tmp_array = cat(3, tmp_array, cell_data);
        wincomp_pop.(ex_mouse_name).(ex_hva) = tmp_array;
        
        % add the cell_data matrix to either the "med" or "lat" arrays in
        % the populaiton database
        switch lower(ex_hva)
            case {'am', 'pm'}
                alt_hva = 'med';
            case {'al', 'lm', 'rl'}
                alt_hva = 'lat';
        end
        tmp_array = wincomp_pop.(ex_mouse_name).(alt_hva);
        tmp_array = cat(3, tmp_array, cell_data);
        wincomp_pop.(ex_mouse_name).(alt_hva) = tmp_array;
        
    end
    
end

% 3) make the full-factorized comparisons of the HVAs.
%hvas = {'LM', 'AL', 'PM', 'AM'};
hvas = {'med', 'lat'};
hva_pair_inds = combinator(numel(hvas), 2, 'c');
pulse_num_for_plot = 3;
hf_p1pn = figure;
hf_tfs = figure;
n_cols = size(hva_pair_inds, 1);
n_rows = numel(freq_tags);

for i_comp = 1:n_cols
    
    hva1 = hvas{hva_pair_inds(i_comp, 1)};
    hva2 = hvas{hva_pair_inds(i_comp, 2)};
    
    [hva1_dat, hva2_dat] = deal([]);
    mouse_names = fieldnames(wincomp_pop);
    for i_mouse = 1:numel(mouse_names)
        
        % do not analyze mice for which paired data aren't available
        has_hva1 = ~isempty(wincomp_pop.(mouse_names{i_mouse}).(hva1));
        has_hva2 = ~isempty(wincomp_pop.(mouse_names{i_mouse}).(hva2));
        if ~has_hva1 || ~has_hva2; continue; end
        
        % add this mouse's data to the population array. Average when there
        % are multiple cells for a single area with in a mouse
        hva1_dat = cat(3, hva1_dat, mean(wincomp_pop.(mouse_names{i_mouse}).(hva1), 3));
        hva2_dat = cat(3, hva2_dat, mean(wincomp_pop.(mouse_names{i_mouse}).(hva2), 3));    
    end
    
    
    % deal with the full plot of all PPRs in sequence
    figure(hf_tfs)
    for i_tf = 1:n_rows
        pltidx = sub2ind([n_cols, n_rows], i_comp, i_tf);
        subplot(n_rows, n_cols, pltidx);
        hold on,
        hva1_clr = hvaPlotColor(hva1);
        errorbar(1:10, mean(hva1_dat(:, i_tf, :), 3),...
                       stderr(hva1_dat(:, i_tf, :), 3),...
                       '-', 'color', hva1_clr);
                   
        hva2_clr = hvaPlotColor(hva2);
        errorbar(1:10, mean(hva2_dat(:, i_tf, :), 3),...
                       stderr(hva2_dat(:, i_tf, :), 3),...
                       '-', 'color', hva2_clr);
        
        if i_tf == 1
            leg_text = {sprintf('%s, n=%d', hva1, size(hva1_dat,3)),...
                        sprintf('%s, n=%d', hva2, size(hva2_dat,3))};
            legend(leg_text, 'location', 'best');
            title(sprintf('%s vs. %s', hva1, hva2));
        end
        if i_comp == 1
            ylabel(sprintf('Norm EPSC \n tf = %s Hz', freq_tags{i_tf}(3:4)));
        end
        xlim([0, 11])
    end
    
    % deal with the pn:p1 ratios
    figure(hf_p1pn)
    for i_tf = 1:n_rows
        pltidx = sub2ind([n_cols, n_rows], i_comp, i_tf);
        subplot(n_rows, n_cols, pltidx);
        hold on,
        plt_mtx = [permute(hva1_dat(pulse_num_for_plot, i_tf, :), [1,3,2]);...
                   permute(hva2_dat(pulse_num_for_plot, i_tf, :), [1,3,2])];
        plot([1,2], plt_mtx, 'k-');
        
        if i_tf == 1
            title(sprintf('%s vs. %s', hva1, hva2));
        end
        if i_comp == 1
            ylabel(sprintf('Norm EPSC \n tf = %s Hz', freq_tags{i_tf}(3:4)));
        end
        xlim([0.5, 2.5])
        axis tight
    end
    
end







%% ESTIMATE TIME CONSTANTS OF SHORT TERM PLASTICITY

MODEL = 'ddff'; % string (eg 'ddff') or number of d and f terms (eg, 3) for FULLFACT
TRAINSET = 'recovery';  % could be 'rit', 'recovery', 'all'
PLOT_TRAINING_DATA = true;
FIT_RECOVERY_PULSE = true;
DATA_FOR_FIT = 'avg'; % can be 'avg', 'norm', 'raw'. 
CONVERT_TO_SMOOTH_P1 = true;
FORCE_RECOVERY = false;

% write an anyonomous helper function to find the training data
clear isTrainingSet
switch TRAINSET
    case 'rit'
        isTrainingSet = @(x) strncmp(x, 'RITv', 4);
    case 'recovery'
        isTrainingSet = @(x) ~strncmp(x, 'RITv', 4);
    case 'all'
        isTrainingSet = @(x) true(size(x));
end


for i_ex = 1:numel(dat)   %  i_ex = 17,46 is a good one for debugging.
    clc
    hf = figure;
    hf.Units = 'Normalized';
    hf.Position = [0.3347    0.0422    0.6639    0.7489];
    hf.Name = sprintf('Mouse %s, site %s, opsin: %s.  Train with: %s, Plot training set: %d',...
        dat{i_ex}.info.mouseName, dat{i_ex}.info.siteNum, dat{i_ex}.info.opsin, TRAINSET, PLOT_TRAINING_DATA);
    chempty = false(1,2);
    for i_ch = 1:2
        
        % determine if there are data to fit
        isvalid = dat{i_ex}.info.HS_is_valid_Vclamp(i_ch);
        if ~isvalid
            chempty(i_ch) = true;
            if  all(chempty)
                close(hf)
            end
            continue
        end
        
        
        % make a mini dataset that's composed of only pOnTimes, EPSC amps,
        % and p1Amps
        [training_data, xval_data] = deal(struct('pOnTimes', {{}}, 'rawAmps', {{}}, 'tdict', {[]}, 'unique_conds_rawAmps', {{}}));
        condnames = fieldnames(dat{i_ex}.expt);
        for i_cond = 1:numel(condnames)
            
            % check to make sure there are data for this condition
            if isempty(dat{i_ex}.expt.(condnames{i_cond}).stats.EPSCamp{i_ch})
                continue
            end
            
            % grab pOnTimes, raw EPSC amps, and tdict
            pOnTimes = dat{i_ex}.expt.(condnames{i_cond}).pOnTimes;
            rawAmps = dat{i_ex}.expt.(condnames{i_cond}).stats.EPSCamp{i_ch};
            ex_tdict = dat{i_ex}.expt.(condnames{i_cond}).tdict;
            istrain = isTrainingSet(condnames{i_cond});
            
            if CONVERT_TO_SMOOTH_P1
                real_trl_nums = dat{i_ex}.expt.(condnames{i_cond}).realTrialNum{i_ch};
                smooth_p1_amps = dat{i_ex}.qc.p1amp_norm{i_ch}(real_trl_nums);
                rawAmps(1,1,:) = smooth_p1_amps;
                rawAmps = bsxfun(@rdivide, rawAmps, rawAmps(1,1,:));
            end
            
            % delete the recovery pulse if need be
            if ~FIT_RECOVERY_PULSE
                ex_tdict(4) = 0; % set the recovery field to zero so that down-stream analysis will not think this is a recovery pulse (and hack off an additional pulse)
                has_recov_pulse = dat{i_ex}.expt.(condnames{i_cond}).tdict(4) > 0;
                if has_recov_pulse
                    pOnTimes(end) = [];
                    rawAmps = rawAmps(1:end-1,1,:);
                end
            end
            
            % force the model to fit a ficticious data point at +30 seconds
            % at a normalized value of 1
            if FORCE_RECOVERY
                if ~FIT_RECOVERY_PULSE; error('ERROR: need to fit recov for force_recovery'); end
                if ~CONVERT_TO_SMOOTH_P1; error('ERROR: need to normalize for force_recovery'); end
                pOnTimes(end+1:end+3) = [30, 35, 40]; % in seconds
                nsweeps = size(rawAmps,3);
                rawAmps(end+1:end+3,:,:) = repmat([1;1;1], 1, 1, nsweeps);
            end
            
            % allocate the data to either the training or xval datasets
            if isTrainingSet(condnames{i_cond})
                training_data.pOnTimes = cat(1, training_data.pOnTimes, pOnTimes);
                training_data.rawAmps = cat(1, training_data.rawAmps, rawAmps);
                training_data.tdict = cat(1, training_data.tdict, ex_tdict);
            else
                xval_data.pOnTimes = cat(1, xval_data.pOnTimes, pOnTimes);
                xval_data.rawAmps = cat(1, xval_data.rawAmps, rawAmps);
                xval_data.tdict = cat(1, xval_data.tdict, ex_tdict);
            end
        end
        
        % determine if there are xval data
        xval_exists = ~isempty(xval_data.tdict);
        
        % Combine across unique conditions. This only makes a difference if
        % there were multiple different recovery conditions within a single
        % TF, AND the analysis is set up to ignore the recovery pulses...
        training_data.unique_conds_tdict = unique(training_data.tdict, 'rows');
        for i_cond = 1:size(training_data.unique_conds_tdict,1)
            cond_inds = ismember(training_data.tdict, training_data.unique_conds_tdict(i_cond,:), 'rows');
            training_data.unique_conds_rawAmps{i_cond} = cat(3, training_data.rawAmps{cond_inds}); % sweeps is 3rd dim
            training_data.unique_conds_pOnTimes{i_cond} = training_data.pOnTimes{find(cond_inds, 1, 'first')};
        end
        if xval_exists
            xval_data.unique_conds_tdict = unique(xval_data.tdict, 'rows');
            for i_cond = 1:size(xval_data.unique_conds_tdict,1)
                cond_inds = ismember(xval_data.tdict, xval_data.unique_conds_tdict(i_cond,:), 'rows');
                xval_data.unique_conds_rawAmps{i_cond} = cat(3, xval_data.rawAmps{cond_inds}); % sweeps is 3rd dim
                xval_data.unique_conds_pOnTimes{i_cond} = xval_data.pOnTimes{find(cond_inds, 1, 'first')};
            end
        end
        
        % average the data sets (possibly for fitting, but also for plotting)
        training_data.unique_conds_xbar = cellfun(@(x) mean(x,3), training_data.unique_conds_rawAmps, 'uniformoutput', false);
        if xval_exists
            xval_data.unique_conds_xbar = cellfun(@(x) mean(x,3), xval_data.unique_conds_rawAmps, 'uniformoutput', false);
        end
        
        % create a normalized version of the xbar arrays (don't normalize the raw data)
        training_data.unique_conds_xbar_norm = cellfun(@(x) x./x(1), training_data.unique_conds_xbar, 'uniformoutput', false);
        if xval_exists
            xval_data.unique_conds_xbar_norm = cellfun(@(x) x./x(1), xval_data.unique_conds_xbar, 'uniformoutput', false);
        end
        
        % now assign data to be fit
        switch lower(DATA_FOR_FIT)
            case 'avg'
                training_data.amps_for_fit = training_data.unique_conds_xbar;
                training_data.ptimes_for_fit = training_data.unique_conds_pOnTimes;
                if xval_exists
                    xval_data.amps_for_fit = xval_data.unique_conds_xbar;
                    xval_data.ptimes_for_fit = xval_data.unique_conds_pOnTimes;
                end
            case 'norm'
                training_data.amps_for_fit = training_data.unique_conds_xbar_norm;
                training_data.ptimes_for_fit = training_data.unique_conds_pOnTimes;
                if xval_exists
                    xval_data.amps_for_fit = xval_data.unique_conds_xbar_norm;
                    xval_data.ptimes_for_fit = xval_data.unique_conds_pOnTimes;
                end
            case 'raw'
                training_data.amps_for_fit = training_data.unique_conds_rawAmps;
                training_data.ptimes_for_fit = training_data.unique_conds_pOnTimes;
                if xval_exists
                    xval_data.amps_for_fit = xval_data.unique_conds_rawAmps;
                    xval_data.ptimes_for_fit = xval_data.unique_conds_pOnTimes;
                end
        end
        
        
        % if there were no training data, then move along,
        if isempty(training_data.amps_for_fit)
            chempty(i_ch) = true;
            if all(chempty)
                close(hf)
                continue
            end
        end
        
        
        % look for instances where there are pOnTimes but no data (could
        % happen if some sweeps get deleted from one HS but not the other.
        l_empty = cellfun(@isempty, training_data.amps_for_fit);
        training_data.amps_for_fit(l_empty) = [];
        training_data.ptimes_for_fit(l_empty) = [];
        
        fit_results = {};
        % fit the data
        if ischar(MODEL)
            [tmp_best_fit_params, fxn_val] = fit_vca_model(training_data.amps_for_fit, training_data.ptimes_for_fit, MODEL);
            fit_results{1}.params = tmp_best_fit_params;
            fit_results{1}.fxn_val = fxn_val;
            fit_results{1}.model = MODEL;
            
        elseif isscalar(MODEL)
            all_models = fullfact([MODEL+1, MODEL+1])-1; % -1 to have models with just d or f. +1 to bring the total back up to the desired number
            l_zeros = sum(all_models, 2)==0;
            all_models(l_zeros,:) = []; % delete models with zero d and zero f terms
            
            for i_mod = 1:size(all_models,1)
                d_terms = repmat('d',1,all_models(i_mod,1));
                f_terms = repmat('f',1,all_models(i_mod,2));
                tmp_model = strcat(d_terms, f_terms);
                
                [tmp_best_fit_params, fxn_val] = fit_vca_model(training_data.amps_for_fit, training_data.ptimes_for_fit, tmp_model);
                fit_results{i_mod}.params = tmp_best_fit_params;
                fit_results{i_mod}.fxn_val = fxn_val;
                fit_results{i_mod}.model = tmp_model;
            end
            
        else
            error('MODEL was neither a scalar or a string')
        end
        
        % predict all the data (for all types of models)
        for i_mod = 1:numel(fit_results)
            [d, dTau, f, fTau] = parse_vca_model_coeffs(fit_results{i_mod}.params, fit_results{i_mod}.model);
            training_data.pred_amps{i_mod} = cellfun(@(x,y) predict_vca_psc(x, d, dTau, f, fTau, mean(y(1,1,:))), training_data.ptimes_for_fit, training_data.amps_for_fit, 'uniformoutput', false);
            if xval_exists
                xval_data.pred_amps{i_mod} = cellfun(@(x,y) predict_vca_psc(x, d, dTau, f, fTau, mean(y(1,1,:))), xval_data.ptimes_for_fit, xval_data.amps_for_fit, 'uniformoutput', false);
            end
        end
        
        %
        % plot the training or cross validation data set, and the prediction
        %
        %%%%%%%%%%%%%%
        figure(hf)
        if PLOT_TRAINING_DATA
            amps_for_plot = training_data.amps_for_fit;
            ptimes_for_plot = training_data.ptimes_for_fit;
            pred_for_plot = training_data.pred_amps;
        else
            if xval_exists
                amps_for_plot = xval_data.amps_for_fit;
                ptimes_for_plot = xval_data.ptimes_for_fit;
                pred_for_plot = xval_data.pred_amps;
            else
                amps_for_plot = [];
                ptimes_for_plot = [];
                pred_for_plot = [];
            end
        end
        
        % make a scatter plot of all predicted and actual PSC amps
        plt_clr = lines(numel(pred_for_plot));
        [crossval_pred, crossval_raw, training_pred, training_raw] = deal([]);
        for i_mod = 1:numel(pred_for_plot)
            hs = [];
            training_pred = cellfun(@(x,y) repmat(x, [1,1,size(y,3)]), training_data.pred_amps{i_mod}, training_data.amps_for_fit, 'uniformoutput', false); % make sure there is one pred for every real amp
            training_pred = cellfun(@(x) x(:), training_pred, 'uniformoutput', false); %  a single col vector for each cond
            training_pred = cat(1, training_pred{:}); % now all conds aggregated into a single col vec
            training_raw = cellfun(@(x) x(:), training_data.amps_for_fit, 'uniformoutput', false);
            training_raw = cat(1, training_raw{:});
            l_ones = training_raw == 1; % valid hits when raw data are normalized
            training_raw(l_ones) = [];
            training_pred(l_ones) = [];
            if xval_exists
                crossval_pred =  cellfun(@(x,y) repmat(x, [1,1,size(y,3)]), xval_data.pred_amps{i_mod}, xval_data.amps_for_fit, 'uniformoutput', false);
                crossval_pred = cellfun(@(x) x(:), crossval_pred, 'uniformoutput', false);
                crossval_pred = cat(1, crossval_pred{:});
                crossval_raw = cellfun(@(x) x(:), xval_data.amps_for_fit, 'uniformoutput', false);
                crossval_raw = cat(1, crossval_raw{:});
                l_ones = crossval_raw == 1;
                crossval_raw(l_ones) = [];
                crossval_pred(l_ones) = [];
            end
            
            % calculate all the fitting errors
            num_params = numel(fit_results{i_mod}.params);
            fit_results{i_mod}.R2_train = get_r2(training_raw, training_pred);
            fit_results{i_mod}.R2_train_adj = get_r2(training_raw, training_pred, num_params);
            fit_results{i_mod}.AICc_train = get_aic(training_raw, training_pred, num_params);
            fit_results{i_mod}.MSE_train = mean((training_raw - training_pred).^2);
            [rho, p] = corr(training_raw, training_pred, 'type', 'Spearman');
            fit_results{i_mod}.corr_train = rho;
            fit_results{i_mod}.p_corr_train = p;
            if xval_exists
                fit_results{i_mod}.R2_crossvalid = get_r2(crossval_raw, crossval_pred);
                [rho, p] = corr(crossval_raw, crossval_pred, 'type', 'Spearman');
                fit_results{i_mod}.corr_crossvalid = rho;
                fit_results{i_mod}.p_corr_crossvalid = p;
            else
                fit_results{i_mod}.R2_crossvalid = nan;
                fit_results{i_mod}.corr_crossvalid = nan;
                fit_results{i_mod}.p_corr_crossvalid = nan;
            end

            
            if i_ch == 1; pltcol=1; else pltcol=4; end
            pltIdx = sub2ind([4, 3], pltcol, 1);
            hs = subplot(3, 4, pltIdx); hold on,
            plot(training_raw, training_pred, 'o', 'color', plt_clr(i_mod,:))
            if xval_exists
                plot(crossval_raw, crossval_pred, '.', 'color', plt_clr(i_mod,:))
            end
            
            maxval = max([hs.XLim, hs.YLim]);
            minval = min([hs.XLim, hs.YLim]);
            plot([minval, maxval], [minval, maxval], 'k--')
            xlabel('raw EPSC amp')
            ylabel('pred amp')
            
            pltIdx = sub2ind([4, 3], pltcol, 2);
            subplot(3,4,pltIdx), hold on,
            resid = training_raw - training_pred;
            histogram(resid, 'FaceColor', plt_clr(i_mod,:))
            plot(mean(resid), 10, 'rv', 'markerfacecolor', 'r')
            xlabel('real-pred')
            if i_mod == numel(pred_for_plot);
                title(sprintf('Best R2 = %.2f', max(cell2mat(structcat(fit_results, 'R2_train')))));
            end
            
            % plot cross validation stuff
            if xval_exists
                pltIdx = sub2ind([4, 3], pltcol, 3);
                subplot(3,4,pltIdx), hold on,
                resid = crossval_raw - crossval_pred;
                histogram(resid, 'FaceColor', plt_clr(i_mod,:));
                plot(mean(resid), 5, 'rv', 'markerfacecolor', 'r')
                xlabel('cross-valid (real-pred)')
                if i_mod == numel(pred_for_plot);
                    [maxval, maxidx] = max(cell2mat(structcat(fit_results, 'corr_crossvalid')));
                    pval = fit_results{maxidx}.p_corr_crossvalid;
                    title(sprintf('Best Rho= %.2f; p= %.2f', maxval, pval));
                end
            end
        end % i_mod scatter plots
        
        % main figures of raw data and predictions
        xlims = [inf -inf];
        ylims = [inf -inf];
        hs = [];
        for i_cond = 1:numel(amps_for_plot)
            
            if isempty(amps_for_plot{i_cond}); continue; end
            
            pltIdx = sub2ind([4, numel(amps_for_plot)], i_ch+1, i_cond);
            hs(i_cond) = subplot(numel(amps_for_plot), 4, pltIdx); hold on,
           
            % plot the model fits
            xx = ptimes_for_plot{i_cond};
            hp = [];
            model_string = {};
            for i_mod = 1:numel(pred_for_plot)
                hp(i_mod) = plot(xx, pred_for_plot{i_mod}{i_cond}, '-', 'color', plt_clr(i_mod,:), 'linewidth', 2);
                model_string{i_mod} = [fit_results{i_mod}.model, ',', blanks(8-numel(fit_results{i_mod}.model)-1) 'R2: ', num2str(fit_results{i_mod}.R2_train, 3)];
            end
            % plot the experimental data
            switch DATA_FOR_FIT
                case {'avg', 'norm'}
                    plot(xx, amps_for_plot{i_cond}, '-k.');
                case 'raw'
                    my_errorbar(xx, mean(amps_for_plot{i_cond},3), stderr(amps_for_plot{i_cond},3), '-k');
            end
            
            %legend(hp, model_string);
            axis tight
            xlims(1) = min([min(xx), xlims(1)]);
            xlims(2) = max([max(xx), xlims(2)]);
            yvals = get(gca, 'ylim');
            ylims(1) = min([yvals(1), ylims(1)]);
            ylims(2) = max([yvals(2), ylims(2)]);
            
        end
        if ~isempty(hs) && sum(hs)>0
            set(hs(hs~=0), 'YLim', ylims) % standardize y axis
            set(hs(hs~=0), 'XLim', xlims)  % standardize x axis
        end
        
        
        % store some parameters in the dat array
        dat{i_ex}.stpfits.trainingSet = TRAINSET;
        dat{i_ex}.stpfits.fitRecovPulse = FIT_RECOVERY_PULSE;
        dat{i_ex}.stpfits.dataforfit = DATA_FOR_FIT;
        dat{i_ex}.stpfits.force_recovery = FORCE_RECOVERY;
        dat{i_ex}.stpfits.convert_to_smooth_p1 = CONVERT_TO_SMOOTH_P1;
        dat{i_ex}.stpfits.fit_results{i_ch} = fit_results; % previously = [d, f, dTau, fTau], now has all models
        dat{i_ex}.stpfits.training_data{i_ch} = training_data;
        dat{i_ex}.stpfits.xval_data{i_ch} = xval_data;
        
    end
    drawnow
    clf
    
    
end



%% STP FITTING PLOTS RE-DO

TIME_IS_SEC = true;
PLOT_TRAINING_DATA = false;

for i_ex = 1:numel(dat)   %  i_ex = 17,46 is a good one for debugging.
    clc
    hf = figure;
    hf.Units = 'Normalized';
    hf.Position = [0.3347    0.0422    0.6639    0.7489];
    hf.Name = sprintf('Mouse %s, site %s, opsin: %s.  Train with: %s, Plot training set: %d',...
        dat{i_ex}.info.mouseName, dat{i_ex}.info.siteNum, dat{i_ex}.info.opsin, TRAINSET, PLOT_TRAINING_DATA);
    chempty = false(1,2);
    for i_ch = 1:2
        
        % determine if there are data to fit
        isvalid = dat{i_ex}.info.HS_is_valid_Vclamp(i_ch);
        if ~isvalid
            chempty(i_ch) = true;
            if  all(chempty)
                close(hf)
            end
            continue
        end
        
        % define training_data, xval_data
        fit_results = dat{i_ex}.stpfits.fit_results{i_ch}; % previously = [d, f, dTau, fTau], now has all models
        training_data = dat{i_ex}.stpfits.training_data{i_ch};
        xval_data = dat{i_ex}.stpfits.xval_data{i_ch};
        
        % determine if there are xval data
        xval_exists = ~isempty(xval_data.tdict);
        if ~PLOT_TRAINING_DATA && ~xval_exists
            continue
        end
        
        %
        % plot the training or cross validation data set, and the prediction
        %
        %%%%%%%%%%%%%%
        figure(hf)
        if PLOT_TRAINING_DATA
            amps_for_plot = training_data.amps_for_fit;
            ptimes_for_plot = training_data.ptimes_for_fit;
            pred_for_plot = training_data.pred_amps;
        else
            amps_for_plot = xval_data.amps_for_fit;
            ptimes_for_plot = xval_data.ptimes_for_fit;
            pred_for_plot = xval_data.pred_amps;
        end
        
        % make a scatter plot of all predicted and actual PSC amps
        plt_clr = lines(numel(pred_for_plot));
        for i_mod = 1:numel(pred_for_plot)
            hs = [];
            training_pred = cellfun(@(x,y) repmat(x, [1,1,size(y,3)]), training_data.pred_amps{i_mod}, training_data.amps_for_fit, 'uniformoutput', false); % make sure there is one pred for every real amp
            training_pred = cellfun(@(x) x(:), training_pred, 'uniformoutput', false); %  a single col vector for each cond
            training_pred = cat(1, training_pred{:}); % now all conds aggregated into a single col vec
            training_raw = cellfun(@(x) x(:), training_data.amps_for_fit, 'uniformoutput', false);
            training_raw = cat(1, training_raw{:});
            l_ones = training_raw == 1; % valid hits when raw data are normalized
            training_raw(l_ones) = [];
            training_pred(l_ones) = [];
            if xval_exists
                crossval_pred =  cellfun(@(x,y) repmat(x, [1,1,size(y,3)]), xval_data.pred_amps{i_mod}, xval_data.amps_for_fit, 'uniformoutput', false);
                crossval_pred = cellfun(@(x) x(:), crossval_pred, 'uniformoutput', false);
                crossval_pred = cat(1, crossval_pred{:});
                crossval_raw = cellfun(@(x) x(:), xval_data.amps_for_fit, 'uniformoutput', false);
                crossval_raw = cat(1, crossval_raw{:});
                l_ones = crossval_raw == 1;
                crossval_raw(l_ones) = [];
                crossval_pred(l_ones) = [];
            end
            
            if i_ch == 1; pltcol=1; else pltcol=4; end
            pltIdx = sub2ind([4, 3], pltcol, 1);
            hs = subplot(3, 4, pltIdx); hold on,
            plot(training_raw, training_pred, 'o', 'color', plt_clr(i_mod,:))
            if xval_exists
                plot(crossval_raw, crossval_pred, '.', 'color', plt_clr(i_mod,:))
            end
            
            maxval = max([hs.XLim, hs.YLim]);
            minval = min([hs.XLim, hs.YLim]);
            plot([minval, maxval], [minval, maxval], 'k--')
            xlabel('raw EPSC amp')
            ylabel('pred amp')
            
            pltIdx = sub2ind([4, 3], pltcol, 2);
            subplot(3,4,pltIdx), hold on,
            resid = training_raw - training_pred;
            histogram(resid, 'FaceColor', plt_clr(i_mod,:))
            plot(mean(resid), 10, 'rv', 'markerfacecolor', 'r')
            xlabel('real-pred')
            if i_mod == numel(pred_for_plot);
                title(sprintf('Best R2 = %.2f', max(cell2mat(structcat(fit_results, 'R2_train')))));
            end
            
            % plot cross validation stuff
            if xval_exists
                pltIdx = sub2ind([4, 3], pltcol, 3);
                subplot(3,4,pltIdx), hold on,
                resid = crossval_raw - crossval_pred;
                histogram(resid, 'FaceColor', plt_clr(i_mod,:));
                plot(mean(resid), 5, 'rv', 'markerfacecolor', 'r')
                xlabel('cross-valid (real-pred)')
                if i_mod == numel(pred_for_plot);
                    [maxval, maxidx] = max(cell2mat(structcat(fit_results, 'corr_crossvalid')));
                    pval = fit_results{maxidx}.p_corr_crossvalid;
                    title(sprintf('Best rho= %.2f p=%.2f', maxval, pval));
                end
            end
        end % i_mod scatter plots
        
        % main figures of raw data and predictions
        xlims = [inf -inf];
        ylims = [inf -inf];
        hs = [];
        for i_cond = 1:numel(amps_for_plot)
            
            if isempty(amps_for_plot{i_cond}); continue; end
            
            pltIdx = sub2ind([4, numel(amps_for_plot)], i_ch+1, i_cond);
            hs(i_cond) = subplot(numel(amps_for_plot), 4, pltIdx); hold on,
            
            % plot the model fits
            xx = ptimes_for_plot{i_cond};
            max_x = 7;
            min_x = 0;
            if ~TIME_IS_SEC
                xx = 1:numel(xx);
                max_x = numel(xx);
                min_x = 1;
            end
            hp = [];
            model_string = {};
            for i_mod = 1:numel(pred_for_plot)
                hp(i_mod) = plot(xx, pred_for_plot{i_mod}{i_cond}, '-', 'color', plt_clr(i_mod,:), 'linewidth', 2);
                model_string{i_mod} = [fit_results{i_mod}.model, ',', blanks(8-numel(fit_results{i_mod}.model)-1) 'R2: ', num2str(fit_results{i_mod}.R2_train, 3)];
            end
            % plot the experimental data
            switch DATA_FOR_FIT
                case {'avg', 'norm'}
                    plot(xx, amps_for_plot{i_cond}, '-k.');
                case 'raw'
                    my_errorbar(xx, mean(amps_for_plot{i_cond},3), stderr(amps_for_plot{i_cond},3), '-k');
            end
            
            %legend(hp, model_string);
            axis tight
            yvals = get(gca, 'ylim');
            ylims(1) = min([yvals(1), ylims(1)]);
            ylims(2) = max([yvals(2), ylims(2)]);
            xlim([min_x,max_x])
            
        end
        if ~isempty(hs) && sum(hs)>0
            set(hs(hs~=0), 'YLim', ylims) % standardize y axis
            %set(hs(hs~=0), 'XLim', xlims)  % standardize x axis
        end
        
        
    end
    drawnow
end

%% MODEL SELECTION 1: BEST FITTING MODELS
close all; clc;

NORMALIZATION = 'probability'; % how to normalize the histos. can be 'probability', or 'count'


% define a set of attributes for each analysis group
% {CellType, Layer,  BrainArea,  OpsinType}
% Brain Area can be: 'AL', 'PM', 'AM', 'LM', 'any', 'med', 'lat'. CASE SENSITIVE
plotgroups = {
    'PY', 'L23', 'any', 'chief';...
    'PY', 'L23', 'any', 'chronos';...
    'all_som', 'L23', 'any', 'chief';...
    'all_pv', 'L23', 'any', 'chief';...
    };
N_groups = size(plotgroups,2);

% initalize the aggregate datasets
error_types = {'R2_train', 'R2_train_adj', 'AICc_train', 'MSE_train'};
N_error_types = numel(error_types);
empty_arrays = repmat({[]}, 1, size(plotgroups, 1)); % should only have N cells, where N = size(plotgroups, 1).
group_data = [];
for i_err = 1:N_error_types
    group_data.(error_types{i_err}) = empty_arrays;
end

% initialze an analysis-wide model order. Start as empty, then fill in on
% the first looop. compare suscessive instances to this template
model_order_template = {};

for i_ex = 1:numel(dat)
    
    for i_ch = 1:2
        
        % check to make sure this neuron was defined
        isvalid = dat{i_ex}.info.HS_is_valid_Vclamp(i_ch);
        if ~isvalid
            continue
        end
        
        % check the attributes
        ch_attribs = {dat{i_ex}.info.cellType{i_ch}, upper(dat{i_ex}.info.brainArea), dat{i_ex}.info.opsin}; % force the brain area to be uppercase
        group_idx = groupMatcher(plotgroups, ch_attribs);
        if sum(group_idx) == 0; continue; end
        
        % get & check the model order for this experiment. Make sure that
        % the order is consistent from run to run
        ex_model_order = cellfun(@(x) x.model, dat{i_ex}.stpfits.fit_results{i_ch}, 'uniformoutput', false);
        if isempty(model_order_template)
            model_order_template = ex_model_order;
            N_models = numel(model_order_template);
        else
            model_matches = cellfun(@(x,y) strcmp(x,y), ex_model_order, model_order_template);
            assert(all(model_matches), 'ERROR: inconsistent model types found');
        end
        
        % concatenate data into the appropriate group for each of the error
        % types. one row for each Neuron, one column for each model type.
        for i_err = 1:N_error_types
            
            tmp_pop_array = group_data.(error_types{i_err}){group_idx}; % grab the existing group data
            
            tmp_ex_array = cellfun(@(x) x.(error_types{i_err}), dat{i_ex}.stpfits.fit_results{i_ch}); % err vals for all models
            tmp_pop_array = cat(1, tmp_pop_array, tmp_ex_array);
            
            group_data.(error_types{i_err}){group_idx} = tmp_pop_array; % update the group data
        end
        
        
    end
end


%
%  HISTOGRAMS OF BEST FITTING MODELS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = figure;
f.Position = [349          55        1152         913];
plt_clrs = parula(N_groups);
for i_err = 1:N_error_types
    
    % loop over group types, make normalized histogram over models
    normalized_histos = nan(N_groups, N_models); % each row is a histogram over model types
    legend_text = {};
    
    for i_group = 1:N_groups
        
        tmp_errs = group_data.(error_types{i_err}){i_group};
        
        % find the best fitting model
        switch error_types{i_err}
            case {'R2_train', 'R2_train_adj'}
                [~, best_model] = max(tmp_errs, [] , 2); % index to best model
            case {'MSE_train', 'AICc_train'}
                [~, best_model] = min(tmp_errs, [] , 2); % index to best model
        end
        
        % calculate the normalized histogram and put it in the population
        % array
        edges = (1:N_models+1)-0.5; % one bin for each model's index value
        counts = histcounts(best_model, edges, 'Normalization', NORMALIZATION);
        normalized_histos(i_group,:) = counts;
        
        legend_text{i_group} = sprintf('%s, %s, %s, %s, N=%d', plotgroups{i_group, 1}, plotgroups{i_group, 2}, plotgroups{i_group, 3}, plotgroups{i_group, 4}, size(group_data.(error_types{i_err}){i_group},1));
    end
    
    h_ax = subplot(N_error_types, 1, i_err);
    set(gca, 'ColorOrder', plt_clrs);
    h_bar = bar(1:N_models, normalized_histos');
    for i_group = 1:N_groups
        h_bar(i_group).FaceColor = plt_clrs(i_group,:);
    end
    h_ax.Box = 'off';
    h_ax.XLim = [0.5, N_models+0.5];
    h_ax.TickDir = 'out';
    h_ax.FontSize = 14;
    h_ax.Title.String = error_types{i_err};
    h_ax.Title.Interpreter = 'none';
    h_ax.YLabel.String = NORMALIZATION;
    h_ax.XTickLabel = model_order_template;
    if i_err == 1
        h_leg = legend(legend_text);
        h_leg.Location = 'Northwest';
        h_leg.Interpreter = 'none';
        h_leg.Box = 'off';
    end
    
end



% 
%  LINE SERIES OF AVERAGE R2 FOR EACH MODEL
%   
% R2 for most complex model should be superior to the others, so I'm
% looking to verify this in order to assert that the models are getting fit
% appropriately by my fitting routines (i.e., enough start points, good
% start points, etc.)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = figure;
f.Position = [349          55        1152         913];
for i_err = 1:N_error_types
    
    % loop over group types, aggregating the mean of the raw values
    legend_text = {};
    h_line = [];
    for i_group = 1:N_groups
        
        tmp_errs = group_data.(error_types{i_err}){i_group}; % one row per neuron, one column per model
        
        % find the best fitting model
        switch error_types{i_err}
            case {'R2_train', 'R2_train_adj', 'MSE_train'}
                % no special treatment yet for these error types.
            case {'AICc_train'}
                tmp_errs = bsxfun(@minus, tmp_errs, min(tmp_errs,[],2));
        end
        
        h_ax = subplot(N_error_types, 1, i_err);
        hold on,
        h_line(i_group) = my_errorbar(1:N_models, mean(tmp_errs,1), stderr(tmp_errs, 1), 'color', plt_clrs(i_group, :), 'linewidth', 2);
        
        legend_text{i_group} = sprintf('%s, %s, %s, %s', plotgroups{i_group, 1}, plotgroups{i_group, 2}, plotgroups{i_group, 3}, plotgroups{i_group, 4});
    end
    
    
    h_ax.Box = 'off';
    h_ax.XLim = [0.5, N_models+0.5];
    h_ax.TickDir = 'out';
    h_ax.FontSize = 14;
    h_ax.Title.String = error_types{i_err};
    h_ax.Title.Interpreter = 'none';
    h_ax.YLabel.String = 'err value';
    h_ax.XTick = 1:N_models;
    h_ax.XTickLabel = model_order_template;
    if i_err == 1
        h_leg = legend(h_line, legend_text);
        h_leg.Location = 'Southeast';
        h_leg.Interpreter = 'none';
        h_leg.Box = 'off';
    end
    
end



%% SUMMARY OF FIT PARAMETERS
clc

MODEL = 'ddff'; % only consider one of the models at a time (if many are present)

% define a set of attributes for each analysis group
% {CellType, Layer,  BrainArea,  OpsinType}
% Brain Area can be: 'AL', 'PM', 'AM', 'LM', 'any', 'med', 'lat'. CASE SENSITIVE
plotgroups = {
    'all_som', 'L23', 'any', 'any';...
    'PY', 'L23', 'lat', 'chief';...
%     'PY', 'L23', 'med', 'chronos';...
%     'PY', 'L23', 'lat', 'chronos';...
%     'all_som', 'any', 'any', 'any';...
%     'all_pv', 'any', 'any', 'any';...
    };

% initalize the aggregate datasets
empty_array = repmat({[]}, 1, size(plotgroups, 1)); 
[groupdata_fit_params, groupdata_training_MSE, groupdata_xval_MSE] = deal(empty_array);


for i_ex = 1:numel(dat)
    
    for i_ch = 1:2
        
        % check to make sure this neuron was defined
        isvalid = dat{i_ex}.info.HS_is_valid_Vclamp(i_ch);
        if ~isvalid
            continue
        end
        
        % check the attributes
        ch_attribs = {dat{i_ex}.info.cellType{i_ch}, upper(dat{i_ex}.info.brainArea), dat{i_ex}.info.opsin}; % force the brain area to be uppercase
        group_idx = groupMatcher(plotgroups, ch_attribs);
        if sum(group_idx) == 0; continue; end
        
        % find the index to the model requested
        tmp_fit_results = dat{i_ex}.stpfits.fit_results{i_ch};
        mod_idx = cellfun(@(x) ~isempty(strcmpi(x.model, MODEL)), tmp_fit_results);
        
        % add data to the appropriate group data array
        assert(numel(tmp_fit_results{mod_idx}.params) == numel(MODEL)*2, 'ERROR: params mismatch')
        groupdata_fit_params{group_idx} = cat(1, groupdata_fit_params{group_idx}, tmp_fit_results{mod_idx}.params);
        
    end
end


% fix the order of the params: ascending order
for i_group = 1:size(plotgroups,1)
    tmp_dat = groupdata_fit_params{i_group};

    for ii = 1:size(tmp_dat,1)
        [d, tau_d, f, tau_f] = parse_vca_model_coeffs(tmp_dat(ii,:), MODEL);
        [~, d_sort_idx] = sort(d, 2, 'ascend');
        [~, f_sort_idx] = sort(f, 2, 'ascend');
        d = d(d_sort_idx);
        tau_d = tau_d(d_sort_idx);
        f = f(f_sort_idx);
        tau_f = tau_f(f_sort_idx);
        tmp_dat(ii,:) = [d, tau_d, f, tau_f];        
    end
    groupdata_fit_params{i_group} = tmp_dat;
end

% make a summary table of the average fit params
T = table;
for i_group = 1:size(plotgroups,1)
    tmp_dat = groupdata_fit_params{i_group};
    avg_params = mean(tmp_dat, 1);
    [d, tau_d, f, tau_f] = parse_vca_model_coeffs(avg_params, MODEL);
    cell_type = plotgroups(i_group, 1);
    opsin_type = plotgroups(i_group, 4);
    brain_area = plotgroups(i_group, 3);
    N_expts = size(tmp_dat,1);
    T = [T; table(cell_type, brain_area, opsin_type, d, tau_d, f, tau_f, N_expts)];
end
T


% line plots of parameter fits
figure, hold on,
for i_group = 1:size(plotgroups,1)
    plt_clr = hvaPlotColor(plotgroups{i_group, 3});
    tmp_dat = groupdata_fit_params{i_group};
    
    xbar = mean(tmp_dat,1);
    sem = stderr(tmp_dat,1);
    %p1 = plot(tmp_dat', '-', 'color', plotcolors{i_group}, 'linewidth', 0.5);
    %for ii = 1:numel(p1); p1(ii).Color(4) = 0.5; end
    %plot(xbar, '-', 'color', plotcolors{i_group}, 'linewidth', 4)
    shadedErrorBar(1:numel(MODEL)*2, xbar, sem, {'color', plt_clr, 'linewidth', 4}, true);   
    %set(gca, 'yscale', 'log')
    %set(gca, 'xtick', [1:size(tmp_dat,2)], 'xticklabel', {'d1', 'd2', 'Tau d1', 'Tau d2', 'f', 'Tau f1'})
    set(gca, 'TickDir', 'out')
end

% error('Does not work with arbitrary model types')
% figure, hold on,
% for i_group = 1:size(plotgroups,1)
%     tmp_dat = groupdata_fit_params{i_group};
%     
%     xbar = mean(tmp_dat,1);
%     sem = stderr(tmp_dat,1);
%     
%     % D terms
%     subplot(1,4,1), hold on,
%     p1 = plot(tmp_dat(:,1:2)', '-', 'color', plotcolors{i_group}, 'linewidth', 0.5);
%     for ii = 1:numel(p1); p1(ii).Color(4) = 0.5; end
%     plot(xbar(1:2), '-', 'color', plotcolors{i_group}, 'linewidth', 4)
%     %set(gca, 'yscale', 'log')
%     set(gca, 'xtick', [1:2], 'xticklabel', {'d1', 'd2'})
%     set(gca, 'TickDir', 'out')
%     
%     % D taus
%     subplot(1,4,2), hold on,
%     p1 = plot(tmp_dat(:,3:4)', '-', 'color', plotcolors{i_group}, 'linewidth', 0.5);
%     for ii = 1:numel(p1); p1(ii).Color(4) = 0.5; end
%     plot(xbar(3:4), '-', 'color', plotcolors{i_group}, 'linewidth', 4)
%     %set(gca, 'yscale', 'log')
%     set(gca, 'xtick', [1:2], 'xticklabel', {'Tau d1', 'Tau d2'})
%     set(gca, 'TickDir', 'out')
%     
%     % F terms
%     subplot(1,4,3), hold on,
%     p1 = plot(ones(size(tmp_dat,1)), tmp_dat(:,5)', 'o', 'markeredgecolor', plotcolors{i_group}, 'linewidth', 0.5);
%     for ii = 1:numel(p1); p1(ii).Color(4) = 0.5; end
%     plot(1, xbar(5), 'o', 'markeredgecolor', plotcolors{i_group}, 'markerfacecolor', plotcolors{i_group}, 'linewidth', 4)
%     %set(gca, 'yscale', 'log')
%     set(gca, 'xtick', [1], 'xticklabel', {'f'})
%     set(gca, 'TickDir', 'out')
%     
%     % F taus
%     subplot(1,4,4), hold on,
%     p1 = plot(ones(size(tmp_dat,1)), tmp_dat(:,6)', 'o', 'markeredgecolor', plotcolors{i_group}, 'linewidth', 0.5);
%     for ii = 1:numel(p1); p1(ii).Color(4) = 0.5; end
%     plot(1, xbar(6), 'o', 'markeredgecolor', plotcolors{i_group}, 'markerfacecolor', plotcolors{i_group}, 'linewidth', 4)
%     %set(gca, 'yscale', 'log')
%     set(gca, 'xtick', [1], 'xticklabel', {'Tau f1'})
%     set(gca, 'TickDir', 'out')
% end


%% PLOT R2 AND CORR FOR TRAIN/TEST


train_stat = 'R2_train'; % 'corr_train' or 'R2_train'
test_stat = 'corr_crossvalid'; % 'corr_crossvalid' or 'R2_crossvalid'

cell_type = 'PY_L23';
opsin = 'chief';

[train_vals, test_vals, p_test] = deal([]);
param_norms = [];

for i_ex = 1:numel(dat)
    for i_ch = 1:2
        
        isvalid = dat{i_ex}.info.HS_is_valid_Vclamp(i_ch);
        if ~isvalid
            continue
        end
        
        cell_type_match = strcmpi(dat{i_ex}.info.cellType{i_ch}, cell_type);
        opsin_match = ~isempty(regexpi(dat{i_ex}.info.opsin, opsin));
        if ~(cell_type_match && opsin_match)
            continue
        end
        
        
        train_vals = cat(1, train_vals, dat{i_ex}.stpfits.fit_results{i_ch}{1}.(train_stat));
        test_vals = cat(1, test_vals, dat{i_ex}.stpfits.fit_results{i_ch}{1}.(test_stat));
        
        if ~isempty(regexpi(test_stat, 'corr'))
            p_test_fld = ['p_', test_stat];
            p_test = cat(1, p_test, dat{i_ex}.stpfits.fit_results{i_ch}{1}.(p_test_fld));
        end
        
        tmp_norm =  norm(dat{i_ex}.stpfits.fit_results{i_ch}{1}.params);
        param_norms = cat(1, param_norms, tmp_norm);
        
    end
end

%figure, hold on,
plot(train_vals, test_vals, 'bo', 'markersize', 6)
plot(train_vals(p_test<0.05), test_vals(p_test<0.05), 'bo', 'markerfacecolor', 'b', 'markersize', 6)
ax_vals = [get(gca, 'xlim'); get(gca, 'ylim')];
ax_min = min(ax_vals(:,1));
ax_max = max(ax_vals(:,2));
set(gca, 'xlim', [ax_min, ax_max], 'ylim', [ax_min, ax_max], 'tickdir', 'out')
xlabel(train_stat, 'fontsize', 14, 'interpreter', 'none')
ylabel(test_stat, 'fontsize', 14, 'interpreter', 'none')
[rho, p] = corr(train_vals, test_vals, 'type', 'Spearman');
title(sprintf('Corr of train and test: %.2f  p: %.2f', rho, p), 'fontsize', 16)

%% PAIRED PULSE PLASTICITY MANIFOLDS (DATA COLLECTION)
% TODO convert raw data used here to the data explicitly used in the fits.
% (including the smoothed version)

%error('Does not generalize to arbitrary fit params')
MODEL = 'ddff';

% loop through the experiments. Pull out the trains data. Ignore the
% recovery train (if present) and aggregate across recovery conditions.
pprpop = [];
pprpop.TFsAllExpts = [];
pprpop.MaxNPulses = 0;
pprpop.smoothManifold_isi_ms = fliplr([1000/60 : 2 : 1000/.1]);
pprpop.smoothManifold_isi_ms = fliplr(logspace(log10(20), log10(100), 150));
pprpop.smoothManifold_numPulses = 10;
for i_ex = 1:numel(dat)
    
    % assume that the model was trained using the recovery trains data
    %assert(strcmpi(dat{i_ex}.stpfits.trainingSet, 'recovery'), 'ERROR: model needs to be fit with the recovery trains')
    
    % store some metadata
    pprpop.info{i_ex} = dat{i_ex}.info;
    
    % pre allocate the outputs in case a channel is not defined
    pprpop.dat{i_ex}.xbar = {[],[]};
    pprpop.tfs{i_ex} = {[],[]};
    pprpop.smoothManifold{i_ex} = {[],[]};
    pprpop.params{i_ex} = {[],[]};
    pprpop.R2{i_ex}= {[],[]};
    
    % loop over channels
    for i_ch = 1:2
        
        % was this channel defined? If not, set the default outputs
        isvalid = dat{i_ex}.info.HS_is_valid_Vclamp(i_ch);
        if ~isvalid
            continue
        end
       
        % pull out the ttype dictionary for the traning data. These will be
        % constrained to be trains or recovery trains, but could have
        % multiple recovery times for a single base train TF
        all_training_conds = dat{i_ex}.stpfits.training_data{i_ch}.unique_conds_tdict;
        unique_tconds = unique(all_training_conds, 'rows');
        l_rit = unique_tconds(:,5) > 0;
        unique_tconds(l_rit,:) = [];
        %assert(~any(unique_tconds(:,5) >0), 'ERROR: found a RIT train')
        
        % identify the unique TF conditions for this experiment, and update the
        % running log of TFs used across all experiments
        uniqueTFs = unique(unique_tconds(:,3));
        tmp = cat(1, pprpop.TFsAllExpts, uniqueTFs);
        pprpop.TFsAllExpts = unique(tmp);
        
        % aggregate data within each unique TF
        for i_tf = 1:numel(uniqueTFs)
            
            catdat = [];
            list_of_conds_at_specific_tf = find(all_training_conds(:,3) == uniqueTFs(i_tf));
            for i_cond = 1:numel(list_of_conds_at_specific_tf)
                cond_idx = list_of_conds_at_specific_tf(i_cond);
                catdat = cat(3, catdat, dat{i_ex}.stpfits.training_data{i_ch}.amps_for_fit{cond_idx});
            end
            
            pprpop.dat{i_ex}.xbar{i_ch}{i_tf} = mean(catdat,3);
            pprpop.tfs{i_ex}{i_ch} = uniqueTFs;
            pprpop.MaxNPulses = max([pprpop.MaxNPulses, size(catdat, 1)]);
            
        end
        
        
        % make a smooth manifold for this neuron.
        tmp_fit_results = dat{i_ex}.stpfits.fit_results{i_ch};
        mod_idx = cellfun(@(x) ~isempty(strcmpi(x.model, MODEL)), tmp_fit_results);
        assert(numel(tmp_fit_results{mod_idx}.params) == numel(MODEL)*2, 'ERROR: params mismatch')
        params = tmp_fit_results{mod_idx}.params;
        isi_ms = pprpop.smoothManifold_isi_ms;
        NumPulses = pprpop.smoothManifold_numPulses;
        
        % solve for P1 amp
        all_p1_amps = cat(3, pprpop.dat{i_ex}.xbar{i_ch}{:});
        A0 = mean(all_p1_amps(1,1,:), 3);
        
        smoothManifold = nan(NumPulses, numel(isi_ms));
        for i_isi = 1:numel(isi_ms)
            tmp_pOntimes_ms = 0 : isi_ms(i_isi) : isi_ms(i_isi)*(NumPulses-1);
            tmp_pOntimes_sec = tmp_pOntimes_ms ./ 1000;
            [d, tau_d, f, tau_f] = parse_vca_model_coeffs(params, MODEL);
            smoothManifold(:,i_isi) = predict_vca_psc(tmp_pOntimes_sec, d, tau_d, f, tau_f, A0);
        end
        pprpop.smoothManifold{i_ex}{i_ch} = smoothManifold;
        pprpop.params{i_ex}{i_ch} = params;
        
    end
    
end



%% PAIRED PULSE PLASTICITY MANIFOLDS (PLOTS)
clc;

% switches for first figure
PLOT_AVG_SMOOTH_MANIFOLD = true;
PLOT_BOOTSTRAP_SMOOTH_MANIFOLD = false;
PLOT_MANIFOLD_OF_AVG_PARAMS = false;
PLOT_MINIFOLD_FROM_FITTED_GRAND_AVERAGE_MANIFOLD = false;
PLOT_MANIFOLD_OF_AVG_RAW_DATA = false;
PLOT_OVERLAY_ALL_MANIFOLDS = false;
PLOT_AVG_RAW_DATA = true;

% add other plots?
PLOT_DATASETS_INDIVIDUALLY = false;

% define a set of attributes for each analysis group
% {CellType, Layer,  BrainArea,  OpsinType}
% Brain Area can be: 'AL', 'PM', 'AM', 'LM', 'any', 'med', 'lat'. CASE SENSITIVE
man_plotgrps = {
    'PY', 'L23', 'PM', 'chief';...
    'PY', 'L23', 'AM', 'chief';...
    'PY', 'L23', 'LM', 'chief';...
    'PY', 'L23', 'AL', 'chief';...
    };

empty_array = repmat({[]}, 1, size(man_plotgrps, 1)); % should only have N cells, where N = size(man_plotgrps, 1). Each cell has a matrix with a cononicalGrid:
[groupdata_raw, groupdata_smooth, groupexpinds, groupchinds, group_params] = deal(empty_array);
canonicalGrid = nan(pprpop.MaxNPulses, numel(pprpop.TFsAllExpts)); % need to get the maxNpulses into the pprpop struct
allTFs = pprpop.TFsAllExpts;

% iterate over the experiments. For each recording channel, determine whata
% the attributes are, and place the data in the correct ploting group.
for i_ex = 1:numel(dat)
    
    for i_ch = 1:2
        
        % check to make sure this neuron was defined
        isvalid = dat{i_ex}.info.HS_is_valid_Vclamp(i_ch);
        if ~isvalid
            continue
        end
        
        % check the attributes
        ch_attribs = {dat{i_ex}.info.cellType{i_ch}, upper(dat{i_ex}.info.brainArea), dat{i_ex}.info.opsin}; % force the brain area to be uppercase
        group_idx = groupMatcher(man_plotgrps, ch_attribs);
        if sum(group_idx) == 0; continue; end
        
        % add data to the appropriate group data array
        tmpgrid = canonicalGrid;
        ch_tfs = pprpop.tfs{i_ex}{i_ch};
        for i_tf = 1:numel(ch_tfs)
            grid_col_idx = allTFs == ch_tfs(i_tf);
            tmpdat = pprpop.dat{i_ex}.xbar{i_ch}{i_tf};
            Npulses = numel(tmpdat);
            tmpgrid(1:Npulses, grid_col_idx) = tmpdat(:);
        end
        
        if all(isnan(tmpgrid(:))); continue; end
        
        % aggregate the data
        groupdata_raw{group_idx} = cat(3, groupdata_raw{group_idx}, tmpgrid);
        groupexpinds{group_idx} = cat(1, groupexpinds{group_idx}, i_ex);
        groupchinds{group_idx} = cat(1, groupchinds{group_idx}, i_ch);
        
        if PLOT_AVG_SMOOTH_MANIFOLD || PLOT_OVERLAY_ALL_MANIFOLDS || PLOT_MANIFOLD_OF_AVG_RAW_DATA || PLOT_BOOTSTRAP_SMOOTH_MANIFOLD
            groupdata_smooth{group_idx} = cat(3, groupdata_smooth{group_idx}, pprpop.smoothManifold{i_ex}{i_ch});
            group_params{group_idx} = cat(1, group_params{group_idx}, pprpop.params{i_ex}{i_ch});
        end
    end
end


% Plot manfold(s)
f = figure;
hold on,
for i_group = 1:numel(groupdata_raw)
    plt_clr = hvaPlotColor(man_plotgrps{i_group, 3});
    
    if PLOT_AVG_RAW_DATA
        grid_average = nanmean(groupdata_raw{i_group},3);
        grid_average = grid_average(1:10, :);
        grid_N = sum(~isnan(groupdata_raw{i_group}),3)
        
        Y = 10:-1:1;
        X = allTFs';
        
        l_nan_tfs = all(isnan(grid_average), 1);
        
        hs = surf(X(:,~l_nan_tfs), Y, flipud(grid_average(:,~l_nan_tfs)));
        hs.EdgeColor = plt_clr;
        hs.EdgeAlpha = 1;
        hs.LineWidth = 1.5;
        hs.FaceColor = plt_clr;
        hs.FaceAlpha = 0.1;
    end
    
    % plot the average smoothManifold
    isi_ms = pprpop.smoothManifold_isi_ms;
    N_pulses_smooth = pprpop.smoothManifold_numPulses;
    X = 1000./isi_ms;
    Y = N_pulses_smooth:-1:1;
    
    if PLOT_AVG_SMOOTH_MANIFOLD
        grid_average = nanmean(groupdata_smooth{i_group}, 3);
        hmod = surf(X,Y, flipud(grid_average));
        hmod.EdgeAlpha = 0;
        hmod.FaceColor = plt_clr;
        hmod.FaceAlpha = 0.6;
    end
    
    if PLOT_BOOTSTRAP_SMOOTH_MANIFOLD
        sample_surfs = groupdata_raw{i_group};
        N_samps = size(sample_surfs,3);
        N_boots = 10000;
        btstrp_idx = randi(N_samps, N_boots, 1);
        grid_average = nanmean(sample_surfs(:,:,btstrp_idx), 3);

        hmod = surf(X,Y, flipud(grid_average));
        hmod.EdgeAlpha = .2;
        hmod.FaceColor = 'k';
        hmod.EdgeColor = 'k';
        hmod.FaceAlpha = 0;
    end
    
    if PLOT_MANIFOLD_OF_AVG_PARAMS
        avg_params = mean(group_params{i_group}, 1);
        [d, tau_d, f, tau_f] = parse_vca_model_coeffs(avg_params, MODEL);
        A0 = 1;
        smoothManifold = nan(N_pulses_smooth, numel(isi_ms));
        for i_isi = 1:numel(isi_ms)
            tmp_pOntimes_ms = 0 : isi_ms(i_isi) : (isi_ms(i_isi)*N_pulses_smooth)-1;
            tmp_pOntimes_sec = tmp_pOntimes_ms ./ 1000;
            smoothManifold(:,i_isi) = predict_vca_psc(tmp_pOntimes_sec, d, tau_d, f, tau_f, A0);
        end
        predmod = surf(X,Y, flipud(smoothManifold));
        predmod.EdgeAlpha = 0;
        predmod.FaceColor = plt_clr;
        predmod.FaceAlpha = 0.2;
    end
    
    if PLOT_OVERLAY_ALL_MANIFOLDS
        N_manifolds = size(groupdata_smooth{i_group}, 3);
        for i_examp = 1:N_manifolds
            hmod = surf(X,Y, flipud(groupdata_smooth{i_group}(:,:,i_examp)));
            hmod.EdgeAlpha = .2;
            hmod.FaceColor = plt_clr;
            hmod.FaceAlpha = 0.5;
        end
    end
    
end
%set(gca, 'zlim', [0.45, 1.4])
set(gca, 'zscale', 'linear', 'view', [177 20])
set(gca, 'Xdir', 'reverse')
zmax = get(gca, 'zlim');
set(gca, 'XGrid', 'on', 'Ygrid', 'on', 'Zgrid', 'on')
set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'off', 'ZMinorGrid', 'off')
xlabel('Temporal Frequency')
ylabel('Pulse Number')
zlabel('norm amp')
set(gca, 'zscale', 'log')

% plot manifold of params fit to grand average
if PLOT_MINIFOLD_FROM_FITTED_GRAND_AVERAGE_MANIFOLD
    for i_group = 1:numel(groupdata_raw)
        plt_clr = hvaPlotColor(man_plotgrps{i_group, 3});
        
        % remake the smooth manifold
        isi_ms = pprpop.smoothManifold_isi_ms;
        N_pulses_smooth = pprpop.smoothManifold_numPulses;
        X = 1000./isi_ms;
        Y = N_pulses_smooth:-1:1;
        [d, tau_d, f, tau_f] = parse_vca_model_coeffs(grand_avg_params{i_group}, MODEL);
        A0 = 1;
        smoothManifold = nan(N_pulses_smooth, numel(isi_ms));
        for i_isi = 1:numel(isi_ms)
            tmp_pOntimes_ms = 0 : isi_ms(i_isi) : (isi_ms(i_isi)*N_pulses_smooth)-1;
            tmp_pOntimes_sec = tmp_pOntimes_ms ./ 1000;
            smoothManifold(:,i_isi) = predict_vca_psc(tmp_pOntimes_sec, d, tau_d, f, tau_f, A0);
        end
        
        %re-plot the smooth manifold
        predmod = surf(X,Y, flipud(smoothManifold));
        predmod.EdgeAlpha = 0.5;
        predmod.EdgeColor = plt_clr;
        predmod.FaceColor = plt_clr;
        predmod.FaceAlpha = 0;
        
    end
end

if PLOT_MANIFOLD_OF_AVG_RAW_DATA
    for i_group = 1:size(man_plotgrps,1)
        
        % check that the plot_group defined here is the same as the
        % plot-group defined by the group_fit_to_grand_avg dataset
        data_matches = all(strcmpi(man_plotgrps(i_group,:), group_fit_to_grand_avg{i_group}.grouptypes));
        assert(data_matches, 'ERROR: plot groups do not match')
        
        isi_ms = pprpop.smoothManifold_isi_ms;
        N_pulses_smooth = pprpop.smoothManifold_numPulses;
        X = 1000./isi_ms;
        Y = N_pulses_smooth:-1:1;
        smoothManifold = group_fit_to_grand_avg{i_group}.smooth_manifold;
        
        %plot the smooth manifold
        predmod = surf(X,Y, flipud(smoothManifold));
        predmod.EdgeAlpha = 0.5;
        predmod.EdgeColor = 'k';
        predmod.FaceColor = 'k';
        predmod.FaceAlpha = 0;
        
    end
    set(gca, 'Xdir', 'reverse')
end

%
% plot manifolds for each dataset individually
% (lots of figures)
%%%%%%%
if PLOT_DATASETS_INDIVIDUALLY
    for i_group = 1:numel(groupdata_raw)
        plt_clr = hvaPlotColor(man_plotgrps{i_group, 3});

        for i_examp = 1:size(groupdata_raw{i_group},3)
            f = figure;
            
            grid_average = groupdata_raw{i_group}(:,:,i_examp);
            
            Y = 1:size(groupdata_raw{i_group},1);
            X = allTFs';
            
            l_nan_tfs = all(isnan(grid_average), 1);
            
            hs = surf(X(:,~l_nan_tfs), Y, flipud(grid_average(:,~l_nan_tfs)));
            hs.EdgeColor = plt_clr;
            hs.EdgeAlpha = 1;
            hs.FaceColor = plt_clr;
            hs.FaceAlpha = 0;
            hs.LineWidth = 1.5;
            
            set(gca, 'zscale', 'linear', 'view', [-43    16])
            set(gca, 'ytick', Y, 'yticklabel', cellfun(@(x) num2str(x), num2cell(rot90(Y,2)), 'uniformoutput', false))
            set(gca, 'xtick', X, 'xticklabel', cellfun(@(x) num2str(x), num2cell(allTFs), 'uniformoutput', false))
            xlabel('Temporal Frequency')
            ylabel('Pulse Number')
            zlabel('norm amp')
            
            i_ex = groupexpinds{i_group}(i_examp);
            i_ch = groupchinds{i_group}(i_examp);
            %f.Name = sprintf('Mouse %s, site %s, HS%d, R2: %.2f', dat{i_ex}.info.mouseName, dat{i_ex}.info.siteNum, i_ch, pprpop.R2{i_ex}{i_ch});
            
            % display the smooth manifold
            if PLOT_AVG_SMOOTH_MANIFOLD
                smoothManifold = pprpop.smoothManifold{i_ex}{i_ch};
                isi_ms = pprpop.smoothManifold_isi_ms;
                N_pulses_smooth = pprpop.smoothManifold_numPulses;
                X = 1000./isi_ms;
                Y = N_pulses_smooth:-1:1;
                hold on,
                hmod = surf(X,Y, flipud(smoothManifold));
                hmod.EdgeAlpha = 0;
                hmod.FaceAlpha = 0.5;
            end
            
            
            
        end
        
    end
end

%% PLOT THE RATIO OF THE MANIFOLDS
close all; clc

% plot the average smoothManifold
isi_ms = pprpop.smoothManifold_isi_ms;
N_pulses_smooth = pprpop.smoothManifold_numPulses;
X = 1000./isi_ms;
Y = N_pulses_smooth:-1:1;

grid_averages = {};
min_val = Inf;
max_val = -Inf;
max_abs_diff = 0;
for i_group = 1:numel(groupdata_raw)
    grid_averages{i_group} = nanmean(groupdata_smooth{i_group}, 3);
    min_val = min(min_val, min(reshape(grid_averages{i_group}(2:end,:), 1, [])));
    max_val = max(max_val, max(reshape(grid_averages{i_group}(2:end,:), 1, [])));
    max_abs_diff = max(abs(min_val-1), abs(max_val-1));
end

plt_x = round(linspace(1, numel(X), 5));
plt_x_labels = round(X(plt_x));

f = figure;
f.Position = [17         508        1391         273];
for i_group = 1:numel(grid_averages)
    subplot(1, numel(grid_averages)+1, i_group)
    manifold = grid_averages{i_group}(2:end, :);
    imagesc(manifold)
    set(gca, 'xtick', plt_x, 'xticklabels', plt_x_labels)
    set(gca, 'ytick', [1:9], 'yticklabel', [2:10]) 
    axis tight
    colorbar
    colormap(brewermap(256, '*RdYlBu'))
%     caxis([1-max_abs_diff, 1+max_abs_diff])
    xlabel('temporal frequency')
    ylabel('pulse number')
    title(sprintf('HVA: %s', man_plotgrps{i_group,3}))
end

ratio = grid_averages{1}(2:end, :) ./ grid_averages{2}(2:end, :);
diffvals = grid_averages{1}(2:end, :) - grid_averages{2}(2:end, :);

y_raw = [1:9];
y_interp = 1:0.01:9;
x_raw = 1000./isi_ms;
x_interp = x_raw;

[X,Y] = meshgrid(x_raw, y_raw);
diffvals = interp2(X, Y, diffvals, x_interp, y_interp(:));
plt_max_val = max(abs(diffvals(:)));

subplot(1, numel(grid_averages)+1, numel(grid_averages)+1)
imagesc(diffvals)
set(gca, 'xtick', plt_x, 'xticklabels', plt_x_labels)
set(gca, 'ytick', [1:100:900], 'yticklabel', [2:10]) 
axis tight
colorbar
% caxis([-plt_max_val, plt_max_val])
xlabel('temporal frequency')
ylabel('pulse number')
title('ratio')




%%  FIT THE "GRAND AVERAGE" MANIFOLDS
% the manifold of the average params is dead wrong. Instead, find the 
% best fitting params of the "average" manifold.
grand_avg_params={};
for i_group = 1:size(man_plotgrps,1)
    
    % calculate the average across manifolds
    avg_manifold = mean(groupdata_smooth{i_group}, 3);
    
    % fit the grand-average manifold with the VCA model
    isi_ms = pprpop.smoothManifold_isi_ms;
    N_pulses_smooth = pprpop.smoothManifold_numPulses;
    X = 1000./isi_ms;
    Y = N_pulses_smooth:-1:1;
    grand_avg_amps = {};
    grand_avg_p_times = {};
    for i_isi = 1:numel(isi_ms)
        tmp_pOntimes_ms = 0 : isi_ms(i_isi) : (isi_ms(i_isi)*N_pulses_smooth)-1;
        grand_avg_p_times{i_isi} = tmp_pOntimes_ms ./ 1000;
        grand_avg_amps{i_isi} = avg_manifold(:,i_isi);
    end
    grand_avg_params{i_group} = fit_vca_model(grand_avg_amps, grand_avg_p_times, MODEL);
end

% SUMMARY TABLE OF GRAND AVERAGE PARAMS
clc
T = table;
for i_group = 1:numel(groupdata_raw)
    avg_params = grand_avg_params{i_group};
    [d, tau_d, f, tau_f] = parse_vca_model_coeffs(avg_params, MODEL);
    cell_type = man_plotgrps(i_group, 1);
    opsin_type = man_plotgrps(i_group, 4);
    brain_area = man_plotgrps(i_group, 3);
    T = [T; table(cell_type, brain_area, opsin_type, d, tau_d, f, tau_f)];
end
T_fit_to_avg_manifold = T


%% FIT THE MODEL TO THE AVERAGE RAW DATA


% initialize the output cell array
empty_struct.amps_for_fit = {};
empty_struct.ptimes_for_fit = {};
empty_struct.model = {};
empty_struct.params = {};
empty_struct.grouptypes = {};
group_fit_to_grand_avg = repmat({empty_struct}, size(man_plotgrps,1), 1);

% aggregate data across plotgroups
for i_ex = 1:numel(dat)

    for i_ch = 1:2
        % check to make sure this neuron was defined
        isvalid = dat{i_ex}.info.HS_is_valid_Vclamp(i_ch);
        if ~isvalid
            continue
        end
        
        % check the attributes
        ch_attribs = {dat{i_ex}.info.cellType{i_ch}, dat{i_ex}.info.brainArea, dat{i_ex}.info.opsin};
        group_idx = groupMatcher(man_plotgrps, ch_attribs);
        if sum(group_idx) == 0; continue; end
        
%         % check to make sure the amplitudes exceeded a certian criterion
%         if regexpi(dat{i_ex}.info.cellType{i_ch}, 'PY_L23')
%             p1amp_avg = nanmean(dat{i_ex}.qc.p1amp{i_ch});
%             if p1amp_avg < 200
%                 warning('culling data on basis of p1amp')
%                 continue
%             end
%         end
        
        % pull out the data used for fitting (from before) and place it in
        % the correct group_data array
        tmp_amps = group_fit_to_grand_avg{group_idx}.amps_for_fit;
        ex_amps = dat{i_ex}.stpfits.training_data{i_ch}.amps_for_fit;
        tmp_amps = cat(1, tmp_amps, ex_amps');
        group_fit_to_grand_avg{group_idx}.amps_for_fit = tmp_amps;
        
        % pull out the p_times used for fitting (from before) and place it in
        % the correct group_data array
        tmp_ptimes = group_fit_to_grand_avg{group_idx}.ptimes_for_fit;
        ex_ptimes = dat{i_ex}.stpfits.training_data{i_ch}.ptimes_for_fit;
        tmp_ptimes = cat(1, tmp_ptimes, ex_ptimes');
        group_fit_to_grand_avg{group_idx}.ptimes_for_fit = tmp_ptimes;

    end
end

% fit the population datasets
for i_group = 1:size(man_plotgrps, 1)
    group_fit_to_grand_avg{i_group}.model = MODEL;
    tmp_amps = group_fit_to_grand_avg{i_group}.amps_for_fit;
    tmp_ptimes = group_fit_to_grand_avg{i_group}.ptimes_for_fit;
    params = fit_vca_model(tmp_amps, tmp_ptimes, MODEL);
    group_fit_to_grand_avg{i_group}.params = params;
    
    % make a smooth surface for plotting
    isi_ms = pprpop.smoothManifold_isi_ms;
    N_pulses_smooth = pprpop.smoothManifold_numPulses;
    X = 1000./isi_ms;
    Y = N_pulses_smooth:-1:1;
    [d, tau_d, f, tau_f] = parse_vca_model_coeffs(params, MODEL);
    A0 = 1;
    smoothManifold = nan(N_pulses_smooth, numel(isi_ms));
    for i_isi = 1:numel(isi_ms)
        tmp_pOntimes_ms = 0 : isi_ms(i_isi) : (isi_ms(i_isi)*N_pulses_smooth)-1;
        tmp_pOntimes_sec = tmp_pOntimes_ms ./ 1000;
        smoothManifold(:,i_isi) = predict_vca_psc(tmp_pOntimes_sec, d, tau_d, f, tau_f, A0);
    end
    group_fit_to_grand_avg{i_group}.smooth_manifold = smoothManifold;
    group_fit_to_grand_avg{i_group}.grouptypes = man_plotgrps(i_group,:);
end


% SUMMARY TABLE OF GRAND AVERAGE PARAMS
clc
T = table;
for i_group = 1:size(man_plotgrps, 1)
    avg_params = group_fit_to_grand_avg{i_group}.params;
    [d, tau_d, f, tau_f] = parse_vca_model_coeffs(avg_params, MODEL);
    cell_type = man_plotgrps(i_group, 1);
    opsin_type = man_plotgrps(i_group, 4);
    brain_area = man_plotgrps(i_group, 3);
    T = [T; table(cell_type, brain_area, opsin_type, d, tau_d, f, tau_f)];
end
T_fit_to_avg_raw = T
writetable(T_fit_to_avg_raw, 'params_from_fit_to_grand_avg_raw_withRecovAndRIT_forcedRecov_newDefaults.csv')


%% STRENGTH OF INPUT ONTO INs

% define a set of attributes for each analysis group. Cell_type_pairs will
% be split according to opsin, layer, and/or brain areas
PLOTVAL = 'p1_ratio'; % 'p1_ratio', 'num_val', 'denom_val' 
cell_type_pairs = {'all_som', 'PY';...  % defines the different groups that will appear on the x-axis
                   'all_pv',  'PY';...
                   };
grp_opsins = {'any_opsin'}; % need to be encapsulated in a cell array {'any_opsin'}
grp_layers = {'L23'};
grp_brain_areas = {'AL', 'LM', 'PM', 'AM'}; % case sensitive, 'any_hva', 'med', 'lat', 'AL', 'LM' etc..

close all; clc

pop_in_strength = [];
for i_ex = 1:numel(dat)
    
    % book keeping
    pop_in_strength{i_ex}.opsin = dat{i_ex}.info.opsin;
    pop_in_strength{i_ex}.cellType = dat{i_ex}.info.cellType;
    pop_in_strength{i_ex}.brainArea = dat{i_ex}.info.brainArea;
    
    % get the data
    for i_ch = 1:2
        
        % double check that the data are present
        
        % pull out p1 amps across all trials
        p1_amps =  dat{i_ex}.qc.p1amp{i_ch};
        assert(~isempty(p1_amps), 'ERROR: could not find the p1 data')
        
        pop_in_strength{i_ex}.p1amps{i_ch} = nanmean(p1_amps);
    end
end

% initialize a group_data structure
bar_groups = [];
for i_grp = 1:size(cell_type_pairs,1)
    fld_name = sprintf('%s_and_%s', cell_type_pairs{i_grp,1}, cell_type_pairs{i_grp,2});
    for i_opsin = 1:numel(grp_opsins)
        for i_layer = 1:numel(grp_layers)
            for i_hva = 1:numel(grp_brain_areas)
                bar_groups.(fld_name).(grp_opsins{i_opsin}).(grp_layers{i_layer}).(grp_brain_areas{i_hva}).p1_ratio = [];
                bar_groups.(fld_name).(grp_opsins{i_opsin}).(grp_layers{i_layer}).(grp_brain_areas{i_hva}).num_val = [];
                bar_groups.(fld_name).(grp_opsins{i_opsin}).(grp_layers{i_layer}).(grp_brain_areas{i_hva}).denom_val = [];
            end
        end
    end
end


% fill in the group_data structure
for i_ex = 1:numel(pop_in_strength)
    
    % retrieve the brain area, and catagorize according the user defined groups
    ex_brain_area = pop_in_strength{i_ex}.brainArea;
    switch lower(ex_brain_area)
        case {'am', 'pm', 'am/pm', 'pm/am'}
            ex_brain_area = [ex_brain_area, ' med'];
        case {'al', 'lm', 'al/lm', 'lm/al'}
            ex_brain_area = [ex_brain_area, ' lat'];
    end
    l_brainarea_is_any = strcmpi(grp_brain_areas, 'any_hva');
    l_brainarea_matches_input = cellfun(@(x) ~isempty(regexp(ex_brain_area, x, 'once')), grp_brain_areas);
    l_brainarea_match = l_brainarea_is_any | l_brainarea_matches_input;
    assert(sum(l_brainarea_match)<=1, 'ERROR: >1 match to brain area');
    if sum(l_brainarea_match)==0
        continue
    end
    ex_brain_area = grp_brain_areas{l_brainarea_match};
    
    % retrieve the opsin, and catagorize according the user defined groups
    ex_opsin = pop_in_strength{i_ex}.opsin;
    l_opsin_is_any = strcmpi(grp_opsins, 'any_opsin');
    l_opsin_matches_input = cellfun(@(x) ~isempty(regexpi(ex_opsin, x)), grp_opsins);
    l_opsin_match = l_opsin_matches_input | l_opsin_is_any;
    assert(sum(l_opsin_match)==1, 'ERROR: opsin was not correctly identified')
    ex_opsin = grp_opsins{l_opsin_match};
    
    % figure out what cell & layer the data are from for each head stage
    idx_flag = cellfun(@(x) regexp(x, '_L'), pop_in_strength{i_ex}.cellType, 'uniformoutput', false);
    ex_cell_type = cellfun(@(x, y) x(1:y-1), pop_in_strength{i_ex}.cellType, idx_flag, 'uniformoutput', false);
    ex_cell_layer = cellfun(@(x, y) x(y+1:end), pop_in_strength{i_ex}.cellType, idx_flag, 'uniformoutput', false);
    
    % error checking
    all_same_layer = all(strcmpi(ex_cell_layer, ex_cell_layer{1}));
    assert(all_same_layer, 'ERROR: recordings were from different layers')
    ex_cell_layer = ex_cell_layer{1};
    
    % modify the cell types to deal with different types of entries in the
    % cell_type_pairs array
    for i_ch = 1:2
        switch lower(ex_cell_type{i_ch})
            case {'pvcre', 'fs'}
                ex_cell_type{i_ch} = [ex_cell_type{i_ch}, ' all_pv'];
            case {'somcre', 'ltsin'}
                ex_cell_type{i_ch} = [ex_cell_type{i_ch}, ' all_som'];
        end
    end
    
    % which plot group (cell_type_pair) do these data contribute to?
    cell_1_matches = cellfun(@(x) ~isempty(regexp(ex_cell_type{1}, x, 'once')), cell_type_pairs);
    cell_2_matches = cellfun(@(x) ~isempty(regexp(ex_cell_type{2}, x, 'once')), cell_type_pairs);
    l_group_match = sum(cell_1_matches | cell_2_matches, 2) == 2;
    assert(sum(l_group_match) == 1, 'ERROR: plot group has not be correctly identified')
    ex_group_name =  sprintf('%s_and_%s', cell_type_pairs{l_group_match,1}, cell_type_pairs{l_group_match,2});
    
    
    % switch the order of the cell type labels acorrding to the pattern
    % above. In so doing, figure out the P1:P1 ratio
    idx_to_cell_1 = cellfun(@(x) ~isempty(regexp(ex_cell_type{1}, x, 'once')), cell_type_pairs(l_group_match,:));
    assert(sum(idx_to_cell_1)==1);
    idx_to_cell_2 = cellfun(@(x) ~isempty(regexp(ex_cell_type{2}, x, 'once')), cell_type_pairs(l_group_match,:));
    assert(sum(idx_to_cell_2)==1);
    assert(sum(idx_to_cell_1 | idx_to_cell_2) == 2)
    
    % assign the P1:P1 ration
    p1_numerator = pop_in_strength{i_ex}.p1amps{idx_to_cell_1};
    p1_denomenator = pop_in_strength{i_ex}.p1amps{idx_to_cell_2};
    p1_ratio = p1_numerator ./ p1_denomenator;
    bar_groups.(ex_group_name).(ex_opsin).(ex_cell_layer).(ex_brain_area).p1_ratio(end+1) = p1_ratio;
    bar_groups.(ex_group_name).(ex_opsin).(ex_cell_layer).(ex_brain_area).num_val(end+1) = p1_numerator;
    bar_groups.(ex_group_name).(ex_opsin).(ex_cell_layer).(ex_brain_area).denom_val(end+1) = p1_denomenator;
    
    % do some verbose error checking
    first_cell_name = pop_in_strength{i_ex}.cellType{idx_to_cell_1};
    second_cell_name = pop_in_strength{i_ex}.cellType{idx_to_cell_2};
    original_pair = sprintf('%s_and_%s', first_cell_name, second_cell_name);
    fprintf('[%s] is [%s]\n', original_pair, ex_group_name) 

end


% now do the plotting
plot_yy = [];
plot_sem = [];
plot_labels = {};
for i_grp = 1:size(cell_type_pairs,1)
    fld_name = sprintf('%s_and_%s', cell_type_pairs{i_grp,1}, cell_type_pairs{i_grp,2});
    for i_opsin = 1:numel(grp_opsins)
        for i_layer = 1:numel(grp_layers)
            for i_hva = 1:numel(grp_brain_areas)
                tmp_data = bar_groups.(fld_name).(grp_opsins{i_opsin}).(grp_layers{i_layer}).(grp_brain_areas{i_hva}).(PLOTVAL);
                assert(~any(isnan(tmp_data)))
                if isempty(tmp_data)
                    xbar_ratios = nan;
                    sem_ratios = nan;
                else
                    xbar_ratios = mean(tmp_data);
                    sem_ratios = stderr(tmp_data);
                end
                tmp_label = sprintf('%s, %s, %s, %s, n=%d', fld_name, grp_opsins{i_opsin}, grp_layers{i_layer}, grp_brain_areas{i_hva}, numel(tmp_data));
                
                plot_yy(end+1) = xbar_ratios;
                plot_sem(end+1) = sem_ratios;
                plot_labels{end+1} = tmp_label;
            end
        end
    end
end

N_plots = numel(plot_yy);
plt_colors = lines(N_plots);
hf = figure; hold on,
for i_plt = 1:N_plots
    errorbar(i_plt, plot_yy(i_plt), plot_sem(i_plt), 's',...
                                                     'markeredgecolor', plt_colors(i_plt,:),...
                                                     'markerfacecolor', plt_colors(i_plt,:),...
                                                     'markersize', 10)
end
plot([0.8, N_plots+.2], [1,1], ':k')
h_leg = legend(plot_labels);
h_leg.Interpreter = 'none';
h_leg.Location = 'best';
legend boxoff
hf.Position = [226         619        1063         342];
set(gca, 'xtick', [1:N_plots], 'xticklabel', [])
if strcmpi(PLOTVAL, 'p1_ratio')
    set(gca, 'yscale', 'log', 'ytick', [0.0313, 0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32])
    ylim([0.0313, 32])
    xlim([0.8, N_plots+0.2])
    ylabel('P1 : P1 ratio')
end


%% INTERNEURON LATENCY ANALYSIS (DATA COLLECTION)

ENFORCE_MIN_EPSC = true;

% this analysis is only reasonable when the INs are recorded simultaneous
% to the PY cells.
assert(strcmp(EXPTTYPE, 'IN_strength'), 'ERROR: not the correct type of experiment for latency analysis')


latency_pop = [];
for i_ex = 1:numel(dat)
    
    
    % both recording channels should have valid Vclamp data. Check this,
    % and then (optionally) exclude files with very small peak epscs.
    assert(all(dat{i_ex}.info.HS_is_valid_Vclamp), 'Error: vclamp not defined for at least one channel')
    if ENFORCE_MIN_EPSC
        p1amp_avg = cellfun(@nanmean, dat{i_ex}.qc.p1amp);
        baseline_noise = cellfun(@nanmean, dat{i_ex}.qc.instNoise);
        assert(all([p1amp_avg ./ baseline_noise] >= 1.25), 'ERROR: some of the recordings are too noisy')
        if [p1amp_avg ./ baseline_noise] < 1.25
            latency_pop.dat{i_ex}.ignore{i_ch} = true;
            warning('culling data on basis of P1 amp')
        end
    end

    %
    % analyze only the first pulse latency (aggregated across all train
    % types). pull out the raw wfs. , average, re-measure the latency. This is a
    % less noisy estimate of latency than the average of single pulse
    % latencys
    %
    %%%%%%%
    latency_pop.dat{i_ex}.p1_latency_ms = {};
    latency_pop.dat{i_ex}.p1_wf = {};
    condnames = fieldnames(dat{i_ex}.expt);
    for i_ch = 1:2
        wfs = cellfun(@(x) dat{i_ex}.expt.(x).raw.snips{i_ch}(1,:,:), condnames, 'uniformoutput', false);
        wfs = cat(3, wfs{:});
        avg_p1_wf = mean(wfs, 3);
        
        N = numel(avg_p1_wf);
        tt = [0:N-1] ./ dat{i_ex}.info.sampRate.vclamp;
        tt = tt - dat{i_ex}.info.pretime.vclamp;
        peak_window = (tt > 500e-6) & (tt < 6e-3);
        minval = min(avg_p1_wf(peak_window));
        eq2min = avg_p1_wf == minval;
        peak_idx = find(eq2min & peak_window, 1, 'first');
        p1_latency_ms = tt(peak_idx) * 1000;
        
        % estimate the EPSCs rise moment
        threshold = std(avg_p1_wf(tt<0)) * -4 + mean(avg_p1_wf(tt<0));
        below_thresh = avg_p1_wf < threshold;
        below_thresh(peak_idx:end) = false;
        below_thresh(tt<0) = false;
        rise_idx = find(below_thresh, 1, 'first');
        if ~isempty(rise_idx)
            p1_latency_rise_ms = tt(rise_idx) * 1000;
        else
            p1_latency_rise_ms = nan;
        end
        
        latency_pop.dat{i_ex}.p1_latency_ms{i_ch} = p1_latency_rise_ms;
        latency_pop.dat{i_ex}.p1_wf{i_ch} = avg_p1_wf;
    end
end


%% INTERNEURON LATENCY ANALYSIS (PLOT RAW WAVEFORMS)

NORM_WF = false;

close all; clc
assert(strcmp(EXPTTYPE, 'IN_strength'), 'ERROR: not the correct type of experiment for latency analysis')

groupdata_wfs.som.py = [];
groupdata_wfs.som.in = [];
groupdata_wfs.pv.in = [];
groupdata_wfs.pv.py = [];
groupdata_wfs.pv.brainarea = {};
groupdata_wfs.som.brainarea = {};

for i_ex = 1:numel(dat)
     
    % figure out which channel had the IN, and what it is.
    ex_in_type = NaN;
    in_ch = NaN;
    for i_ch = 1:2
        cell_type = dat{i_ex}.info.cellType{i_ch};
        if ~strcmpi(cell_type, 'py_l23')
            ex_in_type = cell_type;
            in_ch = i_ch;
        end
    end
    
    % configure the field name according to IN type
    if regexpi(ex_in_type, 'LTSIN|SOMCRE')
        in_fld_name = 'som';
    elseif regexpi(ex_in_type, 'FS|PVCRE')
        in_fld_name = 'pv';
    else
        error('did not identify cell type')
    end
    
    % determine the brain area
    ex_brain_area = dat{i_ex}.info.brainArea;
    if regexpi(ex_brain_area, 'pm|am')
        ex_brain_area = 'med';
    elseif regexpi(ex_brain_area, 'lm|al')
        ex_brain_area = 'lat';
    else
        error('did not identify area')
    end
    groupdata_wfs.(in_fld_name).brainarea = cat(1, groupdata_wfs.(in_fld_name).brainarea, ex_brain_area);
    
    % assign the data to the correct field in the structure
    for i_ch = 1:2
        avg_p1_wf =  latency_pop.dat{i_ex}.p1_wf{i_ch};
        if NORM_WF
            avg_p1_wf = (avg_p1_wf ./ min(avg_p1_wf)) .* -1; % -1 to preserve sign
        end
        
        if i_ch == in_ch
            groupdata_wfs.(in_fld_name).in = cat(1, groupdata_wfs.(in_fld_name).in, avg_p1_wf);
        else
            groupdata_wfs.(in_fld_name).py = cat(1, groupdata_wfs.(in_fld_name).py, avg_p1_wf);
        end
    end
end

% figure of all the data
tt_ms = [0:numel(avg_p1_wf)-1] ./ dat{1}.info.sampRate.vclamp * 1000; 
tt_ms = tt_ms - dat{1}.info.pretime.vclamp;
hf = figure;
hf.Position = [299   430   975   420];
subplot(1,2,1), hold on,
title('PV-TOM cells', 'fontsize', 14)
plot(tt_ms, groupdata_wfs.pv.py', 'k')
plot(tt_ms, -groupdata_wfs.pv.in', 'g')
xlabel('time (ms)', 'fontsize', 14)
ylabel('current', 'fontsize', 14)
subplot(1,2,2), hold on,
title('SOM-TOM cells', 'fontsize', 14)
plot(tt_ms, groupdata_wfs.som.py', 'b')
plot(tt_ms, -groupdata_wfs.som.in', '-', 'color', [255, 140, 0]./255)
xlabel('time (ms)', 'fontsize', 14)
ylabel('current', 'fontsize', 14)

% now plot only the averages
hf = figure;
hf.Position = [299   430   975   420];
subplot(1,2,1), hold on,
title('PV-TOM cells', 'fontsize', 14)
shadedErrorBar(tt_ms, mean(groupdata_wfs.pv.py, 1), stderr(groupdata_wfs.pv.py, 1), {'k', 'linewidth', 3});
shadedErrorBar(tt_ms, -1.*mean(groupdata_wfs.pv.in, 1), stderr(-groupdata_wfs.pv.in, 1), {'g', 'linewidth', 3});
xlabel('time (ms)', 'fontsize', 14)
ylabel('current', 'fontsize', 14)
subplot(1,2,2), hold on,
title('SOM-TOM cells', 'fontsize', 14)
shadedErrorBar(tt_ms, mean(groupdata_wfs.som.py, 1), stderr(groupdata_wfs.som.py, 1), {'k', 'linewidth', 3});
shadedErrorBar(tt_ms, -1.*mean(groupdata_wfs.som.in, 1), stderr(-groupdata_wfs.som.in, 1), {'-', 'color', [255, 140, 0]./255, 'linewidth', 3});
xlabel('time (ms)', 'fontsize', 14)
ylabel('current', 'fontsize', 14)


% now plot averages, but separated by 'med', 'lat' brain areas
hf = figure;
hf.Position = [299   430   975   420];
cell_type_list = {'pv', 'som'};
for i_in = 1:numel(cell_type_list)
    for hva = {'med', 'lat'}
        l_hva = strcmpi(groupdata_wfs.(cell_type_list{i_in}).brainarea, hva{1});
        
        py_avg = mean(groupdata_wfs.(cell_type_list{i_in}).py(l_hva,:), 1);
        py_sem = stderr(groupdata_wfs.(cell_type_list{i_in}).py(l_hva,:), 1);
        in_avg = -1 .* mean(groupdata_wfs.(cell_type_list{i_in}).in(l_hva,:), 1);
        in_sem = stderr(-groupdata_wfs.(cell_type_list{i_in}).in(l_hva,:), 1);
        
        % make the plot
        switch hva{1}
            case 'med'
                lineclr = 'b';
            case 'lat'
                lineclr = 'r';
        end
        subplot(1,2,i_in), hold on,
        shadedErrorBar(tt_ms, py_avg, py_sem, {'-', 'color', lineclr});
        shadedErrorBar(tt_ms, in_avg, in_sem, {'-', 'color', lineclr});
        xlabel('time (ms)', 'fontsize', 14)
        ylabel('current', 'fontsize', 14)
        title(sprintf('%s-tom cells', cell_type_list{i_in}), 'fontsize', 14)
    end
end
%% INTERNEURON LATENCY ANALYSIS (PLOTS)


close all; clc


COMBINE_HVAS = true;
COMBINE_OPSINS = true;
COMBINE_INS = true;

assert(strcmp(EXPTTYPE, 'IN_strength'), 'ERROR: not the correct type of experiment for latency analysis')
assert(COMBINE_OPSINS, 'error: combine opsins needs to be true, otherwise things may break')
assert(COMBINE_INS, 'error: combine interneurons needs to be true, otherwise things may break')

% make a table with the following fields
% {BrainArea, OpsinType, IN_type, PY_latency, IN_latency}
groupdata = {};
row_names = {'BrainArea', 'OpsinType', 'IN_type', 'PY_latency', 'IN_latency'};
template = {'str', 'str', 'str', NaN, NaN};

for i_ex = 1:numel(dat)
    
    ex_template = template;
    ex_brain_area = dat{i_ex}.info.brainArea;
    ex_opsin = dat{i_ex}.info.opsin;
    
    if COMBINE_HVAS
        if regexpi(ex_brain_area, 'pm|am')
            ex_brain_area = 'med';
        elseif regexpi(ex_brain_area, 'lm|al')
            ex_brain_area = 'lat';
        else
            error('did not identify area')
        end
    end
    ex_template{1} = ex_brain_area;
    
    if COMBINE_OPSINS
        if regexpi(ex_opsin, 'chronos')
            ex_opsin = 'chronos';
        elseif regexpi(ex_opsin, 'chief')
            ex_opsin = 'chief';
        else
            error('did not identify opsin')
        end
    end
    ex_template{2} = ex_opsin;
    
    
    
    
    for i_ch = 1:2
        
        cell_type = dat{i_ex}.info.cellType{i_ch};
        if ~strcmpi(cell_type, 'py_l23')
            if COMBINE_INS
                if regexpi(cell_type, 'LTSIN|SOMCRE')
                    cell_type = 'all_som';
                elseif regexpi(cell_type, 'FS|PVCRE')
                    cell_type = 'all_pv';
                else
                    error('did not identify cell type')
                end
            end
            ex_template{3} = cell_type;
            latency_idx = strcmp(row_names, 'IN_latency');
        else
            latency_idx = strcmp(row_names, 'PY_latency');
        end
        
        ex_template{latency_idx} = latency_pop.dat{i_ex}.p1_latency_ms{i_ch};
        
    end
    
    groupdata(i_ex,:) = ex_template;
end

t_latency = cell2table(groupdata, 'VariableNames', row_names);
t_latency.diffval = minus(t_latency.IN_latency, t_latency.PY_latency);

% summary table:
t_xbar = grpstats(t_latency,...
                  {'BrainArea','IN_type'},...
                  {'nanmean'},...
                  'DataVars', {'PY_latency', 'IN_latency', 'diffval'},...
                  'VarNames', {'BrainArea','IN_type','GroupCount','PY_mean','IN_mean', 'diff_mean'})
              
t_sem = grpstats(t_latency,...
                  {'BrainArea','IN_type'},...
                  {'stderr'},...
                  'DataVars', {'PY_latency', 'IN_latency', 'diffval'},...
                  'VarNames', {'BrainArea','IN_type','GroupCount','PY_sem','IN_sem', 'diff_sem'});

unique_cell_types = unique(t_latency.IN_type);
unique_areas = unique(t_latency.BrainArea);
[py_xbar, in_xbar, diff_vals] = deal([]);
[py_sem, in_sem, diff_sem] = deal([]);
for i_area = 1:numel(unique_areas)
    for i_in = 1:numel(unique_cell_types)
        l_area = strcmpi(t_latency.BrainArea, unique_areas{i_area});
        l_in = strcmpi(t_latency.IN_type, unique_cell_types{i_in});
        
        t_grp = t_latency(l_area & l_in, {'PY_latency', 'IN_latency', 'diffval'});
        py_xbar(i_area, i_in) = nanmean(t_grp.PY_latency);
        in_xbar(i_area, i_in) = nanmean(t_grp.IN_latency);
        diff_vals(i_area, i_in) = nanmean(t_grp.diffval);
        
        py_sem(i_area, i_in) = stderr(t_grp.PY_latency);
        in_sem(i_area, i_in) = stderr(t_grp.IN_latency);
        diff_sem(i_area, i_in) = stderr(t_grp.diffval);
    end
end

in_pltclrs = {'g', [255, 140, 0]./255};

hf = figure;
hf.Position = [402         216        1289         368];
subplot(1,3,1)
hold on,
xmax = numel(unique_areas);
for i_in = 1:numel(unique_cell_types)
    errorbar(1:numel(unique_areas), diff_vals(:,i_in), diff_sem(:,i_in),...
        'linewidth', 3, 'color', in_pltclrs{i_in})
end
plot([0.75, xmax+.25], [0 0], '--k')
xlim([0.75, xmax+.25])
hl = legend(unique_cell_types);
hl.Interpreter = 'None';
hl.Location = 'Best';
hl.Box = 'off';
title('diff val: IN - PY')
ylabel('latency (ms)')
xlabel('brain area')
set(gca, 'tickDir', 'out', 'xtick', [1:xmax], 'xticklabel', unique_areas)

subplot(1,3,2)
hold on,
for i_in = 1:numel(unique_cell_types)
    errorbar(1:numel(unique_areas), py_xbar(:,i_in), py_sem(:,i_in),...
        'linewidth', 3, 'color', in_pltclrs{i_in})
end
% ylim([3e-3, 5.8e-3])
xlim([0.75, xmax+.25])
hl = legend(unique_cell_types);
hl.Interpreter = 'None';
hl.Location = 'Best';
hl.Box = 'off';
title('raw latency PY only')
ylabel('latency (ms)')
xlabel('brain area')
set(gca, 'tickDir', 'out', 'xtick', [1:xmax], 'xticklabel', unique_areas)


subplot(1,3,3)
hold on,
for i_in = 1:numel(unique_cell_types)
    errorbar(1:numel(unique_areas), in_xbar(:,i_in), in_sem(:,i_in),...
        'linewidth', 3, 'color', in_pltclrs{i_in})
end
% ylim([3e-3, 5.8e-3])
xlim([0.75, xmax+.25])
hl = legend(unique_cell_types);
hl.Interpreter = 'None';
hl.Location = 'Best';
hl.Box = 'off';
title('raw latency IN only')
ylabel('latency (ms)')
xlabel('brain area')
set(gca, 'tickDir', 'out', 'xtick', [1:xmax], 'xticklabel', unique_areas)



%% INTERNEURON AMPLITUDE ANALYSIS (DATA COLLECTION)

ENFORCE_MIN_EPSC = false;

% this analysis is only reasonable when the INs are recorded simultaneous
% to the PY cells.
assert(strcmp(EXPTTYPE, 'IN_strength'), 'ERROR: not the correct type of experiment for latency analysis')


inamp_pop = [];
inamp_pop.TFsAllExpts = [];
inamp_pop.recoveryTimesAllExpts = [];
for i_ex = 1:numel(dat)
    
    
    % both recording channels should have valid Vclamp data. Check this,
    % and then (optionally) exclude files with very small peak epscs.
    assert(all(dat{i_ex}.info.HS_is_valid_Vclamp), 'Error: vclamp not defined for at least one channel')
    if ENFORCE_MIN_EPSC
        p1amp_avg = cellfun(@nanmean, dat{i_ex}.qc.p1amp);
        baseline_noise = cellfun(@nanmean, dat{i_ex}.qc.instNoise);
        assert(all([p1amp_avg ./ baseline_noise] >= 1.25), 'ERROR: some of the recordings are too noisy')
        if [p1amp_avg ./ baseline_noise] < 1.25
            inamp_pop.dat{i_ex}.ignore{i_ch} = true;
            warning('culling data on basis of P1 amp')
        end
    end
    
    % I'm including the multi-power datasets. I need to flag them and
    % select out a single amplitude to analyze
    condnames = fieldnames(dat{i_ex}.expt);
    powers = cellfun(@(x) dat{i_ex}.expt.(x).tdict(1), condnames);
    n_powers = numel(unique(powers));
    if n_powers > 1
        err_from_2pt7 = powers - 2.7;
        err_from_2pt7(err_from_2pt7 < 0) = Inf;
        assert(~all(isinf(err_from_2pt7)), 'ERROR: all the Voltages are below 2.7')
        idx_to_2pt7 = find(err_from_2pt7 == min(err_from_2pt7));
        
        % keep only the condition where the laser power is closest to 2.7V
        condnames = condnames{idx_to_2pt7};
    end

    %
    % analyze only the first pulse amplitude (aggregated across all train
    % types). pull out the raw wfs. , average, re-measure the amplitude. This is a
    % less noisy estimate of latency than the average of single pulse
    % latencys
    %
    %%%%%%%
    inamp_pop.dat{i_ex}.p1_amp = {};
    for i_ch = 1:2
        if n_powers > 1
            avg_p1_wf = dat{i_ex}.expt.(condnames).raw.snips{i_ch};
        else
            wfs = cellfun(@(x) dat{i_ex}.expt.(x).raw.snips{i_ch}(1,:,:), condnames, 'uniformoutput', false);
            wfs = cat(3, wfs{:});
            avg_p1_wf = mean(wfs, 3);
        end
        
        N = numel(avg_p1_wf);
        tt = [0:N-1] ./ dat{i_ex}.info.sampRate.vclamp;
        tt = tt - dat{i_ex}.info.pretime.vclamp;
        peak_window = (tt > 500e-6) & (tt < 6e-3);
        p1_amp = min(avg_p1_wf(peak_window));
        
        inamp_pop.dat{i_ex}.p1_amp{i_ch} = abs(p1_amp);
    end
end

%% INTERNEURON AMPLITUDE ANALYSIS (PLOTS)

close all; clc


COMBINE_HVAS = true;
COMBINE_OPSINS = true;
COMBINE_INS = true;

assert(strcmp(EXPTTYPE, 'IN_strength'), 'ERROR: not the correct type of experiment for latency analysis')
assert(COMBINE_OPSINS, 'error: combine opsins needs to be true, otherwise things may break')
assert(COMBINE_INS, 'error: combine interneurons needs to be true, otherwise things may break')

% make a table with the following fields
% {BrainArea, OpsinType, IN_type, PY_amp, IN_amp}
groupdata = {};
row_names = {'BrainArea', 'OpsinType', 'IN_type', 'PY_amp', 'IN_amp'};
template = {'str', 'str', 'str', NaN, NaN};

for i_ex = 1:numel(dat)
    
    ex_template = template;
    ex_brain_area = dat{i_ex}.info.brainArea;
    ex_opsin = dat{i_ex}.info.opsin;
    
    if COMBINE_HVAS
        if regexpi(ex_brain_area, 'pm|am')
            ex_brain_area = 'med';
        elseif regexpi(ex_brain_area, 'lm|al')
            ex_brain_area = 'lat';
        else
            error('did not identify area')
        end
    end
    ex_template{1} = ex_brain_area;
    
    if COMBINE_OPSINS
        if regexpi(ex_opsin, 'chronos')
            ex_opsin = 'chronos';
        elseif regexpi(ex_opsin, 'chief')
            ex_opsin = 'chief';
        else
            error('did not identify opsin')
        end
    end
    ex_template{2} = ex_opsin;
    
    
    
    
    for i_ch = 1:2
        
        cell_type = dat{i_ex}.info.cellType{i_ch};
        if ~strcmpi(cell_type, 'py_l23')
            if COMBINE_INS
                if regexpi(cell_type, 'LTSIN|SOMCRE')
                    cell_type = 'all_som';
                elseif regexpi(cell_type, 'FS|PVCRE')
                    cell_type = 'all_pv';
                else
                    error('did not identify cell type')
                end
            end
            ex_template{3} = cell_type;
            amp_idx = strcmp(row_names, 'IN_amp');
        else
            amp_idx = strcmp(row_names, 'PY_amp');
        end
        
        ex_template{amp_idx} = inamp_pop.dat{i_ex}.p1_amp{i_ch};
        
    end
    
    groupdata(i_ex,:) = ex_template;
end

t_amp = cell2table(groupdata, 'VariableNames', row_names);
t_amp.ratioval = t_amp.IN_amp ./  t_amp.PY_amp;

% summary table:
t_xbar = grpstats(t_amp,...
                  {'BrainArea','IN_type'},...
                  {'mean'},...
                  'DataVars', {'PY_amp', 'IN_amp', 'ratioval'},...
                  'VarNames', {'BrainArea','IN_type','GroupCount','PY_mean','IN_mean', 'ratio_mean'})
              
t_sem = grpstats(t_amp,...
                  {'BrainArea','IN_type'},...
                  {'stderr'},...
                  'DataVars', {'PY_amp', 'IN_amp', 'ratioval'},...
                  'VarNames', {'BrainArea','IN_type','GroupCount','PY_sem','IN_sem', 'ratio_sem'});

unique_cell_types = unique(t_amp.IN_type);
unique_areas = unique(t_amp.BrainArea);
[py_xbar, in_xbar, ratio_xbar] = deal([]);
[py_sem, in_sem, ratio_sem] = deal([]);
all_ratio_vals = {};
for i_area = 1:numel(unique_areas)
    for i_in = 1:numel(unique_cell_types)
        l_area = strcmpi(t_amp.BrainArea, unique_areas{i_area});
        l_in = strcmpi(t_amp.IN_type, unique_cell_types{i_in});
        
        t_grp = t_amp(l_area & l_in, {'PY_amp', 'IN_amp', 'ratioval'});
        py_xbar(i_area, i_in) = nanmean(t_grp.PY_amp);
        in_xbar(i_area, i_in) = nanmean(t_grp.IN_amp);
        ratio_xbar(i_area, i_in) = nanmean(t_grp.ratioval);
        all_ratio_vals{i_area, i_in} = t_grp.ratioval;
        
        py_sem(i_area, i_in) = stderr(t_grp.PY_amp);
        in_sem(i_area, i_in) = stderr(t_grp.IN_amp);
        ratio_sem(i_area, i_in) = stderr(t_grp.ratioval);
    end
end

pltclrs.all_pv = [0, 100, 0]./255;
pltclrs.all_som = [255, 140, 0]./255;

% boxplot of ratio values to see the full distribution
hf = figure;
hf.Position = [1001         388         738         368];
for i_in = 1:numel(unique_cell_types)
    subplot(1,numel(unique_cell_types),i_in), hold on
    %pclr = pltclrs.(unique_cell_types{i_in});
    pclr='rbcy';
    scatter_clr = {'r', 'b', 'c', 'y'};
    box_x = [];
    box_g = [];
    for i_hva = 1:numel(unique_areas)
        box_x = cat(1, box_x, all_ratio_vals{i_hva, i_in});
        box_g = cat(1, box_g, repmat(i_hva, numel(all_ratio_vals{i_hva, i_in}), 1));
    end
    boxplot(box_x, box_g, 'colors', pclr, 'symbol', '+');
    for hc1 = get(gca, 'children')'
        for hc2 = hc1.Children'
            if strcmpi(hc2.Type, 'line')
                hc2.LineWidth = 2;
            end
        end
    end
    scatter(box_g, box_x, 'markeredgecolor', 'k', 'linewidth', 1)
    title(sprintf('IN:PY for %s', unique_cell_types{i_in}))
    ylabel('amplitude ratio', 'fontsize', 14)
    xlabel('brain area', 'fontsize', 14)
    set(gca, 'tickDir', 'out', 'xticklabel', unique_areas, 'box', 'off')
    set(gca, 'yscale', 'linear')
end



hf = figure;
hf.Position = [402         216        1289         368];
xmax = numel(unique_areas);
for i_in = 1:numel(unique_cell_types)
    subplot(1,numel(unique_cell_types),i_in), hold on
    errorbar(1:numel(unique_areas), in_xbar(:,i_in), in_sem(:,i_in),...
        'linewidth', 3)
    plot([0.75, xmax+.25], [0 0], '--k')
    
    xlim([0.75, xmax+.25])
    hl = legend(unique_cell_types{i_in});
    hl.Interpreter = 'None';
    hl.Location = 'Best';
    hl.Box = 'off';
    title('EPSCs onto INs')
    ylabel('EPSC (pA)')
    xlabel('brain area')
    set(gca, 'tickDir', 'out', 'xtick', [1:xmax], 'xticklabel', unique_areas)
end


hf = figure; hold on,
hf.Position = [402         216        1289         368];
xmax = numel(unique_areas);
for i_in = 1:numel(unique_cell_types)
    errorbar(1:numel(unique_areas), py_xbar(:,i_in), py_sem(:,i_in),...
        'linewidth', 3)
end
plot([0.75, xmax+.25], [0 0], '--k')
xlim([0.75, xmax+.25])
hl = legend(unique_cell_types);
hl.Interpreter = 'None';
hl.Location = 'Best';
hl.Box = 'off';
title('EPSCs onto PY cells')
ylabel('EPSC (pA)')
xlabel('brain area')
set(gca, 'tickDir', 'out', 'xtick', [1:xmax], 'xticklabel', unique_areas)

% PY and IN data separately on YY axis
PLOT_NORM_TO_HVA1 = true;
hf = figure; hold on,
hf.Position = [304   335   777   368];
xmax = numel(unique_areas);
for i_in = 1:numel(unique_cell_types)
    subplot(1,numel(unique_cell_types),i_in), hold on,
    title('EPSCs normalized to HVA1')
    
    tmp_py_xbar =  py_xbar(:,i_in);
    tmp_in_xbar = in_xbar(:,i_in);
    x = 1:numel(unique_areas);
    
    if PLOT_NORM_TO_HVA1
        
        prcnt_change_py = tmp_py_xbar./tmp_py_xbar(1);
        prcnt_change_in = tmp_in_xbar./tmp_in_xbar(1);
        plot(x, prcnt_change_py, 'k', 'linewidth', 3);
        plot(x, prcnt_change_in, 'color', pltclrs.(unique_cell_types{i_in}), 'linewidth', 3);
        
        
        hax = gca;
        hax.YScale = 'log';
        hax.YLim = [0.5, 2];
        hax.Box = 'off';
        hax.YTick = [0.5, 0.7, 1, 1.4, 2];
        hax.FontSize = 14;
        hax.LineWidth = 1;
        hax.XTick = 1:numel(unique_areas);
        hax.XTickLabels = unique_areas;
        hax.XLim = [0.8, numel(unique_areas)+0.2];
    else
        [hax, hxbar_l, hxbar_r] = plotyy(x, tmp_py_xbar, x, tmp_in_xbar);
        hxbar_l.Color = 'k';
        hxbar_l.LineWidth = 3;
        hxbar_r.Color = pltclrs.(unique_cell_types{i_in});
        hxbar_r.LineWidth = 3;
        set(hax, 'box', 'off',...
            'fontsize', 14,...
            'linewidth', 1,...
            'xtick', 1:numel(unique_areas),...
            'xticklabels', unique_areas,...
            'xlim', [0.8, numel(unique_areas)+0.2])
        hax(1).YColor = 'k';
        hax(2).YColor = pltclrs.(unique_cell_types{i_in});
        
        tmp_py_sem =  py_sem(:,i_in);
        tmp_in_sem = in_sem(:,i_in);
        
    end
end


%% INTERNEURON MULTIPOWER (DATA COLLECTION)

% this analysis is only reasonable when the INs are recorded simultaneous
% to the PY cells.
assert(strcmp(EXPTTYPE, 'IN_powers'), 'ERROR: not the correct type of experiment for latency analysis')


inpow_pop = [];
for i_ex = 1:numel(dat)
    
    
    % both recording channels should have valid Vclamp data. Check this,
    % and then (optionally) exclude files with very small peak epscs.
    assert(all(dat{i_ex}.info.HS_is_valid_Vclamp), 'Error: vclamp not defined for at least one channel')

    %
    % analyze the first pulse amplitude across powers
    %
    %%%%%%%
    inpow_pop.dat{i_ex}.p1_amp = {};
    inpow_pop.dat{i_ex}.laser_V = [];
    condnames = fieldnames(dat{i_ex}.expt);
    ex_laser_powers = cellfun(@(x) dat{i_ex}.expt.(x).tdict(1), condnames);
    l_ex_good_powers = true(size(ex_laser_powers));
    ch_amps = {};
    for i_ch = 1:2
        amps = cellfun(@(x) dat{i_ex}.expt.(x).stats.EPSCamp{i_ch}(1,:,:), condnames, 'uniformoutput', false);
        % hack to remove some sweeps that were deleted in the .abf file
        % post-hoc. Greedy update of exp-wide list
        l_ch_good_powers = cellfun(@(x) ~isempty(x), amps);
        l_ex_good_powers = l_ex_good_powers & l_ch_good_powers;
        ch_amps{i_ch} = amps;
    end
    for i_ch = 1:2
        inpow_pop.dat{i_ex}.p1_amp{i_ch} = cat(1, ch_amps{i_ch}{l_ex_good_powers});
    end
    inpow_pop.dat{i_ex}.laser_V = ex_laser_powers(l_ex_good_powers);
    
    above_10pa = cellfun(@(x) x>10, inpow_pop.dat{i_ex}.p1_amp, 'uniformoutput', false);
    l_gt_10pA = above_10pa{1} & above_10pa{2};
    inpow_pop.dat{i_ex}.p1_amp = cellfun(@(x) x(l_gt_10pA), inpow_pop.dat{i_ex}.p1_amp, 'uniformoutput', false);
    inpow_pop.dat{i_ex}.laser_V = inpow_pop.dat{i_ex}.laser_V(l_gt_10pA);
    
    
    % grab some meta data
    inpow_pop.dat{i_ex}.hva = dat{i_ex}.info.brainArea;
    inpow_pop.dat{i_ex}.cell_type = dat{i_ex}.info.cellType;
end


%% INTERNEURON MULTIPOWER (PLOTS)



PLOT_TYPE = 'ratio'; % could be 'ratio', or 'scatter'
COLOR_BY = 'in';  % could be 'hva', or 'in'

NORM_TO_PY_MAX = false;
COMBINE_HVAS = true;
COMBINE_INS = true;

assert(strcmp(EXPTTYPE, 'IN_powers'), 'ERROR: not the correct type of experiment for multipower analysis')
assert(COMBINE_INS, 'error: combine interneurons needs to be true, otherwise things may break')

figure, hold on,

for i_ex = 1:numel(dat)
    
    tmp_dat = {[], []};
    for i_ch = 1:2
        
        ex_brain_area = inpow_pop.dat{i_ex}.hva;
        if COMBINE_HVAS
            if regexpi(ex_brain_area, 'pm|am')
                ex_brain_area = 'med';
                hva_clr = 'b';
            elseif regexpi(ex_brain_area, 'lm|al')
                ex_brain_area = 'lat';
                hva_clr = 'r';
            else
                error('did not identify area')
            end
        end
        
        cell_type = inpow_pop.dat{i_ex}.cell_type{i_ch};
        if ~strcmpi(cell_type, 'py_l23')
            if COMBINE_INS
                if regexpi(cell_type, 'LTSIN|SOMCRE')
                    cell_type = 'all_som';
                    cell_clr = [250, 100, 0]./255;
                elseif regexpi(cell_type, 'FS|PVCRE')
                    cell_type = 'all_pv';
                    cell_clr = [0, 128, 0]./255;
                else
                    error('did not identify cell type')
                end
            end
            plt_idx = 2; % interneurons are on Y
        else
            plt_idx = 1;  % PY cell on X
        end
       
        tmp_dat{plt_idx} = inpow_pop.dat{i_ex}.p1_amp{i_ch};
    end
    
    if NORM_TO_PY_MAX
        maxval = tmp_dat{1}(end);
        tmp_dat = cellfun(@(x) x./maxval, tmp_dat, 'uniformoutput', false);
    end
    
    
    switch COLOR_BY
        case 'in'
            plt_clr = cell_clr;
        case 'hva'
            plt_clr = hva_clr;
    end
    
    % to plot the raw values
    if strcmpi(PLOT_TYPE, 'scatter')
        plot(tmp_dat{1}, tmp_dat{2}, '-o', 'color', plt_clr, 'linewidth', 2)
    elseif strcmpi(PLOT_TYPE, 'ratio')
        % to plot the ratio values
        ratio_vals = tmp_dat{2} ./ tmp_dat{1};
        plot(inpow_pop.dat{i_ex}.laser_V, ratio_vals, '-o', 'color', plt_clr, 'linewidth', 2)
        set(gca, 'yscale', 'log')
        ylabel('IN : PY ratio', 'fontsize', 14)
        xlabel('Laser Voltage', 'fontsize', 14)
    end
end


%% DISTANCE BETWEEN PY AND IN FOR PAIRED RECORDINGS

xy_diff = nan(numel(dat), 2);
for i_ex = 1:numel(dat)
    % both recording channels should have valid Vclamp data. Check this,
    % and then (optionally) exclude files with very small peak epscs.
    assert(all(dat{i_ex}.info.HS_is_valid_Vclamp), 'Error: vclamp not defined for at least one channel')
    
    % find indicies to PY cell and IN
    idx_py = strcmpi(dat{i_ex}.info.cellType, 'py_l23');
    assert(sum(idx_py) == 1, 'ERROR: did not find the PY cell')
    
    % pull out xy cordinates for IN and PY cell
    xy_py = dat{i_ex}.info.HS_xy_pos{idx_py};
    xy_in = dat{i_ex}.info.HS_xy_pos{~idx_py};
    
    xy_diff(i_ex,:) = xy_py - xy_in;
end

figure
[theta, radius] = cart2pol(xy_diff(:,1), xy_diff(:,2));
hp = polar(theta, log10(radius), 'ko');
hp.MarkerFaceColor = 'k';
ytick = get(gca, 'ytick');
set(gca, 'yticklabel', cellfun(@(x) num2str(10^x, 3), num2cell(ytick), 'uniformoutput', false))

figure
histogram(radius, 15)
xlabel('distance between PY and IN (um)')
ylabel('count')

%% ZUCKER POPULATION ANALYSIS (DATA COLLECTION) WITH TABLES

%  loop over all the datasets. compile:
% celltype, layer, opsin, brainarea, train_freq, p_train_times, p_zucker_times, train_amps, zucker_amps
%
% all times should be relative to the first pulse time. 

[pop_celltype, pop_opsin, pop_brainarea, pop_train_freq, pop_train_times, pop_zucker_times, pop_train_amps, pop_zucker_amps] = deal({});
for i_ex = 1:numel(dat)
    
    for i_ch = 1:2
        
        isvalid = dat{i_ex}.info.HS_is_valid_Vclamp(i_ch);
        if ~isvalid; continue; end
        
        % define the easy things first
        pop_celltype{end+1} = dat{i_ex}.info.cellType{i_ch};
        pop_opsin{end+1} = dat{i_ex}.info.opsin;
        pop_brainarea{end+1} = dat{i_ex}.info.brainArea;
        
        % figure out the train frequency
        condnames = fieldnames(dat{i_ex}.expt);
        trainParams = cellfun(@(x) dat{i_ex}.expt.(x).tdict, condnames, 'uniformoutput', false);
        trainParams = cat(1, trainParams{:});
        assert(numel(unique(trainParams(:,3))) == 1, 'ERROR: too many TFs')
        pop_train_freq{end+1} = unique(trainParams(:,3));
        
        [ch_train_amps, ch_recov_amps, ch_train_times,ch_recov_times] = deal([]);
        % loop through the ex conds and pull out the amps, ptimes
        for cond = condnames'
            cond_ptimes = dat{i_ex}.expt.(cond{1}).pOnTimes;
            cond_ptimes = cond_ptimes - cond_ptimes(1);
            ipis = diff(cond_ptimes);
            ipis = [ipis(1), ipis]; % add back the one that 'diff' stole
            ipis = round(ipis, 6);
            unique_ipis = unique(ipis, 'stable');
            l_train = ipis == unique_ipis(1);
            
            % pull out the amps, norm to the smooth p1amp
            epsc_amps = dat{i_ex}.expt.(cond{1}).stats.EPSCamp{i_ch};
            real_trl_nums = dat{i_ex}.expt.(cond{1}).realTrialNum{i_ch};
            norm_vals = dat{i_ex}.qc.p1amp_norm{i_ch}(real_trl_nums);
            norm_vals = permute(norm_vals, [3,1,2]);
            epsc_amps = bsxfun(@rdivide, epsc_amps, norm_vals);
            epsc_amps = mean(epsc_amps, 3);
            
            % partion everything to the train or recovery groups
            for i_pulse = 1:numel(epsc_amps)
                if l_train(i_pulse)
                    ch_train_amps(end+1) = epsc_amps(i_pulse);
                    ch_train_times(end+1) = cond_ptimes(i_pulse);
                    
                else
                    ch_recov_amps(end+1) = epsc_amps(i_pulse);
                    ch_recov_times(end+1) = cond_ptimes(i_pulse);
                end
            end
        end
        
        % average across duplicates of the same pulse time for the recov
        % and train conditions
        [avg_train_amps, avg_recov_amps, avg_train_times, avg_recov_times] = deal([]);
        unique_train_times = unique(ch_train_times);
        for i_time = unique_train_times
            idx = ch_train_times == i_time;
            avg_train_amps(end+1) = mean(ch_train_amps(idx));
            avg_train_times(end+1) = i_time;
        end
        
        unique_recov_times = unique(ch_recov_times);
        for i_time = unique_recov_times
            idx = ch_recov_times == i_time;
            avg_recov_amps(end+1) = mean(ch_recov_amps(idx));
            avg_recov_times(end+1) = i_time;
        end
        
        % assign the ch data to the pop data
        pop_train_times{end+1} = avg_train_times;
        pop_zucker_times{end+1} = avg_recov_times;
        pop_train_amps{end+1} = avg_train_amps;
        pop_zucker_amps{end+1} = avg_recov_amps;
        
    end
end

% make a table for aggregation functions down the line...


% ------ make a crappy summary figure ----- %

l_plt = [1,2,3,4]

p_times_induction = cat(1, pop_train_times{l_plt});
p_times_induction = unique(p_times_induction, 'rows');

p_times_zucker = cat(1, pop_zucker_times{l_plt});
p_times_zucker = unique(p_times_zucker, 'rows');

train_amps = cat(1, pop_train_amps{l_plt});
zucker_amps = cat(1, pop_zucker_amps{l_plt});

figure, hold on,
plot(p_times_induction, train_amps', 'linewidth', 2)
ax = gca;
ax.ColorOrderIndex = 1;
plot(p_times_zucker, zucker_amps', 'linewidth', 2)



%% CONTROLS: LOOKING FOR RECURRENT OR POLY-SYNAPTIC PSCs

% logic, longer 90-10 times correspond to recurrent activity. Monosynaptic
% events should be ~2 ms to peak, recurrent activity might be ~4. And the
% 90-10 time might be related to the EPSC magnitude.

CELL_TYPE = {'PY', 'L23', 'any', 'any'};  % for use with groupMatcher.m {'celltype', 'layer', 'brainarea', 'opsin'} 

opsin_types = {'chief', 'chronos'};
opsin_counters = [0 0];
p_counter = 1;
opsin_clrs = {'r', 'b'};
f = figure;
subplot(1,2,1); hold on,
subplot(1,2,2); hold on,
for i_ex = 1:numel(dat)
      for i_ch = 1:2
        
        % skip past recording channels that have no data
        isvalid = dat{i_ex}.info.HS_is_valid_Vclamp(i_ch);
        if ~isvalid
            continue
        end
        
        % skip recordings that are not the correct cell type
        ch_attribs = {dat{i_ex}.info.cellType{i_ch}, upper(dat{i_ex}.info.brainArea), dat{i_ex}.info.opsin};
        is_correct_cell_type = groupMatcher(CELL_TYPE, ch_attribs);
        if ~is_correct_cell_type
            continue
        end
        
        conds = fieldnames(dat{i_ex}.expt);
        ex_risetime = [];
        for i_cond = 1:numel(conds)
            cond_risetime = mean(dat{i_ex}.expt.(conds{i_cond}).stats.rise_time{i_ch}(1,1,:), 3);
            ex_risetime = cat(1, ex_risetime, cond_risetime);
        end
        
        % add to the correct plot
        ex_risetime = ex_risetime .* 1000; % now in ms
        opsin_idx = cellfun(@(x) any(regexp(dat{i_ex}.info.opsin, x)), opsin_types);
        assert(sum(opsin_idx)==1, 'ERROR: incorect opsin type')
        opsin_counters = opsin_counters + opsin_idx;
        
        % plot
        subplot(1,2,find(opsin_idx))
        plot(opsin_counters(opsin_idx), ex_risetime, '.', 'color', opsin_clrs{opsin_idx})
        p(p_counter) = plot(opsin_counters(opsin_idx), mean(ex_risetime), 's', 'color', opsin_clrs{opsin_idx}, 'markerfacecolor', opsin_clrs{opsin_idx});
        fsep_idx = find(dat{i_ex}.info.fid.vclamp==filesep, 1, 'last');
        p(p_counter).UserData = sprintf('%s, ch: %d', dat{i_ex}.info.fid.vclamp(fsep_idx+1 : end-4), i_ch);
        p(p_counter).ButtonDownFcn = @(a,~) disp(a.UserData);
        p_counter = p_counter+1;
        axis tight
        xlabel('Cell Number')
        ylabel('EPSC rise time (ms)')
      end
end

ymax = 0;
for i_plt = 1:2;
    subplot(1,2,i_plt)
    vals = get(gca, 'ylim');
    ymax = max([ymax, vals(2)]);
end
for i_plt = 1:2;
    subplot(1,2,i_plt)
    set(gca, 'ylim', [0, ymax]);
end


figure,
subplot(1,2,1); hold on,
subplot(1,2,2); hold on,
for i_ex = 1:numel(dat)
      for i_ch = 1:2
        
        % skip past recording channels that have no data
        isvalid = dat{i_ex}.info.HS_is_valid_Vclamp(i_ch);
        if ~isvalid
            continue
        end
        
        % skip recordings that are not the correct cell type
        ch_attribs = {dat{i_ex}.info.cellType{i_ch}, upper(dat{i_ex}.info.brainArea), dat{i_ex}.info.opsin};
        is_correct_cell_type = groupMatcher(CELL_TYPE, ch_attribs);
        if ~is_correct_cell_type
            continue
        end
        
        conds = fieldnames(dat{i_ex}.expt);
        ex_risetime = [];
        ex_amps = [];
        for i_cond = 1:numel(conds)
            if isempty(dat{i_ex}.expt.(conds{i_cond}).stats.EPSCamp{i_ch})
                continue
            end
            cond_risetime = mean(dat{i_ex}.expt.(conds{i_cond}).stats.rise_time{i_ch}(1,1,:), 3);
            ex_risetime = cat(1, ex_risetime, cond_risetime);
            
            cond_amps =  mean(dat{i_ex}.expt.(conds{i_cond}).stats.EPSCamp{i_ch}(1,1,:), 3);
            ex_amps = cat(1, ex_amps, cond_amps);
        end
        
        % add to the correct plot
        ex_risetime = ex_risetime .* 1000; % now in ms
        opsin_idx = cellfun(@(x) any(regexp(dat{i_ex}.info.opsin, x)), opsin_types);
        assert(sum(opsin_idx)==1, 'ERROR: incorect opsin type')
        
        % plot
        subplot(1,2,find(opsin_idx))
        plot(ex_amps, ex_risetime, '.', 'color', opsin_clrs{opsin_idx})
        plot(mean(ex_amps), mean(ex_risetime), 's', 'color', opsin_clrs{opsin_idx}, 'markerfacecolor', opsin_clrs{opsin_idx})
        xlabel('EPSC amplitude')
        ylabel('EPSC latency to peak')
      end
end

ymax = 0;
for i_plt = 1:2;
    subplot(1,2,i_plt)
    vals = get(gca, 'ylim');
    ymax = max([ymax, vals(2)]);
end
for i_plt = 1:2;
    subplot(1,2,i_plt)
    set(gca, 'ylim', [0, ymax]);
end

%% CONTROLS: VARIABILITY BY OPSIN FOR 1'ST PULSE

close all
ENFORCE_ABOVE_100PA = true;

%
% Plot the trial to trial P1 amplitude expressed as a fractional change
% from the smooth running average
%
figure, hold on,
chronos_std = [];
chief_std = [];
for i_ex = 1:numel(dat)
    
    for i_ch = 1:2
        
        % check to make sure this neuron was defined
        has_Vclamp_data = ~strcmpi(dat{i_ex}.info.fid.vclamp(end-7:end), 'none.abf');
        isvalid = dat{i_ex}.info.HS_is_valid_Vclamp(i_ch);
        if ~isvalid || ~has_Vclamp_data 
            continue
        end
        
        tmp = dat{i_ex}.qc.p1amp{i_ch}(:) - dat{i_ex}.qc.p1amp_norm{i_ch}(:);
        tmp = tmp ./ dat{i_ex}.qc.p1amp_norm{i_ch}(:);
        if all(isnan(tmp))
            continue
        end
        
        if ENFORCE_ABOVE_100PA && (dat{i_ex}.qc.p1amp{i_ch}(1) < 100)
            continue
        end
        
        if strcmpi(dat{i_ex}.info.opsin, 'chronos') || strcmpi(dat{i_ex}.info.opsin, 'chronos_flx')
            plot(tmp, 'r')
            chronos_std = cat(1, chronos_std, nanstd(tmp));
        else
            plot(tmp, 'b')
            chief_std = cat(1, chief_std, nanstd(tmp));
        end
    end
end
title('Fractional Error of P1 Amp')
xlabel('Trial Number')
ylabel('Fractional Error')


allvals = cat(1, chief_std, chronos_std);
edges = linspace(min(allvals).*0.95, max(allvals).*1.05, 20);
f = figure;
f.Position = [440   135   508   663];
subplot(2,1,1), hold on,
hh = histogram(chief_std, edges);
plot(mean(chief_std), max(get(gca, 'ylim'))*.95, 'rv', 'markerfacecolor', 'w', 'linewidth', 2, 'markersize', 10)
hh.FaceColor = 'r';
xlabel('standard dev of fractional error, p1 amps')
ylabel('counts')
title('oChIEF')
subplot(2,1,2), hold on,
hh = histogram(chronos_std, edges);
hh.FaceColor = 'g';
plot(mean(chronos_std), max(get(gca, 'ylim'))*.95, 'gv', 'markerfacecolor', 'w', 'linewidth', 2, 'markersize', 10)
xlabel('standard dev of fractional error, p1 amps')
ylabel('counts')
title('Chronos')


% 
% Plot the raw values of the P1 amplitude to look for differences in
% stability b/w chronos and ChIEF
%
figure, hold on,
all_chronos = {};
all_chief = {};
for i_ex = 1:numel(dat)
    
    for i_ch = 1:2
        
        % check to make sure this neuron was defined
        has_Vclamp_data = ~strcmpi(dat{i_ex}.info.fid.vclamp(end-7:end), 'none.abf');
        isvalid = dat{i_ex}.info.HS_is_valid_Vclamp(i_ch);
        if ~isvalid || ~has_Vclamp_data 
            continue
        end
        
        tmp = dat{i_ex}.qc.p1amp{i_ch}(:);
        if all(isnan(tmp))
            continue
        end
        if ENFORCE_ABOVE_100PA && (tmp(1) < 100)
            continue
        end
        
        tmp = (tmp-tmp(1)) ./ tmp(1) .* 100;
        
            
            
        if strcmpi(dat{i_ex}.info.opsin, 'chronos') || strcmpi(dat{i_ex}.info.opsin, 'chronos_flx')
            %plot(tmp, 'g')
            all_chronos = cat(1, all_chronos, tmp);
        else
            %plot(tmp, 'r')
            all_chief = cat(1, all_chief, tmp);
        end
        
    end
end

% add the averages across expts
maxlength = max([cellfun(@numel, all_chronos) ; cellfun(@numel, all_chief)]);
all_chronos = cellfun(@(x) cat(1, x, nan(maxlength-numel(x), 1)), all_chronos, 'uniformoutput', false);
all_chronos = cat(2, all_chronos{:})';
all_chief = cellfun(@(x) cat(1, x, nan(maxlength-numel(x), 1)), all_chief, 'uniformoutput', false);
all_chief = cat(2, all_chief{:})';

xbar_chronos = nanmean(all_chronos, 1);
sem_chronos = nanstd(all_chronos, [], 1) ./ sqrt(sum(~isnan(all_chronos), 1));
shadedErrorBar(1:maxlength, xbar_chronos, sem_chronos, {'g', 'linewidth', 3});

xbar_chief = nanmean(all_chief, 1);
sem_chief = nanstd(all_chief, [], 1) ./ sqrt(sum(~isnan(all_chief), 1));
shadedErrorBar(1:maxlength, xbar_chief, sem_chief, {'r', 'linewidth', 3});

title('Percent change in P1 amplitude')
xlabel('Trial Number')
ylabel('Percent change from sweep 1')

%% CONTROLS: CO-VARIATION IN P1 INPUTS FOR PAIRED RECORDINGS
clc; close all;

% initalize the output structure
paired_pop = [];
paired_pop.amps = {};
paired_pop.rho_fixed = [];
paired_pop.rho_raw = [];
paired_pop.info = {};

% pull out the data
counter = 1;
Nexpts = numel(dat);
for i_ex = 1:Nexpts
    % is this a paired recording?
    is_paired = all(dat{i_ex}.info.HS_is_valid_Vclamp);
    if ~is_paired; continue; end
    
    % do the p1 amps exist
    if ~isfield(dat{i_ex}.qc, 'p1amp'); continue; end
    
    % cull cases where one of the cells is a SOM cell
    is_som = cellfun(@(x) ~isempty(x), regexpi(dat{i_ex}.info.cellType, 'somcre|ltsin|udfrs'));
    if any(is_som); continue; end
    
    % extract the p1 amps, and back out the non-stationarities using the
    % boxcar filtered estimate of the underlying P1 amps
    raw_p1_amps = cellfun(@(x) permute(x, [3,1,2]), dat{i_ex}.qc.p1amp, 'uniformoutput', false);
    boxcar_p1_amps = cellfun(@(x) permute(x, [2,1]), dat{i_ex}.qc.p1amp_norm, 'uniformoutput', false);
    residual_p1_amps = cellfun(@(x,y) x-y, raw_p1_amps, boxcar_p1_amps, 'uniformoutput', false);
    
    % cull low trial counts
    trl_count = cellfun(@(x) sum(~isnan(x)), residual_p1_amps);
    if any(trl_count < 10); continue; end
    
    % run the correlation, store the value
    [rho_fixed, pval] = corr(residual_p1_amps{1}, residual_p1_amps{2}, 'type', 'spearman', 'rows', 'complete');
    [rho_raw, ~] = corr(raw_p1_amps{1}, raw_p1_amps{2}, 'type', 'spearman', 'rows', 'complete');
    
    % compute the distance between the cells
    cell_distance = norm(dat{i_ex}.info.HS_xy_pos{1} - dat{i_ex}.info.HS_xy_pos{2});
    
    % store some things across experiments
    paired_pop.amps{counter} = residual_p1_amps;
    paired_pop.rho_raw(counter) = rho_raw;
    paired_pop.rho_fixed(counter) = rho_fixed;
    paired_pop.info{counter} = dat{i_ex}.info;
    paired_pop.cell_distance(counter) = cell_distance;
    
    % update the counter
    counter = counter + 1;
end


% histogram of rhow for all pairs
hf = figure;
hold on,
hp1 = cdfplot(paired_pop.rho_fixed);
hp1.LineWidth = 2;
hp1.Color = 'b';
hp1 = cdfplot(paired_pop.rho_raw);
hp1.LineWidth = 2;
hp1.Color = 'r';
legend({'fixed', 'raw'}, 'location', 'best')
xlabel('Spearman''s rho', 'fontsize', 14)
ylabel('Cumulative Density', 'fontsize', 14)


% scatter plot of all x,y pairs
tmp = cat(1, paired_pop.amps{:});
x_zscore = cellfun(@(x) (x-nanmean(x))./nanstd(x), tmp(:,1), 'uniformoutput', false);
y_zscore = cellfun(@(x) (x-nanmean(x))./nanstd(x), tmp(:,2), 'uniformoutput', false);

x_zscore = cat(1, x_zscore{:});
y_zscore = cat(1, y_zscore{:});
[rho, p] = corr(x_zscore, y_zscore, 'type', 'spearman', 'rows', 'complete');

hf = figure;
hp = plot(x_zscore, y_zscore, 'k.');
title(sprintf('Rho: %.3f, p = %.3f', rho, p))
xlabel('Z-score cell 1')
ylabel('Z-score cell 2')


% scatter plot of rho vs. distance
hf = figure;
hp = plot(paired_pop.cell_distance, paired_pop.rho_fixed, 'ko');
[rho, p] = corr(paired_pop.cell_distance(:), paired_pop.rho_fixed(:), 'type', 'spearman', 'rows', 'complete');
title(sprintf('Rho: %.3f, p = %.3f', rho, p))
xlabel('distance between cells (um)')
ylabel('Spearman rho')






%% PHOTOS OF THE RECORDING SITES

close all; clc


NORMTYPE = 'none';


group_data_img = {};
for i_ex = 1:numel(dat)
   
    % add a picture of the slice, only once per figure
    mousename = dat{i_ex}.info.mouseName;
    sitenum = dat{i_ex}.info.siteNum;
    cd([GL_DATPATH, mousename, filesep, 'Other'])
    d = dir;
    d.name;
    
    photoprefix = sprintf('site%s|cell%s', sitenum, sitenum);
    l_img = cellfun(@(x) ~isempty(regexpi(x, photoprefix)), structcat(d, 'name'));
    assert(sum(l_img)<=1, 'ERROR: too many images found')
    if any(l_img)
        img = imread(d(l_img).name);
        
        % store the data
        brainArea = dat{i_ex}.info.brainArea;
        if isfield(group_data_img, brainArea)
            tmp = group_data_img.(brainArea);
            tmp = cat(1, tmp, img);
            group_data_img.(brainArea) = tmp;
        else
            group_data_img.(brainArea) = {img};
        end
    end
    
end

hvas = fieldnames(group_data_img);
for i_hva = 1:numel(hvas)
    f = figure;
    f.Name = hvas{i_hva};
    f.Position = [139          26        1614         961];
    n_plts = numel(group_data_img.(hvas{i_hva}));
    n_cols = 5;
    n_rows = ceil(n_plts ./ n_cols);
    for i_plt = 1:n_plts
        subplot(n_rows, n_cols, i_plt)
        img = group_data_img.(hvas{i_hva}){i_plt};
        minpixval = min(img(:));
        maxpixval = round(0.9 .* max(img(:)));
        imshow(img, [minpixval, maxpixval])
        colormap gray
        set(gca, 'xtick', [], 'ytick', [])
    end
end

