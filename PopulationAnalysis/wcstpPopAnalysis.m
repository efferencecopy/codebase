%% NOTES

% Does Chronos cause artifical depression? Look for SOM cells that should
% evoke facilitation. Look for any instances of good facilitation using
% Chronos.

% Does oChIEF cause wonky STP? Strange facilitation followed by strong
% depression followed by rebound

% Do the oChIEF fits improve if I add an extra facilitation term?


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
end



%%%% DEFINE THE ANALYSIS PARAMS %%%%

params.pretime.vclamp = 0.002;     % seconds before pulse onset
params.posttime.vclamp = 0.015;    % seconds after pulse onset
params.pretime.dcsteps = 0.100;    % seconds before current onset
params.posttime.dcsteps = 0.300;   % seconds after current offset 
params.pretime.iclamp = 0.005;
params.posttime.iclamp = 0.015;    % actually a minimum value, real value stored in the dat structure, and depends on ISI
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% figure out which experiments should be analyzed
wb_path = [GL_DOCUPATH, 'Other_workbooks', filesep, 'wholeCellSTPCellList.xlsx'];
[~, ~, wb_expt] = xlsread(wb_path, 2);
if ~strcmpi(EXPTTYPE, 'all')
    header_idx = strcmpi(EXPTTYPE, wb_expt(1,:));
    assert(sum(header_idx) == 1)
    l_to_analyze = cellfun(@(x) numel(x)==1 && x==1, wb_expt(:, header_idx));
else
    l_to_analyze = ones(size(wb_expt,1), 1);
end
l_to_analyze(1) = []; % delete the header row
expt_idx = find(l_to_analyze);

% load in the workbook that contains all the experimental information
[~,~,wb_info] = xlsread(wb_path, 1);

% generate the header index info
for i_atrib = 1:size(wb_info,2)
    fldname = wb_info{1,i_atrib};
    fldname(isspace(fldname)) = [];
    hidx.(fldname) = i_atrib;
end
    
% now that the header is formed, delete the first row.
wb_info(1,:) = [];

% convert the file names into fully qualified paths so that par-for can run
% without calling a global
dcsteps_fpath = cellfun(@(x,y) strcat(GL_DATPATH, x, filesep, 'Physiology', filesep, y, '.abf'), wb_info(:,hidx.MouseName), wb_info(:, hidx.ABFDCsteps), 'uniformoutput', false);
vclamp_fpath = cellfun(@(x,y) strcat(GL_DATPATH, x, filesep, 'Physiology', filesep, y, '.abf'), wb_info(:,hidx.MouseName), wb_info(:, hidx.ABFOptostimVclamp), 'uniformoutput', false);
iclamp_fpath = cellfun(@(x,y) strcat(GL_DATPATH, x, filesep, 'Physiology', filesep, y, '.abf'), wb_info(:,hidx.MouseName), wb_info(:, hidx.ABFOptostimIclamp), 'uniformoutput', false);
wb_info(:,hidx.ABFDCsteps) = dcsteps_fpath;
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
    pool = parpool(21);
end


parfor i_ex = 1:Nexpts
    dat{i_ex} = wcstp_compile_data(attributes{i_ex}, hidx, params);
end

fprintf('All done importing data\n')

%% QULAITY CONTROL PLOTS

close all

for i_ex = 1:numel(dat)

    f = figure;
    f.Name = sprintf('Mouse %s, site %s', dat{i_ex}.info.mouseName, dat{i_ex}.info.siteNum);
    f.Position = [332  96 1259 665];
    
    for i_ch = 1:2
        
        isvalid = dat{i_ex}.info.HS_is_valid_Vclamp(i_ch);
        if ~isvalid
            continue
        end
        
        % plot the current step data set to help identify cell types
        if isfield(dat{i_ex}, 'dcsteps')
            Vm = dat{i_ex}.dcsteps.Vm_raw{i_ch};
            Icmd = dat{i_ex}.dcsteps.Icmd{i_ch};
            N = size(Vm,2);
            tt = [0:N-1] ./ dat{i_ex}.info.sampRate.dcsteps;
            if i_ch == 1; col = 1; else col = 3; end
            
            if any(Icmd<0)
                pltidx = sub2ind([4,3], col, 3);
                ha = subplot(3,4,pltidx);
                plot(tt, Vm(Icmd<0,:)');
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
            end
        end
        
        % figure out the subplot column number for the following axes
        if i_ch == 1; col = 2; else col = 4; end
        
        % series resistance
        if isfield(dat{i_ex}.qc, 'Rs') && ~all(isnan(dat{i_ex}.qc.Rs{i_ch}))
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
        if isfield(dat{i_ex}.qc, 'Verr') && ~all(isnan(dat{i_ex}.qc.verr{i_ch}))            
            pltidx = sub2ind([4,3],col,2);
            subplot(3,4,pltidx)
            tmp = squeeze(dat{i_ex}.qc.verr{i_ch});
            plot(tmp)
            ylim([min([0,min(tmp(:))]), max(tmp(:))*1.5])
            ylabel('SS Verr (mV)')
            xlabel('trial number')
            
        end
        
        % p1amps
        if isfield(dat{i_ex}.qc, 'p1amp') && ~all(isnan(dat{i_ex}.qc.p1amp{i_ch}))
            pltidx = sub2ind([4,3], col, 3);
            subplot(3,4,pltidx), hold on,
            tmp = squeeze(dat{i_ex}.qc.p1amp{i_ch});
            tmp_norm = squeeze(dat{i_ex}.qc.p1amp_norm{i_ch});
            plot(tmp)
            plot(tmp_norm, 'r')
            ylim([min([0,min(tmp(:))]), max(tmp(:))*1.5])
            ylabel('P1 Amp')
            xlabel('trial number')
            %set(gca, 'yscale', 'log', 'ylim', [0.3333, 3])
        end
            
        
    end
    drawnow
    
end


%% SUMMARY OF PASSIVE PROPERTIES

close all; clc

% define a set of attributes for each lineseries (or manifold) in the plot
% {CellType,  BrainArea,  OpsinType}
% Brain Area can be: 'AL', 'PM', 'AM', 'LM', 'any', 'med', 'lat'. CASE SENSITIVE
plotgroups = {
    'PY_L23',    'any', 'any';...
    'NPVIN',    'any', 'any';...
    };

groupdata.Rin_peak = repmat({[]}, 1, size(plotgroups, 1)); % should only have N cells, where N = size(plotgroups, 1). Each cell has a matrix with a cononicalGrid:
groupdata.Rin_asym = repmat({[]}, 1, size(plotgroups, 1));
% groupdata.Rin_vclamp = repmat({[]}, 1, size(plotgroups, 1)); % a place  holder for when this analysis is up and running
groupdata.Vrest = repmat({[]}, 1, size(plotgroups, 1));
groupdata.Depth = repmat({[]}, 1, size(plotgroups, 1));
groupdata.IVcurve_peak = repmat({{}}, 1, size(plotgroups,1));
groupdata.IVcurve_asym = repmat({{}}, 1, size(plotgroups,1));
groupdata.Ih_sag = repmat({[]}, 1, size(plotgroups,1));
groupdata.starttime = repmat({[]}, 1, size(plotgroups,1));


for i_ex = 1:numel(dat)
    
    for i_ch = 1:2
        
        % check to make sure this neuron was defined
        isvalid = dat{i_ex}.info.HS_is_valid_Vclamp(i_ch);
        if ~isvalid
            continue
        end
        
        % check the attributes
        ch_attribs = {dat{i_ex}.info.cellType{i_ch}, dat{i_ex}.info.brainArea, dat{i_ex}.info.opsin};
        l_nan = cellfun(@(x) all(isnan(x)), ch_attribs);
        ch_attribs(l_nan) = cellfun(@num2str, ch_attribs(l_nan), 'uniformoutput', false);
        
        % is this expt and channel cooresponds to one of the plot_groups in
        % terms of cellType and opsin
        l_cellType_match = cellfun(@(x) ~isempty(regexpi(ch_attribs{1}, x)), plotgroups(:,1)) | strcmpi(plotgroups(:,1), 'any');
        l_opsinMatch = cellfun(@(x) ~isempty(regexpi(ch_attribs{3}, x)), plotgroups(:,3)) | strcmpi(plotgroups(:,3), 'any');
        
        
        % determine if this experiment has a brain area that corresponds to
        % one of the plot_groups. Start by adding a 'medial', 'lateral'
        % assignment to the brain area
        expt_area = ch_attribs{2};
        if ~isempty(regexp(expt_area, 'AM', 'once')) || ~isempty(regexp(expt_area, 'PM', 'once'))
            expt_area = [expt_area, ' med'];
        elseif ~isempty(regexp(expt_area, 'AL', 'once')) || ~isempty(regexp(expt_area, 'LM', 'once'))
            expt_area = [expt_area, ' lat'];
        end
        l_brainArea_match = cellfun(@(x) ~isempty(regexp(expt_area, x, 'once')), plotgroups(:,2)) | strcmpi(plotgroups(:,2), 'any');
        
        group_idx = sum([l_cellType_match, l_brainArea_match, l_opsinMatch], 2) == 3;
        assert(sum(group_idx)<=1, 'ERROR: found too many group indicies')
        if sum(group_idx) == 0; continue; end
        
        % add the statistics to the appropriate cell array
        groupdata.Rin_peak{group_idx} = cat(1, groupdata.Rin_peak{group_idx}, dat{i_ex}.dcsteps.IVpeak.Rin{i_ch});
        groupdata.Rin_asym{group_idx} = cat(1, groupdata.Rin_asym{group_idx}, dat{i_ex}.dcsteps.IVasym.Rin{i_ch});
        groupdata.Vrest{group_idx} = cat(1, groupdata.Vrest{group_idx}, dat{i_ex}.dcsteps.Vrest{i_ch});
        groupdata.Depth{group_idx} = cat(1, groupdata.Depth{group_idx}, dat{i_ex}.info.cellDepth_um(i_ch));
        groupdata.IVcurve_peak{group_idx} = cat(1, groupdata.IVcurve_peak{group_idx}, {dat{i_ex}.dcsteps.IVpeak.raw{i_ch}});
        groupdata.IVcurve_asym{group_idx} = cat(1, groupdata.IVcurve_asym{group_idx}, {dat{i_ex}.dcsteps.IVasym.raw{i_ch}});
        groupdata.Ih_sag{group_idx} = cat(1, groupdata.Ih_sag{group_idx}, dat{i_ex}.dcsteps.Ih_sag{i_ch}(2));
        groupdata.starttime{group_idx} = cat(1, groupdata.starttime{group_idx}, dat{i_ex}.info.fileStartTime_24hrs);
    end
end


%
% plot histograms of Input resistance measured three different ways
%
%%%%%%%%%%%%%%%%%%%%%%%%%
f=figure;
Ngroups = (numel(groupdata.Rin_peak));
groupcolors = lines(Ngroups);
f.Position = [276         333        1276         478];
allNums = cat(1, groupdata.Rin_peak{:},  groupdata.Rin_asym{:});
edges = linspace(min(allNums)-10, max(allNums)+10, 30);
for i_group = 1:numel(groupdata.Rin_asym)
    % current clamp, peak vals
    pltidx = sub2ind([3, Ngroups], 1, i_group);
    subplot(Ngroups, 3, pltidx)
    
    h = histogram(groupdata.Rin_peak{i_group}, edges);
    h.FaceColor = groupcolors(i_group,:);
    xbar = nanmean(groupdata.Rin_peak{i_group});
    hold on,
    plot(xbar, 0.5, 'kv', 'markerfacecolor', 'k', 'markersize', 5)
    legtext = sprintf('%s, %s, %s', plotgroups{i_group,1},plotgroups{i_group, 2},plotgroups{i_group,3} );
    legend(legtext, 'Location', 'northeast')
    legend boxoff
    xlabel('Input R (MOhm)')
    ylabel('counts')
    if i_group == 1
        title('Rin (iClamp, peak)')
    end
    
    % current clamp asym vals
    pltidx = sub2ind([3, Ngroups], 2, i_group);
    subplot(Ngroups, 3, pltidx)
    
    h = histogram(groupdata.Rin_asym{i_group}, edges);
    h.FaceColor = groupcolors(i_group,:);
    xbar = nanmean(groupdata.Rin_asym{i_group});
    hold on,
    plot(xbar, 0.5, 'kv', 'markerfacecolor', 'k', 'markersize', 5)
    xlabel('Input R (MOhm)')
    ylabel('counts')
    if i_group == 1
        title('Rin (iClamp, asym)')
    end
    
    % voltage clamp
end


%
% plot a histograms of Vrest and depth for each group
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=figure;
Ngroups = (numel(groupdata.Rin_peak));
groupcolors = lines(Ngroups);
f.Position = [276         333        1276         478];
for i_group = 1:numel(groupdata.Vrest)
    
    % Vrest
    allNums = cat(1, groupdata.Vrest{:});
    edges = linspace(min(allNums)-10, max(allNums)+10, 30);
    
    pltidx = sub2ind([3, Ngroups], 1, i_group);
    subplot(Ngroups, 3, pltidx)
    
    h = histogram(groupdata.Vrest{i_group}, edges);
    h.FaceColor = groupcolors(i_group,:);
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
    
    pltidx = sub2ind([3, Ngroups], 2, i_group);
    subplot(Ngroups, 3, pltidx)
    
    h = histogram(groupdata.Depth{i_group}, edges);
    h.FaceColor = groupcolors(i_group,:);
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
% NEED TO ADD histogram MEMBRANE TAU TO FIGURE WITH Vrest AND Rm
% 




%
% scatter plots of cell depth vs. Rin, Vrest, tau
%
%%%%%%%%%%%%%%%%%%%%%%%%%
f=figure;
Ngroups = (numel(groupdata.Rin_asym));
groupcolors = lines(Ngroups);
f.Position = [49         305        1333         602];
for i_group = 1:numel(groupdata.Rin_asym)
    
    % cell depth vs. Rin
    allRin = cat(1, groupdata.Rin_asym{:});
    allDepth = cat(1, groupdata.Depth{:});
    pltidx = sub2ind([3, Ngroups], 1, i_group);
    subplot(Ngroups, 3, pltidx)
    
    X = groupdata.Depth{i_group};
    Y = groupdata.Rin_asym{i_group};
    plot(X, Y, 'o', 'color', groupcolors(i_group,:), 'markerfacecolor', groupcolors(i_group,:), 'markersize', 5)
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
    allRin = cat(1, groupdata.Rin_asym{:});
    allTimes = cat(1, groupdata.starttime{:});
    pltidx = sub2ind([3, Ngroups], 2, i_group);
    subplot(Ngroups, 3, pltidx)
    
    X = groupdata.starttime{i_group};
    Y = groupdata.Rin_asym{i_group};
    plot(X, Y, 'o', 'color', groupcolors(i_group,:), 'markerfacecolor', groupcolors(i_group,:), 'markersize', 5)
    xlabel('start time (24hrs)')
    ylabel('Input Resistance')
    if i_group == 1
        title('Rin vs. Incubation Time')
    end
    xlim([min(allTimes).*.95, max(allTimes).*1.05])
    ylim([min(allRin).*.95, max(allRin).*1.05])
    
    % scatter plot of Vrest vs depth
    allVrest = cat(1, groupdata.Vrest{:});
    pltidx = sub2ind([3, Ngroups], 3, i_group);
    subplot(Ngroups, 3, pltidx)
    
    X = groupdata.Depth{i_group};
    Y = groupdata.Vrest{i_group};
    plot(X, Y, 'o', 'color', groupcolors(i_group,:), 'markerfacecolor', groupcolors(i_group,:), 'markersize', 5)
    xlabel('Cell Depth (um)')
    ylabel('Resting Potential (mV)')
    if i_group == 1
        title('Vrest vs. Cell Depth')
    end
    xlim([min(allDepth).*.95, max(allDepth).*1.05])
    ylim([min(allVrest).*.95, max(allVrest).*1.05])
    
end



%
% Line plots of I-V curves, and peak vs. steady state
%
%%%%%%%%%%%%%%%%%%%%%%%%%
f=figure;
Ngroups = (numel(groupdata.IVcurve_peak));
groupcolors = lines(Ngroups);
f.Position = [114 426 1667 456];
legtext = {};

% I-V curves for peak values
subplot(1,3,1), hold on,
%allX = % for linear interpolation
for i_group = 1:numel(groupdata.IVcurve_peak)
    for i_cell = 1:numel(groupdata.IVcurve_peak{i_group})
        X = groupdata.IVcurve_peak{i_group}{i_cell}(:,1);
        Y = groupdata.IVcurve_peak{i_group}{i_cell}(:,2);
        plot(X, Y, '-', 'color', groupcolors(i_group,:))
    end
    %legtext{i_group} = sprintf('%s, %s, %s', plotgroups{i_group,1},plotgroups{i_group, 2},plotgroups{i_group,3} );
end
title('IV curve, peak vals')
xlabel('current injection (pA)')
ylabel('voltage response (mV)')
% legend(legtext, 'Location', 'best')
% legend boxoff

% I-V curves for asym values
subplot(1,3,2), hold on,
for i_group = 1:numel(groupdata.IVcurve_asym)
    for i_cell = 1:numel(groupdata.IVcurve_asym{i_group})
        X = groupdata.IVcurve_asym{i_group}{i_cell}(:,1);
        Y = groupdata.IVcurve_asym{i_group}{i_cell}(:,2);
        plot(X, Y, '-', 'color', groupcolors(i_group,:))
    end
    %legtext{i_group} = sprintf('%s, %s, %s', plotgroups{i_group,1},plotgroups{i_group, 2},plotgroups{i_group,3} );
end
title('IV curve, asym vals')
xlabel('current injection (pA)')
ylabel('voltage response (mV)')
% legend(legtext, 'Location', 'best')
% legend boxoff

% directly compare peak vs. steady state DC response 
subplot(1,3,3), hold on,
for i_group = 1:numel(groupdata.IVcurve_peak)
    for i_cell = 1:numel(groupdata.IVcurve_peak{i_group})
        X = groupdata.IVcurve_peak{i_group}{i_cell}(:,1);
        Y_peak = groupdata.IVcurve_peak{i_group}{i_cell}(:,2);
        Y_asym = groupdata.IVcurve_asym{i_group}{i_cell}(:,2);
        plot(X, Y_peak-Y_asym, '-', 'color', groupcolors(i_group,:))
    end
end
axis tight
title('Peak - Steady State')
xlabel('DC injection (pA)')
ylabel('Difference (mV)')


% compare distributions of Ih sag
f = figure;
f.Position = [524   406   493   550];
allNums = cat(1, groupdata.Ih_sag{:});
edges = linspace(min(allNums)*0.90, max(allNums)*0.11, 30);
for i_group = 1:numel(groupdata.Vrest)
    
    pltidx = sub2ind([1, Ngroups], 1, i_group);
    subplot(Ngroups, 1, pltidx)
    
    h = histogram(groupdata.Ih_sag{i_group}, edges);
    h.FaceColor = groupcolors(i_group,:);
    xbar = nanmean(groupdata.Ih_sag{i_group});
    hold on,
    plot(xbar, 0.5, 'kv', 'markerfacecolor', 'k', 'markersize', 5)
    legtext = sprintf('%s, %s, %s', plotgroups{i_group,1},plotgroups{i_group, 2},plotgroups{i_group,3} );
    legend(legtext, 'Location', 'northeast')
    legend boxoff
    xlabel('Ih sag (mV)')
    ylabel('counts')
    if i_group == 1
        title('Ih Sag')
    end

end

%% PLOT THE RAW VCLAMP WAVEFORMS FOLLOWING EACH PULSE

close all


PLOT_ALL_TRIALS = true;
DEBUG_MEAN = false;
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
        isvalid = dat{i_ex}.info.HS_is_valid_Iclamp(i_ch);
        if ~isvalid
            continue
        end
        
        
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


PLOT_ALL_TRIALS = true;
DEBUG_MEAN = false;
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
            baselines = dat{i_ex}.iclamp.(conds{idx}).raw.bkgndVm{i_ch};
            snips = bsxfun(@plus, snips, baselines);
            
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
            peakVals = -1 .* dat{i_ex}.iclamp.(conds{idx}).stats.EPSPamp{i_ch};
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




%% ESTIMATE TIME CONSTANTS OF SHORT TERM PLASTICITY

TRAINSET = 'recovery';  % could be 'rit', 'recovery', 'all'
PLOTTRAININGDATA = true;
NORMALIZEDATA = true;
FITAVERAGEDATA = true;
FITRECOVERYPULSE = false;


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


for i_ex = 1:numel(dat)
   clc
   hf = figure;
   hf.Units = 'Normalized';
   hf.Position = [0.3347    0.0422    0.6639    0.7489];
   hf.Name = sprintf('Mouse %s, site %s, opsin: %s.  Train with: %s, Plot training set: %d',...
                      dat{i_ex}.info.mouseName, dat{i_ex}.info.siteNum, dat{i_ex}.info.opsin, TRAINSET, PLOTTRAININGDATA);
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
       [pOnTimes, rawAmps, p1Amps, raw_xbar, raw_sem] = deal({});
       condnames = fieldnames(dat{i_ex}.expt);
       l_trainingSet = isTrainingSet(condnames);
       for i_cond = 1:numel(condnames)
           
           % check to make sure there are data for this condition
           if isempty(dat{i_ex}.expt.(condnames{i_cond}).stats.EPSCamp{i_ch})
               pOnTimes{i_cond} = [];
               rawAmps{i_cond} = [];
               raw_sem{i_cond} = [];
               raw_xbar{i_cond} = [];
               p1Amps{i_cond} = [];
               continue
           end
           
           % grab pOnTimes
           pOnTimes{i_cond} = dat{i_ex}.expt.(condnames{i_cond}).pOnTimes;
           rawAmps{i_cond} = dat{i_ex}.expt.(condnames{i_cond}).stats.EPSCamp{i_ch};
           
           % normalize to the p1amp_norm if need be
           trlNums = dat{i_ex}.expt.(condnames{i_cond}).realTrialNum{i_ch};
           normfacts = dat{i_ex}.qc.p1amp_norm{i_ch}(trlNums);
           normfacts = permute(normfacts, [1,3,2]);
           
           % delete the recovery pulse if need be
           if ~FITRECOVERYPULSE
               has_recov_pulse = dat{i_ex}.expt.(condnames{i_cond}).tdict(4) > 0;
               if has_recov_pulse
                   pOnTimes{i_cond}(end) = [];
                   rawAmps{i_cond} = rawAmps{i_cond}(1:end-1,1,:);
               end
           end
           
           
           if NORMALIZEDATA
               
               rawAmps{i_cond} = bsxfun(@rdivide, rawAmps{i_cond}, normfacts);
               raw_sem{i_cond} = stderr(rawAmps{i_cond}, 3);
               raw_xbar{i_cond} = mean(rawAmps{i_cond},3);
               if FITAVERAGEDATA
                   p1Amps{i_cond} = 1;
               else
                   p1Amps{i_cond} = ones(1, size(rawAmps{i_cond},3));
               end
               
           else
               
               % leave 'rawAmps' unchanged, but update these:
               raw_sem{i_cond} = stderr(rawAmps{i_cond}, 3);
               raw_xbar{i_cond} = mean(rawAmps{i_cond},3);
               if FITAVERAGEDATA
                   p1Amps{i_cond} = raw_xbar{i_cond}(1);
               else
                   p1Amps{i_cond} = squeeze(normfacts);
               end
               
           end
           
       end
       
       % allocate the training data
       [training_amps, training_pOnTimes, training_p1Amps] = deal({});
       if FITAVERAGEDATA
           training_amps = raw_xbar(l_trainingSet);
       else
           training_amps = rawAmps(l_trainingSet);
       end
       training_pOnTimes = pOnTimes(l_trainingSet);
       training_p1Amps = p1Amps(l_trainingSet);
       
       
       % if there were no training data, then move along,
       if isempty(training_amps)
           chempty(i_ch) = true;
           if i_ch == 1
               continue
           elseif all(chempty)
               close(hf)
               continue
           end
       end
       NtrainingSets = sum(l_trainingSet);
       
       % look for instances where there are pOnTimes but no data (could
       % happen if some sweeps get deleted from one HS but not the other.
       l_empty = cellfun(@isempty, training_amps);
       training_amps(l_empty) = [];
       training_pOnTimes(l_empty) = [];
       training_p1Amps(l_empty) = [];
       
       % fit the RIT data, but only if the fit params do not already exist
       [d, f, dTau, fTau] = deal([]); %#ok<*ASGLU>
       MAKENEWFITS = true;
       if [isfield(dat{i_ex}, 'stpfits') ...
               && isfield(dat{i_ex}.stpfits, 'modelParams') ...
               && numel(dat{i_ex}.stpfits.modelParams)>= i_ch ...
               && ~isempty(dat{i_ex}.stpfits.modelParams{i_ch})] % end of if conditional
           
           % if you've made it here, then the old params are potentially
           % still appliciable, but I should check to make sure that they
           % were fit using the same technique
           matches = [strcmpi(dat{i_ex}.stpfits.trainingSet, TRAINSET); ...
                      dat{i_ex}.stpfits.normalizeData == NORMALIZEDATA; ...
                      dat{i_ex}.stpfits.fitWithAvg == FITAVERAGEDATA; ...
                      dat{i_ex}.stpfits.fitRecovPulse == FITRECOVERYPULSE];
           
          if all(matches)
              d = dat{i_ex}.stpfits.modelParams{i_ch}(1:2);
              f = dat{i_ex}.stpfits.modelParams{i_ch}(3);
              dTau = dat{i_ex}.stpfits.modelParams{i_ch}(4:5);
              fTau = dat{i_ex}.stpfits.modelParams{i_ch}(6);

              MAKENEWFITS = ralse;
          end
       end
       if MAKENEWFITS
           keyboard
           [d, f, dTau, fTau] = fitTau2STP(training_amps, training_pOnTimes, training_p1Amps, 'multistart');
       end
       
       % predict all the data
       pred = {};
       for i_cond = 1:numel(condnames)
           if ~isempty(rawAmps{i_cond})
               A0 = mean(p1Amps{i_cond});
               pred{i_cond} = predictPSCfromTau(pOnTimes{i_cond}, d, dTau, f, fTau, A0);
           else
               pred{i_cond} = [];
           end
       end
       
       
       % plot the training or cross validation data set, and the prediction
       figure(hf)
       if PLOTTRAININGDATA
           l_condsToPlot = isTrainingSet(condnames);
       else
           l_condsToPlot = ~isTrainingSet(condnames);
       end
       idx_condsToPlot = find(l_condsToPlot);
       xlims = [inf -inf];
       ylims = [inf -inf];
       hs = [];
       for i_cond = 1:numel(idx_condsToPlot)
           
           typeIdx = idx_condsToPlot(i_cond);
           if isempty(raw_xbar{typeIdx}); continue; end
           
           pltIdx = sub2ind([4, numel(idx_condsToPlot)], i_ch+1, i_cond);
           hs(i_cond) = subplot(numel(idx_condsToPlot), 4, pltIdx); hold on,
           
           xx = pOnTimes{typeIdx};
           my_errorbar(xx, raw_xbar{typeIdx}, raw_sem{typeIdx}, 'k');
           plot(xx, pred{typeIdx}, 'r', 'linewidth', 2)
           xlims(1) = min([min(xx), xlims(1)]);
           xlims(2) = max([max(xx), xlims(2)]);
           yvals = get(gca, 'ylim');
           ylims(1) = min([yvals(1), ylims(1)]);
           ylims(2) = max([yvals(2), ylims(2)]);
       end
       if ~isempty(hs) && sum(hs)>0
           set(hs(hs~=0), 'XLim', xlims, 'YLim', ylims)
       end
       
       
       
       % make a scatter plot of all predicted and actual PSC amps
       hs = [];
       training_raw = [];
       training_pred = [];
       crossval_raw = [];
       crossval_pred = [];
       for i_cond = 1:numel(condnames)
           if isTrainingSet(condnames{i_cond}); pltclr = 'k';else pltclr = 'r';end
           if i_ch == 1; pltcol=1; else pltcol=4; end
           
           if isempty(raw_xbar{i_cond}); continue; end
           
           tmp_raw = raw_xbar{i_cond};
           tmp_pred = pred{i_cond};
           
           tmp_raw = tmp_raw(:);
           tmp_pred = tmp_pred(:);
           assert(all(size(tmp_pred) == size(tmp_raw)))
           
           pltIdx = sub2ind([4, 3], pltcol, 1);
           hs = subplot(3, 4, pltIdx); hold on,
           plot(tmp_raw, tmp_pred, '.', 'color', pltclr)
           axis tight
           
           % concatenate the training data
           if isTrainingSet(condnames{i_cond})
               training_raw = cat(1, training_raw, tmp_raw(:));
               training_pred = cat(1, training_pred, tmp_pred(:));
           end
           
           % concatenate the cross-validation data
           if ~isTrainingSet(condnames{i_cond})
               crossval_raw = cat(1, crossval_raw, tmp_raw(:));
               crossval_pred = cat(1, crossval_pred, tmp_pred(:));
           end
       end
       maxval = max([hs.XLim, hs.YLim]);
       minval = min([hs.XLim, hs.YLim]);
       plot([minval, maxval], [minval, maxval], 'k--')
       %hs.XScale = 'log';
       %hs.YScale = 'log';
       xlabel('raw EPSC amp')
       ylabel('pred amp')
       title(sprintf('num RITs fit: %d', NtrainingSets))
       
       pltIdx = sub2ind([4, 3], pltcol, 2);
       subplot(3,4,pltIdx), hold on,
       resid = training_pred - training_raw;
       histogram(resid)
       plot(mean(resid), 10, 'rv', 'markerfacecolor', 'r')
       R2_train = 1 - (sum(resid.^2) ./ sum((training_raw - mean(training_raw)).^2));
       xlabel('pred-real')
       title(sprintf('R2 = %.2f', R2_train));
       
       
       
       % plot cross validation stuff
       R2_crossvalid = [];
       if ~isempty(crossval_raw)
           pltIdx = sub2ind([4, 3], pltcol, 3);
           subplot(3,4,pltIdx), hold on,
           resid = crossval_pred - crossval_raw;
           histogram(resid);
           plot(mean(resid), 5, 'rv', 'markerfacecolor', 'r')
           R2_crossvalid = 1 - (sum(resid.^2) ./ sum((crossval_raw-mean(crossval_raw)).^2));
           xlabel('cross-valid (pred-real)')
           title(sprintf('R2 = %.2f', R2_crossvalid));
       end
       
       
       % store some parameters in the dat array
       dat{i_ex}.stpfits.trainingSet = TRAINSET;
       dat{i_ex}.stpfits.normalizeData = NORMALIZEDATA;
       dat{i_ex}.stpfits.fitWithAvg = FITAVERAGEDATA;
       dat{i_ex}.stpfits.fitRecovPulse = FITRECOVERYPULSE;
       dat{i_ex}.stpfits.modelParams{i_ch} = [d, f, dTau, fTau];
       dat{i_ex}.stpfits.R2.training{i_ch} = R2_train;
       dat{i_ex}.stpfits.R2.crossvalid{i_ch} = R2_crossvalid;

   end
   drawnow
   
   
end

%% iCLAMP POPULATION ANALYSIS (DATA COLLECTION)

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
                iclamp_pop.dat{i_ex}.xbar{i_tf}{i_ch} = [];
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


%% iCLAMP POPULATION ANALYSIS (PLOTTING)

clc; close all

PLOTERRBAR = true;

% define a set of attributes for each lineseries (or manifold) in the plot
% {CellType,  BrainArea,  OpsinType}
% Brain Area can be: 'AL', 'PM', 'AM', 'LM', 'any', 'med', 'lat'. CASE SENSITIVE
plotgroups = {
    'PY',    'PM', 'any';...
    'PY',    'AL', 'any';...
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
        l_nan = cellfun(@(x) all(isnan(x)), ch_attribs);
        ch_attribs(l_nan) = cellfun(@num2str, ch_attribs(l_nan), 'uniformoutput', false);
        
        % is this expt and channel cooresponds to one of the plot_groups in
        % terms of cellType and opsin
        l_cellType_match = cellfun(@(x) ~isempty(regexpi(ch_attribs{1}, x)), plotgroups(:,1)) | strcmpi(plotgroups(:,1), 'any');
        l_opsinMatch = cellfun(@(x) ~isempty(regexpi(ch_attribs{3}, x)), plotgroups(:,3)) | strcmpi(plotgroups(:,3), 'any');
        
        % determine if this experiment has a brain area that corresponds to
        % one of the plot_groups. Start by adding a 'medial', 'lateral'
        % assignment to the brain area
        expt_area = ch_attribs{2};
        if ~isempty(regexp(expt_area, 'AM', 'once')) || ~isempty(regexp(expt_area, 'PM', 'once'))
            expt_area = [expt_area, ' med'];
        elseif ~isempty(regexp(expt_area, 'AL', 'once')) || ~isempty(regexp(expt_area, 'LM', 'once'))
            expt_area = [expt_area, ' lat'];
        end
        l_brainArea_match = cellfun(@(x) ~isempty(regexp(expt_area, x, 'once')), plotgroups(:,2)) | strcmpi(plotgroups(:,2), 'any');
        
        group_idx = sum([l_cellType_match, l_brainArea_match, l_opsinMatch], 2) == 3;
        assert(sum(group_idx)<=1, 'ERROR: found too many group indicies')
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
            groupdata_raw{group_idx}{tf_idx} = cat(1, groupdata_raw{group_idx}{tf_idx}, tmpdat);
        end

    end
end


% plot the data, all line series
f = figure;
f.Units = 'normalized';
f.Position = [0.1, 0.01, 0.3, 0.9];
Ntfs = numel(allTFs);
groupcolors = {'r', 'b', 'g'};
for i_tf = 1:Ntfs
    subplot(Ntfs, 1, i_tf), hold on,
    for i_grp = 1:numel(groupdata_raw)
        tmp = groupdata_raw{i_grp}{i_tf};
        xbar = nanmean(tmp, 1);
        sem = nanstd(tmp, [], 1) ./ sqrt(sum(~isnan(tmp),1));
        if PLOTERRBAR
            shadedErrorBar(1:size(tmp,2), xbar, sem, {'color', groupcolors{i_grp}, 'linewidth', 2});
        else
            plot(tmp', '-', 'color', groupcolors{i_grp})
            %plot(mean(tmp,1), '-', 'color', groupcolors{i_grp}, 'linewidth', 3)
        end
        
    end
end



%% PAIRED PULSE PLASTICITY MANIFOLDS (DATA COLLECTION)

% loop through the experiments. Pull out the trains data. Ignore the
% recovery train (if present) and aggregate across recovery conditions.
pprpop = [];
pprpop.TFsAllExpts = [];
pprpop.MaxNPulses = 0;
for i_ex = 1:numel(dat)
    
    % find the normal trains. Assume the field name is NOT 'ritv'
    condnames = fieldnames(dat{i_ex}.expt);
    l_trains = ~strncmp(condnames, 'RITv', 4);
    if sum(l_trains)==0; continue; end % no trains data
    trainParams = cellfun(@(x) dat{i_ex}.expt.(condnames{x}).tdict, mat2cell(find(l_trains), ones(sum(l_trains),1), 1), 'uniformoutput', false);
    trainParams = cat(1, trainParams{:});
    
    % make sure the pulse amplitude and width were identical across ttypes
    assert(numel(unique(trainParams(:,1)))==1, 'ERROR: more than one pulse amplitude')
    assert(numel(unique(trainParams(:,2)))==1, 'ERROR: more than one pulse width');
    
    % identify the unique TF conditions for this experiment, and update the
    % running log of TFs used across all experiments
    uniqueTFs = unique(trainParams(:,3));
    tmp = cat(1, pprpop.TFsAllExpts, uniqueTFs);
    pprpop.TFsAllExpts = unique(tmp);
    
    % store the stim params for all stim types, which will be useful for
    % indexing later.
    allStimParams = cellfun(@(x) dat{i_ex}.expt.(condnames{x}).tdict, mat2cell((1:numel(l_trains))', ones(numel(l_trains),1), 1), 'uniformoutput', false);
    allStimParams = cat(1, allStimParams{:});
    
    % store some metadata
    pprpop.info{i_ex} = dat{i_ex}.info;
    
    
    % aggregate data within TF conditons (Separately for each
    % recording channel
    for i_tf = 1:numel(uniqueTFs);
        
        tfidx = find(trainParams(:,3) == uniqueTFs(i_tf)); % condnames that contain the train with a particular TF
        
        for i_ch = 1:2
            
            % check to make sure this neuron was defined
            isvalid = dat{i_ex}.info.HS_is_valid_Vclamp(i_ch);
            if ~isvalid
                pprpop.dat{i_ex}.xbar{i_tf}{i_ch} = [];
                pprpop.tfs{i_ex}{i_ch} = [];
                continue
            end
            
            % iterate over the trains with the same freq. 
            catdat = [];
            for i_cond = 1:numel(tfidx)
                
                % pull out the data
                condIdx = ismember(allStimParams, trainParams(tfidx(i_cond),:), 'rows');
                assert(sum(condIdx)==1, 'ERROR, found zero or more than 1 instance of the trial type')
                tmpdat = dat{i_ex}.expt.(condnames{condIdx}).stats.EPSCamp{i_ch}; % [Npulses, 1, Nsweeps]
                if ~isempty(tmpdat)
                    % normalize by the p1Amp_norm
                    realTrlNums = dat{i_ex}.expt.(condnames{condIdx}).realTrialNum{i_ch};
                    p1Amp_norm = dat{i_ex}.qc.p1amp_norm{i_ch}(realTrlNums);
                    assert(~any(isnan(p1Amp_norm)), 'ERROR: scale factor is a nan');
                    p1Amp_norm = permute(p1Amp_norm, [1,3,2]);
                    tmpdat = bsxfun(@rdivide, tmpdat, p1Amp_norm); % [Npulses, 1, Nsweeps]
                    
                    % store in a matrix for averaging later
                    catdat = cat(3, catdat, tmpdat);
                end
            end
            
            % check to make sure there are data for these conditions. Even
            % though this recording channel should be defined (See above),
            % it's possible there are no data due to deletion of single
            % sweeps
            if isempty(catdat)
                pprpop.dat{i_ex}.xbar{i_tf}{i_ch} = [];
                pprpop.tfs{i_ex}{i_ch} = [];
            else
                
                % now average across trails, and re-normalize to the first
                % pulse. Store in the population data structure.
                avg = mean(catdat,3);
                ppr = avg./avg(1);
                
                % delete the last pulse (if it's a recovery pulse, but
                % store the recovery pulse in a different field of the
                % population structure
                isrecovery = trainParams(tfidx(i_cond),4) > 0;
                if isrecovery
                    ppr(end) = [];
                end
                
                pprpop.dat{i_ex}.xbar{i_tf}{i_ch} = ppr;
                pprpop.tfs{i_ex}{i_ch} = uniqueTFs;
                pprpop.MaxNPulses = max([pprpop.MaxNPulses, numel(ppr)]);
                
                
                
                % make a smooth manifold for this neuron.
                % Assume TF = 10 : 50;
                if isfield(dat{i_ex}, 'stpfits')
                    params = dat{i_ex}.stpfits.modelParams{i_ch};
                    isi_ms = fliplr([1000/50 : 1 : 1000/10]);
                    NumPulses = 10;
                    
                    
                    smoothManifold = nan(NumPulses, numel(isi_ms));
                    for i_isi = 1:numel(isi_ms)
                        A0 = 1;
                        tmp_pOntimes_ms = 0 : isi_ms(i_isi) : (isi_ms(i_isi)*NumPulses)-1;
                        tmp_pOntimes_sec = tmp_pOntimes_ms ./ 1000;
                        smoothManifold(:,i_isi) = predictPSCfromTau(tmp_pOntimes_sec, params(1:2), params(4:5), params(3), params(6), A0);
                    end
                    pprpop.smoothManifold{i_ex}{i_ch} = smoothManifold;
                    pprpop.smoothManifold_isi{i_ex}{i_ch} = isi_ms;
                end
                
            end
        end
        
    end
    
end


%% PAIRED PULSE PLASTICITY MANIFOLDS (PLOTS)
clc; close all

PLOT_SMOOTH_MANIFOLD = false;
PLOT_RAW_DATA = true;
PLOT_INDIVIDUAL_DATASETS = false;

% define a set of attributes for each lineseries (or manifold) in the plot
% {CellType,  BrainArea,  OpsinType}
% Brain Area can be: 'AL', 'PM', 'AM', 'LM', 'any', 'med', 'lat'. CASE SENSITIVE
plotgroups = {
    'PY',    'any', 'any';...
    'FS',    'any', 'any';...
    'NPVIN',    'any', 'any';...
    };

groupdata_raw = repmat({[]}, 1, size(plotgroups, 1)); % should only have N cells, where N = size(plotgroups, 1). Each cell has a matrix with a cononicalGrid:
groupdata_smooth =  repmat({[]}, 1, size(plotgroups, 1));
groupexpinds =  repmat({[]}, 1, size(plotgroups, 1));
groupchinds = repmat({[]}, 1, size(plotgroups, 1));
canonicalGrid = nan(pprpop.MaxNPulses, numel(pprpop.TFsAllExpts)); % need to get the maxNpulses into the pprpop struct
allTFs = pprpop.TFsAllExpts;

% iterate over the experiments. For each recording channel, determine what
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
        l_nan = cellfun(@(x) all(isnan(x)), ch_attribs);
        ch_attribs(l_nan) = cellfun(@num2str, ch_attribs(l_nan), 'uniformoutput', false);
        
        % is this expt and channel cooresponds to one of the plot_groups in
        % terms of cellType and opsin
        l_cellType_match = cellfun(@(x) ~isempty(regexpi(ch_attribs{1}, x)), plotgroups(:,1)) | strcmpi(plotgroups(:,1), 'any');
        l_opsinMatch = cellfun(@(x) ~isempty(regexpi(ch_attribs{3}, x)), plotgroups(:,3)) | strcmpi(plotgroups(:,3), 'any');
        
        % determine if this experiment has a brain area that corresponds to
        % one of the plot_groups. Start by adding a 'medial', 'lateral'
        % assignment to the brain area
        expt_area = ch_attribs{2};
        if ~isempty(regexp(expt_area, 'AM', 'once')) || ~isempty(regexp(expt_area, 'PM', 'once'))
            expt_area = [expt_area, ' med'];
        elseif ~isempty(regexp(expt_area, 'AL', 'once')) || ~isempty(regexp(expt_area, 'LM', 'once'))
            expt_area = [expt_area, ' lat'];
        end
        l_brainArea_match = cellfun(@(x) ~isempty(regexp(expt_area, x, 'once')), plotgroups(:,2)) | strcmpi(plotgroups(:,2), 'any');
        
        group_idx = sum([l_cellType_match, l_brainArea_match, l_opsinMatch], 2) == 3;
        assert(sum(group_idx)<=1, 'ERROR: found too many group indicies')
        if sum(group_idx) == 0; continue; end
        
        % add data to the appropriate group data array
        tmpgrid = canonicalGrid;
        ch_tfs = pprpop.tfs{i_ex}{i_ch};
        for i_tf = 1:numel(ch_tfs)
            grid_col_idx = allTFs == ch_tfs(i_tf);
            tmpdat = pprpop.dat{i_ex}.xbar{i_tf}{i_ch};
            Npulses = numel(tmpdat);
            tmpgrid(1:Npulses, grid_col_idx) = tmpdat(:);
        end
        
        if all(isnan(tmpgrid(:))); continue; end
        
        % aggregate the data
        groupdata_raw{group_idx} = cat(3, groupdata_raw{group_idx}, tmpgrid);
        groupexpinds{group_idx} = cat(1, groupexpinds{group_idx}, i_ex);
        groupchinds{group_idx} = cat(1, groupchinds{group_idx}, i_ch);
        
        if PLOT_SMOOTH_MANIFOLD
            groupdata_smooth{group_idx} = cat(3, groupdata_smooth{group_idx}, pprpop.smoothManifold{i_ex}{i_ch});
        end
    end
end


% Plot average manfold. any global (above) that has multiple elements will
% be ploted against eachother, but the other globals will be held constant.
hf = figure;
hold on,
plotcolors = {'r', 'b', 'g'};
for i_group = 1:numel(groupdata_raw)
    
    if PLOT_RAW_DATA
        grid_average = nanmean(groupdata_raw{i_group},3);
        grid_N = sum(~isnan(groupdata_raw{i_group}),3)
        
        Y = 1:size(groupdata_raw{i_group},1);
        X = allTFs';
        
        l_nan_tfs = all(isnan(grid_average), 1);
        
        hs = surf(X(:,~l_nan_tfs), Y, flipud(grid_average(:,~l_nan_tfs)));
        hs.EdgeColor = plotcolors{i_group};
        hs.EdgeAlpha = 1;
        hs.LineWidth = 1.5;
        hs.FaceColor = plotcolors{i_group};
        hs.FaceAlpha = 0;
    end
    
    % plot the average smoothManifold
    if PLOT_SMOOTH_MANIFOLD
        grid_average = mean(groupdata_smooth{i_group}, 3);
        isi_ms = fliplr([1000/50 : 1 : 1000/10]);
        X = 1000./isi_ms;
        Y = 1:size(smoothManifold,1);
        hmod = surf(X,Y, flipud(grid_average));
        hmod.EdgeAlpha = .2;
        hmod.FaceColor = plotcolors{i_group};
        hmod.FaceAlpha = 0.5;
    end
    
end
set(gca, 'zscale', 'log', 'view', [-43    16])
set(gca, 'YTick', 1:10, 'YTickLabel', {'10','9','8','7','6','5','4','3','2','1'})
zmax = get(gca, 'zlim');
set(gca, 'XGrid', 'on', 'Ygrid', 'on', 'Zgrid', 'on')
set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'off', 'ZMinorGrid', 'off')
xlabel('Temporal Frequency')
ylabel('Pulse Number')
zlabel('norm amp')


% plot manifolds for each dataset individually
if PLOT_INDIVIDUAL_DATASETS
    for i_group = 1:numel(groupdata_raw)
        
        for i_examp = 1:size(groupdata_raw{i_group},3)
            f = figure;
            
            grid_average = groupdata_raw{i_group}(:,:,i_examp);
            
            Y = 1:size(groupdata_raw{i_group},1);
            X = allTFs';
            
            l_nan_tfs = all(isnan(grid_average), 1);
            
            hs = surf(X(:,~l_nan_tfs), Y, flipud(grid_average(:,~l_nan_tfs)));
            hs.EdgeColor = plotcolors{i_group};
            hs.EdgeAlpha = 1;
            hs.FaceColor = plotcolors{i_group};
            hs.FaceAlpha = 0;
            hs.LineWidth = 1.5;
            
            set(gca, 'zscale', 'log', 'view', [-43    16])
            set(gca, 'ytick', Y, 'yticklabel', cellfun(@(x) num2str(x), num2cell(rot90(Y,2)), 'uniformoutput', false))
            set(gca, 'xtick', X, 'xticklabel', cellfun(@(x) num2str(x), num2cell(allTFs), 'uniformoutput', false))
            xlabel('Temporal Frequency')
            ylabel('Pulse Number')
            zlabel('norm amp')
            
            i_ex = groupexpinds{i_group}(i_examp);
            i_ch = groupchinds{i_group}(i_examp);
            f.Name = sprintf('Mouse %s, site %s, HS%d', dat{i_ex}.info.mouseName, dat{i_ex}.info.siteNum, i_ch);
            
            % display the smooth manifold
            smoothManifold = pprpop.smoothManifold{i_ex}{i_ch};
            isi_ms = pprpop.smoothManifold_isi{i_ex}{i_ch};
            X = 1000./isi_ms;
            Y = 1:size(smoothManifold,1);
            hold on,
            hmod = surf(X,Y, flipud(smoothManifold));
            hmod.EdgeAlpha = 0;
            hmod.FaceAlpha = 0.5;
            
            
            
        end
        
    end
end


%% POST-TETANIC POTENTIATION (DATA COLLECTION)

% loop through the experiments. Pull out the trains data. Ignore the
% recovery train (if present) and aggregate across recovery conditions.
recovpop = [];
recovpop.TFsAllExpts = [];
recovpop.recoveryTimesAllExpts = [];
for i_ex = 1:numel(dat)
    
    % find the normal trains. Assume the field name is NOT 'ritv'
    condnames = fieldnames(dat{i_ex}.expt);
    l_trains = ~strncmp(condnames, 'RITv', 4);
    if sum(l_trains)==0; continue; end % no trains data
    trainParams = cellfun(@(x) dat{i_ex}.expt.(condnames{x}).tdict, num2cell(find(l_trains)), 'uniformoutput', false);
    trainParams = cat(1, trainParams{:});
    recovpop.trainParams{i_ex} = trainParams;
    
    % make sure the pulse amplitude and width were identical across ttypes
    assert(numel(unique(trainParams(:,1)))==1, 'ERROR: more than one pulse amplitude')
    assert(numel(unique(trainParams(:,2)))==1, 'ERROR: more than one pulse width');
    
    % store the stim params for all stim types, which will be useful for
    % indexing later.
    allStimParams = cellfun(@(x) dat{i_ex}.expt.(condnames{x}).tdict, num2cell((1:numel(l_trains))'), 'uniformoutput', false);
    allStimParams = cat(1, allStimParams{:});
    
    
    % aggregate data across TF conditons and recovery times (Separately for each
    % recording channel)
    
    for i_ch = 1:2
        
        % check to make sure this neuron was defined
        isvalid = dat{i_ex}.info.HS_is_valid_Vclamp(i_ch);
        
        for i_cond = 1:sum(l_trains)
            
            condIdx = ismember(allStimParams, trainParams(i_cond,:), 'rows');
            
            % check to make sure this is a recovery train
            isrecovery = allStimParams(condIdx,4) > 0;
            if ~isvalid || ~isrecovery
                recovpop.dat{i_ex}.recovAmp{i_ch}(i_cond,1) = NaN;
                continue
            end
            
            % pull out the data
            tmpdat = dat{i_ex}.expt.(condnames{condIdx}).stats.EPSCamp{i_ch}; % [Npulses, 1, Nsweeps]
            

            % check to make sure there are data for these conditions. Even
            % though this recording channel should be defined (See above),
            % it's possible there are no data due to deletion of single
            % sweeps
            if isempty(tmpdat)
                recovpop.dat{i_ex}.recovAmp{i_ch}(i_cond,1) = NaN;
                continue
            else
                
                % normalize by the p1Amp_norm
                realTrlNums = dat{i_ex}.expt.(condnames{condIdx}).realTrialNum{i_ch};
                p1Amp_norm = dat{i_ex}.qc.p1amp_norm{i_ch}(realTrlNums);
                assert(~any(isnan(p1Amp_norm)), 'ERROR: scale factor is a nan');
                p1Amp_norm = permute(p1Amp_norm, [1,3,2]);
                tmpdat = bsxfun(@rdivide, tmpdat, p1Amp_norm); % [Npulses, 1, Nsweeps]
                
                
                % now average across trails, and re-normalize to the first
                % pulse. Store in the population data structure.
                avg = mean(tmpdat,3);
                ppr = avg./avg(end-1);
                
                recovpop.dat{i_ex}.recovAmp{i_ch}(i_cond,1) = ppr(end);
                
            end
            
        end
        
    end
    
    %
    % update these things, but only if there were recovery data aggregated
    % for this experiment
    %
    dataPresent = cellfun(@(x) any(~isnan(x)), recovpop.dat{i_ex}.recovAmp);
    if any(dataPresent)
        
        % identify the unique TF conditions for this experiment, and update the
        % running log of TFs used across all experiments. This will be used by
        % the plotting routines to figure out the "canonical grid"
        uniqueTFs = unique(trainParams(:,3));
        tmp = cat(1, recovpop.TFsAllExpts, uniqueTFs);
        recovpop.TFsAllExpts = unique(tmp);
        
        % identify the unique recovery times for this experiment, and update the
        % running log of recov times used across all experiments. This will be used by
        % the plotting routines to figure out the "canonical grid"
        uniqueRecoveryTimes = unique(trainParams(:,4));
        tmp = cat(1, recovpop.recoveryTimesAllExpts, uniqueRecoveryTimes);
        recovpop.recoveryTimesAllExpts = unique(tmp);
        
    end
    
end



%% POST-TETANIC POTENTIATION (PLOTS)

% plot recovery pulse amplitude (normalized) vs. train frequency. Do this
% separately for cell types, brain areas, opsins, etc...
clc


% define a set of attributes for each lineseries (or manifold) in the plot
% {CellType,  BrainArea,  OpsinType}
plotgroups = {
              'PY',    'med', 'chronos';...
              'PY',    'lat', 'chronos';...
              };

groupdata_raw = repmat({[]}, 1, size(plotgroups, 1)); % should only have N cells, where N = size(plotgroups, 1). Each cell has a matrix with a cononicalGrid:
canonicalGrid = nan(numel(recovpop.TFsAllExpts), numel(recovpop.recoveryTimesAllExpts)); % need to get the maxNpulses into the recovpop struct
allTFs = recovpop.TFsAllExpts;
allRecoveryTimes = recovpop.recoveryTimesAllExpts';

% iterate over the experiments. For each recording channel, determine what
% the attributes are, and place the data in the correct ploting group.
for i_ex = 1:numel(dat)
    
    for i_ch = 1:2
        % check to make sure this neuron was defined
        isvalid = dat{i_ex}.info.HS_is_valid_Vclamp(i_ch);
        if ~isvalid
            continue
        end
        
        % check the attributes
        ch_attribs = {dat{i_ex}.info.cellType{i_ch}, dat{i_ex}.info.brainArea, dat{i_ex}.info.opsin};
        l_nan = cellfun(@(x) all(isnan(x)), ch_attribs);
        ch_attribs(l_nan) = cellfun(@num2str, ch_attribs(l_nan), 'uniformoutput', false);
        
        % is this expt and channel cooresponds to one of the plot_groups in
        % terms of cellType and opsin
        l_cellType_match = cellfun(@(x) ~isempty(regexpi(ch_attribs{1}, x)), plotgroups(:,1)) | strcmpi(plotgroups(:,1), 'any');
        l_opsinMatch = cellfun(@(x) ~isempty(regexpi(ch_attribs{3}, x)), plotgroups(:,3)) | strcmpi(plotgroups(:,3), 'any');
        
        
        % determine if this experiment has a brain area that corresponds to
        % one of the plot_groups. Start by adding a 'medial', 'lateral'
        % assignment to the brain area
        expt_area = ch_attribs{2};
        if ~isempty(regexp(expt_area, 'AM', 'once')) || ~isempty(regexp(expt_area, 'PM', 'once'))
            expt_area = [expt_area, ' med'];
        elseif ~isempty(regexp(expt_area, 'AL', 'once')) || ~isempty(regexp(expt_area, 'LM', 'once'))
            expt_area = [expt_area, ' lat'];
        end
        l_brainArea_match = cellfun(@(x) ~isempty(regexp(expt_area, x, 'once')), plotgroups(:,2)) | strcmpi(plotgroups(:,2), 'any');
        
        group_idx = sum([l_cellType_match, l_brainArea_match, l_opsinMatch], 2) == 3;
        assert(sum(group_idx)<=1, 'ERROR: found too many group indicies')
        if sum(group_idx) == 0; continue; end

        
        % add data to the appropriate group data array
        tmpgrid = canonicalGrid;
        for i_cond = 1:size(recovpop.trainParams{i_ex}, 1);
            grid_row_idx = allTFs == recovpop.trainParams{i_ex}(i_cond, 3);
            grid_col_idx = allRecoveryTimes == recovpop.trainParams{i_ex}(i_cond, 4);
            tmpgrid(grid_row_idx, grid_col_idx) = recovpop.dat{i_ex}.recovAmp{i_ch}(i_cond);
        end
        
        if all(isnan(tmpgrid(:))); continue; end
        
        % aggregate the data
        groupdata_raw{group_idx} = cat(3, groupdata_raw{group_idx}, tmpgrid);
        
    end
end




% Plot average lineseries.
hf = figure;
hold on,
plotcolors = [1,0,0;...
              0,0,1;...
              0,1,0];
          
for i_group = 1:numel(groupdata_raw)
    
    grid_average = nanmean(groupdata_raw{i_group},3);
    if isempty(grid_average); continue; end
    
    N = sum(~isnan(groupdata_raw{i_group}),3)
    grid_sem = nanstd(groupdata_raw{i_group},[],3) ./ sqrt(N);
    
    groupcolors = repmat(plotcolors(i_group,:), size(grid_average,1), 1);
    ramp = linspace(0,0.8, size(grid_average,1))';
    l_off = sum(groupcolors, 1)==0;
    groupcolors(:, l_off) = [ramp, ramp];
    
    % loop over TF conds
    for i_tf = 1:size(grid_average, 1)
        l_nan = isnan(grid_average(i_tf,:));
        tmp_amps = grid_average(i_tf, ~l_nan);
        tmp_sem = grid_sem(i_tf, ~l_nan);
        tmp_recovTimes = allRecoveryTimes(~l_nan);
        my_errorbar(tmp_recovTimes, tmp_amps, tmp_sem, '-', 'color', groupcolors(i_tf,:), 'linewidth', 3);
    end
    
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
            plot(tmp, 'g')
            chronos_std = cat(1, chronos_std, nanstd(tmp));
        else
            plot(tmp, 'r')
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

%% CONTROLS: POSITION OF LASER STIMULUS

for i_ex = 1:numel(dat)
   
    % add a picture of the slice, only once per figure
    f = figure;    
    mousename = dat{i_ex}.info.mouseName;
    sitenum = dat{i_ex}.info.siteNum;
    cd([GL_DATPATH, mousename, filesep, 'Other'])
    d = dir;
    d.name;
    
    photoprefix = [mousename, '_site', sitenum];
    FOUND_PHOTO = false;
    for i_photo = 1:numel(d)
        if strncmpi(d(i_photo).name, photoprefix, numel(photoprefix))
            img = double(imread(d(i_photo).name));
            iminfo = imfinfo(d(i_photo).name);
            maxpixval = prctile(img(:), [99.9]);
            normfact = min([maxpixval, 2^iminfo.BitDepth-1]);
            img = (img ./ normfact) .* (2^iminfo.BitDepth-1); % auto adjust LUT
            
            imshow(uint8(img));
            FOUND_PHOTO = true;
        end
    end
    f.Position = [600   118   741   532];
    hold on
    
    if ~FOUND_PHOTO
        close(f)
        continue
    end
    
    % add an icon for the stimulus locations
    uinput_xy = round(ginput(2));
    i_ch = 1;
    DONE = false;
    while ~DONE && i_ch<=2
        
        % look for a valid recording channel and define the coordinate
        % based off that.
        isvalid = dat{i_ex}.info.HS_is_valid_Vclamp(i_ch) || dat{i_ex}.info.HS_is_valid_Iclamp(i_ch);
        if isvalid
            hs_xy = dat{i_ex}.info.HS_xy_pos{i_ch};
            pixperum = pixPerMicron(size(img,1), size(img,2));
            hs_xy_correction = round(hs_xy .* pixperum); %now in pix
            
            % plot the HS point
            plot(uinput_xy(i_ch, 1), uinput_xy(i_ch, 2), 'ko', 'markerfacecolor', 'k', 'markersize', 8)
            
            % define the position of the new stim site and pia site based
            % off the uinput coordinates
            if isfield(dat{i_ex}.info, 'pia_xy_pos')
                pia_xy = round(dat{i_ex}.info.pia_xy_pos.*.80 .* pixperum) +  uinput_xy(i_ch,:);
                plot(pia_xy(1), pia_xy(2), 'o', 'markeredgecolor', 'r', 'markerfacecolor', 'r', 'markersize', 8)
            end
            stim_xy = round(dat{i_ex}.info.stim_xy_pos.*.80 .* pixperum) +  uinput_xy(i_ch,:);
            plot(stim_xy(1), stim_xy(2), 'o', 'markeredgecolor', 'y', 'markerfacecolor', 'y', 'markersize', 8)
            drawnow
        end
        
        i_ch = i_ch + 1;
        
    end
end


