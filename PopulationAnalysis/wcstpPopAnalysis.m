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


parfor i_ex_par = 1:Nexpts
    dat{i_ex_par} = wcstp_compile_data(attributes{i_ex_par}, hidx, params);
    close all
end
fprintf('All done importing data\n')

%% QULAITY CONTROL PLOTS

close all

for i_ex = 1:numel(dat)

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
    'PY', 'L23', 'med', 'any';...
    'PY', 'L23', 'lat', 'any';...
    };

groupdata.Rin_peak = repmat({[]}, 1, size(plotgroups, 1)); % should only have N cells, where N = size(plotgroups, 1). Each cell has a matrix with a cononicalGrid:
groupdata.Rin_asym = repmat({[]}, 1, size(plotgroups, 1));
% groupdata.Rin_vclamp = repmat({[]}, 1, size(plotgroups, 1)); % a place  holder for when this analysis is up and running
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
    end
end


%
% plot histograms of Input resistance measured two different ways
%
%%%%%%%%%%%%%%%%%%%%%%%%%
f=figure;
Ngroups = (numel(groupdata.Rin_peak)); % adding one for "summary" cdf fig
N_plt_rows = Ngroups+1;
groupcolors = lines(Ngroups);
f.Position = [451   187   968   665];
allNums = cat(1, groupdata.Rin_peak{:}, groupdata.Rin_asym{:});
edges = linspace(min(allNums)-1, max(allNums)+1, 30);
for i_group = 1:numel(groupdata.Rin_asym)
    
    %
    % current clamp, peak vals
    %%%%%%%%%%%%%%%%%%
    pltidx = sub2ind([2, N_plt_rows], 1, i_group);
    subplot(N_plt_rows, 2, pltidx)
    
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
    
    % CDF summary
    pltidx = sub2ind([2, N_plt_rows], 1, Ngroups+1);
    subplot(N_plt_rows, 2, pltidx), hold on,
    N = histcounts(groupdata.Rin_peak{i_group}, edges);
    N(end+1) = 0;
    cdf_vals = cumsum(N)./sum(N);
    stairs(edges, cdf_vals, '-', 'color', groupcolors(i_group,:))
    xlabel('Input R (MOhm)')
    ylabel('Proportion')
    
    %
    % current clamp asym vals
    %%%%%%%%%%%%%%%%%%%%%%%%
    pltidx = sub2ind([2, N_plt_rows], 2, i_group);
    subplot(N_plt_rows, 2, pltidx)
    
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
    
    % CDF summary
    pltidx = sub2ind([2, N_plt_rows], 2, Ngroups+1);
    subplot(N_plt_rows, 2, pltidx), hold on,
    N = histcounts(groupdata.Rin_asym{i_group}, edges);
    N(end+1) = 0;
    cdf_vals = cumsum(N)./sum(N);
    stairs(edges, cdf_vals, '-', 'color', groupcolors(i_group,:))
    xlabel('Input R (MOhm)')
    ylabel('Proportion')
end


%
% plot a histograms of Vrest and depth for each group
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=figure;
Ngroups = (numel(groupdata.Rin_peak));
groupcolors = lines(Ngroups);
f.Position = [466   464   835   440];
for i_group = 1:numel(groupdata.Vrest)
    
    % Vrest
    allNums = cat(1, groupdata.Vrest{:});
    edges = linspace(min(allNums)-10, max(allNums)+10, 30);
    
    pltidx = sub2ind([2, Ngroups], 1, i_group);
    subplot(Ngroups, 2, pltidx)
    
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
    
    pltidx = sub2ind([2, Ngroups], 2, i_group);
    subplot(Ngroups, 2, pltidx)
    
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
% scatter plots of cell depth vs. Rin, Vrest, tau
%
%%%%%%%%%%%%%%%%%%%%%%%%%
f=figure;
Ngroups = (numel(groupdata.Rin_asym));
groupcolors = lines(Ngroups);
f.Position = [49         305        1333         602];
for i_group = 1:numel(groupdata.Rin_asym)
    
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
    plot(X, Y, 'o', 'color', groupcolors(i_group,:), 'markerfacecolor', groupcolors(i_group,:), 'markersize', 5);
    shadedErrorBar(x_mod, y_mod, [CI_up ; CI_down], {'color', groupcolors(i_group,:)});
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
    plot(X, Y, 'o', 'color', groupcolors(i_group,:), 'markerfacecolor', groupcolors(i_group,:), 'markersize', 5)
    shadedErrorBar(x_mod, y_mod, [CI_up ; CI_down], {'color', groupcolors(i_group,:)});
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
    plot(X, Y, 'o', 'color', groupcolors(i_group,:), 'markerfacecolor', groupcolors(i_group,:), 'markersize', 5)
    shadedErrorBar(x_mod, y_mod, [CI_up ; CI_down], {'color', groupcolors(i_group,:)});
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
groupcolors = lines(Ngroups);
f.Position = [114 426 1667 456];
legtext = {};
leghand = [];

% I-V curves for peak values
subplot(1,3,1), hold on,
for i_group = 1:numel(groupdata.IVcurve_peak)
    for i_cell = 1:numel(groupdata.IVcurve_peak{i_group})
        X = groupdata.IVcurve_peak{i_group}{i_cell}(:,1);
        Y = groupdata.IVcurve_peak{i_group}{i_cell}(:,2);
        hp = plot(X, Y, '-', 'color', groupcolors(i_group,:));
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
    for i_cell = 1:numel(groupdata.IVcurve_asym{i_group})
        X = groupdata.IVcurve_asym{i_group}{i_cell}(:,1);
        Y = groupdata.IVcurve_asym{i_group}{i_cell}(:,2);
        plot(X, Y, '-', 'color', groupcolors(i_group,:))
    end
end
title('IV curve, asym vals')
xlabel('current injection (pA)')
ylabel('voltage response (mV)')

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


% compare Ih sag
legtext = {};
leghand = [];
f = figure;
f.Position = [384 353 1071 415];
subplot(1,2,1), hold on,
for i_group = 1:numel(groupdata.Ih_sag)
    N = numel(groupdata.Ih_sag{i_group});
    for i_ex = 1:N
        hp = plot(groupdata.Ih_Vm{i_group}{i_ex}, groupdata.Ih_sag{i_group}{i_ex}, 'o-', 'color', groupcolors(i_group,:));
    end
    legtext{i_group} = sprintf('%s, %s, %s', plotgroups{i_group,1},plotgroups{i_group, 2},plotgroups{i_group,3} );
    leghand(i_group) = hp;
end
xlabel('Membrane Potential At Sag Peak (mV)')
ylabel('Sag Amplitude (mV)')
legend(leghand, legtext, 'location', 'best')
legend boxoff

subplot(1,2,2), hold on,
for i_group = 1:numel(groupdata.Ih_asym)
    N = numel(groupdata.Ih_asym{i_group});
    for i_ex = 1:N
        sag_Vm = groupdata.Ih_Vm{i_group}{i_ex} - groupdata.Ih_sag{i_group}{i_ex};
        asym_Vm =  groupdata.Ih_asym{i_group}{i_ex};
        rebound = sag_Vm - asym_Vm;
        plot(groupdata.Ih_Vm{i_group}{i_ex}, rebound, 'o-', 'color', groupcolors(i_group,:))
    end
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

trl_cutoff_num = 10;

% define a set of attributes for each analysis group
% {CellType, Layer,  BrainArea,  OpsinType}
% Brain Area can be: 'AL', 'PM', 'AM', 'LM', 'any', 'med', 'lat'. CASE SENSITIVE
plotgroups = {
    'PY', 'L23', 'med', 'any';...
    'PY', 'L23', 'lat', 'any';...
    %'all_pv', 'any', 'any', 'any';...
    %'all_som', 'any', 'any', 'any';...
    };

groupdata = [];
groupdata.Vthresh = repmat({[]}, 1, size(plotgroups, 1)); % should only have N cells, where N = size(plotgroups, 1).
groupdata.Vreset = repmat({[]}, 1, size(plotgroups, 1));
groupdata.fi = repmat({[]}, 1, size(plotgroups, 1));

Ngroups = (numel(groupdata.Vthresh));
groupcolors = lines(Ngroups);

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
    
    % first spike
    tmp_dat = groupdata.Vthresh_ordered_first{i_group};
    N = sum(~isnan(tmp_dat),1);
    
    l_N_enough = N>= trl_cutoff_num;
    xbar = nanmean(tmp_dat(:,l_N_enough), 1);
    sem = stderr(tmp_dat(:,l_N_enough), 1);
    plt_pA = pA_cmd_bins(l_N_enough);
    
    hp = shadedErrorBar(plt_pA, xbar, sem, {'o-', 'color', groupcolors(i_group,:)});
    leghands(i_group) = hp.mainLine;
    legtext{i_group} =  sprintf('first sipke: %s, %s, %s', plotgroups{i_group,1},plotgroups{i_group, 2},plotgroups{i_group,3} );
    
    % last spike
    tmp_dat = groupdata.Vthresh_ordered_last{i_group};
    N = sum(~isnan(tmp_dat),1);
    
    l_N_enough = N>=trl_cutoff_num;
    xbar = nanmean(tmp_dat(:,l_N_enough), 1);
    sem = stderr(tmp_dat(:,l_N_enough), 1);
    plt_pA = pA_cmd_bins(l_N_enough);
    
    shadedErrorBar(plt_pA, xbar, sem, {'s--', 'color', groupcolors(i_group,:)});
    
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
    
    % first spike
    tmp_dat = groupdata.Vreset_ordered_first{i_group};
    N = sum(~isnan(tmp_dat),1);
    
    l_N_enough = N>= trl_cutoff_num;
    xbar = nanmean(tmp_dat(:,l_N_enough), 1);
    sem = stderr(tmp_dat(:,l_N_enough), 1);
    plt_pA = pA_cmd_bins(l_N_enough);
    
    hp = shadedErrorBar(plt_pA, xbar, sem, {'o-', 'color', groupcolors(i_group,:)});
    leghands(i_group) = hp.mainLine;
    legtext{i_group} =  sprintf('first sipke: %s, %s, %s', plotgroups{i_group,1},plotgroups{i_group, 2},plotgroups{i_group,3} );
   
    
    % last spike
    tmp_dat = groupdata.Vreset_ordered_last{i_group};
    N = sum(~isnan(tmp_dat),1);
    
    l_N_enough = N>=trl_cutoff_num;
    xbar = nanmean(tmp_dat(:,l_N_enough), 1);
    sem = stderr(tmp_dat(:,l_N_enough), 1);
    plt_pA = pA_cmd_bins(l_N_enough);
    
    shadedErrorBar(plt_pA, xbar, sem, {'s--', 'color', groupcolors(i_group,:)});
    
end
xlabel('Current injection (pA)')
ylabel('Avg AHP voltage (mV)')
legend(leghands, legtext, 'location', 'northwest')
legend boxoff

% third, plot the spike rates
figure, hold on,
legtext = {};
leghands = [];
for i_group = 1:Ngroups
    tmp_dat = groupdata.fi_ordered{i_group};
    N = sum(~isnan(tmp_dat), 1);
    
    l_N_enough = N>=trl_cutoff_num;
    xbar = nanmean(tmp_dat(:,l_N_enough), 1);
    sem = stderr(tmp_dat(:,l_N_enough), 1);
    plt_pA = pA_cmd_bins(l_N_enough);
    
    hp = shadedErrorBar(plt_pA, xbar, sem, {'-', 'color', groupcolors(i_group,:)});
    leghands(i_group) = hp.mainLine;
    legtext{i_group} =  sprintf('%s, %s, %s', plotgroups{i_group,1},plotgroups{i_group, 2},plotgroups{i_group,3} );
    
    
end
xlabel('Current injection (pA)')
ylabel('Avg Spike Rate (Hz)')
legend(leghands, legtext, 'location', 'northwest')
legend boxoff
%% PLOT THE RAW VCLAMP WAVEFORMS FOLLOWING EACH PULSE

close all


PLOT_ALL_TRIALS = true;
DEBUG_MEAN = true;
DEBUG_ALL = false;
NORM_TO_SMOOTH_P1 = false;
PLOT_RIT = true;


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
CELL_TYPE = 'PY_L23';

p1pop.opsin = {};
p1pop.avgwf_vclamp = [];
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
            
        end
        
    end    
end

f = figure; hold on,
clrs.chronos = 'b';
clrs.chief = 'r';
for opsin = {'chronos', 'chief'}
    l_opsin = cellfun(@(x) ~isempty(x), regexpi(p1pop.opsin, opsin));
    wf_avg = mean(p1pop.avgwf_vclamp(l_opsin,:), 1);
    wf_sem = stderr(p1pop.avgwf_vclamp(l_opsin, :), 1);
    shadedErrorBar(p1pop.tt.*1000, wf_avg, wf_sem, {'linewidth', 2, 'color', clrs.(opsin{1})});
end
xlabel('time (ms)')
legend('chronos', 'chief')
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
PULSE_NUM = 10; % PnP1

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
    'PY', 'L23', 'med', 'chief';...
    'PY', 'L23', 'lat', 'chief';...
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
groupcolors = lines(size(plotgroups,1));
for i_tf = 1:Ntfs
    subplot(Ntfs, 1, i_tf), hold on,
    for i_grp = 1:numel(groupdata_raw)
        tmp = groupdata_raw{i_grp}{i_tf};
        xbar = nanmean(tmp, 1);
        N = sum(~isnan(tmp),1);
        unique(N)
        sem = nanstd(tmp, [], 1) ./ sqrt(N);
        if PLOTERRBAR
            shadedErrorBar(1:size(tmp,2), xbar, sem, {'color', groupcolors(i_grp,:), 'linewidth', 2});
        else
            plot(tmp', '-', 'color', groupcolors(i_grp,:))
            %plot(mean(tmp,1), '-', 'color', groupcolors{i_grp}, 'linewidth', 3)
        end
    end
    ylim([0,2])
end




%% V-CLAMP POPULATION ANALYSIS (DATA COLLECTION)

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
close all; clc
% plot recovery pulse amplitude (normalized) vs. train frequency. Do this
% separately for cell types, brain areas, opsins, etc...

PLOT_RAW_DATA = false;
PLOT_AVG_DATA = true;
PLOT_MANIFOLD_OF_AVG_RAW_DATA = false;

% define a set of attributes for each analysis group
% {CellType, Layer,  BrainArea,  OpsinType}
% Brain Area can be: 'AL', 'PM', 'AM', 'LM', 'any', 'med', 'lat'. CASE SENSITIVE
plotgroups = {
    %'PY', 'L23', 'any', 'chief';...
    %'PY', 'L23', 'lat', 'chief';...
    %'PY', 'L23', 'AM', 'chronos';...
    %'PY', 'L23', 'any', 'chronos';...
    %'all_pv', 'any', 'any', 'any';...
    'NPVIN', 'any', 'any', 'any';...
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
groupdata.wfs = repmat({template}, size(plotgroups,1), 1); 
groupdata.amps = repmat({template}, size(plotgroups,1), 1);

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
        
        % was the the EPSC amp too small?
        if recovpop.dat{i_ex}.ignore{i_ch}
            continue
        end
        
        % check the attributes
        ch_attribs = {dat{i_ex}.info.cellType{i_ch}, dat{i_ex}.info.brainArea, dat{i_ex}.info.opsin};
        group_idx = groupMatcher(plotgroups, ch_attribs);
        if sum(group_idx) == 0; continue; end
        
        % pull out the waveform data
        %
        % Need to order things by ALL TFS and RECOVs, and not w/r/t an
        % individual data file's lattice
        %
        wf_top_level = recovpop.dat{i_ex}.psc_wfs{i_ch};
        amps_top_level = recovpop.dat{i_ex}.psc_amps{i_ch};
        ex_tfs = cell2mat(wf_top_level{1});
        %if any(ex_tfs==10);
        i_ex    
        group_idx
            ch_attribs
            
        %end
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
Ngroups = size(plotgroups,1);
plotcolors = lines(Ngroups);
legtext = repmat({{}}, 1, Ntfs);
leghand = repmat({[]}, 1, Ntfs);
for i_group = 1:size(plotgroups, 1)
    
    for i_tf = 1:Ntfs
        
        subplot(Ntfs, 1, i_tf), hold on,
        
        % check to make sure there are data for this TF condition
        data_exist = ~isempty(groupdata.wfs{i_group}{i_tf,1});
        if data_exist
            
            
            % concatenate the time series, padding with nans
            N_expts = cellfun(@(x) size(x,1), groupdata.wfs{i_group}(i_tf,:));
            N_expts = max(N_expts);
            all_wfs = cellfun(@(x) cat(2, x, nan(N_expts, 200)), groupdata.wfs{i_group}(i_tf,:), 'uniformoutput', false);
            all_wfs = cat(2, all_wfs{:});
            
            % determine the mean across experiments
            xbar = nanmean(all_wfs, 1);
            sem = stderr(all_wfs, 1);
            
            % mean +/- sem
            if PLOT_AVG_DATA
                hp = shadedErrorBar(1:numel(xbar), xbar, sem, {'color', plotcolors(i_group,:)});
                leghand{i_tf}(end+1) = hp.mainLine;
                legtext{i_tf}{end+1} = sprintf('%s, %s, %s, N=%d', plotgroups{i_group,1}, plotgroups{i_group,2}, plotgroups{i_group,3}, N_expts);
            end
            
            
            % all the data
            if PLOT_RAW_DATA
                plot(1:numel(xbar), all_wfs', '-', 'color', plotcolors(i_group,:))
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
Ngroups = size(plotgroups,1);
plotcolors = lines(Ngroups);
legtext = repmat({{}}, 1, Ntfs);
leghand = repmat({[]}, 1, Ntfs);
for i_group = 1:size(plotgroups, 1)
    
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
            recov_amps = cat(2, groupdata.amps{i_group}{i_tf,2:end});
            
            % determine the mean across experiments
            xbar_trains = nanmean(train_amps, 1);
            sem_trains = stderr(train_amps, 1);
            xx_trains = 1:numel(xbar_trains);
            xbar_recov = nanmean(recov_amps, 1);
            sem_recov = stderr(recov_amps, 1);
            xx_recov = numel(xbar_trains)+1 : numel(xbar_trains)+numel(xbar_recov);
            
            % plot
            hp = errorbar(xx_trains, xbar_trains, sem_trains, 'color', plotcolors(i_group,:));
            errorbar(xx_recov, xbar_recov, sem_recov, 'o', 'markerfacecolor', plotcolors(i_group,:), 'markeredgecolor', plotcolors(i_group,:))
            leghand{i_tf}(end+1) = hp(1);
            legtext{i_tf}{end+1} = sprintf('%s, %s, %s, N=%d', plotgroups{i_group,1}, plotgroups{i_group,2}, plotgroups{i_group,3}, N_expts);
            
            % add predictions from the fits to the averaged raw data
            if PLOT_MANIFOLD_OF_AVG_RAW_DATA
                % check that the plot_group defined here is the same as the
                % plot-group defined by the group_fit_to_grand_avg dataset
                data_matches = all(strcmpi(plotgroups(i_group,:), group_fit_to_grand_avg{i_group}.grouptypes));
                assert(data_matches, 'ERROR: plot groups do not match')
                
                smooth_manifold = group_fit_to_grand_avg{i_group}.smooth_manifold;
                
                train_tf = allTFs(i_tf);
                isi_ms = fliplr([1000/50 : 1 : 1000/10]);
                manifold_freqs = 1000./isi_ms;
                [~, freq_idx] = min(abs(manifold_freqs - train_tf));
                smooth_xx = 1:size(smooth_manifold, 1);
                smooth_yy = smooth_manifold(:, freq_idx);
                
                %plot the smooth manifold
                predmod = plot(smooth_xx,smooth_yy, '--', 'color', plotcolors(i_group,:));
            end
        end
        
    end
end
% indicate the TF, and recovery times, add a legend
figure(f)
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


%% ESTIMATE TIME CONSTANTS OF SHORT TERM PLASTICITY
tic
MODEL = 'ddff'; % string (eg 'ddff') or number of d and f terms (eg, 3) for FULLFACT
TRAINSET = 'all';  % could be 'rit', 'recovery', 'all'
PLOT_TRAINING_DATA = true;
FIT_RECOVERY_PULSE = true;
DATA_FOR_FIT = 'avg'; % can be 'avg', 'norm', 'raw'. 
CONVERT_TO_SMOOTH_P1 = true;
FORCE_RECOVERY = true;

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
            
            % force the model to fit a ficticious data point at +15 seconds
            % at a normalized value of 1
            if FORCE_RECOVERY
                if ~FIT_RECOVERY_PULSE; error('ERROR: need to fit recov for force_recovery'); end
                if ~CONVERT_TO_SMOOTH_P1; error('ERROR: need to normalize for force_recovery'); end
                pOnTimes(end+1:end+3) = [15, 17, 25]; % in seconds
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
            
            % calculate all the fitting errors
            num_params = numel(fit_results{i_mod}.params);
            fit_results{i_mod}.R2_train = get_r2(training_raw, training_pred);
            fit_results{i_mod}.R2_train_adj = get_r2(training_raw, training_pred, num_params);
            fit_results{i_mod}.AICc_train = get_aic(training_raw, training_pred, num_params);
            fit_results{i_mod}.MSE_train = mean((training_raw - training_pred).^2);
            fit_results{i_mod}.R2_crossvalid = NaN; % this needs to be fixed... get_r2(crossval_raw, crossval_pred);
            
            
            
            
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
                    title(sprintf('Best R2 = %.2f', max(cell2mat(structcat(fit_results, 'R2_train')))));
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
        dat{i_ex}.stpfits.fit_results{i_ch} = fit_results; % previously = [d, f, dTau, fTau], now has all models
        dat{i_ex}.stpfits.training_data{i_ch} = training_data;
        dat{i_ex}.stpfits.xval_data{i_ch} = xval_data;
        
    end
    drawnow
    
    
end
toc


%% STP FITTING PLOTS RE-DO

PLOT_TRAINING_DATA = true;

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
                    title(sprintf('Best R2 = %.2f', max(cell2mat(structcat(fit_results, 'R2_train')))));
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
            %xlims(1) = min([min(xx), xlims(1)]);
            %xlims(2) = max([max(xx), xlims(2)]);
            yvals = get(gca, 'ylim');
            ylims(1) = min([yvals(1), ylims(1)]);
            ylims(2) = max([yvals(2), ylims(2)]);
            xlim([0,7])
            
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
    'PY', 'L23', 'med', 'chief';...
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
plotcolors = {'r', 'b', 'g', 'k'};

figure, hold on,
for i_group = 1:size(plotgroups,1)
    tmp_dat = groupdata_fit_params{i_group};
    
    xbar = mean(tmp_dat,1);
    sem = stderr(tmp_dat,1);
    %p1 = plot(tmp_dat', '-', 'color', plotcolors{i_group}, 'linewidth', 0.5);
    %for ii = 1:numel(p1); p1(ii).Color(4) = 0.5; end
    %plot(xbar, '-', 'color', plotcolors{i_group}, 'linewidth', 4)
    shadedErrorBar(1:numel(MODEL)*2, xbar, sem, {'color', plotcolors{i_group}, 'linewidth', 4}, true);   
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
for i_ex = 1:numel(dat)
    
    % assume that the model was trained using the recovery trains data
    %assert(strcmpi(dat{i_ex}.stpfits.trainingSet, 'recovery'), 'ERROR: model needs to be fit with the recovery trains')
    
    % store some metadata
    pprpop.info{i_ex} = dat{i_ex}.info;
    
    % pre allocate the outputs in case a channel is not defined
    pprpop.dat{i_ex}.xbar = {[],[]};
    pprpop.tfs{i_ex} = {[],[]};
    pprpop.smoothManifold{i_ex} = {[],[]};
    pprpop.smoothManifold_isi{i_ex} = {[],[]};
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
        isi_ms = fliplr([1000/50 : 1 : 1000/10]);
        NumPulses = 10;
        
        % solve for P1 amp
        all_p1_amps = cat(3, pprpop.dat{i_ex}.xbar{i_ch}{:});
        A0 = mean(all_p1_amps(1,1,:), 3);
        
        smoothManifold = nan(NumPulses, numel(isi_ms));
        for i_isi = 1:numel(isi_ms)
            tmp_pOntimes_ms = 0 : isi_ms(i_isi) : (isi_ms(i_isi)*NumPulses)-1;
            tmp_pOntimes_sec = tmp_pOntimes_ms ./ 1000;
            smoothManifold(:,i_isi) = predict_vca_psc(tmp_pOntimes_sec, params(1:2), params(3:4), params(5:6), params(7:8), A0);
        end
        pprpop.smoothManifold{i_ex}{i_ch} = smoothManifold;
        pprpop.smoothManifold_isi{i_ex}{i_ch} = isi_ms;
        pprpop.params{i_ex}{i_ch} = params;
        
    end
    
end



%% PAIRED PULSE PLASTICITY MANIFOLDS (PLOTS)
clc;

% switches for first figure
PLOT_AVG_SMOOTH_MANIFOLD = true;
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
    'PY', 'L23', 'LM', 'chronos';...
    'PY', 'L23', 'med', 'chronos';...
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
        
        if PLOT_AVG_SMOOTH_MANIFOLD || PLOT_OVERLAY_ALL_MANIFOLDS || PLOT_MANIFOLD_OF_AVG_RAW_DATA
            groupdata_smooth{group_idx} = cat(3, groupdata_smooth{group_idx}, pprpop.smoothManifold{i_ex}{i_ch});
            group_params{group_idx} = cat(1, group_params{group_idx}, pprpop.params{i_ex}{i_ch});
        end
    end
end


% Plot manfold(s)
f = figure;
hold on,
plotcolors = lines(size(man_plotgrps, 1));
for i_group = 1:numel(groupdata_raw)
    
    if PLOT_AVG_RAW_DATA
        grid_average = nanmean(groupdata_raw{i_group},3);
        grid_average = grid_average(1:10, :);
        grid_N = sum(~isnan(groupdata_raw{i_group}),3)
        
        Y = 1:10;
        X = allTFs';
        
        l_nan_tfs = all(isnan(grid_average), 1);
        
        hs = surf(X(:,~l_nan_tfs), Y, flipud(grid_average(:,~l_nan_tfs)));
        hs.EdgeColor = plotcolors(i_group,:);
        hs.EdgeAlpha = 1;
        hs.LineWidth = 1.5;
        hs.FaceColor = plotcolors(i_group,:);
        hs.FaceAlpha = 0.3;
    end
    
    % plot the average smoothManifold
    isi_ms = fliplr([1000/50 : 1 : 1000/10]);
    X = 1000./isi_ms;
    N_pulses_smooth = size(groupdata_smooth{i_group}, 1);
    Y = 1:N_pulses_smooth;
    
    if PLOT_AVG_SMOOTH_MANIFOLD
        grid_average = nanmean(groupdata_smooth{i_group}, 3);
        hmod = surf(X,Y, flipud(grid_average));
        hmod.EdgeAlpha = 0;
        hmod.FaceColor = plotcolors(i_group,:);
        hmod.FaceAlpha = 0.6;
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
        predmod.FaceColor = plotcolors(i_group,:);
        predmod.FaceAlpha = 0.2;
    end
    
    if PLOT_OVERLAY_ALL_MANIFOLDS
        N_manifolds = size(groupdata_smooth{i_group}, 3);
        for i_examp = 1:N_manifolds
            hmod = surf(X,Y, flipud(groupdata_smooth{i_group}(:,:,i_examp)));
            hmod.EdgeAlpha = .2;
            hmod.FaceColor = plotcolors(i_group,:);
            hmod.FaceAlpha = 0.5;
        end
    end
    
end
%set(gca, 'zlim', [0.45, 1.4])
set(gca, 'zscale', 'linear', 'view', [-43    16])
set(gca, 'YTick', 1:10, 'YTickLabel', {'10','9','8','7','6','5','4','3','2','1'})
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
        % remake the smooth manifold
        isi_ms = fliplr([1000/50 : 1 : 1000/10]);
        X = 1000./isi_ms;
        N_pulses_smooth = size(groupdata_smooth{i_group}, 1);
        Y = 1:N_pulses_smooth;
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
        predmod.EdgeColor = plotcolors(i_group,:);
        predmod.FaceColor = plotcolors(i_group,:);
        predmod.FaceAlpha = 0;
        
    end
end

if PLOT_MANIFOLD_OF_AVG_RAW_DATA
    for i_group = 1:size(man_plotgrps,1)
        
        % check that the plot_group defined here is the same as the
        % plot-group defined by the group_fit_to_grand_avg dataset
        data_matches = all(strcmpi(man_plotgrps(i_group,:), group_fit_to_grand_avg{i_group}.grouptypes));
        assert(data_matches, 'ERROR: plot groups do not match')
        
        isi_ms = fliplr([1000/50 : 1 : 1000/10]);
        X = 1000./isi_ms;
        N_pulses_smooth = size(groupdata_smooth{i_group}, 1);
        Y = 1:N_pulses_smooth;
        smoothManifold = group_fit_to_grand_avg{i_group}.smooth_manifold;
        
        %plot the smooth manifold
        predmod = surf(X,Y, flipud(smoothManifold));
        predmod.EdgeAlpha = 0.5;
        predmod.EdgeColor = 'k';
        predmod.FaceColor = 'k';
        predmod.FaceAlpha = 0;
        
    end
end

%
% plot manifolds for each dataset individually
% (lots of figures)
%%%%%%%
if PLOT_DATASETS_INDIVIDUALLY
    for i_group = 1:numel(groupdata_raw)
        
        for i_examp = 1:size(groupdata_raw{i_group},3)
            f = figure;
            
            grid_average = groupdata_raw{i_group}(:,:,i_examp);
            
            Y = 1:size(groupdata_raw{i_group},1);
            X = allTFs';
            
            l_nan_tfs = all(isnan(grid_average), 1);
            
            hs = surf(X(:,~l_nan_tfs), Y, flipud(grid_average(:,~l_nan_tfs)));
            hs.EdgeColor = plotcolors(i_group,:)
            hs.EdgeAlpha = 1;
            hs.FaceColor = plotcolors(i_group,:);
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
end

%%  FIT THE "GRAND AVERAGE" MANIFOLDS
% the manifold of the average params is dead wrong. Instead, find the 
% best fitting params of the "average" manifold.
grand_avg_params={};
for i_group = 1:size(man_plotgrps,1)
    
    % calculate the average across manifolds
    avg_manifold = mean(groupdata_smooth{i_group}, 3);
    
    % fit the grand-average manifold with the VCA model
    isi_ms = fliplr([1000/50 : 1 : 1000/10]);
    X = 1000./isi_ms;
    N_pulses_smooth = size(groupdata_smooth{i_group}, 1);
    Y = 1:N_pulses_smooth;
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
    isi_ms = fliplr([1000/50 : 1 : 1000/10]);
    X = 1000./isi_ms;
    N_pulses_smooth = size(groupdata_smooth{i_group}, 1);
    Y = 1:N_pulses_smooth;
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
writetable(T_fit_to_avg_raw, 'params_from_fit_to_grand_avg_raw_withRecovAndRIT_forcedRecov.csv')


%% STRENGTH OF INPUT ONTO INs

% define a set of attributes for each analysis group. Cell_type_pairs will
% be split according to opsin, layer, and/or brain areas
PLOTVAL = 'p1_ratio'; % 'p1_ratio', 'num_val', 'denom_val' 
cell_type_pairs = {'all_som', 'PY';...  % defines the different groups that will appear on the x-axis
                   'all_pv',  'PY';...
                   };
grp_opsins = {'any_opsin'}; % need to be encapsulated in a cell array {'any_opsin'}
grp_layers = {'L23'};
grp_brain_areas = {'med', 'lat'}; % case sensitive, 'any_hva', 'med', 'lat', 'AL', 'LM' etc..

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


