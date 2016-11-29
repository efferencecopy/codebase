%% SPECIFY WHICH EXPERIMENTS SHOULD CONTRIBUTE, LOAD THE DATA

fin

GL_SAVE_FIGURE = false;
GL_FIG_DIRECTORY = [GL_DATPATH(1:end-11), filesep, 'wcstp_figures'];


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
    pool = parpool(25);
end


parfor i_ex_par = 1:Nexpts
    dat{i_ex_par} = wcstp_compile_data(attributes{i_ex_par}, hidx, params);
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
        if ~all(cellfun(@isempty, dat{i_ex}.dcsteps.Vm_raw))
            Vm = dat{i_ex}.dcsteps.Vm_raw{i_ch};
            Icmd = dat{i_ex}.dcsteps.Icmd{i_ch};
            N = size(Vm,2);
            tt = ([0:N-1] ./ dat{i_ex}.info.sampRate.dcsteps) - params.pretime.dcsteps;
            if i_ch == 1; col = 1; else col = 3; end
            
            if any(Icmd<0)
                pltidx = sub2ind([4,3], col, 3);
                ha = subplot(3,4,pltidx); hold on,
                neg_cmds = unique(Icmd(Icmd < 0));
                for i_cmd = 1:numel(neg_cmds)
                    l_cmd = Icmd == neg_cmds(i_cmd);
                    plot(tt, mean(Vm(l_cmd,:), 1))
                end
                if isfield(dat{i_ex}.dcsteps, 'Ih_sag')
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
    
    if GL_SAVE_FIGURE
        
        % cd to the correct directory
        cd(GL_FIG_DIRECTORY)
        d = dir;
        mouse_dir_exists = any(strcmpi(structcat(d, 'name'), dat{i_ex}.info.mouseName));
        if ~mouse_dir_exists
            mkdir('./', dat{i_ex}.info.mouseName)
        end
        cd([GL_FIG_DIRECTORY, filesep, dat{i_ex}.info.mouseName])
        
        % save the figure
        print('quality_control.eps', '-depsc')
        
    end
    
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
% plot histograms of Input resistance measured three different ways
%
%%%%%%%%%%%%%%%%%%%%%%%%%
f=figure;
Ngroups = (numel(groupdata.Rin_peak)) + 1; % adding one for "summary" cdf fig
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


% compare Ih sag
f = figure;
f.Position = [384 353 1071 415];
subplot(1,2,1), hold on,
for i_group = 1:numel(groupdata.Ih_sag)
    N = numel(groupdata.Ih_sag{i_group});
    for i_ex = 1:N
        plot(groupdata.Ih_Vm{i_group}{i_ex}, groupdata.Ih_sag{i_group}{i_ex}, 'o-', 'color', groupcolors(i_group,:))
    end
end
xlabel('Membrane Potential At Sag Peak (mV)')
ylabel('Sag Amplitude (mV)')

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

%% PLOT THE RAW VCLAMP WAVEFORMS FOLLOWING EACH PULSE

close all


PLOT_ALL_TRIALS = true;
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
            case {'chronos', 'chronos_flx', 'chronos_gfp', 'chronos_tom'}
                pltclr = 'g';
            otherwise
                pltclr = nan;
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


%% I-CLAMP POPULATION ANALYSIS (PLOTTING)

clc; close all

PLOTERRBAR = true;

% define a set of attributes for each analysis group
% {CellType, Layer,  BrainArea,  OpsinType}
% Brain Area can be: 'AL', 'PM', 'AM', 'LM', 'any', 'med', 'lat'. CASE SENSITIVE
plotgroups = {
    'PY', 'L23', 'AM', 'chief';...
    'PY', 'L23', 'LM', 'chief';...
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
        N = sum(~isnan(tmp),1);
        unique(N)
        sem = nanstd(tmp, [], 1) ./ sqrt(N);
        if PLOTERRBAR
            shadedErrorBar(1:size(tmp,2), xbar, sem, {'color', groupcolors{i_grp}, 'linewidth', 2});
        else
            plot(tmp', '-', 'color', groupcolors{i_grp})
            %plot(mean(tmp,1), '-', 'color', groupcolors{i_grp}, 'linewidth', 3)
        end
        
    end
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
    Ntfs = numel(unique_tfs);
    Nrecov = numel(unique_recov);
    for i_ch = 1:2
        
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

% plot recovery pulse amplitude (normalized) vs. train frequency. Do this
% separately for cell types, brain areas, opsins, etc...
clc
PLOT_RAW_DATA = false;
PLOT_AVG_DATA = true;

% define a set of attributes for each analysis group
% {CellType, Layer,  BrainArea,  OpsinType}
% Brain Area can be: 'AL', 'PM', 'AM', 'LM', 'any', 'med', 'lat'. CASE SENSITIVE
plotgroups = {
    %'all_pv', 'any', 'any', 'chief';...
    %'all_som', 'any', 'any', 'chief';...
    %'all_pv', 'any', 'any', 'chief';...
    'PY', 'L23', 'AM', 'chief';...
    'PY', 'L23', 'LM', 'chief';...
    %'PY', 'L23', 'LM', 'chief';...
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
groupdata.wfs = repmat({template}, size(plotgroups,1), 1); 

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
        ex_tfs = cell2mat(wf_top_level{1});
        wfs_train = cellfun(@(x) mean(x, 3), wf_top_level{2}, 'uniformoutput', false); % average across replicates
        for i_tf = 1:numel(ex_tfs)
            % grab the data from the trains
            tf_wf = wfs_train{i_tf};
            tf_wf = cat(2, tf_wf, nan(size(tf_wf, 1), 25)); % add some nans for separation
            tf_wf = reshape(tf_wf', 1, []); % a row vector
            
            % push the data into the correct slot in the groupdata
            row_idx = allTFs == ex_tfs(i_tf);
            tmp_group = groupdata.wfs{group_idx}{row_idx, 1};
            tmp_group = cat(1, tmp_group, tf_wf); % Ntraces x Ntime
            groupdata.wfs{group_idx}{row_idx, 1} = tmp_group;
            
            % loop over the recovery conditions and deal with those
            ex_recovs = cell2mat(wf_top_level{3}{i_tf, 1});
            wfs_recovs = cellfun(@(x) mean(x, 3), wf_top_level{3}{i_tf,2}, 'uniformoutput', false); % average across replicates
            for i_recov = 1:numel(ex_recovs)
                col_idx = [false, allRecoveryTimes == ex_recovs(i_recov)]; % leading False to offset the "trains" position
                tmp_group = groupdata.wfs{group_idx}{row_idx, col_idx};
                tmp_group = cat(1, tmp_group, wfs_recovs{i_recov});
                groupdata.wfs{group_idx}{row_idx, col_idx} = tmp_group;                
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
plotcolors = {'r', 'b', 'g'};
for i_group = 1:size(plotgroups, 1)
    
    for i_tf = 1:Ntfs
        
        subplot(Ntfs, 1, i_tf), hold on,
        
        % check to make sure there are data for this TF condition
        data_exist = ~isempty(groupdata.wfs{i_group}{i_tf,1});
        
        if data_exist
            
            
            % concatenate the time series, padding with nans
            N_expts = cellfun(@(x) size(x,1), groupdata.wfs{i_group}(i_tf,:));
            N_expts = max(N_expts)
            all_wfs = cellfun(@(x) cat(2, x, nan(N_expts, 200)), groupdata.wfs{i_group}(i_tf,:), 'uniformoutput', false);
            all_wfs = cat(2, all_wfs{:});
            
            % determine the mean across experiments
            xbar = nanmean(all_wfs, 1);
            sem = stderr(all_wfs, 1);
            
            % mean +/- sem
            if PLOT_AVG_DATA
                shadedErrorBar(1:numel(xbar), xbar, sem, {'color', plotcolors{i_group}});
            end
            
            
            % all the data
            if PLOT_RAW_DATA
                plot(1:numel(xbar), all_wfs', '-', 'color', plotcolors{i_group})
            end
            
        end
        
    end
end
% indicate the TF, and recovery times
for i_tf = 1:Ntfs
    subplot(Ntfs, 1, i_tf)
    ylim([-1.5, 0.2])
    plt_train_tf = allTFs(i_tf);
    l_recov_exist = cellfun(@(x) ~isempty(x), groupdata.wfs{i_group}(i_tf,:));
    l_recov_exist(1) = []; % ignore the first array b/c it's not recov pulse
    plt_recov_ms = allRecoveryTimes(l_recov_exist);
    plt_string = sprintf('Train = %d Hz, Recovery (in ms): %s', plt_train_tf, num2str(plt_recov_ms));
    txt_hand = text(750, 0.2, plt_string);
    txt_hand.FontSize = 9;
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
           
           % delete the recovery pulse if need be
           if ~FITRECOVERYPULSE
               has_recov_pulse = dat{i_ex}.expt.(condnames{i_cond}).tdict(4) > 0;
               if has_recov_pulse
                   pOnTimes{i_cond}(end) = [];
                   rawAmps{i_cond} = rawAmps{i_cond}(1:end-1,1,:);
               end
           end
           
           % normalize to the p1amp_norm if need be
           if NORMALIZEDATA
               % norm to first pulse
               rawAmps{i_cond} = bsxfun(@rdivide, rawAmps{i_cond}, rawAmps{i_cond}(1,:,:));
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
                   trlNums = dat{i_ex}.expt.(condnames{i_cond}).realTrialNum{i_ch};
                   smooth_p1 = dat{i_ex}.qc.p1amp_norm{i_ch}(trlNums);
                   smooth_p1 = permute(smooth_p1, [1,3,2]);
                   p1Amps{i_cond} = squeeze(smooth_p1);
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
%        if [isfield(dat{i_ex}, 'stpfits') ...
%                && isfield(dat{i_ex}.stpfits, 'modelParams') ...
%                && numel(dat{i_ex}.stpfits.modelParams)>= i_ch ...
%                && ~isempty(dat{i_ex}.stpfits.modelParams{i_ch})] % end of if conditional
%            
%            % if you've made it here, then the old params are potentially
%            % still appliciable, but I should check to make sure that they
%            % were fit using the same technique
%            matches = [strcmpi(dat{i_ex}.stpfits.trainingSet, TRAINSET); ...
%                       dat{i_ex}.stpfits.normalizeData == NORMALIZEDATA; ...
%                       dat{i_ex}.stpfits.fitWithAvg == FITAVERAGEDATA; ...
%                       dat{i_ex}.stpfits.fitRecovPulse == FITRECOVERYPULSE];
%            
%           if all(matches)
%               d = dat{i_ex}.stpfits.modelParams{i_ch}(1:2);
%               f = dat{i_ex}.stpfits.modelParams{i_ch}(3);
%               dTau = dat{i_ex}.stpfits.modelParams{i_ch}(4:5);
%               fTau = dat{i_ex}.stpfits.modelParams{i_ch}(6);
% 
%               MAKENEWFITS = false;
%           end
%        end
       if MAKENEWFITS
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
           axis tight
       end
       if ~isempty(hs) && sum(hs)>0
           %set(hs(hs~=0), 'XLim', xlims, 'YLim', ylims)
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


%% PAIRED PULSE PLASTICITY MANIFOLDS (DATA COLLECTION)

% loop through the experiments. Pull out the trains data. Ignore the
% recovery train (if present) and aggregate across recovery conditions.
pprpop = [];
pprpop.TFsAllExpts = [];
pprpop.MaxNPulses = 0;
for i_ex = 1:numel(dat)
    
    %if i_ex == 48 || i_ex == 50; continue; end
    
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
                pprpop.smoothManifold{i_ex}{i_ch} = [];
                pprpop.smoothManifold_isi{i_ex}{i_ch} = [];
                pprpop.R2{i_ex}{i_ch} = [];
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
                pprpop.smoothManifold{i_ex}{i_ch} = [];
                pprpop.smoothManifold_isi{i_ex}{i_ch} = [];
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
                    pprpop.R2{i_ex}{i_ch} = dat{i_ex}.stpfits.R2.training{i_ch};
                end
                
            end
        end
        
    end
    
end


%% PAIRED PULSE PLASTICITY MANIFOLDS (PLOTS)
clc;

% switches for first figure
PLOT_AVG_SMOOTH_MANIFOLD = false;
PLOT_OVERLAY_ALL_MANIFOLDS = false;
PLOT_AVG_RAW_DATA = true;

% add other plots?
PLOT_DATASETS_INDIVIDUALLY = false;

% define a set of attributes for each analysis group
% {CellType, Layer,  BrainArea,  OpsinType}
% Brain Area can be: 'AL', 'PM', 'AM', 'LM', 'any', 'med', 'lat'. CASE SENSITIVE
plotgroups = {
    'PY', 'L23', 'AM', 'chief';...
    'PY', 'L23', 'AL', 'chief';...
    };

groupdata_raw = repmat({[]}, 1, size(plotgroups, 1)); % should only have N cells, where N = size(plotgroups, 1). Each cell has a matrix with a cononicalGrid:
groupdata_smooth =  repmat({[]}, 1, size(plotgroups, 1));
groupexpinds =  repmat({[]}, 1, size(plotgroups, 1));
groupchinds = repmat({[]}, 1, size(plotgroups, 1));
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
        group_idx = groupMatcher(plotgroups, ch_attribs);
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
        
        if PLOT_AVG_SMOOTH_MANIFOLD || PLOT_OVERLAY_ALL_MANIFOLDS
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
    
    if PLOT_AVG_RAW_DATA
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
        hmod.EdgeAlpha = .2;
        hmod.FaceColor = plotcolors{i_group};
        hmod.FaceAlpha = 0.5;
    end
    
    if PLOT_OVERLAY_ALL_MANIFOLDS
        N_manifolds = size(groupdata_smooth{i_group}, 3);
        for i_examp = 1:N_manifolds
            hmod = surf(X,Y, flipud(groupdata_smooth{i_group}(:,:,i_examp)));
            hmod.EdgeAlpha = .2;
            hmod.FaceColor = plotcolors{i_group};
            hmod.FaceAlpha = 0.5;
        end
    end
    
end
set(gca, 'zscale', 'linear', 'view', [-43    16], 'zlim', [0.45, 1.4])
set(gca, 'YTick', 1:10, 'YTickLabel', {'10','9','8','7','6','5','4','3','2','1'})
zmax = get(gca, 'zlim');
set(gca, 'XGrid', 'on', 'Ygrid', 'on', 'Zgrid', 'on')
set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'off', 'ZMinorGrid', 'off')
xlabel('Temporal Frequency')
ylabel('Pulse Number')
zlabel('norm amp')


% plot manifolds for each dataset individually
if PLOT_DATASETS_INDIVIDUALLY
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


