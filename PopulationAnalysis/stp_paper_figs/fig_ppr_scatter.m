function fig_ppr_scatter(dat, plotgroups)


MIN_TRL_COUNT = 0; % set to 0 to turn off
MIN_TFS_COUNT = 0; % set to 0 to turn off
XP = 3;
YP = 10;


% loop through the experiments. Pull out the trains data. Ignore the
% recovery train (if present) and aggregate across recovery conditions.
recovpop = [];
recovpop.TFsAllExpts = [];
for i_ex = 1:numel(dat)
    
    recovpop.dat{i_ex}.psc_amps = {}; % will have one cell or each ch
    recovpop.dat{i_ex}.ignore={[],[]};
    
    % check the attributes
    ch_attribs = {'PY_L23', dat{i_ex}.info.brainArea, dat{i_ex}.info.opsin};
    group_idx = groupMatcher(plotgroups, ch_attribs);
    if sum(group_idx) == 0; continue; end
    
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
    if isempty(unique_recov); continue; end % leaving this in for consistency with the other figuress
    
    % make sure the pulse amplitude and width were identical across ttypes
    assert(numel(unique(trainParams(:,1)))==1, 'ERROR: more than one pulse amplitude')
    assert(numel(unique(trainParams(:,2)))==1, 'ERROR: more than one pulse width');
    
    % aggregate data across TF conditons and recovery times (Separately for each
    % recording channel)
    Ntfs = numel(unique_tfs);
    Nrecov = numel(unique_recov);
    for i_ch = 1:2
        
        % remove instances where there is insufficient data to constrain
        % the fit
        insufficient_data = size(dat{i_ex}.qc.p1amp{i_ch}, 3) < MIN_TRL_COUNT;
        insufficient_tfs = Ntfs < MIN_TFS_COUNT;
        if insufficient_data || insufficient_tfs
            recovpop.dat{i_ex}.ignore{i_ch} = true;
        end
        
        
        % initialize the outputs. I need to concatenate similar TFs
        % together, but keep each of their recovery conditions separate
        %
        % recovpop.dat.psc_amps{i_ch} = {{TF}, {psc_amps}}
        recovpop.dat{i_ex}.psc_amps{i_ch} = {num2cell(unique_tfs), repmat({[]}, Ntfs,1)};
         
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
                % first pulse (i.e., deal with non-stationarity). Mean
                % across trials comes later...
                tmp_psc = dat{i_ex}.expt.(condnames_trains{cond_idx}).stats.EPSCamp{i_ch};
                real_trl_nums = dat{i_ex}.expt.(condnames_trains{cond_idx}).realTrialNum{i_ch};
                p1_norm_vals = dat{i_ex}.qc.p1amp_norm{i_ch}(real_trl_nums);
                tmp_psc = bsxfun(@rdivide, tmp_psc, permute(p1_norm_vals, [3,1,2]));
                
                % concatenate 1:10 into 'trains', and the 11th into 'recov'
                recovpop.dat{i_ex}.psc_amps{i_ch}{2}{i_tf} = cat(3, recovpop.dat{i_ex}.psc_amps{i_ch}{2}{i_tf}, tmp_psc(1:10,:,:));
            end
        end
    end
    
    
    
    % identify the unique TF conditions for this experiment, and update the
    % running log of TFs used across all experiments. This will be used by
    % the plotting routines to figure out the "canonical grid"
    tmp = cat(1, recovpop.TFsAllExpts, unique_tfs);
    recovpop.TFsAllExpts = unique(tmp);
    
end



% _________ PLOTS __________________
% plot recovery pulse amplitude (normalized) vs. train frequency. Do this
% separately for cell types, brain areas, opsins, etc...
allTFs = recovpop.TFsAllExpts;
Ntfs = numel(allTFs);

% set up the output arrays
template = repmat({[]}, Ntfs, 1);
groupdata = [];
groupdata.amps = repmat({template}, size(plotgroups,1), 1);

% iterate over the experiments. For each recording channel, determine what
% the attributes are, and place the data in the correct ploting group.
for i_ex = 1:numel(dat)
    
    % skip experiments that have no data (i.e., trains but no recovery)
    if isempty(recovpop.dat{i_ex});
        keyboard;
    end
    
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
        group_idx = groupMatcher(plotgroups, ch_attribs);
        if sum(group_idx) == 0; continue; end
        
        % pull out the data
        amps_top_level = recovpop.dat{i_ex}.psc_amps{i_ch};
        ex_tfs = cell2mat(amps_top_level{1});
        
        amps_train = cellfun(@(x) mean(x, 3), amps_top_level{2}, 'uniformoutput', false); % average across replicates
        amps_train = cellfun(@(x) x./x(1), amps_train, 'uniformoutput', false); % norm to P1
        for i_tf = 1:numel(ex_tfs)
            
            % push the amps data into the correct slot in the groupdata
            tf_amps = amps_train{i_tf};
            row_idx = allTFs == ex_tfs(i_tf);
            tmp_group = groupdata.amps{group_idx}{row_idx, 1};
            tmp_group = cat(1, tmp_group, tf_amps'); % Nneurons x Npulses
            groupdata.amps{group_idx}{row_idx, 1} = tmp_group;
            
        end
        
    end
end



%
%  Scater plots
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = figure; hold on,
f.Units = 'Normalized';
Ngroups = size(plotgroups,1);
facealpha = linspace(1, 0.5, Ntfs);
x_vals = []; % xbar and sem for each tf
y_vals = []; % xbar and sem for each tf
p2_vals = []; % xbar for pulse 2 data
for i_group = 1:Ngroups
    plt_clr = hvaPlotColor(plotgroups{i_group, 3});
    
    for i_tf = 1:Ntfs
        pprs = groupdata.amps{i_group}{i_tf,1};
        
        X = pprs(:, XP);
        Y = pprs(:, YP);
        Z = pprs(:, 2);
        xvals(i_group, :, i_tf) = [mean(X), stderr(X)];
        yvals(i_group, :, i_tf) = [mean(Y), stderr(Y)];
        p2_vals(i_group, :, i_tf) = [mean(Z), stderr(Z)];
        scatter(X, Y, 40, 'o', 'markeredgecolor', 'w',...
                               'markerfacecolor', plt_clr,...
                               'markerfacealpha', facealpha(i_tf))
    end
end
xlabel(sprintf('Pulse %d (normalized)', XP))
ylabel(sprintf('Pulse %d (normalized)', YP))

% scatter plot of just the means
figure, hold on,
for i_group = 1:Ngroups
    
    plt_clr = hvaPlotColor(plotgroups{i_group, 3});
    
    x_xbar = permute(xvals(i_group, 1, :), [3,1,2]);
    x_sem = permute(xvals(i_group, 2, :), [3,1,2]);
    
    y_xbar = permute(yvals(i_group, 1, :), [3,1,2]);
    y_sem = permute(yvals(i_group, 2, :), [3,1,2]);
    
    plot(x_xbar, y_xbar, '-', 'color', plt_clr, 'linewidth', 2)
    
    for i_tf = 1:Ntfs
        xerr = [x_xbar(i_tf)-x_sem(i_tf), x_xbar(i_tf)+x_sem(i_tf)];
        yerr = [y_xbar(i_tf)-y_sem(i_tf), y_xbar(i_tf)+y_sem(i_tf)];
        plot(xerr, [y_xbar(i_tf), y_xbar(i_tf)], '-', 'color', plt_clr, 'linewidth', 2)
        plot([x_xbar(i_tf), x_xbar(i_tf)], yerr, '-', 'color', plt_clr, 'linewidth', 2)
        scatter(x_xbar(i_tf), y_xbar(i_tf), 200, 'markeredgecolor', 'w',...
                               'markerfacecolor', plt_clr,...
                               'markerfacealpha', facealpha(i_tf))
    end
    
%     % if this is for the lateral HVA, re-plot after shifting so that P2s
%     % are equal
%     if strcmpi(plotgroups{i_group, 3}, 'lat')
%         p2_vals_lat = permute(p2_vals(i_group, 1, :), [3,1,2]);
%         idx_med = strcmpi(plotgroups(:, 3), 'med');
%         p2_vals_med = permute(p2_vals(idx_med, 1, :), [3,1,2]);
%         new_x_xbar = x_xbar - p2_vals_lat + p2_vals_med;
%         new_y_xbar = y_xbar - p2_vals_lat + p2_vals_med;
%         plot(new_x_xbar, new_y_xbar, '--', 'color', plt_clr, 'linewidth', 2)
%     end
    
    
end
xlabel(sprintf('Pulse %d (normalized)', XP))
ylabel(sprintf('Pulse %d (normalized)', YP))
