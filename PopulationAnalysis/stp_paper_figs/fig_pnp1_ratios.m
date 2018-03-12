function fig_pnp1_ratios(dat, plotgroups)


MIN_TRL_COUNT = 0; % set to 0 to turn off
MIN_TFS_COUNT = 0; % set to 0 to turn off
MIN_RECOV_COUNT = 0; % set to 0 to turn off


% loop through the experiments. Pull out the trains data. Ignore the
% recovery train (if present) and aggregate across recovery conditions.
recovpop = [];
recovpop.TFsAllExpts = [];
recovpop.recoveryTimesAllExpts = [];
for i_ex = 1:numel(dat)
    
    recovpop.dat{i_ex}.psc_amps = {};
    recovpop.dat{i_ex}.psc_wfs = {};
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
    if isempty(unique_recov); continue; end
    
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
        insufficient_recovs = Nrecov < MIN_RECOV_COUNT;
        if insufficient_data || insufficient_tfs || insufficient_recovs
            recovpop.dat{i_ex}.ignore{i_ch} = true;
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



% _________ PLOTS __________________
% plot recovery pulse amplitude (normalized) vs. train frequency. Do this
% separately for cell types, brain areas, opsins, etc...
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
    if isempty(recovpop.dat{i_ex});
        disp('here')
        continue;
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
        
    end
end


%
% Plots of the normalized amplitudes
%
%%%%%%%%%%%%
f = figure;
f.Units = 'Normalized';
f.Position = [0.0865, 0, 0.3161, 1];
Ngroups = size(plotgroups,1);
legtext = repmat({{}}, 1, Ntfs);
leghand = repmat({[]}, 1, Ntfs);
for i_group = 1:size(plotgroups, 1)
    plt_clr = hvaPlotColor(plotgroups{i_group, 3});
    
    for i_tf = 1:Ntfs
        
        hs = subplot(Ntfs, 1, i_tf); hold on,
        force_consistent_figs(hs, 'ax');
        
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
            hp = my_errorbar(xx_trains, xbar_trains, sem_trains, 'color', plt_clr, 'linewidth', 2);
            my_errorbar(xx_recov, xbar_recov, sem_recov, 'o', 'markerfacecolor', plt_clr, 'markeredgecolor', plt_clr, 'color', plt_clr, 'linewidth', 2);
            leghand{i_tf}(end+1) = hp(1);
            legtext{i_tf}{end+1} = sprintf('%s, %s, %s, N=%d', plotgroups{i_group,1}, plotgroups{i_group,2}, plotgroups{i_group,3}, N_expts);
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
    title(plt_string, 'fontsize', 12)
    legend(leghand{i_tf}, legtext{i_tf}, 'location', 'best', 'interpreter', 'none')
    legend boxoff
    ylim([0.85, 1.6])
end

