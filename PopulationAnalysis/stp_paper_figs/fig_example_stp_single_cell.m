function fig_example_stp_single_cell(dat, mouse, site, i_ch)

for idx = 1:numel(dat)
    is_mouse = strcmpi(dat{idx}.info.mouseName, mouse);
    is_site = str2double(dat{idx}.info.siteNum) == site;
    if is_mouse && is_site
        i_ex = idx;
    end
end


% loop through the experiments. Pull out the trains data. Ignore the
% recovery train (if present) and aggregate across recovery conditions.
MIN_TRL_COUNT = 0; % set to 0 to turn off
MIN_TFS_COUNT = 0; % set to 0 to turn off
MIN_RECOV_COUNT = 0; % set to 0 to turn off
recovpop = [];
recovpop.TFsAllExpts = [];
recovpop.recoveryTimesAllExpts = [];

% find the normal trains. Assume the field name is NOT 'ritv'
condnames = fieldnames(dat{i_ex}.expt);
l_trains = ~strncmp(condnames, 'RITv', 4);
if sum(l_trains)==0
    error('Could not find trains')
end % no trains data

condnames_trains = condnames(l_trains);
trainParams = cellfun(@(x) dat{i_ex}.expt.(x).tdict, condnames_trains, 'uniformoutput', false);
trainParams = cat(1, trainParams{:});
recovpop.trainParams{i_ex} = trainParams;
unique_tfs = unique(trainParams(:,3));
unique_recov = unique(trainParams(:,4));
unique_recov(unique_recov == 0) = [];
if isempty(unique_recov)
    error('Could not find recovery pulses')
end

% make sure the pulse amplitude and width were identical across ttypes
assert(numel(unique(trainParams(:,1)))==1, 'ERROR: more than one pulse amplitude')
assert(numel(unique(trainParams(:,2)))==1, 'ERROR: more than one pulse width');

% aggregate data across TF conditons and recovery times (Separately for each
% recording channel)
recovpop.psc_amps = {};
recovpop.psc_wfs = {};
recovpop.ignore={[],[]};
Ntfs = numel(unique_tfs);
Nrecov = numel(unique_recov);

    
    % remove instances where there is insufficient data to constrain
    % the fit
    insufficient_data = size(dat{i_ex}.qc.p1amp{i_ch}, 3) < MIN_TRL_COUNT;
    insufficient_tfs = Ntfs < MIN_TFS_COUNT;
    insufficient_recovs = Nrecov < MIN_RECOV_COUNT;
    if insufficient_data || insufficient_tfs || insufficient_recovs
        error('Not enough data')
    end
    
    
    % initialize the outputs. I need to concatenate similar TFs
    % together, but keep each of their recovery conditions separate
    %
    % recovpop.dat.psc_amps = {{TF}, {psc_amps}, {{recov_time}, {recov_psc_amp}}}
    % recovpop.dat.psc_wfs = {{TF}, {psc_wfs}, {{recov_time}, {recov_psc_wfs}}}
    recovpop.psc_amps{i_ch} = {num2cell(unique_tfs), repmat({[]}, Ntfs,1), repmat({num2cell(unique_recov), repmat({[]}, Nrecov,1)}, Ntfs,1)};
    recovpop.psc_wfs{i_ch} =  {num2cell(unique_tfs), repmat({[]}, Ntfs,1), repmat({num2cell(unique_recov), repmat({[]}, Nrecov,1)}, Ntfs,1)};
    
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
                error('could not find data')
            end
            
            % pull out the EPSC amplitudes (mean across sweeps), norm to
            % first pulse
            tmp_psc = mean(dat{i_ex}.expt.(condnames_trains{cond_idx}).stats.EPSCamp{i_ch}, 3);
            norm_fact = tmp_psc(1);
            tmp_psc = tmp_psc ./ norm_fact;
            
            % concatenate 1:10 into 'trains', and the 11th into 'recov'
            recovpop.psc_amps{i_ch}{2}{i_tf} = cat(1, recovpop.psc_amps{i_ch}{2}{i_tf}, tmp_psc(1:end-1)');
            recovpop.psc_amps{i_ch}{3}{i_tf,2}{i_recov} = cat(1, recovpop.psc_amps{i_ch}{3}{i_tf,2}{i_recov}, tmp_psc(end));
            
            % pull out the EPSC waveforms (mean across sweeps)
            tmp_wfs = mean(dat{i_ex}.expt.(condnames_trains{cond_idx}).raw.snips{i_ch}, 3);% Npulses x Ntime
            norm_fact = min(tmp_wfs(1,:), [], 2);
            tmp_wfs = tmp_wfs ./ -norm_fact;
            
            % store the waveforms in the population structures
            recovpop.psc_wfs{i_ch}{2}{i_tf} = cat(3, recovpop.psc_wfs{i_ch}{2}{i_tf}, tmp_wfs(1:end-1, :));
            recovpop.psc_wfs{i_ch}{3}{i_tf,2}{i_recov} = cat(3, recovpop.psc_wfs{i_ch}{3}{i_tf,2}{i_recov}, tmp_wfs(end,:));
            
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


% set up the output arrays
%
%  {{trains_(tf=X)}{recov=A, tf=X}{recov=B, tf=X}...} one cell for each recovery pulse, one row for each TF
%
allTFs = recovpop.TFsAllExpts;
Ntfs = numel(allTFs);
allRecoveryTimes = recovpop.recoveryTimesAllExpts';
Nrecovs = numel(allRecoveryTimes);
template = repmat({[]}, Ntfs, Nrecovs+1);
groupdata.wfs = template;
groupdata.amps = template;


% pull out the waveform data
%
% Need to order things by ALL TFS and RECOVs, and not w/r/t an
% individual data file's lattice
%
wf_top_level = recovpop.psc_wfs{i_ch};
amps_top_level = recovpop.psc_amps{i_ch};
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
    tmp_group = groupdata.wfs{row_idx, 1};
    tmp_group = cat(1, tmp_group, tf_wf); % Ntraces x Ntime
    groupdata.wfs{row_idx, 1} = tmp_group;
    
    % push the amps data into the correct slot in the groupdata
    tmp_group = groupdata.amps{row_idx, 1};
    tmp_group = cat(1, tmp_group, tf_amps); % Nneurons x Npulses
    groupdata.amps{row_idx, 1} = tmp_group;
    
    % loop over the recovery conditions and deal with those
    ex_recovs = cell2mat(wf_top_level{3}{i_tf, 1});
    wfs_recovs = cellfun(@(x) mean(x, 3), wf_top_level{3}{i_tf,2}, 'uniformoutput', false); % average across replicates
    amps_recovs = cellfun(@(x) mean(x, 3), amps_top_level{3}{i_tf, 2}, 'uniformoutput', false); % average across replicates. Each recovTime gets a scalar
    for i_recov = 1:numel(ex_recovs)
        col_idx = [false, allRecoveryTimes == ex_recovs(i_recov)]; % leading False to offset the "trains" position
        tmp_group = groupdata.wfs{row_idx, col_idx};
        tmp_group = cat(1, tmp_group, wfs_recovs{i_recov});
        groupdata.wfs{row_idx, col_idx} = tmp_group;
        
        % deal with amps
        tmp_group = groupdata.amps{row_idx, col_idx};
        tmp_group = cat(1, tmp_group, amps_recovs{i_recov});
        groupdata.amps{row_idx, col_idx} = tmp_group;
    end
end






%
% Plots of the waveforms
%
%%%%%%%%%%%%
f = figure;
f.Units = 'Normalized';
f.Position = [0.2889    0.1011    0.3076    0.7722];
f.Name = sprintf('%s, site_%d, ch_%d', mouse, site, i_ch);
f.Color = 'w';
min_val = 0;
for i_tf = 1:Ntfs
    
    subplot(Ntfs, 1, i_tf), hold on,
    force_consistent_figs(gca, 'ax');
    
    % concatenate the time series, padding with nans
    N_expts = cellfun(@(x) size(x,1), groupdata.wfs(i_tf,:));
    N_expts = max(N_expts);
    all_wfs = cellfun(@(x) cat(1, x, nan(N_expts-size(x,1), size(x,2))), groupdata.wfs(i_tf,:), 'uniformoutput', false);
    all_wfs = cellfun(@(x) cat(2, x, nan(N_expts, 200)), all_wfs, 'uniformoutput', false);
    all_wfs = cat(2, all_wfs{:});
    
    % determine the mean across experiments
    xbar = nanmean(all_wfs, 1);
    sem = stderr(all_wfs, 1);
    
    % mean +/- sem
    hp = shadedErrorBar(1:numel(xbar), xbar, sem, {'color', 'k'});
    
    % for the ylim:
    min_val = min([0, min(all_wfs(:))]);
end
% indicate the TF, and recovery times, add a legend
opsin = dat{i_ex}.info.opsin;
cell_type = dat{i_ex}.info.cellType{i_ch};
hva = dat{i_ex}.info.brainArea;
for i_tf = 1:Ntfs
    subplot(Ntfs, 1, i_tf)
    plt_train_tf = allTFs(i_tf);
    l_recov_exist = cellfun(@(x) ~isempty(x), groupdata.wfs(i_tf,:));
    l_recov_exist(1) = []; % ignore the first array b/c it's not recov pulse
    plt_recov_ms = allRecoveryTimes(l_recov_exist);
    plt_string = sprintf('Train = %d Hz, Recovery (in ms): %s, opsin: %s, cell: %s, HVA: %s', plt_train_tf, num2str(plt_recov_ms), opsin, cell_type, hva);
    title(plt_string, 'fontsize', 12)
    %ylim([min_val*1.05, 10])
    axis tight
end







