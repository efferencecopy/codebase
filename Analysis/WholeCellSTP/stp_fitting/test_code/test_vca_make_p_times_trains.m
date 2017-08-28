function [p_times, freq_t_list, recov_t_list] = test_vca_make_p_times_trains(induction_freqs, n_pulses, recov_sec, n_repeats)

% for book keeping
n_recovs = numel(recov_sec);
n_freqs = numel(induction_freqs);
inds = fullfact([n_freqs, n_recovs]);
n_unique_conds = size(inds,1);

% expand according to n_repeats
inds = repmat(inds, n_repeats, 1);
freq_idx(:,1) = inds(:,1);
recov_idx(:,1) = inds(:,2);


% make ISIs for each unique frequency
isis = 1./induction_freqs(freq_idx);
p_nums = repmat(0:n_pulses-1, numel(isis), 1);
p_times = bsxfun(@times, p_nums, isis(:));

% Add each recovery time to each induction frequency
all_recov_sec = recov_sec(recov_idx);
p_times(:, n_pulses+1) = p_times(:,n_pulses) + all_recov_sec(:);

% make the repeats (shuffling within block)
t_start = 1;
for i_block = 1:n_repeats
    t_stop = t_start+n_unique_conds-1;
    shuffle_inds = randperm(n_unique_conds) + t_start - 1;
    
    % blockwise shuffle the p_times, freq_idx, recov_idx
    p_times(t_start:t_stop, :) = p_times(shuffle_inds, :);
    freq_idx(t_start:t_stop) = freq_idx(shuffle_inds);
    recov_idx(t_start:t_stop) = recov_idx(shuffle_inds);
    
    % update the t_start
    t_start = t_stop + 1;
end

% return the p_times and the trial-by-trial definitions of freq and recov
% (for book keeping later)
 freq_t_list = induction_freqs(freq_idx)';
 recov_t_list = recov_sec(recov_idx)';





