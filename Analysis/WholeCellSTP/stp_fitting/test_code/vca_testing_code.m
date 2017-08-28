%% TEST VCA FITTING 

% 1) Generate pTimes (for induction and recovery)

% 1.1) Generate A0 for the model responses (can have non-stationarities,
% and noisy draws from an underlying model)

% 2) Generate model responses to these trains (given A0 and the STP params)

% 3) Fit the vca model to these fake data

% 5) Evaluate the fitted params against the real params.

% 6) Do this for a bunch of different params
%      * add noise to A0
%      * add noise to all pred amps
%      * add non-stationarities (w/wo noise)

%% MAKE A FAKE DATASET
fin

% ----  setup the testing params ------------
model = 'ddff';
induction_freqs = [12, 25, 50];
n_pulses = 10;
recov_sec = [0.3, 1, 5.5];
n_repeats = 3;

A0_init = 250; % pA for first pulse EPSC on sweep1

nsparams.type = 'linear'; % 'none', 'linear'
nsparams.prcnt_rundown = 0;

train_noise_prcnt = 0.15;

n_iters = 25;
%----------------------------------------------

%fake_params = test_vca_make_fake_params(model);
fake_params = [.6, .8, 4, 2, .25, 1, 3, 4];
all_fake_params = [];
all_fitted_params = [];
all_R2 = [];
t_start = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')
for i_iter = 1:n_iters;
    close all
    %
    %  generate the fake data
    %%%%%%
    % make the pulse trains
    [p_times, freq_t_list, recov_t_list] = test_vca_make_p_times_trains(induction_freqs, n_pulses, recov_sec, n_repeats);
    n_sweeps = numel(freq_t_list);
    
    % generate P1 amps, with noise or non-stationarities if necessary
    [p1_amps, noiseless_p1_amps] = test_vca_make_p1_amp(n_sweeps, A0_init, [], nsparams); %ns_params are for non-stationarities
    
    % generate fake data from a synapse that obeys the dynamics specified by
    % the vca model
    %fake_params = test_vca_make_fake_params(model);
    [d, tau_d, f, tau_f] = parse_vca_model_coeffs(fake_params, model);
    fake_data_with_ns = nan(size(p_times));
    for i_swp = 1:n_sweeps
        fake_data_with_ns(i_swp,:) = predict_vca_psc(p_times(i_swp,:), d, tau_d, f, tau_f, p1_amps(i_swp));
    end
    fake_data_no_ns = nan(size(p_times));
    for i_swp = 1:n_sweeps
        fake_data_no_ns(i_swp,:) = predict_vca_psc(p_times(i_swp,:), d, tau_d, f, tau_f, noiseless_p1_amps(i_swp));
    end
    
    noisy_fake_data_with_ns = test_vca_add_noise(fake_data_with_ns, fake_data_no_ns, train_noise_prcnt);
    
    
    
    
    %
    % RECAPITUALATE THD FITTING ROUTINES FROM MY ANALYSIS PIPLINE
    %
    %%%%%%%%%%%%%%
    
    % ---- define the fitting params -----------
    params_fitting.MODEL = model; % could theoretically be different than the one setup above, or could be a scalar
    params_fitting.FIT_RECOVERY_PULSE = true;
    params_fitting.DATA_FOR_FIT = 'avg'; % can be 'avg', 'norm', 'raw'.
    params_fitting.CONVERT_TO_SMOOTH_P1 = true;
    params_fitting.FORCE_RECOVERY = false;
    % -------------------------------------------
    
    fake_p1_amps = noisy_fake_data_with_ns(:,1);
    smooth_p1_amps = test_vca_make_smooth_p1(fake_p1_amps, 6);
    
    fprintf('Fitting iteration %d\n', i_iter)
    [fit_results, training_data] = test_vca_fitting_routines(noisy_fake_data_with_ns,...
        p_times,...
        smooth_p1_amps,...
        freq_t_list,...
        recov_t_list,...
        params_fitting);
    
    
    fit_results = test_vca_plot_comparisons(fit_results, training_data, params_fitting, fake_params);
    
    % aggregate the fake and fitted params to test for bias
    all_fake_params = cat(1, all_fake_params, fake_params);
    all_fitted_params = cat(1, all_fitted_params, fit_results{1}.params);
    all_R2 = cat(1, all_R2, fit_results{1}.R2_train);
end

save('vca_model_testing.mat')
t_finish = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')

%% EXAMINE BIAS IN FITTED PARAMS
close all

n_params = size(all_fake_params, 2);
param_diff_vals = all_fitted_params - all_fake_params;

hf = figure;
hf.Units = 'Normalized';
hf.Position = [0.0583   -0.1044    1.1653    1.0067];
for i_param = 1:n_params
    subplot(4,n_params,i_param);
    hold on,
    hh = histogram(param_diff_vals(:,i_param), 'normalization', 'probability'); %, 'displaystyle', 'stairs');
    hh.EdgeColor = 'w';
    hh.FaceColor = 'k';
    hh.FaceAlpha = 0.8;
    xbar = mean(param_diff_vals(:,i_param));
    s = std(param_diff_vals(:,i_param));
    plot(xbar, 0.05, 'rs', 'markerfacecolor', 'r', 'markersize', 8);
    plot([xbar-s, xbar+s], [0.05, 0.05], '-r', 'linewidth', 3)
    title(sprintf('Param %d', i_param))
    xlabel('fitted - real')
    ylabel('Probability')
    set(gca, 'tickdir', 'out')
    axis tight
    
    subplot(4,n_params,i_param+n_params), hold on,
    plot(all_fitted_params(:,i_param), all_fake_params(:,i_param), 'wo', 'markerfacecolor', 'b', 'markersize', 6);
    maxval = max([all_fitted_params(:,i_param); all_fake_params(:,i_param)]);
    plot([0, maxval], [0 maxval], '--k')
    xlabel('fitted param')
    ylabel('real param')
    set(gca, 'tickdir', 'out')
    axis tight
end

% R2 values for the training data
subplot(4, n_params, (n_params*2)+1)
cdfplot(all_R2)%, 25, 'normalization', 'cdf', 'displaystyle', 'stairs')
title('')
ylabel('CDF')
xlabel('R2 train')

% Line plots of each fit individually. turns out it's not that useful...
subplot(4, n_params, (n_params*2)+2), hold on,
plot([1:n_params], param_diff_vals', 'r', 'linewidth', 0.25)
my_errorbar([1:n_params], mean(param_diff_vals, 1), stderr(param_diff_vals, 1), 'sk', 'linewidth', 1, 'markerfacecolor', 'k');
xlim([0.75, n_params+0.25])
axis tight
ylabel('difference')
xlabel('parameter')


% scatter plots for d's, tau's f's
params_for_scatter = param_diff_vals; % all_fake_params or param_diff_vals or all_fitted_params
subplot(4,n_params, (n_params*2)+3)
plot(params_for_scatter(:,1), params_for_scatter(:,2), 'wo', 'markerfacecolor', 'k', 'linewidth', 0.25, 'markersize', 4)
[rho, p] = corr(params_for_scatter(:,1), params_for_scatter(:,2), 'type', 'spearman');
title(sprintf('r: %.2f p: %.2f', rho, p))
xlabel('d_1')
ylabel('d_2')

subplot(4,n_params, (n_params*2)+4)
plot(params_for_scatter(:,3), params_for_scatter(:,4), 'wo', 'markerfacecolor', 'k', 'linewidth', 0.25, 'markersize', 4)
[rho, p] = corr(params_for_scatter(:,3), params_for_scatter(:,4), 'type', 'spearman');
title(sprintf('r: %.2f p: %.2f', rho, p))
xlabel('Tau d_1')
ylabel('Tau d_2')

subplot(4,n_params, (n_params*2)+5)
plot(params_for_scatter(:,5), params_for_scatter(:,6), 'wo', 'markerfacecolor', 'k', 'linewidth', 0.25, 'markersize', 4)
[rho, p] = corr(params_for_scatter(:,5), params_for_scatter(:,6), 'type', 'spearman');
title(sprintf('r: %.2f p: %.2f', rho, p))
xlabel('f_1')
ylabel('f_2')

subplot(4,n_params, (n_params*2)+6)
plot(params_for_scatter(:,7), params_for_scatter(:,8), 'wo', 'markerfacecolor', 'k', 'linewidth', 0.25, 'markersize', 4)
[rho, p] = corr(params_for_scatter(:,7), params_for_scatter(:,8), 'type', 'spearman');
title(sprintf('r: %.2f p: %.2f', rho, p))
xlabel('Tau f_1')
ylabel('Tau f_2')

% Make some fake trains to test the generalizability of the model fits
% This is useful only when a single set of fake params was used (to test
% bias)

unique_fake_params = unique(all_fake_params, 'rows');

xval_tf = [5, 10, 60];
xval_recov = [0.200];
[p_times_xval, freq_t_list_xval, recov_t_list_xval] = test_vca_make_p_times_trains(xval_tf, n_pulses, xval_recov, 1);
[freq_t_list_xval, sort_idx] = sort(freq_t_list_xval);
p_times_xval = p_times_xval(sort_idx,:);
recov_t_list_xval = recov_t_list_xval(sort_idx);
for i_xval = 1:numel(xval_tf)
    subplot(4, n_params, (n_params*3)+i_xval), hold on,
    
    % make predictions from the fitted params
    pred_data = nan(n_iters, n_pulses+1);
    for i_iter = 1:n_iters
        [d, tau_d, f, tau_f] = parse_vca_model_coeffs(all_fitted_params(i_iter,:), model);
        pred_data(i_iter, :) = predict_vca_psc(p_times_xval(i_xval,:), d, tau_d, f, tau_f, 1);
    end
    
    % make epsc amps from real params
    if size(unique_fake_params, 1)==1
        real_data = nan(1, n_pulses+1);
        [d, tau_d, f, tau_f] = parse_vca_model_coeffs(unique_fake_params, model);
        real_data(1,:) = predict_vca_psc(p_times_xval(i_xval,:), d, tau_d, f, tau_f, 1);
    else
        real_data = nan(n_iters, n_pulses+1);
        for i_iter = 1:n_iters
            [d, tau_d, f, tau_f] = parse_vca_model_coeffs(all_fake_params(i_iter,:), model);
            real_data(i_iter, :) = predict_vca_psc(p_times_xval(i_xval,:), d, tau_d, f, tau_f, 1);
        end
    end
    
    if size(unique_fake_params, 1)==1
        plot(p_times_xval(i_xval,1:n_pulses), pred_data(:, 1:n_pulses)', 'k-')
        plot(p_times_xval(i_xval,n_pulses+1), pred_data(:, n_pulses+1)', 'k.')
        plot(p_times_xval(i_xval,1:n_pulses), real_data(:, 1:n_pulses)', 'r', 'linewidth', 2)
        plot(p_times_xval(i_xval,n_pulses+1), real_data(:, n_pulses+1)', 'sr', 'markerfacecolor', 'r')
        ylabel('norm EPSC')
    else
        dat_for_plot = (pred_data - real_data) ./ real_data;
        std_of_diff = std(dat_for_plot, [], 1);
        mean_of_diff = mean(dat_for_plot, 1);
        herr = shadedErrorBar(p_times_xval(i_xval,1:n_pulses),...
                              mean_of_diff(1:n_pulses),...
                              std_of_diff(1:n_pulses),...
                              {'k-', 'linewidth', 2});
        hrecov = errorbar(p_times_xval(i_xval,n_pulses+1),...
                          mean_of_diff(n_pulses+1),...
                          std_of_diff(n_pulses+1),...
                          'sk', 'linewidth', 2, 'markerfacecolor', 'k');
        plot([0, n_pulses+1], [0 0], ':k')
        ylabel('%diff in norm EPSC')
    end
    isi = min(diff(p_times_xval(i_xval,1:n_pulses)));
    xlim([-isi, p_times_xval(i_xval,n_pulses+1)+isi])
    xlabel('spike time')
    title(sprintf('train: %dHz', xval_tf(i_xval)))
end



% make poisson spike trains and predict them too
% but only if there's a single fake_param set
enforce_simple_cell_envelope = true;
dt = 0.001; % in ms
totalTime = 10; % in seconds
tt = 0:dt:(totalTime-dt);
drift_rate = 2; %Hz
simple_cell_env = sin(2.*pi.*drift_rate.*tt);
simple_cell_env(simple_cell_env<0) = 0;
for i_xval = 1:numel(xval_tf)
    
    % make some isi's that are poisson process
    spkThresh = xval_tf(i_xval) .* dt;
    if enforce_simple_cell_envelope;
        spkThresh = spkThresh .* simple_cell_env;
    end
    randNums = unifrnd(0,1, 1, numel(tt));
    spikeTrain = randNums < spkThresh;
    spikeTimes = tt(spikeTrain);
    
    % make predictions from fitted params
    pred_data = nan(n_iters, numel(spikeTimes));
    for i_iter = 1:n_iters
        [d, tau_d, f, tau_f] = parse_vca_model_coeffs(all_fitted_params(i_iter,:), model);
        pred_data(i_iter, :) = predict_vca_psc(spikeTimes, d, tau_d, f, tau_f, 1);
    end
    
    % make epsc amps from real params
    if size(unique_fake_params, 1)==1
        real_data = nan(1, numel(spikeTimes));
        [d, tau_d, f, tau_f] = parse_vca_model_coeffs(unique_fake_params, model);
        real_data(1,:) = predict_vca_psc(spikeTimes, d, tau_d, f, tau_f, 1);
    else
        real_data = nan(n_iters, numel(spikeTimes));
        for i_iter = 1:n_iters
            [d, tau_d, f, tau_f] = parse_vca_model_coeffs(all_fake_params(i_iter,:), model);
            real_data(i_iter, :) = predict_vca_psc(spikeTimes, d, tau_d, f, tau_f, 1);
        end
    end
    
    % make the subplots
    subplot(4, n_params, (n_params*3)+3+i_xval), hold on,
    if size(unique_fake_params, 1)==1
        plot(spikeTimes, pred_data', 'k-')
        plot(spikeTimes, real_data, 'r', 'linewidth', 2)
        ylabel('norm EPSC')
    else
        dat_for_plot = (pred_data - real_data) ./ real_data;
        std_of_diff = std(dat_for_plot, [], 1);
        mean_of_diff = mean(dat_for_plot, 1);
        herr = shadedErrorBar(spikeTimes, mean_of_diff, std_of_diff, {'k-', 'linewidth', 2});
        plot([1, spikeTimes(end)], [0 0], ':k')
        ylabel('%diff in norm EPSC')
    end
    xlabel('spike time')
    title(sprintf('poiss: %dHz', xval_tf(i_xval)))
    axis tight
    
end



savefig('vca_test_figure')



