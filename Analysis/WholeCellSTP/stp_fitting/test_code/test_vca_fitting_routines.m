function [fit_results, training_data] = test_vca_fitting_routines(noisy_fake_data_with_ns,...
                                                                  p_times,...
                                                                  smooth_p1_amps,...
                                                                  freq_t_list,...
                                                                  recov_t_list,...
                                                                  params_fitting)

% make a mini dataset that's composed of only pOnTimes, EPSC amps,
% and p1Amps. be sure to aggregate across similar trial types.
training_data = struct('pOnTimes', {{}}, 'rawAmps', {{}}, 'tdict', {[]}, 'unique_conds_rawAmps', {{}});
unique_trial_types = unique([freq_t_list, recov_t_list], 'rows');
for i_cond = 1:size(unique_trial_types, 1)
    
    cond_idx = (freq_t_list == unique_trial_types(i_cond, 1)) & (recov_t_list == unique_trial_types(i_cond, 2));
    
    % grab pOnTimes, raw EPSC amps, and tdict
    pOnTimes = p_times(cond_idx, :);
    pOnTimes = pOnTimes(1,:); % only need one copy
    rawAmps = noisy_fake_data_with_ns(cond_idx, :);
    rawAmps = permute(rawAmps, [2,3,1]); % the rest of the code expects [Npulses x 1 x Nsweeps]
    
    if params_fitting.CONVERT_TO_SMOOTH_P1
        new_p1_amps = smooth_p1_amps(cond_idx);
        rawAmps(1,1,:) = new_p1_amps;
        rawAmps = bsxfun(@rdivide, rawAmps, rawAmps(1,1,:));
    end
    
    % delete the recovery pulse if need be
    if ~params_fitting.FIT_RECOVERY_PULSE
        has_recov_pulse = all(unique_trial_types(i_cond, 2)) > 0;
        if has_recov_pulse
            pOnTimes(end) = [];
            rawAmps = rawAmps(1:end-1,1,:);
        end
    end
    
    % force the model to fit a ficticious data point at +15 seconds
    % at a normalized value of 1
    if params_fitting.FORCE_RECOVERY
        if ~params_fitting.FIT_RECOVERY_PULSE; error('ERROR: need to fit recov for force_recovery'); end
        if ~params_fitting.CONVERT_TO_SMOOTH_P1; error('ERROR: need to normalize for force_recovery'); end
        pOnTimes(end+1:end+3) = [15, 17, 25]; % in seconds
        nsweeps = size(rawAmps,3);
        rawAmps(end+1:end+3,:,:) = repmat([1;1;1], 1, 1, nsweeps);
    end
    
    % allocate the data to either the training or xval datasets
    training_data.pOnTimes = cat(1, training_data.pOnTimes, pOnTimes);
    training_data.rawAmps = cat(1, training_data.rawAmps, rawAmps);
end

% average the data sets (possibly for fitting, but also for plotting)
training_data.unique_conds_xbar = cellfun(@(x) mean(x,3), training_data.rawAmps, 'uniformoutput', false);


% create a normalized version of the xbar arrays (don't normalize the raw data)
training_data.unique_conds_xbar_norm = cellfun(@(x) x./x(1), training_data.unique_conds_xbar, 'uniformoutput', false);


% now assign data to be fit
training_data.ptimes_for_fit = training_data.pOnTimes;
switch lower(params_fitting.DATA_FOR_FIT)
    case 'avg'
        training_data.amps_for_fit = training_data.unique_conds_xbar;
    case 'norm'
        training_data.amps_for_fit = training_data.unique_conds_xbar_norm;
    case 'raw'
        training_data.amps_for_fit = training_data.unique_conds_rawAmps;
        
end



fit_results = {};
% fit the data
if ischar(params_fitting.MODEL)
    [tmp_best_fit_params, fxn_val] = fit_vca_model(training_data.amps_for_fit, training_data.ptimes_for_fit, params_fitting.MODEL);
    [d, tau_d, f, tau_f] = parse_vca_model_coeffs(tmp_best_fit_params, params_fitting.MODEL);
    [d, didx] = sort(d);
    tau_d = tau_d(didx);
    [f, fidx] = sort(f);
    tau_f = tau_f(fidx);
    tmp_best_fit_params = [d, tau_d, f, tau_f];
    fit_results{1}.params = tmp_best_fit_params;
    fit_results{1}.fxn_val = fxn_val;
    fit_results{1}.model = params_fitting.MODEL;
    
elseif isscalar(params_fitting.MODEL)
    all_models = fullfact([params_fitting.MODEL+1, params_fitting.MODEL+1])-1; % -1 to have models with just d or f. +1 to bring the total back up to the desired number
    l_zeros = sum(all_models, 2)==0;
    all_models(l_zeros,:) = []; % delete models with zero d and zero f terms
    
    for i_mod = 1:size(all_models,1)
        d_terms = repmat('d',1,all_models(i_mod,1));
        f_terms = repmat('f',1,all_models(i_mod,2));
        tmp_model = strcat(d_terms, f_terms);
        
        [tmp_best_fit_params, fxn_val] = fit_vca_model(training_data.amps_for_fit, training_data.ptimes_for_fit, tmp_model);
        [d, tau_d, f, tau_f] = parse_vca_model_coeffs(tmp_best_fit_params, params_fitting.MODEL);
        [d, didx] = sort(d);
        tau_d = tau_d(didx);
        [f, fidx] = sort(f);
        tau_f = tau_f(fidx);
        tmp_best_fit_params = [d, tau_d, f, tau_f];
        fit_results{i_mod}.params = tmp_best_fit_params;
        fit_results{i_mod}.fxn_val = fxn_val;
        fit_results{i_mod}.model = tmp_model;
    end
    
else
    error('MODEL was neither a scalar or a string')
end




