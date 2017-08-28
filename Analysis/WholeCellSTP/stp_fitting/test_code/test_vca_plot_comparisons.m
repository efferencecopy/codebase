function fit_results = test_vca_plot_comparisons(fit_results, training_data, params, fake_params)

% predict all the data (for all types of models)
for i_mod = 1:numel(fit_results)
    [d, dTau, f, fTau] = parse_vca_model_coeffs(fit_results{i_mod}.params, fit_results{i_mod}.model);
    training_data.pred_amps{i_mod} = cellfun(@(x,y) predict_vca_psc(x, d, dTau, f, fTau, mean(y(1,1,:))), training_data.ptimes_for_fit, training_data.amps_for_fit, 'uniformoutput', false);
end


clc
hf = figure;
hf.Units = 'Normalized';
hf.Position = [0.4109    0.1093    0.3786    0.7500];



% make a scatter plot of all predicted and actual PSC amps
N_mods = numel(training_data.pred_amps);
plt_clr = lines(N_mods);
for i_mod = 1:N_mods
    hs = [];
    training_pred = cellfun(@(x,y) repmat(x, [1,1,size(y,3)]), training_data.pred_amps{i_mod}, training_data.amps_for_fit, 'uniformoutput', false); % make sure there is one pred for every real amp
    training_pred = cellfun(@(x) x(:), training_pred, 'uniformoutput', false); %  a single col vector for each cond
    training_pred = cat(1, training_pred{:}); % now all conds aggregated into a single col vec
    training_raw = cellfun(@(x) x(:), training_data.amps_for_fit, 'uniformoutput', false);
    training_raw = cat(1, training_raw{:});
    l_ones = training_raw == 1; % valid hits when raw data are normalized
    training_raw(l_ones) = [];
    training_pred(l_ones) = [];
    
    % calculate all the fitting errors
    num_params = numel(fit_results{i_mod}.params);
    fit_results{i_mod}.R2_train = get_r2(training_raw, training_pred);
    fit_results{i_mod}.R2_train_adj = get_r2(training_raw, training_pred, num_params);
    fit_results{i_mod}.AICc_train = get_aic(training_raw, training_pred, num_params);
    fit_results{i_mod}.MSE_train = mean((training_raw - training_pred).^2);
    fit_results{i_mod}.R2_crossvalid = NaN; % this needs to be fixed... get_r2(crossval_raw, crossval_pred);
    

    pltIdx = sub2ind([2, 3], 1, 1);
    hs = subplot(3, 2, pltIdx); hold on,
    plot(training_raw, training_pred, 'o', 'color', plt_clr(i_mod,:))
    
    maxval = max([hs.XLim, hs.YLim]);
    minval = min([hs.XLim, hs.YLim]);
    plot([minval, maxval], [minval, maxval], 'k--')
    xlabel('raw EPSC amp')
    ylabel('pred amp')
    
    pltIdx = sub2ind([2, 3], 1, 2);
    subplot(3,2,pltIdx), hold on,
    resid = training_raw - training_pred;
    histogram(resid, 'FaceColor', plt_clr(i_mod,:))
    plot(mean(resid), 10, 'rv', 'markerfacecolor', 'r')
    xlabel('real-pred')
    if i_mod == N_mods;
        title(sprintf('Best R2 = %.2f', max(cell2mat(structcat(fit_results, 'R2_train')))));
    end
end % i_mod scatter plots

% main figures of raw data and predictions
amps_for_plot = training_data.amps_for_fit;
ptimes_for_plot = training_data.ptimes_for_fit;
pred_for_plot = training_data.pred_amps;

xlims = [inf -inf];
ylims = [inf -inf];
hs = [];
for i_cond = 1:numel(amps_for_plot)
    
    pltIdx = sub2ind([2, numel(amps_for_plot)], 2, i_cond);
    hs(i_cond) = subplot(numel(amps_for_plot), 2, pltIdx); hold on,
    
    % plot the model fits
    xx = ptimes_for_plot{i_cond};
    hp = [];
    model_string = {};
    for i_mod = 1:numel(pred_for_plot)
        hp(i_mod) = plot(xx, pred_for_plot{i_mod}{i_cond}, '-', 'color', plt_clr(i_mod,:), 'linewidth', 2);
        model_string{i_mod} = [fit_results{i_mod}.model, ',', blanks(8-numel(fit_results{i_mod}.model)-1) 'R2: ', num2str(fit_results{i_mod}.R2_train, 3)];
    end
    % plot the experimental data
    switch params.DATA_FOR_FIT
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

if ~isscalar(fake_params)
    pltIdx = sub2ind([2, 3], 1, 3);
    subplot(3,2,pltIdx),, hold on,
    plot(fit_results{1}.params, fake_params, 'o')
    maxval = max([fit_results{1}.params(:); fake_params(:)]);
    plot([0, maxval], [0, maxval], 'k')
    xlim([0, maxval])
    ylim([0, maxval])
end
drawnow

