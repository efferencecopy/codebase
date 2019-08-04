function quant_test_pprs(recovpop, groupdata, plotgroups)


anovaTFs = [12, 25, 50];
[all_pnum_labels, all_log_pprs, all_tf_labels, all_raw_pprs] = deal([]);
all_hva_labels = {};
for i_group = 1:size(plotgroups, 1)
    for i_tf = 1:numel(anovaTFs)
        tf_idx = recovpop.TFsAllExpts == anovaTFs(i_tf);
        
        % get the PPRs for the trains only (no recov conds)
        tmp_pprs = log(groupdata.amps{i_group}{tf_idx,1});
        tmp_raw_pprs = (groupdata.amps{i_group}{tf_idx,1});
        
        % make a pulse number index vector
        tmp_pnum_labels = repmat(1:size(tmp_pprs,2), size(tmp_pprs, 1), 1);
        tmp_pnum_labels = tmp_pnum_labels(:);
        
        % make the other grouping variables
        tmp_tf_labels = ones(size(tmp_pnum_labels)) .* anovaTFs(i_tf);
        tmp_hva_labels = repmat({plotgroups{i_group, 3}}, numel(tmp_pnum_labels), 1);
        
        % concatenate the tmp labels into the population vectors
        all_log_pprs = cat(1, all_log_pprs, tmp_pprs(:));
        all_raw_pprs = cat(1, all_raw_pprs, tmp_raw_pprs(:));
        all_pnum_labels = cat(1, all_pnum_labels, tmp_pnum_labels);
        all_tf_labels = cat(1, all_tf_labels, tmp_tf_labels);
        all_hva_labels = cat(1, all_hva_labels, tmp_hva_labels);
    end
end

[~, ~, stats] = anovan(all_log_pprs,...
                        {all_hva_labels, all_tf_labels, all_pnum_labels},...
                        'varnames', {'all_hva_labels', 'all_tf_labels', 'all_pnum_labels'});
if numel(unique(all_pnum_labels)) > 1
    figure
    multcompare(stats, 'Dimension', [3]);
end
if numel(unique(all_tf_labels)) > 1
    figure
    multcompare(stats, 'Dimension', [2]);
end
if size(unique(all_hva_labels), 1) > 1
    figure
    multcompare(stats, 'Dimension', [1]);
end


% p10 LM at 12 Hz is facilitating
fprintf('\n\nLM P10 at 12 Hz is facilitating:\n')
l_p10 = all_pnum_labels == 10;
l_tf12 = all_tf_labels == 12;
l_lm = strcmpi(all_hva_labels, 'lm');
ttest_pprs = all_log_pprs(l_p10 & l_lm & l_tf12);
[h, p] = ttest(ttest_pprs)
mean(ttest_pprs)

% p10 LM at 25 Hz is depressing
fprintf('\n\nLM P10 at 25 Hz is different:\n')
l_p10 = all_pnum_labels == 10;
l_tf25 = all_tf_labels == 25;
l_lm = strcmpi(all_hva_labels, 'lm');
ttest_pprs = all_log_pprs(l_p10 & l_lm & l_tf25);
[h, p] = ttest(ttest_pprs)
mean(ttest_pprs)

% p10 LM at 50 Hz is depressing
fprintf('\n\nLM P10 at 50 Hz is depressing:\n')
l_p10 = all_pnum_labels == 10;
l_tf50 = all_tf_labels == 50;
l_lm = strcmpi(all_hva_labels, 'lm');
ttest_pprs = all_log_pprs(l_p10 & l_lm & l_tf50);
[h, p] = ttest(ttest_pprs)
mean(ttest_pprs)



% p10 AM at 12 Hz is facilitating
fprintf('\n\nAM P10 at 12 Hz is facilitating:\n')
l_p10 = all_pnum_labels == 10;
l_tf12 = all_tf_labels == 12;
l_AM = strcmpi(all_hva_labels, 'AM');
ttest_pprs = all_log_pprs(l_p10 & l_AM & l_tf12);
[h, p] = ttest(ttest_pprs)
mean(ttest_pprs)

% p10 AM at 25 Hz is depressing
fprintf('\n\nAM P10 at 25 Hz is different:\n')
l_p10 = all_pnum_labels == 10;
l_tf25 = all_tf_labels == 25;
l_AM = strcmpi(all_hva_labels, 'AM');
ttest_pprs = all_log_pprs(l_p10 & l_AM & l_tf25);
[h, p] = ttest(ttest_pprs)
mean(ttest_pprs)

% p10 AM at 50 Hz is depressing
fprintf('\n\nAM P10 at 50 Hz is depressing:\n')
l_p10 = all_pnum_labels == 10;
l_tf50 = all_tf_labels == 50;
l_AM = strcmpi(all_hva_labels, 'AM');
ttest_pprs = all_log_pprs(l_p10 & l_AM & l_tf50);
[h, p] = ttest(ttest_pprs)
mean(ttest_pprs)
