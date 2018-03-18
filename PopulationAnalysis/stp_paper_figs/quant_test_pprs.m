function quant_test_pprs(recovpop, groupdata, plotgroups)


anovaTFs = [12, 25, 50];
[all_pnum_labels, all_log_pprs, all_tf_labels, all_hva_labels] = deal([]);
for i_group = 1:size(plotgroups, 1)
    for i_tf = 1:numel(anovaTFs)
        tf_idx = recovpop.TFsAllExpts == anovaTFs(i_tf);
        
        % get the PPRs for the trains only (no recov conds)
        tmp_pprs = log(groupdata.amps{i_group}{tf_idx,1});
        
        % make a pulse number index vector
        tmp_pnum_labels = repmat(1:size(tmp_pprs,2), size(tmp_pprs, 1), 1);
        tmp_pnum_labels = tmp_pnum_labels(:);
        
        % make the other grouping variables
        tmp_tf_labels = ones(size(tmp_pnum_labels)) .* anovaTFs(i_tf);
        tmp_hva_labels = repmat(plotgroups{i_group, 3}, numel(tmp_pnum_labels), 1);
        
        % concatenate the tmp labels into the population vectors
        all_log_pprs = cat(1, all_log_pprs, tmp_pprs(:));
        all_pnum_labels = cat(1, all_pnum_labels, tmp_pnum_labels);
        all_tf_labels = cat(1, all_tf_labels, tmp_tf_labels);
        all_hva_labels = cat(1, all_hva_labels, tmp_hva_labels);
    end
end

[p, tbl, stats] = anovan(all_log_pprs,...
                        {all_hva_labels, all_tf_labels, all_pnum_labels},...
                        'model', 'interaction',...
                        'varnames', {'all_hva_labels', 'all_tf_labels', 'all_pnum_labels'});
figure
multcompare(stats, 'Dimension', [3]);
figure
multcompare(stats, 'Dimension', [2]);
figure
multcompare(stats, 'Dimension', [1]);