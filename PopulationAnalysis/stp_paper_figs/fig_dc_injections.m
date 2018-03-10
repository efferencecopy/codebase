function fig_1_dc_injections(dat, ex_num, ch_num)



% plot the current step data set to help identify cell types
if ~isempty(dat{ex_num}.dcsteps.Vm_raw{ch_num})
    f = figure;
    f.Units = 'Pixels';
    f.Position = [440    88   483   710];
    
    Vm = dat{ex_num}.dcsteps.Vm_raw{ch_num};
    Icmd = dat{ex_num}.dcsteps.Icmd{ch_num};
    N = size(Vm, 2);
    tt = ([0:N-1] ./ dat{ex_num}.info.sampRate.dcsteps) - dat{ex_num}.info.pretime.dcsteps;
    
    
    ha = subplot(2,1,1); hold on,
    unique_cmd = unique(Icmd);
    trl_list = [];
    for i_swp = 1:numel(unique_cmd)
        trl_idx = find(Icmd==unique_cmd(i_swp), 1, 'first');
        trl_list = cat(1, trl_list, trl_idx);
        plot(tt, Vm(trl_idx,:), 'k', 'linewidth', 2)
    end
    ylim([min(Vm(:)).*1.1, max(Vm(:)).*1.1])
    xlim([min(tt(:)), 0.825])
    force_consistent_figs(ha, 'ax');
    xlabel('Seconds')
    ylabel('mV')
    
    ha = subplot(2,1,2); hold on,
    Icmd_ts = zeros(size(tt));
    l_stimon = (tt>0) & (tt<0.700);
    Icmd_ts(l_stimon) = 1;
    Icmd_ts = repmat(Icmd_ts, numel(Icmd), 1);
    Icmd_ts = bsxfun(@times, Icmd_ts, Icmd');
    plot(tt, Icmd_ts(trl_list,:), 'k', 'linewidth', 2)
    ylim([min(Icmd_ts(:)).*1.1, max(Icmd_ts(:)).*1.1])
    xlim([min(tt(:)), 0.825])
    force_consistent_figs(ha, 'ax');
    xlabel('Seconds')
    ylabel('mV')
end


% series resistance
if isfield(dat{ex_num}.qc, 'Rs') && dat{ex_num}.info.HS_is_valid_Vclamp(ch_num)
    f = figure;
    tmp = squeeze(dat{ex_num}.qc.Rs{ch_num});
    plot(tmp)
    ylim([min([0,min(tmp(:))]), max(tmp(:))*1.5])
    force_consistent_figs(gca, 'ax');
    ylabel('Series Resistance (MOhm)')
    xlabel('trial number')
    title(sprintf('Channel %d', ch_num))
end


% p1amps
if isfield(dat{ex_num}.qc, 'p1amp') && dat{ex_num}.info.HS_is_valid_Vclamp(ch_num)
    f = figure; hold on,
    tmp = squeeze(dat{ex_num}.qc.p1amp{ch_num});
    tmp_norm = squeeze(dat{ex_num}.qc.p1amp_norm{ch_num});
    plot(tmp)
    plot(tmp_norm, 'r')
    ylim([min([0,min(tmp(:))]), max(tmp(:))*1.5])
    force_consistent_figs(gca, 'ax');
    ylabel('P1 Amp')
    xlabel('trial number')
end

