function dat = wcstp_compile_data(exinfo, hidx, params)
%
% INPUTS:
%
% exinfo:  a cell array of experiment attributes, one attribute per column.
%          Contains things like file names etc.
% 
% hidx:    a structure, that contains the lookup info for the fields in
%          'exinfo'. 
%
% OUTPUTS:
%
% dat.qc.Rs         -> series resistance for HS1,2 [1xx1xNsweeps]
% dat.qc.p1amp      -> amplitude for first pulse of each sweep for HS1,2 [1xx1xNsweeps]
% dat.qc.verr       -> measured holding potential for each sweep, HS1,2 [1xx1xNsweeps]
% dat.qc.instNoise  -> the noise in the holding current which indicates changes in recording quality
%
%
% dat.dcsteps.Vm_raw       -> raw voltage traces during stimulus on period for HS1,2 {N_pA_Cmd}[Nsweeps x Ntime]
% dat.dcsteps.Icmd         -> magnitude of current injection on a sweep by sweep basis. [N_pA_Cmd x 1]
% dat.dcsteps.Vrest        -> the resting potential across sweeps.
% dat.dcsteps.IVpeak.raw   -> the pA_cmd and peak Vm change [pA_Cmd, Vm_change]
% dat.dcsteps.IVasym.raw   -> the pA_cmd and the steady-state Vm change [pA_Cmd, Vm_change]
% dat.dcsteps.IVpeak.Rin   -> the input resistance, estimated from small current injections
% dat.dcsteps.pA_holding   -> Holding current (to keep a crappy cell ~ -65 mV)
% dat.dcsteps.Ih_sag       -> Difference b/w peak and steady state Vm for fairly large DC hyperpolarizations
% 
%
% dat.expt.[stimType].raw.snips           -> trace surrounding the EPSC for HS1,2 [Npulses x Ntime x Nsweeps]
% dat.expt.[stimType].stats.EPSCamp       -> EPSC amplitude for HS1,2 [Npulses x 1 x Nsweeps], when raw current is inward
% dat.expt.[stimType].stats.IPSCamp       -> IPSC amplitude for HS1,2 [Npulses x 1 x Nsweeps], when raw current is outward
% dat.expt.[stimType].stats.latency       -> time to peak of PSC in seconds [Npulses x 1 x Nsweeps]
% dat.expt.[stimType].realTrialNum        -> the actual trial number for each of the sweeps, this could change between HS b/c trials can get cut from one but not the other.
% dat.expt.[stimType].tdict               -> the trialtype dictionary entry for this condition
% 
% dat.info.opsin
% dat.info.fid.dcsteps
% dat.info.fid.vclamp
% dat.info.mouseName
% dat.info.siteNum
% dat.info.cellType
% dat.info.brainArea
% dat.info.HS_is_valid_Vclamp
% dat.info.HS_is_valid_Iclamp
% dat.info.pretime.vclamp
% dat.info.pretime.dcsteps
% dat.info.posttime.vclamp
% dat.info.posttime.dcsteps
% dat.info.sampRate.vclamp
% dat.info.sampRate.dcsteps
% dat.info.sweepLength.vclamp   -> Total time of sweep. Useful for plotting 



% fill out the info field and get things ready for the analysis
dat.info = defineInfo(exinfo, hidx, params);

try
    % the function is verbose. Tell the user what's going on
    fprintf('Analyzing data from mouse %s, site %s\n', dat.info.mouseName, dat.info.siteNum)
    
    % unpack the DC current injection dataset
    dat = unpack_dc_injections(dat);
    
    % unpack the Vclamp trains dataset(s)
    dat = unpack_vclamp_trains(dat, exinfo, hidx, params.force_epsc);
    
    % unpack the Iclamp trains dataset
    dat = unpack_iclamp_trains(dat, exinfo, hidx);
    
    % get an estimate of the P1 amp smoothed across time
    N = 6;
    dat = smoothP1amp(dat, N);
    
catch
    fprintf('crashed on this cell:\n')
    dat.info
    exinfo
    error('analysis has stopped')
end

end


function info = defineInfo(exinfo, hidx, params)

    info.fid.dcsteps_hs1 = exinfo{hidx.ABFDCstepsHS1};
    info.fid.dcsteps_hs2 = exinfo{hidx.ABFDCstepsHS2};
    info.fid.vclamp = exinfo{hidx.ABFOptostimVclamp};
    info.fid.iclamp = exinfo{hidx.ABFOptostimIclamp};

    info.opsin = exinfo{hidx.opsin};
    info.mouseName = exinfo{hidx.MouseName};
    info.siteNum = exinfo{hidx.Site};
    info.cellType = exinfo([hidx.Celltype1, hidx.Celltype2]); % one for each HS
    info.brainArea = exinfo{hidx.brainarea};
    info.HS_is_valid_Vclamp = [str2double(exinfo{hidx.VC_HS1_valid}), str2double(exinfo{hidx.VC_HS2_valid})];
    info.HS_is_valid_Iclamp = [str2double(exinfo{hidx.IC_HS1_valid}), str2double(exinfo{hidx.IC_HS2_valid})];
    info.laser_power_supply_ND_filt = exinfo{hidx.laserpowsupply};
    
    info.pretime.vclamp = params.pretime.vclamp;
    info.posttime.vclamp = params.posttime.vclamp;
    info.pretime.dcsteps = params.pretime.dcsteps;
    info.posttime.dcsteps = params.posttime.dcsteps;
    info.pretime.iclamp = params.pretime.iclamp;
    info.posttime.iclamp = params.posttime.iclamp;
    
    % stuff that is only available when there was a Vclamp experiment
    if any(info.HS_is_valid_Vclamp)
        if ~isnan(exinfo{hidx.xystim}) % sometime the field is blank (i.e., I forgot to measure this...)
            info.stim_xy_pos = str2num(exinfo{hidx.xystim});
        end
    end
    
    % determine the depth of each recording channel w/r/t the pia
    info.cellDepth_um = [nan, nan];
    info.HS_xy_pos = {[nan,nan],[nan,nan]};
    for i_hs = 1:2
        if ~(info.HS_is_valid_Vclamp(i_hs) || info.HS_is_valid_Iclamp(i_hs)); continue; end
        
        hsname = sprintf('xyHS%d', i_hs);
        
        if any(isnan(exinfo{hidx.(hsname)})) || isempty(exinfo{hidx.(hsname)}); continue; end
        xy_hs = str2num(exinfo{hidx.(hsname)});
        
        if any(isnan(exinfo{hidx.xypia})) || isempty(exinfo{hidx.xypia}); continue; end
        xy_pia = str2num(exinfo{hidx.xypia});
        
        info.cellDepth_um(i_hs) = norm(xy_hs - xy_pia);
        info.HS_xy_pos{i_hs} = xy_hs;
        info.pia_xy_pos = xy_pia;
    end
    
end

function dat = unpack_dc_injections(dat)

    % dat.dcsteps.Vm_raw: raw voltage traces during stimulus on period for HS1,2 [Nsweeps x Ntime]
    % dat.dcsteps.Icmd: magnitude of current injection on a sweep by sweep basis. [Nsweeps x 1]
    
    fprintf('  Unpacking DC current injections\n')
    
    % loop over channels, which may have their own ABF file
    for i_ch = 1:2
        
        % initialize the ouput in case there isn't a vaild datafile, or
        % there were other issues that prevented storing data
        dat.dcsteps.Vm_raw{i_ch} = [];
        dat.dcsteps.Icmd{i_ch} = [];
        
        % make sure there is a valid file.
        hs_fid = sprintf('dcsteps_hs%d', i_ch);
        if strncmpi(dat.info.fid.(hs_fid)(end-8:end), [filesep, 'none.abf'], 9)
            continue
        end
        
        % open an ABF file for HS1, and if there's a different file for HS2
        if i_ch == 1 || (i_ch == 2 && ~strcmpi(dat.info.fid.dcsteps_hs1, dat.info.fid.dcsteps_hs2))
            ax = abfobj(dat.info.fid.(hs_fid));
        end
        
        
        % fill out the info struct
        dat.info.sampRate.dcsteps = ax.head.sampRate;
        dat.info.fileStartTime_24hrs = ax.head.uFileStartTimeMS ./ (60*60*1000);
        
        % figure out which channels were turned on, and set to Iclamp
        Iclamp_idx = strcmpi(ax.head.recChUnits, 'mV');
        secCh = cellfun(@any, regexpi(ax.head.recChNames, '_sec'));
        Iclamp_idx = Iclamp_idx & ~secCh;
        
        % make sure data are present for this recording channel. If not,
        % continue to the next channel
        Iclamp_name_is_HS_name = cellfun(@(x) ~isempty(regexpi(x,sprintf('HS%d_', i_ch))), ax.head.recChNames);
        HSpresent = any(Iclamp_name_is_HS_name & Iclamp_idx);
        assert(HSpresent, 'ERROR: DC ABF file supplied, but data from HS could not be located')
        HSname = sprintf('HS%d_Vm', i_ch);
        
        % extract raw voltage traces and spike times.
        preTime_samps = ceil(dat.info.pretime.dcsteps .* dat.info.sampRate.dcsteps);
        postTime_samps = ceil(dat.info.posttime.dcsteps .* dat.info.sampRate.dcsteps);

        
        % determine the magnitude of each sweep's current injection and
        % make a trial type dictionary
        Iclamp_cmd_name = sprintf('HS%d_secIm', i_ch);
        stimWF = ax.dat(:,ax.idx.(Iclamp_cmd_name),:);
        stimWF = permute(stimWF, [1,3,2]);
        stimWF = bsxfun(@minus, stimWF, mean(stimWF(1:preTime_samps, :),1));
        
        diff_stim_wf = diff(stimWF,1,1);
        threshold = 9;
        l_aboveThreshold = sum((diff_stim_wf > threshold) | (diff_stim_wf < -threshold), 1) > 0;
        diff_stim_wf = diff_stim_wf(:,l_aboveThreshold);
        diff_stim_wf = bsxfun(@rdivide, diff_stim_wf, max(diff_stim_wf,[],1));
        diff_stim_wf = mean(abs(diff_stim_wf),2);
        aboveThreshold = diff_stim_wf > 0.8;
        crossing_on = [false(1, size(aboveThreshold,2)); diff(aboveThreshold, 1, 1) == 1];
        crossing_on = find(crossing_on == 1);
        assert(numel(crossing_on) == 2, 'ERROR: may not have correctly identified the start stop of Iclamp pulse')
        
        idx_on = crossing_on(1);
        idx_off = crossing_on(2);
        
        stimWF = bsxfun(@minus, stimWF, mean(stimWF([idx_on-preTime_samps : idx_on-1], :),1)); % re-baseline now that you know when the stimulus came on
        
        quarter_time = round((idx_off - idx_on) .* 0.20);
        idx_Icmd_middle_on = idx_on + quarter_time;
        idx_Icmd_middle_off = idx_off - quarter_time;
        trl_Iclamp_cmd_pA = mean(stimWF([idx_Icmd_middle_on:idx_Icmd_middle_off], :), 1);
        
        % deal with small amounts of noise that make idential Icmd pulses
        % seem different
        for i_cmd = 1:numel(trl_Iclamp_cmd_pA)
            diffvals = trl_Iclamp_cmd_pA - trl_Iclamp_cmd_pA(i_cmd);
            l_same = abs(diffvals) < 2.5;
            trl_Iclamp_cmd_pA(l_same) = round(mean(trl_Iclamp_cmd_pA(l_same)));
        end
        unique_cmd_pA = unique(trl_Iclamp_cmd_pA);
        assert(~any(diff(unique_cmd_pA)<5), 'ERROR: possible duplicate Iclamp command pA due to noise?')
        
        % store the Iclamp command injection magnitudes
        dat.dcsteps.Icmd{i_ch} = trl_Iclamp_cmd_pA;
        
        % now make a list of indicies that correspond to stim on +/-
        % pre/post time
        stim_idx = (idx_on - preTime_samps) : (idx_off + postTime_samps);
        
        % pull out the Vm data, store it in the output structure along with
        % the resting membrane potential
        raw_Vm = permute(ax.dat(stim_idx, ax.idx.(HSname), :), [3,1,2]);
        baseline = mean(raw_Vm(:,1:preTime_samps), 2);
        delta_baseline = baseline - mean(baseline);
        assert(~any(delta_baseline>7.5), sprintf('ERROR: file has Vm that changed from sweep to sweep.\n fid: %s', dat.info.fid.(hs_fid)))
        raw_Vm = bsxfun(@minus, raw_Vm, delta_baseline); % adjust for small differences in resting Vm from sweep to sweep
        
        dat.dcsteps.Vm_raw{i_ch} = raw_Vm; % [Nsweeps x Ntime]
        dat.dcsteps.Vrest{i_ch} = mean(baseline);
        
        %
        % estimate the input resistance as the slope of the I-V curve
        %%%%%%%%%%%%%%%%
        raw_Vm = permute(ax.dat([idx_on:(idx_off-1)], ax.idx.(HSname), :), [3,1,2]); % subtracting off 1 samps to avoid issues with obo errors with turning off currents
        raw_Vm = bsxfun(@minus, raw_Vm, delta_baseline);
        tt = [0:size(raw_Vm,2)-1] ./ dat.info.sampRate.dcsteps;
        l_window_asym = (tt >= 0.450)&(tt < 0.500);
        assert(any(tt >= 0.499), 'ERROR: current injection was too short')
        steady_state_mV = mean(raw_Vm(:, l_window_asym), 2);
        
        % store the I-V data, but only for sub threshold responses
        spike_threshold = -10;
        l_no_spikes = max(raw_Vm,[],2) < spike_threshold;
        subthresh_cmd_pA = trl_Iclamp_cmd_pA(l_no_spikes);
        subthresh_resp_mV =  steady_state_mV(l_no_spikes);
        unique_subthresh_cmd = unique(subthresh_cmd_pA);
        asym_vals = [];
        for i_cmd = 1:numel(unique_subthresh_cmd)
            t_list = subthresh_cmd_pA == unique_subthresh_cmd(i_cmd);
            asym_vals(i_cmd,:) = [unique_subthresh_cmd(i_cmd),  mean(subthresh_resp_mV(t_list)), stderr(subthresh_resp_mV(t_list))];
        end
        dat.dcsteps.IVasym.raw{i_ch} = asym_vals;
        
        % now find the maximal slope of the I-V curve in a small neighborhood around Vrest
        pA_cmd = asym_vals(:,1);
        mV_asym = asym_vals(:,2);
        l_valid_pA = (pA_cmd ~= 0) & (pA_cmd <= 40) & (pA_cmd >= -65);
        if any(l_valid_pA)
            MOhm = (mV_asym(l_valid_pA) - dat.dcsteps.Vrest{i_ch}) ./ (pA_cmd(l_valid_pA)/1000);
            MOhm = max(abs(MOhm));
            assert(~isinf(MOhm), 'ERROR: Rinput is Inf')
        else
            MOhm = NaN;
        end
        dat.dcsteps.IVasym.Rin{i_ch} = MOhm;
        
        % now only consider +/- going current steps
        MOhm = (mV_asym - dat.dcsteps.Vrest{i_ch}) ./ (pA_cmd/1000);
        MOhm = abs(MOhm);
        l_neg = l_valid_pA & pA_cmd<0;
        l_pos = l_valid_pA & pA_cmd>0;
        dat.dcsteps.Rin_posPA{i_ch} = max(MOhm(l_pos));
        dat.dcsteps.Rin_negPA{i_ch} = max(MOhm(l_neg));
        
        
        % compute the Vm response right after the DC injection starts (as a
        % proxy for the 'peak' response). Average a window 6 to 16 ms
        % after the injection starts.
        tt = [0:size(raw_Vm,2)-1] ./ dat.info.sampRate.dcsteps;
        l_window_peak = (tt >= 0.012)&(tt < 0.025);
        peak_mV = mean(raw_Vm(:, l_window_peak),2);
        subthresh_cmd_pA = trl_Iclamp_cmd_pA(l_no_spikes);
        subthresh_resp_mV =  peak_mV(l_no_spikes);
        unique_subthresh_cmd = unique(subthresh_cmd_pA);
        peak_vals = [];
        for i_cmd = 1:numel(unique_subthresh_cmd)
            t_list = subthresh_cmd_pA == unique_subthresh_cmd(i_cmd);
            peak_vals(i_cmd,:) = [unique_subthresh_cmd(i_cmd),  mean(subthresh_resp_mV(t_list)), stderr(subthresh_resp_mV(t_list))];
        end
        dat.dcsteps.IVpeak.raw{i_ch} = peak_vals;
        
        % Re-estimate input resistance according to the "peak" values
        pA_cmd = peak_vals(:,1);
        mV_peak = peak_vals(:,2);
        l_valid_pA = (pA_cmd ~= 0) & (pA_cmd <= 40) & (pA_cmd >= -65);
        if any(l_valid_pA)
            MOhm = (mV_peak(l_valid_pA) - dat.dcsteps.Vrest{i_ch}) ./ (pA_cmd(l_valid_pA)/1000);
            MOhm = max(abs(MOhm));
            assert(~isinf(MOhm), 'ERROR: Rinput is Inf')
        else
            MOhm = NaN;
        end
        dat.dcsteps.IVpeak.Rin{i_ch} = MOhm;
        
        % compute the membrane time constant for sub-threhold current
        % injections. Do this two ways: by estimation of exponential tau,
        % and by eye-balling 63% decay
        %
        % return [pA_cmd, Vm_at_cmd, tau_fitted, tau_quick]
        
        % iterate over subthreshold current injections:
        raw_Vm_offset = permute(ax.dat(:, ax.idx.(HSname), :), [3,1,2]);
        baseline = mean(raw_Vm_offset(:, idx_on-preTime_samps:idx_on), 2);
        delta_baseline = baseline - mean(baseline);
        assert(~any(delta_baseline>7.5), sprintf('ERROR: file has Vm that changed from sweep to sweep.\n fid: %s', dat.info.fid.(hs_fid)))
        raw_Vm_offset = bsxfun(@minus, raw_Vm_offset, delta_baseline); % adjust for small differences in resting Vm from sweep to sweep
        raw_Vm_offset = raw_Vm_offset(:, idx_off+10:idx_off+postTime_samps); % just getting the decay
        l_valid_pA = (trl_Iclamp_cmd_pA ~= 0) & (trl_Iclamp_cmd_pA <= 40) & (trl_Iclamp_cmd_pA >= -65);
        unique_valid_cmd = unique(trl_Iclamp_cmd_pA(l_valid_pA));
        tau_data = nan(numel(unique_valid_cmd), 3);
        for i_pa = 1:numel(unique_valid_cmd)
            pa_idx = trl_Iclamp_cmd_pA == unique_valid_cmd(i_pa);
            avg_off_Vm = mean(raw_Vm_offset(pa_idx, :), 1);
            mv_at_offset = avg_off_Vm(1);
            tt_off_tau = [0:numel(avg_off_Vm)-1] ./ dat.info.sampRate.dcsteps;
            if unique_valid_cmd(i_pa) < 0;
                avg_off_Vm = -avg_off_Vm;
            end
            [Ao, tau, C, tau_quick_est] = fit_exp1(avg_off_Vm, tt_off_tau);
            tau_data(i_pa, 1) = unique_valid_cmd(i_pa);
            tau_data(i_pa, 2) = mv_at_offset;
            tau_data(i_pa, 3) = tau;
            tau_data(i_pa, 4) = tau_quick_est;
        end
        dat.dcsteps.tau{i_ch} = tau_data;
        
        
        % estimate the amount of Ih at the onset of the negative current
        % injection. Find the "true peak" of the Ih sag instead of using
        % the IVpeak.raw values (which underestimate the Ih sag). Then find
        % the saddle point (the most positive point) between there and the
        % end of the injection. This counteracts any hyperpolarizing
        % current that exists to work against the Ih.
        l_window_peak = (tt >= 0.002)&(tt < 0.050);
        n_samps_true_peak = round(ax.head.sampRate .*125e-6);
        n_samps_smooth_kernel = round(ax.head.sampRate .* 20e-3);
        l_sweeps_for_Ih_analysis = trl_Iclamp_cmd_pA < -400;
        unique_pA_levels = unique(trl_Iclamp_cmd_pA(l_sweeps_for_Ih_analysis));
        
        % initialize the outputs
        dat.dcsteps.Ih_sag.pA{i_ch} = [];
        dat.dcsteps.Ih_sag.peak_Vm{i_ch} = [];
        dat.dcsteps.Ih_sag.sag{i_ch} = [];
        dat.dcsteps.Ih_sag.Vm_asym{i_ch} = [];
        
        for i_pa = 1:numel(unique_pA_levels)
            
            % average sweeps at the same current injection strength
            l_sweeps = trl_Iclamp_cmd_pA == unique_pA_levels(i_pa);
            tmp_raw = mean(raw_Vm(l_sweeps, :), 1);
            
            % average samps for true peak.
            true_peak = min(tmp_raw(l_window_peak));
            peak_idx = find(tmp_raw == true_peak, 1, 'first');
            true_peak = mean( tmp_raw((peak_idx - n_samps_true_peak) : (peak_idx + n_samps_true_peak)) );
            
            % locate the saddle point b/w the true peak and the asympote point (500
            % ms)
            l_window_saddle = (tt > tt(round(n_samps_smooth_kernel))) & (tt <= 0.500); % smoothing affects the waveform, so ignore the first samples
            smooth_raw = tmp_raw;
            smooth_raw(1:peak_idx) = true_peak;
            smooth_raw = filtfilt(ones(1, n_samps_smooth_kernel)./n_samps_smooth_kernel, 1, smooth_raw);
            smooth_raw = smooth_raw(l_window_saddle);
            [~, saddle_idx] = max(smooth_raw(n_samps_true_peak : end-n_samps_true_peak-1)); % watch for edge effects
            saddle_idx = saddle_idx + n_samps_true_peak;
            saddle_value = mean(smooth_raw((saddle_idx - n_samps_true_peak) : (saddle_idx + n_samps_true_peak)));
            
            
            % store the data
            dat.dcsteps.Ih_sag.pA{i_ch}(i_pa) = unique_pA_levels(i_pa);
            dat.dcsteps.Ih_sag.peak_Vm{i_ch}(i_pa) = true_peak;
            dat.dcsteps.Ih_sag.sag{i_ch}(i_pa) = true_peak - saddle_value;
            dat.dcsteps.Ih_sag.Vm_asym{i_ch}(i_pa) = mean(tmp_raw(l_window_asym));
            
        end
        
        % flag instances with holding current
        idx_secCh = ax.idx.(sprintf('HS%d_secIm', i_ch));
        sec_raw_pA = permute(ax.dat(:, idx_secCh, :), [3, 1, 2]);
        sec_raw_pA = sec_raw_pA(:,[idx_on-preTime_samps : idx_on-1]);
        sec_raw_pA = mean(sec_raw_pA, 2);
        dat.dcsteps.pA_holding{i_ch} = [min(sec_raw_pA), max(sec_raw_pA)];
        
        
        %
        % estimate spike f-I relationship. Need to get spike rate, and Icmd
        %%%%%%%%%%
        raw_Vm = permute(ax.dat([idx_on:idx_off], ax.idx.(HSname), :), [3,1,2]); % subtracting off 1 samps to avoid issues with obo errors with turning off currents
        raw_Vm = bsxfun(@minus, raw_Vm, delta_baseline);
        above_threshold_cross_up = [false(size(raw_Vm,1), 1), diff(raw_Vm > spike_threshold, 1, 2)==1];
        n_spikes = sum(above_threshold_cross_up,2);
        total_time = size(raw_Vm,2) ./ dat.info.sampRate.dcsteps;
        spike_rate = n_spikes ./ total_time;
        fi_dat = [];
        for pA = unique(trl_Iclamp_cmd_pA)
            l_trl = trl_Iclamp_cmd_pA == pA;
            fi_dat = cat(1, fi_dat, [pA, mean(spike_rate(l_trl))]);
        end
        dat.dcsteps.fi_curve{i_ch} = fi_dat;
        
        
        %
        % Estimate the spike threshold
        %%%%%%%%%%%%
        dvdt = [zeros(size(raw_Vm,1),1), diff(raw_Vm, 1, 2)];
        max_dvdt = max(dvdt, [], 2);
        dvdt_thresh = 0.04 .* max_dvdt;
        ms1_in_samps = round(0.001 .* dat.info.sampRate.dcsteps);
        thresholds_first_last = nan(size(raw_Vm,1), 2);
        resets_first_last = nan(size(raw_Vm,1), 2);
%         figure
%         set(gcf, 'position', [12         383        1903         420])
%         fig_name = sprintf('  %s, site: %s, ch: %d', dat.info.mouseName, dat.info.siteNum, i_ch);
%         set(gcf, 'name', fig_name)
%         fig_count = 0;
        for i_trl = 1:size(raw_Vm,1)
%             if ~l_no_spikes(i_trl)
%                 fig_count = fig_count + 1;
%                 subplot(1,sum(~l_no_spikes), fig_count), hold on,
%             end
            
            spike_idx = [false, diff(raw_Vm(i_trl,:)>spike_threshold)==1];
            tmp_spike_thresholds = [];
            tmp_spike_reset = [];
            if any(spike_idx)
                trl_dvdt = dvdt(i_trl,:);
                spike_idx = find(spike_idx);
                for i_spk = 1:numel(spike_idx)
                    % get spike thresholds
                    first_samp = max([1, spike_idx(i_spk)-ms1_in_samps]);
                    spk_dvdt = trl_dvdt(first_samp : spike_idx(i_spk));
                    spk_raw_Vm = raw_Vm(i_trl, first_samp : spike_idx(i_spk));
                    thresh_idx = find(spk_dvdt>dvdt_thresh(i_trl), 1, 'first');
                    if isempty(thresh_idx); % slope criterion not met
                        thresh_val = nan;
                    else
                        thresh_val = spk_raw_Vm(thresh_idx);
                    end
                    tmp_spike_thresholds = cat(1, tmp_spike_thresholds, thresh_val);
%                     if ~l_no_spikes(i_trl)
%                         plot(spk_raw_Vm, 'k')
%                         plot(thresh_idx, spk_raw_Vm(thresh_idx), 'ko', 'markerfacecolor', 'k')
%                     end
                    
                    % get spike reset Vm
                    if  (i_spk < numel(spike_idx)) && (numel(spike_idx)>=2)
                        spk_raw_Vm = raw_Vm(i_trl, spike_idx(i_spk) : spike_idx(i_spk+1));
                        [ahp_val, ahp_idx] = min(spk_raw_Vm);
                        tmp_spike_reset = cat(1, tmp_spike_reset, ahp_val);
                        
%                         if ~l_no_spikes(i_trl)
%                             plot(spk_raw_Vm, 'b')
%                             plot(ahp_idx, spk_raw_Vm(ahp_idx), 'bo', 'markerfacecolor','b')
%                         end
                    end
                end
            end
            
            % compute trial wide stats
            if ~isempty(tmp_spike_thresholds)
                thresholds_first_last(i_trl,1) = tmp_spike_thresholds(1);
            end
            if numel(tmp_spike_thresholds)>2
                thresholds_first_last(i_trl,2) = tmp_spike_thresholds(end);
                resets_first_last(i_trl,1) = tmp_spike_reset(1); % need 2 spikes to define an AHP
            end
            if numel(tmp_spike_thresholds)>3
                resets_first_last(i_trl,2) = tmp_spike_reset(end); % need 3 spikes to define "last" AHP
            end
        end
%         drawnow
        
        % aggregate Vthresh and Vreset data across Icmd conditions, keeping
        % track of first/last spike
        Vthresh_dat = [];
        Vreset_dat = [];
        for pA = unique(trl_Iclamp_cmd_pA)
            l_trl = trl_Iclamp_cmd_pA == pA;
            Vthresh_dat = cat(1, Vthresh_dat, [pA, nanmean(thresholds_first_last(l_trl,:),1)]);
            Vreset_dat = cat(1, Vreset_dat, [pA, nanmean(resets_first_last(l_trl,:),1)]);
        end
        dat.dcsteps.Vthresh{i_ch} = Vthresh_dat;
        dat.dcsteps.Vreset{i_ch} = Vreset_dat;
        
    end % for i_ch

end

function dat = unpack_vclamp_trains(dat, exinfo, hidx, FORCE_EPSC)
    %
    % NEED TO DEFINE THE FOLLOWING:
    %
    % dat.expt.[stimType].raw.snips: pA trace surrounding the EPSC for HS1,2 [Npulses x Ntime x Nsweeps]
    % dat.expt.[stimType].stats.EPSCamp: EPSC amplitude for HS1,2 [Npulses x 1 x Nsweeps], when raw current is inward
    % dat.expt.[stimType].stats.IPSCamp: IPSC amplitude for HS1,2 [Npulses x 1 x Nsweeps], when raw current is outward
    % dat.expt.[stimType].realTrialNum: the actual trial number for each of the sweeps
    % dat.expt.[stimType].tdict: the trialtype dictionary entry for this condition
    % dat.expt.[stimType].pOnTimes: the time of each pulse (in seconds). Includes the recovery pulse
    %
    % dat.qc.Rs: series resistance for HS1,2 [1xNtrials]
    % dat.qc.p1amp: amplitude for first pulse of each sweep for HS1,2 [1xNtrials]
    % dat.qc.verr: measured holding potential for each sweep, HS1,2 [1xNtrials]
    % dat.qc.instNoise:  the noise in the holding current which indicates changes in recording quality [1xNtrials]
    
    
    if strncmpi(dat.info.fid.vclamp(end-8:end), [filesep,'none.abf'], 9)
        return
    end
    fprintf('  Unpacking Vclamp trains\n')
    ax = abfobj(dat.info.fid.vclamp);
    
    % fill out the info struct
    dat.info.sampRate.vclamp = ax.head.sampRate;
    dat.info.sweepLength.vclamp = size(ax.dat, 1) ./ ax.head.sampRate;
    
    % figure out which channels were turned on, and set to Vclamp
    Vclamp_idx = strcmpi(ax.head.recChUnits, 'pA'); % units are telegraphed from Multiclamp and more accurate than "names"
    secCh = cellfun(@any, regexpi(ax.head.recChNames, '_sec'));
    Vclamp_idx = Vclamp_idx & ~secCh;
    Vclamp_names = ax.head.recChNames(Vclamp_idx);
    
    
    % extract raw current traces.
    preTime_samps = ceil(dat.info.pretime.vclamp .* dat.info.sampRate.vclamp);
    postTime_samps = ceil(dat.info.posttime.vclamp .* dat.info.sampRate.vclamp);
    
    
    % figure out what kind of stimuli were presented (e.g., simple trains,
    % recovery trains, RITs). Make a tdict
    tdict = outerleave(ax.dat(:,ax.idx.Laser,:), dat.info.sampRate.vclamp, false);
    Nconds = size(tdict.conds,1);

    % iterate over stimulus types aggregating data.
    for i_cond = 1:Nconds
        
        % determine the name of the condition, which will be used as the
        % name of a field in the output structure
        condname = make_condition_name(tdict, i_cond);
        
        % store some stuff from the tdict struct
        dat.expt.(condname).realTrialNum = repmat({find(tdict.trlList == i_cond)}, 1,2);
        dat.expt.(condname).tdict = tdict.conds(i_cond,:);
        
        % determine the pulse on times
        Nsweeps = numel(dat.expt.(condname).realTrialNum{1});
        laser_cmd = ax.dat(:,ax.idx.Laser, dat.expt.(condname).realTrialNum{1});
        laser_cmd = permute(laser_cmd, [3,1,2]);% [Nsweeps, Ntime]
        thresh = tdict.conds(i_cond, 1) .* 0.8;
        aboveThresh = laser_cmd > thresh;
        crossing_on_idx = [false(Nsweeps,1), diff(aboveThresh,1,2)] == 1;
        
        % make sure pulse on times are identical from sweep to sweep
        Npulses = unique(sum(crossing_on_idx, 2));
        assert(numel(Npulses) == 1, 'ERROR: Npulses changes from sweep to sweep for the same Ttype')
        tmp = sum(crossing_on_idx, 1);
        assert(sum(tmp==Nsweeps) == Npulses, 'ERROR: pON idicies are not consistent from sweep to sweep')
        
        % since we now know the stimuli are all the same, use the first
        % sweep as the template.
        crossing_on_idx = find(crossing_on_idx(1,:));
        dat.expt.(condname).pOnTimes = (crossing_on_idx-1) ./ dat.info.sampRate.vclamp; % seconds from sweep start
        
        % loop through channels
        for i_ch = 1:2
            
            % make sure data are present for this recording channel and
            % initialize the outputs
            HSname = sprintf('HS%d_', i_ch);
            HSpresent = strncmp(Vclamp_names, HSname, numel(HSname));
            assert(sum(HSpresent)<=1, 'ERROR: too many matches');
            if ~any(HSpresent)
                dat.expt.(condname).raw.snips{i_ch} = [];
                dat.expt.(condname).stats.EPSCamp{i_ch} = [];
                dat.expt.(condname).stats.IPSCamp{i_ch} = [];
                dat.expt.(condname).stats.latency{i_ch} = [];
                continue
            else
                HSname = Vclamp_names{HSpresent};  % need to modify in cases where Clampex thinks Iclamp but multiclamp set to Vclamp
                dat.expt.(condname).raw.snips{i_ch} = nan(Npulses, preTime_samps+postTime_samps+1, Nsweeps);
           end
            
            
            % iterate through the pulses and store the snippets
            for i_pulse = 1:Npulses
                idx = crossing_on_idx(i_pulse);
                idx = idx-preTime_samps : idx+postTime_samps;
                sweep_idx = dat.expt.(condname).realTrialNum{i_ch};
                tmp = ax.dat(idx, ax.idx.(HSname), sweep_idx); % Nth pulse, all sweeps
                dat.expt.(condname).raw.snips{i_ch}(i_pulse,:,:) = permute(tmp, [2,1,3]);
            end
            
            % deal with deleted sweeps here. things that need to get
            % culled:
            %  dat.expt.(condname).raw.snips{i_ch}
            %  dat.expt.(condname).realTrialNum
            rmSweep_string = exinfo{hidx.(sprintf('VC_rm_swp_HS%d', i_ch))};
            if ~any(isnan(rmSweep_string))
                badSweeps = eval(['[',rmSweep_string,']']);
                [~, l_badSweeps] = intersect(dat.expt.(condname).realTrialNum{i_ch}, badSweeps);
                if ~isempty(l_badSweeps)
                    fprintf('    Deleting %d sweeps from %s site %s ch %d\n',...
                               numel(l_badSweeps), dat.info.mouseName, dat.info.siteNum, i_ch);
                    dat.expt.(condname).realTrialNum{i_ch}(l_badSweeps) = [];
                    dat.expt.(condname).raw.snips{i_ch}(:,:,l_badSweeps) = [];
                end
            end
            
            
            % subtract off the background from all the pulses
            bkgnd = mean(dat.expt.(condname).raw.snips{i_ch}(:,1:preTime_samps,:), 2);
            dat.expt.(condname).raw.snips{i_ch} = bsxfun(@minus, dat.expt.(condname).raw.snips{i_ch}, bkgnd);
            
            
            % analyze the snippets and determine the peak current (EPSC or IPSC)
            nsamps_latency = round(1.75e-3 .* dat.info.sampRate.vclamp);
            nsamps_pulse = round(5.75e-3 .* dat.info.sampRate.vclamp);
            sign_idx = (preTime_samps + nsamps_latency) : (preTime_samps + nsamps_pulse);
            psc_sign = sign(sum(mean(dat.expt.(condname).raw.snips{i_ch}(1,sign_idx,:),3)));
            [peak_pA, peak_tt] = get_peak_psc(dat.expt.(condname).raw.snips{i_ch}, psc_sign, dat, 'vclamp');
            dat.expt.(condname).stats.latency{i_ch} = peak_tt;
            if psc_sign == 1 && ~FORCE_EPSC  %IPSP
                    dat.expt.(condname).stats.EPSCamp{i_ch} = [];
                    dat.expt.(condname).stats.IPSCamp{i_ch} = peak_pA;
            elseif psc_sign == -1 || FORCE_EPSC %EPSC
                    dat.expt.(condname).stats.EPSCamp{i_ch} = peak_pA;
                    dat.expt.(condname).stats.IPSCamp{i_ch} = [];
            elseif isnan(psc_sign)
                    dat.expt.(condname).stats.EPSCamp{i_ch} = [];
                    dat.expt.(condname).stats.IPSCamp{i_ch} = [];
            else
                    error('unexpected psc_sign')
            end
            
            % compute the rise times
            dat.expt.(condname).stats.rise_time{i_ch} = compute_rise_times(dat, psc_sign, i_ch, condname, FORCE_EPSC);

        end

    end
    
    %
    % PULL IN THE QUALITY CONTROL DATA
    %
    qc = ax.getRa;
    for i_ch = 1:2
        nTrialsExpt = size(ax.dat,3);
        dat.qc.Rs{i_ch} = nan(1,1,nTrialsExpt);
        dat.qc.Rinput{i_ch} = nan(1,1,nTrialsExpt);
        dat.qc.verr{i_ch} = nan(1,1,nTrialsExpt);
        dat.qc.p1amp{i_ch} = nan(1,1,nTrialsExpt);
        dat.qc.instNoise{i_ch} = nan(1,1,nTrialsExpt);
        
        
        % make sure data are present for this recording channel and
        % initialize the outputs
        HSname = sprintf('HS%d_', i_ch);
        HSpresent = strncmp(Vclamp_names, HSname, numel(HSname));
        if ~any(HSpresent)
            continue
        end
        
        % extract things from the Ra structure
        HSname = Vclamp_names{HSpresent};  % need to modify in cases where Clampex thinks Iclamp but multiclamp set to Vclamp
        idx = strcmpi(qc.chNames, HSname);
        dat.qc.Rs{i_ch} = qc.Ra(1,idx,:);
        dat.qc.verr{i_ch} = qc.Verr(1,idx,:);
        dat.qc.Rinput{i_ch} = qc.Rinput(1,idx,:);
        
        % dat.qc.Rs, Rinput, and dat.qc.verr contain info for all sweeps. Delete the
        % ones that are not analyzed.
        rmSweep_string = exinfo{hidx.(sprintf('VC_rm_swp_HS%d', i_ch))};
        if ~any(isnan(rmSweep_string))
            l_badSweeps = eval(['[',rmSweep_string,']']);
            dat.qc.Rs{i_ch}(l_badSweeps) = nan;
            dat.qc.Rinput{i_ch}(l_badSweeps) = nan;
            dat.qc.verr{i_ch}(l_badSweeps) = nan;
        end
        
        % Extract the p1 amp and Instrument Noise from the existing data in
        % a for-loop.
        conds = fieldnames(dat.expt);
        for i_cond = 1:Nconds
            real_trl_nums = dat.expt.(conds{i_cond}).realTrialNum{i_ch};
            if ~isempty(real_trl_nums)
                % deal with the p1 amps
                epsc = dat.expt.(conds{i_cond}).stats.EPSCamp{i_ch};
                ipsc = dat.expt.(conds{i_cond}).stats.IPSCamp{i_ch};
                assert(sum([isempty(epsc), isempty(ipsc)])==1, 'ERROR: was not expecting ipsc and epsc to both be defined')
                if ~isempty(epsc)
                    psc = epsc;
                elseif ~isempty(ipsc)
                    psc = ipsc;
                end
                dat.qc.p1amp{i_ch}(real_trl_nums) = psc(1,1,:);
                
                
                % deal with the instrument noise
                snips = dat.expt.(conds{i_cond}).raw.snips{i_ch};
                snips = snips(:, 1:preTime_samps, :);
                instNoise = std(snips, [], 2);
                instNoise  = mean(instNoise, 1);
                dat.qc.instNoise{i_ch}(real_trl_nums) = instNoise;
            end
        end
    end
end


function dat = unpack_iclamp_trains(dat, exinfo, hidx) % added args are for sweep deletion.
    %
    % NEED TO DEFINE THE FOLLOWING:
    %
    % dat.iclamp.[stimType].raw.snips: mV trace surrounding the EPSP for HS1,2 [Npulses x Ntime x Nsweeps]
    % dat.iclamp.[stimType].stats.EPSPamp: EPSP amplitude for HS1,2 [Npulses x 1 x Nsweeps], when raw current is inward
    % dat.iclamp.[stimType].stats.IPSPamp: IPSP amplitude for HS1,2 [Npulses x 1 x Nsweeps], when raw current is outward
    % dat.iclamp.[stimType].realTrialNum: the actual trial number for each of the sweeps
    % dat.iclamp.[stimType].tdict: the trialtype dictionary entry for this condition
    % dat.iclamp.[stimType].pOnTimes: the time of each pulse (in seconds). Includes the recovery pulse
    
    if strncmpi(dat.info.fid.iclamp(end-8:end), [filesep,'none.abf'], 9)
        return
    end
    fprintf('  Unpacking Iclamp trains\n')
    ax = abfobj(dat.info.fid.iclamp);
    
    % fill out the info struct
    dat.info.sampRate.iclamp = ax.head.sampRate;
    dat.info.sweepLength.iclamp = size(ax.dat, 1) ./ ax.head.sampRate;
    
    % figure out which channels were turned on, and set to Iclamp
    Iclamp_idx = strcmpi(ax.head.recChUnits, 'mV'); % units are telegraphed from Multiclamp and more accurate than "names"
    secCh = cellfun(@any, regexpi(ax.head.recChNames, '_sec'));
    Iclamp_idx = Iclamp_idx & ~secCh;
    Iclamp_names = ax.head.recChNames(Iclamp_idx);
    
    
    % figure out what kind of stimuli were presented (e.g., simple trains,
    % recovery trains, RITs). Make a tdict
    tdict = outerleave(ax.dat(:,ax.idx.Laser,:), dat.info.sampRate.iclamp, false);
    Nconds = size(tdict.conds,1);

    % iterate over stimulus types aggregating data.
    for i_cond = 1:Nconds
        
        % determine the name of the condition, which will be used as the
        % name of a field in the output structure
        condname = make_condition_name(tdict, i_cond);
        
        % store some stuff from the tdict struct
        dat.iclamp.(condname).realTrialNum = repmat({find(tdict.trlList == i_cond)}, 1,2);
        dat.iclamp.(condname).tdict = tdict.conds(i_cond,:);
        
        % determine the pulse on times
        Nsweeps = numel(dat.iclamp.(condname).realTrialNum{1});
        laser_cmd = ax.dat(:,ax.idx.Laser, dat.iclamp.(condname).realTrialNum{1});
        laser_cmd = permute(laser_cmd, [3,1,2]);% [Nsweeps, Ntime]
        thresh = tdict.conds(i_cond, 1) .* 0.8;
        aboveThresh = laser_cmd > thresh;
        crossing_on_idx = [false(Nsweeps,1), diff(aboveThresh,1,2)] == 1;
        
        % make sure pulse on times are identical from sweep to sweep
        Npulses = unique(sum(crossing_on_idx, 2));
        assert(numel(Npulses) == 1, 'ERROR: Npulses changes from sweep to sweep for the same Ttype')
        tmp = sum(crossing_on_idx, 1);
        assert(sum(tmp==Nsweeps) == Npulses, 'ERROR: pON idicies are not consistent from sweep to sweep')
        
        % since we now know the stimuli are all the same, use the first
        % sweep as the template.
        crossing_on_idx = find(crossing_on_idx(1,:));
        dat.iclamp.(condname).pOnTimes = (crossing_on_idx-1) ./ dat.info.sampRate.iclamp; % seconds from sweep start
        
        % determine the pre/post time in samples. Allow the post time to
        % float depending on the ISI
        preTime_samps = ceil(dat.info.pretime.iclamp .* dat.info.sampRate.iclamp);
        if dat.iclamp.(condname).tdict(5)==0; % RIT and normal trains get postTimes assigned differently
            ISI_sec = 1 ./ dat.iclamp.(condname).tdict(3);
            postTime_sec = max([dat.info.posttime.iclamp, ISI_sec-0.005]);
        else
            postTime_sec = dat.info.posttime.iclamp;
        end
        dat.iclamp.(condname).posttime_sec = postTime_sec; % store for later analysis functions
        postTime_samps = ceil(postTime_sec .* dat.info.sampRate.iclamp);
        
        % loop through channels
        for i_ch = 1:2
            
            % make sure data are present for this recording channel and
            % initialize the outputs
            HSname = sprintf('HS%d_', i_ch);
            HSpresent = strncmp(Iclamp_names, HSname, numel(HSname));
            assert(sum(HSpresent)<=1, 'ERROR: too many matches');
            if ~any(HSpresent)
                dat.iclamp.(condname).raw.snips{i_ch} = [];
                dat.iclamp.(condname).stats.EPSPamp{i_ch} = [];
                dat.iclamp.(condname).stats.IPSPamp{i_ch} = [];
                continue
            else
                HSname = Iclamp_names{HSpresent};  % need to modify in cases where Clampex thinks Iclamp but multiclamp set to Vclamp
                dat.iclamp.(condname).raw.snips{i_ch} = nan(Npulses, preTime_samps+postTime_samps+1, Nsweeps);
           end
            
            
            % iterate through the pulses and store the snippets
            for i_pulse = 1:Npulses
                idx = crossing_on_idx(i_pulse);
                idx = idx-preTime_samps : idx+postTime_samps;
                sweep_idx = dat.iclamp.(condname).realTrialNum{i_ch};
                tmp = ax.dat(idx, ax.idx.(HSname), sweep_idx); % Nth pulse, all sweeps
                dat.iclamp.(condname).raw.snips{i_ch}(i_pulse,:,:) = permute(tmp, [2,1,3]);
            end
            
            
            
            % subtract off the background from all the pulses. DEFINE THE
            % BACKROUND AS THE TIME BEFORE ONLY THE FIRST PULSE. This
            % differes from the Vclamp experiment, in which each pulse gets
            % baselined.
            baseline_samps = ceil(0.150 .* dat.info.sampRate.iclamp);
            idx = crossing_on_idx(1);
            idx = idx-baseline_samps : idx-1;
            sweep_idx = dat.iclamp.(condname).realTrialNum{i_ch};
            tmp = ax.dat(idx, ax.idx.(HSname), sweep_idx); % baseline data, all sweeps
            bkgnd = mean(tmp, 1);
            dat.iclamp.(condname).raw.snips{i_ch} = bsxfun(@minus, dat.iclamp.(condname).raw.snips{i_ch}, bkgnd);
            
            % store the background Vm for debugging and such.
            dat.iclamp.(condname).raw.bkgndVm{i_ch} = bkgnd;
            
            
            % analyze the snippets and determine the peak current (EPSP or IPSP)
            ch_is_valid = dat.info.HS_is_valid_Iclamp(i_ch);
            psp_sign = sign(sum(mean(dat.iclamp.(condname).raw.snips{i_ch}(1,:,:),3)));
            [peak_pA, peak_tt] = get_peak_psc(dat.iclamp.(condname).raw.snips{i_ch}, psp_sign, dat, 'iclamp');
            dat.iclamp.(condname).stats.latency{i_ch} = peak_tt;
            if ~ch_is_valid
                dat.iclamp.(condname).stats.EPSPamp{i_ch} = [];
                dat.iclamp.(condname).stats.IPSPamp{i_ch} = [];
            elseif psp_sign == -1  %IPSP
                dat.iclamp.(condname).stats.EPSPamp{i_ch} = [];
                dat.iclamp.(condname).stats.IPSPamp{i_ch} = peak_pA;
                warning('found an IPSP')
            elseif psp_sign == 1 %EPSP
                dat.iclamp.(condname).stats.EPSPamp{i_ch} = peak_pA;
                dat.iclamp.(condname).stats.IPSPamp{i_ch} = [];
            elseif isnan(psp_sign) || psp_sign == 0
                dat.iclamp.(condname).stats.EPSPamp{i_ch} = [];
                dat.iclamp.(condname).stats.IPSPamp{i_ch} = [];
                channel_is_defined = str2double(exinfo{hidx.(['IC_HS', num2str(i_ch), '_valid'])});
                assert(channel_is_defined == 0, 'ERROR: Channel is defined but psp sign is not identified')
            end
                      
        end

    end
    
end


function dat = smoothP1amp(dat, N)
    
    if ~(isfield(dat, 'qc') && isfield(dat.qc, 'p1amp'))
        dat.qc.p1amp_norm = [];
        return
    end
    
    for i_ch = 1:2
        
        p1 = squeeze(dat.qc.p1amp{i_ch});
        
        initmean = nanmean(p1(1:5));
        endmean = nanmean(p1(end-4:end));
        tmp = [ones(1,N).*initmean , p1', ones(1,N+2).*endmean];
        normfact = nan(size(tmp));
        for i_swp = 1:(numel(normfact)-N)
            normfact(i_swp) = nanmean(tmp(i_swp:(i_swp+(N-1))));
        end
        normfact(1:N-3)=[];
        normfact(numel(p1)+1:end) = [];
        
        dat.qc.p1amp_norm{i_ch} = normfact;
    end
end

function [peak_pA, peak_tt] = get_peak_psc(snips, psc_sign, dat, METHOD)
    
    % make a time vector
    tt = [0:size(snips,2)-1] ./ dat.info.sampRate.vclamp;
    tt = tt - dat.info.pretime.vclamp;
    switch METHOD
        case 'vclamp'
            peak_window = (tt > 500e-6) & (tt < 6e-3);
        case 'iclamp'
            peak_window = (tt > 500e-6) & (tt < 15e-3);
    end
    
    % define a window around the peak to average the data
    fullwindow = 300e-6 .* dat.info.sampRate.vclamp;
    halfwidth = ceil(fullwindow./2);
    avg_window = -halfwidth : 1 : halfwidth;
    
    % flip the sign if need be, so that we can always use max
    if psc_sign == -1
        snips = -snips;
    end
    
    % iterate through the sweeps and pulses and define the peak, and latency
    [peak_pA, peak_tt] = deal(nan(size(snips,1), 1, size(snips,3)));
    for i_sweep = 1:size(snips,3)
        for i_pulse = 1:size(snips,1)
            tmp_snip = squeeze(snips(i_pulse,:,i_sweep));
            maxval = max(tmp_snip(peak_window));
            eq2max = tmp_snip == maxval;
            peak_idx = find(eq2max & peak_window, 1, 'first');
            
            peak_pA(i_pulse, 1, i_sweep) = mean(tmp_snip(avg_window+peak_idx));
            peak_tt(i_pulse, 1, i_sweep) = tt(peak_idx);
        end
    end
    assert(~any(isnan(peak_pA(:))), 'ERROR: found some nans');
    assert(~any(isnan(peak_tt(:))), 'ERROR: found some nans');
        
end

function rise_times = compute_rise_times(dat, psc_sign, i_ch, condname, FORCE_EPSC)

            % compute the 5-95% time time here. multipy the snips by the
            % psc_sign to make everything positive.
            Npulses = size(dat.expt.(condname).raw.snips{i_ch}, 1);
            Nsweeps = size(dat.expt.(condname).raw.snips{i_ch}, 3);
            rise_times = nan(Npulses,1,Nsweeps); % pre-allocate
            for i_swp = 1:Nsweeps
                for i_pulse = 1:Npulses
                    if psc_sign == 1 && ~FORCE_EPSC %IPSP
                        peak_val = dat.expt.(condname).stats.IPSCamp{i_ch}(i_pulse,:,i_swp);
                    elseif psc_sign == -1  || FORCE_EPSC %EPSC
                        peak_val = dat.expt.(condname).stats.EPSCamp{i_ch}(i_pulse,:,i_swp);
                    else
                        error('unexpected psc_sign')
                    end
                    wf = dat.expt.(condname).raw.snips{i_ch}(i_pulse,:,i_swp);
                    wf = wf .* psc_sign; % make sure all the values are positive.
                    peak_idx = (dat.info.pretime.vclamp + dat.expt.(condname).stats.latency{i_ch}(i_pulse,1,i_swp)) .* dat.info.sampRate.vclamp;
                    peak_idx = round(peak_idx);
                    laser_idx = round(dat.info.pretime.vclamp .* dat.info.sampRate.vclamp);
                    idx_95 = find(wf(laser_idx:peak_idx) <= (peak_val .* 0.95), 1, 'last');
                    idx_05 = find(wf(laser_idx:peak_idx) <= (peak_val .* 0.05), 1, 'last');
                    if ~isempty(idx_95) && ~isempty(idx_05)
                        rise_times(i_pulse,:,i_swp) = (idx_95 - idx_05) ./ dat.info.sampRate.vclamp;
                    else
                        rise_times(i_pulse,:,i_swp) = nan;
                    end
                end
            end
            
end

function [Ao, tau, C, tau_guess] = fit_exp1(vm, tt)
    
    % guess values for Ao, tau, C
    total_decay = vm(1)-vm(end);
    idx_to_63prcnt = find(vm<(vm(1)-(total_decay*0.63)), 1, 'first');
    tau_guess = tt(idx_to_63prcnt);
    C_guess = mean(vm(end-5000:end));
    Ao_guess = total_decay;
    
    % run the fminsearch, minimize sse
    [fitted_params] = fminsearch(@sse_err, [Ao_guess, tau_guess, C_guess]);
    Ao = fitted_params(1);
    tau = fitted_params(2);
    C = fitted_params(3);
    
    figure, hold on
    plot(tt, vm, 'k')
    vm_pred = Ao .* exp(-tt./tau) + C;
    plot(tt, vm_pred, 'r')
    drawnow
    function err = sse_err(params)
        
        fit = params(1) .* exp(-tt./params(2)) + params(3);
        err = sum((vm - fit).^2);
        
    end
end

function name = make_condition_name(tdict, i_cond)
    pamp = round(tdict.conds(i_cond,1), 1) .* 10; % in tens of mV
    pwidth = round(tdict.conds(i_cond,2), 6) .* 1e6; % in us
    tfreq = round(tdict.conds(i_cond,3)); % in Hz
    trecov = round(tdict.conds(i_cond,4), 3) .* 1000; % in ms
    tRIT = tdict.conds(i_cond,5); %  'version number'
    if tRIT ~= 0
        name = sprintf('RITv%d', tRIT);
    elseif trecov ~= 0
        name = sprintf('a%d_w%d_f%d_r%d', pamp, pwidth, tfreq, trecov);
    else
        name = sprintf('a%d_w%d_f%d', pamp, pwidth, tfreq);
    end
end

