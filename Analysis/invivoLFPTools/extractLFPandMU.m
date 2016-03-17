function dat = extractLFPandMU(blk, lfp_ch_idx, sampFreq_lfp, NOISEMETHOD, STIMTYPE, NSX)
   
    lfp_data = blk.ras(:, lfp_ch_idx);
    sampFreq_nsx = blk.sum.(NSX).MetaTags.SamplingFreq;

    % design and implement a notch filter
    lp_cutoff = 60.5;
    hp_cutoff = 59.5;
    order = 2;
    Wn = [hp_cutoff, lp_cutoff] ./ (sampFreq_nsx ./ 2);
    [B_notch, A_notch] = butter(ceil(order/2), Wn, 'stop');
    
    
    % design and implement a band pass filter
    bp_lp_cutoff = 250; % guards against alaising during downsampling
    bp_hp_cutoff = 2;
    order = 6;
    Wn = [bp_hp_cutoff, bp_lp_cutoff] ./ (sampFreq_nsx ./ 2);
    [B_bp, A_bp] = butter(ceil(order/2), Wn);
    assert(bp_lp_cutoff <= sampFreq_lfp./5, 'ERROR: LFP sampRate and filter are incompatable');
    
    parfor i_idx = 1: numel(lfp_data);
        
        % band pass filter
        lfp_data{i_idx} = filtfilt(B_bp, A_bp, lfp_data{i_idx})
        
        % remove 60Hz noise
        switch NOISEMETHOD
            case 'subtract'
                N = numel(lfp_data{i_idx});
                lfp_data{i_idx} = rmhum_2(lfp_data{i_idx}, sampFreq_nsx, 1, N, 60);
            case 'filter'
                lfp_data{i_idx} = filtfilt(B_notch, A_notch, lfp_data{i_idx})
            case 'none'
                % do nothing
        end
        
    end
    
    
    % determine which trials have the different pulse train conditions (in the
    % case that things are interleaved)
    disp(' Getting time locked data for LFP')
    lfpSnips = getTimeLockedLFP(lfp_data, blk, STIMTYPE, NSX);
    clear lfp_data
    
    
    ptypes = fieldnames(lfpSnips);
    i_bad = strcmpi(ptypes, 'preTime') | strcmpi(ptypes, 'postTime') | strcmpi(ptypes, 'pulseOnT_sec');
    ptypes(i_bad) = [];
    
    
    disp(' Down sampling LFP data')
    % downsample
    for i_cond = 1:numel(ptypes);
        
        tmp_data = lfpSnips.(ptypes{i_cond});
        parfor i_ch = 1:numel(tmp_data);
            
            % downsample
            NdownSample = sampFreq_nsx ./ sampFreq_lfp;
            assert(rem(NdownSample,1)==0, 'ERROR: new sampFreq does not divide evenly into old rate')
            tmp_data{i_ch} = downsample(tmp_data{i_ch}', NdownSample)'; % notice the transpose(s)
            
        end
        
        % put things back
        lfpSnips.(ptypes{i_cond}) = tmp_data;
        
    end
    
    dat.lfpsnips = lfpSnips;
    clear lfpSnips;
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Calculate the multi-unit spike times
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(' Preparing multi-unit spike times (parfor loops BP filter)')
    RMS_MULTIPLE = 4.5;
    
    mu_data = blk.ras(:, lfp_ch_idx);
    
    % design and implement a band pass filter
    bp_lp_cutoff = 6000;
    bp_hp_cutoff = 700;
    order = 12;
    Wn = [bp_hp_cutoff, bp_lp_cutoff] ./ (sampFreq_nsx ./ 2);
    [B_mu, A_mu] = butter(ceil(order/2), Wn);
    assert(bp_hp_cutoff>100, 'ERROR: may need notch filter')
    
    parfor i_idx = 1:numel(mu_data)
        mu_data{i_idx} = filtfilt(B_mu, A_mu, mu_data{i_idx});
    end
    
    % sync to the stimulus onset, and baseline subtract
    disp(' Getting time locked data for MU spike times')
    mu_snips = getTimeLockedLFP(mu_data, blk, STIMTYPE, NSX);
    clear mu_data;
    
    
    % identify the different pulse types in play
    ptypes = fieldnames(mu_snips);
    i_bad = strcmpi(ptypes, 'preTime') | strcmpi(ptypes, 'postTime') | strcmpi(ptypes, 'pulseOnT_sec');
    ptypes(i_bad) = [];
    Nptypes = numel(ptypes);
    
    preSamps = round(mu_snips.preTime .* sampFreq_nsx);
    preTime = mu_snips.preTime;
    
    disp(' Extracting spike times')
    spikeTimes = [];
    for i_type = 1:Nptypes
        
        tmp_data = mu_snips.(ptypes{i_type});
        tmp_spiketimes = {};
        tmp_Ntrials = {};
        tmp_spikeTrialID = {};
        
        for i_ch = 1:numel(mu_snips.(ptypes{i_type}))
            
            if isempty(tmp_data{i_ch})
                tmp_spiketimes{i_ch} = [];
                tmp_Ntrials{i_ch} = [];
                
            else
                
                % estimate RMS during the pre-stimulus period
                bkgnd = tmp_data{i_ch}(:,1:preSamps);
                rmsVal = rms(bkgnd,2);
                threshold = rmsVal .* RMS_MULTIPLE;
                
                % full-wave rectify the data, and impose a threshold. Find the
                % crossings and store the spike times.
                aboveThresh = bsxfun(@lt, tmp_data{i_ch}, -threshold);
                diffAboveThresh = [zeros(size(aboveThresh,1),1), diff(aboveThresh, 1, 2)];
                aboveThresh = diffAboveThresh == 1;
                
                % store the spike times (synced to stim onset) but put all the spike into a single basket.
                % I'll store the number of trials, so that I can easily compute
                % PSTHs later
                tt = [0:size(aboveThresh,2)-1] ./ sampFreq_nsx;
                tt = tt - preTime;
                
                [trlID, t_idx]= find(aboveThresh);
                tmp_spiketimes{i_ch} = tt(t_idx);
                tmp_spikeTrialID{i_ch} = trlID;
                tmp_Ntrials{i_ch} = size(aboveThresh, 1);
                
            end
        end
        
        % put the data back where they belong
        spikeTimes.(ptypes{i_type}).tt = tmp_spiketimes;
        spikeTimes.(ptypes{i_type}).trialIDs = tmp_spikeTrialID;
        spikeTimes.(ptypes{i_type}).Ntrials = tmp_Ntrials;
        
    end
    
    
    dat.spikeTimes = spikeTimes;
    dat.spikeTimes.preTime = mu_snips.preTime;
    dat.spikeTimes.postTime = mu_snips.postTime;
    dat.spikeTimes.pulseOnT_sec = mu_snips.pulseOnT_sec;
    
    dat.info.sampFreq_lfp = sampFreq_lfp;
    dat.info.rms_multiple = RMS_MULTIPLE;
    dat.info.nsxtype = NSX; % the version with the LFP continuous data
    dat.info.stimtype = STIMTYPE;% train or sinusoid
    dat.info.noisemethod = NOISEMETHOD;% 'subtract' or 'filter'