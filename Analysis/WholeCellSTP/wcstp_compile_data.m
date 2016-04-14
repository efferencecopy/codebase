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
% dat.qc.Rs: series resistance for HS1,2 [1xNtrials]
% dat.qc.p1amp: amplitude for first pulse of each sweep for HS1,2 [1xNtrials]
% dat.qc.vhold: measured holding potential for each sweep, HS1,2 [1xNtrials]
% 
% dat.dcsteps.Vm_raw: raw voltage traces during stimulus on period for HS1,2 [Nsweeps x Ntime]
% dat.dcsteps.Icmd: magnitude of current injection on a sweep by sweep basis. [Nsweeps x 1]
% 
% dat.expt.[stimType].raw.snips: pA trace surrounding the EPSC for HS1,2 [Npulses x Ntime x Nsweeps]
% dat.expt.[stimType].stats.EPSCamp: EPSC amplitude for HS1,2 [Npulses x 1 x Nsweeps], when raw current is inward
% dat.expt.[stimType].stats.IPSCamp: IPSC amplitude for HS1,2 [Npulses x 1 x Nsweeps], when raw current is outward
% dat.expt.[stimType].stats.latency: time to peak of PSC in seconds [Npulses x 1 x Nsweeps]
% dat.expt.[stimType].realTrialNum: the actual trial number for each of the sweeps
% dat.expt.[stimType].tdict: the trialtype dictionary entry for this condition
% 
% dat.info.opsin
% dat.info.fid.dcsteps
% dat.info.fid.vclamp
% dat.info.mouseName
% dat.info.siteNum
% dat.info.pretime.vclamp
% dat.info.pretime.dcsteps
% dat.info.posttime.vclamp
% dat.info.posttime.dcsteps
% dat.info.sampRate.vclamp
% dat.info.sampRate.iclamp



% fill out the info field and get things ready for the analysis
dat.info = defineInfo(exinfo, hidx, params);

% the function is verbose. Tell the user what's going on
fprintf('Analyzing data from mouse %s, site %s\n', dat.info.mouseName, dat.info.siteNum)

% unpack the DC current injection data set
dat = unpack_dc_injections(dat);


% unpack the Vclamp trains data set(s)
dat = unpack_vclamp_trains(dat);

end


function info = defineInfo(exinfo, hidx, params)
    info.fid.dcsteps = exinfo{hidx.ABFDCsteps};
    info.fid.vclamp = exinfo{hidx.ABFOptostim};

    info.opsin = exinfo{hidx.opsin};
    info.mouseName = exinfo{hidx.MouseName};
    info.siteNum = exinfo{hidx.Site};
    
    info.pretime.vclamp = params.pretime.vclamp;
    info.posttime.vclamp = params.posttime.vclamp;
    info.pretime.dcsteps = params.pretime.dcsteps;
    info.posttime.dcsteps = params.posttime.dcsteps;
end

function dat = unpack_dc_injections(dat)

    % dat.dcsteps.Vm_raw: raw voltage traces during stimulus on period for HS1,2 [Nsweeps x Ntime]
    % dat.dcsteps.Icmd: magnitude of current injection on a sweep by sweep basis. [Nsweeps x 1]

    global GL_DATPATH 
    
    fprintf('  Unpacking DC current injections. Mouse %s, file %s\n', dat.info.mouseName, dat.info.fid.dcsteps)
    fpath = strcat(GL_DATPATH, dat.info.mouseName, filesep, 'Physiology', filesep, dat.info.fid.dcsteps, '.abf');
    ax = abfobj(fpath);
    
    % fill out the info struct
    dat.info.sampRate.iclamp = ax.head.sampRate;
    
    % figure out which channels were turned on, and set to Iclamp
    Iclamp_idx = strcmpi(ax.head.recChUnits, 'mV');
    secCh = cellfun(@any, regexpi(ax.head.recChNames, '_sec'));
    Iclamp_idx = Iclamp_idx & ~secCh;
    Iclamp_names = ax.head.recChNames(Iclamp_idx);
    
    
    % extract raw voltage traces and spike times.
    preTime_samps = ceil(dat.info.pretime.dcsteps .* dat.info.sampRate.iclamp);
    postTime_samps = ceil(dat.info.posttime.dcsteps .* dat.info.sampRate.iclamp);
    for i_ch = 1:2
        
        % make sure data are present for this recording channel. If not,
        % store the data as an empty matrix.
        HSname = sprintf('HS%d_Vm', i_ch);
        HSpresent = strcmp(Iclamp_names, HSname);
        if ~any(HSpresent)
            dat.dcsteps.Vm_raw{i_ch} = [];
            dat.dcsteps.Icmd{i_ch} = [];
            continue
        end
        
        % figure out when the DC injection started and stoped
        Iclamp_name = sprintf('HS%d_Iclamp', i_ch);
        wf = permute(ax.wf(:,ax.idx.(Iclamp_name),:), [1,3,2]); %[Ntime x Nsweeps]
        minval = min(abs(wf), [], 1); % do the thresholding on the abs of the wf
        maxval = max(abs(wf), [], 1);
        
        thresh = (maxval - minval) .* 0.8;
        aboveThresh = bsxfun(@gt, abs(wf), thresh);
        crossing_on = cat(1, false(size(thresh)), diff(aboveThresh, 1, 1) == 1);
        crossing_off = cat(1, false(size(thresh)), diff(aboveThresh, 1, 1) == -1);
        
        % convert the onset/offset indicies to R/C notation
        [ridx_on,~] = find(crossing_on);
        ridx_on = unique(ridx_on);
        [ridx_off,~] = find(crossing_off);
        ridx_off = unique(ridx_off);
        assert(numel([ridx_on, ridx_off])==2, 'ERROR: did not identify the correct number of threshold crossings')
        
        % now make a list of indicies that correspond to stim on +/-
        % pre/post time
        stim_idx = (ridx_on - preTime_samps) : (ridx_off + postTime_samps);
        
        % pull out the data, store it in the output structure
        raw_Vm = permute(ax.dat(stim_idx, ax.idx.(HSname), :), [3,1,2]);
        dat.dcsteps.Vm_raw{i_ch} = raw_Vm; % [Nsweeps x Ntime]
        
        % pull out the Iclamp command data (current injection)
        stim_idx = (ridx_on+10) : (ridx_off-10); % make sure to only include the current injection portion and not the pre/post time
        raw_pA = permute(ax.wf(stim_idx, ax.idx.(Iclamp_name), :), [3, 1, 2]);
        dat.dcsteps.Icmd{i_ch} = mean(raw_pA,2);
        
    end

end

function dat = unpack_vclamp_trains(dat)
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
    % dat.qc.vhold: measured holding potential for each sweep, HS1,2 [1xNtrials]
    
    global GL_DATPATH 

    fprintf('  Unpacking Vclamp trains. Mouse %s, file %s\n', dat.info.mouseName, dat.info.fid.vclamp)
    fpath = strcat(GL_DATPATH, dat.info.mouseName, filesep, 'Physiology', filesep, dat.info.fid.vclamp, '.abf');
    ax = abfobj(fpath);
    
    % fill out the info struct
    dat.info.sampRate.vclamp = ax.head.sampRate;
    
    % figure out which channels were turned on, and set to Vclamp
    Vclamp_idx = strcmpi(ax.head.recChUnits, 'pA'); % units are telegraphed from Multiclamp and more accurate than "names"
    secCh = cellfun(@any, regexpi(ax.head.recChNames, '_sec'));
    Vclamp_idx = Vclamp_idx & ~secCh;
    Vclamp_names = ax.head.recChNames(Vclamp_idx);
    
    
    % extract raw voltage traces and spike times.
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
        dat.expt.(condname).realTrialNum = find(tdict.trlList == i_cond);
        dat.expt.(condname).tdict = tdict.conds(i_cond,:);
        
        % determine the pulse on times
        Nsweeps = numel(dat.expt.(condname).realTrialNum);
        laser_cmd = ax.dat(:,ax.idx.Laser, dat.expt.(condname).realTrialNum);
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
        
        preTime_samps = ceil(dat.info.pretime.vclamp .* dat.info.sampRate.vclamp);
        postTime_samps = ceil(dat.info.posttime.vclamp .* dat.info.sampRate.vclamp);
        
        
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
                dat.expt.(condname).stats.latency{i_ch} = [];continue
            else
                HSname = Vclamp_names{HSpresent};  % need to modify in cases where Clampex thinks Iclamp but multiclamp set to Vclamp
                dat.expt.(condname).raw.snips{i_ch} = nan(Npulses, preTime_samps+postTime_samps+1, Nsweeps);
           end
            
            
            % iterate through the pulses and store the snippets
            for i_pulse = 1:Npulses
                idx = crossing_on_idx(i_pulse);
                idx = idx-preTime_samps : idx+postTime_samps;
                sweep_idx = dat.expt.(condname).realTrialNum;
                tmp = ax.dat(idx, ax.idx.(HSname), sweep_idx); % Nth pulse, all sweeps
                dat.expt.(condname).raw.snips{i_ch}(i_pulse,:,:) = permute(tmp, [2,1,3]);
            end
            
            % subtract off the background from all the pulses
            bkgnd = mean(dat.expt.(condname).raw.snips{i_ch}(:,1:preTime_samps,:), 2);
            dat.expt.(condname).raw.snips{i_ch} = bsxfun(@minus, dat.expt.(condname).raw.snips{i_ch}, bkgnd);
            
            
            % analyze the snippets and determine the peak current (EPSC or IPSC)
            psc_sign = sign(sum(mean(dat.expt.(condname).raw.snips{i_ch}(1,:,:),3)));
            [peak_pA, peak_tt] = get_peak_psc(dat.expt.(condname).raw.snips{i_ch}, psc_sign, dat);
            dat.expt.(condname).stats.latency{i_ch} = peak_tt;
            switch psc_sign
                case 1  %IPSP
                    dat.expt.(condname).stats.EPSCamp{i_ch} = [];
                    dat.expt.(condname).stats.IPSCamp{i_ch} = peak_pA;
                case -1 %EPSC
                    dat.expt.(condname).stats.EPSCamp{i_ch} = peak_pA;
                    dat.expt.(condname).stats.IPSCamp{i_ch} = [];
            end
                      
        end

    end
    
    %
    % PULL IN THE QUALITY CONTROL DATA
    %
    qc = ax.getRa;
    for i_ch = 1:2
        
        % make sure data are present for this recording channel and
        % initialize the outputs
        HSname = sprintf('HS%d_', i_ch);
        HSpresent = strncmp(Vclamp_names, HSname, numel(HSname));
        if ~any(HSpresent)
            dat.qc.Rs{i_ch} = [];
            dat.qc.p1amp{i_ch} = [];
            dat.qc.vhold{i_ch} = [];
            continue
        end
        
        % extract things from the Ra structure
        HSname = Vclamp_names{HSpresent};  % need to modify in cases where Clampex thinks Iclamp but multiclamp set to Vclamp
        idx = strcmpi(qc.chNames, HSname);
        dat.qc.Rs{i_ch} = permute(qc.dat(1,idx,:), [3,1,2]);
        dat.qc.vhold{i_ch} = permute(qc.Verr(1,idx,:), [3,1,2]);
        
        
        % Extract the p1 amp from the existing data in a for-loop
        dat.qc.p1amp{i_ch} = nan(1, size(ax.dat, 3));
        conds = fieldnames(dat.expt);
        for i_cond = 1:Nconds
            real_trl_nums = dat.expt.(conds{i_cond}).realTrialNum;
            
            epsc = dat.expt.(conds{i_cond}).stats.EPSCamp{i_ch};
            ipsc = dat.expt.(conds{i_cond}).stats.IPSCamp{i_ch};
            assert(sum([isempty(epsc), isempty(ipsc)])==1, 'ERROR: was not expecting ipsc and epsc to both be defined')
            if ~isempty(epsc)
                psc = epsc;
            elseif ~isempty(ipsc)
                psc = ipsc;
            end
            
            dat.qc.p1amp{i_ch}(real_trl_nums) = squeeze(psc(1,1,:));
        end
    end
end

function [peak_pA, peak_tt] = get_peak_psc(snips, psc_sign, dat)
    
    % make a time vector
    tt = [0:size(snips,2)-1] ./ dat.info.sampRate.vclamp;
    tt = tt - dat.info.pretime.vclamp;
    peak_window = (tt > 500e-6) & (tt < 6e-3);
    
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

