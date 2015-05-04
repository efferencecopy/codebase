%% WHICH MICE SHOULD CONTRIBUTE? [ALL MICE]

% clear out the workspace
fin

% in = {Mouse Name, Site}

in = {'CH_141124_B', 1;...
    'CH_141215_B', 1;...
    'CH_141215_B', 2;...
    'CH_141215_D', 3;...
    'CH_141215_E', 1;...
    'CH_141215_E', 2;...
    'CH_141215_F', 1;...
    'CH_150105_A', 1;...
    'CH_150105_B', 1;...
    'CH_150105_B', 2;...
    'CH_150105_C', 1;...
    'CH_150105_D', 2;...
    'CH_150112_A', 1;...
    'CH_150112_A', 2;...
    'CH_150112_B', 1;...
    'CH_150112_B', 2;...
    'CH_150112_D', 1;...
    'CH_150112_D', 2;...
    'CH_150112_C', 1;...
    'CH_150112_C', 2;...
    'CH_150119_D', 1;...
    'CH_150119_D', 2};


%% WHICH MICE SHOULD CONTRIBUTE?  [GOOD MICE]

% Anlyze data sets that used 300us pulses, and that have a pure Na+ FV

% clear out the workspace
fin

% in = {Mouse Name, Site}

in = {
        'CH_141215_F', 1;...
        'CH_141215_E', 1;...
        'CH_141215_E', 2;...
        %'CH_150105_A', 1;...  % might get cut
        'CH_150105_B', 1;...
        'CH_150105_B', 2;...
        %'CH_150105_C', 1;...  % might get cut
        %'CH_150112_A', 1;...  % might get cut
        'CH_150112_A', 2;...
        'CH_150112_B', 1;...
        'CH_150112_B', 2;...
        'CH_150112_D', 1;...
        'CH_150112_D', 2;...
        'CH_150112_C', 1;...
        'CH_150119_C', 1;...
        'CH_150119_C', 2;...
        'CH_150119_D', 1;...
        'CH_150119_D', 2;...
        'CH_150119_B', 1;...
        'CH_150119_B', 2;...
        %'CH_150302_C', 1;... % different led powers, good for avg oChIEF waveform, but nothing else
        };


%% WHICH MICE SHOULD CONTRIBUTE?  [GOOD MICE FOR DIFFERENT PULSE AMP EXPTS]


% % clear out the workspace
% fin
% 
% % in = {Mouse Name, Site}
% 
% in = {
%       'CH_150112_C', 2;...
%       'CH_150302_C', 1;...
%       'CH_150302_A', 1;...
%       'CH_150302_D', 1;...
%       };


%% LOOP THOUGH EACH MOUSE AND CREATE THE NECESSARY RAW DATA TRACES


% grab the fiber volley pop excel workbook
fname = [GL_DOCUPATH, 'Other_workbooks', filesep, 'fiberVolleyCellList.xlsx'];
[~,txt, raw] = xlsread(fname);
raw(size(txt,1)+1:end, :) = [];
raw(:,size(txt,2)+1:end) = [];
channelIdx = cellfun(@(x) ~isempty(x), regexpi(raw(1,:), 'CH\d'));
opsinIdx = strcmpi(raw(1,:), 'opsin');
stimSiteIdx = strcmpi(raw(1,:), 'stim site');
clear txt

% do the analysis
Nexpts = size(in,1);
dat = {};
for i_ex = 1:Nexpts;
    
    % figure out what rows in the work book to pay attention to
    l_mouse = cellfun(@(x) ~isempty(x), regexp(raw(:,1), in{i_ex,1}));
    l_site = [false ; cell2mat(raw(2:end,2)) == in{i_ex,2}]; % add a leading 'false' to account for the header row in 'raw'
    l_expt = l_mouse & l_site;
    
    % run the analysis
    [dat{i_ex}, info{i_ex}] = fiberVolleyAnalysis(l_expt, raw, false);
    
    % enter a few other useful things into the structure
    info{i_ex}.mouse =  in{i_ex,1};
    info{i_ex}.opsin = unique(cell2mat(raw(l_expt, opsinIdx)), 'rows');
    info{i_ex}.ignoreChans = unique(cell2mat(raw(l_expt, channelIdx)), 'rows');
    info{i_ex}.stimSite = unique(cell2mat(raw(l_expt, stimSiteIdx)), 'rows');
    if ischar(info{i_ex}.stimSite); % deals with nans
        info{i_ex}.stimSite = str2num(info{i_ex}.stimSite);
    end
end



%% PULL OUT SNIPPETS OF DATA FOR EACH PULSE (ANALYZE THEM LATER)

% PEAK TO PEAK AMP
% INTEGRAL OF THE TRACE
% save the snippets for plotting

prePulseTime = 0.001; % in sec
postPulseTime = 0.009; % in sec

for i_ex = 1:Nexpts
    
    pTypes = fieldnames(dat{i_ex});
    Ntfs = numel(pTypes);
    for i_tf = 1:Ntfs
        
        conds = {'FV_Na', 'nbqx_apv_cd2_ttx', 'synapticTransmission', 'none', 'nbqx_apv', 'nbqx_apv_cd2'};
        for i_cond = 1:numel(conds)
            
            sampRate = info{i_ex}.(pTypes{i_tf}).(conds{i_cond}).sampRate;
            prePulseSamps = ceil(prePulseTime .* sampRate);
            postPulseSamps = ceil(postPulseTime .* sampRate);
            Nsamps = prePulseSamps + postPulseSamps + 1;
            
            Npulses = sum(info{i_ex}.(pTypes{i_tf}).(conds{i_cond}).pulseOn_idx);
            pulseOn_idx = find(info{i_ex}.(pTypes{i_tf}).(conds{i_cond}).pulseOn_idx);
            dat{i_ex}.(pTypes{i_tf}).snips.(conds{i_cond}) = {nan(Npulses, Nsamps), nan(Npulses, Nsamps)};
            for i_pulse = 1:Npulses
                
                snip_idx = pulseOn_idx(i_pulse)-prePulseSamps : 1 : pulseOn_idx(i_pulse)+postPulseSamps;
                
                for i_ch = 1:2;
                    
                    % deal with some exceptions
                    if ~info{i_ex}.ignoreChans(i_ch)
                        continue
                    end
                    
                    % a strange case where HS1 was set to Vclamp, and so
                    % HS2's data got put in the first column...
                    if strcmpi(info{i_ex}.mouse, 'CH_150105_D')
                        if i_ch == 1; error('something went wrong'); end
                        i_ch = 1;
                    end
                    
                    
                    % pull out the snippet, subtract off the baseline and
                    % store it for each pulse in the train
                    snippet_full = dat{i_ex}.(pTypes{i_tf}).(conds{i_cond})(snip_idx ,i_ch);
                    baseline = mean(snippet_full(1:prePulseSamps));
                    snippet_full = snippet_full - baseline;
                    
                    dat{i_ex}.(pTypes{i_tf}).snips.(conds{i_cond}){i_ch}(i_pulse,:) = snippet_full;
                    
                    
                end
            end
            
        end
        
    end
    
end


%% ANALYZE THE SNIPPETS AND CALCULATE VARIOUS STATS

% REMINDER: this is happening in separate cell-script because I will define
% the analysis region based off the average first pulse across TF conds.
% This requires grabbing all the snippets b/4 any analysis can proceed. The
% preceeding cell will need to be run before this one.

% should the analysis window be defined based off the average response
% following the first pulse (across conditions) or should the analysis
% window be unique for each pulse?
FIRSTPULSE = false;

for i_ex = 1:Nexpts
    
    conds = {'FV_Na', 'nbqx_apv_cd2_ttx', 'synapticTransmission'};
    for i_cond = 1:numel(conds)
        
        pTypes = fieldnames(dat{i_ex});
        Ntfs = numel(pTypes);
        
        sampRate = info{i_ex}.(pTypes{1}).(conds{i_cond}).sampRate;
        prePulseSamps = ceil(prePulseTime .* sampRate); % samples prior to pulse onset
        postPulseSamps = ceil(postPulseTime .* sampRate); % samples available after pulse ONSET
        photoDelay= 0.0006; % timeout following pulse offset
        
        for i_ch = 1:2;
            
            % deal with some exceptions
            % exception 1 (common)
            if ~info{i_ex}.ignoreChans(i_ch)
                continue
            end
            
            % exception 2 (uncommon)
            % a strange case where HS1 was set to Vclamp, and so
            % HS2's data got put in the first column...
            if strcmpi(info{i_ex}.mouse, 'CH_150105_D')
                if i_ch == 1; error('something went wrong'); end
                i_ch = 1;
            end
            
            
            %
            % For each 'condition', I need to identify keep points on the
            % raw waveforms (eg., trough time, peak time, etc.). These
            % time points will be estimated from the first pulse in the
            % train (averaged across TF conditions) or on a pulse-by-pulse
            % basis.
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%
            if FIRSTPULSE
                
                % calculate the analysis windows based off the average 1st pulse,
                % which should be the same across TFs (within pharmacology
                % condition and channel).
                firstPulse = nan(Ntfs, prePulseSamps+postPulseSamps+1);
                for i_tf = 1:Ntfs
                    firstPulse(i_tf,:) = dat{i_ex}.(pTypes{i_tf}).snips.(conds{i_cond}){i_ch}(1,:);
                end
                
                firstPulse = mean(firstPulse,1);
                pWidth = info{i_ex}.(pTypes{i_tf}).(conds{i_cond}).pWidth;
                tt = ([0:numel(firstPulse)-1] - prePulseSamps) ./ sampRate; % time=0 is when the LED comes on
                [troughidx, peakidx, takeoff]  = anlyMod_getWFepochs(firstPulse, tt, conds{i_cond}, pWidth, photoDelay);
                
            end
            

            % 
            % Now do the analysis on a pulse by pulse basis. Loop over TF
            % conditions, and pulses. If the important time points haven't
            % already been determined, do so right before the analysis.
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for i_tf = 1:Ntfs
                
                Npulses = sum(info{i_ex}.(pTypes{i_tf}).(conds{i_cond}).pulseOn_idx);
                pOnIdx = find(info{i_ex}.(pTypes{i_tf}).(conds{i_cond}).pulseOn_idx);
                
                for i_pulse = 1:Npulses
                    
                    snippet = dat{i_ex}.(pTypes{i_tf}).snips.(conds{i_cond}){i_ch}(i_pulse,:);
                    tt = ([0:numel(snippet)-1] - prePulseSamps) ./ sampRate; % time=0 is when the LED comes on
                    
                    if ~FIRSTPULSE
                        
                        pWidth = info{i_ex}.(pTypes{i_tf}).(conds{i_cond}).pWidth;
                        [troughidx, peakidx]  = anlyMod_getWFepochs(snippet, tt, conds{i_cond}, pWidth, photoDelay);
                        
                    end
                    
                    % store the peak and trough indicies for plotting later (if desired)
                    dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).trpk_inds{i_ch}(i_pulse,:) = [troughidx, peakidx];
                    
                    
                    % add a few points on either side of the true trough/peak
                    trough_window = troughidx-4: min([troughidx+4, numel(snippet)]);
                    peak_window = peakidx-4: min([peakidx+4, numel(snippet)]);
                    
                    %
                    % store some stats for each pulse
                    %
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%
                    switch conds{i_cond}
                        case 'FV_Na'
                            
                            trough = mean(snippet(trough_window));
                            peak = mean(snippet(peak_window));
                            
                            pk2tr = peak-trough;
                            dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).pk2tr{i_ch}(i_pulse) = pk2tr;
                            dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).diffval{i_ch}(i_pulse) = trough;
                            
                        case 'nbqx_apv_cd2_ttx'
                            
                            trough = mean(snippet(trough_window));
                            dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).diffval{i_ch}(i_pulse) = trough;
                            
                            % fit a single tau to the decay using OLS
                            startVal = snippet(troughidx) .* 0.80;
                            startIdx = find((snippet > startVal) & (tt > tt(troughidx)), 1, 'first');
                            fit_tt = tt(startIdx : end);
                            fit_dat = snippet(startIdx : end);
                            
                            % make sure none of the fit_dat points are
                            % positive because the fitting routine will
                            % assume that all the points are negative
                            critval = -0.005;
                            if any(fit_dat > critval)
                                l_pos = fit_dat > critval;
                                fit_tt = fit_tt(~l_pos);
                                fit_dat = fit_dat(~l_pos);
                            end

                            betas = [fit_tt(:), ones(size(fit_tt(:)))] \ log(abs(fit_dat(:)));
                            
%                             if strcmpi(info{i_ex}.mouse, 'CH_150119_C') && (info{i_ex}.stimSite == 1) && (i_pulse > 1) && (info{i_ex}.(pTypes{i_tf}).(conds{i_cond}).pTF > 70)
%                                 keyboard;
%                             end
                            
                            % store the slope params
                            dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).tau_m{i_ch}(i_pulse) = betas(1);
                            dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).tau_b{i_ch}(i_pulse) = betas(2);
                            dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).tau_ind{i_ch}(i_pulse) = startIdx;
                            
                            
                        case 'synapticTransmission'
                            
                            trough = mean(snippet(trough_window));
                            dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).diffval{i_ch}(i_pulse) = trough;
                            
                            %
                            % fit the slope using OLS regression
                            %%%%%%%%%%%%%%%%
                            synapseDelay = 0.0018;
                            startSlopeVal = snippet(troughidx) .* 0.20;
                            slopeStartIdx = find((tt >= synapseDelay) & (snippet <= startSlopeVal), 1, 'first');
                            stopSlopeVal = trough .* 0.80;
                            slopeStopIdx = find((tt < tt(troughidx)) & (snippet >= stopSlopeVal) , 1, 'last');
                            fit_dat = snippet(slopeStartIdx : slopeStopIdx);
                            fit_tt = tt(slopeStartIdx : slopeStopIdx);
                            
                            %check to make sure the waveform is not
                            %contaminated by two peaks
                            critval = -2 ./ sampRate; % empirically pretty good at catching non-monotonic trends
                            posslope = diff(fit_dat)> critval; 
                            consecPosSlope = conv(double(posslope), [1 1]);
                            if any(consecPosSlope >= 2);
                                posSlopeIdx = find(consecPosSlope == 2, 1, 'first')-2;
                                slopeStopIdx = slopeStartIdx + posSlopeIdx;
                                fit_dat = snippet(slopeStartIdx : slopeStopIdx);
                                fit_tt = tt(slopeStartIdx : slopeStopIdx);
                            end
                            
                            % do the fit
                            betas = [fit_tt(:), ones(numel(fit_tt), 1)] \ fit_dat(:);
                            
                            % store the slope params
                            dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).slope{i_ch}(i_pulse) = betas(1);
                            dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).slope_intercept{i_ch}(i_pulse) = betas(2);
                            dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).slope_inds{i_ch}(i_pulse,:) = [slopeStartIdx, slopeStopIdx];
                            
                    end
                    
                    
                    %
                    % calculate the integral of the LFP signal. This
                    % integral will get adjusted later to reflect the
                    % integral of a noisy signal with no response to the
                    % LED...
                    %%%%%%%%%%%%%%%%%%%%%%
                    pWidth = info{i_ex}.(pTypes{i_tf}).(conds{i_cond}).pWidth;
                    pulse_window = (tt >= (pWidth+photoDelay)) & (tt < 0.008);
                    snippet_pulse = snippet(pulse_window);
                    area = sum(abs(snippet_pulse)) ./ (numel(snippet_pulse) ./sampRate); % adjusts for differences in sample rate between experiments
                    dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).area{i_ch}(i_pulse) = area;
                    
                    
                end % pulses
                
                
                % now calculate the integral of part of the signal
                % preceeding the LED pulses and use this to correct the
                % integral calculated above. The idea is that the integral
                % can be positive even if there is no response to the LED.
                Nsamps = numel(snippet_pulse);
                fullsweep = dat{i_ex}.(pTypes{i_tf}).(conds{i_cond})(:,i_ch);
                firstPulseOn = find(info{i_ex}.(pTypes{i_tf}).(conds{i_cond}).pulseOn_idx, 1, 'first');
                maxstartidx = firstPulseOn-Nsamps-1;
                randStartIdx = unidrnd(maxstartidx, 1e4, 1);
                indexMtx = repmat(0:Nsamps-1, size(randStartIdx,1), 1);
                indexMtx = bsxfun(@plus, randStartIdx, indexMtx);
                shuffleBaseline = fullsweep(indexMtx);
                noiseIntegral = sum(abs(shuffleBaseline), 2) ./ (Nsamps ./ sampRate);
                noiseIntegral = mean(noiseIntegral);
                
                % correct the integral from above, based on the shuffle
                % corrected noiseIntegral.
                tmp = dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).area{i_ch};
                tmp = tmp - noiseIntegral;
                dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).area{i_ch} = tmp;
                
            end % tfs
        end % channels
        
    end % conditions
    
end % expts



%% MAKE SOME PLOTS: ONE SUMMARY PLOT PER RECORDING

close all

conds = {'nbqx_apv_cd2_ttx', 'FV_Na', 'synapticTransmission'};

CHECK_TRIAL_STATS = true;
RESTRICT_TO_STIM_SITE = true;
NORM_TO_PULSE1 = true;

for i_ex = 1:Nexpts
    
    for i_ch = 1:2;
        
        if RESTRICT_TO_STIM_SITE
            if i_ch ~= info{i_ex}.stimSite
                continue
            end
        end
        
        % deal with data and channels that need to be ignored.
        if ~info{i_ex}.ignoreChans(i_ch)
            continue
        else
            % a strange case where HS1 was set to Vclamp, and so
            % HS2's data got put in the first column...
            if strcmpi(info{i_ex}.mouse, 'CH_150105_D')
                if i_ch == 1; error('something went wrong'); end
                i_ch = 1;
            end
        end
        
        hFig = figure;
        set(gcf, 'position', [87 6 1260 799]);
        set(gcf, 'name', sprintf('%s, site %d, chan: %d', info{i_ex}.mouse, info{i_ex}.stimSite, i_ch))
        set(gcf, 'defaulttextinterpreter', 'none')
        s = warning('off', 'MATLAB:uitabgroup:OldVersion');
        hTabGroup = uitabgroup('Parent',hFig);
        
        
        hTabs(i_cond) = uitab('Parent', hTabGroup, 'Title', 'Raw Data');
        hAx(i_cond) = axes('Parent', hTabs(i_cond));
        hold on,
        for i_cond = 1:numel(conds)
            
            %
            % plot the raw (Average) traces
            %
            pTypes = fieldnames(dat{i_ex});
            Ntfs = numel(pTypes);
            Nconds = numel(conds);
            for i_tf = 1:Ntfs
                
                tmp_raw = dat{i_ex}.(pTypes{i_tf}).snips.(conds{i_cond}){i_ch}';
                inds = dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).trpk_inds{i_ch};
                Ntime = size(tmp_raw,1);
                Ncols = size(tmp_raw,2);
                inds_tt = inds;
                inds_idx = bsxfun(@plus, inds, Ntime .* (0:Ncols-1)');
                
                tt = (0:Ntime-1) ./ info{i_ex}.(pTypes{i_tf}).(conds{i_cond}).sampRate;
                tt = (tt - prePulseTime) .* 1000; % in ms.
                
                subplot(Nconds, Ntfs, i_tf+((i_cond-1) * Ntfs)), 
                cmap = colormap('copper');
                cidx = round(linspace(1, size(cmap,1), size(tmp_raw,2)));
                cmap = cmap(cidx,:);
                set(gca, 'colororder', cmap, 'NextPlot', 'replacechildren');
                
                plot(tt, tmp_raw, 'linewidth', 2), hold on,
                
                if CHECK_TRIAL_STATS
                    
                    % check the peak and trough indicies
                    plot(tt(inds(:,1)), tmp_raw(inds_idx(:,1)), 'ro', 'markerfacecolor', 'r')
                    if strcmpi(conds{i_cond}, 'FV_Na')
                        plot(tt(inds(:,2)), tmp_raw(inds_idx(:,2)), 'co', 'markerfacecolor', 'c')
                    end
                    
                    % check slope fitting for the field PSP
                    if strcmpi(conds{i_cond}, 'synapticTransmission')
                        for i_pulse = 1:size(tmp_raw,2);
                            startIdx = dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).slope_inds{i_ch}(i_pulse,1);
                            stopIdx = dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).slope_inds{i_ch}(i_pulse,2);
                            m = dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).slope{i_ch}(i_pulse);
                            m = m./1000; % slope was calculated in seconds, but ploting is in ms.
                            b = dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).slope_intercept{i_ch}(i_pulse);
                            fit_tt = tt(startIdx:stopIdx);
                            fit_vals = m .* fit_tt + b;
                            plot(fit_tt, fit_vals, 'g-')
                        end
                    end
                    
                    
                    % check best fitting decay tau for the opsin current
                    if strcmpi(conds{i_cond}, 'nbqx_apv_cd2_ttx')
                        for i_pulse = 1:size(tmp_raw,2);
                            startIdx = dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).tau_ind{i_ch}(i_pulse);
                            m = dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).tau_m{i_ch}(i_pulse);
                            m = m./1000; % slope was calculated in seconds, but ploting is in ms.
                            b = dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).tau_b{i_ch}(i_pulse);
                            fit_tt = tt(startIdx:end);
                            fit_vals = exp(m .* fit_tt) .* exp(b);
                            fit_vals = -fit_vals; % compensate for the fact that the fit was on abs(rawdata), but opsin current is negative
                            plot(fit_tt, fit_vals, 'm-')
                        end
                    end
                    
                end
                
                axis tight
                title(pTypes{i_tf})
                xlabel('time (ms)')
                if i_tf==1;
                    switch conds{i_cond}
                        case 'nbqx_apv_cd2_ttx'
                            ylabel(sprintf('%s Current', info{i_ex}.opsin))
                        case 'FV_Na'
                            ylabel('Fiber Volley')
                        case 'synapticTransmission'
                            ylabel('Field PSP')
                    end
                    
                end

            end
        end
        
        
        %
        % Now plot the summary statistics
        %
        %%%%%%%%%%%%%%%%%%
        
        
        hTabs(i_cond) = uitab('Parent', hTabGroup, 'Title', 'Summary Stats');
        hAx(i_cond) = axes('Parent', hTabs(i_cond));
        hold on,
        statTypes = {'diffval', 'area', 'pk2tr', 'slope', 'tau_m'};
        Nstats = numel(statTypes);
        for i_cond = 1:Nconds
            
            for i_stat = 1:Nstats
                
                if ~isfield(dat{i_ex}.(pTypes{1}).stats.(conds{i_cond}), statTypes{i_stat})
                    continue
                end
                
                subplot(Nconds, Nstats, i_stat+((i_cond-1)*Nstats)), hold on,
                legtext = {};
                rawvals={};
                tmptfs = [];
                for i_tf = 1:Ntfs
                    tmp = dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).(statTypes{i_stat}){i_ch};
                    if NORM_TO_PULSE1
                        tmp = tmp ./ tmp(1);
                    end
                    rawvals = cat(1, rawvals, tmp);
                    tmptfs = cat(1, tmptfs, info{i_ex}.(pTypes{i_tf}).(conds{i_cond}).pTF);
                end
                [~, order] = sort(tmptfs);
                
                for i_plt = 1:numel(order)
                    idx = order(i_plt);
                    plot(1:numel(rawvals{idx}), rawvals{idx}, 'o-', 'color', cmap(i_plt,:), 'linewidth', 2)
                    legtext = cat(2, legtext, num2str(tmptfs(idx)));
                end
                xlabel('Pulse number')
                xlim([1, max(cellfun(@numel, rawvals))])
                title(sprintf('%s', statTypes{i_stat}));
                if NORM_TO_PULSE1
                    yvals = get(gca, 'ylim');
                    yvals(1) = min([0, yvals(1)]);
                    set(gca, 'ylim', [0, yvals(2)]);
                end
                if all([i_cond, i_stat] == [2,1])
                    legend(legtext, 'location', 'best')
                    legend boxoff
                end
                if i_stat == 1
                    switch conds{i_cond}
                        case 'nbqx_apv_cd2_ttx'
                            ylabel(sprintf('%s current', info{i_ex}.opsin))
                        case 'FV_Na'
                            ylabel('Fiber Volley')
                        case 'synapticTransmission'
                            ylabel('field PSP')
                    end
                end
                
            end
            
            drawnow % force the uitab plot to update in quasi real time
            
        end
        
    end
end

%% POPULATION SUMMARY PLOTS


% plan: loop over opsins. Only consider a single recording channel (distal
% or proximal). Show PP ratio as a function of TF. Show Pn:P1 ratio as a
% function of pulse number for each frequency

STIMSITE = true;  % true => stimsite,  false => distal site


%initialize the outputs
opsinTypes = {'chr2', 'ochief'};
conds = {'nbqx_apv_cd2_ttx', 'FV_Na', 'synapticTransmission'};
statTypes = {'diffval', 'area', 'pk2tr', 'slope', 'tau_m'};
for i_opsin = 1:numel(opsinTypes)
    for i_cond = 1:numel(conds)
        for i_stat = 1:numel(statTypes)
            pop.(opsinTypes{i_opsin}).pnp1.(conds{i_cond}).(statTypes{i_stat}) = {[],{}}; % {{TF}, {Vals}}
            pop.(opsinTypes{i_opsin}).raw.(conds{i_cond}).(statTypes{i_stat}) = {[],{}}; % {{TF}, {Vals}}
        end
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  sort the data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i_ex = 1:numel(dat)
    
    
    % skip cases where the led was not targeted to either of the recording
    % sites
    if isnan(info{i_ex}.stimSite)
        continue
    end
    
    % Determine which recording channel to analyze
    if strcmpi(info{i_ex}.mouse, 'CH_150105_D') % a strange exception that bucks the rules.
        CHANNEL = 1;
    else % all other experiments...
        if STIMSITE
            CHANNEL = info{i_ex}.stimSite;
        else
            if info{i_ex}.stimSite == 1
                CHANNEL = 2;
            elseif info{i_ex}.stimSite == 2
                CHANNEL = 1;
            end
        end
    end
    
    opsin = lower(info{i_ex}.opsin);
    
    pTypes = fieldnames(dat{i_ex});
    Ntfs = numel(pTypes);
    for i_tf = 1:Ntfs
        
        for i_cond = 1:numel(conds);
            
            for i_stat = 1:numel(statTypes)
                
                % skip instances where the data do not exist
                if ~isfield(dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}), statTypes{i_stat})
                    continue
                end
                
                % structure the pnp1 data
                ex_stat = dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).(statTypes{i_stat}){CHANNEL};
                if strcmpi(statTypes{i_stat}, 'tau_m')
                    ex_stat = 1./ex_stat; % tau_m isn't a tau, need to convert to tau
                end
                
                % convert to paired pulse measures by normalizing by the
                % height of the first pulse.
                ex_p1p2 = ex_stat ./ ex_stat(1);
                
                % grab the data field of the 'pop' structure
                tmp_pop_p1p2 = pop.(opsin).pnp1.(conds{i_cond}).(statTypes{i_stat});
                tmp_pop_raw = pop.(opsin).raw.(conds{i_cond}).(statTypes{i_stat});
                tf = info{i_ex}.(pTypes{i_tf}).(conds{i_cond}).pTF;
                
                % is there already an entry for this TF and drug condition?
                alreadyThere = any(tmp_pop_p1p2{1} == tf);
                if alreadyThere
                    idx = tmp_pop_p1p2{1} == tf;
                    
                    % deal with cases where there are different numbers of
                    % pulses
                    nPulsesExisting = size(tmp_pop_p1p2{2}{idx},2);
                    if numel(ex_p1p2) < nPulsesExisting
                        ex_p1p2(1,end+1:nPulsesExisting) = nan;
                        ex_stat(1,end+1:nPulsesExisting) = nan;
                    else
                        tmp_pop_p1p2{2}{idx}(:,end+1:numel(ex_p1p2)) = nan;
                        tmp_pop_raw{2}{idx}(:,end+1:numel(ex_p1p2)) = nan;
                    end
                    
                    tmp_pop_p1p2{2}{idx} = cat(1, tmp_pop_p1p2{2}{idx}, ex_p1p2);
                    tmp_pop_raw{2}{idx} = cat(1, tmp_pop_raw{2}{idx}, ex_stat);
                    
                else
                    tmp_pop_p1p2{1} = cat(1, tmp_pop_p1p2{1}, tf);
                    tmp_pop_p1p2{2} = cat(1, tmp_pop_p1p2{2}, ex_p1p2);
                    
                    tmp_pop_raw{1} = cat(1, tmp_pop_raw{1}, tf);
                    tmp_pop_raw{2} = cat(1, tmp_pop_raw{2}, ex_stat);
                end
                
                % replace the 'pop' structure version with the tmp version
                pop.(opsin).pnp1.(conds{i_cond}).(statTypes{i_stat}) = tmp_pop_p1p2;
                pop.(opsin).raw.(conds{i_cond}).(statTypes{i_stat}) = tmp_pop_raw;
                
            end
        end
    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  plotting routines
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i_opsin = 1:numel(opsinTypes)
    
    figure
    set(gcf, 'name', opsinTypes{i_opsin}, 'defaulttextinterpreter', 'none')
    set(gcf, 'position', [109    31   891   754])
    
    Nconds = numel(conds);
    for i_cond = 1:Nconds;
        
        Nstats = numel(statTypes);
        for i_stat = 1:Nstats
            
            % grab the data
            tmp_dat = pop.(opsinTypes{i_opsin}).pnp1.(conds{i_cond}).(statTypes{i_stat});
            
            % make sure there's actually data there
            if isempty(tmp_dat{1}); continue; end % some things are not in the data set
            
            % now do the plotting
            subplot(Nconds, Nstats, (i_cond-1).*Nstats + i_stat)
            title(statTypes{i_stat})
            
            tfs = tmp_dat{1};
            cmap = colormap('copper');
            cidx = round(linspace(1, size(cmap,1), 6));
            cmap = cmap(cidx,:);
            set(gca, 'colororder', cmap, 'NextPlot', 'replacechildren');
            hold on
            
            legtext = {};
            [~, tforder] = sort(tmp_dat{1}, 1, 'ascend');
            for i_tf = 1:numel(tfs);
                
                tf_idx = tforder(i_tf);
                
                xbar = nanmean(tmp_dat{2}{tf_idx}, 1);
                sem = nanstd(tmp_dat{2}{tf_idx}, [], 1) ./ sqrt(sum(~isnan(tmp_dat{2}{tf_idx}), 1));
                
                %errorbar(1:numel(xbar), xbar, sem, 'color', cmap(i_tf,:), 'linewidth', 2)
                errorbar(1:7, xbar(1:7), sem(1:7), 'color', cmap(i_tf,:), 'linewidth', 2) % only plot the first 7 pulses
                
                legtext = cat(2, legtext, num2str(tfs(tf_idx)));
            end
            axis tight
            xlabel('Pulse number')
            if i_stat==1
                ylabel(sprintf('%s \n Pn:P1 ratio', conds{i_cond}))
            end
            if i_stat==1  && i_cond==2
                legend(legtext, 'location', 'best')
                legend boxoff
            end
            yvals = get(gca, 'ylim');
            yvals(1) = min([0, yvals(1)]);
            set(gca, 'ylim', yvals);
            if yvals(1)<0
                plot([1,numel(xbar)], [0,0] , 'k--', 'linewidth', 2)
            end
            
        end
    end
end



%% OPSIN CURRENT VS. FIBER VOLLEY FOR CHIEF AND CHR2


FV_STAT = 'pk2tr';
OPSIN_STAT = 'area';
STAT_TYPE = 'raw';
TFs = [10,20,40, 60, 100];


cmap = colormap('copper');
cidx = round(linspace(1, size(cmap,1), 6));
cmap = cmap(cidx,:);
figure, hold on,

for i_tf = 1:numel(TFs)
   
    % pull out the data for ChR2:
    tf_idx = pop.chr2.(STAT_TYPE).FV_Na.(FV_STAT){1} == TFs(i_tf);
    if sum(tf_idx == 1);
        chr2_fv_raw = pop.chr2.(STAT_TYPE).FV_Na.(FV_STAT){2}{tf_idx};
        chr2_opsin_raw = pop.chr2.(STAT_TYPE).nbqx_apv_cd2_ttx.(OPSIN_STAT){2}{tf_idx};
    else
        chr2_fv_raw = nan(5,7); % a hack to allow plotting of tf conds for oChIEF that don't exist for ChR2
        chr2_opsin_raw = nan(5,7);
    end
    
    % pull out the data for oChIEF:
    tf_idx = pop.ochief.(STAT_TYPE).FV_Na.(FV_STAT){1} == TFs(i_tf);
    ochief_fv_raw = pop.ochief.(STAT_TYPE).FV_Na.(FV_STAT){2}{tf_idx};
    ochief_opsin_raw = pop.ochief.(STAT_TYPE).nbqx_apv_cd2_ttx.(OPSIN_STAT){2}{tf_idx};
    
    
    for i_pulse = 1:6
        plot(chr2_opsin_raw(:,i_pulse), chr2_fv_raw(:,i_pulse), '+', 'markeredgecolor', cmap(i_pulse,:))
        plot(ochief_opsin_raw(:,i_pulse), ochief_fv_raw(:,i_pulse), 'o', 'markeredgecolor', cmap(i_pulse,:))
    end
    plot(mean(chr2_opsin_raw(:,1:6),1), mean(chr2_fv_raw(:,1:6),1), '-sb', 'markerfacecolor', 'b')
    plot(mean(ochief_opsin_raw(:,1:6),1), mean(ochief_fv_raw(:,1:6),1), '-sr', 'markerfacecolor', 'r')
    
    
    
    
end






%%  SHAPE OF THE FIRST PULSE RESPONSE FOR oChIEF AND ChR2


STIMSITE = true;

conds = {'nbqx_apv_cd2_ttx', 'FV_Na'};
chr2_examp = {[] []};
ochief_examp = {[] []};
for i_ex = 1:Nexpts
    
    
    % skip cases where the led was not targeted to either of the recording
    % sites
    if isnan(info{i_ex}.stimSite)
        continue
    end
    
    % Determine which recording channel to analyze
    if strcmpi(info{i_ex}.mouse, 'CH_150105_D') % a strange exception that bucks the rules.
        CHANNEL = 1;
    else % all other experiments...
        if STIMSITE
            CHANNEL = info{i_ex}.stimSite;
        else
            if info{i_ex}.stimSite == 1
                CHANNEL = 2;
            elseif info{i_ex}.stimSite == 2
                CHANNEL = 1;
            end
        end
    end
    
     
    for i_cond = 1:numel(conds)
        
        pTypes = fieldnames(dat{i_ex});
        Ntfs = numel(pTypes);
        
        sampRate = info{i_ex}.(pTypes{1}).(conds{i_cond}).sampRate;
        prePulseSamps = ceil(prePulseTime .* sampRate); % samples prior to pulse onset
        postPulseSamps = ceil(postPulseTime .* sampRate); % samples available after pulse ONSET
        
        
        firstPulse = nan(Ntfs, prePulseSamps+postPulseSamps+1);
        for i_tf = 1:Ntfs
            
            tmp_dat = dat{i_ex}.(pTypes{i_tf}).snips.(conds{i_cond}){CHANNEL}(1,:);
            switch conds{i_cond}
                case 'FV_Na'
                    normFact = dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).pk2tr{CHANNEL}(1);
                case 'nbqx_apv_cd2_ttx'
                    normFact = dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).diffval{CHANNEL}(1);
            end
            
            firstPulse(i_tf,:) = tmp_dat ./ normFact;
            
        end
        
        % take the average, and normalize
        firstPulse = mean(firstPulse,1);
        
        
        % store the first pulse according to the opsin type, but only
        % if the sampling rate = 20kHz (so that the sampling latice is
        % the same...)
        correctSR = info{i_ex}.(pTypes{i_tf}).(conds{i_cond}).sampRate == 20e3;
        if ~correctSR
            
            oldSampRate = info{i_ex}.(pTypes{i_tf}).(conds{i_cond}).sampRate;
            old_tt = [0 : numel(firstPulse)-1] ./ oldSampRate;
            
            totalTime = numel(firstPulse) ./ oldSampRate;
            newSampRate = 20e3;
            newNumSamps = ceil(totalTime .* newSampRate);
            new_tt = [0 : newNumSamps-1] ./ newSampRate;
            firstPulse = interp1(old_tt, firstPulse, new_tt);
            
        end
        switch lower(info{i_ex}.opsin)
                case 'ochief'
                    ochief_examp{i_cond} = cat(1, ochief_examp{i_cond}, firstPulse);
                case 'chr2'
                    chr2_examp{i_cond} = cat(1, chr2_examp{i_cond}, firstPulse);
        end
        
    end % conditions
    
end % expts


for i_cond = 1:numel(conds)
    
    % the example first 
    tt = (0:size(ochief_examp{i_cond},2)-1)./20e3;
    tt = tt-tt(prePulseSamps);
    tt = tt.*1000;
    
    figure, hold on,
    plot(tt, ochief_examp{i_cond}', 'm');
    plot(tt, chr2_examp{i_cond}', 'c');
    plot(tt, mean(ochief_examp{i_cond}, 1), 'r', 'linewidth', 4);
    plot(tt, mean(chr2_examp{i_cond}, 1), 'b', 'linewidth', 4);
    xlabel('time (msec)')
    ylabel('Normalized LFP amplitude')
    axis tight
    switch conds{i_cond}
        case 'FV_Na'
            title('Average FV for first pulse')
        case 'nbqx_apv_cd2_ttx'
            title('Average opsin current for first pulse')
    end
    
    
    
    if strcmpi(conds{i_cond}, 'FV_Na')
        figure, hold on,
        plot(tt, cumsum(-mean(ochief_examp{i_cond}, 1)), 'r', 'linewidth', 4);
        plot(tt, cumsum(-mean(chr2_examp{i_cond}, 1)), 'b', 'linewidth', 4);
        xlabel('time (msec)')
        ylabel('Normalized LFP amplitude')
        title('Integral of (negative) FV pulse')
        axis tight
    end
    
end




