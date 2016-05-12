%% WHICH MICE SHOULD CONTRIBUTE?  [Main project]

% clear out the workspace
fin


% decide what experiment to run
EXPTTYPE = 1;
BRAINAREA = 'al';
COMBINE_CHIEF = false;
switch EXPTTYPE
    case 1
        EXPTTYPE = 'main expt';
    case 2
        EXPTTYPE = 'interleaved amps';
    case 3
        EXPTTYPE = 'rundown';
    case 4
        EXPTTYPE = 'stim positions';
    case 5
        EXPTTYPE = 'L5 record';
    case 6
        EXPTTYPE = '4-AP';
    case 7
        EXPTTYPE = 'Baclofen';
    case 8
        EXPTTYPE = 'Intracellular';
    case 9
        EXPTTYPE = 'E-stim';
    case 10
        EXPTTYPE = 'all manuscript';
end

% grab the mouse names and sites from the excel workbook.
fname = [GL_DOCUPATH, 'Other_workbooks', filesep, 'fiberVolleyCellList.xlsx'];
[~,~,wb_expt] = xlsread(fname, 2);

wb_expt_nameidx = strcmpi(wb_expt(1,:), 'mouse name');
wb_expt_siteidx = strcmpi(wb_expt(1,:), 'site');
if strcmpi(EXPTTYPE, 'all manuscript')
    exptlistidx = strcmpi(wb_expt(1,:), 'main expt') | strcmpi(wb_expt(1,:), 'interleaved amps');% | strcmpi(wb_expt(1,:), 'Intracellular') | strcmpi(wb_expt(1,:), 'E-stim');
else
    exptlistidx = strcmpi(wb_expt(1,:), EXPTTYPE);
end

wb_expt = wb_expt(2:end, :); % notice that I'm hacking off the header row

% figure out the appropriate expts to analyze
l_expt = cellfun(@(x) isnumeric(x) && x==1, wb_expt(:, exptlistidx));
l_expt = sum(l_expt,2) > 0;
MouseName = wb_expt(l_expt, wb_expt_nameidx);
Site = wb_expt(l_expt, wb_expt_siteidx);

if ~strcmpi(BRAINAREA, 'any')
    % open up the other workbook and extract brain area information
    [~,~, wb_info] = xlsread(fname, 1);
    wb_info_areaidx = find(strcmpi(wb_info(1,:), 'brain area')==1, 1, 'first');
    wb_info_nameidx = find(strcmpi(wb_info(1,:), 'mouse name')==1, 1, 'first');
    wb_info_siteidx = find(strcmpi(wb_info(1,:), 'site') == 1, 1, 'first');
    wb_info = wb_info(2:end,:);
    
    l_correctArea = false(numel(MouseName), 1);
    for i_ex = 1:numel(MouseName)
        l_mouse = regexpi(wb_info(:, wb_info_nameidx), MouseName{i_ex});
        l_mouse = cellfun(@(x) ~isempty(x), l_mouse);
        l_site = cell2mat(wb_info(:, wb_info_siteidx)) == Site{i_ex};
        l_union = l_mouse & l_site;
        
        areamatch = regexpi(wb_info(l_union, wb_info_areaidx), BRAINAREA);
        areamatch = all(cellfun(@(x) ~isempty(x), areamatch));
        l_correctArea(i_ex) = areamatch;
    end
    
    MouseName = MouseName(l_correctArea);
    Site = Site(l_correctArea);
end


in = [MouseName, Site]


%% EXTRACT THE RAW DATA FROM EACH DATA FILE

RMLINENOISE = false;

% grab the fiber volley pop excel workbook
fname = [GL_DOCUPATH, 'Other_workbooks', filesep, 'fiberVolleyCellList.xlsx'];
[~,txt, raw] = xlsread(fname);
raw(size(txt,1)+1:end, :) = [];
raw(:,size(txt,2)+1:end) = [];
channelIdx = cellfun(@(x) ~isempty(x), regexpi(raw(1,:), 'CH\d'));
opsinIdx = strcmpi(raw(1,:), 'opsin');
primaryStimSiteIdx = strcmpi(raw(1,:), 'Primary Stim Site (0,0)');
stimSiteCordIdx = strcmpi(raw(1,:), 'optostim target');
clear txt

% replace the 'file names' with fully qualified paths. I'm hoping that this
% allows for parallel operations
Nfiles = size(raw,1)-1;
idx_mousename = strcmpi(raw(1,:), 'Mouse Name');
idx_fname =  strcmpi(raw(1,:), 'file name');
tmp_fnames = raw(2:end, idx_fname);
tmp_mousenames = raw(2:end, idx_mousename);
prefix = cellfun(@(x,y) [x,y], repmat({GL_DATPATH}, Nfiles, 1), tmp_mousenames, 'uniformoutput', false);
prefix = cellfun(@(x,y) [x,y], prefix, repmat({[filesep, 'Physiology', filesep]}, Nfiles, 1), 'uniformoutput', false);
tmp_fnames = cellfun(@(x,y) [x,y], prefix, tmp_fnames, 'uniformoutput', false);
tmp_fnames = cellfun(@(x,y) [x,y], tmp_fnames, repmat({'.abf'}, Nfiles, 1), 'uniformoutput', false);
raw(2:end, idx_fname) = tmp_fnames;

% combine the chief variants if desired
if COMBINE_CHIEF
    newopsin = raw(:, opsinIdx);
    newopsin = cellfun(@(x) regexprep(x, 'chief_cit', 'chief_all'), newopsin, 'uniformoutput', false);
    newopsin = cellfun(@(x) regexprep(x, 'chief_flx', 'chief_all'), newopsin, 'uniformoutput', false);
    raw(:, opsinIdx) = newopsin;
end


% do the analysis
Nexpts = size(in,1);
dat = {};

pool = gcp('nocreate');
if isempty(pool)
    pool = parpool(min([32, Nexpts]));
end


parfor i_ex = 1:Nexpts;
    
    % figure out what rows in the work book to pay attention to
    l_mouse = cellfun(@(x) ~isempty(x), regexp(raw(:,1), in{i_ex,1}));
    l_site = [false ; cell2mat(raw(2:end,2)) == in{i_ex,2}]; % add a leading 'false' to account for the header row in 'raw'
    l_expt = l_mouse & l_site;
    
    % run the analysis
    [dat{i_ex}, info{i_ex}] = fiberVolleyAnalysis(l_expt, raw, false, RMLINENOISE, EXPTTYPE);
    
    % enter a few other useful things into the structure
    info{i_ex}.mouse =  in{i_ex,1};
    info{i_ex}.opsin = unique(cell2mat(raw(l_expt, opsinIdx)), 'rows');
    info{i_ex}.ignoreChans = unique(cell2mat(raw(l_expt, channelIdx)), 'rows');
    info{i_ex}.stimSite = unique(cell2mat(raw(l_expt, primaryStimSiteIdx)), 'rows');
    if ischar(info{i_ex}.stimSite); % deals with nans
        info{i_ex}.stimSite = str2num(info{i_ex}.stimSite);
    end
    tmp = raw(l_expt, stimSiteCordIdx);
    if ischar(tmp{1})
        tmp = unique(tmp);
        info{i_ex}.optStimCords  = str2num(tmp{1});
    else
        tmp = cell2mat(tmp);
        info{i_ex}.optStimCords  = unique(tmp);
    end
end

disp('All done importing data')





 %% PULL OUT SNIPPETS OF DATA FOR EACH PULSE (ANALYZE THEM LATER)

prePulseTime = 0.001; % in sec
postPulseTime = 0.010; % in sec
Nexpts = size(in,1);

for i_ex = 1:Nexpts
    
    pTypes = fieldnames(dat{i_ex});
    Ntfs = numel(pTypes);
    for i_tf = 1:Ntfs
        
        conds = fieldnames(dat{i_ex}.(pTypes{i_tf}));
        for i_cond = 1:numel(conds)
            
            % check to see if this condition exists
            if ~isfield(info{i_ex}.(pTypes{i_tf}), conds{i_cond})
                continue
            end
            
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
                    else
                        % occasionally CH1 is disabled, and CH2 occupies the
                        % 1st and only column, adjust i_ch to account for these
                        % cases
                        Nchannels = sum(info{i_ex}.ignoreChans);
                        if (Nchannels == 1) && (i_ch == 2)
                            i_ch = 1;
                        end
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

clc

% REMINDER: this is happening in separate cell-script because I will define
% the analysis region based off the average first pulse across TF conds.
% This requires grabbing all the snippets b/4 any analysis can proceed. The
% preceeding cell will need to be run before this one.

% should the analysis window be defined based off the average response
% following the first pulse (across conditions) or should the analysis
% window be unique for each pulse?
FIRSTPULSE = false;

Nexpts = size(in,1);
parfor i_ex = 1:Nexpts
    
    % Determine the pharmacology conditions that are present. Look
    % specifically for a fiber volley, a TTX, and a synapticTransmission
    % condition
    pTypes = fieldnames(dat{i_ex});
    Ntfs = numel(pTypes);
    fldnames = {};
    for i_tf = 1:Ntfs
        fldnames = cat(1, fldnames, fieldnames(info{i_ex}.(pTypes{i_tf})));
    end
    fldnames = unique(fldnames);
    tags = {'FV_', 'ttx', 'synapticTransmission'}; %default line
    
    conds = {};
    for i_t = 1:numel(tags)
        idx = cellfun(@(x) ~isempty(x), regexpi(fldnames, tags{i_t}));
        conds = cat(2, conds, fldnames(idx)');
    end

    for i_cond = 1:numel(conds)  %#ok<*PFTUS>
        

        tmpnames = fieldnames(info{i_ex}.(pTypes{1}));
        sampRate = info{i_ex}.(pTypes{1}).(tmpnames{1}).sampRate; % assuming sampRates are the same across experments within a day
        prePulseSamps = ceil(prePulseTime .* sampRate); % % samples prior to pulse onset
        postPulseSamps = ceil(postPulseTime .* sampRate); % samples available after pulse ONSET
        if strcmpi(conds{i_cond}, 'synapticTransmission')
            photoDelay= 1e-3;
        elseif strcmpi(EXPTTYPE, 'E-stim')
            photoDelay = 300e-6;
        else
            photoDelay= 200e-6; % timeout following pulse offset (in sec)
        end
        anlyWindowInPts = ceil(450e-6 .* sampRate);
        if rem(anlyWindowInPts, 2)==0;
            anlyWindowInPts = anlyWindowInPts+1;
        end
        
        for i_ch = 1:2;
            
            % deal with data and channels that need to be ignored.
            if ~info{i_ex}.ignoreChans(i_ch)
                continue
            else
                % occasionally CH1 is disabled, and CH2 ocupys the
                % 1st and only column, adjust i_ch to account for these
                % cases
                Nchannels = sum(info{i_ex}.ignoreChans);
                if (Nchannels == 1) && (i_ch == 2)
                    i_ch = 1;
                end
            end
            
            
            %
            % For each 'condition', I need to identify key points on the
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
                switch conds{i_cond}
                    case {'FV_Na', 'FV_Na_Ca2_mGluR'}
                        troughTime = 0.00225;
                    otherwise
                        troughTime = 0.005;
                end
                trough_window = (tt >= pWidth+photoDelay) & (tt <= troughTime);
                if mean(firstPulse(trough_window))> 0;
                    direction = 'outward';
                else
                    direction = 'inward';
                end
                [troughidx, peakidx, takeoff]  = anlyMod_getWFepochs(firstPulse, tt, conds{i_cond}, pWidth, photoDelay, direction);
                
            end
            
            
            %
            % Now do the analysis on a pulse by pulse basis. Loop over TF
            % conditions, and pulses. If the important time points haven't
            % already been determined, do so right before the analysis.
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for i_tf = 1:Ntfs   %
                
                % check to see if this condition exists
                if ~isfield(info{i_ex}.(pTypes{i_tf}), conds{i_cond})
                    continue
                end
                
                Npulses = sum(info{i_ex}.(pTypes{i_tf}).(conds{i_cond}).pulseOn_idx);
                pOnIdx = find(info{i_ex}.(pTypes{i_tf}).(conds{i_cond}).pulseOn_idx);
                pWidth = info{i_ex}.(pTypes{i_tf}).(conds{i_cond}).pWidth;
                
                % figure out the direction of the responses
                allSnippets = dat{i_ex}.(pTypes{i_tf}).snips.(conds{i_cond}){i_ch};
                meanSnippet = mean(allSnippets, 1);
                Ntime = numel(meanSnippet);
                tt = ([0:Ntime-1] - prePulseSamps) ./ sampRate;   % time=0 is when the LED comes on
                switch conds{i_cond}
                    case {'FV_Na', 'FV_Na_Ca2_mGluR'}
                        troughTime = 0.00225;
                    otherwise
                        troughTime = 0.005;
                end
                trough_window = (tt >= pWidth+photoDelay) & (tt <= troughTime);
                if mean(meanSnippet(trough_window))> 0;
                    direction = 'outward';
                else
                    direction = 'inward';
                end
                
                
                for i_pulse = 1:Npulses
                    
                    snippet = dat{i_ex}.(pTypes{i_tf}).snips.(conds{i_cond}){i_ch}(i_pulse,:);
                    
                    if ~FIRSTPULSE
                        [troughidx, peakidx]  = anlyMod_getWFepochs(snippet, tt, conds{i_cond}, pWidth, photoDelay, direction, info{i_ex}.opsin);
                    end
                    
                    % store the peak and trough indicies for plotting later (if desired)
                    dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).trpk_inds{i_ch}(i_pulse,:) = [troughidx, peakidx];
                    
                    
                    % add a few points on either side of the true trough/peak
                    halfWindowInPts = round((anlyWindowInPts-1) / 2);% subtract one b/c centering will add one later, round so that it's an integer
                    trough_window = troughidx-halfWindowInPts: min([troughidx+halfWindowInPts, numel(snippet)]);
                    assert(~isempty(trough_window) && numel(trough_window)==anlyWindowInPts, 'ERROR: wrong number of points in trough window')
                    if any(strcmpi(conds{i_cond}, {'FV_Na', 'FV_Na_Ca2_mGluR'}))
                        peak_window = peakidx-halfWindowInPts: min([peakidx+halfWindowInPts, numel(snippet)]);
                        assert(~isempty(peak_window) && numel(peak_window)==anlyWindowInPts, 'ERROR: no data for peak windwo')
                    end
                    
                    
                    %
                    % store some stats for each pulse
                    %
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if any(strcmpi(conds{i_cond}, {'FV_Na', 'FV_Na_Ca2_mGluR'}))
                        
                        trough = mean(snippet(trough_window));
                        peak = mean(snippet(peak_window));
                        
                        pk2tr = peak-trough;
                        dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).pk2tr{i_ch}(i_pulse) = pk2tr;
                        dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).diffval{i_ch}(i_pulse) = trough;
                        
                        latency = tt(troughidx);
                        dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).latency{i_ch}(i_pulse) = latency;
                        
                    elseif strcmpi(conds{i_cond}, 'nbqx_apv_cd2')
                        
                        trough = mean(snippet(trough_window));
                        dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).diffval{i_ch}(i_pulse) = trough;
                        
                        % deal with cases where the opsin current is
                        % outward
                        tmp_snippet = snippet;
                        if strcmpi(direction, 'outward');
                            tmp_snippet = -tmp_snippet;
                        end
                        
                    elseif any(regexpi(conds{i_cond}, 'ttx')) && ~strcmpi(info{i_ex}.opsin, 'estim')
                        
                        trough = mean(snippet(trough_window));
                        dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).diffval{i_ch}(i_pulse) = trough;
                        
                        % deal with cases where the opsin current is
                        % outward
                        tmp_snippet = snippet;
                        if strcmpi(direction, 'outward');
                            tmp_snippet = -tmp_snippet;
                        end
                        
                        
                        % fit a single tau to the decay using OLS
                        startVal = tmp_snippet(troughidx) .* 0.85;
                        startIdx = find((tmp_snippet > startVal) & (tt > tt(troughidx)), 1, 'first');
                        fit_tt = tt(startIdx : end);
                        fit_dat = tmp_snippet(startIdx : end);
                       
                        dataforfit = ~isempty(fit_dat);
                        if dataforfit
                            fit_dat = -fit_dat; % so that it's a decaying exponential
                            
                            N_2ms = round(sampRate .* 0.002);
                            pred = [ones(N_2ms,1), fit_tt(1:N_2ms)'];
                            resp = fit_dat(1:N_2ms)';
                            guess = pred \ resp;
                            guess(1) = exp(guess(1)); % to solve for amplitude of exponential
                            fitopt = fitoptions('exp1');
                            fitopt.StartPoint = guess;
                            fitopt.Lower = [0, -inf]; % constrain amps to be +, Taus to be -.
                            fitopt.Upper = [inf, 0];
                            fitopt.TolFun = 1e-10;
                            fitopt.TolX = 1e-10;
                            bobj = fit(fit_tt(:), fit_dat(:), 'exp1'); % fit model: YY = beta0 .* exp(beta1 * XX)
                        else
                            betas = [nan, nan];
                            startIdx = nan;
                        end
                        
                        % store the slope params
                        dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).tau_yint{i_ch}(i_pulse) = bobj.a;
                        dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).tau_invtc{i_ch}(i_pulse) = bobj.b;
                        dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).tau_ind{i_ch}(i_pulse) = startIdx;
                        
                        % store the latency to the trough
                        latency = tt(troughidx);
                        dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).latency{i_ch}(i_pulse) = latency;
                        
                    elseif any(regexpi(conds{i_cond}, 'ttx')) && strcmpi(info{i_ex}.opsin, 'estim')
                        
                        trough = mean(snippet(trough_window));
                        peak = mean(snippet(peak_window));
                        
                        pk2tr = peak-trough;
                        dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).pk2tr{i_ch}(i_pulse) = pk2tr;
                        dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).diffval{i_ch}(i_pulse) = trough;
                        
                        % store fake slope params
                        dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).tau_yint{i_ch}(i_pulse) = nan;
                        dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).tau_invtc{i_ch}(i_pulse) = nan;
                        dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).tau_ind{i_ch}(i_pulse) = nan;
                        
                        % store latency to the trough
                        latency = tt(troughidx);
                        dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).latency{i_ch}(i_pulse) = latency;
                        
                    elseif strcmpi(conds{i_cond},  'synapticTransmission')
                        
                        trough = mean(snippet(trough_window));
                        dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).diffval{i_ch}(i_pulse) = trough;
                        
                        %
                        % fit the slope using OLS regression
                        %%%%%%%%%%%%%%
                        tmp_snippet = snippet;
                        stopSlopeVal = trough .* 0.80;
                        if strcmpi(direction, 'outward');  % a hack to deal with inwards and outwards fEPSPs
                            tmp_snippet = -tmp_snippet;
                            stopSlopeVal = -trough .* 0.80;
                        end
                        synapseDelay = 0.001;
                        slopeStopIdx = find((tt < tt(troughidx)) & (tmp_snippet >= stopSlopeVal) , 1, 'last');
                        startSlopeVal = tmp_snippet(troughidx) .* 0.12;
                        slopeStartIdx = (tmp_snippet >= startSlopeVal) & (tt < tt(slopeStopIdx));
                        slopeStartIdx = find(slopeStartIdx, 1, 'last');
                        slopeStartIdx = max([slopeStartIdx, find(tt>synapseDelay, 1, 'first')]); % make sure that you don't encroach into the synaptic delay timeout
                        
                        % do the referencing with tmp_snippet, but grab
                        % the actual data from snippet
                        fit_dat = snippet(slopeStartIdx : slopeStopIdx);
                        fit_tt = tt(slopeStartIdx : slopeStopIdx);
                        
                        % find the max slope and center the analysis
                        % region around that
                        slope = [nan, diff(fit_dat)];
                        if strcmpi(direction, 'outward')
                            [~, idx] = max(slope);
                        else
                            [~, idx] = min(slope);
                        end
                        
                        idx = idx-halfWindowInPts : idx+halfWindowInPts;
                        if numel(fit_dat) <= anlyWindowInPts
                            idx = 1:numel(fit_dat);
                        elseif idx(1)<=0
                            idx = 1:anlyWindowInPts;
                        elseif idx(end)>numel(fit_dat)
                            idx = numel(fit_dat)- (anlyWindowInPts-1) : numel(fit_dat);
                        end
                        fit_tt = fit_tt(idx);
                        fit_dat = fit_dat(idx);
                        
                        if numel(fit_dat) > halfWindowInPts
                            % do the fit
                            betas = [fit_tt(:), ones(numel(fit_tt), 1)] \ fit_dat(:);
                            
                        else
                            warning('no data for slope, betas set to nan')
                            betas = [nan, nan];
                        end
                        
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
                    dt = 1./sampRate;
                    area = sum(abs(snippet_pulse)) .* dt; % adjusts for differences in sample rate between experiments
                    dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).area{i_ch}(i_pulse) = area;
                    areaStartIdx = find(tt >= (pWidth+photoDelay), 1, 'first');
                    areaStopIdx = find(tt >= 0.008, 1, 'first');
                    dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).area_inds{i_ch}(i_pulse,:) = [areaStartIdx, areaStopIdx];
                    
                    
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
                dt = 1./sampRate;
                noiseIntegral = sum(abs(shuffleBaseline), 2) .* dt;
                noiseIntegral = mean(noiseIntegral);
                
                % correct the integral from above, based on the shuffle
                % corrected noiseIntegral.
                tmp = dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).area{i_ch};
                tmp = tmp - noiseIntegral;
                dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).area{i_ch} = tmp;
                
                
                % derive a standard deviation prior to stimulus onset. Use
                % this later to do a sanity check on the PPRs for FVs
                numPtsInWindow = round(sampRate .* 0.010); % look at 10ms before stim onset
                startIdx = firstPulseOn - numPtsInWindow - 10;
                sigma = std(fullsweep(startIdx : (firstPulseOn-10)));
                dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).bkgnd_sigma{i_ch} = sigma;
                
                
                
            end % tfs
        end % channels
        
    end % conditions
    
end % expts



%% MAKE SOME PLOTS: ONE SUMMARY PLOT PER RECORDING

close all

CHECK_TRIAL_STATS = true;
RESTRICT_TO_STIM_SITE = true;
NORM_TO_PULSE1 = true;


% define the conditions that will get plotted
switch EXPTTYPE
    case '4-AP'
        conds = {'ttx_cd2', 'ttx_cd2_4AP', 'potassium'};
    case 'Baclofen'
        conds = {'ttx_cd2', 'ttx_cd2_bac10'};
    case 'Intracellular'
        conds = {'nbqx_apv_ttx'};
    case 'E-stim'
        conds = {'nbqx_apv_cd2', 'nbqx_apv_cd2_ttx', 'FV_Na'};
        
    otherwise
        
        conds = {'nbqx_apv_cd2_ttx', 'FV_Na', 'synapticTransmission'};
        %conds = {'FV_Na', 'nbqx_apv_cd2_ttx', 'nbqx_apv_cd2'}; % last two pharm steps with subtraction
        %conds = {'FV_Na', 'nbqx_apv_cd2_ttx', 'synapticTransmission'};
        %conds = {'FV_Na', 'FV_Na_Ca2_mGluR'}; % both %FVs
        %conds = {'FV_Na_Ca2_mGluR', 'nbqx_apv_ttx',  'synapticTransmission'}; % no cadmium
        %conds = {'nbqx_apv_cd2_ttx', 'nbqx_apv_ttx',  'synapticTransmission'}; % both opsin current verisons
end


Nexpts = size(in,1);
for i_ex = 1:Nexpts
    
    for i_ch = 1:2;
        
        if RESTRICT_TO_STIM_SITE && (i_ch ~= info{i_ex}.stimSite)
            continue
        end
        
        % deal with data and channels that need to be ignored.
        if ~info{i_ex}.ignoreChans(i_ch)
            continue
        else
            % occasionally CH1 is disabled, and CH2 occupies the
            % 1st and only column, adjust i_ch to account for these
            % cases
            Nchannels = sum(info{i_ex}.ignoreChans);
            if (Nchannels == 1) && (i_ch == 2)
                i_ch = 1;
            end
        end
        
        % establish the drug conditions, and the pulse-type conditions (TF
        % and pulse AMP)
        pTypes = fieldnames(dat{i_ex});
        Ntfs = numel(pTypes);
        Nconds = numel(conds);
        
        hFig = figure;
        if Ntfs == 1
            set(gcf, 'position', [414 31 329 754]);
        else
            set(gcf, 'position', [9 10 1135 776]);
        end
        set(gcf, 'name', sprintf('%s, site %.1f, chan: %d', info{i_ex}.mouse, in{i_ex, 2}, i_ch))
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
            for i_tf = 1:Ntfs
                
                
                % check to see if this ptype exists
                if ~isfield(dat{i_ex}.(pTypes{i_tf}), 'stats')
                    continue
                end
                
                % check to see if this drug condition exists
                if ~isfield(dat{i_ex}.(pTypes{i_tf}).snips, conds{i_cond})
                    continue
                end
                
                tmp_raw = dat{i_ex}.(pTypes{i_tf}).snips.(conds{i_cond}){i_ch}';
                Ntime = size(tmp_raw,1);
                Ncols = size(tmp_raw,2);
                
                tt = (0:Ntime-1) ./ info{i_ex}.(pTypes{i_tf}).(conds{i_cond}).sampRate;
                tt = (tt - prePulseTime) .* 1000; % in ms.
                
                subplot(Nconds, Ntfs, i_tf+((i_cond-1) * Ntfs), 'align'),
                cmap = gray(max([Ncols, Ntfs])+3);
                cmap = cmap(1:Ncols, :);
                set(gca, 'colororder', cmap, 'NextPlot', 'replacechildren');
                
                plot(tt, tmp_raw, 'linewidth', 2), hold on,
                
                if CHECK_TRIAL_STATS && isfield(dat{i_ex}.(pTypes{i_tf}).stats, conds{i_cond})
                    
                    % check the peak and trough indicies
                    inds = dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).trpk_inds{i_ch};
                    inds_idx = bsxfun(@plus, inds, Ntime .* (0:Ncols-1)');
                    plot(tt(inds(:,1)), tmp_raw(inds_idx(:,1)), 'ro', 'markerfacecolor', 'r')
                    
                    if any(strcmpi(conds{i_cond}, {'FV_Na', 'FV_Na_Ca2_mGluR'})) || strcmpi(info{i_ex}.opsin, 'estim')
                        plot(tt(inds(:,2)), tmp_raw(inds_idx(:,2)), 'co', 'markerfacecolor', 'c')
                    end
                    
                    % check the area start/stop inds
                    inds = dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).area_inds{i_ch};
                    inds_idx = bsxfun(@plus, inds, Ntime .* (0:Ncols-1)');
                    plot(tt(inds), tmp_raw(inds_idx), 'ko', 'markerfacecolor', 'k')
                    
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
                    if any(strcmpi(conds{i_cond}, {'nbqx_apv_cd2_ttx', 'nbqx_apv_ttx'})) && ~strcmpi(info{i_ex}.opsin, 'estim')
                        for i_pulse = 1:size(tmp_raw,2);
                            startIdx = dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).tau_ind{i_ch}(i_pulse);
                            invtc = dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).tau_invtc{i_ch}(i_pulse);
                            invtc = invtc./1000; % slope was calculated in seconds, but ploting is in ms.
                            beta0 = dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).tau_yint{i_ch}(i_pulse);
                            fit_tt = tt(startIdx:end);
                            fit_vals = beta0 .* exp(invtc .* fit_tt);
                            fit_vals = -fit_vals; % compensate for the fact that the fit was on abs(rawdata), but opsin current is negative
                            plot(fit_tt, fit_vals, 'm-')
                        end
                    end
                    
                end
                
                axis tight
                if i_cond == 1
                    mytitle(pTypes{i_tf})
                end
                if i_cond == numel(conds)
                    xlabel('time (ms)')
                end
                if i_tf == 1
                    if any(regexpi(conds{i_cond}, 'ttx'))
                        switch EXPTTYPE
                            case {'4-AP', 'Baclofen'}
                                h = ylabel(conds{i_cond});
                            otherwise
                                h = ylabel(sprintf('%s Current', info{i_ex}.opsin));
                        end
                        set(h, 'Interpreter', 'none')
                    elseif any(regexpi(conds{i_cond}, 'FV_Na'))
                        h = ylabel(conds{i_cond});
                        set(h, 'Interpreter', 'none')
                    elseif any(regexpi(conds{i_cond}, 'synapticTransmission'))
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
        statTypes = {'diffval', 'area', 'pk2tr', 'slope', 'tau_invtc'}; % could be 'tau_invtc', or 'latency'
        Nstats = numel(statTypes);
        cmap = copper(Ntfs);
        for i_cond = 1:Nconds
            
            for i_stat = 1:Nstats
                
                
                subplot(Nconds, Nstats, i_stat+((i_cond-1)*Nstats)), hold on,
                legtext = {};
                rawvals={};
                tmplegend = [];
                noise_calibration = [];
                for i_tf = 1:Ntfs
                    
                    if ~isfield(dat{i_ex}.(pTypes{i_tf}), 'stats')
                        continue
                    elseif ~isfield(dat{i_ex}.(pTypes{i_tf}).stats, conds{i_cond}) %need to have drug condition
                        continue
                    elseif ~isfield(dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}), statTypes{i_stat}) % need to have stat type
                        continue
                    end
                    
                    
                    tmp_raw = dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).(statTypes{i_stat}){i_ch};
                    tmp_sigma = dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).bkgnd_sigma{i_ch};
                    if NORM_TO_PULSE1
                        tmp_sigma = tmp_sigma ./ tmp_raw(1); % needs to be b/4 the next line b/c tmp_raw is destructively modified
                        tmp_raw = tmp_raw ./ tmp_raw(1);
                    end
                    rawvals = cat(1, rawvals, tmp_raw);
                    
                    noise_calibration = cat(1, noise_calibration, tmp_sigma);
                    switch EXPTTYPE
                        case 'interleaved amps'
                            tmplegend = cat(1, tmplegend, info{i_ex}.(pTypes{i_tf}).(conds{i_cond}).pAmp);
                        otherwise
                            tmplegend = cat(1, tmplegend, info{i_ex}.(pTypes{i_tf}).(conds{i_cond}).pTF);
                    end
                end
                
                
                if isempty(rawvals) % there's no data
                    continue
                end
                
                [~, order] = sort(tmplegend);
                for i_plt = 1:numel(order)
                    idx = order(i_plt);
                    plot(1:numel(rawvals{idx}), rawvals{idx}, 'o-', 'color', cmap(i_plt,:), 'linewidth', 2)
                    legtext = cat(2, legtext, num2str(tmplegend(idx)));
                end
                xlabel('Pulse number')
                xlim([1, max(cellfun(@numel, rawvals))])
                mytitle(sprintf('%s', statTypes{i_stat}));
                if NORM_TO_PULSE1
                    yvals = get(gca, 'ylim');
                    yvals(1) = min([0, yvals(1)]);
                    set(gca, 'ylim', [0, yvals(2)]);
                end
                if i_stat == 1
                    legend(legtext, 'location', 'best')
                    legend boxoff
                end
                
                if any(regexpi(conds{i_cond}, 'ttx'))
                    switch EXPTTYPE
                        case {'4-AP', 'Baclofen'}
                            h = ylabel(conds{i_cond});
                        otherwise
                            h = ylabel(sprintf('%s Current', info{i_ex}.opsin));
                    end
                    set(h, 'Interpreter', 'none')
                elseif any(regexpi(conds{i_cond}, 'FV_Na'))
                    h = ylabel(conds{i_cond});
                    set(h, 'Interpreter', 'none')
                elseif any(regexpi(conds{i_cond}, 'synapticTransmission'))
                    ylabel('Field PSP')
                end
                
                
                if CHECK_TRIAL_STATS
                    if strcmpi(statTypes{i_stat}, 'diffval')
                        if any(strcmpi(conds{i_cond}, {'FV_Na', 'FV_Na_Ca2_mGluR', 'nbqx_apv_cd2_ttx', 'nbqx_apv_ttx'}))
                            xvals = get(gca, 'xlim');
                            yvals = [1-noise_calibration, 1+noise_calibration];
                            for i_plt = 1:numel(order)
                                idx = order(i_plt);
                                %plot(xvals, repmat(yvals(i_plt,:), 2, 1), ':', 'color', cmap(i_plt,:), 'linewidth', 1);
                            end
                        end
                    end
                end
                
                
            end
            
            drawnow % force the uitab plot to update in quasi real time
            
        end
        
    end
end

%% POPULATION SUMMARY PLOTS: MAIN EXPERIMENT



STIMSITE = true;  % true => stimsite,  false => distal site
CALCPVALS = false;
NORMALIZE = true;

% only do this for the main experiment (TF and FV)
assert(~strcmpi(EXPTTYPE, 'interleaved amps'), 'ERROR: this anaysis is only for the TF and FV experiments');


% define the conditions that will get plotted
if COMBINE_CHIEF
    opsinTypes = {'chr2', 'chief_all', 'chronos'};
else
    opsinTypes = {'chr2', 'chief_cit','chief_flx', 'chronos'};
end
switch EXPTTYPE
    case '4-AP'
        conds = {'ttx_cd2', 'ttx_cd2_4AP'};
    case 'Baclofen'
        conds = {'ttx_cd2', 'ttx_cd2_bac10'};
    case 'Intracellular'
        conds = {'nbqx_apv_ttx'};
    case 'E-stim'
        conds = {'nbqx_apv_cd2_ttx', 'FV_Na'};
        opsinTypes = {'estim'}; % overwrite the default
    otherwise
        conds = {'FV_Na', 'nbqx_apv_cd2_ttx', 'synapticTransmission'}; % with cadmium
        %     conds = {'FV_Na', 'FV_Na_Ca2_mGluR', 'synapticTransmission'}; % both %FVs
        %     conds = {'FV_Na_Ca2_mGluR', 'nbqx_apv_ttx',  'synapticTransmission'}; % no cadmium
        %     conds = {'nbqx_apv_cd2_ttx', 'nbqx_apv_ttx',  'synapticTransmission'}; % both opsin current verisons
end

%initialize the outputs
statTypes = {'diffval', 'area', 'pk2tr', 'slope', 'latency'}; % can also include 'tau_invtc' or 'latency'
pop = [];
avgtf = [];
for i_opsin = 1:numel(opsinTypes)
    for i_cond = 1:numel(conds)
        for i_stat = 1:numel(statTypes)
            pop.(opsinTypes{i_opsin}).pnp1.(conds{i_cond}).(statTypes{i_stat}) = {[],[],{}}; % {[TF], [pulseAmp], {Vals}}
            pop.(opsinTypes{i_opsin}).raw.(conds{i_cond}).(statTypes{i_stat}) = {[],[],{}}; % {[TF], [pulseAmp], {Vals}}
            avgtf.(opsinTypes{i_opsin}).avg10to40.(conds{i_cond}).(statTypes{i_stat}) = []; % Nexpts x Npulses
            avgtf.(opsinTypes{i_opsin}).avg10to60.(conds{i_cond}).(statTypes{i_stat}) = []; % Nexpts x Npulses
            avgtf.(opsinTypes{i_opsin}).avg10to100.(conds{i_cond}).(statTypes{i_stat}) = [];
        end
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  sort the data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i_ex = 1:numel(dat)
    
    % figure out which (if any) channels should be analyzed
    CHANNEL = info{i_ex}.stimSite;
    if isnan(CHANNEL) % cases where neither recording sites were targeted
        continue
    end
    Nchannels = sum(info{i_ex}.ignoreChans);
    if Nchannels == 1
        if ~STIMSITE
            continue % no other stim site to show...
        elseif CHANNEL == 2 % cases where CH = 2, but only one channel recorded
            CHANNEL = 1;
        end
    else
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
    
    
    
    
    for i_cond = 1:numel(conds);
      
        for i_stat = 1:numel(statTypes)
            
            Nptypes = numel(pTypes);
            avgOverTF.tf = [];
            avgOverTF.vals = [];
            for i_ptype = 1:Nptypes
                
                % skip instances where the drug condition does not exist
                if ~isfield(dat{i_ex}.(pTypes{i_ptype}).stats, conds{i_cond})
                    continue
                end
                
                % skip instances where the stat type doesn't exist
                if ~isfield(dat{i_ex}.(pTypes{i_ptype}).stats.(conds{i_cond}), statTypes{i_stat})
                    continue
                end
                
                % structure the pnp1 data
                ex_stat = dat{i_ex}.(pTypes{i_ptype}).stats.(conds{i_cond}).(statTypes{i_stat}){CHANNEL};
                if strcmpi(statTypes{i_stat}, 'tau_invtc')
                    ex_stat = 1./ex_stat; % tau_invtc isn't a tau, need to convert to tau
                end
                
                % convert to paired pulse measures by normalizing by the
                % height of the first pulse.
                if NORMALIZE
                    ex_p1p2 = ex_stat ./ ex_stat(1);
                else
                    ex_p1p2 = ex_stat;
                end
                
                % grab the data field of the 'pop' structure
                tmp_pop_p1p2 = pop.(opsin).pnp1.(conds{i_cond}).(statTypes{i_stat});
                tmp_pop_raw = pop.(opsin).raw.(conds{i_cond}).(statTypes{i_stat});
                tf = info{i_ex}.(pTypes{i_ptype}).(conds{i_cond}).pTF;
                pAmp = info{i_ex}.(pTypes{i_ptype}).(conds{i_cond}).pAmp;
                
                % is there already an entry for this TF and drug condition?
                l_tf_match = tmp_pop_p1p2{1} == tf;
                l_pAmp_match = tmp_pop_p1p2{2} == pAmp;
                alreadyThere = any(l_tf_match & l_pAmp_match);
                if alreadyThere
                    
                    idx = l_tf_match & l_pAmp_match;
                    
                    % deal with cases where there are different numbers of
                    % pulses
                    nPulsesExisting = size(tmp_pop_p1p2{3}{idx},2);
                    if numel(ex_p1p2) < nPulsesExisting
                        ex_p1p2(1,end+1:nPulsesExisting) = nan;
                        ex_stat(1,end+1:nPulsesExisting) = nan;
                    else
                        tmp_pop_p1p2{3}{idx}(:,end+1:numel(ex_p1p2)) = nan;
                        tmp_pop_raw{3}{idx}(:,end+1:numel(ex_p1p2)) = nan;
                    end
                    
                    tmp_pop_p1p2{3}{idx} = cat(1, tmp_pop_p1p2{3}{idx}, ex_p1p2);
                    tmp_pop_raw{3}{idx} = cat(1, tmp_pop_raw{3}{idx}, ex_stat);
                    
                else
                    tmp_pop_p1p2{1} = cat(1, tmp_pop_p1p2{1}, tf);
                    tmp_pop_p1p2{2} = cat(1, tmp_pop_p1p2{2}, pAmp);
                    tmp_pop_p1p2{3} = cat(1, tmp_pop_p1p2{3}, ex_p1p2);
                    
                    tmp_pop_raw{1} = cat(1, tmp_pop_raw{1}, tf);
                    tmp_pop_raw{2} = cat(1, tmp_pop_raw{2}, pAmp);
                    tmp_pop_raw{3} = cat(1, tmp_pop_raw{3}, ex_stat);
                end
                
                % replace the 'pop' structure version with the tmp version
                pop.(opsin).pnp1.(conds{i_cond}).(statTypes{i_stat}) = tmp_pop_p1p2;
                pop.(opsin).raw.(conds{i_cond}).(statTypes{i_stat}) = tmp_pop_raw;
                
                % deal with averaging across TFs
                if numel(ex_p1p2)<7;
                    ex_p1p2(end+1) = nan;
                end
                avgOverTF.vals = cat(1, avgOverTF.vals, ex_p1p2(1:7));
                avgOverTF.tf = cat(1, avgOverTF.tf, tf);
            end
            
            % average over TFs, and store the data
            l_10to40 = avgOverTF.tf <= 40;
            if any(l_10to40)
                tmpavg = nanmean(avgOverTF.vals(l_10to40,:), 1);
                avgtf.(opsin).avg10to40.(conds{i_cond}).(statTypes{i_stat})(end+1,:) = tmpavg; % Nexpts x Npulses
            end
            
            if any(avgOverTF.tf <= 60);
                l_10to60 = avgOverTF.tf <= 60;
                if any(l_10to60)
                    tmpavg = nanmean(avgOverTF.vals(l_10to60,:), 1);
                    avgtf.(opsin).avg10to60.(conds{i_cond}).(statTypes{i_stat})(end+1,:) = tmpavg; % Nexpts x Npulses
                end
            end
            
            if any(avgOverTF.tf <= 100);
                l_10to100 = avgOverTF.tf <= 100;
                if any(l_10to100)
                    tmpavg = nanmean(avgOverTF.vals(l_10to60,:), 1);
                    avgtf.(opsin).avg10to100.(conds{i_cond}).(statTypes{i_stat})(end+1,:) = tmpavg; % Nexpts x Npulses
                end
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
    opsinTypes{i_opsin}
    
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
            mytitle(statTypes{i_stat})
            
            tfs = tmp_dat{1};
            tfs = unique(tfs);
            Ntfs = numel(tfs);
            pAmps = tmp_dat{2};
            cmap = copper(max([Ntfs, 5])); % a minimum of 5 tfs
            set(gca, 'colororder', cmap, 'NextPlot', 'replacechildren');
            hold on
            
            legtext = {};
            hvec = [];
            pp_dat_allTFs = {};
            for i_tf = 1:Ntfs;
                
                tf_idx = tmp_dat{1} == tfs(i_tf);
                l_7pulses = cellfun(@(x) size(x,2)>=7, tmp_dat{3});
                tf_idx = tf_idx & l_7pulses; % only analyzes experiments where there are >= 7 pulses per train
                
                % only take the first 7 pulses
                trim_dat = cellfun(@(x) x(:,1:7), tmp_dat{3}(tf_idx), 'uniformoutput', false);
                cat_dat = vertcat(trim_dat{:});
                
                xbar = nanmean(cat_dat, 1);
                sem = nanstd(cat_dat, [], 1) ./ sqrt(sum(~isnan(cat_dat), 1));
                if ~isempty(regexpi(conds{i_cond}, 'FV_Na'))
                    nrepeats = sum(~isnan(cat_dat), 1)
                end
                
                [hvec(i_tf),~] = my_errorbar(1:7, xbar, sem, 'color', cmap(i_tf,:), 'linewidth', 2); % only plot the first 7 pulses
                legtext = cat(2, legtext, num2str(tfs(i_tf)));
                
                % store the pp data in an array that I can use to perform
                % inferential stats
                pp_dat_allTFs{i_tf} = cat_dat;
                
            end
            axis tight
            xlabel('Pulse number')
            if i_stat==1
                ylabel(sprintf('%s \n Pn:P1 ratio', conds{i_cond}))
            end
            if i_stat==1  && i_cond==1
                legend(hvec, legtext, 'location', 'best')
                legend boxoff
            end
            yvals = get(gca, 'ylim');
            yvals(1) = min([0, yvals(1)]);
            set(gca, 'ylim', yvals);
            if yvals(1)<0
                plot([1,numel(xbar)], [0,0] , 'k--', 'linewidth', 2)
            end
            
            
            if CALCPVALS
                
                whichOpsin = 'estim';    % 'chronos', 'chr2', 'chief_all', 'chief_cit', 'chief_flx'
                whichCond =  'Fv_Na';   % 'Fv_Na', 'nbqx_apv_cd2_ttx', 'synapticTransmission'
                whichStat = 'latency';  % 'diffval', 'slope', 'latency'
                
                opsinMatch = strcmpi(opsinTypes{i_opsin}, whichOpsin);
                condMatch = strcmpi(conds{i_cond}, whichCond);
                statMatch = strcmpi(statTypes{i_stat}, whichStat);
                
                if all([opsinMatch, condMatch, statMatch])
                    keyboard
                    %
                    %  Do some inferential tests and present a table with the
                    %  results
                    %
                    %%%%%%%%%%%%%%%%%%%%%%%%%%
                    [p2_for_each_tf, p7_for_each_tf, p1_p7_across_tfs, p1_p2_across_tfs] = deal(nan(Ntfs, 1));
                    
                    % Test for significant decreases in PPRs. Look specifically at
                    % P2:P1. Do the test individually for each TF
                    %
                    % Ho -> distribution of PPRs has Xbar = 1;
                    for i_tf = 1:Ntfs
                        tmp = pp_dat_allTFs{i_tf}(:,2);
                        tmp(tmp<0) = 0.01;
                        [~, p2_for_each_tf(i_tf, 1), ~, stattab] = ttest(log(tmp))
                    end
                    
                    % Test for significant decreases in PPRs. Look specifically at
                    % P7:P1. Do the test individually for each TF
                    %
                    % Ho -> distribution of PPRs has Xbar = 1;
                    for i_tf = 1:Ntfs
                        tmp = pp_dat_allTFs{i_tf}(:,7);
                        tmp(tmp<0) = 0.01;
                        [~, p7_for_each_tf(i_tf, 1), ~, stattab] = ttest(log(tmp))
                    end
                    
                    
                    % Test for differences in P7:P1 across TFs. Here I'm asking if
                    % there's a significant effect of TF on the PPR. Do an anova
                    % like test with post-hoc comparisons
                    %
                    % Ho -> distributions of PPRs have the identical median
                    tmp_dat_p2 = [];
                    tmp_dat_p7 = [];
                    tmp_group = [];
                    for i_tf = 1:Ntfs
                        tmp_dat_p2 = cat(1, tmp_dat_p2, pp_dat_allTFs{i_tf}(:,2));
                        tmp_dat_p7 = cat(1, tmp_dat_p7, pp_dat_allTFs{i_tf}(:,7));
                        tmp_group = cat(1, tmp_group, ones(size(pp_dat_allTFs{i_tf}(:,7))).*i_tf);
                    end
                    tmp_dat_p7(tmp_dat_p7<0) = 0.01; % so that the log doesn't bonk
                    [p1_p7_across_tfs(1), tab, stattab] = anova1(log(tmp_dat_p7), tmp_group);
                    multcompare(stattab);
                    
                    
                    % Test for differences in P2:P1 across TFs. Here I'm asking if
                    % there's a significant effect of TF on the PPR. Do an anova
                    % like test with post-hoc comparisons
                    %
                    % Ho -> distributions of PPRs have the identical median
                    tmp_dat_p2(tmp_dat_p2<0) = 0.01; % so that the log doesn't bonk
                    [p1_p2_across_tfs(1), tab, stattab] = anova1(log(tmp_dat_p2), tmp_group);
                    multcompare(stattab);
                    
                    
                    fprintf('  comparisons for %s values for %s  \n',...
                        upper(opsinTypes{i_opsin}), upper(conds{i_cond}))
                    
                    rownames = cellfun(@num2str, mat2cell(tfs, ones(size(tfs))), 'uniformoutput', false);
                    T = table(p2_for_each_tf, p1_p2_across_tfs, p7_for_each_tf, p1_p7_across_tfs, 'RowNames', rownames)
                    fprintf('\n\n')
                end
            end
            
        end
    end
end


%% POPULATION SUMMARY: INTERLEAVED POWERS


% plan: loop over opsins. Only consider a single recording channel (distal
% or proximal). Show PP ratio as a function of TF. Show Pn:P1 ratio as a
% function of pulse number for each frequency

STIMSITE = true;  % true => stimsite,  false => distal site
PULSENUM = 1;
PLOTPPR = false; % if true, PPR = P(PULSENUM) ./ P(1), else the stat is raw P(PULSENUM)
NORM_TO_LED10 = false; % norm data to LED=10V
STATTYPE = 'diffval';
CONDITION = 'FV_Na'; % 'FV_Na', or 'nbqx_apv_cd2_ttx'

% only do this for the main experiment (TF and FV)
assert(strcmpi(EXPTTYPE, 'interleaved amps'), 'ERROR: this anaysis is only for experiments with interleaved amplitudes');

% load in the calibration data to convert from Volts to power denisty
load led_470_cal.mat;
cal_led = cal;
load laser_450_cal.mat;
cal_laser = cal;

% pull out the waveform of the first pulse, keep track of the power (TF is
% irrelevant) so that I can average acorss mice
firstPulse = [];
firstPulse.wf = {}; % waveforms, 1 cell for each mouse, 1 row for each amp
firstPulse.amp = {}; % pAmps, 1 cell for each mouse, 1 row for each amp
firstPulse.tf = {}; % just make sure that only 1 TF was tested per mouse
firstPulse.stats = {}; % for the line series plots
firstPulse.opsin = {}; % 1 opsin per mouse

allPulses = [];
allPulses.rawstats = {};


for i_ex = 1:numel(dat)
    
    
    % figure out which (if any) channels should be analyzed
    CHANNEL = info{i_ex}.stimSite;
    if isnan(CHANNEL) % cases where neither recording sites were targeted
        continue
    end
    Nchannels = sum(info{i_ex}.ignoreChans);
    if Nchannels == 1
        if ~STIMSITE
            continue % no other stim site to show...
        elseif CHANNEL == 2 % cases where CH = 2, but only one channel recorded
            CHANNEL = 1;
        end
    else
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
    
    % get P2:P1 ratios
    ex_tf = [];
    ex_powDensity = [];
    ex_ppr = [];
    
    pTypes = fieldnames(dat{i_ex});
    for i_ptype = 1:numel(pTypes);
        
        if ~isfield(dat{i_ex}.(pTypes{i_ptype}).stats, CONDITION)
            continue
        end
        
        
        % deal with the raw data
        tmp = dat{i_ex}.(pTypes{i_ptype}).stats.(CONDITION).(STATTYPE){CHANNEL};
        if numel(tmp)<7; continue; end;
        if PLOTPPR
            assert(PULSENUM>1, 'ERROR: for PPR, PULSENUM must be > 1')
            amp = tmp(PULSENUM)/tmp(1);
        else
            amp = tmp(PULSENUM);
        end
        
        % convert amplitudes from volts to mW/mm2. Treat the LED
        % differently then the Laser
        cmd_voltage = info{i_ex}.(pTypes{i_ptype}).(CONDITION).pAmp;
        switch info{i_ex}.(pTypes{i_ptype}).(CONDITION).stimtype
            case 'laser'
                mW_mm2 = ppval(cal_laser.pp_mW_mm2, cmd_voltage);
            case 'led'
                mW_mm2 = ppval(cal_led.pp_mW_mm2, cmd_voltage);
        end
        
        % compile the data for a single experiment
        ex_tf = cat(1, ex_tf, info{i_ex}.(pTypes{i_ptype}).(CONDITION).pTF);
        assert(all(ex_tf == 20), 'ERROR: not all TFs equal 20')
        
        % store the representative pulse waveform
        tmp_wf = dat{i_ex}.(pTypes{i_ptype}).snips.(CONDITION){CHANNEL}(PULSENUM,:);
        
        firstPulse.wf{i_ex}(i_ptype,:) = tmp_wf;
        firstPulse.mW_mm2{i_ex}(i_ptype) = mW_mm2; 
        firstPulse.tf{i_ex}(i_ptype) = info{i_ex}.(pTypes{i_ptype}).(CONDITION).pTF; 
        firstPulse.opsin{i_ex} = info{i_ex}.opsin;
        firstPulse.stats{i_ex}(i_ptype) = amp;
        
        % store all the raw pulse stats
        if strcmpi(info{i_ex}.opsin, 'chr2')
            tmp = dat{i_ex}.(pTypes{i_ptype}).stats.(CONDITION).(STATTYPE){CHANNEL};
            allPulses.rawstats{i_ex}(i_ptype,:) = tmp(:)';
            allPulses.mW_mm2{i_ex}(i_ptype) = mW_mm2;
        end
        
    end

end


%
% Plot the shape of the Nth pulse across pAmps, and a line series showing
% the dependence of the stat type on the light power. This shows only the
% Nth pulse
%
opsins = {'chr2'};

for i_opsin = 1:numel(opsins)
    
    idx_opsin = strcmpi(firstPulse.opsin, opsins{i_opsin});
    idx_opsin = find(idx_opsin);
    
    opsin_wfs = [];
    opsin_amps = [];
    opsin_stats = [];
    for i_ex = idx_opsin
        
        tmp_tf = firstPulse.tf{i_ex};
        tmp_pamp = round(firstPulse.mW_mm2{i_ex}, 1);
        tmp_wf = firstPulse.wf{i_ex};
        tmp_stats = firstPulse.stats{i_ex};
        
        % I'm going to normalize to the value at 10Volts. so the data need
        % to contain some 10 volt trials
        if NORM_TO_LED10
            
            assert(~PLOTPPR, 'ERROR: PLOTPPR can not be true')
            
            normidx = tmp_pamp == 25.7; % in mW_mm2
            if ~any(normidx); continue; end
            
            % normalize the wfs
            normval = min(tmp_wf(normidx,:));
            tmp_wf = tmp_wf ./ abs(normval); % preserve the sign
            
            
            % normalize the stats
            normval = tmp_stats(normidx);
            tmp_stats = tmp_stats ./ normval; % do not preserve the sign (all positive)
        end
        
        opsin_wfs = cat(1, opsin_wfs, tmp_wf);
        opsin_amps = cat(1, opsin_amps, tmp_pamp(:));
        opsin_stats = cat(1, opsin_stats, tmp_stats(:));
        
    end
 
    % plot of waveforms
    legtext = {};
    f = figure; hold on,
    f.Name = opsins{i_opsin};
    f.Position = [ 244   307   425   443];
    unique_amps = unique(opsin_amps);
    for i_amp = 1:numel(unique_amps)
        pamp = unique_amps(i_amp);
        l_wfs = opsin_amps == pamp;
        
        tt = ((0:size(opsin_wfs, 2)-1)-prePulseSamps) ./ sampRate;
        tt = tt.*1000;
        plot(tt, mean(opsin_wfs(l_wfs,:), 1), 'linewidth', 2)
        
        legtext = cat(2, legtext, [num2str(pamp) 'mW per mm2']);
    end
    axis tight
    legend(legtext, 'Location', 'best')
    legend boxoff
    xlabel('time (ms)')
    if NORM_TO_LED10
        ylabel('amplitude (norm to 10V LED)')
    else
        ylabel(sprintf('raw amplitudes, pulse=%d', PULSENUM))
    end
    
    % line series plot of stats
    plt_pow = [];
    plt_xbar = [];
    plt_sem = [];
    for i_amp = 1:numel(unique_amps)
        pamp = unique_amps(i_amp);
        l_wfs = opsin_amps == pamp;
        
        plt_pow(i_amp) = pamp;
        plt_xbar(i_amp) = mean(opsin_stats(l_wfs));
        plt_sem(i_amp) = stderr(opsin_stats(l_wfs));
    end
    
    f = figure;
    my_errorbar(plt_pow, plt_xbar, plt_sem, '-ks', 'markersize', 8, 'markerfacecolor', 'k', 'linewidth', 2);
    f.Name = opsins{i_opsin};
    f.Position = [660   172   464   620];
    axis tight
    xlabel('Power (mW per mm2)')
    xlim([min(plt_pow)-1, max(plt_pow)+1])
    box off
    set(gca, 'xscale', 'log')
    if NORM_TO_LED10
        ylabel('amplitude (norm to LED=10V)')
    elseif PLOTPPR
        ylabel(sprintf('paired pulse ratio: P%d:P1', PULSENUM))
        ylim([0 1])
    else
        ylabel(sprintf('raw amplitudes, pulse=%d', PULSENUM))
        set(gca, 'XAxisLocation', 'top', 'xtick', [10, 100])
        ymin = min(get(gca, 'ylim'));
        ylim([ymin*1.05, 0])
    end
    
    
    
    
    
end

%
% Now plot all pulses and all amplitudes. Normalize to the first pulse, and
% average across experiments. Make sure that the pAmps are ordered
% consistently across experiments
%

all_pAmps = cat(1, allPulses.mW_mm2{:});
unique_pAmps = unique(all_pAmps, 'rows');
assert(size(unique_pAmps,1)==1, 'ERROR: pulse amps are out of order')

% grab the data, normalize, and average
all_rawstats = cat(3, allPulses.rawstats{:});
all_normstats = bsxfun(@rdivide, all_rawstats, all_rawstats(:,1,:));

all_avg_normstats = mean(all_normstats, 3);
all_sem_normstats = stderr(all_normstats, 3);

f = figure;
f.Position = [440   285   469   513];
nAmps = numel(unique_pAmps);
nPulses = size(all_avg_normstats, 2);
cmap = winter(nAmps);
set(gca, 'colororder', cmap, 'NextPlot', 'replacechildren');
hold on,
hvec = [];
for i_amp = 1:nAmps
    [hvec(i_amp),~] = my_errorbar(1:nPulses, all_avg_normstats(i_amp,:), all_sem_normstats(i_amp,:), '-', 'color', cmap(i_amp,:), 'linewidth', 3);
end
yvals = get(gca,'ylim');
set(gca, 'ylim', [0, max([yvals(2), 1])]);
xlabel('pulse number')
hy = ylabel(sprintf('normalized %s %s', CONDITION, STATTYPE));
hy.Interpreter = 'none';
legtext = cellfun(@(x) [num2str(x, '%3.1f\n')  ' mW/mm2'], num2cell(unique_pAmps), 'uniformoutput', false);
legend(hvec, legtext, 'location', 'southwest')
legend boxoff

% do some inferential tests on the PPRs
P2_vals = permute(all_normstats(:,2,:), [3, 1, 2]) % dimensionality is [Nexp x Npowers];
[P2_rho_on_means, Pval] = corr(mean(P2_vals, 1)', unique_pAmps(:), 'type', 'spearman')


% is there a differences in the P2:P1 vals across light powers
P2_vals = permute(all_normstats(:,2,:), [3, 1, 2]) % dimensionality is [Nexp x Npowers];
groups = repmat([1:size(P2_vals, 2)], size(P2_vals, 1), 1);
P2_vals = P2_vals(:);
groups = groups(:);
[pval, ~, stattab] = anova1(log(P2_vals), groups)


P7_vals = permute(all_normstats(:,7,:), [3, 1, 2]) % dimensionality is [Nexp x Npowers];
[P7_rho_on_means, Pval] = corr(mean(P7_vals, 1)', unique_pAmps(:), 'type', 'spearman')


% is there a differences in the P7:P1 vals across light powers
P7_vals = permute(all_normstats(:,7,:), [3, 1, 2]) % dimensionality is [Nexp x Npowers];
groups = repmat([1:size(P7_vals, 2)], size(P7_vals, 1), 1);
P7_vals = P7_vals(:);
groups = groups(:);
[pval, ~, stattab] = anova1(log(P7_vals), groups)

%% OPSIN CURRENT VS. FIBER VOLLEY [SPIDER PLOT]

% only do this for the main experiment (TF and FV)
assert(strcmpi(EXPTTYPE, 'main expt'), 'ERROR: this anaysis is only for the TF and FV experiments');


close all

FV_STAT = 'diffval';
OPSIN_STAT = 'diffval';
STAT_TYPE = 'pnp1';  % can be 'pnp1' or 'raw'
PLOT_RAW = true;
PLOT_AVG = false;
PLOT_ERR = false;
REMOVE_OUTLIER = true;
INDIVIDUAL_PLOTS = true;
N_PULSES = 7;
TFs = [10,20,40,60,100];


if ~INDIVIDUAL_PLOTS
    figure, hold on,
    cmap = colormap('copper');
    cidx = round(linspace(1, size(cmap,1), N_PULSES+1));
    cmap = cmap(cidx(1:end-1),:);
end

TF_line_scalar = linspace(0,0.7,numel(TFs));

for i_tf = 1:numel(TFs)
    
    if INDIVIDUAL_PLOTS
        figure, hold on,
        cmap = colormap('copper');
        cidx = round(linspace(1, size(cmap,1), N_PULSES+1));
        cmap = cmap(cidx(1:end-1),:);
    end
    
    % CAH new: nned to deal with sum(idx == 1) > 1, which occurs now b/c
    % there is more than one row for each tf.
    
    % pull out the data for ChR2:
    tf_idx = pop.chr2.(STAT_TYPE).FV_Na.(FV_STAT){1} == TFs(i_tf);
    if sum(tf_idx) >= 1;
        
        trim_dat = cellfun(@(x) x(:,1:7), pop.chr2.(STAT_TYPE).FV_Na.(FV_STAT){3}(tf_idx), 'uniformoutput', false);
        chr2_fv_raw = vertcat(trim_dat{:});
        
        trim_dat = cellfun(@(x) x(:,1:7), pop.chr2.(STAT_TYPE).nbqx_apv_cd2_ttx.(OPSIN_STAT){3}(tf_idx), 'uniformoutput', false);
        chr2_opsin_raw = vertcat(trim_dat{:});
    else
        chr2_fv_raw = nan(5,7); % a hack to allow plotting of tf conds for oChIEF that don't exist for ChR2
        chr2_opsin_raw = nan(5,7);
    end
    
    
    % pull out the data for oChIEF:
    tf_idx = pop.chief_cit.(STAT_TYPE).FV_Na.(FV_STAT){1} == TFs(i_tf);
    if sum(tf_idx) >= 1;
        trim_dat = cellfun(@(x) x(:,1:7), pop.chief_cit.(STAT_TYPE).FV_Na.(FV_STAT){3}(tf_idx), 'uniformoutput', false);
        ochief_fv_raw = vertcat(trim_dat{:});
        
        trim_dat = cellfun(@(x) x(:,1:7), pop.chief_cit.(STAT_TYPE).nbqx_apv_cd2_ttx.(OPSIN_STAT){3}(tf_idx), 'uniformoutput', false);
        ochief_opsin_raw = vertcat(trim_dat{:});
    else
        ochief_fv_raw = nan(5,7); % a hack to allow plotting of tf conds for oChIEF that don't exist for ChR2
        ochief_opsin_raw = nan(5,7);
    end
    if REMOVE_OUTLIER && any(ochief_opsin_raw(:)>1400)
        l_bad = sum(ochief_opsin_raw>1400,2) > 0;
        ochief_opsin_raw(l_bad,:) = [];
        ochief_fv_raw(l_bad,:) = [];
    end
    
    % pull out the data for Chronos:
    tf_idx = pop.chronos.(STAT_TYPE).FV_Na.(FV_STAT){1} == TFs(i_tf);
    if sum(tf_idx) >= 1;
        trim_dat = cellfun(@(x) x(:,1:7), pop.chronos.(STAT_TYPE).FV_Na.(FV_STAT){3}(tf_idx), 'uniformoutput', false);
        chronos_fv_raw = vertcat(trim_dat{:});
        
        trim_dat = cellfun(@(x) x(:,1:7), pop.chronos.(STAT_TYPE).nbqx_apv_cd2_ttx.(OPSIN_STAT){3}(tf_idx), 'uniformoutput', false);
        chronos_opsin_raw = vertcat(trim_dat{:});
    else
        chronos_fv_raw = nan(5,7); % a hack to allow plotting of tf conds for oChIEF that don't exist for ChR2
        chronos_opsin_raw = nan(5,7);
    end
    
    
    ochiefclr = [1, TF_line_scalar(i_tf), TF_line_scalar(i_tf)];
    chr2clr = [TF_line_scalar(i_tf), TF_line_scalar(i_tf), 1];
    chronosclr = [TF_line_scalar(i_tf), 1, TF_line_scalar(i_tf)];
    
    if PLOT_RAW
        for i_pulse = 1:N_PULSES
            plot(chr2_opsin_raw(:,i_pulse), chr2_fv_raw(:,i_pulse), '+', 'markeredgecolor', chr2clr)
            plot(ochief_opsin_raw(:,i_pulse), ochief_fv_raw(:,i_pulse), 'o', 'markeredgecolor', ochiefclr)
            plot(chronos_opsin_raw(:,i_pulse), chronos_fv_raw(:,i_pulse), 'o', 'markeredgecolor', chronosclr)
        end
    end
    
    if PLOT_AVG
        % plot the avereages across conditions
        plot(mean(chr2_opsin_raw(:,1:N_PULSES),1), mean(chr2_fv_raw(:,1:N_PULSES),1), '-s', 'markerfacecolor', chr2clr, 'color', chr2clr, 'markeredgecolor', chr2clr, 'markersize', 10, 'linewidth', 2)
        plot(mean(ochief_opsin_raw(:,1:N_PULSES),1), mean(ochief_fv_raw(:,1:N_PULSES),1), '-s', 'markerfacecolor', ochiefclr, 'color', ochiefclr, 'markeredgecolor', ochiefclr,  'markersize', 10, 'linewidth', 2)
        plot(mean(chronos_opsin_raw(:,1:N_PULSES),1), mean(chronos_fv_raw(:,1:N_PULSES),1), '-s', 'markerfacecolor', chronosclr, 'color', chronosclr, 'markeredgecolor', chronosclr,  'markersize', 10, 'linewidth', 2)
        
        % plot SEM for the various conditions
        if PLOT_ERR
            tmpX = mean(chr2_opsin_raw(:,1:N_PULSES),1);
            tmpY = mean(chr2_fv_raw(:,1:N_PULSES),1);
            SEM_X = nanstd(chr2_opsin_raw(:,1:N_PULSES), [] ,1) ./ sqrt(sum(~isnan(chr2_opsin_raw(:,1:N_PULSES)), 1));
            SEM_Y = nanstd(chr2_fv_raw(:,1:N_PULSES), [] ,1) ./ sqrt(sum(~isnan(chr2_fv_raw(:,1:N_PULSES)), 1));
            plot([tmpX-SEM_X/2; tmpX+SEM_X/2] , [tmpY; tmpY], 'color', chr2clr)
            plot([tmpX ; tmpX], [tmpY-SEM_Y/2; tmpY+SEM_Y/2], 'color', chr2clr)
            
            tmpX = mean(ochief_opsin_raw(:,1:N_PULSES),1);
            tmpY = mean(ochief_fv_raw(:,1:N_PULSES),1);
            SEM_X = nanstd(ochief_opsin_raw(:,1:N_PULSES), [] ,1) ./ sqrt(sum(~isnan(ochief_opsin_raw(:,1:N_PULSES)), 1));
            SEM_Y = nanstd(ochief_fv_raw(:,1:N_PULSES), [] ,1) ./ sqrt(sum(~isnan(ochief_fv_raw(:,1:N_PULSES)), 1));
            plot([tmpX-SEM_X/2; tmpX+SEM_X/2] , [tmpY; tmpY], 'color', ochiefclr)
            plot([tmpX ; tmpX], [tmpY-SEM_Y/2; tmpY+SEM_Y/2], 'color', ochiefclr)
            
            tmpX = mean(chronos_opsin_raw(:,1:N_PULSES),1);
            tmpY = mean(chronos_fv_raw(:,1:N_PULSES),1);
            SEM_X = nanstd(chronos_opsin_raw(:,1:N_PULSES), [] ,1) ./ sqrt(sum(~isnan(chronos_opsin_raw(:,1:N_PULSES)), 1));
            SEM_Y = nanstd(chronos_fv_raw(:,1:N_PULSES), [] ,1) ./ sqrt(sum(~isnan(chronos_fv_raw(:,1:N_PULSES)), 1));
            plot([tmpX-SEM_X/2; tmpX+SEM_X/2] , [tmpY; tmpY], 'color', chronosclr)
            plot([tmpX ; tmpX], [tmpY-SEM_Y/2; tmpY+SEM_Y/2], 'color', chronosclr)
        end
    end
    
    % axis labels and such
    xlabel(sprintf('Opsin current (%s)', OPSIN_STAT));
    ylabel(sprintf('Fiber Volley (%s)', FV_STAT));
    
    if strcmpi(STAT_TYPE, 'pnp1')
        xlim([0, 1.05])
        ylim([0, 1.05])
    end
end





%% OPSIN CURRENT VS. FIBER VOLLEY FOR [FIRST PULSE ONLY]


FV_STAT = 'diffval';
OPSIN_STAT = 'diffval';
STIMSITE = true;
NORMTOMAX = false;

[pop_fv, pop_opsin, pop_pamp] = deal({});

% load in the LED calibration data to convert from Volts to power denisty
load led_470_cal.mat % data now stored in 'cal' struct

figure, hold on,
for i_ex = 1:numel(dat)
     
    % figure out which (if any) channels should be analyzed
    CHANNEL = info{i_ex}.stimSite;
    if isnan(CHANNEL) % cases where neither recording sites were targeted
        continue
    end
    Nchannels = sum(info{i_ex}.ignoreChans);
    if Nchannels == 1
        if ~STIMSITE
            continue % no other stim site to show...
        elseif CHANNEL == 2 % cases where CH = 2, but only one channel recorded
            CHANNEL = 1;
        end
    else
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
    
    
    
    % there could be multiple TFs and pAmps. figure out what the deal is.
    pTypeNames = fieldnames(dat{i_ex});
    ex_pAmps = nan(numel(pTypeNames), 1);
    ex_TFs = nan(numel(pTypeNames), 1);
    for i_pType = 1:numel(pTypeNames);
        if ~isfield(info{i_ex}.(pTypeNames{i_pType}), 'FV_Na')
            continue
        end
        ex_pAmps(i_pType) = info{i_ex}.(pTypeNames{i_pType}).nbqx_apv_cd2_ttx.pAmp;
        ex_TFs(i_pType) = info{i_ex}.(pTypeNames{i_pType}).nbqx_apv_cd2_ttx.pTF;
    end
    
    
    % for each pAmp, average the first pulse stats across available TFs
    pAmps_unique = unique(ex_pAmps);
    [opsin_avg, fv_avg] = deal(nan(1, numel(pAmps_unique)));
    for i_pamp = 1:numel(pAmps_unique);
        
        % find the field names with the TF conds that have a particular
        % pulse amplitude.
        idx = ex_pAmps == pAmps_unique(i_pamp);
        fldnames = pTypeNames(idx);
        
        [opsin_raw, fv_raw] = deal(nan(numel(fldnames), 1));
        for i_fld = 1:numel(fldnames)
            
            if ~isfield(dat{i_ex}.(fldnames{i_fld}).stats, 'FV_Na')
                continue
            end
            
            % get the FV_Na stat
            tmp = dat{i_ex}.(fldnames{i_fld}).stats.FV_Na.(FV_STAT){CHANNEL}(1);
            if strcmpi(FV_STAT, 'diffval')
                tmp = abs(tmp);
            end
            fv_raw(i_fld) = tmp;
            
            % get the opsin stat
            tmp = dat{i_ex}.(fldnames{i_fld}).stats.nbqx_apv_cd2_ttx.(OPSIN_STAT){CHANNEL}(1);
            if strcmpi(OPSIN_STAT, 'diffval')
                tmp = abs(tmp);
            end
            opsin_raw(i_fld) = tmp;
        end
        
        % average across TF conditions
        fv_avg(i_pamp) = mean(fv_raw);
        opsin_avg(i_pamp) = mean(opsin_raw);
        
    end
    
    % Normalize the values if desired. The data should be in order from
    % smallest pulse amp to biggest pulse amp
    if NORMTOMAX
        [~, sortidx] = sort(pAmps_unique);
        assert(all(sortidx == [1:numel(pAmps_unique)]'), 'ERROR: pAmps out of order')
        fv_avg = fv_avg ./ fv_avg(end);
        opsin_avg = opsin_avg ./ opsin_avg(end);
    end
    
    % compile the data across experiments
    pop_fv{i_ex} = fv_avg;
    pop_opsin{i_ex} = opsin_avg;
    pop_pamp{i_ex} =  ppval(cal.pp_mW_mm2, pAmps_unique);
    
    
    
    switch lower(info{i_ex}.opsin)
        case 'chr2'
            p = plot(opsin_avg, fv_avg, 'b-o', 'markerfacecolor', 'b');
        case {'chief_cit', 'chief_all'}
            p = plot(opsin_avg, fv_avg, 'r-o', 'markerfacecolor', 'r');
        case 'chief_flx'
            p = plot(opsin_avg, fv_avg, 'm-o', 'markerfacecolor', 'm');
        case 'chronos'
            p = plot(opsin_avg, fv_avg, 'g-o', 'markerfacecolor', 'g');
    end
    
    % axis labels and such
    xlabel(sprintf('Opsin current (%s)', OPSIN_STAT));
    ylabel(sprintf('Fiber Volley (%s)', FV_STAT));
    bdfxn = @(x,y,z,zz) mytitle(sprintf('%s, site: %d', z, zz));
    set(p, 'buttonDownFcn', {bdfxn, in{i_ex,1}, in{i_ex,2}})
    if NORMTOMAX; ylim([0 1]); xlim([0 1]); end
end


% plot the opsin current and FV against light power
figure
for i_ex = 1:numel(dat)
    
    switch lower(info{i_ex}.opsin)
        case 'chr2'
            pltclr = 'b';
        case {'chief_cit', 'chief_all'}
            pltclr = 'r';
        case 'chief_flx'
            pltclr = 'm';
        case 'chronos'
            pltclr = 'g';
    end
    
    % opsin first
    subplot(1,2,1), hold on
    plot(pop_pamp{i_ex}, pop_opsin{i_ex}, '-o', 'markerfacecolor', pltclr, 'color', pltclr);
    xlabel('Light Power')
    ylabel('Opsin Curent')
    
    
    % FV second
    subplot(1,2,2), hold on
    plot(pop_pamp{i_ex}, pop_fv{i_ex}, '-o', 'markerfacecolor', pltclr, 'color', pltclr);
    xlabel('Light Power')
    ylabel('FV amplitude')
end

opsinTypes = {'chr2', 'chronos', 'chief_all', 'chief_flx', 'chief_cit'};
popPower = [];
for i_opsin = 1:numel(opsinTypes)
    
    l_opsin = strcmpi(opsinTypes{i_opsin}, structcat(info, 'opsin'));
    tmpdat = cat(1,pop_pamp{l_opsin});
    tmpdat(isnan(tmpdat)) = []; % in case there was no data
    popPower.(opsinTypes{i_opsin}).avg = mean(tmpdat);
    popPower.(opsinTypes{i_opsin}).sem = stderr(tmpdat);
    popPower.(opsinTypes{i_opsin}).allPows = tmpdat;
    
    % store the raw opsin currents
    tmpdat = cat(2, pop_opsin{l_opsin})';
    tmpdat(isnan(tmpdat)) = [];
    popPower.(opsinTypes{i_opsin}).allOcurrents = tmpdat;
    
end

% Kruskal-wallis for light power
groupName = {};
groupData = [];
for i_opsin = 1:numel(opsinTypes)
   tmpdata = popPower.(opsinTypes{i_opsin}).allPows;
   groupData = cat(1, groupData, tmpdata);
   groupName = cat(1, groupName, repmat({opsinTypes{i_opsin}}, numel(tmpdata), 1));
end
[P,ANOVATAB,STATS] = kruskalwallis(groupData, groupName, 'on');
multcompare(STATS)

% Anova for opsin current
groupName = {};
groupData = [];
for i_opsin = 1:numel(opsinTypes)
   tmpdata = popPower.(opsinTypes{i_opsin}).allOcurrents;
   groupData = cat(1, groupData, tmpdata);
   groupName = cat(1, groupName, repmat({opsinTypes{i_opsin}}, numel(tmpdata), 1));
end
[P,ANOVATAB,STATS] = anova1(groupData, groupName, 'on');
multcompare(STATS)

%%  SHAPE OF THE FIRST PULSE RESPONSE

% only do this for the main experiment (TF and FV)
assert(any(strcmpi(EXPTTYPE, {'main expt', 'Intracellular', '4-AP'})), 'ERROR: this anaysis is only for the TF and FV experiments');


STIMSITE = true;
PLOTERR = true;

switch EXPTTYPE
    case 'Intracellular'
        conds = {'nbqx_apv_ttx'};
    case '4-AP'
        conds = {'ttx_cd2', 'ttx_cd2_4AP'};
    otherwise
        conds = {'nbqx_apv_cd2_ttx', 'FV_Na'};
end

chr2_examp = {[] []};
chief_cit_examp = {[] []};
chief_flx_examp = {[] []};
chief_all_examp = {[] []};
chronos_examp = {[] []};
Nexpts = size(in,1);

for i_ex = 1:Nexpts
    
    % figure out which (if any) channels should be analyzed
    CHANNEL = info{i_ex}.stimSite;
    if isnan(CHANNEL) % cases where neither recording sites was targeted
        continue
    end
    Nchannels = sum(info{i_ex}.ignoreChans);
    if Nchannels == 1
        if ~STIMSITE
            continue % no other stim site to show...
        elseif CHANNEL == 2 % cases where CH = 2, but only one channel recorded
            CHANNEL = 1;
        end
    else
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
    
    
    % run the analysis
    for i_cond = 1:numel(conds)
        
        pTypes = fieldnames(dat{i_ex});
        Ntfs = numel(pTypes);
        
        if ~isfield(info{i_ex}.(pTypes{1}), conds{i_cond})
            continue
        end
        
        sampRate = info{i_ex}.(pTypes{1}).(conds{i_cond}).sampRate;
        prePulseSamps = ceil(prePulseTime .* sampRate); % samples prior to pulse onset
        postPulseSamps = ceil(postPulseTime .* sampRate); % samples available after pulse ONSET
        
        
        firstPulse = nan(Ntfs, prePulseSamps+postPulseSamps+1);
        for i_tf = 1:Ntfs
            
            if any(strcmpi(pTypes{i_tf}, {'tf100', 'tf60'}))
                continue
            end
            
            tmp_dat = dat{i_ex}.(pTypes{i_tf}).snips.(conds{i_cond}){CHANNEL}(1,:);
            normFact = dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).diffval{CHANNEL}(1);
            normFact = min(tmp_dat);
            normFact = normFact * -1;
            
            firstPulse(i_tf,:) = tmp_dat ./ normFact;
            
        end
        
        % take the average, and normalize
        firstPulse = nanmean(firstPulse,1);
        
        
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
            case 'chief_cit'
                chief_cit_examp{i_cond} = cat(1, chief_cit_examp{i_cond}, firstPulse);
            case 'chief_flx'
                chief_flx_examp{i_cond} = cat(1, chief_flx_examp{i_cond}, firstPulse);
            case 'chief_all'
                chief_all_examp{i_cond} = cat(1, chief_all_examp{i_cond}, firstPulse);
            case 'chr2'
                chr2_examp{i_cond} = cat(1, chr2_examp{i_cond}, firstPulse);
            case 'chronos'
                chronos_examp{i_cond} = cat(1, chronos_examp{i_cond}, firstPulse);
        end
        
    end % conditions
    
end % expts


for i_cond = 1:numel(conds)
    
    % the example first
    tt = (0:size(chr2_examp{i_cond},2)-1)./20e3;
    tt = tt-tt(prePulseSamps);
    tt = tt.*1000;
    
    figure, hold on,
    if ~PLOTERR
%         plot(tt, chief_cit_examp{i_cond}', 'r');
%         plot(tt, chief_flx_examp{i_cond}', 'm');
%         plot(tt, chief_all_examp{i_cond}', 'r');
%         plot(tt, chr2_examp{i_cond}', 'b');
%         plot(tt, chronos_examp{i_cond}', 'g');
        plot(tt, mean(chief_cit_examp{i_cond}, 1), 'r', 'linewidth', 4);
        plot(tt, mean(chief_flx_examp{i_cond}, 1), 'm', 'linewidth', 4);
%        plot(tt, mean(chief_all_examp{i_cond}, 1), 'r', 'linewidth', 4);
        plot(tt, mean(chr2_examp{i_cond}, 1), 'b', 'linewidth', 4);
        plot(tt, mean(chronos_examp{i_cond}, 1), 'g', 'linewidth', 4);
    else
        % chief cit 
        xbar = mean(chief_cit_examp{i_cond}, 1);
        sem = stderr(chief_cit_examp{i_cond}, 1);
        shadedErrorBar(tt, xbar, sem, {'r', 'linewidth', 3}, 1);
        % chief flx 
        xbar = mean(chief_flx_examp{i_cond}, 1);
        sem = stderr(chief_flx_examp{i_cond}, 1);
        shadedErrorBar(tt, xbar, sem, {'m', 'linewidth', 3}, 1);
%         % chief all
%         xbar = mean(chief_all_examp{i_cond}, 1);
%         sem = stderr(chief_all_examp{i_cond}, 1);
%         shadedErrorBar(tt, xbar, sem, {'r', 'linewidth', 3}, 1);
        % chronos 
        xbar = mean(chronos_examp{i_cond}, 1);
        sem = stderr(chronos_examp{i_cond}, 1);
        shadedErrorBar(tt, xbar, sem, {'g', 'linewidth', 3}, 1);
        % chr2 
        xbar = mean(chr2_examp{i_cond}, 1);
        sem = stderr(chr2_examp{i_cond}, 1);
        shadedErrorBar(tt, xbar, sem, {'b', 'linewidth', 3}, 1);
    end
    xlabel('time (msec)')
    ylabel('Normalized amplitude')
    axis tight
    switch conds{i_cond}
        case 'FV_Na'
            mytitle('Average FV for first pulse')
        case 'nbqx_apv_cd2_ttx'
            mytitle('Average opsin current for first pulse')
        case 'nbqx_apv_ttx'
            mytitle('Average intracellular current for first pulse')
    end
    
    
    
    if strcmpi(conds{i_cond}, 'FV_Na')
        figure, hold on,
        %plot(tt, cumsum(-mean(chief_cit_examp{i_cond}, 1)), 'r', 'linewidth', 4);
        %plot(tt, cumsum(-mean(chief_flx_examp{i_cond}, 1)), 'm', 'linewidth', 4);
        plot(tt, cumsum(-mean(chr2_examp{i_cond}, 1)), 'b', 'linewidth', 4);
        plot(tt, cumsum(-mean(chronos_examp{i_cond}, 1)), 'g', 'linewidth', 4);
        xlabel('time (msec)')
        ylabel('Normalized LFP amplitude')
        mytitle('Integral of (negative) FV pulse')
        axis tight
    end
    
end


%% SHAPE OF WAVEFORMS FOR EACH PULSE (NORMALIZED AND UN-NORMALIZED)

% choose a stimulation location:
NORMVALS = false;
PLOTERR = true;
STIMSITE = true;
ENFORCETFS = true; % culls some expts from Chronos that are not interleaved
Npulses = 7;

clc; close all

% prepare the population structure
if COMBINE_CHIEF
    opsinTypes = {'chr2', 'chief_all', 'chronos'};
else
    opsinTypes = {'chr2', 'chief_cit', 'chronos', 'chief_flx'};
end

tfnames = {'tf10_led', 'tf20_led', 'tf40_led', 'tf60_led', 'tf100_led'};
switch EXPTTYPE
    case '4-AP'
        drugconds = {'ttx_cd2', 'ttx_cd2_4AP'};
    case 'Baclofen'
        drugconds = {'ttx_cd2', 'ttx_cd2_bac10'};
    case 'Intracellular'
        drugconds = {'nbqx_apv_ttx'};
    case 'main expt'
        drugconds = {'nbqx_apv_cd2_ttx', 'FV_Na'};
    otherwise
        error('This analysis is not supported for the specified EXPTTYPE')
end

for i_opsin = 1:numel(opsinTypes)
     opsincurrent.(opsinTypes{i_opsin}) = repmat({[]}, numel(tfnames), Npulses, numel(drugconds));
     fvlatency.(opsinTypes{i_opsin}) = repmat({[]}, numel(tfnames), Npulses);
     opsin_amps.(opsinTypes{i_opsin}) = repmat({[]}, numel(tfnames), Npulses);
end


% organize and collect the data
Nexpt = size(in,1);
for a_ex = 1:Nexpt
    
    info{a_ex}.mouse;
    opsin = lower(info{a_ex}.opsin);
    
    % what TFs are present in this file?
    ex_tfs = fieldnames(dat{a_ex});
    has10 = any(strcmpi(ex_tfs, 'tf10_led'));
    has20 = any(strcmpi(ex_tfs, 'tf20_led'));
    has40 = any(strcmpi(ex_tfs, 'tf40_led'));
    has60 = any(strcmpi(ex_tfs, 'tf60_led'));
    has100 = any(strcmpi(ex_tfs, 'tf100_led'));
    
    if ENFORCETFS
        if ~(has40 && has60) && strcmpi(opsin, 'chronos')
            continue
        end
    end
    
    for a_ptype = 1:numel(tfnames)
        
        % this pulse type might not exist. Detect these cases and continue
        tmp_tfcond = tfnames{a_ptype};
        if ~isfield(dat{a_ex}, tmp_tfcond)
            continue
        end
        
        % figure out which (if any) channels should be analyzed
        CHANNEL = info{a_ex}.stimSite;
        if isnan(CHANNEL) % cases where neither recording sites was targeted
            continue
        end
        Nchannels = sum(info{a_ex}.ignoreChans);
        if Nchannels == 1
            if ~STIMSITE
                continue % no other stim site to show...
            elseif CHANNEL == 2 % cases where CH = 2, but only one channel recorded
                CHANNEL = 1;
            end
        else
            if STIMSITE
                CHANNEL = info{a_ex}.stimSite;
            else
                if info{a_ex}.stimSite == 1
                    CHANNEL = 2;
                elseif info{a_ex}.stimSite == 2
                    CHANNEL = 1;
                end
            end
        end
        
        
        % iterate over the pharmacology conditions
        for a_drug = 1:numel(drugconds)
            
            tmp_drugcond = drugconds{a_drug};
            
            if ~isfield(dat{a_ex}.(tmp_tfcond).snips, tmp_drugcond);
                continue
            end
            
            % pull out the relevant data, normalize, and add to the population
            % structure
            snips = dat{a_ex}.(tmp_tfcond).snips.(tmp_drugcond){CHANNEL};
            exptSampRate = info{a_ex}.(tmp_tfcond).(tmp_drugcond).sampRate;
            diffvals = [];
            if ~NORMVALS
                normsnips = snips;
                
                %keep track of the raw diff vals, but only for the opsin
                %currents
                if any(strcmpi(tmp_drugcond, {'nbqx_apv_cd2_ttx', 'nbqx_apv_ttx'}))
                    diffvals = dat{a_ex}.(tmp_tfcond).stats.(tmp_drugcond).diffval{CHANNEL};
                    diffvals = abs(diffvals);                    
                end
                
            elseif NORMVALS
                switch tmp_drugcond
                    case {'nbqx_apv_cd2_ttx', 'nbqx_apv_ttx'}
                        troughidx_col = dat{a_ex}.(tmp_tfcond).stats.(tmp_drugcond).trpk_inds{CHANNEL}(:,1);
                        idx_row = 1:size(snips,1);
                        trough_idx = sub2ind(size(snips), idx_row(:), troughidx_col(:));
                        normvals = snips(trough_idx);
                        normvals = abs(normvals); % maintains the sign of the opsin current
                        normsnips = bsxfun(@rdivide, snips, normvals);
                    case 'FV_Na'
                        troughidx_col = dat{a_ex}.(tmp_tfcond).stats.(tmp_drugcond).trpk_inds{CHANNEL}(:,1);
                        peakidx_col = dat{a_ex}.(tmp_tfcond).stats.(tmp_drugcond).trpk_inds{CHANNEL}(:,2);
                        idx_row = 1:size(snips,1);
                        trough_idx = sub2ind(size(snips), idx_row(:), troughidx_col(:));
                        peak_idx = sub2ind(size(snips), idx_row(:), peakidx_col(:));
                        trough_vals = snips(trough_idx);
                        peak_vals = snips(peak_idx);
                        normvals = peak_vals - trough_vals;
                        normsnips = bsxfun(@rdivide, snips, normvals(:));
                end
                
                % estimate the latency, but only for the FV
                pWidth = info{a_ex}.(tmp_tfcond).(tmp_drugcond).pWidth;
                if strcmpi(opsin, 'chronos')
                    preSamps = ceil((prePulseTime+pWidth/2).*exptSampRate);
                else
                    preSamps = ceil((prePulseTime+pWidth).*exptSampRate);
                end
                tt = [0 : size(normsnips,2)-1] ./ exptSampRate;
                tt = tt - prePulseTime;
                if strcmpi(tmp_drugcond, 'FV_Na')
                    tmp_snips = normsnips(:,preSamps:end);
                    thresh  = normsnips(trough_idx) .* 0.85;
                    aboveThresh = bsxfun(@le, tmp_snips, thresh);
                    cross_down = [nan(size(tmp_snips,1),1), diff(aboveThresh, 1, 2)==1]; % rember that the wf is negative...
                    cross_down = mat2cell(cross_down, ones(size(tmp_snips,1),1), size(tmp_snips,2));
                    latency_idx = cellfun(@(x) find(x==1, 1, 'first'), cross_down);
                    latency_idx = latency_idx + preSamps;
                    latency_tt = tt(latency_idx);
                end
                
            end
            
           % correct for differences in the sampling rate
            correctSR = exptSampRate == 20e3;
            if ~correctSR
                
                old_tt = [0 : size(normsnips,2)-1] ./ exptSampRate;
                
                totalTime = size(normsnips,2) ./ exptSampRate;
                newSampRate = 20e3;
                newNumSamps = ceil(totalTime .* newSampRate);
                new_tt = [0 : newNumSamps-1] ./ newSampRate;
                normsnips = interp1(old_tt(:), normsnips', new_tt(:));
                normsnips = normsnips';
                
            end
            
            
            Npulses_ex = size(snips, 1);
            for a_pulse = 1 : min([Npulses, Npulses_ex])
                
                % store the waveform
                tmp = opsincurrent.(opsin){a_ptype, a_pulse, a_drug};
                tmp = cat(1, tmp, normsnips(a_pulse,:));
                opsincurrent.(opsin){a_ptype, a_pulse, a_drug} = tmp;
                
                % store the FV latencies
                if strcmpi(tmp_drugcond, 'FV_Na') && NORMVALS
                    tmp = fvlatency.(opsin){a_ptype, a_pulse};
                    tmp = cat(1, tmp, latency_tt(a_pulse));
                    fvlatency.(opsin){a_ptype, a_pulse} = tmp;
                end
                
                % store the opsin current diff vals
                if any(strcmpi(tmp_drugcond, {'nbqx_apv_cd2_ttx', 'nbqx_apv_ttx'})) && ~NORMVALS
                    tmp = opsin_amps.(opsin){a_ptype, a_pulse};
                    tmp = cat(1, tmp, diffvals(a_pulse));
                    opsin_amps.(opsin){a_ptype, a_pulse} = tmp;
                end
                
            end
        end
        
    end
    
end



%
% PLOT THE WAVEFORMS
%
cmap = gray(Npulses+3);
cmap = cmap(1:Npulses, :);
for a_opsin = 1:numel(opsinTypes)
    
    % first, determine if there are any data for this drug/opsin combo
    tmpdat = opsincurrent.(opsinTypes{a_opsin});
    nodata = cellfun(@(x) isempty(x), tmpdat(:));
    if all(nodata==1)
        continue
    end
    
    f = figure;
    set(f, 'position', [8 41 1135 730], 'name', opsinTypes{a_opsin})
    
    Ndrugs = numel(drugconds);
    for a_drug = 1:Ndrugs
        
        Ntypes = numel(tfnames);
        for a_ptype = 1:Ntypes
            
            %grab the data
            tmpdat = opsincurrent.(opsinTypes{a_opsin})(a_ptype, :, a_drug);
            
            pltidx = Ntypes*(a_drug-1) + a_ptype;
            subplot(Ndrugs, Ntypes, pltidx, 'align')
            hold on,
            
            if a_drug == 1
                title(tfnames{a_ptype})
            end
            if a_ptype == 1
                h = ylabel(drugconds{a_drug});
                h.Interpreter = 'none';
            end
            if a_drug == Ndrugs
                xlabel('Time (ms)')
            end
            
            
            % iterate over the pulses
            for a_pulse = 1:Npulses
                
                xbar = mean(tmpdat{a_pulse}, 1);
                N = sqrt(size(tmpdat{a_pulse}, 1));
                sem = std(tmpdat{a_pulse}, [], 1) ./ sqrt(N);
                if N==1; sem = nan(size(sem)); end
                
                tt = [0:numel(xbar)-1] ./ 20e3;
                tt = tt - prePulseTime;
                tt = tt.*1000; % in ms
                
                if isempty(xbar)
                    continue
                end
                
                pltclr = cmap(a_pulse,:);
                
                % plot the SEM
                if PLOTERR
                    transparent = true;
                    shadedErrorBar(tt, xbar, sem, {'-','color', pltclr, 'linewidth', 2}, transparent);
                else
                    plot(tt, xbar, '-', 'color', pltclr, 'linewidth', 2);
                end
                
                
            end
            axis tight
        end
    end
end


%
% Plot the latencies
%
if NORMVALS
    cmap = copper(numel(tfnames));
    f=figure;
    f.Position = [99         370        1072         307];
    for i_opsin = 1:numel(opsinTypes)
        
        tmpdat = fvlatency.(opsinTypes{i_opsin});
        xbar = cellfun(@mean, tmpdat);
        sem = cellfun(@stderr, tmpdat);
        
        subplot(1, numel(opsinTypes), i_opsin)
        hold on,
        for i_tf = 1:numel(tfnames);
            my_errorbar(1:numel(xbar(i_tf,:)), xbar(i_tf,:), sem(i_tf,:), 'color', cmap(i_tf,:));
        end
        ylim([0, 2.2e-3])
        xlim([0.5, 7.5])
        title(opsinTypes{i_opsin})
        legend(tfnames, 'location', 'best')
        legend boxoff
        xlabel('pulse number')
        ylabel('FV latency (sec)')
    end
end


%
% Plot the raw opsin currents to compare P1 and P7
%
if ~NORMVALS
    f=figure;
    hold on;
    pltclr = {'b', 'r', 'g', 'm'};
    
    for i_opsin = 1:numel(opsinTypes)
        
        opsinTypes{i_opsin}
        tmpdat = opsin_amps.(opsinTypes{i_opsin});
        xbar = cellfun(@mean, tmpdat)
        sem = cellfun(@stderr, tmpdat)
        
        % plot the first pulse
        x = (1:numel(tfnames)) + (0.2*i_opsin-1);
        my_errorbar(x, xbar(:,1)', sem(:,1)', 's', 'color', pltclr{i_opsin}, 'linewidth', 2, 'markersize', 10, 'markerfacecolor', pltclr{i_opsin});
        
        % plot the last pulse
        x2 = x + 0.1;
        my_errorbar(x2, xbar(:,7)', sem(:,7)', 's', 'color', pltclr{i_opsin}, 'linewidth', 2, 'markersize', 10);
        
        % plot the line segments
        plot([x;x2], xbar(:, [1,7])', 'color', pltclr{i_opsin}, 'linewidth', 2)
        
    end
    ylim([0 0.09])
    set(gcf, 'Position', [13 430 1420 325]);
    set(gca, 'Position', [0.0423 0.1100 0.9352 0.8150]);
end


%% AVERAGE RECORDING WAVEFORM OF OPSIN CURRENT FOR WHOLE SWEEP

% choose a stimulation location:
STIMSITE = true;
PLOTERR = true;
NORMALIZE = false;


% only do this for the main experiment (TF and FV)
correctExpt = any(strcmpi(EXPTTYPE, {'main expt', '4-AP', 'baclofen', 'intracellular'}));
assert(correctExpt, 'ERROR: This analysis is for a different type of experiment');

% prepare the population structure
if COMBINE_CHIEF
    opsinTypes = {'chr2', 'chief_all', 'chronos'};
else
    opsinTypes = {'chr2', 'chief_cit', 'chronos', 'chief_flx'};
end
tfnames = {'tf10_led', 'tf20_led', 'tf40_led', 'tf60_led', 'tf100_led'};
avgopsincurent = [];
for i_opsin = 1:numel(opsinTypes)
    for i_tf = 1:numel(tfnames)
        avgopsincurent.(opsinTypes{i_opsin}).(tfnames{i_tf}).raw = [];
    end
end


Nexpt = size(in,1);
for a_ex = 1:Nexpt
    opsin = lower(info{a_ex}.opsin);
    
    pTypes = fieldnames(dat{a_ex});
    
    for a_ptype = 1:numel(pTypes)
        
        % find the TTX pharm condition
        conds = fieldnames(dat{a_ex}.(pTypes{a_ptype}).snips);
        ttxidx = cellfun(@any, regexpi(conds, 'ttx'));
        
        if ~any(ttxidx);
            continue
        end
        
        
        ttxidx = find(ttxidx);
        for i_ttx = ttxidx'; % so that I can compare multiple K+ drugs
            
            
            % figure out which (if any) channels should be analyzed
            CHANNEL = info{a_ex}.stimSite;
            if isnan(CHANNEL) % cases where neither recording sites was targeted
                continue
            end
            Nchannels = sum(info{a_ex}.ignoreChans);
            if Nchannels == 1
                if ~STIMSITE
                    continue % no other stim site to show...
                elseif CHANNEL == 2 % cases where CH = 2, but only one channel recorded
                    CHANNEL = 1;
                end
            else
                if STIMSITE
                    CHANNEL = info{a_ex}.stimSite;
                else
                    if info{a_ex}.stimSite == 1
                        CHANNEL = 2;
                    elseif info{a_ex}.stimSite == 2
                        CHANNEL = 1;
                    end
                end
            end
            
            % pull out the relevant data, normalize, and add to the population
            % structure
            raw = dat{a_ex}.(pTypes{a_ptype}).(conds{i_ttx})(:,CHANNEL);
            if NORMALIZE
                pulse1_snip = dat{a_ex}.(pTypes{a_ptype}).snips.(conds{i_ttx}){CHANNEL}(1,:);
                troughidx = dat{a_ex}.(pTypes{a_ptype}).stats.(conds{i_ttx}).trpk_inds{CHANNEL}(1,1);
                normval = pulse1_snip(troughidx);
                normval = abs(normval); % this preserves the sign of the original signal
                raw = raw ./ normval;
            end
            
            % figure out how many pulses there were, cut off everything past
            % the 7th pulse. If not 7 pulses, move along
            sampRate = info{a_ex}.(pTypes{a_ptype}).(conds{i_ttx}).sampRate;
            pulseOnIdx = find(info{a_ex}.(pTypes{a_ptype}).(conds{i_ttx}).pulseOn_idx);
            ipi_samps = unique(round(diff(pulseOnIdx)));
            pulseOnIdx = pulseOnIdx(1);
            pulseOffIdx = find(info{a_ex}.(pTypes{a_ptype}).(conds{i_ttx}).pulseOff_idx);
            if numel(pulseOffIdx)<7; continue; end
            pulseOffIdx = pulseOffIdx(7);
            raw = raw(pulseOnIdx - ipi_samps : pulseOffIdx+ipi_samps-5);
            
            % correct for differences in the sampling rate
            correctSR = sampRate == 20e3;
            if ~correctSR
                oldSampRate = sampRate;
                old_tt = [0 : numel(raw)-1] ./ oldSampRate;
                
                totalTime = numel(raw) ./ oldSampRate;
                newSampRate = 20e3;
                newNumSamps = ceil(totalTime .* newSampRate);
                new_tt = [0 : newNumSamps-1] ./ newSampRate;
                raw = interp1(old_tt(:), raw(:), new_tt(:));
            end
            
            
            % add the raw data to the array. Make sure that there is a spot
            % to put the data first
            
            if ~isfield(avgopsincurent.(opsin).(pTypes{a_ptype}).raw, conds{i_ttx})
                avgopsincurent.(opsin).(pTypes{a_ptype}).raw.(conds{i_ttx}) = [];
            end
            tmp =  avgopsincurent.(opsin).(pTypes{a_ptype}).raw.(conds{i_ttx});
                
            raw = raw(:)';
            if ~isempty(tmp)
                raw(size(tmp,2):end) = [];
                raw(end:size(tmp,2)) = nan;
            end
            
            avgopsincurent.(opsin).(pTypes{a_ptype}).raw.(conds{i_ttx}) = cat(1, tmp, raw);
            
        end
    end
    
end


for a_tf = 1:numel(tfnames)

    f = figure; hold on,
    f.Name = tfnames{a_tf};
    
    legtext = {};
    for a_opsin = 1:numel(opsinTypes)
        
        
        
        if isempty(avgopsincurent.(opsinTypes{a_opsin}).(tfnames{a_tf}).raw);
            continue;
        end
        
        drugconds = fieldnames(avgopsincurent.(opsinTypes{a_opsin}).(tfnames{a_tf}).raw);
        
        switch opsinTypes{a_opsin}
            case 'chr2'
                pltclr = [0 0 1];
            case {'chief_cit', 'chief_all'}
                pltclr = [1 0 0];
            case 'chief_flx'
                pltclr = [1 0 .5];
            case 'chronos'
                pltclr = [0 1 0];
        end
        
        pltclr = repmat(pltclr, numel(drugconds), 1);
        pltclr = bsxfun(@times, pltclr, linspace(0.5, 1, numel(drugconds))');
        for a_drug = 1:numel(drugconds)
            
            tmpdat = avgopsincurent.(opsinTypes{a_opsin}).(tfnames{a_tf}).raw.(drugconds{a_drug});
            tmpmean = nanmean(tmpdat, 1);
            N = size(tmpdat,1);
            if N>1
            tmpsem = nanstd(tmpdat, [], 1) ./ sqrt(N);
            else
                tmpsem = nan(size(tmpdat));
            end
            Nsamps = numel(tmpmean);
            tt = ([0:Nsamps-1]-ipi_samps) ./ 20e3;
            
            if PLOTERR
                shadedErrorBar(tt, tmpmean, tmpsem, {'-','color', pltclr(a_drug,:), 'linewidth', 2}, 1);
            else
                plot(tt, tmpdat', '-', 'color', pltclr(a_drug,:), 'linewidth', 2)
            end
            
            tmp = [opsinTypes{a_opsin}, ': ', drugconds{a_drug}];
            legtext = cat(2, legtext, tmp);
            
        end
    end
    if ~isempty(legtext)
        L = legend(legtext);
        L.Location = 'best';
        L.Box = 'off';
        L.Interpreter = 'none';
    end
end



%% TIME CONSTANT ANALYSIS (ONE VS. TWO)

close all

STIMSITE = true;
PLOTALL = true;


tau_prePulseTime = 0; % in sec
tau_postPulseTime = 0.023; % in sec, chosen so that I can average 10 and 20 Hz sweeps.
Nexpts = size(in,1);
poptau = {};
for i_ex = 1:Nexpts

    pTypes = fieldnames(dat{i_ex});
    Ntfs = numel(pTypes);
    tmp_snippet = []; % concatenate across TFs for each channel
    for i_tf = 1:Ntfs
        
        % only allow 10 and 20 Hz to pass through
        if ~any(strcmpi(pTypes{i_tf}, {'tf10_led', 'tf20_led', 'tf_40_led'}))
            continue
        end
        
        conds = fieldnames(dat{i_ex}.(pTypes{i_tf}));
        ttx_idx = cellfun(@(x) ~isempty(regexpi(x, 'ttx')), conds);
        assert(sum(ttx_idx)<=1, 'ERROR: too many ttx conditions')
        if ~any(ttx_idx)
            continue
        end
        i_cond = find(ttx_idx);
        
        sampRate = info{i_ex}.(pTypes{i_tf}).(conds{i_cond}).sampRate;
        prePulseSamps = ceil(tau_prePulseTime .* sampRate);
        postPulseSamps = ceil(tau_postPulseTime .* sampRate);
        pulseOn_idx = find(info{i_ex}.(pTypes{i_tf}).(conds{i_cond}).pulseOn_idx);
        
        i_pulse = 1;
        snip_idx = pulseOn_idx(i_pulse)-prePulseSamps : 1 : pulseOn_idx(i_pulse)+postPulseSamps;
        
        % figure out which (if any) channels should be analyzed
        CHANNEL = info{i_ex}.stimSite;
        if isnan(CHANNEL) % cases where neither recording sites were targeted
            continue
        end
        Nchannels = sum(info{i_ex}.ignoreChans);
        if Nchannels == 1
            if ~STIMSITE
                continue % no other stim site to show...
            elseif CHANNEL == 2 % cases where CH = 2, but only one channel recorded
                CHANNEL = 1;
            end
        else
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
        
        % pull out the snippet, subtract off the baseline and
        % store it for each pulse in the train
        snippet_full = dat{i_ex}.(pTypes{i_tf}).(conds{i_cond})(snip_idx , CHANNEL);
        basetime = round(0.001 .* sampRate);
        baseline = mean(snippet_full(end-basetime : end));
        snippet_full = snippet_full - baseline;
        
        tmp_snippet = cat(1, tmp_snippet, snippet_full');
        
    end
    
    if isempty(tmp_snippet)
        poptau{i_ex}.snips = []; % in case there are no data
        continue
    end
    
    
    % average across 10, 20 Hz
    tmp_snippet = mean(tmp_snippet, 1);
    tt = [0:numel(tmp_snippet)-1] ./ sampRate;
    poptau{i_ex}.snips = tmp_snippet;
    
    % grab the subset of the data to fit
    [troughidx, peakidx]  = anlyMod_getWFepochs(tmp_snippet, tt, conds{i_cond}, 300e-6, 200e-6, 'inward');
    startVal = tmp_snippet(troughidx) .* 0.95;
    startIdx = find((tmp_snippet > startVal) & (tt > tt(troughidx)), 1, 'first');
    fit_tt = tt(startIdx : end);
    fit_dat = tmp_snippet(startIdx : end);
    
    %inital guesses, and 1exp fit
    N_2ms = round(sampRate .* 0.002);
    pred = [ones(N_2ms,1), fit_tt(1:N_2ms)'];
    resp = fit_dat(1:N_2ms)';
    guess = pred \ resp;
    guess(1) = exp(guess(1)); % to solve for amplitude of exponential
    fitopt = fitoptions('exp1');
    fitopt.StartPoint = guess;
    fitopt.Lower = [0, -inf]; % constrain amps to be +, Taus to be -.
    fitopt.Upper = [inf, 0];
    fitopt.TolFun = 1e-10;
    fitopt.TolX = 1e-10;
    [fexp1, gof1] = fit(fit_tt(:), -fit_dat(:), 'exp1', fitopt);
    poptau{i_ex}.tau_1 = 1./fexp1.b;
    poptau{i_ex}.amp_1 = fexp1.a;
    poptau{i_ex}.gof_1 = gof1.rsquare;
    
    % fit a bi-exponential
    fitopt = fitoptions('exp2');
    fitopt.StartPoint = [guess(1)/2, guess(1)*10, guess(1)/2, guess(1)/5];
    fitopt.Lower = [0, -inf, 0, -inf]; % constrain amps to be +, Taus to be -.
    fitopt.Upper = [inf, 0, inf, 0];
    fitopt.TolFun = 1e-10;
    fitopt.TolX = 1e-10;
    [fexp2, gof2] = fit(fit_tt(:), -fit_dat(:), 'exp2', fitopt);
    bothTaus = [1./fexp2.b, 1./fexp2.d];
    bothAmps = [fexp2.a, fexp2.c];
    [~, idx] = sort(abs(bothTaus)); %sort by faster tau first
    poptau{i_ex}.tau_2 = bothTaus(idx);
    poptau{i_ex}.amp_2 = bothAmps(idx);
    poptau{i_ex}.gof_2 = gof2.rsquare;
    
    if PLOTALL
        figure, hold on,
        set(gcf, 'name', num2str(i_ex))
        plot(tt, tmp_snippet);
        fit1 = -1.*(fexp1.a .* exp(fexp1.b .* fit_tt));
        plot(fit_tt, fit1, 'r')
        fit2 = -1.* ( (fexp2.a .* exp(fexp2.b .* fit_tt)) + (fexp2.c .* exp(fexp2.d .* fit_tt)) );
        plot(fit_tt, fit2, 'g')
        legend({'raw', 'exp1', 'exp2'}, 'location', 'best')
        ylabel(info{i_ex}.opsin)
        title(sprintf('gof1: %.2f, gof2: %.2f', gof1.rsquare, gof2.rsquare))
        drawnow
    end
    
    
end

if COMBINE_CHIEF
    opsins = {'chr2', 'chronos', 'chief_all'};
else
    opsins = {'chr2', 'chronos', 'chief_cit', 'chief_flx'};
end

opsintau = [];
for i_opsin = opsins
    opsintau.(i_opsin{1}).tau1 = [];
    opsintau.(i_opsin{1}).gof1 = [];
    opsintau.(i_opsin{1}).gof2 = [];
    opsintau.(i_opsin{1}).amp2 = [];
end

for i_ex = 1:Nexpts
    tmpopsin = lower(info{i_ex}.opsin);
    if ~isempty(poptau{i_ex}.snips) % data were fitted
        opsintau.(tmpopsin).tau1 = cat(1, opsintau.(tmpopsin).tau1, poptau{i_ex}.tau_1);
        opsintau.(tmpopsin).gof1 = cat(1, opsintau.(tmpopsin).gof1, poptau{i_ex}.gof_1);
        opsintau.(tmpopsin).gof2 = cat(1, opsintau.(tmpopsin).gof2, poptau{i_ex}.gof_2);
        opsintau.(tmpopsin).amp2 = cat(1, opsintau.(tmpopsin).amp2, poptau{i_ex}.amp_2);
    end
end


% plot of single exp time constants
figure, hold on
pltclr = {'b', 'g', 'r', 'm'};
for i_op=1:numel(opsins)
    if ~isempty(opsintau.(opsins{i_op}).tau1)
        plot(i_op, -opsintau.(opsins{i_op}).tau1*1000, 'o', 'color', pltclr{i_op})
        plot(i_op, mean(-opsintau.(opsins{i_op}).tau1*1000), 's', 'color', pltclr{i_op}, 'markerfacecolor', pltclr{i_op})
    end
end
ylabel('Time constant (ms)')
xlim([0 numel(opsins)+1])
set(gca, 'xtick', 1:numel(opsins), 'xticklabel', opsins)


% plot of R-squared values for single/double exp fits
figure, hold on
for i_op = 1:numel(opsins)
    if ~isempty(opsintau.(opsins{i_op}).gof2)
        plot([i_op-.2, i_op+.2], [opsintau.(opsins{i_op}).gof1(:),opsintau.(opsins{i_op}).gof2(:)], '-o', 'color', pltclr{i_op})
    end
end
ylabel('R-squared')
xlim([0 numel(opsins)+1])
title('R-square of 1exp on left')
set(gca, 'xtick', 1:numel(opsins), 'xticklabel', opsins)


% plot of Amplitudes for each exp in the DOUBLE exp fit
figure, hold on
for i_op = 1:numel(opsins)
    if ~isempty(opsintau.(opsins{i_op}).amp2)
        plot([i_op-.2, i_op+.2], opsintau.(opsins{i_op}).amp2, '-o', 'color', pltclr{i_op})
    end
end
ylabel('Amplitudes')
xlim([0 numel(opsins)+1])
set(gca, 'xtick', 1:numel(opsins), 'xticklabel', opsins)
title('Fast Tau on left')




%% RUNDOWN: RUNNING AVERAGE OF STATS ACROSS SWEEPS

SWEEPSTOAVG = 7;
DEBUG = false;



% grab the fiber volley pop excel workbook
fname = [GL_DOCUPATH, 'Other_workbooks', filesep, 'fiberVolleyCellList.xlsx'];
[~,txt, raw] = xlsread(fname);
raw(size(txt,1)+1:end, :) = [];
raw(:,size(txt,2)+1:end) = [];
channelIdx = cellfun(@(x) ~isempty(x), regexpi(raw(1,:), 'CH\d'));
opsinIdx = strcmpi(raw(1,:), 'opsin');
primaryStimSiteIdx = strcmpi(raw(1,:), 'Primary Stim Site (0,0)');
channels = raw(:, cellfun(@(x) ~isempty(x), regexpi(raw(1,:), 'CH\d')));
rmsweeps = raw(:, strcmpi(raw(1,:), 'rmSweeps'));


% replace the 'file names' with fully qualified paths. I'm hoping that this
% allows for parallel operations
Nfiles = size(raw,1)-1;
idx_mousename = strcmpi(raw(1,:), 'Mouse Name');
idx_fname =  strcmpi(raw(1,:), 'file name');
tmp_fnames = raw(2:end, idx_fname);
tmp_mousenames = raw(2:end, idx_mousename);
prefix = cellfun(@(x,y) [x,y], repmat({GL_DATPATH}, Nfiles, 1), tmp_mousenames, 'uniformoutput', false);
prefix = cellfun(@(x,y) [x,y], prefix, repmat({[filesep, 'Physiology', filesep]}, Nfiles, 1), 'uniformoutput', false);
tmp_fnames = cellfun(@(x,y) [x,y], prefix, tmp_fnames, 'uniformoutput', false);
tmp_fnames = cellfun(@(x,y) [x,y], tmp_fnames, repmat({'.abf'}, Nfiles, 1), 'uniformoutput', false);
raw(2:end, idx_fname) = tmp_fnames;


% do the analysis
smoothStats = {};
Nexpts = size(in,1);
conds = {'ttx', 'ttx_cd2', 'nbqx_apv', 'nbqx_apv_cd2', 'nbqx_apv_cd2_ttx'};
for i_ex = 1:Nexpts;
    
    % update the user with what's going on
    fprintf('Analyzing mouse %s site %d, file %d of %d\n', in{i_ex, 1}, in{i_ex, 2}, i_ex, Nexpts)
    
    for i_cond = 1:numel(conds)
        
        % figure out what rows in the work book to pay attention to,
        % which will define a particular file
        l_mouse = cellfun(@(x) ~isempty(x), regexp(raw(:,1), in{i_ex,1}));
        l_site = [false ; cell2mat(raw(2:end,2)) == in{i_ex,2}]; % add a leading 'false' to account for the header row in 'raw'
        l_cond = strcmpi(raw(:,5), regexprep(conds{i_cond}, '_', '\+'));
        l_expt = l_mouse & l_site & l_cond;
        
        % only analyze files where the conditions were interleaved.
        % otherwise the ordinal relationship between trials will be wrong
        interleaved = cellfun(@(x) strcmpi(x, 'interleaved'), raw(:,4));
        if any(~interleaved(l_expt))
            continue
        end
        
        % simple error check
        if sum(l_expt) == 0
            continue
        else
            assert(sum(l_expt) <= 1, 'ERROR: too many data files')
        end
        
        % add a few useful peices of information to the data structure
        smoothStats{i_ex}.opsin = raw{l_expt, opsinIdx};
        smoothStats{i_ex}.stimSite = raw{l_expt, primaryStimSiteIdx};
        smoothStats{i_ex}.mouseName = in{i_ex, 1};
        smoothStats{i_ex}.siteNum = in{i_ex, 2};
        
        
        % extract the raw data
        fname = raw{l_expt, 3};
        ax = abfobj(fname);
        ax.head.validChans = cat(1, channels{l_expt,:});
        
        % remove some sweeps if need be
        if ~isnan(rmsweeps{l_expt})
            goodSweeps = true(size(ax.dat,3),1);
            badSweeps = eval(['[',rmsweeps{l_expt},']']);
            goodSweeps(badSweeps) = false;
            ax.dat = ax.dat(:,:,goodSweeps);
            ax.wf = ax.wf(:,:,goodSweeps);
        end
        
        % store the data, grouped by unique conditions. First, I need to figure
        % out the appropriate field name for each condition
        ledChIdx = strcmpi('LED_470', ax.head.recChNames);
        if ~any(ledChIdx)
            ledChIdx = strcmpi('Laser', ax.head.recChNames);
        end
        
        tmpWF = ax.dat(:, ledChIdx, :);
        sampRate = ax.head.sampRate;
        tdict = outerleave(tmpWF, sampRate);
        
        
        % make sure the first LED pulse is identical across trial types
        assert(numel(unique(tdict.conds(:,1)))==1, 'ERROR: too many pAmps')
        assert(all(diff(tdict.conds(:,2)) < 1e-14), 'ERROR: too many pWidths')
        pwidth = tdict.conds(1,2);
        
        % store some useful parameters
        sampRate = ax.head.sampRate;
        prePulseTime = 0.006; % in sec
        postPulseTime = 0.009; % in sec
        prePulseSamps = ceil(prePulseTime .* sampRate); % samples prior to pulse onset
        postPulseSamps = ceil(postPulseTime .* sampRate); % samples available after pulse ONSET
        photoDelay= 150e-6; % timeout following pulse offset
        
        
        % for each channel, run through sweeps and calculate the running
        % average raw data trace. Extract stats for the first pulse
        for i_ch = 1:2;
            
            % deal with data and channels that need to be ignored.
            if ~channels{l_expt, i_ch}
                continue
            else
                % occasionally CH1 is disabled, and CH2 occupies the
                % 1st and only column, adjust i_ch to account for these
                % cases
                recChs = cellfun(@(x) ~isempty(regexpi(x, '_Vm')), ax.head.recChNames);
                Nchannels = sum(recChs);
                if (Nchannels == 1) && (i_ch == 2)
                    i_ch = 1;
                end
            end
            
            
            % figure out the appropriate index to the data. It should have
            % the correct "HSx_" prefix, and have the correct units
            validChans = find([channels{l_expt,:}] == 1);
            whichChan = validChans(i_ch);
            correctUnits = strcmpi('mV', ax.head.recChUnits);
            correctHS = strncmpi(sprintf('HS%d_', whichChan), ax.head.recChNames, 3);
            datChIdx = correctUnits & correctHS;
            assert(sum(datChIdx)==1, 'ERROR: incorrect channel selection')

            
            %
            % Extract the first pulse from all the sweeps, do the running
            % average later
            %
            %%%%%%%%%%%%%%%%%%%%
            
            % grab the raw data and the sampled LED waveforms
            lfp_raw = ax.dat(:, datChIdx, :);
            lfp_raw = permute(lfp_raw, [3,1,2]);
            led = ax.dat(:, ledChIdx, :);
            led = permute(led, [3,1,2]);
            
            
            % make sure the LED pulse comes on at the correct time for
            % each sweep
            tcross = diff(led > 0.05, [], 2)== 1;
            tcross = cat(2, false(size(led, 1),1), tcross);
            firstCross = cellfun(@(x) find(x==1, 1, 'first'), mat2cell(tcross, ones(size(led, 1), 1)));
            assert(numel(unique(firstCross))==1, 'ERROR: first pulse times are inconsistent')
            firstCross = firstCross(1);
            
            % extract the first pulse, and subtract off the baseline
            col_idx = (firstCross-prePulseSamps) : (firstCross+postPulseSamps);
            firstPulses = lfp_raw(:, col_idx);
            bkgnd = mean(firstPulses(:, 1:prePulseSamps), 2);
            firstPulses = bsxfun(@minus, firstPulses, bkgnd);
            
            
            %
            %  Here is where I'll pull out the adjcent sweeps. To start,
            %  pull out all snippets and rearrage them acording to trial
            %  order
            %
            Nsweeps = size(firstPulses, 1);
            smoothStats{i_ex}.(conds{i_cond}).trough{i_ch} = nan(Nsweeps,1);
            for i_swp = 1:(Nsweeps-SWEEPSTOAVG)
                
                % make a list of relavent trials (a sliding window)
                tList = i_swp : 1 : (i_swp+SWEEPSTOAVG-1);
                
                % grab the appropriate sweeps
                tmp = mean(firstPulses(tList,:), 1);
                
                lp_freq = 2500;
                tmp = butterfilt(tmp, lp_freq, sampRate, 'low', 2);
                
                % cacluate the peak and trough idx
                tt = ([0:numel(tmp)-1] - prePulseSamps) ./ sampRate; % time=0 is when the LED comes on
                [troughidx, peakidx]  = anlyMod_getWFepochs(tmp, tt, conds{i_cond}, pwidth, photoDelay, 'inward');
                
                % add a few points on either side of the true trough/peak
                trough_window = troughidx-4: min([troughidx+4, numel(tmp)]);
                assert(~isempty(trough_window) && numel(trough_window)==9, 'ERROR: no data for trough window')
                if strcmpi(conds{i_cond}, 'FV_Na')
                    peak_window = peakidx-4: min([peakidx+4, numel(tmp)]);
                    assert(~isempty(peak_window) && numel(peak_window)==9, 'ERROR: no data for peak window')
                end
                
                
                %
                % store some stats for each condition type
                %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                switch conds{i_cond}
                    case {'nbqx_apv_cd2_ttx', 'nbqx_apv_cd2', 'nbqx_apv', 'ttx', 'ttx_cd2'}
                        
                        trough = mean(tmp(trough_window));
                        smoothStats{i_ex}.(conds{i_cond}).trough{i_ch}(i_swp) = trough;
                    case 'synapticTransmission'
                end
                
                
                if DEBUG
                    cla
                    hold on,
                    plot(tt, tmp);
                    plot(0, 0, 'ro')
                    plot(pwidth, 0, 'ro')
                    plot(pwidth+photoDelay, 0, 'ko')
                    plot(tt(troughidx), trough, 'go')
                    drawnow
                end
            end % sweeps
            
            
            % keep track of the correlation between trial number and the
            % diff val
            all_difvals = smoothStats{i_ex}.(conds{i_cond}).trough{i_ch};
            all_difvals(isnan(all_difvals)) = [];
            [r, p] = corr((1:numel(all_difvals))', abs(all_difvals), 'type', 'spearman');
            smoothStats{i_ex}.(conds{i_cond}).rho(i_ch) = r;
            
            
            % store the (average) waveforms for the first/last 20 trials
            tmp = mean(firstPulses(1:20, :), 1);
            tmp = butterfilt(tmp, lp_freq, sampRate, 'low', 2);
            smoothStats{i_ex}.(conds{i_cond}).first20{i_ch} = tmp;
            
            tmp = mean(firstPulses(end-19:end, :), 1);
            tmp = butterfilt(tmp, lp_freq, sampRate, 'low', 2);
            smoothStats{i_ex}.(conds{i_cond}).last20{i_ch} = tmp;
            
        end % channels
        
    end % conditions
end % expts


%%  PLOTTING ROUTINES FOR RUNDOWN ANALYSIS

close all
NORMVALS = true;
STIMSITE = true;

% clear out the structures that have no data
l_empty = cellfun(@isempty, smoothStats);
smoothStats = smoothStats(~l_empty);
catdat = {};

exptOpsins = structcat(smoothStats, 'opsin');

if COMBINE_CHIEF
    opsins = {'chr2', 'chief_all', 'chronos'};
else
    opsins = {'chr2', 'chief_cit', 'chronos'};
end

for i_cond = 1:numel(conds)
    for i_opsin = 1:numel(opsins);
        pooledRho.(opsins{i_opsin}).(conds{i_cond}) = [];
    end
end


for i_opsin = 1:numel(opsins)
    
    l_opsin = strcmpi(exptOpsins, opsins{i_opsin});
    l_opsin = find(l_opsin);
    
    h(i_opsin) = myfig([1 379 1141 405], opsins{i_opsin});
    for i_ex = 1:numel(l_opsin)
        idx = l_opsin(i_ex);
        
        % Determine which recording channel to analyze
        if STIMSITE
            CHANNEL = smoothStats{idx}.stimSite;
        else
            if smoothStats{idx}.stimSite == 1
                CHANNEL = 2;
            elseif smoothStats{idx}.stimSite == 2
                CHANNEL = 1;
            end
        end

        % plot each of the drug conditions
        for i_cond = 1:numel(conds)
            if ~isfield(smoothStats{idx}, conds{i_cond})
                continue
            end
            
            % fix the CHANNEL in the cases where HS1 was not recorded from
            if CHANNEL == 2 && numel(smoothStats{idx}.(conds{i_cond}).trough)==1
                CHANNEL = 1;
                disp('here')
            end
            
            tmp = smoothStats{idx}.(conds{i_cond}).trough{CHANNEL};
            tmp = abs(tmp);
            if NORMVALS
                tmp = (tmp - tmp(1)) ./ tmp(1);
            end
            subplot(1,numel(conds), i_cond)
            hold on,
            plot(tmp, 'k')
            ylabel('diffval')
            xlabel('trial number')
            mytitle(conds{i_cond})
            
            % compile the pooled rho structure for later plotting
            pooledRho.(opsins{i_opsin}).(conds{i_cond}) = cat(1, ...
                pooledRho.(opsins{i_opsin}).(conds{i_cond}),...
                smoothStats{idx}.(conds{i_cond}).rho(CHANNEL));
            
            % concatenate all the relevant data
            catdat{i_opsin,i_cond}{i_ex} = tmp;
        end
        
    end
end




% plot the mean of each noisy thing
maxY = 0;
for i_opsin = 1:numel(opsins)
    figure(h(i_opsin));
    for i_cond = 1:numel(conds)
        subplot(1,numel(conds), i_cond)
        tmp = catdat{i_opsin, i_cond};
        if isempty(tmp)
            continue
        end
        maxsize = max(cellfun(@numel, tmp));
        newdat = nan(numel(tmp), maxsize);
        for i_ex = 1:numel(tmp)
            newdat(i_ex,1:numel(tmp{i_ex})) = tmp{i_ex};
        end
        plot(nanmean(newdat), 'b', 'linewidth', 3)
        axis tight
        ylims = get(gca, 'ylim');
        maxY = max([maxY, abs(ylims)]);
    end
end
% standardize the axes
for i_opsin = 1:numel(opsins)
    figure(h(i_opsin));
    for i_cond = 1:numel(conds)
        subplot(1,numel(conds), i_cond)
        set(gca, 'ylim', [-maxY, maxY]);
    end
end






%
% plot the waveform for the first/last 20 pulses
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i_opsin = 1:numel(opsins)
    f = figure;
    f.Name = opsins{i_opsin};
    f.Position = [1 37 1152 793];
    
    l_opsin = strcmpi(exptOpsins, opsins{i_opsin});
    l_opsin = find(l_opsin);
    Nexpts = numel(l_opsin);
    Nconds = numel(conds);
    
    for i_ex = 1:Nexpts;
        
        idx = l_opsin(i_ex);
        
        % Determine which recording channel to analyze
        % Determine which recording channel to analyze
        if STIMSITE
            CHANNEL = smoothStats{idx}.stimSite;
        else
            if smoothStats{idx}.stimSite == 1
                CHANNEL = 2;
            elseif smoothStats{idx}.stimSite == 2
                CHANNEL = 1;
            end
        end
        
        
        for i_cond = 1:Nconds
            if ~isfield(smoothStats{idx}, conds{i_cond})
                continue
            end
            
            % fix the CHANNEL in the cases where HS1 was not recorded from
            if CHANNEL == 2 && numel(smoothStats{idx}.(conds{i_cond}).trough)==1
                CHANNEL = 1;
                disp('here')
            end
            
            s = subplot(Nexpts, Nconds, (Nconds.*(i_ex-1) + i_cond) );
            hold on,
            plot(tt.*1000, smoothStats{idx}.(conds{i_cond}).first20{CHANNEL}, 'k')
            plot(tt.*1000, smoothStats{idx}.(conds{i_cond}).last20{CHANNEL}, 'b')
            s.XTick = [];
            axis tight
            s.XLim = [-1, 6];
            if i_ex == 1
                mytitle(conds{i_cond})
            end
            drawnow
            ylabel(smoothStats{idx}.mouseName)
            
        end
    end
end



%% RUNDOWN ANALYSIS: 1ST_5 VS. LAST_5

% this script requires the aggregated data on the first/last 5 sweeps
% (saved .mat files), and also requires importing the stats for cd2 and ttx
% conditions separately.

fin

load 'first5.mat';
dat_first = dat;
info_first = info;

load 'last5.mat';
dat_last = dat;
info_last = info;

clear info dat

STIMSITE = true;
DRUGCOND = 'nbqx_apv_cd2'; 
STATTYPE = 'diffval';


[pop_first, pop_last, pop_pamp] = deal({});

% load in the LED calibration data to convert from Volts to power denisty
load led_470_cal.mat % data now stored in 'cal' struct

figure; hold on,
for i_ex = 1:numel(dat_first)
     
    % figure out which (if any) channels should be analyzed
    CHANNEL = info_first{i_ex}.stimSite;
    if isnan(CHANNEL) % cases where neither recording sites were targeted
        continue
    end
    Nchannels = sum(info_first{i_ex}.ignoreChans);
    if Nchannels == 1
        if ~STIMSITE
            continue % no other stim site to show...
        elseif CHANNEL == 2 % cases where CH = 2, but only one channel recorded
            CHANNEL = 1;
        end
    else
        if STIMSITE
            CHANNEL = info_first{i_ex}.stimSite;
        else
            if info_first{i_ex}.stimSite == 1
                CHANNEL = 2;
            elseif info_first{i_ex}.stimSite == 2
                CHANNEL = 1;
            end
        end
    end
    
    
    
    % there could be multiple TFs and pAmps. figure out what the deal is.
    pTypeNames = fieldnames(dat_first{i_ex});
    ex_pAmps = nan(numel(pTypeNames), 1);
    ex_TFs = nan(numel(pTypeNames), 1);
    for i_pType = 1:numel(pTypeNames);
        if ~isfield(info_first{i_ex}.(pTypeNames{i_pType}), DRUGCOND)
            continue
        end
        ex_pAmps(i_pType) = info_first{i_ex}.(pTypeNames{i_pType}).(DRUGCOND).pAmp;
        ex_TFs(i_pType) = info_first{i_ex}.(pTypeNames{i_pType}).(DRUGCOND).pTF;
    end
    
    
    % for each pAmp, average the first pulse stats across available TFs
    pAmps_unique = unique(ex_pAmps);
    [first_avg, last_avg] = deal(nan(1, numel(pAmps_unique)));
    for i_pamp = 1:numel(pAmps_unique);
        
        % find the field names with the TF conds that have a particular
        % pulse amplitude.
        idx = ex_pAmps == pAmps_unique(i_pamp);
        fldnames = pTypeNames(idx);
        
        [first_raw, last_raw] = deal(nan(numel(fldnames), 1));
        for i_fld = 1:numel(fldnames)
            
            if ~isfield(dat_first{i_ex}.(fldnames{i_fld}).stats, DRUGCOND)
                continue
            end
            
            % get the first5 stat
            tmp = dat_first{i_ex}.(fldnames{i_fld}).stats.(DRUGCOND).(STATTYPE){CHANNEL}(1);
            tmp = abs(tmp);
            first_raw(i_fld) = tmp;
            
            % get the last5 stat
            tmp = dat_last{i_ex}.(fldnames{i_fld}).stats.(DRUGCOND).(STATTYPE){CHANNEL}(1);
            tmp = abs(tmp);
            last_raw(i_fld) = tmp;
        end
        
        % average across TF conditions
        last_avg(i_pamp) = mean(last_raw);
        first_avg(i_pamp) = mean(first_raw);
        
    end
    

    % compile the data across experiments
    pop_first{i_ex} = first_avg;
    pop_last{i_ex} = last_avg;
    pop_pamp{i_ex} =  ppval(cal.pp_mW_mm2, pAmps_unique);
    
    
    
    switch lower(info_first{i_ex}.opsin)
        case 'chr2'
            p = plot(first_avg, last_avg, 'b-o', 'markerfacecolor', 'b');
        case 'chief_cit'
            p = plot(first_avg, last_avg, 'r-o', 'markerfacecolor', 'r');
        case 'chronos'
            p = plot(first_avg, last_avg, 'g-o', 'markerfacecolor', 'g');
    end
    
    % axis labels and such
    xlabel('first 5 amplitude');
    ylabel('last 5 amplitude');
end

% add a unity line
ymax = max(get(gca, 'YLim'));
xmax = max(get(gca, 'XLim'));
maxval = max([xmax, ymax]);
plot([0 maxval], [0 maxval], 'k');
axis square
axis tight

% plot the opsin current and FV against light power
figure, hold on
for i_ex = 1:numel(dat_first)
    
    switch lower(info_first{i_ex}.opsin)
        case 'chr2'
            pltclr = 'b';
        case 'chief_cit'
            pltclr = 'r';
        case 'chronos'
            pltclr = 'g';
    end
    
    % first5 - last5
    plot(pop_pamp{i_ex}, pop_first{i_ex}-pop_last{i_ex}, '-o', 'color', pltclr);
    xlabel('Light Power')
    ylabel('first5 - last5')
end



%% LASER STIM SITE ANALYSIS


% figure out how many unique experiments there were
mouseNames = in(:,1);
sites = cellfun(@floor, in(:,2), 'uniformoutput', false); % use floor to eliminate the position information
comboname = cellfun(@(x,y) [x,num2str(y)], mouseNames, sites, 'uniformoutput', false);
[uniqueExpts, ~, exptID] = unique(comboname);

% initalize the pop structure
siteraw = {};

% iterate over the unique experiments and overlay the PPR metrics for FV,
% opsin current, and fEPSP slope
Nexpts = size(uniqueExpts, 1);
for i_ex = 1:Nexpts
    
    % order the files by proximal to distal (which should be hardcoded in
    % the excel workbook
    inds_mouse = find(exptID == i_ex);
    sites = cat(1, in{inds_mouse,2});
    opsinTag = structcat(info, 'opsin');
    opsinTag = unique(opsinTag(inds_mouse));
    opsinTag = opsinTag{1};
    
    [~, analysis_order] = sort(sites);
    inds_mouse = inds_mouse(analysis_order); % now the list is ordered by proximal to distal
    assert(numel(inds_mouse) == 2,'error: more than 2 stim sites')
    
    f = figure;
    set(f, 'defaulttextinterpreter', 'none', 'units', 'normalized', 'position', [0 0.1435 1 0.6620]);
    set(f, 'name', sprintf('%s, site: %d, opsin: %s', in{inds_mouse(1),1}, floor(in{inds_mouse(1),2}), opsinTag))
    cmap = copper(4);
    cmap_tf = {'tf10_', 'tf20_', 'tf40_', 'tf60_'};
    
    % assume that the primary channel is the one that's in L2/3
    CH_L23 = structcat(info(inds_mouse), 'stimSite');
    CH_L23 = unique(cell2mat(CH_L23));
    assert(numel(CH_L23)==1, 'ERROR: too many ''L23'' sites');
    if CH_L23 == 1
        CH_L5 = 2;
    else
        CH_L5 = 1;
    end
    
    % figure out the total number of tfconds, across stimulation locations
    tfconds = {};
    for i_site = 1:numel(inds_mouse);
        idx = inds_mouse(i_site);
        tfconds = cat(1, tfconds, fieldnames(dat{idx}));
    end
    assert(iscell(tfconds), 'ERROR: tfconds are unreliable')
    tfconds = unique(tfconds);
    
    
    for i_tf = 1:numel(tfconds)
        
        conds = {'nbqx_apv_cd2_ttx', 'FV_Na', 'synapticTransmission'}; % both opsin current verisons
        npltcols = numel(conds)+1;
        for i_cond = 1:numel(conds)
            
            % figure out what stat to use for this condition
            switch conds{i_cond}
                case {'nbqx_apv_cd2_ttx', 'FV_Na'}
                    STAT = 'diffval';
                case 'synapticTransmission'
                    STAT = 'slope';
            end
            
            for i_site = 1:numel(inds_mouse);
                
                idx = inds_mouse(i_site);
                
                
                % the tfconds could be laser or led, so skip a dat{idx} if
                % this site doesn't contain the apropriate stimulation
                % type. Also skip if the pharacology condition is not present
                if ~isfield(dat{idx}, tfconds{i_tf})
                    BLANKPLOT = true;
                    linestyle = '.';
                elseif ~isfield(dat{idx}.(tfconds{i_tf}).stats, conds{i_cond})
                    BLANKPLOT = true;
                    linestyle = '.';
                else
                    
                    BLANKPLOT = false;
                    
                    % figure out the plot colors. The LED should have a
                    % different color as the laser. and I should make sure that
                    % all LED stim is in L2/3, all laser stim is in L5;
                    if i_site == 1; % L2/3
                        ledcheck = regexpi(tfconds{i_tf}, '_led');
                        assert(~isempty(ledcheck), 'ERROR: led should be L2/3')
                        linestyle = ':';
                    elseif i_site == 2; % L5
                        lasercheck = regexpi(tfconds{i_tf}, '_laser');
                        assert(~isempty(lasercheck), 'ERROR: laser should be L5')
                        linestyle = '-';
                    end
                end
                
                
                
                % figure for L2/3 recording channel goes on top, figure for L5 recording channel goes on bottom
                for i_pltrow = 1:2;
                    
                    if i_pltrow == 1
                        CHANNEL = CH_L23;
                    elseif i_pltrow == 2;
                        CHANNEL = CH_L5;
                    end
                    
                    if BLANKPLOT
                        tmp_raw = nan;
                    else
                        if (CHANNEL == 2) && numel(dat{idx}.(tfconds{i_tf}).stats.(conds{i_cond}).(STAT)) == 1
                            error('Could not find the data')
                        end
                        
                        
                        tmp_raw = dat{idx}.(tfconds{i_tf}).stats.(conds{i_cond}).(STAT){CHANNEL};
                        if ~any(regexpi(conds{i_cond}, 'ttx')) % PPR for everything but opsin current
                            tmp_raw = tmp_raw ./ tmp_raw(1);
                        end
                    end
                    
                    
                    % store the data into a population structure for future
                    % analysis:
                    % sitepop{i_ex}{i_stimSite}.(recLayer).(tfcond).(pharmCond)
                    % 
                    % i_stimSite will be a 1 or 2. 1 is in L2/3, 2=L5
                    % recLayer is L2/3 or L5.
                    % tfcond will be LED_yyy or Laser_yyy
                    if ~isnan(tmp_raw)
                        recLayer = {'recL23', 'recL5'};
                        stimSite = {'stimL23' 'stimL5'};
                        recLayer = recLayer{i_pltrow};
                        stimSite = stimSite{i_site};
                        siteraw{i_ex}.(stimSite).(recLayer).(tfconds{i_tf}).(conds{i_cond}) = tmp_raw;
                        siteraw{i_ex}.opsin = opsinTag;
                    end
                    
                    
                    
                    subplot(2, npltcols, (i_pltrow-1)*npltcols + i_cond)
                    hold on,
                    tfcolor_idx = find(strcmpi(tfconds{i_tf}(1:5), cmap_tf));
                    plot(1:numel(tmp_raw), tmp_raw, linestyle, 'color', cmap(tfcolor_idx,:), 'linewidth', 2)
                    axis tight
                    
                    if i_pltrow == 1
                        if any(regexpi(conds{i_cond}, '_ttx'))
                            title(info{inds_mouse(1)}.opsin)
                        else
                            title(conds{i_cond})
                        end
                        ylabel('Record in L2/3')
                    else
                        ylabel('Record in L5')
                    end
                    
                    % adjust the axes
                    switch conds{i_cond}
                        case {'FV_Na', 'synapticTransmission'}
                            yvals = get(gca, 'ylim');
                            set(gca, 'ylim', 1.05.* [0, max([yvals(2), 1])]);
                        otherwise
                    end
                end
                
            end
        end
        
    end
    
    
    % add a picture of the slice, only once per figure
    s = subplot(1, npltcols, npltcols);
    
    mousename = uniqueExpts{i_ex}(1:end-1);
    siteidx =  uniqueExpts{i_ex}(end);
    cd([GL_DATPATH, mousename, filesep, 'Other'])
    d = dir;
    d.name;
    
    photoprefix = [mousename, '_site', num2str(siteidx)];
    for i_photo = 1:numel(d)
        if strncmpi(d(i_photo).name, photoprefix, numel(photoprefix))
            img = double(imread(d(i_photo).name));
            iminfo = imfinfo(d(i_photo).name);
            maxpixval = prctile(img(:), [99.9]);
            normfact = min([maxpixval, 2^iminfo.BitDepth-1]);
            img = (img ./ normfact) .* (2^iminfo.BitDepth-1); % auto adjust LUT
            
            imshow(uint8(img));
            s.Position = [0.72, 0.1, 0.27, 0.85];
            title(sprintf('(0,0) is HS%d', CH_L23))
        end
    end
    
%     
%     % add an icon for the stimulus locations
%     centPos = round(ginput(1));
%     xy = cell2mat(structcat(info(inds_mouse), 'optStimCords'));
%     pixperum = pixPerMicron(size(img,1), size(img,2));
%     xy = round(xy .* pixperum); %now in pix
%     xy = bsxfun(@plus, xy, centPos); % pix relative to neuron
%     hold on,
%     for a = 1:size(xy,1)
%         plot(xy(a,1), xy(a,2), 'o', 'markeredgecolor', cmap(a,:), 'markerfacecolor', cmap(a,:), 'markersize', 10)
%     end
    
    % update plots in quasi realtime
    drawnow
    
end

%% POPULATION SUMMARY FOR L5 STIM ANALYSIS

% plot fEPSP slopes measured in L2/3, but stimulated in L2/3 or L5. Do this
% for each opsin separately. 
USEALLLFP = false;

% initialize the population structure
opsinTypes = {'chronos', 'chief_all'};
tfconds.stimL23 = {'tf10_led', 'tf20_led', 'tf40_led', 'tf60_led'};
tfconds.stimL5 =  {'tf10_laser', 'tf20_laser', 'tf40_laser', 'tf60_laser'};
stimpos = {'stimL23', 'stimL5'};
recpos = 'recL23';
conds = {'synapticTransmission', 'FV_Na'};
sitepop = [];
for i_opsin = 1:numel(opsinTypes);
    for i_site = 1:numel(stimpos)
        for i_tf = 1:numel(tfconds.(stimpos{i_site}));
            for i_cond = 1:numel(conds)
                sitepop.(opsinTypes{i_opsin}).(recpos).(stimpos{i_site}).(tfconds.(stimpos{i_site}){i_tf}).(conds{i_cond}) = [];
            end
        end
    end
end

% compile the data from the siteraw structure into the sitepop structure
for i_ex = 1:numel(siteraw)
    
    tmpopsin = siteraw{i_ex}.opsin;
    
    for i_site = 1:numel(stimpos)
        
        tffields = fieldnames(siteraw{i_ex}.(stimpos{i_site}).(recpos));
        Ntfs = numel(tffields);
        for i_tf = 1:Ntfs
            
            for i_cond = 1:numel(conds)
                
                % make sure that the drug condition exists
                if ~isfield(siteraw{i_ex}.(stimpos{i_site}).(recpos).(tffields{i_tf}), conds{i_cond})
                    continue
                end
                
                % pull out the data from the siteraw structure
                tmpvals = siteraw{i_ex}.(stimpos{i_site}).(recpos).(tffields{i_tf}).(conds{i_cond});
                
                % grab the pop structure
                tmppop = sitepop.(tmpopsin).(recpos).(stimpos{i_site}).(tffields{i_tf}).(conds{i_cond});
                
                % concatenate the data, and put it back into the pop
                % structure
                tmppop = cat(1, tmppop, tmpvals(1:7));
                sitepop.(tmpopsin).(recpos).(stimpos{i_site}).(tffields{i_tf}).(conds{i_cond}) = tmppop;
            end
        end
    end
    
end


if USEALLLFP
    % load the big LFP dataset
    poplfp = load('lfp_pop_combinechief.mat');
end

% ploting routines
pltclr = copper(6);
for i_opsin = 1:numel(opsinTypes)
    
    f = figure;
    f.Units = 'Normalized';
    f.Position = [0.1007    0.1979    0.7109    0.5544];
    f.Name = opsinTypes{i_opsin};
    
    for i_cond = 1:numel(conds)
        
        subplot(1,numel(conds), i_cond), hold on,
        title(conds{i_cond})
        
        lineseries = {'--', '-'};
        %build a sensible legend
        plot(nan, nan, '--k', 'linewidth', 2)
        plot(nan, nan, '-k', 'linewidth', 3)
        legend(stimpos, 'location', 'southeast')
        legend boxoff
        
        for i_site = 1:numel(stimpos)
            
            if USEALLLFP && strcmpi(stimpos{i_site}, 'stimL23')
                
                tmpdat = poplfp.pop.(opsinTypes{i_opsin}).pnp1.(conds{i_cond}).diffval;
                trim_dat = cellfun(@(x) x(:,1:7), tmpdat{3}, 'uniformoutput', false);
                uniqueTFs = unique(tmpdat{1});
                uniqueTFs(uniqueTFs>60)=[];
                Ntfs = numel(uniqueTFs);
                for i_tf = 1:Ntfs
                    
                    l_tf = tmpdat{1} == uniqueTFs(i_tf);
                    catdat = cat(1, trim_dat{l_tf});
                    xbar = mean(catdat,1);
                    sem = stderr(catdat, 1);
                    
                    my_errorbar(1:numel(xbar), xbar, sem, lineseries{i_site}, 'color', pltclr(i_tf,:), 'linewidth', 3);
                    
                end
                
            else
                
                tffields = fieldnames(sitepop.(opsinTypes{i_opsin}).(recpos).(stimpos{i_site}));
                Ntfs = numel(tffields);
                for i_tf = 1:Ntfs
                    tmpvals = sitepop.(opsinTypes{i_opsin}).(recpos).(stimpos{i_site}).(tffields{i_tf}).(conds{i_cond});
                    xbar = mean(tmpvals, 1);
                    sem = stderr(tmpvals,1);
                    
                    my_errorbar(1:numel(xbar), xbar, sem, lineseries{i_site}, 'color', pltclr(i_tf,:), 'linewidth', 3);
                end
                
                
                
            end
        end
        
        
        yvals = get(gca, 'ylim');
        ylim([0, max([yvals(2), 1])])
        
    end
end


%% COMPARE PLASTICITY ACROSS BRAIN AREAS




ERRBARS = true;
NPULSES = 7;
PHARMCONDITION = 'synapticTransmission'; % could be 'synapticTransmission' or 'nbqx_apv_cd2_ttx' or 'FV_Na'
STATTYPE = 'pnp1'; % could be 'raw' or 'pnp1'

% load in data from each area
load('pop_al_stimL23_recL23.mat');
al_pop = pop;

load('pop_pm_stimL23_recL23.mat');
pm_pop = pop;


%
% make a 3 x N plot. one row for each opsin, one column for each TF. Plot
% PPRs as a function of pulse number
%

F = figure;
set(F, 'position', [109    31   891   754])

opsinTypes = {'chr2', 'chief_all', 'chronos'};
for i_opsin = 1:numel(opsinTypes)
    
    % grab the data
    pm_tmp_dat = pm_pop.(opsinTypes{i_opsin}).(STATTYPE).(PHARMCONDITION).diffval;
    al_tmp_dat = al_pop.(opsinTypes{i_opsin}).(STATTYPE).(PHARMCONDITION).diffval;
    
    pm_tfs = pm_tmp_dat{1};
    pm_tfs = unique(pm_tfs);
    al_tfs = al_tmp_dat{1};
    al_tfs = unique(al_tfs);
    unique_tfs = unique([pm_tfs ; al_tfs]);
    Ntfs = numel(unique_tfs);
    
    for i_tf = 1:Ntfs;
        
        % now do the plotting
        subplot(3, Ntfs, (i_opsin-1).*Ntfs + i_tf)
        hold on,
        mytitle(sprintf('%s %dHz', opsinTypes{i_opsin}, unique_tfs(i_tf)))
        
        
        % plot for PM
        tf_idx = pm_tmp_dat{1} == unique_tfs(i_tf);
        l_pulses = cellfun(@(x) size(x,2)>=NPULSES, pm_tmp_dat{3});
        tf_idx = tf_idx & l_pulses; % only analyzes experiments where there are >= NPULSES per train
        trim_dat = cellfun(@(x) x(:,1:NPULSES), pm_tmp_dat{3}(tf_idx), 'uniformoutput', false);
        cat_dat = vertcat(trim_dat{:});
        pm_N = size(cat_dat, 1);
        pm_leg = sprintf('PM N=%d', pm_N);
        
        xbar = nanmean(cat_dat, 1);
        sem = nanstd(cat_dat, [], 1) ./ sqrt(sum(~isnan(cat_dat), 1));
        
        [~, pltclr] = hvaPlotColor('pm');
        if ERRBARS
            h_pm = my_errorbar(1:NPULSES, xbar, sem, '-', 'color', pltclr, 'linewidth', 2);
        else
            h_pm = plot(1:NPULSES, cat_dat', '-', 'color', pltclr, 'linewidth', 2);
        end
        
        
        % plot for AL
        tf_idx = al_tmp_dat{1} == unique_tfs(i_tf);
        l_pulses = cellfun(@(x) size(x,2)>=NPULSES, al_tmp_dat{3});
        tf_idx = tf_idx & l_pulses; % only analyzes experiments where there are >= NPULSES per train
        trim_dat = cellfun(@(x) x(:,1:NPULSES), al_tmp_dat{3}(tf_idx), 'uniformoutput', false);
        cat_dat = vertcat(trim_dat{:});
        al_N = size(cat_dat, 1);
        al_leg = sprintf('AL N=%d', al_N);
        
        xbar = nanmean(cat_dat, 1);
        sem = nanstd(cat_dat, [], 1) ./ sqrt(sum(~isnan(cat_dat), 1));
        
        [~, pltclr] = hvaPlotColor('al');
        if ERRBARS
            h_al = my_errorbar(1:NPULSES, xbar, sem, '-', 'color', pltclr, 'linewidth', 2);
        else
            h_al = plot(1:NPULSES, cat_dat', '-', 'color', pltclr, 'linewidth', 2);
        end
        
        % tidy up
        axis tight
        xlabel('Pulse number')
        if i_tf==1
            ylabel(sprintf('fEPSP slope \n Pn:P1 ratio'))
        end
        yvals = get(gca, 'ylim');
        yvals(1) = min([0, yvals(1)]);
        set(gca, 'ylim', yvals);
        legend([h_pm, h_al], {pm_leg, al_leg}, 'location', 'best')
        legend boxoff
    end
end



%% COMPARE OPSIN CURRENTS (INTRA VS LFP)

fin

STIMSITE = true;
NORMVALS = true;
LOWPOWER = true;
COMBINE_CHIEF = true;

% need to load the data
load('lfp_all_pow_combine_chief.mat');
pop.lfp.dat = dat; clear dat;
pop.lfp.info = info; clear info;

load('intra_lowRa_all_pow_combine_chief.mat'); % or 'intra_lowRa_all_pow_split_chief.mat'
pop.intra.dat = dat; clear dat;
pop.intra.info = info; clear info;


tforder = [10, 20, 40, 60, 100];

% preallocate all the arrays
if COMBINE_CHIEF
    opsins = {'chr2', 'chronos', 'chief_all'};
else
    opsins = {'chr2', 'chief_cit', 'chief_flx', 'chronos'};
end
for i_fld = {'lfp', 'intra'};
    Nexp = numel(pop.(i_fld{1}).dat);
    for i_opsin = opsins;
        catdat.(i_fld{1}).(i_opsin{1}).raw = nan(numel(tforder), 7, Nexp); 
    end
end

for i_fld = {'lfp', 'intra'}
    
    
    Nexp = numel(pop.(i_fld{1}).dat);
    for i_ex = 1:Nexp
        
       ex_opsin = lower(pop.(i_fld{1}).info{i_ex}.opsin);

       % which channel should be analyzed?
       CHANNEL = pop.(i_fld{1}).info{i_ex}.stimSite;
       if isnan(CHANNEL) % cases where neither recording sites were targeted
           continue
       end
       Nchannels = sum(pop.(i_fld{1}).info{i_ex}.ignoreChans);
       if Nchannels == 1
           if ~STIMSITE
               continue % no other stim site to show...
           elseif CHANNEL == 2 % cases where CH = 2, but only one channel recorded
               CHANNEL = 1;
           end
       else
           if ~STIMSITE
               if CHANNEL == 1
                   CHANNEL = 2;
               elseif CHANNEL == 2
                   CHANNEL = 1;
               end
           end
       end
       
       
       % extract the raw data for the TTX condition
       tfnames = fieldnames(pop.(i_fld{1}).dat{i_ex});
       for i_tf = 1:numel(tfnames)
           
           drugconds = fieldnames(pop.(i_fld{1}).info{i_ex}.(tfnames{i_tf}));
           ttxidx = cellfun(@(x) regexpi(x, 'ttx'), drugconds, 'uniformoutput', false);
           ttxidx = cellfun(@(x) ~isempty(x), ttxidx);
           if ~any(ttxidx)
               continue
           else
               assert(sum(ttxidx) == 1, 'ERROR: found multiple TTX indx')
           end
           ttxcond = drugconds{ttxidx};
           
           % allow the user to only anlyze data that came from weak
           % expression intracellular data
           if LOWPOWER
               if strcmpi(i_fld{1}, 'intra')
                   tmp_pamp = pop.(i_fld{1}).info{i_ex}.(tfnames{i_tf}).(ttxcond).pAmp;
                   if tmp_pamp <7; continue; end
               end
           end
           
           tmptf = pop.(i_fld{1}).info{i_ex}.(tfnames{i_tf}).(ttxcond).pTF;
           rowidx = tforder == tmptf;
           
           tmpraw = pop.(i_fld{1}).dat{i_ex}.(tfnames{i_tf}).stats.(ttxcond).diffval{CHANNEL};
           if numel(tmpraw)<7; continue; end
           catdat.(i_fld{1}).(ex_opsin).raw(rowidx, 1:7, i_ex) = tmpraw(1:7);
           
               
       end
    end
end


% normalize the data, and keep track of the xbar and sem for each TF
for i_fld = {'lfp', 'intra'}
    for i_opsin = opsins
        
        % grab the raw data and normalize by the first pulse
        tfval = 40;
        tmpdat = catdat.(i_fld{1}).(i_opsin{1}).raw;
        p1_vals = tmpdat(:,1,:); % all 1st pulses of all TFs and expts;
        tmpdat = bsxfun(@rdivide, tmpdat, p1_vals);
        
        if NORMVALS
            
            % normalize the difference between the P1 and P7 for a particular
            % TF (which depends on the opsin)
            diff_from_1 = 1-tmpdat;
            normfact = mean(diff_from_1(tforder == tfval,5:7,:), 2);
            norm_diffs = bsxfun(@rdivide, diff_from_1, normfact);
            
            
            catdat.(i_fld{1}).(i_opsin{1}).normvals = 1-norm_diffs;
            
        else
            catdat.(i_fld{1}).(i_opsin{1}).normvals = tmpdat;
        end
        
        
    end
end

% plot the results

for i_opsin = opsins
    
    f = figure;
    f.Position = [415   417   348   344];
    f.Name = i_opsin{1};
    hold on,
    
    tmp_lfp = catdat.lfp.(i_opsin{1}).normvals;
    xbar_lfp = nanmean(tmp_lfp, 3);
        
    tmp_intra = catdat.intra.(i_opsin{1}).normvals;
    xbar_intra = nanmean(tmp_intra, 3);
    i_opsin
    N = sum(~isnan(tmp_intra),3)
    
    cmap = copper(numel(tforder));
    for i_tf = 1:numel(tforder)
        plot(xbar_lfp(i_tf,:), '--', 'color', cmap(i_tf,:), 'linewidth', 2)
        plot(xbar_intra(i_tf,:), '-', 'color', cmap(i_tf,:), 'linewidth', 2)
    end
    xlabel('pulse number')
    ylabel('opsin current');
    if ~NORMVALS
        ylim([0 1])
    end
    
end


% 
% 
% % compile the stats
% for i_opsin = 1:numel(opsins)
%     
%     
%     % grab the data
%     tmp_lfp = catdat.lfp.(opsins{i_opsin}).normvals;
%     tmp_intra = catdat.intra.(opsins{i_opsin}).normvals;
%     
%     % pull out pulse 7 data
%     tmp_lfp = tmp_lfp(:,7,:);
%     tmp_lfp = permute(tmp_lfp, [1,3,2]);
%     allnans = sum(isnan(tmp_lfp),1) == size(tmp_lfp,1);
%     tmp_lfp = tmp_lfp(:, ~allnans);
%     
%     tmp_intra = tmp_intra(:,7,:);
%     tmp_intra = permute(tmp_intra, [1,3,2]);
%     allnans = sum(isnan(tmp_intra),1) == size(tmp_intra,1);
%     tmp_intra = tmp_intra(:, ~allnans);
%     
%     % the flx version of chief was not tested at 100Hz for the LFP, cut
%     % these data from the intracellular traces
%     if ~COMBINE_CHIEF && strcmpi(opsins{i_opsin}, 'chief_flx');
%         tmp_intra = tmp_intra(1:4,:);
%         tmp_lfp = tmp_lfp(1:4,:);
%     end
%     
%     % compile the lfp data for the anova
%     group_tf = [];
%     group_method = {};
%     YY = [];
%     
%     TFs = [10 20 40 60 100];
%     Ntf = size(tmp_lfp, 1);
%     Nex = size(tmp_lfp, 2);
%     for i_row = 1:Ntf;
%         for i_col = 1:Nex
%             if ~isnan(tmp_lfp(i_row, i_col))
%                 group_tf = cat(1, group_tf, TFs(i_row));
%                 group_method = cat(1, group_method, 'LFP');
%                 YY = cat(1, YY, tmp_lfp(i_row, i_col));
%             end
%         end
%     end
%     
%     Ntf = size(tmp_intra, 1);
%     Nex = size(tmp_intra, 2);
%     for i_row = 1:Ntf;
%         for i_col = 1:Nex
%             if ~isnan(tmp_intra(i_row, i_col))
%                 group_tf = cat(1, group_tf, TFs(i_row));
%                 group_method = cat(1, group_method, 'Intra');
%                 YY = cat(1, YY, tmp_intra(i_row, i_col));
%             end
%         end
%     end
%     
%     % run the anova
%     [p, tab, stats] = anovan(YY, {group_tf, group_method}, 'model', 2, 'varnames', {'Temp Freq', 'Method'});
%     
% end
% 


%% COMPARE FIBER VOLLEYS (OPTICAL VS. ELECTRICAL)

fin

STIMSITE = true;
COMBINE_CHIEF = true;

% need to load the data
load('lfp_all_pow_combine_chief.mat');
tmp.optical.dat = dat; clear dat;
tmp.optical.info = info; clear info;

load('estim_all_pow.mat'); % intra_all_pow.mat or intra_lowRa_all_pow.mat
tmp.estim.dat = dat; clear dat;
tmp.estim.info = info; clear info;

% concatenate the datasets
pop = cat(2, tmp.estim.dat, tmp.optical.dat);
info = cat(2, tmp.estim.info, tmp.optical.info);


tforder = [10, 20, 40, 60, 100];

% preallocate all the arrays
if COMBINE_CHIEF
    opsins = {'estim', 'chr2', 'chronos', 'chief_all'};
else
    opsins = {'estim', 'chr2', 'chronos', 'chief_cit', 'chief_flx'};
end

Nexp = numel(pop);
for i_opsin = opsins;
    catdat.(i_opsin{1}).raw = nan(numel(tforder), 7, Nexp);
end


for i_ex = 1:Nexp
    
    ex_opsin = lower(info{i_ex}.opsin);
    
    % which channel should be analyzed?
    CHANNEL = info{i_ex}.stimSite;
    if isnan(CHANNEL) % cases where neither recording sites were targeted
        continue
    end
    Nchannels = sum(info{i_ex}.ignoreChans);
    if Nchannels == 1
        if ~STIMSITE
            continue % no other stim site to show...
        elseif CHANNEL == 2 % cases where CH = 2, but only one channel recorded
            CHANNEL = 1;
        end
    else
        if ~STIMSITE
            if CHANNEL == 1
                CHANNEL = 2;
            elseif CHANNEL == 2
                CHANNEL = 1;
            end
        end
    end
    
    
    % extract the raw data for the TTX condition
    tfnames = fieldnames(pop{i_ex});
    for i_tf = 1:numel(tfnames)
        
        drugconds = fieldnames(info{i_ex}.(tfnames{i_tf}));
        idx_FV_Na = strcmpi('FV_Na', drugconds);
        if ~any(idx_FV_Na)
            continue
        else
            assert(sum(idx_FV_Na) == 1, 'ERROR: found multiple TTX indx')
        end
        
        tmptf = info{i_ex}.(tfnames{i_tf}).FV_Na.pTF;
        if any(tmptf > max(tforder))
            continue
        end
        rowidx = tforder == tmptf;
        
        tmpraw = pop{i_ex}.(tfnames{i_tf}).stats.FV_Na.diffval{CHANNEL};
        if numel(tmpraw)<7; continue; end
        catdat.(ex_opsin).raw(rowidx, 1:7, i_ex) = tmpraw(1:7);
    end
end



% normalize the data, and keep track of the xbar and sem for each TF
for i_opsin = opsins
    
    % grab the raw data and normalize by the first pulse
    tmpdat = catdat.(i_opsin{1}).raw;
    p1_vals = tmpdat(:,1,:); % all 1st pulses of all TFs and expts;
    tmpdat = bsxfun(@rdivide, tmpdat, p1_vals);
    catdat.(i_opsin{1}).normvals = tmpdat;
end


% compile the stats
for i_opsin = 2:numel(opsins)
    
    fprintf('  ********  %s   ****\n', opsins{i_opsin})
    
    % grab the data
    tmp_lfp = catdat.(opsins{i_opsin}).normvals;
    tmp_estim = catdat.estim.normvals;
    
    % pull out pulse 7 data
    tmp_lfp = tmp_lfp(:,7,:);
    tmp_lfp = permute(tmp_lfp, [1,3,2]);
    allnans = sum(isnan(tmp_lfp),1) == size(tmp_lfp,1);
    tmp_lfp = tmp_lfp(:, ~allnans);
    
    tmp_estim = tmp_estim(:,7,:);
    tmp_estim = permute(tmp_estim, [1,3,2]);
    allnans = sum(isnan(tmp_estim),1) == size(tmp_estim,1);
    tmp_estim = tmp_estim(:, ~allnans);
    
    % the flx version of chief was not tested at 100Hz for the LFP, cut
    % these data from the intracellular traces
    switch opsins{i_opsin}
        case {'chief_flx', 'chronos'}
            tmp_estim = tmp_estim(1:4,:);
            tmp_lfp = tmp_lfp(1:4,:);
        case 'chr2';
            tmp_estim = tmp_estim(1:3,:);
            tmp_lfp = tmp_lfp(1:3,:);       
    end
    
    % loop over each of the TFs and do a ttest
    pval_p7 = [];
    for i_tf = 1:size(tmp_lfp, 1)
       ttest_estim = tmp_estim(i_tf,:);
       ttest_lfp = tmp_lfp(i_tf, :);
       
       [~, pval_p7, ~, stats] = ttest2(ttest_estim(:), ttest_lfp(:));
       
    end
    
    % look at the differences in means
    estim_avg = nanmean(tmp_estim, 2)
    lfp_avg = nanmean(tmp_lfp, 2)
    
    % compile the lfp data for the anova
    group_tf = [];
    group_method = {};
    YY = [];
    
    TFs = [10 20 40 60 100];
    Ntf = size(tmp_lfp, 1);
    Nex = size(tmp_lfp, 2)
    for i_row = 1:Ntf;
        for i_col = 1:Nex
            if ~isnan(tmp_lfp(i_row, i_col))
                group_tf = cat(1, group_tf, TFs(i_row));
                group_method = cat(1, group_method, 'LFP');
                YY = cat(1, YY, tmp_lfp(i_row, i_col));
            end
        end
    end
    
    Ntf = size(tmp_estim, 1);
    Nex = size(tmp_estim, 2)
    for i_row = 1:Ntf;
        for i_col = 1:Nex
            if ~isnan(tmp_estim(i_row, i_col))
                group_tf = cat(1, group_tf, TFs(i_row));
                group_method = cat(1, group_method, 'Intra');
                YY = cat(1, YY, tmp_estim(i_row, i_col));
            end
        end
    end
    
    % run the anova
    [p, tab, stats] = anovan(YY, {group_tf, group_method}, 'model', 2, 'varnames', {'Temp Freq', 'Method'});
end




%% OPSIN CURRENT FALL OFF WITH DISTANCE

fin

% load in the workbook 
fname = [GL_DOCUPATH, 'Other_workbooks', filesep, 'fiberVolleyCellList.xlsx'];
[~,~,wb_expt] = xlsread(fname, 3);

header = wb_expt(1,:);
for i_head = 1:numel(header)
    tmp = deblank(header{i_head});
    tmp(isspace(tmp)) = [];
    idx.(tmp) = i_head;
end

wb_expt(1,:) = [];
mouseNames = wb_expt(:, idx.MouseName);
sites = cell2mat(wb_expt(:, idx.Site));
fnames = wb_expt(:, idx.FileName);
HS1_loc = cell2mat(wb_expt(:, [idx.HS1X, idx.HS1Y]));
HS2_loc = cell2mat(wb_expt(:, [idx.HS2X, idx.HS2Y]));
led_loc = cell2mat(wb_expt(:, [idx.StimX, idx.StimY]));

name_site = cellfun(@(x,y) [x, '_', num2str(y)], mouseNames, num2cell(sites), 'uniformoutput', false);
unique_expts = unique(name_site);

dat = {};
for i_ex = 1:numel(unique_expts);
    
    l_expt = strcmpi(name_site, unique_expts{i_ex});
    
    expt_fnames = fnames(l_expt)
    expt_HS1 = unique(HS1_loc(l_expt,:), 'rows');
    expt_HS2 = unique(HS2_loc(l_expt,:), 'rows');
    expt_led = led_loc(l_expt,:);
    
    figure
    [hs1stat, hs2stat, hs1dist, hs2dist] = deal(nan(numel(expt_fnames),1));
    for i_ax = 1:numel(expt_fnames)
        
        ax = abfobj(expt_fnames{i_ax});
        
        % find the stim on time
        above_thresh = ax.dat(:, ax.idx.LED_470)' > 0.5;
        crossing_up = [false, diff(above_thresh)==1];
        assert(sum(crossing_up)==1, 'ERROR: too many stim on times')
        
        % pull out the data, baseline subtract, average
        preSamps = ceil(0.050 .* ax.head.sampRate);
        pulse_width_samps = sum(above_thresh);
        photoDelay = 2;
        postSamps = ceil(0.100 .* ax.head.sampRate);
        totalPostSamps = postSamps + photoDelay + pulse_width_samps;
        idx = find(crossing_up)-preSamps : find(crossing_up)+totalPostSamps;
        snips = ax.dat(idx, :, :);
        
        bkgnd = mean(snips(1:preSamps, :, :), 1);
        snips = bsxfun(@minus, snips, bkgnd);
        
        snips = mean(snips, 3);
        
        tt = ((0:size(snips,1)-1)-preSamps) ./ ax.head.sampRate .* 1000;
        anly_window = (tt>=1) & (tt<2);
        
        
        if ~any(isnan(expt_HS1))
            expt_hs1 = snips(:, ax.idx.HS1_Vm);
            expt_hs1 = butterfilt(expt_hs1, 1000, ax.head.sampRate, 'low', 1);
            
            subplot(2,1,1), hold on,
            plot(tt, expt_hs1)
            
            % get the hs1 stat
            hs1stat(i_ax) = mean(expt_hs1(anly_window));
            signval = norm(expt_led(i_ax,:)) - norm(expt_HS1);
            tmp = norm(expt_led(i_ax,:)-expt_HS1);
            if signval<0
                hs1dist(i_ax) = -tmp;
            else
                hs1dist(i_ax) =tmp;
            end
        end
        
        
        if ~any(isnan(expt_HS2))
            expt_hs2 = snips(:, ax.idx.HS2_Vm);
            expt_hs2 = butterfilt(expt_hs2, 1000, ax.head.sampRate, 'low', 1);
            
            subplot(2,1,2), hold on,
            plot(tt, expt_hs2)
            
            % get the HS2 stat
            hs2stat(i_ax) = mean(expt_hs2(anly_window));
            signval = norm(expt_led(i_ax,:)) - norm(expt_HS2);
            tmp = norm(expt_led(i_ax,:)-expt_HS2);
            if signval<0
                hs2dist(i_ax) = -tmp;
            else
                hs2dist(i_ax) =tmp;
            end
        end
        
    end
    
    % store the stats across files within an experiment
    dat{i_ex}.hs1stat = hs1stat;
    dat{i_ex}.hs1dist = hs1dist;
    dat{i_ex}.hs2stat = hs2stat;
    dat{i_ex}.hs2dist = hs2dist;
    drawnow
    
    
end

figure
for i_ex = 1:numel(unique_expts)
    
    
    subplot(2,1,1), hold on,
    
    %hs1
    minx = min(dat{i_ex}.hs1dist);
    maxx = max(dat{i_ex}.hs1dist);
    x = minx:maxx;
    plot(dat{i_ex}.hs1dist, dat{i_ex}.hs1stat, 'ko:')
    if ~any(isnan(minx))
        plot(x, ppval(spline(dat{i_ex}.hs1dist, dat{i_ex}.hs1stat), x), 'k')
    end
    
    %hs2
    minx = min(dat{i_ex}.hs2dist);
    maxx = max(dat{i_ex}.hs2dist);
    x = minx:maxx;
    plot(dat{i_ex}.hs2dist, dat{i_ex}.hs2stat, 'bo:')
    plot(x, ppval(spline(dat{i_ex}.hs2dist, dat{i_ex}.hs2stat), x), 'b')
    
    subplot(2,1,2), hold on,
    plot(dat{i_ex}.hs1dist, i_ex, 'ko')
    plot(dat{i_ex}.hs2dist, i_ex+.5, 'bo')
end


% dump all the stats in the same bucket, bin distances, and make a summary
% figure
USEBOTHCHANS = false;
popstat = [];
popdist = [];
dbin = 100;
edges = -150:dbin:650;
for i_ex = 1:numel(unique_expts)
    if USEBOTHCHANS
        popstat = cat(1, popstat, dat{i_ex}.hs1stat, dat{i_ex}.hs2stat);
        popdist = cat(1, popdist, dat{i_ex}.hs1dist, dat{i_ex}.hs2dist);
    else
        
        % only look at one recording channel, and average duplicate points
        % in each bin
        tmpdist = dat{i_ex}.hs2dist;
        tmpstat = dat{i_ex}.hs2stat;
        
        [~, ex_inds] = histc(tmpdist, edges);
        unique_inds = unique(ex_inds);
        corrected_dist = [];
        corrected_stat = [];
        for i_ind = 1:numel(unique_inds)
            l_ind = ex_inds == unique_inds(i_ind);
            corrected_dist(i_ind) = mean(tmpdist(l_ind));
            corrected_stat(i_ind) = mean(tmpstat(l_ind));
        end
        
        popstat = cat(1, popstat, corrected_stat(:));
        popdist = cat(1, popdist, corrected_dist(:));
        
    end
end

% remove the nans
l_nan = isnan(popstat);
popstat(l_nan) = [];
popdist(l_nan) = [];

% bin distances, and make a summary figure
[~, inds] = histc(popdist, edges);

bin_xbar = [];
bin_sem = [];
for i_bin = 1:numel(edges)
    bin_xbar(i_bin) = mean(popstat(inds==i_bin));
    bin_sem(i_bin) = stderr(popstat(inds==i_bin));
end

figure;
my_errorbar(edges+(dbin/2), bin_xbar, bin_sem);
xlabel('distance from rec site (um)')
ylabel('opsin current amplitude')
xlim([-150 600])

%% SERIES RESISTANCE FOR INTRACELLULAR RECORDINGS


% grab the fiber volley pop excel workbook
fname = [GL_DOCUPATH, 'Other_workbooks', filesep, 'fiberVolleyCellList.xlsx'];
[~,txt, raw] = xlsread(fname);
raw(size(txt,1)+1:end, :) = [];
raw(:,size(txt,2)+1:end) = [];
channelIdx = cellfun(@(x) ~isempty(x), regexpi(raw(1,:), 'CH\d'));
opsinIdx = strcmpi(raw(1,:), 'opsin');
idx_rmSweeps = strcmpi(raw(1,:), 'rmSweeps');
idx_mousename = strcmpi(raw(1,:), 'Mouse Name');
idx_fname =  strcmpi(raw(1,:), 'file name');
primaryStimSiteIdx = strcmpi(raw(1,:), 'Primary Stim Site (0,0)');
drugIdx = strcmpi(raw(1,:), 'drugs');
clear txt

% delete the header row
raw(1,:) = [];

% replace the 'file names' with fully qualified paths. I'm hoping that this
% allows for parallel operations
Nfiles = size(raw,1);
tmp_fnames = raw(:, idx_fname);
tmp_mousenames = raw(:, idx_mousename);
rmsweeps = raw(:, idx_rmSweeps);
prefix = cellfun(@(x,y) [x,y], repmat({GL_DATPATH}, Nfiles, 1), tmp_mousenames, 'uniformoutput', false);
prefix = cellfun(@(x,y) [x,y], prefix, repmat({[filesep, 'Physiology', filesep]}, Nfiles, 1), 'uniformoutput', false);
tmp_fnames = cellfun(@(x,y) [x,y], prefix, tmp_fnames, 'uniformoutput', false);
fid_path = cellfun(@(x,y) [x,y], tmp_fnames, repmat({'.abf'}, Nfiles, 1), 'uniformoutput', false);



% do the analysis
Nexpts = size(in,1);
popRa = nan(Nexpts, 1);
for i_ex = 1:Nexpts;
    
    % figure out what rows in the work book to pay attention to
    l_mouse = cellfun(@(x) ~isempty(x), regexp(raw(:,1), in{i_ex,1}));
    l_site = cell2mat(raw(:,2)) == in{i_ex,2}; 
    l_ttx = cellfun(@(x) ~isempty(regexpi(x, 'ttx')), raw(:, drugIdx));
    l_expt = l_mouse & l_site & l_ttx;
    assert(sum(l_expt)==1, 'ERROR: too many experiments')
    
    % run the analysis
    ax = abfobj(fid_path{l_expt});
    
    % get the series resistance for the proper channel
    ra = ax.getRa;
    channel = raw{l_expt, primaryStimSiteIdx};
    if (channel==2) && (size(ra.dat, 2)==1)
        channel=1;
    end
    ra = squeeze(ra.dat(:,channel,:));
    
    % remove bad trials if necessary    % remove some sweeps if need be
    if any(~isnan(rmsweeps{l_expt}))
        goodSweeps = true(size(ra,1),1);
        badSweeps = eval(['[',rmsweeps{l_expt},']']);
        goodSweeps(badSweeps) = false;
        ra = ra(goodSweeps,:);
    end
    
    popRa(i_ex) = mean(ra);
    
end

figure
histogram(popRa, 10)

average_Ra = mean(popRa)

l_gt20 = popRa > 20;
in(l_gt20, :)


%% JITTER ANALYSIS

% load a data file with cell attached spikes (from Tim).
[d, h, wf] = my_abfload('C:\Users\charlie\Desktop\14d17018.abf');

% select a single spike, reverse the sign to make it look like an LFP
template = -d(100600:101100,:,89);
tt = [0:numel(template)-1] ./ h.sampRate;
plot(tt, template)

% jitter the start time
jittered = nan(500, numel(template));
jit_max = 1e-3 ./ (1./h.sampRate);  % max number of timesteps to shift
for i_iter = 1:size(jittered,1)
    jittertime = unidrnd



%% GET TRIAL COUNTS AND SNR AND GENOTYPES

clc

ISLFP = false; % false only returns age, sex, genotype


mdb = initMouseDB('update', 'notext');


genotype = {};
sex = {};
date_inj = [];
date_record = [];
date_birth = [];
uniqueMice = unique(in(:,1));
Nmice = numel(uniqueMice);
for i_mouse = 1:Nmice
    
    % get genotype and sex
    [~, idx] = mdb.search(uniqueMice{i_mouse});
    genotype{i_mouse, 1} = mdb.mice{idx}.info.strain;
    sex{i_mouse, 1} = upper(mdb.mice{idx}.info.sex(1));
    
    
    % get injection date
    tmpDate = regexpi(uniqueMice{i_mouse}, '_\w*_', 'match');
    tmpDate = tmpDate{1}(2:end-1);
    if str2double(tmpDate(1:2)) <= 12
        % assume mmddyy
        date_inj(i_mouse,1) = datenum(tmpDate, 'mmddyy');
    else
        %assume yymmdd
        date_inj(i_mouse,1) = datenum(tmpDate, 'yymmdd');
    end
    
    % get the date of recording, grab the name of the first file from the
    % first cell
    tmpDate = mdb.mice{idx}.phys.cell(1).file(1).FileName;
    tmpDate = regexpi(tmpDate, '\d{4}_\d{2}_\d{2}', 'match'); % remove the file number
    tmpDate = regexprep(tmpDate{1}, '_', '');
    date_record(i_mouse,1) = datenum(tmpDate, 'yyyymmdd');
    
    % get the date of birth
    tmpDate = mdb.mice{idx}.info.dob;
    date_birth(i_mouse,1) = datenum(tmpDate, 'mm/dd/yyyy');
    
end

if ISLFP
    
    SNR = [];
    nexpts = numel(info);
    minExTrialCount = [];
    for a_ex = 1:nexpts
        fields = fieldnames(info{a_ex});
        
        fieldTcount = inf;
        for a_field = 1:numel(fields)
            
            % make sure you're looking at data and not header info
            if ~isfield(info{a_ex}.(fields{a_field}) , 'nbqx_apv_cd2_ttx')
                continue
            end
            
            conds = {'nbqx_apv_cd2', 'nbqx_apv_cd2_ttx'};
            
            for a_cond = 1:numel(conds)
                
                if ~isfield(info{a_ex}.(fields{a_field}), conds{a_cond});
                    warning('could not find data')
                end
                
                tlist = info{a_ex}.(fields{a_field}).(conds{a_cond}).realTrialNum;
                fieldTcount = min([fieldTcount, numel(tlist)]);
            end
            
            minExTrialCount = cat(1, minExTrialCount, fieldTcount);
            
            % now deal with SNR
            CHANNEL = info{a_ex}.stimSite;
            if ~info{a_ex}.ignoreChans(CHANNEL)
                continue
            else
                % occasionally CH1 is disabled, and CH2 occupies the
                % 1st and only column, adjust i_ch to account for these
                % cases
                Nchannels = sum(info{a_ex}.ignoreChans);
                if (Nchannels == 1) && (CHANNEL == 2)
                    CHANNEL = 1;
                end
            end
            
            noise = dat{a_ex}.(fields{a_field}).stats.FV_Na.bkgnd_sigma{CHANNEL};
            signal = dat{a_ex}.(fields{a_field}).stats.FV_Na.diffval{CHANNEL}(1);
            SNR = cat(1, SNR, abs(signal./noise));
            
            
            
        end
    end
    
    minimumAcrossExperiments = min(minExTrialCount)
    maxAcrossExperiments = max(minExTrialCount)
    meanTrialCount = mean(minExTrialCount)
    stdTrialCount = std(minExTrialCount)
    
    
    minSNR = min(SNR)
    meanSNR = mean(SNR)
    stdSNR = std(SNR)
end




% determine the age at recording and incubation time for each construct
popage = [];
for i_mouse = 1:Nmice
    
    
    % pull out the opsin for this mouse, just select the first experiment
    idx_expt = strcmpi(in(:,1), uniqueMice{i_mouse});
    idx_expt = find(idx_expt, 1, 'first');
    opsin = info{idx_expt}.opsin;
    
    if ~isfield(popage, opsin)
        popage.(opsin) = [];
        popage.(opsin).days_incubate = [];
        popage.(opsin).days_old = [];
    end
   
   incubation = date_record(i_mouse) - date_inj(i_mouse);
   popage.(opsin).days_incubate = cat(1, popage.(opsin).days_incubate, incubation);
   
   daysold = date_record(i_mouse) - date_birth(i_mouse);
   popage.(opsin).days_old = cat(1, popage.(opsin).days_old, daysold);
    
end

% print out the mean age at recording and days incubated
fields = fieldnames(popage);
for opsinType = fields'
    tmpmean = mean(popage.(opsinType{1}).days_old);
    fprintf('\n %s, mean age at recording = %.2f\n', opsinType{1},  tmpmean);
    
    tmpmean = mean(popage.(opsinType{1}).days_incubate);
    fprintf('%s, incubation in days = %.2f\n\n', opsinType{1},  tmpmean);
end

