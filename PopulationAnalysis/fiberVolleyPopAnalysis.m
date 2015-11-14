%% WHICH MICE SHOULD CONTRIBUTE?  [Main project]

% clear out the workspace
fin


% decide what experiment to run
EXPTTYPE = 1;
BRAINAREA = 'al';
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
        EXPTTYPE = 'all manuscript';
end

% grab the mouse names and sites from the excel workbook.
fname = [GL_DOCUPATH, 'Other_workbooks', filesep, 'fiberVolleyCellList.xlsx'];
[~,~,wb_expt] = xlsread(fname, 2);

wb_expt_nameidx = strcmpi(wb_expt(1,:), 'mouse name');
wb_expt_siteidx = strcmpi(wb_expt(1,:), 'site');
if strcmpi(EXPTTYPE, 'all manuscript')
    exptlistidx = strcmpi(wb_expt(1,:), 'main expt') | strcmpi(wb_expt(1,:), 'interleaved amps');
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

% % uncomment these lines for a few files that are useful for code
% % development:
%
% in = {'CH_150112_B', [1];...
%       'CH_141215_E', [2]};

% flag some strange data files where the channels were not properly
% indicated (channel 2 appears in the first and only column...)
EXCEPTIONS = {'EB_150529_A', 1; 'EB_150529_B', 1; 'EB_150630_D', 1};

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
    [dat{i_ex}, info{i_ex}] = fiberVolleyAnalysis(l_expt, raw, false, RMLINENOISE);
    
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
postPulseTime = 0.009; % in sec

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
                    end
                    
                    % strange cases
                    mouseMatch = strcmpi(in{i_ex,1}, EXCEPTIONS(:,1));
                    siteMatch = in{i_ex,2} == vertcat(EXCEPTIONS{:,2});
                    if any(mouseMatch & siteMatch)
                        %  HS2 is the data channel, but since HS1 wasn't
                        %  used, the data are in the first column
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

clc

% REMINDER: this is happening in separate cell-script because I will define
% the analysis region based off the average first pulse across TF conds.
% This requires grabbing all the snippets b/4 any analysis can proceed. The
% preceeding cell will need to be run before this one.

% should the analysis window be defined based off the average response
% following the first pulse (across conditions) or should the analysis
% window be unique for each pulse?
FIRSTPULSE = false;

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
    tags = {'FV_', 'ttx', 'synapticTransmission'};
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
        switch conds{i_cond}
            case 'synapticTransmission'
                photoDelay= 1e-3;
            otherwise
                photoDelay= 300e-6; % timeout following pulse offset (in sec)
        end
        anlyWindowInPts = ceil(450e-6 .* sampRate);
        
        for i_ch = 1:2;
            
            % deal with data and channels that need to be ignored.
            if ~info{i_ex}.ignoreChans(i_ch)
                continue
            else
                
                % strange cases
                mouseMatch = strcmpi(in{i_ex,1}, EXCEPTIONS(:,1));
                siteMatch = in{i_ex,2} == vertcat(EXCEPTIONS{:,2});
                if any(mouseMatch & siteMatch)
                    %  HS2 is the data channel, but since HS1 wasn't
                    %  used, the data are in the first column
                    if i_ch == 1; error('something went wrong'); end
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
                        [troughidx, peakidx]  = anlyMod_getWFepochs(snippet, tt, conds{i_cond}, pWidth, photoDelay, direction);
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
                        
                    elseif any(regexpi(conds{i_cond}, 'ttx'))
                        
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
                            % make sure none of the fit_dat points are
                            % positive because the fitting routine will
                            % assume that all the points are negative
                            critval_slope = log(abs(fit_dat(1))) - 2; % two orders of magnitude
                            if any(log(abs(fit_dat)) < critval_slope)
                                l_zero = log(abs(fit_dat)) <= critval_slope;
                                stopIdx = find(l_zero==1, 1, 'first')-1;
                                fit_tt = fit_tt(1:stopIdx);
                                fit_dat = fit_dat(1:stopIdx);
                                dataforfit = ~isempty(fit_dat);
                            end
                        end
                        
                        
                        if dataforfit
                            betas = [fit_tt(:), ones(size(fit_tt(:)))] \ log(abs(fit_dat(:)));
                        else
                            betas = [nan, nan];
                            startIdx = nan;
                        end
                        
                        % store the slope params
                        dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).tau_m{i_ch}(i_pulse) = betas(1);
                        dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).tau_b{i_ch}(i_pulse) = betas(2);
                        dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).tau_ind{i_ch}(i_pulse) = startIdx;
                        
                        
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

Nexpts = size(in,1);
for i_ex = 1:Nexpts
    
    
    % define the conditions that will get plotted
    switch EXPTTYPE
        case '4-AP'
            conds = {'ttx_cd2', 'ttx_cd2_4AP800', 'ttx_cd2_4AP1800'};
        case 'Baclofen'
            conds = {'ttx_cd2', 'ttx_cd2_bac10'};
        otherwise
            
            conds = {'FV_Na', 'nbqx_apv_cd2_ttx', 'synapticTransmission'}; % with cadmium
            %     conds = {'FV_Na', 'FV_Na_Ca2_mGluR', 'synapticTransmission'}; % both %FVs
            %     conds = {'FV_Na_Ca2_mGluR', 'nbqx_apv_ttx',  'synapticTransmission'}; % no cadmium
            %     conds = {'nbqx_apv_cd2_ttx', 'nbqx_apv_ttx',  'synapticTransmission'}; % both opsin current verisons
    end
    
    
    for i_ch = 1:2;
        
        if RESTRICT_TO_STIM_SITE && (i_ch ~= info{i_ex}.stimSite)
           continue
        end
        
        % deal with data and channels that need to be ignored.
        if ~info{i_ex}.ignoreChans(i_ch)
            continue
        else
            % strange cases
            mouseMatch = strcmpi(in{i_ex,1}, EXCEPTIONS(:,1));
            siteMatch = in{i_ex,2} == vertcat(EXCEPTIONS{:,2});
            if any(mouseMatch & siteMatch)
                % HS2 is the data channel, but since HS1 wasn't
                % used, the data are in the first column
                if i_ch == 1; error('something went wrong'); end
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
                
                % check to see if this condition exists
                if ~isfield(dat{i_ex}.(pTypes{i_tf}), 'stats') || ~isfield(dat{i_ex}.(pTypes{i_tf}).stats, conds{i_cond})
                    continue
                end
                
                tmp_raw = dat{i_ex}.(pTypes{i_tf}).snips.(conds{i_cond}){i_ch}';
                Ntime = size(tmp_raw,1);
                Ncols = size(tmp_raw,2);
                
                tt = (0:Ntime-1) ./ info{i_ex}.(pTypes{i_tf}).(conds{i_cond}).sampRate;
                tt = (tt - prePulseTime) .* 1000; % in ms.
                
                subplot(Nconds, Ntfs, i_tf+((i_cond-1) * Ntfs), 'align'),
                cmap = colormap('copper');
                cidx = round(linspace(1, size(cmap,1), max([Ncols, Ntfs])));
                cmap = cmap(cidx,:);
                set(gca, 'colororder', cmap, 'NextPlot', 'replacechildren');
                
                plot(tt, tmp_raw, 'linewidth', 2), hold on,
                
                if CHECK_TRIAL_STATS
                    
                    % check the peak and trough indicies
                    inds = dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).trpk_inds{i_ch};
                    inds_idx = bsxfun(@plus, inds, Ntime .* (0:Ncols-1)');
                    plot(tt(inds(:,1)), tmp_raw(inds_idx(:,1)), 'ro', 'markerfacecolor', 'r')
                    
                    if any(strcmpi(conds{i_cond}, {'FV_Na', 'FV_Na_Ca2_mGluR'}))
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
                    if any(strcmpi(conds{i_cond}, {'nbqx_apv_cd2_ttx', 'nbqx_apv_ttx'}))
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
        statTypes = {'diffval', 'area', 'pk2tr', 'slope', 'tau_m'};
        Nstats = numel(statTypes);
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
                                plot(xvals, repmat(yvals(i_plt,:), 2, 1), ':', 'color', cmap(i_plt,:), 'linewidth', 1);
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
PVALTYPE = 'diffval';

% only do this for the main experiment (TF and FV)
assert(strcmpi(EXPTTYPE, 'main expt'), 'ERROR: this anaysis is only for the TF and FV experiments');

%initialize the outputs
opsinTypes = {'chr2', 'chief_cit', 'chronos', 'chief_flx'};
conds = {'nbqx_apv_cd2_ttx', 'FV_Na', 'synapticTransmission'};
statTypes = {'diffval', 'area', 'pk2tr', 'slope', 'tau_m'};
for i_opsin = 1:numel(opsinTypes)
    for i_cond = 1:numel(conds)
        for i_stat = 1:numel(statTypes)
            pop.(opsinTypes{i_opsin}).pnp1.(conds{i_cond}).(statTypes{i_stat}) = {[],[],{}}; % {[TF], [pulseAmp], {Vals}}
            pop.(opsinTypes{i_opsin}).raw.(conds{i_cond}).(statTypes{i_stat}) = {[],[],{}}; % {[TF], [pulseAmp], {Vals}}
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
    mouseMatch = strcmpi(in{i_ex,1}, EXCEPTIONS(:,1));
    siteMatch = in{i_ex,2} == vertcat(EXCEPTIONS{:,2});
    if any(mouseMatch & siteMatch) % a strange exception that bucks the rules.
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
    Nptypes = numel(pTypes);
    for i_ptype = 1:Nptypes
        
        for i_cond = 1:numel(conds);
            
            % skip instances where the drug condition does not exist
            if ~isfield(dat{i_ex}.(pTypes{i_ptype}).stats, conds{i_cond})
                continue
            end
            
            for i_stat = 1:numel(statTypes)
                
                % skip instances where the stat type doesn't exist
                if ~isfield(dat{i_ex}.(pTypes{i_ptype}).stats.(conds{i_cond}), statTypes{i_stat})
                    continue
                end
                
                % structure the pnp1 data
                ex_stat = dat{i_ex}.(pTypes{i_ptype}).stats.(conds{i_cond}).(statTypes{i_stat}){CHANNEL};
                if strcmpi(statTypes{i_stat}, 'tau_m')
                    ex_stat = 1./ex_stat; % tau_m isn't a tau, need to convert to tau
                end
                
                % convert to paired pulse measures by normalizing by the
                % height of the first pulse.
                ex_p1p2 = ex_stat ./ ex_stat(1);
                
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
            cmap = colormap('copper');
            cidx = round(linspace(1, size(cmap,1), Ntfs));
            cmap = cmap(cidx,:);
            set(gca, 'colororder', cmap, 'NextPlot', 'replacechildren');
            hold on
            
            legtext = {};
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
                if strcmpi(conds{i_cond}, 'FV_Na')
                    nrepeats = sum(~isnan(cat_dat), 1)
                end
                
                errorbar(1:7, xbar, sem, 'color', cmap(i_tf,:), 'linewidth', 2) % only plot the first 7 pulses
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
            
            
%             if strcmpi(PVALTYPE, statTypes{i_stat}) && ~strcmpi(conds{i_cond}, 'synapticTransmission')
%                 %
%                 %  Do some inferential tests and present a table with the
%                 %  results
%                 %
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%
%                 [p_for_each_tf, p_p7_across_tfs, p_p2_across_tfs] = deal(nan(Ntfs, 1));
%                 
%                 % Test for significant decreases in PPRs. Look specifically at
%                 % P7:P1. Do the test individually for each TF
%                 %
%                 % Ho -> distribution of PPRs has Xbar = 1;
%                 for i_tf = 1:Ntfs
%                     tmp = pp_dat_allTFs{i_tf}(:,7);
%                     p_for_each_tf(i_tf, 1) = signrank(tmp-1);
%                 end
%                 
%                 
%                 % Test for differences in P7:P1 across TFs. Here I'm asking if
%                 % there's a significant effect of TF on the PPR. Do an anova
%                 % like test with post-hoc comparisons
%                 %
%                 % Ho -> distributions of PPRs have the identical median
%                 tmp_dat_p2 = [];
%                 tmp_dat_p7 = [];
%                 tmp_group = [];
%                 for i_tf = 1:Ntfs
%                     tmp_dat_p2 = cat(1, tmp_dat_p2, pp_dat_allTFs{i_tf}(:,2));
%                     tmp_dat_p7 = cat(1, tmp_dat_p7, pp_dat_allTFs{i_tf}(:,7));
%                     tmp_group = cat(1, tmp_group, ones(size(pp_dat_allTFs{i_tf}(:,7))).*i_tf);
%                 end
%                 p_p7_across_tfs(1) = kruskalwallis(tmp_dat_p7, tmp_group, 'off');
%                 
%                 
%                 % Test for differences in P2:P1 across TFs. Here I'm asking if
%                 % there's a significant effect of TF on the PPR. Do an anova
%                 % like test with post-hoc comparisons
%                 %
%                 % Ho -> distributions of PPRs have the identical median
%                 p_p2_across_tfs(1) = kruskalwallis(tmp_dat_p2, tmp_group, 'off');
%                 
%                 
%                 fprintf('  comparisons for %s values for %s  \n',...
%                     upper(opsinTypes{i_opsin}), upper(conds{i_cond}))
%                 
%                 rownames = cellfun(@num2str, mat2cell(tfs, ones(size(tfs))), 'uniformoutput', false);
%                 T = table(p_for_each_tf, p_p7_across_tfs, p_p2_across_tfs, 'RowNames', rownames)
%                 fprintf('\n\n')
%             end
            
        end
    end
end


%% POPULATION SUMMARY: INTERLEAVED POWERS


% plan: loop over opsins. Only consider a single recording channel (distal
% or proximal). Show PP ratio as a function of TF. Show Pn:P1 ratio as a
% function of pulse number for each frequency

STIMSITE = true;  % true => stimsite,  false => distal site
PULSENUM = 3;
STATTYPE = 'diffval';
CONDITION = 'FV_Na';

% only do this for the main experiment (TF and FV)
assert(strcmpi(EXPTTYPE, 'interleaved amps'), 'ERROR: this anaysis is only for experiments with interleaved amplitudes');

% load in the LED calibration data to convert from Volts to power denisty
cd(GL_DATPATH)
cd '../../Calibration files'
load led_cal.mat % data now stored in 'cal' struct


% initialize the figure
f = figure; hold on;

[rho_chr2, rho_ochief, rho_chronos] = deal([]);
for i_ex = 1:numel(dat)
    
    
    % skip cases where the led was not targeted to either of the recording
    % sites
    if isnan(info{i_ex}.stimSite)
        continue
    end
    
    % Determine which recording channel to analyze
    mouseMatch = strcmpi(in{i_ex,1}, EXCEPTIONS(:,1));
    siteMatch = in{i_ex,2} == vertcat(EXCEPTIONS{:,2});
    if any(mouseMatch & siteMatch) % a strange exception that bucks the rules.
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
    
    % get P2:P1 ratios for the FV
    ex_tf = [];
    ex_amps = [];
    ex_ppr = [];
    
    pTypes = fieldnames(dat{i_ex});
    for i_ptype = 1:numel(pTypes);
        
        if ~isfield(dat{i_ex}.(pTypes{i_ptype}).stats, CONDITION)
            continue
        end
        
        tmp = dat{i_ex}.(pTypes{i_ptype}).stats.(CONDITION).(STATTYPE){CHANNEL};
        if numel(tmp)<PULSENUM; continue; end;
        ppr = tmp(PULSENUM)./tmp(1);
        
        ex_ppr = cat(1, ex_ppr, ppr);
        ex_tf = cat(1, ex_tf, info{i_ex}.(pTypes{i_ptype}).none.pTF);
        ex_amps = cat(1, ex_amps, info{i_ex}.(pTypes{i_ptype}).none.pAmp);
        
        % convert amplitudes from volts to mW/mm2
        ex_powDensity = ppval(cal.pp.led_470, ex_amps);
        
        % run a correlation b/w PPR and light power
        rho = corr(ex_ppr, ex_powDensity, 'type', 'Pearson');
    end
    
    %
    % do some plotting if there are data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(ex_ppr)
        
        switch lower(info{i_ex}.opsin)
            case 'chr2'
                pltclr = 'b';
                rho_chr2 = cat(1, rho_chr2, rho);
            case 'chief_cit'
                pltclr = 'r';
                rho_ochief = cat(1, rho_ochief, rho);
            case 'chronos'
                pltclr = 'g';
                rho_chronos = cat(1, rho_chronos, rho);
        end
        
        assert(all(ex_tf == ex_tf(1)), 'ERROR: not all TFs equal')
        if ex_tf(1) == 10
            symbol = 'o';
        elseif ex_tf(1) == 20;
            symbol = 's';
        elseif ex_tf(1) == 50;
            symbol = 'v';
        else
            error('undefined tf')
        end
        
        
        figure(f); hold on,
        plot(ex_powDensity, ex_ppr, ['-', symbol, pltclr])
        
    end
    
end
%set(gca, 'yscale', 'log')
xlabel('Power (mW per mm^2)')
ylabel(sprintf('Paired Pulse Ratio \n (%s)', STATTYPE))
set(gca, 'yscale', 'log')






%% OPSIN CURRENT VS. FIBER VOLLEY [SPIDER PLOT]

% only do this for the main experiment (TF and FV)
assert(strcmpi(EXPTTYPE, 'main expt'), 'ERROR: this anaysis is only for the TF and FV experiments');


close all

FV_STAT = 'diffval';
OPSIN_STAT = 'diffval';
STAT_TYPE = 'raw';  % can be 'pnp1' or 'raw'
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
    tf_idx = pop.ochief.(STAT_TYPE).FV_Na.(FV_STAT){1} == TFs(i_tf);
    if sum(tf_idx) >= 1;
        trim_dat = cellfun(@(x) x(:,1:7), pop.ochief.(STAT_TYPE).FV_Na.(FV_STAT){3}(tf_idx), 'uniformoutput', false);
        ochief_fv_raw = vertcat(trim_dat{:});
        
        trim_dat = cellfun(@(x) x(:,1:7), pop.ochief.(STAT_TYPE).nbqx_apv_cd2_ttx.(OPSIN_STAT){3}(tf_idx), 'uniformoutput', false);
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

figure, hold on,
for i_ex = 1:numel(dat)
    
    
    % Determine which recording channel to analyze
    mouseMatch = strcmpi(in{i_ex,1}, EXCEPTIONS(:,1));
    siteMatch = in{i_ex,2} == vertcat(EXCEPTIONS{:,2});
    if any(mouseMatch & siteMatch) % a strange exception that bucks the rules.
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
    
    
    % there could be multiple TFs and pAmps. figure out what the deal is.
    pTypeNames = fieldnames(dat{i_ex});
    ex_pAmps = nan(numel(pTypeNames), 1);
    ex_TFs = nan(numel(pTypeNames), 1);
    for i_pType = 1:numel(pTypeNames);
        ex_pAmps(i_pType) = info{i_ex}.(pTypeNames{i_pType}).none.pAmp;
        ex_TFs(i_pType) = info{i_ex}.(pTypeNames{i_pType}).none.pTF;
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
    
    
    switch lower(info{i_ex}.opsin)
        case 'chr2'
            p = plot(opsin_avg, fv_avg, 'b-o');
        case 'chief_cit'
            p = plot(opsin_avg, fv_avg, 'r-o');
        case 'chronos'
            p = plot(opsin_avg, fv_avg, 'g-o');
    end
    
    % axis labels and such
    xlabel(sprintf('Opsin current (%s)', OPSIN_STAT));
    ylabel(sprintf('Fiber Volley (%s)', FV_STAT));
    bdfxn = @(x,y,z,zz) mytitle(sprintf('%s, site: %d', z, zz));
    set(p, 'buttonDownFcn', {bdfxn, in{i_ex,1}, in{i_ex,2}})
    if NORMTOMAX; ylim([0 1]); xlim([0 1]); end
end







%%  SHAPE OF THE FIRST PULSE RESPONSE


STIMSITE = true;

conds = {'nbqx_apv_cd2_ttx', 'FV_Na'};
chr2_examp = {[] []};
ochief_examp = {[] []};
chronos_examp = {[] []};
for i_ex = 1:Nexpts
    
    
    % skip cases where the led was not targeted to either of the recording
    % sites
    if isnan(info{i_ex}.stimSite)
        continue
    end
    
    % Determine which recording channel to analyze
    mouseMatch = strcmpi(in{i_ex,1}, EXCEPTIONS(:,1));
    siteMatch = in{i_ex,2} == vertcat(EXCEPTIONS{:,2});
    if any(mouseMatch & siteMatch) % a strange exception that bucks the rules.
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
                    normFact = dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).diffval{CHANNEL}(1);
                    normFact = normFact * -1;
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
            case 'chief_cit'
                ochief_examp{i_cond} = cat(1, ochief_examp{i_cond}, firstPulse);
            case 'chr2'
                chr2_examp{i_cond} = cat(1, chr2_examp{i_cond}, firstPulse);
            case 'chronos'
                chronos_examp{i_cond} = cat(1, chronos_examp{i_cond}, firstPulse);
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
    plot(tt, chronos_examp{i_cond}', 'g');
    plot(tt, mean(ochief_examp{i_cond}, 1), 'r', 'linewidth', 4);
    plot(tt, mean(chr2_examp{i_cond}, 1), 'b', 'linewidth', 4);
    plot(tt, mean(chronos_examp{i_cond}, 1), 'g', 'linewidth', 4);
    xlabel('time (msec)')
    ylabel('Normalized LFP amplitude')
    axis tight
    switch conds{i_cond}
        case 'FV_Na'
            mytitle('Average FV for first pulse')
        case 'nbqx_apv_cd2_ttx'
            mytitle('Average opsin current for first pulse')
    end
    
    
    
    if strcmpi(conds{i_cond}, 'FV_Na')
        figure, hold on,
        plot(tt, cumsum(-mean(ochief_examp{i_cond}, 1)), 'r', 'linewidth', 4);
        plot(tt, cumsum(-mean(chr2_examp{i_cond}, 1)), 'b', 'linewidth', 4);
        plot(tt, cumsum(-mean(chronos_examp{i_cond}, 1)), 'g', 'linewidth', 4);
        xlabel('time (msec)')
        ylabel('Normalized LFP amplitude')
        mytitle('Integral of (negative) FV pulse')
        axis tight
    end
    
end


%% SHAPE OF OPSIN CURRENT FOR EACH PULSE, WITH TAU FITS

% choose a stimulation location:
STIMSITE = true;
PLOTERR = false;
Npulses = 7;

clc; close all

% only do this for the main experiment (TF and FV)
%assert(strcmpi(EXPTTYPE, 'main expt'), 'ERROR: this anaysis is only for the TF and FV experiments');

% prepare the population structure
opsinTypes = {'chr2', 'chief_cit', 'chronos', 'chief_flx'};
tfnames = {'tf10', 'tf20', 'tf40', 'tf60', 'tf100'};
switch EXPTTYPE
    case '4-AP'
        drugconds = {'ttx_cd2', 'ttx_cd2_4AP200' 'ttx_cd2_4AP800', 'ttx_cd2_4AP1800'};
    case 'Baclofen'
        drugconds = {'ttx_cd2', 'ttx_cd2_bac10'};
    case 'otherwise'
        drugconds = {'nbqx_apv_cd2_ttx', 'FV_Na'};
end

for i_opsin = 1:numel(opsinTypes)
     opsincurrent.(opsinTypes{i_opsin}) = repmat({[]}, numel(tfnames), Npulses, numel(drugconds));
end


% organize and collect the data
Nexpt = size(in,1);
for a_ex = 1:Nexpt
    
    opsin = lower(info{a_ex}.opsin);
    for a_ptype = 1:numel(tfnames)
        
        % this pulse type might not exist. Detect these cases and continue
        tmp_tfcond = tfnames{a_ptype};
        if ~isfield(dat{a_ex}, tmp_tfcond)
            continue
        end
        
        % determine which channel to analyze
        % Determine which recording channel to analyze
        mouseMatch = strcmpi(in{a_ex,1}, EXCEPTIONS(:,1));
        siteMatch = in{a_ex,2} == vertcat(EXCEPTIONS{:,2});
        if any(mouseMatch & siteMatch) % a strange exception that bucks the rules.
            CHANNEL = 1;
        else % all other experiments...
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
            troughidx_col = dat{a_ex}.(tmp_tfcond).stats.(tmp_drugcond).trpk_inds{CHANNEL}(:,1);
            troughidx_row = 1:size(snips,1);
            troughidx = sub2ind(size(snips), troughidx_row(:), troughidx_col(:));
            normvals = snips(troughidx);
            normvals = abs(normvals); % maintains the sign of the opsin current
            normsnips = bsxfun(@rdivide, snips, normvals);
            
            
            % correct for differences in the sampling rate
            exptSampRate = info{a_ex}.(tmp_tfcond).(tmp_drugcond).sampRate;
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
                
                tmp = opsincurrent.(opsin){a_ptype, a_pulse, a_drug};
                tmp = cat(1, tmp, normsnips(a_pulse,:));
                opsincurrent.(opsin){a_ptype, a_pulse, a_drug} = tmp;
                
            end
        end
        
    end
    
end

cmap = colormap('copper'); close;
cidx = round(linspace(1, size(cmap,1), Npulses));
cmap = cmap(cidx,:);
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
                sem = std(tmpdat{a_pulse}, [], 1) ./ sqrt(size(tmpdat{a_pulse}, 1));
                tt = [0:numel(xbar)-1] ./ 20e3;
                if isempty(xbar)
                    continue
                end
                
                pltclr = cmap(a_pulse,:);
                
                % plot the SEM
                if PLOTERR
                    shadedErrorBar(tt, xbar, sem, {'-','color', pltclr, 'linewidth', 2})
                else
                    plot(tt, xbar, '-', 'color', pltclr, 'linewidth', 2);
                end
                
                
            end
            axis tight
        end
    end
end

%% AVERAGE RECORDING WAVEFORM OF OPSIN CURRENT FOR WHOLE SWEEP

% choose a stimulation location:
STIMSITE = true;
PLOTERR = false;

clc; close all

% only do this for the main experiment (TF and FV)
correctExpt = any(strcmpi(EXPTTYPE, {'main expt', '4-AP', 'baclofen'}));
assert(correctExpt, 'ERROR: This analysis is for a different type of experiment');

% prepare the population structure
opsinTypes = {'chr2', 'chief_cit', 'chronos', 'chief_flx'};
tfnames = {'tf10', 'tf20', 'tf40', 'tf60', 'tf100'};
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
        % this pulse type might not exist. Detect these cases and continue
        if ~isfield(dat{a_ex}, pTypes{a_ptype})
            continue
        end
        
        
        % find the TTX pharm condition
        conds = fieldnames(dat{a_ex}.(pTypes{a_ptype}).snips);
        ttxidx = cellfun(@any, regexpi(conds, 'ttx'));
        
        if ~any(ttxidx);
            continue
        end
        
        
        ttxidx = find(ttxidx);
        for i_ttx = ttxidx'; % so that I can compare multiple K+ drugs
            
            % determine which channel to analyze
            % Determine which recording channel to analyze
            mouseMatch = strcmpi(in{a_ex,1}, EXCEPTIONS(:,1));
            siteMatch = in{a_ex,2} == vertcat(EXCEPTIONS{:,2});
            if any(mouseMatch & siteMatch) % a strange exception that bucks the rules.
                CHANNEL = 1;
            else % all other experiments...
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
            pulse1_snip = dat{a_ex}.(pTypes{a_ptype}).snips.(conds{i_ttx}){CHANNEL}(1,:);
            troughidx = dat{a_ex}.(pTypes{a_ptype}).stats.(conds{i_ttx}).trpk_inds{CHANNEL}(1,1);
            normval = pulse1_snip(troughidx);
            normval = abs(normval); % this preserves the sign of the original signal
            raw = raw ./ normval;
            
            % figure out how many pulses there were, cut off everything past
            % the 7th pulse. If not 7 pulses, move along
            pulseOnIdx = find(info{a_ex}.(pTypes{a_ptype}).(conds{i_ttx}).pulseOn_idx);
            pulseOnIdx = pulseOnIdx(1);
            pulseOffIdx = find(info{a_ex}.(pTypes{a_ptype}).(conds{i_ttx}).pulseOff_idx);
            if numel(pulseOffIdx)<7; continue; end
            pulseOffIdx = pulseOffIdx(7);
            sampRate = info{a_ex}.(pTypes{a_ptype}).(conds{i_ttx}).sampRate;
            preTimeSec = 0.030;
            preTimeSamps = ceil(preTimeSec .* sampRate);
            raw = raw(pulseOnIdx - preTimeSamps : pulseOffIdx+preTimeSamps);
            
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
% 
% cmap = colormap('parula'); close;%also opens a figure;
% inds = round(linspace(1, 50, numel(tfnames)));
% cmap = cmap(inds,:);
% for a_opsin = 1:numel(opsinTypes)
%     figure, hold on,
%     
%     
%     for a_tf = 1:numel(tfnames)
%         
%         if isempty(avgopsincurent.(opsinTypes{a_opsin}).(tfnames{a_tf}).raw)
%             continue
%         end
%         
%         drugconds = fieldnames(avgopsincurent.(opsinTypes{a_opsin}).(tfnames{a_tf}).raw);
%         tfclr = cmap(a_tf,:);
%         drugcolors = repmat(tfclr, numel(drugconds), 1);
%         drugcolors = bsxfun(@times, drugcolors, linspace(1, 0.5, numel(drugconds))');
%         for a_drug = 1:numel(drugconds)
%             
%             tmpdat = avgopsincurent.(opsinTypes{a_opsin}).(tfnames{a_tf}).raw.(drugconds{a_drug});
%             if isempty(tmpdat); continue; end
%             
%             tmpmean = nanmean(tmpdat, 1);
%             tmpsem = nanstd(tmpdat, [], 1) ./ sqrt(size(tmpdat,1));
%             Nsamps = numel(tmpmean);
%             preTimeSamps = ceil(preTimeSec .* 20e3);
%             tt = ([0:Nsamps-1]-preTimeSamps) ./ 20e3;
%             if PLOTERR
%                 shadedErrorBar(tt, tmpmean, tmpsem, {'-','color', drugcolors(a_drug,:), 'linewidth', 2})
%             else
%                 plot(tt, tmpmean, '-', 'color',drugcolors(a_drug,:), 'linewidth', 2)
%             end
%         end
%         
%     end
% end



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
            case 'chief_cit'
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
            tmpsem = nanstd(tmpdat, [], 1) ./ sqrt(size(tmpdat,1));
            Nsamps = numel(tmpmean);
            preTimeSamps = ceil(preTimeSec .* 20e3);
            tt = ([0:Nsamps-1]-preTimeSamps) ./ 20e3;
            
            if PLOTERR
                shadedErrorBar(tt, tmpmean, tmpsem, {'-','color', pltclr(a_drug,:), 'linewidth', 2})
            else
                plot(tt, tmpmean, '-', 'color', pltclr(a_drug,:), 'linewidth', 2)
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
        tdict = outerleave(ax, ax.idx.LED_470);
        
        
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
                
                % strange cases
                mouseMatch = strcmpi(in{i_ex,1}, EXCEPTIONS(:,1));
                siteMatch = in{i_ex,2} == vertcat(EXCEPTIONS{:,2});
                if any(mouseMatch & siteMatch)
                    %  HS2 is the data channel, but since HS1 wasn't
                    %  used, the data are in the first column
                    if i_ch == 1; error('something went wrong'); end
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
            ledChIdx = strcmpi('LED_470', ax.head.recChNames);
            
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
                [troughidx, peakidx]  = anlyMod_getWFepochs(tmp, tt, conds{i_cond}, pwidth, photoDelay);
                
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
opsins = {'chr2', 'chief_cit', 'chronos'};

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
        mouseMatch = strcmpi(smoothStats{idx}.mouseName, EXCEPTIONS(:,1));
        siteMatch = smoothStats{idx}.siteNum == vertcat(EXCEPTIONS{:,2});
        if any(mouseMatch & siteMatch) % a strange exception that bucks the rules.
            CHANNEL = 1;
        else % all other experiments...
            if STIMSITE
                CHANNEL = smoothStats{idx}.stimSite;
            else
                if smoothStats{idx}.stimSite == 1
                    CHANNEL = 2;
                elseif smoothStats{idx}.stimSite == 2
                    CHANNEL = 1;
                end
            end
        end
        
        % plot each of the drug conditions
        for i_cond = 1:numel(conds)
            if ~isfield(smoothStats{idx}, conds{i_cond})
                continue
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
        mouseMatch = strcmpi(smoothStats{idx}.mouseName, EXCEPTIONS(:,1));
        siteMatch = smoothStats{idx}.siteNum == vertcat(EXCEPTIONS{:,2});
        if any(mouseMatch & siteMatch) % a strange exception that bucks the rules.
            CHANNEL = 1;
        else % all other experiments...
            if STIMSITE
                CHANNEL = smoothStats{idx}.stimSite;
            else
                if smoothStats{idx}.stimSite == 1
                    CHANNEL = 2;
                elseif smoothStats{idx}.stimSite == 2
                    CHANNEL = 1;
                end
            end
        end
        
        
        for i_cond = 1:Nconds
            if ~isfield(smoothStats{idx}, conds{i_cond})
                continue
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
            
        end
    end
end




%% LASER STIM SITE ANALYSIS


% figure out how many unique experiments there were
mouseNames = in(:,1);
sites = cellfun(@floor, in(:,2), 'uniformoutput', false); % use floor to eliminate the position information
comboname = cellfun(@(x,y) [x,num2str(y)], mouseNames, sites, 'uniformoutput', false);
[uniqueExpts, ~, exptID] = unique(comboname);


% iterate over the unique experiments and overlay the PPR metrics for FV,
% opsin current, and fEPSP slope
Nexpts = size(uniqueExpts, 1);
for i_ex = 1:Nexpts
    
    % order the files by proximal to distal (which should be hardcoded in
    % the excel workbook
    inds_mouse = find(exptID == i_ex);
    sites = cat(1, in{inds_mouse,2});
    [~, analysis_order] = sort(sites);
    inds_mouse = inds_mouse(analysis_order); % now the list is ordered by proximal to distal
    
    figure
    set(gcf, 'defaulttextinterpreter', 'none', 'position', [17 278 1353 420]);
    set(gcf, 'name', sprintf('%s, site: %d', in{inds_mouse(1),1}, floor(in{inds_mouse(1),2})))
    cmap = colormap('parula'); clf;%also opens a figure;
    inds = round(linspace(1, 50, numel(inds_mouse)));
    cmap = cmap(inds,:);
    
    % assume that the primary channel is the one that's in L2/3
    CH_proximal = structcat(info(inds_mouse), 'stimSite');
    CH_proximal = unique(cell2mat(CH_proximal));
    assert(numel(CH_proximal)==1, 'ERROR: too many ''proximal'' sites');
    if CH_proximal == 1
        CH_distal = 2;
    else
        CH_distal = 1;
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
        
        % insert tab group here, initialized the x,y coordinates to the
        % stimulus positions
        
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
                
                if ~isfield(dat{idx}.(tfconds{i_tf}).stats, conds{i_cond})
                    continue
                end
                
                % figure for distal channel goes on top, figure for proximal channel goes on bottom
                for i_ch = 1:2;
                    switch i_ch
                        case 1
                            CHANNEL = CH_proximal;
                        case 2
                            CHANNEL = CH_distal;
                    end
                    
                    if (CHANNEL == 2) && numel(dat{idx}.(tfconds{i_tf}).stats.(conds{i_cond}).(STAT)) == 1
                        continue % the data don't exist; that channel was not recorded
                    end
                    
                    
                    tmp_raw = dat{idx}.(tfconds{i_tf}).stats.(conds{i_cond}).(STAT){CHANNEL};
                    tmp_sigma = dat{idx}.(tfconds{i_tf}).stats.(conds{i_cond}).bkgnd_sigma{CHANNEL};
                    if ~any(regexpi(conds{i_cond}, 'ttx'))
                        tmp_sigma = tmp_sigma ./ tmp_raw(1); % needs to be b/4 the next line b/c tmp_raw is destructively modified
                        tmp_raw = tmp_raw ./ tmp_raw(1);
                    end
                    
                    subplot(2, npltcols, (i_ch-1)*npltcols + i_cond)
                    hold on,
                    plot(1:numel(tmp_raw), tmp_raw, 'o-', 'color', cmap(i_site,:), 'linewidth', 2)
                    xvals = get(gca, 'xlim');
                    if any(regexpi(conds{i_cond}, 'FV_Na'))
                        yvals = [1-tmp_sigma, 1+tmp_sigma];
                        plot(xvals, repmat(yvals, 2, 1), ':', 'color', cmap(i_site,:), 'linewidth', 1);
                    end
                    if i_ch == 1
                        if any(regexpi(conds{i_cond}, '_ttx'))
                            title(info{inds_mouse(1)}.opsin)
                        else
                            title(conds{i_cond})
                        end
                    end
                    % adjust the axes
                    switch conds{i_cond}
                        case {'FV_Na', 'synapticTransmission'}
                            yvals = get(gca, 'ylim');
                            set(gca, 'ylim', 1.05.* [0, max([yvals(2), max(tmp_raw)])]);
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
            title(sprintf('(0,0) is HS%d', CH_proximal))
        end
    end
    
    
    % add an icon for the stimulus locations
    centPos = round(ginput(1));
    xy = cell2mat(structcat(info(inds_mouse), 'optStimCords'));
    pixperum = pixPerMicron(size(img,1), size(img,2));
    xy = round(xy .* pixperum); %now in pix
    xy = bsxfun(@plus, xy, centPos); % pix relative to neuron
    hold on,
    for a = 1:size(xy,1)
        plot(xy(a,1), xy(a,2), 'o', 'markeredgecolor', cmap(a,:), 'markerfacecolor', cmap(a,:), 'markersize', 10)
    end
    
    % update plots in quasi realtime
    drawnow
    
end



%% COMPARE PLASTICITY ACROSS BRAIN AREAS

fin


ERRBARS = true;
NPULSES = 7;

% load in data from each area
load('C:\Users\charlie\Dropbox\Duke on Dropbox\oChIEF Methods\pop_al.mat');
al_pop = pop;

load('C:\Users\charlie\Dropbox\Duke on Dropbox\oChIEF Methods\pop_pm.mat');
pm_pop = pop;


%
% make a 3 x N plot. one row for each opsin, one column for each TF. Plot
% PPRs as a function of pulse number
%

F = figure;
set(F, 'position', [109    31   891   754])

opsinTypes = {'chr2', 'chief_cit', 'chronos'};
for i_opsin = 1:numel(opsinTypes)
    
    % grab the data
    pm_tmp_dat = pm_pop.(opsinTypes{i_opsin}).pnp1.synapticTransmission.slope;
    al_tmp_dat = al_pop.(opsinTypes{i_opsin}).pnp1.synapticTransmission.slope;
    
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
            errorbar(1:NPULSES, xbar, sem, '-', 'color', pltclr, 'linewidth', 2)
        else
            plot(1:NPULSES, xbar, '-', 'color', pltclr, 'linewidth', 2)
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
            errorbar(1:NPULSES, xbar, sem, '-', 'color', pltclr, 'linewidth', 2)
        else
            plot(1:NPULSES, xbar, '-', 'color', pltclr, 'linewidth', 2)
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
        legend({pm_leg, al_leg}, 'location', 'best')
        legend boxoff
    end
end











%% GET TRIAL COUNTS AND SNR AND GENOTYPES

clc

mdb = initMouseDB('update', 'notext');
nexpts = numel(info);
minExTrialCount = [];
SNR = [];
genotypes = {};
sex = {};
for a_ex = 1:nexpts
    
    [~, idx] = mdb.search(info{a_ex}.mouse);
    genotypes{a_ex, 1} = mdb.mice{idx}.info.strain;
    sex{a_ex, 1} = mdb.mice{idx}.info.sex;
    
    fields = fieldnames(info{a_ex});
    
    fieldTcount = inf;
    for a_field = 1:numel(fields)
        
        % make sure you're looking at data and not header info
        if ~isfield(info{a_ex}.(fields{a_field}) , 'none')
            continue
        end
        
        conds = {'nbqx_apv_cd2', 'nbqx_apv_cd2_ttx'};
        
        for a_cond = 1:numel(conds)
            
            if ~isfield(info{a_ex}.(fields{a_field}), conds{a_cond});
                error('could not find data')
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
            % strange cases
            mouseMatch = strcmpi(in{a_ex,1}, EXCEPTIONS(:,1));
            siteMatch = in{a_ex,2} == vertcat(EXCEPTIONS{:,2});
            if any(mouseMatch & siteMatch)
                % HS2 is the data channel, but since HS1 wasn't
                % used, the data are in the first column
                if CHANNEL == 1; error('something went wrong'); end
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



