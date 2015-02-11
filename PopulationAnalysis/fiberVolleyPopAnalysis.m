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
    'CH_150112_C', 1};


%% WHICH MICE SHOULD CONTRIBUTE?  [GOOD MICE]

% Anlyze data sets that used 300us pulses, and that have a pure Na+ FV

% clear out the workspace
fin

% in = {Mouse Name, Site}

in = {'CH_141215_F', 1;...
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
    'CH_150119_C', 2};


%% WHICH MICE SHOULD CONTRIBUTE?  [GOOD MICE FOR DIFFERENT PULSE AMP EXPTS]

% Anlyze data sets that used 300us pulses, and that have a pure Na+ FV

% clear out the workspace
fin

% in = {Mouse Name, Site}

in = {'CH_150112_C', 2};


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
        
        conds = {'FV_Na', 'nbqx_apv_cd2_ttx', 'synapticTransmission'};
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


%% ANALYZE THE SNIPPETS: PEAK-TO-PEAK, PEAK-2-BASELINE

% REMINDER: this is happening in separate cell-script because I will define
% the analysis region based off the average first pulse across TF conds.
% This requires grabbing all the snippets b/4 any analysis can proceed. The
% preceeding cell will need to be run before this one.

% should the analysis window be defined based off the average response
% following the first pulse (across conditions) or should the analysis
% window be unique for each pulse?
FIRSTPULSE = true;

for i_ex = 1:Nexpts
    
    conds = {'FV_Na', 'nbqx_apv_cd2_ttx'};
    for i_cond = 1:numel(conds)
        
        pTypes = fieldnames(dat{i_ex});
        Ntfs = numel(pTypes);
        
        sampRate = info{i_ex}.(pTypes{1}).(conds{i_cond}).sampRate;
        prePulseSamps = ceil(prePulseTime .* sampRate); % samples prior to pulse onset
        postPulseSamps = ceil(postPulseTime .* sampRate); % samples available after pulse ONSET
        photoDelaySamps = ceil(0 .* sampRate); % timeout following pulse offset
        
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
            
            
            
            if FIRSTPULSE
                % calculate the analysis windows based off the average 1st pulse,
                % which should be the same across TFs (within pharmacology
                % condition and channel).
                firstPulse = nan(Ntfs, prePulseSamps+postPulseSamps+1);
                for i_tf = 1:Ntfs
                    
                    firstPulse(i_tf,:) = dat{i_ex}.(pTypes{i_tf}).snips.(conds{i_cond}){i_ch}(1,:);
                    
                end
                pWidth = info{i_ex}.(pTypes{i_tf}).(conds{i_cond}).pWidth;
                pWidthSamps = ceil(pWidth .* sampRate);
                
                firstPulse = mean(firstPulse,1);
                
                firstValidPostPulseIdx = prePulseSamps + pWidthSamps + photoDelaySamps;
                [~, troughidx] = min(firstPulse(firstValidPostPulseIdx+1:end)); % only look after the pulse has come on
                troughidx = troughidx + firstValidPostPulseIdx;
                [~, peakidx] = max(firstPulse(firstValidPostPulseIdx+1:end));
                peakidx = peakidx + firstValidPostPulseIdx;
                
                % some error checking for the FV case
                if strcmpi('FV_Na', conds{i_cond})
                    assert(peakidx > troughidx, 'ERROR: negativity does not lead the positivity')
                    assert(peakidx./sampRate < 0.007, 'ERROR: positivity occurs too late');
                    assert(peakidx./sampRate > 250e-6, 'ERROR: positivity occurs too early');
                end
                
                % add a few points on either side of the true trough/peak
                troughidx = troughidx-2:troughidx+2;
                peakidx = peakidx-2:peakidx+2;
            end
            
            
            
            % now do the analysis
            for i_tf = 1:Ntfs
                
                Npulses = sum(info{i_ex}.(pTypes{i_tf}).(conds{i_cond}).pulseOn_idx);
                pOnIdx = find(info{i_ex}.(pTypes{i_tf}).(conds{i_cond}).pulseOn_idx);
                for i_pulse = 1:Npulses
                    
                    snippet = dat{i_ex}.(pTypes{i_tf}).snips.(conds{i_cond}){i_ch}(i_pulse,:);
                    
                    % store some stats for each pulse
                    switch conds{i_cond}
                        case 'FV_Na'
                            
                            trough = mean(snippet(troughidx));
                            peak = mean(snippet(peakidx));
                            
                            pk2tr = peak-trough;
                            dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).pk2tr{i_ch}(i_pulse) = pk2tr;
                            dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).diffval{i_ch}(i_pulse) = trough;
                            
                        case 'nbqx_apv_cd2_ttx'
                            
                            trough = mean(snippet(troughidx));
                            
                            dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).diffval{i_ch}(i_pulse) = trough;
                            
                    end
                    
                    % calculate the integral of the LFP signal. This
                    % integral will get adjusted later to reflect the
                    % integral of a noisy signal with no response to the
                    % LED...
                    pWidth = info{i_ex}.(pTypes{i_tf}).(conds{i_cond}).pWidth;
                    pWidthSamps = ceil(pWidth .* sampRate);
                    firstValidPostPulseIdx = prePulseSamps + pWidthSamps + photoDelaySamps;
                    snippet_pulse = snippet(firstValidPostPulseIdx:end);
                    
                    area = sum(abs(snippet_pulse)) ./ (numel(snippet_pulse) ./sampRate);
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

conds = {'nbqx_apv_cd2_ttx', 'FV_Na'};

for i_ex = 1:Nexpts
    
    for i_ch = 1:2;
        
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
        
        figure
        figName = sprintf('%s: site %d, ch %d', info{i_ex}.mouse, in{i_ex,2}, i_ch);
        set(gcf, 'name', figName,...
            'position', [148    33   813   752],...
            'defaulttextinterpreter', 'none');
        
        for i_cond = 1:numel(conds)
            
            %
            % plot the raw (Average) traces
            %
            pTypes = fieldnames(dat{i_ex});
            Ntfs = numel(pTypes);
            Nplts = max([3, Ntfs]);
            for i_tf = 1:Ntfs
                
                tmp_raw = dat{i_ex}.(pTypes{i_tf}).snips.(conds{i_cond}){i_ch};
                
                tt = (0:size(tmp_raw,2)-1) ./ info{i_ex}.(pTypes{i_tf}).(conds{i_cond}).sampRate;
                tt = (tt - prePulseTime) ./ 1000;
                
                pltidx = (i_cond-1)*Nplts + i_tf;
                subplot(3, Nplts, pltidx)
                cmap = colormap('copper');
                cidx = round(linspace(1, size(cmap,1), size(tmp_raw,1)));
                cmap = cmap(cidx,:);
                set(gca, 'colororder', cmap, 'NextPlot', 'replacechildren');
                
                plot(tt, tmp_raw', 'linewidth', 2)
                axis tight
                title(pTypes{i_tf})
                xlabel('time (ms)')
                if i_tf==1;
                    switch conds{i_cond}
                        case 'nbqx_apv_cd2_ttx'
                            ylabel(sprintf('%s Current', info{i_ex}.opsin))
                        case 'FV_Na'
                            ylabel('Fiber Volley')
                    end
                    
                end
            end
        end
        
        
        %
        % plot the summary stat (pk2tr or diff from baseline) as a function
        % of pulse number
        %
        
        % diff val for opsin current
        subplot(3,Nplts, Nplts*2+1), hold on,
        for i_tf = 1:Ntfs
            
            diffval = dat{i_ex}.(pTypes{i_tf}).stats.nbqx_apv_cd2_ttx.diffval{i_ch};
            diffval = abs(diffval);
            plot(1:numel(diffval), diffval, 'o-', 'color', cmap(i_tf,:), 'linewidth', 2)
        end
        xlabel('Pulse number')
        ylabel(sprintf('%s diffVal', info{i_ex}.opsin))
        xlim([1, numel(diffval)])
        yvals = get(gca, 'ylim');
        yvals(1) = min([0, yvals(1)]);
        set(gca, 'ylim', [0, yvals(2)]);
        if yvals(1)<0
            plot([1,numel(diffval)], [0,0] , 'k--', 'linewidth', 2)
        end
        
        
        % Fiber volley peak to trough
        subplot(3,Nplts, Nplts*2+2), hold on,
        for i_tf = 1:Ntfs
            
            pk2tr = dat{i_ex}.(pTypes{i_tf}).stats.FV_Na.pk2tr{i_ch};
            plot(1:numel(pk2tr), pk2tr, 'o-', 'color', cmap(i_tf,:), 'linewidth', 2)
        end
        xlabel('Pulse number')
        ylabel('Fiber Volley pk2tr')
        xlim([1, numel(pk2tr)])
        yvals = get(gca, 'ylim');
        yvals(1) = min([0, yvals(1)]);
        set(gca, 'ylim', yvals);
        if yvals(1)<0
            plot([1,numel(diffval)], [0,0] , 'k--', 'linewidth', 2)
        end
        
        % fiber volley integral
        subplot(3,Nplts, Nplts*2+3), hold on,
        for i_tf = 1:Ntfs
            
            area = dat{i_ex}.(pTypes{i_tf}).stats.FV_Na.area{i_ch};
            plot(1:numel(area), area, 'o-', 'color', cmap(i_tf,:), 'linewidth', 2)
        end
        xlabel('Pulse number')
        ylabel('FV Area (mV/sec)')
        xlim([1, numel(area)])
        yvals = get(gca, 'ylim');
        yvals(1) = min([0, yvals(1)]);
        set(gca, 'ylim', yvals);
        if yvals(1)<0
            plot([1,numel(area)], [0,0] , 'k--', 'linewidth', 2)
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
conds = {'nbqx_apv_cd2_ttx', 'FV_Na'};
stattype = {'diffval', 'pk2tr', 'area'};
for i_opsin = 1:numel(opsinTypes)
    for i_cond = 1:numel(conds)
        for i_stat = 1:numel(stattype)
            pop.(opsinTypes{i_opsin}).pnp1.(conds{i_cond}).(stattype{i_stat}) = {[],{}}; % {{TF}, {Vals}}
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
            
            stattype = fieldnames(dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}));
            
            for i_stat = 1:numel(stattype)
                
                % structure the pnp1 data
                tmp_stat = dat{i_ex}.(pTypes{i_tf}).stats.(conds{i_cond}).(stattype{i_stat}){CHANNEL};
                tmp_pnp1 = tmp_stat ./ tmp_stat(1);
                
                % grab the data field of the 'pop' structure
                tmp_pop = pop.(opsin).pnp1.(conds{i_cond}).(stattype{i_stat});
                tf = info{i_ex}.(pTypes{i_tf}).(conds{i_cond}).pTF;
                
                % is there already an entry for this TF and drug condition?
                alreadyThere = any(tmp_pop{1} == tf);
                if alreadyThere
                    idx = tmp_pop{1} == tf;
                    
                    % deal with cases where there are different numbers of
                    % pulses
                    nPulsesExisting = size(tmp_pop{2}{idx},2);
                    if numel(tmp_pnp1) < nPulsesExisting
                        tmp_pnp1(1,end+1:nPulsesExisting) = nan;
                    else
                        tmp_pop{2}{idx}(:,end+1:numel(tmp_pnp1)) = nan;
                    end
                    
                    tmp_pop{2}{idx} = cat(1, tmp_pop{2}{idx}, tmp_pnp1);
                else
                    tmp_pop{1} = cat(1, tmp_pop{1}, tf);
                    tmp_pop{2} = cat(1, tmp_pop{2}, tmp_pnp1);
                end
                
                % replace the 'pop' structure version with the tmp version
                pop.(opsin).pnp1.(conds{i_cond}).(stattype{i_stat}) = tmp_pop;
                
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
    
    
    for i_cond = 1:numel(conds);
        
        stattype = fieldnames(pop.(opsinTypes{i_opsin}).pnp1.(conds{i_cond}));
        
        for i_stat = 1:numel(stattype)
            
            % grab the data
            tmp_dat = pop.(opsinTypes{i_opsin}).pnp1.(conds{i_cond}).(stattype{i_stat});
            
            % make sure there's actually data there
            if isempty(tmp_dat{1}); continue; end % some things are not in the data set
            
            % now do the plotting
            subplot(2, numel(stattype), (i_cond-1).*numel(stattype) + i_stat)
            title(stattype{i_stat})
            
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
                
                errorbar(1:numel(xbar), xbar, sem, 'color', cmap(i_tf,:), 'linewidth', 2)
                
                legtext = cat(2, legtext, num2str(tfs(tf_idx)));
            end
            axis tight
            xlabel('Pulse number')
            if i_stat==1
                ylabel(sprintf('%s \n Pn:P1 ratio', conds{i_cond}))
            else
                ylabel('Pn:P1 ratio')
            end
            if i_stat==1  && i_cond==1
                legend(legtext, 'location', 'southeast')
                legend boxoff
            end
            yvals = get(gca, 'ylim');
            yvals(1) = min([0, yvals(1)]);
            set(gca, 'ylim', yvals);
            if yvals(1)<0
                plot([1,numel(diffval)], [0,0] , 'k--', 'linewidth', 2)
            end
            
        end
    end
end




%%  ESTIMATE THE SLOPE OF THE FIELD EPSP


STIMSITE = true;  % true => stimsite,  false => distal site

figure, hold on,
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
    
    
    conds = {'synapticTransmission'}; % 'for loop' of one for now, but room to grow...
    for i_cond = 1:numel(conds)
        
        pTypes = fieldnames(dat{i_ex});
        Ntfs = numel(pTypes);
        
        sampRate = info{i_ex}.(pTypes{1}).(conds{i_cond}).sampRate;
        prePulseSamps = ceil(prePulseTime .* sampRate); % samples prior to pulse onset
        postPulseSamps = ceil(postPulseTime .* sampRate); % samples available after pulse ONSET
        photoDelaySamps = ceil(0 .* sampRate); % 500us timeout following pulse offset
        
        
        
        firstPulse = nan(Ntfs, prePulseSamps+postPulseSamps+1);
        for i_tf = 1:Ntfs
            
            firstPulse(i_tf,:) = dat{i_ex}.(pTypes{i_tf}).snips.(conds{i_cond}){CHANNEL}(1,:);
            
        end
        firstPulse = mean(firstPulse,1);
        firstValidPostPulseIdx = prePulseSamps + pWidthSamps + photoDelaySamps;
        
        
        switch conds{i_cond}
            case 'synapticTransmission'
                
                % Analyze a time window that's approx 2.5 to 4 ms from
                % LED pulse onset
                windowStartIdx = prePulseSamps + ceil(0.0025 .* sampRate);
                windowStopIdx = prePulseSamps + ceil(0.004 .* sampRate);
                avgVal = mean(firstPulse(windowStartIdx:windowStopIdx));
                
                p1stats{i_ex}.(conds{i_cond}).avgval(1, CHANNEL) = abs(avgVal);
        end
        
        
        % provisional plotting (delete later)
        fPSP = p1stats{i_ex}.(conds{i_cond}).avgval(1, CHANNEL);
        pk2tr = dat{i_ex}.(pTypes{1}).stats.FV_Na.pk2tr{CHANNEL};
        pnp1 = pk2tr(end)./pk2tr(1);
        switch lower(info{i_ex}.opsin)
            case 'ochief'
                plot(fPSP, pnp1, 'ro')
            case 'chr2'
                plot(fPSP, pnp1, 'bo')
        end
        ylabel('diffval')
        xlabel('fPSP')
        
        
    end % conditions
    
end % expts



%%  SHAPE OF THE FIRST FIBER VOLLEY PULSE FOR oChIEF AND ChR2

close all
STIMSITE = true;

chr2_examp = [];
ochief_examp = [];
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
    
    
    conds = {'FV_Na'}; % for loop of one for now, but room to grow...
    for i_cond = 1:numel(conds)
        
        pTypes = fieldnames(dat{i_ex});
        Ntfs = numel(pTypes);
        
        sampRate = info{i_ex}.(pTypes{1}).(conds{i_cond}).sampRate;
        prePulseSamps = ceil(prePulseTime .* sampRate); % samples prior to pulse onset
        postPulseSamps = ceil(postPulseTime .* sampRate); % samples available after pulse ONSET
        photoDelaySamps = ceil(0 .* sampRate); % 500us timeout following pulse offset
        
        
        firstPulse = nan(Ntfs, prePulseSamps+postPulseSamps+1);
        for i_tf = 1:Ntfs
            
            firstPulse(i_tf,:) = dat{i_ex}.(pTypes{i_tf}).snips.(conds{i_cond}){CHANNEL}(1,:);
            
        end
        firstPulse = mean(firstPulse,1);
        normFact = dat{i_ex}.(pTypes{i_tf}).stats.FV_Na.pk2tr{CHANNEL}(1);
        firstPulse = firstPulse ./ normFact;
        
        
        % store the first pulse according to the opsin type, but only
        % if the sampling rate = 20kHz (so that the sampling latice is
        % the same...)
        correctSR = info{i_ex}.(pTypes{i_tf}).(conds{i_cond}).sampRate == 20e3;
        correctCond = strcmpi(conds{i_cond}, 'FV_Na');
        if correctSR && correctCond
            switch lower(info{i_ex}.opsin)
                case 'ochief'
                    ochief_examp = cat(1, ochief_examp, firstPulse);
                case 'chr2'
                    chr2_examp = cat(1, chr2_examp, firstPulse);
            end
        end
        
        
    end % conditions
    
end % expts


% the example first pulse FV
tt = (0:size(ochief_examp,2)-1)./20e3;
tt = tt-tt(prePulseSamps);
tt = tt.*1000;
figure, hold on,
plot(tt, ochief_examp', 'm');
plot(tt, chr2_examp', 'c');
plot(tt, mean(ochief_examp, 1), 'r', 'linewidth', 4);
plot(tt, mean(chr2_examp, 1), 'b', 'linewidth', 4);
xlabel('time (msec)')
ylabel('Normalized LFP amplitude')
title('Average FV for first pulse')
axis tight



figure, hold on,
plot(tt, cumsum(-mean(ochief_examp, 1)), 'r', 'linewidth', 4);
plot(tt, cumsum(-mean(chr2_examp, 1)), 'b', 'linewidth', 4);
xlabel('time (msec)')
ylabel('Normalized LFP amplitude')
title('Integral of (negative) FV pulse')
axis tight






