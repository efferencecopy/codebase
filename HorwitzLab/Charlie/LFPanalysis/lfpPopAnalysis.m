%% SETTING THINGS UP

% to do
% 2) Change xcorr to only 32 comparisons, yoked to the distance between each
%    etrode and the "closest" etrode to the stimulus.
% 
% 3) Calculate "mean" response locked to stim onset, and plot as a function
%    etrode distance to stimulus.
%
% 4) ID bad etrodes on the basis of %change

fin

% expts = {{'A051611001'},...
%          {'A051611002'}}; % bank?, position = 'center'

% expts = {{'F032613001'},...
%          {'F032613002'}}; % bank = 'b', position = 'ul'

% expts = {{'F032813002', 'F032813004', 'F032813006', 'F032813008'},...
%          {'F032813003', 'F032813005', 'F032813007'}}; % bank = c, position = 'up'

expts = {{'F032913001', 'F032913003', 'F032913005'},...
         {'F032913002', 'F032913004', 'F032913006'}}; % bank = c, position = 'right'




% params for the lfp average
global params
params.expts = expts;
params.freqBand = [5 20];
params.pretime = 0.500;
params.posttime = 6.00;
params.zscore = false;
params.norm2bkgnd = false;
params.absval = true;
params.trig = 'stimon';
params.bank = 'b';
params.stimPos = 'right';
params.orderDist = 'arbitrary';

%% PSTH OF ALL CHANNELS SYNCED TO STIM ONSET

clear,
global params

for ex = 1:numel(cellfun(@iscell, params.expts))
    
    % open the file
    tmp = {};
    if numel(params.expts{ex}) == 1
        WN = wnobj(params.expts{ex});
    elseif numel(params.expts{ex}) > 1
        for i = 1:numel(params.expts{ex})
            tmp{i} = wnobj(params.expts{ex}{i});
        end
        WN = strocat(tmp);
    end
    
    % a few preliminary definitions
    switch params.trig
        case 'stimon'
            TRIG = WN.idx.stim_on;
    end    
    
    % determine the color direction
    if all(WN.trial(:, WN.idx.noise_type) == 3)
        colordir = 'L-M';
    elseif all(all([WN.trial(:, [WN.idx.sigma1, WN.idx.sigma2])==0, WN.trial(:, WN.idx.sigma3)>0]))
        colordir = 'S-iso';
    else
        %error('unknown color')
        colordir = 'unknown';
    end
    
    lfpchan = find(WN.idx.lfp);
    timeVec = -params.pretime:1/WN.sum.analog.storeRates{1}:params.posttime;
    lfpaverages = nan(numel(lfpchan), numel(timeVec));
    avgVar = nan(numel(lfpchan),1);
    for chan = lfpchan
        
        % calculate the band limited power. 'blplfp' takes care of the oob
        % issues in the raw lfp.
        blp = blplfp(WN, chan, params.freqBand);
        
        % now find the average variance for each electrode
        rawvar = cellfun(@var, blp);
        l_valid = ~WN.checkTrials(chan);
        avgVar((chan == lfpchan),1) = mean(rawvar(l_valid));
        
        
        % zscore if desired
        if params.zscore
            blp = cellfun(@zscore, blp, 'uniformoutput', 0);
        end
        
        % rectify if desired
        if params.absval
            blp = cellfun(@abs, blp, 'uniformoutput', 0);
        end
        
        % pull out the time points surrounding the triger event
        anlgRate = WN.sum.analog.storeRates{chan==find(WN.idx.lfp)};
        preSamples = params.pretime * anlgRate;
        postSamples = params.posttime * anlgRate;
        psth = nan(numel(blp), preSamples+postSamples+1);
        for trl = 1:numel(blp)
            trigTime = WN.trial(trl, TRIG);
            anlgStart = WN.ras{trl, WN.idx.anlgStart};
            idxAtTrig = round((trigTime-anlgStart)*anlgRate + 1);
            if (idxAtTrig-preSamples)>=1 && (idxAtTrig+postSamples)<=numel(blp{trl});
                psth(trl,:) = blp{trl}(idxAtTrig-preSamples:idxAtTrig+postSamples);
            end
        end
        
        
        % normalize to the bkgnd (pre trig) resp
        if params.norm2bkgnd
            bkgnd = mean(psth(:, 1:preSamples),2);
            psth = bsxfun(@rdivide, psth, bkgnd);
        end
        
        % average the PSTHs
        lfpaverages((chan == lfpchan), :) = nanmean(psth,1);
        
        % do a t-test to see 
        pre = lfpaverages((chan == lfpchan), 1:preSamples);
        post = lfpaverages((chan == lfpchan), preSamples+1:end);
        P((chan == lfpchan), 1) = ranksum(pre, post);
    end
    
    % order the figures by the distance b/w each etrode and the stimulus.
    switch params.orderDist
        case 'none'
            % do nothing
        case 'arbitrary'
            [idx, relDist] = reorderPinDistance(params.stimPos, params.bank);
            lfpaverages = lfpaverages(idx,:);
            avgVar = avgVar(idx);
            relDist = relDist(idx);
    end
    
    
    figure, hold on,
    set(gcf, 'name', colordir)
    imagesc(lfpaverages)
    [~, idx] = min(abs(timeVec));
    plot([idx, idx], [1,numel(lfpchan)], 'k')
    xlim([0 numel(timeVec)]);
    ylim([0.5 numel(lfpchan)+0.5]);
    set(gca, 'xticklabel', [get(gca, 'xtick')./WN.sum.analog.storeRates{1}-params.pretime])
    set(gcf, 'position', [262   535   710   294]);
    colormap(pmkmp(250, 'CubicL'))
    colorbar
    set(gca, 'fontsize', 16)
    xlabel('Time (sec)')
    ylabel('Channel')
    if params.zscore
        title('PSTH on Zscores')
    else
        title('PSTH amp relative to bkgnd')
    end
    
    
    figure, hold on,
    set(gcf, 'name', colordir)
    pltClr = colormap(pmkmp(size(lfpaverages,1), 'CubicL'));
    for a = 1:size(lfpaverages,1);
        plot(timeVec, lfpaverages(a,:), '-', 'linewidth', 2, 'color', pltClr(a,:))
    end
    colorbar
    
    figure
    set(gcf, 'name', colordir)
    plot(timeVec, mean(lfpaverages,1))
    xlabel('Time (sec)')
    ylabel('Average LFP across channels')
    
    figure
    set(gcf, 'name', colordir)
    plot(relDist, avgVar, 'k.')
    xlabel('Relative pin distance')
    ylabel('Average variance')
end


%% CORRELATION BETWEEN CHANNELS


clear, 
global params %#ok<*REDEF>

MAXXCTIME = 0.030; % in seconds
for ex = 1:numel(cellfun(@iscell, params.expts))
    
    % open the file(s)
    tmp = {};
    if numel(params.expts{ex}) == 1
        WN = wnobj(params.expts{ex});
    elseif numel(params.expts{ex}) > 1
        for i = 1:numel(params.expts{ex})
            tmp{i} = wnobj(params.expts{ex}{i});
        end
        WN = strocat(tmp);
    end
    
    
    % for each trial, compute the xcorr for all possible combinations of
    % channels. Do this for all the available data, but at the end, only
    % average across correllelegrames that have enough data
    nTrials = size(WN.trial,1);
    lfpchan = find(WN.idx.lfp);
    LFPsamprate = WN.sum.analog.storeRates{lfpchan(1)};
    % iterate over the channels, making a mega array of blp for each channel
    % and trial
    blLFP = {};
    for chan = lfpchan
        blLFP(:, chan==lfpchan) = blplfp(WN, chan, params.freqBand);
    end
    
    % estimate the variance of each LFP signal and average within channels
    varLFP{ex} = cellfun(@var, blLFP);
    
    nSamps = MAXXCTIME.*LFPsamprate;
    nSampsForXcorr = 2.*nSamps+1;
    xcorrMtx = zeros(nSampsForXcorr, sum(WN.idx.lfp).^2);
    N_validTrials = zeros(nSampsForXcorr, sum(WN.idx.lfp).^2);
    for trl = 1:nTrials
        
        % determine how many time samples the xcorr will have
        tvec = [0:numel(WN.ras{trl, find(WN.idx.lfp>0,1,'first')})] ./ LFPsamprate;
        tvec = tvec + WN.ras{trl, WN.idx.anlgStart};
        [~, stimOnIdx] = min(abs(tvec - WN.trial(trl, WN.idx.stim_on)));
        stimOffTime = (WN.trial(trl, WN.idx.num_frames) ./ WN.sum.exptParams.framerate) + WN.trial(trl, WN.idx.stim_on);
        [~, stimOffIdx] = min(abs(tvec - stimOffTime));
        
        % skip trials where there isn't enough data
        if (stimOffIdx - stimOnIdx) < nSamps
            continue
        elseif ~rem(trl, 10)
            fprintf('Processing trial <%d> of <%d>\n', trl, nTrials);
        end
        
        % preallocate space for the xcorr
        nTimeSamps = stimOffIdx-stimOnIdx+1;
        trlblp_mtx = nan(nTimeSamps, numel(lfpchan));
        
        % make the corrMtx
        for chan = 1:sum(WN.idx.lfp);
            trlblp_mtx(:,chan) = blLFP{trl,chan}(stimOnIdx : stimOffIdx);
        end
        
        % run the xcorr for this trial
        trl_xcorr = xcorr(trlblp_mtx, nSamps, 'coeff');
        xcorrMtx = nansum(cat(3, trl_xcorr, xcorrMtx), 3);
        N_validTrials = sum(cat(3, N_validTrials, ~isnan(trl_xcorr)), 3);
    end
    
    
    % divide by the number of valid trials
    xcorrMtx = xcorrMtx ./ N_validTrials;
    lags = (-nSamps : nSamps) / LFPsamprate;
    dists = distance_between_pins(pins_rel2abs(params.bank, 1:32)); % going to assume the order of AD channels maps directly to pin number
    
    
    % now binning "distances"
    nbins = 11;
    bins = linspace(min(dists)-5*eps, max(dists)+5*eps, nbins);
    xcorrAvgOverDist{ex} = zeros(numel(bins)-1,numel(lags));
    N = [];
    for a = 1:(numel(bins)-1)
        idx = (dists'>=bins(a)) & (dists'<bins(a+1));
        N = [N; idx];
        xcorrAvgOverDist{ex}(a,:) = mean(xcorrMtx(:, idx), 2); % the xcorrs now go across rows!!!
    end
    
end


% plot each channel individually
for ex = 1:numel(cellfun(@iscell, params.expts))
    figure;
    s = surf(lags, bins(1:end-1), xcorrAvgOverDist{ex});
    xlabel('Lag (sec)'); ylabel('Pin distance'); zlabel('Correlation');
    colormap(pmkmp(250, 'CubicL'))
    set(s, 'edgealpha', 0.2, 'linewidth', 0.2)
end

% plot the difference
if numel(cellfun(@iscell, params.expts)) == 2
    
    %
    % plot the difference in the correllograms
    %
    %%%%%%%%%%%%%%%%%%%%%
    figure;
    difference = xcorrAvgOverDist{1} - xcorrAvgOverDist{2};
    s = surf(lags, bins(1:end-1), difference);
    xlabel('Lag (sec)'); ylabel('Pin distance'); zlabel('Difference in correlation');
    colormap(pmkmp(250, 'CubicL'))
    set(s, 'edgealpha', 0.2, 'linewidth', 0.2)
    
    %
    % plot a comparison of the variance of each channel, ordered by the
    % distance of each channel to the center of the stimulus.
    %
    %%%%%%%%%%%%%%%%%%%%
    
    xbar_var = cellfun(@(x) nanmean(x,1), varLFP, 'uniformoutput', 0);
    xbar_var = cat(1,xbar_var{:})';
    SEM_var = cellfun(@(x) nanstd(x,[],1)./sqrt(size(x,1)), varLFP, 'uniformoutput', 0);
    SEM_var = cat(1, SEM_var{:})';
    
    % reorder the data
    [idx, relDist] = reorderPinDistance(params.stimPos, params.bank);
    xbar_var = xbar_var(idx, :);
    SEM_var = SEM_var(idx, :);
    relDist = relDist(idx);
    
    % do some plotting
    figure
    subplot(1,2,1), hold on,
    pltClr = colormap(pmkmp(size(xbar_var,1), 'CubicL'));
    for a = 1:size(xbar_var,1)
        plot(xbar_var(a,1), xbar_var(a,2), '.', 'color', pltClr(a,:));
        plot(repmat(xbar_var(a,1), 2, 1), [xbar_var(a,2)+SEM_var(a,2);xbar_var(a,2)-SEM_var(a,2)], '-', 'color', pltClr(a,:))
        plot([xbar_var(a,1)+SEM_var(a,1);xbar_var(a,1)-SEM_var(a,1)], [xbar_var(a,2);xbar_var(a,2)], '-', 'color', pltClr(a,:))
    end
    xlims = [min(xbar_var(:,1)).*.9, max(xbar_var(:,1)).*1.1];
    ylims = [min(xbar_var(:,2)).*.9, max(xbar_var(:,2)).*1.1];
    minval = min([xlims,ylims]);
    maxval = max([xlims,ylims]);
    plot([minval maxval], [minval maxval], 'k')
    xlabel('Variance of etrodes in expt 1')
    ylabel('Variance of etrodes in expt 2')
    colorbar
    
    subplot(1,2,2)
    hist(xbar_var(:,1)-xbar_var(:,2))
    xlabel('diff in variance')
    ylabel('count')
    [h,p] = ttest(xbar_var(:,1)-xbar_var(:,2));
    title(sprintf('H: %d, p = %g', h, p));
    
    %
    % Variance as a function of distance to stimulus
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure, 
    subplot(1,2,1), hold on,
    plot(relDist, xbar_var(:,1), 'k.');
    plot(relDist, xbar_var(:,2), 'b.');
    xlabel('etrode distance to stimulus')
    ylabel('variance of BL_LFP')
    subplot(1,2,2)
    plot(xbar_var(:,1)-xbar_var(:,2), 'k.')
    
    
end


%% MAP OF THE ARRAY ELECTRODES

array_layout = flipud(reshape(1:100,[10 10])');
map = zeros(size(array_layout));
bank = ['a', 'b', 'c'];
for b = 1:3;
    out_pins = pins_rel2abs(bank(b), 1:32);
    for a = 1:numel(out_pins);
        map(array_layout(:)==out_pins(a)) = b;
    end
end

figure
imagesc(map)
colorbar



%% MAP OF RF LOCATIONS FOR EACH BANK

bank_c = [14 10 15 9  11 10 11 11 14 13 10 10 9 10;...
          10 13 9  11 11 9  10 10 13 10 13 11 7 10]';


figure, hold on,
plot(bank_c(:,1), bank_c(:,2), 'bo', 'markerfacecolor', 'b')
xlim([0 20])
ylim([0 20])




