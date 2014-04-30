function anlyMod_pulseTrains_stimLoc(params)

% An analysis module that gets called by 'invitroAnalysisOverview' and can
% be referenced in the 'physiology_notes.m' script
%
% This module plots the average PSC following a set of LED light pulses
% (either a single pulse or trains). Data from different stimulus locations
% are plotted on the same set of axes. The function will also plot
% normalized PSCs (normalized to the first PSC). 


%
% PLOT THE AVERAGE TRACE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
avgCurrent = {};
tt = {}; % unique tvec for each file in case there are diffs in acquisition length
for a = 1:numel(params.files)
    idx_Im = eval(['params.ax{a}.idx.',params.validCh,'Im']);
    idx_SecVm = regexpi(params.ax{a}.head.recChUnits, 'mV', 'match');
    idx_SecVm = cellfun(@(x) ~isempty(x), idx_SecVm);
    
    % check the Vhold
    if sum(idx_SecVm) == 1;
        actVhold = params.ax{a}.dat(1:100, idx_SecVm, :);
        actVhold = mean(mean(actVhold, 3));
        if (actVhold-params.vHold(a)) > 1; error('Vhold is wrong'); end
    else
        warning('Holding potential could not be verified on file %s', params.files{a})
    end
    
    % pull out the data
    raw = mean(params.ax{a}.dat(:,idx_Im,:),3);
    raw = raw - mean(raw(1:100));
    avgCurrent{a} = raw;
    
    % change the tt (time vector) so that time=0 is the onset of the first
    % pulse
    idx_pulseOn = params.ax{a}.threshold(0.5, params.ax{a}.idx.LEDcmd_470, size(params.ax{a}.dat,3), 'u');
    idx_pulseOn = find(idx_pulseOn, 1, 'first');
    tt{a} = params.ax{a}.tt-params.ax{a}.tt(idx_pulseOn);

end


f = figure; map = colormap('jet'); close(f);
clrIdx = round(linspace(1,size(map,1), numel(params.files)));
figure, hold on,
set(gcf, 'position', [1621 181 1020 589], 'name', sprintf('%s cell %d', params.mouse, params.cellNum))
for a = 1:numel(params.files)
    plot(tt{a},  avgCurrent{a}, 'color', map(clrIdx(a),:), 'linewidth', 2)
end
hold off,
legend(params.legTxt, 'location', 'southeast')
xlabel('Time (sec)')
ylabel('Baseline subtracted current (pA)')
drawnow



%
% PLOT THE NORMALIZED PSC AMPLITUDES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% iterate over the average traces, pull out the data between each pulse,
% and calculate the magnitude of the PSC
for i = 1:numel(params.files)
    
    trigChIdx = params.ax{i}.idx.LEDcmd_470;
    
    % find the pulse offsets. Using the off sets so that I don't include
    % the light onset/offset artifact when looking for the max/min current
    thresh = 0.02;
    sweep = size(params.ax{i}.wf,3);
    [~, cross_time] = params.ax{i}.threshold(thresh, trigChIdx, sweep, 'd');
    
    % find the peaks
    isi = mean(diff(cross_time));
    window = .5e-3; % seconds on either side of the peak
    ptsPerWindow = ceil(window .* params.ax{i}.head.sampRate);
    
    for a = 1:numel(cross_time)
        
        timeStart = cross_time(a) + 1e-3; % add a little to avoid the LED artifact
        timeEnd = timeStart + isi - 2e-3; % subtract a little to avoid the next pulse's LED artifact
        idx = (params.ax{i}.tt >= timeStart) & (params.ax{i}.tt < timeEnd);
        
        % find the max
        vals = avgCurrent{i}(idx);
        [~, minIdx] = min(vals);
        
        % compute the average
        idx = [minIdx-ptsPerWindow : minIdx+ptsPerWindow];
        amp{i}(a) = mean(vals(idx));
        
        % subtract off the decaying current from this PSC from the average
        % trace so that subsequent currents are judged from the correct
        % baseline.
        params.subtractexp = true;
        if params.subtractexp
            fittedCurve = fitexp(vals)
        end

    end
end


% plot a summary of the pulse train data (the peak inward current), both
% normalized and un-normalized. The un-normalized version serves as a
% sanity check for the automated routine to calculate PSC. Only plot when
% there was more than one pulse
if numel(cross_time)>1
    figure, hold on,
    for a = 1:numel(amp)
        plot(1:numel(amp{a}), amp{a}, 'o-', 'color', map(clrIdx(a),:), 'linewidth', 2, 'markerfacecolor', map(clrIdx(a),:))
    end
    hold off
    legend(params.legTxt, 'location', 'southeast')
    set(gcf, 'name', sprintf('%s cell %d: TF = %.3f', params.mouse, params.cellNum, 1/isi), 'position', [2186 296 683 491])
    set(gca, 'xtick', [1:numel(amp{a})], 'xlim', [0.75 numel(amp{a})+.25])
    xlabel('Pulse Number')
    ylabel('Raw PSC (pA)')
    
    figure, hold on,
    for a = 1:numel(amp)
        plot(1:numel(amp{a}), amp{a}./amp{a}(1), 'o-', 'color', map(clrIdx(a),:), 'linewidth', 2, 'markerfacecolor', map(clrIdx(a),:))
    end
    hold off
    legend(params.legTxt, 'location', 'northeast')
    set(gcf, 'name', sprintf('%s cell %d: TF = %.3f', params.mouse, params.cellNum, 1/isi), 'position', [1471 276 683 491])
    set(gca, 'xtick', [1:numel(amp{a})], 'xlim', [0.75 numel(amp{a})+.25])
    xlabel('Pulse Number')
    ylabel('Normalized PSC')
end





