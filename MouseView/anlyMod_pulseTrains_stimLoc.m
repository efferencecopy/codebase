function params = anlyMod_pulseTrains_stimLoc(params)

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
    
    % change the tt (time vector) so that time=0 is the onset of the first
    % pulse
    index = [params.ax{a}.idx.LEDcmd_470, size(params.ax{a}.dat,3)];
    idx_pulseOn = params.ax{a}.threshold(0.1, index, 'u');
    idx_pulseOn = find(idx_pulseOn, 1, 'first');
    tt{a} = params.ax{a}.tt-params.ax{a}.tt(idx_pulseOn);

    
    % pull out the data. baseline subtract everything. define baseline as
    % the time directly prior to the first pulse.
    raw = mean(params.ax{a}.dat(:,idx_Im,:),3);
    bkgndTimeInPts = params.ax{a}.head.sampRate .* 0.100; % 100 msec.
    raw = raw - mean(raw((idx_pulseOn-bkgndTimeInPts):idx_pulseOn));
    avgCurrent{a} = raw;
    
    if ~isempty(params.filter)
       avgCurrent{a} = butterfilt(avgCurrent{a}, params.filter, params.ax{a}.head.sampRate, 'low');
    end
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
fprintf('**** Calculating PSCs for %d files:\n', numel(params.files));
for i = 1:numel(params.files)
    
    fprintf(' Analyzing file %d: %s \n', i,  params.files{i})
    trigChIdx = params.ax{i}.idx.LEDcmd_470;
    
    % find the pulse offsets. Using the off sets so that I don't include
    % the light onset/offset artifact when looking for the max/min current
    thresh = 0.02;
    sweep = size(params.ax{i}.wf,3);
    index = [trigChIdx, sweep];
    [~, cross_time] = params.ax{i}.threshold(thresh, index, 'd');
    
    % find the peaks
    isi = mean(diff(cross_time));
    window = .5e-3; % seconds on either side of the peak
    ptsPerWindow = ceil(window .* params.ax{i}.head.sampRate);
    params.TF(i) = 1/isi;
    
    tmp_avg = avgCurrent{i};
    Npulses = numel(cross_time);
    for a = 1:Npulses
        
        timeStart = cross_time(a) + 1e-3; % add a little to avoid the LED artifact
        timeEnd = timeStart + isi - 2e-3; % subtract a little to avoid the next pulse's LED artifact
        idx_pulse = (params.ax{i}.tt >= timeStart) & (params.ax{i}.tt < timeEnd);
        
        % find the max
        vals = tmp_avg(idx_pulse);
        [~, minIdx] = min(vals);
        
        % check for out of bounds errors
        if minIdx<ptsPerWindow
            minIdx = ptsPerWindow+1;
            warning('No valid start point located')
        end
        
        % compute the average
        idx = [minIdx-ptsPerWindow : minIdx+ptsPerWindow];
        amp{i}(a) = mean(vals(idx));
        
        % subtract off the decaying current from this PSC from the average
        % trace so that subsequent currents are judged from the correct
        % baseline. In order for this to work, the trace should be baseline
        % subtracted from directly before the first pulse. This is because
        % the exponential fits are constrained to asympote at zero.
        if params.subtractexp && Npulses>1
            [fitfun, startidx] = fitexp(vals, 0, params);
            
            % define the new time "zero", which needs to be identical to
            % the one used in fitexp.m
            idx_pulseOn = find(idx_pulse, 1, 'first');
            idx_newZero = idx_pulseOn + startidx;
            
            % make a 'dummy' xx vec, which is the domain over which to evaluate the
            % exponential. It will be shifted with respect to the real
            % domain, and zero padded on the left. Zero pad after
            % evaluating the fit fun.
            xx = 1:(numel(tmp_avg)-idx_newZero+1);
            
            % eval the fitfun and subtract if from the original
            fitted = fitfun(xx);
            fitted = [zeros(1,idx_newZero-1), fitted];
            
            % update the average trace to reflect the subtracted waveform
            tmp_avg = tmp_avg - fitted(:);
            
            if isfield(params, 'debug') && params.debug
                % plot the results for debugging
                if ~exist('dbFig', 'var')
                    dbFig = figure;
                end
                figure(dbFig);
                subplot(Npulses, 2, a+(a-1))
                hold on,
                plot(tmp_avg, 'b')
                plot(fitted, 'r', 'linewidth', 2)
                xmin = find(params.ax{i}.tt == cross_time(1))-100;
                xmax = find(params.ax{i}.tt == cross_time(end))+(isi .* params.ax{i}.head.sampRate)+100;
                xlim([xmin, xmax])
                hold off
                
                subplot(Npulses, 2, a.*2)
                hold on,
                plot(avgCurrent{i}, 'k', 'linewidth', 2)
                plot(tmp_avg, 'r')
                xlim([xmin, xmax])
                hold off
            end
        end

    end
    clear dbFig % so that each iter gets its own fig
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


% as a final step, package the PSC amplitudes into the output argument
params.pscamp = amp;


