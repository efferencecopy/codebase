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
for a = 1:size(params.files ,1)
    idx_Im = eval(['params.ax{a}.idx.',params.validCh,'Im']);
    idx_Vclamp = eval(['params.ax{a}.idx.',params.validCh,'Vclamp']);
    idx_SecVm = eval(['params.ax{a}.idx.',params.validCh,'secVm']);
    
    % check the Vhold
    actVhold = params.ax{a}.dat(1:100, idx_SecVm, :);
    actVhold = mean(mean(actVhold, 3));
    if (actVhold-params.vHold(a)) > 1; error('Vhold is wrong'); end
    
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


map = colormap('jet'); close;
clrIdx = round(linspace(1,size(map,1), size(params.files,1)));
figure, hold on,
set(gcf, 'position', [212 51 1067 758])
for a = 1:size(params.files,1)
    plot(tt{a},  avgCurrent{a}, 'color', map(clrIdx(a),:), 'linewidth', 2)
end
leg = params.legTxt;
legend(leg')
xlabel('Time (sec)')
ylabel('Baseline subtracted current (pA)')



%
% PLOT THE NORMALIZED PSC AMPLITUDES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% iterate over the average traces, pull out the data between each pulse,
% and calculate the magnitude of the PSC
for i = a:numel(params.files)
    
    trigChIdx = params.ax{i}.idx.LEDcmd_470;
    dataChIdx = eval(['params.ax{a}.idx.',params.validCh,'Im']);
    
    % find the pulse onsets
    thresh = 0.02;
    sweep = size(params.ax{i}.wf,3);
    [~, cross_time] = params.ax{i}.threshold(thresh, trigChIdx, sweep, 'u');
    
    % find the peaks
    isi = mean(diff(cross_time)) - 0.010;
    window = .5e-3; % seconds on either side of the peak
    ptsPerWindow = ceil(window .* params.ax{i}.head.sampRate);
    for a = 1:numel(cross_time)
        
        timeStart = cross_time(a);
        timeEnd = timeStart + isi;
        idx = (params.ax{i}.tt >= timeStart) & (params.ax{i}.tt < timeEnd);
        
        % find the max
        vals = avgCurrent{i}(idx);
        [~, minIdx] = min(vals);
        
        % compute the average
        idx = [minIdx-ptsPerWindow : minIdx+ptsPerWindow];
        amp = mean(vals(idx));
        
        % deal with the background
        timeEnd = cross_time(a);
        timeStart = timeEnd - window;
        idx = (ax.tt >= timeStart) & (ax.tt < timeEnd);
        bkgnd = mean(avg(idx, dataChIdx));
        psc{fl}(a) = amp - bkgnd;
    end
end




