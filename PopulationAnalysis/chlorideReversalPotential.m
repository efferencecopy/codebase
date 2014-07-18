%% CHLORIDE REVERSAL POTENTIAL FOR K-GLUCONATE 'BUSTED' INTERNAL MADE ON 11/12/13

fin

mdb = initMouseDB();

in = {'EB_012414_A', '2014_02_13_0002', -60, 1.3306, 1.3310;...
      'EB_012414_A', '2014_02_13_0014', -60, 1.3296, 1.3302;...
      'EB_012414_A', '2014_02_13_0029', -60, 1.3320, 1.3331;...
      'EB_012414_A', '2014_02_13_0038', -60, 1.3294, 1.3304;...
      'EB_012414_A', '2014_02_13_0039', -70, 1.3302, 1.3319;...
      'EB_012214_C', '2014_02_11_0009', -60, 1.3333, 1.3359;...
      'EB_012214_A', '2014_02_09_0002', -60, 1.3313, 1.3333};
  
  
% roll through the data files, plot the raw data, and extract the
% parameters that are important

for a = 1:size(in,1)
    
    
    fpath = findfile(in{a,2}, [GL_DATPATH, in{a,1}], '.abf');
    [dat, si, h] = abfload(fpath);
    
    N = size(dat, 1);
    tt = (0:N-1) .* (si * 1e-6);
    
    % plot the data
    window = [in{a,4}-0.030, in{a,5}+0.050];
    idx = (tt >= window(1)) & (tt < window(2));
    figure
    plot(tt(idx), permute(dat(idx,1,:), [1,3,2]))
    axis tight
    xlabel('time')
    ylabel('Current')
    title([in{a,1}, ' ', in{a,2}])
    
    % extract the average current in the analysis window
     window = [in{a,4}, in{a,5}]; % seconds from sweep start
     idx = (tt >= window(1)) & (tt < window(2));
     tmp = permute(dat(idx,1,:), [1,3,2]);
     Im_peak = mean(tmp, 1);
     
     % deal with the baseline response
     window = [1.32, 1.3275]; % seconds from sweep start
     idx = (tt >= window(1)) & (tt < window(2));
     tmp = permute(dat(idx,1,:), [1,3,2]);
     Im(:,a) =Im_peak - mean(tmp, 1);
    
     % customize the voltage relative to the holding potential used in the
     % experiments
     Vcmd(:,a) = linspace(-40,20,13) + in{a,3};
    
end


figure, hold on,
plot(Vcmd, Im, 'linewidth', 3)
plot([-110 -70], [0 0], 'k', 'linewidth', 2)
legend({in{:,2}}, 'location', 'northwest')
xlabel('Voltage (mv)')
ylabel('Current (pA)')




%% CHLORIDE REVERSAL POTENTIAL FOR Cs-GLUCONATE INTERNAL MADE ON 07/03/14

fin

% THE DATA: < MOUSE NAME, DATA PREFIX, DATA FILE NUMBER, ANALYSIS WINDOW, PAIR NUMBER >
in = {'CH_061614_D', '2014_07_09_', 3,       (37:90), 1;... % 18 minutes, holding at -60
      'CH_071514_A', '2014_07_15_', (7:37),  (30:50), 2;... % minimum of 4 min, max of 12 min, Vhold = -60, only negative pulses
      'CH_071514_A', '2014_07_15_', (53:57), (32:54), 3;... % Vhold = -60, only negative pulses
      'CH_071614_A', '2014_07_16_', (1:8),   (33:49), 4;... % 1 to 2 minutes of dialysis, Vhold = -40
      'CH_071614_A', '2014_07_16_', (9:28),  (31:53), 4;... % 7 to 13 minutes, Vhold = -60
      'CH_071614_A', '2014_07_16_', (37:52), (41:80), 5;... % zero to 8 minutes, Vhold = -30
      'CH_071614_A', '2014_07_16_', (56:74), (38:54), 5};   % 15 to 32 minutes, Vhold = 0

for ex = 1:size(in,1);
    fprintf('now working on experiment %d\n', ex)
    
    % make a list of ABF files to concatenate
    fNames = {};
    for a = 1:numel(in{ex,3})
        suffix = num2str(in{ex,3}(a));
        nZerosNeeded = 4-numel(suffix);
        suffix = [repmat('0',1,nZerosNeeded), suffix];
        fNames{a} = [in{ex,2}, suffix];
    end
    
    % concatenate the files
    catdim = 4;
    ax = abfcat(catdim, fNames);
    
    % figure out which channel has the threshold crossing
    rec_mV = strcmpi(ax.head.recChUnits, 'mv');
    sec_ch = cellfun(@(x) ~isempty(x), regexpi(ax.head.recChNames, 'sec'));
    Iclamp_ch = rec_mV & ~sec_ch;
    assert(sum(Iclamp_ch)==1, 'Error: too many Iclamp channels')
    Iclamp_HS = ax.head.recChNames{Iclamp_ch}(1:3);
    thresh_ch = cellfun(@(x) ~isempty(x), regexpi(ax.head.DACchNames, Iclamp_HS));
    thresh_ch = find(thresh_ch);
    
    % figure out which channel has the voltage clamp measurement
    rec_pA = strcmpi(ax.head.recChUnits, 'pa');
    vclamp_ch = rec_pA & ~sec_ch;
    vclamp_monitor = rec_mV & sec_ch;
    
    % designate the amount of time pre and post pulse to consider
    baselinePts = round(ax.head.sampRate .* .050);
    analysisPts = round(ax.head.sampRate .* .050);
    
    %ax.dat is (time x channels x sweeps x files). Iterate over files
    %within each sweep and create a distribution of IPCSs for each sweep.
    %Because each sweep is a different holding potential, this will build
    %up a distrbution of currents for each Vhold.
    avgCurrentPerVHold = nan(size(ax.dat,3), baselinePts+analysisPts+1);
    avgHoldingPotential = nan(1, size(ax.dat,3));
    for swp = 1:size(ax.dat, 3);
        
        rawData = nan(size(ax.dat,4), baselinePts+analysisPts+1);
        rawVhold = nan(size(ax.dat,4), 1);
        for fl = 1:size(ax.dat,4)
            
            % find the time of the presynaptic spike.
            idx = threshold(ax, 5, [thresh_ch, swp, fl], 'u');
            idx = find(idx);
            assert(numel(idx)==1, 'ERROR: too many threshold crossings');
            
            % grab the data (pA)
            raw = ax.dat([idx-baselinePts:idx+analysisPts], vclamp_ch, swp, fl);
            baseline = mean(raw(1:baselinePts));
            raw = raw-baseline;
            raw = butterfilt(raw, 500, ax.head.sampRate, 'low');
            rawData(fl,:) = raw;
            
            % grab the actual membrane potential
            rawVhold(fl) = mean(ax.dat([idx-baselinePts:idx+analysisPts], vclamp_monitor, swp, fl));
        end
        
        avgCurrentPerVHold(swp, :) = mean(rawData, 1);
        avgHoldingPotential(swp) = mean(rawVhold);
    end
    
    % determine the magnitude of the current
    [PSC_avg, PSC_max] = deal([]);
    window = in{ex, 4};
    
    [~, maxIdx] = max(abs(avgCurrentPerVHold(:,[baselinePts+window])), [],  2);
    maxIdx = maxIdx+baselinePts+window(1);
    for i = 1:numel(maxIdx);
        PSC_max(i) = avgCurrentPerVHold(i, maxIdx(i));
        PSC_avg(i) = mean(avgCurrentPerVHold(i,[maxIdx(i)-10 : maxIdx+10]), 2);
    end
    
    
    
    % plot the data
    figure
    subplot(1,2,1), hold on,
    tt = ([idx-baselinePts:idx+analysisPts] - idx) ./ ax.head.sampRate .* 1000;
    plot(tt, avgCurrentPerVHold')
    plot(tt([baselinePts+window]), avgCurrentPerVHold(:,[baselinePts+window])', '.');
    xlabel('Time (ms)');
    ylabel('current (pA)')
    hold off
    subplot(1,2,2), hold on,
    plot(avgHoldingPotential, PSC_avg, '-ko')
    plot(avgHoldingPotential, PSC_max, '-bo')
    xlabel('Voltage (mV')
    ylabel('Current (pA)')
    hold off
    drawnow
    
    % combine the population data
    popdat{ex}.psc_avg = PSC_avg;
    popdat{ex}.psc_max = PSC_max;
    popdat{ex}.vhold = avgHoldingPotential;
    
end


% present all the data on the same axis;
pairs = [in{:, 5}]';
clrmap = colormap('jet'); close;
inds = round(linspace(1, size(clrmap,1), numel(unique(pairs))));
clrmap = clrmap(inds, :);

figure, hold on,
for a = 1:numel(popdat)
    plot(popdat{a}.vhold, popdat{a}.psc_avg, '-o',...
         'color', clrmap(pairs(a), :),...
         'markeredgecolor', clrmap(pairs(a), :),...
         'linewidth', 3)
end
plot([-100 -40], [0 0], 'k', 'linewidth', 2)
xlabel('Holding potential (mV)')
ylabel('IPSC magnitude (pA)')






  
  