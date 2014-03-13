%% TO DO 

% automate the plot windows based off the thresh crossings
%
% make a method for psp, specify 'i', or 'e' for polarity
%
% Do the exponential fit
%
% Scale the epsp and fit the inhibition.
%
% make abfcat more descriptive and diagnostic. Failures are a puzzle...

%% CH_020314_B CELL #1

% might be in V1?!?!?
% some facilitation at 20 hz

fin

% load in the data.
% <abf name(s)>, <elim sweeps>, <vhold>, <pulse freq>, <led voltage>
files = {
    {'2014_02_20_0015', '2014_02_20_0016'}, {}, -85, 5;
    {'2014_02_20_0017', '2014_02_20_0018'}, {}, -85, 20;
    {'2014_02_20_0019'}, {}, -85, 40;
    };


%% CH_020314_B CELL #2

% definitely not in V1
% band pass for the P1/P2 comparison, otherwise, depressing


fin

% load in the data.
% <abf name(s)>, <elim sweeps>, <vhold>, <pulse freq>, <led voltage>
files = {
    {'2014_02_20_0025'}, {}, -85, 5;
    {'2014_02_20_0026'}, {}, -85, 40;
    {'2014_02_20_0027'}, {}, -85, 60;
    {'2014_02_20_0032'}, {}, -85, 2;
    {'2014_02_20_0033'}, {}, -85, 20;
    };


%% CH_020314_B CELL #3

% definitely not in V1.
% mostly depressing, but different lamina than cell #2, even though it was
% in the same cortical area.

fin

% load in the data.
% <abf name(s)>, <elim sweeps>, <vhold>, <pulse freq>, <led voltage>
files = {
    {'2014_02_20_0052'}, {}, -85, 40;
    {'2014_02_20_0057'}, {}, -85, 5;
    {'2014_02_20_0049', '2014_02_20_0059'}, {}, -85, 20;
    };


% alternate file for 20 Hz, with slight larger amp LED
% {'2014_02_20_0046'}, {}, -85, 20;

% or just look at _0059 instead of _0049 and _0059


%% CH_020314_A CELL #1

fin

% possibly V1

% load in the data.
% <abf name(s)>, <elim sweeps>, <vhold>, <pulse freq>, <led voltage>
files = {
    {'2014_02_19_0004'}, {}, -85, 5;
    {'2014_02_19_0005'}, {}, -85, 20;
    };


%% CH_020314_A CELL #2

fin

% almost certainly V1

% load in the data.
% <abf name(s)>, <elim sweeps>, <vhold>, <pulse freq>, <led voltage>
files = {
    {'2014_02_19_0016'}, {}, -85, 5;
    {'2014_02_19_0017'}, {}, -85, 40;
    };



%% CH_020314_C CELL #3

fin

% unknown area...
% strongly depressing


% load in the data.
% <abf name(s)>, <elim sweeps>, <vhold>, <pulse freq>, <led voltage>
files = {
    {'2014_02_18_0029'}, {}, -85, 10;
    {'2014_02_18_0030'}, {}, -85, 20;
    {'2014_02_18_0032'}, {}, -85, 20;
    {'2014_02_18_0033'}, {}, -85, 5;
    {'2014_02_18_0038'}, {}, -85, 20;
    {'2014_02_18_0039'}, {}, -85, 10;
    {'2014_02_18_0040'}, {}, -85, 40;
    {'2014_02_18_0041'}, {}, -85, 5;
    };


%% CH_020314_C CELL #4

fin

% unknown area
% depressing

% load in the data.
% <abf name(s)>, <elim sweeps>, <vhold>, <pulse freq>, <led voltage>
files = {
    {'2014_02_18_0044'}, {}, -85, 5;
    {'2014_02_18_0045'}, {}, -85, 20;
    {'2014_02_18_0046'}, {}, -85, 10;
    {'2014_02_18_0047'}, {}, -85, 40;
    };



%% CH_020314_D CELL #1.... FS cell in V1?

fin

% Possibly in V1
% Facilitation at 20 and 40 Hz

% load in the data.
% <abf name(s)>, <elim sweeps>, <vhold>, <pulse freq>, <led voltage>
files = {
    {'2014_02_22_0005'}, {}, -85, 5;
    {'2014_02_22_0006'}, {}, -85, 20;
    {'2014_02_22_0010'}, {}, -80, 40; % could also be _0010 or _0009
    };



%% CH_020314_D CELL #2

fin

% definintely in a HOA
% P1/P2 close to one, but depressing after that


% load in the data.
% <abf name(s)>, <elim sweeps>, <vhold>, <pulse freq>, <led voltage>
files_85 = {    
    {'2014_02_22_0022'}, {}, -85, 20, 1.05;    
    {'2014_02_22_0024'}, {}, -85, 5, 1.05;
    {'2014_02_22_0027'}, {}, -85, 40, 1.05;
    };

files_45 = {    
    {'2014_02_22_0023'}, {}, -45, 20, 1.05;    
    {'2014_02_22_0025'}, {}, -45, 5, 1.05;
    {'2014_02_22_0026'}, {}, -45, 40, 1.05;
    };

files = files_85;


%% CH_020314_D CELL #3

fin

% load in the data.
% <abf name(s)>, <elim sweeps>, <vhold>, <pulse freq>, <led voltage>
files_85 = {
    {'2014_02_22_0031'}, {}, -85, 5, 2.2;
    {'2014_02_22_0033'}, {}, -85, 20, 2.2;    
    {'2014_02_22_0036'}, {}, -85, 40, 2.2;    
    {'2014_02_22_0037'}, {}, -85, 60, 2.2;
    };


files_45 = {
    {'2014_02_22_0032'}, {}, -45, 5, 2.2;
    {'2014_02_22_0034'}, {}, -45, 20, 2.2;    
    {'2014_02_22_0035'}, {}, -45, 40, 2.2;    
    {'2014_02_22_0038'}, {}, -45, 60, 2.2;
    };

files = files_85;

%% %% CRUNCH THE NUMBERS (E/I ratio)

% OUTLINE
%
% 1) import the raw data, from -85 and -45. 
%
% 2) Identify the pulse times, plot the -45 and -85 on top of one another
% to make sure that the timing is corrent.
%
% 3) Do some error checking to make sure that pulse freq, Vhold, and pulse
% times are the same between the two files
%
% 4) By subtraction, pull out the purely inhibitory CURRENT from the -45 mv
% file
%
% 5) Find the (Baseline corrected) peak current (ex/inhib). Back out the
% baseline either by (1) fitting an exponential, or (2) subtracting the
% previous decaying trace off (or the residuals after subtraction for
% pulses 3:end)
%
% 6) Turn the "currents" into conductances, and then normalize all the
% values to the first pulse. Plot normalized Ge and Gi as a function of
% pulse number.
%
% 7) Sumarize the figure (1st/2nd, or 1st/last) and plot that value as a
% function of FREQUENCY across all the frequencies tested



% lump the different conditions together and sort them out later
files = cat(1, files_45, files_85); 

% loop through all the files
dat = [];
for a = 1:size(files,1)
    
    % import the data
    ax = abfcat(files{a,1})
    
    % id the data and triger channels
    dataChIdx = ax.idx.HS2_Im;
    trigChIdx = ax.idx.LEDcmd_470;
    
    % compute the average across sweeps
    avg = mean(ax.dat, 3);
    
    % find the pulse times and pull out the relevant epoch. Save the pulse
    % period, a time vector (adjust the tvec so that t=0 is the first
    % pulse)
    thresh = 0.02;
    sweep = 1;
    [cross_idx, cross_time] = ax.threshold(thresh, trigChIdx, sweep, 'u');
    cross_idx = find(cross_idx);
    ipi = unique(diff(cross_idx));
    epoch_idx = cross_idx(1)-ipi(1) : cross_idx(end)+ipi(1);
    epoch.Im = avg(epoch_idx);
    epoch.tvec = ax.tt(epoch_idx) - cross_time(1); % make the first pulse t=0;
    epoch.tcross = cross_time - cross_time(1);
    
    %
    % do some error checking:
    %
    
    % Vclamp
    expectedVhold = files{a, 3};
    actualVhold = ax.dat(epoch_idx, ax.idx.HS2_secVm, :);
    correctHold = all((actualVhold(:) > expectedVhold - 1) & (actualVhold(:) < expectedVhold + 1)); % +/- 1 mV during epoch
    assert(correctHold, 'ERROR: Vclamp tolerance exceeded');    
    
    % pulse frequency
    expectedPulseFreq = files{a, 4};
    actualPulseFreq = 1./(cross_time(2) - cross_time(1));
    assert(abs((expectedPulseFreq-actualPulseFreq))<1, 'ERROR: Pulse frequency incorrectly specified');
    
    % store the data
    dat{a,1} = expectedPulseFreq;
    dat{a,2} = expectedVhold;
    dat{a,3} = epoch;
end



% Now do the subtraction...
TFs = [dat{:,1}]';
Vhold = [dat{:,2}]';
uniqueTFs = unique(TFs);
for a = 1:numel(uniqueTFs)
   % find the corresponding -85 and -45 files 
    l_TF = TFs == uniqueTFs(a);
    l_85 = Vhold == -85;
    l_45 = Vhold == -45;
    
    % get ready for the subtraction
    epsc = dat{l_TF&l_85, 3}; % pure espc
    mpsc = dat{l_TF&l_45, 3}; % mixed current psc
    
    % do the subtraction
    ipsc = subtractEPSC(epsc, mpsc);
end




%% CRUNCH THE NUMBERS
dat = [];
for fl = 1:size(files,1)
    
    % condition the inputs
    ax = abfcat(files{fl,1});
    
    % average the data across sweeps.
    dataChIdx = ax.idx.HS2_Im;
    trigChIdx = ax.idx.LEDcmd_470;
    
    % compute the average across sweeps
    avg = mean(ax.dat, 3);
    
    switch files{fl, 4}
        case 2
            t_idx = ax.tt>1.05 & ax.tt<3.2;
        case 5
            t_idx = ax.tt>1.05 & ax.tt<2;
        case 10
            t_idx = ax.tt>1.05 & ax.tt<2;
        case 20
            t_idx = ax.tt>1.05 & ax.tt<1.39;
        case 40
            t_idx = ax.tt>1.05 & ax.tt<1.25;
        case 60
            t_idx = ax.tt>1.06 & ax.tt<1.2;
    end
    
    
    figure
    s1 = subplot(3,1,1);
    plot(ax.tt(t_idx), avg(t_idx, dataChIdx))
    axis tight
    ylabel('current (pa)')
    title(ax.name)
    s2 = subplot(3,1,2);
    plot(ax.tt(t_idx), ax.wf(t_idx, trigChIdx, 1))
    axis tight
    ylabel('voltage')
    xlabel('time (ms)')
    
    
    % find the pulse onsets
    thresh = 0.02;
    sweep = 1;
    [cross_idx, cross_time] = ax.threshold(thresh, trigChIdx, sweep, 'u');
    
    % make sure the threshold crossing are correct
    subplot(3,1,1), hold on,
    plot(cross_time, avg(cross_idx, dataChIdx), 'mo', 'markerfacecolor', 'm')
    
    % find the peaks
    isi = mean(diff(cross_time)) - 0.010;
    window = .5e-3; % seconds on either side of the peak
    ptsPerWindow = ceil(window .* ax.head.sampRate);
    for a = 1:numel(cross_time)
        
        timeStart = cross_time(a);
        timeEnd = timeStart + isi;
        idx = (ax.tt >= timeStart) & (ax.tt < timeEnd);
        
        % find the max
        vals = avg(idx, dataChIdx);
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
    
    
    subplot(3,1,3)
    plot(psc{fl}./psc{fl}(1), 'o-')
    set(gcf, 'position', [53    20   560   791], 'name', [num2str(files{fl, 4}), ' Hz'])
    
    % fot comparision with all the other TFs
    TFs(fl) = files{fl,4};
    TF_leg{fl} = num2str(TFs(fl));
end % files


% all the data together
colors = {'k', 'r', 'b', 'g', 'c', 'm'};
figure, hold on,
for a = 1:numel(TFs)
    plot(psc{a}./psc{a}(1), 'color', colors{a}, 'linewidth', 3)
end
set(gca, 'yscale', 'log')
set(gca, 'fontSize', 26)
xlabel('Pulse number')
ylabel('PN / P1')
legend(TF_leg)



% p1/p2 acrss TF
[~, idx] = sort(TFs);
p1p2 = cellfun(@(x) x(2)/x(1), psc);
p1p2 = p1p2(idx);
p1pend = cellfun(@(x) x(end)/x(1), psc)
p1pend = p1pend(idx);

figure, hold on,
plot(TFs(idx), p1p2, 'k', 'linewidth', 3)
plot(TFs(idx), p1pend, 'b', 'linewidth', 3)
set(gca, 'xscale', 'log', 'yscale', 'log')
ylim([min([p1p2, p1pend]).*.9, max([p1p2, p1pend]).*1.1])
xlim([min(TFs).*.9, max(TFs).*1.1])
set(gca, 'xtick', [5, 10, 20, 40])
set(gca, 'fontSize', 26)


    

