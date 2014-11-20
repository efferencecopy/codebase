%% EB_060914_C site 2  CHR2

% TTX alone looks inwards.

fin

% % 5 Hz data
% ax.control = abfobj('2014_06_25_0008');
% ax.synapticBlockers = abfobj('2014_06_25_0010');
% ax.ttx = abfobj('2014_06_25_0012');

% 20 Hz data
ax.control = abfobj('2014_06_25_0007');
ax.synapticBlockers = abfobj('2014_06_25_0009');
ax.ttx = abfobj('2014_06_25_0011');
chIdx = ax.control.idx.HS2_Im;

%% EB_060314_A site 1  CHR2

% AT ROOM TEMP!!!!

fin

% 20 Hz data
ax.control = abfobj('2014_06_17_0000');
ax.synapticBlockers = abfobj('2014_06_17_0004');
ax.ttx = abfobj('2014_06_17_0009');

% % 40 Hz data
% ax.control = abfobj('2014_06_17_0001');
% ax.synapticBlockers = abfobj('2014_06_17_0005');
% ax.ttx = abfobj('2014_06_17_0007');
% xlims = [1.05 1.23];


%% EB_060314_A site 2  CHR2

% the only normal ACSF condition is 40 Hz, so the 20 and 5 Hz stuff (Aside
% from the fiber volley) will be junk.


fin

% % 5 Hz data
ax.control = abfobj('2014_06_17_0010');
ax.synapticBlockers = abfobj('2014_06_17_0013');
ax.ttx = abfobj('2014_06_17_0014');
xlims = [1.05 2.05];
chIdx = ax.control.idx.HS2_Im;


% 20 Hz data
% ax.control = abfobj('2014_06_17_0010');
% ax.synapticBlockers = abfobj('2014_06_17_0012');
% ax.ttx = abfobj('2014_06_17_0015');
% xlims = [1.05 1.34]


% % 40 Hz data
% ax.control = abfobj('2014_06_17_0010');
% ax.synapticBlockers = abfobj('2014_06_17_0011');
% ax.ttx = abfobj('2014_06_17_0016');
% xlims = [1.05 1.23];


%% EB_051914 site 1 CHR2
%
% shows strong depression of synaptic transmission

fin

% 20 Hz FS open
ax.control = abfobj('2014_06_12_0006');
ax.synapticBlockers = abfobj('2014_06_12_0007');
ax.ttx = abfobj('2014_06_12_0009');
chIdx = ax.control.idx.HS2_Im;
% 
% % 20 Hz FS closed
% ax.control = abfobj('2014_06_12_0005');
% ax.synapticBlockers = abfobj('2014_06_12_0008');
% ax.ttx = abfobj('2014_06_12_0010');


%% CH_091114_B Site 1 CHIEF

% 20 Hz FS open
ax.control = abfobj('2014_09_30_0012');
ax.control.head.DACchNames{1} = 'HS1_Vclamp'; %mis named signal
ax.synapticBlockers = abfobj('2014_09_30_0013');
ax.synapticBlockers.head.DACchNames{1} = 'HS1_Vclamp'; %mis named signal
ax.ttx = abfobj('2014_09_30_0015');
ax.ttx.head.DACchNames{1} = 'HS1_Vclamp'; %mis named signal
chIdx = ax.control.idx.HS2_Im;



%% CH_091114_C Site 1 CHIEF

% 20 Hz FS open
ax.control = abfobj('2014_10_02_0001');
ax.control.head.DACchNames{1} = 'HS1_Vclamp'; %mis named signal
ax.synapticBlockers = abfobj('2014_10_02_0002');
ax.synapticBlockers.head.DACchNames{1} = 'HS1_Vclamp'; %mis named signal
ax.ttx = abfobj('2014_10_02_0004');
ax.ttx.head.DACchNames{1} = 'HS1_Vclamp'; %mis named signal
chIdx = ax.control.idx.HS2_Im;

% 40 Hz FS open (the code below won't work for this data set b/c the number
% of samples differ between the control (20hz) and the other 40 Hz expts.)
ax.control = abfobj('2014_10_02_0001');
ax.control.head.DACchNames{1} = 'HS1_Vclamp'; %mis named signal
ax.synapticBlockers = abfobj('2014_10_02_0003');
ax.synapticBlockers.head.DACchNames{1} = 'HS1_Vclamp'; %mis named signal
ax.ttx = abfobj('2014_10_02_0005');
ax.ttx.head.DACchNames{1} = 'HS1_Vclamp'; %mis named signal
chIdx = ax.control.idx.HS2_Im;


%% AK_141103_A Site 2 & 3

fin

% 20 Hz FS open
ax.control = abfobj('2014_11_18_0006');
ax.cadmium = abfobj('2014_11_18_0007');
ax.synapticBlockers = abfobj('2014_11_18_0008');
ax.ttx = abfobj('2014_11_18_0009');
chIdx = ax.control.idx.HS2_Im;


%% Do the analysis

close all

sampFreq = ax.control.head.sampRate;
tt = (0:size(ax.control.dat, 1)-1) ./ sampFreq .* 1000;
filter.Freqs = 400;
filter.Type = 'low';
filter.windowSize = sampFreq .* 0.300;
filter.stepSize = sampFreq .* 0.100;

% figure out when the pulses came on:
thresh = 0.05;
index = [ax.control.idx.LEDcmd_470, 1];
crossings = ax.control.threshold(thresh, index, 'u');
pulseOnset = find(crossings==1, 1, 'first');
tt = tt-tt(pulseOnset);


% mean subtract and filter. store the data in a new structure called "raw"
raw = [];
exptCond = {'control', 'synapticBlockers', 'ttx', 'cadmium'};
for a = 1:numel(exptCond)
    
    if ~isfield(ax, exptCond{a});
        continue
    end
    
    % grab the raw data and baseline subtract
    tmp = ax.(exptCond{a}).dat(:, chIdx,:); % grab the raw data
    tmp = squeeze(tmp);
    tmp = bsxfun(@minus, tmp, mean(tmp(pulseOnset-201:pulseOnset-1, :),1));
    
    % locally detrend
    trendless = locdetrend(tmp, sampFreq);
    
    % filter out the high frequency stuff
    filtered = butterfilt(trendless, filter.Freqs, sampFreq, filter.Type, 1);
    
    % take the mean
    raw.(exptCond{a}) = mean(filtered,2);
    
end





glutamateOnly = raw.control - raw.synapticBlockers;
fiberVolley = raw.synapticBlockers - raw.ttx;



figure
set(gcf, 'position', [28 296 1375 506])
plot(tt, raw.control, 'k', 'linewidth', 3)
%xlim(xlims)
xlabel('time (sec)')
ylabel('LFP amplitude')
title('Control Condition')
axis tight


figure
set(gcf, 'position', [28 296 1375 506])
plot(tt, raw.synapticBlockers, 'k', 'linewidth', 3)
xlabel('time (sec)')
ylabel('LFP amplitude')
title('With Synaptic Blockers')
axis tight


figure
set(gcf, 'position', [28 296 1375 506])
plot(tt, raw.ttx, 'k', 'linewidth', 3)
xlabel('time (sec)')
ylabel('LFP amplitude')
title('With TTX')
axis tight


figure
set(gcf, 'position', [28 296 1375 506])
plot(tt, glutamateOnly, 'k', 'linewidth', 3)
xlabel('time (sec)')
ylabel('LFP amplitude')
title('LFP due to synaptic transmission')
axis tight


figure, hold on,
set(gcf, 'position', [28 296 1375 506])
plot(tt, fiberVolley, 'k', 'linewidth', 3)
plot(tt(crossings), zeros(1, sum(crossings)), 'or', 'markerfacecolor', 'r')
xlabel('time (sec)')
ylabel('LFP amplitude')
title('LFP due to fiber volley')
axis tight




