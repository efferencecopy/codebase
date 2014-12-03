%% EB_060914_C site 2  CHR2

% TTX alone looks inwards.

fin

% 5 Hz data
ax.control = abfobj('2014_06_25_0008');
ax.synapticBlockers = abfobj('2014_06_25_0010');
ax.ttx = abfobj('2014_06_25_0012');

% % 20 Hz data
% ax.control = abfobj('2014_06_25_0007');
% ax.synapticBlockers = abfobj('2014_06_25_0009');
% ax.ttx = abfobj('2014_06_25_0011');
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




%% AK_141103_A Site 1

fin

% 5 Hz FS open
ax.control = abfobj('2014_11_18_0001');

ax.cadmium = abfobj('2014_11_18_0002');
ax.cadmium = ax.cadmium.removeSweeps([1:20]); % CdCl comes in on sweep 6, stable at sweep 20

ax.synapticBlockers = abfobj('2014_11_18_0003'); % includes CdCl

ax.ttx = abfobj('2014_11_18_0004');

chIdx = ax.control.idx.HS2_Im;



%% AK_141103_A Site 2 & 3

% notes:
%
% 20 Hz FS open
%
% be sure to look at both channels

fin


ax.control = abfcat(3, {'2014_11_18_0005', '2014_11_18_0006'});

ax.cadmium = abfobj('2014_11_18_0007');
ax.cadmium = ax.cadmium.removeSweeps([1:20]); % CdCl comes in on sweep 6, stable at sweep 20

ax.synapticBlockers = abfobj('2014_11_18_0008'); % includes CdCl

ax.ttx = abfobj('2014_11_18_0009');

chIdx = ax.control.idx.HS2_Im; 



%%

%
%
%

fin

% 20 Hz
ax.control = abfobj('2014_11_21_0000');
ax.synapticBlockers = abfobj('2014_11_21_0003'); % includes CdCl
ax.ttx = abfobj('2014_11_21_0005'); % includes CdCl
% 
% % 40 Hz
% ax.control = abfobj('2014_11_21_0001');
% ax.synapticBlockers = abfobj('2014_11_21_0004'); % includes CdCl
% ax.ttx = abfobj('2014_11_21_0006'); % includes CdCl

% look at both channels
chIdx = ax.control.idx.HS2_Im; 

%% Do the analysis

close all

sampFreq = ax.control.head.sampRate;
tt = (0:size(ax.control.dat, 1)-1) ./ sampFreq .* 1000;

filter.lp_freqs = 2000; % lower than 2e3 is bad b/c it cuts our important events...


% figure out when the pulses came on:
thresh = 0.05;
index = [ax.control.idx.LEDcmd_470, 1];
crossings = ax.control.threshold(thresh, index, 'u');
pulseOnset = find(crossings==1, 1, 'first');

% figure out if there was a test pulse (which is an artifact from using a
% Vclamp protocol)
thresh = -0.001; % in mV
index = [ax.control.idx.HS2_Vclamp, 1];
crossings = ax.control.threshold(thresh, index, 'u');
Rs_testPulse_off = find(crossings==1, 1, 'first');

% mean subtract and filter. store the data in a new structure called "trace"
trace = [];
exptCond = {'control', 'synapticBlockers', 'ttx', 'cadmium'};
for a = 1:numel(exptCond)
    
    if ~isfield(ax, exptCond{a});
        continueRs_testPulse_off
    end
    
    
    
    % grab the raw data
    tmp = ax.(exptCond{a}).dat(:, chIdx,:); % grab the raw data
    tmp = squeeze(tmp);
    
    
    % baseline subtract
    tmp = bsxfun(@minus, tmp, mean(tmp(pulseOnset-201:pulseOnset-1, :),1));
    
    
    % filter out the high frequency stuff
    filtered = butterfilt(tmp, filter.lp_freqs, sampFreq, 'low', 1);
    
    % take the mean
    average = mean(filtered,2);
    
    % try to reduce 60 cycle noise
    lines = 60;
    winStart = Rs_testPulse_off + (sampFreq * 0.150);
    winEnd = pulseOnset;
    average= rmhum(average, sampFreq, winStart, winEnd,  lines, 1);
    
    
    trace.(exptCond{a}) = average;
    
end





trace.glutamateOnly = trace.control - trace.synapticBlockers;
trace.fiberVolley = trace.synapticBlockers - trace.ttx;

%
% plot the results
%
exptCond = {exptCond{:}, 'glutamateOnly', 'fiberVolley'};
for a = 1:numel(exptCond)
    
    if ~isfield(trace, exptCond{a})
        continue
    end
    
    figure, hold on,
    set(gcf, 'position', [28 296 1375 506], 'name', exptCond{a})
    tt = tt-tt(pulseOnset);
    plot(tt, trace.(exptCond{a}), 'k', 'linewidth', 3)
    plot(tt(crossings), zeros(1,sum(crossings)), 'ro', 'markerfacecolor', 'r')
    xlabel('time (ms)')
    ylabel('LFP amplitude')
    title(exptCond{a})
    axis tight
    
end




