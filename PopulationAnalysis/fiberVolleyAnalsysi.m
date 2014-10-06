%% EB_060914_C site 2 CHR2

% TTX alone looks inwards.

fin

% % 5 Hz data
% ax_control = abfobj('2014_06_25_0008');
% ax_synapticBlockers = abfobj('2014_06_25_0010');
% ax_ttx = abfobj('2014_06_25_0012');

% 20 Hz data
ax_control = abfobj('2014_06_25_0007');
ax_synapticBlockers = abfobj('2014_06_25_0009');
ax_ttx = abfobj('2014_06_25_0011');


%% EB_060314_A site 1 CHR2

% AT ROOM TEMP!!!!

fin

% 20 Hz data
ax_control = abfobj('2014_06_17_0000');
ax_synapticBlockers = abfobj('2014_06_17_0004');
ax_ttx = abfobj('2014_06_17_0009');

% % 40 Hz data
% ax_control = abfobj('2014_06_17_0001');
% ax_synapticBlockers = abfobj('2014_06_17_0005');
% ax_ttx = abfobj('2014_06_17_0007');
% xlims = [1.05 1.23];


%% EB_060314_A site 2 CHR2

% the only normal ACSF condition is 40 Hz, so the 20 and 5 Hz stuff (Aside
% from the fiber volley) will be junk.


fin

% % 5 Hz data
% ax_control = abfobj('2014_06_17_0010');
% ax_synapticBlockers = abfobj('2014_06_17_0013');
% ax_ttx = abfobj('2014_06_17_0014');
% xlims = [1.05 2.05];


% 20 Hz data
ax_control = abfobj('2014_06_17_0010');
ax_synapticBlockers = abfobj('2014_06_17_0012');
ax_ttx = abfobj('2014_06_17_0015');
xlims = [1.05 1.34]


% % 40 Hz data
% ax_control = abfobj('2014_06_17_0010');
% ax_synapticBlockers = abfobj('2014_06_17_0011');
% ax_ttx = abfobj('2014_06_17_0016');
% xlims = [1.05 1.23];


%% EB_051914 site 1 CHR2
%
% shows strong depression of synaptic transmission

fin

% 20 Hz FS open
ax_control = abfobj('2014_06_12_0006');
ax_synapticBlockers = abfobj('2014_06_12_0007');
ax_ttx = abfobj('2014_06_12_0009');

% 
% % 20 Hz FS closed
% ax_control = abfobj('2014_06_12_0005');
% ax_synapticBlockers = abfobj('2014_06_12_0008');
% ax_ttx = abfobj('2014_06_12_0010');



%% CH_091114_C site 1 ChIEF

% 20 Hz FS open
ax_control = abfobj('2014_06_12_0006');
ax_synapticBlockers = abfobj('2014_06_12_0007');
ax_ttx = abfobj('2014_06_12_0009');

% 40 Hz FS open
ax_control = abfobj('2014_06_12_0006');
ax_synapticBlockers = abfobj('2014_06_12_0007');
ax_ttx = abfobj('2014_06_12_0009');



%% Do the analysis

close all

sampFreq = ax_control.head.sampRate;
tt = (0:size(ax_control.dat, 1)-1) ./ sampFreq;
filterFreqs = 500;
filterType = 'low';

% notch filter
CONDITION = false;
d = designfilt('bandstopiir','FilterOrder',6, ...
    'HalfPowerFrequency1',5,'HalfPowerFrequency2',6, ...
    'DesignMethod','butter','SampleRate',sampFreq);

% control data
tmp = permute(ax_control.dat(:,1,:), [3,1,2]);
baseline = mean(tmp(:, 3000:4000), 2);
tmp = bsxfun(@minus, tmp, baseline);
tmp = butterfilt(tmp, filterFreqs, sampFreq, filterType, 2);
if CONDITION
    for a = 1:size(tmp, 1)
        tmp(a,:) = filtfilt(d, tmp(a,:));
    end
end
trace_control = mean(tmp, 1);

% synaptic blockers data
tmp = permute(ax_synapticBlockers.dat(:,1,:), [3,1,2]);
baseline = mean(tmp(:, 3000:4000), 2);
tmp = bsxfun(@minus, tmp, baseline);
tmp = butterfilt(tmp, filterFreqs, sampFreq, filterType, 2);
if CONDITION
    for a = 1:size(tmp, 1)
        tmp(a,:) = filtfilt(d, tmp(a,:));
    end
end
trace_synapticBlockers = mean(tmp, 1);

% ttx data
tmp = permute(ax_ttx.dat(:,1,:), [3,1,2]);
baseline = mean(tmp(:, 3000:4000), 2);
tmp = bsxfun(@minus, tmp, baseline);
tmp = butterfilt(tmp, filterFreqs, sampFreq, filterType, 2);
if CONDITION
    for a = 1:size(tmp, 1)
        tmp(a,:) = filtfilt(d, tmp(a,:));
    end
end
trace_ttx = mean(tmp, 1);


glutamateOnly = trace_control - trace_synapticBlockers;
fiberVolley = trace_synapticBlockers - trace_ttx;



figure
set(gcf, 'position', [28 296 1375 506])
plot(tt, trace_control, 'k', 'linewidth', 3)
%xlim(xlims)
xlabel('time (sec)')
ylabel('LFP amplitude')
title('Control Condition')


figure
set(gcf, 'position', [28 296 1375 506])
plot(tt, trace_synapticBlockers, 'k', 'linewidth', 3)
%xlim(xlims)
xlabel('time (sec)')
ylabel('LFP amplitude')
title('With Synaptic Blockers')


figure
set(gcf, 'position', [28 296 1375 506])
plot(tt, trace_ttx, 'k', 'linewidth', 3)
%xlim(xlims)
xlabel('time (sec)')
ylabel('LFP amplitude')
title('With TTX')

figure
set(gcf, 'position', [28 296 1375 506])
plot(tt, glutamateOnly, 'k', 'linewidth', 3)
%xlim(xlims)
xlabel('time (sec)')
ylabel('LFP amplitude')
title('LFP due to synaptic transmission')

figure
set(gcf, 'position', [28 296 1375 506])
plot(tt, fiberVolley, 'k', 'linewidth', 3)
%xlim(xlims)
xlabel('time (sec)')
ylabel('LFP amplitude')
title('LFP due to fiber volley')





