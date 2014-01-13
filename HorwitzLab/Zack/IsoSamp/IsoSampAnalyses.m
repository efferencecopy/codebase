%%
% stro = nex2stro();
stro = stro3;
trialFields = stro.sum.trialFields(1,:);
stimon_t = stro.trial(:,strcmp(trialFields, 'stimon_t'));
stimoff_t = stro.trial(:,strcmp(trialFields, 'stimoff_t'));
anlgstart_t = stro.ras(:,strcmp(stro.sum.rasterCells, 'anlgStartTime'));
spikes = stro.ras(:,strcmp(stro.sum.rasterCells, 'sig001a'));
sf = stro.trial(:,strcmp(trialFields, 'sf'));
num_stims = min(stro.sum.exptParams.num_stims, size(stro.trial, 1));

lms = [stro.trial(:,strcmp(trialFields, 'stim_l')) ...
       stro.trial(:,strcmp(trialFields, 'stim_m')) ...
       stro.trial(:,strcmp(trialFields, 'stim_s'))] ./ 100; % stro1,2 don't require this
   
[unqLms,~,unqIdxs] = unique(lms, 'rows');
   
binwidth = .005;
stim_dur = mode(stimoff_t - stimon_t);
offset = [-0.1 stim_dur + 0.1]; % 100 ms window
bins = offset(1):binwidth:offset(2);
PSTH = zeros(size(bins));
stim_resp = zeros(num_stims, 1);
num_trials = size(spikes, 1);
theseSpikes = cell(size(spikes));

for trialIdx = 1:num_trials
    tempSpikes = spikes{trialIdx} - stimon_t(trialIdx);
    tempSpikes = tempSpikes(tempSpikes >= offset(1) ...
                              & tempSpikes <= offset(2));
    theseSpikes{trialIdx} = tempSpikes;
end

figure(); axes(); hold on;
for stimIdx = 1:num_stims
    Lstims = unqIdxs == stimIdx;
    if sum(Lstims) > 0
        spikeTimes = cat(1, theseSpikes{Lstims});
        stim_resp(stimIdx) = numel(spikeTimes);
        plot([spikeTimes spikeTimes]', ...
             [zeros(stim_resp(stimIdx), 1) .5 * ones(stim_resp(stimIdx), 1)]' + stimIdx, 'k-');
        PSTH = PSTH + hist(spikeTimes, bins);
    end
end

set(gca, 'XLim', offset, 'YTick', [], 'YLim', [-4 num_stims + 5]);
plot([0 0], get(gca, 'YLim'), 'm-', ...
     [stim_dur stim_dur], get(gca, 'YLim'), 'm-');

PSTH = PSTH / binwidth / num_trials;

figure(); axes(); hold on;
plot(bins, PSTH, 'k-', 'linewidth', 2);
set(gca, 'YLim', [0 10 * ceil(max(PSTH/10))]);
set(gca, 'Xlim', offset);
plot([0 0], get(gca, 'YLim'), 'm-', ...
     [stim_dur stim_dur], get(gca, 'YLim'), 'm-');
 
figure(); axes(); hold on;
% scatter3(unqLms(:,1), unqLms(:,2), unqLms(:,3), 8, stim_resp/max(stim_resp), 'o', 'filled');
cmap = jet(64);
for i = 1:num_stims
    plot3(unqLms(i,1), unqLms(i,2), unqLms(i,3), 'o', 'markersize', 5, ...
        'markerfacecolor', cmap(ceil(63*stim_resp(i)/max(stim_resp))+1,:), ...
        'markeredgecolor', 'none');
end
xlabel('L'); ylabel('M'); zlabel('S');

%%
% Greg's attempts to find a "preferred direction"
% [b,bint] = regress(stim_resp, [ones(num_stims, 1) abs(unqLms)]);
% plot3([0 b(2)/1e4], [0 b(3)/1e4], [0 b(4)/1e4], 'm-', 'linewidth', 3);
% plot3([0 -b(2)/1e4], [0 -b(3)/1e4], [0 -b(4)/1e4], 'm-', 'linewidth', 3);
maxidx = find(max(stim_resp) == stim_resp);
options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-6,'TolX',10^-6);
[tmpconeweights, tmpfval, exitflag] = fminsearch(@(x) linmodfiterr(unqLms, stim_resp, x), unqLms(maxidx,:),options);
plot3([0 tmpconeweights(1)],[0 tmpconeweights(2)],[0 tmpconeweights(3)],'k-');

figure; axes; hold on;
projs = unqLms*tmpconeweights';
L = projs<0;
plot(abs(projs(L)),stim_resp(L),'ko');
plot(abs(projs(~L)),stim_resp(~L),'ro');
