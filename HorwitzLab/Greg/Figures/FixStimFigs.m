% FixStim figures for talks
%
% Section 1)
% Rasters and PSTHs

%
%%
% Section 1
% Rasters and PSTH
%stro = nex2stro(findfile('S081111018.nex'));
stro = nex2stro(findfile('S081111026.nex'));
stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
fpoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fpoff_t'));
stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimoff_t'));
optfreq = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'elec_stimfreq'));
targ_shown = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_shown'));
stim_type = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_type'));

offset = [-.15 .4];  % pre and post time
spikechan = strcmp(stro.sum.rasterCells,'sig001a');
spikes = stro.ras(:,spikechan);
binwidth = .005;
bins = offset(1):binwidth:offset(2);
PSTH = zeros(1,length(bins));
L = targ_shown == 0 & stim_type == 2; % No target, optostim
trlidxs = find(L);

figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);

% Panel 1: Rasters
axes('position',[2 6 3 2]); hold on;
for counter = 1:sum(L)
    trlidx = trlidxs(counter);
    tmpspikes = spikes{trlidx}-stimon_t(trlidx);
    tmpspikes(tmpspikes < offset(1) | tmpspikes > offset(2)) = [];
    nspikestot = length(tmpspikes);
    plot([tmpspikes tmpspikes]',[zeros(nspikestot,1) .7*ones(nspikestot,1)]'+counter,'k-');
    PSTH = PSTH + hist(tmpspikes, bins);
end
set(gca,'Xlim',offset,'Ylim',[0 counter+1],'XTick',[],'YTick',[])

% Panel 2: PSTHs
axes('position',[2 4 3 1]); hold on; % PSTHs
firingrate = PSTH/sum(L)/binwidth;
plot(bins,firingrate,'k-','LineWidth',2)
set(gca,'XLim', offset,'Ytick',[0 floor(max(firingrate))],'YLim',[0 max(firingrate)+5],'XTick',[]);

% Panel 3: Pulse sequence
axes('position',[2 3 3 .5]); hold on;
dur = mode(stimoff_t-stimon_t);
secspercycle = 1/unique(optfreq(optfreq > 0));
transitions = 0:secspercycle/2:dur;
x = [transitions; transitions];
x = [offset(1); x(:); offset(2)];
y = [repmat([0 1],1,length(transitions)/2) 0];
y = [y;y];
if (length(x(:)) == length(y(:)))
    plot(x,y(:)','k-','linewidth',2);
end
set(gca,'XLim', offset,'Ytick',[]);