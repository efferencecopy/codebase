%%
NTfilename = 'S030510003';
NT = nex2stro(findfile(NTfilename));
out = NTpreprocess(NT, .4, 1);
out(out(:,end) == 1,:) = [];  % Getting rid of the OOGamuts
spikeidx = strcmp(NT.sum.rasterCells(1,:), getSpikenum(NT));
coloridxs = NT.trial(:,strcmp(NT.sum.trialFields(1,:),'coloridx'));
stimon_t = NT.trial(:,strcmp(NT.sum.trialFields(1,:),'stim_on'));
stimoff_t = NT.trial(:,strcmp(NT.sum.trialFields(1,:),'stim_off'));
lms = NT.trial(:,[find(strcmp(NT.sum.trialFields(1,:),'lcont'))...
            find(strcmp(NT.sum.trialFields(1,:),'mcont'))...
            find(strcmp(NT.sum.trialFields(1,:),'scont'))]);

for i = 1:size(NT.trial,1)
    spiketimes = NT.ras{i,spikeidx};
    nspikes(i,:) = sum(spiketimes > stimon_t(i)+NT.sum.exptParams.latency/1000 && spiketimes < stimoff_t(i));
end

uniquecoloridxs = out(:,1);
uniquecolordirs = out(:,2:4);

figure; subplot(2,1,1); hold on;
data = [];

for i = uniquecoloridxs'
    L = coloridxs == i;
    contrast = abs(uniquecolordirs(uniquecoloridxs == i,:)*lms(L,:)');
    data(length(data)+1).contrast = contrast';
    data(length(data)).response = nspikes(L);
    data(length(data)).lms = uniquecolordirs(uniquecoloridxs == i,:);
    
    h = plot(contrast,nspikes(L),'ko');
    set(h,'MarkerFaceColor',unifrnd(0,1,3,1));
end