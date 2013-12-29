% For Jonathan Pillow
filename = 'K033111002';
stro = nex2stro(findfile(filename));

spikeidx = strcmp(stro.sum.rasterCells(1,:),getSpikenum(stro,'first'));
Lbldone = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'bldone'));
coloridxs = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'coloridx'));
levels = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'levelidx'));
reversals = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'reversal'));
stepsize = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stepsize'));

lmsidxs = [find(strcmp(stro.sum.trialFields(1,:),'lcont'))...
            find(strcmp(stro.sum.trialFields(1,:),'mcont'))...
            find(strcmp(stro.sum.trialFields(1,:),'scont'))];

fpacq_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fp_acq'));
stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_on'));
stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_off'));
eot_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'eot'));

lms = stro.trial(:,lmsidxs);
lat = stro.sum.exptParams.latency;
mono = nanmax([0 stro.sum.exptParams.monopolar]);

Lbaseline = ~levels & all(lms == zeros(size(lms)),2);
Lprobe = ~levels & ~all(lms == zeros(size(lms)),2);

nspikes = zeros(size(stro.trial,1),1);
for i = 1:size(stro.trial,1)
    spiketimes = stro.ras{i,spikeidx};
    nspikes(i) = sum(spiketimes > stimon_t(i)+stro.sum.exptParams.latency/1000 & spiketimes < stimoff_t(i));
end

L = ~Lprobe & ~Lbaseline;

data = [];
data.filename = filename;
data.staircaseidx = coloridxs(L);
data.lms = lms(L,:);
data.nspikes = nspikes(L);
data.stimon_t = stimon_t(L);
