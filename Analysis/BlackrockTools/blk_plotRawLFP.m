function blk_plotRawLFP(stro, nsxtype, lp_cutoff, plotall)

% find the LFP channels
searchString = ['chan\w*', nsxtype];
LFPchannels = cellfun(@(x) ~isempty(x), regexpi(stro.sum.rasterFields, searchString));
LFPchannels = find(LFPchannels);


dataOnTidx = strcmpi(stro.sum.trialFields, [nsxtype,'_start_sec']);
stimonTidx = stro.sum.idx.stimon;
sampRate_nsx = stro.sum.(nsxtype).MetaTags.SamplingFreq;
preSamps = round(0.100 .* sampRate_nsx);
postSamps = round(0.200.* sampRate_nsx);
Ntrials = size(stro.trial, 1);
figure, hold on,
nplts = ceil(sqrt(numel(LFPchannels)));
for i_ch = 1:numel(LFPchannels)
    
    chdata = [];
    for i_trl = 1:Ntrials
        
        raw = stro.ras{i_trl, LFPchannels(i_ch)};
        raw = butterfilt(raw, lp_cutoff, sampRate_nsx, 'low');
        tstart = stro.trial(i_trl, dataOnTidx);
        tt = (0:numel(raw)-1) ./ sampRate_nsx;
        tt = tt + tstart; % now in seconds from beginning of data file
        
        stimont = stro.trial(i_trl, stimonTidx);
        tt = tt - stimont; % now relative to stim onset
        [~, zeroidx] = min(abs(tt-0));
        
        bkgnd = mean(raw(1:zeroidx));
        raw = raw - bkgnd;
        
        idx = [zeroidx - preSamps : zeroidx + postSamps];
        raw = raw(idx);
        chdata = cat(1, chdata, raw);
        
    end
    
    subplot(nplts, nplts, i_ch)
    tt = tt(idx).*1000;
    if plotall
        plot(tt, chdata)
    else
        plot(tt, mean(chdata,1))
    end
    axis tight
    xlim([-10, 95])
    title(sprintf('CH: %d', i_ch));
    
end



