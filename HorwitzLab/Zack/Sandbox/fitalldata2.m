%% brute force ugliness
cellcount = 1;
initguesses = combn(linspace(-250,250,20),3); % all possible combinations: 20^3 of them
modelouts = nan(size(initguesses,1),4,14);
errorouts = nan(size(initguesses,1),14);
for s = {{'K062909007' 'K062909008'}, {'K070209006' 'K070209007'}, ...
        {'K071709005' 'K071709006'}, {'K071709007' 'K071709008'}, ...
        {'K080709005' 'K080709006'}, {'K082709003' 'K082709004'}, ...
        {'K101409002' 'K101409003'}, {'K101509004' 'K101509005'}, ...
        {'K102809002' 'K102809003'}, {'K010110003' 'K010110005'}, ...
        {'S012710002' 'S012710003'}, {'S030510002' 'S030510003'}, ...
        {'S070910003' 'S070910004'}, {'S091310001' 'S091310003'}}
    
    stro = nex2stro(findfile(s{1}{2}));
    
    spikeidx = strcmp(stro.sum.rasterCells(1,:),getSpikenum(stro, 'first'));
    levels = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'levelidx'));
    stepsize = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stepsize'));
    
    lmsidxs = [find(strcmp(stro.sum.trialFields(1,:),'lcont')) ...
        find(strcmp(stro.sum.trialFields(1,:),'mcont')) ...
        find(strcmp(stro.sum.trialFields(1,:),'scont'))];
    
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_on'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_off'));
    
    lms = stro.trial(:,lmsidxs);
    
    Lbaseline = ~levels & all(lms == zeros(size(lms)),2);
    Lprobe = ~levels & ~all(lms == zeros(size(lms)),2);
    Lreplay = (levels == max(levels)) & ~stepsize;
    
    spikecounts = zeros(size(stro.trial,1),1);
    for i = 1:size(stro.trial,1)
        spiketimes = stro.ras{i,spikeidx};
        spikecounts(i) = sum(spiketimes > stimon_t(i)+stro.sum.exptParams.latency/1000 & spiketimes < stimoff_t(i));
    end
    
    opt = 'baseline';
    err = 'fish';
    
    Laccept = ~Lreplay & ~Lbaseline & ~Lprobe & stepsize;
    newlms = lms(Laccept,:);
    newsc = spikecounts(Laccept);
    
    bline = nan;
    optsbase = optimset('MaxFunEvals',5e3,'MaxIter',5e3,'TolFun',1e-9,'TolX',1e-9,'Display','off');
    
    for i = 1:size(modelouts,1)
        initguess = initguesses(i,:);
        switch opt
            case 'abs'
                [xout errorouts(i,cellcount)] = ...
                    fminsearch(@(x) modelfun(x,newlms,newsc,bline,opt,err),initguess(:),optsbase);
                modelouts(i,1:3,cellcount) = xout';
            case 'exponent'
                options = optimset(optsbase,'LargeScale','off');
                [xout errorouts(i,cellcount)] = ...
                    fminunc(@(x) modelfun(x,newlms,newsc,bline,opt,err),[initguess(:);2],options);
                modelouts(i,:,cellcount) = xout';
            case 'baseline'
                options = optimset(optsbase,'Algorithm','sqp');
                [xout errorouts(i,cellcount)] = ...
                    fmincon(@(x) modelfun(x,newlms,newsc,bline,opt,err),[initguess(:);0],...
                    [],[],[],[],[-inf(3,1);0],[inf(3,1);max(newsc)],[],options);
                modelouts(i,:,cellcount) = xout';
        end
    end
    
    fprintf('done with cell %d\n', cellcount);
    cellcount = cellcount + 1;
end