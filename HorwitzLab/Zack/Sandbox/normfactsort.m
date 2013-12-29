%% normfactsort.m
datapath = 'N:\NexFiles';
dirstruct = subdir([datapath filesep '*.nex']);
numfiles = length(dirstruct);
counter = 0; normfacts = cell(numfiles,2);
for i = 1:numfiles    
    try
        stro = nex2stro(dirstruct(i).name);
    catch exception
        wh = findall(0, 'tag', 'TMWWaitbar');
        if ~isempty(wh), delete(wh); end
        continue
    end
    
    if stro.sum.paradigmID ~= 103, continue; end
    
    prelimData = NTpreprocess(stro,0,inf,inf); %prelim filter
    if ~isempty(prelimData)
        prelimData(prelimData(:,end) == 1,:) = [];
        if length(unique(prelimData(:,1))) < 20 % at least 20 color dirs
            continue
        end
    end
    
    stringeData = NTpreprocess(stro,0.4,1,0); %stringent filter
    if ~isempty(stringeData)
        if length(stringeData(:,end)) - sum(stringeData(:,end) == 1) < 10 % 10 nOOG points
            continue
        end
    else
        continue
    end

    sf = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'sf'));
    if sf(1) == 0 || isempty(sf)
        disp('No sf');
        continue
    end
    
    load('T_cones_smj10');
    
    slashidxs = strfind(dirstruct(i).name, '\');
    NTfilename = dirstruct(i).name(slashidxs(end)+1:end-4);
    
    lmsidxs = [find(strcmp(stro.sum.trialFields(1,:),'lcont')) ...
        find(strcmp(stro.sum.trialFields(1,:),'mcont')) ...
        find(strcmp(stro.sum.trialFields(1,:),'scont'))];
    
    fundamentals = stro.sum.exptParams.fundamentals;
    try
        fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
    catch exception
        continue
    end
    
    mon_spd = stro.sum.exptParams.mon_spd;
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    mon_spd = SplineRaw((380:4:780)', mon_spd, (380:5:780)');
    M = fundamentals'*mon_spd;
    M10 = T_cones_smj10*mon_spd;
    
    stro.trial(:,lmsidxs) = stro.trial(:,lmsidxs)*100*(M10/M)';
    lms = stro.trial(:,lmsidxs);
    
    Loog = logical(stringeData(:,end));
    NTscaled = stringeData(:,2:4) .* repmat(stringeData(:,5),1,3);
    
    spikeidx = strcmp(stro.sum.rasterCells(1,:),getSpikenum(stro, 'first'));
    levels = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'levelidx'));
    stepsize = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stepsize'));
    
    fpacq_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fp_acq'));
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_on'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_off'));
    
    Lbaseline = ~levels & all(lms == zeros(size(lms)),2);
    Lprobe = ~levels & ~all(lms == zeros(size(lms)),2);
    Lreplay = (levels == max(levels)) & ~stepsize;
    
    spikecounts = zeros(size(stro.trial,1),1);
    for j = 1:size(stro.trial,1)
        spiketimes = stro.ras{j,spikeidx};
        spikecounts(j) = sum(spiketimes > stimon_t(j)+stro.sum.exptParams.latency/1000 & spiketimes < stimoff_t(j));
    end
    
    [pparams,~,~,~,xformmat] = NTsurfacefit(NTscaled,Loog);
    coneweights = xformmat*pparams;
    
    Laccept = ~Lreplay & ~Lbaseline & ~Lprobe & stepsize;
    lms = lms(Laccept,:);
    spikecounts = spikecounts(Laccept);
    
    DTNTsfs = [0.5008 0.9919 1.9839 3.9677];
    err = abs(DTNTsfs-sf(1));
    whichsf = err == min(err);
    params = [1.1817	2.8475	2.8556	3.0997;
              0.1203	0.0147	-0.0035	0.0012;
              -0.0052	-0.0033	0.0072	0.0126;
              0.1445	0.2288	0.1570	0.0491;
              -0.0367	-0.5428	0.2218	-0.2422;
              -0.3620	-0.2225	0.4693	0.2491;
              0.1134	0.0706	0.0446  0.0476;
              0.7796	1.2509	-0.7186	-0.3305;
              -0.6624	-1.2128	0.5492	-0.2676;
              -0.0508	-0.0515	0.0766	-0.0105];
    params = params(:,whichsf);
    mechs = reshape(params(2:end),3,3);
    
    [th,phi,distPlane] = cart2sph(coneweights(1),coneweights(2),coneweights(3));
    wholesum = sum(abs([cos(phi)*cos(th) cos(phi)*sin(th) sin(phi)]*mechs).^params(1),2);
    distDTSurf = (1/wholesum)^(1/params(1));

    normfact = distDTSurf./distPlane;

    counter = counter + 1;
    normfacts{counter,1} = NTfilename;
    normfacts{counter,2} = normfact;
end
normfacts(counter+1:end,:) = [];
% save('normfacts0617.mat', 'normfacts');