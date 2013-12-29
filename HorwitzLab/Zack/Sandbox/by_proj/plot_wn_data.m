%%
% Section 10
%
% STAs and PC1s for a bunch of cells.

[fnames, spikeIdx] = fnamesFromTxt2('C:\NO BACKUP\NexFiles\nexfilelists\Greg\WhiteNoise\BYcandidates.txt');
%fnames = {{'K042109005' 'K042109009' 'K042109012'}};
%spikeIdx = 1;
nExpts = size(fnames, 1);
h1 = figure;
h2 = figure;
spike_channel_names = {'sig001a' 'sig001b'};
for a = 1:nExpts;
    stro = {};
    for i = 1:size(fnames{a},2)
        tmpstro = nex2stro(findfile(char(fnames{a}(i))));
        if (isempty(stro))
            stro = tmpstro;
        else
            stro = strocat(stro, tmpstro);
        end
    end
    maxT = 8;
    framerate = stro.sum.exptParams.framerate;
    nstixperside = stro.sum.exptParams.nstixperside;
    ntrials = length(stro.sum.absTrialNum);
    stimonidx = find(strcmp(stro.sum.trialFields(1,:),'stim_on'));
    stimoffidx = find(strcmp(stro.sum.trialFields(1,:),'all_off'));
    nframesidx = find(strcmp(stro.sum.trialFields(1,:),'num_frames'));
    noisetypeidx = find(strcmp(stro.sum.trialFields(1,:),'noise_type'));
    sigmaidxs = strmatch('sigma',stro.sum.trialFields(1,:));

    L = stro.trial(:,noisetypeidx) == 1;
    stro.ras(~L,:) = [];
    stro.trial(~L,:) = [];

    out = getWhtnsStats(stro,maxT,'STCOVmex', {nstixperside^2, 3, maxT}, spike_channel_names{spikeIdx(a)});
    STAs = out{1};
    STCs = out{2};

    energy = sum(STAs.^2);
    whichframe = find(energy == max(energy));
    subSTA = STAs(:,whichframe);
    P = eye(length(subSTA))-subSTA/(subSTA'*subSTA)*subSTA';

    STC = reshape(STCs(:,whichframe),sqrt(size(STCs,1)),sqrt(size(STCs,1)));
    [PC1,d] = eigs(P*STC*P',1);

    % Normalizing images
    template = reshape([1:nstixperside^2],nstixperside,nstixperside);
    edgepixels = [template(:,1); template(1,[2:end-1])'; template(end,[2:end-1])'; template(:,end)];
    edgepixelidxs = [edgepixels; edgepixels+nstixperside^2; edgepixels+2*(nstixperside^2)];
    PCelements = PC1(edgepixelidxs,:,:);
    PCsds = std(PCelements);    % One std calculated per PC
    PC1 = PC1.*repmat(std(STAs(:,1))./PCsds,[size(PC1,1),1,1]);

    rowidxs = reshape([1:3*nstixperside^2],[nstixperside^2 3]);
    maxes = []; mins = [];
    imagevectors = [STAs(:,whichframe), PC1];
    for i = 1:3
        maxes = [maxes; max(max(imagevectors(rowidxs(:,i),:)))];
        mins = [mins; min(min(imagevectors(rowidxs(:,i),:)))];
    end
    potentialnormfactors = [(1-[.5; .5; .5]-eps)./(maxes-[.5; .5; .5]); (-[.5; .5; .5]+eps)./(mins-[.5; .5; .5])];
    % 'eps' in above line is a kludge that is required for avoiding
    % out of bounds errors.
    potentialnormfactors(potentialnormfactors < 0) = []; % if min > mu or max < mu
    normfactor = min(potentialnormfactors);
    muvect = reshape(repmat([.5 .5 .5],nstixperside^2,1),nstixperside^2*3,1);

    % Plotting
    figure(h1);
    subplot(ceil(sqrt(nExpts)),ceil(sqrt(nExpts)),a);
    STA = normfactor*(STAs(:,whichframe)-muvect)+muvect;
    STA = reshape(STA,[nstixperside nstixperside 3]);
    image(STA);
    set(gca,'XTick',[],'YTick',[]); axis square;

    figure(h2);
    subplot(ceil(sqrt(nExpts)),ceil(sqrt(nExpts)),a);
    v = normfactor*(PC1-muvect)+muvect;
    v = reshape(v,[nstixperside nstixperside 3]);
    image(v);
    set(gca,'XTick',[],'YTick',[]); axis square;

    drawnow;
end