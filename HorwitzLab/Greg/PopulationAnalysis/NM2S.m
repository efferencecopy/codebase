% NM2S 
% Batch behavioral analyses
%
% Contents:
%
% Section 1) relationship between performance (as a function of contrast) 
% and the occurance of fixational saccades around the time of the stimulus
% presentation.
%

%%%%
% Section 1)

datapath = 'C:\PlexonData\Kali Psychophysics';
filenames = fnamesFromTxt('C:\NexFiles\Kali\NM2S.1.txt',0);
data = [];
binwidth = .1;
timebins = [-.3:binwidth:.3];  % time relative to discriminanda on (sec)
nbins = length(timebins);
for i = 1:size(filenames,1)
    stro = nex2stro([datapath,'\',filenames(i,:),'.nex']);
    % Getting all the relevant stuff out of the stro file
    ntrials = size(stro.trial,1);
    fpacqidx = find(strcmp(stro.sum.trialFields(1,:),'fp_acq'));
    stimonidx = find(strcmp(stro.sum.trialFields(1,:),'stim_on'));
    fpoffidx = find(strcmp(stro.sum.trialFields(1,:),'fp_off'));
    saccadeidx = find(strcmp(stro.sum.trialFields(1,:),'saccade'));
    fpacq_t = stro.trial(:,fpacqidx);
    stimon_t = stro.trial(:,stimonidx);
    framerate = stro.sum.exptParams.framerate;
    stimoff_t = stimon_t+stro.sum.exptParams.nframes_discriminanda/framerate;
    fpoff_t = stro.trial(:,fpoffidx);
    saccade_t = stro.trial(:,saccadeidx);
    correctidx = find(strcmp(stro.sum.trialFields(1,:),'correct'));
    correct = stro.trial(:,correctidx);
    nonmatchsideidx = find(strcmp(stro.sum.trialFields(1,:),'nonmatchside'));
    nonmatchside = stro.trial(:,nonmatchsideidx);
    deltaccs = stro.trial(:,find(strcmp(stro.sum.trialFields(1,:),'T1_Scc')));
    deltaccs(nonmatchside == 0) = stro.trial(nonmatchside == 0,find(strcmp(stro.sum.trialFields(1,:),'T2_Scc')));
    
    
    tmpdata = nan*ones(ntrials, 2+nbins);
    
    % Now getting the fixational saccade parameters
    sacstats = getSacData(stro);
    delete(gcf);
    for j = 1:ntrials
        st = sacstats.starttimes{j}-stimon_t(j);
        [n,x] = hist(st(st>timebins(1)-binwidth/2 & st<timebins(end)+binwidth/2) ,timebins);
        tmpdata(j,1:nbins) = n>0;
        tmpdata(j,nbins+1) = deltaccs(j);
        tmpdata(j,nbins+2) = correct(j);
    end
    data = [data; tmpdata];
    size(data)
end

% Now looking at the data.
figure;axes; hold on;
for i = 1:length(timebins)
    r = corrcoef(data(:,[i end]));
    plot(timebins(i),r(1,2),'ko');
end

uniquecontrasts = unique(data(:,end-1));
figure; 
for whichtimebin = 1:length(timebins) 
    subplot(3,3,whichtimebin); hold on;
    for i = 1:length(uniquecontrasts)
        L = data(:,end-1) == uniquecontrasts(i);
        r = corrcoef(data(L,[whichtimebin end]));
        plot(uniquecontrasts(i),r(1,2),'ko');
    end
    set(gca,'YLim',[-.3 .3]);
end