%%
addpath('C:\Documents and Settings\horwitzlab\Desktop\MatlabCode\Analysis\FixStim');
stro = nex2stro();
fpoffidx = strcmp(stro.sum.trialFields(1,:), 'fp_off');
targidx = strcmp(stro.sum.trialFields(1,:), 'targ');
stidx = strcmp(stro.sum.trialFields(1,:), 'stim_type');
ttinfo = [stro.trial(:,targidx) stro.trial(:,stidx)];
% build logical array with columns denoting task type
ttall = [ttinfo(:,1) == 1 & ttinfo(:,2) == 0, ...
         ttinfo(:,1) == 1 & ttinfo(:,2) == 1, ...
         ttinfo(:,1) == 1 & ttinfo(:,2) == 2, ...
         ttinfo(:,1) == 0 & ttinfo(:,2) == 1, ...
         ttinfo(:,1) == 0 & ttinfo(:,2) == 2, ...
         ttinfo(:,1) == 0 & ttinfo(:,2) == 0];
     
if stro.sum.analog.storeRates{1} ~= stro.sum.analog.storeRates{2}
    keyboard
end

% looking at saccades aligned on fp off
sampleperiod = 1/stro.sum.analog.storeRates{1};
fpoffalign = floor((stro.trial(:,fpoffidx) - [stro.ras{:,3}]')/sampleperiod);
nsamples = cellfun(@length, stro.ras(:,1));
ntrials = size(stro.ras,1);
X = nan(max(nsamples),ntrials);
YHE = nan(max(nsamples),ntrials);
YVE = nan(max(nsamples),ntrials);
for i = 1:ntrials
    alignedtime = (-(fpoffalign(i)-1):(nsamples(i)-fpoffalign(i)))*sampleperiod;
    X(:,i) = [alignedtime nan(1,max(nsamples)-length(alignedtime))];
    YHE(:,i) = [stro.ras{i,1}' nan(1,max(nsamples)-length(alignedtime))];
    YVE(:,i) = [stro.ras{i,2}' nan(1,max(nsamples)-length(alignedtime))];
end

% plot like this, where # is the task type
% figure; plot(X(:,ttall(:,#)), YHE(:,ttall(:,#)))
