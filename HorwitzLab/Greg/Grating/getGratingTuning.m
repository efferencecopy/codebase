function out = getGratingTuning(stro, spikenum)
% function out = getGratingTuning(stro)
%
% Takes a stro file from the grating paradigm (#150) and returns a
% structure full of tuning parameters.
% GDLH 4/27

if (nargin == 1)
    spikenum = strcmp(stro.sum.rasterCells(1,:),getSpikenum(stro));
end

orients = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'orient'));
sfs = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'sf'));
diams = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'diam'));
protocols = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'protocol'));
framerate = stro.sum.exptParams.framerate;
nframes = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'nframes'));
stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_on'));
stimoff_t= stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_off'));
Lcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'lcont'));
Mcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'mcont'));
Scc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'scont'));
spikerates = [];
baselines = [];  baseline_t = 0.25;
for i = 1:size(stro.trial,1)
    spiketimes = stro.ras{i,spikenum};
    nspikes = sum(spiketimes > stimon_t(i) & spiketimes < stimoff_t(i));
    spikerates = [spikerates; nspikes./(stimoff_t(i)-stimon_t(i))];
    nspikes = sum(spiketimes > stimon_t(i)-baseline_t & spiketimes < stimon_t(i));
    baselines = [baselines; nspikes./baseline_t];
end
out.baselines = [mean(baselines) std(baselines) size(stro.trial,1)];

lastorienttrial = find(protocols == 1,1,'last');
if (all(protocols(1:lastorienttrial) == 1))
    firstorienttrial = 1;
else
    firstorienttrial = find(protocols(1:lastorienttrial) ~= 1,1,'last')+1;
end
lastSFtrial = find(protocols == 2,1,'last');
if (all(protocols(1:lastSFtrial) == 1))
    firstSFtrial = 1;
else
    firstSFtrial = find(protocols(1:lastSFtrial) ~= 2,1,'last')+1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orientation
% Getting the preferred orientations
if (any(protocols == 1))
    trlidxs = [firstorienttrial:lastorienttrial];
    x = orients(trlidxs);
    y = spikerates(trlidxs);
    Ltmp = x == min(x);
    y = [y; y(Ltmp)];
    x = [x; x(Ltmp)+2*pi];
    pp = csape(x,y,'periodic');
    xx = linspace(0,2*pi,100);
    fit = ppval(pp,xx);
    out.orient.preforient = xx(find(fit == max(fit),1));
    
    % Getting means, standard deviations, and n's for orientation tuning.
    out.orient.stim = unique(x);
    for j = out.orient.stim'
        out.orient.resp(j == out.orient.stim,1) = mean(y(x==j));
        out.orient.resp(j == out.orient.stim,2) = std(y(x==j));
        out.orient.resp(j == out.orient.stim,3) = sum(x==j);
    end
else
    out.orient.preforient = nan;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SFs
% Getting the preferred SF
if (any(protocols == 2))
    trlidxs = [firstSFtrial:lastSFtrial];
    x = sfs(trlidxs);
    y = spikerates(trlidxs);
    pp = csape(x,y,'variational');
    xx = linspace(min(sfs),max(sfs),100);
    fit = ppval(pp,xx);
    out.sf.prefSF = xx(find(fit == max(fit),1));
    
    % Getting means, standard deviations, and n's for spatial frequency tuning.
    out.sf.stim = unique(x);
    for j = out.sf.stim'
        out.sf.resp(j == out.sf.stim,1) = mean(y(x==j));
        out.sf.resp(j == out.sf.stim,2) = std(y(x==j));
        out.sf.resp(j == out.sf.stim,3) = sum(x==j);
    end
else
    out.sf.prefSF = nan;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Looking at area summation curve
if (any(protocols == 5))
    Lprotocol = protocols == 5;
    x = diams(Lprotocol);
    y = spikerates(Lprotocol);
    pp = csape(x,y,'not-a-knot');
    xx = linspace(min(diams),max(diams),100);
    fit = ppval(pp,xx);
    out.areasummation.prefsize = xx(find(fit == max(fit),1));
    
    % Getting means, standard deviations, and n's for area summation curve.
    out.areasummation.stim = unique(x);
    for j = out.areasummation.stim'
        out.areasummation.resp(j == out.areasummation.stim,1) = mean(y(x==j));
        out.areasummation.resp(j == out.areasummation.stim,2) = std(y(x==j));
        out.areasummation.resp(j == out.areasummation.stim,3) = sum(x==j);
    end
else
    out.areasummation.prefsize = nan;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Color
% Using a new model fitting procedure  8/20/09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L4 = protocols == 4;
errorflag = 0;
if sum(L4) > 8
    colordirections = [Lcc(L4) Mcc(L4) Scc(L4)];
    uniquecolordirections = unique(colordirections,'rows');
    normresps = []; rawresps = []; sds = []; ns = [];
    for i = 1:size(uniquecolordirections,1)
        cdir = uniquecolordirections(i,:);
        L = Lcc == cdir(1) &  Mcc == cdir(2) & Scc == cdir(3);
        normresps = [normresps; mean(spikerates(L&L4))];   % NOT normalizing by cone contrast
        rawresps = [rawresps;mean(spikerates(L&L4))];
        sds = [sds; std(spikerates(L&L4))];
        ns = [ns; sum(L&L4)];
    end
    responses = spikerates(L4);
  %  maxrespidx = find(responses == max(responses),1);
  %  initguess = colordirections(maxrespidx,:);
  %  initguess = initguess.*(responses(maxrespidx)./norm(initguess));
  % Using each tested color direction as an initial guess.
    options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-6,'TolX',10^-6);
    minfval = Inf; conweights = [nan nan nan];
    for i = 1:size(uniquecolordirections,1)
        initguess = uniquecolordirections(i,:);
        [tmpconeweights, tmpfval, exitflag] = fminsearch(@(x) linmodfiterr(colordirections, responses, x), initguess, options);
        if (tmpfval < minfval)
            %disp([num2str(tmpconeweights),'   ',num2str(minfval-tmpfval)]);
            minfval = tmpfval;
            coneweights = tmpconeweights;
        end
    end
    %return the raw data (always) and the model fit (if appropriate)
    out.color.colors = uniquecolordirections;
    out.color.colresp = [rawresps sds ns];
    if (exitflag)
        out.color.prefcolor = coneweights;
    else
        out.color.prefcolor = [nan nan nan];
    end
else
    out.color.prefcolor = [nan nan nan];
    out.color.colors = nan;
    out.color.colresp = nan;
end


%%%%%%%%%%%%%
% F1/F0
% A nan can mean that there were no conditions in which the firing rate
% exceeeded MINFRTHRESH.
MINFRTHRESH = 10;
out.modulationratio = nan; % gets changed if there's enough data to compute
if (any(L4))
    binwidth = .001;
    start_t = stimon_t(L4);
    end_t = stimoff_t(L4);
    minstimdur = min(end_t-start_t);
    bins = [0:binwidth:minstimdur];  % Bad assuming stimulus dur
    colordirections = [Lcc(L4) Mcc(L4) Scc(L4)];
    uniquecolordirections = unique(colordirections,'rows');
    tmp =[];
    for i = 1:size(uniquecolordirections,1)
        L = L4 &...
            Lcc == uniquecolordirections(i,1) & ...
            Mcc == uniquecolordirections(i,2) &...
            Scc == uniquecolordirections(i,3);
        start_t = stimon_t(L);
        spiketimes = stro.ras(L,1);
        totalspikes = [];
        for j = 1:sum(L)
            totalspikes = [totalspikes; spiketimes{j}-start_t(j)];
        end
        psth = histc(totalspikes',bins)/binwidth/sum(L);

        k = unique(stro.trial(L,strcmp(stro.sum.trialFields(1,:), 'tf')));
        if (length(k) > 1)
            error('Too many temporal frequencies!');
        end
        basis0 = ones(length(psth),1);
        basis1 = exp(-2*pi*sqrt(-1)*k*[0:length(psth)-1]/length(psth))';
        F0 = abs((psth-mean(baselines))*basis0);
        F1 = 2*abs((psth-mean(baselines))*basis1);
        % Not sure what the logic for the the "2*" in the line above
        % but it makes a half-wave rectified sinewave have a modulation
        % ratio of pi/2, which is what it's supposed to have.
        tmp = [tmp; F1 F0];
    end
    L = tmp(:,2)/length(psth) > MINFRTHRESH;  % Threshold on mean firing rate
    if any(L)
        out.modulationratio =  mean(tmp(L,1)./tmp(L,2));
        out.modulationratio =  geomean(tmp(L,1)./tmp(L,2));
    else
        out.modulationratio = nan;
    end
end
    
    
    % if any(strcmp('AD11',stro.sum.rasterCells))  % i.e, eye position saved
%     sacstats = getSacData(stro);
%     close;
%     
%     bindur = .0025;  % The precision of Plexon?
%     
%     % Setting up an L vector of trials to analyze
%     % Only looking at conditions in protocol 4 with > 10 sp/sec
%     tfs = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'tf'));
%     L = protocols == 4 & spikerates-mean(baselines) > 10;
%     trlidxs = find(L);
%     
%     % Now doing the analysis on a trial by trial basis
%     F0 = [];
%     F1 = [];
%     epochdurs = [];
%     for i = trlidxs'
%         startsac_t = sacstats.starttimes{i};
%         endsac_t = sacstats.endtimes{i}+.1; % cutting out from sac start to 100 ms after sac end
%         L = startsac_t < stimon_t(i) | endsac_t > stimoff_t(i);
%         startsac_t(L) = [];
%         endsac_t(L) = [];
%         % Counterintuitively, we want to start epochs at the *ends*
%         % of saccades and end them at the *beginings* of saccades.
%         start_t = [stimon_t(i); endsac_t];
%         end_t = [startsac_t;stimoff_t(i)];
%         
%         % getting intersaccade times
%         epochdur = end_t-start_t;
%         
%         for j = 1:length(epochdur)
%             if (epochdur(j) < 1/tfs(i))  % Not taking short epochs
%                 continue;
%             end
%             psth = histc(stro.ras{i},[start_t(j):bindur:end_t(j)]);
%             psth = psth/bindur;  % To get back to sp/sec (which is the units of "baseline").
%             k = (end_t(j)-start_t(j))/(1/tfs(i)); % k is the number of cycles in the snippet
%             basis1 = exp(-2*pi*sqrt(-1)*k*[0:length(psth)-1]/length(psth))';
%             F0(length(F0)+1) = abs(sum(psth-mean(baselines)));
%             F1(length(F1)+1) = abs(sum((psth-mean(baselines))'*basis1));
%             epochdurs(length(epochdurs)+1) = epochdur(j);
%         end
%     end
%     
%     %only compute the modulation ratio if there is sufficient data.
%     %Otherwise, return a nan (initialized to a nan above).
%     if ~any(isempty([F0, F1, epochdurs]))
%         %out.modulationratio = nanmean(F1./F0);
%         out.modulationratio = (epochdurs./sum(epochdurs))*(F1./F0)';
%     end
%     
% end

%%%%%%%%%%%%%
% Eccentricity
x = stro.sum.exptParams.rf_x ./ 10; %represented in the stro file in 10ths of degrees?
y = stro.sum.exptParams.rf_y ./ 10;
out.eccentricity = sqrt(x.^2 + y.^2);

%%%%%%%%%%%%%
% Contrast Response Function (if present). Returns a structure.
% xxx.crf.rates is a 2x1 cell array. Each row of the array represents a
% distinct color direction. Entries are the firing rates for each contrast.
% xxx.colors gives an index into xxx.crf to determine which row is which
% color.
l_crf = protocols == 7;
tColors = [Lcc, Mcc, Scc];
crfColors = unique(sign(tColors(l_crf,:)), 'rows');
crfColors(sum(abs(crfColors), 2) == 0, :) = []; %remove the zero contrast condition
out.crf.colors = crfColors;
for clr = 1:size(crfColors, 1);
    l_color = ismember(sign(tColors), crfColors(clr,:), 'rows');
    tmp = unique(tColors(l_color & l_crf, :), 'rows');
    contrasts = zeros(size(tmp,1)+1, size(tmp,2));
    contrasts(2:end,:) = tmp; %by indexing 2:end I'm explictly adding the zero contrast condition
    out.crf.norms{clr, 1} = sqrt(sum(contrasts.^2, 2));
    for cntrst = 1:size(contrasts,1)
        l_cntrst = softEq(contrasts(cntrst, :), tColors, 7, 'rows');
        if cntrst == 1 %ie the zero contrast condition
            l_trials = l_crf & l_cntrst;
        else
            l_trials = (l_crf & l_color) & l_cntrst;
        end
        out.crf.rates{clr,1}{cntrst} = spikerates(l_trials);
    end
end




