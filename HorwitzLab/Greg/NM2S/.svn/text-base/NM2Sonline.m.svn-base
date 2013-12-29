function p = NM2Sonline()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Online analysis stuff for NM2S paradigm
% GDLH 1/26/09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               
% Some initializations
global udpCom;  % Needed by "DealWithMessage()" and associated subfunctions
    udpCom.sock = [];
    udpCom.port = 6665;
    udpCom.rexip = '192.168.1.120';

SetUpFigure;
C = NM2SCodes;        % Script that defines a bunch of constants
p = InitPStruct(C.EOTCD); % Structure of the events from each trial
s = InitPlex();     % Link to Plexon datastream
InitStatsStruct();  % Structure that holds the statistics

[sock, Success] = pnetStart(6665);  % UDP communication with REX
stopnow = 0;    % Flag to break out of the endless while loop (hit ESC to set to 1)

% The main loop
while (~stopnow)
    stopnow = CheckForESCKey();
 %    if stopnow
 %        keyboard
 %    end
    % Check for a message from REX
    msgSize = pnet(udpCom.sock, 'readpacket', 200, 'noblock');
    if(msgSize)
        stopnow = DealWithMessage(msgSize);
    end

    [n, eventList] = PL_GetTS(s);
    if (n > 0)
       p = ProcessEventList(p, eventList);
    end
    
    if (~p.processnowflag)
        continue;
    end
    gotfulltrial = GetTrialParams(p);
    if (~gotfulltrial)
        p = CleanUpEvents(p); % Remove aborted trials
        continue;
    end
    IntegrateTrial;
    PlotFig;
    p = CleanUpEvents(p);
end
PL_Close(s);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Statistics that get updated every trial
function InitStatsStruct()
    a = get(gcf,'UserData');
    
    a.samplel = nan;
    a.samplem = nan;
    a.samples = nan;
    a.t1l = nan;
    a.t1m = nan;
    a.t1s = nan;
    a.t2l = nan;
    a.t2m = nan;
    a.t2s = nan;
    a.grad1l = nan;
    a.grad1m = nan;
    a.grad1s = nan;
    a.grad2l = nan;
    a.grad2m = nan;
    a.grad2s = nan;
    a.correct = nan;
    a.mstim = nan;
    a.storage.trialidentifiers = [];
    a.storage.correct = [];
    a.storage.ntrials = [];
    a.trialidentstrings = {'nonmatchl', 'nonmatchm', 'nonmatchs', 'gradl', 'gradm', 'grads', 'mstim'};
    
    set(gcf,'UserData',a);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Getting the trial parameters from the ecodes
function gotfulltrial = GetTrialParams(p)
    a = get(gcf,'UserData');
    C = NM2SCodes;
    if (any(p.events == C.ABORTCD))
        gotfulltrial = 0;
        return;
    end
    try
        a.samplel = GetVal(C.SAMPLELCD, p.events, 'float');
        a.samplem = GetVal(C.SAMPLEMCD, p.events, 'float');
        a.samples = GetVal(C.SAMPLESCD, p.events, 'float');
        a.t1l = GetVal(C.T1LCD, p.events, 'float');
        a.t1m = GetVal(C.T1MCD, p.events, 'float');
        a.t1s = GetVal(C.T1SCD, p.events, 'float');
        a.t2l = GetVal(C.T2LCD, p.events, 'float');
        a.t2m = GetVal(C.T2MCD, p.events, 'float');
        a.t2s = GetVal(C.T2SCD, p.events, 'float');
        a.grad1l = GetVal(C.GRAD1LCD, p.events, 'float');
        a.grad1m = GetVal(C.GRAD1MCD, p.events, 'float');
        a.grad1s = GetVal(C.GRAD1SCD, p.events, 'float');
        a.grad2l = GetVal(C.GRAD2LCD, p.events, 'float');
        a.grad2m = GetVal(C.GRAD2MCD, p.events, 'float');
        a.grad2s = GetVal(C.GRAD2SCD, p.events, 'float');
        a.nonmatchside = GetVal(C.NONMATCHSIDECD, p.events);
        a.correct = GetVal(C.CORRECTCD, p.events);
        a.mstim = GetVal(C.MSTIMCD, p.events);
        if (a.nonmatchside == 0)
            a.nonmatchl = a.t2l-a.samplel;
            a.nonmatchm = a.t2m-a.samplem;
            a.nonmatchs = a.t2s-a.samples;
            a.gradl = a.grad2l; a.gradm = a.grad2m; a.grads = a.grad2s;
        else
            a.nonmatchl = a.t1l-a.samplel;
            a.nonmatchm = a.t1m-a.samplem;
            a.nonmatchs = a.t1s-a.samples;
            a.gradl = a.grad1l; a.gradm = a.grad1m; a.grads = a.grad1s;
        end
        stimon_t = p.times(find(p.events == C.STIMONCD,1));
        stimoff_t = p.times(find(p.events == C.STIMOFFCD,1));
        spiketimes = p.spikes{1}-stimon_t;
        a.spiketimes = spiketimes(spiketimes > 0);
        gotfulltrial = 1;
    catch
        keyboard
        gotfulltrial = 0;
    end
    set(gcf,'UserData',a);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extracting the parameters that define the current trial,
% setting up a new row in the 'storage.trialidentifiers' field
% if necessary, and updating the corresponding number correct
% and total number of trials.
function IntegrateTrial()
    a = get(gcf,'UserData');
    foundit = 0;
    
    % Extracting the parameters that define this trial
    currenttrialident = nan*ones(1,length(a.trialidentstrings));
    for i = 1:length(a.trialidentstrings)
        currenttrialident(i) = eval(['a.',char(a.trialidentstrings{i})]);
    end
  
    if (isempty(a.storage.trialidentifiers))
        disp('First trial');
        a.storage.trialidentifiers = currenttrialident;
        a.storage.correct = a.correct;
        a.storage.spikes = {[a.spiketimes; nan]};
        a.storage.ntrials = 1;
    else
        for i = 1:size(a.storage.trialidentifiers,1)
            if (all(currenttrialident == a.storage.trialidentifiers(i,:)))
               % disp('Found trial type');
                a.storage.correct(i) = a.storage.correct(i)+a.correct;
                a.storage.ntrials(i) = a.storage.ntrials(i)+1;
                a.storage.spikes(i) = {[a.storage.spikes{i}; a.spiketimes; nan]};
                foundit = 1;
                break;
            end
        end
        if (~foundit)
          %  disp('New trial type');
            a.storage.trialidentifiers(i+1,:) = currenttrialident;
            a.storage.correct(i+1) = a.correct;
            a.storage.ntrials(i+1) = 1;
            a.storage.spikes(i+1) = {[a.spiketimes; nan]};
        end
    end
    set(gcf,'UserData',a);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UDP communication functions below...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stopnow = DealWithMessage(msgSize)
    global udpCom;
    stopnow = 0;
    
    message = pnet(udpCom.sock, 'read', msgSize, 'char');
    if (strncmp(message,'return',6))
        a = dbstack;  % Check whether called from another function or from command line 
        if (~strcmp(a(end).name, mfilename))
            stopnow = 1;
        end
    end
    try
        eval(message);
    catch
        fprintf('Unknown message: %s\n', message);
        error = lasterror;
        disp(error.message);
        disp(error.identifier);
        disp(error.stack); 
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting functions below...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SetUpFigure()
    figure(1);
    a = get(gcf,'UserData');  % In case there's already things there (eg from whitenoise)
    set(gcf,'position',[300 200 750 500])
    set(gcf,'DefaultAxesUnits','pixels')
    clf;
    a.axeshandles.behavior = axes('position',[50 50 300 400]);
    a.axeshandles.neuron = axes('position',[400 50 300 400]);
    a.uicontrols.starttime = uicontrol('style','Edit','string',50,'Position',[400 450 40 20]);
    uicontrol('style','text','string','start','Position',[400 470 40 15]);
    a.uicontrols.stoptime = uicontrol('style','Edit','string',100,'Position',[460 450 40 20]);
    uicontrol('style','text','string','stop','Position',[460 470 40 15]);
    a.uicontrols.reset = uicontrol('style','pushbutton','Callback',@ResetCallback, 'string','RESET','Position',[520 450 60 20]);
    set(gcf,'UserData',a);
    drawnow;
end

function PlotFig()
    a = get(gcf,'UserData');
    
    colors = [.8 .8 0; .5 .5 .5; 0 0 1];
    linestyles = {'-',':'};
   
    nonmatchlidx = find(strcmp('nonmatchl',a.trialidentstrings));
    nonmatchmidx = find(strcmp('nonmatchm',a.trialidentstrings));
    nonmatchsidx = find(strcmp('nonmatchs',a.trialidentstrings));
    gradlidx = find(strcmp('gradl',a.trialidentstrings));
    gradmidx = find(strcmp('gradm',a.trialidentstrings));
    gradsidx = find(strcmp('grads',a.trialidentstrings));
    mstimidx = find(strcmp('mstim',a.trialidentstrings));

    start_t = str2double(get(a.uicontrols.starttime,'String'));
    stop_t = str2double(get(a.uicontrols.stoptime,'String'));
    
    ident = a.storage.trialidentifiers;
    ccs = ident(:, [nonmatchlidx nonmatchmidx nonmatchsidx]);
    [u,s,v] = svd(ccs);
    colordir = v(:,1);
    if (colordir(3) < 0)
        colordir = -colordir;
    end
    if (colordir(1) < 0)
        colordir = -colordir;
    end
    
    grads = ident(:,[gradlidx gradmidx gradsidx]);
    uniquegrads = unique(grads,'rows');
    
    [u,s,v] = svd(grads);
    gradcolordir = v(:,1);
    if (gradcolordir(3) < 0)
        gradcolordir = -gradcolordir;
    end
    if (gradcolordir(1) < 0)
        gradcolordir = -gradcolordir;
    end
    gradcolors = sort(grads*gradcolordir);
    
    ntrialtypes = size(ident,1);
    if (ntrialtypes > 1)
        % First extracting spike stats
        mn = nan*ones(ntrialtypes,1);
        sem = nan*ones(ntrialtypes,1);
        for i = 1:ntrialtypes
            tmp = a.storage.spikes{i};
            nanidxs = find(isnan(tmp));
            trialvect = zeros(size(tmp));
            trialvect(nanidxs) = 1;
            trialvect = cumsum(trialvect);
            L = tmp < start_t | tmp > stop_t | isnan(tmp);
            tmp(L) = [];
            trialvect(L) = [];
            mn(i) = length(tmp)/a.storage.ntrials(i)/(stop_t-start_t)*1000;
            if (~isempty(trialvect))
                counts = histc(trialvect,[0:1:trialvect(end)]);
                sem(i) = sqrt((var(counts))/a.storage.ntrials(i))/(stop_t-start_t)*1000;
            else
                sem(i) = nan;
            end
        end
        
        axes(a.axeshandles.behavior);
        cla; hold on;
        for i = 0:1 % mstim
            for j = 1:size(uniquegrads,1)
                Lgrad = all(grads == repmat(uniquegrads(j,:),size(grads,1),1),2);
                L = ident(:,mstimidx) == i & Lgrad;
                if (any(L))
                    [nmcolors, idxs] = sort(ccs(L,:)*colordir);
                    correct = a.storage.correct(L);
                    correct = correct(idxs)';
                    ntrials = a.storage.ntrials(L);
                    ntrials = ntrials(idxs)';  % Ugly
                    pcor = correct./ntrials;
                    h = errorbar(nmcolors, pcor, sqrt((pcor.*(1-pcor))./ntrials),'Linewidth',2);
                    set(h,'LineStyle',char(linestyles{i+1}));
                    tmp = uniquegrads(j,:)*10+.5;
                    tmp(tmp>1) = 1;
                    tmp(tmp<0) = 0;
                    set(h,'Color',tmp);
                end
            end
        end
        axes(a.axeshandles.neuron);  % Somewhat redundant with above code
        cla; hold on;
        for j = 1:size(uniquegrads,1)
        	L = all(grads == repmat(uniquegrads(j,:),size(grads,1),1),2);
            if (any(L))
                [nmcolors, idxs] = sort(ccs(L,:)*colordir);
                resp = mn(L);
                resp = resp(idxs)';
                semresp = sem(L);
                semresp = semresp(idxs)';
                h = errorbar(nmcolors, resp, semresp, 'Linewidth',2);
                tmp = uniquegrads(j,:)*10+.5;
                tmp(tmp>1) = 1;
                tmp(tmp<0) = 0;
                set(h,'Color',tmp);
            end
        end
        set(gca,'YLim',[0 1.1*max([mn; 1])]);
    end
end

%%%%%%%%%%%%%%%%%
% Callbacks
%%%%%%%%%%%%%%%%%

function ResetCallback(h, ev)
    a = get(gcf,'UserData');
    cla(a.axeshandles.behavior);
    cla(a.axeshandles.neuron);
    a.storage.trialidentifiers = [];
    set(gcf,'UserData',a);
end