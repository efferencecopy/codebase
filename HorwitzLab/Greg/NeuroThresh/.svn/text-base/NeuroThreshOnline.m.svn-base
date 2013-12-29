function NeuroThreshOnline()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Online analysis stuff for the NeuroThresh paradigm.
% Based on GratingOnline.m, this function will keep a 
% running tally of previous trials and spike rates and
% based on a set of rules will pick the stimulus for 
% the next trial
%
% Started 6/11/09 by GDLH
%
% Adding contrast scaled replay conditions
% GDLH 12/7/10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global Variables
global gl;

global udpCom;  % Needed by "DealWithMessage()" and associated subfunctions
    udpCom.sock = [];
    udpCom.port = 6665;
    udpCom.rexip = '192.168.1.120';

% Some initializations
s = InitPlex();     % Link to Plexon datastream
C = NeuroThreshCodes;        % Script that defines a bunch of constants
p = InitPStruct(0, C.FIXXCD); % Structure of the events from each trial
InitGLStruct();  % Global data structure
[sock, Success] = pnetStart(6665);  % UDP communication with REX
SetUpFigure;
PlotPSTH(nan, nan);
PickNewColorDirs(nan);  % Setting the persistent "levelcount" variable to 0.
spiketally = [];
stimon_t = [];
stimoff_t = [];
abortflag = 0;
plxrec = 0;

codestruct = {...
        C.FPONCD,@fponfn;...
        C.STIMONCD,@stimonfn;...
        C.ABORTCD,@abortfn;...
        C.STIMOFFCD,@stimofffn;...
        C.CORRECTCD,@correctfn;...
        C.EOTCD,@eotfn;...
        };

% The main loop
while (~gl.timetoleave)
    gl.timetoleave = CheckForESCKey();
    % Check for outgoing messages
    % Careful: codes don't get sent/dropped when PLXREC is off
    if (~isempty(gl.outmsg) && plxrec)
        eval(gl.outmsg{1});
        gl.outmsg(1) = [];
    end
    
    % Check for a message from REX
    msgSize = pnet(udpCom.sock, 'readpacket', 200, 'noblock');
    if(msgSize)
        DealWithMessage(msgSize);
    end
    [n, eventList] = PL_GetTS(s);
    if (n > 0)
       p = ProcessEventList(p, eventList);
    end
    % Getting the header
    if (~gl.gotheader)
        p = GetHeader(p);
    end
    % Accumulating spikes
    if (any(p.spikes{1}) && ~isempty(stimon_t) && ~abortflag)
        spiketally = [spiketally; p.spikes{1}];
        p.spikes{1} = []; 
    end
    % Dealing with each code, one by one
    for i = 1:size(codestruct,1)
        if (any(p.events == codestruct{i,1}))
            feval(codestruct{i,2})
        end
    end
    p = CleanUpEvents(p);
end % big while loop

    %%%
    % Nested functions.
    % Because these are nested, we can rely on inheriting all the relevant
    % variables.
    % fpon
    function fponfn
        p.lastprocessed_t = p.times(find(p.events == C.FPONCD,1));
    end
    % stimon
    function stimonfn
        stimon_t = p.times(find(p.events == C.STIMONCD,1));
        spiketally = [];
        abortflag = 0;
        plxrec = 1;
        ClearStatusText();
        p.lastprocessed_t = p.times(find(p.events == C.STIMONCD,1));
    end
    % abort
    function abortfn
        abortflag = 1;
        if (isempty(gl.trialspecs))
            return
        end
        if (length(gl.trialspecs(gl.trialspecidx).norms) ==...
            length(gl.trialspecs(gl.trialspecidx).spikerates)+1)
            gl.trialspecs(gl.trialspecidx).norms(end) = [];
        end
        ClearStatusText();
        spiketally = [];
        p.lastprocessed_t = p.times(find(p.events == C.ABORTCD,1,'last'));
    end
    % stimoff
    function stimofffn
        stimoff_t = p.times(find(p.events == C.STIMOFFCD,1));
        p.lastprocessed_t = stimoff_t;
    end
    % correct
    function correctfn
        correct_t = p.times(find(p.events == C.CORRECTCD,1));
        a = get(gcf,'UserData');
        thresh = get(a.uicontrols.threshslider,'value');
        % Skipping the first few msecs as determined by PSTHslider
        t_offset = get(a.uicontrols.PSTHslider,'value');
        dur = stimoff_t-stimon_t-t_offset;
        nspikes = sum(spiketally<stimoff_t & spiketally>stimon_t+t_offset);
        spikerates = [gl.trialspecs(gl.trialspecidx).spikerates nspikes/dur*1000];
        gl.trialspecs(gl.trialspecidx).spikerates = spikerates;
        stepsize = gl.trialspecs(gl.trialspecidx).stepsize;
        norms = gl.trialspecs(gl.trialspecidx).norms;
        scale = gl.initparams.stepsizescale;
        nrev = gl.initparams.nreversals;
        if (stepsize ~= 0)
            if (length(spikerates) >= 1)
                s2 = spikerates(end)-thresh;
                if (s2 > 0)
                    set(a.uicontrols.statustext,'String','Response','BackgroundColor',[1 .5 .5])
                end
            end
            if (length(spikerates) >= 2)
                s1 = spikerates(end-1)-thresh;
                if (sign(s1) ~= sign(s2))  % A reversal has just occurred
                    stepsize = stepsize*scale;
                    % eps in below line avoids a round off error
                    if (stepsize-eps <= gl.initparams.stepsize*(scale^nrev))
                        gl.trialspecs(gl.trialspecidx).done= 1;
                        gl.trialspecs(gl.trialspecidx).measuredthreshold = norms(end);
                        disp(['Done with idx ',num2str(gl.trialspecidx)]);
                    end
                    gl.trialspecs(gl.trialspecidx).stepsize = stepsize;
                end
            end
        else
            if (strcmp(gl.trialspecs(gl.trialspecidx).type,'stationarity check'))
                set(a.uicontrols.statustext,'String','Stat.Chk.')
                gl.trialspecs(gl.trialspecidx).done = 1;
            elseif (strcmp(gl.trialspecs(gl.trialspecidx).type,'replay'))
                set(a.uicontrols.statustext,'String','Replay')
                % Ugly hard coding below.  Make this user definable. GH
                if (length(gl.trialspecs(gl.trialspecidx).spikerates) == 5)
                    gl.trialspecs(gl.trialspecidx).done = 1;
                end
            end
        end
        p.lastprocessed_t = correct_t; 
    end
    %eot
    function eotfn
        eot_t = p.times(find(p.events == C.EOTCD,1));
        if (~abortflag && ~gl.bldone)
            PlotFRhist();
            PlotPSTH(spiketally, stimon_t);
        end
        plxrec = 0;
        p.lastprocessed_t = eot_t;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UDP communication functions below...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DealWithMessage(msgSize)
    global udpCom;
    
    message = pnet(udpCom.sock, 'read', msgSize, 'char');
    if (strncmp(message,'return',6))
        a = dbstack;  % Check whether called from another function or from command line 
        if (~strcmp(a(end).name, mfilename))
            gl.timetoleave = 1;
        end
    end
    try
        if (strcmp(message,'sendToRex(udpCom, gaborparams, ''double'', message);'))
           gaborparams = PickGaborParams;
        end
        if (strcmp(message,'sendToRex(udpCom, trialparams, ''double'', message);'))
            trialparams = GetTrialParams;
        end
        if (strcmp(message,'sendToRex(udpCom, monopolar, ''integer'', message);'))
            b = get(gcf,'UserData');
            monopolar = get(b.uicontrols.monopolarcheckbox,'Value');
            set(b.uicontrols.monopolarcheckbox,'Enable','off');
        end
        eval(message);
    catch
        fprintf('Unknown message: %s\n', message);
        error = lasterror;
        disp(error.message);
        disp(error.identifier);
        for i = 1:length(error.stack)
            disp(error.stack(i));
        end 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
function InitGLStruct()
    global gl;
    
    gl.initparams.LMSdir = [];   % Initial 3 color dirs (sent from REX)
    gl.initparams.preflms = [];  % color direction for baseline and stationarity check trials
    gl.initparams.preflmsnorm = []; % contrast for stationarity check trials
    gl.initparams.stepsize = [];
    gl.initparams.stepsizescale = [];
    gl.initparams.nreversals = [];
    gl.initparams.stimdur = 0;
    gl.initparams.nnewdirs = 0;
    gl.initparams.linpredtol = 0;
    gl.initparams.oogscale = 0.5;
    gl.disp.monSpd = [];
    gl.disp.fundamentals = [];
    gl.disp.bkgndrgb = [];
    gl.disp.msperframe = 0;
    gl.trialspecs = []; % A structure containing all the information about the trial types
    gl.trialspecidx = []; % index into gl.trialspecs indicating the current trial
    gl.gotheader = 0;
    gl.timetoleave = 0;  % Flag when to leave this .m file
    gl.bldone = 0;  % Finished baseline measurement
    gl.outmsg = {}; % Queue of messages to be sent to REX
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Getting the header data from the ecodes
% Destructively modifies the p structure to update p.lastprocessed_t.
function p = GetHeader(p)
    global gl;
    
    C = NeuroThreshCodes;

    % Creating a cell array that we iterate over to collect the
    % header parameters.
    headerstruct = {...
        [C.MONSPDCD C.MONSPDCD],'double',@monspdfn;...
        [C.FUNDAMENTALSCD C.FUNDAMENTALSCD],'double',@fundamentalsfn;...
        C.FRAMERATECD,'double',@frameratefn;...
        [C.INITLMS1CD C.INITLMS1CD],'float',@initlmsfn;...
        [C.INITLMS2CD C.INITLMS2CD],'float',@initlmsfn;...
        [C.INITLMS3CD C.INITLMS3CD],'float',@initlmsfn;...
        C.INITSTEPSIZECD,'float',@stepsizefn;...
        C.STEPSIZESCALECD,'float',@stepsizescalefn;...
        C.NREVERSALSCD,'int',@nreversalsfn;...
        C.NNEWDIRSCD,'int',@nnewdirsfn;
        C.LINPREDTOLCD,'float',@linpredtolfn;
        C.OOGSCALECD,'float',@oogscalefn;
        [C.BKGNDRGBCD C.BKGNDRGBCD],'double',@bkgndrgbfn;
        [C.PREFLMSCD C.PREFLMSCD],'float',@preflmsfn;
        C.NFRAMESPLATCD,'int',@nframesplatfn;
        C.NFRAMESRAMPCD,'int',@nframesrampfn;
        C.MONOPOLARCD,'int',@monopolarfn;
        };

    for i = 1:size(headerstruct,1)
        [out, p] = GetValsIfPossible(headerstruct{i,1}, p, headerstruct{i,2});
        if (~isempty(out))
            feval(headerstruct{i,3},out)
        end
    end
    % Checking to see if we have the entire header
    % bkgndrgb is assumed to be the last thing dropped.
    if ~isempty(gl.disp.bkgndrgb)
        gl.gotheader = 1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Nested functions for processing header params
    function monspdfn(out)
        disp('Got monSPD');
        gl.disp.monSpd = reshape(out, length(out)/3,3);
    end
    function fundamentalsfn(out)
        disp('Got fundamentals');
        gl.disp.fundamentals = reshape(out, length(out)/3,3);
        wldelta = (780-380+2)/size(gl.disp.monSpd,1);
        P_device = SplineSpd([380:wldelta:780]',gl.disp.monSpd,[380:5:780]');
        gl.disp.M = gl.disp.fundamentals'*P_device;
        gl.disp.invM = inv(gl.disp.M);
    end
    function preflmsfn(out)
        gl.initparams.preflms = out'./norm(out);
    end
    function frameratefn(out)
        gl.disp.msperframe = 1000/out;
    end
    function initlmsfn(out)
        if (isempty(gl.initparams.LMSdir))
            gl.initparams.LMSdir = out';
        else
            L = all(gl.initparams.LMSdir == repmat(out',size(gl.initparams.LMSdir,1),1),2);
            if (~all(out == 0) && ~any(L))
                gl.initparams.LMSdir = [gl.initparams.LMSdir; out'];
            end
        end
    end
    function stepsizefn(out)
        gl.initparams.stepsize = out;
    end
    function stepsizescalefn(out)
        gl.initparams.stepsizescale = out;
    end
    function nreversalsfn(out)
        gl.initparams.nreversals = out;
    end
    function bkgndrgbfn(out)
        gl.disp.bkgndrgb = out;
    end
    function nnewdirsfn(out)
        gl.initparams.nnewdirs = out;
    end
    function linpredtolfn(out)
        gl.initparams.linpredtol = out;
    end
    function oogscalefn(out)
        gl.initparams.oogscale = out;
    end
    function nframesplatfn(out)
        gl.initparams.stimdur = gl.initparams.stimdur+out;
    end
    function nframesrampfn(out)
        gl.initparams.stimdur = gl.initparams.stimdur+2*out;
    end
    function monopolarfn(out)
        gl.monopolar = out;
    end
end

%%%%%
% Preparing and sending trial parameters when requested to by REX
% Called by REX after rewoff (so we know the just-past trial was correct)
% and before droptrialparams (so that the numbers returned by this function
% make it into the data file).

function out = GetTrialParams()
    global gl;

    nlev = PickNewColorDirs('?');  % query for which level of iteration we're on - a persistent variable in PickNewColorDirs
    
    % Ugly hack below to avoid the problem that as soon as the user clicks the
    % BLdone button the level advances, even if we're showing a probe or
    % blanks trial (trial types 1 and 2 which occur only before bldone).
    if (gl.trialspecidx) <=2
        nlev = 0;
    end
    a = get(gcf,'UserData');
    thresh = get(a.uicontrols.threshslider,'value');
    spikerates = gl.trialspecs(gl.trialspecidx).spikerates;
    stepsize = gl.trialspecs(gl.trialspecidx).stepsize;
    s1 = 0; s2 = 0;
    if (stepsize ~= 0)
        if (length(spikerates) >= 1)
            s2 = spikerates(end)-thresh;
        end
        if (length(spikerates) >= 2)
            s1 = spikerates(end-1)-thresh;
        end
    end
    out = [gl.trialspecidx (sign(s2)-sign(s1))/2 gl.trialspecs(gl.trialspecidx).stepsize s2>0 nlev];
end

%%%%%%%%%%%%%%%%%%%%%%
% Adding a new trial type to the gl.trialspecs structure.
% Parentvertices: each *row* is an LMS triplet, not each column
% Parentoutofgamut: 1x3 boolean vector
% Type: 'standard','baseline','stationarity check','replay'
function SetUpTrialSpecs(colordir, predthresh, stepsize, type, parentvertices, parentoog, active)
    global gl;
    

    if (nargin < 7)
        active = 1;
    end
    if (nargin < 6)
        parentoog = [];
    end
    if (nargin < 5)
        parentvertices = [];
    end
    if (nargin < 4)
        type = 'standard';
    end
    if (nargin < 3)
        stepsize = gl.initparams.stepsize;
    end

    idx = length(gl.trialspecs)+1;
    gl.trialspecs(idx).colordir = colordir;
    gl.trialspecs(idx).predictedthreshold = predthresh;
    gl.trialspecs(idx).stepsize = stepsize;
    gl.trialspecs(idx).type = {type};
    gl.trialspecs(idx).norms = [];
    gl.trialspecs(idx).spikerates = [];
    gl.trialspecs(idx).done = 0;  % For selecting each trial 
    gl.trialspecs(idx).active = active;   % For selecting color directions 
    gl.trialspecs(idx).parentvertices = parentvertices;
    gl.trialspecs(idx).parentoutofgamut = parentoog;
    gl.trialspecs(idx).measuredthreshold = [];
    gl.trialspecs(idx).outofgamut = 0;

    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Randomly picks a valid color direction, finds the next contrast
% and returns the cone contrast triplet at the peak of the Gabor.
% Calls SetUpTrialSpecs if it hasn't yet been called.
% Setting stepsize to 0 keeps contrast fixed (for baseline measurements,
% stationarity checks).
function out = PickGaborParams()
    global gl
    
    a = get(gcf,'UserData');
    thresh = get(a.uicontrols.threshslider,'Value');
    gamutviolation = 0;
    
    if (~isfield(gl.trialspecs, 'done'))
        SetUpTrialSpecs([0 0 0], 0, 0,'baseline');     % By default start with zero contrast
        [~,scalar] = gamutCheck(gl.initparams.preflms, gl.disp.bkgndrgb, gl.disp.M, 'both');
        SetUpTrialSpecs(gl.initparams.preflms, get(a.uicontrols.prefcontslider,'Value')*scalar, 0, 'baseline');  % and a non-zero contrast
    end
    idxs = find([gl.trialspecs.active]&~[gl.trialspecs.done]);
    if (isempty(idxs))
        PickNewColorDirs();
        idxs = find([gl.trialspecs.active]&~[gl.trialspecs.done]);
    end
    gl.trialspecidx = idxs(unidrnd(length(idxs)));
    
    if (isempty(gl.trialspecs(gl.trialspecidx).norms))
        norm = gl.trialspecs(gl.trialspecidx).predictedthreshold;
    else
        if (gl.trialspecs(gl.trialspecidx).spikerates(end) > thresh)
            % Decrementing contrast
            norm = gl.trialspecs(gl.trialspecidx).norms(end).*(1-gl.trialspecs(gl.trialspecidx).stepsize);
        else  % Incrementing contrast
            norm = gl.trialspecs(gl.trialspecidx).norms(end).*(1+gl.trialspecs(gl.trialspecidx).stepsize);
            gamutviolation = GamutViolation(norm, gl.trialspecs(gl.trialspecidx).colordir);
            if (gamutviolation)
                if (gl.bldone == 1)
                    disp(['Done with idx ',num2str(gl.trialspecidx),' (out of gamut)']);
                    gl.trialspecs(gl.trialspecidx).done = 1;
                    gl.trialspecs(gl.trialspecidx).measuredthreshold =...
                        gl.trialspecs(gl.trialspecidx).norms(end);
                    gl.trialspecs(gl.trialspecidx).outofgamut = 1;
                    out = PickGaborParams;
                elseif (gl.bldone == 0)
                    disp('Error: Out of gamut but bl not done yet!');
                    keyboard
                end
            end
        end
    end

    if (~gamutviolation)
        out = norm*gl.trialspecs(gl.trialspecidx).colordir;
        % Because we're appending to gl.trialspecs().norms here we have 
        % to strip elements off gl.trialspecs().norms on aborted trials.
        gl.trialspecs(gl.trialspecidx).norms = ...
                [gl.trialspecs(gl.trialspecidx).norms, norm];
    end
    
    % Checking for gamut violations.
    function out = GamutViolation(norm, colordir)
        bkgndlms = gl.disp.M*gl.disp.bkgndrgb;
        rgbs = [gl.disp.invM*(bkgndlms.*(1+norm*colordir'));
                gl.disp.invM*(bkgndlms.*(1-norm*colordir'))];
        if (any(rgbs > 1) || any(rgbs < 0))
            out = 1;
            set(a.uicontrols.statustext,'String','Out of gamut','BackgroundColor',[.5 .5 1]);
        else
            out = 0;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Called once all the threshold have been determined for the current set of
% color directions.  Optional argument is for querying 'levelcount' 
% (or resetting it to 0) which is how many levels of iteration we are at.
% 0: baseline
% 1: initial color direction (from REX)
% 2: first 4 triangles
% 3+: additional triangles
%
% Function of the 'active' and 'done' fields in gl.trialspecs(i):
% Active  Done
% ------  ----
%   0      0  Color directions waiting the the queue to be probed
%   0      1  Completed color directions
%   1      0  Color directions that are currently being probed
%   1      1  Color directions that have been completes, but haven't yet
%   been broken up into subtriangles.
%   So the progression is: (0,0) -> (1,0) -> (1,1) -> (0,1)
%
% Here is the basic flow of events (for level 3 and above)
% 1) Active color directions are extracted (all should also be "done" or
% this function would not be called).
% 2) Color directions satisfing the local linearity tolerance are
% inactivated.
% 3) Color directions NOT satisfing the local linearity tolerance are
% used as the vertices for three new triangles (the centers of these new
% triangles are ~active and ~done).
% 4) All ~active and ~done color directions are pooled, and a subset are
% activated for probing.

function out = PickNewColorDirs(in)
    global gl
    persistent levelcount
    a = get(gcf,'UserData');
    out = [];
    
    if (nargin)
        if (strcmp(in,'?'))
            out = levelcount;
        elseif (isnan(in))
            levelcount = 0;
        end
        return
    end
    
    nlev = str2double(get(a.uicontrols.nlevels,'String'));
    if (levelcount >= nlev)  % A little ugly, this code is duplicated below.
        if (get(a.uicontrols.replaycheckbox,'Value'))
            SetUpReplay();
        else
            SetUpTrialSpecs([100 0 0],1, 0); % send stop signal (L contrast > 1)
        end
        disp('Incrementing level count for replays');
        levelcount = levelcount + 1;
        return;
    end
    Lact = logical([gl.trialspecs.active]);
    Lprobe = logical(strcmp([gl.trialspecs.type],'stationarity check'));
    % First time through
    if (~gl.bldone)  % in BlDoneCallback: PickNewColorDirs is called before gl.done is set to 1.
        for i = 1:size(gl.initparams.LMSdir)
            v = gl.initparams.LMSdir(i,:);
            if (get(a.uicontrols.polarityflip,'Value'))
                v = -v;
            end
            SetUpTrialSpecs(v/norm(v), norm(v), gl.initparams.stepsize);
            if (gl.monopolar)
                SetUpTrialSpecs(-v/norm(v), norm(v), gl.initparams.stepsize);
            end
        end
        % Stationarity check trials
        SetUpTrialSpecs(gl.initparams.preflms, gl.initparams.preflmsnorm, 0, 'stationarity check');
    % Second time through (first time using threshold triangles to pick
    % color directions).  Logic: if there are no parents, this must be the
    % first time through.  No triangles omitted from this round even if all
    % parents are OOG.
    elseif (all(isempty([gl.trialspecs(logical(Lact)).parentvertices])))
        [gl.trialspecs(logical(Lact)).active] = deal(0); % deactivating initial color directions.
        lmsmat = reshape([gl.trialspecs(Lact&~Lprobe).colordir],3,sum(Lact&~Lprobe))';
        parentoog = [gl.trialspecs(Lact&~Lprobe).outofgamut];
        threshs  = [gl.trialspecs(Lact&~Lprobe).measuredthreshold];
        % For the purposes of selecting vertices, decrease the threshold for 
        % color directions that went outside of the gamut.
        % This prevents points that are out of gamut from
        % "sucking" the next point out of the gamut too.  Do this only when
        % calculating 'v', the new color direction, not when storing
        % the threshold point in parentvertices (otherwise we keep on
        % decrementing it!).  Line below, (v = mean...), maps parentoog = 0
        % to parentvertices and parentoog = 1 to parentvertices*oogscale.
       if (~gl.monopolar)  % Bipolar stimuli
            signmat = [(fullfact([2,2])-1.5)*2 ones(4,1)];  % assumes three color dirs    
            % each *row* is an LMS triplet
            for i = 1:size(signmat,1)
                parentvertices = diag(signmat(i,:).*threshs)*lmsmat;
                v = mean(parentvertices.*repmat((gl.initparams.oogscale-1)*parentoog'+1,1,3));
                SetUpTrialSpecs(v/norm(v), norm(v), gl.initparams.stepsize, 'standard', parentvertices, parentoog);
            end
        else  % Monopolar stimuli
            signmat = (fullfact([2,2,2])-1.5)*2;
            selmat = false(size(signmat,1), size(signmat,2),2);
            selmat(:,:,1) = signmat == 1;
            selmat(:,:,2) = ~selmat(:,:,1);
            selmat = reshape(permute(selmat,[3 2 1]), [size(signmat,2)*2 size(signmat,1)])';
            % Assumes that complimentary color directions are in pairs
            % (e.g. A, A', B, B',...) 
            for i = 1:size(selmat,1)
                parentvertices = diag(selmat(i,:).*threshs)*lmsmat;
                parentvertices = parentvertices(selmat(i,:),:);
                poog = parentoog(selmat(i,:));
                v = mean(parentvertices.*repmat((gl.initparams.oogscale-1)*poog'+1,1,3));
                SetUpTrialSpecs(v/norm(v), norm(v), gl.initparams.stepsize, 'standard', parentvertices, poog);
            end        
        end
        
        % Stationarity check trial
        SetUpTrialSpecs(gl.initparams.preflms, gl.initparams.preflmsnorm, 0, 'stationarity check');
    else
        % Third time through or later
        % Inactivating active stationarity check trials
        [gl.trialspecs(Lact&Lprobe).active] = deal(0);
        % First, eliminating locally linear color directions
        for i = find(Lact&~Lprobe)
            threshratio = gl.trialspecs(i).measuredthreshold/gl.trialspecs(i).predictedthreshold;
            if (abs(log(threshratio)) < abs(log(1+gl.initparams.linpredtol)))
                disp([num2str(i),': consistent with linear model']);
                gl.trialspecs(i).active = 0;
                Lact(i) = 0;
            end
        end
        % Second, subdividing triangles in which the added interior point 
        % does not lie near the plane (is not locally linear).
        for i = find(Lact&~Lprobe)
            grandparentvertices = gl.trialspecs(i).parentvertices;
            grandparentoog = gl.trialspecs(i).parentoutofgamut;
            for j = 1:3
                parentvertices = grandparentvertices;
                parentoog = grandparentoog;
                parentvertices(j,:) = gl.trialspecs(i).colordir.*gl.trialspecs(i).measuredthreshold;
                parentoog(j) = gl.trialspecs(i).outofgamut;
                if (~all(parentoog))
                    v = mean(parentvertices.*repmat((gl.initparams.oogscale-1)*parentoog'+1,1,3));
                    SetUpTrialSpecs(v/norm(v), norm(v), gl.initparams.stepsize, 'standard', parentvertices, parentoog, 0);
                else
                    disp(['Skipping a triangle because of OOG: (', num2str(i),', ',num2str(j),')']);
                end
            end
            gl.trialspecs(i).active = 0;  % inactivating subdivided triangle
        end
        % Third, creating a list of ~active, ~done triangles that we will pick from
        % sanity check: should be zero
        if any([gl.trialspecs.active])
           disp('Error: some color directions are still active.');
           keyboard;
        end
        L = logical(~[gl.trialspecs.done]);
        if (sum(L) == 0)
            if (get(a.uicontrols.replaycheckbox,'Value'))
                SetUpReplay();
            else
                SetUpTrialSpecs([100 0 0],1, 0); % send stop signal (L contrast > 1)
            end
            disp('Incrementing level count for replays');
            levelcount = levelcount + 1;
            return;
        end
        idxs = find(L);
        disp('Color directions in queue');
        idxs
        if (gl.initparams.nnewdirs < length(idxs)) % More triangles in queue than we can show in one round
            idxs = idxs(floor([0:gl.initparams.nnewdirs-1]*(length(idxs)/gl.initparams.nnewdirs))+1);
            % Taking a broad assortment of triangles throughout the list
            % (oldest ones are at the top, newest ones are at the bottom).
        end
        disp('Color directions to be shown');
        idxs
        [gl.trialspecs(idxs).active] = deal(1);
        % Stationarity check trial
        SetUpTrialSpecs(gl.initparams.preflms, gl.initparams.preflmsnorm, 0 ,'stationarity check');
    end
    levelcount = levelcount + 1;
    set(a.uicontrols.currentleveltextbox,'String',num2str(levelcount));
end

% For the moment, hardcoding multipliers of 1 and 2
function SetUpReplay()
    global gl;
    a = get(gcf,'UserData');
    
    disp('doing replay experiment');
    set(a.uicontrols.replaycheckbox,'Value',0,'Enable','off');
    [gl.trialspecs.active] = deal(0);
    Loog = [gl.trialspecs.outofgamut];
    Luntested = ~[gl.trialspecs.done];
    Lmeasuredthresh = (~isempty([gl.trialspecs.measuredthreshold]));
    Lstandard = strcmp([gl.trialspecs.type],'standard');
    L = Lmeasuredthresh & ~Loog & Lstandard & ~Luntested;
    if any(L)
        for i = find(L)
            for multiplier = 1:2  % GH: Change this hardcoding
                cc = gl.trialspecs(i).colordir*gl.trialspecs(i).measuredthreshold*multiplier;
                if (gamutCheck(cc, gl.disp.bkgndrgb, gl.disp.M, 'both'))
                    SetUpTrialSpecs(gl.trialspecs(i).colordir, gl.trialspecs(i).measuredthreshold*multiplier,0, 'replay');
                end
            end
        end
    else
        SetUpTrialSpecs([100 0 0],1, 0);  % No directions to test
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting & analysis functions below...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SetUpFigure()

    figure(1);
    a = get(gcf,'UserData');  % In case there's already things there (eg from whitenoise)
    if (isfield(a,'uicontrols'))
        staticcontrols = {'threshslider','threshval', 0;...
            'PSTHslider','PSTHval',0;...
            'prefcontslider','contval', 0;...
            'monopolarcheckbox','monoval', 0;...
            'replaycheckbox','replay', 0};
        for i = 1:size(staticcontrols,1)
            if (isfield(a.uicontrols,staticcontrols{i,1}))
                if (ishandle(eval(['a.uicontrols.',staticcontrols{i,1}])))
                    eval([staticcontrols{i,2},' = get(a.uicontrols.',staticcontrols{i,1},',''Value'');']);
                else
                    eval([staticcontrols{i,2},' = ',num2str(staticcontrols{i,3}),';']);
                end
            else
                eval([staticcontrols{i,2},' = ',num2str(staticcontrols{i,3}),';']);
            end
        end
    else
        threshval = 0;
        contval = 0;
        PSTHval = 0;
        monoval = 0;
        replay = 0;
    end
        
    set(gcf,'position',[200 200 900 400])
    set(gcf,'DefaultAxesUnits','pixels')
    clf;
    a.axes.h1 = axes('position',[50 50 250 150]);
    a.axes.h2 = axes('position',[350 50 250 150]);
    a.uicontrols.threshslider = uicontrol('style','slider','Min',0,'Max',1.2*threshval+.05,'Callback',@ThreshSliderCallback,'SliderStep',[.01 .1],'Value',threshval,'Position',[35 210 280 15]);
    a.uicontrols.prefcontslider = uicontrol('style','slider','Min',0,'Max',.9,'Callback',@PrefContSliderCallback,'SliderStep',[.01 .1],'Value',contval,'Position',[35 260 280 15]);
    a.uicontrols.mkachrom = uicontrol('style','checkbox','Min',0,'Max',1,'CallBack',@MkAchromCallback,'Position',[200 360 15 15]);
    a.uicontrols.labels.mkchrom = uicontrol('style','text','string','Achrom','Position',[220 360 45 15]);
    a.uicontrols.polarityflip = uicontrol('style','checkbox','Min',0,'Max',1,'CallBack',@PolarityFlipCallback,'Position',[200 340 15 15]);
    a.uicontrols.labels.polarityflip = uicontrol('style','text','string','Polarity','Position',[220 340 45 15]);
    a.uicontrols.monopolarcheckbox = uicontrol('style','checkbox','Min',0,'Max',1,'Position',[200 320 15 15],'Value',monoval);
    a.uicontrols.labels.monopolarcheckbox = uicontrol('style','text','string','Monopolar','Position',[220 320 45 15]);
    a.uicontrols.replaycheckbox = uicontrol('style','checkbox','Min',0,'Max',1,'Position',[200 300 15 15],'Value',replay);
    a.uicontrols.labels.replaycheckbox = uicontrol('style','text','string','Replay','Position',[220 300 45 15]);
    a.uicontrols.bldone = uicontrol('style','pushbutton','string','Got Baseline','Callback',@BlDoneCallback,'Position',[130 230 80 25]);
    a.uicontrols.clearspikes = uicontrol('style','pushbutton','string','Reset','Callback',@ClearSpikesCallback,'Position',[50 230 80 25]);
    a.uicontrols.PSTHslider = uicontrol('style','slider','Min',0,'Max',500,'Callback',@PSTHSliderCallback,'SliderStep',[.01 .1],'Value',PSTHval,'Position',[335 210 280 15]);
    a.uicontrols.statustext = uicontrol('style','text','Position',[35 300 90 20]);
    a.uicontrols.nlevels = uicontrol('style','edit','Position',[130 300 30 20],'String','10');
    a.uicontrols.currentleveltextbox = uicontrol('style','text','Position',[130 330 30 20],'String',num2str(0));
    uicontrol('style','text','Position',[35 330 90 20],'String','Current level:');
    set(gcf,'UserData',a);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Callbacks for uicontrols
    % Callback executed when user manipulates threshold slider
    function ThreshSliderCallback(ev,h)
        a = get(gcf,'UserData');
        axes(a.axes.h1);
        b = get(a.axes.h1,'Children');
        L = strcmp(get(b,'tag'),'thresh');
        if (any(L))
            val = get(a.uicontrols.threshslider,'Value');
            set(b(L),'Xdata',[val val]);
        end
    end

    % Moves the latency bar in axes h2
    function PSTHSliderCallback(ev,h)
        a = get(gcf,'UserData');
        axes(a.axes.h2);
        b = get(a.axes.h2,'Children');
        L = strcmp(get(b,'tag'),'PSTH');
        if (any(L))
            val = get(a.uicontrols.PSTHslider,'Value');
            set(b(L),'Xdata',[val val]);
        end
    end

    % Reset plots
    function ClearSpikesCallback(ev,h)
       global gl
       
       [gl.trialspecs.spikerates] = deal([]);
       [gl.trialspecs.norms] = deal([]);
       set(a.uicontrols.threshslider,'Value',0);
       set(a.uicontrols.prefcontslider,'Value',.03);
       set(a.uicontrols.PSTHslider,'Value',.0);
       PlotPSTH(nan, nan);
    end

    % Once threshold and latency parameters have been set, user clicks
    % the bldone button whereupon new color directions are picked (given by
    % REX in header), threshold and latency values are packed into messages
    % for REX, and uicontrols are disabled.  Also the contrast of the probe
    % trials (adjusted with the slider) gets preserved in the global variable
    % gl.initparams.preflmsnorm for use in stationarity check trials.
    function BlDoneCallback(ev,h)
        global gl;
        
        [gl.trialspecs.done] = deal(1);
        [gl.trialspecs.active] = deal(0);
        
        n = length(gl.trialspecs);
        L = reshape([gl.trialspecs.colordir],3,n) == repmat([0 0 0]',1,n);
        idx = find(~all(L),1);
        if (~isempty(idx) && ~isempty(gl.trialspecs(idx).norms))
            [~,scalar] = gamutCheck(gl.trialspecs(idx).colordir, gl.disp.bkgndrgb, gl.disp.M, 'both');
            gl.initparams.preflmsnorm = get(a.uicontrols.prefcontslider,'Value')*scalar;
        end
        
        PickNewColorDirs;
        a = get(gcf,'UserData');
        axes(a.axes.h1);
        set(a.uicontrols.threshslider,'Enable','off');
        set(a.uicontrols.bldone,'Enable','off');
        set(a.uicontrols.prefcontslider,'Enable','off');
        set(a.uicontrols.PSTHslider,'Enable','off');
        set(a.uicontrols.clearspikes,'Enable','off');
        set(a.uicontrols.polarityflip,'Enable','off');
        set(a.uicontrols.mkachrom,'Enable','off');
        val = get(a.uicontrols.threshslider,'Value');
        gl.outmsg = [gl.outmsg; 'sendToRex(udpCom,',num2str(val),',''double'', ''Threshold'');'];
        val = get(a.uicontrols.PSTHslider,'Value');
        gl.outmsg = [gl.outmsg; 'sendToRex(udpCom,',num2str(val),',''double'', ''Latency'');'];
        gl.bldone = 1;
    end

    % Changes the contrast of the 'probe' stimulus used during the baseline
    % measurement.
    function PrefContSliderCallback(ev,h)
        global gl
        if (~isfield(gl.trialspecs,'colordir'))
            return
        end
        a = get(gcf,'UserData');
        val = get(a.uicontrols.prefcontslider,'Value');
        n = length(gl.trialspecs);
        L = reshape([gl.trialspecs.colordir],3,n) == repmat([0 0 0]',1,n);
        idx = find(~all(L),1);
        if (~isempty(idx) && ~isempty(gl.trialspecs(idx).norms))
            [~,scalar] = gamutCheck(gl.trialspecs(idx).colordir, gl.disp.bkgndrgb, gl.disp.M, 'both');
            gl.trialspecs(idx).norms(end) = val*scalar;
        end
    end
    
    % Toggle probe stimulus between the "preferred color" as determined
    % by the grating paradigm and achromatic.
    function MkAchromCallback(ev,h)
        global gl
        
        if (~isfield(gl.trialspecs,'colordir'))
            return
        end
        n = length(gl.trialspecs);
        L = reshape([gl.trialspecs.colordir],3,n) == repmat([0 0 0]',1,n);
        idx = find(~all(L),1);
        if (~isempty(idx) && ~isempty(gl.trialspecs(idx).norms))
            a = get(gcf,'UserData');
            if (get(a.uicontrols.mkachrom,'Value'))
                gl.trialspecs(idx).colordir = [1 1 1]./sqrt(3);
            else
                gl.trialspecs(idx).colordir = gl.initparams.preflms;
            end
            PrefContSliderCallback;
        end
    end

    %  Flip the phase of the stimuli shown during the baseline measurement.
    %  This will hopefully permit a better estimate of the response latency for arbitrary phase.
    %  User now controls polarity (in monopolar case) or 0/180 deg phase (in
    %  bipolar case)with a radio button.
    % Are we losing the identity of which trials were collected with which
    % polarity?
    function PolarityFlipCallback(~,h)
       global gl
        
        if (~isfield(gl.trialspecs,'colordir'))
            return
        end
        n = length(gl.trialspecs);
        L = reshape([gl.trialspecs.colordir],3,n) == repmat([0 0 0]',1,n);
        idx = find(~all(L),1);
        if (~isempty(idx) && ~isempty(gl.trialspecs(idx).norms))
            gl.trialspecs(idx).colordir = -gl.trialspecs(idx).colordir;
            PrefContSliderCallback
        end
    end
end

% Clear the status text box after each trial
function ClearStatusText()
    a = get(gcf,'UserData');
    set(a.uicontrols.statustext,'String',[],'BackgroundColor',[0.92549 0.913725 0.847059]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting functions
% Plotting a histogram of firing rates to allow the user to set a threshold
function PlotFRhist()
    global gl
    
    MINBINS = 5;
    a = get(gcf,'UserData');
    axes(a.axes.h1);
    maxrate = max([1.5*max([gl.trialspecs.spikerates]) MINBINS]);
    bins = linspace(0,maxrate,20);
    actidxs = find(~[gl.trialspecs.done]);
    n = zeros(length(gl.trialspecs), length(bins));
    for i = 1:length(actidxs)
        spikerates = gl.trialspecs(actidxs(i)).spikerates;
        [n(i,:),x] = hist(spikerates,bins);
    end
    bar(x,n','stacked');
    ylabel('count');
    xlabel('spike rate (sp/sec)');
    threshval = get(a.uicontrols.threshslider,'Value');
    set(gca,'Xlim',[-.5 max(maxrate, threshval)]);
    set(a.uicontrols.threshslider,'Min',0,'Max',max(maxrate, threshval));
    hold on;
    val = get(a.uicontrols.threshslider,'Value');
    hthreshline = plot([val val],[0 max(n(:))],'k-');
    set(hthreshline,'tag','thresh'); % So that ThreshSliderCallback knows which axis child to manipulate
    hold off;
end

% Plotting a PSTH to assist the user in selecting a latency
function PlotPSTH(spiketally, stimon_t)
    global gl
    persistent psth
    persistent trialcounter
    
    if (isnan(spiketally))
        psth = [];
        trialcounter = 0;
    end
    if (isempty(spiketally))
        return;
    end
    a = get(gcf,'UserData');
    axes(a.axes.h2);
    
    binwidth = 20;
    bins = [0:binwidth:gl.initparams.stimdur*gl.disp.msperframe];
    [n,x] = histc(spiketally-stimon_t, bins);
    if (length(spiketally) < 2)  % annoying matlab transpose issue
        n = n';
    end
    if (isempty(psth))
        psth = n;
        trialcounter = 1;
    else
        psth = psth+n;
        trialcounter = trialcounter+1;
    end
    plot(bins, (psth./trialcounter)*(1000/binwidth));
    set(gca,'XLim',[0 500]);
    hold on;
    val = get(a.uicontrols.PSTHslider,'Value');
    hpsthline = plot([val val],[0 max(psth./trialcounter)*(1000/binwidth)],'k-');
    set(hpsthline,'tag','PSTH'); % So that PSTHSliderCallback knows which axis child to manipulate
    hold off;
end