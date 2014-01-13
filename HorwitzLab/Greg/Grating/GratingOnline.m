function GratingOnline()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Online analysis stuff for the grating paradigm.
% Unlike WhiteNoiseOnline.m, this function keeps going 
% around the main loop as quickly as possible without
% waiting for any particular code to be dropped.  This
% will be important if we ever want to do anything in
% "real time" (ie within a trial).
%
% Started 3/3/08 by GDLH
%
% list of protocols:
% 0) Null protocol (all done)
% 1) Orientation (move to 2 or 5)
% 2) Spatial frequency (move to 1 or 5)
% 3) Phase tuning a la Williams and Shapley (move to 0)
% 4) Classical color tuning (move to 0)
% 5) Aperture size (move to 4)
% 6) Preferred grating a bunch of times (move to 0)
% 7) Contrast response functions with Gabor patch (a la DTspot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 8/10/12 Added a radio button to make color tuning (protocol 4)
% optional. Currently, protocols 3,6,and 7 are unused (and there's
% no way to get to them). Also increased the contrast of the achromatic 
% stimulus from .4 to .6. GDLH
%
% 7/4/13 GDLH Adding a new protocol (8)
% Joint orientation/spatial frequency open loop
% measurement (primarily for use with the Utah array). 
% Only accessible via the the "protocol" user menu from REX.

% Global Variables
global gl;

global udpCom;  % Needed by "DealWithMessage()" and associated subfunctions
    udpCom.sock = [];
    udpCom.port = 6665;
    udpCom.rexip = '192.168.1.120';

% Some initializations
SetUpFigure;
C = GratingCodes;        % Script that defines a bunch of constants
TS = TrialSpecCodes;    % Another script that defines a bunch of constants
p = InitPStruct(0, C.FIXXCD); % Structure of the events from each trial
s = InitPlex();     % Link to Plexon datastream
InitGLStruct();  % Structure that holds the statistics
[sock, Success] = pnetStart(6665);  % UDP communication with REX
currentspiketally = [];
stimon_t = [];
fix_t = [];

% The main loop
while (~gl.timetoleave)
    gl.timetoleave = CheckForESCKey();

    % Check for a message from REX
    msgSize = pnet(udpCom.sock, 'readpacket', 200, 'noblock');
    if(msgSize)
        DealWithMessage(msgSize);
    end
    [n, eventList] = PL_GetTS(s);
    if (n > 0)
       p = ProcessEventList(p, eventList);
    end
    if (~gl.gotheader)
        p = GetHeader(p);  % "SetUpTrials" called from here if possible
    end
    if (any(p.events == C.FPACQCD))
        fix_t = p.times(find(p.events == C.FPACQCD,1));
        abortflag = 0;
    end
    if (any(p.events == C.STIMONCD))
        stimon_t = p.times(find(p.events == C.STIMONCD,1));
        nspikes = sum(currentspiketally > 0); % Spike times are relative to fix_t
        gl.trialspecs{gl.paramsidx,TS.BASELINE} =...
            [gl.trialspecs{gl.paramsidx,TS.BASELINE}; nspikes/(stimon_t-fix_t)*1000];
        p.lastprocessed_t = p.times(find(p.events == C.STIMONCD,1));
    end
    if (any(p.spikes{1}) && ~isempty(fix_t) && ~abortflag)
        currentspiketally = [currentspiketally; p.spikes{1}-fix_t];
        PlotRaster(currentspiketally);
        p.spikes{1} = []; 
    end
    if (any(p.events == C.ABORTCD))
        abortflag = 1;
        currentspiketally = [];
        ClearRaster;
        p.lastprocessed_t = p.times(find(p.events == C.ABORTCD,1,'last'));
    end
    if (any(p.events == C.STIMOFFCD))
        stimoff_t = p.times(find(p.events == C.STIMOFFCD,1));
        stimdur = stimoff_t-stimon_t;
        nspikes = sum(currentspiketally > stimon_t-fix_t);
        gl.trialspecs{gl.paramsidx,TS.SPIKES} = [gl.trialspecs{gl.paramsidx,TS.SPIKES}; nspikes/stimdur*1000];
        currentspiketally = [];
        stimon_t = [];
        fix_t = [];
        if (gl.protocol == 4)
            PlotColorTuning;
        elseif (gl.protocol == 7)
            PlotContrastResponseFunctions;
        else
            PlotTuningCurve;
        end
        ClearRaster;
        UpdateTrialCountdown;
        p.lastprocessed_t = stimoff_t; 
    end
    p = CleanUpEvents(p);
end % big while loop
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UDP communication functions below...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DealWithMessage(msgSize)
    global gl;
    global udpCom;
    
    message = pnet(udpCom.sock, 'read', msgSize, 'char');
    if (strncmp(message,'return',6))
        a = dbstack;  % Check whether called from another function or from command line 
        if (~strcmp(a(end).name, mfilename))
            gl.timetoleave = 1;
        end
    end
    try
        if (strcmp(message,'sendToRex(udpCom, gl.gratingparams, ''double'', message);'))
            PickGratingParams;
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
    if (~isfield(gl,'preforient'))
        gl.preforient = nan;
    end
    if (~isfield(gl,'prefsf'))
        gl.prefsf = nan;
    end
    if (~isfield(gl,'prefdiam'))
        gl.prefdiam = nan;
    end
    if (~isfield(gl,'prefccdir'))
        gl.prefccdir = [0 0 0];
    end
    gl.splineargmax = [];
    gl.disp.gammaTable = [];
    gl.disp.monSpd = [];
    gl.disp.fundamentals = [];
    gl.disp.bkgndrgb = [];    
    gl.ntrlspercond = [];
    gl.protocol = [];
    gl.gratingparamcols = []; % which columns of gl.trialspecs contain the grating parameters
    gl.gratingparams = []; % For sending grtaing parameters to REX.  Updated every trial.
    gl.trialspecs = {}; % A structure containing all the grating parameters used in a given
                        % protocol.  Updated each time the protocol is
                        % changed.
    gl.paramsidx = []; % index into gl.trialspecs indicating the current trial
                       % (from which gl.gratingparams is extracted)
    
    gl.gotheader = 0;
    gl.timetoleave = 0;  % Flag when to leave this .m file
    gl.stimon_t = [];
    gl.prefisolumcolors = [0 0 0];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Getting the header data from the ecodes
% Destructively modifies the p structure to update p.lastprocessed_t.
% Calls SetUpTrials if the complete header has been obtained.  
function p = GetHeader(p)
    global gl;
    
    C = GratingCodes;
    
    [out, p] = GetValsIfPossible([C.GAMMATABLECD C.GAMMATABLECD], p);
    if ~isempty(out)
        disp('Got gammatable');
        gl.disp.gammaTable = reshape(out, length(out)/3,3);
        gl.disp.invgammaTable = InvertGamma(gl.disp.gammaTable,1);
    end

    [out, p] = GetValsIfPossible([C.MONSPDCD C.MONSPDCD], p);
    if ~isempty(out)
        disp('Got monSPD');
        gl.disp.monSpd = reshape(out, length(out)/3,3);
        
    end

    [out, p] = GetValsIfPossible([C.FUNDAMENTALSCD C.FUNDAMENTALSCD], p);
    if ~isempty(out)
        disp('Got fundamentals');
        gl.disp.fundamentals = reshape(out, length(out)/3,3);
        if (size(gl.disp.monSpd, 1) == 81)
            P_device = SplineSpd([380:4:780]',gl.disp.monSpd,[380:5:780]');
        else
            P_device = SplineSpd([380:2:780]',gl.disp.monSpd,[380:5:780]');
        end
        gl.disp.M = gl.disp.fundamentals'*P_device;
        gl.disp.invM = inv(gl.disp.M);
    end

    [out, p] = GetValsIfPossible(C.FRAMERATECD, p, 'double');
    if (~isempty(out))
        disp('Got framerate');
        gl.disp.msperframe = 1000/out;
    end
    
    [out, p] = GetValsIfPossible([C.BKGNDRGBCD C.BKGNDRGBCD], p, 'double');
    if ~isempty(out)
        disp('Got background rgb');
        gl.disp.bkgndrgb = out;
    end

    headerstruct = {...
	C.INITPROTOCOLCD, 'gl.protocol',1,'int';...
    C.NTRLSPERCONDCD, 'gl.ntrlspercond',1,'int';...
    C.INITDIAMCD,'gl.DEFAULTDIAM',100,'int';...
    C.INITSFCD,'gl.DEFAULTSF',100,'int';...
    C.INITTFCD,'gl.DEFAULTTF',100,'int';...
    C.INITORIENTCD,'gl.DEFAULTORIENT',180/pi,'int';...
    C.INITPHASECD,'gl.DEFAULTPHASE',100,'int';...
    C.INITLCONTCD,'gl.DEFAULTCONECONTRASTS(1)',100,'int';...
    C.INITMCONTCD,'gl.DEFAULTCONECONTRASTS(2)',100,'int';...
    C.INITSCONTCD,'gl.DEFAULTCONECONTRASTS(3)',100,'int';...
    C.INITNCYCLESCD,'gl.DEFAULTNCYCLES',1,'int';...
    C.NCYCLESPLATCD,'gl.NCYCLESPLATCD',1,'float';...
    C.NCYCLESRAMPCD,'gl.NCYCLESRAMPCD',1,'float';...
    C.NCONTRASTSCD,'gl.NCONTRASTS',1,'int';...
    C.NGABORREPSCD,'gl.NGABORREPS',1,'int';...
    C.MAXCONTSCALEFACTCD,'gl.MAXCONTSCALEFACTOR',1,'float';...
    C.MINCONTSCALEFACTCD,'gl.MINCONTSCALEFACTOR',1,'float'};

    for i = 1:size(headerstruct,1)
        [out, p] = GetValsIfPossible(headerstruct{i,1}, p, headerstruct{i,4});
        if ~isempty(out)
            eval([headerstruct{i,2},' = out/',num2str(headerstruct{i,3}),';']);
        end
    end

    % Checking to see if we have the entire header
    if ~isempty(gl.disp.bkgndrgb) && ~isempty(gl.disp.monSpd) && ~isempty(gl.MINCONTSCALEFACTOR) && ~isempty(gl.disp.fundamentals)
        gl.gotheader = 1;    % For now assuming if we got this, we got everything
        SetUpTrials(gl.protocol);
    end
end

%%%%%%%%%%%%%%%%%%%%%%
% Setting up the gl structure according to the protocol specified by
% by REX. 
%
% gl.ntrlspercond = 0 means: "do an ANOVA to determine whether we should do another
% block of trials or (if p < 0.05) move to the next protocol".  If gl.ntrlspercond == 0 
% then set gl.trialspecs{:,1} = 1 so that we do one trial in each condition
% before doing the ANOVA again.
function SetUpTrials(protocol)
    global gl;
    
    a = get(gcf,'UserData');
    TS = TrialSpecCodes;
    % Defining "constants"
    gl.gratingparamcols = find(~ismember(fieldnames(TS),{'TRIALCOUNT','SPIKES','BASELINE'}));  % Columns of gl.trialspecs containing the grating parameters
    gl.trialspecs = {};
    gl.protocol = protocol;
    
    switch (protocol)
        case 1,
            orients = [0:pi/4:7*pi/4];
            ntrialtypes = length(orients);
        case 2,
%            sfs = [0.5 1 2 4];
            sfs = [0.5000    0.8409    1.4142    2.3784    4.0000];
            ntrialtypes = length(sfs);
        case 3,
            phases = [0:pi/4:7*pi/4];
            ntrialtypes = length(phases);
        case 4,
            colordirections = [.09 .09 0;
                .09 -.09 0;
                0 0 0.6364;
                .1273 0 0;
                0 .1273 0;
                0.0636   -0.0636    -0.4500;
                0.0636   -0.0636    0.4500;
                0.0636   0.0636    0.4500;
                0.0636  0.0636    -0.4500;
                ];
            colordirections = colordirections*.9; % Ugly hack to keep things in gamut on rig 1
            ntrialtypes = size(colordirections,1);
        case 5,
            apsizes = linspace(0.2,2,5);
            ntrialtypes = length(apsizes);
        case 6,
            ntrialtypes = 1;
        case 7,
            load sCSF_Fits_Kali.mat
            LvMidx = find(ismember(sign(fp.colors), [1 -1 0],'rows'));
            Sidx = find(ismember(sign(fp.colors), [0 0 1],'rows'));
            if (isnan(gl.prefsf))
                LvMsensitivity = polyval(fp.polyCoeff(LvMidx,:), gl.DEFAULTSF);
                Ssensitivity = polyval(fp.polyCoeff(Sidx,:), gl.DEFAULTSF);
            else
                LvMsensitivity = polyval(fp.polyCoeff(LvMidx,:), gl.prefsf);
                Ssensitivity = polyval(fp.polyCoeff(Sidx,:), gl.prefsf);
            end
            LM = gl.MAXCONTSCALEFACTOR .* (1./LvMsensitivity) .* fp.colors(LvMidx, :);
            [h,scale] = gamutCheck(LM, gl.disp.bkgndrgb, gl.disp.M, 'both');
            if (h)
                LM = LM.*scale;
            end
            S = gl.MAXCONTSCALEFACTOR .* (1./Ssensitivity) .* fp.colors(Sidx, :);
            [h,scale] = gamutCheck(S, gl.disp.bkgndrgb, gl.disp.M, 'both');
            if (h)
                S = S.*scale;
            end
            scaleVector = logspace(log10(gl.MINCONTSCALEFACTOR), 0, gl.NCONTRASTS);
            colordirections = [0 0 0; scaleVector'*LM; scaleVector'*S];
            ntrialtypes = size(colordirections,1);
        case 8,
            orients = [0:pi/4:7*pi/4];
            sfs = [0.5000 1 2];
            ntrialtypes = length(orients)*length(sfs);
        case 2,
%            sfs = [0.5 1 2 4];
            ntrialtypes = length(sfs);
        otherwise
            ntrialtypes = 0;
    end

    % Setting up default values
    for i = 1:ntrialtypes
        gl.trialspecs{i,TS.TRIALCOUNT} = max([gl.ntrlspercond 1]);  % No fewer than 1 trial/cond
        gl.trialspecs{i,TS.SPIKES} = [];  % spike counts
        if (~isnan(gl.prefdiam))
            gl.trialspecs{i,TS.DIAM} = gl.prefdiam;    % diam (deg)
        else
            gl.trialspecs{i,TS.DIAM} = gl.DEFAULTDIAM;    % diam (deg)
        end
        if (~isnan(gl.prefsf))
            gl.trialspecs{i,TS.SF} = gl.prefsf;    % diam (deg)
        else
            gl.trialspecs{i,TS.SF} = gl.DEFAULTSF;     % SF (cycles/deg)
        end
        gl.trialspecs{i,TS.TF} = gl.DEFAULTTF;     % TF (cycles/sec)
        if ~isempty(gl.prefccdir) & any(gl.prefccdir)
            gl.trialspecs{i,TS.CC} = gl.prefccdir;  % cone contrasts            
        else
            gl.trialspecs{i,TS.CC} = gl.DEFAULTCONECONTRASTS;  % cone contrasts
        end
        if (~isnan(gl.preforient))
            gl.trialspecs{i,TS.ORIENT} = gl.preforient;  % orientation (rad)
        else
            gl.trialspecs{i,TS.ORIENT} = gl.DEFAULTORIENT;  % orientation (rad)
        end
        gl.trialspecs{i,TS.PHI} = gl.DEFAULTPHASE;  % phase (rad)
        gl.trialspecs{i,TS.NCYCLES} = gl.DEFAULTNCYCLES;  % number of cycles
        gl.trialspecs{i,TS.STIMTYPE} = 0;  % 0 = grating, 1 = gabor
        gl.trialspecs{i,TS.PROTOCOL} = protocol;  % which protocol
    end
    
    if (protocol == 1)  % 8-point orientation tuning curve
        disp('Setting up trial parameters for orientation tuning');
        for i = 1:ntrialtypes
            gl.trialspecs{i,TS.ORIENT} = orients(i);  % orientation (rad)
        end
    elseif (protocol == 2)  % 5-point spatial frequency tuning curve
        disp('Setting up trial parameters for SF tuning');
        for i = 1:ntrialtypes
            gl.trialspecs{i,TS.SF} = sfs(i);     % SF (cycles/deg)
        end
    elseif (protocol == 3)  % 8-phase tuning curve
        disp('Setting up trial parameters for phase tuning');
        for i = 1:ntrialtypes
            gl.trialspecs{i,TS.PHI} = phases(i); % phase (rad)
            gl.trialspecs{i,TS.TF} = 0;  % TF has to equal 0, that's how REX knows to flash the stimulus
            gl.trialspecs{i,TS.NCYCLES} = 0;  % number of cycles (irrelevant)
        end
    elseif (protocol == 4)  % Color directions
        disp('Setting up trial parameters for color tuning');
        for i = 1:ntrialtypes
            gl.trialspecs{i,TS.CC} = colordirections(i,:);  % cone contrasts
            gl.trialspecs{i,TS.TRIALCOUNT} = gl.trialspecs{i,TS.TRIALCOUNT};
            % No longer doubling the number of color protocol trials
        end
    elseif (protocol == 5)  % Aperture sizes
        disp('Setting up trial parameters for aperture size tuning');
        for i = 1:ntrialtypes
            gl.trialspecs{i,TS.DIAM} = apsizes(i);  % aperture diameters (deg)
        end
    elseif (protocol == 7)  % Charlie's gabor contrast-response functions
        disp('Setting up trial parameters for contrast-response functions');
        for i = 1:ntrialtypes
            gl.trialspecs{i,TS.CC} = colordirections(i,:);  % cone contrasts
            gl.trialspecs{i,TS.STIMTYPE} = 1;
            gl.trialspecs{i,TS.TRIALCOUNT} = gl.NGABORREPS;
            gl.trialspecs{i,TS.NCYCLES} = round(gl.NCYCLESPLATCD+2*gl.NCYCLESRAMPCD);
        end
        set(a.uicontrols.reset,'Enable','off');        
    elseif (protocol == 8)  % 8-point orientation, 3-point SF joint tuning curve
        disp('Setting up trial parameters for brute force joint orientation/SF tuning');
        i = 1;
        for j = 1:length(orients)
            for k = 1:length(sfs)
                gl.trialspecs{i,TS.ORIENT} = orients(j);  % orientation (rad)
                gl.trialspecs{i,TS.SF} = sfs(k);     % SF (cycles/deg)
                i=i+1;
            end
        end
    elseif (protocol == 0)
        disp('Early termination request received');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Called when we've completed a protocol, this function figures out which
% protocol to run next.  
function protocol = PickNewProtocol
    global gl;
    uicontrols = getfield(get(gcf,'UserData'),'uicontrols');

    
    if (gl.protocol == 1) % orientation
        preforient = gl.splineargmax;
        if (isnan(gl.preforient))
            gl.preforient = preforient;
            lastpreforient = gl.DEFAULTORIENT;
        else
            lastpreforient = gl.preforient;
        end
        disp(['new preferred orientation is ',num2str(preforient*180/pi)]);
        gl.preforient = preforient;
        if (isnan(gl.prefsf) || mod(abs(preforient-lastpreforient),2*pi) > 20*pi/180) % 20 deg criterion
            protocol = 2;
            disp(['old preferred orientation was ',num2str(lastpreforient*180/pi)]);
        else
            protocol = 5;
        end
    elseif (gl.protocol == 2)  % spatial frequency
        prefsf = gl.splineargmax;
        if (isnan(gl.prefsf))
            gl.prefsf = prefsf;
            lastprefsf = gl.DEFAULTSF;
        else
            lastprefsf = gl.prefsf;
        end
        disp(['new preferred sf is ',num2str(prefsf)]);
        gl.prefsf = prefsf;
        if (isnan(gl.preforient) || (abs(prefsf-lastprefsf) > 1))
            protocol = 1;
            disp(['old preferred sf was ',num2str(lastprefsf)]);
        else
            protocol = 5;
        end
    elseif (gl.protocol == 3)  % Flashed random spatial phases
        protocol = 0;  % Done
    elseif (gl.protocol == 4)  % Colors
        [gl.prefisolumcolors, gl.prefccdir] = computePrefColors();
        protocol = 0;
    elseif (gl.protocol == 5)  % aperture sizes
        gl.prefdiam = gl.splineargmax;
        disp(['new preferred aperture size is ',num2str(gl.prefdiam)]);
        if (get(uicontrols.colorradio,'Value'))
            protocol = 4;
        else
            protocol = 0;
        end
    elseif (gl.protocol == 6)  % Preferred grating many times to get F1/F0
        protocol = 0;   
    elseif (gl.protocol == 7)  % Contrast response with gabors
        protocol = 0;
    elseif (gl.protocol == 8)  % Joint orientation/SF
        protocol = 0;
    elseif (gl.protocol == 0)
%         disp('Gratings is Finished:')
%         disp('Preferred Orientation =');
%         gl.preforient
%         disp('Preferred Spatial Frequency =')
%         gl.prefsf
%         disp('Preferred Aperture Size =')
%         gl.prefdiam
        protocol = 0;
    end
    disp(['Switching to protocol ',num2str(protocol)]);
    gl.splineargmax = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading the appropriate parameters into gl.gratingparams
% Order of grating parameters:
% [diam, sf, tf, lcont, mcont, scont, orient, phase, ncycles, stimtype, protocol].
% See bottom of code.
%
% This function gets called whenever REX sends a message
% requesting gl.gratingparams.
function PickGratingParams()
    global gl;
    global udpCom;
    TS = TrialSpecCodes;
    
    if (~gl.gotheader)
        keyboard
        error('We have not yet received the header yet!');
    end

    if (all([gl.trialspecs{:,TS.TRIALCOUNT}] == 0) && gl.ntrlspercond == 0)
        % If we're out of trials and REX has specified that we should do an ANOVA
        % to figure out whether it's time to advance to the next protocol.
        p = anova1([gl.trialspecs{:,TS.SPIKES}],[],'off')
        if (isnan(p) || p > 0.05) % Do another block 
            [gl.trialspecs{:,TS.TRIALCOUNT}] = deal(1);
        end
    end
    if (isempty(gl.trialspecs) || all([gl.trialspecs{:,TS.TRIALCOUNT}] == 0))
        % Out of trials either because the ANOVA was significant or because
        % we're using a fixed number of trials (as opposed to ANOVA).
        protocol = PickNewProtocol;
        if (protocol == 0)   % we're done, send preferred parameters
            gl.gratingparams = [gl.prefdiam; gl.prefsf; 0; 0; 0; 0; gl.preforient; 0; 0; 0];
            gl.gratingparams(isnan(gl.gratingparams)) = 0;  % Can't send nans for some reason, so sending 0 instead
            return;
        else
            SetUpTrials(protocol);
        end
    end

    % Picking a randomly chosen trial
    ntrialsleft = [gl.trialspecs{:,TS.TRIALCOUNT}]';
    choices = find(ntrialsleft == max(ntrialsleft));
    gl.paramsidx = choices(unidrnd(length(choices)));
    gl.gratingparams = [gl.trialspecs{gl.paramsidx, gl.gratingparamcols}];
    gl.trialspecs{gl.paramsidx,TS.TRIALCOUNT} = gl.trialspecs{gl.paramsidx,TS.TRIALCOUNT}-1;
    
    %alert Rex that the next trial's params have been chosen
    pnet(udpCom.sock, 'write', 'paramSetupComplete>> >>');
    pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting & analysis functions below...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SetUpFigure()
    figure(1);
    a = get(gcf,'UserData');  % In case there's already things there (eg from whitenoise)
    set(gcf,'position',[300 200 600 400])
    set(gcf,'DefaultAxesUnits','pixels')
    clf;
        
    h1 = axes('position',[50 50 300 200]);
    ylabel('sp/sec');
    
    h2 = axes('position',[50 300 300 50]);
    set(h2,'Ylim',[-1 1],'Xlim',[0 1000],'NextPlot','add');
    xlabel('Time (ms)');
    
    for i = 1:3
        h3(i) = axes('position',[50+(i-1)*120 100 80 80],'Visible','off',...
            'XTick',[],'YTick',[],'XColor',[1 1 1],'YColor',[1 1 1]);
    end
    
    % UICONTROLS
    uicontrol('style','text','Position',[380 40 80 20],'String','Trial countdown:');
    a.uicontrols.ntrialstextbox = uicontrol('style','text','Position',[460 40 30 20],'String',num2str(0));
    uicontrol('style','text','Position',[380 65 80 20],'String','Protocol number:');
    a.uicontrols.protocoltextbox = uicontrol('style','text','Position',[460 65 30 20],'String',num2str(0));
  %  uicontrol('style','text','Position',[500 65 80 20],'String','Include Color?');
  %  a.uicontrols.gaborradio = uicontrol('style','radio','Position',[535 50 10 10],'Value',0);
    uicontrol('style','text','Position',[500 65 80 20],'String','Include Color?');
    a.uicontrols.colorradio = uicontrol('style','radio','Position',[535 50 10 10],'Value',1);
    a.uicontrols.changecol = uicontrol('style','popupmenu','string','Achrom|LvsM|S|L-M+S|L-M-S','Callback',@ChangeCol,'Position',[370 330 80 20]);
    a.uicontrols.sftextbox = uicontrol('style','text','Position',[460 350 80 12],'String','SF');
    a.uicontrols.changesf = uicontrol('style','slider','Min',0.5,'Max',4,'Callback',@ChangeSF,'SliderStep',[.01 .1],'Value',2,'Position',[460 330 80 20]);
    a.uicontrols.diamtextbox = uicontrol('style','text','Position',[460 300 80 12],'String','Diam');
    a.uicontrols.changediam = uicontrol('style','slider','Min',0.2,'Max',2,'Callback',@ChangeDiam,'SliderStep',[.01 .1],'Value',1,'Position',[460 280 80 20]);
    a.uicontrols.reset = uicontrol('style','pushbutton','String','Reset','Position',[370 350 80 20],'Callback',@ResetCallback);
    
    a.axeshandles.tuningcurve = h1;
    a.axeshandles.raster = h2;
    a.axeshandles.colorspace = h3;
    set(gcf,'UserData',a);
    
    function ChangeCol(ev,h)
        global gl
        TS = TrialSpecCodes;
        
        uicontrols = getfield(get(gcf,'UserData'),'uicontrols');
        whichcolidx = get(uicontrols.changecol,'Value');
        whichcolnames = get(uicontrols.changecol,'String');
        whichcolname = whichcolnames(whichcolidx,:);
        switch deblank(whichcolname)
            case 'Achrom'
                cc = [.6 .6 .6];
            case 'LvsM'
                cc = [.09 -.09 0];
            case 'S'
                cc = [0 0 .6364];
            case 'L-M-S'
                cc = [0.0636 -0.0636 -0.45];
            case 'L-M+S'
                cc = [0.0636 -0.0636 0.45];
            otherwise
                disp('Unknown color choice');
                keyboard
        end
        cc = cc*.9;
        for i = 1:size(gl.trialspecs,1)
            gl.trialspecs{i,TS.CC} = cc;
        end
        gl.DEFAULTCONECONTRASTS = cc;
        gl.prefccdir = cc;
   end
   function ChangeSF(ev,h)
        global gl
        TS = TrialSpecCodes;
        
        uicontrols = getfield(get(gcf,'UserData'),'uicontrols');
        sf = get(uicontrols.changesf,'Value');
        for i = 1:size(gl.trialspecs,1)
            gl.trialspecs{i,TS.SF} = sf;
        end
        gl.DEFAULTSF = sf;
        gl.prefsf = sf;
   end

   function ChangeDiam(ev,h)
        global gl
        TS = TrialSpecCodes;
        
        uicontrols = getfield(get(gcf,'UserData'),'uicontrols');
        diam = get(uicontrols.changediam,'Value');
        for i = 1:size(gl.trialspecs,1)
            gl.trialspecs{i,TS.DIAM} = diam;
        end
        gl.DEFAULTDIAM = diam;
        gl.prefdiam = diam;
   end

    function ResetCallback(ev,h)
        global gl
        
        gl.preforient = nan;
        gl.prefsf = nan;
        gl.prefdiam = nan;
        gl.prefccdir = [0 0 0]; 
        SetUpTrials(1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots a tuning curve with a spline fit
% and updates gl.splineargmax (for efficiency).
% Replotting every single point each time we come through here (not adding
% to the plot).
function PlotTuningCurve()
    global gl;
    TS = TrialSpecCodes;
    
    spikerates = [];
    baselines = [];
    orientations  = [];
    sfs  = [];
    diams = [];
    phases = [];
    fit = [];
    x = [];
    
   for i = 1:size(gl.trialspecs,1)
        spikerate = gl.trialspecs{i,TS.SPIKES};
        spikerates = [spikerates; spikerate];
        orientations = [orientations; repmat(gl.trialspecs{i,TS.ORIENT}, length(spikerate),1)];
        sfs = [sfs; repmat(gl.trialspecs{i,TS.SF}, length(spikerate),1)];
        diams = [diams; repmat(gl.trialspecs{i,TS.DIAM}, length(spikerate),1)];
        phases = [phases; repmat(gl.trialspecs{i,TS.PHI}, length(spikerate),1)];
        baselines = [baselines; gl.trialspecs{i,TS.BASELINE}];
    end
    meanbaseline = mean(baselines);
    
    allorientations = unique([gl.trialspecs{:,TS.ORIENT}]);
    allsfs = unique([gl.trialspecs{:,TS.SF}]);
    alldiams = unique([gl.trialspecs{:,TS.DIAM}]);
    allphases = unique([gl.trialspecs{:,TS.PHI}]);

    % Setting up axes
    axesh = getfield(get(gcf,'UserData'),'axeshandles');
    set(axesh.tuningcurve,'Visible','on');
    set(axesh.colorspace,'Visible','off');
    for i = 1:length(axesh.colorspace)
        delete(get(axesh.colorspace(i),'Children'));
    end
    axes(axesh.tuningcurve);
    h = get(axesh.tuningcurve,'Children');
    delete(h);
   
    hold on;
    
    if (length(allorientations) > 2)
        set(gca,'XTick',[0:90:360],'XTickLabel',num2str([0:90:360]'),'XLim',[0 360],'XScale','linear');
        xlabel('orientation (deg)');
        parametervals = orientations;
        L = parametervals == min(parametervals);
        parametervals = [parametervals; parametervals(L)+2*pi];  % computing in rads, plotting in deg
        spikerates = [spikerates; spikerates(L)];
        splineopt = 'periodic';
        plotcoef = 180/pi;
    elseif (length(allsfs) > 2)
        set(gca,'XTick',unique([gl.trialspecs{:,TS.SF}]),'XTickLabel',num2str(unique([gl.trialspecs{:,TS.SF}]')),'XLim',[0 8],'XScale','log');
        xlabel('sf (cycles/deg)');
        parametervals = sfs;
        splineopt = 'variational';
        plotcoef = 1;
    elseif (length(alldiams) > 2)
        set(gca,'XTick',unique([gl.trialspecs{:,TS.DIAM}]),'XTickLabel',num2str(unique([gl.trialspecs{:,TS.DIAM}]')),'XLim',[0 3],'XScale','linear');
        xlabel('diam (deg)');
        parametervals = diams;
        splineopt = 'not-a-knot';
        plotcoef = 1;
    elseif (length(allphases) > 2)
        set(gca,'XTick',[0:90:360],'XTickLabel',num2str([0:90:360]'),'XLim',[0 360],'XScale','linear');
        xlabel('phase (deg)');
        parametervals = phases;
        means = nan*ones(length(allphases), 1);
        for i = 1:length(allphases)
            means(i) = mean(spikerates(find(phases == allphases(i))));
        end
        L = parametervals == min(parametervals);
        parametervals = [parametervals; parametervals(L)+2*pi];
        spikerates = [spikerates; spikerates(L)];
        splineopt = 'periodic';
        plotcoef = 180/pi;  % computing in rads, plotting in deg
      %  keyboard
    else
        return
    end
    
    if (length(unique(parametervals)) > 2)  % Compute the spline
        x = linspace(min(parametervals), max(parametervals),100);
        pp = csape(parametervals,spikerates,splineopt);
        fit = ppval(pp,x);
        gl.splineargmax = min(x(fit == max(fit))); % min to avoid ties
        plot(plotcoef*x,fit,'b-');
        plot(plotcoef*gl.splineargmax,max(fit),'m*','MarkerSize',5);
        plot(plotcoef*[min(parametervals) max(parametervals)], [meanbaseline meanbaseline], 'k:');
        % Updating the global indicating maximum of the tuning function
        if (length(allphases) > 2)  % Calculate F1/F0
            basis1 = exp(-2*pi*sqrt(-1)*[0:length(allphases)-1]/length(allphases))';
            F1 = abs(sum((means-meanbaseline)'*basis1));
            F0 = abs(sum(means-meanbaseline));
            text(allphases(1)+20,meanbaseline+.1*max(means),['F1/F0 = ',num2str(F1/F0)]);
        end
    end

    plot(plotcoef*parametervals, spikerates,'k.');
    if (length(parametervals) > 1)
        set(gca,'XLim',[min(parametervals) max(parametervals)]*plotcoef);
    end
    hold off;
end

function PlotColorTuning
    global gl
    TS = TrialSpecCodes;

    % Setting up axes
    axesh = getfield(get(gcf,'UserData'),'axeshandles');

    spikerates = [];
    baselines = [];
    ccs = [];
    for i = 1:size(gl.trialspecs,1)
        if (gl.trialspecs{i,TS.PROTOCOL} == 4)
            spikerate = gl.trialspecs{i,TS.SPIKES};
            spikerates = [spikerates; spikerate];
            ccs = [ccs; repmat(gl.trialspecs{i,TS.CC}, length(spikerate),1)];
            baselines = [baselines; gl.trialspecs{i,TS.BASELINE}];
        end
    end
    meanbaseline = mean(baselines);
 
    set(axesh.tuningcurve,'Visible','off');
    set(axesh.colorspace,'Visible','on');
    delete(get(axesh.tuningcurve,'Children'));
    for i = 1:length(axesh.colorspace)
        delete(get(axesh.colorspace(i),'Children'));
    end

    % Plotting
    basis = [1 1 0; 1 -1 0; 0 0 1];
    DKLs = (basis*ccs')';
    coordinates = sign(DKLs);
    peakfr = max([spikerates; .1]);  % .1 in case all spikerates = 0
    labels = {'L+M','L-M';'S','L-M';'S','L+M'};
    permmat = [3 1 2; 1 3 2; 2 3 1];
    for i = 1:3
        axes(axesh.colorspace(i));
        set(gca,'XLim',[-peakfr peakfr],'YLim',[-peakfr peakfr]);
        hold on;
        plot([0 0],[-peakfr peakfr],'k-'); plot([-peakfr peakfr],[0 0],'k-');
        % polar(linspace(0,2*pi,100), repmat(mean(baselines),1,100));
        text(0,1.5*peakfr,labels{i,1},'HorizontalAlignment','center');
        text(1.1*peakfr,.2*peakfr,labels{i,2});
        L = coordinates(:,permmat(i,1)) == 0;
        theta = atan2(coordinates(:,permmat(i,2)),coordinates(:,permmat(i,3)));
        [x,y] = pol2cart(theta(L), spikerates(L));
        plot(x, y,'k.')
        plot(-x,-y,'k.')
    end
end

function PlotContrastResponseFunctions
    global gl
    TS = TrialSpecCodes;

    % Setting up axes
    axesh = getfield(get(gcf,'UserData'),'axeshandles');

    spikerates = [];
    baselines = [];
    ccs = [];
    for i = 1:size(gl.trialspecs,1)
        if (gl.trialspecs{i,TS.PROTOCOL} == 7)
            spikerate = gl.trialspecs{i,TS.SPIKES};
            spikerates = [spikerates; spikerate];
            ccs = [ccs; repmat(gl.trialspecs{i,TS.CC}, length(spikerate),1)];
            baselines = [baselines; gl.trialspecs{i,TS.BASELINE}];
        end
    end
    meanbaseline = mean(baselines);
 
    set(axesh.tuningcurve,'Visible','on');
    set(axesh.colorspace,'Visible','off');
    delete(get(axesh.tuningcurve,'Children'));
    for i = 1:length(axesh.colorspace)
        delete(get(axesh.colorspace(i),'Children'));
    end
    axes(axesh.tuningcurve); hold on;
    set(gca,'XLim',[0 1],'XTick',[0 1],'XTickLabel',[0 1],'XScale','log');
    xlabel('Normalized contrast');
    Lblank = sum(abs(ccs),2) == 0;
    if (any(Lblank))
        baseline = mean(spikerates(Lblank));
        plot([0.0625 1],[baseline baseline],'k-');  % Ugly hardcoding.
    end
    if (any(~Lblank))
        L(:,1) = ccs*[1 -1 0]' ~= 0;
        L(:,2) = ccs*[0 0 1]' ~= 0;
        plotcols = {'red','blue'};
        for i = 1:2
            x = sum(ccs(L(:,i),:).^2,2);
            x = x./max(x);
            y = spikerates(L(:,i));
            plot(x,y,'k.','MarkerEdgeColor',char(plotcols{i}));
            uniquecontrasts = unique(x);
            meanresponses = [];
            for j = 1:length(uniquecontrasts)
                meanresponses(j) = mean(y(x == uniquecontrasts(j)));
            end
            plot(uniquecontrasts, meanresponses,'-','Color',char(plotcols{i}));
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the prefered cardinal and intermediate colors
% Also the color direction (in or out of isoluminant plane) that 
% evoked the greatest response.
function [prefisolumcolors, prefccdir] = computePrefColors()
    global gl
    TS = TrialSpecCodes;


    spikerates = [];
    baselines = [];
    ccs = [];
    for i = 1:size(gl.trialspecs,1)
        if (gl.trialspecs{i,TS.PROTOCOL} == 4)
            spikerate = gl.trialspecs{i,TS.SPIKES};
            spikerates = [spikerates; spikerate];
            ccs = [ccs; repmat(gl.trialspecs{i,TS.CC}, length(spikerate),1)];
            baselines = [baselines; gl.trialspecs{i,TS.BASELINE}];
        end
    end
 %   meanbaseline = mean(baselines);
    prefccdir = [];
    
    uniqueccs = unique(ccs,'rows');
    for a = 1:size(uniqueccs,1)
        tList = ismember(ccs, uniqueccs(a,:), 'rows');
        allRates(a,1) = mean(spikerates(tList)); % -baselines(tList)); % Commented out for now
    end
    prefccdir = uniqueccs(find(allRates==max(allRates),1),:);
    
    % Now getting prefisolumdir
    cardinals = [1 -1 0; 0 0 1];
    intermediates = [1 -1 -1; 1 -1 1];
    trialColors = sign(ccs);
    for a = 1:size(cardinals, 1);
        [tList, dummy] = ismember(trialColors, cardinals(a,:), 'rows');
        cardinalRates(a,1) = mean(spikerates(tList)-baselines(tList));
        cardinalCCS(a,:) = unique(ccs(tList,:), 'rows');
    end
    for a = 1:size(intermediates, 1);
        [tList, dummy] = ismember(trialColors, intermediates(a,:), 'rows');
        intermediateRates(a,1) = mean(spikerates(tList)-baselines(tList));
        intermediateCCS(a,:) = unique(ccs(tList,:),'rows');
    end
    
    %determine the perfered cardinal and intermediate colors. Send them to
    %rex as a vector where the first triplet is the "true" prefered color
    %(on the basis of maximum firing rate)
    [val, cardInd] = max(cardinalRates);
    [val, intInd] = max(intermediateRates);
    if cardinalRates(cardInd) > intermediateRates(intInd)
        prefisolumcolors(1,1:3) = cardinalCCS(cardInd,:);
        prefisolumcolors(1,4:6) = intermediateCCS(intInd,:);
    else
        prefisolumcolors(1,1:3) = intermediateCCS(intInd,:);
        prefisolumcolors(1,4:6) = cardinalCCS(cardInd,:);
    end
    disp(['prefcolor:', num2str(prefccdir)])
end


%%%%%%%%%%%%%%%%
% Raster functions below
function PlotRaster(spiketally)
    axesh = getfield(get(gcf,'UserData'),'axeshandles');
    axes(axesh.raster);
    plot([spiketally, spiketally]',[-.5; .5]*ones(1,length(spiketally)),'k-');
end

function ClearRaster
    axesh = getfield(get(gcf,'UserData'),'axeshandles');
    h = get(axesh.raster,'Children');  % Should be able to do this with cla but can't
    delete(h);
end

function UpdateTrialCountdown
    global gl;
    uicontrols = getfield(get(gcf,'UserData'),'uicontrols');
    TS = TrialSpecCodes;
    ntrials = sum([gl.trialspecs{:,TS.TRIALCOUNT}]);
    set(uicontrols.ntrialstextbox,'String',num2str(ntrials));
    set(uicontrols.protocoltextbox,'String',num2str(gl.trialspecs{gl.paramsidx,TS.PROTOCOL}));
end
% Here are the columns of the gl.trialspecs array:
% 1: Number of trials in this condition more to do
% 2: Driven spike counts
% 3: Baseline spike counts
% 4: Diameter
% 5: SF
% 6: TF
% 7: Cone contrasts
% 8: Orientation
% 9: Spatial phase
% 10: Number of cycles
% 11: Protocol number

function out = TrialSpecCodes
    out.TRIALCOUNT = 1;
    out.SPIKES = 2;
    out.BASELINE = 3;
    out.DIAM = 4;
    out.SF = 5;
    out.TF = 6;
    out.CC = 7;
    out.ORIENT = 8;
    out.PHI = 9;
    out.NCYCLES = 10;
    out.STIMTYPE = 11;    
    out.PROTOCOL = 12;
end