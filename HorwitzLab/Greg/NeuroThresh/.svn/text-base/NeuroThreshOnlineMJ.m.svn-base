function NeuroThreshOnlineMJ()

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
% 
% Setting things up for Mijung Park 
% GDLH 7/22/13
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

savefilename = ['MJ',num2str(now),'.mat',];
save(savefilename, 'gl');

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
        disp('abort');
        abortflag = 1;
        if (isempty(gl.trialspecs))
            return
        end
        % This might need to be fixed

        if (isfield(gl.MJ.totData,'r'))
            if (size(gl.MJ.totData.x,1) == size(gl.MJ.totData.r,1)+1)
                gl.MJ.totData.x(end,:) = [];
                size(gl.MJ.totData.x)
                size(gl.MJ.totData.r)
            end
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
        disp('correct');
        correct_t = p.times(find(p.events == C.CORRECTCD,1));
        a = get(gcf,'UserData');
        % Skipping the first few msecs as determined by PSTHslider
        t_offset = get(a.uicontrols.PSTHslider,'value');
        dur = stimoff_t-stimon_t-t_offset;
        if (isempty(spiketally))
            nspikes = 0;
        else
            nspikes = sum(spiketally<stimoff_t & spiketally>stimon_t+t_offset);
        end
        if (~isempty(gl.trialspecidx))
            nalg = length(gl.algorithms);
            
            if (mod(gl.trialspecidx,nalg) == 0)
                trialtype = nalg;
            else
                trialtype = mod(gl.trialspecidx,nalg);
            end
            whichalg = gl.algorithms(trialtype);
            if (strcmp(whichalg,'Adaptive'))  % GP adapative algorithm
                disp('updating r');
                gl.MJ.totData.r(end+1,1) = nspikes;

              %  size(gl.MJ.totData.x)
              %  size(gl.MJ.totData.r)
                
                if (size(gl.MJ.totData.r,1) == size(gl.MJ.totData.x,1))
                    if (rem(length(gl.MJ.totData.r), 10)==0) % GH this 10 should not be hardcoded
                        gl.MJ.totData.prs = estimateHyperparameters(gl.MJ.totData, gl.MJ.support, gl.MJ.maxsupport, gl.MJ.minsupport, gl.MJ.whichkernel);
                    end
                else
                    disp('x and r do not match');
              %      keyboard
                end
            end
            gl.trialspecidx = gl.trialspecidx + 1;
        elseif (gl.bldone) % We'll do this exactly once after "Got Baseline" button is clicked
            disp('setting trialspecidx to 1');
            gl.trialspecidx = 1;
        end
        p.lastprocessed_t = correct_t; 
    end
    %eot
    function eotfn
        disp('eot');
        eot_t = p.times(find(p.events == C.EOTCD,1));
        if (~abortflag && ~gl.bldone)
            PlotPSTH(spiketally, stimon_t);
        end
        nalg = length(gl.algorithms);
        realtrialidx = gl.trialspecidx-1;
        if (mod(realtrialidx,nalg) == 0)  % -1 becasue EOT is called after correctfn
            trialtype = nalg;
        else
            trialtype = mod(realtrialidx,nalg);
        end
        whichalg = gl.algorithms(trialtype);
        if (strcmp(whichalg,'Adaptive') & size(gl.MJ.totData.r,1) == size(gl.MJ.totData.x,1))
            if (isfield(gl.MJ.totData,'lambFinal'))  % 1= we into the adaptive phase, 0 = still in LH phase
                PlotGPsurface(gl.MJ.totData);
            end
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
             return; % Outdated 
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
    gl.initparams.stimdur = 0;
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
    gl.algorithms = {};
    gl.MJ = [];
    gl.MJ.totData.x = [];  % substructure for mijung's stuff
    gl.MJ.totData.r = [];  % substructure for mijung's stuff
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Getting the header data from the ecodes
% Destructively modifies the p structure to update p.lastprocessed_t.
function p = GetHeader(p)
    global gl;
    
    C = NeuroThreshCodes;

    % Creating a cell array that we iterate over to collect the
    % header parameters.
    % GH: add new header parameters? Which algorithm is being used?

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
        C.NPGRIDCD,'int',@npgridfn;
        C.WHICHKERNELCD,'int',@whichkernelfn;
        C.RANDSTIMSELCD,'int',@randstimselfn;
        C.LHSTIMSELCD,'int',@lhstimselfn;
        C.ADAPTSTIMSELCD,'int',@adaptstimselfn;
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
        % ---------------------------------
        % Setting up all of MiJung's stuff
        % Just hardcoding a bunch of stuff
        % ---------------------------------
        maxsupport_x = 1;
        minsupport_x = -1;
        maxsupport_y = 1;
        minsupport_y = -1;
        maxsupport_z = 0;
        minsupport_z = 1;
        gl.MJ.maxsupport = max([maxsupport_x maxsupport_y maxsupport_z]);
        gl.MJ.minsupport = min([minsupport_x minsupport_y minsupport_z]);
        numInitData = 3*gl.MJ.np_grid; % using these initial points, my algorithm will estimate hyperparameters, and then start adaptive sampling.
        gl.MJ.rgbscalefactors = min([gl.disp.bkgndrgb'; 1-gl.disp.bkgndrgb']);
        % initial points from LH design
        for i = 1:2
            if (i == 1)
                x_lh = lhsdesign(numInitData, 3);
                gl.MJ.totData.x = zeros(numInitData, 3);
                gl.MJ.totData.x(:,1) = bsxfun(@plus, bsxfun(@times, (maxsupport_x - minsupport_x), x_lh(:,1)), minsupport_x);
                gl.MJ.totData.x(:,2) = bsxfun(@plus, bsxfun(@times, (maxsupport_y - minsupport_y), x_lh(:,2)), minsupport_y);
                gl.MJ.totData.x(:,3) = bsxfun(@plus, bsxfun(@times, (maxsupport_z - minsupport_z), x_lh(:,3)), minsupport_z);
                gl.MJ.LHrgbs30 = gl.MJ.totData.x.*repmat(gl.MJ.rgbscalefactors,numInitData,1)+repmat(gl.disp.bkgndrgb',numInitData,1);
            else
                x_lh = lhsdesign(1000, 3); 
                x = zeros(1000, 3);
                x(:,1) = bsxfun(@plus, bsxfun(@times, (maxsupport_x - minsupport_x), x_lh(:,1)), minsupport_x);
                x(:,2) = bsxfun(@plus, bsxfun(@times, (maxsupport_y - minsupport_y), x_lh(:,2)), minsupport_y);
                x(:,3) = bsxfun(@plus, bsxfun(@times, (maxsupport_z - minsupport_z), x_lh(:,3)), minsupport_z);
                gl.MJ.LHrgbs1000 = x.*repmat(gl.MJ.rgbscalefactors,1000,1)+repmat(gl.disp.bkgndrgb',1000,1);
            end
        end
                
        stimrange_x = linspace(minsupport_x, maxsupport_x, gl.MJ.np_grid);
        stimrange_y = linspace(minsupport_y, maxsupport_y, gl.MJ.np_grid);
        stimrange_z = linspace(minsupport_z, maxsupport_z, gl.MJ.np_grid);
        [support_x,support_y,support_z] = meshgrid(stimrange_x,stimrange_y,stimrange_z);
        gl.MJ.support = [support_x(:),support_y(:),support_z(:)];
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
        gl.initparams.preflms = out;
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
    function npgridfn(out)
        gl.MJ.np_grid = out;
    end
    function whichkernelfn(out)
        % whichKernel = 1;  for non-symmetric stimulus space
        % whichKernel = 2;  for symmetric space (symmetric through the origin)
        gl.MJ.whichkernel = out;
    end
    function randstimselfn(out)
        if (out)
            gl.algorithms{length(gl.algorithms)+1}= 'Random';
        end
    end
    function lhstimselfn(out)
        if (out)
            gl.algorithms{length(gl.algorithms)+1}= 'LH';
        end
    end
    function adaptstimselfn(out)
        if (out)
            gl.algorithms{length(gl.algorithms)+1}= 'Adaptive';
        end
    end
end

%%%%%
% Preparing and sending trial parameters when requested to by REX
% Called by REX after rewoff (so we know the just-past trial was correct)
% and before droptrialparams (so that the numbers returned by this function
% make it into the data file).
function out = GetTrialParams()  % GH send up other parameters?
    global gl
    % By this point gl.trialspecidx has *already* been incremented! (Happens
    % in the correct trial state. Need to deduct '1'. Ugly.
    if (~isempty(gl.trialspecs))
        nalg = length(gl.algorithms);
        realtrialidx = gl.trialspecidx-1;        
        if (mod(realtrialidx,nalg) == 0)
            trialtype = nalg;
        else
            trialtype = mod(realtrialidx,nalg);
        end
    else
        trialtype = 0;
    end
    % Could use these fields to drop alternative representation of x,y,z
    out = [trialtype 1 1 1 1 1 1];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Randomly picks a valid color direction, finds the next contrast
% and returns the cone contrast triplet at the peak of the Gabor.
function out = PickGaborParams()
    global gl

    a = get(gcf,'UserData');
    if (isempty(gl.trialspecidx)) % If still in the latency measuring part
        if (~get(a.uicontrols.mkachrom,'Value') & ~isempty(gl.initparams.preflms))
            out = gl.initparams.preflms;
        else
            out = [.5 .5 .5]; % arbitrary
        end
        pause(.5);
    else
        % Mijung's stuff
    	nalg = length(gl.algorithms);
        if (mod(gl.trialspecidx,nalg) == 0)
            whichalg = gl.algorithms(nalg);
        else
            whichalg = gl.algorithms(mod(gl.trialspecidx,nalg));
        end
        % disp(char(whichalg));
        if (strcmp(whichalg,'Random'))  % Uniform random RGB
            disp('Random');
            rgb = unifrnd(-1,1,1,3).*gl.MJ.rgbscalefactors+gl.disp.bkgndrgb';
            pause(.5);

        elseif (strcmp(whichalg,'LH'))  % Latin hypercube with 1000 points
            disp('LH');
            rgb = gl.MJ.LHrgbs1000(ceil(gl.trialspecidx/nalg),:);
            pause(.5);

        elseif (strcmp(whichalg,'Adaptive'))  % GP adapative algorithm
            if (gl.trialspecidx <= nalg*3*gl.MJ.np_grid); % '3' in this line is for the 3-D stimulus space
                % initial points from LH design
                disp('initial LH');
                rgb = gl.MJ.LHrgbs30(ceil(gl.trialspecidx/nalg),:);
                pause(.5);
            else % we're done with LH and have moved on to GP part
                % Stuff for adaptive algorithm
                
%                 if (rem(gl.trialspecidx, 10)==0)&(gl.trialspecidx>nalg*3*gl.MJ.np_grid)
                disp('adaptive');
                try
            %    [nextX, totData] = computeNextStim_ALalgorithm_NewModel_varyingP(gl.MJ.totData.x, gl.MJ.totData.r, 3, gl.MJ.support, 3*gl.MJ.np_grid, gl.MJ.maxsupport, gl.MJ.minsupport, gl.MJ.whichkernel);
            
                size(gl.MJ.totData.x)
                size(gl.MJ.totData.r)
                
                [nextX, totData] = computeNextStim_ALalgorithm_NewModel_varyingP(gl.MJ.totData.x, gl.MJ.totData.r, gl.MJ.totData.prs, 3, gl.MJ.support,  3*gl.MJ.np_grid, gl.MJ.whichkernel);
                catch
                    keyboard
                end
                rgb = nextX.*gl.MJ.rgbscalefactors+gl.disp.bkgndrgb';
                totData.x(end+1,:) = nextX;
                gl.MJ.totData = totData; 
            end
        else
            error('Unknown algorithm');
        end
        
        
        %    if (strcmp(whichalg,'Adaptive')&(gl.trialspecidx <= nalg*3*gl.MJ.np_grid))
        if (~isempty(gl.MJ.totData.r) & (strcmp(whichalg,'Adaptive')))
            try
                disp([gl.trialspecidx gl.MJ.totData.x(length(gl.MJ.totData.r),:) gl.MJ.totData.r(length(gl.MJ.totData.r),:)]);
            catch
                keyboard
            end
        else
            disp(gl.trialspecidx);
        end
        tmp = gl.disp.M*rgb';  % Taking random rgbs and turning them to lms
        bkgndlms = gl.disp.M*gl.disp.bkgndrgb;
        out = (tmp-bkgndlms)./bkgndlms;
        gl.trialspecs(gl.trialspecidx).cc = out;
        % end of debugging code
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting & analysis functions below...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SetUpFigure()

    figure(1);
    a = get(gcf,'UserData');  % In case there's already things there (eg from whitenoise)
    if (isfield(a,'uicontrols'))
        staticcontrols = {'PSTHslider','PSTHval',0;...
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
        PSTHval = 0;
        monoval = 0;
        replay = 0;
    end
        
    set(gcf,'position',[200 200 900 400])
    set(gcf,'DefaultAxesUnits','pixels')
    clf;
    a.axes.h1 = axes('position',[50 50 250 150]);
    a.axes.h2 = axes('position',[350 50 250 150]);
    a.uicontrols.mkachrom = uicontrol('style','checkbox','Min',0,'Max',1,'Position',[200 360 15 15]);
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
       
       [gl.trialspecs.nspikes] = deal([]);
       [gl.trialspecs.cc] = deal([]);
       set(a.uicontrols.PSTHslider,'Value',.0);
       PlotPSTH(nan, nan);
    end

    % Once threshold and latency parameters have been set, user clicks
    % the bldone button whereupon new color directions are picked,
    % latency value is packed into messages for REX, and uicontrols 
    % are disabled. 
    function BlDoneCallback(ev,h)
        global gl;

        a = get(gcf,'UserData');
        axes(a.axes.h1);
        set(a.uicontrols.bldone,'Enable','off');
        set(a.uicontrols.PSTHslider,'Enable','off');
        set(a.uicontrols.clearspikes,'Enable','off');
        set(a.uicontrols.polarityflip,'Enable','off');
        set(a.uicontrols.mkachrom,'Enable','off');
        val = get(a.uicontrols.PSTHslider,'Value');
        gl.outmsg = [gl.outmsg; 'sendToRex(udpCom,',num2str(val),',''double'', ''Latency'');'];
        gl.bldone = 1;
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
        L = reshape([gl.trialspecs.cc],3,n) == repmat([0 0 0]',1,n);
        idx = find(~all(L),1);
        if (~isempty(idx) && ~isempty(gl.trialspecs(idx).cc))
            gl.trialspecs(idx).cc = -gl.trialspecs(idx).cc;
        end
    end
end

% Clear the status text box after each trial
function ClearStatusText()
    a = get(gcf,'UserData');
    set(a.uicontrols.statustext,'String',[],'BackgroundColor',[0.92549 0.913725 0.847059]);
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

% Plotting Mijung's surface
% Is this code supposed to destructively modify "totData" because it
% doesn't the way's it's implemented now. GH
function PlotGPsurface(totData)
    global gl
%     if (size(totData.x,1) == size(totData.r,1)+1)
%         disp('Using Mijung''s "suggestion"');
%         totData.x = totData.x(1:end-1,:); % This was Mijung's idea
%     end
%     
%     minsupport_x = -1;
%     maxsupport_x = 1;
%     minsupport_y = -1;
%     maxsupport_y = 1;
%     
%     stimrange_x = linspace(minsupport_x, maxsupport_x, gl.MJ.np_grid);
%     stimrange_y = linspace(minsupport_y, maxsupport_y, gl.MJ.np_grid);
%     [aa, bb] = meshgrid(stimrange_x, stimrange_y);
%     
%     totData.support = [aa(:) bb(:) zeros(length(aa(:)), 1)];
%     if gl.MJ.whichkernel == 2
%         [totData.norm_mat_Kstar, totData.norm_mat_Kstar2] = form_normMat(totData.support, totData.x);
%         totData.Kstar = totData.prs(2)*exp(-.5/totData.prs(3).*totData.norm_mat_Kstar) + totData.prs(2)*exp(-.5/totData.prs(3).*totData.norm_mat_Kstar2);
%         [totData.norm_mat_support, totData.norm_mat_support2] = form_normMat(totData.support, totData.support);  % squared distance
%     else
%         totData.norm_mat_Kstar = form_normMat(totData.support, totData.x);
%         totData.Kstar = totData.prs(2)*exp(-.5/totData.prs(3).*totData.norm_mat_Kstar);
%         totData.norm_mat_support = form_normMat(totData.support, totData.support);  % squared distance
%     end
%     [mean_lambda_LM] = makePrediction_nse(totData.prs, totData, totData.aFinal, totData.WFinal, totData.sqrtLFinal);
     a = get(gcf,'UserData');
     axes(a.axes.h1);
%     contour(aa, bb, reshape(mean_lambda_LM, gl.MJ.np_grid, []), 5); title('LM');
    plot(gl.MJ.totData.lambFinal)
    
end


% -----------------------------
% Support functions for Mijung's stuff
% -----------------------------
function [neglogev, hmapR, a, W, sqrtL]  = computeFmap_nse(datastruct)

% given H, find hmap
h0 = datastruct.hinit;
thresh = 0.05;
diffIn_h = 1;
maxIter = 50;
count = 1;
while (diffIn_h>thresh)&&(count<=maxIter)
    
    [diffIn_h, hmapR, a, W, sqrtL, obj] = objFun_NewtonsMethod_nse(h0, datastruct);
    if count ==1
        datastruct.obj = obj;
    end
    %     diffInF
    % check if the obj is increasing
    if count >1
        if obj < datastruct.obj
            h0 = (hmapR+h0)/2; % if obj is decreasing, choose a smaller step
        else
            datastruct.obj = obj;
            h0 = hmapR;
        end
    end
    count = count+1;
    
    
end

% 3. given fmap, compute neglogev
neglogev = computeLogevidence_nse(hmapR, W, datastruct.muf, a, datastruct.r, datastruct.pest);
end
function  [prs, hmapFinal, aFinal, WFinal, sqrtLFinal, neglogev0, detH] = computeFmapAndUpdateTheta_nse(prs0, datastruct)

threshEvDiff = 0.1;

count = 1;
maxIter = 10;

% 1. given theta0, form K
[neglogev0, hmapR0, a0, W0, sqrtL0]  = computeFmap_nse(datastruct);
datastruct.L = sqrtL0.^2;
datastruct.Lminit = datastruct.L*hmapR0 + a0;

datastruct.neglogev = neglogev0;

while(count<maxIter)
    
    % 0. from theta0, optimize hyperparameters using analytic form
    [neglogev1, prs, hmapFinal, aFinal, WFinal, sqrtLFinal, Hfinal, detH] = updateThetaGivenL_nse(prs0, datastruct);
    datastruct.hinit =  hmapFinal;
    datastruct.H = Hfinal; % H = K + nsevar I
    datastruct.muf = abs(prs(1));
    datastruct.pest = abs(prs(5));
    
    % 1. given theta0, form H
    [neglogev0, hmapR0, a0, W0, sqrtL0]  = computeFmap_nse(datastruct);
    
    evDiff = datastruct.neglogev - neglogev0; % should be positive, if it goes to the right direction
    
    if (abs(evDiff))<= threshEvDiff
        %         neglog = datastruct.neglog0;
        return;
    elseif evDiff>0
        prs0 = prs;
        datastruct.neglogev = neglogev0;
        datastruct.L = sqrtL0.^2;
        datastruct.Lminit = datastruct.L*hmapR0 + a0;
    else
        prs0 = (prs0+prs)./2;
        %         datastruct.finit = fmapR0;
    end
    
    count = count +1;
    
end
end
function neglogev = computeLogevidence_nse(hmap, W, muf, a, r, p)

g = logexp1(hmap, p);
g(g<1e-6) = 1e-6; 
g(isinf(g)) = 1e6;

logev = r'*log(g) - sum(g) - logdetns(W) - 0.5*(hmap - muf)'*a;

neglogev = - logev;
end
function [neglogev, hmap, a, W, sqrtL, H]  = computeMarginal_usingAnalyticForm_nse(prs, datastruct)

% 1. given theta, form H
muf = abs(prs(1));
param = abs(prs(2:3));
nsevar = abs(prs(4));
p = abs(prs(5));

norm_mat = datastruct.norm_mat;

if datastruct.whichKernel == 2
    norm_mat2 = datastruct.norm_mat2;
    K = param(1)*exp(-.5/param(2).*norm_mat) + param(1)*exp(-.5/param(2).*norm_mat2);
else
    K = param(1)*exp(-.5/param(2).*norm_mat);
end

nsevarI = nsevar*eye(size(norm_mat)); 
H = K + nsevarI;

% 2. given K, find fmap
Lm_init = datastruct.Lminit;
sqrtL = sqrt(datastruct.L);

B = eye(size(sqrtL)) + sqrtL*H*sqrtL;

W = chol(B, 'lower');

M = W\sqrtL;
sqrtLinvBsqrtL = M'*M;
a = Lm_init - sqrtLinvBsqrtL*(H*Lm_init + muf);

hmap = H*a + muf;

% 3. given fmap, compute neglogev
neglogev = computeLogevidence_nse(hmap, W, muf, a, datastruct.r, p);
end

function [nextX, totData] = computeNextStim_ALalgorithm_NewModel_varyingP(x, r, prs, stimDim, support,  numInitData, whichKernel)

%% to compute next stimulus using active learning algorithm under a new model for overdispersed data
% we do not fix nonlinearity here; we estimate the order of nonlinearity
% from the data. (p is the order of nonlinearity)

% --------- input -----------
% 1. x, r: all the stimuli/response pairs
% 2. prs: hyperparameters
% 3. stimDim: stimulus dimension
% 4. support: grid of points where we estimate tuning curve
% 5. numbInit: initial number of data points
% 5. whichKernel: the new kernel for symmetric space (#2) or squared distance kernel (#1)

% --------- output -----------
% 1. nextX: selected stimulus to present next
% 2. totData: datastructure that has important quantities
%    totData.x: all the stimuli presented
%    totData.r: all the responses measured
%    totData.support: grid of points where we estimate tuning curve
%    totData.ndim: stimulus dimension
%    totData.norm_mat_support: matrix with squared distance of support
%    totData.norm_mat: matrix with squared distance of x
%    totData.norm_mat_Kstar: matrix with squared distance between x and support
%    totData.hinit: initial values of hmap
%    totData.g: nonlinearity (exponential)
%    totData.ginv: inverse nonlinearity (log)
%    totData.nstim: number of stimuli presented
%    totData.prs: estimated hyperparameters
%    totData.lambFinal: estimated tuning curve 

%%

global datastruct;

% to avoid numerical problem with r=0
thTrial = length(r);
r(r==0) = 0.1; 

%% construct data structure

if thTrial ==numInitData
    
    datastruct.x = x;
    datastruct.r = r;
    datastruct.support = support;
    datastruct.ndim = stimDim;
    datastruct.whichKernel = whichKernel;
    
    g = @(a, p) (log(exp(a)+1)).^p;
    ginv = @(a, p) log(exp(a.^(1./p))-1); % here r is non-negative
    
    datastruct.g = g;
    datastruct.ginv = ginv;
    
    datastruct.pinit = 1;
    datastruct.hinit = datastruct.ginv(datastruct.r+0.1, datastruct.pinit);
    
    if whichKernel==2
        % if it's symmetric kernel
        [datastruct.norm_mat_support, datastruct.norm_mat_support2] = form_normMat(support, support);
        [datastruct.norm_mat, datastruct.norm_mat2] = form_normMat(datastruct.x, datastruct.x);  
        [datastruct.norm_mat_Kstar, datastruct.norm_mat_Kstar2] = form_normMat(datastruct.support, datastruct.x);
    else
        % kernel with squared distance 
        datastruct.norm_mat_support = form_normMat(support, support);
        datastruct.norm_mat = form_normMat(datastruct.x, datastruct.x);  
        datastruct.norm_mat_Kstar = form_normMat(datastruct.support, datastruct.x);
    end
       
else
    
    datastruct.x = x;
    datastruct.r = r;
    
    try
            
        if whichKernel==2
            % if it's symmetric kernel
            
            [sqrDist_new, sqrDist_new2] = form_normMat(datastruct.x(end,:), datastruct.x);
            datastruct.norm_mat = [[datastruct.norm_mat; sqrDist_new(1:end-1)] sqrDist_new'];
            datastruct.norm_mat2 = [[datastruct.norm_mat2; sqrDist_new2(1:end-1)] sqrDist_new2'];

            [normMat_Kstar_new, normMat_Kstar_new2] = form_normMat(datastruct.support, datastruct.x(end,:));
            datastruct.norm_mat_Kstar = [datastruct.norm_mat_Kstar normMat_Kstar_new];
            datastruct.norm_mat_Kstar2 = [datastruct.norm_mat_Kstar2 normMat_Kstar_new2];
            
        else
                        
            sqrDist_new = form_normMat(datastruct.x(end,:), datastruct.x);
            datastruct.norm_mat = [[datastruct.norm_mat; sqrDist_new(1:end-1)] sqrDist_new'];
            
            normMat_Kstar_new = form_normMat(datastruct.support, datastruct.x(end,:));
            datastruct.norm_mat_Kstar = [datastruct.norm_mat_Kstar normMat_Kstar_new];
        end
        
        datastruct.hinit = [datastruct.hinit; datastruct.ginv(r(end)+0.1, datastruct.pest)];
        
    catch
        
        disp('wow'); % if it breaks, then generate the distance matrix from scratch

        if whichKernel==2
            % if it's symmetric kernel
            [datastruct.norm_mat, datastruct.norm_mat2] = form_normMat(datastruct.x, datastruct.x);
            [datastruct.norm_mat_Kstar, datastruct.norm_mat_Kstar2] = form_normMat(datastruct.support, datastruct.x);
        else
            % kernel with squared distance
            datastruct.norm_mat = form_normMat(datastruct.x, datastruct.x);
            datastruct.norm_mat_Kstar = form_normMat(datastruct.support, datastruct.x);
        end
        datastruct.hinit = datastruct.ginv(datastruct.r+0.1, datastruct.pest);
    end     
        
end

% update number of stimulus 
datastruct.nstim = length(datastruct.r);

%% update fmap given hyperparameters

% update fmap given hyperparameters

if whichKernel==2
    K = prs(2)*exp(-.5/abs(prs(3)).*datastruct.norm_mat)+prs(2)*exp(-.5/prs(3).*datastruct.norm_mat2);
else
    K = prs(2)*exp(-.5/abs(prs(3)).*datastruct.norm_mat);
end
datastruct.H = K + abs(prs(4))*eye(size(datastruct.norm_mat));
datastruct.muf = abs(prs(1));
datastruct.pest = prs(5); % p_est

[neglogev, hmapFinal, aFinal, WFinal, sqrtLFinal]  = computeFmap_nse(datastruct);

datastruct.neglogev = neglogev;
datastruct.hinit = hmapFinal;

%% select next stimulus

if whichKernel ==2
    datastruct.Kstar = prs(2)*exp(-.5/prs(3).*datastruct.norm_mat_Kstar) + prs(2)*exp(-.5/prs(3).*datastruct.norm_mat_Kstar2); 
else
    datastruct.Kstar = prs(2)*exp(-.5/prs(3).*datastruct.norm_mat_Kstar); 
end

[mean_lambda, nextX, idxNext, var_lambda] = makePrediction_nse(prs, datastruct, aFinal, WFinal, sqrtLFinal);

datastruct.idxNext = idxNext;
datastruct.nextX = nextX;
datastruct.lambFinal = mean_lambda;
datastruct.varlamb = var_lambda;
datastruct.aFinal = aFinal;
datastruct.WFinal = WFinal;
datastruct.sqrtLFinal = sqrtLFinal;
datastruct.prs = prs; 

totData = datastruct;
end

%  function [nextX, totData] = computeNextStim_ALalgorithm_NewModel_varyingP(x, r, prs, stimDim, support,  numInitData, whichKernel)
% 
% %% to compute next stimulus using active learning algorithm under a new model for overdispersed data
% % we do not fix nonlinearity here; we estimate the order of nonlinearity
% % from the data. (p is the order of nonlinearity)
% 
% % --------- input -----------
% % 1. x, r: all the stimuli/response pairs
% % 2. prs: hyperparameters
% % 3. stimDim: stimulus dimension
% % 4. support: grid of points where we estimate tuning curve
% % 5. numbInit: initial number of data points
% % 5. whichKernel: the new kernel for symmetric space (#2) or squared distance kernel (#1)
% 
% % --------- output -----------
% % 1. nextX: selected stimulus to present next
% % 2. totData: datastructure that has important quantities
% %    totData.x: all the stimuli presented
% %    totData.r: all the responses measured
% %    totData.support: grid of points where we estimate tuning curve
% %    totData.ndim: stimulus dimension
% %    totData.norm_mat_support: matrix with squared distance of support
% %    totData.norm_mat: matrix with squared distance of x
% %    totData.norm_mat_Kstar: matrix with squared distance between x and support
% %    totData.hinit: initial values of hmap
% %    totData.g: nonlinearity (exponential)
% %    totData.ginv: inverse nonlinearity (log)
% %    totData.nstim: number of stimuli presented
% %    totData.prs: estimated hyperparameters
% %    totData.lambFinal: estimated tuning curve 
% 
% %%
% 
% global datastruct;
% 
% % to avoid numerical problem with r=0
% thTrial = length(r);
% r(r==0) = 0.1; 
% 
% %% construct data structure
% 
% if thTrial ==numInitData
%    
% else
%     
%     datastruct.x = x;
%     datastruct.r = r;
%     
%     try
%             
%         if whichKernel==2
%             % if it's symmetric kernel
%             
%             [sqrDist_new, sqrDist_new2] = form_normMat(datastruct.x(end,:), datastruct.x);
%             datastruct.norm_mat = [[datastruct.norm_mat; sqrDist_new(1:end-1)] sqrDist_new'];
%             datastruct.norm_mat2 = [[datastruct.norm_mat2; sqrDist_new2(1:end-1)] sqrDist_new2'];
% 
%             [normMat_Kstar_new, normMat_Kstar_new2] = form_normMat(datastruct.support, datastruct.x(end,:));
%             datastruct.norm_mat_Kstar = [datastruct.norm_mat_Kstar normMat_Kstar_new];
%             datastruct.norm_mat_Kstar2 = [datastruct.norm_mat_Kstar2 normMat_Kstar_new2];
%             
%         else
%                         
%             sqrDist_new = form_normMat(datastruct.x(end,:), datastruct.x);
%             datastruct.norm_mat = [[datastruct.norm_mat; sqrDist_new(1:end-1)] sqrDist_new'];
%             
%             normMat_Kstar_new = form_normMat(datastruct.support, datastruct.x(end,:));
%             datastruct.norm_mat_Kstar = [datastruct.norm_mat_Kstar normMat_Kstar_new];
%         end
%         
%         datastruct.hinit = [datastruct.hinit; datastruct.ginv(r(end)+0.1, datastruct.pest)];
%         
%     catch
%         
%         disp('wow'); % if it breaks, then generate the distance matrix from scratch
% 
%         if whichKernel==2
%             % if it's symmetric kernel
%             [datastruct.norm_mat, datastruct.norm_mat2] = form_normMat(datastruct.x, datastruct.x);
%             [datastruct.norm_mat_Kstar, datastruct.norm_mat_Kstar2] = form_normMat(datastruct.support, datastruct.x);
%         else
%             % kernel with squared distance
%             datastruct.norm_mat = form_normMat(datastruct.x, datastruct.x);
%             datastruct.norm_mat_Kstar = form_normMat(datastruct.support, datastruct.x);
%         end
%         datastruct.hinit = datastruct.ginv(datastruct.r+0.1, datastruct.pest);
%     end     
%         
% end
% 
% % update number of stimulus 
% datastruct.nstim = length(datastruct.r);
% 
% %% update fmap given hyperparameters
% 
% % update fmap given hyperparameters
% 
% if whichKernel==2
%     K = prs(2)*exp(-.5/abs(prs(3)).*datastruct.norm_mat)+prs(2)*exp(-.5/prs(3).*datastruct.norm_mat2);
% else
%     K = prs(2)*exp(-.5/abs(prs(3)).*datastruct.norm_mat);
% end
% datastruct.H = K + abs(prs(4))*eye(size(datastruct.norm_mat));
% datastruct.muf = abs(prs(1));
% datastruct.pest = prs(5); % p_est
% 
% [neglogev, hmapFinal, aFinal, WFinal, sqrtLFinal]  = computeFmap_nse(datastruct);
% 
% datastruct.neglogev = neglogev;
% datastruct.hinit = hmapFinal;
% 
% %% select next stimulus
% 
% if whichKernel ==2
%     datastruct.Kstar = prs(2)*exp(-.5/prs(3).*datastruct.norm_mat_Kstar) + prs(2)*exp(-.5/prs(3).*datastruct.norm_mat_Kstar2); 
% else
%     datastruct.Kstar = prs(2)*exp(-.5/prs(3).*datastruct.norm_mat_Kstar); 
% end
% 
% [mean_lambda, nextX, idxNext, var_lambda] = makePrediction_nse(prs, datastruct, aFinal, WFinal, sqrtLFinal);
% 
% datastruct.idxNext = idxNext;
% datastruct.nextX = nextX;
% datastruct.lambFinal = mean_lambda;
% datastruct.varlamb = var_lambda;
% datastruct.aFinal = aFinal;
% datastruct.WFinal = WFinal;
% datastruct.sqrtLFinal = sqrtLFinal;
% datastruct.prs = prs; 
% 
% totData = datastruct;
% end

function hyparamEstimates = estimateHyperparameters(datastruct, support, maxsupport, minsupport, whichKernel)

% ----------------------------------------------
% The stuff below was cut and pasted from the begining of
% computeNexStim_AL....  I'mk not sure what it does. GDLH
% if (~isfield(datastruct,'ginv')) % first time through?  GH
%     datastruct.support = support;
%     datastruct.whichKernel = whichKernel;

g = @(a, p) (log(exp(a)+1)).^p;
ginv = @(a, p) log(exp(a.^(1./p))-1); % here r is non-negative

datastruct.g = g;
datastruct.ginv = ginv;

% end
% ----------------------------------------------

% optimize hyperparameters with analytic form
datastruct.pinit = 1;
datastruct.hinit = datastruct.ginv(datastruct.r+0.1, datastruct.pinit);

nsevar_init = 1;
ovrscl_1 = datastruct.ginv(mean(datastruct.r), datastruct.pinit); % overall scale
lngthscl_1 = (maxsupport-minsupport)/2; % variance
prs0 = [datastruct.ginv(mean(datastruct.r), datastruct.pinit); ovrscl_1; lngthscl_1; nsevar_init; datastruct.pinit];

% kernel matrix

if whichKernel==2
    % if it's symmetric kernel
    [datastruct.norm_mat_support, datastruct.norm_mat_support2] = form_normMat(support, support);
    [datastruct.norm_mat, datastruct.norm_mat2] = form_normMat(datastruct.x, datastruct.x);
    K = prs0(2)*exp(-.5/prs0(3).*datastruct.norm_mat)+prs0(2)*exp(-.5/prs0(3).*datastruct.norm_mat2);
else
    % kernel with squared distance
    datastruct.norm_mat_support = form_normMat(support, support);
    datastruct.norm_mat = form_normMat(datastruct.x, datastruct.x);
    K = prs0(2)*exp(-.5/prs0(3).*datastruct.norm_mat);
end

datastruct.H = K + prs0(4)*eye(size(datastruct.norm_mat));
datastruct.muf = prs0(1);
datastruct.pest = prs0(5);
datastruct.whichKernel = whichKernel; 

hyparamEstimates = computeFmapAndUpdateTheta_nse(prs0, datastruct);

end
function [norm_mat_subtract, norm_mat_sum] = form_normMat(x1, x2)

% compute the l-2 norm
x_dim = length(x1(1,:));
nsamps_x1 = size(x1, 1);
nsamps_x2 = size(x2, 1);

xn_i = repmat(x1, 1, nsamps_x2);
xn = reshape(xn_i', x_dim, []);
xn = xn';
xm = repmat(x2, nsamps_x1, 1);

if nargout ==1
    
    norm_vec = sum(abs(xn-xm).^2, 2);
    norm_mat_tr = reshape(norm_vec, nsamps_x2, nsamps_x1);
    norm_mat_subtract = norm_mat_tr';
    
else

    norm_vec1 = sum(abs(xn-xm).^2, 2);
    norm_mat_tr1 = reshape(norm_vec1, nsamps_x2, nsamps_x1);
    norm_mat_subtract = norm_mat_tr1';

        
    norm_vec2 = sum(abs(xn+xm).^2, 2); 
    norm_mat_tr2 = reshape(norm_vec2, nsamps_x2, nsamps_x1);
    norm_mat_sum = norm_mat_tr2';
    
end

end
function fn_val = GaussHermite(func, npt, varargin) 
%
%       HERMITE INTEGRATION, I.E. EVALUATION OF THE INTEGRAL OF
%       exp(-x**2)*f(x) FROM -INFINITY TO +INFINITY
%
%       NPT     = INPUT, NO. OF POINTS AT WHICH f(x) IS TO BE EVALUATED.
%                 NPT MUST BE ONE OF 2, 4, 6, 8, 10, 12, 16, 20.
%       FUNC    = INPUT, NAME OF THE USER'S FUNCTION TO SUPPLY VALUES
%                 OF F(X).
%       varargin contains any extra parameters for FUNC
%
% Based on hermite.f90 from Alan Miller's Fortran site:
% http://users.bigpond.net.au/amiller/
% Converted from Fortran90 using f2matlab + manual editing
% Code rewritten with MATLAB constructs
% S Bocquet 28 July 2008 Tested with MATLAB Version 7.6.0.324 (R2008a)
%************************************************************************

persistent npts ipos wt xpt 

%     Local variables

if isempty(xpt), xpt(1:39) =[ 0.707106781186548,0.524647623275290, 1.650680123885785,0.436077411927617, 1.335849074013697, 2.350604973674492,0.381186990207322, 1.157193712446780, 1.981656756695843,2.930637420257244,0.342901327223705, 1.036610829789514, 1.756683649299882,2.532731674232790, 3.436159118837738,0.314240376254359, 0.947788391240164, 1.597682635152605,2.279507080501060, 3.020637025120890, 3.889724897869782,0.27348104613815,  0.82295144914466, 1.38025853919888,1.95178799091625,  2.54620215784748, 3.17699916197996,3.86944790486012,  4.68873893930582,0.2453407083009,   0.7374737285454, 1.2340762153953,1.7385377121166,   2.2549740020893, 2.7888060584281,3.3478545673832,   3.9447640401156, 4.6036824495507,5.3874808900112 ]; end;
if isempty(wt), wt(1:39) =[ 8.862269254528d-1,8.049140900055d-1, 8.131283544725d-2,7.246295952244d-1, 1.570673203229d-1, 4.530009905509d-3,6.611470125582d-1, 2.078023258149d-1, 1.707798300741d-2,1.996040722114d-4,6.108626337353d-1, 2.401386110823d-1, 3.387439445548d-2,1.343645746781d-3, 7.640432855233d-6,5.701352362625d-1, 2.604923102642d-1, 5.160798561588d-2,3.905390584629d-3, 8.573687043588d-5, 2.658551684356d-7,5.079294790166d-1, 2.806474585285d-1, 8.381004139899d-2,1.288031153551d-2, 9.322840086242d-4, 2.711860092538d-5,2.320980844865d-7, 2.654807474011d-10,4.622436696006d-1, 2.866755053628d-1, 1.090172060200d-1,2.481052088746d-2, 3.243773342238d-3, 2.283386360163d-4,7.802556478532d-6, 1.086069370769d-7, 4.399340992273d-10,2.229393645534d-13 ]; end;
if isempty(npts), npts(1:8) = int8([ 2, 4, 6, 8, 10, 12, 16, 20 ]); end;
if isempty(ipos), ipos(1:9) = int8([ 1, 2, 4, 7, 11, 16, 22, 30, 40 ]); end;

%       CHECK FOR PERMISSIBLE VALUE OF NPT.

j = find(npt==npts);
if isempty(j)
  error('GaussHermite:InvalidNpt','Invalid number of points: npt must be 2,4,6,8,10,12,16 or 20')
end

%       EVALUATE SUM OF WT(I) * FUNC(X(I))

i = ipos(j):ipos(j+1)-1;
fn_val = dot(wt(i),(func(xpt(i),varargin{:}) + func(-xpt(i),varargin{:})));

return

end %function GaussHermite
function x = logdetns(A)
% LOGDET - computes the log-determinant of a matrix A using Cholesky or LU
% factorization
%
% LOGDET
%
% x = logdet(A);
%
% This is faster and more stable than using log(det(A))
%
% Input:
%     A NxN - A must be sqaure, positive SYMMETRIC and semi-definite
%     (Chol assumes A is symmetric)

% if checksym(A)  % crude symmetry check
% 
%     x = 2*sum(log(diag(chol(A))));
% 
% else
    
    x = sum(log(abs(diag(lu(A)))));

% end

% ----------------------------
% function z = checksym(A)
% % Crude (but fast) algorithm for checking symmetry, just by looking at the
% % first column and row
% 
% z = isequal(A(:,1),A(1,:)');
end
function [f,df,ddf] = logexp1(x, pow)
%  [f,df,ddf] = logexp1(x);
%
%  Implements the nonlinearity:  
%     f(x) = log(1+exp(x)).^pow;
%  Where pow = 1;
%  plus first and second derivatives
%


% pow = 1;

f0 = log(1+exp(x));
f = f0.^pow;

if nargout > 1
    df = pow*f0.^(pow-1).*exp(x)./(1+exp(x));
end
if nargout > 2
    ddf = pow*f0.^(pow-1).*exp(x)./(1+exp(x)).^2 + pow*(pow-1)*f0.^(pow-2).*(exp(x)./(1+exp(x))).^2;
end
end
function [mean_lambda, xNext, idxNext, var_lambda, predictiveMean] = makePrediction_nse(prs, datastruct, aFinal, WFinal, sqrtLFinal)

muf = prs(1);
param = prs(2:3);
nv = prs(4);

Kstar = datastruct.Kstar;

norm_mat_support = datastruct.norm_mat_support;

if datastruct.whichKernel == 2
    norm_mat_support2 = datastruct.norm_mat_support2;
    Kstarstar = param(1)*exp(-.5/param(2).*norm_mat_support) + param(1)*exp(-.5/param(2).*norm_mat_support2);
else
    Kstarstar = param(1)*exp(-.5/param(2).*norm_mat_support);
end

try
predictiveMean = muf + Kstar*aFinal; % of f
catch
    keyboard
end

M = WFinal\(sqrtLFinal*Kstar');

predictiveCov = Kstarstar - M'*M; % of f

predictiveVar = diag(predictiveCov);

%%

xx = length(predictiveMean);
mean_lambda = zeros(xx, 1);
var_lambda = zeros(xx, 1); 

for i = 1: xx
    [mean_lambda(i), var_lambda(i)] = meanvar_t_GaussHermiteQuadrature(nv, predictiveMean(i), predictiveVar(i),datastruct.g, prs(5));
end

%% sampling based on variance reduction

% how many same points do we tolerate?
tol_man = 3; 

if size(unique(datastruct.x(end:-1:end-tol_man,:), 'rows'), 1)==1
    
    % if the previous three stimuli are the same, we will
    % randomly jump to another stimulus.
    idxNext = floor(rand*length(datastruct.support))+1;
    xNext = datastruct.support(idxNext,:);
    
else
    
    idx = find(var_lambda == max(var_lambda));
    idxNext = idx(floor(rand*length(idx))+1);
    xNext = datastruct.support(idxNext,:);
    
end

end
function [mt, vt] = meanvar_t_GaussHermiteQuadrature(nsevar, muf, varf, g, p)

% compute the mean and variance of the estimate of t
% using Gauss-Hermite Quadrature, given nsevar, muf, varf, and g
% g = log(exp + 1)^p;


% 1. GHQ
n = 10;
xi = roots([1024 0 -23040 0 161280 0 -403200 0 302400 0 -30240]);
wi =  (2.^9*factorial(10)*sqrt(pi))./(10^2*(512*xi.^9 - 9216*xi.^7 + 48384*xi.^5 - 80640*xi.^3 + 30240*xi).^2);

% n = 5;
% xi = roots([32 0 -160 0 120 0]);
% wi = (2.^4*factorial(5)*sqrt(pi))./(25*(16*(xi.^4)-48.*(xi.^2)+12).^2);

% mean
mt = sum(wi.*g(muf+sqrt(2*(nsevar+varf))*xi, p))/sqrt(pi);

% var
fi = sqrt(2*varf).*xi + muf;
h = zeros(n,1);
for i=1:n
    h(i) = sum(wi.*g(fi(i)+sqrt(2*nsevar).*xi, p))/sqrt(pi);
end

vt = sum(wi.*(h.^2))/sqrt(pi) - mt.^2;
end
function [differece_in_h, hmap, a, W, sqrtL, L] = objFun_NewtonsMethod_nse(h0, datastruct)

H = datastruct.H;

r = datastruct.r;
p = datastruct.pest; 

[g,dg,ddg] = logexp1(h0, p);
diagL = r./(g.^2).*dg.^2 - r./g.*ddg + ddg;
diagL(diagL<0) = 1e-6; % to avoid numerical problem.
diagL(isnan(diagL)) = 1e-6; % in case g=0, diagL=nan.
L = diag(diagL);

muf = datastruct.muf;

sqrtL = sqrt(L);

B = eye(size(L)) + sqrtL*H*sqrtL;

W = chol(B, 'lower');

b = L*h0 + r./(g+1e-6).*dg - dg;

M = W\sqrtL;
sqrtLinvBsqrtL = M'*M;
a = b - sqrtLinvBsqrtL*(H*b + muf);

hmap = H*a + muf;

differece_in_h = norm(hmap - h0);


end
function [neglogev, prs, hmapFinal, aFinal, WFinal, sqrtLFinal, Hfinal, detH] = updateThetaGivenL_nse(prs0, datastruct)

% set bounds on estimated parameters
% these parameters are for f, i.e., log-lambda 
muf = [1e-3, 1e3]; % mean of f 
alpha = [1e-3, 1e3]; % overall scale of f
gamma = [1e-3, 1e3]; % smoothness scale of f
nsevar = [1e-3, 1e2]; % nsevar
p = [1, 5]; % p for log(exp+1).^p

LB = [muf(1); alpha(1); gamma(1); nsevar(1); p(1)];
UB = [muf(2); alpha(2); gamma(2); nsevar(2); p(2)]; 

lfun= @(pp)computeMarginal_usingAnalyticForm_nse(pp, datastruct);

opts = optimset('Display','iter','TolFun',1e-8, 'algorithm', 'active-set', 'maxIter', 3*1e3, 'MaxFunEvals', 5*1e3, 'TolX', 1e-8);

% ------ Optimize evidence --------------------------------------
[prs, fval, exitflag, output, lambda, grad, hessian] = fmincon(lfun,prs0,[],[],[],[],LB,UB,[],opts);

[neglogev, hmapFinal, aFinal, WFinal, sqrtLFinal, Hfinal] = lfun(prs);
detH = 1/det(hessian); 
end

