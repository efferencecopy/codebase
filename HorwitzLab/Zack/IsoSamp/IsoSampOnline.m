%#ok<*DEFNU>
function IsoSampOnline()
    global gl udpCom

    % assign gl.<monkey name> model params
    loadMonkeyData(load('IsoSamp/isosurf_params.mat'));

    gl.receivedHdr = false;
    gl.genStims = false;
    gl.monkey = 0;
    gl.model = zeros(10, 1);
    gl.plotlims = zeros(3, 1);
    gl.sf = 0;
    gl.plexRecOn = false;
    gl.stims = zeros(3000, 3);
    gl.nextStim = zeros(1, 3);
    gl.stimIdx = 1; % increments when a CORRECTCD drops

    udpCom.sock = [];
    udpCom.port = 6665;
    udpCom.rexip = '192.168.1.120';
    udpCom.plexip = '192.168.1.121';
    udpCom.mac = '192.168.1.122';

    all_codes = IsoSampCodes(); % defined in a separate m file

    p = InitPStruct(0, all_codes.HDRCOMPLETECD); % HDRCOMPLETECD dropped last

    socketOpen = 0;
    while ~socketOpen
        [udpCom.sock, socketOpen] = pnetStart(udpCom.port);
    end
    plxServer = InitPlex();

    bounceOut = false;
    while ~bounceOut
        if CheckForESCKey() || dealWithMsgs(udpCom.sock)
            bounceOut = true;
        end

        [n, eventList] = PL_GetTS(plxServer);
        if n
            p = ProcessEventList(p, eventList);
            if ~gl.receivedHdr && any(p.events == p.headercode)
%                 p = parseHeader(p, codes);
                p = RemoveOldEvents(p, p.headercode);
                gl.receivedHdr = true;
            end
            p = parseTrialCodes(all_codes, p);
        end
        p = CleanUpEvents(p);
    end
    PL_Close(plxServer);
end

% function p = parseHeader(p, codes)
% end

function loadMonkeyData(monkey_struct)
    global gl
    gl.monkeynames = fieldnames(monkey_struct);

    for s = gl.monkeynames'
        gl.(char(s)) = monkey_struct.(char(s));
    end
end

function p = parseTrialCodes(all_codes, p)
    global gl

    code_func_map = {
        all_codes.FPONCD @fFpon;
        all_codes.CORRECTCD @fCorrect;
        all_codes.EOTCD @fEot;
    };

    for codenum = 1:size(code_func_map,1)
        if any(p.events == code_func_map{codenum,1})
            feval(code_func_map{codenum,2});
            codeidx = find(p.events == code_func_map{codenum,1}, 1, 'last');
            p.lastprocessed_t = p.times(codeidx);
        end
    end

    function fFpon
        gl.plexRecOn = true;
    end
%     function fFpoff
%     end
    function fCorrect
        gl.stimIdx = gl.stimIdx + 1;
    end
%     function fAbort
%     end
    function fEot
        gl.plexRecOn = false;
    end
end

function num_stims = genStimSet(monkey, sf, target_nstims)
    global gl

    gl.monkey = monkey;
    gl.stimIdx = 1;

    num_stims = 0; % if zero is returned REX pauses

    try
        [gl.model, gl.plotlims, gl.keepratio] = ...
            setModelParams(monkey, sf, target_nstims);
    catch E
        if strcmp(E.identifier, 'MATLAB:badsubscript') || ...
                strcmp(E.identifier, 'MATLAB:nonExistentField')
            logger('Problem setting up the monkey''s parameters\n');
            return
        else
            rethrow(E);
        end
    end

    n_meshpts = 70;

    [x,y,z] = meshgrid( ...
        linspace(-gl.plotlims(1), gl.plotlims(1), n_meshpts), ...
        linspace(-gl.plotlims(2), gl.plotlims(2), n_meshpts), ...
        linspace(-gl.plotlims(3), gl.plotlims(3), n_meshpts));
    xyz = [x(:) y(:) z(:)];

    v = sum(abs(xyz * reshape(gl.model(2:end), [3 3])) .^ gl.model(1), 2);
    fv = isosurface(x, y, z, reshape(v, size(x)), 1);

    p0 = min(fv.vertices(:,1:3));
    p1 = max(fv.vertices(:,1:3));

    try
        [node,~,face] = surf2mesh(fv.vertices, fv.faces, p0, p1, ...
            gl.keepratio, 5, [], []);
    catch E
        if strcmp(E.identifier, 'MATLAB:unassignedOutputs') || ...
                strcmp(E.message, 'node file is missing!')
            logger('No stimuli generated\n');
            return
        else
            rethrow(E);
        end
    end

    num_stims = size(face, 1);

    logger('Generated %d stimuli!\n', num_stims);

    gl.stims(num_stims+1:end) = [];

    facerows = num2cell(face(:,1:3), 2);
    gl.stims = cell2mat(cellfun(@(x) mean(node(x,:), 1), ...
        facerows, 'uniformoutput', 0));
    gl.stims = gl.stims(randperm(size(gl.stims, 1)),:); % / 100 on slave
end

function [model, plotlims, keepratio] = setModelParams(monkeyidx, sf, target_nstims)
    global gl

    % find the model that corresponds to the requested spatial frequency
    errs = abs(sf - gl.(gl.monkeynames{monkeyidx}).sfs);
    sfidx = find(errs == min(errs), 1);

    % find the nearest keepratio that maps to the requested number of stims
    keepratios = linspace(.1, .0001, 100);
    errs = abs( ...
        target_nstims - gl.(gl.monkeynames{monkeyidx}).kr_stims_map(:,sfidx));

    keepratio = keepratios(errs == min(errs));
    keepratio = keepratio(1);
	model = gl.(gl.monkeynames{monkeyidx}).model(:,sfidx);
	plotlims = gl.(gl.monkeynames{monkeyidx}).plotlims(:,sfidx);
end

function allDone = dealWithMsgs(socket)
    allDone = false;
    msgSize = pnet(socket, 'readpacket', 250, 'noblock');
    if ~msgSize, return; end

    message = pnet(socket, 'read', msgSize, 'char');
    if ~isempty(message)
        [message,allDone] = parseMsg(message);
        evalMsg(message);
    else
        allDone = true;
    end
end

% The following function does 2 things. First, it handles the case when the
% Online program was asked to return but was called from the command
% line (not through MasterOnline). Second, if 'message' begins with
% 'sendToRex', the 2nd argument to the sendToRex function is parsed out.
% Once we have that string, we can do any necessary processing related to
% it before it is sent to Rex.
% This function may be unnecessary (apart from the 'return' check--you need
% that) so you can get rid of it.
function [message,doneFlag] = parseMsg(message)
    doneFlag = false;
    
    matches = strtrim(regexp(message, ',', 'split'));
    if strncmp(message, 'return', 6)
        stk = dbstack();
        if ~strcmp(stk(end).name, mfilename)
            doneFlag = true;
        end
    elseif matches{2}(1) ~= '[' % assuming '[a b c ...]' won't require more processing
        if strncmp(matches{1}, 'sendToRex', 9)
            % the regex below finds the proper variable or function name
            % (minus the parens and any function arguments).
            % 'thisFunction(23,51,2.2)' -> 'thisFunction'
            process_me = regexp(matches{2}, ...
                '[A-Za-z][.A-Za-z0-9_]*', 'match', 'once');
            if strcmp(process_me, 'gl.nextStim')
                setNextStim();
                % elseif strcmp(process_me, 'someExampleVariable')
            end
        end
    end
end

function setNextStim()
    global gl
    if gl.stimIdx > size(gl.stims, 1)
        gl.stims = gl.stims(randperm(size(gl.stims, 1)),:);
        gl.stimIdx = 1;
    end
    gl.nextStim = gl.stims(gl.stimIdx,:);
end

function evalMsg(message)
    global gl udpCom %#ok<NUSED> eval needs access to gl and udpCom
    try
        eval(message);
    catch exception
        fprintf('Trouble with message: "%s"\n', message);
        disp(getReport(exception));
    end
end
