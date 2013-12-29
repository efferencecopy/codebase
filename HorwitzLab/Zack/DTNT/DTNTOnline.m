function DTNTOnline
global gl

MAXTRIALSPECS = 500;

udpCom.sock = [];
udpCom.port = 6665;
udpCom.mac = '192.168.1.122';
udpCom.rexip = '192.168.1.120';
udpCom.plexip = '192.168.1.121';

stim.color_dirs = zeros(1,9);
stim.threshold_guesses = zeros(1,3);

socketOpen = 0;
while ~socketOpen
    [udpCom.sock, socketOpen] = pnetStart(udpCom.port);
end

[filelist,filelistpath,~] = uigetfile('N:\NexFiles\nexfilelists\Greg\DTNT\*.txt', ...
    'Please select the file list or press Cancel to start a new experiment');
if ~isequal(filelist,0)
    parsed_lines = open_parse_filenames([filelistpath filesep filelist]);
    
    if numel(parsed_lines) < 2
        error('DTNTOnline:TooFewLines', ...
            'The text file must contain at least 2 lines (''sf:#'' header and filename)');
    end
    
    % get the sfs and make the sf (column 1) -> filenames (column 2) map
    [sf_header_idx,listed_sfs] = lines_with_sf_header(parsed_lines);
    sf_filename_map = cell(length(listed_sfs), 2);
    sf_filename_map(:,1) = num2cell(listed_sfs);
    
    % Get the indices of the filenames belonging to every 'sf:' and stick
    % them in the map
    file_starts = sf_header_idx(:)' + 1;
    file_stops = [file_starts(2:end)-2 length(parsed_lines)];
    
    for i = 1:length(listed_sfs)
        sf_filename_map{i,2} = parsed_lines(file_starts(i):file_stops(i));
    end
    
    % Get rid of the invalid NEX filenames in the map (prints invalid ones
    % to command window)
    valid_nexfiles = cellfun(@isvalidnexfilename, sf_filename_map(:,2), 'unif', 0);
    sf_filename_map(:,2) = cellfun(@(x,y) {x(y)}, sf_filename_map(:,2), valid_nexfiles);
    sf_filename_map(cellfun(@isempty, sf_filename_map(:,2)),:) = [];
    
    % Consolidate all filenames of the same sf and sort the map in
    % ascending spatial frequency
    sf_filename_map = consolidate_sort_map(sf_filename_map, listed_sfs);
    
    if size(sf_filename_map, 1) > 1
        [sfidx,okay] = listdlg('PromptString', 'Select a spatial frequency:', ...
            'SelectionMode', 'single', 'ListString', cellstr(num2str([sf_filename_map{:,1}]')));
        if ~okay, return; end
    else
        sfidx = 1;
    end
    
    if isfield(gl, 'dtnt') && ~isempty(gl.dtnt.filenames) ...
            && gl.dtnt.nfiles <= numel(sf_filename_map{sfidx,2}) ...
            && all(strcmp(gl.dtnt.filenames, sf_filename_map{sfidx,2}(1:gl.dtnt.nfiles)))
        gl.dtnt.filenames = sf_filename_map{sfidx,2};

        for filename = gl.dtnt.filenames(gl.dtnt.nfiles+1:end)'
            if size(gl.dtnt.ts, 1) - gl.dtnt.nspecs <= 20
                gl.dtnt.ts = resize_trialspecs_byfactor(2); % resize trialspecs by 2x
            end
            
            stro = nex2stro(findfile(filename{:}));
            [threshs, color_dirs, ~, trajectories] = DTquestUnpackGH(stro, 'mode');            
            gl.dtnt.nfiles = gl.dtnt.nfiles + 1;
            
            setup_next_rounds(threshs, color_dirs, trajectories, filename{:});
        end
    else % either no files have been loaded, or the previously loaded files don't match what was just parsed
        % reload everything
        gl.dtnt.ts = init_dtnt_struct(MAXTRIALSPECS);
        gl.dtnt.filenames = sf_filename_map{sfidx,2};
        
        for filename = gl.dtnt.filenames'
            stro = nex2stro(findfile(filename{:}));
            [threshs, color_dirs, ~, trajectories] = DTquestUnpackGH(stro, 'mode');
            
            gl.dtnt.nfiles = gl.dtnt.nfiles + 1;
            
            if gl.dtnt.nfiles == 1
                setup_init_rounds(threshs, color_dirs, trajectories);
            else
                setup_next_rounds(threshs, color_dirs, trajectories, filename{:});
            end
        end
    end
    
    INITCONTRASTFACTOR = 4; % multiplier for initial contrast (not starting in the plane of
    % the parent triangle, otherwise contrast would be too low).

    candidates = cellfun(@isempty, {gl.dtnt.ts(1:gl.dtnt.nspecs).measuredthreshold});
    fprintf('\n%d color directions currently in the queue\n', sum(candidates));

    idxs = find(candidates);
    if sum(candidates) > 3  % can show maximum three color directions at a time
        idxs = idxs(floor((0:2)*(length(idxs)/3))+1);
    end
    
    % Need to pare the list down to three entries
    % Printing out the color directions and thresholds to use
    for idx = idxs
        fprintf('LMS = (% .6f,% .6f,% .6f) SCALE = %.6f\n', gl.dtnt.ts(idx).colordir, ...
            INITCONTRASTFACTOR * gl.dtnt.ts(idx).predictedthreshold);
    end

    stim.color_dirs = [gl.dtnt.ts(idxs).colordir zeros(1,3*(3-length(idxs)))];
    stim.threshold_guesses = INITCONTRASTFACTOR*[gl.dtnt.ts(idxs).predictedthreshold zeros(1,3-length(idxs))]; %#ok<STRNU>
else % This is the first time through (user clicked "cancel" on the UI)
    gl.dtnt.ts = init_dtnt_struct(MAXTRIALSPECS);
    stim.color_dirs = [1 1 1 1 -1 0 0 0 1];
    stim.threshold_guesses = [6 20 80];
    disp('Using default color directions and threshold guesses')
    disp(stim.color_dirs);
    disp(stim.threshold_guesses);
end

% MAIN LOOP
bounceOut = false;
while ~bounceOut
    if CheckForESCKey() || dealWithMsgs(udpCom.sock)
        bounceOut = true;
    end
end

    % nested functions
    function allDone = dealWithMsgs(socket)
        allDone = false;
        msgSize = pnet(socket, 'readpacket', 250, 'noblock');
        if ~msgSize, return; end
        
        message = pnet(socket, 'read', msgSize, 'char');
        if ~isempty(message)
            evalMsg(message);
        else
            allDone = true;
        end
    end

    function evalMsg(message)
        try
            eval(message);
        catch exception
            fprintf('Trouble with message: "%s"\n', message);
            disp(getReport(exception));
        end
    end
end
