function DTNTOnline
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

[filelist,filelistpath,~] = uigetfile('N:\NexFiles\nexfilelists\Greg\DTNT\*.txt','Please select the file list');
if ~isequal(filelist,0)
    % parse every line while ignoring comments
    fid = fopen([filelistpath filesep filelist]);
    listData = textscan(fid,'%s','CommentStyle','%','Delimiter','\n');
    % NOTE: textscan will ignore lines with a comment character (even if it's after the filename)!
    fclose(fid);
    listData = strtrim(listData{1});
    
    badEntries = cellfun('isempty',...
        regexp(listData,'(^[A-Z]{1,2}(0[1-9]|1[012])(0[1-9]|[12][0-9]|3[01])\d{5}(\.\d)?(\.nex)?$)|(^sf:)'));
    
    blankLines = cellfun('length',listData) == 0;
    if sum(badEntries & ~blankLines) ~= 0
        warning('The following lines will be ignored:');
        disp(listData(badEntries & ~blankLines));
    end
    listData(badEntries | blankLines) = [];
    
    sfHeadingIdxs = find(~cellfun('isempty',strfind(listData,'sf:')));
    if ~isempty(sfHeadingIdxs)
        listSFs = textscan([listData{sfHeadingIdxs}],'sf:%f');
        listSFs = listSFs{1};
        
        % sort file names by SF
        filenamesBySF = cell(length(listSFs),2);
        for i = 1:length(listSFs)
            if i ~= length(listSFs)
                fns = listData(sfHeadingIdxs(i)+1:sfHeadingIdxs(i+1)-1);
            else
                fns = listData(sfHeadingIdxs(i)+1:end);
            end
            filenamesBySF{i,1} = fns;
            filenamesBySF{i,2} = listSFs(i);
        end
    end
    % sort list by sf
    sfNums = cell2mat(filenamesBySF(:,2));
    [sfNums,sortedOrder] = sort(sfNums);
    filenamesBySF = filenamesBySF(sortedOrder,:);
    
    % find the indices of the dupes (based on repval from Mathworks FEX)
    [B,~,J] = unique(sfNums);
    I = 1:length(J);
    I2 = find(diff(J) == 0);
    B2 = unique([I2(:)' I2(:)'+1]);
    [~,~,IR] = unique(J(B2));
    if ~isempty(IR)
        POS = I(B2);
        idxsToDelete = [];
        for i = unique(IR)'
            idxsDupes = POS(IR == i);
            temp = filenamesBySF(idxsDupes,1);
            filenamesBySF{idxsDupes(1),1} = vertcat(temp{:});
            idxsToDelete = [idxsToDelete; idxsDupes(2:end)];
        end
        filenamesBySF(idxsToDelete,:) = []; % trim them out
    end
    
    LINPREDTOL = 0.3;
    INITCONTRASTFACTOR = 4; % multiplier for initial contrast (not starting in the plane of
    % the parent triangle, otherwise contrast would be too low).
    if size(filenamesBySF, 1) > 1
        [sfidx,okay] = listdlg('PromptString','Select a spatial frequency:',...
            'SelectionMode','single','ListString',cellstr(num2str([filenamesBySF{:,2}]')));
        if ~okay, return; end
    else
        sfidx = 1;
    end
    
    % Loading a list of data files, one by one
%     assumedDataPath = '';
    filenames = filenamesBySF{sfidx,1};
    for fileidx = 1:size(filenames,1)
%         if fileidx == 1
%             firstFile = findfile(char(filenames{fileidx}));
%             if isempty(firstFile)
%                 warning('DTNT:fnf','File ''%s'' wasn''t found',char(filenames{fileidx}));
%                 continue
%             else
%                 [assumedDataPath,~,~] = fileparts(firstFile);
%             end
%             currFile = firstFile;
%         else
%             nextFile = [assumedDataPath filesep char(filenames{fileidx})];
%             if ~exist(nextFile,'file')
%                 nextFile = findfile(char(filenames{fileidx}));
%                 if isempty(nextFile)
%                     warning('DTNT:fnf','File ''%s'' wasn''t found',char(filenames{fileidx}));
%                     continue
%                 else
%                     [assumedDataPath,~,~] = fileparts(nextFile);
%                 end
%             end
%             currFile = nextFile;
%         end
        stro = nex2stro(findfile(filenames{fileidx}));
        [thresholds, colorDirs, sfs, QuestTrajectories] = DTquestUnpackGH(stro, 'mode');
        %     figure; set(gcf,'name',char(filenames{fileidx}));
        %     for i = 1:length(thresholds)
        %         subplot(2,2,i);
        %         plot(QuestTrajectories{i},'k.');
        %         title(num2str(colorDirs(i,:)./norm(colorDirs(i,:))));
        %     end
        
        threshs = thresholds';
        lmsmat = mkbasis(colorDirs')';
        
        mon_spd = stro.sum.exptParams.mon_spect;
        mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
        M = reshape(stro.sum.exptParams.m_mtx,[3 3]);
        
        % Each row is a color direction for consistency with NeuroThreshOnline.m
        if (fileidx == 1)
            % Making entries into the trialspecs array: round 1, (3 color dirs)
            for i = 1:size(lmsmat,1)
                trialspecs(i).colordir = lmsmat(i,:)./norm(lmsmat(i,:));
                trialspecs(i).predictedthreshold = [];
                trialspecs(i).parentvertices = [];
                trialspecs(i).parentOOGs = [];
                trialspecs(i).measuredthreshold = threshs(i);
                if (threshs(i) == max(QuestTrajectories{i}))
                    trialspecs(i).OOG = 1;
                else
                    trialspecs(i).OOG = 0;
                end
                % trialspecs(i).M = M;
                % trialspecs(i).bkgndrgb = [stro.sum.exptParams.bkgnd_r stro.sum.exptParams.bkgnd_g stro.sum.exptParams.bkgnd_b];
            end
            
            % Making entries into the trialspecs array: round 2, (4 color dirs)
            signmat = [(fullfact([2,2])-1.5)*2 ones(4,1)];  % assumes three color dirs
            for i = 1:size(signmat,1)
                idx = length(trialspecs)+1;
                parentvertices = diag(signmat(i,:).*threshs)*lmsmat;
                v = mean(parentvertices);
                trialspecs(idx).colordir = v./norm(v);
                trialspecs(idx).predictedthreshold = norm(v);
                trialspecs(idx).parentvertices = parentvertices;
                trialspecs(idx).parentOOGs =  [trialspecs(1:size(lmsmat,1)).OOG];
                trialspecs(idx).measuredthreshold = [];
                trialspecs(idx).M = nan*ones(3);
            end
        else  % end first data file (setting up round 1 and 2)
            
            colordirs = reshape([trialspecs.colordir],3,length(trialspecs))';
            for i = 1:size(lmsmat,1)
                L = logical(softEq((lmsmat(i,:)./norm(lmsmat(i,:)))*colordirs', 1,10));
                % Error checking
                if (sum(L) == 0)
                    disp(['cannot find color direction in trialspecs: ',num2str(lmsmat(i,:))]);
                    continue
                end
                if (sum(L) > 1)
                    error('Multiple identical color directions detected.');
                    keyboard;
                end
                if (~isempty(trialspecs(L).measuredthreshold))
                    error('Already measured a threshold in this direction');
                end
                trialspecs(L).measuredthreshold = threshs(i);
                if (threshs(i) == max(QuestTrajectories{i}))
                    trialspecs(L).OOG = 1;
                else
                    trialspecs(L).OOG = 0;
                end
                threshratio = trialspecs(L).measuredthreshold./trialspecs(L).predictedthreshold;
                if (abs(log(threshratio)) > abs(log(1+LINPREDTOL)))  % making new nodes (directions to search)
                    grandparentvertices = trialspecs(L).parentvertices;
                    grandparentOOGs = trialspecs(L).parentOOGs;
                    
                    for j = 1:3
                        parentOOGs = grandparentOOGs;
                        parentOOGs(j) = trialspecs(L).OOG;
                        if (~all(parentOOGs))
                            idx = length(trialspecs)+1;
                            parentvertices = grandparentvertices;
                            parentvertices(j,:) = trialspecs(L).colordir.*trialspecs(L).measuredthreshold;
                            v = mean(parentvertices);
                            trialspecs(idx).colordir = v./norm(v);
                            trialspecs(idx).predictedthreshold = norm(v);
                            trialspecs(idx).parentvertices = parentvertices;
                            trialspecs(idx).parentOOGs = parentOOGs;
                            trialspecs(idx).measuredthreshold = [];
                        end
                    end
                else
                    disp(['File #',num2str(fileidx),': Threshold in direction ',num2str(lmsmat(i,:)),' is consistent with linear prediction']);
                end
                if (abs(log(threshratio)) > log(3))
                    fprintf('Warning: aberrant threshold point: %s\n',[char(filenames{fileidx}),' (',num2str(trialspecs(L).colordir),')']);
                end
            end
        end
    end % Datafile loop
    candidates = zeros(length(trialspecs),1);
    for i = 1:length(trialspecs)
        if (isempty(trialspecs(i).measuredthreshold))
            candidates(i) = 1;
        end
    end
    disp(['Number of color directions currently in queue: ',num2str(sum(candidates))]);
    idxs = find(candidates);
    if (sum(candidates) > 3)  % can show maximum three color directions at a time
        idxs = idxs(floor([0:2]*(length(idxs)/3))+1);
    end
    % Need to pare the list down to three entries
    % Printing out the color directions and thresholds to use
    for i = idxs'
        disp(['LMS = ',num2str(trialspecs(i).colordir),' SCALE = ',num2str(INITCONTRASTFACTOR*trialspecs(i).predictedthreshold)]);
    end
    stim.color_dirs = [trialspecs(idxs).colordir];
    stim.threshold_guesses = INITCONTRASTFACTOR*[trialspecs(idxs).predictedthreshold];
else % This is the first time through (user clicked "cancel" on the UI)
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