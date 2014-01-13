% Makes and plots to one figure for every two function calls
function DTNTdualplot
persistent times_plotted
if isempty(times_plotted), times_plotted = 0; end
%%
% Section 1
% Reading a text file of DTNT filenames and creating an Nx2 matrix.
% The first element in each row is a cell array, each element of which
% is a filename. The second element in each row is just a scalar which is
% the spatial frequency used for those data files. This is a convenient way
% to hold onto the data since many color directions were tested but only a 
% a few SFs were tested.
%
% This script is smart in the sense that it finds entries in the textfile
% that begin with 'sf:' and uses whatever comes next to group the files.
if ispc
    filelistpath = 'N:\NexFiles\nexfilelists\Greg\DTNT';
else
    filelistpath = '/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/DTNT';
end

[filelist,filelistpath,~] = uigetfile([filelistpath filesep '*.txt'], 'Please select the DTNT file list');
if isequal(filelist,0), return; end

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
    disp('Warning -- the following lines will be ignored:');
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

%%
% Section 2
% Extracting thresholds from the datafiles collected with a user-entered
% spatial frequency. Creates a "trialspecs" data structure, just like in
% NeuroThreshOnline, that keeps track of which color directions have been
% tested and which are in the queue, ready to be tested. This structure
% allows us to check for errors in the color directions that were tested
% and also lets us see how many (and which) color directions matched the
% linear prediction.

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

clear trialspecs;
load ('T_cones_smj10');
load ('T_cones_smj');

% Loading a list of data files, one by one
filenames = filenamesBySF{sfidx,1};
for fileidx = 1:size(filenames,1)
    stro = nex2stro(findfile(filenames{fileidx}));
    [thresholds, colorDirs, sfs, QuestTrajectories] = DTquestUnpackGH(stro, 'mode');

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
            trialspecs(i).parentOOG = [];
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
           % trialspecs(i).bkgndrgb = [nan nan nan];
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
            if (threshs(i) == QuestTrajectories{i}(1))
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
%%
% Section 3
% Displaying the data acquired so far
SHOWDONETRIANGLES = 0;
SHOWELLIPSOID = 0;

scaled = []; LOOG = [];
for i = 1:length(trialspecs)
    if (~isempty(trialspecs(i).measuredthreshold))
        scaled = [scaled; trialspecs(i).colordir.*trialspecs(i).measuredthreshold];
        LOOG = logical([LOOG; trialspecs(i).OOG]);
    end
end

figure(999+floor(times_plotted/2)); hold on;
if ~mod(times_plotted,2)
    point_color = 'k';
else
    point_color = 'b';
end
h = [];
for i = 1:size(scaled,1)
    if (~LOOG(i))
        plot3(scaled(i,1),scaled(i,2),scaled(i,3),[point_color 'o'],'markerfacecolor',point_color);
        plot3(-scaled(i,1),-scaled(i,2),-scaled(i,3),[point_color 'o'],'markerfacecolor',point_color);
    else
        h(1) = plot3(scaled(i,1)*[-1 1],scaled(i,2)*[-1 1],scaled(i,3)*[-1 1],'r-');
    end
end
xlabel('L contrast');
ylabel('M contrast');
zlabel('S contrast');

% Plotting 'done' triangles
if (SHOWDONETRIANGLES)
    for i = 1:length(trialspecs)
        if (~isempty(trialspecs(i).measuredthreshold))
            threshratio = trialspecs(i).measuredthreshold./trialspecs(i).predictedthreshold;
            if (abs(log(threshratio)) <= abs(log(1+LINPREDTOL)))
                h(1) = patch(trialspecs(i).parentvertices(:,1),trialspecs(i).parentvertices(:,2),trialspecs(i).parentvertices(:,3),'green');
                h(2) = patch(-trialspecs(i).parentvertices(:,1),-trialspecs(i).parentvertices(:,2),-trialspecs(i).parentvertices(:,3),'green');
                set(h,'FaceAlpha',.2);
            end
        end
    end
end

if (SHOWELLIPSOID)
    [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, zeros(size(scaled,1),1))
    A = [quadparams(1) quadparams(4) quadparams(5);...
        quadparams(4) quadparams(2) quadparams(6);...
        quadparams(5) quadparams(6) quadparams(3)];   [evecs, evals] = eig(A);
    radii = sqrt(1./diag(evals));
    [x, y, z] = ellipsoid(0,0,0,radii(1),radii(2), radii(3));
    [newxyz] = [x(:) y(:) z(:)]*evecs'/xformmat;
    
    h = surf(reshape(newxyz(:,1),size(x)), reshape(newxyz(:,2),size(y)), reshape(newxyz(:,3),size(z)));
    set(h,'FaceAlpha',.2,'FaceColor','green','Edgealpha',0);
    camlight;
    lighting phong;
end

times_plotted = times_plotted + 1;
