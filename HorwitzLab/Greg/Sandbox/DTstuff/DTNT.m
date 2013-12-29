%%
% Contents
%
% 1) Reading a text file of filenames and organizing these filenames by spatial
% frequency.
%
% 2) Extracting detection thresholds from the files in the data structure
% generated in the previous section.
%
% 3) Displaying the data as thresholds in 3-D cone contrast space.
%
% 4) Finding the next set of color directions to test
%
% 5) Checking to see whether a common set of cone fundmantals was used for
% every file in a filelist (This checks out. GDLH 5/26/12)


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

[filelist,filelistpath,~] = uigetfile('N:\NexFiles\*.txt','Please select the file list');
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
[sfidx,okay] = listdlg('PromptString','Select a spatial frequency:',...
    'SelectionMode','single','ListString',cellstr(num2str([filenamesBySF{:,2}]')));
if ~okay, return; end

clear trialspecs;
load ('T_cones_smj10');
load ('T_cones_smj');

% Loading a list of data files, one by one
assumedDataPath = '';
filenames = filenamesBySF{sfidx,1};
for fileidx = 1:size(filenames,1)
    if fileidx == 1
        firstFile = findfile(char(filenames{fileidx}));
        if isempty(firstFile)
            warning('DTNT:fnf','File ''%s'' wasn''t found',char(filenames{fileidx}));
            continue
        else
            [assumedDataPath,~,~] = fileparts(firstFile);
        end
        currFile = firstFile;
    else
        nextFile = [assumedDataPath filesep char(filenames{fileidx})];
        if ~exist(nextFile,'file')
            nextFile = findfile(char(filenames{fileidx}));
            if isempty(nextFile)
                warning('DTNT:fnf','File ''%s'' wasn''t found',char(filenames{fileidx}));
                continue
            else
                [assumedDataPath,~,~] = fileparts(nextFile);
            end
        end
        currFile = nextFile;
    end
    stro = nex2stro(currFile);
    [thresholds, colorDirs, sfs, QuestTrajectories] = DTquestUnpackGH(stro, 'mode');
    figure; set(gcf,'name',char(filenames{fileidx}));
    for i = 1:length(thresholds)
        subplot(2,2,i);
        plot(QuestTrajectories{i},'k.');
        title(num2str(colorDirs(i,:)./norm(colorDirs(i,:))));
    end

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
            trialspecs(i).measuredthreshold = threshs(i);
            trialspecs(i).M = M;
            trialspecs(i).bkgndrgb = [stro.sum.exptParams.bkgnd_r stro.sum.exptParams.bkgnd_g stro.sum.exptParams.bkgnd_b];
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
            trialspecs(idx).measuredthreshold = [];
            trialspecs(idx).M = nan*ones(3);
            trialspecs(i).bkgndrgb = [nan nan nan];
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
            threshratio = trialspecs(L).measuredthreshold./trialspecs(L).predictedthreshold;
            if (abs(log(threshratio)) > abs(log(1+LINPREDTOL)))
                grandparentvertices = trialspecs(L).parentvertices;
                for j = 1:3
                    idx = length(trialspecs)+1;
                    parentvertices = grandparentvertices;
                    parentvertices(j,:) = trialspecs(L).colordir.*trialspecs(L).measuredthreshold;
                    v = mean(parentvertices);
                    trialspecs(idx).colordir = v./norm(v);
                    trialspecs(idx).predictedthreshold = norm(v);
                    trialspecs(idx).parentvertices = parentvertices;
                    trialspecs(idx).measuredthreshold = [];
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
% Displaying the data acqured so far
SHOWDONETRIANGLES = 0;
SHOWELLIPSOID = 0;

scaled = [];
for i = 1:length(trialspecs)
    if (~isempty(trialspecs(i).measuredthreshold))
        scaled = [scaled; trialspecs(i).colordir.*trialspecs(i).measuredthreshold];
    end
end

figure('Position',[140 250 400 400]);
axes('Position',[.25 .25 .5 .5]); hold on;
cmap = jet(size(scaled,1));
h = [];
for i = 1:size(scaled,1)
    h(1) = plot3(scaled(i,1),scaled(i,2),scaled(i,3),'o');
    h(2) = plot3(-scaled(i,1),-scaled(i,2),-scaled(i,3),'o');
    set(h,'MarkerFaceColor',cmap(i,:),'MarkerEdgeColor',cmap(i,:));
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

if sfs < 3
    set(gca,'XLim',[-3 3],'YLim',[-3 3],'Zlim',[-8 8],'view',[140 20])
else
    set(gca,'XLim',[-6 6],'YLim',[-6 6],'Zlim',[-30 30],'view',[140 20])
end

%%
% Section 4
% Finding the next set of color directions
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
%%
% Section 5
% Running through all the files in a filelist and making sure that the 
% set of cone fundamentals used is consistent. Changes in the calibration
% structure are fine but changes in the cone fundamentals used are not.
% Everything checked out fine. GDLH 5/26/12

load ('T_cones_smj10');
load ('T_cones_smj');

fnames = [];
for i = 1:size(filenamesBySF,1)
    fnames = [fnames; filenamesBySF{i,1}];
end

whichfunds = [];
for i = 1:size(fnames,1)
    stro = nex2stro(findfile(fnames(i,:)));
    mon_spd = stro.sum.exptParams.mon_spect;
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    M = reshape(stro.sum.exptParams.m_mtx,[3 3]);
    M10 = T_cones_smj10*mon_spd;
    M2 = T_cones_smj*mon_spd;
    bkgndrgb = [stro.sum.exptParams.bkgnd_r stro.sum.exptParams.bkgnd_g stro.sum.exptParams.bkgnd_b];
    if (all(M(:) == M10(:)))
        whichfunds(i) = 10;
    elseif (all(M(:) == M2(:)))
        whichfunds(i) = 2;
    else
         error('Not sure which fundamentals we''re using!');
    end
end
if (all(whichfunds == whichfunds(1)))
    disp([num2str(whichfunds(1)),' deg funds used in all experiments']);
else
    error('Different fundamentals used in different experiments');
end

%%
% Section 6
% Plotting figures that can be compared directly with Cole et al. 1993
% Dependent on sections 1 and 2, above.

scaled = [];
for i = 1:length(trialspecs)
    if (~isempty(trialspecs(i).measuredthreshold))
        scaled = [scaled; trialspecs(i).colordir.*trialspecs(i).measuredthreshold];
    end
end

figure;
% L/M plane
subplot(2,2,1); hold on;
plot(scaled(:,1)/100,scaled(:,2)/100,'k.')
plot(-scaled(:,1)/100,-scaled(:,2)/100,'k.')
axis square; xlabel('\DeltaL/L'); ylabel('\DeltaM/M');
set(gca,'Xlim',max(scaled(:,1)/100)*[-1.1 1.1],'Ylim',max(scaled(:,2)/100)*[-1.1 1.1])

% L-M/S plane
subplot(2,2,2); hold on;
lvmcontrast = scaled(:,[1 2])*[1/sqrt(2) -1/sqrt(2)]'/100;
plot(lvmcontrast,scaled(:,3)/100,'k.')
plot(-lvmcontrast,-scaled(:,3)/100,'k.')
axis square; xlabel('(-.71, .71, 0)'); ylabel('\DeltaS/S');
set(gca,'Xlim',max(lvmcontrast)*[-1.1 1.1],'Ylim',max(scaled(:,3)/100)*[-1.1 1.1])

% L+M/S plane
subplot(2,2,3); hold on;
lmcontrast = scaled(:,[1 2])*[1/sqrt(2) 1/sqrt(2)]'/100;
plot(lmcontrast,scaled(:,3)/100,'k.')
plot(-lmcontrast,-scaled(:,3)/100,'k.')
axis square; xlabel('(.71, .71, 0)'); ylabel('\DeltaS/S');
set(gca,'Xlim',max(lvmcontrast)*[-1.1 1.1],'Ylim',max(scaled(:,3)/100)*[-1.1 1.1])

%%
% Section 7
% Fitting surfaces (a la Cole et al 1993) to the data from each spatial
% frequency.  Using LS fit of ellipsoid as initial guess. Can we do better?
% Dependent on sections 1 and 2, above.

scaled = [];
for i = 1:length(trialspecs)
    if (~isempty(trialspecs(i).measuredthreshold))
        scaled = [scaled; trialspecs(i).colordir.*trialspecs(i).measuredthreshold];
    end
end

% Doing the fitting
scaled_dup = [scaled; -scaled];
% THIS COULD BE IMPROVED BY CONVERTING TO POLAR COORDINATES AND RESCALING r to LOG(r) 
% Using an ellipsoid fit to start
D = [scaled_dup(:,1) .* scaled_dup(:,1),...
    scaled_dup(:,2) .* scaled_dup(:,2),...
    scaled_dup(:,3) .* scaled_dup(:,3),...
    2*scaled_dup(:,1) .* scaled_dup(:,2),...
    2*scaled_dup(:,1) .* scaled_dup(:,3),...
    2*scaled_dup(:,2) .* scaled_dup(:,3)];
lssoln = (D' * D) \(D' * ones(size(scaled_dup,1),1));
A = [lssoln(1) lssoln(4) lssoln(5);...
    lssoln(4) lssoln(2) lssoln(6);...
    lssoln(5) lssoln(6) lssoln(3)];
[evecs, evals] = eig(A);
initparams = [2; reshape(evecs*sqrt(evals),9,1)];
options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-5,'Display','none');
[fpar,fv] = fminsearch(@(x) colefiterr(x,scaled_dup,0),initparams, options);

% Plotting the datapoints
figure; axes; hold on;
plot3(scaled(:,1),scaled(:,2),scaled(:,3),'k.')
plot3(-scaled(:,1),-scaled(:,2),-scaled(:,3),'k.')

% Plotting the fit
plotlims = max(abs(scaled));
[x,y,z] = meshgrid(linspace(-plotlims(1),plotlims(1),150),linspace(-plotlims(2),plotlims(2),150),linspace(-plotlims(3),plotlims(3),50));
tmp = [x(:) y(:) z(:)];
%   v = sum(abs(tmp *reshape(initparams(2:end),3,3)).^initparams(1),2);
v = sum(abs(tmp *reshape(fpar(2:end),3,3)).^fpar(1),2);

fv = isosurface(x,y,z,reshape(v,size(x,1),size(x,2),size(x,3)),1);
h = patch(fv,'EdgeAlpha',0,'FaceAlpha',.2);
plot3(scaled(:,1),scaled(:,2),scaled(:,3),'k.');
xlabel('\DeltaL/L'); ylabel('\DeltaM/M'); zlabel('\DeltaS/S');