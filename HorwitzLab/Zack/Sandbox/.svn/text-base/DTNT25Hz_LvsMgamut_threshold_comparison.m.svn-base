% This is a hacked together script to answer the question (in the negative
% hopefully) "Is it still possible that at 25 Hz monkeys' are equally
% sensitive to L+M and L-M?"
%
% The purpose of this script is to superimpose the 25 Hz detection thresholds
% over the L-M gamut sample points
%%
% Section 1
if ispc
    filelistpath = 'C:\NO BACKUP\NexFiles\nexfilelists\Greg\DTNT';
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

% extract valid, in-gamut L+M thresholds
colordirs = reshape([trialspecs.colordir],3,[])';
has_thresh = ~cellfun(@isempty, {trialspecs.measuredthreshold})';
trialspecs = trialspecs(has_thresh);
colordirs = colordirs(has_thresh,:);
scaled = bsxfun(@times,colordirs,[trialspecs.measuredthreshold]');
LOOG = [trialspecs.OOG]';

%%
% Section 4
% drawing the extent of the gamut in L-M

% generate color dirs
% define the color directions.... just putting points on a sphere for now
nColors = 2000;
tmp = ceil(sqrt(nColors));
az = linspace(0, 2*pi, tmp);
el = linspace(0, pi/2, tmp);
inds = fullfact([tmp, tmp]);
dirs_pol = [az(inds(:,1))', el(inds(:,2))'];
[x,y,z] = sph2cart(dirs_pol(:,1), dirs_pol(:,2), ones(size(dirs_pol,1),1));
%z = z.*3; % scale the s-cone axis
colorDirs = [x,y,z];
norms = sqrt(sum(colorDirs.^2,2));
colorDirs = bsxfun(@rdivide, colorDirs, norms);
colorDirs(abs(colorDirs)<1e-14) = 0;
colorDirs = [colorDirs; 1 -1 0; -1 1 0; 0 0 -1; 0 0 1];
colorDirs = unique(colorDirs, 'rows'); % remove the duplicates
nColors = size(colorDirs,1);

lvsm_idxs = find(sign(colorDirs(:,1)) .* sign(colorDirs(:,2)) < 0);

load Dell4BitsCal.mat
s = load('T_cones_smj10');
fns = fieldnames(s);
fundamentals = s.(fns{1});
fundWavelengthSpacing = s.(fns{2});
calData = cals{end};
ptb.bkgndRGB = round(255 .* calData.bgColor); %bkgnd voltages discretized b/w 0&255
ptb.bkgndrgb = [calData.gammaTable(ptb.bkgndRGB(1)+1, 1), calData.gammaTable(ptb.bkgndRGB(2)+1, 2), ...
    calData.gammaTable(ptb.bkgndRGB(3)+1, 3)]; %add one to create and index b/w 1 & 256

ptb.maxdac = 2^16-1;
ptb.gammaTable = calData.gammaTable;
ptb.monSpd = SplineSpd(calData.S_device, calData.P_device, fundWavelengthSpacing);
ptb.M = fundamentals * ptb.monSpd;
ptb.invM = inv(ptb.M);
ptb.bkgndlms = ptb.M * ptb.bkgndrgb';
ptb.invGamma = InvertGammaTable(calData.gammaInput, ptb.gammaTable, ptb.maxdac+1);

for c = 1:size(colorDirs, 1);
    guess = 0;
    incSize = 0.1; %1 tenth of a %CC
    while(1)
        tmpLMS = colorDirs(c,:) .* guess;
        
        %peak of the gabor:
        peakrgb = ptb.invM*((1+tmpLMS./100)' .* ptb.bkgndlms);
        peakRGB = round(ptb.maxdac .* peakrgb) + 1;
        peakInGamut = all(peakRGB < (ptb.maxdac+1)) && all(peakRGB > 0);
        
        %trough of the gabor:
        troughrgb = ptb.invM*((1-tmpLMS./100)' .* ptb.bkgndlms);
        troughRGB = round(ptb.maxdac .* troughrgb) + 1;
        troughInGamut = all(troughRGB < (ptb.maxdac+1)) && all(troughRGB > 0);
        
        if (peakInGamut && troughInGamut)
            guess = guess + incSize;
        else
            maxContrast(c) = guess - incSize;
            break
        end
    end
end

gamut_in_cc = bsxfun(@times,colorDirs,maxContrast');
lvsm_gamut = gamut_in_cc(lvsm_idxs,1:2);

%% plotting
figure();
ax_main = axes(); hold on;
set(gca,'XLim',[-50 50],'YLim',[-50 50],'View',[0 90]);

scatter(ax_main, lvsm_gamut(:,1),lvsm_gamut(:,2), 15, 'fill', 'markerface', 'r', 'markeredge', 'r');
scatter([scaled(~LOOG,1); -scaled(~LOOG,1)],[scaled(~LOOG,2); -scaled(~LOOG,2)], 25, ...
    'fill', 'markerface', 'k', 'markeredge', 'k');

title(filelist(1:end-12)); legend('L-M gamut','thresholds'); legend boxoff;
axis square
