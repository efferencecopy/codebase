% A few population analysis for DT spot data
%
% Section 1: Looking at superimposed Quest functions for a variety of runs.
%
% Section 2: Comparing psychometric functions for T1 and T2 presentations
% (is the monkey better at detecting stimuli on one side over the other?)
%
% Section 3: Finding (Quest) experiments with unusually high thresholds.
%
% Section 4: Reading in a list of DTNT files (Excel?) from multiple observers (also
% different display devices) and displaying them (or fits to them).
% This will be useful for comparing humans to monkeys or for comparing
% ViewPixx to Dell 4.
%
% Section 5: Nested test of Cole model to get a low dimensional description
% of the isoresponse surface (e.g. a plane or a tube)
%
% Section 6: Comparing ViewPixx to Dell4 DTNT threshold surfaces using
% different sets of cone fundamentals
%%
% Section 1
% Looking at superimposed Quest functions
filelist = 'N:\NexFiles\nexfilelists\Greg\SednaQuestLrg.txt';
filenames = fnamesFromTxt(filelist);
defaultcolors = [1 1 1; 1 -1 0; 0 0 1];
defaultcolors = mkbasis(defaultcolors');
colorlabels = {'Ach','L-M','S'};
defaultsfs = [3.2238 0.8893 0.2504];
figure;
for a = 1:size(filenames,1)
    stro = nex2stro(findfile(filenames(a,:)));
    stro = DTfilterquesttrials(stro,0,nan);
    if isempty(stro.sum.analog.sigid) || isempty(stro.trial)
        continue;
    end
    [thresholds, colorDirs, sfs, Q] = DTquestUnpackGH(stro, 'mode', nan);
    for i = 1:size(defaultcolors,2)
        for j = 1:length(defaultsfs)
            clridx = find(abs(mkbasis(colorDirs')'*defaultcolors(:,i)) > .99);
            sfidx = find(abs(sfs-defaultsfs(j)) < 0.001);
            if (~isempty(clridx) & ~isempty(sfidx))
                subplot(3,3,sub2ind([3 3],i,j)); hold on;
                plot(Q{clridx,sfidx},'k-');
            end
        end
    end
end
for i =1:3
    for j = 1:3
        subplot(3,3,sub2ind([3 3],i,j));
        title([colorlabels{i},' ',num2str(defaultsfs(j))]);
    end
end

%%
% Section 2
% Is the monkey better at detecting stimuli on one side over the other?

filelist = 'N:\NexFiles\nexfilelists\Greg\tmp.txt';
filenames = fnamesFromTxt(filelist);
defaultcolors = [1 1 1; 1 -1 0; 0 0 1];
defaultcolors = mkbasis(defaultcolors');
colorlabels = {'Ach','L-M','S'};
defaultlambdas = [16 58 206];
data = zeros(3,3,2,2)
for a = 1:size(filenames,1)
    stro = nex2stro(findfile(filenames(a,:)));
  %  stro = DTfilterquesttrials(stro,0,nan);
    if isempty(stro.sum.analog.sigid) || isempty(stro.trial)
        continue;
    end
    
    flashside = sign(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_x')));
    color_dir = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'color_dir'));
    thresh = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'quest_thresh'));
    lambda = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'gabor_lambda'));
    correct = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));

  %  [thresholds, colorDirs, sfs, Q] = DTquestUnpackGH(stro, 'mode', nan);
    for i = 1:size(defaultcolors,2)
        for j = 1:length(defaultlambdas)
            clridx = find(abs(mkbasis(colorDirs')'*defaultcolors(:,i)) > .99);
            if (~isempty(clridx))
                L = lambda == defaultlambdas(j) & color_dir == clridx;
                Lside = flashside == -1;   % T1 presentations
                a = sum(correct(L & Lside));
                b = sum(correct(L & ~Lside));
                c = sum(~correct(L & Lside));
                d = sum(~correct(L & ~Lside));
                data(i,j,:,:) = squeeze(data(i,j,:,:))+[a b; c d];
            end
        end
    end
end

%       T1  T2
% corr   a  b
% inc    c  d

% First index into data is color, second is spatial frequency
for i = 1:3
    for j = 1:3
        x = squeeze(data(i,j,:,:))
        subplot(3,3,sub2ind([3 3],i,j));
        imagesc(x); 
        colormap(gray);
        set(gca,'XTick',[],'YTick',[]);
        m1 = sum(x,1);
        m2 = sum(x,2);
        tot = sum(x(:));
        expected = (m2./sum(m1))*(m1./sum(m2))*tot;
        chi2stat = sum(sum((expected-x).^2./expected))
        p = 1-chi2cdf(chi2stat, 1);
        title(['p = ',num2str(p)]);
        if (j == 3)
            xlabel(colorlabels{i});
        end
        if (i == 1)
             ylabel(num2str(defaultlambdas(j)));
        end
    end
end
% Within the figure (across the subplots):
% Each row is a different spatial period (high sf on left, low sf on right)
% Each column is a different color (Ach, L-M, S)
% Within each subplot
% Each column is a different stimulus location (T1 and T2)
% The top row is "correct" the bottom row is "incorrect"

%%
% Section 3
% Finding Quest files with unusually large (or small) thresholds.  What
% happened?  There appears to be little consistency between conditions. 
% It's not that there are a few experiments during which behavior was
% uniformly terrible.
filelist = 'N:\NexFiles\nexfilelists\Greg\KaliQuestLrg.txt';
filenames = fnamesFromTxt(filelist);
defaultcolors = [1 1 1; 1 -1 0; 0 0 1];
defaultcolors = mkbasis(defaultcolors');
colorlabels = {'Ach','L-M','S'};
defaultsfs = [0.2504 0.8893 3.2238];
data = [];
for a = 1:size(filenames,1)
    stro = nex2stro(findfile(filenames(a,:)));
    stro = DTfilterquesttrials(stro,0,nan);
    if isempty(stro.sum.analog.sigid) || isempty(stro.trial)
        continue;
    end
    [thresholds, colorDirs, sfs, Q] = DTquestUnpackGH(stro, 'mode', nan);
    tmp = nan*ones(size(defaultcolors,2),length(defaultsfs));
    for i = 1:size(defaultcolors,2)
        for j = 1:length(defaultsfs)
            clridx = find(abs(mkbasis(colorDirs')'*defaultcolors(:,i)) > .99);
            sfidx = find(abs(sfs-defaultsfs(j)) < 0.001);
            if (~isempty(clridx) & ~isempty(sfidx))
                tmp(i,j) = thresholds(clridx,sfidx);
            end
        end
    end
    data = cat(3,data,tmp);
end

figure;
for i = 1:size(data,1)*size(data,2)
    subplot(3,3,i);
    [sfidx,clridx] = ind2sub([size(data,1),size(data,2)],i);
    plot(squeeze(data(sfidx,clridx,:)));
    set(gca,'YScale','log','XLim',[1 size(data,3)]);
    set(gca,'YLim',[.01 100]);
end

%%
% Section 4
% DTNT surface fits for multiple observers.
% Also comparing data collected on multiple display devices
% (e.g. Dell4 and ViewPixx) to make sure there is consistency across these.
if (ismac)
    listpath = '/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/DTNT';
else
    listpath = 'N:\NexFiles\nexfilelists\Greg\DTNT';
end

[~,txt,raw]=xlsread([listpath,filesep,'DTNT15Hz']);
subjects = unique(txt(:,1));
displays = unique(txt(:,2));
filenames = txt(:,3);
data = [];
for subjectidx = 1:length(subjects)
    Lsub = strcmp(subjects{subjectidx},txt(:,1));
    for fileidx = find(Lsub)'
        whichdisplay = find(strcmp(displays,txt(fileidx,2)));
        stro = nex2stro(findfile(filenames{fileidx}));
        [thresholds, colordirs, sfs, QuestTrajectories] = DTquestUnpackGH(stro, 'mode');   
        threshs = thresholds';
        lmsmat = mkbasis(colordirs')';
        repmat(thresholds,1,3).*colordirs;
        Loog = [];
        for j = 1:length(QuestTrajectories)
<<<<<<< .mine
            Loog(j,1) = QuestTrajectories{j}(end) == max(QuestTrajectories{j});
%            Loog(j,1) = QuestTrajectories{j}(end) >= .9*max(QuestTrajectories{j});

=======
            Loog(j,1) = QuestTrajectories{j}(end) == max(QuestTrajectories{j});
%            Loog(j,1) = QuestTrajectories{j}(end) >= .9*max(QuestTrajectories{j});
>>>>>>> .r2960
        end
        data = [data; repmat(subjectidx,3,1) repmat(whichdisplay,3,1) repmat(thresholds,1,3).*colordirs Loog];
    end
end
% Columns: 1)Subject, 2)Display, 3)L, 4)M, 5)S, 6)Loog
% Getting the gamut of the ViewPixx for plotting
load Dell4BitsCal.mat
%load ViewPixx.mat
load T_cones_smj10
fundamentals = T_cones_smj10;
fundWavelengthSpacing = S_cones_smj10;
calData = cals{end};
ptb.bkgndRGB = round(255 .* calData.bgColor); %bkgnd voltages discretized b/w 0&255
ptb.bkgndrgb = [calData.gammaTable(ptb.bkgndRGB(1)+1, 1), calData.gammaTable(ptb.bkgndRGB(2)+1, 2), ...
    calData.gammaTable(ptb.bkgndRGB(3)+1, 3)]; %add one to create and index b/w 1 & 256

ptb.monSpd = SplineSpd(calData.S_device, calData.P_device, fundWavelengthSpacing);
ptb.M = fundamentals * ptb.monSpd;
ptb.bkgndlms = ptb.M * ptb.bkgndrgb';

rgb_cube_verts = fullfact([2 2 2]) - 1;
lms_cube_verts = ptb.M * rgb_cube_verts';

cc_cube_verts = bsxfun(@rdivide, bsxfun(@minus, lms_cube_verts, ptb.bkgndlms), ptb.bkgndlms)';
% K = convhull(cc_cube_verts);
% T == NaN -> points outside the gamut
lims = max(abs(data(:,[3 4 5])));
npts = 75;
[xx yy zz] = meshgrid(linspace(-lims(1),lims(1),npts),...
    linspace(-lims(2),lims(2),npts),...
    linspace(-lims(3),lims(3),npts));
xformedxyz = [xx(:) yy(:) zz(:)];
PERCGAMUT = 90; % For fitting, how near to gamut edge to get (in each direction)
T = tsearchn(cc_cube_verts*PERCGAMUT, delaunayn(cc_cube_verts*PERCGAMUT), xformedxyz);
xformedxyz(isnan(T),:) = nan;
T = tsearchn(cc_cube_verts*PERCGAMUT, delaunayn(cc_cube_verts*PERCGAMUT), -xformedxyz);
xformedxyz(isnan(T),:) = nan;
%plot3(xformedxyz(~isnan(T),1),xformedxyz(~isnan(T),2),xformedxyz(~isnan(T),3),'k,') % Check

%%
% Plotting data from one or more observers
figure; axes; hold on; camlight; lighting phong; axis square;
%whichsubjects = {'Zack','Greg','Leah'};
OVERLAYPLOTS = 1;
whichsubjects = {'Freya','Nut','Sedna'};
<<<<<<< .mine
<<<<<<< .mine
whichsubjects = {'Freya'};
whichdisplay = 'Dell 4';
Ldisplay = data(:,2) == find(strcmp(whichdisplay,displays));
if (~OVERLAYPLOTS)
    figure; axes; hold on; camlight; lighting phong; axis square;
end
=======
whichsubjects = {'Nut'};
=======
whichsubjects = {'Sedna'};
>>>>>>> .r2962
DISPLAYNAME = 'Dell 4';

>>>>>>> .r2960
colors = hsv(length(whichsubjects));
Loog = logical(data(:,6) == 1);
<<<<<<< .mine

FIT = 3; % 1 = plane, 2 = quad, 3 = superquadric, 0 = nothing
=======
FIT = 1; % 1 = plane, 2 = quad, 3 = superquadric, 0 = nothing
>>>>>>> .r2960
hs = [];
figure; axes; hold on;
for i = 1:length(whichsubjects)
    subjectidx = find(strcmp(subjects, whichsubjects(i)));
<<<<<<< .mine
    L = logical(data(:,1) == subjectidx) & Ldisplay;
=======
    displayidx = find(strcmp(DISPLAYNAME, displays));
    if (~isempty(displayidx))
        Ldisplay = logical(displayidx == data(:,2));
    else
        Ldisplay = true(size(data,1),1);
    end
    L = logical(data(:,1) == subjectidx);
>>>>>>> .r2960
    L = L & Ldisplay;
    
    for j = [-1, 1]
        h = plot3(j*data(L&~Loog,3),j*data(L&~Loog,4),j*data(L&~Loog,5),'ko');
        set(h,'MarkerFaceColor',colors(i,:),'MarkerSize',4,'MarkerEdgeColor','none');
        plot3([zeros(sum(L&Loog),1) j*data(L&Loog,3)]',[zeros(sum(L&Loog),1) j*data(L&Loog,4)]',[zeros(sum(L&Loog),1) j*data(L&Loog,5)]','k-','Color',colors(i,:));
    end
    hs = [hs;subjectidx h];
    % Fitting data
    [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(data(L,[3 4 5]), Loog(L));
    planeparams=(planeparams'*xformmat'); % adjusting planeparams
    A = [quadparams(1) quadparams(4) quadparams(5);... % adjusting quadparams
        quadparams(4) quadparams(2) quadparams(6);...
        quadparams(5) quadparams(6) quadparams(3)];
    B = xformmat*A*xformmat';
    quadparams = [B(1,1) B(2,2) B(3,3) B(1,2) B(1,3) B(2,3)]';
    variables = [xformedxyz(:,1).^2 xformedxyz(:,2).^2 xformedxyz(:,3).^2 2*xformedxyz(:,1).*xformedxyz(:,2) 2*xformedxyz(:,1).*xformedxyz(:,3) 2*xformedxyz(:,2).*xformedxyz(:,3)];
    if (FIT == 1)
        coefficients = [planeparams.*planeparams planeparams(1)*planeparams(2) planeparams(1)*planeparams(3) planeparams(2)*planeparams(3)]';
        fr = variables*coefficients;
    elseif (FIT == 2)
        fr = variables*quadparams;
    elseif (FIT == 3)
        [evecs, evals] = eig(B);
        initparams = [2; reshape(evecs*sqrt(evals),9,1)];
        options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-5,'Display','none');
        [fpar,fv] = fminsearch(@(x) colefiterr(x,data(L,[3 4 5]),Loog(L),0),initparams, options);
        fr = sum(abs(xformedxyz *reshape(fpar(2:end),3,3)).^fpar(1),2);        
    end
    if (FIT)
        surfstruct = isosurface(xx,yy,zz,reshape(fr,size(xx)), 1);
        p = patch(surfstruct);
        set(p,'EdgeColor', 'none', 'FaceAlpha',.5,'FaceColor',colors(i,:),'Edgealpha',0);
    end
end
legend(hs(:,2), subjects(hs(:,1)));
xlabel('\DeltaL/L'); ylabel('\DeltaM/M'); zlabel('\DeltaS/S');
%set(gca,'Ylim',[-20 20]);

%%
% For a single subject, comparing CRT and ViewPixx
% need "data" from section 4
whichsubject = 'Freya';
colors = [1 0 0; 0 .5 0; 0 0 1];
figure; axes; hold on;
Lsub = logical(data(:,1) == find(ismember(subjects, whichsubject)));
Loog = logical(data(:,6) == 1);
hs = [];

% Setting up monitor gamut for function fit plotting
load Dell4BitsCal.mat
%load ViewPixx.mat
load T_cones_smj10
fundamentals = T_cones_smj10;
fundWavelengthSpacing = S_cones_smj10;
calData = cals{end};
ptb.bkgndRGB = round(255 .* calData.bgColor); %bkgnd voltages discretized b/w 0&255
ptb.bkgndrgb = [calData.gammaTable(ptb.bkgndRGB(1)+1, 1), calData.gammaTable(ptb.bkgndRGB(2)+1, 2), ...
    calData.gammaTable(ptb.bkgndRGB(3)+1, 3)]; %add one to create and index b/w 1 & 256

ptb.monSpd = SplineSpd(calData.S_device, calData.P_device, fundWavelengthSpacing);
ptb.M = fundamentals * ptb.monSpd;
ptb.bkgndlms = ptb.M * ptb.bkgndrgb';
rgb_cube_verts = fullfact([2 2 2]) - 1;
lms_cube_verts = ptb.M * rgb_cube_verts';
cc_cube_verts = bsxfun(@rdivide, bsxfun(@minus, lms_cube_verts, ptb.bkgndlms), ptb.bkgndlms)';
% K = convhull(cc_cube_verts);
% T == NaN -> points outside the gamut
lims = max(abs(data(:,[3 4 5])));
npts = 75;
[xx yy zz] = meshgrid(linspace(-lims(1),lims(1),npts),...
    linspace(-lims(2),lims(2),npts),...
    linspace(-lims(3),lims(3),npts));
xformedxyz = [xx(:) yy(:) zz(:)];
PERCGAMUT = 90; % For fitting, how near to gamut edge to get (in each direction)
T = tsearchn(cc_cube_verts*PERCGAMUT, delaunayn(cc_cube_verts*PERCGAMUT), xformedxyz);
xformedxyz(isnan(T),:) = nan;
T = tsearchn(cc_cube_verts*PERCGAMUT, delaunayn(cc_cube_verts*PERCGAMUT), -xformedxyz);
xformedxyz(isnan(T),:) = nan;

for i = 1:length(displays)
    Ldisplay = logical(i == data(:,2));
    L = Ldisplay&Lsub;
    for j = [-1, 1]
        h = plot3(j*data(L&~Loog,3),j*data(L&~Loog,4),j*data(L&~Loog,5),'ko');
        set(h,'MarkerFaceColor',colors(i,:),'MarkerSize',5,'MarkerEdgeColor','none');
        plot3([zeros(sum(L&Loog),1) j*data(L&Loog,3)]',[zeros(sum(L&Loog),1) j*data(L&Loog,4)]',[zeros(sum(L&Loog),1) j*data(L&Loog,5)]','k-','Color',colors(i,:));
    end
    [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(data(L,[3 4 5]), Loog(L));
    planeparams=(planeparams'*xformmat'); % adjusting planeparams
    A = [quadparams(1) quadparams(4) quadparams(5);... % adjusting quadparams
        quadparams(4) quadparams(2) quadparams(6);...
        quadparams(5) quadparams(6) quadparams(3)];
    B = xformmat*A*xformmat';
    quadparams = [B(1,1) B(2,2) B(3,3) B(1,2) B(1,3) B(2,3)]';
    variables = [xformedxyz(:,1).^2 xformedxyz(:,2).^2 xformedxyz(:,3).^2 2*xformedxyz(:,1).*xformedxyz(:,2) 2*xformedxyz(:,1).*xformedxyz(:,3) 2*xformedxyz(:,2).*xformedxyz(:,3)];
    
    [evecs, evals] = eig(B);
    initparams = [2; reshape(evecs*sqrt(evals),9,1)];
    options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-5,'Display','none');
    [fpar,fv] = fminsearch(@(x) colefiterr(x,data(L,[3 4 5]),Loog(L),0),initparams, options);
    fr = sum(abs(xformedxyz *reshape(fpar(2:end),3,3)).^fpar(1),2);
    
    surfstruct = isosurface(xx,yy,zz,reshape(fr,size(xx)), 1);
    p = patch(surfstruct);
    set(p,'EdgeColor', 'none', 'FaceAlpha',.5,'FaceColor',colors(i,:),'Edgealpha',0);
    
    hs = [hs; h];
end
if (length(hs) == length(displays))
    legend(hs,displays);
end
set(gca,'View',[-141 -38]);
title(whichsubject)
xlabel('\DeltaL/L'); ylabel('\DeltaM/M'); zlabel('\DeltaS/S');
%%
% For a single subject 
% Plot DTMacPig thresholds on the same axes (as Section 3) for comparison
% (as magenta stars)

<<<<<<< .mine
<<<<<<< .mine
whichsubject = 'Freya';
=======
whichsubject = 'Leah';
=======
whichsubject = 'Freya';
>>>>>>> .r2962
DISPLAYNAME = 'Both';
>>>>>>> .r2960
figure; axes; hold on;
displayidx = find(strcmp(DISPLAYNAME, displays));
if (~isempty(displayidx))
    Ldisplay = logical(displayidx == data(:,2));
else
    Ldisplay = true(size(data,1),1);
end
L = logical(data(:,1) == find(ismember(subjects, whichsubject)));
L = L & Ldisplay;
Loog = logical(data(:,6) == 1);
hs = [];
for j = [-1, 1]
    h = plot3(j*data(L&~Loog,3),j*data(L&~Loog,4),j*data(L&~Loog,5),'ko');
    set(h,'MarkerFaceColor','red','MarkerSize',5,'MarkerEdgeColor','none');
    plot3([zeros(sum(L&Loog),1) j*data(L&Loog,3)]',[zeros(sum(L&Loog),1) j*data(L&Loog,4)]',[zeros(sum(L&Loog),1) j*data(L&Loog,5)]','r-');
end
title(whichsubject)

% Plotting a 2-D Cole fit

% First fitting a plane and quadratic for initial guesses
[planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(data(L,[3 4 5]), Loog(L));
planeparams=(planeparams'*xformmat'); % adjusting planeparams
A = [quadparams(1) quadparams(4) quadparams(5);... % adjusting quadparams
    quadparams(4) quadparams(2) quadparams(6);...
    quadparams(5) quadparams(6) quadparams(3)];
B = xformmat*A*xformmat';
quadparams = [B(1,1) B(2,2) B(3,3) B(1,2) B(1,3) B(2,3)]';
[evecs, evals] = eig(B);
initparams = [2; reshape(evecs*sqrt(evals),9,1)];
options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-5,'Display','none');
[fpar,fv] = fminsearch(@(x) colefiterr(x,data(L,[3 4 5]),Loog(L),0),initparams, options);
fpar'
fv
[fpar, fval] = fitDetectionSurface(data(L,[3 4 5]), sqrt(sum(data(L,[3 4 5]).^2,2)),'pso',1,Loog(L))

% Plotting
npts = 75;
[xx yy zz] = meshgrid(linspace(-lims(1),lims(1),npts),...
    linspace(-lims(2),lims(2),npts),...
    linspace(-lims(3),lims(3),npts));
xformedxyz = [xx(:) yy(:) zz(:)];
PERCGAMUT = 90; % For fitting, how near to gamut edge to get (in each direction)
T = tsearchn(cc_cube_verts*PERCGAMUT, delaunayn(cc_cube_verts*PERCGAMUT), xformedxyz);
xformedxyz(isnan(T),:) = nan;
T = tsearchn(cc_cube_verts*PERCGAMUT, delaunayn(cc_cube_verts*PERCGAMUT), -xformedxyz);
xformedxyz(isnan(T),:) = nan;
fr = sum(abs(xformedxyz *reshape(fpar(2:end),3,3)).^fpar(1),2);        
%fr = abs(xformedxyz*planeparams');
surfstruct = isosurface(xx,yy,zz,reshape(fr,size(xx)), 1);
p = patch(surfstruct);
set(p,'EdgeColor', 'none', 'FaceAlpha',.5,'FaceColor','red','Edgealpha',0);


if ispc
    flpmacpig = 'N:\NexFiles\nexfilelists\Greg\DTEM';
else
    flpmacpig = '/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/DTEM';
end

[flmacpig,flpmacpig,~] = uigetfile([flpmacpig filesep '*.txt'], ['Please select the MacPig for ',whichsubject]);
if ~isequal(flmacpig,0)    
    macpig_filenames = fnamesFromTxt2([flpmacpig filesep flmacpig]);
    macpig_data = [];
    for file = macpig_filenames'
        stro = nex2stro(findfile(char(file{:})));
        if stro.sum.exptParams.rf_x == 50 % only get the 5 degree eccentric thresholds
            [threshes, color_dirs] = DTquestUnpackGH(stro, 'mode');
            macpig_scaled = bsxfun(@times, color_dirs, threshes./sqrt(sum(color_dirs.^2, 2))); % make color_dirs unit vectors and multiply by thresholds
            blueidx = find(color_dirs(:,3) == max(color_dirs(:,3))); % blue has the greatest S-cone component
            greenidx = find(color_dirs(:,2) == max(color_dirs(:,2)));% green has the greatest M-cone component
            macpig_data = [macpig_data; permute([macpig_scaled(greenidx,:); macpig_scaled(blueidx,:)],[3 2 1])];
        end
    end
end
c = {'green','blue'};
for i = 1:size(macpig_data,3) % looping over green, blue
    h(1) = plot3(macpig_data(:,1,i), macpig_data(:,2,i), macpig_data(:,3,i), 'mp', 'MarkerSize', 10,'MarkerFaceColor',c{i},'MarkerEdgeColor',c{i});
    h(2) = plot3(-macpig_data(:,1,i), -macpig_data(:,2,i), -macpig_data(:,3,i), 'mp', 'MarkerSize', 10,'MarkerFaceColor',c{i},'MarkerEdgeColor',c{i});
end


% Plotting the macpig data on the fitted surface and conditioning so that
% you can see that the magpig data line up with the DTNT isodetection
% surfaces
unitvect_g = squeeze(macpig_data(1,:,1));
unitvect_g = unitvect_g./norm(unitvect_g);
unitvect_b = squeeze(macpig_data(1,:,2));
unitvect_b = unitvect_b./norm(unitvect_b);
nullvect = null([unitvect_b; unitvect_g]); % orthogonal to b and g. Should point into the page of the plot.

v = surfstruct.vertices;
nvproj = v*nullvect;
TOL = .25;
Lproj = abs(nvproj)<TOL;
figure; axes; hold on;
plot3(v(Lproj,1),v(Lproj,2),v(Lproj,3),'k.')

view(nullvect);
c = {'green','blue'};
for i = 1:size(macpig_data,3) % looping over green, blue
    h(1) = plot3(macpig_data(:,1,i), macpig_data(:,2,i), macpig_data(:,3,i), 'mp', 'MarkerSize', 10,'MarkerFaceColor',c{i},'MarkerEdgeColor',c{i});
    h(2) = plot3(-macpig_data(:,1,i), -macpig_data(:,2,i), -macpig_data(:,3,i), 'mp', 'MarkerSize', 10,'MarkerFaceColor',c{i},'MarkerEdgeColor',c{i});
end
title(whichsubject);

%%
% Section 5
% Fitting the Cole et al. model to isodetection data and building it up bit by bit. 
% Can we get away with fewer than three mechanisms? 

if (ismac)
    listpath = '/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/DTNT';
else
    listpath = 'N:\NexFiles\nexfilelists\Greg\DTNT';
end
filelist = 'DTNT15Hz.txt';
filenames = fnamesFromTxt2([listpath,filesep,filelist]);
data = [];
for a = 1:size(filenames,1)
    if (length(char(filenames{a})) > 6)   % Only taking real files
        stro = nex2stro(findfile(filenames{a}));
        [thresholds, colordirs, sfs, QuestTrajectories] = DTquestUnpackGH(stro, 'mode');
        threshs = thresholds';
        lmsmat = mkbasis(colordirs')';
        repmat(thresholds,1,3).*colordirs;
        Loog = [];
        for j = 1:length(QuestTrajectories)
            Loog(j,1) = QuestTrajectories{j}(end) == max(QuestTrajectories{j});
        end
        data = [data; repmat(thresholds,1,3).*colordirs Loog];
    end
end
% First fitting a plane and quadratic for initial guesses
Loog = data(:,end);
[planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(data(:,[1 2 3]), Loog);
planeparams=(planeparams'*xformmat'); % adjusting planeparams
A = [quadparams(1) quadparams(4) quadparams(5);... % adjusting quadparams
    quadparams(4) quadparams(2) quadparams(6);...
    quadparams(5) quadparams(6) quadparams(3)];
B = xformmat*A*xformmat';
quadparams = [B(1,1) B(2,2) B(3,3) B(1,2) B(1,3) B(2,3)]';
% Estimated L cone ratio 
% planeparams(1)./(planeparams(1)+planeparams(2))
[evecs, evals] = eig(B);
initparams = [2; reshape(evecs*real(sqrt(evals)),9,1)];
options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-5,'Display','none');
lms = data(:,[1 2 3]);

% Fitting reduced Cole model (1-D)
initparams = planeparams;
[fpar1,fv1] = fminsearch(@(x) colefiterr(x,lms,Loog,0),initparams, options);

% Fitting reduced Cole model (2-D) 
initparams =  [2; planeparams'; reshape(evecs*sqrt(evals),9,1)];
[fpar2(:,1),fv2(1)] = fminsearch(@(x) colefiterr(x,lms,Loog,0),initparams(1:7), options);
initparams =  [2; planeparams'; 0;0;0];
[fpar2(:,2),fv2(2)] = fminsearch(@(x) colefiterr(x,lms,Loog,0),initparams(1:7), options);
initparams = [2; reshape(evecs*real(sqrt(evals)),9,1)]; % Seems to be best
[fpar2(:,3),fv2(3)] = fminsearch(@(x) colefiterr(x,lms,Loog,0),initparams(1:7), options);
fv2
fpar2 = fpar2(:,min(fv2) == fv2);
fv2 = min(fv2);

% Fitting reduced Cole model (3-D)
initparams = [fpar2; 0; 0; 0];
[fpar3(:,1),fv3(1)] = fminsearch(@(x) colefiterr(x,lms,Loog,0),initparams, options);
initparams = [fpar2(1); reshape(evecs*real(sqrt(evals)),9,1)];
[fpar3(:,2),fv3(2)] = fminsearch(@(x) colefiterr(x,lms,Loog,0),initparams, options);
fv3
fpar3 = fpar3(:,min(fv3) == fv3);
fv3 = min(fv3);


%chi squared tests
% This could be a problem because we're not really fitting by maximum
% likelihood.
% Numerator is smaller model
%stat = -2*log(exp(-fv1)./exp(-fv2));
%stat = -2*log(exp(-fv1) - log(exp(-fv2));
%stat = -2*((-fv1) - (-fv2));
%stat = -2*(-fv1+fv2);

stat = 2*(fv1-fv2)
p = 1-chi2cdf(stat,length(fpar2)-length(fpar1))

stat = 2*(fv2-fv3)
p = 1-chi2cdf(stat,length(fpar3)-length(fpar2))

%%
% Section 6
% Trying to represent DTNT surfaces in cone contrast spaces with whichever
% set of cone fundamentals the user wants to use.
if (ismac)
    listpath = '/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/DTNT';
else
    listpath = 'N:\NexFiles\nexfilelists\Greg\DTNT';
end
load('T_cones_smj10.mat'); % These are the standard fundamentals used with DTNT
load('T_cones_synthgh2'); % These are the cone fundamentals we're converting to
newfunds = T_cones_synthgh2;
[~,txt,raw]=xlsread([listpath,filesep,'DTNT15Hz']);
subjects = unique(txt(:,1));
displays = unique(txt(:,2));
filenames = txt(:,3);
data = [];
for subjectidx = 1:length(subjects)
    Lsub = strcmp(subjects{subjectidx},txt(:,1));
    for fileidx = find(Lsub)'
        whichdisplay = find(strcmp(displays,txt(fileidx,2)));
        stro = nex2stro(findfile(filenames{fileidx}));
        monspd = reshape(stro.sum.exptParams.mon_spect, length(stro.sum.exptParams.mon_spect)/3,3);
        bkgndrgb = [stro.sum.exptParams.bkgnd_r stro.sum.exptParams.bkgnd_g stro.sum.exptParams.bkgnd_b]; % intensities
        computedM = T_cones_smj10*monspd;
        stroM = reshape(stro.sum.exptParams.m_mtx,3,3);
        if (any(any(~isequal(stroM, computedM))))
            error('M in stro file differs from computed M');
        end
        bkgndspect = monspd*bkgndrgb';
        bkgndoldlms = stroM*bkgndrgb';  % background LMS in terms of the old fundamentals
        bkgndnewlms = newfunds'*bkgndspect;  % background LMS in terms of the new fundamentals
        
        [thresholds, colordirs, sfs, QuestTrajectories] = DTquestUnpackGH(stro, 'mode');
        oldccthresh = repmat(thresholds,1,3).*colordirs./100;
        newlmscc = ConvertConeContrastBasis(stroM, newfunds'*monspd, bkgndrgb, oldccthresh);
        
        Loog = [];
        for j = 1:length(QuestTrajectories)
            Loog(j,1) = QuestTrajectories{j}(end) == max(QuestTrajectories{j});
%            Loog(j,1) = QuestTrajectories{j}(end) >= .9*max(QuestTrajectories{j});
        end
        data = [data; repmat(subjectidx,3,1) repmat(whichdisplay,3,1) newlmscc.*100  Loog];
    end
end
% Columns: 1)Subject, 2)Display, 3)L, 4)M, 5)S, 6)Loog
