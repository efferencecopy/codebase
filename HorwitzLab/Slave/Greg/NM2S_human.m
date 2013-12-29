% Psychophysical task to probe color contrast effects.
% Do we base judgements of color similarity on cone excitation
% ratios?  Experiment is set up as an non-match to sample 
% where the sample is always the same.
%
% Parameters that define a trial
% 1) side with target
% 2) gradient type
% 3) distractor type
%
% GDLH 10/27/08

MONDIST = 100; % cm
gl.screenWidthcm = 40; %  from REX
stimx = 6;
stimy = 0;
ngrads = 3;
nextrastim = 2; % number of extra stim past max grad (x2 for + and -)
gradpos = 5;
stimsizeinpix = 20;
samplecc = [-.1 -.1 -.1];

%gradlimscc = [0 0 0; 0 0 0];
gradlimscc = [0 0 0.15; 0 0 0];
%gradlimscc = [.0025 -.0025 -.005; 0 0 0];

load ('T_cones_smj.mat');
fundamentals = T_cones_smj;
load('/Monitor Calibration/Monitor data/Dell 4/Dell4BitsCal');
cal = cals{end};
P_device = SplineSpd(SToWls(cal.S_device), cal.P_device, SToWls([380 5 81]));
M = fundamentals*P_device;
invM = inv(M);
bkgndRGB = round(255*cal.bgColor)';
gl.cal.gammaTable = cal.gammaTable;
gl.cal.invgammaTable = InvertGamma(gl.cal.gammaTable, 1);
gl.cal.monSpd = cal.P_device;
gl.mondistcm = MONDIST;

bkgndrgb = [gl.cal.gammaTable(bkgndRGB(1)+1,1); ...
            gl.cal.gammaTable(bkgndRGB(2)+1,2); ...
            gl.cal.gammaTable(bkgndRGB(3)+1,3)];

bkgndlms = M*bkgndrgb;
% Setting up colors for stimuli
samplergb = inv(M)*(bkgndlms'.*(1+samplecc))';
gradmagnitudes = linspace(-1,1,ngrads);
gradccs = gradmagnitudes'*gradlimscc(1,:);
gradlmss = repmat(bkgndlms',ngrads,1).*(1+gradccs);
distractorlmss = gradlmss.*(1+repmat(samplecc,ngrads,1));  % Cone excitations for the distractors
if (ngrads > 1)
    deltalmss = distractorlmss(2,:) - distractorlmss(1,:);
else
    deltalmss = [0 0 bkgndlms(3)/10];  % Cone excitation differences between neighboring distractors
end                         % Ugly code!  GDLH 2/2/10

for i = 1:nextrastim
    distractorlmss = [distractorlmss(1,:)-deltalmss; distractorlmss; distractorlmss(end,:)+deltalmss];
end
ndistractors = ngrads+2*nextrastim;

% Only reinitialize if we haven't collected any data yet
if (~exist('trialtypes') || all(trialtypes(:,5) == 0))
    trialtypes = fullfact([2,ngrads,ndistractors]);
    trialtypes(:,[4 5]) = 0;
end

% Opening the screen
gl.windowPtr = Screen('OpenWindow',0, 255*cal.bgColor);
[screenwidthpix, screenheightpix]  = Screen('WindowSize', gl.windowPtr);
gl.screenCenterXpix = screenwidthpix/2;
gl.screenCenterYpix = screenheightpix/2;
pixpercm = screenwidthpix/gl.screenWidthcm;
cmperdeg = gl.screenWidthcm/(2*atan2(gl.screenWidthcm/2, gl.mondistcm)*180/pi);
gl.pixperdeg = pixpercm*cmperdeg;  % Small pixels

clut = repmat(linspace(0,1,256),3,1)';
Screen('LoadNormalizedGammaTable', 0, clut);

% Creating drawing widows for the targets
targpos = [0 0; stimx stimy; -stimx stimy];
for j = 1:3
    drawrect = round([gl.screenCenterXpix+(targpos(j,1)*gl.pixperdeg)-stimsizeinpix ...
                      gl.screenCenterYpix-(targpos(j,2)*gl.pixperdeg)-stimsizeinpix ...
                      gl.screenCenterXpix+(targpos(j,1)*gl.pixperdeg)+stimsizeinpix...
                      gl.screenCenterYpix-(targpos(j,2)*gl.pixperdeg)+stimsizeinpix]);

    if(rem(drawrect(1), 2)) %if the rectangle starts on an odd pixel
           drawrect(1) = drawrect(1) - 1;
           drawrect(3) = drawrect(3) - 1;
    end
    drawrects(j,:) = drawrect;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Making textures that we want to keep around for the entire expt.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Texture for background screen without the gradient
im = zeros(screenheightpix, screenwidthpix/2, 3);
for plane = 1:3
    tmp = bkgndrgb(plane);
    tmp = round(tmp*size(gl.cal.invgammaTable,1)-1)+1;
    tmp = gl.cal.invgammaTable(tmp, plane);
    tmp = round(tmp*(size(gl.cal.invgammaTable,1)-1));
    im(:,:,plane) = repmat(tmp, screenheightpix, screenwidthpix/2);
end
bkgndtex=Screen('MakeTexture', gl.windowPtr, TranslateToColourMode(im,1));

% Texture for the sample
im = zeros(2*stimsizeinpix, stimsizeinpix, 3);
for plane = 1:3
    tmp = samplergb(plane);
    tmp = round(tmp*size(gl.cal.invgammaTable,1)-1)+1;
    tmp = gl.cal.invgammaTable(tmp, plane);
    tmp = round(tmp*(size(gl.cal.invgammaTable,1)-1));
    im(:,:,plane) = repmat(tmp, 2*stimsizeinpix, stimsizeinpix);
end
sampletex=Screen('MakeTexture', gl.windowPtr, TranslateToColourMode(im,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Top of the big loop
pause();   % hit return to start
trialorder = randperm(size(trialtypes,1));
for i = trialorder
    whichside = trialtypes(i,1);    
    whichgrad = trialtypes(i,2);
    whichdist = trialtypes(i,3);

    % Making background screen with the gradient
    im = zeros(screenheightpix, screenwidthpix/2, 3);
    gradendpointsinpix = round([gl.screenCenterXpix-(gradpos*gl.pixperdeg),...
                            gl.screenCenterXpix+(gradpos*gl.pixperdeg)]./2);
    for cone = 1:3
        tmp = [repmat(gradccs(whichgrad,cone),1,gradendpointsinpix(1)),...
                      linspace(gradccs(whichgrad,cone),gradlimscc(2,cone),gradendpointsinpix(2)-gradendpointsinpix(1)),...
                      repmat(0,1,screenwidthpix/2-gradendpointsinpix(2))];
        gradientlms(cone,:) = bkgndlms(cone)*(1+tmp);
    end
    if (whichside == 2); % nonmatch on right
        gradientlms = fliplr(gradientlms);
    end
    gradientrgb = inv(M)*gradientlms;
    for plane = 1:3
        tmp = gradientrgb(plane,:);
        tmp = round(tmp*size(gl.cal.invgammaTable,1)-1)+1;
        tmp = gl.cal.invgammaTable(tmp, plane);
        tmp = round(tmp*(size(gl.cal.invgammaTable,1)-1));
        im(:,:,plane) = repmat(tmp', screenheightpix, 1);
    end
    bkgndgradtex=Screen('MakeTexture', gl.windowPtr, TranslateToColourMode(im,1));

    % Making the distractor
    distractorrgb = inv(M)*distractorlmss(whichdist,:)';
    im = zeros(2*stimsizeinpix, stimsizeinpix, 3);
    for plane = 1:3
        tmp = distractorrgb(plane);
        tmp = round(tmp*size(gl.cal.invgammaTable,1)-1)+1;
        tmp = gl.cal.invgammaTable(tmp, plane);
        tmp = round(tmp*(size(gl.cal.invgammaTable,1)-1));
        im(:,:,plane) = repmat(tmp, 2*stimsizeinpix, stimsizeinpix);
    end
    distractortex=Screen('MakeTexture', gl.windowPtr, TranslateToColourMode(im,1));

    % drawing
    Screen('DrawTexture',gl.windowPtr,bkgndtex,[],[],[],0);
    Screen('Flip',gl.windowPtr,0,0);
    pause(1);
    Screen('DrawTexture',gl.windowPtr,bkgndtex,[],[],[],0);
    Screen('DrawTexture',gl.windowPtr,sampletex,[],drawrects(1,:),[],0);
    Screen('Flip',gl.windowPtr,0,0);
    pause(1);
    Screen('DrawTexture',gl.windowPtr,bkgndgradtex,[],[],[],0);
    if (whichside == 1)  % distractor on left
        Screen('DrawTexture',gl.windowPtr,sampletex,[],drawrects(2,:),[],0);
        Screen('DrawTexture',gl.windowPtr,distractortex,[],drawrects(3,:),[],0);
    else
        Screen('DrawTexture',gl.windowPtr,sampletex,[],drawrects(3,:),[],0);
        Screen('DrawTexture',gl.windowPtr,distractortex,[],drawrects(2,:),[],0);
    end
    Screen('Flip',gl.windowPtr,0,0); 
    pause(0.1);
    Screen('DrawTexture',gl.windowPtr,bkgndtex,[],[],[],0);
    Screen('Close',bkgndgradtex);
    Screen('Close',distractortex);
    Screen('Flip',gl.windowPtr,0,0);
    FlushEvents;
    ch = GetChar();
    if ((ch == '1') && whichside == 1) || ((ch == '0') && whichside == 2)
        trialtypes(i,5) = trialtypes(i,5)+1;
        trialtypes(i,4) = trialtypes(i,4)+1;
    elseif ((ch == '1') && whichside == 2) || ((ch == '0') && whichside == 1)
        trialtypes(i,4) = trialtypes(i,4)+1;
    elseif ((ch ~= '1') && (ch ~= '0'))
        break
        
    end
end
Screen('Close',bkgndtex);
Screen('CloseAll');

%%
% Some analysis
for i = 1:3
    levels = unique(trialtypes(:,i));
    for j = levels'
        L = trialtypes(:,i) == j;
        n = sum(trialtypes(L,4));
        x = sum(trialtypes(L,5));
        pctcor = x/n
    end
end

pctcor = nan*ones(max(trialtypes(:,2)),max(trialtypes(:,3)));
se = nan*ones(max(trialtypes(:,2)),max(trialtypes(:,3)));
for grad = unique(trialtypes(:,2))'
    for dist = unique(trialtypes(:,3))'
        L = trialtypes(:,2) == grad & trialtypes(:,3) == dist;
        n = sum(trialtypes(L,4));
        x = sum(trialtypes(L,5));
        pctcor(grad,dist) = x/n;
        se(grad,dist) = sqrt((x/n)*(1-x/n))./sqrt(n);
    end
end
imagesc(pctcor);
colormap(gray);

%%
distractordeltas = distractorlmss-repmat((samplecc.*bkgndlms'),size(distractorlmss,1),1);
tmp = repmat(bkgndlms',size(distractorlmss,1),1);

distractorccs = (distractordeltas-tmp)./tmp;
figure; axes; hold on;
plot(distractorccs(:,3),100*pctcor','LineWidth',2);
legend(num2str(gradccs(:,3)));
set(gca,'Xtick',distractorccs(:,3));
errorbar(repmat(distractorccs(:,3),1,3),100*pctcor', 100*se')

xlabel('S-cone contrast');
ylabel('% Chosen');

save trialtypes.dat trialtypes -ascii  % Saves the table "trialtypes" as a text document

