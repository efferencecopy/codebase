% S-coneMatch.m
%
% This script is going to measure the variability in matches made between
% S-cone isolating stimuli, or stimuli with more or less LM excitation
% added to them.  The observer has control over the S-cone excitation of
% the match to match the appearance of a test stimulus.  The theory that
% I'm testing is that Weber's law applies to pure S-cone matches (and maybe
% to matches that consist of some S-cone and some L and M cone excitation),
% but in the special case in which LM and S all covary such that they are
% achromatic (pure radiance modulations from a EEW background) Weber's law
% goes away: i.e. we are very sensitive to deviations from the achromatic 
% locus and this sensitivity is independent of the mean photoreceptor
% excitations.

% Loading calibration data
CALPATH = '/Monitor Calibration/Monitor data/Dell 4';
CALFILE = 'Dell4BitsCal';
CALIDX = 1;
load ([CALPATH,'/',CALFILE]);
cal = cals{CALIDX};

% Expt parameters
NBLOCKS = 1;
NCONTRASTS = 1;

% Display parameters
NPIXPERSIDE = 0;
XPOS = 0;
YPOS = 40;
CONTRASTS = linspace(.1,.5, NCONTRASTS);
COLORDIRS = [0 0 1; .9 1.1 1];  % S cone and achromatic
MAXDAC = 65535;

% Data structure
data = {};
for i = 1:NCONTRASTS
    data{i}.lms = [];
    data{i}.s = [];
    data{i}.contrast = CONTRASTS(i);
end

% Setting up the NormalizedGammaTable
bkgndrgb = cal.bgColor;
oldGammaTable = Screen('ReadNormalizedGammaTable',0);
invGammaTable = InvertGamma(cal, 1);

% Preparing the M matrix
load ('T_cones_smj.mat');
fundamentals = T_cones_smj;
cal = cals{CALIDX};
P_device = SplineSpd(SToWls(cal.S_device), cal.P_device, SToWls([380 5 81]));
M = fundamentals*P_device;
% Ambient is already built into cal.bgColor which is voltage
% bkgndrgb is in intensity
bkgndrgb = [cal.gammaTable(round(256*cal.bgColor(1)),1);...
            cal.gammaTable(round(256*cal.bgColor(2)),2);...
            cal.gammaTable(round(256*cal.bgColor(3)),3)];

% Opening a window
w = Screen('OpenWindow',0,round(255*cal.bgColor));
%w = Screen('OpenWindow',0,[0 0 0]);
HideCursor;
[winWidth, winHeight] = Screen('WindowSize',w);

winCenterX = winWidth/2;
winCenterY = winHeight/2;
drawrect = [-floor(NPIXPERSIDE/2) -floor(NPIXPERSIDE/2) ceil(NPIXPERSIDE/2) ceil(NPIXPERSIDE/2)];
drawrect1 = drawrect + [winCenterX winCenterY+YPOS winCenterX winCenterY+YPOS];
if (drawrect1(1)/2 ~= floor(drawrect1(1)/2))  % ensure that the image starts on an even pixel
    drawrect1(1) = drawrect1(1) + 1;
    drawrect1(3) = drawrect1(3) + 1;
end
drawrect2 = drawrect + [winCenterX winCenterY-YPOS winCenterX winCenterY-YPOS];
if (drawrect2(1)/2 ~= floor(drawrect2(1)/2))  % ensure that the image starts on an even pixel
    drawrect2(1) = drawrect2(1) + 1;
    drawrect2(3) = drawrect2(3) + 1;
end

trialtypes = fullfact([NCONTRASTS,2]);
for block = 1:NBLOCKS
    trialidxs = randperm(size(trialtypes,1));
    for trial = trialidxs
        contrast = CONTRASTS(trialtypes(trial,1));
        dir = trialtypes(trial,2);
        
        % Figuring out where the adjustable patch goes
        whereowhere = unidrnd(2)-1;
        if (whereowhere)
            testdrawrect = drawrect1;
            matchdrawrect = drawrect2;
        else
            testdrawrect = drawrect2;
            matchdrawrect = drawrect1;
        end
        
        %RGBs for test
        lmsbkgnd = M*bkgndrgb;
        lmstest = lmsbkgnd.*(1+COLORDIRS(dir, :)'*contrast);
        rgbtest = inv(M)*lmstest;
        RGBtest = rgb2RGB(rgbtest, invGammaTable);
        testimage = repmat(permute(MAXDAC*RGBtest,[3 2 1]),NPIXPERSIDE,NPIXPERSIDE/2);
        testimage = round(TranslateToColourMode(testimage,1));
        
        MINVAL = -.12; MAXVAL = .12;
        RANDVAL = unifrnd(MINVAL, MAXVAL);
        % Getting user input
        while(1)
            HideCursor;
            [x,y,buttons] = GetMouse(w);  
            if (any(buttons))
                while (any(buttons)) % wait for release
  	                [x,y,buttons] = GetMouse;
                end
                break 
            end
            val = ((x+1)/winWidth)*(MAXVAL-MINVAL)+MINVAL+RANDVAL;
        
            lmsmatch = lmsbkgnd.*(1+COLORDIRS(dir, :)'*contrast);
            lmsmatch(3) = lmsmatch(3)*(1+val);
            rgbmatch = inv(M)*lmsmatch;
            RGBmatch = rgb2RGB(rgbmatch, invGammaTable);
            
            matchimage = repmat(permute(MAXDAC*RGBmatch,[3 2 1]),NPIXPERSIDE,NPIXPERSIDE/2);
            matchimage = round(TranslateToColourMode(matchimage,1));
            matchtex=Screen('MakeTexture', w, matchimage);
            
            Screen('DrawTexture',w,matchtex,[],matchdrawrect,[],0);
            testtex=Screen('MakeTexture', w, testimage);
            Screen('DrawTexture',w,testtex,[],testdrawrect,[],0);
            Screen('Flip',w);
            Screen('Close',testtex);
            Screen('Close',matchtex);
            [keyIsDown, secs, keyCode] = KbCheck;
            if (keyIsDown)
                if (find(keyCode) == 41) % Escape
                   Screen('Close',w);  % This will result in an error
                   break;              % which is ugly but at least is stops the program
                end
            end   
        end
        matchcontrast = (lmsmatch(3)-lmsbkgnd(3))/lmsbkgnd(3);
        if (dir == 1)
           data{trialtypes(trial,1)}.s = [data{trialtypes(trial,1)}.s, matchcontrast];
        elseif (dir == 2)
           data{trialtypes(trial,1)}.lms = [data{trialtypes(trial,1)}.lms, matchcontrast];
        end
    end
end
Screen('Close',w);
ShowCursor;

% Data analysis
figure; axes; hold on;
means = [];
stds = [];
for i = 1:NCONTRASTS
    plot(data{i}.contrast,data{i}.s-data{i}.contrast,'b.');
    plot(data{i}.contrast,data{i}.lms-data{i}.contrast,'k.');
    means = [means; mean(data{i}.s), mean(data{i}.lms)];
    stds = [stds; std(data{i}.s), std(data{i}.lms)]; 
end
xlabel('S cone contrast');
ylabel('delta S cone contrast');

% And another figure
figure; subplot(2,1,1); hold on;
plot(CONTRASTS,stds(:,1),'k-');
plot(CONTRASTS,stds(:,2),'m-');
legend({'S only','LMS'});
ylabel('stdev S cone contrast');
subplot(2,1,2); hold on;
plot(CONTRASTS,means(:,1)-CONTRASTS','k-');
plot(CONTRASTS,means(:,2)-CONTRASTS','m-');
ylabel('mean S cone contrast');
xlabel('S cone contrast');

% And another plot
tmpdata = [];
for i = 1:NCONTRASTS
    tmpdata(:,i,1) = data{i}.s'-data{i}.contrast;
    tmpdata(:,i,2) = data{i}.lms'-data{i}.contrast;
end
figure; axes; hold on;
plot(prctile(squeeze(tmpdata(:,:,1)),75),'b:');
plot(prctile(squeeze(tmpdata(:,:,1)),25),'b:');
plot(median(squeeze(tmpdata(:,:,1))),'b-');
plot(prctile(squeeze(tmpdata(:,:,2)),75),'k:');
plot(prctile(squeeze(tmpdata(:,:,2)),25),'k:');
plot(median(squeeze(tmpdata(:,:,2))),'k-');
