% Script to present a chromatic (S-cone isolating?) test stimulus on a
% Gaussian white noise luminance pedestal.  The idea is to make a classification
% image to see the temporal dynamics of the luminance facilitation.

% Loading calibration data
CALPATH = '/Monitor Calibration/Monitor data/Dell 4';
CALFILE = 'Dell4BitsCal';
CALIDX = 1;
load ([CALPATH,'/',CALFILE]);
cal = cals{CALIDX};

% Display parameters
NTRIALS = 500;
NPIXPERSIDE = 30;
NFRAMESPED = 50;
NFRAMESTEST = 10;
NFRAMETESTSTART = NFRAMESPED-NFRAMESTEST-5;
NOISEAMP = 0.05;  % standard deviation for luminance noise process
SIGNALCONTRAST = 0.075;  % cone contrast
SIGNALLMS = [0 0 1];
MAXDAC = 65535;
FPXY = [-200 0];
FPPIX = 6;

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
HideCursor;
[winWidth, winHeight] = Screen('WindowSize',w);

winCenterX = winWidth/2;
winCenterY = winHeight/2;

% The stimulus drawing window
drawrect = [-floor(NPIXPERSIDE/2) -floor(NPIXPERSIDE/2) ceil(NPIXPERSIDE/2) ceil(NPIXPERSIDE/2)];
drawrect = drawrect + [winCenterX winCenterY winCenterX winCenterY];
if (drawrect(1)/2 ~= floor(drawrect(1)/2))  % ensure that the image starts on an even pixel
    drawrect(1) = drawrect(1) + 1;
    drawrect(3) = drawrect(3) + 1;
end

% The fixation pint drawing window
fprect = [-floor(FPPIX/2) -floor(FPPIX/2) ceil(FPPIX/2) ceil(FPPIX/2)];
fprect = fprect+[winCenterX+FPXY(1) winCenterY+FPXY(2) winCenterX+FPXY(1) winCenterY+FPXY(2)];
if (fprect(1)/2 ~= floor(fprect(1)/2))  % ensure that the image starts on an even pixel
    fprect(1) = fprect(1) + 1;
    fprect(3) = fprect(3) + 1;
end
fpimage = zeros(FPPIX,FPPIX,3);

% Preparing the signal vector
bkgndlms = M*bkgndrgb;
siglms = bkgndlms.*SIGNALCONTRAST.*SIGNALLMS'./norm(SIGNALLMS);
sigrgb = inv(M)*siglms;

if ~exist('hits')
    hits = []; 
    misses = []; 
    FAs = []; 
    CRs = [];
end
pause();  % need to hit return to get the first trial
for trial = 1:NTRIALS
    % First the noise
    NoiseTraj = normrnd(0, NOISEAMP, 1, NFRAMESPED);
    Noisergb = bkgndrgb*NoiseTraj;

    sigpresent = unidrnd(2)-1;
    SigTraj = zeros(1,NFRAMESPED);
    if (sigpresent)
        SigTraj(NFRAMETESTSTART:NFRAMETESTSTART+NFRAMESTEST) = 1;
    end
    Sigrgb = sigrgb*SigTraj;
    rgb = repmat(bkgndrgb,1,NFRAMESPED)+Sigrgb+Noisergb;
    rgb(:,end) = bkgndrgb;
    
    RGBs = [];
    for i = 1:NFRAMESPED
        RGBs(i,:) = [invGammaTable(round(rgb(1,i)*65535)+1,1),...
                     invGammaTable(round(rgb(2,i)*65535)+1,2),...
                     invGammaTable(round(rgb(3,i)*65535)+1,3)];
    end

    % Actually doing the drawing
    for framecount = 1:NFRAMESPED
        image = repmat(permute(MAXDAC*RGBs(framecount,:),[3 1 2]),NPIXPERSIDE,NPIXPERSIDE/2);
        image = round(TranslateToColourMode(image,1));
        tex=Screen('MakeTexture', w, image);
        Screen('DrawTexture',w,tex,[],drawrect,[],0);
        % drawing the fix point
        tex=Screen('MakeTexture', w, fpimage);
        Screen('DrawTexture',w,tex,[],fprect,[],0);
        Screen('Flip',w);
        Screen('Close',tex);
    end
    FlushEvents;
    ch = GetChar();
    if ((ch == '1') & sigpresent)
        hits = [hits; NoiseTraj];
    elseif ((ch == '0') & sigpresent)
        misses = [misses; NoiseTraj];
    elseif((ch == '1') & ~sigpresent)
        FAs = [FAs; NoiseTraj];
    elseif((ch == '0') & ~sigpresent)
        CRs = [CRs; NoiseTraj];
    end
end

hz = Screen('FrameRate',w);
Screen('Close',w);
s = input('Save data?','s');
if (~isempty(s) & s(1) == 'y' | s(1) == 'Y')
%     %use a save dialog here: these paths are out-of-date. zlb 01/2012
%     filename = input('filename: ','s');
%     savepath = '/Users/horwitzlab/MatlabCode/PsychophysicsData';
%     eval(['save ',[savepath,'/',filename],' hits misses FAs CRs']);
end

ShowCursor;
nhits = size(hits,1);
nmisses = size(misses,1);
nFAs = size(FAs,1);
nCRs = size(CRs,1);

trialtypes = {'hits','FAs','misses','CRs'};
figure;
tvect = [0:NFRAMESPED-1]/hz*1000;  % time in ms
tvect = tvect-tvect(NFRAMETESTSTART-1);
for i = 1:length(trialtypes)
    eval(['n=size(',trialtypes{i},',1);']);
    eval(['m=mean(',trialtypes{i},');']);
    eval(['v=var(',trialtypes{i},');']);
    ci_m = norminv([.025 .975],0,NOISEAMP/sqrt(n));
    ci_v = (NOISEAMP^2)/n*chi2inv([.025 .975],n);
    subplot(2,4,i); hold on;
    plot(tvect, m);
    plot([tvect(NFRAMETESTSTART-1) tvect(NFRAMETESTSTART+NFRAMESTEST-1)],[0 0],'m-','Linewidth',3);
    plot([tvect(1) tvect(end)], [ci_m(1) ci_m(1)],'k:');
    plot([tvect(1) tvect(end)], [ci_m(2) ci_m(2)],'k:');
    title([trialtypes{i},': ',num2str(n)]);
    set(gca,'Xlim',[tvect(1) tvect(end)]);
    subplot(2,4,i+4); hold on;
    plot(tvect, v);
    plot([tvect(1) tvect(end)], [ci_v(1) ci_v(1)],'k:');
    plot([tvect(1) tvect(end)], [ci_v(2) ci_v(2)],'k:');
    plot([tvect(NFRAMETESTSTART-1) tvect(NFRAMETESTSTART+NFRAMESTEST-1)],[NOISEAMP^2 NOISEAMP^2],'m-','Linewidth',3);
    set(gca,'Xlim',[tvect(1) tvect(end)]);
end

nhits = size(hits,1);
nmisses = size(misses,1);
nFAs = size(FAs,1);
nCRs = size(CRs,1);
pFAs = nFAs./(nFAs+nCRs);
phits = nhits./(nhits+nmisses);
z1 = norminv(1-pFAs,0,1);
z2 = norminv(1-phits,0,1);
dprime = z1-z2

[v,d] = eig(cov(hits));
