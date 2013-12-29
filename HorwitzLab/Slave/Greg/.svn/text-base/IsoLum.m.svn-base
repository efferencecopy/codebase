% This is going to a script for adjusting L to M contrast ratios until we
% attain isoluminance.  I'm mostly just doing this to make sure I can do it
% (and to provide a demo).  I'm going to use the quadrature motion stimulus
% of Cavanaugh and Anstis (1983).
%
% GDLH 6/9/07

% Constants
NPIXROWS = 50;
NPIXCOLS = 20;
NCYCLES = 3;
NITER = 5;
CONSTANTGUN = 3;
VARIABLEGUN = 1;
CONSTANTGUNVAL = 1;

% Loading cal file so we know if we're using the Bits++ or not
CALPATH = '/Monitor Calibration/Monitor data/Dell 4';
CALFILE = 'Dell4BitsCal';
CALIDX = 1;
load ([CALPATH,'/',CALFILE]);
cal = cals{CALIDX};

if (strfind(cal.describe.comment,'Colour'))
    COLOURMODE = 1;
    MAXDAC = 65535;
else
    COLOURMODE = 0;
    MAXDAC = 255;
end
MAXVAL = round(1*MAXDAC);
MINVAL = round(.4*MAXDAC);

% Preparing the M matrix and the inverse gamma table
load ('T_cones_smj.mat');
fundamentals = T_cones_smj;
cal = cals{CALIDX};
P_device = SplineSpd(SToWls(cal.S_device), cal.P_device, SToWls([380 5 81]));
M = fundamentals*P_device;
rgbAmbient = FindModelWeights(cal.P_ambient, cal.P_device);
invGammaTable = InvertGamma(cal, 1);
gammaTable = zeros(MAXDAC+1,3);
if (COLOURMODE)
    for i = 1:3
        gammaTable(:,i) = interp1([0:255], cal.gammaTable(:,i), linspace(0,255,65536))';
    end
else
    gammaTable = cal.gammaTable;
end

% calculating initial RGBs
bkgndrgb = [cal.gammaTable(round(256*cal.bgColor(1)),1);...
            cal.gammaTable(round(256*cal.bgColor(2)),2);...
            cal.gammaTable(round(256*cal.bgColor(3)),3)];

rgbs = [1.5*bkgndrgb'; 1.5*bkgndrgb'];
RGBs = zeros(4,3);
for i = 1:2
    RGBs(i,:) = round(MAXDAC*[invGammaTable(round(rgbs(i,1)*65535)+1,1);...
                              invGammaTable(round(rgbs(i,2)*65535)+1,2);...
                              invGammaTable(round(rgbs(i,3)*65535)+1,3)]);
end

RGBs(3,CONSTANTGUN) = MAXDAC*invGammaTable((CONSTANTGUNVAL*65535)+1,CONSTANTGUN);
RGBs(4,VARIABLEGUN) = MAXDAC;

rgb1mat = repmat(permute(RGBs(1,:),[3 1 2]),NPIXROWS, floor(NPIXCOLS/2));
rgb2mat = repmat(permute(RGBs(2,:),[3 1 2]),NPIXROWS, floor(NPIXCOLS/2));
frame1 = repmat([rgb1mat, rgb2mat], 1,NCYCLES);
frame3 = repmat([rgb2mat, rgb1mat], 1,NCYCLES);
if (COLOURMODE)
    frame1 = TranslateToColourMode(frame1);
    frame3 = TranslateToColourMode(frame3);
end

% Opening the window
w = Screen('OpenWindow',0,255*cal.bgColor);
HideCursor;
[winWidth, winHeight] = Screen('WindowSize',w);
Screen('LoadNormalizedGammaTable',0,repmat(linspace(0,1,256)',1,3));

% Main loop
betas = [];
for iter = 1:NITER
    i = 0;
    while(1)
        [x,y,buttons] = GetMouse(w);  
        if (any(buttons))
           while (any(buttons)) % wait for release
  	           [x,y,buttons] = GetMouse;
           end
           break 
        end
        val = round(((x+1)/winWidth)*(MAXVAL-MINVAL)+MINVAL);
        RGBs(4,VARIABLEGUN) = val;
        if (COLOURMODE)
            highbyte = floor(RGBs(3,:)/256);
            lowbyte = rem(RGBs(3,:),256);
            rgb1mat = repmat(permute([highbyte; lowbyte],[3 1 2]),NPIXROWS, floor(NPIXCOLS/2));
            highbyte = floor(RGBs(4,:)/256);
            lowbyte = rem(RGBs(4,:),256);
            rgb2mat = repmat(permute([highbyte; lowbyte],[3 1 2]),NPIXROWS, floor(NPIXCOLS/2));
        else
            rgb1mat = repmat(permute(RGBs(3,:),[3 1 2]),NPIXROWS, NPIXCOLS);
            rgb2mat = repmat(permute(RGBs(4,:),[3 1 2]),NPIXROWS, NPIXCOLS);
        end
        frame2 = [rgb1mat(:, [1:floor(NPIXCOLS/2)],:),...
                  repmat([rgb2mat, rgb1mat], 1,NCYCLES-1),...
                  rgb2mat,...
                  rgb1mat(:,[1:end-ceil(NPIXCOLS/2)],:)];
        frame4 = [rgb2mat(:, [1:floor(NPIXCOLS/2)],:),...
                  repmat([rgb1mat, rgb2mat], 1,NCYCLES-1),...
                  rgb1mat,...
                  rgb2mat(:,[1:end-ceil(NPIXCOLS/2)],:)];
        frame = mod(i,4)+1;
        tex=Screen('MakeTexture', w, eval(['frame',num2str(frame)]));
        Screen('DrawTexture',w,tex,[],[],[],0);
        Screen('Flip',w);
        %pause(.01);
        Screen('Close',tex);
        i = i +1;
    end
    int1 = gammaTable(RGBs(3,CONSTANTGUN)+1,CONSTANTGUN);
    int2 = gammaTable(RGBs(4,VARIABLEGUN)+1,VARIABLEGUN);
    LMS1 = M(:,CONSTANTGUN)*int1;
    LMS2 = M(:,VARIABLEGUN)*int2;
    beta = (LMS2(1)-LMS1(1))/(LMS1(2)-LMS2(2));
    betas = [betas; beta];
    SetMouse(unidrnd(winWidth+1),1,w);
    ShowCursor;
    HideCursor;
end
Screen('Close',w);
ShowCursor;
figure;
hist(betas);
mean(1./betas)

