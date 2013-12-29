% Section 1
% Measuring cone contrasts of cone-isolating stimuli created
% with 2 deg fundamentals and analyzed with 10 deg fundamentals.
%
% Section 2
% Finding LMS values for use with  DTspot (which defaults to the SMJ 2
% degree fundamentals) 
% that can be entered in REX that will provide an S-cone isolation stimulus
% 
% according to the 10 deg fundamentals
%
% Section 3
% Trying to analyze a DTspot file with parameters appropriate for testing
% the difference between 2 and 10 deg S-cone fundamentals
% 
% Section 4
% Comparing performance (and microsaccadic suppression) across potential
% S-cone isolation color directions at high and low temporal frequency.
%
% Section 5
% Finding color directions in DTspot that you can use to modulate the green
% and blue phosphors individually.
%
% Section 6
% Analyzing macular pigment DTspot data (blue vs green phosphor flicker)
%
% Section 6.1
% Comparing green to blue threshold ratios at 15 and 25 Hz. Is there a
% difference?
%
% Section 7
% Computing luminance for blue and green phosphors from synthetic cone
% fundamentals. Ability to manipulate macular pigment density, lens
% density, self screening, L to M ratio.
%
% Section 7.1
% Creating a set of cone fundamentals from the Baylor/Schnapf suction 
% electrode recordings. This may have worked well. (3/2/12)
%
% Section 7.2
% Creating a bunch of synthetic cone fundamentals based on the Stockman and
% Sharpe 10 degree fundamentals.
%
% Section 7.25
% Creaing cone fundamentals starting with the Baylor et al. data.
%
% Section 7.3
% Do the synthetic fundamentals based on the Baylor/Boettner data eliminate
% the S-cone signals from L+M Neurothresh cells? This is the section where
% T_cones_synthgh1 is created.
% the S-cone signals from L+M Neurothresh cells? 
%
% Section 8
% Trying to get from the Baylor rod measurements to T_rods (the human
% scotoptic sensitivity function). If monkey and human lens densities are
% different, will DT_scot be sufficient to discriminate humans from
% monkeys.
% 
% Section 9 
% Comparing the Stockman, MacLeod fundamentals, Stockman and Sharpe, and the
% "synthetic" ones.
%
% Section 10 
% Comparing the predictions of monkey luminance from measurements from the
% literature.
%
% Section 11
% Generating a set of transmittance values for the Kodak 1 log unit neutral
% density filter. Tranmittance values from 400 to 700 nm are from the
% published numbers. Values from 380 to 400 and from 700 to 780 are 
% extrapolations based on a graph that accompany the numbers on the PDF of
% page 182 in the Kodak filter book.
%
% Section 12
% Playing with the Smith and Pokorny lens density as a function of age
% model.
%
% Section 13
% Getting the spectra of RGBs at scotopic threshold
%
% Section 14
% What blue/green ratio is predicted from the synthetic cone fundamentals
% (with an artificially low lens density) and realistic amounts of macular
% pigment (e.g. 0.35 at 460 nm in fovea)?
%%
% Section 1
% If we make lights based on the 2 deg fundamentals, but the 10 deg
% fundamentals are the correct ones, how badly have we blown it?
load('C:\Documents and Settings\horwitzlab\Desktop\MatlabCode\Slave\Monitor Calibration\Monitor data\Dell 4\Dell4BitsCal.mat')
cal = cals{end};
load 'T_cones_smj10'
load 'T_cones_smj'
P_device = SplineSpd(SToWls(cal.S_device), cal.P_device, SToWls([380 5 81]));
M1 = T_cones_smj * P_device;
M2 = T_cones_smj10 * P_device;

bkgndRGB = round(cal.bgColor * 255) +1;
bkgndrgb  = [cal.gammaTable(bkgndRGB(1),1), cal.gammaTable(bkgndRGB(2),2), cal.gammaTable(bkgndRGB(3),3)];


colors = [.09 -.09 0;
    0 0 0.6364;
    0.0636   -0.0636    -0.4500;
    0.0636   -0.0636    0.4500;
    ];


stimrgb = inv(M1)*(repmat(bkgndlms1,1,size(colors,1)).*(1+colors'));

% check
bkgndlms1 = M1 * bkgndrgb';
twodegcc = (M1*stimrgb-repmat(bkgndlms1,1,4))./repmat(bkgndlms1,1,4)
% good

bkgndlms2 = M2 * bkgndrgb';
tendegcc = (M2*stimrgb-repmat(bkgndlms2,1,4))./repmat(bkgndlms2,1,4)

xy1 = [1 -1 0; 0 0 1]*twodegcc;
xy2 = [1 -1 0; 0 0 1]*tendegcc;

plot(xy2(1,:), xy2(2,:),'m*')

%%
% Section 2
% Finding LMS values for DTspot that can be entered in REX that will 
% provide an S-cone isolation stimulus according to the 10 deg fundamentals
load('C:\Documents and Settings\horwitzlab\Desktop\MatlabCode\Slave\Monitor Calibration\Monitor data\Dell 4\Dell4BitsCal.mat')
cal = cals{end};
load 'T_cones_smj10'
load 'T_cones_smj'
P_device = SplineSpd(SToWls(cal.S_device), cal.P_device, SToWls([380 5 81]));
M1 = T_cones_smj * P_device;
M2 = T_cones_smj10 * P_device;

bkgndRGB = round(cal.bgColor * 255) +1;
bkgndrgb  = [cal.gammaTable(bkgndRGB(1),1), cal.gammaTable(bkgndRGB(2),2), cal.gammaTable(bkgndRGB(3),3)];
bkgndlms1 = M1 * bkgndrgb';
bkgndlms2 = M2 * bkgndrgb';
% What we are using for and S-cone isolating direction
stimrgb1 = inv(M1)*(bkgndlms1.*(1+[0 0 1]'));
% What we want
stimrgb2 = inv(M2)*(bkgndlms2.*(1+[0 0 1]'));
stimrgb3 = 2*stimrgb2-stimrgb1;  % Even *less* macular pigment than 10°

% Now I need to find some LMSs which, when I put them through M1,
% give me stimrgb2
LMS = ((M1 * stimrgb2)./bkgndlms1)-1;

% Checking
stimrgb0 = inv(M1)*(bkgndlms1.*(1+LMS));

% Now I need to find some LMSs which, when I put them through M1,
% give me stimrgb3 (even les macular pigment than 10°)
LMS = ((M1 * stimrgb3)./bkgndlms1)-1;



figure; axes; hold on;
plot(stimrgb1,'k-');
plot(stimrgb2,'b-');
plot(stimrgb0,'m-');

%%%%%
% OK, now doing it in reverse
% Here's the RGBs that we'll use
stimrgb1 = inv(M1)*(bkgndlms1.*(1+LMS));
%Here's the LMS that corresponds to, using the 10 degree fundamentals
stimlms2 = M2*stimrgb1;
% And here's the cone contrast of that stimulus, again using the 10 deg
% fundamentals
(stimlms2-bkgndlms2)./bkgndlms2

% Good.


%%
% Section 3
% Trying to analyze a DTspot file with parameters appropriate for testing
% the difference between 2 and 10 deg S-cone fundamentals
stro = nex2stro;
unique(stro.trial(:,21))  % Should be 58 (spatial period)

load('/Users/greghorwitz/Desktop/MatlabCode/Slave/Monitor Calibration/Monitor data/Dell 4/Dell4BitsCal.mat')
cal = cals{end};
load 'T_cones_smj10'
load 'T_cones_smj'
P_device = SplineSpd(SToWls(cal.S_device), cal.P_device, SToWls([380 5 81]));
M1 = T_cones_smj * P_device;
M2 = T_cones_smj10 * P_device;

% Looking at the RGBs
exptColors = reshape(stro.sum.exptParams.RF_colors, 3, 3)';
exptColors = exptColors.*repmat(sign(sum(sign(exptColors+eps),2)),1,size(exptColors,2)); %standardizing color directions
noColor = sum(abs(exptColors), 2) == 0;
exptColors(noColor,:) = [];
whichcol = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'color_dir'));
norms = sqrt(sum(exptColors.^2, 2));
exptColors = exptColors./repmat(norms, 1, 3);
bkgndrgb = [stro.sum.exptParams.bkgnd_r, stro.sum.exptParams.bkgnd_g, stro.sum.exptParams.bkgnd_b];
M = reshape(stro.sum.exptParams.m_mtx, 3, 3);
bkgndlms = M * bkgndrgb';
x = 0:255; %the normal range of the gamma look up table
xx = linspace(0, 255, 2^16); %the desired quantization of the gammaTable
g1 = reshape(stro.sum.exptParams.gamma_table, 256, 3);
gammaTable = [spline(x, g1(:,1), xx)', spline(x, g1(:,2), xx)', spline(x, g1(:,3), xx)'];

RGB = stro.trial(:,14:16)+1;
rgb = [gammaTable(RGB(:,1), 1), gammaTable(RGB(:,2), 2), gammaTable(RGB(:,3), 3)];
LMS = (M*rgb')-repmat(bkgndlms,1,size(rgb,1));
ccs = LMS' ./ repmat(bkgndlms',size(rgb,1),1);

figure; axes; hold on;
plot3(rgb(whichcol == 1,1),rgb(whichcol == 1,2),rgb(whichcol == 1,3),'m.');
plot3(rgb(whichcol == 2,1),rgb(whichcol == 2,2),rgb(whichcol == 2,3),'k.');

% Plugging in the 10 degree fundamentals
bkgndlms2 = M2 * bkgndrgb';
LMS2 = M2*rgb';
ccs2 = (LMS2'- repmat(bkgndlms2',size(rgb,1),1))./ repmat(bkgndlms2',size(rgb,1),1);

figure;
subplot(3,2,1)
plot(ccs(whichcol == 1,:))
subplot(3,2,2)
plot(ccs2(whichcol == 1,:))
subplot(3,2,3)
plot(ccs(whichcol == 2,:))
subplot(3,2,4)
plot(ccs2(whichcol == 2,:))
subplot(3,2,5)
plot(ccs(whichcol == 3,:))
subplot(3,2,6)
plot(ccs2(whichcol == 3,:))
%%
% Section 4
% Comparing performance (and microsaccadic suppression) across potential
% S-cone isolation color directions at high and low temporal frequency.

filenames = fnamesFromTxt();
offset = [.167 -.167];
allcolors = [];
data = [];
for a = 1:size(filenames,1)
    stro = nex2stro(findfile(filenames(a,:)));
    sacstats = getSacData(stro);
    close;
    ntrials = size(stro.trial,1);
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_on'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_off'));
    correct = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));

    [thresholds, colorDirs, sfs] = DTquestUnpack(stro,'mode');
    % close;
    speed = stro.trial(:,strcmp(stro.sum.trialFields(1,:), 'gabor_speed'));
    whichcol = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'color_dir'));

    colors = reshape(stro.sum.exptParams.RF_colors,[3 3]);
    if (isempty(allcolors))
        allcolors = colors;
    else
        if (~all(colors == allcolors))
            error('Disagreement in colors across experiments')
        end
    end
    
    % Looping over trials
    saccadeoccurred = nan*ones(ntrials,1);
    for j = 1:ntrials
        st = sacstats.starttimes{j};
        Lsac = (st > stimon_t(j)+offset(1)) & (st < stimoff_t(j)+offset(2));
        if (any(Lsac))
            saccadeoccurred(j) = 1;
        else
            saccadeoccurred(j) = 0;
        end
    end
    trialtypes = unique(whichcol);
    tmpdata = zeros(size(trialtypes,1),1);
    for j = 1:size(trialtypes,1)
        L = whichcol == trialtypes(j);
        r = corrcoef(correct(L), saccadeoccurred(L));
        tmpdata(j) = r(1,2);  % Could be Nan if one of the two vectors has var = 0
    end
    
    
    data = [data; thresholds', tmpdata', speed(1)];
end
speeds = unique(data(:,7));

% All possible threshold ratios (3 Hz/not 3 Hz, within color directions)
% Too bad these things aren't independent
figure; axes; hold on;
for i = 1:size(allcolors,2)
    L = data(:,7) == 3;
    a = find(L);
    b = find(~L);
    designmat = fullfact([sum(L), sum(~L)]);
    for j = 1:size(designmat,1)
        plot(i,data(a(designmat(j,1)),i)./data(b(designmat(j,2)),i),'k.')
    end
end

% Now looking at paired threshold ratios.  Lame, but whatever.
L = data(:,7) == 3;
a = find(L);
b = find(~L);
a = a(1:min(length(a),length(b)))
b = b(1:min(length(a),length(b)))

ratios = data(a,[1 2 3])./data(b,[1 2 3]);
figure; subplot(2,1,1); hold on;
plot(ratios(:,1),'b.-')
plot(ratios(:,2),'y.-')
legend({'2°','10°'});
[h,p] = ttest(ratios(:,1)-ratios(:,2))
subplot(2,1,2);
hist(ratios(:,1)-ratios(:,2))
title(['p = ',num2str(p)]);


% Getting summary statistics on the behavioral performance
out = [];
for i = 1:size(allcolors,2)
    for j = 1:length(speeds)
        spd = speeds(j);
        L = data(:,7) == spd;
        out(j,:,i) = [mean(data(L,i)) std(data(L,i)) sum(L)];
    end
end

% Making a plot of the threshold ratios
% Using the formulas from Kamerund 1978
figure; axes; hold on;
for i = 1:size(out,3)
    muX = out(1,1,i);
    sigmaX = out(1,2,i)./sqrt(out(1,3,i));
    muY = out(2,1,i);
    sigmaY = out(2,2,i)./sqrt(out(2,3,i));
    
    z = linspace(0.2,.6,1000);
    w = (sigmaY/sigmaX)*z;  % w is a scaled version of z
    s = 1./sqrt(w.^2+1);
    k = (muX/sigmaX.*w+muY/sigmaY).*s.^2;
    M = -.5*(muY/sigmaY.*w-muX/sigmaX).^2.*s.^2;
    Q = k.*s.*sqrt(2*pi).*(1-2.*normcdf(-k./s))+(2*s.^2.*exp(-k.^2./(2*s.^2)));
    
    fg = 1/(2*pi).*Q.*exp(M);
    fz = (sigmaY/sigmaX)*fg;
   % plot(z,fz./max(fz),'m-');
    
    % CDF
    Fz = cumsum(fz);
    lastidx = find(diff(Fz) == 0,1,'first');
    if (~isempty(lastidx))
        Fz = Fz(1:lastidx);
        z = z(1:lastidx);
    end
    Fz = Fz./Fz(end);
    CI = interp1(Fz,z,[.05 .95]);
    plot(i,muX/muY,'k*')
    plot([i,i],CI);
end
set(gca,'Xlim',[0 4])


% Looking at the microsaccadic suppression effect
figure;
for j = 1:length(speeds)
    spd = speeds(j);
    L = data(:,7) == spd;
    subplot(length(speeds),1,j)
    boxplot(data(L,[4 5 6]));
    hold on;
    title(['speed: ',num2str(spd)]);
    for i = 1:3
        [h,p] = ttest(data(L,i+3));
        if (h)
            plot(i,0,'m*');
        end
    end
end

xlabel('Color direction');
ylabel('Microsaccade effect');
% Comparing microsaccade effects between S-cone isolating
% directions at 3 Hz.
L = data(:,7) == 3;
[h,p] = ttest2(data(L,4),data(L,5))

%%
% Section 5
% Scripts for estimating macular pigmentation from DTspot data
% Finding color directions in DTspot that you can use to modulate the green
% and blue phosphors individually.

load('C:\Documents and Settings\ghorwitz\Desktop\MatlabCode\Slave\Monitor Calibration\Monitor data\Dell 4\Dell4BitsCal.mat')
cal = cals{end};
load 'T_cones_smj10'  % The default cone fundamentals
P_device = SplineSpd(SToWls(cal.S_device), cal.P_device, SToWls([380 5 81]));
M = T_cones_smj10 * P_device;

bkgndRGB = round(cal.bgColor * 255) +1;
bkgndrgb  = [cal.gammaTable(bkgndRGB(1),1), cal.gammaTable(bkgndRGB(2),2), cal.gammaTable(bkgndRGB(3),3)];

bluergb = [bkgndrgb(1) bkgndrgb(2) 1];
greenrgb = [bkgndrgb(1) 1 bkgndrgb(3)];
redrgb = [1 bkgndrgb(2) bkgndrgb(3)];

bluelms = M*bluergb';
greenlms = M*greenrgb';
redlms = M*redrgb';
bkgndlms = M * bkgndrgb';
% Cone contrasts
redcc = (redlms-bkgndlms)./bkgndlms;
greencc = (greenlms-bkgndlms)./bkgndlms;
bluecc = (bluelms-bkgndlms)./bkgndlms;

% Checking
% stimrgb = inv(M)*(bkgndlms.*(1+bluecc))
% stimrgb = inv(M)*(bkgndlms.*(1+greencc))
% stimrgb = inv(M)*(bkgndlms.*(1+redcc))
%bkgndrgb
% Good.

% making them unit vectors for no good reason
% Don't worry about these numbers not matching the ones above
% the ones above
%redunit = redcc./norm(redcc)
%greenunit = greencc./norm(greencc)
%blueunit = bluecc./norm(bluecc)
%%
% Section 6
% Analyzing macular pigment DTspot data (blue vs green phosphor flicker)
filenames = fnamesFromTxt2();
data = [];
% Expressing thresholds in luminance contrast
% Warning: 1:1 L to M ratio assumed below
load T_xyz1931;
Vlambda = T_xyz1931(2,:);
for fileidx = 1:size(filenames,1)
    stro = nex2stro(findfile(char(filenames{fileidx,:})));
    [thresh, colorDirs, sfs] = DTquestUnpack(stro, 'mode');
    normcolordirs= colorDirs./repmat(sqrt(sum(colorDirs.^2,2)),1,3);
    guns = reshape(stro.sum.exptParams.mon_spect,81,3);
    % close(gcf); % for when you just want the summary figures
    bkgndrgb = [stro.sum.exptParams.bkgnd_r stro.sum.exptParams.bkgnd_g stro.sum.exptParams.bkgnd_b];
    M = reshape(stro.sum.exptParams.m_mtx,3,3);
    % Quantifying thresholds in 2 deg luminance contrast units
    bkgndlum = bkgndrgb*guns'*Vlambda';
    bkgndlms = M*bkgndrgb';
    threshrgb = inv(M)*(repmat(bkgndlms',2,1).*(1+repmat(thresh/100,[1 3]).*normcolordirs([1 2],:)))';
%    thresholds = (threshrgb'*guns'*Vlambda')./bkgndlum;
    deltagun = threshrgb-repmat(bkgndrgb',1,2)
    % Quantifying thresholds in 2 deg luminance units (should convert to
    % cd/m^2?)
    thresholds = [];
    for i = 1:size(deltagun,2)
        thresholds(i) =  (guns*deltagun(:,i))'*Vlambda';
        thresholds(i) =   thresholds(i)./bkgndlum; % dividing by bkgndlum to make the 
        % units of the absolute thresholds eaasier to think about. This is
        % irrelevant to the blum/glum ratio since bkgndlum will cancel.
    end
    set(gcf,'Name',stro.sum.fileName);
    ecc = stro.sum.exptParams.rf_x;
    data = [data; abs(ecc/10) thresholds];
%    data = [data; ecc/10 thresh'];
end

% Overall performance
figure;
subplot(2,1,1);
[r,p] = corr([1:size(data,1)]',data(:,2),'Type','Pearson');
plot(data(:,2),'k.')
ylabel('Green threshold');
title(['r = ',num2str(r,3),' p = ',num2str(p,3)]);
subplot(2,1,2);
plot(data(:,3),'k.')
ylabel('Blue threshold');
[r,p] = corr([1:size(data,1)]',data(:,3),'Type','Pearson');
title(['r = ',num2str(r,3),' p = ',num2str(p,3)]);
xlabel('File idx');

figure;
subplot(2,1,1); hold on;
plot(data(:,1),data(:,2)*100,'g.')
plot(data(:,1),data(:,3)*100,'b.')
%plot(data(:,1),data(:,2),'g.')
%plot(data(:,1),data(:,3),'b.')
ylabel({'Detection threshold','(% contrast)'});
xlabel('Eccentricity (deg)');

pp_g  = csaps(data(:,1),log10(data(:,2)*100),.5);
pp_b  = csaps(data(:,1),log10(data(:,3)*100),.5);
%pp_g  = csaps(data(:,1),log10(data(:,2)),.5);
%pp_b  = csaps(data(:,1),log10(data(:,3)),.5);
x = linspace(0,8,100);
plot(x,10.^(fnval(x,pp_g)),'g-','Linewidth',2)
plot(x,10.^(fnval(x,pp_b)),'b-','Linewidth',2)
% Doing the fits on the logs doesn't seem to make any difference at all.

%regcoef_g = regress(data(:,2),[data(:,1), ones(size(data,1),1)]);
%regcoef_b = regress(data(:,3),[data(:,1), ones(size(data,1),1)]);
%plot([0 8],[0 8]*regcoef_g(1)+regcoef_g(2),'g-');
%plot([0 8],[0 8]*regcoef_b(1)+regcoef_b(2),'b-');

subplot(2,1,2); hold on;
%plot(data(:,1),data(:,3)./data(:,2),'k.');
plot(data(:,1),data(:,3)./data(:,2),'k.');
%[r,p] = corrcoef ([data(:,1),data(:,3)./data(:,2)])
[r,p] = corrcoef ([data(:,1),data(:,3)./data(:,2)])
ylabel('Bthresh/Gthresh'); xlabel('Eccentricity (deg)');
x = linspace(0,8,100);
plot(x,10.^fnval(x,pp_b)./10.^fnval(x,pp_g),'k-','LineWidth',2)
%plot(linspace(0,8,100),(linspace(0,8,100)*regcoef_b(1)+regcoef_b(2))./(linspace(0,8,100)*regcoef_g(1)+regcoef_g(2)),':');

% % testing "gun isolation"
% M = reshape(stro.sum.exptParams.m_mtx,3,3);
% bkgndrgb = [stro.sum.exptParams.bkgnd_r stro.sum.exptParams.bkgnd_g stro.sum.exptP
%ecc = input('new eccentricity (in deg)? ');
%disp(['Set t1 and t2 dist to ',num2str(20/ecc)]);
%disp(['Set RF_X to ',num2str(ecc*10)]);

% % testing "gun isolation"
% M = reshape(stro.sum.exptParams.m_mtx,3,3);
% bkgndrgb = [stro.sum.exptParams.bkgnd_r stro.sum.exptParams.bkgnd_g stro.sum.exptParams.bkgnd_b];
% colordirs = reshape(stro.sum.exptParams.RF_colors,3,3);
% bkgndlms = M*bkgndrgb';
% stimrgb1 = inv(M)*(bkgndlms.*(colordirs(:,1)+1));
% stimrgb2 = inv(M)*(bkgndlms.*(colordirs(:,2)+1));
% [stimrgb1 stimrgb2 bkgndrgb']


%%
% Section 6.1
% Does flicker frequency influence blue to green threshold ratio?

filenames = fnamesFromTxt2();
data = [];
load T_xyz1931;
Vlambda = T_xyz1931(2,:);
for fileidx = 1:size(filenames,1)
    whichfiles = filenames{fileidx};
    for i = 1:length(whichfiles)
        stro = nex2stro(findfile(char(whichfiles{i})));
        TFidx = find(strcmp(stro.sum.trialFields(1,:),'gabor_speed'));
        [thresh, colorDirs, sfs] = DTquestUnpack(stro, 'mode');
        % close(gcf); % for when you just want the summary figures
        normcolordirs= colorDirs./repmat(sqrt(sum(colorDirs.^2,2)),1,3);
        guns = reshape(stro.sum.exptParams.mon_spect,81,3);
        bkgndrgb = [stro.sum.exptParams.bkgnd_r stro.sum.exptParams.bkgnd_g stro.sum.exptParams.bkgnd_b];
        M = reshape(stro.sum.exptParams.m_mtx,3,3);
        bkgndlms = M*bkgndrgb';
        bkgndlum = sum(bkgndlms([1 2]));
        % Quantifying thresholds in 2 deg luminance contrast units
        bkgndlum = bkgndrgb*guns'*Vlambda';
        threshrgb = inv(M)*(repmat(bkgndlms',2,1).*(1+repmat(thresh/100,[1 3]).*normcolordirs([1 2],:)))';
        thresholds = (threshrgb'*guns'*Vlambda')./bkgndlum;
        
        set(gcf,'Name',stro.sum.fileName);
        ecc = stro.sum.exptParams.rf_x;
        data = [data; ecc/10 thresholds' stro.trial(1,TFidx)];
    end
end

figure; axes; hold on;
h = [];
for i = [15 25]
    L = data(:,4) == i;
    if (i == 15)
        h(1) = plot(data(L,1),(data(L,3)-1)./(data(L,2)-1),'ko','MarkerSize',5);
    else
        h(2) = plot(data(L,1),(data(L,3)-1)./(data(L,2)-1),'m*','MarkerSize',5); 
    end    
end
legend(h,{'15 Hz','25 Hz'})
% A little regression analysis
whichecc = unique(data(:,1))';
p = []; y =[]; X = [];
for ecc = whichecc
    Lecc = data(:,1) == ecc;
    LTF15 = data(:,4) == 15;
    ratios15 = (data(Lecc&LTF15,3)-1)./(data(Lecc&LTF15,2)-1);
    LTF25 = data(:,4) == 25;
    ratios25 = (data(Lecc&LTF25,3)-1)./(data(Lecc&LTF25,2)-1);
    [h,p(ecc == whichecc)] = ttest(ratios15,ratios25);
    n15 = length(ratios15);
    n25 = length(ratios25);
    y = [y; ratios15; ratios25];
    X = [X; [repmat(ecc,n15+n25,1), [repmat(15,n15,1);repmat(25,n25,1)] ones(n15+n25,1)]];
end
% At which eccentricities is there a difference between 15 and 25 Hz? 
whichecc(p<0.05)
% y = b1*eccentricity+b2*frequency+b3;
[b,bint,r,rint,stats] = regress(y,X);
disp(['ecc'; 'TF ';'int'])
disp(bint)

[p,antab] = anovan(y,{X(:,1),X(:,2)},'model','interaction','varnames',strvcat('Ecc','TF'));

%%
% Section 7
% Computing luminance from cone fundamentals from which macular pigment
% (and lens?) have been removed (adjusted). Continuing from the previous
% cell. 
load ('den_mac_ws');
den_mac_ws(4) = 0.0425; %  minor tweak made by SMJ 1993: density of 0.0425 @ 395 nm
% dens_lens_smj: This is SMJ Table 7 *1.16 (small pupil adjustment)
% CVRL says "These values are for a small pupil, for an open pupil divide by 1.16."
den_lens_smj =...  %
[380 3.2;... % linear extrapolation
385  2.8;... % linear extrapolation
390 2.4;...
395	1.985;...
400	1.62;...
405	1.308;...
410	1.05;...
415	0.84;...
420	0.675;...
425	0.557;...
430	0.468;...
435	0.393;...
440	0.335;...
445	0.29;...
450	0.26;...
455	0.24;...
460	0.225;...
465	0.215;...
470	0.203;...
475	0.191;...
480	0.18;...
485	0.168;...
490	0.162;...
495	0.151;...
500	0.145;...
505	0.133;...
510	0.128;...
515	0.122;...
520	0.116;...
525	0.11;...
530	0.104;...
535	0.099;...
540	0.093;...
545	0.087;...
550	0.081;...
555	0.075;...
560	0.07;...
565	0.064;...
570	0.058;...
575	0.052;...
580	0.046;...
585	0.041;...
590	0.036;...
595	0.031;...
600	0.028;...
605	0.024;...
610	0.021;...
615	0.017;...
620	0.014;...
625	0.012;...
630	0.009;...
635	0.007;...
640	0.005;...
645	0.003;...
650	0.002;...
655	0.001;...
[[660:5:780]' zeros(25,1)]];

% Some constants
macpig = den_mac_ws;
lens = den_lens_smj(:,2)./1.16;  % Dividing by 1.16 to make it "open pupil". Thats what SMJ did.
macpigtransmittance = 1./(10.^(macpig*.28));
lenstransmittance = 1./(10.^(lens*1.28));
opticaldensity = 0.3;
LVLAMBDA = 1; % L contribution to vlambda

% Computing predictions from vlambda
% Just using the default fundamentals for now.
% Assuming lum is L+M
%M = reshape(stro.sum.exptParams.m_mtx,3,3);
spds = reshape(stro.sum.exptParams.mon_spect,81,3);

load('T_cones_smj')
M = T_cones_smj*spds;
bkgndlms = M*bkgndrgb';
rgbdirs=  inv(M)*(repmat(bkgndlms',2,1).*(1+colorDirs))';
gpred = [LVLAMBDA 1]*M([1 2],2);
bpred = [LVLAMBDA 1]*M([1 2],3);
glum = Vlambda*guns(:,2);
blum = Vlambda*guns(:,3);
synthpred = (gpred./bpred)./(glum./blum)
subplot(2,1,2);
plot([0 8], synthpred*[1 1],'k:');
text(max(data(:,1)),synthpred, '2 deg');

load('T_cones_smj10')
M = T_cones_smj10*spds;
bkgndlms = M*bkgndrgb';
rgbdirs=  inv(M)*(repmat(bkgndlms',2,1).*(1+colorDirs))';
gpred = [LVLAMBDA 1]*M([1 2],2);
bpred = [LVLAMBDA 1]*M([1 2],3);
synthpred = (gpred./bpred)./(glum./blum)
subplot(2,1,2);
plot([0 8], synthpred*[1 1],'k:');
text(max(data(:,1)),synthpred, '10 deg','FontSize',9)

% ****************
% Undoing preretinal filters (a la SJM 1993)
% ****************
fund10 = T_cones_smj10./repmat(max(T_cones_smj10,[],2),1,size(T_cones_smj10,2));
absorptance = fund10'./repmat(macpigtransmittance,1,3)./repmat(lenstransmittance,1,3);
absorptance = absorptance./repmat(max(absorptance),81,1);
actionspectra = -log10(1-absorptance*(1-10^-opticaldensity));
actionspectra = actionspectra/opticaldensity;

% Now creating a NEW set of fundamentals starting with the action spectra
% Change parameters here
macpigtransmittance = 1./(10.^(macpig*0));
%lenstransmittance = 1./(10.^(lens*1.22/lens(5)));  % lens = 1.22 at 400 nm as per Cottaris 2003
%lenstransmittance = 1./(10.^(lens*1));  % default
lenstransmittance = 1./(10.^(lens/lens(5)));  % lens = 1 at 400 nm as per Norren
% lens = 1 at 400 (.71 multiplier), macpig*0 works OK...
% lens multipliers of .7 to 1.2 are in the reasonable range (Norren)
opticaldensity = 0.3; % Greater optical density *decreases* the B to G lum thresh ratio

funds = 1-10.^(-actionspectra.*opticaldensity);
funds = funds.*repmat(macpigtransmittance,1,3).*repmat(lenstransmittance,1,3);
funds = funds./repmat(max(funds),81,1);
%funds = funds.*repmat(macpigtransmittance,1,3).*repmat(lenstransmittance,1,3);

% Sanity check. Find should equal fund10 if preretinal screening 
% parameters are set identically in the removal and re-addition processes.
%figure; axes; hold on;
%plot(fund10','k-');
%plot(funds,'m-');

M = funds'*spds;
gpred = [LVLAMBDA 1]*M([1 2],2);
bpred = [LVLAMBDA 1]*M([1 2],3);
synthpred = (gpred./bpred)./(glum./blum)
h(1)= plot([0 8], synthpred*[1 1],'k:');
h(2)= text(max(data(:,1)),synthpred, 'reduced lens','FontSize',9)

% Adding a bit of macpig = 0.5 GDLH
lenstransmittance = 1./(10.^(lens/lens(5)));  % lens = 1 at 400 nm as per Norren
macpigtransmittance = 1./(10.^(macpig*.5));
funds = 1-10.^(-actionspectra.*opticaldensity);
funds = funds.*repmat(macpigtransmittance,1,3).*repmat(lenstransmittance,1,3);
funds = funds./repmat(max(funds),81,1);
M = funds'*spds;
gpred = [LVLAMBDA 1]*M([1 2],2);
bpred = [LVLAMBDA 1]*M([1 2],3);
synthpred = (gpred./bpred)./(glum./blum)
h(1)= plot([0 8], synthpred*[1 1],'k:');
h(2)= text(max(data(:,1)),synthpred, ['MPpeak = ',num2str(log10(1./min(macpigtransmittance))),' red. lens'],'FontSize',9)


% Parameters that can be manipulated:
% Macular pigmentation
% Photopigment optical density
% Lens pigment density
% L to M ratio in luminance
%%
% Section 7.1
% Comparing the "action spectra" that I get from backing the pretinal
% filters out of the 10 deg fundamentals to those tabulated by Baylor et
% al. This is also a test to see whether I know what I'm doing when I 
% do the conversion from cone fundamentals to action spectra.
% The Baylor S-cone action spectrum is more sensitive to long wavelengths
% than the action spectrum derived from the SMJ fundamentals.

load ('den_mac_ws');
den_mac_ws(4) = 0.0425; %  minor tweak made by SMJ 1993: density of 0.0425 @ 395 nm
% dens_lens_smj: This is SMJ Table 7 *1.16 (small pupil adjustment)
% CVRL says "These values are for a small pupil, for an open pupil divide by 1.16."
den_lens_smj =...  %
[380 3.3;... % extrapolation to smooth action spectra
385  2.8;... % extrapolation to smooth action spectra
390 2.4;...
395	1.985;...
400	1.62;...
405	1.308;...
410	1.05;...
415	0.84;...
420	0.675;...
425	0.557;...
430	0.468;...
435	0.393;...
440	0.335;...
445	0.29;...
450	0.26;...
455	0.24;...
460	0.225;...
465	0.215;...
470	0.203;...
475	0.191;...
480	0.18;...
485	0.168;...
490	0.162;...
495	0.151;...
500	0.145;...
505	0.133;...
510	0.128;...
515	0.122;...
520	0.116;...
525	0.11;...
530	0.104;...
535	0.099;...
540	0.093;...
545	0.087;...
550	0.081;...
555	0.075;...
560	0.07;...
565	0.064;...
570	0.058;...
575	0.052;...
580	0.046;...
585	0.041;...
590	0.036;...
595	0.031;...
600	0.028;...
605	0.024;...
610	0.021;...
615	0.017;...
620	0.014;...
625	0.012;...
630	0.009;...
635	0.007;...
640	0.005;...
645	0.003;...
650	0.002;...
655	0.001;...
[[660:5:780]' zeros(25,1)]];

% Some constants
macpig = den_mac_ws;
lens = den_lens_smj(:,2)./1.16;  % Dividing by 1.16 to make it "open pupil". Thats what SMJ did.
macpigtransmittance = 1./(10.^(macpig*.28));
lenstransmittance = 1./(10.^(lens*1.28));
opticaldensity = 0.3;

load('T_cones_smj10')
% Undoing preretinal filters (a la SJM 1993)
fund10 = T_cones_smj10./repmat(max(T_cones_smj10,[],2),1,size(T_cones_smj10,2));
absorptance = fund10'./repmat(macpigtransmittance,1,3)./repmat(lenstransmittance,1,3);
absorptance = absorptance./repmat(max(absorptance),81,1);
actionspectra = -log10(1-absorptance*(1-10^-opticaldensity));
actionspectra_smj = actionspectra/opticaldensity;

% Comparing to the Stockman and Sharpe absorbances 

% Parametric fits to the Baylor data
a = [-5.2734 -87.403 1228.4 -3346.3 -5070.3 30881 -31607];
lmm = [561 531 430];
S = [];
for j = 1:length(lmm)
    tmp = log10((1./([380:5:780]/1000))*(lmm(j)/561));
    for i = 1:length(a)
        S(j,:,i) = a(i)*(tmp.^(i-1));
    end
end
actionspectra_baylor = sum(S,3);

% need to extrapolate the S-cone fundamental because the polynomial fit
% goes up again outside of the range of the data.
Lin = [380:5:780] >= 535 & [380:5:780] <= 575;
Lout = [380:5:780] > 575;
actionspectra_baylor(3,Lout) = interp1([535:5:575],actionspectra_baylor(3,Lin),[580:5:780],'linear','extrap')

%figure; axes; hold on;
%plot([380:5:780],actionspectra_smj);
%plot([380:5:780],10.^actionspectra_baylor,':')
%set(gca,'YScale','log','XLim',[380 780]);
%set(gca,'Ylim',[10^-9 1]);
%title('Action Spectra from SMJ (solid) and Baylor data (dotted)');

% -----------------------
% Now trying monkey data from Boettner 1967 (direct measurements)
% Are these a good set of fundamentals for monkeys? 
% Lens by itself is too little; cornea in addition works better
% Using Boet_all works well.
% -----------------------
load('T_cones_synthgh1');
% Let's see what the predicted B/G ratio is
M = T_cones_synthgh1'*guns;
gpred = [1 1]*M([1 2],2);
bpred = [1 1]*M([1 2],3);
load T_xyz1931;
Vlambda = T_xyz1931(2,:);
glum = Vlambda*guns(:,2);
blum = Vlambda*guns(:,3);

synthpred = (gpred./bpred)./(glum./blum)  % ratio of how much the guns hit the putative detection 
% mechanism divided by how much they hit a vlambda mechanism. This
% caluculation gives the same answer if we'd used (gpred./bpred) and 
% normalized gun intensities (background subtracted) for thresholds. 
% Suggests that mac pig ranges from a peak intensity of 0 to 0.2.
% So there appears to be no problem with units (threshold normalization).
h(1)= plot([0 8], synthpred*[1 1],'k:');
h(2)= text(max(data(:,1)),synthpred, 'Synthetic funds','FontSize',9)

%%
% Section 7.2
% Making a new prediction based on a modifed version of the SS fundamentals
% (with lower lens density). Also preparing the SS macular pigment and lens
% densities to be on the 380:5:780 lattice and saving them. This is also
% where we create T_cones_synthgh2.
% % 
%   load('den_lens_ssf');
%   load('macss_5.csv');
% %  den_lens_ss = 10.^(interp1([390:1:830],log10(den_lens_ssf),[380:5:780],'linear','extrap'))';
% ------------- Below didn't work. It led to weird looking action spectra -----------
% %  den_mac_ss = 10.^(interp1(macss_5(:,1),log10(macss_5(:,2)),[380:5:780],'linear','extrap'))';
%   den_lens_ss = interp1([390:1:830],den_lens_ssf,[380:5:780],'linear','extrap')';
% ------------- Don't interpolate in log space I guess.  -----------
%   den_mac_ss = interp1(macss_5(:,1),macss_5(:,2),[380:5:780],'linear','extrap')';
%   den_lens_ss(isnan(den_lens_ss)) = 0;
%   den_mac_ss(isnan(den_mac_ss)) = 0;
%   
%   save('den_lens_ss.mat','den_lens_ss');
%   save('den_mac_ss.mat','den_mac_ss');

load den_lens_ss; 
load den_mac_ss;

% Putting the SS 2 degs funds on the 380:5:780 sampling lattice
%load('T_cones_ss2');
%T_cones_ss2 = 10.^(interp1([390:1:830]',log10(T_cones_ss2'),[380:5:780],'linear','extrap'));
%T_cones_ss2(isnan(T_cones_ss2)) = 0;
%save T_cones_myss2 T_cones_ss2

% % Comparing lens densities
boet_wls = [320 340 360 380 400 420 440 460 480 500 550 600 650 700 750 800];
boet_lens = -log10([5 1 .5 .5 2 35 71 82 84.5 86.5 87.5 89 89.5 90 90.5 90.5]/100);
boet_lens = interp1(boet_wls,boet_lens,[380:5:780],'linear')';
boet_all = -log10([.8 .3 .2 .2 1 18 38.5 47 50.5 53.5 58.5 63.5 65.5 67.5 69.5 70.5]/100);
boet_all = interp1(boet_wls,boet_all,[380:5:780],'linear')';
boet_all = boet_all-boet_all(end); % shifted to make more similar to other lens spectra

figure; axes; hold on;
plot([380:5:780], boet_all,'b.-');
plot([380:5:780], den_lens_ss./den_lens_ss(5),'m.-');

% Removing macular pigment and lens from SS 10 deg fundamentals
load('T_cones_ss10');
T_cones_ss10 = 10.^(interp1([390:1:830]',log10(T_cones_ss10'),[380:5:780],'linear','extrap'));
T_cones_ss10(isnan(T_cones_ss10)) = 0;
%save T_cones_myss10 T_cones_ss10

macpigtransmittance = 1./(10.^(den_mac_ss./max(den_mac_ss)*0.095));
lenstransmittance = 1./(10.^(den_lens_ss));
absorptance = T_cones_ss10./repmat(macpigtransmittance,1,3)./repmat(lenstransmittance,1,3);
absorptance = absorptance./repmat(max(absorptance),81,1);

% Adding back a bit of lens density to see if we can make a reasonable
% prediction for the blue to green threshold ratio.
%lenstransmittance = 1./(10.^(boet_all./boet_all(5)));  % density of 1 at 400 nm
%lenstransmittance = 1./(10.^boet_all);  % way too high
lenstransmittance = 1./(10.^(den_lens_ss./den_lens_ss(5)*1)); % density of 1 at 400 nm
plot([380:5:780],-log10(lenstransmittance),'b.-');
synthfunds2 = absorptance.*repmat(lenstransmittance,1,3);
synthfunds2 = synthfunds2./repmat(max(synthfunds2),81,1);

% Finding the prediction of the gpred/bpred from these new fundamentals
load T_xyz1931;
Vlambda = T_xyz1931(2,:);
stro = nex2stro(findfile('S100511013'));
spds = reshape(stro.sum.exptParams.mon_spect,81,3);
bkgndrgb = [stro.sum.exptParams.bkgnd_r stro.sum.exptParams.bkgnd_g stro.sum.exptParams.bkgnd_b];
M = synthfunds2'*spds;
gpred = [1 1]*M([1 2],2);
bpred = [1 1]*M([1 2],3);
glum = Vlambda*spds(:,2);
blum = Vlambda*spds(:,3);
synthpred = (gpred./bpred)./(glum./blum)

T_cones_synthgh2= synthfunds2;
%save T_cones_synthgh2.mat T_cones_synthgh2;

% .. and a third set of fundamentals. SS 10 deg funds + boettner_all set to
% a density of 1 at 400 nm. 
lenstransmittance = 1./(10.^(boet_all./boet_all(5)));  % density of 1 at 400 nm
synthfunds3 = absorptance.*repmat(lenstransmittance,1,3);
synthfunds3 = synthfunds3./repmat(max(synthfunds3),81,1);
T_cones_synthgh3= synthfunds3;
save T_cones_synthgh3.mat T_cones_synthgh3;
M = synthfunds3'*spds;
gpred = [1 1]*M([1 2],2);
bpred = [1 1]*M([1 2],3);
synthpred = (gpred./bpred)./(glum./blum)

T_cones_synthgh3= synthfunds3;
%save T_cones_synthgh3.mat T_cones_synthgh3;

% ... and a fourth set of fundamentals with relatively little lens pigment
% (a la synthpred2) and S-cone action spectra nudged over a bit to longer
% wavelengths, consistent with Baylor et al.
tmpabsorp = absorptance;
tmpabsorp(:,3) = interp1([385:5:780]',tmpabsorp([1:end-1],3),[380:5:780]','linear','extrap');
lenstransmittance = 1./(10.^(den_lens_ss./den_lens_ss(5)*1)); % density of 1 at 400 nm
synthfunds4 = tmpabsorp.*repmat(lenstransmittance,1,3);
synthfunds4 = synthfunds4./repmat(max(synthfunds4),81,1);
T_cones_synthgh4= synthfunds4;
save T_cones_synthgh4.mat T_cones_synthgh4;

% ... and a fifth set using the Pokorny age-compensated lens density
age = 3;
den_lens_pokorny = [400 .600 1.00;
                     410 .510 .583;
                     420 .433 .300;
                     430 .377 .116;
                     440 .327 .033;
                     450 .295 .005;
                     460 .267 0;
                     470 .233 0;
                     480 .207 0;
                     490 .187 0;
                     500 .167 0;
                     510 .147 0;
                     520 .133 0;
                     530 .120 0;
                     540 .107 0;
                     550 .093 0;
                     560 .080 0;
                     570 .067 0;
                     580 .053 0;
                     590 .040 0;
                     600 .033 0;
                     610 .027 0;
                     620 .020 0;
                     630 .013 0;
                     640 .007 0;
                     650 0 0];

Tl = den_lens_pokorny(:,2)*(1+.02*(age-32))+den_lens_pokorny(:,3);
plot(den_lens_pokorny(:,1),Tl,'r*');
Tl = interp1(den_lens_pokorny(:,1),Tl,[380:5:780],'linear','extrap')';
Tl(Tl < 0) = 0;
lenstransmittance = 1./(10.^Tl)
synthfunds5 = absorptance.*repmat(lenstransmittance,1,3);
synthfunds5 = synthfunds5./repmat(max(synthfunds5),81,1);
T_cones_synthgh5= synthfunds5;
save T_cones_synthgh5.mat T_cones_synthgh5;

%% 
%Section 7.25
% Making synthetic cone fundamentals starting with the Baylor data.

load den_lens_ss; 
load den_mac_ss;
macpig = den_mac_ss;
lens = den_lens_ss;

a = [-5.2734 -87.403 1228.4 -3346.3 -5070.3 30881 -31607];
lmm = [561 531 430];
S = [];
for j = 1:length(lmm)
    tmp = log10((1./([380:5:780]/1000))*(lmm(j)/561));
    for i = 1:length(a)
        S(j,:,i) = a(i)*(tmp.^(i-1));
    end
end
actionspectra_baylor = sum(S,3);

% need to extrapolate the S-cone fundamental because the polynomial fit
% goes up again outside of the range of the data.
Lin = [380:5:780] >= 535 & [380:5:780] <= 575;
Lout = [380:5:780] > 575;
actionspectra_baylor(3,Lout) = interp1([535:5:575],actionspectra_baylor(3,Lin),[580:5:780],'linear','extrap')
actionspectra_baylor = 10.^actionspectra_baylor'; % getting rid of log.

macpigtransmittance = 1./(10.^(macpig./max(macpig)*.1));
lenstransmittance = 1./(10.^(lens*1));  % default
lenstransmittance = 1./(10.^(lens./lens(5)*1)); % density of 1 at 400 nm
opticaldensity = 0.3;
funds = 1-10.^(-actionspectra_baylor.*opticaldensity);
funds = funds.*repmat(macpigtransmittance,1,3).*repmat(lenstransmittance,1,3);
T_cones_synthgh6 = funds./repmat(max(funds),81,1);
save T_cones_synthgh6.mat T_cones_synthgh6;


%% 
% Outdated stuff for predicting the blue to green threshold ratio
% based on the SMJ-derived action spectra with various lens and macular
% pigment densities.

% % Adding optical density and lens density
% opticaldensity = .3;
% macpigtransmittance = 1./(10.^(macpig./max(macpig)*.4));
% lenstransmittance = 1./(10.^(lens*1.28));
% %lenstransmittance = 1./(10.^(lens/lens(5)));  % lens = 1 at 400 nm as per Norren
% plot(actionspectra_baylor.*opticaldensity)
% funds = 1-10.^(-(10.^actionspectra_baylor').*opticaldensity);
% funds = funds.*repmat(macpigtransmittance,1,3).*repmat(lenstransmittance,1,3);
% funds_baylor = funds./repmat(max(funds),81,1);
% figure; axes; hold on;
% tmph = plot([380:5:780],funds_baylor);
% h = tmph(1);
% % Comparing to the SMJ 10 deg fundamentals
% tmph = plot([380:5:780],T_cones_smj10',':');
% h(2) = tmph(1);
% legend(h,'Baylor','SMJ10');
% 
% % something around .49 would be good for Kali. .45 for Sedna
% % Looks like optical density = 0.3 and lens * 1.28 using Baylor
% % action spectra works well to predict monkey's data (mac pig goes from 0
% % to .4 as ecc. goes from 7 to 0.5)
% % According to Boetner, lens at 400 nm of 1.7 to 1 is OK.
% % But his transmitance curve looks very different from what I'm using
% (this should be action spectra?)
% load('T_log10coneabsorbance_ss.mat');
% plot([390:1:830],10.^T_log10coneabsorbance_ss','--')
% set(gca,'YScale','log','XLim',[380 780]);
% set(gca,'Ylim',[10^-9 1]);
% Making some cone fundmentals from the S&S absorption spectra

% Creating a set of cone fundamentals from the Baylor suction recordings
% What B/G ratio does this predict?

%load('succones.csv');
%figure; axes; hold on;
%plot(succones(:,1),10.^succones(:,[2:end]),'o');
%plot(succones(:,1),succones(:,[2:end]),'o');
% The action spectra from the SMJ fundamentals are too narrow. To make them
% sufficiently wide, you have to assume an unrealistically low optical
% pigment density when converting the SMJ fundamentals to action spectra (like 0.1).


%%
% Section 7.3
% Do these new (synth) cone fundamentals eliminate the S-cone tilt of the L+M
% planes measured with NeuroThresh? No. Not unless we add lots of MP.

% % Setting up the new cone fundamentals
% % load('succones.csv');
% load ('den_mac_ws'); % unused
% macpig = den_mac_ws; % unused
% a = [-5.2734 -87.403 1228.4 -3346.3 -5070.3 30881 -31607];
% lmm = [561 531 430];
% S = [];
% for j = 1:length(lmm)
%    tmp = log10((1./([380:5:780]/1000))*(lmm(j)/561));
%    for i = 1:length(a)
%        S(j,:,i) = a(i)*(tmp.^(i-1));
%    end
% end
% actionspectra_baylor = sum(S,3);
% % need to extrapolate the S-cone fundamental because the polynomial fit
% % goes up again outside of the range of the data.
% Lin = [380:5:780] >= 535 & [380:5:780] <= 575;
% Lout = [380:5:780] > 575;
% actionspectra_baylor(3,Lout) = interp1([535:5:575],actionspectra_baylor(3,Lin),[580:5:780],'linear','extrap');
% boet_wls = [320 340 360 380 400 420 440 460 480 500 550 600 650 700 750 800];
% boet_all = -log10([.8 .3 .2 .2 1 18 38.5 47 50.5 53.5 58.5 63.5 65.5 67.5 69.5 70.5]/100);
% boet_all = interp1(boet_wls,boet_all,[380:5:780],'linear')';
% %boet_all = boet_all-boet_all(end);
% opticaldensity = .3;
% macpigtransmittance = 1./(10.^(macpig./max(macpig)*0));  % .4 is the value from Snodderly
% alltransmittance = 1./(10.^(boet_all));
% funds = 1-10.^(-(10.^actionspectra_baylor').*opticaldensity);
% funds = funds.*repmat(macpigtransmittance,1,3).*repmat(alltransmittance,1,3);
% T_cones_synthgh1 = funds./repmat(max(funds),81,1);
% save T_cones_synthgh1.mat T_cones_synthgh1

% Above code is how we calculate T_cones_synthgh1
load('T_cones_synthgh1');
funds = T_cones_synthgh1;

% Gathering up the NT data
data = [];
[fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\NTlum.txt');
for cellcounter = 1:size(fnames,1)
    filename = findfile(char(fnames{cellcounter}));
    paradigmID = getparadigmID(filename);
    if paradigmID == 103
        NT = nex2stro(filename);
    end
    out = NTpreprocess(NT,.4,2);
    fundamentals = NT.sum.exptParams.fundamentals;
    fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
    mon_spd = NT.sum.exptParams.mon_spd;
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
    M = fundamentals'*mon_spd;
    Loog = logical(out(:,7));
    ecc = sqrt((NT.sum.exptParams.rf_x/10).^2 + (NT.sum.exptParams.rf_y/10).^2);
    data(cellcounter).M = M;
    data(cellcounter).data = out;
    data(cellcounter).bkgndrgb = NT.sum.exptParams.bkgndrgb;
    data(cellcounter).ecc = ecc;
    data(cellcounter).monspd = mon_spd;
end

for i = 1:length(data)
    cc_original = data(i).data(:,[2:4]).*repmat(data(i).data(:,5),1,3);
    Loog = data(i).data(:,end);
    M = funds'*mon_spd;
    scaled = ConvertConeContrastBasis(data(i).M, M, data(i).bkgndrgb, cc_original);
    [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);
    coneweights = (xformmat*planeparams)';
    if (sign(coneweights(1)) == -1)
        coneweights = -coneweights;
    end
    data(i).coneweights = coneweights;
end
allconeweights = cat(1,data(:).coneweights);
norms = sqrt(sum(allconeweights.*allconeweights,2));
normconeweights = allconeweights./repmat(norms,1,3);
[h,p] = ttest(normconeweights(:,3))

%%
% Section 8
% Trying to get from Baylor rod data to the human scotopic luminous
% sensitivity functions. Or vice versa. Something weird about about the
% lens density at very short wavelengths.
wls = [398 420 440 459 479 501 520 539 560 579 599 619 640 659 681 700 720];
s = [-.697 -.436 -.314 -.146 -.017 0 -.151 -.356 -.749 -1.22 -1.755 -2.312 -3.093 -3.743 -4.503 -5.147 -5.657];
load T_rods;
load den_lens_ws;
%lens = den_lens_smj(:,2);
lens = den_lens_ws;
opticaldensity = .35;
lensdensityat400 = 1.7;
lenstransmittance = 1./(10.^(lens*(lensdensityat400./lens(5))));
absorptance = T_rods'./lenstransmittance;
absorptance = absorptance./max(absorptance);
actionspectra = -log10(1-absorptance*(1-10^-opticaldensity));
actionspectra = actionspectra/opticaldensity;
figure; axes;
plot([380:5:780],log10(actionspectra));
hold on;
plot(wls,s,'ko');

% Predicting monkey rod fundamentals
rodactionspectra = interp1(wls,s,[380:5:780],'linear','extrap');
boet_wls = [320 340 360 380 400 420 440 460 480 500 550 600 650 700 750 800];
boet_all = -log10([.8 .3 .2 .2 1 18 38.5 47 50.5 53.5 58.5 63.5 65.5 67.5 69.5 70.5]/100);
boet_all = interp1(boet_wls,boet_all,[380:5:780],'spline')';
opticaldensity = .3;
alltransmittance = 1./(10.^(boet_all));
fund = 1-10.^(-(10.^rodactionspectra').*opticaldensity);
fund = fund.*alltransmittance;
fund = fund./repmat(max(fund),81,1);
figure; axes; hold on;
plot([380:5:780],fund,'b-');
plot([380:5:780],T_rods,'g-');

% Predicting monkey and human performance on DT_scot
load('Dell4BitsCal.mat')
cal = cals{end};
spds = SplineSpd(SToWls(cal.S_device), cal.P_device, SToWls([380 5 81]));

humanpreds = T_rods*spds;
monkeypreds = fund'*spds;
figure; axes; hold on;
plot(humanpreds,'k-'); plot(monkeypreds,'r-');

% Loading in a DT_scot datafile
stro = nex2stro(findfile('G031312001'));
all_data = stro.trial(:,[10:12 14]); % r, g, b, color_dir

for gun = 1:3
    
end

%%
% Section 9
% Finding the differences between the SMJ 10 deg funds, the SS 10 deg funds,
% and a set of synthetic fundamentals.
load ('den_mac_ws');
den_mac_ws(4) = 0.0425; %  minor tweak made by SMJ 1993: density of 0.0425 @ 395 nm
macpig = den_mac_ws./max(den_mac_ws);

load('T_cones_synthgh1');
load('T_cones_smj10');
load('T_cones_ss10');
T_cones_ss10 = interp1([390:1:830],T_cones_ss10',[380:5:780],'linear','extrap')';

% adding a tiny bit of macular pigment to the synth fundamentals
%macpigtransmittance = 1./(10.^(macpig*.05));
%T_cones_synthgh1 = T_cones_synthgh1.*repmat(macpigtransmittance,1,3);
%T_cones_synthgh1 = T_cones_synthgh1./repmat(max(T_cones_synthgh1),81,1);

% removing macular pigment from SMJ fundamentals
macpigtransmittance = 1./(10.^(macpig*.14));
T_cones_smj10 = T_cones_smj10'./repmat(macpigtransmittance,1,3);
T_cones_smj10 = T_cones_smj10./repmat(max(T_cones_smj10),81,1);

% removing macular pigment from SS fundamentals
macpigtransmittance = 1./(10.^(macpig*.095));
T_cones_ss10 = T_cones_ss10'./repmat(macpigtransmittance,1,3);
T_cones_ss10 = T_cones_ss10./repmat(max(T_cones_ss10),81,1);

figure; 
for conetype = 1:3
    subplot(3,1,conetype); hold on;
    plot([380:5:780],T_cones_synthgh1(:,conetype),'k-');
    plot([380:5:780],T_cones_smj10(:,conetype),'b-');
    plot([380:5:780],T_cones_ss10(:,conetype),'r-');
end
% SMJ and SS are very similar. Synth is quite different.

% Where are the differences in the SMJ and the synth funds?
figure; axes; hold on;
plot([380:5:780],T_cones_smj10-T_cones_synthgh1);
% T_cones_synthgh1 appears to have too much macular pigment. Weird.
% Uncomment the above to see it.

%%
% Section 10
% Comparing estimates of monkey photopic spectral sensitivity from the
% literature to sums of L and M cone fundamentals.
load T_cones_synthgh2;
load T_cones_synthgh1;
load T_cones_ss2;
load T_cones_ss10;
load T_cones_myss10;
load schrier_1966;
%load jacobs_1997;
load de_valois_1974;
load van_norren_1971;
h = []; g = [];

schrier1966(:,2) = 10.^schrier1966(:,2);
van_norren_1971(:,2) = 10.^(van_norren_1971(:,2)-max(van_norren_1971(:,2)));
de_valois_1974(:,2) = 10.^de_valois_1974(:,2);

figure; axes; hold on;
% 10 deg funds 2:1 (L:M)
pred_ss10_2 = 2*T_cones_ss10(:,1)+T_cones_ss10(:,2);
pred_ss10_2 = pred_ss10_2./max(pred_ss10_2);
h(1)=plot([380:5:780],pred_ss10_2,'r-','LineWidth',2);

% 10 deg funds 1:1 (L:M)
pred_ss10_1 = T_cones_ss10(:,1)+T_cones_ss10(:,2);
pred_ss10_1 = pred_ss10_1./max(pred_ss10_1);
h(2)=plot([380:5:780],pred_ss10_1,'k-','LineWidth',2);

% Lowered lens pigmentation
pred_synth2 = T_cones_synthgh2(:,1)+T_cones_synthgh2(:,2);
pred_synth2 = pred_synth2./max(pred_synth2);
h(3)=plot([380:5:780],pred_synth2,'b-','LineWidth',2);

% % Boettner lens pigmentation
% pred_synth1 = T_cones_synthgh1(:,1)+T_cones_synthgh1(:,2);
% pred_synth1 = pred_synth1./max(pred_synth1);
% plot([380:5:780],pred_synth1,'m-');

set(gca,'Ylim',[10^-3 2]);
set(gca,'YScale','log');

wls = [380:5:780];
wls_thresh = 540;
Lwhich_wls = logical(wls>=wls_thresh);
mymodel = @(scale,x) (scale*interp1(wls(Lwhich_wls),log10(pred_ss10_1(Lwhich_wls)),x,'linear'));

L = van_norren_1971(:,1) >= wls_thresh;
scale = nlinfit(van_norren_1971(L,1),log10(van_norren_1971(L,2)),mymodel,1);
g(1)=plot(van_norren_1971(:,1),scale*van_norren_1971(:,2),'kv','MarkerFaceColor','black');

L = schrier1966(:,1) >= wls_thresh;
scale = nlinfit(schrier1966(L,1),log10(schrier1966(L,2)),mymodel,1);
g(2)=plot(schrier1966(:,1),scale*schrier1966(:,2),'ko','MarkerFaceColor','black');

L = de_valois_1974(:,1) >= wls_thresh;
scale = nlinfit(de_valois_1974(L,1),log10(de_valois_1974(L,2)),mymodel,1);
g(3)=plot(de_valois_1974(:,1),scale*de_valois_1974(:,2),'k^','MarkerFaceColor','black');

xlabel('Wavelength (nm)','FontSize',14);
ylabel('Relative sensitivity','FontSize',14);
legend([h,g],{'2L:M','1L:M','1L:M (lens=1@400nm)','Van Norren 1971','Schrier 1966','DeValois 1974'},'location','south');

%%
% Section 11
% synthesizing a transmittance spectrum for the Kodak neutral density
% filters from the published data.
% "dig" = digitized (from the graph), "trans" = transcribed (from text)

load wrattennd1.mat
S = [380:5:780];
figure; axes; hold on;
plot(wratten_nd1_dig(:,1),wratten_nd1_dig(:,2),'r-')
plot(wratten_nd1_trans(:,1),wratten_nd1_trans(:,2),'b-');
y = nan*ones(size(S));
% First, where we have numerical data
L = S>=wratten_nd1_trans(1,1) & S<=wratten_nd1_trans(end,1);
y(L) = interp1(wratten_nd1_trans(:,1),wratten_nd1_trans(:,2), S(L));
% Now the short wavelengths
L = S<wratten_nd1_trans(1,1);
y(L) = interp1(wratten_nd1_dig(:,1),wratten_nd1_dig(:,2), S(L),'linear','extrap');
% Now the long wavelengths
L = S>wratten_nd1_trans(end,1);
y(L) = interp1(wratten_nd1_dig(:,1),wratten_nd1_dig(:,2), S(L),'linear','extrap');
wratten_nd1 = [S' y'/100];
save wrattennd1 wratten_nd1_dig wratten_nd1_trans wratten_nd1;

%%
% Section 12
% Playing around with the Pokorny et al. (1987) lens model
% The age only changes the shallow part, long wavelength part of the 
% absorption function - doesn't change the density at 400 nm

den_lens_pokorny = [400 .600 1.00;
                     410 .510 .583;
                     420 .433 .300;
                     430 .377 .116;
                     440 .327 .033;
                     450 .295 .005;
                     460 .267 0;
                     470 .233 0;
                     480 .207 0;
                     490 .187 0;
                     500 .167 0;
                     510 .147 0;
                     520 .133 0;
                     530 .120 0;
                     540 .107 0;
                     550 .093 0;
                     560 .080 0;
                     570 .067 0;
                     580 .053 0;
                     590 .040 0;
                     600 .033 0;
                     610 .027 0;
                     620 .020 0;
                     630 .013 0;
                     640 .007 0;
                     650 0 0];

figure; axes; hold on;
age = 5;
Tl = den_lens_pokorny(:,2)*(1+.02*(age-32))+den_lens_pokorny(:,3);
plot(den_lens_pokorny(:,1),Tl,'r*');
Tl = interp1(den_lens_pokorny(:,1),Tl,[380:5:780],'linear','extrap')';
Tl(Tl < 0) = 0;
plot([380:5:780],Tl,'r-');

load den_lens_ws;
plot([380:5:780],den_lens_ws,'r-')
plot([380:5:780],den_lens_ws./den_lens_ws(5),'r:')
%boet_wls = [320 340 360 380 400 420 440 460 480 500 550 600 650 700 750 800];
%boet_all = -log10([.8 .3 .2 .2 1 18 38.5 47 50.5 53.5 58.5 63.5 65.5 67.5 69.5 70.5]/100);
%boet_all = interp1(boet_wls,boet_all,[380:5:780],'linear')';


%% 
% Section 13
% Presenting and measuring the spectra of threshold-level stimuli
% Need to input RGBs by hand
nWavelengthSamples = 101;
S = [380 4 nWavelengthSamples];  % Native PR650 settings
boxSize = 200;
Xoffset = 0;
Yoffset = 0;
bkgndrgb = [0 0 0];

disp('Turn on Photometer and hit <return>');
disp('Then focus on white square and hit <return> again');
input('');
CMCheckInit();

[window,screenRect] = Screen(0,'OpenWindow',bkgndrgb);
Screen('LoadNormalizedGammaTable',window,linspace(0,1,256)'*ones(1,3));

boxRect = [0 0 boxSize boxSize];
boxRect = CenterRect(boxRect,screenRect);
boxRect = boxRect + [Xoffset -1*Yoffset Xoffset -1*Yoffset];
Screen('FillRect',window,[25 25 25],boxRect);
Screen('Flip',window);

input('');
pause(2);  % Time to leave room

spdMat = zeros(size(RGBs,1), nWavelengthSamples);
for i = 1:size(RGBs,1)
    Screen('FillRect',window,255*RGBs(i,:),boxRect);
    Screen('Flip',window);
    spd = PR650measspd(S);
    spdMat(i,:) = spd';
end
sca();
PR650close;

save junkout spd RGBs

% Now seeing how well these match the predictions from human scotopic Vlambda
load('T_rods');
%T_rods = SplineSpd([380:5:780]',T_rods',[380:4:780]');
load wrattennd1
filter = wratten_nd1(:,2);
spds = SplineSpd([380:4:780]',spdMat,[380:5:780]');

