%%
wls = [398 420 440 459 479 501 520 539 560 579 599 619 640 659 681 700 720];
s = [-.697 -.436 -.314 -.146 -.017 0 -.151 -.356 -.749 -1.22 -1.755 -2.312 -3.093 -3.743 -4.503 -5.147 -5.657];
load T_rods;
load den_lens_ss;  % best choice?
lens = den_lens_ss;
opticaldensity = .4;
lensdensityat400 = lens(5);
lenstransmittance = 1./(10.^(lens*(lensdensityat400./lens(5))));
absorptance = T_rods'./lenstransmittance;
absorptance = absorptance./max(absorptance);
%actionspectra = -log10(1-absorptance*(1-10^-opticaldensity));
%actionspectra = actionspectra/opticaldensity;

% Predicting monkey rod fundamentals
rodactionspectra = interp1(wls,s,[380:5:780],'spline','extrap')';
% boet_wls = [320 340 360 380 400 420 440 460 480 500 550 600 650 700 750 800];
% boet_all = -log10([.8 .3 .2 .2 1 18 38.5 47 50.5 53.5 58.5 63.5 65.5 67.5 69.5 70.5]/100);
% boet_all = interp1(boet_wls,boet_all,[380:5:780],'spline')';
opticaldensity = .35;
% alltransmittance = 1./(10.^(boet_all));

%lensdensityat400 = 1;
%lenstransmittance = 1./(10.^(lens*(lensdensityat400./lens(5))));

fund = 1-10.^(-(10.^rodactionspectra).*opticaldensity);
%fund = 1-10.^(-actionspectra.*opticaldensity);
%fund = fund.*lenstransmittance;
fund = fund./repmat(max(fund),81,1);
absorptance = fund';
%figure; axes; hold on; plot(fund,'k:'); plot(T_rods,'k-'); legend('Monkey','human');
%%
% Generating prediction

lensdensityat400 = lens(5);
lenstransmittance = 1./(10.^(lens*(lensdensityat400./lens(5))));
absorptance = T_rods'./lenstransmittance;
absorptance = absorptance./max(absorptance);

lensdensityat400 = 1;
lenstransmittance = 1./(10.^(lens*(lensdensityat400./lens(5))));
fund = absorptance.*lenstransmittance;
fund = fund./repmat(max(fund),81,1);
figure; axes; hold on; plot(fund,'k:'); plot(T_rods,'k-'); legend('Monkey','human');

%%
load('filtertests.mat');
load('T_rods.mat');
transmittance = newtest3{2,2} ./ newtest3{1,2};
transmittance = SplineRaw([380 2 201], transmittance, [380 5 81]);

% transmittance = ones(81,1);

stro = nex2stro;
colordir_idx = strcmp(stro.sum.trialFields(1,:), 'color_dir');
questmode_idx = strcmp(stro.sum.trialFields(1,:), 'quest_mode');
thresholds = zeros(3,1);

mon_spd = reshape(stro.sum.exptParams.mon_spd, 101, 3);
mon_spd = SplineSpd([380 4 101], mon_spd, [380 5 81]);

figure;
for gun = 1:3
    L = stro.trial(:, colordir_idx) == gun;
    all_modes = stro.trial(L, questmode_idx);
    subplot(3,1,gun);
    plot(all_modes,'b.-');
    thresholds(gun) = all_modes(end);
end

new_mon_spd = mon_spd .* repmat(transmittance .^ 6, 1, 3);
sensitivity = 1 ./ thresholds;
theory = T_rods * new_mon_spd;
monkey_theory = fund' * new_mon_spd;
monkey_scale = monkey_theory' \ sensitivity;
scalefactor = theory' \ sensitivity;

figure; subplot(2,1,1); plot(1:3, theory*scalefactor, 'k-');
hold on; plot(1:3, sensitivity, 'r.');
set(gca, 'xlim', [.5 3.5], 'xtick', 1:3, 'xticklabel', {'R' 'G' 'B'});
title('Human solid - monkey dashed');
plot(1:3, monkey_theory*monkey_scale, 'k:');
set(gcf,'Name', stro.sum.fileName);

%%
% A bit of population analysis for an individual observer
%[fnames, spikenums] = fnamesFromTxt2('/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/DTEM/SednaDTscot.txt');

if (ispc)
    [fnames, spikenums] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\DTEM\GregDTscot.txt');
else
    [fnames, spikenums] = fnamesFromTxt2('/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/DTEM/GregDTscot.txt');    
end
data = [];
for i = 1:length(fnames)
    stro = nex2stro(findfile(fnames{i}));
    
    r_idx = strcmp(stro.sum.trialFields(1,:), 'stim_r');
    g_idx = strcmp(stro.sum.trialFields(1,:), 'stim_g');
    b_idx = strcmp(stro.sum.trialFields(1,:), 'stim_b');
    rgbs = stro.trial(:,r_idx|g_idx|b_idx);
    questmode_idx = strcmp(stro.sum.trialFields(1,:), 'quest_mode');
    
    thresholds = nan*ones(3,1);
    whichgunmat = reshape(stro.sum.exptParams.stim_colors,3,3);
    
    mon_spd = reshape(stro.sum.exptParams.mon_spd, length(stro.sum.exptParams.mon_spd)/3, 3);
    if (size(mon_spd,1) == 101)
        mon_spd = SplineSpd([380 4 101], mon_spd, [380 5 81]);
    else
        mon_spd = SplineSpd([380 2 201], mon_spd, [380 5 81]);
    end
    
    figure;
    for gun = 1:3
        L = logical(rgbs(:,gun));
        if (any(L))
            all_modes = stro.trial(L, questmode_idx);
            subplot(3,1,gun);
            plot(all_modes,'b.-');
            thresholds(gun) = all_modes(end);
        end
    end
    data = [data; thresholds'];
end
figure; axes; hold on;
plot(log10(data(:,1)),'r.');
plot(log10(data(:,2)),'g.');
plot(log10(data(:,3)),'b.');
ylabel('Threshold'); xlabel('Session');

%%
% Doing population analysis again
% Making it optional to spline the spds or the T_rods
% This demonstrates that the difference between splining the spds
% is equivalent to splining T_rods.

SPLINESPDS = 1;
load wrattennd1
if (ispc)
    [fnames, spikenums] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\DTEM\GregDTscot.txt');
else
    [fnames, spikenums] = fnamesFromTxt2('/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/DTEM/GregDTscot.txt');    
end
data = []; dotprods = [];
for i = 1:length(fnames)
    stro = nex2stro(findfile(fnames{i}));
    colordir_idx = strcmp(stro.sum.trialFields(1,:), 'color_dir');
    questmode_idx = strcmp(stro.sum.trialFields(1,:), 'quest_mode');
    thresholds = zeros(3,1);
    
    mon_spd = reshape(stro.sum.exptParams.mon_spd, 101, 3);
    if (SPLINESPDS)
        mon_spd = SplineSpd([380 4 101], mon_spd, [380 5 81]);
        wratten = wratten_nd1(:,2);
        rods = T_rods;
    else % spline T_rods
        wratten = SplineSpd([380 5 81], wratten_nd1(:,2), [380 4 101]);
        rods = SplineSpd([380 5 81], T_rods', [380 4 101])';
    end
    
    figure;
    for gun = 1:3
        L = stro.trial(:, colordir_idx) == gun;
        all_modes = stro.trial(L, questmode_idx);
        subplot(3,1,gun);
        plot(all_modes,'b.-');
        thresholds(gun) = all_modes(end);
    end
    data = [data; thresholds'];
    threshspd = repmat(thresholds',size(mon_spd,1),1).*repmat(wratten,1,size(mon_spd,2)).*mon_spd;
    dotprods = [dotprods; rods*threshspd];
end
figure; axes; hold on;
plot(log10(data(:,1)),'r.');
plot(log10(data(:,2)),'g.');
plot(log10(data(:,3)),'b.');
ylabel('Threshold'); xlabel('Session');

% Large dotproduct means low sensitivity
figure; axes; hold on;
plot(log10(dotprods(:,1)),'r.');
plot(log10(dotprods(:,2)),'g.');
plot(log10(dotprods(:,3)),'b.');
ylabel('log dotprod'); xlabel('Session');


%%
% Trying a few lens densities
% Weirdly low threshold to red. What gives?
% Sliding theoretical predictions up to match mean log sensitivities
% Red is mean, green is vlambda prediction.
load ('den_mac_ws');
den_mac_ws(4) = 0.0425;
macpig = den_mac_ws;
macpigdensityat360 = 0;

% Making new monkey predictions
lensdensityat400 = 1;
lenstransmittance = 1./(10.^(lens*(lensdensityat400./lens(5))));
opticaldensity = .4 % .4 is the estimate from Baylor et al.
absorptance = 1-10.^(-(10.^rodactionspectra).*opticaldensity);
absorptance = absorptance./max(absorptance);

fund = absorptance.*lenstransmittance;
fund = fund./repmat(max(fund),81,1);

figure; subplot(2,1,1); hold on;
set(gca,'Color',[.5 .5 .5]);
plot(1./data','k.');
set(gca,'Yscale','log');
% 
% lensdensities = linspace(0,3,20); plotcols = jet(length(lensdensities));
% SSE = [];
% for i = 1:length(lensdensities)
%     lensdensityat400 = lensdensities(i);
%     lenstransmittance = 1./(10.^(lens*(lensdensityat400./lens(5))));
%     mptransmittance = 1./(10.^(macpig*(macpigdensityat360./max(macpig))));
%     fund = absorptance.*lenstransmittance;
%     fund = fund.*mptransmittance;
%     fund = fund./repmat(max(fund),81,1);
%     theory = fund' * new_mon_spd;    
%     corrections = mean(log10(1./data))-mean(log10(theory));
%     scalefactor = 10.^mean(corrections);
%     h = plot(theory.*scalefactor);
%     set(h,'color',plotcols(i,:));
%     err = repmat(log10(theory*scalefactor),size(data,1),1)-log10(1./data);
%     SSE(i,:) = sum(err.^2);
% end
% prediction based on vlambda
theory = T_rods * new_mon_spd;
corrections = mean(log10(1./data))-mean(log10(theory));
scalefactor = 10.^mean(corrections);
plot(theory.*scalefactor,'g:','LineWidth',2);
theory = fund' * new_mon_spd;
corrections = mean(log10(1./data))-mean(log10(theory));
scalefactor = 10.^mean(corrections);
plot(theory.*scalefactor,'b:','LineWidth',2);

% Plotting the geometric means
plot(geomean(1./data),'mo','MarkerSize',6,'MarkerFaceColor','magenta');  

% subplot(2,1,2); hold on;
% plot(lensdensities, sum(SSE,2),'k*-');
% ylabel('Error'); xlabel('lens density');
% err = repmat(log10(theory*scalefactor),size(data,1),1)-log10(1./data);
% plot(lens(5),sum(err(:).^2),'g*');

% Looking at residuals from the Vlambda prediction
theory = T_rods * new_mon_spd;
logsens = log10(1./data);
meancorrection = mean(logsens(:))-mean(log10(theory));
pred = log10(theory) + meancorrection;
resid = logsens-repmat(pred,size(data,1),1);
set(gca,'Yscale','log','Xlim',[.5 3.5]);

subplot(2,1,2); set(gca,'Color',[.5 .5 .5]); hold on; 
plot(10.^resid','k.');
plot(10.^mean(resid),'mo','MarkerSize',6,'MarkerFaceColor','magenta');  
plot([1 2 3],[1 1 1],'b:','LineWidth',2);

%Looking at residuals from a synthetic "monkey" prediction
theory = fund' * new_mon_spd;
logsens = log10(1./data);
meancorrection = mean(logsens(:))-mean(log10(theory));
monkeypred = log10(theory) + meancorrection;
plot([1 2 3],10.^(monkeypred-pred),'g:','linewidth',2)
set(gca,'Yscale','log','Xlim',[.5 3.5]);

% Calculating some F-statistics
% I basically want to do a 1-way ANOVA on the residuals
% Is there significant structure in the residuals once we
% project onto the prediction?
% 
% num = sum((size(resid,1)*mean(resid)).^2);
% num = num./2;
% den = sum(sum((resid - repmat(mean(resid),size(resid,1),1)).^2));
% %den = sum(sum((resid.^2)));
% den = den./(numel(resid)-3);
% F = num/den
% 1-fcdf(F,1,numel(resid))
resid = logsens-repmat(pred,size(data,1),1);
anova1(resid)

% Now monkey pred
resid = logsens-repmat(monkeypred,size(data,1),1);
anova1(resid)

% After projecting onto the vlambda prediction (which soaks up a lot of
% variance) is there a significant component of the Y vector in the
% direction of the component of monkeypred that is orthogonal to pred.
projontopred = (pred*monkeypred')/sum(pred.^2)*pred;
monkeypredorth = monkeypred-projontopred; % This is the component of monkeypred that 
% is orthogonal to pred.
% The question is: is there a significant projection of the residual vector
% onto monkeypredorth?
resid = logsens-repmat(pred,size(data,1),1);
resid = resid(:);


% Now projecting the residual Y vector onto the monkeypredorth vector
X = reshape(repmat(monkeypredorth,size(data,1),1),numel(data),1);
proj = (resid'*(X./norm(X)))*resid;
num = sum(proj.^2);
% getting the squared length of the error vector
denom = (sum(resid.^2)-num)./(length(resid)-2);
F = num./denom
p = 1-fcdf(F,1,length(resid)-2) % one df for mean, one for pred

% Partial correlation
% -------------------
x = logsens(:);
y = reshape(repmat(monkeypred,size(data,1),1),numel(data),1);
z = reshape(repmat(pred,size(data,1),1),numel(data),1);
%Using the formula
r = corrcoef([x,y,z])
partialcor = (r(1,2)-r(1,3)*r(2,3))/(sqrt(1-r(1,3).^2)*sqrt(1-r(2,3).^2));
disp(['Partial correlation (orthogonal to vlambda prediction): ',num2str(partialcor)])

% From first principles
xprime = x-((x'*z)/sum(z.^2))*z;
yprime = y-((y'*z)/sum(z.^2))*z;
xprime'*z
yprime'*z % round off error?
corrcoef([xprime,yprime])
(xprime./norm(xprime))'*(yprime./norm(yprime))

%%
% Looking at each set of measurements in 3-D space. Human and monkey
% predictions as lines? Loading data from all subjects.

if(ispc)
    filelistpath = 'N:\NexFiles\nexfilelists\Greg\DTEM';
else
   filelistpath = '/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/DTEM';
end

listnames = {'GregDTscot.txt','LeahDTscot.txt','FreyaDTscot.txt','SednaDTscot.txt' };
subjectnames = {'Human G','Human L','Monkey F','Monkey S'};
%listnames = {'GregDTscot.txt'};
%subjectnames = {'Human G'};

load('T_rods');
data = [];
for i = 1:length(listnames)
    [fnames, spikenums] = fnamesFromTxt2([filelistpath,filesep,listnames{i}]);
    for j = 1:length(fnames)
        stro = nex2stro(findfile(char(fnames{j})));
        
        colordir_idx = strcmp(stro.sum.trialFields(1,:), 'color_dir');
        questmode_idx = strcmp(stro.sum.trialFields(1,:), 'quest_mode');
        thresholds = zeros(3,1);
        
        for gun = 1:3
            L = stro.trial(:, colordir_idx) == gun;
            all_modes = stro.trial(L, questmode_idx);
            thresholds(gun) = all_modes(end);
        end
        data = [data; thresholds' i];
    end
end
% % A bit of data culling...
% L = [];
% sdthreshold = 2;
% for i = 1:length(listnames)
%     Lobs = data(:,end) == i;
%     logthresh = log10(data(Lobs,[1 2 3]));
%  %   L = [L; logthresh<repmat(mean(logthresh)-sdthreshold*std(logthresh),sum(Lobs),1) |...
%  %          logthresh>repmat(mean(logthresh)+sdthreshold*std(logthresh),sum(Lobs),1)];
%     L = [L; logthresh < repmat(prctile(logthresh,5),sum(Lobs),1) | ...
%             logthresh > repmat(prctile(logthresh,95),sum(Lobs),1)];
% end
% trimmeddata = data(:,[1 2 3]);
% trimmeddata(logical(L)) = nan;
% trimmeddata(:,end+1) = data(:,end);
trimmeddata = data;

% Reds are humans, greens are monkeys
Lhuman = strncmp('Human',subjectnames,5);
cols = zeros(length(listnames),3);
cols(Lhuman,1) = linspace(1,.5,sum(Lhuman));
cols(~Lhuman,2) = linspace(1,.5,sum(~Lhuman)); 

figure; axes; hold on;
xlabel('Red');ylabel('Green');zlabel('Blue');
for i = 1:length(listnames)
    L = trimmeddata(:,end) ==i;
    h = plot3(log10(trimmeddata(L,1)),log10(trimmeddata(L,2)),log10(trimmeddata(L,3)),'o')
    set(h,'MarkerFaceColor',cols(i,:),'MarkerEdgeColor','none');
    mn = nanmean(log10(trimmeddata(L,[1 2 3])));
    S2 = nancov(log10(trimmeddata(L,[1 2 3])));
    [v,d] = eig(S2);
    h = plot3(mn(1),mn(2),mn(3),'h');
    set(h,'MarkerFaceColor',cols(i,:),'MarkerEdgeColor','black','MarkerSize',14)
    % Confidence ellipsoid
    sem = sqrt(diag(d)./sum(L));
    [x, y, z] = ellipsoid(0,0,0,1.96*sem(1),1.96*sem(2),1.96*sem(3));
    newxyz = [x(:) y(:) z(:)]*v+repmat(mn,length(x(:)),1);
    h = surf(reshape(newxyz(:,1),size(x)), reshape(newxyz(:,2),size(y)), reshape(newxyz(:,3),size(z)));
    set(h,'FaceAlpha',.2,'FaceColor',cols(i,:),'Edgealpha',0);
end

% Getting the monitor spectra
load('Dell4blackbkgnd.mat')
cal = cals{end};
mon_spd = SplineSpd(cal.S_device, cal.P_device, [380:5:780]');
%load('filtertests.mat');
%transmittance = newtest3{2,2} ./ newtest3{1,2};
%transmittance = SplineRaw([380 2 201], transmittance, [380 5 81]);
load('wrattennd1.mat');
transmittance = wratten_nd1(:,2);
new_mon_spd = mon_spd .* repmat(transmittance .^ 6, 1, 3);

% Plotting a line for the human prediction
humanpred =  log10(1./(T_rods*new_mon_spd));
roughcenter = humanpred-mean(humanpred)+mean(mn);
plot3(.5*[-1 1]+roughcenter(1),.5*[-1 1]+roughcenter(2),.5*[-1 1]+roughcenter(3),'r-','Linewidth',2);

% Plotting a line for the monkey prediction
wls = [398 420 440 459 479 501 520 539 560 579 599 619 640 659 681 700 720];
s = [-.697 -.436 -.314 -.146 -.017 0 -.151 -.356 -.749 -1.22 -1.755 -2.312 -3.093 -3.743 -4.503 -5.147 -5.657];
opticaldensity = .3;
rodactionspectra = interp1(wls,s,[380:5:780],'linear','extrap');
%boet_wls = [320 340 360 380 400 420 440 460 480 500 550 600 650 700 750 800];
%boet_all = -log10([.8 .3 .2 .2 1 18 38.5 47 50.5 53.5 58.5 63.5 65.5 67.5 69.5 70.5]/100);
%boet_all = interp1(boet_wls,boet_all,[380:5:780],'spline')';opticaldensity = .35;
%alltransmittance = 1./(10.^(boet_all));
load den_lens_ws;
lens = den_lens_ws;
lensdensityat400 = 1;
lenstransmittance = 1./(10.^(lens*(lensdensityat400./lens(5))));
fund = 1-10.^(-(10.^rodactionspectra').*opticaldensity);
fund = fund.*lenstransmittance;  % option to use lenstransmittance or alltransmittance here
fund = fund'./repmat(max(fund),1,81); % transposing here
monkeypred = log10(1./(fund*new_mon_spd));
roughcenter = monkeypred-mean(monkeypred)+mean(mn);
plot3(.5*[-1 1]+roughcenter(1),.5*[-1 1]+roughcenter(2),.5*[-1 1]+roughcenter(3),'g-','Linewidth',2);

% Plotting thresholds over time
for i = 1:max(data(:,end))
    figure; axes; hold on;
    L = data(:,end) == i;
    plot(trimmeddata(L,1),'r.');
    plot(trimmeddata(L,2),'g.');
    plot(trimmeddata(L,3),'b.');
    set(gca,'Xlim',[1 sum(L)],'Yscale','log');
    title(subjectnames{i});
end
    

% Projecting everything orthogonal to the human prediction.
% I think we need to take logs *after* projection, not before.
humanpred =  1./(T_rods*new_mon_spd); humanpred = humanpred./norm(humanpred);
monkeypred =  1./(fund*new_mon_spd); monkeypred = monkeypred./norm(monkeypred); 
projmags = data(:,[1 2 3])*humanpred';
projs = repmat(humanpred,size(projmags,1),1) .*repmat(projmags,1,3);
orthprojs = data(:,[1 2 3]) - projs; % projections orthogonal to the human prediction
monkorthhuman = monkeypred - humanpred*(monkeypred*humanpred');
monkorthhuman = monkorthhuman./norm(monkorthhuman);
% p = orthprojs*monkorthhuman'; % Projection onto the component of the monkey prediction vector that is orthogonal to the human
monkprojmags = orthprojs*monkorthhuman'; % projection onto monkeypred, orthogonal to humanpred 
residual = orthprojs - repmat(monkorthhuman,size(monkprojmags,1),1) .*repmat(monkprojmags,1,3);
L = data(:,end) < 4;
plot(abs(monkprojmags(L)),sqrt(sum(residual(L).^2,2)),'r.');
plot(abs(monkprojmags(~L)),sqrt(sum(residual(~L).^2,2)),'g.');

%%
% An analysis based just on green and blue ratio
bgthreshratio = data(:,2)./data(:,3);
figure; 
subplot(2,1,1); hold on;
hist(bgthreshratio,20);
plot(geomean(bgthreshratio),1.5,'kv','MarkerSize',9,'MarkerFaceColor','black')
% prediction based on vlambda
theory = T_rods * new_mon_spd;
plot(theory(3)./theory(2),0,'r*')
xlabel('b/g threshold ratio');
lensdensities = linspace(0,2,10);
plotcols = hot(length(lensdensities));

% prediction based on synthetics
for i = 1:length(lensdensities)
    lensdensityat400 = lensdensities(i);
    lenstransmittance = 1./(10.^(lens*(lensdensityat400./lens(5))));
    fund = absorptance.*lenstransmittance;
  %  fund = fund.*mptransmittance;
    fund = fund./repmat(max(fund),81,1);
    theory = fund' * new_mon_spd;    
    h = plot(theory(3)/theory(2),.1,'mo');
    set(h,'markerFaceColor',plotcols(i,:));
end

subplot(2,1,2); hold on;
hist(log10(bgthreshratio),20);
plot(mean(log10(bgthreshratio)),1.5,'kv','MarkerSize',9,'MarkerFaceColor','black')
% prediction based on vlambda
theory = T_rods * new_mon_spd;
plot(log10(theory(3)./theory(2)),0,'r*')
xlabel('log b/g threshold ratio');

% prediction based on synthetics
for i = 1:length(lensdensities)
    lensdensityat400 = lensdensities(i);
    lenstransmittance = 1./(10.^(lens*(lensdensityat400./lens(5))));
    fund = absorptance.*lenstransmittance;
 %   fund = fund.*mptransmittance;
    fund = fund./repmat(max(fund),81,1);
    theory = fund' * new_mon_spd;    
    h = plot(log10(theory(3)/theory(2)),.1,'mo');
    set(h,'markerFaceColor',plotcols(i,:));
end


%% Trying to calculate cd/m^2, trolands, etc from the calibration
% information in the stro file.

stro = nex2stro(findfile('G052912008','/Volumes/NO BACKUP/NexFiles/Greg/Psychophysics'));
tmp_spd = reshape(stro.sum.exptParams.mon_spd,101,3);
load('T_xyz1931')
vlambda = T_xyz1931(2,:);
[mon_spd] = SplineSpd([380:4:780]', tmp_spd, [380:5:780]')
lum = 680*vlambda*mon_spd;
% Luminance at the midpoint of the monitor gamut
lum*[.5 .5 .5]'  % that looks about right

% Now getting scotopic trolands for each gun

load('T_rods')
load('filtertests.mat');
transmittance = newtest3{2,2} ./ newtest3{1,2};
transmittance = SplineRaw([380 2 201], transmittance, [380 5 81]);
new_mon_spd = mon_spd .* repmat(transmittance .^ 6, 1, 3);
scotlum = 1700*T_rods*new_mon_spd;

% Now at threshold
colordir_idx = strcmp(stro.sum.trialFields(1,:), 'color_dir');
questmode_idx = strcmp(stro.sum.trialFields(1,:), 'quest_mode');
thresholds = zeros(3,1);
for gun = 1:3
    L = stro.trial(:, colordir_idx) == gun;
    all_modes = stro.trial(L, questmode_idx);
    thresholds(gun) = all_modes(end);
end
scotlumatthresh = scotlum.*thresholds';
pupildiam = 5% % an esitmate
pupilarea = pi*(pupildiam/2)^2; % mm^2
pupilarea = 50; % from W&S
scottrolands = scotlumatthresh*pupilarea


load ('Dell4BitsCal');
cal = cals{end};
plot(cal.P_ambient)
ambient = SplineSpd([380:4:780]', cal.P_ambient, [380:5:780]') .* (transmittance .^ 6);
bkgndscottrol = 1700*T_rods*ambient*pupilarea  % scotopic trolands of background
contrast = scottrolands./bkgndscottrol

% Can I get to quanta/sec?
energy = EnergyToQuanta([380:5:780]',new_mon_spd);
% from a distance of 100 cm (1000 mm), this light is spread over 1000^2 mm
% Lets say a pupil is 50 mm^2
anglefraction = (50^2./1000^2);
energy*anglefraction
% stimulus is 6*6 = 36 cm ^2 = .06*.06 = 0.0036 m^2
stimsize = 0.0036;
sum(energy*stimsize*anglefraction)


