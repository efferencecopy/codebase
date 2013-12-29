% A bit of population analysis. Loading DTspot_scot data files, calculating
% mean thresholds, and measuring them with the PR650. Comparint the results
% to the predictions from the calibration structure.

whichphotometer = 'PR705';
NFILTERS = 1; % number of ND filters used in measurement
nreps = 10;

load('Dell4blackbkgnd.mat');
cal = cals{end};
gammaTable = cal.gammaTable;
invgamma = InvertGamma(gammaTable, 1); %  1 = colour mode
NGAMMASTEPS = size(invgamma,1);

[fnames, spikenums] = fnamesFromTxt2('/Volumes/128.95.153.12/NexFiles/nexfilelists/Greg/DTEM/GregDTscot.txt');
%[fnames, spikenums] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\DTEM\GregDTscot.txt');
data = [];
for i = 1:length(fnames)
    stro = nex2stro(findfile(fnames{i},'/Volumes/128.95.153.12/NexFiles/Greg/Psychophysics'));
    
    colordir_idx = strcmp(stro.sum.trialFields(1,:), 'color_dir');
    questmode_idx = strcmp(stro.sum.trialFields(1,:), 'quest_mode');
    thresholds = zeros(3,1);
    
    mon_spd = reshape(stro.sum.exptParams.mon_spd, 101, 3);
    mon_spd = SplineSpd([380 4 101], mon_spd, [380 5 81]);
    
%     figure;
    for gun = 1:3
        L = stro.trial(:, colordir_idx) == gun;
        all_modes = stro.trial(L, questmode_idx);
%         subplot(3,1,gun);
%         plot(all_modes,'b.-');
        thresholds(gun) = all_modes(end);
    end
    data = [data; thresholds'];
end
means = geomean(data);
rgbs = round(means*(NGAMMASTEPS-1))+1;
RGBs = diag(round((NGAMMASTEPS-1) * ...
            [invgamma(rgbs(1),1) ...
            invgamma(rgbs(2),2) ...
            invgamma(rgbs(3),3)]));
        
        
RGBidxs = repmat([1:size(RGBs,1)],1,nreps);
orderidxs = RGBidxs(randperm(length(RGBidxs)));

% Now presenting and measuring stimuli
% Setting stuff up
if (strcmp(whichphotometer, 'PR705'))
    nWavelengthSamples = 201;
    S = [380 2 nWavelengthSamples];  % Native PR705 settings
else
    nWavelengthSamples = 101;
    S = [380 4 nWavelengthSamples];  % Native PR650 settings
end
boxSize = 200;
Xoffset = 0;
Yoffset = 0;
bkgndrgb = [0 0 0];

disp('Turn on Photometer and hit <return>');
disp('Then focus on white square and hit <return> again');
input('');
if (strcmp(whichphotometer, 'PR705'))
    CMCheckInit(6);
else
    CMCheckInit(1);
end

[window,screenRect] = Screen(0,'OpenWindow',bkgndrgb);
Screen('LoadNormalizedGammaTable',window,linspace(0,1,256)'*ones(1,3));

boxRect = [0 0 boxSize boxSize];
boxRect = CenterRect(boxRect,screenRect);
boxRect = boxRect + [Xoffset -1*Yoffset Xoffset -1*Yoffset];
Screen('FillRect',window,[5 5 5],boxRect);
Screen('Flip',window);

input('');
pause(2);  % Time to leave room

spdMat = zeros(size(RGBs,1), nWavelengthSamples);
for i = 1:length(orderidxs)
    im = ones(boxSize, boxSize / 2, 3);
    for gun = 1:3
        im(:,:,gun) = im(:,:,gun) * RGBs(orderidxs(i),gun);
    end
    tex = Screen('MakeTexture', window, TranslateToColourMode(im, 1));
    Screen('DrawTexture', window, tex, [], boxRect, [], 0);
    Screen('Flip',window);
    if (strcmp(whichphotometer, 'PR705'))
        spd = PR705measspd(S);
    else
        spd = PR650measspd(S);
    end
    spdMat(orderidxs(i),:) =  spdMat(orderidxs(i),:)+spd';
    Screen('Close', tex);
end
spdMat = spdMat./nreps;
sca();
if (strcmp(whichphotometer, 'PR705'))
    PR705close;
else
    PR650close;
end
save junkout spdMat RGBs NFILTERS

load('filtertests.mat');
load('T_rods.mat');
load('wrattennd1.mat');
transmittance1 = newtest3{2,2} ./ newtest3{1,2};
transmittance1 = SplineRaw(S, transmittance1, [380 5 81]);

transmittance2 = wratten_nd1(:,2);

spdMat = SplineSpd(S, spdMat', [380 5 81]);
spdMat_filter1 = spdMat .* repmat(transmittance1 .^ (6-NFILTERS), 1, 3);
spdMat_filter2 = spdMat .* repmat(transmittance2 .^ (6-NFILTERS), 1, 3);

% The two estimates of the filter (empirical, using the fiber optic lamp
% and fromthe Kodak specs) are very similar. GDLH 8/18/12.

% figure();
% subplot(121); hold on;
% plot(spdMat_filter1(:,1), 'color', 'r');
% plot(spdMat_filter1(:,2), 'color', 'g');
% plot(spdMat_filter1(:,3), 'color', 'b');
% plot(spdMat_filter1(:,4), 'color', 'k');
% 
% 
% subplot(122); hold on;
% plot(spdMat_filter2(:,1), 'color', 'r');
% plot(spdMat_filter2(:,2), 'color', 'g');
% plot(spdMat_filter2(:,3), 'color', 'b');
% plot(spdMat_filter2(:,4), 'color', 'k');

scotlum_filter1 = T_rods * spdMat_filter1;
scotlum_filter2 = T_rods * spdMat_filter2;

% Now computing the predictions based on the calibration structure
P_device = SplineSpd([380 4 101], cal.P_device, [380 5 81]);
P_ambient = SplineSpd([380 4 101], cal.P_ambient, [380 5 81]);

threshspect = P_device.*repmat(means,size(P_device,1),1) + repmat(P_ambient,1,3);
threshspect = threshspect .*repmat(transmittance2 .^ 6, 1,3);
% Adding ambient just introductes a vertical shift to everything

scotlum_pred = T_rods * threshspect;

figure(); hold on;
plot(scotlum_filter1(1:3), 'color', 'r');
plot(scotlum_filter2(1:3), 'color', 'b');
plot(scotlum_pred, 'color', 'k');

