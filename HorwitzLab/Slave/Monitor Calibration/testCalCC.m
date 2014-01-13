function [requested, spd] = testCalCC(nRepeats)

ccs = [0.005 0.01 0.015 0.02 0.03 0.04 0.05 0.07 0.09 0.11];
l_ccs = [ccs(:), zeros(length(ccs), 2)];
m_ccs = [zeros(length(ccs), 1), ccs(:), zeros(length(ccs), 1)];
ccs = [0.05 0.01 0.015 0.02 0.03 0.04 0.05 0.07 0.09 0.11 0.15 0.19 0.22 0.28 0.31];
s_ccs = [zeros(length(ccs), 2), ccs(:)];
no_ccs = [0 0 0];
requested = repmat([l_ccs; m_ccs; s_ccs; no_ccs], nRepeats, 1);
idx = randperm(size(requested,1));
requested = requested(idx,:);


%set up the PR650
S = [380, 4, 101];  % Native PR650 settings
spd = nan(101, size(requested,1));
CMCheckInit(1, 'USA19H')
input('Press <Enter> to continue');

%set up the monitor
maxdac = (2^16)-1;
load('T_cones_smj.mat');
load('/Monitor Calibration/Monitor data/Dell 2/Dell2BitsCal.mat');
cal = cals{end};
fundamentals = T_cones_smj;
gammaTable = cal.gammaTable;
bkgndRGB = round(255 .* cal.bgColor)+1; %bkgnd voltages discretized b/w 1&256
bkgndrgb = [gammaTable(bkgndRGB(1),1), gammaTable(bkgndRGB(2),2), gammaTable(bkgndRGB(3),3)];
monSpd = SplineSpd(cal.S_device, cal.P_device, S_cones_smj);
M = fundamentals * monSpd;
invGamma = InvertGamma(cal, 1);
bkgndlms = M * bkgndrgb';
wind = Screen('OpenWindow', 0, bkgndRGB);
textRect = [150, 300, 350, 500];

spd = nan(size(requested,1), 101);
for a = 1:size(requested, 1)
   img = nan(100, 100, 3);
   tCCS = requested(a,:);
   tlms = ((1+tCCS) .* bkgndlms');
   trgb = M \ tlms(:); %left divide instead of inv(M)
   trgb = round(maxdac .* trgb) + 1; %index b/w 1 and 2^16
   tRGB = [invGamma(trgb(1), 1), invGamma(trgb(2), 2), invGamma(trgb(3), 3)];
   tRGB = round(maxdac .* tRGB); %voltage b/w 0 & 2^16
   img(:,:,1) = tRGB(1);
   img(:,:,2) = tRGB(2);
   img(:,:,3) = tRGB(3);
   tex = Screen('MakeTexture', wind, TranslateToColourModeMex(img));
   Screen('DrawTexture', wind, tex, [], textRect, [], 0);
   Screen('Flip', wind);
   tSpd = PR650measspd(S);
   spd(a,:) = tSpd;
   pause(2)
end
Screen('CloseAll')



% analyze the data, start by finding the zero contrast conditions
fundamentals = T_cones_smj';
l_zero = ismember(requested, [0 0 0], 'rows');
spd_zero = spd(l_zero, :);
spd_zero_avg = mean(spd_zero, 1);
bkgndspd = SplineSpd([380, 4, 101], spd_zero_avg(:), S_cones_smj);
bkgndlms = bkgndspd' * fundamentals;
for a = 1:size(spd,1)
    tmp_spd = SplineSpd([380, 4, 101], spd(a,:)', S_cones_smj);
    tmp_lms = tmp_spd' * fundamentals;
    presented(a,:) = (tmp_lms - bkgndlms) ./ bkgndlms;
end

%cycle through the cone isolations and plot (for each cone) the actual vs.
%predicted
signs = sign(requested);
types = [1,0,0; 0,1,0; 0,0,1];
titles = {'L iso', 'M iso', 'S iso'};
colors = {'r.', 'g.', 'b.'};

for a = 1:3
    l_cone = ismember(signs, types(a,:), 'rows');
    pres_cone = presented(l_cone,:);
    req_cone = requested(l_cone,:);
    uniqueRequests = unique(req_cone, 'rows');
    avgs = [];
    stds = [];
    for i = 1:size(uniqueRequests,1)
        l_type = ismember(req_cone, uniqueRequests(i,:), 'rows');
        pres_type = pres_cone(l_type,:);
        avgs(:,i) = mean(pres_type,1)';
        stds(:,i) = std(pres_type, [], 1)';
    end
        
    figure
    for i = 1:3
    subplot(1,3,i)
    errorbar(uniqueRequests(:,a), avgs(i,:), stds(i,:), colors{i});
    ymin = min(min(pres_cone(:,i)), 0);
    ylim([ymin, max(pres_cone(:,a))])
    xlim([0, max(uniqueRequests(:,a))])
    if i == 2
        title(titles{a})
    end
    end
    
    subplot(1,3,a)
    hold on,
    plot([0, max(uniqueRequests(:,a))], [0, max(uniqueRequests(:,a))], 'k')
end



