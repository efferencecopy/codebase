%%
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
%%
load('filtertests.mat');
load('T_rods.mat');
transmittance = newtest3{2,2} ./ newtest3{1,2};
transmittance = SplineRaw([380 2 201], transmittance, [380 5 81]);

% transmittance = ones(81,1);

stro = nex2stro(findfile('Z032312002.nex'));
% stro = nex2stro(findfile('Z022812003.nex'));
colordir_idx = strcmp(stro.sum.trialFields(1,:), 'color_dir');
questmode_idx = strcmp(stro.sum.trialFields(1,:), 'quest_mode');
thresholds = zeros(3,1);

mon_spd = reshape(stro.sum.exptParams.mon_spd, 101, 3);
mon_spd = SplineSpd([380 4 101], mon_spd, [380 5 81]);

for gun = 1:3
    L = stro.trial(:, colordir_idx) == gun;
    all_modes = stro.trial(L, questmode_idx);
    thresholds(gun) = all_modes(end);
end

new_mon_spd = mon_spd .* repmat(transmittance .^ 6, 1, 3);
sensitivity = 1 ./ thresholds;
theory = T_rods * new_mon_spd;
monkey_theory = fund' * new_mon_spd;
monkey_scale = monkey_theory' \ sensitivity;
scalefactor = theory' \ sensitivity;

figure; plot(1:3, theory*scalefactor, 'k.');
hold on; plot(1:3, sensitivity, 'r.');
set(gca, 'xlim', [.5 3.5], 'xtick', 1:3, 'xticklabel', {'R' 'G' 'B'});

figure; plot(1:3, monkey_theory*monkey_scale, 'k.');
hold on; plot(1:3, sensitivity, 'r.');
set(gca, 'xlim', [.5 3.5], 'xtick', 1:3, 'xticklabel', {'R' 'G' 'B'});