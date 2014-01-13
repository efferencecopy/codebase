
%% import the data
habitTxt = 'N:\NexFiles\Charlie\Apollo\text files\habit_lmplane.txt';
noHabitTxt = 'N:\NexFiles\Charlie\Apollo\text files\quest_lmplane.txt';
fitMethod = 'mode';
perfRange = [];
nTrials = 20;
[hab_colors, hab_sfs, hab_data] = questBatchProcess(habitTxt, fitMethod, nTrials, perfRange);
[nohab_colors, nohab_sfs, nohab_data] = questBatchProcess(noHabitTxt, fitMethod, nTrials, perfRange);



%% plot detection contours
%calculate the mean sensitivity
hab_thresh = nanmean(hab_data, 3);
hab_SEM = nanstd(hab_data, [], 3)./sqrt(sum(~isnan(hab_data),3));
nohab_thresh = nanmean(nohab_data, 3);
nohab_SEM = nanstd(nohab_data, [],3)./sqrt(sum(~isnan(nohab_data),3));

%first plot sensitivities in the presence of the habituator
thetas = atan(hab_colors(:,2)./hab_colors(:,1));
Lthresh = hab_thresh .* cos(thetas);
Mthresh = hab_thresh .* sin(thetas);
Lsem_lo = (hab_thresh-hab_SEM) .* cos(thetas);
Lsem_hi = (hab_thresh+hab_SEM) .* cos(thetas);
Msem_lo = (hab_thresh-hab_SEM) .* sin(thetas);
Msem_hi = (hab_thresh+hab_SEM) .* sin(thetas);
f = figure;
hold on,
plot(Lthresh, Mthresh, 'ko', 'markerfacecolor', 'k')
plot(-Lthresh, -Mthresh, 'ko', 'markerfacecolor', 'k')
plot([Lsem_lo(:), Lsem_hi(:)]', [Msem_lo(:), Msem_hi(:)]', 'k-')
plot([-Lsem_lo(:), -Lsem_hi(:)]', [-Msem_lo(:), -Msem_hi(:)]', 'k-')

%next plot sensitivities in the absence of the habituator
thetas = atan(nohab_colors(:,2)./nohab_colors(:,1));
Lthresh = nohab_thresh .* cos(thetas);
Mthresh = nohab_thresh .* sin(thetas);
Lsem_lo = (nohab_thresh-nohab_SEM) .* cos(thetas);
Lsem_hi = (nohab_thresh+nohab_SEM) .* cos(thetas);
Msem_lo = (nohab_thresh-nohab_SEM) .* sin(thetas);
Msem_hi = (nohab_thresh+nohab_SEM) .* sin(thetas);
plot(Lthresh, Mthresh, 'ro', 'markerfacecolor', 'r')
plot(-Lthresh, -Mthresh, 'ro', 'markerfacecolor', 'r')
plot([Lsem_lo(:), Lsem_hi(:)]', [Msem_lo(:), Msem_hi(:)]', 'r-')
plot([-Lsem_lo(:), -Lsem_hi(:)]', [-Msem_lo(:), -Msem_hi(:)]', 'r-')
hold off

%add a line of unity slope
hold on,
lvmintercept = 0.75;
lumintercept = 2.25;
x = -2:1:2;
y1 = x+lvmintercept;
y2 = x-lvmintercept;
y3 = -x + lumintercept;
y4 = -x - lumintercept;
plot([x(:), x(:)],[y1(:), y2(:)], 'r:')
plot([x(:), x(:)],[y3(:), y4(:)], 'k:')
% plot([-2 2], [0 0], 'k')
% plot([0 0], [-2 2], 'k')
xlim([-2 2])
ylim([-2 2])
axis square





