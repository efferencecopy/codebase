
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Check to see how much chromatic aberration affects 
% Stimuli in the L/M plane (generated on a CRT).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Setting up...
% Get a white noise data file for the monitor spectra, etc.
stro = nex2stro('N:\NexFiles\Greg\Kali\2009\K022409001.nex');  

% Getting the fundamentals, monitor SPDs, and gamma table
fundamentals = stro.sum.exptParams.fundamentals;
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3])';
mon_spd = stro.sum.exptParams.mon_spd;
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals*mon_spd;
wls = [380:5:780];
gammaTable = stro.sum.exptParams.gamma_table;
gammaTable = reshape(gammaTable, length(gammaTable)/3, 3);

% Getting the background rgb
ridx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_r'));
gidx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_g'));
bidx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_b'));
bkgndRGB = [mode(stro.trial(:,ridx)), mode(stro.trial(:,gidx)), mode(stro.trial(:,bidx))];
bkgndrgb = [gammaTable(bkgndRGB(1)+1,1); gammaTable(bkgndRGB(2)+1,2); gammaTable(bkgndRGB(3)+1,3)];
bkgndlms = M*bkgndrgb;

% Loading the Wandell and Marimont model OTF
if (exist('C:\Documents and Settings\horwitzlab\Desktop\MatlabCode\Analysis\Sandbox\Greg\ChromAb'))
    cd 'C:\Documents and Settings\horwitzlab\Desktop\MatlabCode\Analysis\Sandbox\Greg\ChromAb';
else
    cd 'C:\Matlab Code\Analysis\Sandbox\Greg\ChromAb'
end
load otf;

% Zeroing out everything that exceeds the maximum wavelength considered
% in the model.
L = wls<=max(wave);
wls = wls(L);
mon_spd = mon_spd(L,:);
fundamentals = fundamentals(:,L);

% Check to see how much chromatic aberration affects 
% Stimuli in the L/M plane (generated on a CRT).

nsteps = 21;  % cone contrast steps
cclim = 0.04;
ccs = linspace(-cclim, cclim, nsteps);
data = {};
for SF = 1:6
    conecontrasts = zeros(nsteps,nsteps,3);
    filter = otf(:,sampleSf == SF);
    filter = filter(ismember(wave, wls));
    for idx1 = 1:nsteps;
        for idx2 = 1:nsteps;
            Lcc = ccs(idx1);
            Mcc = ccs(idx2);
            Scc = 0;
            % First the peak of the stimulus
            stimlms = bkgndlms.*[1+Lcc 1+Mcc 1+Scc]';
            stimrgb = inv(M)*stimlms;
            stimspectrum = mon_spd*stimrgb;
            peakwithchromablms = fundamentals'*(stimspectrum.*filter);
            % Then the mean
            stimspectrum = mon_spd*bkgndrgb;
            meanwithchromablms = fundamentals'*(stimspectrum.*filter);
            conecontrasts(idx1, idx2, :) = (peakwithchromablms-meanwithchromablms)./meanwithchromablms;
        end
    end
    n = length(data)+1;
    data{n}.conecontrasts = conecontrasts;
end

intendedconecontrasts = zeros(nsteps,nsteps,3);
intendedconecontrasts(:,:,1) = repmat(ccs',1,nsteps);
intendedconecontrasts(:,:,2) = intendedconecontrasts(:,:,1)';

% Doing the plotting
SF = 3.2;
figure; axes; hold on;
for i = 1:nsteps
    for j = 1:nsteps
        Li = intendedconecontrasts(i,1,1);
        Mi = intendedconecontrasts(1,j,2);
        La = data{SF}.conecontrasts(i,j,1);
        Ma = data{SF}.conecontrasts(i,j,2);
        h1 = plot(Li,Mi,'k.');
        plot([Li La],[Mi Ma],'b-');
        h2 = plot(La,Ma,'r.');
    end
end
xlabel('L-cone contrast');
ylabel('M-cone contrast');
title(['SF = ',num2str(SF),' cyc/deg']);
axis square;
legend([h1; h2],'Intended','Actual','Location','EastOutside');
% An L-M stmiulus has a greater "L" component than one thinks (by a tiny amount)

%%
% A quantiative analysis of detection thresholds in the L/M plane as a function f
% of spatial frequency.  How much do we expect them to change due to chromatic aberration?
% Set SFMULTIPLIER to a value above 1 to artificially increase the spatial
% frequency.
filelist = 'C:\NO BACKUP\NexFiles\Sedna\textFiles\dtControlExpts.txt';
fundamentalsfile = 'N:\toolbox\Psychtoolbox\PsychColorimetricData\PsychColorimetricMatFiles\T_cones_smj';
SFMULTIPLIER = 1;

% Getting the psychophysical data to be modelled
oldfigs = get(0,'Children');
detectionContours(filelist,'mode');
newfigs = get(0,'Children');
for i = newfigs'
    tmp = get(i,'UserData');
    if ~isempty(tmp);
        psychdata = tmp;
        break;
    end
end
psychdata.mean = nanmean(psychdata.data,3);
delete(setdiff(newfigs, oldfigs));

% Getting fundamentals, spectra, gamma functions
[fnames, spikeIdx] = fnamesFromTxt2(filelist);
stro = nex2stro(findfile(char(fnames{end})));
mon_spd = stro.sum.exptParams.mon_spect;
mon_spd = reshape(mon_spd, length(mon_spd)/3, 3);
M = reshape(stro.sum.exptParams.m_mtx,3,3);

fundamentals = load(fundamentalsfile);
tmp = fieldnames(fundamentals);
fundamentals = getfield(fundamentals,tmp{1});

% Sanity check
Mprime = fundamentals*mon_spd;
if ~all(M==Mprime)
    error('Mismatch between m_mtx in stro file and one calculated here');
end

% Getting the background rgb and lms
bkgndrgb = [stro.sum.exptParams.bkgnd_r; stro.sum.exptParams.bkgnd_g; stro.sum.exptParams.bkgnd_b];
bkgndlms = M*bkgndrgb;

% Loading the Wandell and Marimont model OTF
if (exist('C:\Documents and Settings\horwitzlab\Desktop\MatlabCode\Analysis\Sandbox\Greg\ChromAb'))
    cd 'C:\Documents and Settings\horwitzlab\Desktop\MatlabCode\Analysis\Sandbox\Greg\ChromAb';
else
    cd 'C:\Matlab Code\Analysis\Sandbox\Greg\ChromAb'
end
load otf;
wls = [380:5:780];
L = wls<=max(wave);
wls = wls(L);
mon_spd = mon_spd(L,:);
fundamentals = fundamentals(:,L);

% Simulations
for SFidx = 1:length(psychdata.sfs)
    SF = SFMULTIPLIER*psychdata.sfs(SFidx);
    filter1 = interp2(sampleSf,wave,otf,SF,wls);
    filter2 = interp2(sampleSf,wave,otf,SF,wls);
    for coloridx = 1:size(psychdata.colors,1)
        threshold = psychdata.mean(coloridx,SFidx)/100;
        cc = psychdata.colors(coloridx,:)*threshold;

        % First the peak of the stimulus
        stimlms = bkgndlms.*[1+cc]';
        stimrgb = inv(M)*stimlms;
        stimspectrum = mon_spd*stimrgb;
        peakwithchromablms = fundamentals*(stimspectrum.*filter1);
        % Then the mean
        stimspectrum = mon_spd*bkgndrgb;
        meanwithchromablms = fundamentals*(stimspectrum.*filter2); % Do I really want to include filter here?
       % meanwithchromablms = fundamentals*(stimspectrum);
        actualccs(coloridx, SFidx, :) = cc;
        predccs(coloridx, SFidx, :) = (peakwithchromablms-meanwithchromablms)./meanwithchromablms;
    end
end

for i = 1:size(predccs,2)
    figure; axes; hold on;
    plot(predccs(:,i,1), predccs(:,i,2),'m*');
    plot(actualccs(:,i,1), actualccs(:,i,2),'b*');
    plot([predccs(:,i,1) actualccs(:,i,1)]',[predccs(:,i,2) actualccs(:,i,2)]','r-');
    plot(-predccs(:,i,1), -predccs(:,i,2),'m*');
    plot(-actualccs(:,i,1), -actualccs(:,i,2),'b*');
    plot([-predccs(:,i,1) -actualccs(:,i,1)]',[-predccs(:,i,2) -actualccs(:,i,2)]','r-');
    title(psychdata.sfs(i));
end

%%
% Now looking to see how L and M-cone contrasts vary with position in the 
% "isoluminant plane".  Parametrically varying L-M and S contrasts, and 
% seeing how L and M-cone contrast (which should ideally be equal and
% opposite) depend on the stimulus.

nsteps = 21;  % cone contrast steps
cclim = 0.04;
LvsMccs = linspace(-cclim, cclim, nsteps);
Sccs = 5*linspace(-cclim, cclim, nsteps);
data = {};
for SF = 1:10
    conecontrasts = zeros(nsteps,nsteps,3);
    filter = otf(:,sampleSf == SF);
    filter = filter(ismember(wave, wls));
    for idx1 = 1:nsteps;
        for idx2 = 1:nsteps;
            Lcc = sign(LvsMccs(idx1))*sqrt((LvsMccs(idx1)^2)/2);
            Mcc = -Lcc;
            Scc = Sccs(idx2);
            % First the peak of the stimulus
            stimlms = bkgndlms.*[1+Lcc 1+Mcc 1+Scc]';
            stimrgb = inv(M)*stimlms;
            stimspectrum = mon_spd*stimrgb;
            peakwithchromablms = fundamentals*(stimspectrum.*filter);
            % Then the mean
            stimspectrum = mon_spd*bkgndrgb;
            meanwithchromablms = fundamentals*(stimspectrum.*filter);
            conecontrasts(idx2,idx1,:) = (peakwithchromablms-meanwithchromablms)./meanwithchromablms;
            % S is on the rows, so it's the 'y' variable     
        end
    end
    n = length(data)+1;
    data{n}.conecontrasts = conecontrasts;
end
% concontrasts dimensions: (1)L-M contrast (2)S contrast (3)LMS (data)

% Doing the plotting
SF = 3;
figure; hold on;
subplot(2,2,1);
surf(LvsMccs,Sccs,data{SF}.conecontrasts(:,:,1));
title('Actual L-contrast');
subplot(2,2,2);
surf(LvsMccs,Sccs,data{SF}.conecontrasts(:,:,2));
title('Actual M-contrast');
subplot(2,2,3);
surf(LvsMccs,Sccs,data{SF}.conecontrasts(:,:,1)+data{SF}.conecontrasts(:,:,2));
title('L+M (artifact)');
subplot(2,2,4);
surf(LvsMccs,Sccs,data{SF}.conecontrasts(:,:,1)-data{SF}.conecontrasts(:,:,2));
title('L-M');
for i = 1:4
    subplot(2,2,i);
    axis square;
    xlabel('Intended L-M');
    ylabel('Intended S');
    set(gca,'Zlim',[-cclim cclim]);
    colormap(gray);
end
set(gcf,'name',['SF = ',num2str(SF),' cpd']);

[a,b] = gradient(data{SF}.conecontrasts(:,:,1)+data{SF}.conecontrasts(:,:,2));
disp(['L+M gradient is: ',num2str([a(1) b(1)])]);

%%
% A quantiative analysis of detection thresholds in the isoluminant plane as a function of
% of spatial frequency.
% Set SFMULTIPLIER to a value above 1 to artificially increase the spatial
% frequency.

filelist = 'N:\NexFiles\Charlie\Kali\text files\questCSFdata.txt'; %kali
fundamentalsfile = 'N:\toolbox\Psychtoolbox\PsychColorimetricData\PsychColorimetricMatFiles\T_cones_smj';
SFMULTIPLIER = 1;

% Getting the psychophysical data to be modelled
oldfigs = get(0,'Children');
[quest_colors, quese_sfs, quest_data] = questBatchProcess(filelist, 'mode', 10, []);
detectionContours(quest_colors, quese_sfs, quest_data);
newfigs = get(0,'Children');
for i = newfigs'
    tmp = get(i,'UserData');
    if ~isempty(tmp);
        psychdata = tmp;
        break;
    end
end
psychdata.mean = nanmean(psychdata.data,3);
delete(setdiff(newfigs, oldfigs));

% Getting fundamentals, spectra, gamma functions
[fnames, spikeIdx] = fnamesFromTxt2(filelist);
stro = nex2stro(findfile(char(fnames{end})));
mon_spd = stro.sum.exptParams.mon_spect;
mon_spd = reshape(mon_spd, length(mon_spd)/3, 3);
M = reshape(stro.sum.exptParams.m_mtx,3,3);

fundamentals = load(fundamentalsfile);
tmp = fieldnames(fundamentals);
fundamentals = getfield(fundamentals,tmp{1});

% Sanity check
Mprime = fundamentals*mon_spd;
if ~all(M==Mprime)
    error('Mismatch between m_mtx in stro file and one calculated here');
end

% Getting the background rgb and lms
bkgndrgb = [stro.sum.exptParams.bkgnd_r; stro.sum.exptParams.bkgnd_g; stro.sum.exptParams.bkgnd_b];
bkgndlms = M*bkgndrgb;

% Loading the Wandell and Marimont model OTF
if (exist('C:\Documents and Settings\horwitzlab\Desktop\MatlabCode\Analysis\Sandbox\Greg\ChromAb'))
    cd 'C:\Documents and Settings\horwitzlab\Desktop\MatlabCode\Analysis\Sandbox\Greg\ChromAb';
else
    cd 'C:\Matlab Code\Analysis\Sandbox\Greg\ChromAb'
end
load otf;
wls = [380:5:780];
L = wls<=max(wave);
wls = wls(L);
mon_spd = mon_spd(L,:);
fundamentals = fundamentals(:,L);

% Simulations
for SFidx = 1:length(psychdata.sfs)
    SF = SFMULTIPLIER*psychdata.sfs(SFidx);
    filter = interp2(sampleSf,wave,otf,SF,wls);
    for coloridx = 1:size(psychdata.colors,1)
        threshold = psychdata.mean(coloridx,SFidx)/100;
        cc = psychdata.colors(coloridx,:)*threshold;

        % First the peak of the stimulus
        stimlms = bkgndlms.*[1+cc]';
        stimrgb = inv(M)*stimlms;
        stimspectrum = mon_spd*stimrgb;
        peakwithchromablms = fundamentals*(stimspectrum.*filter);
        % Then the mean
        stimspectrum = mon_spd*bkgndrgb;
        meanwithchromablms = fundamentals*(stimspectrum.*filter);
        % meanwithchromablms = fundamentals*(stimspectrum);
        actualccs(coloridx, SFidx, :) = cc;
        predccs(coloridx, SFidx, :) = (peakwithchromablms-meanwithchromablms)./meanwithchromablms;
    end
end

predlvsm = predccs(:,:,1)-predccs(:,:,2);
preds = predccs(:,:,3);
actuallvsm = actualccs(:,:,1)-actualccs(:,:,2);
actuals = actualccs(:,:,3);
for i = 1:size(predccs,2)
    figure; axes; hold on;
    plot(predlvsm(:,i), preds(:,i),'m*');
    plot(actuallvsm(:,i), actuals(:,i),'b*');
    plot([predlvsm(:,i) actuallvsm(:,i)]',[preds(:,i) actuals(:,i)]','r-');
    plot(-predlvsm(:,i), -preds(:,i),'m*');
    plot(-actuallvsm(:,i), -actuals(:,i),'b*');
    plot(-[predlvsm(:,i) actuallvsm(:,i)]',-[preds(:,i) actuals(:,i)]','r-');
    title(psychdata.sfs(i));
    axis equal
end
%%
% Quantifying the amount of luminance contrast (lum = 2L+M) produced by chromatic
% aberration in the intermediate color directions.
% Also for S-cone isolating direction at 3.2 cpd.  Might the effects of microsaccades
% on detection in this condition be due to chromatic aberration?

stro = nex2stro('N:\NexFiles\Greg\Sedna\S051810003.nex');  

% Getting the fundamentals, monitor SPDs, and gamma table
fundamentals = load ('T_cones_smj');
fundamentals = fundamentals.T_cones_smj;
mon_spd = stro.sum.exptParams.mon_spect;
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
M = fundamentals*mon_spd;
wls = [380:5:780];
gammaTable = stro.sum.exptParams.gamma_table;
gammaTable = reshape(gammaTable, length(gammaTable)/3, 3);

% Getting the background rgb
bkgndrgb = [stro.sum.exptParams.bkgnd_r stro.sum.exptParams.bkgnd_g stro.sum.exptParams.bkgnd_b]';
bkgndlms = M*bkgndrgb;

% Loading the Wandell and Marimont model OTF
if (exist('C:\Documents and Settings\horwitzlab\Desktop\MatlabCode\Analysis\Sandbox\Greg\ChromAb'))
    cd 'C:\Documents and Settings\horwitzlab\Desktop\MatlabCode\Analysis\Sandbox\Greg\ChromAb';
else
    cd 'C:\Matlab Code\Analysis\Sandbox\Greg\ChromAb'
end
load otf;

% Zeroing out everything that exceeds the maximum wavelength considered
% in the model.
L = wls<=max(wave);
wls = wls(L);
mon_spd = mon_spd(L,:);
fundamentals = fundamentals(:,L);

[thresholds, colorDirs, sfs] = DTquestUnpackGH(stro, 'mode')
filter = otf(:,sampleSf == ceil(sfs(1)));
filter = otf(:,sampleSf == 10);
filter = filter(ismember(wave, wls));
lum = [];
vlambda = [2 1 0]*fundamentals;  % Approximate
for i = 1:size(colorDirs,1)
    c = colorDirs(i,:)./norm(colorDirs(i,:));  % make it into a unit vector
    cc = thresholds(i,1)*c;
  
    stimlms = bkgndlms.*[1+cc(1) 1+cc(2) 1+cc(3)]';
    stimrgb = inv(M)*stimlms;
    stimspectrum = mon_spd*stimrgb;
    peakwithchromab = vlambda*(stimspectrum.*filter);
    peakwochromab = vlambda*stimspectrum;
    
    stimspectrum = mon_spd*bkgndrgb;
    meanwithchromab = vlambda*(stimspectrum.*filter);
    meanwochromab = vlambda*stimspectrum;
    a = peakwithchromab-peakwochromab;
    b = meanwithchromab-meanwochromab;
    
    %lum(i) = (peakwithchromab-meanwithchromab)./meanwithchromab;
    %lum(i) = lum(i) - (peakwochromab-meanwochromab)./meanwochromab;
    % Above, correction for 'intended' luminance artifact (because [1 -1 0] has a non-zero projection onto [2 1 0])
    %lum(i) = (a-b)/b
    lum(i) = (a-b)./mean([meanwithchromab, meanwochromab]);
end
% Lum should be in the same units as in DTquestUnpackGH.m (vector norms in
% cone contrast space).

[colorDirs, lum']
% The L-M+S direction has the largest luminance artifacts from chromatic
% aberration, but they only amount to < 1% contrast at detection
% threshold.  Is this sufficient to halve detection threshold?

%%
% This isn't code relating to chromatic aberration, but I'm curious to see
% the directions of the rod null in the isoluminant plane.

load T_rods.mat

% Get a white noise data file for the monitor spectra, etc.
stro = nex2stro('N:\NexFiles\Kali\K022409001.nex');  

% Getting the fundamentals, monitor SPDs, and gamma table
fundamentals = stro.sum.exptParams.fundamentals;
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = stro.sum.exptParams.mon_spd;
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;
wls = [380:5:780];
gammaTable = stro.sum.exptParams.gamma_table;
gammaTable = reshape(gammaTable, length(gammaTable)/3, 3);

% Getting the background rgb
ridx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_r'));
gidx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_g'));
bidx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_b'));
bkgndRGB = [mode(stro.trial(:,ridx)), mode(stro.trial(:,gidx)), mode(stro.trial(:,bidx))];
bkgndrgb = [gammaTable(bkgndRGB(1)+1,1); gammaTable(bkgndRGB(2)+1,2); gammaTable(bkgndRGB(3)+1,3)];
bkgndlms = M*bkgndrgb;

% Looking around in the isoluminant plane and determining
% the rod quantal catch rate for each stimulus.
nsteps = 21;  % cone contrast steps
cclim = 0.04;
LvsMccs = linspace(-cclim, cclim, nsteps);
Sccs = 5*linspace(-cclim, cclim, nsteps);
rodexcitations = zeros(nsteps,nsteps);
for idx1 = 1:nsteps;
    for idx2 = 1:nsteps;
        Lcc = LvsMccs(idx1)/2;
        Mcc = -Lcc;
        Scc = Sccs(idx2);
        % First the peak of the stimulus
        stimlms = bkgndlms.*[1+Lcc 1+Mcc 1+Scc]';
        stimrgb = inv(M)*stimlms;
        stimspectrum = mon_spd*stimrgb;
        rodexcitations(idx2,idx1) = T_rods*stimspectrum;
        % S is on the rows, so it's the 'y' variable
    end
end
figure;
axes; hold on;
surf(LvsMccs,Sccs,rodexcitations);
set(gca,'View',[0 90]);
ylabel('S contrast');
xlabel('L-M contrast');
title('Rod excitation');

[a,b] = gradient(rodexcitations);
disp(['gradient is: ',num2str([a(1) b(1)])]);
disp(['null direction is: ',num2str([-b(1) a(1)])]);



