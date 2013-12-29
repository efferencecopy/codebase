% Code to do analyses and make figures for the SFN2011 poster on the
% orientation of L+M planes in different color spaces (2 vs 10 degree
% fundmantals)
%
% Contents:
% 1) Trying to recreate the SMJ 2 degree fundamentals from the 10
% degree ones.
%
% 2) Plotting cone weights as a function of macular pigment density
% (starting with the 10 degree SMJ cone fundamentals).
%
% 2.1) Looking for an effect of eccentricity on S-cone weight measured with
% the two degree fundamentals
%
% 2.2) Looking for an effect of eccentricity on the macular pigment density
% that predicts zero S-cone weight. Kind of the inverse analysis from 2.1
% above.
%
% 3) Plotting S-cone weight as a function of macular pigment density and
% lens density (jointly).
%
% 4) Orientations of ellipsoid long axes as a function of macular pigment
% density 
%
% 5) Looking at isodetection surfaces using different amounts of macular
% pigment. Wasn't very informative.
%
% 6) How much does the change in pigment self-screening change S-cone
% weight to planar L+M cells?
%
% 7) A couple of views of isoresponse planes from an L+M neuron. Planes
% calculated with different assumed macular pigment densities. POSTER READY.
%
% 7.1) A couple of view of isoresponse planes from an L+M neuron. 2-D plots
% (aL+bM vs. S). POSTER READY.
%
% 8) Plot of mean S-cone weight as a function of peak macular pigment. Also
% A few histograms of S-cone weights.  Planar L+M cells. 
% Based on section(2). POSTER READY.
%
% 9) S-cone component of largest eigenvector of ellipsoids across a
% population of cells. Analogous to section 8, based on section 4. POSTER
% READY.
%
% 10) A couple of views of isoresponse ellipsoids. POSTER READY.
%
% 11) Eccentricity regression. Based on section 2.2. OBSOLETE
%
% 12) What combinations of macular pigment density and lens density
% gives a green to blue luminance ratio of ~3.6? This is roughly 
% what the monkeys show at large eccentricities. Analysis
%
% 13) POSTER READY figure. blum/glum macular pigment density estimation.
% IN PROGRESS
%
% 14) A horizontally oriented gabor for the psychophysics methods
%
% 15) Example Gabors from the best-fit plane from the example neuron
%% Section 1)
% Trying to recreate the SMJ 2 degree fundamentals from the 10 degree ones
load ('T_cones_smj');
load ('T_cones_smj10');
load ('den_mac_ws');
load ('den_lens_ws');

% dens_lens_smj:  From cvrl.org This is *not* the same as SMJ Table 7 
% nor is it the same at the SMJ Table 7 *1.16 (small pupil adjustment)
% But it's close to the latter. CVRL says "These values are for a small
% pupil, for an open pupil divide by 1.16."
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

% Van Norren and Vos (This came from the original 1947 paper)
den_lens_vnv = [2.61 2.31 2.02 1.73 1.45 1.19 0.96 0.76 0.59 0.47 0.37 0.31 0.27 0.24...
    0.225 0.205 0.195 0.185 0.175 0.165 0.155 0.145 0.14 0.13 0.125 0.115 0.11 0.105...
    0.1 0.095 0.09 0.085 0.08 0.075 0.07 0.065 0.06 0.055 0.05 0.045 0.04 0.035...
    0.031 0.027 0.024 0.021 0.018 0.015 0.012 0.01 0.008 0.006 0.004 0.003 0.002 0.001 zeros(1,25)]';



fund10 = T_cones_smj10'./repmat(max(T_cones_smj10'),81,1);
fund2 = T_cones_smj'./repmat(max(T_cones_smj'),81,1);

macpig = den_mac_ws;
lens = den_lens_vnv;

% 
% macpigtransmittance = 1./(10.^(macpig*macpigweight));
% lenstransmittance = 1./(10.^(lens*lensweight));
% 
% synthfund = fund10.*repmat(macpigtransmittance,1,3).*repmat(lenstransmittance,1,3);
% synthfund = synthfund./repmat(max(synthfund),81,1);
% 
% figure; axes; hold on; 
% plot(synthfund);
% plot(T_cones_smj','k.');
% 
% % Recreating the figure from Knowles and Dartnall p.61
% % Note, their y-axis has a typo!
% figure;
% Iinc = 1;
% Itr = linspace(0,1,200);
% plot(log10(Iinc./Itr),(Iinc-Itr)./Iinc,'k.')

%---------------
% Starting the with 10 deg fundamentals are converting to "action spectra"
macpigtransmittance = 1./(10.^(macpig*.28));
lenstransmittance = 1./(10.^(lens*1.28));
opticaldensity = 0.3;

% Undoing preretinal filters
absorptance = fund10./repmat(macpigtransmittance,1,3)./repmat(lenstransmittance,1,3); % yikes!  Exceeds 1
absorptance = absorptance./repmat(max(absorptance),81,1);

% These two lines below are equivalent
%actionspectra = log10(1./(1-fund10.*opticaldensity));
%actionspectra= -log10(1-fund10.*opticaldensity);
actionspectra = -log10(1-absorptance*(1-10^-opticaldensity));
actionspectra = actionspectra/opticaldensity;

% Undoing the manipulations exactly
funds = 1-10.^(-actionspectra.*opticaldensity);
funds = funds.*repmat(macpigtransmittance,1,3).*repmat(lenstransmittance,1,3);
funds = funds./repmat(max(funds),81,1);
figure; axes; hold on;
plot(funds-fund10)
% Good, exactly the same!

% Now going from 10 to 2 degree fundamentals as best I can
macpigtransmittance2 =  1./(10.^(macpig*.70));
lenstransmittance = 1./(10.^(lens*1.28));
opticaldensity2 = 0.4;

synthfunds2 = 1-10.^(-actionspectra.*opticaldensity2);
funds = 1-10.^(-actionspectra.*opticaldensity2);
funds = funds.*repmat(macpigtransmittance2,1,3).*repmat(lenstransmittance,1,3);
funds = funds./repmat(max(funds),81,1);
figure; axes; hold on;
plot(funds-fund2)
set(gca,'YLim',[-0.01 0.01])

%%
% Section 2
% Cone weights of L+M neurons as a function of added macular pigment
% density.
load ('T_cones_smj10');
load ('den_mac_ws');
macpig = den_mac_ws;  % Peaks at 0.5
mult10deg = 0.28; % from SMJ paper (Don't change, used below)
mult2deg = 0.7; % from SMJ paper
macpigtransmittance = 1./(10.^(macpig.*mult10deg));
fund10 = T_cones_smj10'./repmat(max(T_cones_smj10'),81,1);
fundnomacpig = fund10./repmat(macpigtransmittance,1,3);
fundnomacpig = fundnomacpig./repmat(max(fundnomacpig),81,1);

% First gathering all the data
data = [];
[fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\NTlum2.txt');
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
    
    M10 = T_cones_smj10*mon_spd;
    Loog = logical(out(:,7));
    ecc = sqrt((NT.sum.exptParams.rf_x/10).^2 + (NT.sum.exptParams.rf_y/10).^2);
    data(cellcounter).M = M;
    data(cellcounter).data = out;
    data(cellcounter).bkgndrgb = NT.sum.exptParams.bkgndrgb;
    data(cellcounter).ecc = ecc;
    data(cellcounter).monspd = mon_spd;
end

% Now, computing cone weights as a function of macular pigment density
macpigmults = linspace(0,mult2deg,10);
for i = 1:length(data)
    i
    cc_original = data(i).data(:,[2:4]).*repmat(data(i).data(:,5),1,3);
    Loog = data(i).data(:,end);
    allconeweights = [];
    for j = 1:length(macpigmults)
        j
        macpigtransmittance = 1./(10.^(macpig.*macpigmults(j)));
        tmpfund = fundnomacpig.*repmat(macpigtransmittance,1,3);
        tmpfund = tmpfund./repmat(max(tmpfund),81,1);
        M = tmpfund'*mon_spd;
        scaled = ConvertConeContrastBasis(data(i).M, M, data(i).bkgndrgb, cc_original);
        [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);
        coneweights = (xformmat*planeparams)';
        if (sign(coneweights(1)) == -1)
            coneweights = -coneweights;
        end
        allconeweights = [allconeweights; coneweights];
    end
    data(i).coneweights = allconeweights;
end

% DKL plot
figure;
for j = 1:length(macpigmults)
    subplot(ceil(sqrt(length(macpigmults))),ceil(sqrt(length(macpigmults))),j);
    axis square;
    plot([0 1],[1 0],'k-');
    hold on;
    title(num2str(macpigmults(j)));
    for i = 1:length(data)
        normconeweights = data(i).coneweights(j,:)./sum(abs(data(i).coneweights(j,:)));
        if (normconeweights(3) > 0)
            h = plot(normconeweights(1),normconeweights(2),'ko','MarkerFaceColor','black');
        else
            h = plot(normconeweights(1),normconeweights(2),'ko','MarkerFaceColor','white');            
        end
    end
end

% Hist of S-cone weights
figure;
allconeweights = cat(3,data(:).coneweights);
bins = [-.4:.05:.4];
mn = []; stdev = [];
for j = 1:length(macpigmults)
    subplot(ceil(sqrt(length(macpigmults))),ceil(sqrt(length(macpigmults))),j);
    hold on;
    tmpconeweights = squeeze(allconeweights(j,:,:));
    normconeweights = tmpconeweights./repmat(sum(abs(tmpconeweights)),3,1);
    hist(normconeweights(3,:),bins);
    mn(j) = mean(normconeweights(3,:));
    stdev(j) = std(normconeweights(3,:));
    plot(mn(j),0,'m^');
    [~,p] = ttest(normconeweights(3,:));
    title([num2str(mn(j)),' p = ',num2str(p,1)]);
end

% At what macular pigment density do we have zero S-cone weight?
figure; axes; hold on;
errorbar(macpigmults,mn, stdev./sqrt(size(normconeweights,2)));
plot(macpigmults,mn,'k.');
xlabel('macular pigment density multiplier');
ylabel('Mean S-cone weight');
b = regress(mn', [macpigmults' ones(length(mn),1)])
plot([macpigmults(1) macpigmults(end)],b(1).*[macpigmults(1) macpigmults(end)]+b(2),'k:');
title(['Macular pigment density at 0 S-cone weight: ',num2str(-b(2)/b(1))]);
plot([-.1 -b(2)/b(1)],[0 0],'k-','Linewidth',2);
h(1) = plot([-b(2)/b(1) -b(2)/b(1)],[0 -.02],'k-','Linewidth',2);

plot([-.1 mult10deg],[b(1)*mult10deg+b(2) b(1)*mult10deg+b(2)],'g-');
h(2) = plot([mult10deg mult10deg],[b(1)*mult10deg+b(2) -.02],'g-');

plot([-.1 mult2deg],[b(1)*mult2deg+b(2) b(1)*mult2deg+b(2)],'r-');
h(3) = plot([mult2deg mult2deg],[b(1)*mult2deg+b(2) -.02],'r-');
legend(h,{'Zero S-cone weight','10 degree','2 degree'});


% Proportion L-cones (ignoring S-cones)
figure;
allconeweights = cat(3,data(:).coneweights);
bins = [0:.1:1];
mn = [];
for j = 1:length(macpigmults)
    subplot(ceil(sqrt(length(macpigmults))),ceil(sqrt(length(macpigmults))),j);
    hold on;
    tmpconeweights = squeeze(allconeweights(j,:,:));
   % normconeweights = tmpconeweights./repmat(sum(abs(tmpconeweights)),3,1);
    lconeprop = abs(tmpconeweights(1,:))./sum(abs(tmpconeweights([1 2],:)));
    hist(lconeprop,bins);
    mn = mean(lconeprop);
    plot(mn,0,'m^');
    [~,p] = ttest(lconeprop-.5);
    title([num2str(mn),' p = ',num2str(p,1)]);
    set(gca,'XLim',[0 1]);
end

%%
% Section 2.1
% Apparent S-cone weight as a function of eccentricity
% Using 2-degree fundamentals shows no significant relationship, but the
% right trend (more S-cone weight at higher eccentricities where there is
% less macular pigment). OBSOLETE

% First gathering all the data
load ('T_cones_smj');
load ('T_cones_smj10');
load ('den_mac_ws');
data = [];
[fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\NTlum2.txt');
for cellcounter = 1:size(fnames,1)
    filename = findfile(char(fnames{cellcounter}));
    paradigmID = getparadigmID(filename);
    if paradigmID == 103
        NT = nex2stro(filename);
    end
    out = NTpreprocess(NT,0,Inf);
    fundamentals = NT.sum.exptParams.fundamentals;
    fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
    mon_spd = NT.sum.exptParams.mon_spd;
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
    M = fundamentals'*mon_spd;
    M2 = T_cones_smj*mon_spd;
    scaled = ConvertConeContrastBasis(M, M2, NT.sum.exptParams.bkgndrgb, out(:,[2:4]).*repmat(out(:,5),1,3));
    Loog = logical(out(:,7));
    ecc = sqrt((NT.sum.exptParams.rf_x/10).^2 + (NT.sum.exptParams.rf_y/10).^2);
    data(cellcounter).M = M;
    data(cellcounter).scaled = scaled;
    data(cellcounter).Loog = Loog;
    data(cellcounter).bkgndrgb = NT.sum.exptParams.bkgndrgb;
    data(cellcounter).ecc = ecc;
    data(cellcounter).monspd = mon_spd;
    fn = char(fnames{cellcounter});
    data(cellcounter).monk = fn(1);
end

for i = 1:length(data)
    i
    [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(data(i).scaled, data(i).Loog);
    coneweights = (xformmat*planeparams)';
    if (sign(coneweights(1)) == -1)
        coneweights = -coneweights;
    end
    data(i).deg2coneweights = coneweights';
end

iskali = logical([data.monk]=='K');
ecc = [data.ecc];
coneweights = [data.deg2coneweights];
normconeweights = coneweights./repmat(sum(abs(coneweights)),3,1);
figure; axes; hold on;
plot(ecc(iskali),normconeweights(3,iskali),'ko','MarkerFaceColor','black');
plot(ecc(~iskali),normconeweights(3,~iskali),'ko','MarkerFaceColor','red');
[b,bint,r,rint,stats] = regress(normconeweights(3,:)', [ecc' ones(length(ecc),1)]);
plot([min(ecc) max(ecc)],b(1)*[min(ecc) max(ecc)]+b(2),'g:');
xlabel('Eccentricity (deg)');
ylabel('S-cone weight');
title('Calculated with the SMJ 2 degree cone fundamentals');
[r,p] = corr([ecc',normconeweights(3,:)'],'Type','Spearman')

% Now I need to come up with some kind of reaonable prediction.
% I'm going to start with the assumption of 0 S-cone pigment at 10 degrees
% and some reasonable level at 2 degrees.
fund2 = T_cones_smj'./repmat(max(T_cones_smj'),81,1);
fund10 = T_cones_smj10'./repmat(max(T_cones_smj10'),81,1);
macpigmult10deg = 0.28; % from SMJ paper (Don't change, used below)
macpigtransmittance = 1./(10.^(den_mac_ws.*macpigmult10deg));
fundnomacpig = fund10./repmat(macpigtransmittance,1,3);
fundnomacpig = fundnomacpig./repmat(max(fundnomacpig),81,1);


% let's pretend the neurons are tuned for L+M exactly
% Finding preferred color direction in gun space
Mnomacpig = fundnomacpig'*mon_spd;
M10 = fund10'*mon_spd;
M2 = fund2'*mon_spd;
% % Now finding an RGB that isolates an L+M (contrast) mechanism
% lms2bkgnd = M2*NT.sum.exptParams.bkgndrgb;
% % RGBs needed to provide a [1 1 0] contrast stimulus 
% lumisorgb =inv(M2)*([2*lms2bkgnd(1); 2*lms2bkgnd(2); lms2bkgnd(3)]);
% % Taking that light and seeing how much it hits the cones (without macular
% % pigment)
% lmsbkgnd = Mnomacpig*NT.sum.exptParams.bkgndrgb;
% lms = Mnomacpig*lumisorgb;
% cc2 = (lms-lmsbkgnd)./lmsbkgnd;
% % cc2(3)./sum(abs(cc2))
% % cc2 = cc2./max(cc2)
% % % Sheesh, the artifact seems incredibly small

% % Imagine a cell that has cone weights [1 1 0]
% Leftmost "M" is the gun weights assuming some degree of macular pigmentation
% Rightmost "M"  converts these gun weights back to cone weights, but 
% assumes a different degree of macular pigmentation
newconeweights = inv(M2')*Mnomacpig'*[1 1 0]'
newconeweights = newconeweights./sum(abs(newconeweights));
h(1) = plot([10 0],[newconeweights(3) 0],'m:');
% Line of predicted S-cone weights assuming an L+M cell and zero macular
% pigment at 10 degrees and the two degree fundamentals (a multiplier of
% 0.7) at 0 degree eccentricity.

%%
% Section 2.2
% What macular pigment density would give zero S-cone weight?

load ('T_cones_smj10');
load ('den_mac_ws');
% downloaded from cvrl.org 
den_lens_smj = [3.2 2.8 2.4 1.985 1.62 1.308 1.05 0.84 0.675 0.557 0.468 0.393 0.335 0.29 0.26 0.24 0.225 0.215 0.203 0.191 0.18 0.168 0.162 0.151 0.145 0.133...
        0.128 0.122 0.116 0.11 0.104 0.099 0.093 0.087 0.081 0.075 0.07 0.064 0.058 0.052 0.046 0.041 0.036 0.031 0.028 0.024 0.021 0.017 0.014 0.012 0.009 0.007...
        0.005 0.003 0.002 0.001 zeros(1,25)]';
% First two numbers of den_lens_smj are extrapolations.

lens = den_lens_smj./1.16; % Dividing by 1.16 to make it open pupil 
macpig = den_mac_ws;  % Peaks at 0.5
macpigmult10deg = 0.28; % from SMJ paper (Don't change, used below)
macpigmult2deg = 0.7; % from SMJ paper
macpigtransmittance = 1./(10.^(macpig.*macpigmult10deg));
fund10 = T_cones_smj10'./repmat(max(T_cones_smj10'),81,1);
fundnopig = fund10./repmat(macpigtransmittance,1,3);
fundnopig = fundnopig./repmat(max(fundnopig),81,1);

% using a new (synthetic) set of cone fundamentals
load ('T_cones_synthgh1');
fundnopig = T_cones_synthgh1;

% First gathering all the data
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
    
    M10 = T_cones_smj10*mon_spd;
    Loog = logical(out(:,7));
    ecc = sqrt((NT.sum.exptParams.rf_x/10).^2 + (NT.sum.exptParams.rf_y/10).^2);
    data(cellcounter).M = M;
    data(cellcounter).data = out;
    data(cellcounter).bkgndrgb = NT.sum.exptParams.bkgndrgb;
    data(cellcounter).ecc = ecc;
    data(cellcounter).monspd = mon_spd;
    fn = char(fnames{cellcounter});
    data(cellcounter).monk = fn(1);
end

% Now, computing cone weights as a function of macular pigment density
% This takes a while to run
macpigmults = linspace(-.5,1,20);
for i = 1:length(data)
    cc_original = data(i).data(:,[2:4]).*repmat(data(i).data(:,5),1,3);
    Loog = data(i).data(:,end);
    allconeweights = [];
    for j = 1:length(macpigmults)
        macpigtransmittance = 1./(10.^(macpig.*macpigmults(j)));
        tmpfunds = fundnopig.*repmat(macpigtransmittance,1,3);
        tmpfunds = tmpfunds./repmat(max(tmpfunds),81,1);
        M = tmpfunds'*mon_spd;
        scaled = ConvertConeContrastBasis(data(i).M, M, data(i).bkgndrgb, cc_original);
        [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);
        coneweights = (xformmat*planeparams)';
        if (sign(coneweights(1)) == -1)
            coneweights = -coneweights;
        end
        allconeweights = [allconeweights; coneweights];
    end
    data(i).coneweights = allconeweights;
end

predictedmacpig = [];
for i = 1:length(data)
    predictedmacpig(i) = interp1(data(i).coneweights(:,3),macpigmults,0,'linear')
end

iskali = logical([data.monk]=='K');
ecc = [data.ecc]
figure; axes; hold on;
plot(ecc(iskali), predictedmacpig(iskali)*max(macpig),'ko','MarkerFaceColor','black','MarkerSize',7);
plot(ecc(~iskali), predictedmacpig(~iskali)*max(macpig),'ro','MarkerFaceColor','red','MarkerSize',7);

xlabel('eccentricity (deg)');
ylabel('predicted peak macular pigment density (460 nm)');
L = ~isnan(predictedmacpig) & predictedmacpig > min(macpigmults) & predictedmacpig < max(macpigmults);
cordata = [[data.ecc]' predictedmacpig'];
cordata = cordata(L,:);
[r,p] = corr(cordata,'type','Pearson');
title(['r = ',num2str(r(1,2)),' p = ',num2str(p(1,2))]);

[b,bint,r,rint,stats] = regress([predictedmacpig(L)*max(macpig)]', [ecc(L)' ones(sum(L),1)]);
h = plot([min(ecc) max(ecc)],b(1)*[min(ecc) max(ecc)]+b(2),'k:');

set(h,'LineWidth',2,'LineStyle',':');
plot([2 7],[1 1]*macpigmult2deg*max(macpig),':')
plot([2 7],[1 1]*macpigmult10deg*max(macpig),':')
plot([2 7],[0 0],':')
set(gca,'YLim',[-.2 .4]);

%% 
% Section 3
% S-cone weight as a function of macular pigment density and lens density

load ('T_cones_smj10');
load ('den_mac_ws');
%den_lens_vnv = [2.61 2.31 2.02 1.73 1.45 1.19 0.96 0.76 0.59 0.47 0.37 0.31 0.27 0.24...
%    0.225 0.205 0.195 0.185 0.175 0.165 0.155 0.145 0.14 0.13 0.125 0.115 0.11 0.105...
%    0.1 0.095 0.09 0.085 0.08 0.075 0.07 0.065 0.06 0.055 0.05 0.045 0.04 0.035...
%    0.031 0.027 0.024 0.021 0.018 0.015 0.012 0.01 0.008 0.006 0.004 0.003 0.002 0.001 zeros(1,25)]';
% downloaded from cvrl.org 
den_lens_smj = [3.2 2.8 2.4 1.985 1.62 1.308 1.05 0.84 0.675 0.557 0.468 0.393 0.335 0.29 0.26 0.24 0.225 0.215 0.203 0.191 0.18 0.168 0.162 0.151 0.145 0.133...
        0.128 0.122 0.116 0.11 0.104 0.099 0.093 0.087 0.081 0.075 0.07 0.064 0.058 0.052 0.046 0.041 0.036 0.031 0.028 0.024 0.021 0.017 0.014 0.012 0.009 0.007...
        0.005 0.003 0.002 0.001 zeros(1,25)]';
% First two numbers of den_lens_smj are extrapolations. Might want to do
% this better

lens = den_lens_smj./1.16; % Dividing by 1.16 to make it open pupil 
macpig = den_mac_ws;  % Peaks at 0.5
macpigmult10deg = 0.28; % from SMJ paper (Don't change, used below)
macpigmult2deg = 0.7; % from SMJ paper
lensmult = 1.28; % for SMJ paper
macpigtransmittance = 1./(10.^(macpig.*macpigmult10deg));
lenstransmittance = 1./(10.^(lens.*lensmult));
fund10 = T_cones_smj10'./repmat(max(T_cones_smj10'),81,1);
fundnopig = fund10./repmat(macpigtransmittance,1,3)./repmat(lenstransmittance,1,3);
fundnopig = fundnopig./repmat(max(fundnopig),81,1);
opticaldensity = 0.3;
actionspectra = -log10(1-fundnopig*(1-10^-opticaldensity));
actionspectra = actionspectra/opticaldensity;

% Making new fundamentals
% Here we can change the optical density
opticaldensity = 0.2;
fundnopig = 1-10.^(-actionspectra.*opticaldensity);
fundnopig = fundnopig./repmat(max(fundnopig),81,1);

% First gathering all the data
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
    
    M10 = T_cones_smj10*mon_spd;
    Loog = logical(out(:,7));
    ecc = sqrt((NT.sum.exptParams.rf_x/10).^2 + (NT.sum.exptParams.rf_y/10).^2);
    data(cellcounter).M = M;
    data(cellcounter).data = out;
    data(cellcounter).bkgndrgb = NT.sum.exptParams.bkgndrgb;
    data(cellcounter).ecc = ecc;
    data(cellcounter).monspd = mon_spd;
end

% Now, computing cone weights as a function of macular pigment 
% and lens density

macpigmults = linspace(0,macpigmult2deg,8);
lensmults = linspace(.3,1.5,8);  % 1.28 is SMJ value

for i = 1:length(data)
    cc_original = data(i).data(:,[2:4]).*repmat(data(i).data(:,5),1,3);
    Loog = data(i).data(:,end);
    allconeweights = [];
    for j = 1:length(macpigmults)
        for k = 1:length(lensmults)
            macpigtransmittance = 1./(10.^(macpig.*macpigmults(j)));
            lenstransmittance = 1./(10.^(lens.*lensmults(k)));
            tmpfunds = fundnopig.*repmat(macpigtransmittance,1,3).*repmat(lenstransmittance,1,3);
            tmpfunds = tmpfunds./repmat(max(tmpfunds),81,1);
            M = tmpfunds'*mon_spd;
            scaled = ConvertConeContrastBasis(data(i).M, M, data(i).bkgndrgb, cc_original);
            [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);
            coneweights = (xformmat*planeparams)';
            if (sign(coneweights(1)) == -1)
                coneweights = -coneweights;
            end
            allconeweights = [allconeweights; coneweights];
        end
    end
    data(i).coneweights = allconeweights;
end
allconeweights = cat(3,data(:).coneweights);

% Distribution of normalized S-cone weights
bins = [-.4:.05:.4];
norm_s = squeeze(allconeweights(:,3,:)./sum(abs(allconeweights(:,:,:)),2));
norm_s_mat = reshape(norm_s,length(lensmults),length(macpigmults),length(data));
% First column is 1st macpig, second column is 2nd macpig, ...
figure;
for i = 1:length(macpigmults)
    for j = 1:length(lensmults)
        subplot(length(lensmults),length(macpigmults),(j-1)*length(lensmults)+i);
        hist(squeeze(norm_s_mat(i,j,:)),bins)
        mn = mean(squeeze(norm_s_mat(i,j,:)));
        [~,p] = ttest(squeeze(norm_s_mat(i,j,:)));
        if (p< 0.05);
            title(['mean = ',num2str(mn),'*']);
        else
            title(['mean = ',num2str(mn)]);            
        end
    end
    drawnow
end
% Mac pig on the rows lens of the columns

figure; 
contour(abs(mean(norm_s_mat,3)),20);
set(gca,'ytick',1:length(lensmults),'ytickLabel',num2str(lensmults',2))
set(gca,'xtick',1:length(macpigmults),'xtickLabel', num2str(macpigmults',2))
ylabel('Lens multiplier'); xlabel('macular pigment multiplier')
axis xy;

% To minimize S-cone weight we want to multiply the lens density by 0.75
% and the macular pigment density by 0.35
% macpig = 0, lens = 0.5 works pretty well too.

%%
% Section 4
% The long axis orientation of ellipsoidal isoresponse surfaces as a
% function of macular pigment density.  Pretty flat. Some cells are
% peculiarly oriented?

load ('T_cones_smj10');
load ('den_mac_ws');
macpig = den_mac_ws;  % Peaks at 0.5
mult10deg = 0.28; % from SMJ paper (Don't change, used below)
mult2deg = 0.7; % from SMJ paper
macpigtransmittance = 1./(10.^(macpig.*mult10deg));
fund10 = T_cones_smj10'./repmat(max(T_cones_smj10'),81,1);
fundnomacpig = fund10./repmat(macpigtransmittance,1,3);
fundnomacpig = fundnomacpig./repmat(max(fundnomacpig),81,1);

% First gathering all the data
data = [];
[fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\pancolor.txt');
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
    
    M10 = T_cones_smj10*mon_spd;
    Loog = logical(out(:,7));
    ecc = sqrt((NT.sum.exptParams.rf_x/10).^2 + (NT.sum.exptParams.rf_y/10).^2);
    data(cellcounter).M = M;
    data(cellcounter).data = out;
    data(cellcounter).bkgndrgb = NT.sum.exptParams.bkgndrgb;
    data(cellcounter).ecc = ecc;
    data(cellcounter).monspd = mon_spd;
end

% Now, computing cone weights as a function of macular pigment density
macpigmults = linspace(0,mult2deg,10);
for i = 1:length(data)
    i
    cc_original = data(i).data(:,[2:4]).*repmat(data(i).data(:,5),1,3);
    Loog = data(i).data(:,end);
    allconeweights = [];
    for j = 1:length(macpigmults)
        j
        macpigtransmittance = 1./(10.^(macpig.*macpigmults(j)));
        tmpfund = fundnomacpig.*repmat(macpigtransmittance,1,3);
        tmpfund = tmpfund./repmat(max(tmpfund),81,1);
        M = tmpfund'*mon_spd;
        scaled = ConvertConeContrastBasis(data(i).M, M, data(i).bkgndrgb, cc_original);
        [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);
        A = [quadparams(1) quadparams(4) quadparams(5);...
            quadparams(4) quadparams(2) quadparams(6);...
            quadparams(5) quadparams(6) quadparams(3)];
        B = xformmat*A*xformmat';
        [evecs, evals] = eig(B);
        [evals,idxs] = sort(diag(evals));
        evecs = evecs(:,idxs);
        coneweights = evecs(:,1);  % Small eigenvalue = long axis
        if (coneweights(3) < 0)  % Positive S-cone weights
            coneweights = -coneweights;
        end
         allconeweights = [allconeweights; coneweights'];
    end
    data(i).coneweights = allconeweights;
end

allconeweights = cat(3,data(:).coneweights);

% Hist of S-cone weights
figure;
allconeweights = cat(3,data(:).coneweights);
bins = [0:.05:1];
mn = []; stdev = [];
for j = 1:length(macpigmults)
    subplot(ceil(sqrt(length(macpigmults))),ceil(sqrt(length(macpigmults))),j);
    hold on;
    tmpconeweights = squeeze(allconeweights(j,:,:));
    normconeweights = tmpconeweights./repmat(sum(abs(tmpconeweights)),3,1);
    hist(normconeweights(3,:),bins);
    mn(j) = mean(normconeweights(3,:));
    stdev(j) = std(normconeweights(3,:));
    plot(mn(j),0,'m^');
end

% At what macular pigment density do we have zero S-cone weight?
figure; axes; hold on;
errorbar(macpigmults,mn, stdev./sqrt(size(normconeweights,2)));
plot(macpigmults,mn,'k.');
xlabel('macular pigment density multiplier');
ylabel('Mean S-cone weight');

%%
% Section 5
% IsoDETECTION surfaces in cone contrast space with diffent fundamentals
% (different amounts of macular pigment).  
% Adding macular pigment doesn't change things very much. To the extent
% that it does, the long axis becomes more aligned with the S-cone axis
% when more macular pigment was added.

% Preparing fundamentals
load ('T_cones_smj10');
load ('den_mac_ws');
macpig = den_mac_ws;  % Peaks at 0.5
mult10deg = 0.28; % from SMJ paper (Don't change, used below)
mult2deg = 0.7; % from SMJ paper
macpigtransmittance = 1./(10.^(macpig.*mult10deg));
fund10 = T_cones_smj10'./repmat(max(T_cones_smj10'),81,1);
fundnomacpig = fund10./repmat(macpigtransmittance,1,3);
fundnomacpig = fundnomacpig./repmat(max(fundnomacpig),81,1);

% Loading data from DTNT files
[fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\DTEM\SednaDTNT.txt');
collectionofMs = nan(3,3,length(fnames));
for fileidx = 1:size(fnames,1)    
    filename = findfile(char(fnames{fileidx}));
    stro = nex2stro(filename);
    [DTNTthresholds, DTNTcolorDirs, DTNTsf] = DTquestUnpackGH(stro, 'mode');
    DTNTcolorDirs = mkbasis(DTNTcolorDirs');  % Each color is a unit vector in the tested direction
    data(fileidx).colordirs = DTNTcolorDirs;
    data(fileidx).sf = DTNTsf;
    data(fileidx).thresholds = DTNTthresholds;
    data(fileidx).M = reshape(stro.sum.exptParams.m_mtx,3,3);
    data(fileidx).spd = reshape(stro.sum.exptParams.mon_spect,81,3);
    data(fileidx).bkgndrgb = [stro.sum.exptParams.bkgnd_r stro.sum.exptParams.bkgnd_g stro.sum.exptParams.bkgnd_b];
end % Datafile loop

% Was the same calibration file used in all experiments?
if (all(all(min(cat(3,data.M),[],3) ==  max(cat(3,data.M),[],3))))
    disp('All M matrices are the same');
else
    disp('All M matrices are NOT the same');
    disp('The rest of the script will nee editing to deal with this');
    keyboard;
end

% Make sure all the spectr look the same
figure; plot([data.spd]);

% converting cone contrasts to different spaces
sfs = unique([data.sf]);
macpigmults = linspace(0,mult2deg,10);
for sfidx = 1:length(sfs)
   L = [data.sf] == sfs(sfidx);
   unitvects = [data(L).colordirs]';
   thresholds = [data(L).thresholds];
   scaled = unitvects.*repmat(thresholds(:),1,3);
   figure; axes; hold on;
   for j = 1:length(macpigmults)
       macpigtransmittance = 1./(10.^(macpig.*macpigmults(j)));
       tmpfund = fundnomacpig.*repmat(macpigtransmittance,1,3);
       tmpfund = tmpfund./repmat(max(tmpfund),81,1);
       Mnew = tmpfund'*data(1).spd;
       tmpscaled = ConvertConeContrastBasis(data(1).M, Mnew, data(1).bkgndrgb, scaled);
       
       % Trying a pilot plot
       subplot(ceil(sqrt(length(macpigmults))), ceil(sqrt(length(macpigmults))), j);
       hold on;
     %  plot3(scaled(:,1),scaled(:,2),scaled(:,3),'k.');
     %  plot3(-scaled(:,1),-scaled(:,2),-scaled(:,3),'k.');
     %  plot3(tmpscaled(:,1),tmpscaled(:,2),tmpscaled(:,3),'g.');
     %  plot3(-tmpscaled(:,1),-tmpscaled(:,2),-tmpscaled(:,3),'g.');
       

       % Using an ellipsoid fit to start
       scaledsym = [tmpscaled; -tmpscaled];
       D = [scaledsym(:,1) .* scaledsym(:,1),...
           scaledsym(:,2) .* scaledsym(:,2),...
           scaledsym(:,3) .* scaledsym(:,3),...
           2*scaledsym(:,1) .* scaledsym(:,2),...
           2*scaledsym(:,1) .* scaledsym(:,3),...
           2*scaledsym(:,2) .* scaledsym(:,3)];
       lssoln = (D' * D) \(D' * ones(size(scaledsym,1),1));
       A = [lssoln(1) lssoln(4) lssoln(5);...
           lssoln(4) lssoln(2) lssoln(6);...
           lssoln(5) lssoln(6) lssoln(3)];
       [evecs, evals] = eig(A);
       initparams = [2; reshape(evecs*sqrt(evals),9,1)];
       options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-5,'Display','none');
       [fpar,fv] = fminsearch(@(x) colefiterr(x,scaledsym,0),initparams, options)
       
       [x,y,z] = meshgrid(linspace(-6,6,50),linspace(-6,6,50),linspace(-15,15,50));
       tmp = [x(:) y(:) z(:)];
       v = sum(abs(tmp *reshape(fpar(2:end),3,3)).^fpar(1),2);
       fv = isosurface(x,y,z,reshape(v,size(x,1),size(x,2),size(x,3)),1);
       h = patch(fv);
       set(h,'EdgeColor','yellow','FaceColor','yellow');
       if (sfs(sfidx) < .9)
          set(gca,'XLim', [-3 3],'YLim',[-3 3],'ZLim',[-8 8]);
       end
       drawnow
       
       % S-cone component of most S-cone dominated mechanism      
       mechs = mkbasis(reshape(fpar(2:end),3,3))
       title(num2str(max(abs(mechs(3,:)))));
       
       % S-cone component of longest axis (smallest eigenvector)
       [v,d] = eig(A);
       title(num2str(v(3,1)));
   end
end


%%
% Section 6
% Looking at the effects of pigment self-screening on estimated S-cone
% weights. This has no effect on S-cone weights.

load ('T_cones_smj10');
load ('den_mac_ws');
macpig = den_mac_ws;  % Peaks at 0.5
macpigtransmittance = 1./(10.^(macpig.*0.28));
lenstransmittance = 1./(10.^(den_lens_smj(:,2)));

fund10 = T_cones_smj10'./repmat(max(T_cones_smj10'),81,1);
% Undoing macular pigment
fundnomacpig = fund10./repmat(macpigtransmittance,1,3);
% Undoing lens pigment
fundnopig = fundnomacpig./repmat(lenstransmittance,1,3);
fundnopig = fundnopig./repmat(max(fundnopig),81,1);

% First gathering all the data
data = [];
[fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\NTlum2.txt');
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
    
    M10 = T_cones_smj10*mon_spd;
    Loog = logical(out(:,7));
    ecc = sqrt((NT.sum.exptParams.rf_x/10).^2 + (NT.sum.exptParams.rf_y/10).^2);
    data(cellcounter).M = M;
    data(cellcounter).data = out;
    data(cellcounter).bkgndrgb = NT.sum.exptParams.bkgndrgb;
    data(cellcounter).ecc = ecc;
    data(cellcounter).monspd = mon_spd;
end

% Now, computing cone weights as a function of pigment self screening
opticaldensities = linspace(.25,.45,10);
for i = 1:length(data)
    i
    cc_original = data(i).data(:,[2:4]).*repmat(data(i).data(:,5),1,3);
    Loog = data(i).data(:,end);
    allconeweights = [];
    for j = 1:length(opticaldensities)
        j
        actionspectra = -log10(1-fundnopig*(1-10^-opticaldensities(j)));
        % Reapplying the lens and macular pigments
        tmpfund = actionspectra.*repmat(macpigtransmittance,1,3).*repmat(lenstransmittance,1,3);
        tmpfund = tmpfund./repmat(max(tmpfund),81,1);
        M = tmpfund'*mon_spd;
        scaled = ConvertConeContrastBasis(data(i).M, M, data(i).bkgndrgb, cc_original);
        [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);
        coneweights = (xformmat*planeparams)';
        if (sign(coneweights(1)) == -1)
            coneweights = -coneweights;
        end
        allconeweights = [allconeweights; coneweights];
    end
    data(i).coneweights = allconeweights;
end


% DKL plot
figure;
for j = 1:length(opticaldensities)
    subplot(ceil(sqrt(length(opticaldensities))),ceil(sqrt(length(opticaldensities))),j);
    axis square;
    plot([0 1],[1 0],'k-');
    hold on;
    title(num2str(opticaldensities(j)));
    for i = 1:length(data)
        normconeweights = data(i).coneweights(j,:)./sum(abs(data(i).coneweights(j,:)));
        if (normconeweights(3) > 0)
            h = plot(normconeweights(1),normconeweights(2),'ko','MarkerFaceColor','black','MarkerSize',4);
        else
            h = plot(normconeweights(1),normconeweights(2),'ko','MarkerFaceColor','red','MarkerSize',4);            
        end
    end
end


% Hist of S-cone weights
figure;
allconeweights = cat(3,data(:).coneweights);
bins = [-.4:.05:.4];
mn = []; stdev = [];
for j = 1:length(opticaldensities)
    subplot(ceil(sqrt(length(opticaldensities))),ceil(sqrt(length(opticaldensities))),j);
    hold on;
    tmpconeweights = squeeze(allconeweights(j,:,:));
    normconeweights = tmpconeweights./repmat(sum(abs(tmpconeweights)),3,1);
    hist(normconeweights(3,:),bins);
    mn(j) = mean(normconeweights(3,:));
    stdev(j) = std(normconeweights(3,:));
    plot(mn(j),0,'m^');
    [~,p] = ttest(normconeweights(3,:));
    title([num2str(mn(j)),' p = ',num2str(p,1)]);
end


%%
% Section 7) A few views of an example neuron with a planar isoresponse
% surface.
%filename = 'K071309003.nex';  % planar luminance cell looks good
filename = 'K082109009.nex';  % planar luminance cell looks good
filename = 'K082709004.nex'; view = [-4 24]; % planar luminance cell looks good
filename = 'K010110005.nex'; view = [17 26];  % planar luminance cell looks good

load ('T_cones_smj10');
load ('den_mac_ws');
mult2deg = 0.7; % from SMJ paper
macpigmults = linspace(0,mult2deg,10);
macpigmults = [-.2 .15 .7];
macpig = den_mac_ws;  % Peaks at 0.5

% Undoing macular pigment
macpigtransmittance = 1./(10.^(macpig.*0.28));
fundnomacpig = T_cones_smj10'./repmat(macpigtransmittance,1,3);
fundnomacpig = fundnomacpig./repmat(max(fundnomacpig),81,1);

% Getting data ready
stro = nex2stro(findfile(filename));
lmsidxs = [find(strcmp(stro.sum.trialFields(1,:),'lcont'))...
    find(strcmp(stro.sum.trialFields(1,:),'mcont'))...
    find(strcmp(stro.sum.trialFields(1,:),'scont'))];
% -------------------------------
% Converting cone contrasts in nex file to 10 deg fundmentals.
% Must go through excitations first - can't just transform contrasts.
% -------------------------------
fundamentals = stro.sum.exptParams.fundamentals;
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = stro.sum.exptParams.mon_spd;
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
Moriginal = fundamentals'*mon_spd;
bkgndrgb = stro.sum.exptParams.bkgndrgb;
out = NTpreprocess(stro,0,Inf);
scaled = out(:,[2 3 4]).*repmat(out(:,5),1,3);
Loog = logical(out(:,7));
wls = [380:5:780];
XLIM=.25; YLIM=.25; ZLIM=.5;

figure;
for i = 1:length(macpigmults)     
    % making new fundamentals
    macpigtransmittance = 1./(10.^(macpig.*macpigmults(i)));
    tmpfunds = fundnomacpig.*repmat(macpigtransmittance,1,3);
    tmpfunds = tmpfunds./repmat(max(tmpfunds),81,1);
    Mnew = tmpfunds'*mon_spd;
    tmpscaled = ConvertConeContrastBasis(Moriginal, Mnew, bkgndrgb, scaled);
    % Moving to a DKL-ish space
    tmpscaled(:,[1 2]) = tmpscaled(:,[1 2])*mkbasis([1 1; 1 -1]);
    Loog = logical(out(:,end));
    [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(tmpscaled, Loog);
    coneweights = (xformmat*planeparams)'
    
    subplot(2,2,i);
    title(['Macpig density at 460 nm: ',num2str(macpig(wls==460).*macpigmults(i))]);
    hold on;
    axis square; 
    xlabel('L+M'); set(get(gca,'XLabel'),'Color',[0 0 0]);
    ylabel('L-M'); set(get(gca,'YLabel'),'Color',[0 0 0]);
    zlabel('S'); set(get(gca,'ZLabel'),'Color',[0 0 0]);
    h=plot3(tmpscaled(~Loog,1),tmpscaled(~Loog,2),tmpscaled(~Loog,3),'ko');
    set(h,'Markersize',3,'Markerfacecolor','black');
    h=plot3(-tmpscaled(~Loog,1),-tmpscaled(~Loog,2),-tmpscaled(~Loog,3),'ko');
    set(h,'Markersize',3,'Markerfacecolor','black');
    %     % Plotting OOG rays (squeezing them inside the axis limits)
    plotlims = [XLIM YLIM ZLIM];
    normfacts = max(abs(tmpscaled(Loog,:)./repmat(plotlims,sum(Loog),1)),[],2);
    normfacts = max(normfacts,1);
    h=plot3([zeros(sum(Loog),1) tmpscaled(Loog,1)./normfacts]',[zeros(sum(Loog),1) tmpscaled(Loog,2)./normfacts]',[zeros(sum(Loog),1) tmpscaled(Loog,3)./normfacts]','-');
    set(h,'Color',[.8 .8 .8]);
    h=plot3([zeros(sum(Loog),1) -tmpscaled(Loog,1)./normfacts]',[zeros(sum(Loog),1) -tmpscaled(Loog,2)./normfacts]',[zeros(sum(Loog),1) -tmpscaled(Loog,3)./normfacts]','-');
    set(h,'Color',[.8 .8 .8]);
    
    % Adjusting the axes
    set(gca,'XLim',XLIM*[-1 1],'Ylim',YLIM*[-1 1],'Zlim',ZLIM*[-1 1])
    set(gca,'View',view);
    for whichaxes = {'XTick','YTick','ZTick'}
        tmp = get(gca,char(whichaxes)); set(gca,char(whichaxes),[tmp(1) 0 tmp(end)]);
    end
    xticks = get(gca,'XTick');
    set(gca,'FontSize',8)
    %set(gca,'XTick',[],'YTick',[],'ZTick',[])
    
    % Plotting plane fits
    DKLweights = planeparams'*xformmat';
    [x,y,z] = meshgrid(linspace(-plotlims(1),plotlims(1),50),linspace(-plotlims(2),plotlims(2),50),linspace(-plotlims(3),plotlims(3),50));
    v = abs(x.*DKLweights(1)+y.*DKLweights(2)+z.*DKLweights(3));
    fv = isosurface(x,y,z,v,1);
    
    % added by zack
    vertCounts = hist(fv.faces(:), length(fv.vertices)); % vertex count
    edgeVerts = [find(vertCounts==1) find(vertCounts==2) find(vertCounts==3)]; % get vertices that were only used 1-3 times
    
    % separate the vertices by plane
    p1 = []; p2 = [];
    for j = 1:length(edgeVerts)
        if DKLweights*fv.vertices(edgeVerts(j),:)' > 0 % the LHS will be ±1
            p1 = [p1; edgeVerts(j)];
        else
            p2 = [p2; edgeVerts(j)];
        end
    end
    % get the center point of both planes
    center1 = mean(fv.vertices(p1,:));
    center2 = mean(fv.vertices(p2,:));
    
    % thetas for each edge point about the center (per plane)
    theta1 = atan2(fv.vertices(p1,2)-center1(2),fv.vertices(p1,1)-center1(1));
    theta2 = atan2(fv.vertices(p2,2)-center2(2),fv.vertices(p2,1)-center2(1));
    
    % sort 'em
    [~,sortorder1] = sort(theta1);
    [~,sortorder2] = sort(theta2);
    
    ppoints1 = fv.vertices(p1(sortorder1),:);
    ppoints2 = fv.vertices(p2(sortorder2),:);
    
    % plot 'em
    h = patch(ppoints1(:,1),ppoints1(:,2),ppoints1(:,3),[0 .75 0]);
    set(h,'FaceAlpha',0.5,'EdgeAlpha',0,'EdgeColor','none','LineWidth',1);
    h = patch(ppoints2(:,1),ppoints2(:,2),ppoints2(:,3),[0 .75 0]);
    set(h,'FaceAlpha',0.5,'EdgeAlpha',0,'EdgeColor','none','LineWidth',1);
end
text(1,1,0,filename);
set(gcf,'Renderer','painters');  
% Import to Illustrator and make the green parts transparent

%%
% Section 7.1) A few views of an example neuron with a planar isoresponse
% surface.  2-D plots.
filename = 'K071309003.nex';  % planar luminance cell looks good
filename = 'K082109009.nex';  % planar luminance cell looks good
%filename = 'K082709004.nex';  % planar luminance cell looks good
%filename = 'K010110005.nex'; % planar luminance cell looks good
XLIM = .15; YLIM = .7;

load ('T_cones_smj10');
load ('den_mac_ws');
mult2deg = 0.7; % from SMJ paper
macpigmults = linspace(0,mult2deg,10);
macpigmults = [-.2 .13 .7];
macpig = den_mac_ws;  % Peaks at 0.5

% Undoing macular pigment
macpigtransmittance = 1./(10.^(macpig.*0.28));
fundnomacpig = T_cones_smj10'./repmat(macpigtransmittance,1,3);
fundnomacpig = fundnomacpig./repmat(max(fundnomacpig),81,1);

% Getting data ready
stro = nex2stro(findfile(filename));
lmsidxs = [find(strcmp(stro.sum.trialFields(1,:),'lcont'))...
    find(strcmp(stro.sum.trialFields(1,:),'mcont'))...
    find(strcmp(stro.sum.trialFields(1,:),'scont'))];
% -------------------------------
% Converting cone contrasts in nex file to 10 deg fundmentals.
% Must go through excitations first - can't just transform contrasts.
% -------------------------------
fundamentals = stro.sum.exptParams.fundamentals;
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = stro.sum.exptParams.mon_spd;
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
Moriginal = fundamentals'*mon_spd;
bkgndrgb = stro.sum.exptParams.bkgndrgb;
out = NTpreprocess(stro,0,Inf);
scaled = out(:,[2 3 4]).*repmat(out(:,5),1,3);
Loog = logical(out(:,7));
wls = [380:5:780];

figure;
for i = 1:length(macpigmults)     
    % making new fundamentals
    macpigtransmittance = 1./(10.^(macpig.*macpigmults(i)));
    tmpfunds = fundnomacpig.*repmat(macpigtransmittance,1,3);
    tmpfunds = tmpfunds./repmat(max(tmpfunds),81,1);
    Mnew = tmpfunds'*mon_spd;
    tmpscaled = ConvertConeContrastBasis(Moriginal, Mnew, bkgndrgb, scaled);
%    % Moving to a DKL-ish space
%    tmpscaled(:,[1 2]) = tmpscaled(:,[1 2])*mkbasis([1 1; 1 -1]);
    Loog = logical(out(:,end));
    [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(tmpscaled, Loog);
    coneweights = (xformmat*planeparams)'
    if (sign(coneweights(1)) == -1)
        coneweights = -coneweights;
    end
    coneweights./sum(abs(coneweights))
    % Destructively changing tmpscaled here
    lmcoefs = mkbasis(coneweights([1 2]));
    xyscaled = tmpscaled(:,[1 2])*lmcoefs';
    xyscaled(:,2) = tmpscaled(:,3);
    
    subplot(2,2,i);
    title(['Macpig density at 460 nm: ',num2str(macpig(wls==460).*macpigmults(i))]);
    hold on;
    axis square;

    % Plotting the data points
    h = plot(xyscaled(~Loog,1),xyscaled(~Loog,2),'ko');
    set(h,'Markersize',3,'Markerfacecolor','black');
    h = plot(-xyscaled(~Loog,1),-xyscaled(~Loog,2),'ko');
    set(h,'Markersize',3,'Markerfacecolor','black');
    xlabel([num2str(lmcoefs(1)./sum(lmcoefs),2),'L+',num2str(lmcoefs(2)/sum(lmcoefs),2),'M']);
    set(get(gca,'XLabel'),'Color',[0 0 0]);
    ylabel('S'); set(get(gca,'YLabel'),'Color',[0 0 0]);

    % Plotting the OOG directions
    normfact = 1;
    h = plot([zeros(sum(Loog),1) xyscaled(Loog,1)./normfact]',[zeros(sum(Loog),1) xyscaled(Loog,2)./normfact]','-');
    set(h,'Color',[.8 .8 .8]);
    h = plot([zeros(sum(Loog),1) -xyscaled(Loog,1)./normfact]',[zeros(sum(Loog),1) -xyscaled(Loog,2)./normfact]','-');
    set(h,'Color',[.8 .8 .8]);
    
    % Adjusting the axes
    set(gca,'XLim',XLIM*[-1 1],'Ylim',YLIM*[-1 1]);
    
    % Plotting plane fits (which are actually lines in this projection)
    x = mean(abs(xyscaled))+[-.05 .05];
    y = (-norm(coneweights([1 2]))/coneweights(3))*x+(1/coneweights(3));
    h = plot(x,y,'g-','LineWidth',3);
    h = plot(-x,-y,'g-','LineWidth',3);
    
end
text(1,1,0,filename);
set(gcf,'Renderer','painters');  

%%
% Section 8
% Cone weights of L+M neurons as a function of added macular pigment
% density.  Based largely on section 2.
% Removing macular pigment from the SMJ 10 deg cones using .28 * the WS
% macular pigment template (which peaks at 0.5 at 460 nm).

load ('T_cones_smj10');
load ('den_mac_ws');
wls = [380:5:780];
macpig = den_mac_ws;  % Peaks at 0.5
mult10deg = 0.28; % from SMJ paper (Don't change, used below)
mult2deg = 0.7; % from SMJ paper
macpigtransmittance = 1./(10.^(macpig.*mult10deg));
fundnomacpig = T_cones_smj10'./repmat(macpigtransmittance,1,3);
fundnomacpig = fundnomacpig./repmat(max(fundnomacpig),81,1);

% First gathering all the data
data = [];
[fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\NTlum2.txt');
for cellcounter = 1:size(fnames,1)
    filename = findfile(char(fnames{cellcounter}));
    paradigmID = getparadigmID(filename);
    if paradigmID == 103
        NT = nex2stro(filename);
    end
    out = NTpreprocess(NT,0,Inf);
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

% Now, computing cone weights as a function of macular pigment density
macpigmults = linspace(0,mult2deg,10);
for i = 1:length(data)
    i
    cc_original = data(i).data(:,[2:4]).*repmat(data(i).data(:,5),1,3);
    Loog = data(i).data(:,end);
    allconeweights = [];
    for j = 1:length(macpigmults)
        j
        macpigtransmittance = 1./(10.^(macpig.*macpigmults(j)));
        tmpfund = fundnomacpig.*repmat(macpigtransmittance,1,3);
        tmpfund = tmpfund./repmat(max(tmpfund),81,1);
        M = tmpfund'*mon_spd;
        scaled = ConvertConeContrastBasis(data(i).M, M, data(i).bkgndrgb, cc_original);
        [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);
        coneweights = (xformmat*planeparams)';
        if (sign(coneweights(1)) == -1)
            coneweights = -coneweights;
        end
        allconeweights = [allconeweights; coneweights];
    end
    data(i).coneweights = allconeweights;
end

% Computing means and standard deviations
allconeweights = cat(3,data(:).coneweights);
mn = []; stdev = [];
for j = 1:length(macpigmults)
    tmpconeweights = squeeze(allconeweights(j,:,:));
    normconeweights = tmpconeweights./repmat(sum(abs(tmpconeweights)),3,1);
    mn(j) = mean(normconeweights(3,:));
    stdev(j) = std(normconeweights(3,:));
end

% Doing the plotting
% At what macular pigment density do we have zero S-cone weight?
figure; axes; hold on;
h = errorbar(macpigmults*macpig(wls==460),mn, stdev./sqrt(size(normconeweights,2)));
set(h,'LineWidth',2,'Color','black');
h = plot(macpigmults*macpig(wls==460),mn,'ko','LineWidth',1,'MarkerFaceColor','black');
xlabel('macular pigment density at 460 nm');
ylabel('Mean S-cone weight');
b = regress(mn', [macpigmults'*macpig(wls==460) ones(length(mn),1)]);

plot(mult2deg*macpig(wls==460),b(1)*mult2deg*macpig(wls==460)+b(2),'r*');
plot(mult10deg*macpig(wls==460),b(1)*mult10deg*macpig(wls==460)+b(2),'g*');
plot(-b(2)/b(1),0,'y*');
set(gca,'XLim',[-.05 0.4],'YLim',[-.03 0.06])

% Hist of S-cone weights
figure;
allconeweights = cat(3,data(:).coneweights);
bins = [-.1:.01:.1];
maxcount = 0;
for j = 1:length(macpigmults)
    subplot(ceil(sqrt(length(macpigmults))),ceil(sqrt(length(macpigmults))),j);
    hold on;
    tmpconeweights = squeeze(allconeweights(j,:,:));
    normconeweights = tmpconeweights./repmat(sum(abs(tmpconeweights)),3,1);
    [n,x] = hist(normconeweights(3,:),bins);
    maxcount = max(maxcount, max(n));
    bar(x,n,'FaceColor','black');
    [~,p] = ttest(normconeweights(3,:));
    title([' p = ',num2str(p,1)]);
end
for j = 1:length(macpigmults)
    subplot(ceil(sqrt(length(macpigmults))),ceil(sqrt(length(macpigmults))),j);
    set(gca,'Xlim',[bins(1) bins(end)],'Ylim',[0 maxcount+1]);
    plot(mn(j),maxcount,'kv','MarkerFaceColor','magenta');
end


%%
% Section 9
% The long axis orientation of ellipsoidal isoresponse surfaces as a
% function of macular pigment density.  POSTER READY. IN PROGRESS.

load ('T_cones_smj10');
load ('den_mac_ws');
macpig = den_mac_ws;  % Peaks at 0.5
mult10deg = 0.28; % from SMJ paper (Don't change, used below)
mult2deg = 0.7; % from SMJ paper
macpigtransmittance = 1./(10.^(macpig.*mult10deg));
fund10 = T_cones_smj10'./repmat(max(T_cones_smj10'),81,1);
fundnomacpig = fund10./repmat(macpigtransmittance,1,3);
fundnomacpig = fundnomacpig./repmat(max(fundnomacpig),81,1);
wls = [380:5:780];

% First gathering all the data
data = [];
[fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\pancolor.txt');
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
    
    M10 = T_cones_smj10*mon_spd;
    Loog = logical(out(:,7));
    ecc = sqrt((NT.sum.exptParams.rf_x/10).^2 + (NT.sum.exptParams.rf_y/10).^2);
    data(cellcounter).M = M;
    data(cellcounter).data = out;
    data(cellcounter).bkgndrgb = NT.sum.exptParams.bkgndrgb;
    data(cellcounter).ecc = ecc;
    data(cellcounter).monspd = mon_spd;
end

% Now, computing cone weights as a function of macular pigment density
macpigmults = linspace(0,mult2deg,10);
for i = 1:length(data)
    i
    cc_original = data(i).data(:,[2:4]).*repmat(data(i).data(:,5),1,3);
    Loog = data(i).data(:,end);
    allconeweights = [];
    for j = 1:length(macpigmults)
        j
        macpigtransmittance = 1./(10.^(macpig.*macpigmults(j)));
        tmpfund = fundnomacpig.*repmat(macpigtransmittance,1,3);
        tmpfund = tmpfund./repmat(max(tmpfund),81,1);
        M = tmpfund'*mon_spd;
        scaled = ConvertConeContrastBasis(data(i).M, M, data(i).bkgndrgb, cc_original);
        [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);
        A = [quadparams(1) quadparams(4) quadparams(5);...
            quadparams(4) quadparams(2) quadparams(6);...
            quadparams(5) quadparams(6) quadparams(3)];
        B = xformmat*A*xformmat';
        [evecs, evals] = eig(B);
        [evals,idxs] = sort(diag(evals));
        evecs = evecs(:,idxs);
        coneweights = evecs(:,1);  % Small eigenvalue = long axis
        if (coneweights(3) < 0)  % Positive S-cone weights
            coneweights = -coneweights;
        end
         allconeweights = [allconeweights; coneweights'];
    end
    data(i).coneweights = allconeweights;
end

% Computing means and standard deviations
allconeweights = cat(3,data(:).coneweights);
mn = []; stdev = [];
for j = 1:length(macpigmults)
    tmpconeweights = squeeze(allconeweights(j,:,:));
    normconeweights = tmpconeweights./repmat(sum(abs(tmpconeweights)),3,1);
    mn(j) = mean(normconeweights(3,:));
    stdev(j) = std(normconeweights(3,:));
end

% Doing the plotting
figure; axes; hold on;
h = errorbar(macpigmults*macpig(wls==460),mn, stdev./sqrt(size(normconeweights,2)));
set(h,'LineWidth',2,'Color','black');
h = plot(macpigmults*macpig(wls==460),mn,'ko','LineWidth',1,'MarkerFaceColor','black');
xlabel('macular pigment density at 460 nm');
ylabel('Mean S-cone weight');

plot(mult2deg*macpig(wls==460),.78,'r*');
plot(mult10deg*macpig(wls==460),.78,'g*');
set(gca,'XLim',[-.05 .4],'Ylim',[.78 .88]);


%%
% Section 10
% A few views of ellipsoids calculated with different sets of fundamentals.
PLOTSURF = 1;
PLOTDATA = 1;
PLOTAXIS = 1;
filename = 'S073010002.nex';  % Pancolor Good example
filename = 'S041310005.nex';  % ill conditioned B!!
%filename = 'S031610003.nex';  % Good one.
filename = 'S041510003.nex'

load ('T_cones_smj10');
load ('den_mac_ws');
mult2deg = 0.7; % from SMJ paper
macpigmults = linspace(0,mult2deg,10);
macpigmults = [0 .15 .7];
macpigmults = [-1 .15 3];

macpig = den_mac_ws;  % Peaks at 0.5

% Undoing macular pigment
macpigtransmittance = 1./(10.^(macpig.*0.28));
fundnomacpig = T_cones_smj10'./repmat(macpigtransmittance,1,3);
fundnomacpig = fundnomacpig./repmat(max(fundnomacpig),81,1);

% Getting data ready
stro = nex2stro(findfile(filename));
lmsidxs = [find(strcmp(stro.sum.trialFields(1,:),'lcont'))...
    find(strcmp(stro.sum.trialFields(1,:),'mcont'))...
    find(strcmp(stro.sum.trialFields(1,:),'scont'))];
% -------------------------------
% Converting cone contrasts in nex file to 10 deg fundmentals.
% Must go through excitations first - can't just transform contrasts.
% -------------------------------
fundamentals = stro.sum.exptParams.fundamentals;
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = stro.sum.exptParams.mon_spd;
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
Moriginal = fundamentals'*mon_spd;
bkgndrgb = stro.sum.exptParams.bkgndrgb;
out = NTpreprocess(stro,0,Inf);
scaled = out(:,[2 3 4]).*repmat(out(:,5),1,3);
Loog = logical(out(:,7));
wls = [380:5:780];

figure;
for i = 1:length(macpigmults)     
    % making new fundamentals
    macpigtransmittance = 1./(10.^(macpig.*macpigmults(i)));
    tmpfunds = fundnomacpig.*repmat(macpigtransmittance,1,3);
    tmpfunds = tmpfunds./repmat(max(tmpfunds),81,1);
    Mnew = tmpfunds'*mon_spd;
    tmpscaled = ConvertConeContrastBasis(Moriginal, Mnew, bkgndrgb, scaled);
    Loog = logical(out(:,end));
    [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(tmpscaled, Loog);

    A = [quadparams(1) quadparams(4) quadparams(5);...
        quadparams(4) quadparams(2) quadparams(6);...
        quadparams(5) quadparams(6) quadparams(3)];
    B = xformmat*A*xformmat';
    cond(B)
    [evecs, evals] = eig(B);
    [evals,idxs] = sort(diag(evals));
    evecs = evecs(:,idxs);
    coneweights = evecs(:,1)  % Small eigenvalue = long axis
    if (coneweights(1) < 0)  % Positive L-cone weights
        coneweights = -coneweights;
    end

    % Setting x axis to the linear combination of l,m coneweights
    % that makes the maximum projection of long axis visible
    % y is orthogonal to x.
    lmcoefs = mkbasis(coneweights([1 2]));
    xmat = [lmcoefs' 0; [lmcoefs(2) -lmcoefs(1) 0]; 0 0 1];
    xyscaled = tmpscaled*xmat';
    B = xmat*B*xmat';
    coneweights = xmat*coneweights;
    
    subplot(2,2,i);
    title(['Macpig density at 460 nm: ',num2str(macpig(wls==460).*macpigmults(i))]);
    hold on;
    axis square;
    set(gca,'View',[0,0]);

    % Plotting the data points
    if (PLOTDATA)
        h = plot3(xyscaled(~Loog,1),xyscaled(~Loog,2),xyscaled(~Loog,3),'ko');
        set(h,'Markersize',3,'Markerfacecolor','black');
        h = plot3(-xyscaled(~Loog,1),-xyscaled(~Loog,2),-xyscaled(~Loog,3),'ko');
        set(h,'Markersize',3,'Markerfacecolor','black');
        xlabel([num2str(lmcoefs(1)./sqrt(norm(lmcoefs)),2),'L+',num2str(lmcoefs(2)/sqrt(norm(lmcoefs)),2),'M']);
        set(get(gca,'XLabel'),'Color',[0 0 0]);
        ylabel('S'); set(get(gca,'YLabel'),'Color',[0 0 0]);
    end
    
    % Adjusting the axes
  %  LIM = ceil(max(abs(tmpscaled(:,3)))*10)/10;
    LIM = max(abs(xyscaled(:,3)));
    set(gca,'XLim',LIM*[-1 1],'Ylim',LIM*[-1 1])
    
    % Plot surfaces
    if (PLOTSURF)
        [xx yy zz] = meshgrid(linspace(-LIM,LIM,100),...
            linspace(-LIM,LIM,100),...
            linspace(-LIM,LIM,100));
        xformedxyz = [xx(:) yy(:) zz(:)];
        
        variables = [xformedxyz(:,1).^2 xformedxyz(:,2).^2 xformedxyz(:,3).^2 2*xformedxyz(:,1).*xformedxyz(:,2) 2*xformedxyz(:,1).*xformedxyz(:,3) 2*xformedxyz(:,2).*xformedxyz(:,3)];
        coefficients = [B(1,1) B(2,2) B(3,3) B(1,2) B(1,3) B(2,3)]';
        fr = variables*coefficients;
        fv = isosurface(xx,yy,zz,reshape(fr,size(xx)), 1);
        h = patch(fv);
        set(h,'FaceAlpha',0.5,'EdgeAlpha',0.5);
        set(h,'FaceVertexCData',repmat([0 .75 0],size(fv.vertices,1),1))
        set(h,'CDataMapping','direct');
        set(h,'FaceColor','interp','EdgeColor','none');
        lighting gouraud
        material metal
        camlight headlight
    end

    if (PLOTAXIS)
        % Plotting the major axis
        plot3(coneweights(1)*LIM*[-1 1],coneweights(2)*LIM*[-1 1],coneweights(3)*LIM*[-1 1],'k-','LineWidth',4);
        text(LIM/2,0,num2str(coneweights(3)));
    end
end
text(1,1,0,filename);
if (PLOTSURF)
% Import to Illustrator and make the green parts transparent
    set(gcf,'Renderer','openGL');
    set(gcf,'InvertHardCopy','off')
%    print -dtiff -r600 -cmyk junk1
else
    set(gcf,'Renderer','painters');  
end

%%
% Section 11
% Scatterplot of "S-cone weight" calculated using the 2 degree fundamentals
% and retinal eccentricity.

% First gathering all the data
load ('T_cones_smj');
load ('T_cones_smj10');
load ('den_mac_ws');
data = [];
[fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\NTlum2.txt');
for cellcounter = 1:size(fnames,1)
    filename = findfile(char(fnames{cellcounter}));
    paradigmID = getparadigmID(filename);
    if paradigmID == 103
        NT = nex2stro(filename);
    end
    out = NTpreprocess(NT,0,Inf);
    fundamentals = NT.sum.exptParams.fundamentals;
    fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
    mon_spd = NT.sum.exptParams.mon_spd;
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
    M = fundamentals'*mon_spd;
    M2 = T_cones_smj*mon_spd;
    scaled = ConvertConeContrastBasis(M, M2, NT.sum.exptParams.bkgndrgb, out(:,[2:4]).*repmat(out(:,5),1,3));
    Loog = logical(out(:,7));
    ecc = sqrt((NT.sum.exptParams.rf_x/10).^2 + (NT.sum.exptParams.rf_y/10).^2);
    data(cellcounter).M = M;
    data(cellcounter).scaled = scaled;
    data(cellcounter).Loog = Loog;
    data(cellcounter).bkgndrgb = NT.sum.exptParams.bkgndrgb;
    data(cellcounter).ecc = ecc;
    data(cellcounter).monspd = mon_spd;
    fn = char(fnames{cellcounter});
    data(cellcounter).monk = fn(1);
end

for i = 1:length(data)
    i
    [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(data(i).scaled, data(i).Loog);
    coneweights = (xformmat*planeparams)';
    if (sign(coneweights(1)) == -1)
        coneweights = -coneweights;
    end
    data(i).deg2coneweights = coneweights';
end

iskali = logical([data.monk]=='K');
ecc = [data.ecc];
coneweights = [data.deg2coneweights];
normconeweights = coneweights./repmat(sum(abs(coneweights)),3,1);
figure; axes; hold on;
plot(ecc(iskali),normconeweights(3,iskali),'ko','MarkerFaceColor','black','MarkerSize',7);
plot(ecc(~iskali),normconeweights(3,~iskali),'ro','MarkerFaceColor','red','MarkerSize',7);
[b,bint,r,rint,stats] = regress(normconeweights(3,:)', [ecc' ones(length(ecc),1)]);
plot([min(ecc) max(ecc)],b(1)*[min(ecc) max(ecc)]+b(2),'k:','Linewidth',2);
xlabel('Eccentricity (deg)');
ylabel({'S-cone weight','Calculated with the SMJ 2 degree cone fundamentals'});
[r,p] = corr([ecc',normconeweights(3,:)'],'Type','Spearman')
set(gca,'XTick',2:7);
title(['Spearman''s r = ',num2str(r(1,2),3),' p = ',num2str(p(1,2),3)]);

% Now I need to come up with some kind of reaonable prediction.
% I'm going to start with the assumption of 0 S-cone pigment at 10 degrees
% and some reasonable level at 2 degrees.
fund2 = T_cones_smj'./repmat(max(T_cones_smj'),81,1);
fund10 = T_cones_smj10'./repmat(max(T_cones_smj10'),81,1);
macpigmult10deg = 0.28; % from SMJ paper (Don't change, used below)
macpigtransmittance = 1./(10.^(den_mac_ws.*macpigmult10deg));
fundnomacpig = fund10./repmat(macpigtransmittance,1,3);
fundnomacpig = fundnomacpig./repmat(max(fundnomacpig),81,1);

% let's pretend the neurons are tuned for L+M exactly
% Finding preferred color direction in gun space
Mnomacpig = fundnomacpig'*mon_spd;
M10 = fund10'*mon_spd;
M2 = fund2'*mon_spd;
% Imagine a cell that has cone weights [1 1 0]
% Leftmost "M" is the gun weights assuming some degree of macular pigmentation
% Rightmost "M"  converts these gun weights back to cone weights, but 
% assumes a different degree of macular pigmentation
newconeweights = inv(M2')*Mnomacpig'*[1 1 0]'
newconeweights = newconeweights./sum(abs(newconeweights));
slope = newconeweights(3)/10;
h(1) = plot([min(ecc) max(ecc)],slope*[min(ecc) max(ecc)],'m:','Linewidth',2);
% Line of predicted S-cone weights assuming an L+M cell and zero macular
% pigment at 10 degrees and the two degree fundamentals (a multiplier of
% 0.7) at 0 degree eccentricity.


%%
% Section 12
% What are the predicted blue to green luminance ratios for different macular 
% pigment and lens densities?  For comparison with section 3.
stro = nex2stro('N:\NexFiles\Greg\Sedna\2010\S120310003.nex'); % just for the mon spds
mon_spd = stro.sum.exptParams.mon_spd;
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
mult2deg = 0.7; % from SMJ paper
LVLAMBDA = 1;
load ('T_cones_smj10');
load ('den_mac_ws');
%den_lens_vnv = [2.61 2.31 2.02 1.73 1.45 1.19 0.96 0.76 0.59 0.47 0.37 0.31 0.27 0.24...
%    0.225 0.205 0.195 0.185 0.175 0.165 0.155 0.145 0.14 0.13 0.125 0.115 0.11 0.105...
%    0.1 0.095 0.09 0.085 0.08 0.075 0.07 0.065 0.06 0.055 0.05 0.045 0.04 0.035...
%    0.031 0.027 0.024 0.021 0.018 0.015 0.012 0.01 0.008 0.006 0.004 0.003 0.002 0.001 zeros(1,25)]';
% downloaded from cvrl.org 
den_lens_smj = [3.2 2.8 2.4 1.985 1.62 1.308 1.05 0.84 0.675 0.557 0.468 0.393 0.335 0.29 0.26 0.24 0.225 0.215 0.203 0.191 0.18 0.168 0.162 0.151 0.145 0.133...
        0.128 0.122 0.116 0.11 0.104 0.099 0.093 0.087 0.081 0.075 0.07 0.064 0.058 0.052 0.046 0.041 0.036 0.031 0.028 0.024 0.021 0.017 0.014 0.012 0.009 0.007...
        0.005 0.003 0.002 0.001 zeros(1,25)]';
% First two numbers of den_lens_smj are extrapolations.

lens = den_lens_smj./1.16; % Dividing by 1.16 to make it open pupil 
macpig = den_mac_ws;  % Peaks at 0.5
macpigmult10deg = 0.28; % from SMJ paper (Don't change, used below)
macpigmult2deg = 0.7; % from SMJ paper
lensmult = 1.28; % for SMJ paper
macpigtransmittance = 1./(10.^(macpig.*macpigmult10deg));
lenstransmittance = 1./(10.^(lens.*lensmult));
fund10 = T_cones_smj10'./repmat(max(T_cones_smj10'),81,1);
fundnopig = fund10./repmat(macpigtransmittance,1,3)./repmat(lenstransmittance,1,3);
fundnopig = fundnopig./repmat(max(fundnopig),81,1);

macpigmults = linspace(0,macpigmult2deg,8);
lensmults = linspace(.3,1.5,8);  % 1.28 is SMJ value
data = [];
for j = 1:length(macpigmults)
    for k = 1:length(lensmults)
        macpigtransmittance = 1./(10.^(macpig.*macpigmults(j)));
        lenstransmittance = 1./(10.^(lens.*lensmults(k)));
        tmpfunds = fundnopig.*repmat(macpigtransmittance,1,3).*repmat(lenstransmittance,1,3);
        tmpfunds = tmpfunds./repmat(max(tmpfunds),81,1);
        M = tmpfunds'*mon_spd;
        gpred = [LVLAMBDA 1]*M([1 2],2);
        bpred = [LVLAMBDA 1]*M([1 2],3);
        data(j,k) = gpred/bpred;
    end
end
figure;
contour(abs(data-3.6),[0 .1 .2 .5 1]);
set(gca,'XTicklabel',num2str(macpigmults'),'YTicklabel',num2str(lensmults'));
xlabel('mac pig mult'); ylabel('lens mult');

%%
% Section 13
% Estimating macular pigment density from DTspot - blum vs glum as a
% function of retina eccentricity.

for i = 1:3
    if (i == 1)
        filelistname = 'C:\NO BACKUP\NexFiles\nexfilelists\Greg\DTEM\GregMacPig.txt'
    elseif (i == 2)
        filelistname = 'C:\NO BACKUP\NexFiles\nexfilelists\Greg\DTEM\KaliMacPig.txt'
    else
        filelistname = 'C:\NO BACKUP\NexFiles\nexfilelists\Greg\DTEM\SednaMacPig.txt'
    end
    
    filenames = fnamesFromTxt2(filelistname);
    data = [];
    for fileidx = 1:size(filenames,1)
        stro = nex2stro(findfile(char(filenames{fileidx,:})));
        [thresholds, colorDirs, sfs] = DTquestUnpack(stro, 'weibull');
        close;
        ecc = stro.sum.exptParams.rf_x;
        data = [data; ecc/10 thresholds'];
    end
    
    figure('Position',[305    79   329   627]);
    subplot(2,1,1); hold on;
    plot(data(:,1),data(:,2),'go','MarkerFaceColor','green','MarkerSize',5);
    plot(data(:,1),data(:,3),'bo','MarkerFaceColor','blue','MarkerSize',5);
    ylabel({'Detection threshold','(arb. units)'});
    regcoef_g = regress(data(:,2),[data(:,1), ones(size(data,1),1)]);
    regcoef_b = regress(data(:,3),[data(:,1), ones(size(data,1),1)]);
    plot([0 8],[0 8]*regcoef_g(1)+regcoef_g(2),'g-','Linewidth',2);
    plot([0 8],[0 8]*regcoef_b(1)+regcoef_b(2),'b-','Linewidth',2);
    set(gca,'XLim',[0 8],'YLim',[0 60]);
    
    subplot(2,1,2); hold on;
    plot(data(:,1),data(:,3)./data(:,2),'ko','MarkerFaceColor','black','MarkerSize',5);
    [r,p] = corrcoef ([data(:,1),data(:,3)./data(:,2)])
    ylabel('Bthresh/Gthresh'); xlabel('Eccentricity (deg)');
    plot(linspace(0,8,100),(linspace(0,8,100)*regcoef_b(1)+regcoef_b(2))./(linspace(0,8,100)*regcoef_g(1)+regcoef_g(2)),'k-','Linewidth',2);
    title(['r = ',num2str(r(1,2),3),' p = ',num2str(p(1,2),3)]);
    
    % Predictions on blum/glum
    LVLAMBDA = 1; % L contribution to vlambda
    
    % Computing predictions from vlambda
    % Just using the default fundamentals for now.
    % Assuming lum is L+M
    %M = reshape(stro.sum.exptParams.m_mtx,3,3);
    load('T_cones_smj')
    spds = reshape(stro.sum.exptParams.mon_spect,81,3);
    M = T_cones_smj*spds;
    gpred = [LVLAMBDA 1]*M([1 2],2);
    bpred = [LVLAMBDA 1]*M([1 2],3);
    subplot(2,1,2);
    plot([0 8], gpred/bpred*[1 1],'k:');
    text(max(data(:,1)),gpred/bpred, '2 deg');
    
    load('T_cones_smj10')
    spds = reshape(stro.sum.exptParams.mon_spect,81,3)
    M = T_cones_smj10*spds;
    gpred = [LVLAMBDA 1]*M([1 2],2);
    bpred = [LVLAMBDA 1]*M([1 2],3);
    subplot(2,1,2);
    plot([0 8], gpred/bpred*[1 1],'k:');
    text(max(data(:,1)),gpred/bpred, '10 deg','FontSize',9)
    
    % Undoing preretinal filters (a la SJM 1993)
    load ('den_mac_ws');
    macpig = den_mac_ws;
    macpigtransmittance = 1./(10.^(macpig*.28));
    funds = T_cones_smj10./repmat(max(T_cones_smj10,[],2),1,size(T_cones_smj10,2));
    funds = funds./repmat(macpigtransmittance,1,3)';  % removing macular pigment
    funds = funds./repmat(max(funds')',1,81);
    
    M = funds*spds;
    gpred = [LVLAMBDA 1]*M([1 2],2);
    bpred = [LVLAMBDA 1]*M([1 2],3);
    h(1)= plot([0 8], gpred/bpred*[1 1],'k:');
    h(2)= text(max(data(:,1)),gpred/bpred, 'No mac. pig.','FontSize',9)
    set(gca,'XLim',[0 8],'YLim',[2 9]);
    set(gcf,'PaperPositionMode','auto');
    
    % removing some lens pigment too.
    den_lens_smj = [3.2 2.8 2.4 1.985 1.62 1.308 1.05 0.84 0.675 0.557 0.468 0.393 0.335 0.29 0.26 0.24 0.225 0.215 0.203 0.191 0.18 0.168 0.162 0.151 0.145 0.133...
        0.128 0.122 0.116 0.11 0.104 0.099 0.093 0.087 0.081 0.075 0.07 0.064 0.058 0.052 0.046 0.041 0.036 0.031 0.028 0.024 0.021 0.017 0.014 0.012 0.009 0.007...
        0.005 0.003 0.002 0.001 zeros(1,25)]';
    % First two numbers of den_lens_smj are extrapolations.
    lens = den_lens_smj./1.16; % Dividing by 1.16 to make it open pupil
    lenstransmittance = 1./(10.^(lens/lens(5)));  % lens = 1 at 400 nm

    funds = funds./repmat(lenstransmittance,1,3)';  % removing lens pigment
    funds = funds./repmat(max(funds')',1,81);
    M = funds*spds;
    gpred = [LVLAMBDA 1]*M([1 2],2);
    bpred = [LVLAMBDA 1]*M([1 2],3);
    h(1)= plot([0 8], gpred/bpred*[1 1],'k:');
    h(2)= text(max(data(:,1)),gpred/bpred, 'No mac. pig., less lens','FontSize',9)

end

%%
% Section 14
% A horzontally oriented Gabor for the psychophysics methods
bkgndrgb = [.5 .5 .5];
gaborrgb = [.5 1 .5; .5 .5 1];
figure;
for i = 1:2
    subplot(2,2,i);
    im = DrawGaborEdge(bkgndrgb, gaborrgb(i,:)-bkgndrgb, [0 0 0], 0, 1.5, 1, 1, 0, 0, 0, 0, 0, .999, 45);
    image(im);
        
    set(gca,'XTick',[],'YTick',[],'visible','off');
    axis square;
end

%%
% Section 15
% A few exaple Gabors
filename = 'S041310002.nex';
stro = nex2stro(findfile(filename));
lmsidxs = [find(strcmp(stro.sum.trialFields(1,:),'lcont'))...
    find(strcmp(stro.sum.trialFields(1,:),'mcont'))...
    find(strcmp(stro.sum.trialFields(1,:),'scont'))];
ntrials = size(stro.trial(:,lmsidxs),1);
% -------------------------------
% Converting cone contrasts in nex file to 10 deg fundmentals.
% Must go through excitations first - can't just transform contrasts.
% -------------------------------
fundamentals = stro.sum.exptParams.fundamentals;
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = stro.sum.exptParams.mon_spd;
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;
load ('T_cones_smj10');
M10 = T_cones_smj10*mon_spd;
stro.trial(:,lmsidxs) = ConvertConeContrastBasis(M, M10, stro.sum.exptParams.bkgndrgb, stro.trial(:,lmsidxs));
% -------------------------------
out = NTpreprocess(stro,0,Inf);  % Getting the termination points
scaled = out(:,[2 3 4]).*repmat(out(:,5),1,3);
[planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);
coneweights = (planeparams'*xformmat');
nstim = size(scaled,1);
bkgndlms= M*[.5 .5 .5]';
for i = 1:nstim
    figure; axes;
    gaborrgb = inv(M)*(bkgndlms.*(1+scaled(i,:)/2)');
    im = DrawGaborEdge([.5 .5 .5], gaborrgb'-[.5 .5 .5], [0 0 0], 0, 1.5, 1, 1, 0, 0, 0, 0, 0, .999, 45);
    image(im);
end