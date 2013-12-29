% Code for looking at NeuroThresh datafiles off line.
stro = nex2stro;
spikeidx = strcmp(stro.sum.rasterCells(1,:),getSpikenum(stro,'first'));
Lbldone = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'bldone'));
coloridxs = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'coloridx'));
levels = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'levelidx'));
reversals = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'reversal'));
stepsize = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stepsize'));

lmsidxs = [find(strcmp(stro.sum.trialFields(1,:),'lcont'))...
            find(strcmp(stro.sum.trialFields(1,:),'mcont'))...
            find(strcmp(stro.sum.trialFields(1,:),'scont'))];

fpacq_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fp_acq'));
stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_on'));
stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_off'));
eot_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'eot'));

lms = stro.trial(:,lmsidxs);
lat = stro.sum.exptParams.latency;
mono = nanmax([0 stro.sum.exptParams.monopolar]);

Lbaseline = ~levels & all(lms == zeros(size(lms)),2);
Lprobe = ~levels & ~all(lms == zeros(size(lms)),2);
Lreplay = (levels == max(levels)) & ~stepsize;

spikerates = zeros(size(stro.trial,1),1);
for i = 1:size(stro.trial,1)
    spiketimes = stro.ras{i,spikeidx};
    nspikes = sum(spiketimes > stimon_t(i)+stro.sum.exptParams.latency/1000 & spiketimes < stimoff_t(i));
    spikerates(i) = nspikes./(stimoff_t(i)-stimon_t(i)-stro.sum.exptParams.latency/1000);
end

[uniquecoloridxs,tmp1,tmp2] = unique(coloridxs);
uniquecolordirs = mkbasis(lms(tmp1,:)')';

% getting M matrix etc.

fundamentals = stro.sum.exptParams.fundamentals;
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = stro.sum.exptParams.mon_spd;
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;

% Calculating a "data" matrix which is a summary of the data
% to be used in subsequent cells of the script
data = [];
for i = uniquecoloridxs'
    L = coloridxs == i;
    rev = reversals(logical(L));
    if (any(stepsize(L) == 0))
        continue;
    end
    lev = unique(levels(L));
    [u,s,v] = svd(lms(L,:));
    unitvector = v(:,1);
    contrasts = lms(L,:)*unitvector;
    if(sum(contrasts) < 0)
        contrasts = -contrasts;
        unitvector = -unitvector;
    end
    ss = stepsize(L);
    GV1 = GamutViolation(contrasts(end)*(1+ss(end)), unitvector, M, stro.sum.exptParams.bkgndrgb);
    GV2 = GamutViolation(contrasts(end)*(1-ss(end)), unitvector, M, stro.sum.exptParams.bkgndrgb);
    thresh = contrasts(end);
    %if (ss(end) > min(stepsize(stepsize > 0)))  % didn't finish for some reason
     if (sum(abs(rev)) < stro.sum.exptParams.nreversals)
        if (GV1)
            outofgamut = 1;
            thresh = thresh*(1+ss(end));
        elseif (GV2)
            outofgamut = 1;
            thresh = thresh*(1-ss(end));
        else
            outofgamut = nan;  % incomplete color direction
        end
    else
        outofgamut = 0;
    end
    if ~isnan(outofgamut)
        data = [data; i unitvector' thresh lev outofgamut];
    end
end
uniquecoloridxs = data(:,1);  % Ugly coding
data(:,1) = [];  % Ugly coding
Laccept = logical(ones(size(data,1),1));

%%
% Looking at the baseline measurements
figure;
bins = linspace(0,ceil(max(spikerates)),20);
binwidth = bins(2)-bins(1);
ntrialsbaseline = sum(Lbaseline);
ntrialsstim = sum(Lprobe);

subplot(2,2,1);
hist(spikerates(Lbaseline),bins); hold on;
ylabel('count');
xlabel('sp/sec');
title('Baseline trials');
set(gca,'XLim',[bins(1)-binwidth/2 bins(end)+binwidth/2]);
plot(stro.sum.exptParams.threshold, 0, 'm*');

% rasters
subplot(2,2,2); hold on;
trialcounter = 0;
stot = stimoff_t-stimon_t;
eott = eot_t-stimon_t;
for i = find(Lbaseline)'
    sp = stro.ras{i,1}-stimon_t(i);
    plot([sp sp]',[zeros(length(sp),1) .5*ones(length(sp),1)]'+trialcounter,'k-');
    plot([0 0],[0 1]+trialcounter,'m-','Linewidth',2);
    plot([stot(i) stot(i)],[0 1]+trialcounter,'m-','Linewidth',2);
    plot([eott(i) eott(i)],[0 1]+trialcounter,'m-','Linewidth',2);
    trialcounter = trialcounter+1;
end
set(gca,'XLim',[-.5 1],'YLim',[0 max(ntrialsbaseline,ntrialsstim)])

% Driven responses
subplot(2,2,3);
hist(spikerates(Lprobe),bins); hold on;
ylabel('count');
xlabel('sp/sec');
title('Stimulus trials');
set(gca,'XLim',[bins(1)-binwidth/2 bins(end)+binwidth/2]);
plot(stro.sum.exptParams.threshold, 0, 'm*');

% rasters
subplot(2,2,4); hold on;
trialcounter = 0;
stot = stimoff_t-stimon_t;
eott = eot_t-stimon_t;
for i = find(Lprobe)'
    sp = stro.ras{i,1}-stimon_t(i);
    plot([sp sp]',[zeros(length(sp),1) .5*ones(length(sp),1)]'+trialcounter,'k-');
    plot([0 0],[0 1]+trialcounter,'m-','Linewidth',2);
    plot([stot(i) stot(i)],[0 1]+trialcounter,'m-','Linewidth',2);
    plot([eott(i) eott(i)],[0 1]+trialcounter,'m-','Linewidth',2);
    trialcounter = trialcounter+1;
end
set(gca,'XLim',[-.5 1],'YLim',[0 max(ntrialsbaseline,ntrialsstim)])

% Comparing the two spatial phases of stimuli in the probe trials
% First finding the most common probe contrast
% 
lmstmp = stro.trial(Lprobe,lmsidxs);
a = unique(lmstmp,'rows');
counter = [];
for i = 1:size(a,1)
   phase1 = sum(all(repmat(a(i,:),size(lmstmp,1),1) == lmstmp,2));
   phase2 = sum(all(repmat(a(i,:),size(lmstmp,1),1) == -lmstmp,2));
   counter = [counter; phase1+phase2]; 
end
    
probecontrast = a(find(counter == max(counter),1),:);  % This logic doesn't work!
L1 = all(lmstmp == repmat(probecontrast,size(lmstmp,1),1),2);
L2 = all(lmstmp == -repmat(probecontrast,size(lmstmp,1),1),2);
srates = spikerates(Lprobe);
if (any(L1) & any(L2))
    if (ttest2(srates(L1),srates(L2)))
        pctchange = 100*abs(mean(srates(L1))-mean(srates(L2)))/max(mean(srates(L1))+mean(srates(L2)))
        title([num2str(pctchange,2),'% effect of phase!']);
    else
        title('No effect of phase');
    end
end


%%
% Looking for non-stationarities 
L = logical(~isnan(fpacq_t) & ~Lreplay);
prestimfr = [];
for i = find(L)'
    spiketimes = stro.ras{i,spikeidx};
    nspikes = sum(spiketimes > fpacq_t(i) & spiketimes < stimon_t(i));
    prestimfr = [prestimfr; nspikes./(stimon_t(i)-fpacq_t(i))];
end
figure;
subplot(3,1,1); hold on;
plot(prestimfr);
plot([1 length(prestimfr)],[stro.sum.exptParams.threshold stro.sum.exptParams.threshold],'k:');
ylabel('Prestim spike rate');
xlabel('Time (fp on''s)');
[b, bint] = regress(prestimfr, [1:sum(L); ones(1,sum(L))]');
if (sign(bint(1,1)) == sign(bint(1,2)))
    title('ALERT!  BASELINE FIRING RATE CHANGES');
else
    title('Non-significant change in baseline firing rate');
end

stot = stimon_t-fpacq_t;
subplot(3,1,2); hold on;
trialcounter = 1;
for j = find(L')
    sp = stro.ras{j,1}-fpacq_t(j);
    h = plot([sp sp]',[-.3*ones(length(sp),1) .3*ones(length(sp),1)]'+trialcounter,'k-');
    plot([0 0],[-.5 .5]+trialcounter,'m-','Linewidth',2);
    plot([stot(j) stot(j)],[-.5 .5]+trialcounter,'m-','Linewidth',2);
    trialcounter = trialcounter+1;
end
set(gca,'YLim',[-.5 trialcounter],'XLim',[-.1 max(stimoff_t-fpacq_t)+.1]);

% Stationarity check trials
L = levels > 0 & stepsize == 0 & ~Lreplay;
subplot(3,1,3); hold on;
trialcounter = 1;
stationaritychkfr = [];
for j = find(L')
    sp = stro.ras{j,spikeidx}-stimon_t(j);
    h = plot([sp sp]',[-.3*ones(length(sp),1) .3*ones(length(sp),1)]'+trialcounter,'k-');
    plot([0 0],[-.5 .5]+trialcounter,'m-','Linewidth',2);
    plot([stimoff_t(j) stimoff_t(j)]-stimon_t(j),[-.5 .5]+trialcounter,'m-','Linewidth',2);
    nspikes = sum(sp > 0 & sp < stimoff_t(j)-stimon_t(j));
    stationaritychkfr = [stationaritychkfr; nspikes./(stimoff_t(j)-stimon_t(j))];
    trialcounter = trialcounter+1;
end
set(gca,'XLim',[-.5 1]);
[b, bint] = regress(stationaritychkfr, [1:sum(L); ones(1,sum(L))]');
if (sign(bint(1,1)) == sign(bint(1,2)))
    title('Stationarity check trials - ALERT: significant trend');
else 
    title('Stationarity check trials - No significant trend');
end
%%
% Looking at threshold estimates in 3-D cone contrast space.
COLORDIMFACTOR = 1.2;
PAUSE = 0;
PLOTGAMUT = 0;
ADJUSTAXES = 2;  % 0, 1, or 2
PLOTMESH = 0;
SHOWACHROM = 0;
SHOWDONETRIANGLES = 1;

if (ADJUSTAXES == 1)
    scalemat = diag(1./sqrt(var([data(Laccept,[1 2 3]); -data(Laccept,[1 2 3])])));
elseif (ADJUSTAXES == 2)
    [v,d] = eig(cov([data(Laccept,[1 2 3]); -data(Laccept,[1 2 3])]));
    scalemat = v*1/sqrt(d);
else
    scalemat = eye(3);
end

figure; axes; hold on; set(gca,'Color',[.7 .7 .7]);
colors = colormap(jet(6));
plot3([-.1 .1],[0 0], [0 0],'k-');
plot3([0 0], [-.1 .1], [0 0],'k-');
plot3([0 0], [0 0],[-.1 .1],'k-');
if (ADJUSTAXES ~= 2)
    xlabel('L');
    ylabel('M');
    zlabel('S');
end

if (SHOWACHROM)
    plot3([-.2 .2],[-.2 .2],[-.2 .2], 'k-','LineWidth',2)
end

if (PLOTGAMUT && ~ADJUSTAXES)
    a = [0 0 0; 0 0 1; 0 1 1; 0 1 0; 0 0 0]';
    for i = 0:1
        a(1,:) = repmat(i,1,size(a,2));
        for j = 1:3
            a = [a([2,3],:); a(1,:)];
            b = M*a;
            bkgndlms = M*stro.sum.exptParams.bkgndrgb;
            cc = (b-repmat(bkgndlms,1,size(b,2)))./repmat(bkgndlms,1,size(b,2));
            plot3(cc(1,:),cc(2,:),cc(3,:),'k:')
        end
    end
end

if (max(data(:,5)) > size(colors,1))
    colors = [colors; repmat(colors(end,:), max(data(:,5))-size(colors,1),1)];
end
% red, green, blue, pink, aqua
tmp = data(Laccept,[1:3]) .*repmat(data(Laccept,4), 1,3);
%tmp = tmp./repmat(scalefactors, sum(Laccept),1);
tmp = tmp*scalemat;

uniquelevs = unique(data(Laccept,5))';
Loog = logical(data(Laccept,6));
for i  = uniquelevs
    for j = 0:1  % Out of gamut
        L = logical(data(Laccept,5) == i & Loog == j);
        if(any(L))
            plot3(tmp(L,1), tmp(L,2), tmp(L,3),'.','Color',colors(i,:)+j*((1-colors(i,:))/COLORDIMFACTOR));
            if (~mono)
                plot3(-tmp(L,1), -tmp(L,2), -tmp(L,3),'.','Color',colors(i,:)+j*((1-colors(i,:))/COLORDIMFACTOR));
            end
            if (j)
                plot3([zeros(sum(L),1) tmp(L,1)]', [zeros(sum(L),1) tmp(L,2)]', [zeros(sum(L),1) tmp(L,3)]','-','Color',colors(i,:)+j*((1-colors(i,:))/COLORDIMFACTOR));
                if (~mono)
                    plot3([zeros(sum(L),1) -tmp(L,1)]', [zeros(sum(L),1) -tmp(L,2)]', [zeros(sum(L),1) -tmp(L,3)]','-','Color',colors(i,:)+j*((1-colors(i,:))/COLORDIMFACTOR));
                end
            end
        end
    end
    if (PAUSE)
        if (SHOWDONETRIANGLES)
            if (exist('patchhandles'))
                delete(patchhandles);
            end
        end
        pause;
    end
end
if (PLOTMESH)
    for i = 0:1
        if (i == 1 && any(Loog))
            tmp(~Loog,:) = -tmp(~Loog,:);
        end
        T = DelaunayTri(tmp(~Loog,1),tmp(~Loog,2),tmp(~Loog,3));
        [tri, Xb] = convexHull(T);
        htm = trisurf(tri,tmp(~Loog,1),tmp(~Loog,2),tmp(~Loog,3),'FaceColor', 'green', 'FaceAlpha', 0.5, 'EdgeAlpha',.01);
    end
end

%%
% Here's an idea - filter trials based on (rank) correlation between
% contrast and response.
RTHRESH = .4;  % threshold on rank contrast-response correlation.
STBLTHRESH = 1; % Threshold on trjectory violations (positions of peaks/troughs relative to isresponse estimate)
TRAJVIOLTHRESH = 0; % Threshold on trajectory violations (actual vs. predicted based on off-line sorted spikes)

tmpdata = [];
for i = uniquecoloridxs'
    L = coloridxs == i;
        
    [u,s,v] = svd(lms(L,:));
    unitvector = v(:,1);
    if (sum(sign(unitvector)) < 0)
        unitvector = -unitvector;
    end
    contrasts = lms(L,:)*unitvector;
    if (sum(contrasts < 0))
       contrasts = -contrasts; 
    end
    contrastpeaks = contrasts(reversals(L) == 1);
    contrasttroughs = contrasts(reversals(L) == -1);
    
    ss = stepsize(L);
    GV1 = GamutViolation(contrasts(end)*(1+ss(end)), unitvector, M, stro.sum.exptParams.bkgndrgb);
    GV2 = GamutViolation(contrasts(end)*(1-ss(end)), unitvector, M, stro.sum.exptParams.bkgndrgb);
%     if (GV1)
%         tmpdata(3) = 1;
%         thresh = thresh*(1+ss(end));
%     elseif (GV2)
%         tmpdata(3) = 1;
%         thresh = thresh*(1-ss(end));
%     end
    if (GV1 | GV2)
        outofgamut = 1;
    else
        outofgamut = 0;
    end
    
    Loverthresh =  spikerates(L) > stro.sum.exptParams.threshold;
    skipit = 0;
    if (length(contrasts) > 1)
        Lviolation = xor(Loverthresh(1:end-1)', sign(diff(contrasts))'== -1);
    % Contrast trajectory not in agreement with that calculated from
    % post-hoc sorted spikes
    else
        Lviolation = nan;
    end
    
    if (sum(stepsize(L)) == 0)
        skipit = 1;
    end
    if (ss(end) > min(stepsize(stepsize > 0)) & ~outofgamut) % didn't finish in this color dir
        skipit = 1;
        keyboard
    end
    if (~skipit)
        r = corr([spikerates(L), contrasts],'type','Spearman');
        tmpdata = [tmpdata;  r(1,2), sum(contrastpeaks < contrasts(end)), sum(contrasttroughs > contrasts(end)), sum(Lviolation), outofgamut];
    end
end
% Columns:
% 1) Correlation between spikerate and contrast
% 2) Number of contrast peaks that exceed final contrast threshold estimate
% 3) Number of contrast troughs smaller than final contrast threshold estimate
% 4) Number of contrasts at which step direction disagrees with the predicted
% 5) Out of gamut

% Modify the Laccept vector which we use as a filter
% No longer destructively modifying the data matrix
Lfilter = logical(tmpdata(:,1) < RTHRESH & tmpdata(:,5) == 0); % Cor. filter if not oog
Lfilter = Lfilter | tmpdata(:,2) > STBLTHRESH;  % getting rid on unstable trajectories
Lfilter = Lfilter | tmpdata(:,4) > TRAJVIOLTHRESH;  % getting rid of trajectory violations
Laccept = ~Lfilter;

figure;
subplot(2,1,1); hold on;
Loog = logical(data(:,6));
plot(tmpdata(~Loog,2),tmpdata(~Loog,1),'k.');
plot(tmpdata(Loog,2),tmpdata(Loog,1),'y.');
plot(tmpdata(~Laccept,2),tmpdata(~Laccept,1),'rx');

plot([STBLTHRESH+.5 STBLTHRESH+.5],[-1 1],'k:');
plot([-1 4],[RTHRESH RTHRESH],'k:');
set(gca,'Xlim',[-1 4]);
subplot(2,1,2);
hist(tmpdata(~Loog,2),20);
corrcoef(tmpdata(~Loog,2),tmpdata(~Loog,1))
if (any(tmpdata,4))
    disp('Trajectory violation detected.');
end


%%
% An alternative representation: threshold as a function of 
% phi and theta (azimuth and elevation).
% Alternatively (this works too):
%    elev = atan2(projs(1),abs(projs(2))+abs(projs(3)));
%    h = polar(az,pi/2-elev,'ko');

figure; axes;hold on;
basis = mkbasis([1 1 0; 1 -1 0; 0 0 1]')';
for i = find(Laccept')
    projs = basis*data(i,[1 2 3])';
    if (projs(1) < 0)
        projs = -projs;
    end
    az = atan2(projs(3)/5, projs(2)); % L-M on the X-axis
    elev = acos(projs(1)./sqrt(sum(projs.^2)));
    h = polar(az,elev,'ko');
    c = 2./data(i,4)+2;
    set(h,'MarkerSize',c,'MarkerFaceColor',[1 1 1])
    if (data(i,6))  % OOG
       set(h,'MarkerFaceColor',[1 0 0]);
    end
end

x = linspace(0,2*pi,200);
plot(pi/2*cos(x),pi/2*sin(x),'k-');
axis square
xlabel('L-M');
ylabel('S');
set(gcf,'Name',stro.sum.fileName);


%%
% Fitting a pair of planes to the staircase termination points.
% 3/16/11 updating to use the new version of NTsurfacefit (that calls
% surfacefiterr4.m)

scaled = data(Laccept,[1:3]) .*repmat(data(Laccept,4), 1,3);
Loog = logical(data(Laccept,6));
[planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);
xyz = scaled*xformmat;

% What type of quadric?
A = [quadparams(1) quadparams(4) quadparams(5);...
    quadparams(4) quadparams(2) quadparams(6);...
    quadparams(5) quadparams(6) quadparams(3)];
d = eig(A);
if (sum(d<0) == 1)
    disp('hyperboloid of one sheet');
    plotlim = max(abs(xyz(:)))/1.1;
elseif (sum(d<0) == 2)
    disp('hyperboloid of two sheets');
    plotlim = max(abs(xyz(:)))/3;
else %(all(d>0))
    disp('Ellipsoid');
    plotlim = max(abs(xyz(:)))*1.1;
end

% Plotting the points
figure; axes; hold on;
plot3(xyz(~Loog,1),xyz(~Loog,2),xyz(~Loog,3),'k.');
plot3(-xyz(~Loog,1),-xyz(~Loog,2),-xyz(~Loog,3),'k.');
plot3(xyz(Loog,1),xyz(Loog,2),xyz(Loog,3),'y.');
plot3(-xyz(Loog,1),-xyz(Loog,2),-xyz(Loog,3),'y.');

% Plotting the planar fits
tmp1 = linspace(-mean(abs(xyz(:,1))),mean(abs(xyz(:,1))),2);
tmp2 = linspace(-mean(abs(xyz(:,2))),mean(abs(xyz(:,2))),2);
[xx yy] = meshgrid(tmp1,tmp2);
planez = -(1+planeparams(1)*xx+planeparams(2)*yy)/planeparams(3); % plane 1
h = surf(xx,yy,planez);
set(h,'FaceAlpha',.2);
planez = (1-planeparams(1)*xx-planeparams(2)*yy)/planeparams(3); % plane 2
h = surf(xx,yy,planez);
set(h,'FaceAlpha',.2);

% Plotting quadratic fits 
% using isosurface.m - crude, but effective.
tmplattice = linspace(-plotlim,plotlim,20);
[xx yy zz] = meshgrid(tmplattice, tmplattice,tmplattice);
variables = [xx(:).^2 yy(:).^2 zz(:).^2 2*xx(:).*yy(:) 2*xx(:).*zz(:) 2*yy(:).*zz(:)];
fr = variables*quadparams;
p = patch(isosurface(xx,yy,zz,reshape(fr,size(xx)), 1));
set(p, 'FaceColor', 'blue', 'EdgeColor', 'none');
set(p,'FaceAlpha',.2,'FaceColor','green','Edgealpha',0);
axis vis3d;
camlight;
lighting phong;
axis square;



%%
% Continuing from above, comparing to grating data.
GT = nex2stro();
GTstruct = getGratingTuning(GT);
colors = GTstruct.color.colors;
rotcols =  colors*xformmat;
planepredresp = stro.sum.exptParams.threshold*abs(rotcols*planeparams); % strong linearity assumption
trueresp = GTstruct.color.colresp(:,1);
%Outdated code (uses old "simple quadratic")
%A = (quadparams(4).*rotcols(:,1).^2+quadparams(5).*rotcols(:,2).^2+quadparams(6).*rotcols(:,1).*rotcols(:,2))./(r.^2);
%B =(quadparams(1).*rotcols(:,1)+quadparams(2).*rotcols(:,2)+quadparams(3).*rotcols(:,3))./r;
%quadpredr1 = 2./(-B+sqrt(B.^2-4.*A));
%quadpredr2 = 2./(-B-sqrt(B.^2-4.*A));
%quadpredr1(real(quadpredr1) ~= quadpredr1) = Inf;
%quadpredr2(real(quadpredr2) ~= quadpredr2) = Inf;
%quadpredresp = stro.sum.exptParams.threshold.*(r./min([abs(quadpredr1), abs(quadpredr2)],[],2));

[th,ph,r] = cart2sph(rotcols(:,1),rotcols(:,2),rotcols(:,3));
predr2 = 1./(quadparams(1).*(cos(ph).*cos(th)).^2 +...
    quadparams(2).*(cos(ph).*sin(th)).^2 +...
    quadparams(3).*sin(ph).^2 + ...
    2*quadparams(4).*cos(ph).*cos(th).*cos(ph).*sin(th) +...
    2*quadparams(5).*cos(ph).*cos(th).*sin(ph) +...
    2*quadparams(6).*cos(ph).*sin(th).*sin(ph));
predr2(predr2<0) = Inf;
quadpredresp = stro.sum.exptParams.threshold.*(r./sqrt(predr2));

ylabels = {'Linear','Quadratic'};
figure;
for j = 1:2
    subplot(3,1,j); hold on;
    for i = 1:length(trueresp)
        if (j == 1)
            h = plot(trueresp(i), planepredresp(i), 'ko','Markersize',12);
        elseif (j == 2)
            h = plot(trueresp(i), quadpredresp(i), 'ko','Markersize',12);
        else % j == 3
            h = plot(trueresp(i), nppredresp(i), 'ko','Markersize',12);
        end
        col = GTstruct.color.colors(i,:).*[1 1 .2];
        col = col./(2*max(abs(col))) + .5;
        set(h,'MarkerFaceColor',col,'ButtonDownFcn',['disp(num2str(GTstruct.color.colors(',num2str(i),',:)))']);
    end
    xlabel('Response'); ylabel([ylabels{j},' prediction']);
end
subplot(3,1,1);
title(stro.sum.fileName);
corrcoef([trueresp, planepredresp, quadpredresp])
% Order: Actual responses, Linear, Quadratic


%%
% Continuing from above.
% Plotting fitted surfaces in cone contrast space (and also plotting the
% locations of the grating stimuli).
% Needs to be changed
colors = GTstruct.color.colors;

% First plotting the positions of the grating stimuli
figure; axes; hold on;
for i = 1:size(colors,1)
   h1 = plot3(colors(i,1),colors(i,2),colors(i,3),'o','MarkerFaceColor',colors(i,:)/2+.5, 'Color',colors(i,:)/2+.5); 
   h2 = plot3(-colors(i,1),-colors(i,2),-colors(i,3),'o','MarkerFaceColor',-colors(i,:)/2+.5,'Color',-colors(i,:)/2+.5); 
   set(h1,'MarkerSize',GTstruct.color.colresp(i,1)./max(GTstruct.color.colresp(:,1))*14+1);
   set(h2,'MarkerSize',GTstruct.color.colresp(i,1)./max(GTstruct.color.colresp(:,1))*14+1); 
end

% Now plotting surfaces
plotlims = 2*max(abs(colors));
[xx yy zz] = meshgrid(linspace(-plotlims(1),plotlims(1),20),...
                      linspace(-plotlims(2),plotlims(2),20),...
                      linspace(-plotlims(3),plotlims(3),20));
xformedxyz = [xx(:) yy(:) zz(:)]*xformmat;
variables = [xformedxyz(:,1).^2 xformedxyz(:,2).^2 xformedxyz(:,3).^2 2*xformedxyz(:,1).*xformedxyz(:,2) 2*xformedxyz(:,1).*xformedxyz(:,3) 2*xformedxyz(:,2).*xformedxyz(:,3)];
for PLANE = 0:1
    if (PLANE)
        coefficients = [planeparams'.*planeparams' planeparams(1)*planeparams(2) planeparams(1)*planeparams(3) planeparams(2)*planeparams(3)]';
        patchcolor = 'green';
    else
        coefficients = quadparams;
        patchcolor = 'blue';
    end
    
    fr = variables*coefficients;
    surfstruct = isosurface(xx,yy,zz,reshape(fr,size(xx)), 1);
  %  surfstruct.vertices = surfstruct.vertices/xformmat;
    
    p = patch(surfstruct);
    set(p, 'FaceColor', 'blue', 'EdgeColor', 'none');
    set(p,'FaceAlpha',.2,'FaceColor',patchcolor,'Edgealpha',0);
    axis vis3d;
    camlight;
    lighting phong;
    axis square;
end

%% 
% Continuing from above, bootstrapping test comparing linear and quadratic fits to isoresponse
% surface.

nbootiter = 100;

scaled = data(Laccept,[1:3]) .*repmat(data(Laccept,4), 1,3);
Loog = logical(data(Laccept,6));
[planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);
xyz = scaled*xformmat;

[th,ph,r] = cart2sph(xyz(~Loog,1),xyz(~Loog,2),xyz(~Loog,3));
planer = -1./(planeparams(1).*cos(ph).*cos(th)+planeparams(2).*cos(ph).*sin(th)+planeparams(3).*sin(ph));
planer = abs(planer);  % this ok?
resid = log(r)-log(planer);
    
nulldist = nan*ones(nbootiter,1);
wait_h = waitbar(0,'Bootstrapping...');
disp(['Bootstrapping.  PlaneSSE = ',num2str(planeSSE),' QuadSSE = ',num2str(quadSSE)]);
options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-6,'TolX',10^-6);
SSEs = [];
for j = 1:nbootiter
    waitbar(j/nbootiter, wait_h);
    tmpresid = exp(resid(unidrnd(length(resid),[1 length(resid)])));
       
    tmpx = xyz(:,1);
    tmpy = xyz(:,2);
    tmpz = xyz(:,3);
        
    tmpx(~Loog) = planer.*tmpresid.*cos(ph).*cos(th);
    tmpy(~Loog) = planer.*tmpresid.*cos(ph).*sin(th);
    tmpz(~Loog) = planer.*tmpresid.*sin(ph);
        
    [~, SSE1, exitflag1] = fminsearch(@(x) surfacefiterr2([tmpx tmpy tmpz],x, Loog),planeparams, options);
    initguess = [planeparams'.*planeparams' planeparams(1)*planeparams(2) planeparams(1)*planeparams(3) planeparams(2)*planeparams(3)];
    [tmp, SSE2, exitflag2] = fminsearch(@(x) surfacefiterr4([tmpx tmpy tmpz],x, Loog),initguess, options);
    if (~exitflag1)
        disp('Bonked on a plane fit');
    end
    if (~exitflag2)
        disp('Bonked on a quadratic fit');
    end
    if (~exitflag1 | ~exitflag2)
        j = j-1;
    end
    nulldist(j)= SSE1./SSE2;
    SSEs(j,:) = [SSE1 SSE2];
end
close(wait_h);
p = sum(nulldist>planeSSE/quadSSE)./nbootiter;
disp(['p-value is ',num2str(p),' using ',num2str(nbootiter),' iterations']);


%%
% Looking at data from the replay control experiment (Which doesn't
% actually exist yet).
% GDLH: Not including Laccept right now - need to inlude it?

scaled = data(:,[1:3]) .*repmat(data(:,4), 1,3);  % include Laccept later?

tmp = []; clear h;
figure; axes; hold on;
for i = 1:size(scaled,1)
    Lsamecolordir = softEq(mkbasis(lms')'*data(i,[1 2 3])',1); % includes Lonsurface
    if (any(Lsamecolordir&Lreplay))
        lmsincoldir = unique(lms(Lsamecolordir&Lreplay,:),'rows')
        [~,idx] = sort(sum(lmsincoldir.^2,2));
        lmsincoldir = lmsincoldir(idx,:);
        for j = 1:size(lmsincoldir,1)
            Lendstair = all(softEq(repmat(lmsincoldir(j,:),size(lms,1),1), lms,6),2) & Lreplay;
            sum(Lendstair)
            sr = spikerates(Lendstair&Lreplay);
            [lms(Lendstair,:) sr]

            h(1) = plot3(lmsincoldir(j,1),lmsincoldir(j,2),lmsincoldir(j,3),'ko','MarkerSize',mean(sr));
            h(2) = plot3(-lmsincoldir(j,1),-lmsincoldir(j,2),-lmsincoldir(j,3),'ko','MarkerSize',mean(sr));
            if (j == 2)
                set(h,'MarkerfaceColor','red');
                plot3([lmsincoldir(1,1) lmsincoldir(2,1)],...
                    [lmsincoldir(1,2) lmsincoldir(2,2)],...
                    [lmsincoldir(1,3) lmsincoldir(2,3)]);
                plot3(-[lmsincoldir(1,1) lmsincoldir(2,1)],...
                    -[lmsincoldir(1,2) lmsincoldir(2,2)],...
                    -[lmsincoldir(1,3) lmsincoldir(2,3)]);
            end
            if (j == 1 & size(lmsincoldir,1) > 1)
               set(h,'MarkerfaceColor','black');
            end  
            tmp = [tmp; i j sum(Lendstair) mean(sr) std(sr) lmsincoldir(j,:)];
        end
    end
end

size(tmp)
sum(tmp(:,3)) == sum(Lreplay)  % Sanity check
%hist(tmp(:,2),15)  % Bimodal?

bins = linspace(min(tmp(:,4)),max(tmp(:,4)),15);
L = tmp(:,2) == 2;
[n1,x] = hist(tmp(L,4),bins);
[n2,x] = hist(tmp(~L,4),bins);
figure; axes; hold on;
bar(bins, n1+n2,'red');
bar(bins, n2,'black');

figure; axes; hold on; counter = 0; 
for i = unique(tmp(:,1))'
    Lcoldir = tmp(:,1) == i;
    if (sum(Lcoldir) > 1)
        for j = 1:2
            L = Lcoldir & tmp(:,2) == j;
            if (sum(L) > 1)
                disp('Problem')
                keyboard;
            end
            h(1) = plot(counter, tmp(L,4),'ko','MarkerFaceColor','black');
            h(2) = plot([counter counter], tmp(L,4) + [-tmp(L,5) tmp(L,5)]./sqrt(tmp(L,3)),'k-'); 
            if (j == 2)
                set(h,'MarkerFaceColor','red','Color','red');
            end
        end
        counter = counter + 1;
    end
end
%%
% % OK, now fitting surfaces
% scaled = data(Laccept,[1:3]) .*repmat(data(Laccept,4), 1,3);
% mn = mean(scaled(~Loog,:));
% 
% pcs = princomp(scaled(~Loog,:));
% initxyz = -sign(mn*pcs(:,3))* pcs(:,3)./norm(mn*pcs(:,3));
% options = optimset('MaxFunEvals',50000,'MaxIter',50000);
% % Ignoring OOG points
% [planeparams_naive, planeSSE, exitflag, out4] = fminsearch(@(x) surfacefiterr2(scaled(~Loog,:),x),initxyz, options);
% % Now taking OOG points into account
% %[planeparams, planeSSE, exitflag, out4] = fminsearch(@(x) surfacefiterr(scaled,x,Loog),planeparams_naive, options);
% planeparams = planeparams_naive;
% planepredresp = stro.sum.exptParams.threshold*abs(GTstruct.color.colors*planeparams_naive);
% trueresp = GTstruct.color.colresp(:,1);
% 
% % Now getting the quadratic prediction
% % Trying to deal with OOG points.
% %[quadparams, quadSSE, exitflag, out4] = fminsearch(@(x) surfacefiterr(scaled(~Loog,:),x),[planeparams_naive;0;0;0;0;0;0], options);
% %quadSSE
% [quadparams, quadSSE, exitflag, out4] = fminsearch(@(x) surfacefiterr(scaled,x,Loog),[planeparams_naive;0;0;0;0;0;0], options);
% quadSSE
% if (~exitflag)
%     keyboard
% end
% for j = 1:length(trueresp)
%     x = GTstruct.color.colors(j,1);
%     y = GTstruct.color.colors(j,2);
%     z = GTstruct.color.colors(j,3);
%     a = quadparams(1);
%     b = quadparams(2);
%     c = quadparams(3);
%     d = quadparams(4);
%     e = quadparams(5);
%     f = quadparams(6);
%     g = quadparams(7);
%     h = quadparams(8);
%     i = quadparams(9);
% 
%     r = sqrt(x.^2+y.^2+z.^2);
%     A = (d*x.^2+e*y.^2+f*z.^2+g.*x.*y+h.*y.*z+i.*x.*z)./(r.^2);
%     B = (a*x+b*y+c*z)./r;
%     predr1 = 2./(-B+sqrt(B.^2-4.*A));
%     predr2 = 2./(-B-sqrt(B.^2-4.*A));
%     if (~isreal(predr1) & ~isreal(predr2))
%         predr = Inf;
%     else
%         predr = min(abs([predr1, predr2]),[],2);
%     end
%     quadpredresp(j) = stro.sum.exptParams.threshold*norm([x;y;z])/abs(predr);
% end
% 
% % Now a nonparametric estimate.
% % Taking the midpoint through the isoresponse volume along particular radii 
% T = delaunay3(scaled(~Loog,1),scaled(~Loog,2),scaled(~Loog,3));
% k = convhulln(scaled(~Loog,[1 2 3]));
% trueresp = GTstruct.color.colresp(:,1);
% nppredresp = [];
% LIM = 1;
% for i = 1:length(trueresp)
%     c = GTstruct.color.colors(i,:);
%     v = linspace(-LIM,LIM,100000);
%     test = inhull(v'*c,scaled(~Loog,[1 2 3]),k,10^-2);
%     if any(test)
%         nppredresp(i) = stro.sum.exptParams.threshold*norm(c)/min(abs(v(test)));
%     else
%         nppredresp(i) = 0;
%     end
% end
% 
% ylabels = {'Linear','Quadratic','Nonparametric'};
% figure;
% for j = 1:3
%     subplot(3,1,j); hold on;
%     for i = 1:length(trueresp)
%         if (j == 1)
%             h = plot(trueresp(i), planepredresp(i), 'ko','Markersize',12);
%         elseif (j == 2)
%             h = plot(trueresp(i), quadpredresp(i), 'ko','Markersize',12);
%         else % j == 3
%             h = plot(trueresp(i), nppredresp(i), 'ko','Markersize',12);
%         end
%         col = GTstruct.color.colors(i,:).*[1 1 .2];
%         col = col./(2*max(abs(col))) + .5;
%         set(h,'MarkerFaceColor',col,'ButtonDownFcn',['disp(num2str(GTstruct.color.colors(',num2str(i),',:)))']);
%     end
%     xlabel('Response'); ylabel([ylabels{j},' prediction']);
% end
% subplot(3,1,1);
% title(stro.sum.fileName);
% corrcoef([planepredresp, quadpredresp', nppredresp', trueresp])
% 
% %%
% % Come here after plotting data in 3-D cone space.
% for i = 1:size(GTstruct.color.colors,1)
%     a = [GTstruct.color.colors(i,1),GTstruct.color.colors(i,2),GTstruct.color.colors(i,3)];
%     b = GTstruct.color.colresp(i)./max(GTstruct.color.colresp(:,1));
%     plot3([-a(1) a(1)],[-a(2) a(2)],[-a(3) a(3)],'m-*','Linewidth',2*b+1);
% end


% %%
% % Residual distances from origin to newly fitted plane
% % % We know the mean vector lies on the best fit plane, we know how long it
% % % is, and we know the nearest distance from the origin to the plane is in
% % % the direction of the normal vector.
% th = acos(data(:,2));
% ph = atan2(data(:,3),data(:,1));
% L = logical(ph<0);
% ph(L) = -ph(L);
% ph(~L) = 2*pi-ph(~L);
% 
% predr = -1./(pointertoplane(1).*sin(th).*cos(ph)+pointertoplane(2).*cos(th)-pointertoplane(3).*sin(th).*sin(ph));
% resid = data(:,4)./predr;
% 
% figure;
% subplot(2,1,1); 
% hist(resid(~Loog),20);
% xlabel('Proportional residual (actual/pred)');
% ylabel('Count');
% % Negative residuals indicate points that never should have hit the plane
% % (Unless they were going backwards - e.g. we made the point assignment
% % wrong).
% 
% [a,b,c] = princomp(scaled*basis(:,[2 3]));
% xformed = b./repmat(sqrt(c'),size(b,1),1);
% 
% [xi,yi] = meshgrid(linspace(min(xformed(:,1)),max(xformed(:,1)),100),...
%     linspace(min(xformed(:,1)),max(xformed(:,1)),100));
% zi = griddata(xformed(~Loog,1), xformed(~Loog,2), resid(~Loog),xi,yi);
% 
% subplot(2,1,2); hold on;
% surf(xi,yi,zi); 
% h = plot3(xformed(~Loog,1),xformed(~Loog,2),resid(~Loog),'r.','Markersize',20);
% plot3([zeros(1,sum(Loog));xformed(Loog,1)'],...
%     [zeros(1,sum(Loog));xformed(Loog,2)'],...
%     [zeros(1,sum(Loog));resid(Loog)'],'k-');
% plot3([zeros(1,sum(Loog&Lshouldhavehitplane));xformed(Loog&Lshouldhavehitplane,1)'],...
%     [zeros(1,sum(Loog&Lshouldhavehitplane));xformed(Loog&Lshouldhavehitplane,2)'],...
%     [zeros(1,sum(Loog&Lshouldhavehitplane));resid(Loog&Lshouldhavehitplane)'],'r-');
% set(gca,'XLim',[min(xformed(~Loog,1)) max(xformed(~Loog,1))]);
% set(gca,'YLim',[min(xformed(~Loog,2)) max(xformed(~Loog,2))]);
% set(gca,'ZLim',[min(resid(~Loog)) max(resid(~Loog))]);
% xlabel('Basis 1');ylabel('Basis 2');zlabel('proportional resid');
% 
% 
% 
% %%
% %%%%%%
% % Looking at the progression of stimulus contrasts
% % This will hopefully be the start of a set of routines that will filter
% % individual threshold estimates based on the strangeness of the trajectory
% % they terminated. 
% stot = stimoff_t-stimon_t;
% for i = uniquecoloridxs'
%     figure
%     L = coloridxs == i;
%     [u,s,v] = svd(lms(L,:));
%     unitvector = v(:,1);
%     if (sum(sign(unitvector)) < 0)
%         unitvector = -unitvector;
%     end
%     contrasts = lms(L,:)*unitvector;
%     if (sum(contrasts < 0))
%        contrasts = -contrasts; 
%     end
%     contrastpeaks = contrasts(reversals(L) == 1);
%     contrasttroughs = contrasts(reversals(L) == -1);
% 
%     % Now the actual plotting
%     subplot(3,2,1); hold on;
%     plot(contrasts,'ko-','MarkerSize',4,'MarkerFaceColor',uniquecolordirs(i,:)/2+.5);
%     set(gca,'XLim',[0 length(contrasts)],'XTick',[]);
%     
%     % Normal orientation rasters
%     subplot(3,2,2); hold on;
%     trialcounter = 1;
%     for j = find(L')
%         sp = stro.ras{j,1}-stimon_t(j);
%         h = plot([sp sp]',[-.3*ones(length(sp),1) .3*ones(length(sp),1)]'+trialcounter,'k-');
%         plot([0 0],[-.5 .5]+trialcounter,'m-','Linewidth',2);
%         plot([stro.sum.exptParams.latency stro.sum.exptParams.latency]/1000,[-.5 .5]+trialcounter,'m:','Linewidth',2);
%         plot([stot(j) stot(i)],[-.5 .5]+trialcounter,'m-','Linewidth',2);
%         trialcounter = trialcounter+1;
%     end
%     set(gca,'YLim',[-.5 trialcounter],'XLim',[-.1 max(stimoff_t-stimon_t)+.1]);
% 
%     % Contrast as a function of trial number
%     subplot(3,2,3); hold on;
%     plot(spikerates(L),'ko-','MarkerSize',3,'Color',uniquecolordirs(i,:)/4+.5);
%     plot([0 length(contrasts)], repmat(stro.sum.exptParams.threshold,1,2),'k:');
%     set(gca,'XLim',[0 length(contrasts)],'XTick',[]);
%     Loverthresh =  spikerates(L) > stro.sum.exptParams.threshold;
%     Lviolation = xor(Loverthresh(1:end-1)', sign(diff(contrasts))'== -1);
%     if any(Lviolation) % contrast trajectory differs from expected.  Probably due to bad isolation
%         plot(find(Lviolation)+1,0,'rx');
%     end
%     
%     subplot(3,2,5); hold on;
%     trialcounter = 1;
%     for j = find(L')
%         sp = stro.ras{j,1}-stimon_t(j);
%         h = plot([-.3*ones(length(sp),1) .3*ones(length(sp),1)]'+trialcounter, [sp sp]','k-');
%         plot([-.5 .5]+trialcounter,[0 0],'m-','Linewidth',2);
%         plot([-.5 .5]+trialcounter, [stro.sum.exptParams.latency stro.sum.exptParams.latency]/1000,'m:','Linewidth',2);
%         plot([-.5 .5]+trialcounter,[stot(j) stot(i)],'m-','Linewidth',2);
%         trialcounter = trialcounter+1;
%     end
%     set(gca,'XLim',[-.5 trialcounter],'YLim',[-.1 max(stimoff_t-stimon_t)+.1]);
% 
%     subplot(3,2,4); hold on;
%     hist(contrasts);
%     if (~isempty(contrastpeaks))
%         plot(contrastpeaks,3,'m^');
%         plot(contrasttroughs,3,'bv');
%     end
%     plot(contrasts(end),3,'g*');
% 
%     subplot(3,2,6); hold on;
%     plot(contrasts, spikerates(L),'k*');
%     plot([contrasts(end), contrasts(end)],[0 max(spikerates(L))],'k:');
%     plot([min(contrasts) max(contrasts)],[stro.sum.exptParams.threshold stro.sum.exptParams.threshold],'k:');
%     set(gcf,'Name',['condition ',num2str(i)]);
%   
%     %sum(contrastpeaks < contrasts(end))/length(contrastpeaks)
%     %sum(contrasttroughs > contrasts(end))/length(contrasttroughs)
% 
%     pause
% end
% 
% 
% %%
% % Fitting a linear model by PCA
% 
% scaled = data(:,1:3) .*repmat(data(:,4), 1,3);
% figure; axes; hold on;
% plot3(scaled(~Loog,1),scaled(~Loog,2), scaled(~Loog,3),'k.');
% plot3(-scaled(~Loog,1),-scaled(~Loog,2),-scaled(~Loog,3),'b.');
% plot3([zeros(sum(Loog),1) scaled(Loog,1)]', [zeros(sum(Loog),1) scaled(Loog,2)]', [zeros(sum(Loog),1) scaled(Loog,3)]','y.-');
% plot3([zeros(sum(Loog),1) -scaled(Loog,1)]', [zeros(sum(Loog),1) -scaled(Loog,2)]', [zeros(sum(Loog),1) -scaled(Loog,3)]','y.-');
% xlabel('L');
% ylabel('M');
% zlabel('S');
% 
% [coeffs, score, roots]= princomp(scaled(~Loog,:));
% mn = mean(scaled(~Loog,:));
% x = linspace(0,2*pi,100)';
% X = [cos(x)/5 sin(x)/5];
% y = X*coeffs(:,[1:2])';
% linmod = plot3(y(:,1)+mn(1), y(:,2)+mn(2), y(:,3)+mn(3),'k:');
% title(num2str(coeffs(:,3)));
% 
% %%
% % residuals in the PC1/PC2 plane
% 
% figure;
% mn = mean(data(:,[1 2 3]));
% xformed = (data(:,[1 2 3])-repmat(mn,size(data,1),1))*coeffs;
% %xformed = xformed./repmat(sqrt(roots)',size(xformed,1),1);
% % var(xformed(~Loog,:))  % Sanity check
% 
% [xi,yi] = meshgrid(linspace(min(xformed(:,1)),max(xformed(:,1)),50),...
%     linspace(min(xformed(:,2)),max(xformed(:,2)),50));
% zi = griddata(xformed(~Loog,1), xformed(~Loog,2), xformed(~Loog,3),xi,yi,'linear');
% surf(xi,yi,zi); hold on;
% 
% plot3(xformed(~Loog,1),xformed(~Loog,2),xformed(~Loog,3),'r.');
% plot3([zeros(1,sum(Loog));xformed(Loog,1)'],...
%     [zeros(1,sum(Loog));xformed(Loog,2)'],...
%     [zeros(1,sum(Loog));xformed(Loog,3)'],'k-');
% 
% xlabel(['PC 1: ',num2str(coeffs(:,1)',2)]);
% ylabel(['PC 2: ',num2str(coeffs(:,2)',2)]);
% zlabel(['PC 3: ',num2str(coeffs(:,3)',2)]);
% 
% set(gca,'XLim',[min(xformed(~Loog,1)) max(xformed(~Loog,1))].*[.9 1.1]);
% set(gca,'YLim',[min(xformed(~Loog,2)) max(xformed(~Loog,2))].*[.9 1.1]);
% set(gca,'ZLim',[min(xformed(~Loog,3)) max(xformed(~Loog,3))].*[.9 1.1]);
% 
% %%
% % Residual distances from origin to plane
% tmp = data(Laccept,[1:3]) .*repmat(data(Laccept,4), 1,3);
% mn = mean(tmp);
% mnamp = norm(mn);
% % % We know the mean vector lies on the best fit plane, we know how long it
% % % is, and we know the nearest distance from the origin to the plane is in
% % % the direction of the normal vector.
% costheta = (mn./mnamp)*coeffs(:,3);  % angle from mean to normal vector
% disttoplane = mnamp*costheta;
% pointertoplane = coeffs(:,3)*disttoplane;
% 
% figure; axes; hold on;
% plot3(tmp(:,1),tmp(:,2),tmp(:,3),'k.');
% plot3([0 mn(1)],[0 mn(2)],[0 mn(3)],'g*-');
% plot3([0 pointertoplane(1)],[0 pointertoplane(2)],[0 pointertoplane(3)],'m*-');
% plot3([zeros(1,sum(Loog));tmp(Loog,1)'],...
%     [zeros(1,sum(Loog));tmp(Loog,2)'],...
%     [zeros(1,sum(Loog));tmp(Loog,3)'],'y-');
% 
% costheta = mkbasis(tmp')'*coeffs(:,3);
% pred = disttoplane./costheta;
% resid = data(:,4)./pred;
% figure;
% subplot(2,1,1); 
% hist(resid(~Loog),20);
% xlabel('Proportional residual (actual/pred)');
% ylabel('Count');
% % Negative residuals indicate points that never should have hit the plane
% % (Unless they were going backwards - e.g. we made the point assignment
% % wrong).
% 
% xformed = tmp*coeffs*diag(1./sqrt(roots));
% xformed = xformed.*repmat(sign(xformed(:,3)), 1,3);
% [xi,yi] = meshgrid(linspace(min(xformed(:,1)),max(xformed(:,1)),100),...
%     linspace(min(xformed(:,2)),max(xformed(:,2)),100));
% zi = griddata(xformed(~Loog,1), xformed(~Loog,2), resid(~Loog),xi,yi);
% 
% subplot(2,1,2); hold on;
% surf(xi,yi,zi); 
% h = plot3(xformed(~Loog,1),xformed(~Loog,2),resid(~Loog),'r.');
% plot3([zeros(1,sum(Loog));xformed(Loog,1)'],...
%     [zeros(1,sum(Loog));xformed(Loog,2)'],...
%     [zeros(1,sum(Loog));resid(Loog)'],'k-');
% set(gca,'XLim',[min(xformed(~Loog,1)) max(xformed(~Loog,1))]);
% set(gca,'YLim',[min(xformed(~Loog,2)) max(xformed(~Loog,2))]);
% set(gca,'ZLim',[min(resid(~Loog)) max(resid(~Loog))]);
% xlabel('PC1');ylabel('PC2');zlabel('proportional resid');
% 
% 
% %%
% % making sure things are getting dropped correctly
% 
% overthresh = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'overthresh'));
% figure; axes;
% for i = uniquecoloridxs'
%     cla;
%     hold on
%     L = logical(coloridxs == i);
%     plot(lms(L,:));
%     plot(reversal(L)*.8,'k*');
%     plot(stepsize(L),'m:');
%     plot(overthresh(L),'m*');
%     pause
% end
% 
% 
% %%
% % Contrast-response functions
% % Sort of out-dated
% type = 'linear';
% 
% nspikes = zeros(size(stro.trial,1),1);
% for i = 1:size(stro.trial,1)
%     spiketimes = stro.ras{i,spikeidx};
%     nspikes(i) = sum(spiketimes > stimon_t(i) & spiketimes < stimoff_t(i));
% end
% 
% uniquecoloridxs = unique(coloridxs);
% uniquecoloridxs([1:2]) = [];  % hack to avoid baseline measurements.
% 
% blnspikes = nspikes(Lbaseline&~Lprobe);
% 
% data = [];
% for i = uniquecoloridxs'
%     figure
%     L = coloridxs == i;
%     [u,s,v] = svd(lms(L,:));
%     unitvector = v(:,1);
%     if (sum(sign(unitvector)) < 0)
%         unitvector = -unitvector;
%     end
%     contrasts = lms(L,:)*unitvector;
%     if(sum(contrasts) < 0)
%         contrasts = -contrasts;
%     end
%     axes; hold on;
%     plot(contrasts, spikerates(L),'*','Color',uniquecolordirs(i,:)/2+.5);
%     title(num2str(unitvector','%.1f %.1f %.1f'));
%        
%     % Including baseline blank trials (0 contrast)
%     con = [zeros(sum(Lbaseline&~Lprobe),1); contrasts];
%     nsp = [blnspikes; nspikes(L)];
% 
%     [baseline, thresh, slope, success] = contrastrespfit(con, nsp, type);
%     [baseline, thresh, slope, success]
% 
%     plot([0 thresh],[baseline baseline],'k:');
%     if (strcmp(type,'linear'))
%         plot([thresh max(contrasts)],[0 max(contrasts)-thresh]*slope+baseline,'k:');
%         fittedthresh = (stro.sum.exptParams.threshold-baseline)/slope+thresh;
%     elseif(strcmp(type,'quadratic'))
%         x = linspace(thresh, max(contrasts),100);
%         y = baseline+(slope*(x-thresh).^2);
%         plot(x,y,'k:');
%         fittedthresh = sqrt((stro.sum.exptParams.threshold-baseline)/slope)+thresh;
%     end
%     data = [data; fittedthresh];
%     
%     xlabel(thresh);
%     set(gcf,'name',stro.sum.fileName);
% end
% figure;
% bar(data);
% set(gcf,'name',stro.sum.fileName);
% %%
% % Comparing two ways of estimating threshold
% % 1) parametric function fit
% % 2) last contrast used
% % Goes through a list of files.
% % The moral of this story appears to be that just using the last point in
% % the staircase does a better job (in terms of being more consistent across
% % repeats) than fitting a function.
% 
% datapath = 'C:\NO BACKUP\NexFiles\Kali\';
% filenames = {'K061909003.nex','K061909007.nex','K061909008.nex','K061909010.nex'};
% filenames = {'K062209003.nex','K062209004.nex','K062209005.nex','K062209006.nex'};
% 
% type = 'linear';
% data = [];
% for i = 1:length(filenames)
%     stro = nex2stro([datapath, filenames{i}]);
%     spikeidx = strcmp(stro.sum.rasterCells(1,:),getSpikenum(stro));
%     Lbldone = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'bldone'));
%     coloridxs = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'coloridx'));
%     contrasts = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'contrast'));
%     lmsidxs = [find(strcmp(stro.sum.trialFields(1,:),'lcont'))...
%                 find(strcmp(stro.sum.trialFields(1,:),'mcont'))...
%                 find(strcmp(stro.sum.trialFields(1,:),'scont'))];
%         
%     stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_on'));
%     stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_off'));
%     eot_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'eot'));
% 
%     lms = stro.trial(:,lmsidxs);
%     lat = stro.sum.exptParams.latency;
%     Lbaseline = ~Lbldone & all(lms == zeros(size(lms)),2);
%     Lprobe = ~Lbldone & ~all(lms == zeros(size(lms)),2);
%     nspikes = zeros(size(stro.trial,1),1);
%     for i = 1:size(stro.trial,1)
%         spiketimes = stro.ras{i,spikeidx};
%         nspikes(i) = sum(spiketimes > stimon_t(i) & spiketimes < stimoff_t(i));
%     end
% 
%     uniquecoloridxs = unique(coloridxs);
%     uniquecoloridxs = uniquecoloridxs([3:5]);  % hack to avoid baseline and 2nd order measurements.
%     blnspikes = nspikes(Lbaseline&~Lprobe);
%     tmp =[];
%     for i = uniquecoloridxs'
%         L = coloridxs == i;
%         [u,s,v] = svd(lms(L,:));
%         unitvector = v(:,1);
%         if (sum(sign(unitvector)) < 0)
%             unitvector = -unitvector;
%         end
%         contrasts = lms(L,:)*unitvector;
%         if(sum(contrasts) < 0)
%             contrasts = -contrasts;
%         end
%         con = [zeros(sum(Lbaseline&~Lprobe),1); contrasts];
%         nsp = [blnspikes; nspikes(L)];
%         [baseline, thresh, slope, success] = contrastrespfit(con, nsp, type);
%         if (strcmp(type,'linear'))
%             fittedthresh = (stro.sum.exptParams.threshold-baseline)/slope+thresh;
%             fittedthresh = thresh;
%         else
%             fittedthresh = sqrt((stro.sum.exptParams.threshold-baseline)/slope)+thresh;
%                         fittedthresh = thresh;
%         end
%         tmp = [tmp; fittedthresh contrasts(end)];
%     end
%     data = cat(3,data,tmp);
% end
% 
% c = [1 2; 1 3; 2 3];
% figure;
% for i = 1:size(c,1)
%     fittedratios = data(c(i,1),1,:)./data(c(i,2),1,:);
%     endratios = data(c(i,1),2,:)./data(c(i,2),2,:);
%     subplot(2,2,i); hold on;
%     [n1,x] = hist([fittedratios endratios]);
%     [n2,x] = hist(endratios,x);
%     bar(x,n1);
%     hold on;
%     h = bar(x,n2);
%     set(h,'FaceColor','red');
%     title(num2str([std(fittedratios)/mean(fittedratios) std(endratios)/mean(fittedratios)]));
% end
% subplot(2,2,4); hold on;
% plot(squeeze(data(:,1,:)),'b.-');
% plot(squeeze(data(:,2,:)),'r.-');
% set(gca,'YScale','log');
%%
% Grouping points into two clusters (prior to fitting a linear model)
% First, finding a rotation such that the mean (ranked) elevation of the OOG points is
% minimal and mean ranked elevation of the NOOG points is maximal.
% This is good because it only depends on the angles, not on the
% thresholds.
% But it's bad because it discretizes points into OOG and NOOG and doesn't take into 
% % account the dimensions of the gamut.
% % It's a search over a two parameter space (azimuth and elevation).
% % Could be made better if we had a way of linearly transforming the data
% % such that the distribution of (3D) angles was close to uniform.
% 
% USEPCA = 0;
% 
% Loog = data(Laccept,6) == 1;
% azs = linspace(-pi,pi,60);
% elevs = linspace(-pi/2, pi/2, 60);
% oogmat = zeros(60,60);
% noogmat = zeros(60,60);
% tmp = data(Laccept,[1:3]);
% if (~USEPCA)
%     xformed = tmp;
% else
%     [coeffs_all, score, roots_all]= princomp([tmp; -tmp]);
%     coeffs_all = coeffs_all*diag(1./sqrt(roots_all));
%     xformed = tmp*coeffs_all;
% end
% 
% for i = 1:length(azs)
%     for j = 1:length(elevs)
%         az = azs(i);
%         elev = elevs(j);
%         rotmat1 = [cos(az) -sin(az) 0; sin(az) cos(az) 0; 0 0 1];
%         rotmat2 = [1 0 0; 0 cos(elev) -sin(elev); 0 sin(elev) cos(elev)];
%         tmp = rotmat1*[xformed; -xformed]';
%         out = (rotmat2*tmp)';
%         out = mkbasis(out')';
%         tmpoog = out(logical([Loog; Loog]),3);
%         tmpnoog = out(logical([~Loog; ~Loog]),3);
%         tmpoog(tmpoog < 0) = [];
%         tmpnoog(tmpnoog < 0) = [];
%         
%         tmp = tiedrank([tmpoog; tmpnoog]);
%         oogmat(i,j) = sum(tmp(1:length(tmpoog)));
%         noogmat(i,j) = sum(tmp(length(tmpoog)+1:end));
%         %     oogmat(i,j) = acos(mean(tmpoog)); % they're already unit vectors
%         %     noogmat(i,j) = acos(mean(tmpnoog)); % they're already unit vectors
%     end
% end
% ratio = (noogmat./oogmat);
% 
% figure;
% surf(azs, elevs, ratio);
% [i,j] = ind2sub(size(ratio),find(ratio == max(ratio(:)),1));
% bestaz = azs(i);
% bestelev = elevs(j);
% 
% % Now just giving it a try
% rotmat1 = [cos(bestaz) -sin(bestaz) 0; sin(bestaz) cos(bestaz) 0; 0 0 1];
% rotmat2 = [1 0 0; 0 cos(bestelev) -sin(bestelev); 0 sin(bestelev) cos(bestelev)];
% normalvect = inv(rotmat2*rotmat1)*[0 0 1]'; % could be in PCA space.
% if (USEPCA)
%     normalvect = inv(coeffs_all)*normalvect; % in LMS space
% end
% normalvect = normalvect./norm(normalvect);
% out = rotmat2*rotmat1*xformed';
% figure; axes; hold on;
% for i = 1:size(out,2)
%     projs = sign(out(3,i))*out(:,i); % Only looking at 1 hemisphere
%     az = atan2(projs(2), projs(1)); 
%     elev = acos(projs(3)./sqrt(sum(projs.^2)));
%     h = polar(az,elev,'ko');
%     set(h,'MarkerFaceColor',[1 1 1])
%     if (Loog(i))
%         set(h,'MarkerFaceColor',[1 0 0]);
%     end
% end
% x = linspace(0,2*pi,200);
% plot(pi/2*cos(x),pi/2*sin(x),'k-');
% axis square;
% title(['Looking in direction: ',num2str(normalvect',3)]);
% 
% Lgrp1 = logical(data(:,[1 2 3])*normalvect > 0);  % normalvect is in LMS space
% % Destructively modifying data so that all of the points are from 
% % a single cluster.  Flip the signs to get the other cluster.
% data(Lgrp1,[1 2 3]) = -data(Lgrp1,[1 2 3]);
% 
% % Plotting in 3-D space to see how good the clusters look
% tmp = data(Laccept,[1:3]) .*repmat(data(Laccept,4), 1,3);
% figure; axes; hold on;
% plot3(tmp(~Loog,1), tmp(~Loog,2), tmp(~Loog,3),'k.');
% plot3(-tmp(~Loog,1), -tmp(~Loog,2), -tmp(~Loog,3),'b.');
% plot3([-.1 .1],[0 0], [0 0],'k-');
% plot3([0 0], [-.1 .1], [0 0],'k-');
% plot3([0 0], [0 0],[-.1 .1],'k-');
% xlabel('L');
% ylabel('M');
% zlabel('S');
% plot3([zeros(sum(Loog),1) tmp(Loog,1)]', [zeros(sum(Loog),1) tmp(Loog,2)]', [zeros(sum(Loog),1) tmp(Loog,3)]','y-');
% plot3([zeros(sum(Loog),1) -tmp(Loog,1)]', [zeros(sum(Loog),1) -tmp(Loog,2)]', [zeros(sum(Loog),1) -tmp(Loog,3)]','y-');
% plot3([0 normalvect(1)],[0 normalvect(2)],[0 normalvect(3)],'r-')
% 

% %%
% % 3-D plotting in normalized post-receptoral axes (or cone coordinates)
% 
% % First, coming up with an initial guess
% PLOTMESH = 0;
% COORDS = 'CONES';
% COORDS = 'DKL';
% 
% figure; axes; hold on;
% scaled = data(Laccept,[1:3]) .*repmat(data(Laccept,4), 1,3);
% mn = mean(scaled(~Loog,:))';
% [pcs, score, d] = princomp(scaled(~Loog,:));
% initxyz = -sign(mn'*pcs(:,3))* pcs(:,3)./norm(mn'*pcs(:,3));
% 
% % errfun = @(lightvects, mechvect) sum((sqrt(sum(lightvects.^2,2))-(lightvects*mechvect)).^2);
% Loog = data(Laccept,6) == 1;
% if (strcmp('CONES',COORDS))
%     basis = [1 0 0; 0 1 0; 0 0 1];
%     xlabel('L');
%     ylabel('M');
%     zlabel('S');
% elseif (strcmp('DKL',COORDS))
%     basis = [1 -1 0; 0 0 1; 1 1 0];
%     xlabel('L-M');
%     ylabel('S');
%     zlabel('L+M');
% end
% 
% tmp = data(Laccept,[1:3]) .*repmat(data(Laccept,4), 1,3);
% xformed = tmp*basis;
% normfacts = sqrt(var([xformed(~Loog,:); -xformed(~Loog,:)]));
% xformed = xformed./repmat(normfacts, size(xformed,1),1);
% 
% plot3([-1 1],[0 0], [0 0],'k-');
% plot3([0 0], [-1 1], [0 0],'k-');
% plot3([0 0], [0 0],[-1 1],'k-');
% 
% plot3(xformed(~Loog,1),xformed(~Loog,2),xformed(~Loog,3),'k.')
% plot3(-xformed(~Loog,1),-xformed(~Loog,2),-xformed(~Loog,3),'k.')
% plot3(xformed(Loog,1),xformed(Loog,2),xformed(Loog,3),'y.')
% plot3(-xformed(Loog,1),-xformed(Loog,2),-xformed(Loog,3),'y.')
% plot3([zeros(sum(Loog),1) xformed(Loog,1)]', [zeros(sum(Loog),1) xformed(Loog,2)]', [zeros(sum(Loog),1) xformed(Loog,3)]','y-');
% plot3([zeros(sum(Loog),1) -xformed(Loog,1)]', [zeros(sum(Loog),1) -xformed(Loog,2)]', [zeros(sum(Loog),1) -xformed(Loog,3)]','y-');
% 
% set(gcf,'Name',stro.sum.fileName);
% 
% set(gca,'XLim',[-max(abs(xformed(~Loog,1))) max(abs(xformed(~Loog,1)))].*[1.2 1.2]);
% set(gca,'YLim',[-max(abs(xformed(~Loog,2))) max(abs(xformed(~Loog,2)))].*[1.2 1.2]);
% set(gca,'ZLim',[-max(abs(xformed(~Loog,3))) max(abs(xformed(~Loog,3)))].*[1.2 1.2]);
% 
% if (PLOTMESH)
%     T = delaunay3(xformed(~Loog,1),xformed(~Loog,2),xformed(~Loog,3));
%     htm = tetramesh(T,[xformed(~Loog,1),xformed(~Loog,2),xformed(~Loog,3)], ones(length(T),1),'FaceAlpha', .1);
%     set(htm,'EdgeAlpha',.01)
% end
% 
%%
% 3-D plotting in PCA space.
% 
% Loog = data(Laccept,6) == 1;
% [coeffs, score, roots]= princomp([tmp(~Loog,:); -tmp(~Loog,:)]);
% 
% tmp = data(Laccept,[1:3]) .*repmat(data(Laccept,4), 1,3);
% xformed = tmp*coeffs;
% xformed = xformed./repmat(sqrt(roots+eps)',size(xformed,1),1);
% 
% figure; axes; hold on;
% plot3([-1 1],[0 0], [0 0],'k-');
% plot3([0 0], [-1 1], [0 0],'k-');
% plot3([0 0], [0 0],[-1 1],'k-');
% 
% plot3(xformed(~Loog,1),xformed(~Loog,2),xformed(~Loog,3),'k.')
% plot3(-xformed(~Loog,1),-xformed(~Loog,2),-xformed(~Loog,3),'k.')
% plot3(xformed(Loog,1),xformed(Loog,2),xformed(Loog,3),'y.')
% plot3(-xformed(Loog,1),-xformed(Loog,2),-xformed(Loog,3),'y.')
% plot3([zeros(sum(Loog),1) xformed(Loog,1)]', [zeros(sum(Loog),1) xformed(Loog,2)]', [zeros(sum(Loog),1) xformed(Loog,3)]','y-');
% plot3([zeros(sum(Loog),1) -xformed(Loog,1)]', [zeros(sum(Loog),1) -xformed(Loog,2)]', [zeros(sum(Loog),1) -xformed(Loog,3)]','y-');
% 
% xlabel(['PC 1: ',num2str(coeffs(:,1)',2),':   ',num2str(roots(1)./sum(roots))]);
% ylabel(['PC 2: ',num2str(coeffs(:,2)',2),':   ',num2str(roots(2)./sum(roots))]);
% zlabel(['PC 3: ',num2str(coeffs(:,3)',2),':   ',num2str(roots(3)./sum(roots))]);
% set(gcf,'Name',stro.sum.fileName);

%%
% Just hacking around trying to find parents
% a = data(7,[1 2 3])
% az = atan2(a(:,3),a(:,2));
% elev = acos(a(:,1)./sqrt(sum(a.^2)));
% 
% az = atan2(data(:,3),data(:,2));
% elev = acos(data(:,1)./sqrt(sum(data(:,[1 2 3]).^2,2)));

initialdirs = data(data(:,5) == 1,[1 2 3]);
% if (~stro.sum.exptParams.monopolar)
%     initialdirs = [initialdirs; -initialdirs];
% end
% Dumb coding.  Just finding sets of inital dirs that work together
b = [];
for i = 1:size(initialdirs,1)
    for j = i:size(initialdirs,1)
        for k = j:size(initialdirs,1)
            a = [initialdirs(i,:)',initialdirs(j,:)',initialdirs(k,:)'];
            if (rank(a) == 3)
                b = [b; i j k 0]
            end
        end
    end
end
% Below doesn't work yet
% if (~stro.sum.exptParams.monopolar) % round 2 is special - always 4 directions
%     tmp = fullfact([2 2 2])*2-3;
%     Loog = logical(data([1 2 3],6));
%     for i = 1:size(tmp,1)
%         norms = data([1 2 3],4);
%         norms(Loog) = norms(Loog).*stro.sum.exptParams.oogscale;
%        
%         newdir = mean(data([1 2 3],[1 2 3]).*repmat(norms.*tmp(i,:)',1,3));
%         newdir = newdir./norm(newdir);
%         
%         idx = find(all(abs(data(:,[1 2 3])-repmat(newdir,size(data,1),1)) < 0.001,2));
%         if (~isempty(idx))
%             disp('got one');
%             b = [b; b(i,[1 2]) idx];
%             b = [b; b(i,[1 3]) idx];
%             b = [b; b(i,[2 3]) idx];
%         end
%     end
% end
% columns of b are:
% vertex 1, 2, 3 (indices into data), parent triangle (index into b)
alreadyfound = [];
done = 0;
while (done == 0)
    done = 1;
    for i = 1:size(b,1)
        Loog = logical(data(b(i,[1 2 3]),6));
        norms = data(b(i,[1 2 3]),4);
        norms(Loog) = norms(Loog).*stro.sum.exptParams.oogscale;
        
        newdir = mean(data(b(i,[1 2 3]),[1 2 3]).*repmat(norms,1,3));
        newdir = newdir./norm(newdir);
        
        idx = find(all(abs(data(:,[1 2 3])-repmat(newdir,size(data,1),1)) < 0.001,2));
        if (~isempty(idx))
            if (~any(idx == alreadyfound))
                disp('got one');
                done = 0;
                b = [b; b(i,[1 2]) idx i];
                b = [b; b(i,[1 3]) idx i];
                b = [b; b(i,[2 3]) idx i];
                alreadyfound = [alreadyfound idx];
            end
        end
    end
end

% Making sure we're not missing any points...
tmp=unique(b(:,[1 2 3]));
if (any(diff(tmp)) ~= 1)
    error('Missing a point on the interior of data');
end
if (tmp(end) ~= size(data,1))
    error('Missing a point at the end of data');
end

h = [];
for i = 1:size(b,1)
    vertices = data(b(i,[1 2 3]),[1 2 3]).*repmat(data(b(i,[1 2 3]),4),1,3);
    h(i) = patch(vertices(:,1),vertices(:,2),vertices(:,3),1);
    if all(data(b(i,[1 2 3]),6))
        set(h(i),'FaceAlpha',1,'EdgeAlpha',0.01);
    else
        set(h(i),'FaceAlpha',.2,'FaceColor',[1 .5 .5]);
    end
    if (b(i,4) > 0)
        set(h(b(i,4)),'Visible','off');
    end
    pause
end
%%
% Playing around with ellipsoid fitting

scaled = data(Laccept,[1:3]) .*repmat(data(Laccept,4), 1,3);
Loog = logical(data(Laccept,6));
xyz = scaled(~Loog,:);

% First fitting by LS (on coefficients)
D = [xyz(:,1) .* xyz(:,1),...
    xyz(:,2) .* xyz(:,2),...
    xyz(:,3) .* xyz(:,3),...
    2*xyz(:,1) .* xyz(:,2),...
    2*xyz(:,1) .* xyz(:,3),...
    2*xyz(:,2) .* xyz(:,3)];
v = (D' * D) \(D' * ones(size(xyz,1),1));
A = [v(1) v(4) v(5);...
    v(4) v(2) v(6);...
    v(5) v(6) v(3)];
[evecs, evals] = eig(A);
radii = sqrt(1./diag(evals));

 % draw data
figure; axes; hold on;
plot3(xyz(:,1),xyz(:,2),xyz(:,3),'r.' );
plot3(-xyz(:,1),-xyz(:,2),-xyz(:,3),'r.' );

%draw fit
%[x, y, z] = ellipsoid(0,0,0,radii(1),radii(2), radii(3));
%[newxyz] = [x(:) y(:) z(:)]*evecs'
%h = surf(reshape(newxyz(:,1),size(x)), reshape(newxyz(:,2),size(y)), reshape(newxyz(:,3),size(z)));
%set(h,'FaceAlpha',.2,'FaceColor','yellow','Edgealpha',0);
%axis vis3d;
%camlight;
%lighting phong;

% Now fitting by minimizing the sum of the differences between log(r) and log(predr)
% using LS as an initial guess
options = optimset('MaxFunEvals',50000,'MaxIter',50,'TolFun',10^-6,'TolX',10^-6);
[tmp, SSE1, exitflag1] = fminsearch(@(x) surfacefiterr3(scaled,x),v, options);
A = [tmp(1) tmp(4) tmp(5);...
    tmp(4) tmp(2) tmp(6);...
    tmp(5) tmp(6) tmp(3)];
[evecs, evals] = eig(A);
radii = sqrt(1./diag(evals));
[x, y, z] = ellipsoid(0,0,0,radii(1),radii(2), radii(3));
[newxyz] = [x(:) y(:) z(:)]*evecs'
h = surf(reshape(newxyz(:,1),size(x)), reshape(newxyz(:,2),size(y)), reshape(newxyz(:,3),size(z)));
set(h,'FaceAlpha',.2,'FaceColor','green','Edgealpha',0);
axis vis3d;
camlight;
lighting phong;


