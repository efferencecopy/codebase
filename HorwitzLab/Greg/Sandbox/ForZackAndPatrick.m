% For Zack

%%
% Copied from IsoSampAnalyses
stro = nex2stro;
trialFields = stro.sum.trialFields(1,:);
stimon_t = stro.trial(:,strcmp(trialFields, 'stimon_t'));
stimoff_t = stro.trial(:,strcmp(trialFields, 'stimoff_t'));
anlgstart_t = stro.ras(:,strcmp(stro.sum.rasterCells, 'anlgStartTime'));
spikes = stro.ras(:,strcmp(stro.sum.rasterCells, 'sig001a'));
sf = stro.trial(:,strcmp(trialFields, 'sf'));
num_stims = min(stro.sum.exptParams.num_stims, size(stro.trial, 1));

lms = [stro.trial(:,strcmp(trialFields, 'stim_l')) ...
       stro.trial(:,strcmp(trialFields, 'stim_m')) ...
       stro.trial(:,strcmp(trialFields, 'stim_s'))] ./ 100; % stro1,2 don't require this
   
[unqLms,~,unqIdxs] = unique(lms, 'rows');
   
binwidth = .005;
stim_dur = mode(stimoff_t - stimon_t);
offset = [-0.1 stim_dur + 0.1]; % 100 ms window
bins = offset(1):binwidth:offset(2);
PSTH = zeros(size(bins));
stim_resp = zeros(num_stims, 1);
num_trials = size(spikes, 1);
theseSpikes = cell(size(spikes));

for trialIdx = 1:num_trials
    tempSpikes = spikes{trialIdx} - stimon_t(trialIdx);
    tempSpikes = tempSpikes(tempSpikes >= offset(1) ...
                              & tempSpikes <= offset(2));
    theseSpikes{trialIdx} = tempSpikes;
end

figure(); axes(); hold on;
for stimIdx = 1:num_stims
    Lstims = unqIdxs == stimIdx;
    if sum(Lstims) > 0
        spikeTimes = cat(1, theseSpikes{Lstims});
        stim_resp(stimIdx) = numel(spikeTimes);
        plot([spikeTimes spikeTimes]', ...
             [zeros(stim_resp(stimIdx), 1) .5 * ones(stim_resp(stimIdx), 1)]' + stimIdx, 'k-');
        PSTH = PSTH + hist(spikeTimes, bins);
    end
end

set(gca, 'XLim', offset, 'YTick', [], 'YLim', [-4 num_stims + 5]);
plot([0 0], get(gca, 'YLim'), 'm-', ...
     [stim_dur stim_dur], get(gca, 'YLim'), 'm-');

PSTH = PSTH / binwidth / num_trials;

figure(); axes(); hold on;
plot(bins, PSTH, 'k-', 'linewidth', 2);
set(gca, 'YLim', [0 10 * ceil(max(PSTH/10))]);
set(gca, 'Xlim', offset);
plot([0 0], get(gca, 'YLim'), 'm-', ...
     [stim_dur stim_dur], get(gca, 'YLim'), 'm-');
 
figure(); axes(); hold on;
% scatter3(unqLms(:,1), unqLms(:,2), unqLms(:,3), 8, stim_resp/max(stim_resp), 'o', 'filled');
cmap = jet(64);
for i = 1:num_stims
    plot3(unqLms(i,1), unqLms(i,2), unqLms(i,3), 'o', 'markersize', 5, ...
        'markerfacecolor', cmap(ceil(63*stim_resp(i)/max(stim_resp))+1,:), ...
        'markeredgecolor', 'none');
end
xlabel('L'); ylabel('M'); zlabel('S');
%%
% Attempts to find a "preferred direction"
% [b,bint] = regress(stim_resp, [ones(num_stims, 1) abs(unqLms)]);
% plot3([0 b(2)/1e4], [0 b(3)/1e4], [0 b(4)/1e4], 'm-', 'linewidth', 3);
% plot3([0 -b(2)/1e4], [0 -b(3)/1e4], [0 -b(4)/1e4], 'm-', 'linewidth', 3);
maxidx = find(max(stim_resp) == stim_resp);
options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-6,'TolX',10^-6);
[tmpconeweights, tmpfval, exitflag] = fminsearch(@(x) linmodfiterr(unqLms, stim_resp, x), unqLms(maxidx,:),options);
%plot3([0 tmpconeweights(1)],[0 tmpconeweights(2)],[0 tmpconeweights(3)],'k-');

figure; axes; hold on;
projs = unqLms*tmpconeweights';
L = projs<0;
plot(abs(projs(L)),stim_resp(L),'ko');
plot(abs(projs(~L)),stim_resp(~L),'ro');
%%
% A shaded surface plot
k = convhulln(unqLms);
figure; axes; hold on;
patch('Vertices',unqLms,'Faces',k,'FaceVertexCData',repmat(stim_resp,1,3)/max(stim_resp(:)),'FaceColor','interp');
xlabel('L'); ylabel('M'); zlabel('S');
h = colorbar;
colormap(gray);
set(h,'YTick',[0 .5 1],'YTickLabel',[0 max(stim_resp(:))/2 max(stim_resp(:))]);
set(get(h,'YLabel'),'String','Response');
%%
% Plotting response as a function of projection onto some vector
% interactive. Point, click, hit return.
figure;

for i =1:100
    subplot(2,1,1);hold on;
    plot(cos(linspace(0,2*pi,100)),sin(linspace(0,2*pi,100)),'-')
    axis square;
    pause

    whichpt = get(gca,'CurrentPoint')
    whichpt = whichpt(1,[1 2]);
    if (norm(whichpt) > 1);
        whichpt = whichpt./norm(whichpt);
    end
    theta = acos(whichpt(1))
    phi = asin(whichpt(2))
    plot(whichpt(1),whichpt(2),'m*');
    [x,y,z] = sph2cart(theta,phi,1);
    vect = [x y z]'
    projs = unqLms*vect;
    
    subplot(2,1,2);
    plot(projs,stim_resp,'k.');
    title(num2str(vect'));
    set(gca,'XLim',[-.4 .4]);
end

%%
% A couple of analyses for Patrick's experiment

stro = nex2stro;
Lcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'Lcc'));
Mcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'Mcc'));
fpacq_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fpacq_t'));
stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fpoff_t'));
direction = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'Theta'));
ntrials = size(stro.trial,1);
LM = Lcc+Mcc;
LvM = (Lcc-Mcc);

bins = linspace(0,2,50);
PSTH = zeros(1,length(bins));
% First PSTHs to find good values for "offset"
for i = 1:ntrials
    spiketimes = stro.ras{i,1}-stimon_t(i);
    PSTH = PSTH +histc(spiketimes',bins);
end

% Getting baseline firing rate
dt = zeros(ntrials,1);
baselinerates = zeros(ntrials,1);
for i = 1:ntrials
    nsp =  sum(stro.ras{i,1} > fpacq_t(i)+.05 & stro.ras{i,1} < stimon_t(i));
    baselinerates(i) = nsp./(stimon_t(i) - fpacq_t(i) -.05);
    dt(i) = stimon_t(i) - fpacq_t(i);
end

figure; axes; hold on;
bar(bins,PSTH);
gauss = @(x,params)(params(2)+params(1)*exp((-(x-params(3)).^2)/params(4).^2));
bguess = mean(PSTH(bins<.25));
aguess = max(PSTH)-bguess;
muguess = bins(find(PSTH == max(PSTH),1,'first'));
sigmaguess = .1; % terrible. Need something better.
%plot(bins,gauss(bins,[aguess,bguess,muguess,sigmaguess]),'m.')
fittedparams = fminsearch(@(params)sum((gauss(bins,params)-PSTH).^2),[aguess,bguess,muguess,sigmaguess]);
plot(bins,gauss(bins,fittedparams),'b-.','LineWidth',3)
offset = [fittedparams(3)-fittedparams(4) fittedparams(3)+fittedparams(4)];
plot(offset,[fittedparams(1) fittedparams(1)],'k-','linewidth',3);

% Now counting up spikes in a window
ntrials = size(stro.trial,1);
nspikes = zeros(ntrials,1);
spikerates = zeros(ntrials,1);
for i = 1:ntrials
    nspikes(i) = sum(stro.ras{i,1} > stimon_t(i)+offset(1) & stro.ras{i,1}<stimon_t(i)+offset(2));
    dt = offset(2)-offset(1);
    spikerates(i) = nspikes(i)./dt;
end
%%
% Projecting onto a few unit vectors in the LM plane and 
% Fitting the 1-D function (response vs. projection) with 
% regular old 1-D cubic smoothing spline.
% nthetas = 1000;
% thetas = linspace(0,pi,nthetas);
% data = nan*ones(nthetas,1);
% figure; axes; hold on;
% smoothingparam = .9;
% for i = 1:length(thetas);
%     unitvect = [cos(thetas(i)) sin(thetas(i))]';
%     projs = [LM LvM]*unitvect;
%     lastwarn('');
%     pp  = csaps(projs,nspikes,smoothingparam);
%     plot(projs,nspikes,'k.');
%     fnplt(pp);
%     drawnow;
%     cla
%     data(i) = sum((fnval(pp,projs)-nspikes).^2);
%     if ~isempty(lastwarn)
%         data(i) = nan;
%     end
% end
% figure; plot(thetas,data);
% preftheta = thetas(data == min(data))
% % Taking a look at the residuals from a 1-D fit
% unitvect = [cos(preftheta) sin(preftheta)]';
% projs = [LM LvM]*unitvect;
% pp  = csaps(projs,nspikes,smoothingparam);
% predfr1D = fnval(pp,projs);
% figure; 
% plot(predfr1D,nspikes-predfr1D,'k.');
% xlabel('fitted value'); ylabel('residual');
% [x,y]= meshgrid(linspace(-1,1,50),linspace(-1,1,50));
% projs = [x(:) y(:)]*unitvect;
% predfr2D = fnval(pp,projs);
% figure; axes; hold on; axis square;
% surf(x,y,reshape(predfr2D,size(x)));
% plot3(LM,LvM,nspikes,'k.');
% LMconeweights = 1/sqrt(2)*[1 1;1 -1]*unitvect;
% title(['L weight = ',num2str(LMconeweights(1),1),' M weight = ',num2str(LMconeweights(2),1)]);
%%
% Trying a thinplate spline on the sqrt-transformed spikecounts.
% 
% figure; axes; hold on;
% 
% Ldir = logical(direction == 1);
% if sum(Ldir) == 0; Ldir = true(length(Ldir),1); end
% out2d = tpaps([LM(Ldir), LvM(Ldir)]', nspikes(Ldir)');
% plot3(LM(Ldir),LvM(Ldir),nspikes(Ldir),'k.');
% h = fnplt(out2d);
% surf(h{1},h{2},h{3})
% 
%out2d = tpaps([LM, LvM]', sqrt(nspikes'));
%plot3(LM,LvM,sqrt(nspikes),'k.');
%h = fnplt(out2d);
%surf(h{1},h{2},h{3}.^2)
%axis square;

%out2d_noxform = tpaps([LM, LvM]', nspikes');
%h = fnplt(out2d_noxform);
%h =surf(h{1},h{2},h{3},repmat(5,size(h{1})))
%set(h,'FaceAlpha',.7)

% The square root transformation looks like it's pretty subtle.

% Now some 1-D splines
%outx = tpaps([LM, repmat(mean(LvM),size(LvM,1),1)]', sqrt(nspikes)',1);
%outx = csaps(LM, sqrt(nspikes));
%outy = csaps(LvM, sqrt(nspikes));
%figure;
%subplot(2,1,1);
%fnplt(outx);
%subplot(2,1,2);
%fnplt(outy);
%%
% Bootstrapping test to determine whether the 2-D fit is actually helping at all.
% nboot = 2000;
% PARAMETRIC = 2; % 0=Nonparametric, 1=Poisson, 2=Gaussian
% unitvect = [cos(preftheta) sin(preftheta)]';
% projs = [LM LvM]*unitvect;
% %projs = [LM LvM]*[unitvect(2) unitvect(1)]';
% pp1D  = csaps(projs,nspikes,smoothingparam);
% predfr1D = fnval(pp1D,projs);
% resid = predfr1D-nspikes;
% SS1D = sum(resid.^2);
% SS2D = sum((fnval(out2d,[LM LvM]')-nspikes').^2);
% [u,a,b] = unique([LM LvM],'rows');
% SS = zeros(nboot,2);
% for i = 1:nboot
%     i
%     tmpdata = zeros(size(u,1),1);
%     for j = 1:max(b)    
%         L = b == j;
%         r = resid(L);
%         if (PARAMETRIC == 0)
%             tmpdata(L) = predfr1D(L) + r(unidrnd(length(r),length(r),1)); % Nonparametric
%         elseif (PARAMETRIC == 1)
%             tmpdata(L) = poissrnd(predfr1D(L)); % Parametric (Poisson)
%         else
%             tmpdata(L) = predfr1D(L) + normrnd(mean(r),std(r),size(r)); % Parametric (Gaussian)
%         end
%     end
%     pp  = csaps(projs,tmpdata,smoothingparam);
%     predfr = fnval(pp,projs);
%     SS(i,1) = sum((predfr-tmpdata).^2);
%     out = tpaps([LM, LvM]', tmpdata');
%     SS(i,2) = sum((fnval(out,[LM LvM]')-tmpdata').^2);
% end
% SSratio = SS(:,1)./SS(:,2);  % Large values mean that TP spline fits way better
% figure; axes; hold on;
% hist(SSratio);
% plot(SS1D./SS2D,0,'m*');
% p = sum((SS1D./SS2D)<SSratio)/nboot;
% title(['p = ',num2str(p)]);

%%
% Fitting individual Naka-Rushton contrast-response functions.
% This only works for cells that were tested on a radial grid .
% Treating 0 indepndently from 180 (so, unipolar contrast-response functions)
% Fitting spike *rates* with Poisson error. Is this OK?

thetas = atan2(LvM, LM);
bins = linspace(-pi,pi,1000); % way finer than we expect
n = histc(thetas,bins);
uniquethetas = bins(n>(max(n)/2));
wedgewidth = min(diff(uniquethetas));
fittedparams = [];
figure; axes; hold on;
options = optimset;
options = optimset(options,'Diagnostics','off','Display','on');
options = optimset(options,'LargeScale','off','MaxFunEvals',2000);
vlb = [0  0.001 0.001 0];
vub = [1000 1000 100 100];
for i = 1:length(uniquethetas)
    unitvect = [cos(uniquethetas(i)) sin(uniquethetas(i))]';
    projs = [LM LvM]*unitvect;
    Ltheta = thetas > uniquethetas(i)-wedgewidth/2 &...
             thetas < uniquethetas(i)+wedgewidth/2;
    
    Ldir = direction == 1;
    if sum(Ldir) == 0; Ldir = true(length(Ldir),1); end
    %plot(projs(Ltheta&Ldir),spikerates(Ltheta&Ldir),'k.')
    %plot(projs(Ltheta&~Ldir),spikerates(Ltheta&~Ldir),'r.')
   
    currentnspikes = spikerates(Ltheta&Ldir);
    currentprojs = projs(Ltheta&Ldir);
    currentnspikes = [currentnspikes; baselinerates];
    currentprojs = [currentprojs; zeros(length(baselinerates),1)];
    plot(currentprojs,currentnspikes,'r.')
    
    topfr = mean(currentnspikes(currentprojs == max(currentprojs)));
    bottomfr = mean(baselinerates);
    % Trying to get a decent estimate of sigma
    %tmpspline = csaps(currentprojs,currentnspikes,.999);
    %fnplt(tmpspline)
   % sigmaguess = fnval(tmpspline,max(currentnspikes))-fnval(tmpspline,min(currentnspikes));
     sigmaguess = 0.1;
 %    plot(sigmaguess, max(currentnspikes)-min(currentnspikes),'m*');
    params0 = [topfr, sigmaguess, 2, bottomfr];  % need to constrain B across color directions
    f = fmincon('FitNakaRushtonFunJPW',params0,[],[],[],[],vlb,vub,[],options,currentprojs,currentnspikes); 
    
    fittedparams(i,:) = f;
    x = linspace(0,max(currentprojs),100);
    y = ComputeNakaRushtonJPW(f,x);
    plot(x,y,'b-');
    set (gca,'YLim',[0 max(spikerates)],'Xlim',[0 max(currentprojs)]);
    title(num2str(uniquethetas(i)));
    drawnow; 
    title(num2str(unitvect));
    pause;
    cla;
end

figure;
paramnames = {'Amp','sigma','exp'};
for i = 1:size(fittedparams,2)
    subplot(3,1,i);
    polar(uniquethetas(1:end),fittedparams(:,i)','k.-')
    title(paramnames{i});
end

figure; axes; hold on;
plot(fittedparams(:,1),fittedparams(:,2),'k.');
xlabel('Amp'); ylabel('sigma');
% Positive correlation between amp and sigma

%%
% As above but now projecting all the data onto each direction
% Starts with L-M

thetas = atan2(LvM, LM);
bins = linspace(-pi/2,pi/2,1000); % way finer than we expect
%[x,y] = pol2cart(thetas,ones(length(thetas),1)); % endpoints of unit vectors
n = histc(thetas,bins);
uniquethetas = bins(n>(max(n)/2));
wedgewidth = min(diff(uniquethetas));
fittedparams = [];
figure; axes; hold on;
options = optimset;
options = optimset(options,'Diagnostics','off','Display','on');
options = optimset(options,'LargeScale','off','MaxFunEvals',2000);
vlb = [0  0 0.001 0.001 1 1 0];
vub = [1000 1000 1 1 3 3 100];
errs = [];
for i = 1:length(uniquethetas)
    i
    unitvect = [cos(uniquethetas(i)) sin(uniquethetas(i))]';
    projs = [LM LvM]*unitvect;
    Ldir = direction >= 1;
    if sum(Ldir) == 0; Ldir = true(length(Ldir),1); end
    currentnspikes = spikerates(Ldir);
    currentprojs = projs(Ldir);
    plot(currentprojs,currentnspikes,'k.')
    
    topfr = max(spikerates);
    bottomfr = mean(baselinerates);
    
    sigmaguess = 0.1;
    params0 = [topfr, topfr, sigmaguess, sigmaguess, 2, 2, bottomfr];  % need to constrain B across color directions
    [f,errs(i)] = fmincon('FitNakaRushtonFunJPW',params0,[],[],[],[],vlb,vub,[],options,currentprojs,currentnspikes,'asymmetric');
    
    fittedparams(i,:) = f;
    x = linspace(min(currentprojs),max(currentprojs),100);
    y = ComputeNakaRushtonJPW(f,x,'asymmetric');
    plot(x,y,'b-');
    set (gca,'YLim',[0 max(spikerates)],'Xlim',[min(currentprojs) max(currentprojs)]);
    title(num2str(uniquethetas(i)));
    drawnow; 
    title(num2str(unitvect));
    pause; 
    cla;
end

figure;
plot(uniquethetas*180/pi, errs,'bo-')
set(gca,'Xtick',round(uniquethetas*180/pi));
ylabel('-llik');
xlabel('color dir (deg)');

%figure;
%paramnames = {'Amp','sigma','exp'};
%for i = 1:size(fittedparams,2)
%    subplot(3,1,i);
%    polar(uniquethetas(1:end),fittedparams(:,i)','k.-')
%    title(paramnames{i});
%end

%%
% For a trench cell, trying to predict points off the 
% cardinal axes from contrast-response functions fit to data on the
% cardinal axes.

% Setting stuff up
vlb = [0    0 0.001 0.001 1 1 0];
vub = [200 200    1     1 3 3 100];
options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-6,'TolX',10^-6);

bottomfr = 0;
topfr = max(nspikes);

sigmaguess = 0.01;
Ldir = logical(direction == 1);
if sum(Ldir) == 0; Ldir = true(length(Ldir),1); end
maxpred = []; f = [];
titles = {'L+M','L-M'};
for i = 1:2
    if (i == 1)
        L = Ldir & LvM == 0;
        projs = LM(L);
    else
        L = Ldir & LM == 0;
        projs = LvM(L);
    end
    
    topfr = max(nspikes(L));
    params0 = [topfr, topfr, sigmaguess, sigmaguess, 2, 2, bottomfr];  % need to constrain B across color directions
    f(i,:) = fmincon('FitNakaRushtonFunJPW',params0,[],[],[],[],vlb,vub,[],options,projs,nspikes(L),'asymmetric');
    
    figure; axes; hold on;
    plot(projs,nspikes(L),'ko');
    x = linspace(min(projs),max(projs),100);
    y = ComputeNakaRushtonJPW(f(i,:),x,'asymmetric');
    plot(x,y,'b-');
    maxpred(i) = max(y);
    title(titles{i});
end

% Making predictions on the basis of the single direction with more "action".
i = find(maxpred == max(maxpred));

% making predictions
coords = [LM(~L&Ldir) LvM(~L&Ldir)];
projs = coords(:,i);
pred = ComputeNakaRushtonJPW(f(i,:),projs,'asymmetric');
figure; subplot(2,1,1); hold on;
plot(pred,nspikes(~L&Ldir),'.');
plot([0 max(pred)],[0 max(pred)],'k-')
xlabel('pred'); ylabel('actual response');
subplot(2,1,2); hold on;
resid = nspikes(~L&Ldir)-pred;
plot(pred,resid,'.');
plot([min(pred) max(pred)],[0 0],'k-')
xlabel('pred'); ylabel('residual');

% residuals in the plane
figure; hax(1) = subplot(1,2,1); hold on;
plot3(LM(~L&Ldir), LvM(~L&Ldir), nspikes(~L&Ldir)-pred,'k.')
out2d = tpaps([LM(~L&Ldir), LvM(~L&Ldir)]', nspikes(~L&Ldir)'-pred');
h = fnplt(out2d);
surf(h{1},h{2},h{3});
f = sum(pred)-(log(pred')*nspikes(~L&Ldir));
axis tight; axis square; title(['residuals: ',num2str(round(f))]);

hax(2) = subplot(1,2,2); hold on;
plot3(LM(Ldir), LvM(Ldir), nspikes(Ldir),'r.')
out2d = tpaps([LM(Ldir), LvM(Ldir)]', nspikes(Ldir)');
h = fnplt(out2d);
surf(h{1},h{2},h{3});
axis tight; axis square; title('TPS to raw data');

xlims = [get(hax(1),'Xlim'); get(hax(2),'Xlim')];
set(hax,'XLim',[max(xlims(:,1)), min(xlims(:,2))]);
ylims = [get(hax(1),'Ylim'); get(hax(2),'Ylim')];
set(hax,'YLim',[max(ylims(:,1)), min(ylims(:,2))]);
zlims = [get(hax(1),'Zlim'); get(hax(2),'Zlim')];
set(hax,'ZLim',[min(zlims(:,1)), max(zlims(:,2))]);

%%
% Generating an initial guess for the cell below (asymmetric)
% This doesn't work very well yet.

vlb = [0    0 0.001 0.001 1 1 0];
vub = [200 200    1     1 6 6 100];
options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-6,'TolX',10^-6);

sigmaguess = 0.1;
Ldir = logical(direction == 1);
if sum(Ldir) == 0; Ldir = true(length(Ldir),1); end
maxpred = []; f = [];
titles = {'L+M','L-M'};
for i = 1:2
    if (i == 1)
        L = Ldir & LvM == 0;
        projs = LM(L);
    else
        L = Ldir & LM == 0;
        projs = LvM(L);
    end
    
    topfr = max(nspikes(L));
    if (topfr == 0)
        f(i,:) = [0 0 1 1 2 2 0]
    else
        params0 = [topfr, topfr, sigmaguess, sigmaguess, 2, 2, 0];  % need to constrain B across color directions
        f(i,:) = fmincon('FitNakaRushtonFunJPW',params0,[],[],[],[],vlb,vub,[],options,projs,nspikes(L),'asymmetric');
    end
    pred = ComputeNakaRushtonJPW(f(i,:),projs,'asymmetric');
    subplot(2,1,i); hold on;
    plot(pred,nspikes(L),'.');
    plot([0 max(pred)],[0 max(pred)],'k-');
    xlabel('pred'); ylabel('actual response');
end
a = mean(f(1,[1 2]));
b = mean(f(2,[1 2]));
%params0(1) = max(max(f([1 2],[1 2])));
params0 = max(nspikes);
params0(2) = a./sum([a b]); %this doesn't work so well
params0(2) = .5;
params0(3) = .25;
params0(4) = 2; 
params0(5) = (f(1,2)./f(1,4))./(f(1,1)./f(1,3)); if (isnan(params0(5))) params0(5) = 1; end;
params0(6) = (f(2,2)./f(2,4))./(f(2,1)./f(2,3)); if (isnan(params0(6))) params0(6) = 1; end;%
%params0(5) = 1;
%params0(6) = 1;
params0(7) = mean(f(:,end))

%%
% Trying to fit a 2-D Naka-Rushton to a GridLMPlane dataset.
% Need to modify model to give different Rmaxes in different directions.
% Next project: Need to come up with good initial guesses

params0 = [];
params0(1) = max(nspikes); % A
params0(2) = .5; % w
params0(3) = .5; % c50
params0(4) = 2; % exponent
params0(5) = 1; % k1 (asymmetry in dir 1)
params0(6) = 1; % k2 (asymmetry in dir 2)
params0(7) = 0; % baseline

%%
% Fitting the "surface1" model to all the data, or just to the data along
% the L+M and L-M axes.

vlb = [10   0 0.0001 2  0  0  0];
vub = [1000 0     1  2 10 10 20];
options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-6,'TolX',10^-6,'Display','on');

Ldir = logical(direction == 1);
if sum(Ldir) == 0; Ldir = true(length(Ldir),1); end

f0 = fmincon('FitNakaRushtonFunJPW',params0,[],[],[],[],vlb,vub,[],options,[LM(Ldir) LvM(Ldir)],nspikes(Ldir),'surface1');
[x,y] = meshgrid(linspace(min(LM),max(LM),50), linspace(min(LvM),max(LvM),50));
surface = ComputeNakaRushtonJPW(f0,[x(:) y(:)],'surface1');
figure; axes; hold on;
surf(x,y,reshape(surface,size(x)))
plot3(LM(Ldir), LvM(Ldir), nspikes(Ldir),'ko');
title(['Fit to all the data (n = ',num2str(sum(Ldir)),')']); 

% Let's try fitting to the L+M and L-M directions only
% Minor cheat - using result from one analysis to feed the other
L = LM == 0 | LvM == 0;
f1 = fmincon('FitNakaRushtonFunJPW',f0,[],[],[],[],vlb,vub,[],options,[LM(L&Ldir) LvM(L&Ldir)],nspikes(L&Ldir),'surface1');
[x,y] = meshgrid(linspace(min(LM),max(LM),50), linspace(min(LvM),max(LvM),50));
surface = ComputeNakaRushtonJPW(f1,[x(:) y(:)],'surface1');
figure; axes; hold on;
surf(x,y,reshape(surface,size(x)))
plot3(LM(Ldir), LvM(Ldir), nspikes(Ldir),'ko');
title(['Fit to L+M and L-M only (n = ',num2str(sum(L&Ldir)),')']);

% Looking at residuals, first for one model, then the other.
% Predicting data points (~L) because we're looking at the off-cardinal data points
for i = 1:2
    if i == 1
        pred = ComputeNakaRushtonJPW(f0,[LM(~L&Ldir) LvM(~L&Ldir)],'surface1');
    else
        pred = ComputeNakaRushtonJPW(f1,[LM(~L&Ldir) LvM(~L&Ldir)],'surface1');
    end
    
    % making predictions
    figure; subplot(2,1,1); hold on;
    plot(pred,nspikes(~L&Ldir),'.');
    plot([0 max(pred)],[0 max(pred)],'k-')
    xlabel('pred'); ylabel('actual response');
    subplot(2,1,2); hold on;
    resid = nspikes(~L&Ldir)-pred;
    plot(pred,resid,'.');
    plot([min(pred) max(pred)],[0 0],'k-')
    xlabel('pred'); ylabel('residual');
    
    % residuals in the plane
    figure; hax(1) = subplot(1,2,1); hold on;
    plot3(LM(~L&Ldir), LvM(~L&Ldir), nspikes(~L&Ldir)-pred,'k.')
    out2d = tpaps([LM(~L&Ldir), LvM(~L&Ldir)]', nspikes(~L&Ldir)'-pred');
    h = fnplt(out2d);
    surf(h{1},h{2},h{3});
    f = sum(pred)-(log(pred')*nspikes(~L&Ldir));
    axis tight; axis square; title(['residuals: ',num2str(round(f))]);
    
    hax(2) = subplot(1,2,2); hold on;
    plot3(LM(Ldir), LvM(Ldir), nspikes(Ldir),'r.')
    out2d = tpaps([LM(Ldir), LvM(Ldir)]', nspikes(Ldir)');
    h = fnplt(out2d);
    surf(h{1},h{2},h{3});
    axis tight; axis square; title('TPS to raw data');
    
    xlims = [get(hax(1),'Xlim'); get(hax(2),'Xlim')];
    set(hax,'XLim',[max(xlims(:,1)), min(xlims(:,2))]);
    ylims = [get(hax(1),'Ylim'); get(hax(2),'Ylim')];
    set(hax,'YLim',[max(ylims(:,1)), min(ylims(:,2))]);
    zlims = [get(hax(1),'Zlim'); get(hax(2),'Zlim')];
    set(hax,'ZLim',[min(zlims(:,1)), max(zlims(:,2))]);
end
%%
% 

%%
% ---------------------------------------
% ---------------------------------------
% Everything below this point is obsolete
% ---------------------------------------
% ---------------------------------------
%
%%
% Projecting onto a few unit vectors in the LM plane and 
% Fitting the 1-D function (response vs. projection) with a sum
% of weighted Gaussians.
nthetas = 100;
thetas = linspace(0,pi,nthetas);
data = nan*ones(nthetas,1);
nknots = 20;  sigmascale = nknots/70;
X = zeros(length(LM),1);
for i = 1:length(thetas);
    unitvect = [cos(thetas(i)) sin(thetas(i))]';
    projs = [LM LvM]*unitvect;
 %   knots = linspace(min(projs),max(projs),nknots);
    knots = linspace(.5*min(projs),2*max(projs),nknots);
    sigma = range(projs)*sigmascale;
    X = exp((-(repmat(projs,1,size(knots,2))-repmat(knots,size(projs,1),1)).^2)./sigma.^2);
    %for j = 1:nknots-4
    %    sp = spmak(knots(j:j+4),1);
    %    X(:,j) = spval(sp,projs); 
    %end
    [b,dev,stats] = glmfit(X,nspikes,'poisson','link','identity','constant','off');
    data(i) = dev;
end
figure; subplot(2,1,1); hold on;
plot(thetas,data,'.-')
set(gca,'YScale','log')
subplot(2,1,2); hold on;
unitvect = [cos(thetas(data == min(data))) sin(thetas(data == min(data)))]';
projs = [LM LvM]*unitvect;
sigma = range(projs)*sigmascale;
%knots = linspace(min(projs),max(projs),nknots);
knots = linspace(.5*min(projs),2*max(projs),nknots);
%clear y
%for j = 1:nknots-4
%    sp = spmak(knots(j:j+4),1);
%    y(:,j) = spval(sp,projs);
%end
X = exp((-(repmat(projs,1,size(knots,2))-repmat(knots,size(projs,1),1)).^2)./sigma.^2);
[b,dev,stats] = glmfit(X,nspikes,'poisson','link','identity','constant','off');
x = linspace(min(projs),max(projs),200)';
plot(projs,nspikes,'k.');
%clear y
%for j = 1:nknots-4
%    sp = spmak(knots(j:j+4),1);
%    y(:,j) = spval(sp,x);
%end
y = exp((-(repmat(x,1,size(knots,2))-repmat(knots,size(x,1),1)).^2)./sigma.^2);
plot(x,y*b,'y-','LineWidth',3);
title(num2str(unitvect));

% Now seeing whether there is any structure at all in the orthogonal axis
orthvect = [unitvect(2) -unitvect(1)]';
projs = [LM LvM]*orthvect;
figure; axes; hold on;
plot(projs,nspikes,'k.');
[b,dev,stats] = glmfit(projs,nspikes,'poisson','link','identity');
title(['p = ',num2str(stats.p(2))])
plot([min(projs) max(projs)],b(2)*[min(projs) max(projs)]+b(1),'k-')
% Hmm.. still some structure.

% Are there any projections with no linear structure? 
% This analysis is sketchy because many of the fits do not converge
ps = zeros(length(thetas),1);
for i = 1:length(thetas);
    unitvect = [cos(thetas(i)) sin(thetas(i))]';
    projs = [LM LvM]*unitvect;
    [b,dev,stats] = glmfit(projs,nspikes,'poisson','link','identity');
    ps(i) = stats.p(2);
end
orthvect = [cos(thetas(ps == max(ps))) sin(thetas(ps == max(ps)))]';
projs = [LM LvM]*orthvect;
figure; axes; hold on;
plot(projs,nspikes,'k.');
[b,dev,stats] = glmfit(projs,nspikes,'poisson','link','identity');
title(['p = ',num2str(stats.p(2))])
plot([min(projs) max(projs)],b(2)*[min(projs) max(projs)]+b(1),'k-')

figure; axes; hold on;
plot3(LM, LvM, nspikes,'k.');
plot([0 unitvect(1)],[0 unitvect(2)],'m-','LineWidth',2)
plot([0 orthvect(1)],[0 orthvect(2)],'g-','LineWidth',2)

% Projecting onto the L-M axis
%projs = [LM LvM]*orthvect;
%figure; axes; hold on;
%plot(projs,nspikes,'k.');
%[b,dev,stats] = glmfit(projs,nspikes,'poisson','link','identity');
%stats.p
%plot([min(projs) max(projs)],b(2)*[min(projs) max(projs)]+b(1),'k-')

% How good a job do we do representing the surface with a single direction?

%%
% How about just fitting a glm using L+M and L-M?
X = [LM, LM.^2, LvM, LvM.^2, LM.*LvM];
%X = [Lcc, Lcc.^2, Mcc, Mcc.^2, Lcc.*Mcc];

[b,dev,stats] = glmfit(X,nspikes,'poisson','link','identity','constant','on');
stats.t


%%
% Polynomial regression with a bunch of rotations of the data
% Looking for a projection where including the orthogonal axis
% with linear and quadratic terms does not increase significance.
thetas = linspace(-pi/4,pi/4,100);
data = [];
for i = 1:length(thetas);
    unitvect = [cos(thetas(i)) sin(thetas(i))]';
    orthvect = [unitvect(2) -unitvect(1)]';
    projs1 = [LM LvM]*unitvect;
    projs2 = [LM LvM]*orthvect;
    [b1,dev1,stats1] = glmfit([projs1 projs1.^2],nspikes,'poisson','link','identity');
    [b2,dev2,stats2] = glmfit([projs2 projs2.^2],nspikes,'poisson','link','identity');
    [b3,dev3,stats3] = glmfit([projs1 projs1.^2 projs2 projs2.^2 projs1.*projs2],nspikes,'poisson','link','identity');
    overall_p =  1-chi2cdf((dev1-dev3),3);
    data = [data; b1' b2' b3' overall_p];
end

% Trying to find a projection with little structure
idx = find(data(:,end) == max(data(:,end)));
unitvect = [cos(thetas(idx)) sin(thetas(idx))]'; % proj onto this vector has *least* structure
orthvect = [unitvect(2) -unitvect(1)]';
projs1 = [LM LvM]*orthvect;
projs2 = [LM LvM]*unitvect;
[x,y] = meshgrid(linspace(-1,1,30),linspace(-1,1,30));
x = x(:); y = y(:);
z1 = [ones(size(x,1),1) x x.^2 y y.^2 x.*y]*data(idx,[7:12])';
z2 = [ones(size(x,1),1) x x.^2]*data(idx,[1:3])';
figure; axes; hold on;
surf(reshape(x,30,30),reshape(y,30,30),reshape(z1,30,30));
plot3(projs2,projs1,nspikes,'k.');
title(['p = ',num2str(data(idx,end))]);
axis square;

%% 
% Looking for non-stationarity
figure; axes; hold on;
plot(nspikes,'k.')
[b,dev,stats] = glmfit([1:length(nspikes)],nspikes,'poisson','link','identity');
plot([1 length(nspikes)],[1 length(nspikes)]*b(2)+b(1),'m-','LineWidth',5);
title(['p=',num2str(stats.p(2))]);