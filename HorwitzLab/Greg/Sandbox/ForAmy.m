% A potential background

x = linspace(0,10,500);
f = logspace(-1,.6,500);
a = linspace(1,.3,500);
a = linspace(1,1,500);
y = a.*cos(x.*f);
%plot(x,y);

imagesc(y'*y);
colormap(gray);
axis image; axis xy;
set(gca,'XTick',[],'YTick',[]);

%%
% Another potential background (a radial grating)

freq = 18;
DISKRAD = 8000;
DISKPOS1 = [ 781  500]
DISKPOS2 = [ 400  200]
DISKCOL = .6;
BKAMP = 1;


figure; axes;
set(gcf,'Position',[270 44 1171 740]);

set(gca,'Units','pixels');
pos = get(gca,'OuterPosition')

[x,y] = meshgrid([1:pos(3)],[1:pos(4)]);
centerpoint = [pos(3) pos(4)]*2/3;

theta = atan2(y-centerpoint(2),x-centerpoint(1));
colormap(gray(255))
im = BKAMP*cos(freq*theta)/2+.5;

disks = ((x-DISKPOS1(1)).^2+(y-DISKPOS1(2)).^2 < DISKRAD |...
        (x-DISKPOS2(1)).^2+(y-DISKPOS2(2)).^2 < DISKRAD);
im = im.*(1-disks);
    
image(255*max(cat(3,DISKCOL*disks,im),[],3))
axis xy;

%%
% Analysis of gratings data
GT = nex2stro;
orients = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'orient'));
sfs = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'sf'));
unq_sfs = unique(sfs');
unq_orients = unique(orients');
stimon_t = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'stim_on'));
stimoff_t= GT.trial(:,strcmp(GT.sum.trialFields(1,:),'stim_off'));
spikeidxs = strncmp(GT.sum.rasterCells,'sig',3);
hs = [];
figure;
for spikeidx = find(spikeidxs)
    spikerates = [];
    baselines = [];  baseline_t = 0.5;
    for i = 1:size(GT.trial,1)
        spiketimes = GT.ras{i,spikeidx};
        nspikes = sum(spiketimes > stimon_t(i) & spiketimes < stimoff_t(i));
        spikerates = [spikerates; nspikes./(stimoff_t(i)-stimon_t(i))];
        nspikes = sum(spiketimes > stimon_t(i)-baseline_t & spiketimes < stimon_t(i));
        baselines = [baselines; nspikes./baseline_t];
    end
    
    frmap = zeros(length(unq_orients),length(unq_sfs));
    semap = zeros(length(unq_orients),length(unq_sfs));
    for j = 1:length(unq_sfs)
        for k = 1:length(unq_orients)
            L = sfs == unq_sfs(j) & orients == unq_orients(k);
            frmap(k,j) = mean(spikerates(L))
            semap(k,j) = std(spikerates(L))./sqrt(sum(L));
        end
    end
    %subplot(1,sum(spikeidxs),spikeidx); hold on;
    figure; axes; hold on;
    plot(unq_orients*180/pi,frmap,'Linewidth',2);
    plot(unq_orients*180/pi,frmap+semap,':');
    plot(unq_orients*180/pi,frmap-semap,':');
    legend(num2str(unq_sfs'));
    xlabel('orientation (deg)');
    ylabel('firing rate (sp/s)');
    title(GT.sum.rasterCells{spikeidx});
    set(gca,'Xlim',[0 315]);
end

%%
% What magnitude of artifact do we expect from assuming a fixed
% proportion between cm on the computer screen and degrees of visual angle?
% Can thins explain Why Amy needed to adjust the monitor for each cell to
% get the size-tuning curves (on the gray background) to be aligned?

% Bottom line - this effect is very small. I don't think that this could be
% the entire explanation for the discrepansies between the size-tuning
% curves measured on the greay background.

% Shortest distance from eye to screen
VIEWINGDIST = 50; % cm 
% How far the eye has to rotate to obtain an eccentric fixation point
ECCPOS = 20; % deg

SIZEATFP = tan(1*pi/180)*VIEWINGDIST; % distance on screen (in cm) subtended by 1 deg straight ahead
SIZEATECCPOS = tan((ECCPOS+1)*pi/180)*VIEWINGDIST - tan((ECCPOS)*pi/180)*VIEWINGDIST; % distance on screen (in cm) subtended by 1 deg off to the side

% Difference in screen distance travelled by a 1 deg rotation
SIZEATECCPOS-SIZEATFP  % in cm

% Now for a fixed distance on the screen, what is the difference in angle?
DVAATFP = atan2(1,VIEWINGDIST)*180/pi
% Where (in cm) on the monitor is the stimulus that is ECCPOS dva from the
% fp?
ECCPOSCM = tan(ECCPOS*pi/180)*VIEWINGDIST;
DVAATECCPOS = (atan2(ECCPOSCM+1, VIEWINGDIST)-atan2(ECCPOSCM, VIEWINGDIST))*180/pi;
DVAATFP-DVAATECCPOS
% This is the difference (in DVA) betwen the angular subtense of a 1 cm
% stimulus at the FP and at the eccentric location

%%
% Eye movement analysis testing ground
stro = nex2stro(findfile('A081713006.nex'));
offset = [0 .1]; % in sec
hepidx = strcmp(stro.sum.rasterCells(1,:),'AD11');
vepidx = strcmp(stro.sum.rasterCells(1,:),'AD12');
starttimeidx = strcmp(stro.sum.rasterCells(1,:),'anlgStartTime');
stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_on'));
stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_off'));
fpx = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fp_x'));
fpy = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fp_y'));
uniquefppos = unique([fpx fpy],'rows');
bkgnd = unique(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'bkgnd')));
samplerate = stro.sum.analog.storeRates{1};

ntrials = size(stro.ras,1);
% figure; axes; hold on;
data = nan*ones(ntrials,ceil(max(stimoff_t+offset(2)-stimon_t-offset(1))*samplerate),2);  % Eye position snippet
for i = 1:ntrials
    h = stro.ras{i,hepidx}.*4096/400;
    v = stro.ras{i,vepidx}.*4096/400;
    
    % 4096 A/D levels = 10 V
    % 1 degree = 40 A/D levels (according to REX)
    % (4096 levels/10 V) * (1 degree/40 levels) = 4096/400.
    % The key is that both REX and PLEXON use 12 bit A/D boards configured
    % for +/- 5 V.
    
   	e1t = stro.ras{i,starttimeidx};
    neyesamp = size(h,1);
    x = linspace(e1t,neyesamp/samplerate+e1t,neyesamp);
    L = x >=stimon_t(i)+offset(1) & x <=stimoff_t(i)+offset(2);
 %   plot(x(L)-stimon_t(i),h(L),'r-');
 %   plot(x(L)-stimon_t(i),v(L),'g-');
    data(i,1:sum(L),1) = h(L);
    data(i,1:sum(L),2) = v(L);
end

figure;
x = linspace(offset(1), offset(1)+2*(size(data,2)-1), size(data,2));
for i = 1:2
    for j = 1:2
        subplot(2,2,j+2*(i-1)); hold on;
        L = fpx == uniquefppos(i,1) & fpy == uniquefppos(i,2);
        plot(x,squeeze(data(L,:,j))','r-');
        plot(x,mean(squeeze(data(L,:,j))),'k-','LineWidth',2);
        if (j == 1)
            title(['Horz ',num2str(bkgnd)]);
        else
            title(['Vert ',num2str(bkgnd)]);
        end
        set(gca,'Xlim',[x(1) x(end)]);
        % Lnotnan = ~any(isnan(squeeze(data(L,:,1))));
        % [coeff, score, latent] = princomp(squeeze(data(L,Lnotnan,1)));
        % plot(coeff(:,1));
    end
end

%%
% Collecting a few files and seeing if I can eke out any difference in eye
% position between the corridor and the gray backgrounds.

filelist = {'F092412001.nex','F092412002.nex','F092412003.nex','F092412004.nex','F092412005.nex'};  % Size tuning shift, standard geometry
% filelist = {'F111412001.nex','F111412002.nex','F111412003.nex','F111412004.nex','F111412005.nex','F111412006.nex'}; % Size tuning shift, "eye movement control"
%filelist = {'A080313003.nex','A080313004.nex','A080313006.nex','A080313007.nex','A080313008.nex'};  % Size tuning shift, standard geometry
%filelist = {'A081113003.nex','A081113004.nex','A081113005.nex','A081113006.nex','A081113007.nex','A081113008.nex','A081113009','A081113010','A081113011'};  % Size tuning shift, "eye movement control"
%filelist = {'A082213001.nex','A082213002.nex','A082213003.nex','A082213004.nex','A082213005.nex','A082213006.nex','A082213007.nex'}; % standard geometry - shift in EP but not in tuning curve?
%filelist = {'A082113001.nex','A082113002.nex','A082113003.nex','A082113004.nex'};
%filelist = {'F092512001.nex','F092512002.nex','F092512003.nex','F092512004.nex'};
filelist = {'F090912001.nex','F090912002.nex','F090912003.nex','F090912004.nex','F090912005.nex','F090912006.nex'}% Fix move

offset = [.05 .05]; % in sec
granddata = [];
grandidx = [];
grandspikes = [];
for fileidx = 1:length(filelist)
    stro = nex2stro(findfile(char(filelist{fileidx})));
    if (stro.sum.paradigmID ~= 102)
        continue
    end
    spikeidxs = find(strncmp(stro.sum.rasterCells,'sig',3));
    hepidx = strcmp(stro.sum.rasterCells(1,:),'AD11');
    vepidx = strcmp(stro.sum.rasterCells(1,:),'AD12');
    starttimeidx = strcmp(stro.sum.rasterCells(1,:),'anlgStartTime');
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_on'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_off'));
    e1t = [stro.ras{:,starttimeidx}]';
    stimir = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_ir'));
    fpx = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fp_x'));
    fpy = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fp_y'));
    uniquefppos = unique([fpx fpy],'rows');  % These are guaranteed to be sorted
    bkgnd = unique(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'bkgnd')));
    samplerate = stro.sum.analog.storeRates{1};
    ntrials = size(stro.ras,1);
    % figure; axes; hold on;
    data = nan*ones(ntrials,ceil(max(stimoff_t+offset(2)-stimon_t-offset(1))*samplerate),2);  % Eye position snippet
    spikecounts = nan*ones(ntrials,length(spikeidxs));
    for i = 1:ntrials
        h = stro.ras{i,hepidx}.*4096/400;
        v = stro.ras{i,vepidx}.*4096/400;
        if (any(isnan(h)))
            keyboard
        end
        i
        % 4096 A/D levels = 10 V
        % 1 degree = 40 A/D levels (according to REX)
        % (4096 levels/10 V) * (1 degree/40 levels) = 4096/400.
        % The key is that both REX and PLEXON use 12 bit A/D boards configured
        % for +/- 5 V.
        
        neyesamp = size(h,1);
        x = linspace(e1t(i),neyesamp/samplerate+e1t(i),neyesamp);
        L = x >=stimon_t(i)+offset(1) & x <=stimoff_t(i)+offset(2);
        data(i,1:sum(L),1) = h(L);
        data(i,1:sum(L),2) = v(L);
        for j = 1:length(spikeidxs)
            spikecounts(i,j) = sum (stro.ras{i,spikeidxs(j)} > stimon_t(i)+offset(1) & stro.ras{i,spikeidxs(j)} < stimoff_t(i) + offset(2));
        end
    end
    grandidx = [grandidx; [fpx fpy (fpx == uniquefppos(1,1) & fpy == uniquefppos(1,2)) stro.trial(:,strcmp(stro.sum.trialFields(1,:),'bkgnd')) stimir repmat(fileidx,length(fpx),1)]];
    granddata = cat(1,granddata, data);
    grandspikes = cat(1,grandspikes, spikecounts);
end
% Columns of grandidx:
% 1) FPX   2) FPY   3) FPlocidx  4) bkgndidx   5) stimir   6) fileidx
% bkgndidx = 0 means corridor
%%
% Continued from above
% Just analysis here on down - no more disk acccess
if (size(grandspikes,2) > 1)
    whichspike = input(['Which spike? (n=',num2str(size(grandspikes,2)),')']);
else
    whichspike = 1;
end
% Now doing some analyses
Lfp = grandidx(:,3) == 1; % 1 = near
Lcorr = grandidx(:,4) == 0; % 0 means corridor background, 2 means gray.

% Generating a size tuning curve
figure; axes; hold on;
ringrad = unique(grandidx(:, 5));
for i = 0:1 % FP position
    for j = 0:1 % On the corridor?
        tmp = [];
        for k = 1:length(ringrad)
            Lbase = Lfp == i & Lcorr == j;
            L = grandidx(:,5) == ringrad(k);
            tmp(k,1) = mean(grandspikes(L&Lbase,whichspike));
            tmp(k,2) = std(grandspikes(L&Lbase,whichspike));
        end
       % plot([ringrad ringrad]'/10,[tmp(:,1)+tmp(:,2) tmp(:,1)-tmp(:,2)]','k-');
        h = plot(ringrad/10,tmp(:,1),'k-o','LineWidth',2);
        if (i == 1) % If in the near part of the corridor
           set(h,'Color','magenta') 
        end
        if (j == 0)  % if on gray bkgnd
            set(h,'LineStyle','--'); 
        end
    end
end
xlabel('ring radius (°)');
ylabel('spike count');
title('Dashed = gray background, pink = near');


% Let's just start with average fixation position in X and Y
% And pos vs time
avgpos = squeeze(nanmean(granddata,2));
x = linspace(offset(1),size(granddata,2)/samplerate,size(granddata,2));
figure; 
subplot(2,2,1); hold on;
plot(avgpos(Lfp&~Lcorr,1),avgpos(Lfp&~Lcorr,2),'o','Color',[.75 .75 .75]) % Gray means "on the corridor"
plot(avgpos(Lfp&Lcorr,1),avgpos(Lfp&Lcorr,2),'ko');
    plot([mean(avgpos(Lfp&~Lcorr,1)), mean(avgpos(Lfp&~Lcorr,1))+.01*stro.sum.exptParams.rf_x],...
         [mean(avgpos(Lfp&~Lcorr,2)), mean(avgpos(Lfp&~Lcorr,2))+.01*stro.sum.exptParams.rf_y],'b-','LineWidth',2);
plot([mean(avgpos(Lfp&~Lcorr,1)), mean(avgpos(Lfp&Lcorr,1))], [mean(avgpos(Lfp&~Lcorr,2)), mean(avgpos(Lfp&Lcorr,2))],'r-','LineWidth',2)
title('black = corridor, gray = gray bkgnd','Color','magenta');
axis square;
subplot(2,2,2); hold on;
plot(avgpos(~Lfp&~Lcorr,1),avgpos(~Lfp&~Lcorr,2),'o','Color',[.75 .75 .75])
plot(avgpos(~Lfp&Lcorr,1),avgpos(~Lfp&Lcorr,2),'ko')
plot([mean(avgpos(~Lfp&~Lcorr,1)), mean(avgpos(~Lfp&~Lcorr,1))+.01*stro.sum.exptParams.rf_x],...
    [mean(avgpos(~Lfp&~Lcorr,2)), mean(avgpos(~Lfp&~Lcorr,2))+.01*stro.sum.exptParams.rf_y],'b-','LineWidth',2);
plot([mean(avgpos(~Lfp&~Lcorr,1)), mean(avgpos(~Lfp&Lcorr,1))], [mean(avgpos(~Lfp&~Lcorr,2)), mean(avgpos(~Lfp&Lcorr,2))],'r-','LineWidth',2)
title('black = corridor, gray = gray bkgnd','Color','black');
axis square;

% Scaling the size of the mean fixation position to show the spike count

subplot(2,2,3); hold on;
for i = find(Lfp&Lcorr)'
    plot(avgpos(i,1),avgpos(i,2),'o','MarkerFaceColor','black','Color','black','MarkerSize',grandspikes(i,whichspike)*10./max(grandspikes(:,whichspike))+1);
end
for i = find(Lfp&~Lcorr)'
    plot(avgpos(i,1),avgpos(i,2),'o','MarkerFaceColor',[.65 .65 .65],'Color',[.65 .65 .65],'MarkerSize',grandspikes(i,whichspike)*10./max(grandspikes(:,whichspike))+1); % Gray means "on the gray bkgnd"
end
axis square;


subplot(2,2,4); hold on;
for i = find(~Lfp&Lcorr)'
    plot(avgpos(i,1),avgpos(i,2),'o','MarkerFaceColor','black','Color','black','MarkerSize',grandspikes(i,whichspike)*10./max(grandspikes(:,whichspike))+1);
end
for i = find(~Lfp&~Lcorr)'
    plot(avgpos(i,1),avgpos(i,2),'o','MarkerFaceColor',[.65 .65 .65],'Color',[.65 .65 .65],'MarkerSize',grandspikes(i,whichspike)*10./max(grandspikes(:,whichspike))+1); % Gray means "on the gray bkgnd"
end
axis square;

% This time course analysis didn't turn out to be very useful
% % commenting it out for now.
% subplot(2,2,3); hold on;
% plot(x,nanmean(granddata(Lfp&~Lcorr,:,1)),'r-');
% plot(x,nanmean(granddata(Lfp&Lcorr,:,1)),'m-');
% plot(x,nanmean(granddata(Lfp&~Lcorr,:,2)),'g-');
% plot(x,nanmean(granddata(Lfp&Lcorr,:,2)),'c-');
% subplot(2,2,4); hold on;
% plot(x,nanmean(granddata(~Lfp&~Lcorr,:,1)),'r-');
% plot(x,nanmean(granddata(~Lfp&Lcorr,:,1)),'m-');
% plot(x,nanmean(granddata(~Lfp&~Lcorr,:,2)),'g-');
% plot(x,nanmean(granddata(~Lfp&Lcorr,:,2)),'c-');

% Displaying some text
fprintf([filelist{1},' : ',filelist{end},'   Spike #',num2str(whichspike),'\n'])
disp([]);
% Permutation test
niter = 2000;
tmpdata = zeros(niter+1,2);
for fploc = 1:-1:0
    for i = 1:niter+1
        tmp = avgpos(Lfp == fploc,:);
        L = Lcorr(Lfp == fploc);
        if (i > 1)
            L = L(randperm(length(L)));
        end
        v = [mean(tmp(~L,1))-mean(tmp(L,1)), mean(tmp(~L,2))-mean(tmp(L,2))];
        tmpdata(i,:) = v;
    end
    dist2 = sum(tmpdata.^2,2);
    p = sum(dist2(2:end)>dist2(1))./niter;
    if (fploc == 1)
        fprintf(['Near: ']);
    else
        fprintf(['Far: ']);
    end
    disp(['amp: ',num2str(sqrt(dist2(1))),' angle: ',num2str(atan2(tmpdata(1,2),tmpdata(1,1))*180/pi), ' p: ',num2str(p)])
end

% Regression models
% Doing this on a fploc by fploc basis. Doesn't makes sense to model 
% both fixation point locations together
for fploc = 1:-1:0 % first is near, second is far.
    [p,t,stats] = anovan(grandspikes(Lfp == fploc,whichspike),[grandidx(Lfp == fploc,5) avgpos(Lfp == fploc,:)],'varnames',{'stimir','EPx','EPy'},'model','interaction','continuous',[2 3],'display','off');
    p = anovan(stats.resid, grandidx(Lfp == fploc,[4 5]),'varnames',{'bkgnidx','stimir'},'model','interaction','display','off');
    if (fploc == 1)
        disp(['Near: ',num2str(p')]);
    else
        disp(['Far: ',num2str(p')]);    
    end
    [p,t,stats]  = anovan(grandspikes(Lfp == fploc,whichspike),[grandidx(Lfp == fploc,[4 5]) avgpos(Lfp == fploc,:)],'varnames',{'bkgndidx','stimir','EPx','EPy'},'model','full','continuous',[3 4],'display','off');
    if (fploc == 1)
        disp('Near:');
        disp(t(:,[1 3 5 7]));
    else
        disp('Far:');
        disp(t(:,[1 3 5 7]))
    end
    disp(' ');
end

%p = anovan(grandspikes(:,whichspike),grandidx(:,[3 4 5]),'varnames',{'FPloc','bkgnidx','stimir'},'model','interaction');
% Spike count as a fucntion of Predictors: 3) FPlocidx  4) bkgndidx  5) stimir 


%%
% Permutation tests on every neighboring pair of fp locations
niter = 2000;
metadata =[];
for k = 0:1 % 1 = corridor background
    for i = 1:length(uniquefpx)
        Lx = grandidx(:,1) == uniquefpx(i);
        uniquefpy = unique(grandidx(Lx,2));
        for j = 1:length(uniquefpy)-1
            Ly = grandidx(:,2) == uniquefpy(j);
            sample1 = [zeros(sum(Lx&Ly&Lcorr==k),1) avgpos(Lx&Ly&Lcorr==k,1),avgpos(Lx&Ly&Lcorr==k,2)];
            Ly = grandidx(:,2) == uniquefpy(j+1);
            sample2 = [ones(sum(Lx&Ly&Lcorr==k),1) avgpos(Lx&Ly&Lcorr==k,1),avgpos(Lx&Ly&Lcorr==k,2)];
            sample = [sample1; sample2];
            data=nan*ones(niter,2);
            for iter = 1:niter
                if (iter == 1)
                    L = logical(sample(:,1));
                else
                    L = logical(sample(randperm(size(sample,1)),1));
                end
                data(iter,:) = [mean(sample(L,2)-sample(~L,2)) mean(sample(L,3)-sample(~L,3))];
            end
            stat = sqrt(sum(data.^2,2))
            p = sum(stat>=stat(1))./niter;
            metadata = [metadata; p data(1,:) (uniquefpy(j+1)-uniquefpy(j))/10]
        end
    end
end


%%
% Identical to the above cell but now working on a bigger filelist

% Note regarding eye position units
% 4096 A/D levels = 10 V
% 1 degree = 40 A/D levels (according to REX)
% (4096 levels/10 V) * (1 degree/40 levels) = 4096/400.
% The key is that both REX and PLEXON use 12 bit A/D boards configured
% for +/- 5 V.

[fnames, spikenums] = fnamesFromTxt2('/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/AmyFreya.txt');
daylist = [];
currentdate = [];
for i = 1:length(fnames)
    tmp = char(fnames{i});
    if (i == 1)
        daylist = 1;
        currentdate = str2num(tmp(2:end-7));
    else
        if (str2num(tmp(2:end-7)) == currentdate)
            daylist(i) = daylist(i-1);
        else
            daylist(i) = daylist(i-1)+1;
            currentdate = str2num(tmp(2:end-7));
        end
    end
end

offset = [.05 .05]; % in sec
metadata = nan*ones(daylist(end),7);
for daycounter = 1:daylist(end)
    granddata = [];
    grandidx = [];
    for fileidx = min(find(daylist==daycounter)):max(find(daylist==daycounter))
        stro = nex2stro(findfile(char(fnames{fileidx})));
        if (stro.sum.paradigmID ~= 102)
            continue
        end
        hepidx = strcmp(stro.sum.rasterCells(1,:),'AD11');
        vepidx = strcmp(stro.sum.rasterCells(1,:),'AD12');
        starttimeidx = strcmp(stro.sum.rasterCells(1,:),'anlgStartTime');
        stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_on'));
        stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_off'));
        stimx_near = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimx_near'));
        fpx = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fp_x'));
        fpy = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fp_y'));
        uniquefppos = unique([fpx fpy],'rows');  % These are guaranteed to be sorted
        bkgnd = unique(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'bkgnd')));
        samplerate = stro.sum.analog.storeRates{1};
        ntrials = size(stro.ras,1);
        % figure; axes; hold on;
        data = nan*ones(ntrials,ceil(max(stimoff_t+offset(2)-stimon_t-offset(1))*samplerate),2);  % Eye position snippet
        for i = 1:ntrials
            h = stro.ras{i,hepidx}.*4096/400;
            v = stro.ras{i,vepidx}.*4096/400; 
            e1t = stro.ras{i,starttimeidx};
            neyesamp = size(h,1);
            x = linspace(e1t,neyesamp/samplerate+e1t,neyesamp);
            L = x >=stimon_t(i)+offset(1) & x <=stimoff_t(i)+offset(2);
            data(i,1:sum(L),1) = h(L);
            data(i,1:sum(L),2) = v(L);
        end
        grandidx = [grandidx; [fpx fpy (fpx == uniquefppos(1,1) & fpy == uniquefppos(1,2)) stro.trial(:,strcmp(stro.sum.trialFields(1,:),'bkgnd')) repmat(fileidx,length(fpx),1)]];
        granddata = cat(1,granddata, data);
    end
    
    % Now doing some analyses
    Lfp = grandidx(:,3) == 1;  % 1 = lower position (near)
    Lcorr = grandidx(:,4) == 0;  % 0 means corridor. 2 means gray. I double checked.
    
    % Let's just start with average fixation position in X and Y
    % And pos vs time
    avgpos = squeeze(nanmean(granddata,2));
    x = linspace(offset(1),size(granddata,2)/samplerate,size(granddata,2));
    figure;
    subplot(2,2,1); hold on;
    plot(avgpos(Lfp&~Lcorr,1),avgpos(Lfp&~Lcorr,2),'ro');
    plot(avgpos(Lfp&Lcorr,1),avgpos(Lfp&Lcorr,2),'go'); % Green means "on the corridor"
    plot([mean(avgpos(Lfp&~Lcorr,1)), mean(avgpos(Lfp&~Lcorr,1))+.01*stro.sum.exptParams.rf_x],...
         [mean(avgpos(Lfp&~Lcorr,2)), mean(avgpos(Lfp&~Lcorr,2))+.01*stro.sum.exptParams.rf_y],'b-','LineWidth',2);
    plot([mean(avgpos(Lfp&~Lcorr,1)), mean(avgpos(Lfp&Lcorr,1))], [mean(avgpos(Lfp&~Lcorr,2)), mean(avgpos(Lfp&Lcorr,2))],'k-','LineWidth',2);
    axis square;
    
    subplot(2,2,2); hold on;
    plot(avgpos(~Lfp&~Lcorr,1),avgpos(~Lfp&~Lcorr,2),'ro')
    plot(avgpos(~Lfp&Lcorr,1),avgpos(~Lfp&Lcorr,2),'go')
    plot([mean(avgpos(~Lfp&~Lcorr,1)), mean(avgpos(~Lfp&~Lcorr,1))+.01*stro.sum.exptParams.rf_x],...
         [mean(avgpos(~Lfp&~Lcorr,2)), mean(avgpos(~Lfp&~Lcorr,2))+.01*stro.sum.exptParams.rf_y],'b-','LineWidth',2);
    plot([mean(avgpos(~Lfp&~Lcorr,1)), mean(avgpos(~Lfp&Lcorr,1))], [mean(avgpos(~Lfp&~Lcorr,2)), mean(avgpos(~Lfp&Lcorr,2))],'k-','LineWidth',2)
    axis square;
    
    subplot(2,2,3); hold on;
    tmpmn = nanmean(nanmean(granddata(Lfp,:,1)));
    plot(x,nanmean(granddata(Lfp&~Lcorr,:,1))-tmpmn,'r-');
    plot(x,nanmean(granddata(Lfp&Lcorr,:,1))-tmpmn,'m-');
   
    tmpmn = nanmean(nanmean(granddata(Lfp,:,2)));
    plot(x,nanmean(granddata(Lfp&~Lcorr,:,2))-tmpmn,'g-');
    plot(x,nanmean(granddata(Lfp&Lcorr,:,2))-tmpmn,'c-');
    
    subplot(2,2,4); hold on;
    tmpmn = nanmean(nanmean(granddata(~Lfp,:,1)));
    plot(x,nanmean(granddata(~Lfp&~Lcorr,:,1))-tmpmn,'r-');
    plot(x,nanmean(granddata(~Lfp&Lcorr,:,1))-tmpmn,'m-');
    tmpmn = nanmean(nanmean(granddata(~Lfp,:,2)));
    plot(x,nanmean(granddata(~Lfp&~Lcorr,:,2))-tmpmn,'g-');
    plot(x,nanmean(granddata(~Lfp&Lcorr,:,2))-tmpmn,'c-');
    set(gcf,'Name', char(fnames{fileidx}));

    % Permutation test
    niter = 2000;
    tmpdata = zeros(niter+1,2);
    for fploc = 1:-1:0
        for i = 1:niter+1
            tmp = avgpos(Lfp == fploc,:);
            L = Lcorr(Lfp == fploc);
            if (i > 1)
                L = L(randperm(length(L)));
            end
            v = [mean(tmp(L,1))-mean(tmp(~L,1)), mean(tmp(L,2))-mean(tmp(~L,2))];
            tmpdata(i,:) = v;
        end
        dist2 = sum(tmpdata.^2,2);
        p = sum(dist2(2:end)>dist2(1))./niter;
        disp(char(fnames{fileidx}));
        disp(['amp: ',num2str(sqrt(dist2(1))),' angle: ',num2str(atan2(tmpdata(1,2),tmpdata(1,1))*180/pi), ' p: ',num2str(p)])
        disp('------------');
        metadata(daycounter, (1:3)+3*fploc) = [sqrt(dist2(1)) atan2(tmpdata(1,2),tmpdata(1,1))*180/pi p];  % far 1:3, near 4:6
    end
    stro.sum.exptParams.eye_pos_control
    if (isnan(stro.sum.exptParams.eye_pos_control) | stro.sum.exptParams.eye_pos_control == 0)
        metadata(daycounter, 7) = 0;
    else
        metadata(daycounter, 7) = 1;
    end
end
% Angles are from average saccade endpoint on gray to average saccade
% endpoint on corridor. Angular histograms.  
% Angle is from gray EP centroid to corridor EP centroid
figure;
for control = 0:1
    L = metadata(:,7 ) == control;
    subplot(2,2,2+2*control);
    rose(metadata(L,2)*pi/180);
    subplot(2,2,1+2*control);
    rose(metadata(L,5)*pi/180);
end

% Compass plots
figure;
for control = 0:1
    L = metadata(:,7 ) == control;
    subplot(2,2,2+2*control); polar(0,.2,'w'); hold on; 
    x = metadata(L,1).*cos(metadata(L,2)*pi/180);
    y = metadata(L,1).*sin(metadata(L,2)*pi/180);
    compass(x,y);
    subplot(2,2,1+2*control); polar(0,.2,'w'); hold on; 
    x = metadata(L,4).*cos(metadata(L,5)*pi/180);
    y = metadata(L,4).*sin(metadata(L,5)*pi/180);
    compass(x,y);
end

%%
% Putting together some data for Scott. For each neuron pulling out the 
% mean eye positions (on the gray and corridor backgrounds) and the
% direction of the RF. Also grabbing the ring sizes so we can see how much
% the local curvature changes.

[fnames, spikenums] = fnamesFromTxt2('/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/AmyFreya.txt');
daylist = [];
currentdate = [];
for i = 1:length(fnames)
    tmp = char(fnames{i});
    if (i == 1)
        daylist = 1;
        currentdate = str2num(tmp(2:end-7));
    else
        if (str2num(tmp(2:end-7)) == currentdate)
            daylist(i) = daylist(i-1);
        else
            daylist(i) = daylist(i-1)+1;
            currentdate = str2num(tmp(2:end-7));
        end
    end
end

offset = [.05 .05]; % in sec
metadata = [];
for daycounter = 1:daylist(end)
    granddata = [];
    grandidx = [];
    for fileidx = min(find(daylist==daycounter)):max(find(daylist==daycounter))
        stro = nex2stro(findfile(char(fnames{fileidx})));
        if (stro.sum.paradigmID ~= 102)
            continue
        end
        hepidx = strcmp(stro.sum.rasterCells(1,:),'AD11');
        vepidx = strcmp(stro.sum.rasterCells(1,:),'AD12');
        starttimeidx = strcmp(stro.sum.rasterCells(1,:),'anlgStartTime');
        stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_on'));
        stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_off'));
       % stimx_near = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimx_near'));
        stimir = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_ir'));
        stimor = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_or'));

        fpx = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fp_x'));
        fpy = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fp_y'));
        uniquefppos = unique([fpx fpy],'rows');  % These are guaranteed to be sorted
        bkgnd = unique(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'bkgnd')));
        samplerate = stro.sum.analog.storeRates{1};
        ntrials = size(stro.ras,1);
        % figure; axes; hold on;
        data = nan*ones(ntrials,ceil(max(stimoff_t+offset(2)-stimon_t-offset(1))*samplerate),2);  % Eye position snippet
        for i = 1:ntrials
            h = stro.ras{i,hepidx}.*4096/400;
            v = stro.ras{i,vepidx}.*4096/400; 
            e1t = stro.ras{i,starttimeidx};
            neyesamp = size(h,1);
            x = linspace(e1t,neyesamp/samplerate+e1t,neyesamp);
            L = x >=stimon_t(i)+offset(1) & x <=stimoff_t(i)+offset(2);
            data(i,1:sum(L),1) = h(L);
            data(i,1:sum(L),2) = v(L);
        end
        rfx = stro.sum.exptParams.rf_x;
        rfy = stro.sum.exptParams.rf_y;
        Lfp = (fpx == uniquefppos(1,1) & fpy == uniquefppos(1,2));
        bkgnd = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'bkgnd'));
        grandidx = [grandidx; [fpx fpy Lfp bkgnd repmat([rfx rfy],length(Lfp),1) stimir stimor repmat(fileidx,length(fpx),1)]];
        granddata = cat(1,granddata, data);
    end
    
    % Now doing some analyses
    Lfp = grandidx(:,3) == 1;  % 1 = lower position (near)
    Lcorr = grandidx(:,4) == 0;  % 0 means corridor. 2 means gray. I double checked.
    
    % Let's just start with average fixation position in X and Y
    % And pos vs time
    avgpos = squeeze(nanmean(granddata,2));
    Lfp = grandidx(:,3) == 1;  % 1 = lower position (near)
    Lcorr = grandidx(:,4) == 0;  % 0 means corridor. 2 means gray. I double checked.
    fpnear = [unique(grandidx(Lfp,1)) unique(grandidx(Lfp,2))]/10;
    fpfar = [unique(grandidx(~Lfp,1)) unique(grandidx(~Lfp,2))]/10;
    rfpos = [unique(grandidx(:,5)) unique(grandidx(:,6))]/10;
    minringir = min(grandidx(:,7))/10;
    maxringir = min(grandidx(:,8))/10;
    if (isnan(stro.sum.exptParams.eye_pos_control) | stro.sum.exptParams.eye_pos_control == 0)
        whichgeometry = 0;
    else
        whichgeometry = 1;
    end
    
    metadata(daycounter,:) = [fpnear mean(avgpos(Lfp&~Lcorr,:)) mean(avgpos(Lfp&Lcorr,:)) fpfar mean(avgpos(~Lfp&~Lcorr,:)) mean(avgpos(~Lfp&Lcorr,:)) rfpos whichgeometry]
end
% columns
% 1) near FPx
% 2) near FPy
% 3) near EPx (gray)
% 4) near EPy (gray)
% 5) near EPx (corridor)
% 6) near EPy (corridor)
% 7) far FPx
% 8) far FPy
% 9) far EPx (gray)
% 10) far EPy (gray)
% 11) far EPx (corridor)
% 12) far EPy (corridor)
% 13) RFx
% 14) RFy
% 15) which geometry 0 = standard, 1 = eye movement control

nearEPdispX = metadata(:,5)-metadata(:,3);
nearEPdispY = metadata(:,6)-metadata(:,4);
farEPdispX = metadata(:,11)-metadata(:,9);
farEPdispY = metadata(:,12)-metadata(:,10);
nearEPdisplacement = sqrt(nearEPdispX.^2 + nearEPdispY.^2);
farEPdisplacement = sqrt(farEPdispX.^2 + farEPdispY.^2);
% computing a few statistics to quantify the magnitude of the eye position
% displacement vector
unitvectortowardsRF = mkbasis(metadata(:,[13 14])')';
% projection onto vector pointing toward RF
nearproj = sum([nearEPdispX nearEPdispY].*unitvectortowardsRF,2);
farproj = sum([farEPdispX farEPdispY].*unitvectortowardsRF,2);

% How much does the orientation of the tangent to the ring change
% due to changes in eye position?
% These are RF positions relative to the fixation point location
rfneargray = [metadata(:,13)+metadata(:,3)-metadata(:,1) metadata(:,14)+metadata(:,4)-metadata(:,2)]
rfnearcorridor = [metadata(:,13)+metadata(:,5)-metadata(:,1) metadata(:,14)+metadata(:,6)-metadata(:,2)]
thetaneargray = (pi-(pi/2+abs(atan(rfneargray(:,2)./rfneargray(:,1)))))*180/pi; % tangent to ring (deg)
thetanearcorridor = (pi-(pi/2+abs(atan(rfnearcorridor(:,2)./rfnearcorridor(:,1)))))*180/pi;
% Let see hw much the ring tangent changes with the mean shift in EP
thetaneargray-thetanearcorridor

% Again but now for the far fixation point
% These are RF positions relative to the fixation point location
rffargray = [metadata(:,13)+metadata(:,9)-metadata(:,7) metadata(:,14)+metadata(:,10)-metadata(:,8)]
rffarcorridor = [metadata(:,13)+metadata(:,11)-metadata(:,7) metadata(:,14)+metadata(:,12)-metadata(:,8)]
thetafargray = (pi-(pi/2+abs(atan(rffargray(:,2)./rffargray(:,1)))))*180/pi; % tangent to ring (deg)
thetafarcorridor = (pi-(pi/2+abs(atan(rffarcorridor(:,2)./rffarcorridor(:,1)))))*180/pi;
% Let see hw much the ring tangent changes with the mean shift in EP
thetafargray-thetafarcorridor

%EPstats = [nearEPdisplacement farEPdisplacement nearproj farproj abs(thetaneargray-thetanearcorridor) abs(thetafargray-thetafarcorridor)];

nearXgray = metadata(:,3);
nearYgray = metadata(:,4);
nearXcorridor = metadata(:,5);
nearYcorridor = metadata(:,6);
farXgray = metadata(:,9);
farYgray = metadata(:,10);
farXcorridor = metadata(:,11);
farYcorridor = metadata(:,12);

% Below, corridor relative to gray
v1 = [nearXcorridor-nearXgray nearYcorridor-nearYgray];
v2 = [farXcorridor-farXgray farYcorridor-farYgray];



EPstats = sum([v1-v2].*unitvectortowardsRF,2)
% Negative numbers mean that the monkey is looking down at the back and up at
% the front of the corridor. Positive numbers are consistent with the shift in the size
% tuning curves we've been seeing. Most of the numbers are negative because
% the monkey's behavior tends to work against the predicted shift.
% nearproj-farproj 
% This is the same as the thing above.

% Making a plot of eye position displacements between gray and corridor
% (separately for near and far and the two stimulus geometries - concentric
% and eccentric)
Lepcontrol = metadata(:,15) == 1;
figure;
AXLIMS = .2;
for i = 1:4
    subplot(2,2,i);
    if (i == 1)
        tmp = v1(~Lepcontrol,:);
    elseif (i == 2)
        tmp = v2(~Lepcontrol,:);
    elseif (i == 3)
        tmp = v1(Lepcontrol,:);
    elseif (i == 4)
        tmp = v2(Lepcontrol,:);
    end
    plot([zeros(size(tmp,1)) tmp(:,1)]',[zeros(size(tmp,1)), tmp(:,2)]','k-')
    set(gca,'Xlim',[-1 1]*AXLIMS,'Ylim',[-1 1]*AXLIMS); axis square;
    % A Rayleigh test
    unitvects = tmp./repmat(sqrt(sum(tmp.^2,2)),1,2)
    xy = mean(unitvects);
    r = norm(xy);
    n = size(unitvects,1)
    R = n*r;
    p = exp(sqrt(1+4*n+4*(n^2-R^2))-(1+2*n));
    title(num2str(p));
end

% A two-sample test to determine whether the angular distribution of mean eye
% position displacements depends on the stimulus geometry (concentric vs
% eccentric).

[th1,r1] = cart2pol(v1(~Lepcontrol,1), v1(~Lepcontrol,2));
[th2,r2] = cart2pol(v1(Lepcontrol,1), v1(Lepcontrol,2));
[~,p] = WatsonWheelerTest(th1,th2)  % near position comparison
[~,p] = ttest2(r1,r2)

[th1,r1] = cart2pol(v2(~Lepcontrol,1), v2(~Lepcontrol,2));
[th2,r2] = cart2pol(v2(Lepcontrol,1), v2(Lepcontrol,2));
[~,p] = WatsonWheelerTest(th1,th2) % far position comparison
[~,p] = ttest2(r1,r2)
[mean(r1) mean(r2)]  % For Freya eye displacements are little bigger in the
% standard version of the task than the eye position control version.


%save('ForScott','EPstats');
% % Projecting onto a few different unit vectors
% rs = []; ps = [];
% for i = 0:pi/10:2*pi
%     EPstats = [v1-v2]*[cos(i); sin(i)];
%     [r,p] = corrcoef([shift EPstats])
%     rs(end+1) = r(1,2);
%     ps(end+1) = p(1,2);
% end
% 

%%% Plot for Scott
Lepcontrol = metadata(:,15) == 1;
figure;
AXLIMS = .3;
for i = 1:4
    for j = 1:2
        if (j == 1)
            %  Below, all near (first gray then corridor)
            v1 = [metadata(:,3)-metadata(:,1) metadata(:,4)-metadata(:,2)];
            v2 = [metadata(:,5)-metadata(:,1) metadata(:,6)-metadata(:,2)];
        else
            % Below, all far (first gray then corridor) 
            v1 = [metadata(:,9)-metadata(:,7) metadata(:,10)-metadata(:,8)];
            v2 = [metadata(:,11)-metadata(:,7) metadata(:,12)-metadata(:,8)];
        end
        subplot(2,2,i); hold on;
        if (i == 1)
            tmp = v1(~Lepcontrol,:);
        elseif (i == 2)
            tmp = v2(~Lepcontrol,:);
        elseif (i == 3)
            tmp = v1(Lepcontrol,:);
        elseif (i == 4)
            tmp = v2(Lepcontrol,:);
        end
        h = plot(tmp(:,1),tmp(:,2),'ko');
        if (j == 1)
            set(h,'MarkerFaceColor','blue');
        else
            set(h,'MarkerFaceColor','red');    
        end
        set(gca,'Xlim',[-1 1]*AXLIMS,'Ylim',[-1 1]*AXLIMS); axis square;
        plot(0,0,'k+','LineWidth',2)
    end
end


% 

%% Looking at fix move data

filelist = {'F090912001.nex','F090912002.nex','F090912003.nex','F090912004.nex','F090912005.nex','F090912006.nex'}% Fix move

offset = [.05 .05]; % in sec
granddata = [];
grandidx = [];
grandspikes = [];
for fileidx = 1:length(filelist)
    stro = nex2stro(findfile(char(filelist{fileidx})));
    if (stro.sum.paradigmID ~= 102)
        continue
    end
    spikeidxs = find(strncmp(stro.sum.rasterCells,'sig',3));
    hepidx = strcmp(stro.sum.rasterCells(1,:),'AD11');
    vepidx = strcmp(stro.sum.rasterCells(1,:),'AD12');
    starttimeidx = strcmp(stro.sum.rasterCells(1,:),'anlgStartTime');
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_on'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_off'));
    e1t = [stro.ras{:,starttimeidx}]';
    stimir = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_ir'));
    fpx = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fp_x'));
    fpy = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fp_y'));
    uniquefppos = unique([fpx fpy],'rows');  % These are guaranteed to be sorted
    bkgnd = unique(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'bkgnd')));
    samplerate = stro.sum.analog.storeRates{1};
    ntrials = size(stro.ras,1);
    % figure; axes; hold on;
    data = nan*ones(ntrials,ceil(max(stimoff_t+offset(2)-stimon_t-offset(1))*samplerate),2);  % Eye position snippet
    spikecounts = nan*ones(ntrials,length(spikeidxs));
    for i = 1:ntrials
        h = stro.ras{i,hepidx}.*4096/400;
        v = stro.ras{i,vepidx}.*4096/400;
        if (any(isnan(h)))
            keyboard
        end
        i
        % 4096 A/D levels = 10 V
        % 1 degree = 40 A/D levels (according to REX)
        % (4096 levels/10 V) * (1 degree/40 levels) = 4096/400.
        % The key is that both REX and PLEXON use 12 bit A/D boards configured
        % for +/- 5 V.
        
        neyesamp = size(h,1);
        x = linspace(e1t(i),neyesamp/samplerate+e1t(i),neyesamp);
        L = x >=stimon_t(i)+offset(1) & x <=stimoff_t(i)+offset(2);
        data(i,1:sum(L),1) = h(L);
        data(i,1:sum(L),2) = v(L);
        for j = 1:length(spikeidxs)
            spikecounts(i,j) = sum (stro.ras{i,spikeidxs(j)} > stimon_t(i)+offset(1) & stro.ras{i,spikeidxs(j)} < stimoff_t(i) + offset(2));
        end
    end
    grandidx = [grandidx; [fpx fpy (fpx == uniquefppos(1,1) & fpy == uniquefppos(1,2)) stro.trial(:,strcmp(stro.sum.trialFields(1,:),'bkgnd')) stimir repmat(fileidx,length(fpx),1)]];
    granddata = cat(1,granddata, data);
    grandspikes = cat(1,grandspikes, spikecounts);
end

Lcorr = grandidx(:,4) == 0; % 0 means corridor background, 2 means gray.
avgpos = squeeze(nanmean(granddata,2));
fp_pos = unique(grandidx(:,[1 2]),'rows');
uniquefpx = unique(fp_pos(:,1));
cmap = colormap(jet(size(fp_pos(:,2),1)/2));
data = nan*ones(2,2,2,size(cmap,1));
figure;
for k = 0:1 % 1 = corridor background
    for i = 1:length(uniquefpx)
        subplot(2,2,i+k*2); hold on;
        Lx = grandidx(:,1) == uniquefpx(i);
        uniquefpy = unique(grandidx(Lx,2));        
        for j = 1:length(uniquefpy)
            Ly = grandidx(:,2) == uniquefpy(j);
            h = plot(avgpos(Lx&Ly&Lcorr==k,1),avgpos(Lx&Ly&Lcorr==k,2),'o');
            set(h,'Color',cmap(j,:));
            centroid = [mean(avgpos(Lx&Ly&Lcorr==k,1)), mean(avgpos(Lx&Ly&Lcorr==k,2))];
            h = plot(centroid(1),centroid(2),'p','MarkerSize',12); 
            set(h,'MarkerFaceColor',cmap(j,:),'MarkerEdgeColor','black');
            data(:,k+1,i,j) = centroid;
        end
        axis equal
        if (k)
            title('corridor');
        else
            title('gray');
        end
    end
end
%
% Y position as a function of y FP position
figure;
for k = 0:1 % 1 = corridor background
    for i = 1:length(uniquefpx)
        subplot(2,2,i+k*2); hold on;
        Lx = grandidx(:,1) == uniquefpx(i);
        uniquefpy = unique(grandidx(Lx,2));
        for j = 1:length(uniquefpy)
            Ly = grandidx(:,2) == uniquefpy(j);
            plot(uniquefpy(j)./10,avgpos(Lx&Ly&Lcorr==k,2),'ko');
        end
        plot([uniquefpy(1) uniquefpy(end)]/10,[uniquefpy(1) uniquefpy(end)]/10,'g-','LineWidth',2);
        axis square;
        set(gca,'Xlim',[uniquefpy(1)-5 uniquefpy(end)+5]/10);
        set(gca,'Ylim',[uniquefpy(1)-5 uniquefpy(end)+5]/10);
        if (k)
            title('corridor');
        else
            title('gray');
        end
    end
end

%%
% Population analysis of fix move data

[fnames, spikenums] = fnamesFromTxt2('/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/AmyFreyaFixMove.txt');
daylist = [];
currentdate = [];
for i = 1:length(fnames)
    tmp = char(fnames{i});
    if (i == 1)
        daylist = 1;
        currentdate = str2num(tmp(2:end-7));
    else
        if (str2num(tmp(2:end-7)) == currentdate)
            daylist(i) = daylist(i-1);
        else
            daylist(i) = daylist(i-1)+1;
            currentdate = str2num(tmp(2:end-7));
        end
    end
end

offset = [.05 .05]; % in sec
metadata = [];
for daycounter = 1:daylist(end)
    granddata = [];
    grandidx = [];
    for fileidx = min(find(daylist==daycounter)):max(find(daylist==daycounter))
        stro = nex2stro(findfile(char(fnames{fileidx})));
        if (stro.sum.paradigmID ~= 102)
            continue
        end
        spikeidxs = find(strncmp(stro.sum.rasterCells,'sig',3));
        hepidx = strcmp(stro.sum.rasterCells(1,:),'AD11');
        vepidx = strcmp(stro.sum.rasterCells(1,:),'AD12');
        starttimeidx = strcmp(stro.sum.rasterCells(1,:),'anlgStartTime');
        stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_on'));
        stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_off'));
        e1t = [stro.ras{:,starttimeidx}]';
        stimir = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_ir'));
        fpx = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fp_x'));
        fpy = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fp_y'));
        uniquefppos = unique([fpx fpy],'rows');  % These are guaranteed to be sorted
        bkgnd = unique(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'bkgnd')));
        samplerate = stro.sum.analog.storeRates{1};
        ntrials = size(stro.ras,1);
        % figure; axes; hold on;
        data = nan*ones(ntrials,ceil(max(stimoff_t+offset(2)-stimon_t-offset(1))*samplerate),2);  % Eye position snippet
        spikecounts = nan*ones(ntrials,length(spikeidxs));
        for i = 1:ntrials
            h = stro.ras{i,hepidx}.*4096/400;
            v = stro.ras{i,vepidx}.*4096/400;
            if (any(isnan(h)))
                keyboard
            end
            
            % 4096 A/D levels = 10 V
            % 1 degree = 40 A/D levels (according to REX)
            % (4096 levels/10 V) * (1 degree/40 levels) = 4096/400.
            % The key is that both REX and PLEXON use 12 bit A/D boards configured
            % for +/- 5 V.
            
            neyesamp = size(h,1);
            x = linspace(e1t(i),neyesamp/samplerate+e1t(i),neyesamp);
            L = x >=stimon_t(i)+offset(1) & x <=stimoff_t(i)+offset(2);
            data(i,1:sum(L),1) = h(L);
            data(i,1:sum(L),2) = v(L);
            for j = 1:length(spikeidxs)
                spikecounts(i,j) = sum (stro.ras{i,spikeidxs(j)} > stimon_t(i)+offset(1) & stro.ras{i,spikeidxs(j)} < stimoff_t(i) + offset(2));
            end
        end
        grandidx = [grandidx; [fpx fpy (fpx == uniquefppos(1,1) & fpy == uniquefppos(1,2)) stro.trial(:,strcmp(stro.sum.trialFields(1,:),'bkgnd')) stimir repmat(fileidx,length(fpx),1)]];
        if (~isempty(granddata) & size(granddata,2) > size(data,2))
            data = [data, nan*ones(size(data,1),size(granddata,2)-size(data,2),size(data,3))];
        end
        if (~isempty(granddata) & size(granddata,2) < size(data,2))
            granddata = [granddata, nan*ones(size(granddata,1),size(data,2)-size(granddata,2),size(granddata,3))];
        end
        granddata = cat(1,granddata, data);
    end

    Lcorr = grandidx(:,4) == 0; % 0 means corridor background, 2 means gray.
    avgpos = squeeze(nanmean(granddata,2));
    fp_pos = unique(grandidx(:,[1 2]),'rows');
    uniquefpx = unique(fp_pos(:,1));

    % Y position as a function of y FP position
    figure;
    for k = 0:1 % 1 = corridor background
        for i = 1:length(uniquefpx)
            subplot(2,2,i+k*2); hold on;
            Lx = grandidx(:,1) == uniquefpx(i);
            uniquefpy = unique(grandidx(Lx,2));
            for j = 1:length(uniquefpy)
                Ly = grandidx(:,2) == uniquefpy(j);
                plot(uniquefpy(j)./10,avgpos(Lx&Ly&Lcorr==k,2),'ko');
            end
            plot([uniquefpy(1) uniquefpy(end)]/10,[uniquefpy(1) uniquefpy(end)]/10,'g-','LineWidth',2);
            axis square;
            set(gca,'Xlim',[uniquefpy(1)-5 uniquefpy(end)+5]/10);
            set(gca,'Ylim',[uniquefpy(1)-5 uniquefpy(end)+5]/10);
            if (k)
                title('corridor');
            else
                title('gray');
            end
        end
    end
    
    % Permutation tests on every neighboring pair of fp locations
    niter = 2000;
    for k = 0:1 % 1 = corridor background
        for i = 1:length(uniquefpx)
            Lx = grandidx(:,1) == uniquefpx(i);
            uniquefpy = unique(grandidx(Lx,2));
            for j = 1:length(uniquefpy)-1
                Ly = grandidx(:,2) == uniquefpy(j);
                sample1 = [zeros(sum(Lx&Ly&Lcorr==k),1) avgpos(Lx&Ly&Lcorr==k,1),avgpos(Lx&Ly&Lcorr==k,2)];
                Ly = grandidx(:,2) == uniquefpy(j+1);
                sample2 = [ones(sum(Lx&Ly&Lcorr==k),1) avgpos(Lx&Ly&Lcorr==k,1),avgpos(Lx&Ly&Lcorr==k,2)];
                sample = [sample1; sample2];
                data=nan*ones(niter,2);
                for iter = 1:niter
                    if (iter == 1)
                        L = logical(sample(:,1));
                    else
                        L = logical(sample(randperm(size(sample,1)),1));
                    end
                    data(iter,:) = [mean(sample(L,2))-mean(sample(~L,2)) mean(sample(L,3))-mean(sample(~L,3))];
                end
                stat = sqrt(sum(data.^2,2))
                p = sum(stat>=stat(1))./niter;
                metadata = [metadata; p data(1,:) (uniquefpy(j+1)-uniquefpy(j))/10 i k daycounter]
            end
        end
    end
end

% How many fixation point locations used on each day in the "front" and the
% "back" positions. Seems like almost always 7 verically offset fixation 
% point locations (except one file). This gives 6 comparisons per file.
nfpcomparisons = [];
Lcor = logical(metadata(:,6));
for i = 1:metadata(end,end)
    Lday = metadata(:,end) == i;
    nfpcomparisons(i,1) = sum(Lday&~Lcor&metadata(:,5)==1); % near position, gray
    nfpcomparisons(i,2) = sum(Lday&~Lcor&metadata(:,5)==2); % far position, gray
    nfpcomparisons(i,3) = sum(Lday&Lcor&metadata(:,5)==1); % near position, corridor
    nfpcomparisons(i,4) = sum(Lday&Lcor&metadata(:,5)==2); % far position, corridor
end

% How many times do we reject the null hypothesis of same eye position distribution
% as a function of separation between fixation point y offsets?

uniquefpy = unique(metadata(:,4));
for i = 1:length(uniquefpy)
    L = metadata(:,4) == uniquefpy(i);
    uniquefpy(i)
    [sum(metadata(L&Lcor,1) < 0.05) sum(L&Lcor)]
    [sum(metadata(L&~Lcor,1) < 0.05) sum(L&~Lcor)]
end


figure; subplot(2,2,1);
hist(metadata(:,1),30);
subplot(2,2,2);
hist(metadata(:,3)-metadata(:,4),30);
subplot(2,2,3);
hist(metadata(:,3)./metadata(:,4),30);
geomean(metadata(:,3)./metadata(:,4))
geomean(metadata(Lcor,3)./metadata(Lcor,4))
geomean(metadata(~Lcor,3)./metadata(~Lcor,4))
[h,p] = ttest2(metadata(Lcor,3)./metadata(Lcor,4),metadata(~Lcor,3)./metadata(~Lcor,4))

% Making a plot where each mean eye position displacement is represented by
% a vector. Only looking at the 0.1 degree comparisons. Only looking at
% corridor data.
L = metadata(:,6) == 1; % Corridor background only
L = L & metadata(:,4) == .2; % Magnitude of fp displacements to consider
Lfar = metadata(:,5) == 2;
AXLIMS = .3;
figure;
subplot(2,2,1);
plot([zeros(sum(L&~Lfar),1) metadata(L&~Lfar,2)]',[zeros(sum(L&~Lfar),1)  metadata(L&~Lfar,3)]','k-');
set(gca,'Xlim',[-1 1]*AXLIMS,'Ylim',[-1 1]*AXLIMS); axis square;
subplot(2,2,2);
plot([zeros(sum(L&Lfar),1) metadata(L&Lfar,2)]',[zeros(sum(L&Lfar),1)  metadata(L&Lfar,3)]','k-');
set(gca,'Xlim',[-1 1]*AXLIMS,'Ylim',[-1 1]*AXLIMS); axis square;

%%
% Microsaccade analysis


[fnames, spikenums] = fnamesFromTxt2('/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/AmyApollo.txt');
daylist = [];
currentdate = [];
for i = 1:length(fnames)
    tmp = char(fnames{i});
    if (i == 1)
        daylist = 1;
        currentdate = str2num(tmp(2:end-7));
    else
        if (str2num(tmp(2:end-7)) == currentdate)
            daylist(i) = daylist(i-1);
        else
            daylist(i) = daylist(i-1)+1;
            currentdate = str2num(tmp(2:end-7));
        end
    end
end

offset = [.05 .05]; % in sec
granddata = [];
grandidx = [];
for daycounter = 1:daylist(end)
    for fileidx = min(find(daylist==daycounter)):max(find(daylist==daycounter))
        stro = nex2stro(findfile(char(fnames{fileidx})));
        if (stro.sum.paradigmID ~= 102)
            continue
        end
        sacstats = getSacData(stro); close;
        stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_on'));
        stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_off'));
        stimx_near = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimx_near'));
        fpx = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fp_x'));
        fpy = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fp_y'));
        uniquefppos = unique([fpx fpy],'rows');  % These are guaranteed to be sorted
        bkgnd = unique(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'bkgnd')));
        if (isnan(stro.sum.exptParams.eye_pos_control) | stro.sum.exptParams.eye_pos_control == 0)
            whichgeometry = 0;
        else
            whichgeometry = 1;
        end
        
        ntrials = size(stro.ras,1);
        for i = 1:ntrials
            [stimon_t(i)+offset(1) sacstats.starttimes{i}' stimoff_t(i)+offset(2)]
           
            L = sacstats.starttimes{i} > stimon_t(i)+offset(1) & sacstats.starttimes{i} < stimoff_t(i)+offset(2);
            if any(L)
                disp('found one!');
                for j = find(L)
                    j
                    grandidx = [grandidx; fpx(i) fpy(i) (fpx(i) == uniquefppos(1,1) & fpy(i) == uniquefppos(1,2)) stro.trial(i,strcmp(stro.sum.trialFields(1,:),'bkgnd')) whichgeometry daycounter];
                    granddata = [granddata; stimon_t(i) stimoff_t(i) sacstats.starttimes{i}(j) sacstats.directions{i}(j) sacstats.amplitudes{i}(j)];
                end
            end
        end
    end
end
% Indices into grandidx
% 1) fpx
% 2) fpy
% 3) near or far fixation point locations (1 = near)
% 4) What kind of background (0 = corridor, 2 = gray)
% 5) which geometry (0 = concentric, 1 = eccentric stimulus)
% 6) index for day on which files were recorded

% Indices into granddata
% 1) stim on time
% 2) stim off time
% 3) saccade start time
% 4) saccade direction
% 5) saccade amplitude

% First, amplitudes
figure;
subplotcounter = 1;
metadata = [];
EPCONTROL = 0;  % 0 = regular geometry, 1 = eye positon control.
for i = 1:-1:0 % 1 = near, 0 = far
    for j = [2 0] % 2 = gray, 0 = corridor
        L = grandidx(:,3) == i & grandidx(:,4) == j;
        L = L & grandidx(:,5) == EPCONTROL;  
        sum(L)
        subplot(2,2,subplotcounter);
        hist(granddata(L,5));
        if (subplotcounter == 1)
            title('near, gray');
        elseif (subplotcounter == 2)
            title('near, corridor');
        elseif (subplotcounter == 3)
            title('far, gray');
        elseif (subplotcounter == 4)
            title('far, corridor');
        end
        ylabel('count'); xlabel('µsac amp.');
        metadata = [metadata; granddata(L,5) repmat([i j],sum(L),1)];
        subplotcounter = subplotcounter+1;
    end
end
p =anovan(metadata(:,1),metadata(:,[2 3]),'model','interaction','varnames',{'near/far','bkgnd'});


% Second, directions
figure;
subplotcounter = 1;
metadata = [];
for i = 1:-1:0 % 1 = near, 0 = far
    for j = [2 0] % 2 = gray, 0 = corridor
        L = grandidx(:,3) == i & grandidx(:,4) == j;
        L = L & grandidx(:,5) == EPCONTROL;
        sum(L)
        subplot(2,2,subplotcounter);
        [t,r] = rose(granddata(L,4));
        h = polar(t,r,'k-');
        if (subplotcounter == 1)
            title('near, gray');
        elseif (subplotcounter == 2)
            title('near, corridor');
        elseif (subplotcounter == 3)
            title('far, gray');
        elseif (subplotcounter == 4)
            title('far, corridor');
        end
        metadata = [metadata; granddata(L,4) repmat([i j],sum(L),1)];
        subplotcounter = subplotcounter+1;
    end
end

% Third, directions and amplitudes jointly
figure;
subplotcounter = 1;
for i = 1:-1:0 % 1 = near, 0 = far
    for j = [2 0] % 2 = gray, 0 = corridor
        L = grandidx(:,3) == i & grandidx(:,4) == j;
        L = L & grandidx(:,5) == EPCONTROL;
        subplot(2,2,subplotcounter);
        [x,y] = pol2cart(granddata(L,4), granddata(L,5))
        compass(x,y);

        if (subplotcounter == 1)
            title('near, gray');
        elseif (subplotcounter == 2)
            title('near, corridor');
        elseif (subplotcounter == 3)
            title('far, gray');
        elseif (subplotcounter == 4)
            title('far, corridor');
        end
        subplotcounter = subplotcounter+1;
    end
end

% histogram of when microsaccades occur
subplotcounter = 1;
metadata =[];
for i = 1:-1:0 % 1 = near, 0 = far
    for j = [2 0] % 2 = gray, 0 = corridor
        L = grandidx(:,3) == i & grandidx(:,4) == j;
        L = L & grandidx(:,5) == EPCONTROL;
        subplot(2,2,subplotcounter);
        hist(granddata(L,3)-granddata(L,1),[0:.01:.3]);
        set(gca,'Xlim',[0 .3]);
        mean(granddata(L,3)-granddata(L,1))
        
        if (subplotcounter == 1)
            title('near, gray');
        elseif (subplotcounter == 2)
            title('near, corridor');
        elseif (subplotcounter == 3)
            title('far, gray');
        elseif (subplotcounter == 4)
            title('far, corridor');
        end
        xlabel('time from stimon (s)');
        ylabel('µ saccade count');
        
        subplotcounter = subplotcounter+1;
        metadata = [metadata; granddata(L,3)-granddata(L,1) repmat([i j],sum(L),1)];
    end
end
p =anovan(metadata(:,1),metadata(:,[2 3]),'model','interaction','varnames',{'near/far','bkgnd'});
mean(metadata(metadata(:,3) == 2),1) % gray
mean(metadata(metadata(:,3) == 0),1) % corridor