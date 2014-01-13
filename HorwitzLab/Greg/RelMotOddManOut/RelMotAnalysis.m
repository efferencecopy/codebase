stro = nex2stro;
targx = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_x'))/10;
targy = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_y'))/10;
correct = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
hepidx = strcmp(stro.sum.rasterCells(1,:),'AD11');
vepidx = strcmp(stro.sum.rasterCells(1,:),'AD12');
ntrials = size(correct,1);

figure; axes; hold on;
for i = 1:ntrials
    if (correct(i))
        plot(stro.ras{i,hepidx},stro.ras{i,vepidx},'g-');
    else
        plot(stro.ras{i,hepidx},stro.ras{i,vepidx},'r-');
    end
end

figure;
for i = 1:ntrials
    h = stro.ras{i,hepidx}*4096/400;
    v = stro.ras{i,vepidx}*4096/400;
    subplot(ceil(sqrt(ntrials)),ceil(sqrt(ntrials)),i);
    hold on;
    if (correct(i))
        plot(h,v,'g-');
        plot(targx(i),targy(i),'g*');
    else
        plot(h,v,'r-');
        plot(targx(i),targy(i),'r*');
    end
    set(gca,'Xlim',1.2*[min(targx) max(targx)])
    set(gca,'Ylim',1.2*[min(targy) max(targy)])
end