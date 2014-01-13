%% sort filenames based on their sf
filenamesbysf = cell(100,4);
counter = zeros(1,4);
DTNTsfs = [0.5008 0.9919 1.9839 3.9677];
for i = 1:length(fnamesthatpassed)
    stro = nex2stro(fnamesthatpassed{i});
    sf = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'sf'));
    if ~sf(1), continue; end
    err = abs(DTNTsfs-sf(1));
    thissf = err == min(err);
    counter(thissf) = counter(thissf) + 1;
    filenamesbysf{counter(thissf),thissf} = fnamesthatpassed{i};
end
% trim excess
filenamesbysf(max(counter)+1:end,:) = [];
% save('filenamesbysf0628.mat', 'filenamesbysf');

%% store planar fits
storeModelFits = true;
load('filenamesbysf0628'); load('T_cones_smj10');
if storeModelFits
    allStoredModelFits = cell(size(filenamesbysf,1),4);
end

for sfidx = 1:4 % sfidx of [1 2 3 4] -> [.5 1 2 4] cpd
    counter = 1;
    len = sum(~cellfun('isempty',filenamesbysf(:,sfidx)));
    for i = 1:len        
        NTfilename = filenamesbysf{i,sfidx};
        stro = nex2stro(findfile(NTfilename));

        lmsidxs = [find(strcmp(stro.sum.trialFields(1,:),'lcont')) ...
            find(strcmp(stro.sum.trialFields(1,:),'mcont')) ...
            find(strcmp(stro.sum.trialFields(1,:),'scont'))];

        fundamentals = stro.sum.exptParams.fundamentals;
        fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
        mon_spd = stro.sum.exptParams.mon_spd;
        mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
        mon_spd = SplineRaw((380:4:780)', mon_spd, (380:5:780)');
        M = fundamentals'*mon_spd;
        M10 = T_cones_smj10*mon_spd;

        stro.trial(:,lmsidxs) = stro.trial(:,lmsidxs)*100*(M10/M)';
        lms = stro.trial(:,lmsidxs);

        ppOut = NTpreprocess(stro,.4,1,0);
        Loog = logical(ppOut(:,end));
        NTscaled = ppOut(:,2:4) .* repmat(ppOut(:,5),1,3);

        spikeidx = strcmp(stro.sum.rasterCells(1,:),getSpikenum(stro, 'first'));
        levels = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'levelidx'));
        stepsize = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stepsize'));

        fpacq_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fp_acq'));
        stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_on'));
        stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_off'));

        Lbaseline = ~levels & all(lms == zeros(size(lms)),2);
        Lprobe = ~levels & ~all(lms == zeros(size(lms)),2);
        Lreplay = (levels == max(levels)) & ~stepsize;

        spikecounts = zeros(size(stro.trial,1),1);
        for j = 1:size(stro.trial,1)
            spiketimes = stro.ras{j,spikeidx};
            spikecounts(j) = sum(spiketimes > stimon_t(j)+stro.sum.exptParams.latency/1000 & spiketimes < stimoff_t(j));
        end

        [pparams,~,~,~,xformmat] = NTsurfacefit(NTscaled,Loog);
        coneweights = xformmat*pparams;

        Laccept = ~Lreplay & ~Lbaseline & ~Lprobe & stepsize;
        lms = lms(Laccept,:);
        spikecounts = spikecounts(Laccept);

        L = ~isnan(fpacq_t) & ~Lreplay;
        prestimsc = [];
        for j = find(L)'
            spiketimes = stro.ras{j,spikeidx};
            nspikes = sum(spiketimes > fpacq_t(j) & spiketimes < stimon_t(j));
            prestimsc = [prestimsc; nspikes./(stimon_t(j)-fpacq_t(j)) .* ...
                (stimoff_t(j)-(stimon_t(j)+stro.sum.exptParams.latency/1000))];
            % this makes the time interval the same as spikecounts
        end
        bline = mean(prestimsc);

        optsbase = optimset('MaxFunEvals',5e3,'MaxIter',5e3,'TolFun',1e-9,'TolX',1e-9,'Display','off');
        options = optimset(optsbase,'FinDiffType','central','Algorithm','sqp');

        cweights_0 = coneweights./norm(coneweights);
        b1_0 = norm(coneweights)*stro.sum.exptParams.threshold;
        b2_0 = 0;
        initguess = [cweights_0; b1_0; b2_0; bline];
        lb        = [cweights_0; zeros(2,1); bline];
        ub        = [cweights_0; inf(2,1)  ; bline];
        lmsfit = fmincon(@(x) modelfun(x,lms,spikecounts,[],[],'linplusquad','fish'), ...
            initguess,[],[],[],[],lb,ub,[],options);
        
        DTNTsfs = [0.5008 0.9919 1.9839 3.9677];
        params = [1.1817 2.8475 2.8556 3.0997;
                  0.1203 0.0147 -0.0035 0.0012;
                  -0.0052 -0.0033 0.0072 0.0126;
                  0.1445 0.2288 0.1570 0.0491;
                  -0.0367 -0.5428 0.2218 -0.2422;
                  -0.3620 -0.2225 0.4693 0.2491;
                  0.1134 0.0706 0.0446 0.0476;
                  0.7796 1.2509 -0.7186 -0.3305;
                  -0.6624 -1.2128 0.5492 -0.2676;
                  -0.0508 -0.0515 0.0766 -0.0105];
        params = params(:,sfidx);
        mechs = reshape(params(2:end),3,3);

        % Trying it in a few color directions now
        [xs ys zs] = sphere(200);
        orgsize = size(xs);
        [th,phi,~] = cart2sph(xs,ys,zs); % r is always 1
        wholesums = sum(abs([cos(phi(:)).*cos(th(:)) cos(phi(:)).*sin(th(:)) sin(phi(:))]*mechs).^params(1),2);
        newrs = (1./wholesums).^(1/params(1)); newrs = reshape(newrs,orgsize);
        [x y z] = sph2cart(th,phi,newrs);
        xyz = [x(:) y(:) z(:)];
        estresp = max(eps, lmsfit(4)*abs(xyz*lmsfit(1:3)) + lmsfit(5)*(xyz*lmsfit(1:3)).^2 + lmsfit(6));
        estresp = reshape(estresp,orgsize);
        [th,phi,distPlane] = cart2sph(coneweights(1),coneweights(2),coneweights(3));
        wholesum = sum(abs([cos(phi)*cos(th) cos(phi)*sin(th) sin(phi)]*mechs).^params(1),2);
        distDTSurf = (1/wholesum)^(1/params(1));

        normfact = distDTSurf./distPlane;

        if storeModelFits
            allStoredModelFits{counter,sfidx} = {NTfilename normfact lmsfit(:)};
            counter = counter + 1;
        else
            figure; axes; hold on;
            h = surf(x,y,z,estresp./max(estresp(:)));
            set(h,'EdgeColor','none','FaceAlpha',.2,'EdgeAlpha',.2);

            % Plotting the planes
            [rotmat,~] = qr(lmsfit(1:3));
            % rotmat = MakeOrtho([coneweights, normrnd(0,1,3,2)]);
            rotmat = rotmat(:,[2 3 1]);
            corners = rotmat*[-1 1 1 -1; 1 1 -1 -1; 1 1 1 1];
            color = [.5 .5 .5];
            h1 = patch(corners(1,:),corners(2,:),corners(3,:),color);
            set(h1,'CDataMapping','direct','FaceColor',color,'EdgeColor',color,'EdgeAlpha',.2,'FaceAlpha',.2);
            h2 = patch(-corners(1,:),-corners(2,:),-corners(3,:),color);
            set(h2,'CDataMapping','direct','FaceColor',color,'EdgeColor',color,'EdgeAlpha',.2,'FaceAlpha',.2);
            cb = colorbar('YTick',[0 .25 .5 .75 1],'YTickLabel',num2str(max(estresp(:))*[0 .25 .5 .75 1]',2));
            title([NTfilename ' - ' num2str(DTNTsfs(sfidx),2),' cpd, nf = ' num2str(normfact)]);

            ispts = intersectPlaneSurf(corners(:,1),lmsfit(1:3),x,y,z);
            hold on, fill3(-ispts(:,1),-ispts(:,2),-ispts(:,3),[.5 .5 .5]);
            fill3(ispts(:,1),ispts(:,2),ispts(:,3),[.5 .5 .5]);
        end
    end
end
if storeModelFits
    save('allStoredModelFits.mat','allStoredModelFits');
end

%% now that we have stored planar fits, get and sum all estimated responses
load('allStoredModelFits');
sfidx = 1;
len = sum(~cellfun('isempty',allStoredModelFits(:,sfidx)));

params = [1.1817   2.8475   2.8556   3.0997;
    0.1203   0.0147  -0.0035   0.0012;
    -0.0052  -0.0033   0.0072   0.0126;
    0.1445   0.2288   0.1570   0.0491;
    -0.0367  -0.5428   0.2218  -0.2422;
    -0.3620  -0.2225   0.4693   0.2491;
    0.1134   0.0706   0.0446   0.0476;
    0.7796   1.2509  -0.7186  -0.3305;
    -0.6624  -1.2128   0.5492  -0.2676;
    -0.0508  -0.0515   0.0766  -0.0105];
params = params(:,sfidx);
mechs = reshape(params(2:end),3,3);

[xs ys zs] = sphere(200);
orgsize = size(xs);
[th,phi,~] = cart2sph(xs,ys,zs); % r is always 1
wholesums = sum(abs([cos(phi(:)).*cos(th(:)) cos(phi(:)).*sin(th(:)) sin(phi(:))]*mechs).^params(1),2);
newrs = (1./wholesums).^(1/params(1));
newrs = reshape(newrs,orgsize);
[x y z] = sph2cart(th,phi,newrs);
xyz = [x(:) y(:) z(:)];

currStoredFits = allStoredModelFits(:,sfidx);
allEstResps = zeros([orgsize len]);
for i = 1:len
    currFit = currStoredFits{i}{3};
    if ~isempty(currFit)
        estresps = max(eps, currFit(4)*abs(xyz*currFit(1:3)) + currFit(5)*(xyz*currFit(1:3)).^2 + currFit(6));
        %         estresps = estresps/max(estresps);
        allEstResps(:,:,i) = reshape(estresps,orgsize);
    end
end

summedResps = sum(allEstResps,3);
figure; axes; hold on;
h = surf(x,y,z,summedResps);
set(h,'EdgeColor','none');
colorbar; caxis([0 max(summedResps(:))]);
title(sprintf('Summed responses %0.2f cpd', DTNTsfs(sfidx)));

%% add code to plot individual responses (like I do for quads below)?


%% get quad fits
storeModelFitsQuad = true;
load('filenamesbysf0628');
if storeModelFitsQuad
    allStoredModelFitsQuad = cell(size(filenamesbysf,1),4);
end
load('T_cones_smj10');
% sfidx = 1; % sfidx of [1 2 3 4] -> [.5 1 2 4] cpd
for sfidx = 1:4
    counter = 1;
    len = sum(~cellfun('isempty', filenamesbysf(:,sfidx)));
    for i = 1:len
        NTfilename = filenamesbysf{i,sfidx};
        stro = nex2stro(findfile(NTfilename));
        
        if ~isempty(strfind(NTfilename,'K070309004'))
            sf = .5;
        else
            sf = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'sf'));
            if ~sf(1)
                disp('No sf');
                continue;
            end
        end
        
        lmsidxs = [find(strcmp(stro.sum.trialFields(1,:),'lcont')) ...
            find(strcmp(stro.sum.trialFields(1,:),'mcont')) ...
            find(strcmp(stro.sum.trialFields(1,:),'scont'))];
        
        fundamentals = stro.sum.exptParams.fundamentals;
        fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
        mon_spd = stro.sum.exptParams.mon_spd;
        mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
        mon_spd = SplineRaw((380:4:780)', mon_spd, (380:5:780)');
        M = fundamentals'*mon_spd;
        M10 = T_cones_smj10*mon_spd;
        
        stro.trial(:,lmsidxs) = stro.trial(:,lmsidxs)*100*(M10/M)';
        lms = stro.trial(:,lmsidxs);
        
        ppOut = NTpreprocess(stro,.4,1,0);
        Loog = logical(ppOut(:,end));
        NTscaled = ppOut(:,2:4) .* repmat(ppOut(:,5),1,3);
        
        if sum(~Loog) < 10
            error('not enough in-gamut termination points');
        end
        
        spikeidx = strcmp(stro.sum.rasterCells(1,:),getSpikenum(stro, 'first'));
        levels = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'levelidx'));
        stepsize = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stepsize'));
        
        fpacq_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fp_acq'));
        stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_on'));
        stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_off'));
        
        Lbaseline = ~levels & all(lms == zeros(size(lms)),2);
        Lprobe = ~levels & ~all(lms == zeros(size(lms)),2);
        Lreplay = (levels == max(levels)) & ~stepsize;
        
        spikecounts = zeros(size(stro.trial,1),1);
        for j = 1:size(stro.trial,1)
            spiketimes = stro.ras{j,spikeidx};
            spikecounts(j) = sum(spiketimes > stimon_t(j)+stro.sum.exptParams.latency/1000 & spiketimes < stimoff_t(j));
        end
        
        L = ~isnan(fpacq_t) & ~Lreplay;
        prestimsc = [];
        for j = find(L)'
            spiketimes = stro.ras{j,spikeidx};
            nspikes = sum(spiketimes > fpacq_t(j) & spiketimes < stimon_t(j));
            prestimsc = [prestimsc; nspikes./(stimon_t(j)-fpacq_t(j)) .* ...
                (stimoff_t(j)-(stimon_t(j)+stro.sum.exptParams.latency/1000))];
            % this makes the time interval the same as spikecounts
        end
        bline = mean(prestimsc);
        
        [~,~,quadparams,~,xformmat] = NTsurfacefit(NTscaled,Loog);
        
        Laccept = ~Lreplay & ~Lbaseline & ~Lprobe & stepsize;
        spikecounts = spikecounts(Laccept);
        lms = lms(Laccept,:);
        
        A = [quadparams(1) quadparams(4) quadparams(5); ...
            quadparams(4) quadparams(2) quadparams(6); ...
            quadparams(5) quadparams(6) quadparams(3)];
        quadparams = xformmat*A*xformmat';
        
        optsbase = optimset('MaxFunEvals',5e3,'MaxIter',5e3,'TolFun',1e-9,'TolX',1e-9,'Display','off');
        options = optimset(optsbase,'FinDiffType','central','Algorithm','sqp');
        
        [th,ph,r] = cart2sph(lms(:,1),lms(:,2),lms(:,3));
        
        predr2 = 1./(quadparams(1).*(cos(ph(:)).*cos(th(:))).^2 +...
            quadparams(2).*(cos(ph(:)).*sin(th(:))).^2 +...
            quadparams(3).*sin(ph(:)).^2 + ...
            2*quadparams(4).*cos(ph(:)).*cos(th(:)).*cos(ph(:)).*sin(th(:)) +...
            2*quadparams(5).*cos(ph(:)).*cos(th(:)).*sin(ph(:)) +...
            2*quadparams(6).*cos(ph(:)).*sin(th(:)).*sin(ph(:)));
        predr2(predr2<0) = inf;
        
        b1_0 = stro.sum.exptParams.threshold - bline;
        b2_0 = 0;
        initguess = [b1_0; b2_0; bline];
        lb        = [zeros(2,1); bline];
        ub        = [inf(2,1)  ; bline];
        newfit = fmincon(@(x) modelfun(x,lms,spikecounts,r,predr2,'newguy','fish'), ...
            initguess,[],[],[],[],lb,ub,[],options);
        
        DTNTsfs = [0.5008 0.9919 1.9839 3.9677];
        err = abs(DTNTsfs-sf(1));
        whichsf = err == min(err);
        params = [1.1817   2.8475   2.8556   3.0997;
            0.1203   0.0147  -0.0035   0.0012;
            -0.0052  -0.0033   0.0072   0.0126;
            0.1445   0.2288   0.1570   0.0491;
            -0.0367  -0.5428   0.2218  -0.2422;
            -0.3620  -0.2225   0.4693   0.2491;
            0.1134   0.0706   0.0446   0.0476;
            0.7796   1.2509  -0.7186  -0.3305;
            -0.6624  -1.2128   0.5492  -0.2676;
            -0.0508  -0.0515   0.0766  -0.0105];
        params = params(:,whichsf);
        mechs = reshape(params(2:end),3,3);
        
        [xs ys zs] = sphere(200);
        orgsize = size(xs);
        [th,ph,~] = cart2sph(xs,ys,zs); % r is always 1
        
        predr2dt = 1./(quadparams(1).*(cos(ph(:)).*cos(th(:))).^2 +...
            quadparams(2).*(cos(ph(:)).*sin(th(:))).^2 +...
            quadparams(3).*sin(ph(:)).^2 + ...
            2*quadparams(4).*cos(ph(:)).*cos(th(:)).*cos(ph(:)).*sin(th(:)) +...
            2*quadparams(5).*cos(ph(:)).*cos(th(:)).*sin(ph(:)) +...
            2*quadparams(6).*cos(ph(:)).*sin(th(:)).*sin(ph(:)));
        predr2dt(predr2dt<0) = inf;
        
        wholesums = sum(abs([cos(ph(:)).*cos(th(:)) cos(ph(:)).*sin(th(:)) sin(ph(:))]*mechs).^params(1),2);
        newrs = (1./wholesums).^(1/params(1)); newrs = reshape(newrs,orgsize);
        
        estresp = max(eps, newfit(1)*(newrs(:)./sqrt(predr2dt)) + newfit(2)*(newrs(:)./sqrt(predr2dt)) + newfit(3));
        if storeModelFitsQuad
            allStoredModelFitsQuad{counter,sfidx} = {NTfilename; ...
                                               reshape(estresp,orgsize); ...
                                               newfit(:);
                                               [spikecounts r./sqrt(predr2)]};
            counter = counter+1;
        end
    end
end
if storeModelFitsQuad
    save('allStoredModelFitsQuad.mat','allStoredModelFitsQuad');
end

%% use the stored quad fits for summed responses
load('allStoredModelFitsQuad.mat');
DTNTsfs = [0.5008 0.9919 1.9839 3.9677];
allparams = [1.1817   2.8475   2.8556   3.0997;
    0.1203   0.0147  -0.0035   0.0012;
    -0.0052  -0.0033   0.0072   0.0126;
    0.1445   0.2288   0.1570   0.0491;
    -0.0367  -0.5428   0.2218  -0.2422;
    -0.3620  -0.2225   0.4693   0.2491;
    0.1134   0.0706   0.0446   0.0476;
    0.7796   1.2509  -0.7186  -0.3305;
    -0.6624  -1.2128   0.5492  -0.2676;
    -0.0508  -0.0515   0.0766  -0.0105];

for sfidx = 1:4
    len = sum(~cellfun('isempty', allStoredModelFitsQuad(:,sfidx)));
    summedresp = zeros(201,201);
    for cellnum = 1:len
        currcell = allStoredModelFitsQuad{cellnum,sfidx};
        summedresp = summedresp + currcell{2}; % idx 2 has the responses
    end
    params = allparams(:,sfidx);
    mechs = reshape(params(2:end),3,3);
    
    [xs,ys,zs] = sphere(200);
    orgsize = size(xs);
    [th,phi,~] = cart2sph(xs,ys,zs);
    wholesums = sum(abs([cos(phi(:)).*cos(th(:)) cos(phi(:)).*sin(th(:)) sin(phi(:))]*mechs).^params(1),2);
    newrs = (1./wholesums).^(1/params(1));
    newrs = reshape(newrs,orgsize);
    [x y z] = sph2cart(th,phi,newrs);
    
    figure; axes; hold on;
    h = surf(x,y,z,summedresp);
    set(h,'EdgeColor','none');
    colorbar; caxis([0 max(summedresp(:))]);
    title(sprintf('Summed responses %0.2f cpd', DTNTsfs(sfidx)));
end

%% use the stored quad fits for individual responses
SINGLECOLORBAR = false; % have one colorbar for all figures?
load('allStoredModelFitsQuad.mat');
DTNTsfs = [0.5008 0.9919 1.9839 3.9677];
allparams = [1.1817   2.8475   2.8556   3.0997;
    0.1203   0.0147  -0.0035   0.0012;
    -0.0052  -0.0033   0.0072   0.0126;
    0.1445   0.2288   0.1570   0.0491;
    -0.0367  -0.5428   0.2218  -0.2422;
    -0.3620  -0.2225   0.4693   0.2491;
    0.1134   0.0706   0.0446   0.0476;
    0.7796   1.2509  -0.7186  -0.3305;
    -0.6624  -1.2128   0.5492  -0.2676;
    -0.0508  -0.0515   0.0766  -0.0105];

for sfidx = 1 % change this as desired--note: will run out of memory if multiple sfidxs are displayed
    params = allparams(:,sfidx);
    mechs = reshape(params(2:end),3,3);
    [x,y,z] = sphere(200);
    orgsize = size(x);
    [th,phi,~] = cart2sph(x,y,z);
    wholesums = sum(abs([cos(phi(:)).*cos(th(:)) cos(phi(:)).*sin(th(:)) sin(phi(:))]*mechs).^params(1),2);
    newrs = (1./wholesums).^(1/params(1));
    newrs = reshape(newrs,orgsize);
    [x y z] = sph2cart(th,phi,newrs);
    
    len = sum(~cellfun('isempty', allStoredModelFitsQuad(:,sfidx)));
    maxresp = max(cellfun(@(x) max(max(x{2})), allStoredModelFitsQuad(1:len,sfidx)));
    for cellnum = 1:len
        currcell = allStoredModelFitsQuad{cellnum,sfidx};
        estresp = currcell{2}; % idx 2 has the responses
        
        figure; subplot(211); hold on;
        h = surf(x,y,z,estresp);
        set(h,'EdgeColor','none');
        if SINGLECOLORBAR
            colorbar; caxis([0 maxresp]);
        else
            colorbar; caxis([0 max(estresp(:))]);
        end
        title(sprintf('%s %0.2g cpd\n(b1, b2, base) = (%0.2g, %0.2g, %0.2g)', currcell{1}, DTNTsfs(sfidx), currcell{3}));
        subplot(212); hold on;
        data = currcell{4}; currfit = currcell{3};
        Xls = linspace(min(data(:,2)),max(data(:,2)),100);
        plot(data(:,2), data(:,1), 'k.', Xls, currfit(1)*Xls + currfit(2)*Xls + currfit(3), 'b');
    end
end