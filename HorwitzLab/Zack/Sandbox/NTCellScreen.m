%% Put the calculations into a structure
datapath = 'N:\NexFiles';
dirstruct = subdir([datapath filesep '*.nex']);

if exist('okay103files', 'var')
    fidxs = ismember({dirstruct.name}, okay103files)';
    dirstruct = dirstruct(fidxs);
end

len = length(dirstruct);
screen103S(1,len) = struct('one', [], 'two', [], 'three', [], 'four', [], ...
    'five', [], 'six', [], 'seven', [], 'filename', []);

for i = 1:len
    stro = nex2stro(dirstruct(i).name);
    
    out1 = NTpreprocess(stro,0,inf,inf); %prelim filter
    if ~isempty(out1)
        nprelim = length(out1(:,end));
        nprelimoog = sum(out1(:,end) == 1);
        screen103S(i).one = nprelim-nprelimoog;
        screen103S(i).two = nprelimoog;
    else
        screen103S(i).one = NaN;
        screen103S(i).two = NaN;
    end
    
    out2 = NTpreprocess(stro,0.4,1,0); %stringent filter
    if ~isempty(out2)
        nstringe = length(out2(:,end));
        nstringeoog = sum(out2(:,end) == 1);
        screen103S(i).three = nstringe-nstringeoog;
        screen103S(i).four = nstringeoog;
    else
        screen103S(i).three = NaN;
        screen103S(i).four = NaN;
    end
    
    spikeidx = strcmp(stro.sum.rasterCells(1,:), getSpikenum(stro, 'first'));
    levels = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'levelidx'));
    stepsize = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stepsize'));
    
    lmsidxs = [find(strcmp(stro.sum.trialFields(1,:),'lcont')) ...
        find(strcmp(stro.sum.trialFields(1,:),'mcont')) ...
        find(strcmp(stro.sum.trialFields(1,:),'scont'))];
    
    fpacq_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fp_acq'));
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_on'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_off'));
    
    if any(isempty([levels stepsize fpacq_t stimon_t stimoff_t]))
        continue
    end
    
    lms = stro.trial(:,lmsidxs);
    
    try
        Lreplay = (levels == max(levels)) & ~stepsize;
        Lstation = levels > 0 & stepsize == 0 & ~Lreplay;
    catch exception
        fprintf('%s There must have been NaN''s in ''levels''\n', dirstruct(i).name);
        continue;
    end
    
    % In what proportion of trials does the baseline activity exceed threshold?
    L = ~isnan(fpacq_t) & ~Lreplay;
    if ~any(L), continue; end
    
    prestimfr = [];
    for j = find(L)'
        spiketimes = stro.ras{j,spikeidx};
        nspikes = sum(spiketimes > fpacq_t(j) & spiketimes < stimon_t(j));
        prestimfr = [prestimfr; nspikes./(stimon_t(j)-fpacq_t(j))];
    end
    abovethresh = prestimfr > stro.sum.exptParams.threshold;
    screen103S(i).five = sum(abovethresh)/numel(prestimfr);
    
    % regression predictions of the first and last prestimulus and probe
    % trial responses
    spikerates = zeros(size(stro.trial,1),1);
    for j = 1:size(stro.trial,1)
        spiketimes = stro.ras{j,spikeidx};
        nspikes = sum(spiketimes > stimon_t(j)+stro.sum.exptParams.latency/1000 & spiketimes < stimoff_t(j));
        spikerates(j) = nspikes./(stimoff_t(j)-stimon_t(j)-stro.sum.exptParams.latency/1000);
    end
    
    Xps = [(1:sum(L))' ones(sum(L),1)];
    bhatps = Xps\prestimfr;
    yhatps = Xps*bhatps;
    
    Xpt = [(1:sum(Lstation))' ones(sum(Lstation),1)];
    bhatpt = Xpt\spikerates(Lstation);
    yhatpt = Xpt*bhatpt;
    
    if ~numel(yhatps) || ~numel(yhatpt)
        continue
    end
    
    if ~sum(spikerates(Lstation))
        fprintf('cell %d - There weren''t any spikes during the probe trials\n', i);
    end
    
    screen103S(i).six = [yhatps(1) yhatps(end)];
    screen103S(i).seven = [yhatpt(1) yhatpt(end)];
    screen103S(i).filename = dirstruct(i).name;
end

%% Graph the calculations from the structure // outdated!
load('screen103S.mat')

len = length(screen103S);
vals = [([screen103S.one] + [screen103S.two])' ([screen103S.four] ./ [screen103S.two])' ...
    (([screen103S.three] + [screen103S.four]) ./ ([screen103S.one] + [screen103S.two]))' ...
    [screen103S.five]' [screen103S.six]' [screen103S.seven]'];

fnames = {screen103S.filename}';

minvals = min(vals); maxvals = max(vals);

figure;
subplot(2,3,1); hold on;
for i = 1:len
    if isnan(vals(i,1))
        h = plot(i, 0, 'ko', 'MarkerSize', 2, 'MarkerFaceColor', 'k');
    else
        h = plot(i, vals(i,1), 'ko', 'MarkerSize', 2, 'MarkerFaceColor', 'k');
    end
    set(h, 'ButtonDownFcn', ['mybuttondown(fnames{' num2str(i) '}, vals(' num2str(i) ',:))']);
end
title('# points after prelim');
tick = get(gca, 'ytick'); tick = (tick(2)-tick(1))*0.1;
set(gca, 'ylim', [0-tick maxvals(1)+tick]);

subplot(2,3,2); hold on;
for i = 1:len
    if isnan(vals(i,2))
        h = plot(i, 1.03, 'ro', 'MarkerSize', 2, 'MarkerFaceColor', 'r');
    else
        h = plot(i, vals(i,2), 'ko', 'MarkerSize', 2, 'MarkerFaceColor', 'k');
    end
    set(h, 'ButtonDownFcn', ['mybuttondown(fnames{' num2str(i) '}, vals(' num2str(i) ',:))']);
end
title('OOG stringent / OOG prelim (red NaN)');
tick = get(gca, 'ytick'); tick = (tick(2)-tick(1))*0.1;
set(gca, 'ylim', [minvals(2)-tick 1+max(.1,tick)]);

subplot(2,3,3); hold on;
for i = 1:len
    if isnan(vals(i,3))
        h = plot(i, 1.03, 'ro', 'MarkerSize', 2, 'MarkerFaceColor', 'r');
    else
        h = plot(i, vals(i,3), 'ko', 'MarkerSize', 2, 'MarkerFaceColor', 'k');
    end
    set(h, 'ButtonDownFcn', ['mybuttondown(fnames{' num2str(i) '}, vals(' num2str(i) ',:))']);
end
title('# stringent points / # prelim points');
tick = get(gca, 'ytick'); tick = (tick(2)-tick(1))*0.1;
set(gca, 'ylim', [minvals(3)-tick 1+max(.1,tick)]);

subplot(2,3,4); hold on;
for i = 1:len
    h = plot(i, vals(i,4), 'ko', 'MarkerSize', 2, 'MarkerFaceColor', 'k');
    set(h, 'ButtonDownFcn', ['mybuttondown(fnames{' num2str(i) '}, vals(' num2str(i) ',:))']);
end
title('Prop. where background > threshold');
tick = get(gca, 'ytick'); tick = (tick(2)-tick(1))*0.1;
set(gca, 'ylim', [minvals(4)-tick maxvals(4)+tick]);

subplot(2,3,5); hold on;
for i = 1:len
    if isnan(vals(i,5))
        h = plot(i, 1, 'ro', 'MarkerSize', 2, 'MarkerFaceColor', 'r');
    elseif vals(i,5) == -1 %1 trial
        h = plot(i, 1, 'go', 'MarkerSize', 2, 'MarkerFaceColor', 'g');
    else
        h = plot(i, vals(i,5), 'ko', 'MarkerSize', 2, 'MarkerFaceColor', 'k');
    end
    set(h, 'ButtonDownFcn', ['mybuttondown(fnames{' num2str(i) '}, vals(' num2str(i) ',:))']);
end
title('p-val probe FRs (green <2 trials)');
tick = get(gca, 'ytick'); tick = (tick(2)-tick(1))*0.1;
tmpcol = vals(:,5);
set(gca, 'ylim', [min(tmpcol(tmpcol ~= -1))-tick max(tmpcol(tmpcol ~= -1))+tick]);

subplot(2,3,6); hold on;
for i = 1:len
    if isnan(vals(i,6))
        h = plot(i, 1, 'ro', 'MarkerSize', 2, 'MarkerFaceColor', 'r');
    elseif vals(i,6) == -1 %1 trial
        h = plot(i, 1, 'go', 'MarkerSize', 2, 'MarkerFaceColor', 'g');
    else
        h = plot(i, vals(i,6), 'ko', 'MarkerSize', 2, 'MarkerFaceColor', 'k');
    end
    set(h, 'ButtonDownFcn', ['mybuttondown(fnames{' num2str(i) '}, vals(' num2str(i) ',:))']);
end
title('p-val baseline FRs');
tick = get(gca, 'ytick'); tick = (tick(2)-tick(1))*0.1;
tmpcol = vals(:,6);
set(gca, 'ylim', [min(tmpcol(tmpcol ~= -1))-tick max(tmpcol(tmpcol ~= -1))+tick]);

%% Final culling
datapath = 'N:\NexFiles';
dirstruct = subdir([datapath filesep '*.nex']);

if exist('okay103files', 'var')
    fidxs = ismember({dirstruct.name}, okay103files)';
    dirstruct = dirstruct(fidxs);
end

fnamesthatpassed = {}; pvals = [];
for i = 1:length(dirstruct)
    try
        stro = nex2stro(dirstruct(i).name);
    catch exception
        wh = findall(0, 'tag', 'TMWWaitbar');
        if ~isempty(wh), delete(wh); end
        continue
    end
    
    out1 = NTpreprocess(stro,0,inf,inf); %prelim filter
    if ~isempty(out1)
        out1(out1(:,end) == 1,:) = [];
        if length(unique(out1(:,1))) < 20 % at least 20 color dirs
            continue
        end
    end
    
    out2 = NTpreprocess(stro,0.4,1,0); %stringent filter
    if ~isempty(out2)
        if length(out2(:,end)) - sum(out2(:,end) == 1) < 10 % 10 nOOG points
            continue
        end
    else
        continue
    end
    
    %     spikeidx = strcmp(stro.sum.rasterCells(1,:), getSpikenum(stro, 'first'));
    %     levels = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'levelidx'));
    %     stepsize = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stepsize'));
    %
    %     lmsidxs = [find(strcmp(stro.sum.trialFields(1,:),'lcont')) ...
    %         find(strcmp(stro.sum.trialFields(1,:),'mcont')) ...
    %         find(strcmp(stro.sum.trialFields(1,:),'scont'))];
    %
    %     fpacq_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fp_acq'));
    %     stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_on'));
    %     stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_off'));
    %
    %     Lreplay = (levels == max(levels)) & ~stepsize;
    %     Lstation = levels > 0 & stepsize == 0 & ~Lreplay;
    
    %     % Relaxing this restrictions for now
    %     % In what proportion of trials does the baseline activity exceed threshold?
    %     L = ~isnan(fpacq_t) & ~Lreplay;
    %     prestimfr = [];
    %     for j = find(L)'
    %         spiketimes = stro.ras{j,spikeidx};
    %         nspikes = sum(spiketimes > fpacq_t(j) & spiketimes < stimon_t(j));
    %         prestimfr = [prestimfr; nspikes./(stimon_t(j)-fpacq_t(j))];
    %     end
    %
    %     abovethresh = prestimfr > stro.sum.exptParams.threshold;
    %     if sum(abovethresh)/numel(prestimfr) > 0.05
    %         continue
    %     end
    
    %     spikerates = zeros(size(stro.trial,1),1);
    %     for j = 1:size(stro.trial,1)
    %         spiketimes = stro.ras{j,spikeidx};
    %         nspikes = sum(spiketimes > stimon_t(j)+stro.sum.exptParams.latency/1000 & spiketimes < stimoff_t(j));
    %         spikerates(j) = nspikes./(stimoff_t(j)-stimon_t(j)-stro.sum.exptParams.latency/1000);
    %     end
    
    % tests here on the regression predictions, using Lstation
    
    Loog = logical(out2(:,end));
    NTscaled = out2(:,2:4) .* repmat(out2(:,5),1,3);
    [~,~,quadparams,~,~] = NTsurfacefit(NTscaled,Loog);
    %     p = plnquadbootstrap(NTscaled*xformmat,Loog,pparams,psse,qsse,2000);
    
    %     if p < 0.01, continue; end
    
    A = [quadparams(1) quadparams(4) quadparams(5);...
        quadparams(4) quadparams(2) quadparams(6);...
        quadparams(5) quadparams(6) quadparams(3)];
    [~,evals] = eig(A);
    
    if all(diag(evals)) && (sum(diag(evals)<0) == 1 || sum(diag(evals)<0) == 2)
        fnamesthatpassed = [fnamesthatpassed; {dirstruct(i).name}];
    end
end

%%
% here we're going to try the comparing curvatures idea to differentiate
% between planar or quadric surfaces
load('okay103files.mat'); pleaseout = {};
for i = 1:length(okay103files)
    fname = okay103files{i}(end-13:end-4);
    try
        stro = nex2stro(okay103files{i});
    catch exception
        wh = findall(0, 'tag', 'TMWWaitbar');
        if ~isempty(wh), delete(wh); end
        continue
    end
    
    out1 = NTpreprocess(stro,0,inf,inf); %prelim filter
    if isempty(out1), continue; end
    out1(out1(:,end) == 1,:) = []; % cull OOGs
    if length(unique(out1(:,1))) < 10 % at least 10 color dirs
        continue
    end
    
    out2 = NTpreprocess(stro,0.4,1,0); %stringent filter
    if isempty(out2) || length(out2(:,end)) - sum(out2(:,end) == 1) < 10 % at least 10 nOOGs
        continue
    end
    
    Loog = logical(out2(:,end));
    NTscaled = out2(:,2:4) .* repmat(out2(:,5),1,3);
    [~,~,quadparams,~,xformmat] = NTsurfacefit(NTscaled,Loog);
    xyz = NTscaled*xformmat;
    
    A = [quadparams(1) quadparams(4) quadparams(5);...
        quadparams(4) quadparams(2) quadparams(6);...
        quadparams(5) quadparams(6) quadparams(3)];
    [evecs,evals] = eig(A);
    
    if ~softEq(det(evecs),1,5) % det(evecs) = 1 for a rotation matrix
        evecs = evecs(:,[2 1 3]);
        evals = evecs'*A*evecs;
    end
    
    evals = diag(evals);
    a2 = 1/evals(1); b2 = 1/evals(2); c2 = 1/evals(3);
    absevals = sort(abs(evals));
    types = {'hyperboloid 2 sheets' 'hyperboloid 1 sheet' 'ellipsoid'};
    
    typenum = sum(evals>0);
    switch typenum
        case 1
            typestr = types{1};
            lim = max(abs(xyz(:)))/3;
        case 2
            typestr = types{2};
            lim = max(abs(xyz(:)))/1.1;
        case 3
            typestr = types{3};
            lim = max(abs(xyz(:)))*1.1;
    end
    
    figure; axes; hold on;
    plot3(xyz(~Loog,1),xyz(~Loog,2),xyz(~Loog,3),'k.');
    plot3(-xyz(~Loog,1),-xyz(~Loog,2),-xyz(~Loog,3),'k.');
    plot3(xyz(Loog,1),xyz(Loog,2),xyz(Loog,3),'y.');
    plot3(-xyz(Loog,1),-xyz(Loog,2),-xyz(Loog,3),'y.');
    
    tmplattice = linspace(-lim,lim,25);
    [xx yy zz] = meshgrid(tmplattice);
    variables = [xx(:).^2 yy(:).^2 zz(:).^2 2*xx(:).*yy(:) 2*xx(:).*zz(:) 2*yy(:).*zz(:)];
    fr = reshape(variables*quadparams,size(xx));
    
    [faces,vertices] = isosurface(xx,yy,zz,fr,1);
    transvars = evecs'*vertices';
    
    S = sign(evals);
    K = prod(S) ./ (a2*b2*c2*(S(1)^2*transvars(1,:).^2/a2^2 + S(2)^2*transvars(2,:).^2/b2^2 + S(3)^2*transvars(3,:).^2/c2^2)).^2;
    H = (-S(1)^2*b2*c2*(S(3)*b2 + S(2)*c2)*transvars(1,:).^2 + S(2)^2*a2*c2*(S(3)*a2 + S(1)*c2)*transvars(2,:).^2 + S(3)^2*a2*b2*(S(2)*a2 + S(1)*b2)*transvars(3,:).^2) ...
        ./ (2*a2^2*b2^2*c2^2*(S(1)^2*transvars(1,:).^2/a2^2 + S(2)^2*transvars(2,:).^2/b2^2 + S(3)^2*transvars(3,:).^2/c2^2).^(3/2));
    chi = [H + sqrt(H.^2 - K); H - sqrt(H.^2 - K)]; % damn you imaginary numbers
%     chi(imag(chi)~=0) = nan;
%     devplane = sum(chi.^2,1);
    
    p = patch('faces',faces,'vertices',vertices,'facevertexcdata',K','facecolor','interp','facealpha',1,'edgealpha',0);
    isonormals(xx,yy,zz,fr,p);
    view(3); axis equal;
    title(sprintf('%s\n%s\n%0.3g %0.3g', fname, typestr, min(K), max(K)));
    
%     close(gcf);
%     pleaseout{end+1} = {fname;types{typenum};{tmplattice;fr;H';K';chi'}};
end

%     switch typenum
%         case 1 % 2 sheets
%             K = c2^3 ./ (c2^2 - (a2 + c2)*transvars(3,:).^2).^2;
% %             K = a2^3*b2^3*c2 ./ (a2^2*(b2 + c2)*transvars(2,:).^2 + b2^4*(a2^2 + (a2 + c2)*transvars(1,:).^2)).^2;
%         case 2 % 1 sheet
%             K = -c2^3 ./ (c2^2 + (a2 + c2)*transvars(3,:).^2).^2;
% %             K = -a2^3*b2^3*c2 ./ (b2^2*(a2^2 - (a2 + c2)*transvars(1,:).^2) -a2^2*(b2+c2)*transvars(2,:).^2).^2;
%         case 3 % ellipsoid
%             K = a2*b2^3*c2^3 ./ (c2^2*b2^2 + c2^2*(a2 - b2)*transvars(2,:).^2 + b2^2*(a2-c2)*transvars(3,:).^2).^2;
%     end

%     [K,~,~] = zcurvature(fr-1);
%     fvc = isosurface(xx,yy,zz,fr,1,K);
%     p = patch(fvc,'facecolor','interp','facealpha',1,'edgealpha',0);