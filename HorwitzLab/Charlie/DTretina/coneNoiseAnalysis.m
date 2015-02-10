%% IMPORT THE NECESSARY DATA FILE(S)

%clear out the workspace
fin



presdir = pwd;
switch whoami
    case 'hass_mbp' % charlie's laptop
        cd '/Users/charliehass/LabStuff/Huskies/DTcones/Data/';
    case 'nuke'
        cd 'C:\Users\charlie\Desktop\Local Data\DTcones\Data';
end
[fname, fpath] = uigetfile('Pick a text file');
load([fpath, fname]);
cd(presdir);

% print any notes that are attached to the file
if ~isempty(params.notes)
    fprintf('************** NOTES FOR FILE ********************\n\n')
    fprintf(params.notes)
    fprintf('\n\n')
end


% if the simulated expt was based off a V1 DT stro file, load that too...
if strcmpi(params.runType, 'DTV1')
    DT = dtobj(params.DTV1_fname);
    [monk, cell, expt] = DTunpack(DT, 1);
end

% identify the presence/absence of monte-carlo type trials
MONTECARLO = size(idlob.resp,3) > 0

% identify the presence/absence of V1 data
V1DATA = exist('monk', 'var')

% Calculate the area under the ROC curve, and the response as a function of
% contrast (if need be, e.g., obsMethods that require LDA), and fit the
% neurometric function
[idlob, cones] = coneNoiseROC(params, idlob, cones, gab);

%% CONE MOSAIC CONTRAST RESPONSE FUNCTIONS

%
% FIGURE 1
%
% CRFs plotted separately for each color condition. The monte-carlo and
% analytic estimates are overlayed (when possible)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nRows = 4;
nColors = size(idlob.resp,1);
nFigs = ceil(nColors ./ nRows^2);
contrastToPlot = vertcat(gab.contrasts{:}); % also used by some of the other figures below
contrastToPlot(:,1) = contrastToPlot(:,2) .* .20;
idx = 1;
for f = 1:nFigs
    figure;
    set(gcf, 'position', [274     5   983   824])
    
    if nColors-idx < 16;
        nPlots = (nColors-idx)+1;
    else
        nPlots = nRows^2;
    end
    
    for a = 1:nPlots
        subplot(nRows,nRows, a)
        hold on,
        
        % unpack the monte carlo trials
        if MONTECARLO
            switch params.obsMethod
                case 'obsMethod_all'
                    tmpResp = permute(idlob.resp(idx,:,:), [3,2,1]);
                case 'obsMethod_noClrEqSpace'
                    % Do LDA on the 3D clouds of data points
                    tmpResp = nan(size(idlob.resp,3),nContrasts);
                    for cnt = 1:nContrasts
                        noise = cat(2, idlob.resp{idx,1,:})'; % <nTrials x 3>
                        mu_noise = mean(noise,1);
                        cov_noise = cov(noise);
                        
                        sig = cat(2, idlob.resp{idx,cnt,:})';
                        mu_sig = mean(sig,1);
                        cov_sig = cov(sig);
                        
                        % the "within groups" and "between groups" scatter
                        S_within = (cov_noise) + (cov_sig);
                        S_between = (mu_noise - mu_sig)' * (mu_noise - mu_sig); % the outer product
                        
                        % calculate the discriminant vector
                        [vec, val] = eig(S_within \ S_between); % inv(S_within) * S_between
                        [~, maxEigVal] = max(diag(val));
                        W_mc = vec(:,maxEigVal);
                        if W_mc(1)< 0;
                            W_mc = -W_mc; % standardize the representation of eigenvectors
                        end
                        tmpResp(:,cnt) = sig * W_mc;
                    end
            end
            
            % plot the monte carlo data
            plot(contrastToPlot(idx,:), tmpResp, 'bo', 'markerfacecolor', 'b')
            plot(contrastToPlot(idx,:), mean(tmpResp), 'b')
        end
        
        % now deal with the analytic solution
        errorbar(contrastToPlot(idx,:), idlob.analyticMean(idx,:), sqrt(idlob.analyticVar(idx,:)), '-m.', 'markersize', 10, 'linewidth', 2)
        
        % standardize the axes
        set(gca, 'xscale', 'log')
        axis tight
        text_y = max(get(gca, 'ylim')) .* .9;
        text_x = min(get(gca, 'xlim')) .* 1.15;
        text(text_x, text_y, sprintf('[%s]', num2str(gab.colorDirs(idx,:),3)))
        hold off
        
        % update the counter
        idx = idx + 1;
    end
    
    %add labels
    subplot(nRows,nRows, 1)
    ylabel('Ideal Observer Response')
    subplot(nRows, nRows, a)
    xlabel('Contrast')
end

%
% FIGURE 2
%
% CRFs plotted on top of each other. CRFs should be exponential (on a
% log plot), and should be identical when (1) the bkgnd is equated and (2)
% the mosaics are equated
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure, hold on, % CRFs plotted together for all color dirs
if MONTECARLO
    plot(contrastToPlot', mean(idlob.resp, 3)', '-', 'linewidth', 2)
else
    plot(contrastToPlot', idlob.analyticMean', '-', 'linewidth', 2)
end
set(gca, 'xscale', 'log')
xlabel('Contrast')
ylabel('Mean response')
axis tight
hold off


%
% FIGURE 3
%
% Variance vs. contrast. The model's assumption is that the variance of the
% cones is independent of contrast. This assumption should extend to the
% ideal observer responses. Also, the analytic solution to the ideal
% observer variance should match that of the monte-carlo approach.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if MONTECARLO
    figure, hold on % Var vs. Contrast for all the color dirs.
    SEV = var(idlob.resp,[],3) .* (sqrt(2./(size(idlob.resp,3)-1)));
    errorbar(contrastToPlot', var(idlob.resp, [], 3)', SEV', 'o-', 'linewidth', 2)
    plot(contrastToPlot', idlob.analyticVar', '--', 'linewidth', 2)
    set(gca, 'xscale', 'log')
    xlabel('Contrast')
    ylabel('Sample Variance')
    axis tight
    hold off
end



%% NERUOMETRIC FUNCTIONS FOR THE CONE MOSAIC
close all


%
% FIGURE 1
%
% Neurometric (cones and V1) and psychometric function plotted separately
% for each color condition.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nRows = 4;
nColors = size(idlob.resp,1);
nFigs = ceil(nColors ./ nRows^2);
idx = 1;
for f = 1:nFigs
    figure;
    set(gcf, 'position', [274     5   983   824])
    
    if nColors == 1;
        nPlots = 1;
    elseif nColors-idx < 16;
        nPlots = (nColors-idx)+1;
    else
        nPlots = nRows^2;
    end
    
    for a = 1:nPlots
        subplot(nRows, nRows, a), hold on,
        tmp_norms = gab.contrasts{idx};
        xx = logspace(log10(tmp_norms(2).*.3), log10(tmp_norms(end).*1.06), 150);
        
        % plot the model data for the cone mosaic monte-carlo trials
        if MONTECARLO
            mod_cones = 1-0.5.*exp(-(xx./cones.alpha(idx)).^cones.beta(idx));
            plot(xx, mod_cones, 'r')
        end
        
        % plot the model data for the analytic solution
        mod_cones_analytic = 1-0.5.*exp(-(xx./cones.alpha_analytic(idx)).^cones.beta_analytic(idx));
        plot(xx, mod_cones_analytic, 'm')
        plot(cones.alpha_analytic(idx), 0.816, 'ms', 'markerfacecolor', 'c')
        
        % plot the model data for V1
        if V1DATA
            mod_monk = 1-0.5.*exp(-(xx./monk.alpha(idx)).^monk.beta(idx));
            mod_v1 = 1-0.5.*exp(-(xx./cell.alpha(idx)).^cell.beta(idx));
            plot(xx, mod_monk, 'k')
            plot(xx, mod_v1, 'b')
        end
        
        % plot the raw data
        if V1DATA
            plot(tmp_norms, idlob.roc{idx}, 'r.');
            plot(expt.norms{idx}, monk.performance{idx}, 'k.')
            plot(expt.norms{idx}, cell.roc{idx}, 'b.')
            plot(tmp_norms, idlob.roc_analytic{idx}, 'm.')
        else
            plot(tmp_norms, idlob.roc{idx}, 'r.');
            plot(tmp_norms, idlob.roc_analytic{idx}, 'm.')
        end
        
        %adjust the axes
        set(gca, 'xscale', 'log')
        xlim([tmp_norms(2).*.96, tmp_norms(end).*1.05])
        ylim([0.4, 1.02])
        title(sprintf('Color: [%s]', num2str(gab.colorDirs(idx,:),2)))
        
        % increment the counter
        idx = idx+1;
    end
    
    % add labels
    subplot(nRows,nRows, 1)
    ylabel('Ideal Observer Response')
    subplot(nRows, nRows, a)
    xlabel('Contrast')
end



%% ISO-DETECTION SURFACES FOR MONKEY AND RETINA
clc; close all;

% Set the parameters of the analysis
observer = 'sedna';         % Kali_DTNT_0713.mat or Sedna_DTNT_0713.mat
sfIdx = logical([0;1;0;0])  % which SF should be analyzed? [0.5 1 2 4] for kali, [0.5, 1, 2, 4] for sedna
viewSetting = 'isolum';   % could also be 'lm plane', 's vs l+m', 'isolum'
plottype = '3D';

% load the data
clear fpar_monkey
switch observer
    case 'kali'
        load Kali_DTNT_072113.mat % monkey data
        load Kali_fpar_072113.mat % fits to monkey surfaces
    case 'sedna'
        load Sedna_DTNT_072113.mat % monkey data
        load Sedna_fpar_072113.mat % fits to monkey surfaces
end


% divide the behavioral data by 100 so that % is b/w 0&1.
monkeyAlphas = cellfun(@(x,y) x./y, dtnt.alpha',  repmat({100}, numel(dtnt.alpha),1), 'uniformoutput', 0);

% which SF condition should be analyzed? The behvioral data has many SFs,
% but the model data is only at one SF.
sfs = unique(cat(1,dtnt.sfs{:}))
sfs(sfIdx)
warning('hard coding the SF')

% make a logical vector to weed out conditions (if need be)
l_questBeta = dtnt.l_beta > 0;
l_bkgndrgb = dtnt.l_bkgndrgb > 0;

% pull out the spatial frequncy condition
l_sfs = cat(1,dtnt.sfs{:}) == sfs(sfIdx);

% grab the appropriate files
l_valid = [l_sfs & l_questBeta & l_bkgndrgb];
colors = cat(1, dtnt.colorDirs{l_valid});
alphas = cat(1, monkeyAlphas{l_valid});

% l_goodQuestThresh = cat(1, dtnt.validThresh{l_valid});
% colors = colors(l_goodQuestThresh,:);
% alphas = alphas(l_goodQuestThresh);
size(colors)

% find the surface fits or the monkey raw data. Only compute them denovo if
% they don't already exist
if exist('fpar_monkey', 'var')
    fpar_monkey = fpar_monkey{sfIdx}(:);
else
    
    % default init params on the monkey
    [fpar_monkey_default, fval_monkey_default] = fitDetectionSurface(colors, alphas, 'ellipsoid');
    
    % init params from kali at 0.5 cpd (on the monkey)
    initparams_kali = [1.15; -15.82; 38.29; -14.36; -17.29; 6.62; -8.74; -81.82; 68.40; 9.36];  % optional input argument for threshSurfPlot for 0.5cpd
    [fpar_monkey_kali, fval_monkey_kali] = fitDetectionSurface(colors, alphas, initparams_kali);
    
    % particle swarm on the monkey
    npsoiters = 2;
    [fpar_monkey_pso, fval_monkey_pso] = fitDetectionSurface(colors, alphas, 'pso', npsoiters);
    
    
    % find the best method
    method = {'default'; 'kali'; 'pso'};
    fvals = [fval_monkey_default; fval_monkey_kali; fval_monkey_pso]
    fpars = [fpar_monkey_default'; fpar_monkey_kali'; fpar_monkey_pso];
    [~, idx] = min(fvals);
    bestMethod = method{idx}
    fpar_monkey = fpars(idx,:);
end

% fit the cone noise model's data
fpar_cones = fitDetectionSurface(gab.colorDirs, cones.alpha_analytic, 'ellipsoid');


% plot the raw data
figure;
set(gcf, 'position', [56 5 1377 801])
subplot(1,2,1); % behavioral data
title(sprintf('Spatial Frequency = %.2f cpd', sfs(sfIdx)))
threshSurfPlot(colors, alphas, viewSetting, plottype, fpar_monkey)

subplot(1,2,2)% cone noise data
title(sprintf('Spatial Frequency = %.2f cpd', gab.sf))
threshSurfPlot(gab.colorDirs, cones.alpha_analytic, viewSetting, plottype, fpar_cones)




%
% determine the aspect ratio of each ellipsoid
%
%%%%%%%%%%%
colors = [1./sqrt(2), -1./sqrt(2), 0;...
          1./sqrt(2), 1./sqrt(2), 0;...
          0, 0, 1];
cardAlpha_monkey = coleThresh(reshape(fpar_monkey(2:end),3,3), fpar_monkey(1), colors);
aspect_monkey = cardAlpha_monkey(3) ./ cardAlpha_monkey(1)
cardAlpha_cones = coleThresh(reshape(fpar_cones(2:end),3,3), fpar_cones(1), colors);
aspect_cones = cardAlpha_cones(3) ./ cardAlpha_cones(1)

%
% make a map of the ratio b/w Monkey and Retina. Where the ratio is large,
% the monkey does worse than the retina....
%
%%%%%%

[X,Y,Z] = sphere(400);
colors = [X(:), Y(:), Z(:)];
colors = bsxfun(@rdivide, colors, sqrt(sum(colors.^2, 2)));
monkeyThresh = coleThresh(reshape(fpar_monkey(2:end),3,3), fpar_monkey(1), colors);
conesThresh = coleThresh(reshape(fpar_cones(2:end),3,3), fpar_cones(1), colors);
C = monkeyThresh ./ conesThresh;
C = log10(C);
C = reshape(C, size(X,1), size(X,2)); % for plotting
minVal = min(C(:));
minVec = [X(C(:) == minVal), Y(C(:) == minVal), Z(C(:) == minVal)] .* 1.5;
if size(minVec,1) ==1; minVec = [minVec; -minVec]; end
maxVal = max(C(:));
maxVec = [X(C(:) == maxVal), Y(C(:) == maxVal), Z(C(:) == maxVal)] .* 1.5;
if size(maxVec,1) ==1; maxVec = [maxVec; -maxVec]; end




figure, hold on,
set(gcf, 'position', [256 5 1150 824]);
hand = surf(X,Y,Z,C);
set(hand, 'edgealpha', 0);
plot3([-1.4, 1.4], [-1.4 1.4], [0,0], '-', 'linewidth', 18, 'color', 'k') % the L+M direction
%plot3([-1, 1], [-1 1], [-1 1], '-', 'linewidth', 5, 'color', 'k')           % the achromatic direction
plot3([0, 0], [0 0], [-1.4 1.4], '-', 'linewidth', 18, 'color', 'b')         % S-iso color direction
plot3([1.2, -1.2], [-1.2 1.2], [0 0], '-', 'linewidth', 18, 'color', 'r')            % L-M color direction
%plot3(minVec(:,1), minVec(:,2), minVec(:,3), '--', 'linewidth', 5, 'color', 'm')
plot3(maxVec(:,1), maxVec(:,2), maxVec(:,3), '--', 'linewidth', 18, 'color', [.3 .3 .3])
xlabel('\DeltaL/L'); ylabel('\DeltaM/M'); zlabel('\DeltaS/S')
xlim([-1.3, 1.3]); ylim([-1.3, 1.3]); zlim([-1.3, 1.3]);
set(gca, 'view', [172.2 10.161],...
    'visible', 'off')
axis square
maps = {'Edge', 'CubicL', 'CubicYF', 'IsoL'};
colormap(pmkmp(250, maps{2}))
colorbar

fprintf('Max Vec:   %.2f, %.2f, %.2f \n', maxVec(1,1)/1.5, maxVec(1,2)/1.5, maxVec(1,3)/1.5)
fprintf('Min Vec:   %.2f, %.2f, %.2f \n', minVec(1,1)/1.5, minVec(1,2)/1.5, minVec(1,3)/1.5)


%% ISO-DETECTION RATIO AS A FUNCTION OF SPATIAL FREQ

clc; close all;

% Set the parameters of the analysis
observer = 'kali';              % 'kali' or 'sedna'
theseSFToBootstrap = [1 1 1 1]; % generate bootstrap estimates for only a subset of the SFs 
nstraps = 5000;


% load the data
clear fpar_monkey fpar_bootstrap existingBootstraps existingParams % makes auto detection of existing params work
switch observer
    case 'kali'
        load Kali_DTNT_072113.mat % monkey data
        load Kali_fpar_072113.mat % fits to monkey surfaces 050514 for more straps or 072113 for fewer
        load Kali_boot_072113.mat % bootstraps for kali 050514 or 072113
        existingParams = fpar_monkey; % 'fpar_monkey' is a variable that I use below, but I don't want to overwrite it...
        existingBootstraps = fpar_bootstrap;
    case 'sedna'
        load Sedna_DTNT_072113.mat % monkey data
%         load Sedna_fpar_072113.mat % fits to monkey surfaces
%         load Sedna_boot_072113.mat % bootstraps for sedna
%         existingParams = fpar_monkey;
%         existingBootstraps = fpar_bootstrap;
end


% divide the behavioral data by 100 so that % is b/w 0&1.
monkeyAlphas = cellfun(@(x,y) x./y, dtnt.alpha',  repmat({100}, numel(dtnt.alpha),1), 'uniformoutput', 0);

    
fpar_monkey = {};
for a = 1:4;
    clc
    fprintf('Starting SF %d\n', a)
    
    % which SF condition should be analyzed? The behvioral data has many SFs,
    % but the model data is only at one SF.
    sfIdx = false(1,4);
    sfIdx(a) = true
    sfs = unique(cat(1,dtnt.sfs{:}))
    sfs(sfIdx)
    
    % make a logical vector to weed out conditions (if need be)
    l_questBeta = dtnt.l_beta > 0;
    l_bkgndrgb = dtnt.l_bkgndrgb > 0;
    
    % pull out the spatial frequncy condition
    l_sfs = cat(1,dtnt.sfs{:}) == sfs(sfIdx);
    
    % grab the appropriate files
    l_valid = [l_sfs & l_questBeta & l_bkgndrgb];
    colors = cat(1, dtnt.colorDirs{l_valid});
    alphas = cat(1, monkeyAlphas{l_valid});
    
    % l_goodQuestThresh = cat(1, dtnt.validThresh{l_valid});
    % colors = colors(l_goodQuestThresh,:);
    % alphas = alphas(l_goodQuestThresh);
    size(colors)
    
    
    % only fit the raw data if the cole surface fits do not already exist
    if exist('existingParams', 'var')
        fpar_monkey{a} = existingParams{a};
    else
        
        %%%%%%%%%%
        % default init params on the monkey
        [fpar_monkey_default, fval_monkey_default] = fitDetectionSurface(colors, alphas, 'ellipsoid');
        
        % init params from kali at 0.5 cpd (on the monkey)
        initparams_kali = [1.15; -15.82; 38.29; -14.36; -17.29; 6.62; -8.74; -81.82; 68.40; 9.36];  % optional input argument for threshSurfPlot for 0.5cpd
        [fpar_monkey_kali, fval_monkey_kali] = fitDetectionSurface(colors, alphas, initparams_kali);
        
        % particle swarm on the monkey
        npsoiters = 45;
        [fpar_monkey_pso, fval_monkey_pso] = fitDetectionSurface(colors, alphas, 'pso', npsoiters);
        
        
        % find the best method
        method = {'default'; 'kali'; 'pso'};
        fvals = [fval_monkey_default; fval_monkey_kali; fval_monkey_pso]
        fpars = [fpar_monkey_default'; fpar_monkey_kali'; fpar_monkey_pso];
        [~, idx] = min(fvals);
        fpar_monkey{a} = fpars(idx,:);
        ['winner was ', method{idx}]
    end
    
    
    % generate a bunch of bootstrap versions of the fits (if need be). Start
    % by assuming that the best fitting surface is the real surface and then
    % add residuals randomly
    if theseSFToBootstrap(a)
        if exist('existingBootstraps', 'var')
            fpar_bootstrap{a} = existingBootstraps{a};
        else
            fpar_bootstrap{a} = bootstrapDetectSurf(colors, alphas, fpar_monkey{a}, nstraps);
        end
    end
end


% fit the cone noise model's data
fpar_cones = fitDetectionSurface(gab.colorDirs, cones.alpha_analytic, 'ellipsoid');


% now extract the monkey thresholds and the IDR for a handful of color dirs
colors = [1./sqrt(2), -1./sqrt(2), 0;...
          1./sqrt(2), 1./sqrt(2), 0;...
          1./sqrt(3), 1./sqrt(3), 1./sqrt(3);...
          0, 0, 1;...
          0.1386, -0.1386, -0.9806;...
          0.1386, -0.1386,  0.9806];


% calculate the values for each color (IDR and threshold)
[t_monk, t_cones] = deal([]);
t_bootstrap = {};
for a = 1:4
    
    % first the monkey's thresholds
    mech = reshape(fpar_monkey{a}(2:end),3,3);
    beta = fpar_monkey{a}(1);
    t_monk(:,a) = coleThresh(mech, beta, colors); % <nColors x nSF>
    
    % next for the cone model
    mech = reshape(fpar_cones(2:end),3,3);
    beta = fpar_cones(1);
    t_cones(:,a) = coleThresh(mech, beta, colors);
    
    % grab the estimates based on the boot strap <{nSF}(nColors x
    % nBootstraps)>
    if theseSFToBootstrap(a)

        cell_colors = repmat({colors}, numel(fpar_bootstrap{a}), 1);
        cell_mech = cellfun(@(x) reshape(x(2:end),3,3), fpar_bootstrap{a}, 'uniformoutput', false);
        cell_beta = cellfun(@(x) x(1), fpar_bootstrap{a}, 'uniformoutput', false);
        t_bootstrap{a} = cellfun(@(x,y,z) coleThresh(x, y, z), cell_mech, cell_beta, cell_colors', 'uniformoutput', false);

    end
    
end


% determine the errorbars for each color/sf
sem_thresh = nan(size(colors, 1), numel(sfs)); % <nColors x nSF >
for a = 1:numel(sfs)
    if theseSFToBootstrap(a)
        t_bt = cat(2, t_bootstrap{a}{:}); % each row is a color
        sem_thresh(:,a) = std(t_bt, [], 2); % don't need to divide by sqrt(N-1) b/c this is the distribution of sample values 
    end
end


% plot spatial contrast sensitivity
mkClr = {'r', 'k', 'k', 'b', 'c', 'm'};
faceClr = {'r', 'none', 'k', 'b', 'c', 'm'};
lineSty = {'-', '--', '-', '-', '-', '-',};
plotSF = [0.5 1 2 4];
figure, hold on,
for a = 1:size(colors,1)
    plot(plotSF, t_monk(a,:), [lineSty{a}, 'o', mkClr{a}], 'markerfacecolor', faceClr{a})
end
for a = 1:numel(plotSF);
    if theseSFToBootstrap(a)
        t_bt = cat(2, t_bootstrap{a}{:}); % each row is a color
        for clr = 1:size(colors,1)
            %plot(plotSF(a), t_bt(clr,:), ['o', mkClr{clr}], 'markerfacecolor', faceClr{clr})
            errorbar(plotSF(a), t_monk(clr,a), sem_thresh(clr, a), ['o', mkClr{clr}], 'markerfacecolor', faceClr{clr})
        end
    end
end
set(gca, 'xscale', 'log', 'yscale', 'log')
title('Spatial Contrast thresholds')
xlabel('Spatial Frequency')
ylabel('threshold')
axis tight
xlim([0.47 4.1])


% plot IDR as a function of SF
figure, hold on,
for a = 1:size(colors,1)
    plot([0.5 1 2 4], t_monk(a,:)./t_cones(a,:), [lineSty{a}, 'o', mkClr{a}], 'markerfacecolor', faceClr{a})
end
for a = 1:numel(plotSF);
    if theseSFToBootstrap(a)
        t_bt = cat(2, t_bootstrap{a}{:}); % each row is a color
        for clr = 1:size(colors,1)
            IDR_boot = t_bt(clr,:)./t_cones(clr,a);
            %plot(plotSF(a), t_bt(clr,:)./t_cones(clr,a), ['o', mkClr{clr}], 'markerfacecolor', faceClr{clr})
            errorbar(plotSF(a), t_monk(clr,a)./t_cones(clr,a), std(IDR_boot), ['o', mkClr{clr}], 'markerfacecolor', faceClr{clr})
        end
    end
end
set(gca, 'xscale', 'log', 'yscale', 'log')
title('IDR as a fxn or SF')
xlabel('Spatial Frequency')
ylabel('monk thresh / cone thresh')
axis tight
xlim([0.47 4.1])



%% SIGNIFICANCE TEST ON SPECIFIC IDR CONDITIONS

clc; close all;

% Set the parameters of the analysis
observer = 'sedna';              % 'kali' or 'sedna'
sfIdx = logical([0,1,0,0]); % generate bootstrap estimates for only a subset of the SFs 


% load the data
clear fpar_monkey fpar_bootstrap
switch observer
    case 'kali'
        load Kali_DTNT_072113.mat % monkey data
        load Kali_fpar_072113.mat % fits to monkey surfaces
        load Kali_boot_072113.mat % bootstraps for kali
    case 'sedna'
        load Sedna_DTNT_072113.mat % monkey data
        load Sedna_fpar_072113.mat % fits to monkey surfaces
        load Sedna_boot_072113.mat % bootstraps for sedna
end


% divide the behavioral data by 100 so that % is b/w 0&1.
monkeyAlphas = cellfun(@(x,y) x./y, dtnt.alpha',  repmat({100}, numel(dtnt.alpha),1), 'uniformoutput', 0);


% which SF condition should be analyzed? The behvioral data has many SFs,
% but the model data is only at one SF.
sfs = unique(cat(1,dtnt.sfs{:}))
sfs(sfIdx)
l_sfs = cat(1,dtnt.sfs{:}) == sfs(sfIdx);

colors = cat(1, dtnt.colorDirs{l_sfs});
alphas = cat(1, monkeyAlphas{l_sfs});
size(colors)


% fit the cone noise model's data
fpar_cones = fitDetectionSurface(gab.colorDirs, cones.alpha_analytic, 'ellipsoid');

% define values for the cole surfaces for the monkey.
fpar_monkey = fpar_monkey{sfIdx};
fpar_bootstrap = fpar_bootstrap{sfIdx};

% now extract the monkey thresholds and the IDR for a handful of color dirs
% test_colors = [1 1 0;...
%                1 1 1];
           
test_colors = [0.0636 -0.0636 0.45;...
               0.0636 -0.0636 -0.45];



% calculate the distribution IDR values using the distribution of
% bootstraped monkey thresholds

% first the monkey's bootstrap thresholds
cell_colors = repmat({test_colors}, numel(fpar_bootstrap), 1);
cell_mech = cellfun(@(x) reshape(x(2:end),3,3), fpar_bootstrap, 'uniformoutput', false);
cell_beta = cellfun(@(x) x(1), fpar_bootstrap, 'uniformoutput', false);
t_boot = cellfun(@(x,y,z) coleThresh(x, y, z), cell_mech, cell_beta, cell_colors', 'uniformoutput', false);
t_boot = cat(2, t_boot{:});

% next for the cone model
mech = reshape(fpar_cones(2:end),3,3);
beta = fpar_cones(1);
t_cones = coleThresh(mech, beta, test_colors);

% for the monkey
mech = reshape(fpar_monkey(2:end),3,3);
beta = fpar_monkey(1);
t_monk = coleThresh(mech, beta, test_colors);
real_IDRs = t_monk./t_cones
    
IDR_boot = bsxfun(@rdivide, t_boot, t_cones);
diffVals = (IDR_boot(1,:) - IDR_boot(2,:));
alpha = 0.05;
limits = prctile(diffVals, 100.*[alpha/2, 1-(alpha/2)])
figure, hold on,
hist(diffVals, 40)
plot([limits;limits], [0,400])
h = ~((limits(1)<0) && (limits(2)>0))

%% SPHERICAL HISTOGRAM OF BEST AND WORST IDR


clc; close all;

% Set the parameters of the analysis
observer = 'sedna';              % 'kali' or 'sedna'
sfIdx = logical([1,0,0,0]); % generate bootstrap estimates for only a subset of the SFs 


% load the data
clear fpar_monkey fpar_bootstrap
switch observer
    case 'kali'
        load Kali_boot_072113.mat % bootstraps for kali
    case 'sedna'
        load Sedna_boot_072113.mat % bootstraps for sedna
end


% fit the cone noise model's data
fpar_cones = fitDetectionSurface(gab.colorDirs, cones.alpha_analytic, 'ellipsoid');


% define values for the cole surfaces for the monkey.
fpar_bootstrap = fpar_bootstrap{sfIdx};

% select some colors to use to define the IDR values
[X,Y,Z] = sphere(400);
colors = [X(:), Y(:), Z(:)];
colors = bsxfun(@rdivide, colors, sqrt(sum(colors.^2, 2)));

% do this once
conesThresh = coleThresh(reshape(fpar_cones(2:end),3,3), fpar_cones(1), colors);

[maxIDR_idx, minIDR_idx] = deal(zeros(numel(X(:)), numel(fpar_bootstrap)));
for a = 1:numel(fpar_bootstrap);
    mech = reshape(fpar_bootstrap{a}(2:end),3,3);
    beta = fpar_bootstrap{a}(1);
    monkeyThresh = coleThresh(mech, beta, colors);
    
    IDR = monkeyThresh ./ conesThresh;
    
    % find the min IDR and the max IDR, but do a hack to make sure they're
    % always on the same side as each other 
    idx = find(IDR == min(IDR(:)), 1, 'first');
    minVec = [colors(idx,1), colors(idx,2), colors(idx,3)];
    if minVec(1)<0; minVec = -minVec; end
    [~, newIdx] = max(colors * minVec(:));
    minIDR_idx(newIdx, a) = 1;
    

    idx = find(IDR == max(IDR(:)), 1, 'first');
    maxVec = [colors(idx,1), colors(idx,2), colors(idx,3)];
    if maxVec(1)<0; maxVec = -maxVec; end
    [~, newIdx] = max(colors * maxVec(:));
    maxIDR_idx(newIdx,a) = 1;
    
    if rem(a,100)==0; disp(a);end
end

figure, hold on,
minMask = sum(minIDR_idx,2);
minMask = reshape(minMask, size(X,1), size(X,2));
set(gcf, 'position', [256 5 1150 824]);
hand = surf(X,Y,Z,minMask);
set(hand, 'edgealpha', 0);
axis equal

figure, hold on,
maxMask = sum(maxIDR_idx,2);
maxMask = reshape(maxMask, size(X,1), size(X,2));
set(gcf, 'position', [256 5 1150 824]);
hand = surf(X,Y,Z,maxMask);
set(hand, 'edgealpha', 0);
axis equal


%% COMPARING RETINA, V1, AND BEHAVIOR: THRESHOLD RATIOS AND NEUROMETRIC THRESHOLDS

fin

%
% (1) load in the batch data file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd /Users/charliehass/LabStuff/Huskies/DTcones/Data/
[fname, fpath] = uigetfile();
load([fpath, fname])

%iterate over expts and compile the necessary data
for a = 1:length(out.dat);
    rawV1TRs(a,:) = out.dat(a).c.alpha(:)' ./ out.dat(a).m.alpha(:)'; %transpose so that colors go across rows.
    rawV1NTs(a,:) = out.dat(a).c.alpha(:)';
    rawPTs(a,:) = out.dat(a).m.alpha(:)';
    
    rawRetTRs(a,:) = ret.dat(a).alpha(:)' ./ out.dat(a).m.alpha(:)';
    rawRetNTs(a,:) = ret.dat(a).alpha(:)';
    
    prefCards(a,:) = out.dat(a).prefCard;
    prefInts(a,:) = out.dat(a).prefInt;
    prefIsolum(a,:) = out.dat(a).prefIsolum;
    monkInitials(a,1) = out.fnames{a}{1}(1);
end

%retrieve the indices to elements in the 'errors' matrix.
lt16DTtrialsInd = strmatch('<16DTtrials', out.errorTypes);
orientMismatchInd = strmatch('orientMismatch', out.errorTypes);
posMismatchInd = strmatch('posMismatch', out.errorTypes);
commonExclusions = sum(out.errors(:,[lt16DTtrialsInd, orientMismatchInd, posMismatchInd]),2)>0;


%
% TRs and NTs by color
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%cycle through the expts and pull out the indicies to the TRs for each
%color direction
[sTR.ret,lvmTR.ret,swmTR.ret,swlTR.ret] = deal([]);
[sNT.ret,lvmNT.ret,swmNT.ret,swlNT.ret] = deal([]);
[sTR.v1,lvmTR.v1,swmTR.v1,swlTR.v1] = deal([]);
[sNT.v1,lvmNT.v1,swmNT.v1,swlNT.v1] = deal([]);
[sPT, lvmPT, swmPT, swlPT] = deal([]);
for a = find(~commonExclusions)';
    
    tmpColor = sign(out.dat(a).expt.standColors);
    
    %deal with the cardinal TR
    if ismember([0 0 1], tmpColor, 'rows')
        [~,idx] = ismember([0 0 1], tmpColor, 'rows');
        if ~isnan(rawV1TRs(a,idx'))
            sTR.ret(end+1,1) = rawRetTRs(a,idx');
            sNT.ret(end+1,1) = rawRetNTs(a,idx');
            sTR.v1(end+1,1) = rawV1TRs(a,idx');
            sNT.v1(end+1,1) = rawV1NTs(a,idx');
            sPT(end+1,1) = rawPTs(a,idx');
        end
    elseif ismember([1 -1 0], tmpColor, 'rows')
        [~,idx] = ismember([1 -1 0], tmpColor, 'rows');
        if ~isnan(rawV1TRs(a,idx'))
            lvmTR.ret(end+1,1) = rawRetTRs(a,idx');
            lvmNT.ret(end+1,1) = rawRetNTs(a,idx');
            lvmTR.v1(end+1,1) = rawV1TRs(a,idx');
            lvmNT.v1(end+1,1) = rawV1NTs(a,idx');
            lvmPT(end+1,1) = rawPTs(a,idx');
        end
    end
    
    %deal with the intermediate TR
    if ismember([1 -1 1], tmpColor, 'rows')
        [~,idx] = ismember([1 -1 1], tmpColor, 'rows');
        if ~isnan(rawV1TRs(a,idx'))
            swlTR.ret(end+1,1) = rawRetTRs(a,idx');
            swlNT.ret(end+1,1) = rawRetNTs(a,idx');
            swlTR.v1(end+1,1) = rawV1TRs(a,idx');
            swlNT.v1(end+1,1) = rawV1NTs(a,idx');
            swlPT(end+1,1) = rawPTs(a,idx');
        end
    elseif ismember([1 -1 -1], tmpColor, 'rows')
        [~,idx] = ismember([1 -1 -1], tmpColor, 'rows');
        if ~isnan(rawV1TRs(a,idx'))
            swmTR.ret(end+1,1) = rawRetTRs(a,idx');
            swmNT.ret(end+1,1) = rawRetNTs(a,idx');
            swmTR.v1(end+1,1) = rawV1TRs(a,idx');
            swmNT.v1(end+1,1) = rawV1NTs(a,idx');
            swmPT(end+1,1) = rawPTs(a,idx');
        end
    end
end

%
% PLOTTING ALL THE NT'S AND PT'S
% FROM MONKEY, V1, AND MODEL
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
set(gcf, 'position', [193     8   206   821])

subplot(5,1,1) %monk NTs
monkPT = [sPT(:); lvmPT(:); swmPT(:); swlPT(:)];
monkGroup = char(repmat('Siso', length(sPT),1), repmat('LvM', length(lvmPT),1), repmat('SwM', length(swmPT),1), repmat('SwL', length(swlPT),1));
boxplot(monkPT, monkGroup, 'symbol', '')
ylim([0, 0.2])
title('Monkey PTs')

subplot(5,1,2) % V1 NTs
v1data = [sNT.v1(:); lvmNT.v1(:); swmNT.v1(:); swlNT.v1(:)];
v1group = char(repmat('Siso', length(sNT.v1),1), repmat('LvM', length(lvmNT.v1),1), repmat('SwM', length(swmNT.v1),1), repmat('SwL', length(swlNT.v1),1));
boxplot(v1data, v1group, 'symbol', '')
ylim([0, 0.2])
title('V1 NTs')


subplot(5,1,3) % model NTs
retdata = [sNT.ret(:); lvmNT.ret(:); swmNT.ret(:); swlNT.ret(:)];
retgroup = char(repmat('Siso', length(sNT.ret),1), repmat('LvM', length(lvmNT.ret),1), repmat('SwM', length(swmNT.ret),1), repmat('SwL', length(swlNT.ret),1));
boxplot(retdata, retgroup)
title('Retina NTs')


subplot(5,1,4) % Model to Monk TRs
retdata = [sTR.ret(:); lvmTR.ret(:); swmTR.ret(:); swlTR.ret(:)];
retgroup = char(repmat('Siso', length(sTR.ret),1), repmat('LvM', length(lvmTR.ret),1), repmat('SwM', length(swmTR.ret),1), repmat('SwL', length(swlTR.ret),1));
boxplot(retdata, retgroup)
set(gca, 'yscale', 'log')
title('Retina to Monkey TRs')

subplot(5,1,5) % Model to V1 TRs
v1NTs = [sNT.v1(:); lvmNT.v1(:); swmNT.v1(:); swlNT.v1(:)];
modNTs = [sNT.ret(:); lvmNT.ret(:); swmNT.ret(:); swlNT.ret(:)];
ratio = modNTs./v1NTs;
group = char(repmat('Siso', length(sNT.ret),1), repmat('LvM', length(lvmNT.ret),1), repmat('SwM', length(swmNT.ret),1), repmat('SwL', length(swlNT.ret),1));
boxplot(ratio, group, 'symbol', '')
set(gca, 'yscale', 'log')
title('Retina to V1 TRs')


%
% 2D AND 3D SCATTER PLOTS OF NTs AND PTs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v1_nt = {sNT.v1, lvmNT.v1, swmNT.v1, swlNT.v1};
v1_xbar = cellfun(@mean, v1_nt);
v1_sem = cellfun(@(x) nanstd(x)./sqrt(numel(x)), v1_nt);

ret_nt = {sNT.ret, lvmNT.ret, swmNT.ret, swlNT.ret};
ret_xbar = cellfun(@mean, ret_nt);
ret_sem = cellfun(@(x) nanstd(x)./sqrt(numel(x)), ret_nt);

monk_pt = {sPT, lvmPT, swmPT, swlPT};
monk_xbar = cellfun(@mean, monk_pt);
monk_sem = cellfun(@(x) nanstd(x)./sqrt(numel(x)), monk_pt);

colororder = {'b', 'r', [255./256, 140./256, 0], 'm'};

f = figure; hold on
for a = 1:4;
    for d = 1:numel(v1_nt{a})
        plot3(monk_pt{a}(d), ret_nt{a}(d), v1_nt{a}(d), '.', 'color', colororder{a})
    end
    
    p = plot3(monk_xbar(a), ret_xbar(a), v1_xbar(a), 's');
    set(p, 'markerfacecolor', colororder{a},...
        'markeredgecolor', colororder{a},...
        'markersize', 25)
end
set(gca, 'view', [14 18],...
         'xscale', 'log',...
         'yscale', 'log',...
         'zscale', 'log',...
         'xlim', [1e-2, 1],...
         'ylim', [1e-3 1e-1],...
         'zlim', [5e-3 5e-1],...
         'tickdir', 'out',...
         'linewidth', 2);
xlabel('Monkey');
ylabel('Cones');
zlabel('V1');
axis square

%
% ORTHOGONAL REGRESSION IN 3D  
%
%%%%%%%%%%%%%%%%%%%%
v1NTs = [sNT.v1(:); lvmNT.v1(:); swmNT.v1(:); swlNT.v1(:)];
modNTs = [sNT.ret(:); lvmNT.ret(:); swmNT.ret(:); swlNT.ret(:)];
monkPT = [sPT(:); lvmPT(:); swmPT(:); swlPT(:)];
dat = [monkPT, modNTs, v1NTs];
log_dat = log10(dat);
var_log_dat = var(log_dat,[],1)
tmp = bsxfun(@minus, log_dat, mean(log_dat,1));
covMtx = (tmp' * tmp) ./ (size(log_dat,1)-1);
[vec, val] = eig(covMtx);
[~, idx] = max(diag(val));
PC1 = abs(vec(:,idx));
norms = linspace(-2, 1, 100);
xyz = bsxfun(@times, PC1, norms)';
xyz = bsxfun(@plus, xyz, mean(log_dat,1));
xyz = 10.^xyz;
figure(f); hold on,
plot3(xyz(:,1), xyz(:,2), xyz(:,3), '-k', 'linewidth', 3)


% orth regression for just the means
log_dat = log10([monk_xbar(:), ret_xbar(:), v1_xbar(:)]);
[vec, val] = eig(cov(log_dat));
[~, idx] = max(diag(val));
PC1 = abs(vec(:,idx));
norms = linspace(-2, 1, 100);
xyz = bsxfun(@times, PC1, norms)';
xyz = bsxfun(@plus, xyz, mean(log_dat,1));
diffVal = xyz(end,:)-xyz(1,:)
xyz = 10.^xyz;
figure(f); hold on,
plot3(xyz(:,1), xyz(:,2), xyz(:,3), '-c', 'linewidth', 3)


%
% ORTHOGONAL REGRESSION IN 2D
%
%%%%%%%%%%%%%%%%%%%%
for ax = 1:3;
    if ax == 1
        dat = [ret_xbar(:), monk_xbar(:)];
        semX = ret_sem;
        semY = monk_sem;
    elseif ax == 2
        dat = [ret_xbar(:), v1_xbar(:)];
        semX = ret_sem;
        semY = v1_sem;
    elseif ax == 3
        dat = [monk_xbar(:), v1_xbar(:)];
        semX = monk_sem;
        semY = v1_sem;
    end
    covMtx = cov(dat);
    [vec, val] = eig(covMtx);
    [~, idx] = max(diag(val));
    PC1 = abs(vec(:,idx));
    norms = linspace(-0.2, 0.2, 100);
    xy = bsxfun(@times, PC1, norms)';
    xy = bsxfun(@plus, xy, mean(dat,1));
    diffVal = xy(end,:)-xy(1,:);
    m = diffVal(2)./diffVal(1);
    
    % plot the regression
    f = figure; hold on
    plot(xy(:,1), xy(:,2), 'k', 'linewidth', 3)
    
    % the data points for the mean
    for a = 1:4;
        p = plot(dat(a,1), dat(a,2), 's');
        set(p, 'markerfacecolor', colororder{a},...
               'markeredgecolor', colororder{a},...
               'markersize', 25)
    end
    
    % error bars
    for a = 1:4;
        p = plot([dat(a,1)-semX(a); dat(a,1)+semX(a)], [dat(a,2)', dat(a,2)'], '-');
        set(p, 'color', colororder{a}, 'linewidth', 2);
        
        p = plot([dat(a,1)', dat(a,1)'], [dat(a,2)-semY(a); dat(a,2)+semY(a)], '-');
        set(p, 'color', colororder{a}, 'linewidth', 2);
    end
    
    set(gca, 'xlim', [0, max(dat(:,1)).*1.15],...
             'ylim', [0, max(dat(:,2)).*1.15],...
             'tickdir', 'out',...
             'linewidth', 2);
    if ax == 1 || ax == 2;
        xlabel('Cones');
    else
        xlabel('Monkey')
    end
    if ax == 1
        ylabel('Monkey');
    elseif ax == 2 || ax == 3
        ylabel('V1')
    end
    axis square
    title(sprintf('Slope: %.3f', m));
end



%% COMPARING RETINA, V1, AND BEHAVIOR: Old code
fin

%
% (1) load in the batch data file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd /Users/charliehass/LabStuff/Huskies/DTcones/Data/
[fname, fpath] = uigetfile();
load([fpath, fname])

%iterate over expts and compile the necessary data
for a = 1:length(out.dat);
    rawV1TRs(a,:) = out.dat(a).c.alpha(:)' ./ out.dat(a).m.alpha(:)'; %transpose so that colors go across rows.
    rawV1NTs(a,:) = out.dat(a).c.alpha(:)';
    rawPTs(a,:) = out.dat(a).m.alpha(:)';
    
    rawRetTRs(a,:) = ret.dat(a).alpha(:)' ./ out.dat(a).m.alpha(:)';
    rawRetNTs(a,:) = ret.dat(a).alpha(:)';
    
    prefCards(a,:) = out.dat(a).prefCard;
    prefInts(a,:) = out.dat(a).prefInt;
    prefIsolum(a,:) = out.dat(a).prefIsolum;
    monkInitials(a,1) = out.fnames{a}{1}(1);
    
    
    
    %     figure
    %     subplot(1,2,1)
    %     plot(ret.dat(a).norms{1}, ret.dat(a).roc_analytic{1}, 'b.')
    %     set(gca, 'xscale', 'log', 'ylim', [.47 1.02])
    %     subplot(1,2,2)
    %     plot(ret.dat(a).norms{2}, ret.dat(a).roc_analytic{2}, 'b.')
    %     set(gca, 'xscale', 'log', 'ylim', [.47 1.02])
    %     pause(1)
    %     close
end

%retrieve the indices to elements in the 'errors' matrix.
lt16DTtrialsInd = strmatch('<16DTtrials', out.errorTypes);
orientMismatchInd = strmatch('orientMismatch', out.errorTypes);
posMismatchInd = strmatch('posMismatch', out.errorTypes);
commonExclusions = sum(out.errors(:,[lt16DTtrialsInd, orientMismatchInd, posMismatchInd]),2)>0;
l_sedna = monkInitials == 'S';
l_kali = monkInitials == 'K';


%
% (2) Histograms of V1TRs and RetTRs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
set(gcf, 'position', [92         390        1301         439]);
subplot(1,2,1) % TRs from V1
v1TRs = rawV1TRs(~commonExclusions,:);
v1TRs = min(v1TRs,[],2);
valForNan = max(v1TRs).*2;
v1TRs(isnan(v1TRs)) = valForNan;
edges = logspace(log10(min(v1TRs).*.95), log10(max(v1TRs).*1.1), 25);
counts = histc(v1TRs, edges);
b = bar(edges, counts, 'type', 'histc');
set(gca, 'xscale', 'log', 'tickDir', 'out', 'fontsize', 14)
ch = get(gca, 'children');
set(ch(1), 'Visible', 'off')
title('V1 TRs')

subplot(1,2,2)
retTRs = rawRetTRs(~commonExclusions,:);
retTRs = min(retTRs,[],2);
valForNan = max(retTRs).*2;
retTRs(isnan(retTRs)) = valForNan;
edges = logspace(log10(min(retTRs).*.98), log10(max(retTRs).*1.05), 20);
counts = histc(retTRs, edges);
b = bar(edges, counts, 'type', 'histc');
set(gca, 'xscale', 'log', 'tickDir', 'out', 'fontsize', 14)
ch = get(gca, 'children');
set(ch(1), 'Visible', 'off')
title('Cone Mosaic TRs')

%
% (3) scatter plot of V1TRs and RetTRs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
v1TRs = rawV1TRs(~commonExclusions,:);
retTRs = rawRetTRs(~commonExclusions,:);
plot(v1TRs(:), retTRs(:), '.')
xlabel('V1 TRs')
ylabel('Retina TRs')

%
% (4) TRs and NTs by color
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%cycle through the expts and pull out the indicies to the TRs for each
%color direction
[sTR.ret,lvmTR.ret,swmTR.ret,swlTR.ret] = deal([]);
[sNT.ret,lvmNT.ret,swmNT.ret,swlNT.ret] = deal([]);
[sTR.v1,lvmTR.v1,swmTR.v1,swlTR.v1] = deal([]);
[sNT.v1,lvmNT.v1,swmNT.v1,swlNT.v1] = deal([]);
[sPT, lvmPT, swmPT, swlPT] = deal([]);
for a = find(~commonExclusions)';
    
    tmpColor = sign(out.dat(a).expt.standColors);
    
    %deal with the cardinal TR
    if ismember([0 0 1], tmpColor, 'rows')
        [~,idx] = ismember([0 0 1], tmpColor, 'rows');
        if ~isnan(rawV1TRs(a,idx'))
            sTR.ret(end+1,1) = rawRetTRs(a,idx');
            sNT.ret(end+1,1) = rawRetNTs(a,idx');
            sTR.v1(end+1,1) = rawV1TRs(a,idx');
            sNT.v1(end+1,1) = rawV1NTs(a,idx');
            sPT(end+1,1) = rawPTs(a,idx');
        end
    elseif ismember([1 -1 0], tmpColor, 'rows')
        [~,idx] = ismember([1 -1 0], tmpColor, 'rows');
        if ~isnan(rawV1TRs(a,idx'))
            lvmTR.ret(end+1,1) = rawRetTRs(a,idx');
            lvmNT.ret(end+1,1) = rawRetNTs(a,idx');
            lvmTR.v1(end+1,1) = rawV1TRs(a,idx');
            lvmNT.v1(end+1,1) = rawV1NTs(a,idx');
            lvmPT(end+1,1) = rawPTs(a,idx');
        end
    end
    
    %deal with the intermediate TR
    if ismember([1 -1 1], tmpColor, 'rows')
        [~,idx] = ismember([1 -1 1], tmpColor, 'rows');
        if ~isnan(rawV1TRs(a,idx'))
            swlTR.ret(end+1,1) = rawRetTRs(a,idx');
            swlNT.ret(end+1,1) = rawRetNTs(a,idx');
            swlTR.v1(end+1,1) = rawV1TRs(a,idx');
            swlNT.v1(end+1,1) = rawV1NTs(a,idx');
            swlPT(end+1,1) = rawPTs(a,idx');
        end
    elseif ismember([1 -1 -1], tmpColor, 'rows')
        [~,idx] = ismember([1 -1 -1], tmpColor, 'rows');
        if ~isnan(rawV1TRs(a,idx'))
            swmTR.ret(end+1,1) = rawRetTRs(a,idx');
            swmNT.ret(end+1,1) = rawRetNTs(a,idx');
            swmTR.v1(end+1,1) = rawV1TRs(a,idx');
            swmNT.v1(end+1,1) = rawV1NTs(a,idx');
            swmPT(end+1,1) = rawPTs(a,idx');
        end
    end
end

figure
set(gcf, 'position', [95    32   272   774])

subplot(5,1,1) %monk NTs
monkPT = [sPT(:); lvmPT(:); swmPT(:); swlPT(:)];
monkGroup = char(repmat('Siso', length(sPT),1), repmat('LvM', length(lvmPT),1), repmat('SwM', length(swmPT),1), repmat('SwL', length(swlPT),1));
boxplot(monkPT, monkGroup, 'symbol', '')
ylim([0, 0.2])
title('Monkey PTs')

subplot(5,1,2) % V1 NTs
v1data = [sNT.v1(:); lvmNT.v1(:); swmNT.v1(:); swlNT.v1(:)];
v1group = char(repmat('Siso', length(sNT.v1),1), repmat('LvM', length(lvmNT.v1),1), repmat('SwM', length(swmNT.v1),1), repmat('SwL', length(swlNT.v1),1));
boxplot(v1data, v1group, 'symbol', '')
ylim([0, 0.2])
title('V1 NTs')


subplot(5,1,3) % model NTs
retdata = [sNT.ret(:); lvmNT.ret(:); swmNT.ret(:); swlNT.ret(:)];
retgroup = char(repmat('Siso', length(sNT.ret),1), repmat('LvM', length(lvmNT.ret),1), repmat('SwM', length(swmNT.ret),1), repmat('SwL', length(swlNT.ret),1));
boxplot(retdata, retgroup)
title('Retina NTs')


subplot(5,1,4) % Model to Monk TRs
retdata = [sTR.ret(:); lvmTR.ret(:); swmTR.ret(:); swlTR.ret(:)];
retgroup = char(repmat('Siso', length(sTR.ret),1), repmat('LvM', length(lvmTR.ret),1), repmat('SwM', length(swmTR.ret),1), repmat('SwL', length(swlTR.ret),1));
boxplot(retdata, retgroup)
set(gca, 'yscale', 'log')
title('Retina to Monkey TRs')

subplot(5,1,5) % V1 to model TRs
v1NTs = [sNT.v1(:); lvmNT.v1(:); swmNT.v1(:); swlNT.v1(:)];
modNTs = [sNT.ret(:); lvmNT.ret(:); swmNT.ret(:); swlNT.ret(:)];
ratio = v1NTs./modNTs;
group = char(repmat('Siso', length(sNT.ret),1), repmat('LvM', length(lvmNT.ret),1), repmat('SwM', length(swmNT.ret),1), repmat('SwL', length(swlNT.ret),1));
boxplot(ratio, group, 'symbol', '')
set(gca, 'yscale', 'log')
title('V1 to RetinaTRs')


% determine if there's a main effect of color on the retTRs. Here I'll test
% the log(TR) to try to make a fair comparison
retdata = [sTR.ret(:); lvmTR.ret(:); swmTR.ret(:); swlTR.ret(:)];
retgroup = char(repmat('Siso', length(sTR.ret),1), repmat('LvM', length(lvmTR.ret),1), repmat('SwM', length(swmTR.ret),1), repmat('SwL', length(swlTR.ret),1));
[~, ~, stats] = anova1(log(retdata), retgroup, 'off');
figure; multcompare(stats)
title('Comparison of Retinal TRs')


%
% Look at the ratio of V1NT to RetNT, first separated by color
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
retdata = [sNT.ret(:); lvmNT.ret(:); swmNT.ret(:); swlNT.ret(:)];
v1data = [sNT.v1(:); lvmNT.v1(:); swmNT.v1(:); swlNT.v1(:)];
tmpdata = v1data ./ retdata;
group = char(repmat('Siso', length(sTR.ret),1), repmat('LvM', length(lvmTR.ret),1), repmat('SwM', length(swmTR.ret),1), repmat('SwL', length(swlTR.ret),1));
boxplot(tmpdata, group)
set(gca, 'yscale', 'log')
title('V1 to Retina NT Ratio')

% determine if there is a main effect of color on the V1 to Retina NT ratio
[~, ~, stats] = anova1(log(tmpdata), group, 'off');
figure; multcompare(stats)
title('Comparison of V1 to Retina NT Ratio')




%
% Look at the ratio of V1NT to RetNT, now separated by CSI
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~, idx] = min(rawV1TRs,[],2);
V12RetRatio = nan(numel(idx),1);
for a = 1:numel(idx)
    V12RetRatio(a) = rawV1NTs(a,idx(a)) ./ rawRetNTs(a,idx(a));
end
valForNan = max(V12RetRatio).*1.1;
l_nan = isnan(V12RetRatio);
V12RetRatio(l_nan) = valForNan;
rawCSI = cat(1,out.dat(:).csi);

figure
set(gcf, 'position', [163         466        1124         353]);
subplot(1,2,1)
set(gca, 'fontsize', 18)
plot(rawCSI(~commonExclusions), V12RetRatio(~commonExclusions), '.', 'markersize', 20)
set(gca, 'xscale', 'log', 'yscale', 'log')
xlim([0.01, 100])
ylim([1.5, 40])
[r,p] = corr(rawCSI(~commonExclusions), V12RetRatio(~commonExclusions), 'type', 'spearman');
title(sprintf('CSI vs V1:Ret NT Ratio, p = %.3f', p))
xlabel('CSI')
ylabel('V1:Ret NT Ratio')

subplot(1,2,2)
set(gca, 'fontsize', 18)
minV1TRs = min(rawV1TRs,[],2);
valForNan = max(minV1TRs) .* 1.1;
l_nan = isnan(minV1TRs);
minV1TRs(l_nan) = valForNan;
plot(rawCSI(~commonExclusions), minV1TRs(~commonExclusions), '.', 'markersize', 20)
[r,p] = corr(rawCSI(~commonExclusions), minV1TRs(~commonExclusions), 'type', 'spearman');
set(gca, 'xscale', 'log', 'yscale', 'log')
xlim([ 0.01, 100])
ylim([0.48, 7.2])
title(sprintf('CSI vs V1TR, p = %.3f', p))
xlabel('CSI')
ylabel('V1 TR')


%% NEW CODE FOR LINEAR DISCRIMINANT ANALYSIS
clc, close all

MAKEPLOT = true;

% trying to determing the variance in each cone isolating direction. Is the
% variance and covarience of each data cloud the same? Are the covariences
% zero?
nColors = size(idlob.resp,1);
nContrasts = size(idlob.resp,2);
nTrials = size(idlob.resp,3);

pltIdx = 0;
f2 = figure;
hands = f2;
if MONTECARLO; f1 = figure; hands = [hands, f1]; end
set(hands, 'position', [127 74 1257 755])

rocVals_mc = nan(nColors,nContrasts);
rocVals_anly = nan(nColors,nContrasts);
for a = 1:nColors
    for b = 1:nContrasts
        
        if MONTECARLO % if there are monte carlo trials
            
            % this is a hack to make the 6 wt fxn analysis work
            % with the existing code from the 3 wt fxn analysis.
            dims = size(idlob.resp{a,b,1});
            catdim = find(dims == 1);
            needTranspose = catdim == 2;
            
            % Do LDA on the 3D clouds of data points
            noise = cat(catdim, idlob.resp{a,1,:}); % <nTrials x nCones (3or6)>
            if needTranspose; noise = noise'; end
            mu_noise = mean(noise,1);
            cov_noise = cov(noise);
            
            sig = cat(catdim, idlob.resp{a,b,:});
            if needTranspose; sig = sig'; end
            mu_sig = mean(sig,1);
            cov_sig = cov(sig);
            
            S_within = (cov_noise) + (cov_sig);
            S_between = (mu_noise - mu_sig)' * (mu_noise - mu_sig); % the outer product
            
            [vec, val] = eig(S_within \ S_between); % inv(S_within) * S_between
            [~, maxEigVal] = max(diag(val));
            W_mc = vec(:,maxEigVal);
            [~,maxidx] = max(abs(W_mc));
            if W_mc(maxidx)< 0;
                W_mc = -W_mc; % standardize the representation of eigenvectors
            end
            
            newSig_mc = sig * W_mc;
            newNoise_mc = noise * W_mc;
            
            
            auc_mc = roc(newNoise_mc, newSig_mc);
            if mean(newSig_mc) < mean(newNoise_mc);
                auc_mc = 1-auc_mc;
            end
            rocVals_mc(a,b) = auc_mc;
        end
        
        %
        % the analytic solution
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~params.enableScones
            error('Must enable S-cones')
        end
        
        numcones = numel(idlob.analyticMean{a,1});
        
        % assemble the data. enforce  a particular dimensionality
        % (row, or column vectors)...
        mu_noise = idlob.analyticMean{a,1};
        mu_noise = mu_noise(:)'; % needs to be a row vector
        tmp = idlob.analyticVar{a,1}; % turned into a column vector in the line below
        cov_noise = eye(numcones) .* repmat(tmp(:),1,numcones);
        
        mu_sig = idlob.analyticMean{a,b};
        mu_sig = mu_sig(:)'; % needs to be a row vector
        tmp = idlob.analyticVar{a,b}; % turned into a column vector in the line below
        cov_sig = eye(numcones) .* repmat(tmp(:),1,numcones);
        
        
        S_within = (cov_noise) + (cov_sig);
        S_between = (mu_noise - mu_sig)' * (mu_noise - mu_sig); % the outer product
        
        [vec, val] = eig(S_within \ S_between); % inv(S_within) * S_between
        [~, maxEigVal] = max(diag(val));
        W_anly = vec(:,maxEigVal); %assuming the first column corresponds to the largest eigenvalue
        [~,maxidx] = max(abs(W_anly));
        if W_anly(maxidx)< 0;
            W_anly = -W_anly; % standardize the representation of eigenvectors
        end
        
        newSig_mu = mu_sig * W_anly;
        newSig_SD = sqrt(diag(cov_sig)' * W_anly.^2);
        signal_icdf = norminv([0.001, 0.999], newSig_mu, newSig_SD); % the upper and lower portions of the signal distribution
        
        newNoise_mu = mu_noise * W_anly;
        newNoise_SD = sqrt(diag(cov_noise)' * W_anly.^2);
        noise_icdf = norminv([0.001, 0.999], newNoise_mu, newNoise_SD); % the upper and lower portions of the noise distribution
        
        lowVal = min([noise_icdf(1), signal_icdf(1)]);
        highVal = max([noise_icdf(2), signal_icdf(2)]);
        xx = linspace(lowVal, highVal, 500); % a domain that covers both distributions
        xx = [-inf, xx, inf]; %manually add the p=0 and p=1 condition
        pFA = 1-normcdf(xx, newNoise_mu, newNoise_SD);
        pHit = 1 - normcdf(xx, newSig_mu, newSig_SD);
        auc_anly = -trapz(pFA, pHit);
        if newSig_mu < newNoise_mu;
            auc_anly = 1-auc_anly;
        end
        rocVals_anly(a,b) = auc_anly;
        
        
        %
        % do some plotting
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if MAKEPLOT
            pltIdx = pltIdx+1;
            if MONTECARLO
                figure(f1)
                sc = max(sqrt(sum(sig.^2,2)));
                subplot(nColors, nContrasts, pltIdx), hold on,
                plot3(noise(:,1), noise(:,2), noise(:,3), 'k.')
                plot3(sig(:,1), sig(:,2), sig(:,3), 'r.')
                plot3([-W_mc(1);W_mc(1)].*sc, [-W_mc(2);W_mc(2)].*sc, [-W_mc(3);W_mc(3)].*sc, 'b-', 'linewidth', 3)
                plot3([-W_anly(1);W_anly(1)].*sc, [-W_anly(2);W_anly(2)].*sc, [-W_anly(3);W_anly(3)].*sc, 'g-', 'linewidth', 3)
                plot3(0,0,0, 'gs', 'markerfacecolor', 'g', 'markersize', 3)
                xlabel('Lcc'); ylabel('Mcc'); zlabel('Scc')
                set(gca, 'xtick', [0], 'xticklabel', '0', 'ytick', [0], 'yticklabel', '0','ztick', [0], 'zticklabel', '0')
                set(gca, 'view', [-43 14])
                hold off
            end
            
            figure(f2)
            subplot(nColors, nContrasts, pltIdx), hold on,
            if MONTECARLO
                minval = min([newNoise_mu-newNoise_SD.*3;newSig_mu-newSig_SD.*3]);
                maxval = max([newNoise_mu+newNoise_SD.*3;newSig_mu+newSig_SD.*3]);
                edges = linspace(minval, maxval, 30);
                N = numel(newSig_mc);
                delta = edges(2)-edges(1);
                counts_sig = histc(newSig_mc, edges)./N;
                counts_noise = histc(newNoise_mc, edges)./N;
                bar(edges, counts_noise, 'type', 'histc')
                b2 = bar(edges, counts_sig, 'type', 'histc');
                set(b2, 'edgeColor', 'k', 'faceColor', 'r', 'faceAlpha', 0.5, 'edgealpha', 0.3)
            end
            plot(edges, normpdf(edges, newNoise_mu, newNoise_SD).*delta, 'b', 'linewidth', 2)
            plot(edges, normpdf(edges, newSig_mu, newSig_SD).*delta, 'r', 'linewidth', 2)
            set(gca, 'linewidth', 2)
            axis tight
            
            if MONTECARLO
                title(sprintf('MC: %.3f  ANLY: %.3f', auc_mc, auc_anly));
            else
                title(sprintf('MC: %.3f  ANLY: %.3f', auc_mc, auc_anly));
            end
            hold off
            
        end
        
    end
end

%% TEMPORAL FREQUENCY AND APOLLO G VS B GUN DETECTION

fin
modelData = '/Users/charliehass/LabStuff/Huskies/DTcones/HighTF_Simulation_noScones.txt';
behavioralData = '/Users/charliehass/LabStuff/Huskies/DTcones/ApolloMacPig_CH.txt';
USE_CC_THRESHOLD = true;

% determine what parameters change across all these data files.
% p = paramsCheck(behavioralData);

% The results of the paramsCheck demonstrate that the files are identical
% execpt for:
% RF pos X = [5 10 15 20 25 30 35 40 50 60 70] => in 1/10 dva
% Colors directions are either G or R gun iso.
% Speed is either 15 or 25 Hz

% unpack the behavioral data and reproduce Greg's fig of G&R gun threholds
% with eccentricity (include TF)
eccentricities = 5:5:70; % some of these not sampled by the monkey or model...
gunColors = [0.62671      0.77467     0.084372;...
    0.1067      0.19341       0.9753];
[monk_15.thresh, monk_25.thresh] = deal([]);
[monk_15.colors, monk_25.colors] = deal([]);
monkFiles = fnamesFromTxt2(behavioralData);

missing = [];
for a = 1: numel(monkFiles)
    
    % extract the data
    fpath = findfile(monkFiles{a});
    if isempty(fpath);
        missing = [missing, a];
        continue
    end
    DT = dtobj(fpath);
    [t, c, sfs] = DTquestUnpack(DT, 'mode'); close % close the figure that gets printed by default
    t = t ./ 100; % so that %CC is b/w 0 and 1;
    c = bsxfun(@rdivide, c, sqrt(sum(c.^2,2))); % color dirs as vector norms
    
    % convert to lum contrast
    Mmtx = reshape(DT.sum.exptParams.m_mtx, 3, 3);
    monspect = reshape(DT.sum.exptParams.mon_spect, [], 3);
    bkgnd_rgb = [DT.sum.exptParams.bkgnd_r, DT.sum.exptParams.bkgnd_g, DT.sum.exptParams.bkgnd_b];
    for clr = 1:numel(t);
        cnt_LMS = t(clr) .* c(clr,:);
        t(clr) = luminanceContrast(cnt_LMS, Mmtx, bkgnd_rgb, monspect);
                
        % Trying to get things into CC units. 
        if USE_CC_THRESHOLD
            t(clr) = norm(cnt_LMS);
            warning('using cone contrast units')
        end
    end
    
    % parse the threshold data <nColors x nEccentricites>
    tmpThresh = nan(size(c,1), numel(eccentricities));
    idx = eccentricities == sqrt(DT.sum.exptParams.rf_x.^2 + DT.sum.exptParams.rf_y.^2);
    if ~any(idx); error('Eccentricity not found'); end
    tmpThresh(:, idx) = t;
    
    % determine which temporal frequency was used
    switch unique(DT.trial(:,DT.idx.driftRate))
        case 15
            
            % assign the color directions
            if isempty(monk_15.colors);
                monk_15.colors = c;
            else % check the order of colors
                if ~all(c(:) == monk_15.colors(:)); error('colors did not match'); end
            end
            
            % concatenate the threshold measurements
            monk_15.thresh = cat(3, monk_15.thresh, tmpThresh);
            
        case 25
            
            % assign the color directions
            if isempty(monk_25.colors);
                monk_25.colors = c;
            else % check the order of colors
                if ~all(c(:) == monk_25.colors(:)); error('colors did not match'); end
            end
            
            % concatenate the threshold measurements
            monk_25.thresh = cat(3, monk_25.thresh, tmpThresh);
    end
end

% check that the colors for the 15 and 25Hz data are the same
if ~all(monk_15.colors(:) == monk_25.colors(:)) || any(abs(gunColors(:) - monk_15.colors(:)) > 1e5)
    error('The colors are different between the 15 and 25 Hz data and the indicated gunColors')
end


%
% iterate over the model data (with NO S-cones)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[sim_15.thresh, sim_25.thresh] = deal([]);
retinaFiles = fnamesFromTxt2(modelData);

for a = 1:numel(retinaFiles)
    fpath = findfile(retinaFiles{a}, '/Users/charliehass/LabStuff/Huskies/DTcones/Data', '.mat');
    load(fpath)
    [gab.nSd, gab.driftRate, gab.rf_x, gab.rf_y, gab.sd, gab.sf, gab.theta]
    [idlob, cones] = coneNoiseROC(params, idlob, cones, gab); % do the ROC analysis and fit neurometric fxns
    
    % pull out some calibration info
    load(params.monCalFile)
    monspect = reshape(cal.monSpect, [],3);
    bkgnd_rgb = mon.bkgndrgb;
    
    
    % determine threshold in the Green and Blue gun directions
    l_ggun = [strcmpi(gab.colorDirs(1,:), 'ggun'); strcmpi(gab.colorDirs(2,:), 'ggun')];
    l_bgun = [strcmpi(gab.colorDirs(1,:), 'bgun'); strcmpi(gab.colorDirs(2,:), 'bgun')];
    if ~any(l_ggun) && ~any(l_bgun); error('could not find the data'); end
    GBthreshEstimates = [cones.alpha_analytic(l_ggun) ; cones.alpha_analytic(l_bgun)]; % in "gun contrast" units
    gunIsoDirs = [0 1 0; 0 0 1];
    for clr = 1:2;
        cnt_RGB = GBthreshEstimates(clr) .* gunIsoDirs(clr,:);
        GBthreshEstimates(clr) = luminanceContrast([], [], bkgnd_rgb, monspect, cnt_RGB);
        
        % Trying to get things into CC units. 
        if USE_CC_THRESHOLD
            bkgnd_lms = mon.bkgndlms_Rstar;
            rgbAtThresh = (cnt_RGB+1) .* bkgnd_rgb;
            lmsAtThresh = mon.rgb2Rstar * rgbAtThresh(:);
            GBthreshEstimates(clr) = norm((lmsAtThresh - bkgnd_lms)./bkgnd_lms);
        end
    end
    
    
    % allocate the data according to the eccentricity of the stimulus
    tmpEcc = sqrt(gab.rf_x^2 + gab.rf_y^2);
    idx = eccentricities == tmpEcc;
    if ~any(idx); error('Eccentricity not found'); end
    tmpThresh = nan(size(gunColors,1), numel(eccentricities));
    tmpThresh(:, idx) = GBthreshEstimates;
    
    switch gab.driftRate
        case 15
            sim_15.thresh = cat(3, sim_15.thresh, tmpThresh);
        case 25
            sim_25.thresh = cat(3, sim_25.thresh, tmpThresh);
    end
    
end


%
% Threshold vs. eccentricity for B and G gun lights. 15 and 25 Hz
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltEcc = eccentricities ./ 10;
pltColors = [0 0.5 0; 0 0 1; 0.55 0.8 0.55; 0.8 0.8 1];

figure,
subplot(1,2,1), hold on,
set(gca, 'fontsize', 16)
tmpData = [nanmean(monk_15.thresh,3); nanmean(monk_25.thresh,3)];
for a = 1:size(tmpData,1)
    inds = ~isnan(tmpData(a,:));
    plot(pltEcc(inds), tmpData(a,inds), 'linewidth', 3, 'color', pltColors(a,:))
end
title('Monkey''s Behavioral thresholds')
xlabel('Eccentricity (dva)')
ylabel('Threshold (Lum contrast)')
xlim([0 7.5])
ylim([10^(-2.1) 10^(-0.6)])
diff(log10(get(gca, 'ylim')))
set(gca, 'yscale', 'log')

subplot(1,2,2), hold on,
set(gca, 'fontsize', 16)
tmpData = [nanmean(sim_15.thresh,3); nanmean(sim_25.thresh,3)];
for a = 1:size(tmpData,1)
    inds = ~isnan(tmpData(a,:));
    plot(pltEcc(inds), tmpData(a,inds), '--', 'linewidth', 3, 'color', pltColors(a,:))
end
title('Model Thresholds')
xlabel('Eccentricity (dva)')
ylabel('Threshold (Lum contrast)')
xlim([0 7.5])
ylim([10^(-2.75) 10^(-1.25)])
diff(log10(get(gca, 'ylim')))
set(gca, 'yscale', 'log')


%
% Monkey to model ratio
%
%%%%%%%%%%%%%%%%%%%
tmp_sim = [nanmean(sim_15.thresh,3); nanmean(sim_25.thresh,3)];
tmp_monk = [nanmean(monk_15.thresh,3); nanmean(monk_25.thresh,3)];
tmp_dat = tmp_monk ./ tmp_sim;
%tmp_dat = log10(tmp_dat);
figure, hold on,
for a = 1:size(tmp_sim,1)
    inds = ~isnan(tmp_dat(a,:));
    plot(pltEcc(inds), tmp_dat(a,inds), 'linewidth', 3, 'color', pltColors(a,:))
end
set(gca, 'yscale', 'log')
ylim([1, 5])
xlim([0, 7.5])
ylabel('Monkey to Model ratio')
xlabel('Eccentricity (deg)')



%% UNPACKING A DTNT RUN OF EXPERIMENTS AND SAVING COLORS & PSY-THRESH

% This script will unpack a DTNT run of experiments and extract the
% necessary information to (1) replicate the results on the DTcones
% simulation and (2) plot isodetection surfaces
fin
filename = 'sednaDTNT_ch.txt'; % kaliDTNT_ch or sednaDTNT_ch

switch whoami
    case 'hass_mbp'
        DTNT_txtfile = '~/LabStuff/Huskies/DTcones/';
    case 'nuke'
end


fnames = fnamesFromTxt2([DTNT_txtfile, filename]);

idx = 0;
minTrials = inf;
for a = 1:numel(fnames)
    if strncmp('sf:', fnames{a}, 3)
        continue %skip non-data file txt lines that zack uses for something unknown to me
    end
    
    fprintf('Unpacking file %d of %d: <%s>\n', a, numel(fnames), fnames{a}{1})

    if ismac
        path = findfile(fnames{a}{1}, '/Volumes/NO BACKUP/NexFiles/');
        DT = dtobj(path);
    else
        DT = dtobj(fnames{a}{1});
    end
    
    [tmpAlpha, tmpColor] = DTquestUnpack(DT, 'mode'); close(gcf)        %thresholds in CC b/w 0 & 100%
    [tmpBadThresh, ~] = DTquestUnpack(DT, 'mode', 15, 0.10); close(gcf)  % will return a NaN for bad estimates
    
    
    for clr = 1:size(tmpColor,1)
        tmp = sum(DT.trial(:, DT.idx.colorDir) == clr);
        minTrials = min([minTrials, tmp])
    end
        
    
    
    % extract the relavent parameters
    idx = idx + 1;
    dtnt.fname{idx} = fnames{a};
    dtnt.colorDirs{idx} = bsxfun(@rdivide, tmpColor, sqrt(sum(tmpColor.^2, 2)));
    dtnt.rf_x(idx) = DT.sum.exptParams.rf_x;
    dtnt.rf_y(idx) = DT.sum.exptParams.rf_y;
    dtnt.sigma(idx) = unique(DT.trial(:, DT.idx.gaborSigma));
    dtnt.nSD(idx) = DT.sum.exptParams.flash_size;
    dtnt.theta(idx) = unique(DT.trial(:, DT.idx.gaborTheta));
    dtnt.sfs{idx} = 1 ./ (unique(DT.trial(:, DT.idx.gaborLambda)) ./ DT.sum.exptParams.pixperdeg);
    dtnt.gamma(idx) = unique(DT.trial(:, DT.idx.gaborGamma));
    dtnt.length(idx) = DT.sum.exptParams.flash_length ./ 1000; % in seconds
    dtnt.speed(idx) = unique(DT.trial(:, DT.idx.driftRate));
    dtnt.alpha{idx} = tmpAlpha;
    dtnt.questBeta(idx) = DT.sum.exptParams.quest_beta;
    dtnt.questSigma(idx) = DT.sum.exptParams.quest_sigma;
    dtnt.validThresh{idx} = ~isnan(tmpBadThresh);
    
    
    % keep track of the monitor calibration information
    dtnt.monSpect{idx} = DT.sum.exptParams.mon_spect;
    dtnt.bkgndrgb{idx} = [DT.sum.exptParams.bkgnd_r, DT.sum.exptParams.bkgnd_g, DT.sum.exptParams.bkgnd_b];
    dtnt.pixperdeg(idx) = DT.sum.exptParams.pixperdeg;
    dtnt.frameRate(idx) = DT.sum.exptParams.frame_rate;
    
    % update the list of valid data files.
    dtnt.fnames{idx} = fnames{a};
end

% monitor spectra (sedna's can differ slightly over the data set)
monspd = cat(2, dtnt.monSpect{:});
unique_monspd = unique(monspd', 'rows');
dtnt.l_monspd = zeros(size(monspd,2),1);
for a = 1:size(unique_monspd,1)
    dtnt.l_monspd(ismember(monspd',unique_monspd(a,:), 'rows')) = a;
end
dtnt.monSpect = unique_monspd;

% bkgndrgb  (sedna's can differ slightly over the data set)
bkgndrgb = cat(1, dtnt.bkgndrgb{:});
unique_bkgndrgb = unique(bkgndrgb, 'rows');
dtnt.l_bkgndrgb = zeros(size(bkgndrgb, 1), 1);
for a = 1:size(unique_bkgndrgb,1)
    dtnt.l_bkgndrgb(ismember(bkgndrgb,unique_bkgndrgb(a,:), 'rows')) = a;
end
dtnt.bkgndrgb = unique_bkgndrgb;


% quest beta. (sedna's can be different from expt to expt).
unique_beta = unique(dtnt.questBeta);
dtnt.l_beta = zeros(numel(dtnt.questBeta), 1);
for a = 1:numel(unique_beta)
    dtnt.l_beta(dtnt.questBeta == unique_beta(a)) = a;
end
dtnt.questBeta = unique_beta;

% only select one parameter value for each of these parameters
if any(cellfun(@(x) numel(unique(x)), {dtnt.pixperdeg, dtnt.frameRate, dtnt.rf_x, dtnt.rf_y, dtnt.sigma, dtnt.nSD, dtnt.theta, dtnt.gamma, dtnt.length, dtnt.speed, dtnt.questSigma}) ~= 1); error('Duplicate parameter sets'); end
dtnt.pixperdeg = unique(dtnt.pixperdeg);
dtnt.frameRate = unique(dtnt.frameRate);
dtnt.rf_x = unique(dtnt.rf_x);
dtnt.rf_y = unique(dtnt.rf_y);
dtnt.sigma = unique(dtnt.sigma);
dtnt.nSD = unique(dtnt.nSD);
dtnt.theta = unique(dtnt.theta);
dtnt.gamma = unique(dtnt.gamma);
dtnt.length = unique(dtnt.length);
dtnt.speed = unique(dtnt.speed);
dtnt.questSigma = unique(dtnt.questSigma);


%% TEMPORAL CONTRAST SENSITIVITY (MONKEYS)

% the temporal frequencies are [0 0.2 1 5 25 37.5]

fin

NUMTRIALS = 15;
PERFRANGE = 0.10; %
RG = 1;          % column that has the L-M data
LUM = 2;         % column that has the LUM data
humanFile = '/Users/charliehass/LabStuff/Huskies/NexFiles/Charlie/CharliePsychophysics/Text Files/charlieHighTF.txt';
monkeyFile = '/Users/charliehass/LabStuff/Huskies/NexFiles/Charlie/Nut/NutDTspot_cosyne.txt';
modelFiles = '/Users/charliehass/LabStuff/Huskies/DTcones/tCSF_monk_vs_model.txt';


fnames = fnamesFromTxt2(monkeyFile);
Nexpt = numel(fnames);
thresh_mode = nan(Nexpt, 2);
thresh_alpha = nan(Nexpt, 2);
sf = nan(Nexpt,1);
tf = nan(Nexpt,1);
for a = 1:numel(fnames)
    
    DT = dtobj(fnames{a}{1});
    disp(fnames{a});
    
    [ex_thresh_mode, ex_color, ex_Sfs] = DTquestUnpack(DT, 'mode', NUMTRIALS, PERFRANGE);
    set(gcf, 'name', sprintf('TF = %.3f, Name = %s', unique(DT.trial(:,DT.idx.driftRate)), fnames{a}{1}));
    [ex_thresh_alpha] = DTquestUnpack(DT, 'weibull'); close;
    
    
    % keep track of spatial and temporal frequency
    sf(a) = 1./(unique(DT.trial(:,DT.idx.gaborLambda))./DT.sum.exptParams.pixperdeg);
    
    % find the indicies to the lvm and lum colors
    lvmidx = ismember(sign(ex_color), [1 -1 0], 'rows');
    lumidx = ismember(sign(ex_color), [1 1 0], 'rows');
    if ~any(lvmidx); disp('no l-m data'); end
    if ~any(lumidx); disp('no lum data'); end
    tf(a) = unique(DT.trial(:,DT.idx.driftRate));
    
    if ~isempty(PERFRANGE)
        % look for OOG issues
        if ex_thresh_mode(lvmidx) > 12;
            ex_thresh_mode(lvmidx) = nan;
            ex_thresh_alpha(lvmidx) = nan;
            children = get(gcf, 'children');
            if find(lvmidx);idx = 2;else idx = 4; end
            set(children(idx), 'yColor', 'r', 'linewidth', 2)
            
            % look for performance issues
        elseif sum(isnan(ex_thresh_mode))>0
            nanIdx = isnan(ex_thresh_mode);
            for i = find(nanIdx)'
                children = get(gcf, 'children');
                if i == 1; idx = 3; else idx = 1; end
                set(children(idx), 'yColor', 'b', 'linewidth', 2)
            end
            
            %otherwise close the existing fig.
        else
            close
        end
    end
    close;
    
    % sort the thresholds based on color
    thresh_mode(a,RG) = ex_thresh_mode(lvmidx);
    thresh_alpha(a,RG) = ex_thresh_alpha(lvmidx);
    thresh_mode(a,LUM) = ex_thresh_mode(lumidx);
    thresh_alpha(a,LUM) = ex_thresh_alpha(lumidx);
end


% error checking for SF
if numel(unique(sf)) ~= 1;
    error('There was more than one spatial frequency tested')
end

% now cycle through the temporal frequencies, and calculate the mean
% sensitivity for each color direction.
uniqueTF = unique(tf);
sen_xbar_mode = nan(numel(uniqueTF),2);
sen_xbar_alpha = nan(numel(uniqueTF),2);
sen_sem_mode  = nan(numel(uniqueTF),2);
sen_sem_alpha  = nan(numel(uniqueTF),2);
h = figure;
for a = 1:numel(uniqueTF);
    l_tf = tf == uniqueTF(a);
    sum(~isnan(thresh_mode(l_tf,:)))
    
    % once for the mode
    tmp_sen = 1./(thresh_mode(l_tf, :)./100);
    sen_xbar_mode(a,:) = nanmean(tmp_sen, 1);
    sen_sem_mode(a,:) = nanstd(tmp_sen, 1) ./ sqrt(sum(~isnan(tmp_sen),1));
    
    % again for the weibull
    tmp_sen = 1./(thresh_alpha(l_tf, :)./100);
    sen_xbar_alpha(a,:) = nanmean(tmp_sen, 1);
    sen_sem_alpha(a,:) = nanstd(tmp_sen, 1) ./ sqrt(sum(~isnan(tmp_sen),1));
    
    subplot(numel(uniqueTF),1, a), hold on,
    plot(thresh_mode(l_tf, RG), '-r.')
    plot(thresh_mode(l_tf, LUM), '-k.')
    ylabel(sprintf('TF = %.2f', uniqueTF(a)))
end


%
% Iterate over the model's thresholds
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
simThresh = [];
simTF = [];
retinaFiles = fnamesFromTxt2(modelFiles);

for a = 1:numel(retinaFiles)
    
    fpath = findfile(retinaFiles{a}, '/Users/charliehass/LabStuff/Huskies/DTcones/Data', '.mat');
    load(fpath)
    [gab.nSd, gab.driftRate, gab.rf_x, gab.rf_y, gab.sd, gab.sf, gab.theta]
    [idlob, cones] = coneNoiseROC(params, idlob, cones, gab); % do the ROC analysis and fit neurometric fxns
    
    % pull out some calibration info
    load(params.monCalFile)
    monspect = reshape(cal.monSpect, [],3);
    bkgnd_rgb = mon.bkgndrgb;
    
    % allocate the data according to the eccentricity of the stimulus
    simThresh = cat(1, simThresh, cones.alpha_analytic);
    simTF = cat(1, simTF, gab.driftRate);
end

%sort the model data into ascending TF order
[simTF, idx] = sort(simTF, 'ascend');
simThresh = simThresh(idx);
simSensitivity = 1./(simThresh);
simTF(simTF == 0) = 0.01;

% plot the tCSF
pltTFs = uniqueTF;
pltTFs(pltTFs == 0) = 0.01;
figure, hold on,
errorbar(pltTFs, sen_xbar_mode(:,RG), sen_sem_mode(:,RG), '-r.', 'linewidth', 4)
errorbar(pltTFs, sen_xbar_mode(:,LUM), sen_sem_mode(:,LUM), '-k.', 'linewidth', 4)
errorbar(pltTFs, sen_xbar_alpha(:,RG), sen_sem_alpha(:,RG), '--r.', 'linewidth', 3)
errorbar(pltTFs, sen_xbar_alpha(:,LUM), sen_sem_alpha(:,LUM), '--k.', 'linewidth', 3)
plot(simTF, simSensitivity, '-m', 'linewidth', 3);
plot(simTF, 1./(ones(size(simTF)).*.126), 'r:')
hold off
set(gca, 'xscale', 'log', 'yscale', 'log', 'box', 'off')
set(gca, 'fontname', 'helvetica', 'fontsize', 30, 'linewidth', 1)
set(gca, 'xlim', [min(pltTFs).*.7, max(pltTFs).*1.4])
xlabel('Temporal Frequency (Hz)')
ylabel('Sensitivity')


% plot the ratio of model to monkey/human
figure, hold on,
plot(pltTFs, simSensitivity ./ sen_xbar_mode(:,RG), '-ro', 'linewidth', 3)
plot(pltTFs, simSensitivity ./ sen_xbar_mode(:,LUM), '-ko', 'linewidth', 3)
set(gca, 'xscale', 'log', 'yscale', 'log', 'box', 'off')
set(gca, 'fontname', 'helvetica', 'fontsize', 16)
xlim([min(pltTFs).*.7, max(pltTFs).*1.4])
xlabel('Temporal Frequency (Hz)')
ylabel('Model to Monkey Ratio')




%
% FIGURE FOR COSYNE TALK
%
%%%%%%%%%%%%%%%%%%%%%%%%%

% only plot data with TF < 37.5
figure,
l_LT37 = pltTFs < 37;
ax = axes; hold on,
myerrorbar(pltTFs(l_LT37), sen_xbar_mode(l_LT37,RG), sen_sem_mode(l_LT37,RG), sen_sem_mode(l_LT37,RG),...
    'linewidth', 4, 'marker', 'o', 'markersize', 10, 'markerfacecolor', 'r', 'markeredgecolor', 'r', 'color', 'r');
myerrorbar(pltTFs(l_LT37), sen_xbar_mode(l_LT37,LUM), sen_sem_mode(l_LT37,LUM), sen_sem_mode(l_LT37,LUM),...
    'linewidth', 4, 'marker', 'o',  'markersize', 10, 'markerfacecolor', 'k', 'markeredgecolor', 'k', 'color', 'k');
[h_yy, ~, h_line2] = plotyy(nan, nan, simTF(l_LT37), simSensitivity(l_LT37));
hold off

% Things common to both Y axes
set([h_yy, ax],...
    'xscale', 'log',...
    'yscale', 'log',...
    'box', 'off', ...
    'fontname', 'helvetica',...
    'fontsize', 30,...
    'linewidth', 1,...
    'xlim', [min(pltTFs).*.7, max(pltTFs).*1.4])

% Things for only the R yaxis
set(h_yy(2),...
    'ycolor', 'b',...
    'linewidth', 2,...
    'ylim', [48, 480],...
    'ytick', linspace(47, 470, 10))

set(h_line2,...
    'color', 'b',...
    'linewidth', 4,...
    'marker', 'o',...
    'markersize', 10,...
    'markerfacecolor', 'b',...
    'markeredgecolor', 'b');

% Things for only the L y axis
set([ax, h_yy(1)],...
    'ylim', [10, 100],...
    'ycolor', 'k',...
    'linewidth', 2,...
    'ytick', linspace(10, 100, 10))

% Figure labels
xlabel('Temporal Frequency (Hz)')
axes(ax), ylabel('Sensitivity (monkey)')
axes(h_yy(2)), ylabel('Sensitivity (cone model)')
diff(log10(get(h_yy(1), 'ylim')))
diff(log10(get(h_yy(2), 'ylim')))



%% TEMPORAL CONTRAST SENSITIVITY (MODEL, SIG VS NOISE)

fin

prefix ='/Users/charliehass/LabStuff/DTcones/';
txtFiles = {
%             {'tCSF_sig_vs_noise_fullModel.txt',             ' full model',              'k'};...
%             {'tCSF_sig_vs_noise_IRFdelta.txt',              ' IRF delta',               'y'};...
%             {'tCSF_sig_vs_noise_flatPS.txt',                ' Flat PS',                 'c'};...
%             {'tCSF_sig_vs_noise_fullModel_noLatency.txt',   ' no latency',              'm'};...
            {'coneRate_825_monRate_825.txt',                ' m=825, c=825',            'b'};...
            {'coneRate_825_monRate_825_phaseInvariant.txt', ' m=825, c=825 6D',         'b--'};...
            {'coneRate_825_monRate_825_filteredWtFxn.txt',  ' m=825, c=825 fltWtFxn',   'b:'};...
            {'coneRate_825_monRate_75.txt',                 ' m=75, c=825',             'g'};...
            {'coneRate_825_monRate_75_phaseInvariant.txt',  ' m=75, c=825 6D',          'g--'};...
            {'coneRate_825_monRate_75_filteredWtFxn.txt',   ' m=75, c=825 fltWtFxn',    'g:'};...
            };



   
% Iterate over the model's thresholds
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ex = 1:numel(txtFiles)
    
    retFiles = fnamesFromTxt2([prefix, txtFiles{ex}{1}]);
    
    for a = 1:numel(retFiles)
        
        fpath = findfile(retFiles{a}, '/Users/charliehass/LabStuff/DTcones/Data', '.mat');
        
        load(fpath)
        fprintf(params.notes)
        [gab.nSd, gab.driftRate, gab.rf_x, gab.rf_y, gab.sd, gab.sf, gab.theta]
        [idlob, cones] = coneNoiseROC(params, idlob, cones, gab); % do the ROC analysis and fit neurometric fxns
        
        % allocate the data according to the eccentricity of the stimulus
        modThresh{ex}(a) = cones.alpha_analytic;
        modTF{ex}(a) = gab.driftRate;
        %
        %         figure, hold on,
        %         plot(gab.contrasts{1}, idlob.roc_analytic{1}, '.k-')
        %         plot(cones.alpha_analytic, 0.816, 'bs')
        %         set(gca, 'xscale', 'log')
        %         title(sprintf('tf: %.3f', gab.driftRate))
        
        
    end
    
    
    %sort the model data into ascending TF order
    [modTF{ex}, idx] = sort(modTF{ex}, 'ascend');
    modThresh{ex} = modThresh{ex}(idx);
    modSensitivity{ex} = 1./(modThresh{ex});
    modTF{ex}(modTF{ex} == 0) = 0.01;
    
end


figure, hold on,
for ex = 1:numel(txtFiles)
    tmp = modSensitivity{ex};
    %tmp = tmp ./ tmp(1); %normalize by the value at DC
    plot(modTF{ex}, tmp, txtFiles{ex}{3}, 'linewidth', 3)
    legtxt{ex} = txtFiles{ex}{2};
end
legend(legtxt)
set(gca, 'xscale', 'log',...
         'yscale', 'log',...
         'xlim', [1 43],...
         'ylim', [0.03 1.3])

%%%
%%% THIS CODE IS A BIT BUSTED. IT'S SUPPOSED TO COMAPRE TCSF FROM DIFFERENT
%%% CONDITIONS, BUT THERE'S NO DYNAMIC ALLOCATION OF THE INDEX TO THE "FULL
%%% MODEL" DATA SET, AND IT CRASES IF THE NUMBER OF TFS IS DIFFERENT ACROSS
%%% CONDITIONS...
%%%
%
% figure, hold on,
% fullModel = modSensitivity{1};
% fullModel = fullModel ./ fullModel(1); %normalize by the value at DC
% for ex = 1:numel(txtFiles)
%     tmp = modSensitivity{ex};
%     tmp = tmp ./ tmp(1); %normalize by the value at DC
%     plot(modTF{ex}(1:35), tmp(1:35)./fullModel(1:35), '-', 'color', clrOrder{ex}, 'linewidth', 2)
% end
% %legend('Full Model', 'Delta IRF', 'Noise PS Flat', 'Full Model No Latency')
% set(gca, 'xscale', 'log',...
%     'yscale', 'log',...
%     'xlim', [0.3 43],...
%     'ylim', [0.02 1.3])

%% ABS THRESH AS A FUNCTION OF NUM CONES

fin
txtFilePrefix = '/Users/charliehass/LabStuff/DTcones/';
exptFile = {'absThreshVsNumCones_dark.txt',...
    'absThreshVsNumCones_5000.txt',...
    'absThreshVsNumCones_10000.txt'};

for ex = 1:numel(exptFile)
    fnames = fnamesFromTxt2([txtFilePrefix, exptFile{ex}])
    
    alpha{ex} = nan(numel(fnames),1);
    numCones{ex} = nan(numel(fnames),1);
    for a = 1:numel(fnames)
        
        fpath = findfile(fnames{a}, '/Users/charliehass/LabStuff/DTcones/Data', '.mat');
        load(fpath)
        [idlob, cones] = coneNoiseROC(params, idlob, cones, gab); % do the ROC analysis and fit neurometric fxns
        
        % compile the data
        alpha{ex}(a) = cones.alpha_analytic;
        numCones{ex}(a) = params.Ncones;
    end
end

figure, hold on,
plot(numCones{1}, alpha{1}, 'k', 'linewidth', 3)
plot(numCones{2}, alpha{2}, 'k', 'linewidth', 3, 'color', [.4 .4 .4])
plot(numCones{3}, alpha{3}, 'k', 'linewidth', 3, 'color', [.7 .7 .7])
plot([14,14], [2,5], 'b', 'linewidth', 5) % add koenig's estimates
set(gca, 'box', 'off',...
    'xscale', 'log',...
    'yscale', 'log',...
    'ylim', [0.5 100],...
    'fontsize', 14)
xlabel('Number of cones')
ylabel('Threshold (R*/cone)')
legend('dark adapted', '5000 R*', '10000 R*', 'Koenig estimate')
xlim([min(numCones{3}).*.9 max(numCones{3}).*1.1])
