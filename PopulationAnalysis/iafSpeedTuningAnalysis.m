%% LOAD THE DATA


cd('C:\Users\charlie\Desktop')
[filename, path] = uigetfile({'*.mat',...
                               'Related Files (*.mat)'},...
                               'Select the Related MATLAB file');
                           
loadpath = [path, filesep, filename];
load(loadpath)


% initialize the array that agregates data across VAs
num_vas = numel(udat.ROI.VisArea);
for i_va = 1:num_vas
    pixelMatrix{i_va} = {};
end

%%%%% AvgImage
for i_va = 1:num_vas % for every Visual Area 'available' in AvgImage
    
    if any(isnan(udat.ROI.RCidx{i_va}))
        continue % don't bother analyzing ROIs that haven't been defined
    end
    
    % determine the row and column coordinates of the pixels within an ROI
    Npixels = size(udat.ROI.RCidx{i_va}, 1); % Determine # of pixels in one selected area
    columns = udat.ROI.RCidx{i_va}(:,2)';
    rows = udat.ROI.RCidx{i_va}(:,1)';
        
    % convert R/C indicies to linear indicies. Note:
    % size(udat.final_img{1}) == size(udat.AvgTrialON{1}(:,:,1))
    linIdx_roi = sub2ind(size(udat.final_img{1}), rows, columns); %Placement of 'important' pixels for the first frame
    
    % linIdx is the index of each ROI pixel on the 1'st frame. Now we need
    % to slove for the linear index of all pixel across all frames. Do this
    % separately for the ON and OFF movies
    NframesON = size(udat.AvgTrialON{1}, 3);
    NframesOFF =  size(udat.AvgTrialOFF{1}, 3);
 
    % Use the linear index to extract the dfof data from the processed data
    Nttypes = size(udat.ttypes, 1); %Nttypes = Total # of Stimulus Types (Trial Types)
    for i_ttypes = 1 : Nttypes
        
        % first grab the dfof data for the ON period. Reshape the array to
        % [Nframes x NPixels]
        roi_on = udat.AvgTrialON{i_ttypes};
        roi_on = reshape(roi_on, [], NframesON); % npix x nframes
        roi_on = roi_on(linIdx_roi,:);
        roi_on = roi_on'; % transpose to make dims = [nframes x npix]
        
        % now grab the dfof data for the OFF period. Reshape the array to
        % [Nframes x NPixels]
        roi_off = udat.AvgTrialOFF{i_ttypes};
        roi_off = reshape(roi_off, [], NframesOFF);
        roi_off = roi_off(linIdx_roi,:);
        roi_off = roi_off';
        
        % concatenate the on and off portions
        pixelMatrix{i_va}{i_ttypes} = cat(1, roi_on, roi_off);
    end
    
end



%% PLOT THE ROIs AS A SANITY CHECK

figure; hold on,
plotimg = cat(3, udat.final_img{:});
plotimg = nanmean(plotimg, 3);
plotimg = plotimg + (abs(min(plotimg(:))));
plotimg = plotimg ./ (max(plotimg(:)).*1.2);
plotimg = repmat(plotimg, [1,1,3]);
for i_roi = 1:numel(udat.ROI.RCidx)
    
    rc = udat.ROI.RCidx{i_roi};
    if ~isempty(rc) && ~any(isnan(rc(:)))
        linind = sub2ind(size(plotimg(:,:,1)), rc(:,1), rc(:,2));
        plotimg(linind) = 0;
    end
end

h_img = imshow(plotimg);


%% ANALYZE THE dFoF AND PLOT AVERAGE TIMESERIES

MAKEPLOT = true;

Nttypes = size(udat.ttypes, 1);

popmean.on = {};
popmean.off = {};

% now the main plotting routines
for i_va = 1:num_vas
    
    popmean.on{i_va} = [];
    popmean.off{i_va} = [];
    
    % Make sure that there was an ROI defined, and that there is data
    % avaliable
    if isempty(pixelMatrix{i_va});
        continue
    end
    
    for i_ttypes = 1: Nttypes
        
        % calculate the average dFoF during the On and Off periods
        meanON = pixelMatrix{i_va}{i_ttypes}(1:NframesON,:);
        meanON = nanmean(meanON(:)); % identical to taking nanmean across time first, then nanmean across pix.
        popmean.on{i_va}(i_ttypes) = meanON;
        
        meanOFF = pixelMatrix{i_va}{i_ttypes}(NframesON+1:end,:);
        meanOFF = nanmean(meanOFF(:)); % identical to taking nanmean across time first, then nanmean across pix.
        popmean.off{i_va}(i_ttypes) = meanOFF;
        
    end
    
end


%
% PLOT THE AVERAGE TIMESERIES
%
if MAKEPLOT
    
    fig = figure;
    fig.Units = 'normalized';
    fig.Position = [0.0757    0.0422    0.8187    0.8344];
    figName = regexpi(filename, '(\w\d+)', 'match');
    figName = regexpi(figName{1}, '(\d)', 'match');
    figName = ['k', figName{end-1}, figName{end}];
    set(gcf, 'name', sprintf('%s', figName), 'numbertitle', 'off')
    ymax = -inf;
    ymin = inf;
    hs = [];
    
    % add the titles to the first row
    for i_ttypes = 1: Nttypes
        subplot(num_vas, Nttypes, i_ttypes)
        switch udat.experimentType
            case 'sfvstf'
                title(sprintf('SF: %0.2f, TF: %d', udat.ttypes(i_ttypes, 1), udat.ttypes(i_ttypes, 2)), 'FontSize', 10);
            case 'size'
                title(sprintf('Size: %d dva', udat.ttypes(i_ttypes)))
        end
    end
    
    % now the main plotting routines
    for i_va = 1:num_vas
        
        % Make sure that there was an ROI defined
        if isempty(pixelMatrix{i_va});
            continue
        end
        
        for i_ttypes = 1: Nttypes
            
            % Plot the ROI pixels according to stimulus types
            pltIdx = (i_va-1)*Nttypes + i_ttypes;
            hs(end+1) = subplot(num_vas, Nttypes, pltIdx); hold on
            
            Nframes = size(pixelMatrix{i_va}{i_ttypes}, 1);
            tt = [0:Nframes-1] .* (1/udat.frameRate);
            plot(tt, nanmean(pixelMatrix{i_va}{i_ttypes},2))
            axis tight
            
            ttON = (NframesON).* (1/udat.frameRate);
            plot([ttON  ttON], ylim, 'LineStyle', ':', 'Color', 'k') % Mark 'Stimulus Offset'
            
            % add a line for the mean across time
            meanON = popmean.on{i_va}(i_ttypes);
            plot([0, ttON], [meanON, meanON], 'k-', 'linewidth', 2)
            
            meanOFF = popmean.off{i_va}(i_ttypes);
            plot([ttON, tt(end)], [meanOFF, meanOFF], 'k-', 'linewidth', 2)
            
            if i_va == num_vas
                xlabel('Time (sec)')
            else
                set(gca, 'xticklabel', [])
            end
            if i_ttypes == 1
                ylabel(sprintf('%s dF/F', udat.ROI.VisArea{i_va}))
            else
                %set(gca, 'yticklabel', [])
            end
            
            % update the ymin/max
            ylims = get(gca, 'ylim');
            ymin = min([ymin, ylims(1)]);
            ymax = max([ymax, ylims(2)]);
        end
        
    end
end
set(hs, 'ylim', [ymin, ymax]) % standardize the axes


%% JOINT SF-TF TUNING AND 1D SPEED TUNING

assert(strcmpi(udat.experimentType, 'sfvstf'), 'ERROR: not a speed tuning expt')


sf_cpd = udat.ttypes(:,strcmpi(udat.text, 'tGratingSpatialFreqCPD'));
tf_cps = udat.ttypes(:,strcmpi(udat.text, 'tGratingTemporalFreqCPS'));
speed_dps = tf_cps ./ sf_cpd;
unique_tf = unique(tf_cps);
unique_sf = unique(sf_cpd);
unique_speed = unique(speed_dps);

popmean.sftf = {};
popmean.speed = repmat({repmat({[]}, 1, numel(unique_speed))}, 1, num_vas); %ordered according to unique_speed

for i_va = 1:num_vas
    if isempty(popmean.on{i_va})
        popmean.sftf{i_va} = [];
        continue
    end
    
    popmean.sftf{i_va} = nan(numel(unique_sf), numel(unique_tf));
    for i_ttype = 1:numel(popmean.on{i_va})
        
        % grab the data
        tmpdat = popmean.on{i_va}(i_ttype);
        
        % store the data for the SF/TF joint tuning matrix
        ridx = unique_sf == sf_cpd(i_ttype);
        cidx = unique_tf == tf_cps(i_ttype);
        popmean.sftf{i_va}(ridx, cidx) = tmpdat;
        
        % store data for the 1D speed tuning functions
        speed_idx = unique_speed == speed_dps(i_ttype);
        popmean.speed{i_va}{speed_idx}(end+1) = tmpdat;
        
    end
    
    % filp the matrix ud so that the SF values go from little to big from
    % the bottom left corner
    popmean.sftf{i_va} = flipud(popmean.sftf{i_va});
    
end


% plot the joint SF-TF tuning matricies
f = figure;
f.Units = 'normalized';
f.Position = [0.1507    0.0578    0.6729    0.8022];
plotdims = ceil(sqrt(num_vas));
for i_va = 1:num_vas
    
    if isempty(popmean.on{i_va})
        continue
    end    
    
    % plot
    subplot(plotdims, plotdims, i_va)
    imagesc(popmean.sftf{i_va});% ./ max(abs(popdat.sftf{i_va}(:))))
    colormap gray; colorbar
    h = gca;
    h.YTick = [];
    h.XTick = [];
    title(udat.ROI.VisArea{i_va}, 'fontsize', 10)
    axis square
end


% plot the 1D speed tuning function
f = figure;
f.Units = 'normalized';
f.Position = [0.1507    0.0578    0.6729    0.8022];
plotdims = ceil(sqrt(num_vas));
for i_va = 1:num_vas
    
    if isempty(popmean.on{i_va})
        continue
    end    
    
    % grab data
    tmpdat = popmean.speed{i_va};
    
    % plot the raw data points
    subplot(plotdims, plotdims, i_va); hold on,
    for i_speed = 1:numel(tmpdat)
        plot(unique_speed(i_speed), tmpdat{i_speed}, 'ko')
    end
    
    % plot the average
    avgdat = cellfun(@nanmean, tmpdat);
    plot(unique_speed, avgdat, 'k', 'linewidth', 2);
    h = gca;
    h.XScale = 'log';
    title(udat.ROI.VisArea{i_va}, 'fontsize', 10)
    xlabel('Speed (dps)')
    ylabel('dFoF ON')
    
end

% now plot them all on the same axis
f = figure; hold on,
hva_color = lines(num_vas);
legtext = {};
for i_hva = 1:num_vas
    
    if isempty(popmean.on{i_hva})
        continue
    end    
    
    % grab data
    tmpdat = popmean.speed{i_hva};
    
    % plot the average
    avgdat = cellfun(@nanmean, tmpdat);
    plot(unique_speed, avgdat, '-', 'color', hva_color(i_hva, :), 'linewidth', 2);
    
    legtext{end+1} = udat.ROI.VisArea{i_hva};
end
h = gca;
h.XScale = 'log';
xlabel('Speed (dps)')
ylabel('dFoF ON')
legend(legtext)

%% SIZE TUNING

assert(strcmpi(udat.text, 'tGratingDiameterDeg') && strcmpi(udat.experimentType, 'size'), 'ERROR: not a size tuning expt')

f = figure;
f.Units = 'normalized';
f.Position = [0.0703    0.5204    0.8740    0.2222];
for i_va = 1:num_vas
    
    if isempty(popmean.on{i_va})
        continue
    end    
    
    % grab data
    tmpdat = popmean.on{i_va};
    
    % plot the raw data points
    subplot(1, num_vas, i_va); hold on,
    plot(udat.ttypes, tmpdat, '-ko')
    
    h = gca;
    title(udat.ROI.VisArea{i_va}, 'fontsize', 10)
    xlabel('Size (dva)')
    ylabel('Fon')
    
end


%% Surround Suppression

assert(strcmpi(udat.text{1}, 'tAnnulusGratingDiameterDeg') && strcmpi(udat.experimentType, 'annulus'), 'ERROR: not a size tuning expt')

f = figure;
f.Units = 'normalized';
f.Position = [0.0703    0.5204    0.8740    0.2222];
for i_va = 1:num_vas
    
    if isempty(popmean.on{i_va})
        continue
    end    
    
    % grab data
    tmpdat = popmean.on{i_va};
    
    % plot the raw data points
    subplot(1, num_vas, i_va); hold on,
    plot(udat.ttypes(:,1), tmpdat, '-ko')
    
    h = gca;
    title(udat.ROI.VisArea{i_va}, 'fontsize', 10)
    xlabel('Annulus Size (dva)')
    ylabel('Fon')
    
end



