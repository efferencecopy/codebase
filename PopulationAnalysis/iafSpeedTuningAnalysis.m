%% LOAD THE DATA

fin

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
    columns = udat.ROI.RCidx{i_va}(:,1)';
    rows = udat.ROI.RCidx{i_va}(:,2)';
        
    % convert R/C indicies to linear indicies. Note:
    % size(udat.final_img{1}) == size(udat.preProcessed_ON{1}(:,:,1))
    linIdx = sub2ind(size(udat.final_img{1}), rows, columns); %Placement of 'important' pixels for the first frame
    
    % linIdx is the index of each ROI pixel on the 1'st frame. Now we need
    % to slove for the linear index of all pixel across all frames.
    Nframes = size(udat.preProcessed_ON{1}, 3);
    linIdx = repmat(linIdx, Nframes, 1); % lists "placement" Nframes times (not taking into account change in frame); rows = # frames & columns = Npixels

    Nelements = numel(udat.final_img{1}); % Determines the total number of elements in the matrix
    frameAdjuster = repmat(Nelements, 1, Npixels);  % length of vector = Npixels
    frameAdjuster = bsxfun(@times, repmat(frameAdjuster, Nframes, 1), [0:1:Nframes-1]');

    newIdx = linIdx + frameAdjuster; % applies adjuster to the "placement"
    newIdx = newIdx(:); % lists every idx# (in a column) SANITY CHECK: rows of newIdx should = Npixels*Nframes

    % Use the linear index to extract the dfof data from the processed data
    Nttypes = size(udat.ttypes, 1); %Nttypes = Total # of Stimulus Types (Trial Types)
    for i_ttypes = 1 : Nttypes
        
        % first grab the dfof data for the ON period. Reshape the array to
        % [Nframes x NPixels]
        roi_on = udat.preProcessed_ON{i_ttypes}(newIdx);
        roi_on = reshape(roi_on, Nframes, []);
        
        % now grab the dfof data for the OFF period. Reshape the array to
        % [Nframes x NPixels]
        roi_off = udat.preProcessed_OFF{i_ttypes}(newIdx);
        roi_off = reshape(roi_off, Nframes, []);
        
        % concatenate the on and off portions
        pixelMatrix{i_va}{i_ttypes} = cat(1, roi_on, roi_off);
    end
    
end



%% PLOT THE ROIs AS A SANITY CHECK

figure; hold on,
plotimg = cat(3, udat.final_img{:});
plotimg = mean(plotimg, 3);
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


%% PLOT: PIXELS VS. TIME
Nttypes = size(udat.ttypes, 1);
NframesON = size(udat.preProcessed_ON{1}, 3);

fig = figure;
fig.Position = [189    26   668   751];

figName = regexpi(filename, '(\w\d+)', 'match');
figName = regexpi(figName{1}, '(\d)', 'match');
figName = ['k', figName{end-1}, figName{end}];
set(gcf, 'name', sprintf('%s', figName), 'numbertitle', 'off')

% add the titles to the first row
for i_ttypes = 1: Nttypes
    subplot(num_vas, Nttypes, i_ttypes)
    title(sprintf('SF: %0.2f, TF: %d', udat.ttypes(i_ttypes, 1), udat.ttypes(i_ttypes, 2)), 'FontSize', 10);
end

popdat.on = {};
popdat.off = {};

% now the main plotting routines
for i_va = 1:num_vas
    
    popdat.on{i_va} = [];
    popdat.off{i_va} = [];
    
    % Make sure that there was an ROI defined, and that there is data
    % avaliable
    if isempty(pixelMatrix{i_va});
        continue
    end
    
    for i_ttypes = 1: Nttypes
        
        % Plot the ROI pixels according to stimulus types
        pltIdx = (i_va-1)*Nttypes + i_ttypes;
        subplot(num_vas, Nttypes, pltIdx); hold on
        
        Nframes = size(pixelMatrix{i_va}{i_ttypes}, 1);
        tt = [0:Nframes-1] .* (1/udat.frameRate);
        plot(tt, mean(pixelMatrix{i_va}{i_ttypes}, 2))
        axis tight
        %ylim([-0.007 0.007]);        
        
        ttON = (NframesON-1).* (1/udat.frameRate);
        line([ttON  ttON], ylim, 'LineStyle', ':', 'Color', 'k') % Mark 'Stimulus Offset'
        
        % add a line for the mean across time
        meanON = pixelMatrix{i_va}{i_ttypes}(1:NframesON,:);
        meanON = mean(meanON(:)); % identical to taking mean across time first, then mean across pix.
        plot([0, ttON], [meanON, meanON], 'k-', 'linewidth', 2)
        
        meanOFF = pixelMatrix{i_va}{i_ttypes}(NframesON+1:end,:);
        meanOFF = mean(meanOFF(:)); % identical to taking mean across time first, then mean across pix.
        plot([ttON, tt(end)], [meanOFF, meanOFF], 'k-', 'linewidth', 2)
        
        
        popdat.on{i_va}(i_ttypes) = meanON;
        popdat.off{i_va}(i_ttypes) = meanOFF;
        
        
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
    end
    
end


%% MAKE JOINT SF-TF TUNING PLOTS

sf_cpd = udat.ttypes(:,strcmpi(udat.text, 'tGratingSpatialFreqCPD'));
tf_cps = udat.ttypes(:,strcmpi(udat.text, 'tGratingTemporalFreqCPS'));
unique_tf = unique(tf_cps);
unique_sf = unique(sf_cpd);

popdat.speedMtx = {};
figure
set(gcf, 'position', [259    47   969   711])
plotdims = ceil(sqrt(num_vas));
for i_va = 1:num_vas
    if isempty(popdat.on{i_va})
        popdat.speedMtx{i_va} = [];
        continue
    end
    
    popdat.speedMtx{i_va} = nan(numel(unique_sf), numel(unique_tf));
    for i_ttype = 1:numel(popdat.on{i_va})
        
        ridx = unique_sf == sf_cpd(i_ttype);
        cidx = unique_tf == tf_cps(i_ttype);
        
        popdat.speedMtx{i_va}(ridx, cidx) = popdat.on{i_va}(i_ttype) - popdat.off{i_va}(i_ttype);
        
    end
    
    % filp the matrix ud so that the SF values go from little to big from
    % the bottom left corner
    popdat.speedMtx{i_va} = flipud(popdat.speedMtx{i_va});
    
    % plot
    subplot(plotdims,plotdims, i_va)
    imagesc(popdat.speedMtx{i_va});% ./ max(abs(popdat.speedMtx{i_va}(:))))
    colormap gray; colorbar
    h = gca;
    h.YTick = [];
    h.XTick = [];
    title(udat.ROI.VisArea{i_va}, 'fontsize', 10)
    axis square
end

