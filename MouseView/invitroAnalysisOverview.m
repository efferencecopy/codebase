function params = invitroAnalysisOverview(params)

% specify useful globals
global GL_DATPATH


%
% LOAD THE EXPERIMENTAL FILES IF PRESENT
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax = loadFromList(params.files);



% define a color direction for each data file
f = figure; map = colormap('jet'); close(f);
clrIdx = round(linspace(1,size(map,1), numel(ax))); % colors for various plots


%
%  IMAGE OF SLICE (if present)
%
%%%%%%%%%%%%%%%%%%%%
if ~isempty(params.photo)% don't plot the figure when the physiology_notes script is auto-running
    
    % plot the photo
    photoPath = findfile(params.photo, [GL_DATPATH, params.mouse], '.jpg');
    img = imread(photoPath);
    figure
    imshow(img);
    set(gcf, 'name', sprintf('%s cell %d', params.mouse, params.cellNum))
    set(gcf, 'position', [582    17   847   598]);
    drawnow
    
    % highlight the stimulation locations
    if isfield(params, 'stimLoc') && size(params.stimLoc, 1) > 0
        centPos = round(ginput(1));
        if ~isempty(centPos) % the user can press return quickly to avoid this part, which will result in an empty vector.
            stimPoints = params.stimLoc;
            pixperum = pixPerMicron(size(img,1), size(img,2));
            stimPoints = round(stimPoints .* pixperum); %now in pix
            stimPoints = bsxfun(@plus, stimPoints, centPos); % pix relative to neuron
            hold on,
            for a = 1:size(stimPoints,1)
                plot(stimPoints(a,1), stimPoints(a,2), 'o', 'markeredgecolor', map(clrIdx(a),:), 'markerfacecolor', map(clrIdx(a),:), 'markersize', 10)
            end
            drawnow
        end
    end
    
end


%
% DC STEPS if present
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(params, 'DCsteps') && ~isempty(params.DCsteps)
    ax_dc = abfobj(params.DCsteps);
    ax_dc.quickPlot
    set(gcf, 'name', sprintf('%s cell %d', params.mouse, params.cellNum))
    set(gcf, 'position', [12     0   560   420])
    drawnow
end


%
%  REMOVE UNWANTED SWEEPS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(params, 'skipSweeps')
    for a = 1:numel(ax)
        if (numel(params.skipSweeps) >= a) && ~isempty(params.skipSweeps{a});
            disp('removing sweeps')
            ax{a} = ax{a}.removeSweeps(params.skipSweeps{a});
        end
    end
end



%
% PLOT THE SERIES RESISTANCE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if numel(params.files)>0
    fhand = figure; hold on,
    set(gcf, 'position', [15    23   560   651]);
    
    idx = 1;
    for a = 1:numel(ax)
        % pull out the series resistance
        Ra = ax{a}.getRa('quick');
        access = permute(Ra.dat, [3,2,1]);
        Verr = permute(Ra.Verr, [3,2,1]);
        if isempty(access)
            continue
        end
        
        nCh = numel(Ra.chNames);
        pltLoc = [1,2;
                  3,4];
        for ch = 1:nCh
            % plot the access values
            hand(pltLoc(ch,1)) = subplot(nCh, 2, pltLoc(ch,1)); hold on
            xx = idx:(idx+size(access,1)-1);
            clr = map(clrIdx(a),:);
            plot(xx(:), access(:,ch), '-ko', 'markerfacecolor', clr, 'markeredgecolor', clr, 'linewidth', 1.5, 'markersize', 5)
            t = title(sprintf('Channel: %s', Ra.chNames{ch}));
            set(t, 'Interpreter', 'none')
            
            
            % plot the vclamp err values
            hand(pltLoc(ch,1)) = subplot(nCh, 2, pltLoc(ch,2)); hold on
            plot(xx(:), Verr(:,ch), '-ko', 'markerfacecolor', clr, 'markeredgecolor', clr, 'linewidth', 1.5, 'markersize', 5)
            t = title(sprintf('Channel: %s', Ra.chNames{ch}));
            set(t, 'Interpreter', 'none')
        end
        
        % update the index
        idx = xx(end)+1;
    end
    
    % tidy up.
    if exist('hand', 'var') && ~isempty(hand)
        for ch = 1:nCh.*2
            subplot(nCh,2,ch)
            axis tight
            ymax = get(gca, 'ylim');
            ylim([0, ymax(2).*1.05])
            xlabel('Sweep Number')
            if any([1,3] == ch)
                ylabel('Series Resistance (MOhms)')
            else
                ylabel('Vclamp Error')
            end
            
        end
    else
        close(fhand) % the figure opened earlier is not useful
    end
end

%
% PERFORM ALL ADDITIONAL ANALYSES
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.ax = ax; % package the raw data so that I don't have to load the abf files multiple times.
for a = 1:numel(params.fxns)
    params = feval(params.fxns{a}, params);
end












