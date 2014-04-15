function invitroAnalysisOverview(params)

% specify useful globals
global GL_DATPATH

% define a few other things
map = colormap('jet'); close;
clrIdx = round(linspace(1,size(map,1), size(params.files,1))); % colors for various plots



%
%  IMAGE OF SLICE (if present)
%
%%%%%%%%%%%%%%%%%%%%
if ~isempty(params.photo)
    
    % plot the photo
    photoPath = findfile(params.photo, [GL_DATPATH, params.mouse], '.jpg');
    img = imread(photoPath);
    figure
    imshow(img);
    
    % highlight the stimulation locations
    if size(params.stimLoc, 1) > 0
        centPos = round(ginput(1));
        if ~isempty(centPos) % the user can press return quickly to avoid this part, which will result in an empty vector.
            stimPoints = params.stimLoc;
            pixperum = pixPerMicron(size(img,1), size(img,2));
            stimPoints = round(stimPoints .* pixperum); %now in pix
            stimPoints = bsxfun(@plus, stimPoints, centPos); % pix relative to neuron
            hold on,
            for a = 1:size(stimPoints,1)
                plot(stimPoints(a,1), stimPoints(a,2), 'o', 'markeredgecolor', map(clrIdx(a),:), 'markerfacecolor', map(clrIdx(a),:))
            end
            drawnow
        end
    end
    
end


%
% DC STEPS if present
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(params.DCsteps)
    ax_dc = abfobj(params.DCsteps);
    ax_dc.quickPlot
    drawnow
end


%
% LOAD THE REMAINING FILES IF PRESENT
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax = {};
for a = 1:numel(params.files)
    disp(params.files{a})
    ax{a} = abfobj(params.files{a});
    
    % remove unwanted sweeps
    if numel(params.skipSweeps) >= numel(params.files)
        disp('removing sweeps')
        ax{a} = ax{a}.removeSweeps(params.skipSweeps{a});
    end        
end


%
% PLOT THE SERIES RESISTANCE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if numel(params.files)>0
    figure, hold on,
    idx = 1;
    for a = 1:numel(ax)
        % pull out the series resistance
        Ra = ax{a}.getRa('quick');
        Ra = permute(Ra, [3,2,1]);
        idx_Im = eval(['ax{a}.idx.',params.validCh,'Im']);
        Ra = Ra(:,idx_Im);
        
        % plot the values
        xx = idx:(idx+numel(Ra)-1);
        clr = map(clrIdx(a),:);
        plot(xx(:), Ra(:), '-ko', 'markerfacecolor', clr, 'markeredgecolor', clr, 'linewidth', 1.5, 'markersize', 5)
        
        % update the index
        idx = xx(end)+1;
    end
    
    % tidy up. 
    axis tight
    ymax = get(gca, 'ylim');
    ylim([0, ymax(2).*1.05])
    xlabel('Sweep Number')
    ylabel('Series Resistance (MOhms)')
end


%
% PERFORM ALL ADDITIONAL ANALYSES
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.ax = ax; % package the raw data so that I don't have to load the abf files multiple times.
for a = 1:numel(params.fxns)
    feval(params.fxns{a}, params)
end












