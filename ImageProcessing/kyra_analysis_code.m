% LOAD 'udat'
[filename, path] = uigetfile({'*.mat',...
                               'Related Files (*.mat)'},...
                               'Select the Related MATLAB file',...
                               'S:\Data');
loadpath = [path, filesep, filename];
load (loadpath)

% Display choice
    if isequal(filename, 0)
        errordlg(sprintf(' File not found.\n User selected "Cancel", a file MUST\n be selected for code to progress'), 'File Error')
    end


% CREATE MATRIX (SELECTED-PIXELS ACROSS TIME) 
% Preallocate space
TotalNVAs = numel(udat.ROI.VisArea);
for i_TotalNVAs = 1: TotalNVAs
    pixelMatrix.AvgImage{i_TotalNVAs} = {};
end

%%%%% AvgImage
for i_VAs = 1: TotalNVAs % for every Visual Area 'available' in AvgImage
    
    if find(isempty (udat.ROI.AvgImage{i_VAs}))
        udat.ROI.AvgImage{i_VAs} = {};
        pixelMatrix.AvgImage{i_VAs} = udat.ROI.AvgImage{i_VAs};
    else
        Npixels = size(udat.ROI.AvgImage{i_VAs}, 1); % Determine # of pixels in one selected area
            columns = udat.ROI.AvgImage{i_VAs}(:,1)';
            rows = udat.ROI.AvgImage{i_VAs}(:,2)';
        
        linIdx = sub2ind(size(udat.final_img{1}), rows, columns); %Placement of 'important' pixels for the first frame
            % Note: size(udat.final_img{1}) == size(udat.Avg.TrialON{1}(:,:,1))
        Nframes = size(udat.preProcessed_ON{1}, 3) + size(udat.preProcessed_OFF{1}, 3);
        linIdx = repmat(linIdx, Nframes, 1); % lists "placement" Nframes times (not taking into account change in frame); rows = # frames & columns = Npixels
        
        Nelements = prod(size(udat.final_img{1})); % Determines the total number of elements in the matrix
        adjuster = repmat(Nelements, 1, Npixels);  % length of vector = Npixels
        adjuster = bsxfun(@times, repmat(adjuster, Nframes, 1), [0:1:Nframes-1]');
        
        newIdx = linIdx + adjuster; % applies adjuster to the "placement"
        newIdx = newIdx(:); % lists every idx# (in a column) SANITY CHECK: rows of newIdx should = Npixels*Nframes
    
        Nttypes = numel(udat.preProcessed_ON); %Nttypes = Total # of Stimulus Types (Trial Types)
        for i_ttypes = 1 : Nttypes
            AvgTrial{i_ttypes} = cat(3, udat.preProcessed_ON{i_ttypes}, udat.preProcessed_OFF{i_ttypes});
            
            ans = AvgTrial{i_ttypes}(newIdx); % lists the pixel value for each position specified by "placement" list
            pixelMatrix.AvgImage{i_VAs}{i_ttypes} = reshape(ans, Nframes, []); % reshapes the list (columnar) into a matrix pixels across time (rows = Nframes, columns = pixels)
        end
    end
end

% PLOT: PIXELS VS. TIME
TotalNVAs = numel(udat.ROI.VisArea);
Nttypes = numel(udat.preProcessed_ON);
NframesON = size(udat.preProcessed_ON{1}, 3);

for i_TotalNVAs = 1:TotalNVAs
    
    if find(~isempty(pixelMatrix.AvgImage{i_TotalNVAs})); % For the cell arrays that contain something (selected ROIs)
        
        fig = figure;
        fig.Position = [203  199  1451  365];
        
        figName = regexpi(filename, '(\w\d+)', 'match');
        figName = regexpi(figName{1}, '(\d)', 'match');
        figName = ['k', figName{end-1}, figName{end}];
            set(gcf, 'name', sprintf('%s (%s)', figName, udat.ROI.VisArea{i_TotalNVAs}), 'numbertitle', 'off')
        
        for i_ttypes = 1: Nttypes
            % Plot the ROI pixels according to stimulus types
            subplot(1, Nttypes, i_ttypes)
            
            tt = [0:Nframes-1] .* (1/udat.frameRate);
            plot(tt, pixelMatrix.AvgImage{i_TotalNVAs}{i_ttypes})
            
            title(sprintf('\n%s\n(SF: %0.2f,  TF: %d)\n', udat.ROI.VisArea{i_TotalNVAs}, udat.ttypes(i_ttypes), udat.ttypes(i_ttypes+Nttypes)), 'FontSize', 10);
            xlabel(sprintf('Time (sec)', udat.frameRate), 'FontSize', 8)
            ylabel('dF/F', 'FontSize', 8)
            axis([0  tt(end)  -0.05  0.15])
            
            ttON = (NframesON-1).* (1/udat.frameRate);
            line([ttON  ttON], ylim, 'LineStyle', ':', 'Color', 'k') % Mark 'Stimulus Offset'
        end
    end
end


%%% NOTATION %%%
% 'pixelMatrix' includes all the pixel matrices (where rows = Nframes and columns = Npixels of ROI selected)
% 'pixelMatrix.AvgImage' refers to pixel matrices specific to the ROIs selected in the AvgImage
% 'pixelMatrix.AvgImage{?}' (outer cell) refers to the specific ROI selected (i.e. 'A', 'AL'...)
% 'pixelMatrix.AvgImage{}{?}' (inner cell) refers to the stimulus applies to the AvgTrialON section