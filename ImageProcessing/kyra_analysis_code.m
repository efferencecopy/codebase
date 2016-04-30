%% LOAD THE DATA
fin

% Define some parameters for the analysis
dataSource = 'preProcessed_';  % Can be: 'AvgTrial' or 'preProcessed_'


dataSourceOn = [dataSource, 'ON'];
dataSourceOff = [dataSource, 'OFF'];

% Open dialogue box for selecting the processed ROI file
[filename, path] = uigetfile({'*.mat', 'Related Files (*.mat)'},...
                              'Select the Processed ROI file',...
                              'MultiSelect', 'on',...
                              'S:\Data');
loadpath = [path, filesep, filename];
    
clc % Clear the command window so that the user can see the feedback. 

data = {};
% If no file is selected,
if isequal(filename, 0)
    % Present an error message (Code cannot proceed without data)
    errordlg(sprintf(' File not found.\n User selected "Cancel", a file MUST\n be selected for code to progress'),...
             'File Error')                

% If file(s) are selected,
% Determine how many files & deal with accordingly
elseif size(loadpath, 1) == 1 && numel(loadpath) >= 80;
    data{1} = load (loadpath);
    fprintf ('File Selected:\n%s\n\n', filename)

else numel(loadpath) >= 1;
    fprintf('Files Selected:\n')
    % For more than one file,
    % Save the data of each file into a different 'data' cell array
    for i_files = 3: numel(loadpath)
        data{i_files-2} = load ([loadpath{1}, loadpath{i_files}]);
      
        fprintf('%s \n', filename{i_files-2})
        
    end
    fprintf('\n')
end

% Force the analysis to use a "standardized" set of ROIs
% udat.ROI = standardROIs
    
% Determine the # of Imaging Days (according to to the number of files selected)
Ndays = numel(data); 
% Initialize the 'pixelMatrix' variable
for i_day = 1:Ndays
    % If there is more than one data file,
    % Save the Date (when Imaging Data was collected)
    if numel(data) > 1
        Date = {};
        Date{i_day} = regexpi(filename{i_day}, '(\d\d\d\d\d\d)', 'match');
        if Date{i_day}{:}(3) == '0'
            Date{i_day} = [Date{i_day}{:}(4), '/', Date{i_day}{:}(5:6), '/', Date{i_day}{:}(1:2)]; % Format the date
        else
            Date{i_day} = [Date{i_day}{:}(3:4), '/', Date{i_day}{:}(5:6), '/', Date{i_day}{:}(1:2)]; % Format the date
        end
        
        data{i_day}.udat.imgDate = Date{i_day};
    end

    NVAs = numel(data{i_day}.udat.ROI.VisArea);
    for i_VAs = 1: NVAs
        pixelMatrix{i_day}{i_VAs} = {};
    end
end

% Set a 'help' variable: 
% Command Window 'notifications' do not repeat unecessarily
Option = true;
% Error messages only show up once
err = true;

for i_day = 1: Ndays
    
    % Determine the total # of Visual Areas
    NVAs = numel(data{i_day}.udat.ROI.VisArea);
    for i_VAs = 1: NVAs
        
        if any(isnan(data{i_day}.udat.ROI.RCidx{i_VAs}))
            % Don't bother analyzing ROIs that haven't been defined by User
            continue 
        end
        
        % Determine the # of pixels in a selected ROI
        Npixels = size(data{i_day}.udat.ROI.RCidx{i_VAs}, 1); 
        % Determine the 'row' and 'column' coordinates of each pixel (within a ROI)
        columns = data{i_day}.udat.ROI.RCidx{i_VAs}(:,2)';
        rows = data{i_day}.udat.ROI.RCidx{i_VAs}(:,1)';
        
        %%% Convert row/column index to a linear index:
        % Determine the linear index for pixels within the first frame ("placement")
        linIdx = sub2ind(size(data{i_day}.udat.final_img{1}), rows, columns); 
        
        % Take into account that Nframes may not be equal for ON and OFF periods.
        NframesON{i_day} = size(data{i_day}.udat.(dataSourceOn){1}, 3);
        NframesOFF{i_day} = size(data{i_day}.udat.(dataSourceOff){1}, 3);
        
        % Calculate SF & TF for Tuning function
        SF{i_day} = data{i_day}.udat.ttypes(:,strcmpi(data{i_day}.udat.text, 'tGratingSpatialFreqCPD'));
        TF{i_day} = data{i_day}.udat.ttypes(:,strcmpi(data{i_day}.udat.text, 'tGratingTemporalFreqCPS'));

        % Number of Frames must be equal across selected files...
        if NframesON{1} ~= NframesON{i_day} || NframesOFF{1} ~= NframesOFF{i_day}
                error(sprintf(' Number of frames do not match-up!\n Code cannot proceed.'),...
                         'File Mismatch')
        end    
        % Stimulus Types (Number of Stim) must be equal across selected files...
        if numel(SF{1}) ~= numel(SF{i_day}) || numel(TF{1}) ~= numel(TF{i_day});
                errordlg(sprintf('ERROR: Stimulus Types do not match-up! Code cannot proceed.'),...
                         'File Mismatch')
        end
        
        % Lists "placement" Nframes times (Does not tkae into account change in frame)
        % [rows, columns] = [# of frames, # of pixels]
        linIdxON = repmat(linIdx, NframesON{1}, 1); 
        linIdxOFF = repmat(linIdx, NframesOFF{1}, 1);
        
        % Determines the total number of elements in the matrix
        Nelements = numel(data{i_day}.udat.final_img{1}); 
        
        adjuster = repmat(Nelements, 1, Npixels);  % Length of vector = Npixels
        adjusterON = bsxfun(@times, repmat(adjuster, NframesON{1}, 1), [0 : 1 : NframesON{1}-1]');
        adjusterOFF = bsxfun(@times, repmat(adjuster, NframesOFF{1}, 1), [0 : 1 : NframesOFF{1}-1]');
        
        % NewIdx applies the adjuster to the "placement" (linear index)
        newIdxON = linIdxON + adjusterON; 
        newIdxOFF = linIdxOFF +adjusterOFF;
        % Lists every Index # (in a column) 
        % Sanity Check: # of rows should equal Npixels*Nframes
        newIdxON = newIdxON(:); 
        newIdxOFF = newIdxOFF(:);
        %%%
        
        %%%%% Use the linear index to extract the Dfof data from the processed data
        % Nttypes = Total # of Stimulus Types (Trial Types)
        Nttypes = size(data{i_day}.udat.ttypes, 1); 
        for i_ttypes = 1 : Nttypes
            
            % Dfof data for the ON period. [NframesON x NPixels]
            ROI.on = data{i_day}.udat.(dataSourceOn){i_ttypes}(newIdxON); % Lists the pixel value for each position specified by "placement" list
            ROI.on = reshape(ROI.on, NframesON{1}, []); % Reshapes the list (columnar) into a matrix of pixels across time [NframesON, Npixels]
            
            % Dfof data for the OFF period. [NframesOFF x NPixels]
            ROI.off = data{i_day}.udat.(dataSourceOff){i_ttypes}(newIdxOFF);
            ROI.off = reshape(ROI.off, NframesOFF{1}, []);
            
            % Concatenate the ON and OFF periods.
            pixelMatrix{i_day}{i_VAs}{i_ttypes} = cat(1, ROI.on, ROI.off);
        end
        %%%%%      
    end
    
    % Only display the different days if there are multiple days to display
    if Ndays > 1
        
        % HELPON: Only display once per session (as noted above)
        if Option;
            fprintf('Loading from:\n(%d days total)\n', Ndays)
        end
        Option = false;
        
        fprintf('\n  Day %s\n', data{i_day}.udat.imgDate)
        % Display the # of frames (Sanity Check)
        fprintf('   Number of ON Frames: %d  (%d sec)\n   Number of OFF Frames: %d (%d sec)\n',...
            NframesON{i_day}, NframesON{i_day}.* (1/data{i_day}.udat.frameRate),...
            NframesOFF{i_day}, NframesOFF{i_day}.*(1/data{i_day}.udat.frameRate))
    end    
end


fprintf('\nData has been loaded.\n\n')


%% SANITY CHECK: DISPLAY SELECTED-ROIs

figName = regexpi(filename{1}, '(\w\d+)', 'match');
figName = regexpi(figName{1}, '(\d)', 'match');
figName = ['k', figName{end-1}, figName{end}];

figure; hold on,
set(gcf, 'Name', sprintf('%s Selected-ROIs', figName), 'NumberTitle', 'off',...
    'Units', 'normalized', 'Position', [0.0117  0.1388  0.9797  0.7087]);

for i_day = 1: Ndays
    
    plotimg = cat(3, data{i_day}.udat.final_img{:});
    plotimg = mean(plotimg, 3);
    plotimg = plotimg + (abs(min(plotimg(:))));
    plotimg = plotimg ./ (max(plotimg(:)).*1.2);
    plotimg = repmat(plotimg, [1,1,3]);
    
    for i_roi = 1:numel(data{i_day}.udat.ROI.RCidx)
        
        rc = data{i_day}.udat.ROI.RCidx{i_roi};
        if ~isempty(rc) && ~any(isnan(rc(:)))
            linind = sub2ind(size(plotimg(:,:,1)), rc(:,1), rc(:,2));
            plotimg(linind) = 0;
        end
    end
    
    subplot(ceil(Ndays/5), 5, i_day)
    h_img = imshow(plotimg);
    
    title(data{i_day}.udat.imgDate, 'FontSize', 10);
end

%% PLOT: PIXELS VS. TIME

% Initialize variables
for i_day = 1:Ndays
    meanON{i_day} = [];
    meanOFF{i_day} = [];
end

if Ndays > 1
    Option = true;    
end

NVAs = numel(data{i_day}.udat.ROI.VisArea);
for i_VAs = 1:NVAs
    
    % Include a legend for EVERY figure presented,
    % But display legend only ONCE per figure
    repeat = true;
    
    if isempty(pixelMatrix{i_day}{i_VAs});
        % Don't bother plotting ROIs that haven't been defined
        continue
    end
    
    figure;
    if Ndays > 1
        figName = regexpi(filename{1}, '(\w\d+)', 'match');
    else
        figName = regexpi(filename, '(\w\d+)', 'match');
    end
    figName = regexpi(figName{1}, '(\d)', 'match');
    figName = ['k', figName{end-1}, figName{end}];
    set(gcf, 'Name', sprintf('%s (%s)', figName, data{i_day}.udat.ROI.VisArea{i_VAs}), 'numbertitle', 'off');
    
    Nttypes = numel(data{i_day}.udat.(dataSourceOn));
    for i_ttypes = 1: Nttypes
        
        % Change figure formatting if there are a lot of Stimulus Types to display
        if Nttypes > 5
            set(gcf, 'Units', 'normalized');
            set(gcf, 'Position', [0.2969 -0.2375 1.0688 1.0913]);
            subplot (round(sqrt(Nttypes)), round(sqrt(Nttypes)), i_ttypes)
        elseif Nttypes > 3
            set(gcf, 'Position', [95 139 1683 334]);
            subplot(1, Nttypes, i_ttypes)
        else
            set(gcf, 'Position', [203  199  1451  365]);
            subplot(1, Nttypes, i_ttypes)
        end
        
        % Plot the ROI pixels according to stimulus types
        Nframes = size(pixelMatrix{i_day}{i_VAs}{i_ttypes}, 1);
        tt = [0:Nframes-1] .* (1/data{i_day}.udat.frameRate);
        
        if Ndays > 1
            % Do not offer the option of plotting the Raw Data
            hold on;
            for i_day = 1: Ndays
                % **NOTE: Code breaks if Nttypes is not the same across days
                plot(tt, mean(pixelMatrix{i_day}{i_VAs}{i_ttypes}, 2), 'LineWidth', 2); % Plot the mean dF/F over time on the same axis
                axis([0  tt(end)  -0.0035  0.0035])
            end
            hold off;
            hold on;
            
            % Include legend if there is more than one day of data
            if repeat
                Legend = cell(length(pixelMatrix),1);
                for numdays = 1: length(pixelMatrix)
                    Legend{numdays} = strcat(data{numdays}.udat.imgDate);
                end
                legend(Legend, 'Location', 'none', 'FontSize', 8,...
                       'Units', 'normalized', 'Position', [0.02  0.91  0.07  0.07])
                repeat = false;
            end
            
        else
            if Option
                choice = questdlg('Plot the Raw Data?',...
                    '::Question','Yes','No','Yes');
            end
            Option = false;
            
            switch choice
                case 'Yes'
                    % Plot ALL dF/F values over time
                    plot(tt, pixelMatrix{i_day}{i_VAs}{i_ttypes}); hold on 
                    % Plot the mean dF/F over time
                    plot(tt, mean(pixelMatrix{i_day}{i_VAs}{i_ttypes}, 2), 'k--', 'LineWidth', 1.5); 
                    hold on
                    axis([0  tt(end)  -0.0065  0.0065])
                case 'No'
                    % Plot the mean dF/F over time
                    plot(tt, mean(pixelMatrix{i_day}{i_VAs}{i_ttypes}, 2), 'b-', 'LineWidth', 2.5); 
                    hold on 
                    axis([0  tt(end)  -0.0035  0.0035])
            end
        end
        
        % Change figure formatting if there are alot of stim-types to display
        if Nttypes > 5
            title(sprintf('\n(SF: %0.2g,  TF: %0.2g)\n', data{i_day}.udat.ttypes(i_ttypes), data{i_day}.udat.ttypes(i_ttypes+Nttypes)), 'FontSize', 10);
        else
            title(sprintf('\n%s\n(SF: %0.2f,  TF: %d)\n', data{i_day}.udat.ROI.VisArea{i_VAs}, data{i_day}.udat.ttypes(i_ttypes), data{i_day}.udat.ttypes(i_ttypes+Nttypes)), 'FontSize', 10);
        end
        
        xlabel(sprintf('Time (sec)', data{i_day}.udat.frameRate), 'FontSize', 8)
        ylabel('\DeltaF/F', 'FontSize', 8)
        
        ttON = (NframesON{i_day}).* (1/data{i_day}.udat.frameRate);
        line([ttON  ttON], ylim, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1.1) % Mark 'Stimulus Offset'
        
        %%%%% Find mean Dfof across time
        for i_day = 1: Ndays
            % Select only the ON or OFF frames of ROI
            day_meanON{i_day} = pixelMatrix{i_day}{i_VAs}{i_ttypes}(1:NframesON{i_day},:);
            day_meanOFF{i_day} = pixelMatrix{i_day}{i_VAs}{i_ttypes}(NframesON{i_day}+1:Nframes,:);
            
            day_meanON{i_day} = mean(mean(day_meanON{i_day}, 2));
            day_meanOFF{i_day} = mean(mean(day_meanOFF{i_day}, 2));
            
            % Save the mean dF\f  for future analysis/graphs
            ROI.day_meanON{i_day}{i_VAs}{i_ttypes} = day_meanON{i_day};
            ROI.day_meanOFF{i_day}{i_VAs}{i_ttypes} = day_meanOFF{i_day};
        end
        
        meanON = {mean(cat(3, day_meanON{:}), 3)}; % Average across days
        plot([0, ttON], [meanON{:}, meanON{:}], 'k-', 'LineWidth', 1.2); hold on
        
        meanOFF = {mean(cat(3, day_meanOFF{:}), 3)};
        plot([ttON, tt(end)], [meanOFF{:}, meanOFF{:}], 'k-', 'LineWidth', 1.2); hold on
        %%%%%
        
        % Save the mean(across days) dF\f  for future analysis/graphs
        ROI.meanON{i_VAs}{i_ttypes} = meanON;
        ROI.meanOFF{i_VAs}{i_ttypes} = meanOFF;
    end
    fprintf('  Plotting dF/F over time for %s\n', data{i_day}.udat.ROI.VisArea{i_VAs})
end

fprintf('\n')


%% SF-TF TUNING FUNCTION

SF = cat(1, SF{:});
TF = cat(1, TF{:});

uniqueSF = unique(SF);
uniqueTF = unique(TF);

fig = figure;
fig.Units = 'normalized';
fig.Position = [0.1789 0.0688 0.7000 0.7925];
set(gcf, 'name', sprintf('%s Speed Tuning', figName), 'numbertitle', 'off')

ROI.Tuning = {};
NVAs = numel(data{i_day}.udat.ROI.VisArea);
for i_VAs = 1:NVAs
    if isempty(ROI.meanON{i_VAs})
        % Leave ROIs that haven't been defined empty
        ROI.Tuning{i_VAs} = [];
    continue
    end
   
        ROI.Tuning{i_VAs} = nan(numel(uniqueSF), numel(uniqueTF));
        
        % Number of Stimulus types per Visual Area
        NStim = numel(ROI.meanON{i_VAs});
        for i_Stim = 1: NStim
            
            Ridx = uniqueSF == SF(i_Stim);
            Cidx = uniqueTF == TF(i_Stim);
            
            ROI.Tuning{i_VAs}(Ridx, Cidx) = ROI.meanON{i_VAs}{i_Stim}{:} - ROI.meanOFF{i_VAs}{i_Stim}{:};
        end
        
        % Plot
        subplot (round(sqrt(NVAs)), round(sqrt(NVAs)), i_VAs)
       
        imagesc(flipud(ROI.Tuning{i_VAs}))
        colormap copper;
        c = colorbar; c.Label.String = 'Avg. Response (\DeltaF/F)';
        
        title(data{i_day}.udat.ROI.VisArea{i_VAs}, 'FontSize', 11);
        xlabel('TF (Hz)', 'FontSize', 9)
        ylabel('SF (cpd)', 'FontSize', 9)
        h = gca;
        h.XTick = [];
        h.YTick =[];
        set(h,'YDir','normal') % (with flipud) makes SF go from smallest to largest; bottom to top
        
        axis square
        fprintf('Plotting dF/F Response Tuning for %s\n', data{i_day}.udat.ROI.VisArea{i_VAs})
   
end

%% Plot meanON vs Speed

Speed = TF ./ SF; %(degrees of visual angle per sec)
uniqueSpeed = unique(Speed);

% Initialize 'ROI.Speed' & 'ROI.meanSpeed'
ROI.Speed = {};
ROI.meanSpeed = {};
for i_day = 1:Ndays
    ROI.Speed{i_day} = repmat({repmat({[]}, 1, numel(uniqueSpeed))}, 1, NVAs); 
    % Ordered according to uniqueSpeed
end

fig = figure;
fig.Units = 'normalized';
fig.Position = [0.1789 0.0688 0.7000 0.7925];
colormap(subplot (round(sqrt(NVAs)), round(sqrt(NVAs)), i_VAs), copper);

set(gcf, 'name', sprintf('%s Response Strength vs. Speed', figName), 'numbertitle', 'off')

repeat = true;
for i_VAs = 1:NVAs
    
    if isempty(ROI.meanON{i_VAs})
        % Don't bother plotting ROIs that haven't been defined
        continue
    end
   
    % Determine ROI Speed/meanSpeed
    for i_day = 1: Ndays
        NStim = numel(ROI.meanON{i_VAs}); % Number of Stimulus types per Visual Area
        for i_Stim = 1: NStim
            SpeedIdx = uniqueSpeed == Speed(i_Stim);
            
            ROI.Speed{i_day}{i_VAs}{SpeedIdx}(end+1) = ROI.day_meanON{i_day}{i_VAs}{i_Stim} - ROI.day_meanOFF{i_day}{i_VAs}{i_Stim};    
            ROI.meanSpeed{i_VAs}{SpeedIdx} = ROI.meanON{i_VAs}{i_Stim};
        end
    end
    
    subplot(round(sqrt(NVAs)), round(sqrt(NVAs)), i_VAs);
    bfline = {};
        hold on;
        % Plot the 'best-fit' line
        for i_day = 1: Ndays
            bfline{i_day} = plot(uniqueSpeed, cellfun(@mean, ROI.Speed{i_day}{i_VAs}), 'LineWidth', 1.5); 
        end
       
        hold on; 
        % If there is more than one day of data
        if Ndays >1
            % Include a legend (if there is more than one day of data)
            if repeat
                Legend = cell(length(ROI.Speed),1);
                for numdays = 1: length(ROI.Speed)
                    Legend{numdays} = strcat(data{numdays}.udat.imgDate);
                end
                legend(Legend, 'Location', 'none', 'FontSize', 8,...
                       'Units', 'normalized', 'Position', [0.02  0.90  0.07  0.07])
                repeat = false;
            end
            
            % Plot the mean across days (if there is more than one day of data)
            for i_Speed = 1: numel(ROI.Speed{i_day}{i_VAs})
                AvgSpeed(i_VAs, i_Speed) = ROI.meanSpeed{i_VAs}{i_Speed}{:};
            end
                plot(uniqueSpeed, AvgSpeed(i_VAs, :)', 'k-', 'LineWidth', 2)

            % Make each day's data a dotted line (if there is more than one day of data)
            for i_day = 1:Ndays
                set(bfline{i_day}, 'LineStyle', '-.', 'LineWidth', 1.3);
            end
        end
        
        % Plot the Raw Data Points if there are not that many Stimulus Types
        % NOTE: The # of data points is equal to the # of stim-type
        if Nttypes <= 7      
            hold on;
            for i_day = 1:Ndays
                for i_Speed = 1: numel(ROI.Speed{i_day}{i_VAs})
                    raw_data = plot(uniqueSpeed(i_Speed), ROI.Speed{i_day}{i_VAs}{i_Speed}, '.', 'MarkerSize', 17.5); 
                    
                    % Match data-point colors to corresponding line
                    raw_data.Color = bfline{i_day}.Color;
                end
            end
            hold off;
        end
    
    % Formatting Details
    title(sprintf('\n%s', data{i_day}.udat.ROI.VisArea{i_VAs}), 'FontSize', 11);
    xlabel(sprintf('Speed (%c/sec)', char(176)), 'FontSize', 9)
    ylabel('\DeltaF/F ON', 'FontSize', 9)
    set(gca,'XScale', 'log', 'Ylim', [-.004  .003])
end