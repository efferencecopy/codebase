function [bwCell_vol_numbered] = interactiveCellFill(fovImg_vol)

 
% 
% EXAMPLE:       bwCell_vol_numbered = interactiveCellFill(FOVIMG)
% 
% Interactively add/del cells to a mask by clicking on them.
% 
% Saves mask at each step to the UserData field of the figure.If this
% crashes you can pull it out manually:
%    ud = get(gcf,'UserData');
%    recoveredMask = ud.bwMask;
% 
% Hints for using: 
% * Contrast threshold (fract of mean) should be between 0.9
% and 1.05, usually between 0.95 and 1.0.  Higher values mean more
% restrictive/smaller masks, and avoid 'spillover'. I usually use 0.95 for
% dim cells and 1.0 for bright cells. Dilation disk size in pix: 1-3,
% depending on cell sizes.  For cell radii of ~10 pix (512pix,360um FOV), I
% use 3 for 3-5 pix radii, I use 1
% 
% * If two cells are nearby, you can often avoid spillover by clicking the
% bright cell before the dim cell; likewise you can sometimes use bright
% stuff between two cells to segment them, then delete it afterwards.
%
% * If you're not getting a big enough mask for a cell whose center is
% bright, try clicking a bit off center on the dimmer edge.




%
% DEFAULT PARAMETERS
%
brightset=1;


if ischar(fovImg_vol) % if filename
    % read tif file info
    info = imfinfo(fovImg_vol,'tif');
    W=info(1).Width;
    H=info(1).Height;
    nFrames=length(info);
    tmp_im= imread(fovImg_vol,1);
    data_class=class(tmp_im(1));
    
    
    
    %preallocate color stacks
    Red_im=zeros(H,W,nFrames,data_class);
    Green_im=zeros(H,W,nFrames,data_class);
    Blue_im=zeros(H,W,nFrames,data_class);
    
    
    %read frames, separete color chennels
    for frame=1:nFrames
        tmp_im= imread(fovImg_vol,frame);
        Red_im(:,:,frame)=tmp_im(:,:,1);
        Green_im(:,:,frame)=tmp_im(:,:,2);
        Blue_im(:,:,frame)=tmp_im(:,:,3);
    end
end
%fovImg_vol=(Red_im+Green_im+Blue_im)/3;
fovImg_vol=(Green_im);


%zproj:
X1 = 20;
Y1 = 20;
Z1 = 3;
X = X1;
Y = Y1;


%Zproj = squeeze(max(fovImg_vol,[],3));
Zproj = squeeze(fovImg_vol(:,:,Z1));
Xproj = squeeze(fovImg_vol(X1,:,:));
Yproj = squeeze(fovImg_vol(:,Y1,:));

% cah param
Im_contrast = 1; %default is 1, scale min->max of volume

IMG_CONTRAST_INV = 1 - Im_contrast;
MINMAX0 = [min(min(min(fovImg_vol))) max(max(max(fovImg_vol)))];
MINMAX = [(MINMAX0(1) + IMG_CONTRAST_INV*diff(MINMAX0)) ((MINMAX0(2) - IMG_CONTRAST_INV*diff(MINMAX0)))];


%cludgy way to keep scaling on images the same across planes..
fovImg_view = [Xproj' zeros(size(Yproj,2),size(Xproj,2))'; Zproj Yproj];

fovImg_view(1,1) = MINMAX(1);
fovImg_view(1,2) = MINMAX(2);
fovImg_view(find(fovImg_view<MINMAX(1))) = MINMAX(1);
fovImg_view(find(fovImg_view>MINMAX(2))) = MINMAX(2);


fovImg = Zproj;



[nRows0, nCols0, nZ] = size(fovImg_vol);

[nRows, nCols] = size(fovImg_view);
X_skip = nRows - nRows0;


bwCell = false([nRows0,nCols0]);
bwCell_view = false([nRows,nCols]);
bwCell_vol = false([nRows0,nCols0, nZ]);
bwCell_vol_numbered = false([nRows0,nCols0, nZ]);




labCell = bwlabel(bwCell,8);

%%% measure existing cells
stats = regionprops(labCell, 'Area'); % cah this could change
cellAreas = cat(1, stats.Area);

bwvol = bwlabeln(bwCell_vol,6);
nCells = max(max(max(bwvol)));

% cah 2 params:
zoomToAxis = [1 nCols 1 nRows] + 0.5+[-1 1 -1 1];
cellRadius = 1;

%%% auto-computer the cell area dilation disk radius, only if requested
if isempty(cellRadius)
    if nCells < 5
        error('Too few cells to measure cell radius (%d): ', nCells);
    else
        cellRadius = sqrt(mean(cellAreas)/pi);
    end
end
cellRadius = round(cellRadius);

% draw initial figure
figH = figure;

imH = subDrawShaded(fovImg_view, bwCell_view, zoomToAxis, X, Y + X_skip, Z1, nRows0, nCols0, nZ, MINMAX, bwCell_vol);
axH = get(imH, 'Parent');
set(axH, 'Position', [0.06 0.06 0.88 0.88])

%title(sprintf('\\bf%s\\rm: add objects to mask', mfilename));
% set up UI
figPos = get(figH, 'Position');
figColor = get(figH, 'Color');
set(axH, 'Units', 'normalized');
axPos = get(axH, 'Position');
utlH = uicontrol('Style', 'text', ...
    'String', {'Contrast threshold:', ' fract of mean'}, ...
    'Position', [5 figPos(4)-40, 80, 40], ...
    'BackgroundColor', figColor, ...
    'HorizontalAlignment', 'left');
utH = uicontrol('Style', 'edit', ...
    'String', '0.95', ...   % default cThresh
    'Units', 'pixels', ...
    'Position', [5 figPos(4)-70, 60, 30]);

udlH = uicontrol('Style', 'text', ...
    'String', {'Dilation disk:', ' size in pix'}, ...
    'Position', [5 figPos(4)-130, 80, 30], ...
    'BackgroundColor', figColor, ...
    'HorizontalAlignment', 'left');
udH = uicontrol('Style', 'edit', ...
    'String', '1', ...    % default diskR
    'Units', 'pixels', ...
    'Position', [5 figPos(4)-160, 60, 30]);
cdH = uicontrol('Style', 'text', ...
    'String', { 'Cell radius:', ...
    ' (initial disk)', ...
    sprintf(' %5.2fpix', cellRadius) }, ...
    'Units', 'pixels', ...
    'Position', [5 figPos(4)-260, 60, 60], ...
    'BackgroundColor', figColor, ...
    'HorizontalAlignment', 'left');

% title str
tStr = { [mfilename ': Click on a cell region to add'], ...
    '       Shift-click to del, RET to finish, Z to undo' };
tHeight = (1-axPos(4))-axPos(2);
tlH = uicontrol('Style', 'text', ...
    'String', tStr, ...
    'Units', 'normalized', ...
    'Position', [0.25 1-tHeight, 0.5, tHeight*0.8], ...
    'BackgroundColor', figColor, ...
    'HorizontalAlignment', 'left');


%%% iterate: get a point, add it, display it
bwCurr = bwCell;
bwCell_vol(:,:,Z1) = bwCurr;


bright_change=0;

Zproj = squeeze(max(bwCell_vol,[],3));
Xproj = squeeze(bwCell_vol(X1,:,:));
Yproj = squeeze(bwCell_vol(:,Y1,:));

bwCurr_view = [Xproj' zeros(size(Yproj,2),size(Xproj,2))'; Zproj Yproj];


nActions = 0;
nTotal = nCells;
saveTotal = {};

Y_OLD = Y1;
X_OLD = X1;
SHIFTZ = 0; %tag for shifting Z on this iteration
while 1  % till broken out of
    % interactively get clicks
    %   SHIFTZ = 0; %tag for shifting Z on this iteration
    if SHIFTZ == 0
        [X, Y0, selectionType] = getAPoint(gca);
        
        Y = Y0 - X_skip;
        
        if isnan(X)
            key = lower(Y0);
            if isempty(key)
                key = 'NaN';
            end
            if isnumeric(key)
                SHIFTZ = 1;
                shiftz1 = key;
                key='NaN';
            end
            switch key
                case char(13) % return
                    break;  % out of loop, done
                    continue
                case 'q' %quit
                    %prune the values
                    tmp = bwCell_vol_numbered;
                    %first, find all nonzeros values:
                    tmp2 = sort(nonzeros(tmp));
                    tmp2(find(diff([0; tmp2])==0)) = [];
                    %now replace these numbers with ordered numbers from
                    %1:length(tmp2)
                    bwCell_vol_numbered2 = zeros(size(bwCell_vol_numbered));
                    for count5 = 1:length(tmp2)
                        ind = find(bwCell_vol_numbered == tmp2(count5));
                        bwCell_vol_numbered2(ind) = count5;
                    end
                    bwCell_vol_numbered = bwCell_vol_numbered2;
                    
                    break;
                case 'zzz'
                    % undo
                    if nActions == 0
                        printf(1, '** Cannot undo: no cells added yet\n');
                        continue;
                    end
                    tmp = bwSave{nActions};
                    bwCurr = squeeze(tmp(:,:,Z1));
                    %nTotal = saveTotal{nActions};
                    nActions = nActions-1;
                    
                    % draw the new version
                    max_temp=max(max(fovImg_view));
                    fovImg_view_new=fovImg_view.*brightset;
                    fovImg_view_new(find(fovImg_view_new>max_temp)) = max_temp;
                    subDrawShaded(fovImg_view_new, bwCurr, zoomToAxis, X, Y + X_skip, Z1, nRows0, nCols0, nZ, MINMAX, bwCell_vol);
                    fprintf(1, 'Undo!  %d cells total now\n', nTotal);
                    continue
                case 'z'
                    %move up
                    SHIFTZ = 1;
                    shiftz1 = +1;
                case 'x'
                    %move up
                    SHIFTZ = 1;
                    shiftz1 = -1;
                case 'a'
                    brightset=brightset.*1.1;
                    brightset=min([brightset,12]);
                    max_temp=max(max(fovImg_view));
                    fovImg_view_new=fovImg_view.*brightset;
                    fovImg_view_new(find(fovImg_view_new>max_temp)) = max_temp;
                    subDrawShaded(fovImg_view_new, bwCurr_view, zoomToAxis, X, Y + X_skip, Z1, nRows0, nCols0, nZ, MINMAX, bwCell_vol);
                    continue
                case 's'
                    brightset=brightset./1.1;
                    brightset=max([brightset,1]);
                    max_temp=max(max(fovImg_view));
                    fovImg_view_new=fovImg_view.*brightset;
                    fovImg_view_new(find(fovImg_view_new>max_temp)) = max_temp;
                    subDrawShaded(fovImg_view_new, bwCurr_view, zoomToAxis, X, Y + X_skip, Z1, nRows0, nCols0, nZ, MINMAX, bwCell_vol);
                    continue
            end
            continue
        end
    end
    
    
    if bright_change==1
    end
    
    %%% validate XY point
    X = round(X);
    Y = round(Y);
    %   if X <= 0 || X >= nCols || Y <= 0 || Y >= nRows
    if X <= 0 || X >= nCols || Y <= (-1*X_skip) || Y >= nRows0
        % clicked outside axes, repeat
        fprintf(1, 'Click was outside image axes, try again\n');
        continue;
    end
    
    %shift Z axis if click is on XZ or YZ section
    if (((Y > (-1*X_skip)) && (Y <= 0)) || ((X > nCols0) && (X <= nCols))) || SHIFTZ == 1
        if SHIFTZ == 0
            if ((Y > (-1*X_skip)) && (Y <= 0))
                Z1 = nZ + 1 - floor(Y0);
            elseif ((X > nCols0) && (X <= nCols))
                Z1 = nZ + 1 - (X - nCols0);
                Z1 =  (X - nCols0);
            end
        elseif SHIFTZ == 1
            %move up or down
            Z1a = Z1 + shiftz1;
            if Z1a <= nZ && Z1a > 0
                Z1 = Z1a;
            end
            SHIFTZ = 0;
        end
        
        
        
        Y = Y_OLD;
        X = X_OLD;
        
        
        %bwCell_vol(:,:,Z1) = bwCurr;
        %sum(sum(bwCurr))
        %        Zproj = squeeze(max(bwCell_vol,[],3));
        Zproj = squeeze(bwCell_vol(:,:,Z1));
        Xproj = squeeze(bwCell_vol(Y_OLD,:,:));
        Yproj = squeeze(bwCell_vol(:,X_OLD,:));
        bwCurr_view = [flipud(Xproj') zeros(size(Yproj,2),size(Xproj,2))'; Zproj Yproj];
        
        %    Zproj = squeeze(max(fovImg_vol,[],3));
        Zproj = squeeze(fovImg_vol(:,:,Z1));
        Xproj = squeeze(fovImg_vol(Y_OLD,:,:));
        Yproj = squeeze(fovImg_vol(:,X_OLD,:));
        fovImg_view = [flipud(Xproj') zeros(size(Yproj,2),size(Xproj,2))'; Zproj Yproj];
        fovImg = Zproj;
        
        %cludgy way to keep scaling on images the same across planes..
        fovImg_view(1,1) = MINMAX(1);
        fovImg_view(1,2) = MINMAX(2);
        
        fovImg_view(find(fovImg_view<MINMAX(1))) = MINMAX(1);
        fovImg_view(find(fovImg_view>MINMAX(2))) = MINMAX(2);
        
        % draw the new version
        max_temp=max(max(fovImg_view));
        fovImg_view_new=fovImg_view.*brightset;
        fovImg_view_new(find(fovImg_view_new>max_temp)) = max_temp;
        subDrawShaded(fovImg_view_new, bwCurr_view, zoomToAxis, X, Y + X_skip, Z1, nRows0, nCols0, nZ, MINMAX, bwCell_vol);
        fprintf(1, ['New Z position set: ',num2str(Z1),'\n']);
        continue;
    end
    
    
    %%% get ui data
    diskR = str2double(get(udH, 'String'));
    if isnan(diskR)
        fprintf(1, 'Error reading disk radius: %s, try again\n', ...
            get(udH, 'String'));
        continue
    end
    if (diskR < 1) || diskR > 50
        fprintf(1, 'Disk too small or too big, try again\n');
        continue
    elseif ~iswithintol(round(diskR), diskR, 10^2*eps)
        fprintf(1, 'Disk radius must be an integer, try again\n');
        continue
    end
    
    
    cThresh = str2double(get(utH, 'String'));
    if isnan(cThresh)
        fprintf(1, 'Error reading threshold: %s, try again\n', ...
            get(utH, 'String'));
        continue
    end
    if cThresh <= 0 || cThresh >= 100
        fprintf(1, 'Threshold too small or too big, try again\n');
        continue
    end
    
    
    
    %%% what kind of mouse click?
    oldTotal = nTotal;
    switch lower(selectionType)
        case 'extend'    % shift-left: delete
            if bwCurr(Y,X) ~= 1
                
                fprintf(1, '** Trying to delete, not on a cell, try again\n');
                
            else
                
                ind = find(bwCell_vol_numbered == bwCell_vol_numbered(Y,X,Z1));
                
                bwCell_vol_numbered(ind) = 0; % = bwCell_vol_numbered + nActions*bwNew_vol0;
                bwCell_vol = bwCell_vol_numbered>0;
                bwCurr = squeeze(bwCell_vol(:,:,Z1));
                
                
                nTotal = nTotal-1;
                fprintf(1, 'Deleted object, %d total remain\n', nTotal);
            end
            
        case 'normal'    % left-click: add
            % make sure not in a cell
            bwCurr(Y,X)
            if bwCurr(Y,X) == 1
                fprintf(1, '** Trying to add in existing cell region, try again\n');
            else
                bwNew_vol0 = zeros(size(bwCell_vol));
                bwNew_vol = bwCell_vol;
                [bwNew, bwCell, cellMean] = subAddCell(bwCurr, fovImg, X, Y, cellRadius, cThresh, diskR);
                
            end
            
            %step up and down in Z
            for step_Z = [1 -1]
                %            step_Z
                STOP = 0;
                if step_Z == 1
                    Z1_use = Z1;
                elseif step_Z == -1
                    Z1_use = Z1 - 1;
                end
                %Z1_use = Z1_use + step_Z;
                while STOP ~= 1 && (Z1_use > 0) && (Z1_use <= nZ)
                    fovImg_USE = squeeze(fovImg_vol(:,:,Z1_use));
                    bwCurr_USE = squeeze(bwNew_vol(:,:,Z1_use));
                    [bwNew, bwCell] = subAddCell(bwCurr_USE, fovImg_USE, X, Y,cellRadius, cThresh, diskR, cellMean);
                    
                    %check the plane above/below to make sure no touching..
                    Z1_use_tmp = Z1_use + step_Z;
                    
                    %first, check if this plane exists
                    if (Z1_use_tmp > 0) && (Z1_use_tmp <= nZ)
                        bwCell_abovebelow =  squeeze(bwNew_vol(:,:,Z1_use_tmp));
                        %check for direct contact between two cells (diagonal
                        %through Z doesn't get checked, could dilate to do
                        %this)
                        tmp = sum(sum(bwCell_abovebelow.*bwCell));
                        if tmp > 0
                            STOP = 1;
                        end
                    end
                    
                    
                    if sum(sum(bwCell)) > 0 && STOP == 0
                        bwNew_vol(:,:,Z1_use) = bwNew;
                        bwNew_vol0(:,:,Z1_use) = bwCell;
                        Z1_use = Z1_use + step_Z;
                    else
                        
                        STOP = 1;
                    end
                    
                end
            end
            
            %                bwCell_vol_numbered = bwCell_vol_numbered + (nActions+1)*bwNew_vol0;
            bwCell_vol_numbered = bwCell_vol_numbered + (nTotal+1)*bwNew_vol0;
            bwCell_vol = bwCell_vol_numbered>0;
            bwCurr = squeeze(bwCell_vol(:,:,Z1));
            
            
            
            
            nTotal = nTotal+1;
            fprintf(1, 'Added object #%d: %d pix\n', ...
                nTotal, sum(bwNew_vol0(:)));
            
        otherwise
            % other type of click, just continue without comment
            %keyboard
            fprintf(1, 'Unrecognized click occurred (matlab bug?) %s, %g %g\n', ...
                selectionType, X, Y);
            continue
    end
    
    %%% save old mask and update
    bwSave{nActions+1} = bwCell_vol_numbered;
    saveTotal{nActions+1} = oldTotal;
    nActions = nActions+1;
    
    X_OLD = X;
    Y_OLD = Y;
    
    Zproj = squeeze(bwCell_vol(:,:,Z1));
    Xproj = squeeze(bwCell_vol(Y,:,:));
    Yproj = squeeze(bwCell_vol(:,X,:));
    bwCurr_view = [flipud(Xproj') zeros(size(Yproj,2),size(Xproj,2))'; Zproj Yproj];
    
    Zproj = squeeze(fovImg_vol(:,:,Z1));
    Xproj = squeeze(fovImg_vol(Y,:,:));
    Yproj = squeeze(fovImg_vol(:,X,:));
    fovImg_view = [flipud(Xproj') zeros(size(Yproj,2),size(Xproj,2))'; Zproj Yproj];
    fovImg = Zproj;
    fovImg_view(1,1) = MINMAX(1);
    fovImg_view(1,2) = MINMAX(2);
    fovImg_view(find(fovImg_view<MINMAX(1))) = MINMAX(1);
    fovImg_view(find(fovImg_view>MINMAX(2))) = MINMAX(2);
    
    
    % draw the new version
    max_temp=max(max(fovImg_view));
    fovImg_view_new=fovImg_view.*brightset;
    fovImg_view_new(find(fovImg_view_new>max_temp)) = max_temp;
    subDrawShaded(fovImg_view_new, bwCurr_view, zoomToAxis, X, Y + X_skip, Z1, nRows0, nCols0, nZ, MINMAX, bwCell_vol);
end


fprintf(1, '%s: done, %d objects total).\n', mfilename, nTotal);
end % main function


%
%  SUPPORTING FUNCTIONS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function imH = subDrawShaded(fovImg, bw, zoomToAxis, Xuse, Yuse, Zuse, nRows0, nCols0, nZ, MINMAX, bwCell_vol)


    shadeimg = imShade(fovImg, bw);
    cla;
    imH = imagesc(shadeimg,MINMAX);
    set(gca, 'XGrid', 'on', ...
        'YGrid', 'on', ...
        'Visible', 'on', ...
        'YDir', 'reverse', ...
        'DataAspectRatio', [1 1 1]);
    axis(zoomToAxis);

    line(Xuse*ones(1,2),zoomToAxis(3:4),'Color','w');
    line(zoomToAxis(3:4),Yuse*ones(1,2),'Color','w');
    Zuse2 = (nZ + 1 - Zuse);
    Zuse3 = nCols0 + (nZ + 1 - Zuse2);
    line(Zuse3*ones(1,2),zoomToAxis(3:4),'Color','k');
    line(zoomToAxis(1:2),Zuse2*ones(1,2),'Color','k');


    drawnow;


    ud.fovImg = fovImg;
    ud.bwMask = bw;
    ud.bwMask_vol = bwCell_vol;
    set(gcf, 'UserData', ud);
end

function [bwNew, bwCell, cellMean] = subAddCell(bwCurr, fovImg, X, Y, cellRadius, cThresh, diskR, cellMean)

    % Real cell-adding logic is here

    [nRows, nCols] = size(bwCurr);

    %%%% set up dilation disks:  use 4-connected regions for all
    %se = strel('disk',1,8);  % for cells-to-avoid
    se = strel('square',3);  % for cells-to-avoid: all 9 pix in a square;
    % avoid diagonally-connected cells.
    se2 = strel('disk',round(cellRadius),4);  % for region to find mean over
    seJunk = strel('disk', max(round(cellRadius/4), 1), 4);  % remove thin junk
    seExpand = strel('disk', diskR, 4);  % expand thresholded region


    % add a disk around each point, non-overlapping with adj cells
    tempmask = false(nRows, nCols);
    dilateorg = imdilate(bwCurr,se);
    tempmask(Y, X) = 1;
    tempmask = imdilate(tempmask,se2);
    tempmask = tempmask & ~dilateorg;

    % fill region around disk of similar intensity, combine with disk
    if nargin < 8
        cellMean = mean(fovImg(tempmask == 1),1);
    end

    if fovImg(Y, X) > cellMean.*cThresh
        allMeanBw = fovImg >= cellMean.*cThresh;  % threshold by intensity
        %Npts = sum(sum(allMeanBw))
        connMeanBw = bwselect(allMeanBw &~dilateorg, X, Y, 4);
        connMeanBw = connMeanBw |tempmask & ~dilateorg;

        % erode then dilate filled to remove sharp things
        erMean = imerode(connMeanBw, seJunk);
        dilateMean = imdilate(erMean, seJunk);
        dilateMean = imdilate(dilateMean, seExpand); % (thresh is conservative)
        bwCell = dilateMean & ~dilateorg;

        bwNew = bwCurr | bwCell;
    else
        bwCell = zeros(size(bwCurr));
        bwNew = bwCurr;
    end
end


