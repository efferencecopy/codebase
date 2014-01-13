% Figures for a NeuroThresh paper
%
% Section 1) A couple of views of an example data set (termination
% points for a single cell and optional staircases in gray).
%
% Section 1.1) As above, but only two views per cell.
%
% Section 1.2) As above, but only one view per cell + 3-D isosurfaces.
%
% Section 2) Distribution of baseline firing rates and thresholds.
% 
% Section 2.1) Distribution of baseline firing rates and thresholds
% (version 2)
%
% Section 3) Distribution of p-values from bootstrap test. (Or conventional
% F-test)
%
% Section 4) Distribution of cone weights for all cells in a Lennie-style
% plot with planar cells highlighted.
%
% Section 5) Example Gabor.
%
% Section 6) Single cell examples, predicting grating responses from
% NeuroThresh data.
%
% Section 7) Scatterplot of correlation coefficients between responses to
% colored gratings and planar/quadratic predictions.
%
% Section 8) A pretty plot of an example data set (grating responses +
% isoresponse surfaces)
%
% Section 9) Breakdown of the different types of cells and their principal
% axes
%
% Section 9.1) Principal axes of quadratic fits + randomization control + 
% omitting the first round data points
%
% Section 10) Effects of changing the threshold.  Version 1. Transition
% matrix
%
% Section 10.1) Effects of changing the threshold.  Version 2. F-tests.
%
% Section 10.2) Effects of changing threshold. Version 3. Superimposed data
% points with the optimal scaling parameter.
%
% Section 11) Modulation ratios
%
% Section 12) Finding intertrial intervals and RF eccentricites. 
% (not really a figure)
%
% Section 13) Patrick Weller's model simulating the experiment on a cell
% with a linear contrast-response function and an ellipsoidal isoresponse
% surface.
%
% Section 14) Hierarchical model to create nonplanar isoresponse surfaces.
%
% Section 15) Is basing triangles on 2 deg cone-contrast the same as basing
% triangles on (a) 10 deg cone contrast or (b) RGB intensities? 

%% Section 0: Execute this first
%examplefilenames = {'S041310002','K070309004','K082609009'};
examplefilenames = {'K042811002','K070309004','K082609009'};

%% Section 1) Single cell examples.  Termination points only.
LARGEAXISWIDTH = 1.3;
SMALLAXISWIDTH = 1;
PLOTPLANES = 1;
PLOTQUAD = 1;

figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
colordef white;

%filename = 'K082509007.nex';  % planar L-M+S cell
%filename = 'S041310002.nex';  % planar luminance cell
%filename = 'K072009003.nex';  % Luminance funnel cell
%filename = 'K070309004.nex';  % Funnel intermediate color
%filename = 'K082609010.nex';  % pan color
for col = 1:length(examplefilenames)
    filename = examplefilenames{col};
    % Getting best viewangles
    if (strcmp(filename, 'S041310002'))
        customview = [161 22];
        xlim = [-.16 .16]; ylim = [-.16 .16]; zlim = [-1 1];
    elseif (strcmp(filename, 'K070309004'))
        customview = [85 3];
        xlim = [-.07 .07]; ylim = [-.09 .09]; zlim = [-1 1];
    elseif (strcmp(filename, 'K082609010'))
        customview = [-75 9];
        xlim = [-.09 .09]; ylim = [-.09 .09]; zlim = [-.5 .5];
        xlim = [-.07 .07]; ylim = [-.09 .09]; zlim = [-1 1];
    elseif (strcmp(filename, 'K042811002'))
        customview = [-75 9];
        xlim = [-.09 .09]; ylim = [-.09 .09]; zlim = [-.5 .5];
    else
        customview = [-102 -13];
    end
    views = [0 0; 90 0; 0 90; customview];
    
    stro = nex2stro(findfile(filename));
    lmsidxs = [find(strcmp(stro.sum.trialFields(1,:),'lcont'))...
        find(strcmp(stro.sum.trialFields(1,:),'mcont'))...
        find(strcmp(stro.sum.trialFields(1,:),'scont'))];
    ntrials = size(stro.trial(:,lmsidxs),1);
    % -------------------------------
    % Converting cone contrasts in nex file to 10 deg fundmentals.
    % Must go through excitations first - can't just transform contrasts.
    % -------------------------------
    fundamentals = stro.sum.exptParams.fundamentals;
    fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
    mon_spd = stro.sum.exptParams.mon_spd;
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
    M = fundamentals'*mon_spd;
    load ('T_cones_smj10');
    M10 = T_cones_smj10*mon_spd;
    stro.trial(:,lmsidxs) = ConvertConeContrastBasis(M, M10, stro.sum.exptParams.bkgndrgb, stro.trial(:,lmsidxs));
   % -------------------------------
    out = NTpreprocess(stro,0,Inf);  % Getting the termination points
    scaled = out(:,[2 3 4]).*repmat(out(:,5),1,3);
    % Moving to a DKL-ish space
    scaled(:,[1 2]) = scaled(:,[1 2])*mkbasis([1 1; 1 -1]);
    Loog = logical(out(:,end));
    if (PLOTPLANES || PLOTQUAD)
        [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);
    end
    spikeidx = strcmp(stro.sum.rasterCells(1,:),'sig001a');
    % ---------------------------
    % Plotting thresholds
    % ---------------------------
    for i = 1:4
        if (i == 4)
            axes('Position',[.5+(col-1)*LARGEAXISWIDTH*1.7 5.5 LARGEAXISWIDTH LARGEAXISWIDTH]);
            axis vis3d;
            markersize = 3;
        else
            axes('Position',[.7+(col-1)*LARGEAXISWIDTH*1.7 (i-.5)*SMALLAXISWIDTH*1.5 SMALLAXISWIDTH SMALLAXISWIDTH]);
            markersize = 2;
        end
        hold on;
        axis square; set(gca,'View',views(i,:));
        xlabel('L+M'); set(get(gca,'XLabel'),'Color',[0 0 0]);
        ylabel('L-M'); set(get(gca,'YLabel'),'Color',[0 0 0]);
        zlabel('S'); set(get(gca,'ZLabel'),'Color',[0 0 0]);
        h=plot3(scaled(~Loog,1),scaled(~Loog,2),scaled(~Loog,3),'ko');
        set(h,'Markersize',markersize,'Markerfacecolor','black');
        h=plot3(-scaled(~Loog,1),-scaled(~Loog,2),-scaled(~Loog,3),'ko');
        set(h,'Markersize',markersize,'Markerfacecolor','black');
        
        % Plotting OOG rays (squeezing them inside the axis limits)
        plotlims = [xlim(2) ylim(2) zlim(2)];
        normfacts = max(abs(scaled(Loog,:)./repmat(plotlims,sum(Loog),1)),[],2);
        normfacts = max(normfacts,1);
        h=plot3([zeros(sum(Loog),1) scaled(Loog,1)./normfacts]',[zeros(sum(Loog),1) scaled(Loog,2)./normfacts]',[zeros(sum(Loog),1) scaled(Loog,3)./normfacts]','-');
        set(h,'Color',[.8 .8 .8]);
        h=plot3([zeros(sum(Loog),1) -scaled(Loog,1)./normfacts]',[zeros(sum(Loog),1) -scaled(Loog,2)./normfacts]',[zeros(sum(Loog),1) -scaled(Loog,3)./normfacts]','-');
        set(h,'Color',[.8 .8 .8]);
    
        % Adjusting the axes
        set(gca,'XLim',xlim,'Ylim',ylim,'Zlim',zlim)
        for whichaxes = {'XTick','YTick','ZTick'}
            tmp = get(gca,char(whichaxes)); set(gca,char(whichaxes),[tmp(1) 0 tmp(end)]);
        end
         xticks = get(gca,'XTick');
         set(gca,'FontSize',8)
        %set(gca,'XTick',[],'YTick',[],'ZTick',[])
       
        if (PLOTPLANES)
            % Plotting plane fits
            DKLweights = planeparams'*xformmat';
            [x,y,z] = meshgrid(linspace(-plotlims(1),plotlims(1),50),linspace(-plotlims(2),plotlims(2),50),linspace(-plotlims(3),plotlims(3),50));
            v = abs(x.*DKLweights(1)+y.*DKLweights(2)+z.*DKLweights(3));
            fv = isosurface(x,y,z,v,1);
            
            % added by zack
            vertCounts = hist(fv.faces(:), length(fv.vertices)); % vertex count
            edgeVerts = [find(vertCounts==1) find(vertCounts==2) find(vertCounts==3)]; % get vertices that were only used 1-3 times
            
            % separate the vertices by plane
            p1 = []; p2 = [];
            for j = 1:length(edgeVerts)
                if DKLweights*fv.vertices(edgeVerts(j),:)' > 0 % the LHS will be ±1
                    p1 = [p1; edgeVerts(j)];
                else
                    p2 = [p2; edgeVerts(j)];
                end
            end
            % get the center point of both planes
            center1 = mean(fv.vertices(p1,:));
            center2 = mean(fv.vertices(p2,:));
            
            % thetas for each edge point about the center (per plane)
            theta1 = atan2(fv.vertices(p1,2)-center1(2),fv.vertices(p1,1)-center1(1));
            theta2 = atan2(fv.vertices(p2,2)-center2(2),fv.vertices(p2,1)-center2(1));
            
            % sort 'em
            [~,sortorder1] = sort(theta1);
            [~,sortorder2] = sort(theta2);
            
            ppoints1 = fv.vertices(p1(sortorder1),:);
            ppoints2 = fv.vertices(p2(sortorder2),:);
            
            % plot 'em
            h = patch(ppoints1(:,1),ppoints1(:,2),ppoints1(:,3),[0 1 0]);
            set(h,'FaceAlpha',0.5,'EdgeAlpha',0,'EdgeColor','none','LineWidth',1);
            h = patch(ppoints2(:,1),ppoints2(:,2),ppoints2(:,3),[0 1 0]);
            set(h,'FaceAlpha',0.5,'EdgeAlpha',0,'EdgeColor','none','LineWidth',1);
        end
            
        if (PLOTQUAD)
            % Plotting quadratic fits
            A =  [quadparams(1) quadparams(4) quadparams(5);...
                  quadparams(4) quadparams(2) quadparams(6);...
                  quadparams(5) quadparams(6) quadparams(3)];
            % Transforming quadparams into un-whitened space
            B = xformmat*A*xformmat';
            [x,y,z] = meshgrid(linspace(-plotlims(1),plotlims(1),30),linspace(-plotlims(2),plotlims(2),30),linspace(-plotlims(3),plotlims(3),30)); 
            variables = [x(:).^2 y(:).^2 z(:).^2 2*x(:).*y(:) 2*x(:).*z(:) 2*y(:).*z(:)];
            v = variables*[B(1,1) B(2,2) B(3,3) B(1,2) B(1,3) B(2,3)]';
            fv = isosurface(x,y,z,reshape(v,size(x)),1);
            h = patch(fv);
            set(h,'FaceColor',[.4 .7 .4],'EdgeAlpha',.5,'FaceAlpha',.5);
            set(h,'FaceColor','green','EdgeColor','none','LineWidth',1);
        end
    end
end
set(gcf,'Renderer','painters');

% Use the magic wand tool in Illustrator 
% to select all the green stuff and make it transparent.

%%
% Section 1.1
% Three example cells, two views of each.

LARGEAXISWIDTH = 1.3;

figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
colordef white;

markersize = 3;
axislabels = {'L+M','L-M','S'};

for col = 1:length(examplefilenames)
    filename = examplefilenames{col};
    % Getting best viewangles
    if (strcmp(filename, 'S041310002'))
        customview = [173 16; 41 18];
        plotlims = [.16; .16; 1];
    elseif (strcmp(filename, 'K070309004'))
        customview = [80 7;-136 27];
        plotlims = [.1; .1; 1];
    elseif (strcmp(filename, 'K082609009'))
        customview = [-115 11; -196 17];
        plotlims = [.09; .09; .5];
    elseif (strcmp(filename, 'K042811002'))
        customview = [173 16; 59 6];
        plotlims = [.16; .16; 1];
    else
        customview = [-102 -13];
    end
    
    stro = nex2stro(findfile(filename));
    lmsidxs = [find(strcmp(stro.sum.trialFields(1,:),'lcont'))...
        find(strcmp(stro.sum.trialFields(1,:),'mcont'))...
        find(strcmp(stro.sum.trialFields(1,:),'scont'))];
    ntrials = size(stro.trial(:,lmsidxs),1);
    % -------------------------------
    % Converting cone contrasts in nex file to 10 deg fundmentals.
    % Must go through excitations first - can't just transform contrasts.
    % -------------------------------
    fundamentals = stro.sum.exptParams.fundamentals;
    fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
    mon_spd = stro.sum.exptParams.mon_spd;
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
    M = fundamentals'*mon_spd;
    load ('T_cones_smj10');
    M10 = T_cones_smj10*mon_spd;
    stro.trial(:,lmsidxs) = ConvertConeContrastBasis(M, M10, stro.sum.exptParams.bkgndrgb, stro.trial(:,lmsidxs));
   % -------------------------------
    out = NTpreprocess(stro,0,Inf);  % Getting the termination points
    scaled = out(:,[2 3 4]).*repmat(out(:,5),1,3);
    % Moving to a DKL-ish space
    scaled(:,[1 2]) = scaled(:,[1 2])*mkbasis([1 1; 1 -1]);
    Loog = logical(out(:,end));
    [planeparams, ~, quadparams, ~, xformmat] = NTsurfacefit(scaled, Loog);
    spikeidx = strcmp(stro.sum.rasterCells(1,:),'sig001a');
    
    % ---------------------------
    % Plotting thresholds
    % ---------------------------
    for i = 1:2
        if i == 1
            colorder = [1 2 3];
        else
            colorder = [3 1 2];
            colorder = [2 3 1];
        end
        axes('Position',[1+(col-1)*LARGEAXISWIDTH*1.7 5.5-(i-1)*2*LARGEAXISWIDTH LARGEAXISWIDTH LARGEAXISWIDTH]);
        axis vis3d;
        hold on;
        axis square; set(gca,'View',customview(i,:));
        xlabel(axislabels(colorder(1))); set(get(gca,'XLabel'),'Color',[0 0 0]);
        ylabel(axislabels(colorder(2))); set(get(gca,'YLabel'),'Color',[0 0 0]);
        zlabel(axislabels(colorder(3))); set(get(gca,'ZLabel'),'Color',[0 0 0]);
        h=plot3(scaled(~Loog,colorder(1)),scaled(~Loog,colorder(2)),scaled(~Loog,colorder(3)),'ko');
        set(h,'Markersize',markersize,'Markerfacecolor','black');
        h=plot3(-scaled(~Loog,colorder(1)),-scaled(~Loog,colorder(2)),-scaled(~Loog,colorder(3)),'ko');
        set(h,'Markersize',markersize,'Markerfacecolor','black');
        
        % Plotting OOG rays (squeezing them inside the axis limits)
        normfacts = max(abs(scaled(Loog,:)./repmat(plotlims',sum(Loog),1)),[],2);
        normfacts = max(normfacts,1);
        h=plot3([zeros(sum(Loog),1) scaled(Loog,colorder(1))./normfacts]',[zeros(sum(Loog),1) scaled(Loog,colorder(2))./normfacts]',[zeros(sum(Loog),1) scaled(Loog,colorder(3))./normfacts]','-');
        set(h,'Color',[.8 .8 .8]);
        h=plot3([zeros(sum(Loog),1) -scaled(Loog,colorder(1))./normfacts]',[zeros(sum(Loog),1) -scaled(Loog,colorder(2))./normfacts]',[zeros(sum(Loog),1) -scaled(Loog,colorder(3))./normfacts]','-');
        set(h,'Color',[.8 .8 .8]);
        
        % Adjusting the axes
        set(gca,'XLim',[-1 1]*plotlims(colorder(1)),'Ylim',[-1 1]*plotlims(colorder(2)),'Zlim',[-1 1]*plotlims(colorder(3)))
        for whichaxes = {'XTick','YTick','ZTick'}
            tmp = get(gca,char(whichaxes)); set(gca,char(whichaxes),[tmp(1) 0 tmp(end)]);
        end
        xticks = get(gca,'XTick');
        set(gca,'FontSize',8)
        
        % Plotting plane fits
        DKLweights = planeparams'*xformmat';
        DKLweights = DKLweights(colorder);
        [x,y,z] = meshgrid(linspace(-plotlims(colorder(1)),plotlims(colorder(1)),50),linspace(-plotlims(colorder(2)),plotlims(colorder(2)),50),linspace(-plotlims(colorder(3)),plotlims(colorder(3)),50));
        v = abs(x.*DKLweights(1)+y.*DKLweights(2)+z.*DKLweights(3));
        fv = isosurface(x,y,z,v,1);
        
        % added by zack
        vertCounts = hist(fv.faces(:), length(fv.vertices)); % vertex count
        edgeVerts = [find(vertCounts==1) find(vertCounts==2) find(vertCounts==3)]; % get vertices that were only used 1-3 times
        
        % separate the vertices by plane
            p1 = []; p2 = [];
            for j = 1:length(edgeVerts)
                if DKLweights*fv.vertices(edgeVerts(j),:)' > 0 % the LHS will be ±1
                    p1 = [p1; edgeVerts(j)];
                else
                    p2 = [p2; edgeVerts(j)];
                end
            end
            % get the center point of both planes
            center1 = mean(fv.vertices(p1,:));
            center2 = mean(fv.vertices(p2,:));
            
            % thetas for each edge point about the center (per plane)
            theta1 = atan2(fv.vertices(p1,2)-center1(2),fv.vertices(p1,1)-center1(1));
            theta2 = atan2(fv.vertices(p2,2)-center2(2),fv.vertices(p2,1)-center2(1));
            
            % sort 'em
            [~,sortorder1] = sort(theta1);
            [~,sortorder2] = sort(theta2);
            
            ppoints1 = fv.vertices(p1(sortorder1),:);
            ppoints2 = fv.vertices(p2(sortorder2),:);
            
            % plot 'em
            h = patch(ppoints1(:,1),ppoints1(:,2),ppoints1(:,3),[0 1 0]);
            set(h,'FaceAlpha',1,'EdgeAlpha',0,'EdgeColor','none','LineWidth',1);
            h = patch(ppoints2(:,1),ppoints2(:,2),ppoints2(:,3),[0 1 0]);
            set(h,'FaceAlpha',1,'EdgeAlpha',0,'EdgeColor','none','LineWidth',1);
        end
    end    
%set(gcf,'Renderer','painters');

% Use the magic wand tool in Illustrator 
% to select all the green stuff and make it transparent.


%%
% Section 1.2: 3-D fits+data points for the three example cells
% Plot the surfaces only (in the same size as the plot with the data)
% Manually superimpose (in Illustrator) the data (saved as postscript)
% as the surface saved as TIFF. Need to use OpenGL renderer.

LARGEAXISWIDTH = 1.3;
PLOTSURF = 0;
PLOTPLANES = 1;
PLOTDATA = 1;
figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
colordef white;

markersize = 4;
axislabels = {'L+M','L-M','S'};

for col = 1:length(examplefilenames)
    filename = examplefilenames{col};
    % Getting best viewangles
    if (strcmp(filename, 'S041310002'))
        customview = [173 16; 41 18];
        plotlims = [.16; .16; 1];
    elseif (strcmp(filename, 'K070309004'))
        customview = [80 7;-136 27];
        plotlims = [.1; .1; 1];
    elseif (strcmp(filename, 'K082609009'))
        customview = [-115 11; -196 17];
        plotlims = [.09; .09; .5];
    elseif (strcmp(filename, 'K042811002'))
        customview = [173 16; 59 6];
        plotlims = [.16; .16; 1];
    else
        customview = [-102 -13];
    end
    
    stro = nex2stro(findfile(filename));
    lmsidxs = [find(strcmp(stro.sum.trialFields(1,:),'lcont'))...
        find(strcmp(stro.sum.trialFields(1,:),'mcont'))...
        find(strcmp(stro.sum.trialFields(1,:),'scont'))];
    ntrials = size(stro.trial(:,lmsidxs),1);
    % -------------------------------
    % Converting cone contrasts in nex file to 10 deg fundmentals.
    % Must go through excitations first - can't just transform contrasts.
    % -------------------------------
    fundamentals = stro.sum.exptParams.fundamentals;
    fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
    mon_spd = stro.sum.exptParams.mon_spd;
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
    M = fundamentals'*mon_spd;
    load ('T_cones_smj10');
    M10 = T_cones_smj10*mon_spd;
    stro.trial(:,lmsidxs) = ConvertConeContrastBasis(M, M10, stro.sum.exptParams.bkgndrgb, stro.trial(:,lmsidxs));
   % -------------------------------
    out = NTpreprocess(stro,0,Inf);  % Getting the termination points
    scaled = out(:,[2 3 4]).*repmat(out(:,5),1,3);
    % Moving to a DKL-ish space
    scaled(:,[1 2]) = scaled(:,[1 2])*mkbasis([1 1; 1 -1]);
    Loog = logical(out(:,end));
    [planeparams, ~, quadparams, ~, xformmat] = NTsurfacefit(scaled, Loog);
    spikeidx = strcmp(stro.sum.rasterCells(1,:),'sig001a');
    
    % ---------------------------
    % Plotting thresholds
    % ---------------------------
    for i = 1:2
        if i == 1
            colorder = [1 2 3];
        else
            colorder = [2 3 1];
        end
        axes('Position',[1+(i-1)*1.5*LARGEAXISWIDTH 6.5-(col-1)*1.5*LARGEAXISWIDTH LARGEAXISWIDTH LARGEAXISWIDTH]);
        axis vis3d;
        hold on;
        axis square; set(gca,'View',customview(i,:));
        axis square; set(gca,'View',customview(i,:));
        xlabel(axislabels(colorder(1))); set(get(gca,'XLabel'),'Color',[0 0 0]);
        ylabel(axislabels(colorder(2))); set(get(gca,'YLabel'),'Color',[0 0 0]);
        zlabel(axislabels(colorder(3))); set(get(gca,'ZLabel'),'Color',[0 0 0]);
        
        if (PLOTDATA)
            h=plot3(scaled(~Loog,colorder(1)),scaled(~Loog,colorder(2)),scaled(~Loog,colorder(3)),'ko');
            set(h,'Markersize',markersize,'Markerfacecolor','black');
            h=plot3(-scaled(~Loog,colorder(1)),-scaled(~Loog,colorder(2)),-scaled(~Loog,colorder(3)),'ko');
            set(h,'Markersize',markersize,'Markerfacecolor','black');
            % Plotting OOG rays (squeezing them inside the axis limits)
            normfacts = max(abs(scaled(Loog,:)./repmat(plotlims',sum(Loog),1)),[],2);
            normfacts = max(normfacts,1);
            h=plot3([zeros(sum(Loog),1) scaled(Loog,colorder(1))./normfacts]',[zeros(sum(Loog),1) scaled(Loog,colorder(2))./normfacts]',[zeros(sum(Loog),1) scaled(Loog,colorder(3))./normfacts]','-');
            set(h,'Color',[.7 .7 .7]);
            h=plot3([zeros(sum(Loog),1) -scaled(Loog,colorder(1))./normfacts]',[zeros(sum(Loog),1) -scaled(Loog,colorder(2))./normfacts]',[zeros(sum(Loog),1) -scaled(Loog,colorder(3))./normfacts]','-');
            set(h,'Color',[.7 .7 .7]);     
        end
        % Adjusting the axes
        set(gca,'XLim',[-1 1]*plotlims(colorder(1)),'Ylim',[-1 1]*plotlims(colorder(2)),'Zlim',[-1 1]*plotlims(colorder(3)))
        for whichaxes = {'XTick','YTick','ZTick'}
            tmp = get(gca,char(whichaxes)); set(gca,char(whichaxes),[tmp(1) 0 tmp(end)]);
        end
        xticks = get(gca,'XTick');
        set(gca,'FontSize',8)
        
        if (PLOTSURF)
            % Plotting surface fits
      %      if all(colorder == [1 2 3])
                A = [quadparams(1) quadparams(4) quadparams(5);...
                    quadparams(4) quadparams(2) quadparams(6);...
                    quadparams(5) quadparams(6) quadparams(3)];
      %      elseif all(colorder == [2 3 1])
      %          A = [quadparams(2) quadparams(6) quadparams(4);...
      %              quadparams(6) quadparams(3) quadparams(5);...
      %              quadparams(4) quadparams(5) quadparams(1)];
      %      end
            
       %     B =  xformmat(colorder,:)*A*xformmat(colorder,:)';
             B =  xformmat*A*xformmat';
            
            % Plot surfaces
            [xx yy zz] = meshgrid(linspace(-plotlims(1)*.8,plotlims(1)*.8,100),...
                linspace(-plotlims(2)*.8,plotlims(2)*.8,100),...
                linspace(-plotlims(3)*.8,plotlims(3)*.8,100));
            xformedxyz = [xx(:) yy(:) zz(:)];
            xformedxyz = xformedxyz(:,colorder);
            
            variables = [xformedxyz(:,1).^2 xformedxyz(:,2).^2 xformedxyz(:,3).^2 2*xformedxyz(:,1).*xformedxyz(:,2) 2*xformedxyz(:,1).*xformedxyz(:,3) 2*xformedxyz(:,2).*xformedxyz(:,3)];
            coefficients = [B(colorder(1),colorder(1)) B(colorder(2),colorder(2)) B(colorder(3),colorder(3)),...
                            B(colorder(1),colorder(2)) B(colorder(1),colorder(3)) B(colorder(2),colorder(3))]';
            
            surfcolor = [0 1 0];
            fr = variables*coefficients;
            if all(colorder == [1 2 3]);
                fv = isosurface(xx,yy,zz,reshape(fr,size(xx)), 1);
            else
                fv = isosurface(yy,zz,xx,reshape(fr,size(xx)), 1);
            end
            h = patch(fv);
            %        set(h,'FaceAlpha',0.5,'EdgeAlpha',0.5);
            set(h,'FaceVertexCData',repmat(surfcolor,size(fv.vertices,1),1))
            set(h,'CDataMapping','direct');
            set(h,'FaceColor','interp','EdgeColor','none');
            lighting gouraud
            set(h,'SpecularColorReflectance',.1);
            set(h,'SpecularExponent',.2,'SpecularStrength',.2);
            set(h,'DiffuseStrength',.8,'AmbientStrength',.2);
            
            if (strcmp(filename, 'S041310002'))
                camlight(0,0);
                camlight(0,0);
                light('position',[0 1 1])
            elseif (strcmp(filename, 'K070309004'))
                camlight(45,0);
               % camlight(-45,0);
                light('position',[0 1 1])
            elseif (strcmp(filename, 'K082609009'))
                camlight(45,0);
                camlight(-45,0);
            elseif (strcmp(filename, 'K042811002'))
                camlight(0,0);
                camlight(0,0);
                light('position',[0 1 1])
            end
        end
        if (PLOTPLANES)
            % Plotting plane fits
            DKLweights = planeparams'*xformmat';
            DKLweights = DKLweights(colorder);
            [x,y,z] = meshgrid(linspace(-plotlims(colorder(1)),plotlims(colorder(1)),50),linspace(-plotlims(colorder(2)),plotlims(colorder(2)),50),linspace(-plotlims(colorder(3)),plotlims(colorder(3)),50));
            v = abs(x.*DKLweights(1)+y.*DKLweights(2)+z.*DKLweights(3));
            fv = isosurface(x,y,z,v,1);
            % added by zack
            vertCounts = hist(fv.faces(:), length(fv.vertices)); % vertex count
            edgeVerts = [find(vertCounts==1) find(vertCounts==2) find(vertCounts==3)]; % get vertices that were only used 1-3 times
            
            % separate the vertices by plane
            p1 = []; p2 = [];
            for j = 1:length(edgeVerts)
                if DKLweights*fv.vertices(edgeVerts(j),:)' > 0 % the LHS will be ±1
                    p1 = [p1; edgeVerts(j)];
                else
                    p2 = [p2; edgeVerts(j)];
                end
            end
            % get the center point of both planes
            center1 = mean(fv.vertices(p1,:));
            center2 = mean(fv.vertices(p2,:));
            
            % thetas for each edge point about the center (per plane)
            theta1 = atan2(fv.vertices(p1,2)-center1(2),fv.vertices(p1,1)-center1(1));
            theta2 = atan2(fv.vertices(p2,2)-center2(2),fv.vertices(p2,1)-center2(1));
            
            % sort 'em
            [~,sortorder1] = sort(theta1);
            [~,sortorder2] = sort(theta2);
            
            ppoints1 = fv.vertices(p1(sortorder1),:);
            ppoints2 = fv.vertices(p2(sortorder2),:);
            
            % plot 'em
            h = patch(ppoints1(:,1),ppoints1(:,2),ppoints1(:,3),[0 1 0]);
       %     set(h,'FaceAlpha',0.5,'EdgeAlpha',0,'EdgeColor','none','LineWidth',1);
            h = patch(ppoints2(:,1),ppoints2(:,2),ppoints2(:,3),[0 1 0]);
       %     set(h,'FaceAlpha',0.5,'EdgeAlpha',0,'EdgeColor','none','LineWidth',1);
        end
    end 
end
% Copy the green stuff in photoshop and paste it behind the data in
% illustrator

set(gcf,'Renderer','openGL');
set(gcf,'InvertHardCopy','off')
print -dtiff -r600 -cmyk junk2

%% Section 2) Distribution of baseline firing rates and thresholds.
PLOTREDLINES = 0;

[fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\Gratings&NT2.txt');
data = [];
for cellcounter = 1:size(fnames,1)
    NT = {}; GT = {};
    for i = 1:size(fnames{cellcounter},2)
        filename = findfile(char(fnames{cellcounter}(i)));
        paradigmID = getparadigmID(filename);
        if paradigmID == 150
            GT = nex2stro(filename);
            GTstruct = getGratingTuning(GT,1);
        else
            NT = nex2stro(filename);
        end
    end
    
    fpacq_t = NT.trial(:,strcmp(NT.sum.trialFields(1,:),'fp_acq'));
    stimon_t = NT.trial(:,strcmp(NT.sum.trialFields(1,:),'stim_on'));
    L = logical(~isnan(fpacq_t));
    prestimfr = [];
    time = [];
    for i = find(L)'
        spiketimes = NT.ras{i,spikeIdx(cellcounter)};
        nspikes = sum(spiketimes > fpacq_t(i) & spiketimes < stimon_t(i));
        prestimfr = [prestimfr; nspikes./(stimon_t(i)-fpacq_t(i))];
        time = [time; stimon_t(i)-fpacq_t(i)];
    end
    GTmaxresp = max([GTstruct.areasummation.resp(:,1); GTstruct.color.colresp(:,1)]);
    
    data = [data; mean(prestimfr) std(prestimfr) NT.sum.exptParams.threshold mean(time) GTmaxresp];
end

% Plotting (linear)
figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
colordef white;
axes('position',[3,3,3,3]); hold on

% Red lines
if (PLOTREDLINES)
    h = plot([data(:,2) data(:,2)]',[data(:,1) data(:,4)]','Color',[1 .8 .8],'LineWidth',1);  % faint red lines
    upperbound = ceil(max(data(:,4))/10)*10;
    upperbound = 150; 
    tickint = 50;
else
    upperbound = ceil(max(data(:,3))/10)*10;
    tickint = 10;
end
set(gca,'XTick',[0:tickint:upperbound],'YTick',[0:tickint:upperbound]);
[n1,x1] = hist(data(:,1),linspace(0,upperbound,20))
[n2,x2] = hist(data(:,3),linspace(0,upperbound,20))
binwidth = x2(2)-x2(1);
h = plot(data(:,3),data(:,1),'ko','LineWidth',.01)
set(h,'MarkerSize',5,'MarkerFaceColor','black');
set(gca,'Xlim',[0-binwidth/2 upperbound],'Ylim',[0-binwidth/2 upperbound]);
set(gca,'XTickLabel',[],'YTickLabel',[]);
plot([0 upperbound],[0 upperbound],'k-.');
% +/- 1 SD
for i = 1:size(data,1)
    plot([data(i,3) data(i,3)],data(i,1)+[data(i,2) -data(i,2)],'k-','LineWidth',1)
end

% Marginal histogram of thresholds
axes('position',[3,1.5,3,1]);
h = bar(x2,n2);
set(h,'FaceColor','black');
set(gca,'Xlim',[0-binwidth/2 upperbound]);
set(gca,'XTick',[0:tickint:upperbound]);
xlabel('Target rate (sp/sec)','FontSize',14);

% Marginal histogram of baselines
axes('position',[1.5,3,1,3]);
if (PLOTREDLINES)
    hold on;
    bins = linspace(0,upperbound,20);
    binwidth = bins(2)-bins(1);
    [n3,x3] = hist(data(:,4),[bins bins(end)+binwidth])
    h = barh(x3(1:end-1),n3(1:end-1)+n1);  % Want stacked histograms, so adding n1
    set(h,'FaceColor',[1 .8 .8]);
end
h = barh(x1,n1);
set(h,'FaceColor','black');
set(gca,'Ylim',[0-binwidth/2 upperbound]);
set(gca,'YTick',[0:tickint:upperbound]);
set(gca,'XLim',[0 ceil(max(n3(1:end-1)+n1)/10)*10])
ylabel('Baseline rate (sp/sec)','FontSize',14);

sum(data(:,4) > upperbound)  % How many red bars go off the top of the plot?

% -------------
% Plotting (log)
% This didn't work out that well.  All the points get bunched up in a
% different way.
% -------------

% lowerbound = 0.001;
% upperbound = ceil(max(data(:,2))/10)*10;
% moddata = [max(data(:,1),lowerbound),data(:,2)];
% figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
% set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
% set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
% set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
% colordef white;
% 
% [n1,x1] = hist(moddata(:,1),logspace(log10(lowerbound),log10(upperbound),20))
% [n2,x2] = hist(moddata(:,2),logspace(log10(lowerbound),log10(upperbound),20))
% binwidth = x2(2)-x2(1);
% 
% axes('position',[3,3,3,3]); hold on
% h = loglog(moddata(:,2),moddata(:,1),'ko')
% set(h,'MarkerSize',3,'MarkerFaceColor','black');
% set(gca,'Xlim',[lowerbound upperbound],'Ylim',[lowerbound upperbound]);
% set(gca,'XScale','log','YScale','log');
% plot([1 upperbound],[1 upperbound],'k-.');
% 


%% Section 2.1) distribution of baselines, target FRs and maximums
% X axis is some arbitrary cell index (sorted)

[fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\Gratings&NT2.txt');
[fnamesMult, spikeIdxMult] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\NTmultiples.txt')
data = [];
shortfn = cell(1,size(fnames,1))
for cellcounter = 1:size(fnames,1)
    NT = {}; GT = {};
    for i = 1:size(fnames{cellcounter},2)
        filename = findfile(char(fnames{cellcounter}(i)));
        paradigmID = getparadigmID(filename);
        if paradigmID == 150
            GT = nex2stro(filename);
            GTstruct = getGratingTuning(GT,1);
        else
            NT = nex2stro(filename);
        end
    end
    
    fpacq_t = NT.trial(:,strcmp(NT.sum.trialFields(1,:),'fp_acq'));
    stimon_t = NT.trial(:,strcmp(NT.sum.trialFields(1,:),'stim_on'));
    L = logical(~isnan(fpacq_t));
    prestimfr = [];
    for i = find(L)'
        spiketimes = NT.ras{i,1};
        nspikes = sum(spiketimes > fpacq_t(i) & spiketimes < stimon_t(i));
        prestimfr = [prestimfr; nspikes./(stimon_t(i)-fpacq_t(i))];
    end
    GTmaxresp = max([GTstruct.areasummation.resp(:,1); GTstruct.color.colresp(:,1)]);
    shortfn{cellcounter} = NT.sum.fileName(find(NT.sum.fileName == filesep,1,'last')+1:find(NT.sum.fileName == '.',1,'last')-1); 
    data = [data; mean(prestimfr) NT.sum.exptParams.threshold GTmaxresp];
end

% Getting the multiples
multiplesdata = [];
fnameidxs = zeros(size(fnamesMult,1),1);
for i = 1:size(fnamesMult,1)
    L = 0;
    for j = 1:2
        L = L | strcmp(char(fnamesMult{i}(j)), shortfn);
    end
    if (any(L))
        fnameidxs(i) = find(L);
    else
        keyboard
        char(fnamesMult{i}(j))
    end

    for j = 1:2
        NT = nex2stro(findfile(char(fnamesMult{i}(j))));
        fpacq_t = NT.trial(:,strcmp(NT.sum.trialFields(1,:),'fp_acq'));
        stimon_t = NT.trial(:,strcmp(NT.sum.trialFields(1,:),'stim_on'));
        whichtrials = logical(~isnan(fpacq_t));
        prestimfr = [];
        for k = find(whichtrials)'
            spiketimes = NT.ras{k,spikeIdx(cellcounter)};
            nspikes = sum(spiketimes > fpacq_t(k) & spiketimes < stimon_t(k));
            prestimfr = [prestimfr; nspikes./(stimon_t(k)-fpacq_t(k))];
        end
        multiplesdata = [multiplesdata; mean(prestimfr) NT.sum.exptParams.threshold fnameidxs(i)]
    end
end

% Plotting
figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
colordef white;
axes('position',[3,3,3,3]); hold on
[~,sortedidxs] = sort(data(:,2));
counter = 1;
for i = sortedidxs'
    plot([counter counter], [data(i,1) data(i,3)],'k-');
    plot(counter, data(i,2),'m*');
    L = multiplesdata(:,end) == i;
    for j = find(L)'
        plot(counter, multiplesdata(j,2),'m*');
        keyboard
    end
    counter = counter + 1;
end
ylabel('Firing rate (sp/sec)','FontSize',14);
xlabel('Cell number','FontSize',14);
set(gca,'XTick',[]);

%% Section 3) Distribution of p-values from bootstrap test of planarity
BOOTSTRAP = 0;
if BOOTSTRAP
    nbootiter = 200;
else 
    nbootiter = 0;
end
ps = [];
surfacetypes = [];
[fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\Gratings&NT2.txt');
for cellcounter = 1:size(fnames,1)
    NT = {};
    for i = 1:size(fnames{cellcounter},2)
        filename = findfile(char(fnames{cellcounter}(i)));
        paradigmID = getparadigmID(filename);
        if paradigmID == 103
            NT = nex2stro(filename);
        else
            continue;
        end
    end
    
    out = NTpreprocess(NT,0,Inf);
    
    scaled = out(:,[2:4]).*repmat(out(:,5), 1,3);
    Loog = logical(out(:,7));
    [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);
    
    % What type of quadric?
    A = [quadparams(1) quadparams(4) quadparams(5);...
        quadparams(4) quadparams(2) quadparams(6);...
        quadparams(5) quadparams(6) quadparams(3)];
    [evecs,evals] = eig(A);
    [evals,i] = sort(diag(evals),1,'ascend');
    WHICHSURF = sum(evals<0);
    
    xyz = scaled*xformmat;
    % Finding r's that lie on the best fit plane.
    [th,ph,r] = cart2sph(xyz(:,1),xyz(:,2),xyz(:,3));
 	planer = abs(-1./(planeparams(1).*cos(ph).*cos(th)+planeparams(2).*cos(ph).*sin(th)+planeparams(3).*sin(ph)));
    resid = log(r)-log(planer);
    
    
    if (BOOTSTRAP)
        nulldist = nan*ones(nbootiter,1);
        wait_h = waitbar(0,'Bootstrapping...');
        disp(['Bootstrapping.  PlaneSSE = ',num2str(planeSSE),' QuadSSE = ',num2str(quadSSE)]);
        SSEs = [];
        for j = 1:nbootiter
            waitbar(j/nbootiter, wait_h);
            tmpresid = exp(resid(unidrnd(length(resid),[1 length(resid)])));
            
            tmpx = xyz(:,1);
            tmpy = xyz(:,2);
            tmpz = xyz(:,3);
            
            tmpx(~Loog) = planer(~Loog).*tmpresid(~Loog).*cos(ph(~Loog)).*cos(th(~Loog));
            tmpy(~Loog) = planer(~Loog).*tmpresid(~Loog).*cos(ph(~Loog)).*sin(th(~Loog));
            tmpz(~Loog) = planer(~Loog).*tmpresid(~Loog).*sin(ph(~Loog));
            
            options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-6,'TolX',10^-6);
            [tmp, SSE1, exitflag1] = fminsearch(@(x) surfacefiterr2([tmpx tmpy tmpz],x, Loog),planeparams, options);
            initguess = [tmp'.*tmp' tmp(1)*tmp(2) tmp(1)*tmp(3) tmp(2)*tmp(3)];
            [tmp, SSE2, exitflag2] = fminsearch(@(x) surfacefiterr4([tmpx tmpy tmpz],x, Loog),initguess, options);
            if (~exitflag1)
                disp('Bonked on a plane fit');
            end
            if (~exitflag2)
                disp('Bonked on a quadratic fit');
            end
            if (~exitflag1 | ~exitflag2)
                j = j-1;
                keyboard
            end
            
            % Calling NTsurfacefit on every iteration is painfully slow.  Probably not necessary.
            % [planeparams, SSE1, quadparams, SSE2, xformmat] = NTsurfacefit([tmpx tmpy tmpz], Loog);
            nulldist(j)= SSE1./SSE2;
            SSEs(j,:) = [SSE1 SSE2];
        end
        close(wait_h);
        ps = [ps; sum(nulldist>planeSSE/quadSSE)./nbootiter]
    else
        F = ((planeSSE-quadSSE)/3)/(quadSSE/(sum(~Loog)-6))
        pf = 1-fcdf(F,3,sum(~Loog)-3);
        ps = [ps; pf]
    end
    surfacetypes = [surfacetypes; WHICHSURF];
end

figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
colordef white;

axes('position',[2 2 5 2]);
hold on;
[n,x] = hist(ps,linspace(0,1,20));
binwidth = x(2)-x(1); 
h = bar(x,n);
set(h,'FaceColor','black');
set(gca,'Xlim',[-binwidth/2 1+binwidth/2]);
xlabel('Count','FontSize',14);
ylabel('P-value','FontSize',14);

% Arrows for example cells - Yuck, all the good data have p = 0
% yaxlims = get(gca,'Ylim');
% for i = 1:length(examplefilenames)
%     examplecellidx = find(ismember([fnames{:}],examplefilenames{i}))/2
%     % sanity check
%    	plot(ps(examplecellidx),yaxlims(2),'m*')
%     text(ps(examplecellidx),yaxlims(2),num2str(i),'FontWeight','bold');
% end

% Getting pseudo-Rsquare
for i = 1:length(examplefilenames)
    filename = findfile(examplefilenames{i});
    NT = nex2stro(filename);
    out = NTpreprocess(NT,0,Inf);
    
    scaled = out(:,[2:4]).*repmat(out(:,5), 1,3);
    Loog = logical(out(:,7));
    [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);
    
    [th,ph,r] = cart2sph(scaled(:,1),scaled(:,2),scaled(:,3));
    SST = sum((log(r)-mean(log(r))).^2);
    PR2_quad = 1-quadSSE/SST
    PR2_plane = 1-planeSSE/SST
end

%% Section 4) Cone weights for all cells with planar cells highlighted

GRAY = [.7 .7 .7];
data = [];
[planefnames] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\plane0.01.txt');
[fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\Gratings&NT2.txt');
for cellcounter = 1:size(fnames,1)
    NT = {}; 
    for i = 1:size(fnames{cellcounter},2)
        filename = findfile(char(fnames{cellcounter}(i)));
        paradigmID = getparadigmID(filename);
        if paradigmID == 103
            NT = nex2stro(filename);
        else
            continue;
        end
    end

    % -------------------------------
    % Converting cone contrasts in nex file to 10 deg fundmentals.
    % Must go through excitations first - can't just transform contrasts.
    % -------------------------------
    lmsidxs = [find(strcmp(NT.sum.trialFields(1,:),'lcont'))...
        find(strcmp(NT.sum.trialFields(1,:),'mcont'))...
        find(strcmp(NT.sum.trialFields(1,:),'scont'))];
    fundamentals = NT.sum.exptParams.fundamentals;
    fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
    mon_spd = NT.sum.exptParams.mon_spd;
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
    M = fundamentals'*mon_spd;
    load ('T_cones_smj10');
    M10 = T_cones_smj10*mon_spd;
    NT.trial(:,lmsidxs) = ConvertConeContrastBasis(M, M10, NT.sum.exptParams.bkgndrgb, NT.trial(:,lmsidxs));
   % -------------------------------

    out = NTpreprocess(NT,0,Inf);  % .4, 2  or 0, Inf
    scaled = out(:,[2:4]).*repmat(out(:,5), 1,3);
    Loog = logical(out(:,7));
    [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);
    coneweights = (planeparams'*xformmat')./sum(abs(planeparams'*xformmat'));
    if (coneweights(2) < 0)   % convention: positive M-cone weight
        coneweights = -coneweights;
    end
    shortfn = NT.sum.fileName(find(NT.sum.fileName == filesep,1,'last')+1:find(NT.sum.fileName == '.',1,'last')-1);
    data = [data; coneweights ismember(shortfn,[planefnames{:}]')];
end

figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
colordef white;

axes('position',[2 2 4*sqrt(2) 4]);
hold on;
LnegS = data(:,3) < 0;
Lplane = logical(data(:,4));
plot(data(~Lplane&~LnegS,1),data(~Lplane&~LnegS,2),'o','MarkerSize',7,'MarkerFaceColor',GRAY,'MarkerEdgeColor',GRAY);
plot(data(~Lplane&LnegS,1),data(~Lplane&LnegS,2),'o','MarkerSize',7,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',GRAY);
plot(data(Lplane&~LnegS,1),data(Lplane&~LnegS,2),'o','MarkerSize',7,'MarkerFaceColor','black','MarkerEdgeColor','black');
plot(data(Lplane&LnegS,1),data(Lplane&LnegS,2),'o','MarkerSize',7,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor','black');

plot([-1 1 0 -1],[0 0 1 0],'k-');
xlabel('L-cone weight','FontSize',14);
ylabel('M-cone weight','FontSize',14);

% Some stats
Lnonop = data(:,1)+data(:,2) > .9;
Lop = data(:,1) < 0 & abs(data(:,1)-data(:,2)) > .9; 
%plot(data(Lnonop,1),data(Lnonop,2),'ro')
%plot(data(Lop,1),data(Lop,2),'go')

%       Lplane ~Lplane
% L+M     a       b
% L-M     c       d

a = sum(Lnonop & data(:,4));
b = sum(Lnonop & ~data(:,4));
c = sum(Lop & data(:,4));
d = sum(Lop & ~data(:,4));

disp(['prop L+M that are planar: ',num2str(a./(a+b))])
disp(['prop L-M that are planar: ',num2str(c./(c+d))])
[h, p] = fisherExact(a,b,c,d,0.05)


for i = 1:length(examplefilenames)
    examplecellidx = find(ismember([fnames{:}],examplefilenames{i}))/2;
    % sanity check
    [fnames{examplecellidx}(2), examplefilenames(i)]
    h = text(data(examplecellidx,1),data(examplecellidx,2),num2str(i),'FontWeight','bold');
    set(h,'Color','red');
end

% Where are the S-(L+M) cells?
L = data(:,1)+data(:,2)> .5 & data(:,3) < -.1;
plot(data(L,1),data(L,2),'m^');

%% Section 4.1: Revised version of cone weight plot highlighting
% ellipsoidal cells

GRAY1 = [.7 .7 0];
GRAY2 = [.7 .7 .7];
data = [];
[planefnames] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\plane0.01.txt');
%[pancolorfnames] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\pancolor2.txt');
[fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\Gratings&NT2.txt');
for cellcounter = 1:size(fnames,1)
    NT = {}; 
    for i = 1:size(fnames{cellcounter},2)
        filename = findfile(char(fnames{cellcounter}(i)));
        paradigmID = getparadigmID(filename);
        if paradigmID == 103
            NT = nex2stro(filename);
        else
            continue;
        end
    end

    % -------------------------------
    % Converting cone contrasts in nex file to 10 deg fundmentals.
    % Must go through excitations first - can't just transform contrasts.
    % -------------------------------
    lmsidxs = [find(strcmp(NT.sum.trialFields(1,:),'lcont'))...
        find(strcmp(NT.sum.trialFields(1,:),'mcont'))...
        find(strcmp(NT.sum.trialFields(1,:),'scont'))];
    fundamentals = NT.sum.exptParams.fundamentals;
    fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
    mon_spd = NT.sum.exptParams.mon_spd;
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
    M = fundamentals'*mon_spd;
    load ('T_cones_smj10');
    M10 = T_cones_smj10*mon_spd;
    NT.trial(:,lmsidxs) = ConvertConeContrastBasis(M, M10, NT.sum.exptParams.bkgndrgb, NT.trial(:,lmsidxs));
   % -------------------------------

    out = NTpreprocess(NT,0,Inf);  % .4, 2  or 0, Inf
    scaled = out(:,[2:4]).*repmat(out(:,5), 1,3);
    Loog = logical(out(:,7));
    [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);
    coneweights = (planeparams'*xformmat')./sum(abs(planeparams'*xformmat'));
    if (coneweights(2) < 0)   % convention: positive M-cone weight
        coneweights = -coneweights;
    end
    shortfn = NT.sum.fileName(find(NT.sum.fileName == filesep,1,'last')+1:find(NT.sum.fileName == '.',1,'last')-1);
    A = [quadparams(1) quadparams(4) quadparams(5);...
        quadparams(4) quadparams(2) quadparams(6);...
        quadparams(5) quadparams(6) quadparams(3)];
    isellipsoid = all(eig(A)>0) & ~ismember(shortfn,[planefnames{:}]');
    
    data = [data; coneweights isellipsoid ismember(shortfn,[planefnames{:}]')];
end

figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
colordef white;

axes('position',[2 2 4*sqrt(2) 4]);
hold on;
LnegS = data(:,3) < 0;
Lellipsoid = logical(data(:,4));
Lplane = logical(data(:,5));


plot(data(Lellipsoid&~LnegS,1),data(Lellipsoid&~LnegS,2),'o','MarkerSize',7,'MarkerFaceColor',GRAY1,'MarkerEdgeColor',GRAY1);
plot(data(Lellipsoid&LnegS,1),data(Lellipsoid&LnegS,2),'o','MarkerSize',7,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',GRAY1);

plot(data(~Lplane&~Lellipsoid&~LnegS,1),data(~Lplane&~Lellipsoid&~LnegS,2),'o','MarkerSize',7,'MarkerFaceColor',GRAY2,'MarkerEdgeColor',GRAY2);
plot(data(~Lplane&~Lellipsoid&LnegS,1),data(~Lplane&~Lellipsoid&LnegS,2),'o','MarkerSize',7,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',GRAY2);

plot(data(Lplane&~LnegS,1),data(Lplane&~LnegS,2),'o','MarkerSize',7,'MarkerFaceColor','black','MarkerEdgeColor','black');
plot(data(Lplane&LnegS,1),data(Lplane&LnegS,2),'o','MarkerSize',7,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor','black');

plot([-1 1 0 -1],[0 0 1 0],'k-');
xlabel('L-cone weight','FontSize',14);
ylabel('M-cone weight','FontSize',14);

% Any ellipsoids that weren't significantly nonplanar?
any(data(:,end)+data(:,end-1) == 2)
% No, whew.

for i = 1:length(examplefilenames)
    examplecellidx = find(ismember([fnames{:}],examplefilenames{i}))/2;
    % sanity check
    [fnames{examplecellidx}(2), examplefilenames(i)]
    h = text(data(examplecellidx,1),data(examplecellidx,2),num2str(i),'FontWeight','bold');
    set(h,'Color','red');
end

% Looking for S-(L+M) cells
L_s = abs(data(:,3)) > .5;
L_s = abs(data(:,3)) > .3;
sum(L_s)

%% Section 5: Example Gabors for the color space introduction plot

bkgndrgb = [.5 .5 .5];
colordirection = [.05 -.05 .2];  

M = [0.0608    0.1219    0.0175;
    0.0220    0.1266    0.0257;
    0.0019    0.0095    0.0976];
figure;
bkgndlms = M*bkgndrgb';
for i = 1:2
    colordirection = -colordirection;
    stimlms = (colordirection+1).*repmat(bkgndlms',size(colordirection,1),1);
    gaborrgb = (inv(M)*stimlms')';
    subplot(2,2,i);
    im = DrawGaborEdge(bkgndrgb, gaborrgb(j,:)-bkgndrgb, [0 0 0], pi/2, 1.5, 1, 1, 0, 0, 0, 0, 0, .999, 45);
    image(im);
        
    set(gca,'XTick',[],'YTick',[],'visible','off');
    axis square;
end

%% Section 6: Example grating responses in color space, isoresponse
% plane/quadratic surfaces.  Correlations between actual and predicted
% responses
% Move symbol legend to the top of the 3-d plot
% make subplots bigger
% Move large panels fartehr apart in a more principled way
% Move "r=XX" to title? 

LARGEAXISWIDTH = 1.3;
SMALLAXISWIDTH = 1;
COLMARGIN = 1;
PLOTPLANES = 1;
PLOTQUAD = 1;
DKLspace = 1;
PLOTSCATTERPLOTS = 1;
BITMAP = 0;

[fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\Gratings&NT2.txt');

% Setting up the figure
figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
colordef white;

for whichcell = 1:length(examplefilenames)
    filename = examplefilenames{whichcell};
    examplecellidx = find(ismember([fnames{:}],filename))/2;
    GT = nex2stro(findfile(fnames{examplecellidx}(1))); % Assumes GT preceed NT in text file
    GTstruct = getGratingTuning(GT);
    NT = nex2stro(findfile(fnames{examplecellidx}(2)));
    
    lmsidxs = [find(strcmp(NT.sum.trialFields(1,:),'lcont'))...
        find(strcmp(NT.sum.trialFields(1,:),'mcont'))...
        find(strcmp(NT.sum.trialFields(1,:),'scont'))];
    ntrials = size(NT.trial(:,lmsidxs),1);
    
    % Getting best viewangles
    if (strcmp(filename, 'S041310002'))
        viewangle = [161 22];
    elseif (strcmp(filename, 'K070309004'))
        viewangle = [81 7];
    elseif (strcmp(filename, 'K082609009'))
        viewangle = [-115 11];
    elseif (strcmp(filename, 'K042811002'))
        viewangle =  [173 16];
    else
        viewangle = [45 45];
    end
    
    % -------------------------------
    % Converting cone contrasts in nex file to 10 deg fundmentals.
    % Must go through excitations first - can't just transform contrasts.
    % -------------------------------
    fundamentals = NT.sum.exptParams.fundamentals;
    fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
    mon_spd = NT.sum.exptParams.mon_spd;
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
    M = fundamentals'*mon_spd;
    load ('T_cones_smj10');
    M10 = T_cones_smj10*mon_spd;
    NT.trial(:,lmsidxs) = ConvertConeContrastBasis(M, M10, NT.sum.exptParams.bkgndrgb, NT.trial(:,lmsidxs));
    % --------------------------------
    % Converting GT cone contrasts to 10 degree fundamentals
    % Safe to assume the M matrix and bkgndlms is the same for gratings and NT
    GTstruct.color.colors = ConvertConeContrastBasis(M, M10, NT.sum.exptParams.bkgndrgb, GTstruct.color.colors);
    % -------------------------------
    
    out = NTpreprocess(NT,0,Inf);  % Getting the termination points
    scaled = out(:,[2 3 4]).*repmat(out(:,5),1,3);
    if (DKLspace)
        scaled(:,[1 2]) = scaled(:,[1 2])*mkbasis([1 1; 1 -1]);
    end
    Loog = logical(out(:,end));
    [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);
    A = [quadparams(1) quadparams(4) quadparams(5);...
        quadparams(4) quadparams(2) quadparams(6);...
        quadparams(5) quadparams(6) quadparams(3)];
    coneweights = planeparams'*xformmat';
    B = xformmat*A*xformmat';
    
    % Plotting gratings data
    axes('position',[1+(whichcell-1)*(LARGEAXISWIDTH+COLMARGIN) 6 LARGEAXISWIDTH LARGEAXISWIDTH]); hold on;
    GTcols = GTstruct.color.colors;
    if (DKLspace)
        GTcols(:,[1 2]) = GTcols(:,[1 2])*mkbasis([1 1; 1 -1]);
    end
    maxresp = max(GTstruct.color.colresp(:,1));
    minresp = min(GTstruct.color.colresp(:,1));
    maxmarkersize = 12;
    minmarkersize = 2;
    b1 = (maxmarkersize-minmarkersize)./(maxresp-minresp);
    b0 = maxmarkersize-(maxresp*b1);
    clear h;
    for i = 1:size(GTstruct.color.colors,1)
        h(1) = plot3(GTcols(i,1),GTcols(i,2),GTcols(i,3),'ko');
        h(2) = plot3(-GTcols(i,1),-GTcols(i,2),-GTcols(i,3),'ko');
        set(h,'MarkerSize',GTstruct.color.colresp(i,1)*b1+b0);
        set(h,'MarkerfaceColor','black');  % stopgap
    end
    if (DKLspace)
        xlabel('L+M'); ylabel('L-M');
    else
        xlabel('L'); ylabel('M');
    end
    zlabel('S');
    axis vis3d;
    set(gca,'View',viewangle,'Xlim',[-.2 .2],'Ylim',[-.2 .2],'Zlim',[-1 1]);
    
    plotlim = max(abs(GTcols));
    % Plot surfaces
    [xx yy zz] = meshgrid(linspace(-plotlim(1),plotlim(1),100),...
        linspace(-plotlim(2),plotlim(2),100),...
        linspace(-plotlim(3),plotlim(3),100));
    xformedxyz = [xx(:) yy(:) zz(:)];
    
    variables = [xformedxyz(:,1).^2 xformedxyz(:,2).^2 xformedxyz(:,3).^2 2*xformedxyz(:,1).*xformedxyz(:,2) 2*xformedxyz(:,1).*xformedxyz(:,3) 2*xformedxyz(:,2).*xformedxyz(:,3)];
    for i = 1:2
        if (i == 1 & PLOTPLANES)
            coefficients = [coneweights.*coneweights coneweights(1)*coneweights(2) coneweights(1)*coneweights(3) coneweights(2)*coneweights(3)]';
            surfcolor = [.5 .5 .8];
        elseif (i == 2 & PLOTQUAD)
            coefficients = [B(1,1) B(2,2) B(3,3) B(1,2) B(1,3) B(2,3)]';
            surfcolor = [.2 .75 .2];
        end
        fr = variables*coefficients;
        fv = isosurface(xx,yy,zz,reshape(fr,size(xx)), 1);

        if (BITMAP)
            h = patch(fv);
            set(h,'FaceVertexCData',repmat(surfcolor,size(fv.vertices,1),1))
            set(h,'CDataMapping','direct');
            set(h,'FaceColor','interp','EdgeColor','none');
            lighting gouraud
            set(h,'SpecularColorReflectance',.1);
            set(h,'SpecularExponent',.2,'SpecularStrength',.2);
            set(h,'DiffuseStrength',.8,'AmbientStrength',.1);
            
            if (strcmp(filename, 'S041310002'))
                camlight(0,0);
            elseif (strcmp(filename, 'K070309004'))
                camlight(45,0);
                camlight(-45,0);
            elseif (strcmp(filename, 'K082609009'))
                camlight(45,0);
                camlight(-45,0);
            end
        end
    end
    
    % Prettying up the axes for the main plot
    for whichaxes = {'XTick','YTick','ZTick'}
        tmp = get(gca,char(whichaxes)); set(gca,char(whichaxes),[tmp(1) 0 tmp(end)]);
    end
    xticks = get(gca,'XTick');
    set(gca,'FontSize',8);

    % Axes for symbols showing sizes
    axes('position',[1+(whichcell-1)*(LARGEAXISWIDTH+COLMARGIN) 6.5+LARGEAXISWIDTH LARGEAXISWIDTH .2]); hold on;
    set(gca,'YTick',[],'Box','on','TickLength',[0 0]);
    nsymbols = 3;
    resps = linspace(minresp, maxresp,nsymbols);
    for i = 1:nsymbols
        plot(resps(i),0,'ko','MarkerSize',resps(i)*b1+b0,'MarkerFaceColor','black');
    end
    set(gca,'Xlim',[resps(1)-(resps(end)-resps(1))/5 resps(end)+(resps(end)-resps(1))/5],'XTick',round(resps),'XAxisLocation','top');
    %ylabel('Response (sp/sec)','FontSize',14);
    
    if (PLOTSCATTERPLOTS)
        % Getting the predictions
        planepredresp = NT.sum.exptParams.threshold*abs(GTcols*coneweights');
        trueresp = GTstruct.color.colresp(:,1);
        
        % Finding predictions from quadratic model
        [th,ph,r] = cart2sph(GTcols(:,1),GTcols(:,2),GTcols(:,3));
        predr2 = 1./(B(1,1).*(cos(ph).*cos(th)).^2 +...
            B(2,2).*(cos(ph).*sin(th)).^2 +...
            B(3,3).*sin(ph).^2 + ...
            2*B(1,2).*cos(ph).*cos(th).*cos(ph).*sin(th) +...
            2*B(1,3).*cos(ph).*cos(th).*sin(ph) +...
            2*B(2,3).*cos(ph).*sin(th).*sin(ph));
        
        predr2(predr2<0) = Inf;
        quadpredresp = NT.sum.exptParams.threshold.*(r./sqrt(predr2));
        subplotlims = [0 max([trueresp; planepredresp;quadpredresp])*1.05]
        
        axes('position',[1.2+(whichcell-1)*(LARGEAXISWIDTH+COLMARGIN) 4 SMALLAXISWIDTH SMALLAXISWIDTH]); hold on;
        plot(trueresp,planepredresp,'ko','MarkerSize',6,'markerFaceColor','black');
        set(gca,'Ylim',subplotlims,'Xlim',subplotlims)
        if (whichcell == 1)
            ylabel({'Planar','prediction (sp/sec)'},'FontSize',10);
        end
        r = corrcoef([trueresp,planepredresp]);
        title(['r = ',num2str(r(1,2),2)],'FontSize',8);
        
        axes('position',[1.2+(whichcell-1)*(LARGEAXISWIDTH+COLMARGIN) 2 SMALLAXISWIDTH SMALLAXISWIDTH]); hold on;
        plot(trueresp,quadpredresp,'ko','MarkerSize',6,'markerFaceColor','black');
        set(gca,'Ylim',subplotlims,'Xlim',subplotlims)
        if (whichcell == 1)
            ylabel({'Quadratic','prediction (sp/sec)'},'FontSize',10);
        end
        xlabel('Response (sp/sec)','FontSize',10);
        r = corrcoef([trueresp,quadpredresp]);
        title(['r = ',num2str(r(1,2),2)],'FontSize',8);
    end
end

if (BITMAP)
    set(gcf,'InvertHardCopy','off')
    print -dtiff -r600 -cmyk junk
end
%% Section 7) Scatterplot of correlation coefficients
% The results of this anlysis depends slightly on which set of fundamentals
% are used (Both Pearson and Spearman) GDLH 4/27/11 

[fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\Gratings&NT2.txt');
[planefnames] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\plane0.01.txt');
USE10DEG = 1;
CORTYPE = 'Pearson';
data = [];
for cellcounter = 1:size(fnames,1)
    NT = {}; GT = {};
    for i = 1:size(fnames{cellcounter},2)
        filename = findfile(char(fnames{cellcounter}(i)));
        paradigmID = getparadigmID(filename);
        if paradigmID == 150
            GT = nex2stro(filename);
            GTstruct = getGratingTuning(GT,1);
        else
            NT = nex2stro(filename);
        end
    end
    if (USE10DEG)
        % -------------------------------
        % Converting cone contrasts in nex file to 10 deg fundmentals.
        % Must go through excitations first - can't just transform contrasts.
        % -------------------------------
        lmsidxs = [find(strcmp(NT.sum.trialFields(1,:),'lcont'))...
            find(strcmp(NT.sum.trialFields(1,:),'mcont'))...
            find(strcmp(NT.sum.trialFields(1,:),'scont'))];   
        fundamentals = NT.sum.exptParams.fundamentals;
        fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
        mon_spd = NT.sum.exptParams.mon_spd;
        mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
        mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
        M = fundamentals'*mon_spd;
        load ('T_cones_smj10');
        M10 = T_cones_smj10*mon_spd;
        NT.trial(:,lmsidxs) = ConvertConeContrastBasis(M, M10, NT.sum.exptParams.bkgndrgb, NT.trial(:,lmsidxs));
        % --------------------------------
        % Converting GT cone contrasts to 10 degree fundamentals
        % Safe to assume the M matrix and bkgndlms is the same for gratings and NT
        GTstruct.color.colors = ConvertConeContrastBasis(M, M10, NT.sum.exptParams.bkgndrgb, GTstruct.color.colors);
        % ------------------------------- 
    end
    out = NTpreprocess(NT,0,Inf);  % Getting the termination points
    scaled = out(:,[2 3 4]).*repmat(out(:,5),1,3);

    Loog = logical(out(:,end));
    [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);
    A = [quadparams(1) quadparams(4) quadparams(5);...
        quadparams(4) quadparams(2) quadparams(6);...
        quadparams(5) quadparams(6) quadparams(3)];
    coneweights = planeparams'*xformmat';
    B = xformmat*A*xformmat';
    GTcols = GTstruct.color.colors;

    % Getting the predictions
    planepredresp = NT.sum.exptParams.threshold*abs(GTcols*coneweights');
    trueresp = GTstruct.color.colresp(:,1);
    
    % Finding predictions from quadratic model
    [th,ph,r] = cart2sph(GTcols(:,1),GTcols(:,2),GTcols(:,3));
    predr2 = 1./(B(1,1).*(cos(ph).*cos(th)).^2 +...
        B(2,2).*(cos(ph).*sin(th)).^2 +...
        B(3,3).*sin(ph).^2 + ...
        2*B(1,2).*cos(ph).*cos(th).*cos(ph).*sin(th) +...
        2*B(1,3).*cos(ph).*cos(th).*sin(ph) +...
        2*B(2,3).*cos(ph).*sin(th).*sin(ph));
    
    predr2(predr2<0) = Inf;
    quadpredresp = NT.sum.exptParams.threshold.*(r./sqrt(predr2));
    
    r = corr([trueresp,planepredresp,quadpredresp],'type',CORTYPE);
    
    [evecs,evals] = eig(A);
    [evals,i] = sort(diag(evals),1,'ascend');
    WHICHSURF = sum(evals<0);
    shortfn = NT.sum.fileName(find(NT.sum.fileName == filesep,1,'last')+1:find(NT.sum.fileName == '.',1,'last')-1);
    Lplane = ismember(shortfn,[planefnames{:}]');
    data = [data; r(1,2) r(1,3) WHICHSURF Lplane];
end

% Setting up the figure
figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
colordef white;

% % Mariginal histogram (plane)
% axes('position',[2 7.5 3 1]); hold on;
% bins = linspace(-1,1,20);
% binwidth = bins(2)-bins(1);
% [n,x] = hist(data(:,1),bins);
% h = bar(bins,n);
% set(h,'FaceColor','black');
% set(gca,'Xlim',[-1-binwidth/2 1+binwidth/2]);
% plot(nanmedian(data(:,1)),max(n)*1.1,'wv','MarkerFaceColor','black');
% set(gca,'XTick',[-1 -.5 0 .5 1]);
% 
% % Mariginal histogram (quad)
% axes('position',[5.5 4 1 3]); hold on;
% [n,x] = hist(data(:,2),bins);
% h = barh(bins,n);
% set(h,'FaceColor','black');
% set(gca,'Ylim',[-1-binwidth/2 1+binwidth/2]);
% plot(max(n)*1.1,nanmedian(data(:,2)),'w<','MarkerFaceColor','black');
% set(gca,'YTick',[-1 -.5 0 .5 1]);

% Plotting gratings data
axes('position',[2 2 3 3]); hold on;
plot(data(:,1),data(:,2),'ko','MarkerFaceColor','black','MarkerSize',4);
plot([-1 1],[-1 1],'k.-');

Lplane = logical(data(:,4));
Lellipse = logical(data(:,3) == 0);
Lhyper1 = logical(data(:,3) == 1);
Lhyper2 = logical(data(:,3) == 2);
plot(data(Lellipse & ~Lplane,1),data(Lellipse & ~Lplane,2),'ko','MarkerFaceColor',[.7 .7 0],'MarkerSize',4,'MarkerEdgeColor',[.7 .7 0]);
%plot(data(Lellipse&~Lplane,1),data(Lellipse&~Lplane,2),'ko','MarkerFaceColor','black','MarkerSize',4);
%plot(data(Lhyper1&~Lplane,1),data(Lhyper1&~Lplane,2),'ko','MarkerFaceColor','magenta','MarkerSize',4);
%plot(data(Lhyper2&~Lplane,1),data(Lhyper2&~Lplane,2),'ko','MarkerFaceColor','green','MarkerSize',4);

set(gca,'XTick',[-1 -.5 0 .5 1],'YTick',[-1 -.5 0 .5 1]);
set(gca,'Xlim',[-1 1],'Ylim',[-1 1]);
%set(gca,'Xlim',[-1-binwidth/2 1+binwidth/2],'Ylim',[-1-binwidth/2 1+binwidth/2]);
xlabel('Correlation with plane predictions','FontSize',12);
ylabel('Correlation with quadratic predictions','FontSize',12);

% Marginal histogram (difference)
% Have to rotate this thing manually in Illustrator I think
difference = data(:,1)-data(:,2);
lim = max(abs(difference));
[n,x] = hist(difference,linspace(-lim,lim,15));
axespos = get(gca,'Position');
axes('position',[2 6 axespos(3)*lim/sqrt(2) 1]); hold on;
h = bar(x,n);
set(h,'FaceColor','black');
plot(nanmedian(difference),max(n)*.9,'kv','MarkerFaceColor','white')
ylabel('Count','FontSize',12);

% A little analysis
signrank(difference)
%[h,p] = ttest(difference(Lplane))
%[h,p] = ttest(difference(Lellipse&~Lplane))
%[h,p] = ttest(difference(Lhyper1&~Lplane))
%[h,p] = ttest(difference(Lhyper2&~Lplane))
signrank(difference(Lellipse))
signrank(difference(Lhyper1))
signrank(difference(Lhyper2))
anova1(difference(~Lplane),data(~Lplane,3))
anova1(data(~Lplane,2),data(~Lplane,3))
% Huge differences in the quality of the quadratic predictions
% depending on the shape of the quadratic surface - Ellipsoids are the
% worst.


%% Section 8) Pretty plot of an example data set (grating responses +
% isoresponse surfaces)
filename = 'K070309004';
INCLUDEGRATINGDATA = 1;
ALIGNTOPCS = 1;
%examplefilenames = {'S041310002','K070309004','K082609010'};

NT = nex2stro(findfile(filename));
out = NTpreprocess(NT,0,Inf);  % Getting the termination points
scaled = out(:,[2 3 4]).*repmat(out(:,5),1,3);

Loog = logical(out(:,end));
[planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);

fundamentals = NT.sum.exptParams.fundamentals;
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = NT.sum.exptParams.mon_spd;
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;

% Plotting quadratic fits
A =  [quadparams(1) quadparams(4) quadparams(5);...
    quadparams(4) quadparams(2) quadparams(6);...
    quadparams(5) quadparams(6) quadparams(3)];
if (ALIGNTOPCS)
    % Transforming quadparams into un-whitened space
    B = xformmat*A*xformmat';
    plotlims  = max(abs(scaled(~Loog,:)))*1.1;
else
    % or tranforming quadparams into principal axes space
    B = A;
    [v,d] = eig(B);
    d = [-2 0 0; 0 1 0; 0 0 7];
    B = (v*sqrt(d))'*B*(v*sqrt(d));
    plotlims = [2.1 4.1 5.1];
end

npts = 150;
[x,y,z] = meshgrid(linspace(-plotlims(1),plotlims(1),npts),linspace(-plotlims(2),plotlims(2),npts),linspace(-plotlims(3),plotlims(3),npts));
variables = [x(:).^2 y(:).^2 z(:).^2 2*x(:).*y(:) 2*x(:).*z(:) 2*y(:).*z(:)];
v = variables*[B(1,1) B(2,2) B(3,3) B(1,2) B(1,3) B(2,3)]';
fv = isosurface(x,y,z,reshape(v,size(x)),1);
% trying to color the surface
% fv.vertices is in cone contrasts
nvert = size(fv.vertices,1);
bkgndrgb = [.5 .5 .5];
bkgndlms = M*bkgndrgb';
lms = fv.vertices+repmat(bkgndlms',nvert,1);

rgbs = (M*lms')';
multfactors = 3./max(abs(rgbs));

fv.facevertexcdata = rgbs.*repmat(multfactors,nvert,1)+.5;

% Doing the plotting
figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
axes('position',[0 0 8.5 11]); set(gca,'XTick',[],'YTick',[]); hold on;
h = patch(fv);
set(h,'FaceColor','interp')
set(h,'EdgeColor','none');
lighting gouraud

set(h,'SpecularColorReflectance',.4);
set(h,'SpecularExponent',10,'SpecularStrength',.1);
set(h,'DiffuseStrength',.3,'AmbientStrength',.3);
set(gca,'View',[139 2]);
%camlight(0,0);
camlight(30,40);
axis vis3d;
set(gcf,'Color',[0 0 0])
set(gca,'Color',[0 0 0],'Visible','off');


if (INCLUDEGRATINGDATA)
    % Getting grating info
    if (strcmp(filename,'K070309004'));
        GT = nex2stro(findfile('K070309003.2'));
        GTstruct = getGratingTuning(GT);
        GTcols = GTstruct.color.colors;
        maxresp = max(GTstruct.color.colresp(:,1));
        minresp = min(GTstruct.color.colresp(:,1));
        maxsize = .05;
        minsize = .005;
        b1 = (maxsize-minsize)./(maxresp-minresp);
        b0 = maxsize-(maxresp*b1);
        [x,y,z] = sphere(150);
        axislim = [get(gca,'XLim') get(gca,'YLim') get(gca,'ZLim')];
        axisnormfactors = [axislim(1),axislim(3),axislim(5)]./axislim(1);      
        
        colormap('copper');
        for i = 1:size(GTstruct.color.colors,1)
            a = (GTstruct.color.colresp(i,1)*b1+b0).*axisnormfactors;
            h(1) = surf(a(1)*x+GTcols(i,1),a(2)*y+GTcols(i,2),a(3)*z+GTcols(i,3));
            h(2) = surf(a(1)*x-GTcols(i,1),a(2)*y-GTcols(i,2),a(3)*z-GTcols(i,3));
            set(h,'EdgeColor','none','FaceColor','flat');
            set(h,'Cdata',30*ones(size(x)),'CDataMapping','direct');
        end
    end
end
set(gcf,'InvertHardCopy','off')
print -dtiff -r300 -cmyk junk

%% Section 9) Distribution of different types of cells and their principal axes
SKIPINITCOLDIRS = 0;
[fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\Gratings&NT2.txt');
[planefnames] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\plane0.01.txt');
AXESWIDTH = 1.7;
Lplane = []; eigenvectors = []; eigenvalues = [];
for cellcounter = 1:size(fnames,1)
    NT = {};
    for i = 1:size(fnames{cellcounter},2)
        filename = findfile(char(fnames{cellcounter}(i)));
        paradigmID = getparadigmID(filename);
        if paradigmID == 103
            NT = nex2stro(filename);
        end
    end
    lmsidxs = [find(strcmp(NT.sum.trialFields(1,:),'lcont'))...
        find(strcmp(NT.sum.trialFields(1,:),'mcont'))...
        find(strcmp(NT.sum.trialFields(1,:),'scont'))];
    fundamentals = NT.sum.exptParams.fundamentals;
    fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
    mon_spd = NT.sum.exptParams.mon_spd;
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
    M = fundamentals'*mon_spd;
    load ('T_cones_smj10');
    M10 = T_cones_smj10*mon_spd;
    NT.trial(:,lmsidxs) = ConvertConeContrastBasis(M, M10, NT.sum.exptParams.bkgndrgb, NT.trial(:,lmsidxs));
    
    if(SKIPINITCOLDIRS)
        L = NT.trial(:,strcmp(NT.sum.trialFields(1,:),'levelidx')) < 2;
        NT.trial(L,:) = [];
        NT.ras(L,:) = [];
    end
    
    out = NTpreprocess(NT,0,Inf);  % Getting the termination points
    scaled = out(:,[2 3 4]).*repmat(out(:,5),1,3);
    Loog = logical(out(:,end));
 
    [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);
    A = [quadparams(1) quadparams(4) quadparams(5);...
         quadparams(4) quadparams(2) quadparams(6);...
         quadparams(5) quadparams(6) quadparams(3)];
    % Transforming quadparams into un-whitened space
    [evecs,evals] = eig(xformmat*A*xformmat')
    [evals,i] = sort(diag(evals),1,'ascend');
    evecs = evecs(:,i);
    
    eigenvectors = cat(3,eigenvectors,evecs);
    eigenvalues = cat(2,eigenvalues,evals);
    shortfn = NT.sum.fileName(find(NT.sum.fileName == filesep,1,'last')+1:find(NT.sum.fileName == '.',1,'last')-1);
    Lplane = cat(2,Lplane, ismember(shortfn,[planefnames{:}]'))
end
surfacetypes = sum(eigenvalues<0);

% Setting up the figure
figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);

% Piechart
axes('position',[1.5 5 AXESWIDTH AXESWIDTH]);
%piedata = [sum(Lplane) sum(~Lplane&surfacetypes == 0) sum(~Lplane&surfacetypes == 1) sum(~Lplane&surfacetypes == 2)];
colormap([.75 .25 .25; 1 .5 .5; .25 .75 .25; .5 1 .5; .25 .25 .75; .5 .5 1]);
piedata = [sum(Lplane&surfacetypes == 0), sum(~Lplane&surfacetypes == 0)
    sum(Lplane&surfacetypes == 1), sum(~Lplane&surfacetypes == 1)
    sum(Lplane&surfacetypes == 2), sum(~Lplane&surfacetypes == 2)]; 

%piedata = [sum(surfacetypes == 0) sum(surfacetypes == 1) sum(surfacetypes == 2)];
pie(piedata);
%h = legend({'Plane','Ellipsoid','Hyperboloid 1','Hyperboloid 2'});
h = legend({'Ellipsoid','Hyperboloid 1','Hyperboloid 2'});
pos = get(h,'position');
set(h, 'position',[pos(1)-1.3 pos(2)+.7 pos(3)*.9 pos(4)*.9])
set(h, 'fontsize',8);

% % Barchart showing the relationship between surfacetype and planarity
% types = zeros(2,3);
% for i = 0:2
%     types(1,i+1) = sum(surfacetypes == i & Lplane);
%     types(2,i+1) = sum(surfacetypes == i & ~Lplane);
% end
% axes('position',[1.5 3 AXESWIDTH AXESWIDTH/sqrt(2)]);
% bar(types')
% set(gca,'Xticklabel',{'Ellipsoid','Hyper 1','Hyper 2'});
% set(gca,'FontSize',7);

% Principal axes for ellipsoidal cells
evalratiothresh = 5;
for whichsurf = 0:2; % 0 = ellipsoid, 1 = 1 sheet, 2 = 2 sheets
    out = [];
    for i = 1:size(eigenvalues,2)
        if (surfacetypes(i) == whichsurf)
            out = [out; reshape(eigenvectors(:,:,i),1,9)];
            ratio12 = max(abs(sqrt(eigenvalues(1,i))./sqrt(eigenvalues(2,i))),sqrt(abs(eigenvalues(2,i))./sqrt(eigenvalues(1,i))));
            ratio23 = max(abs(sqrt(eigenvalues(2,i))./sqrt(eigenvalues(3,i))),sqrt(abs(eigenvalues(3,i))./sqrt(eigenvalues(2,i))));
            if (ratio12 < evalratiothresh)
                out(end,[1:6]) = nan;
            end
            if (ratio23 < evalratiothresh)
                out(end,[4:9]) = nan;
            end
        end
    end
    axes('position',[4 1+whichsurf*2 AXESWIDTH AXESWIDTH/sqrt(2)]); hold on;
    symbols = {'gs','rd','bv'};
    titles = {'Ellipsoid','Hyperboloid 1 sheet','Hyperboloid 2 sheets'};
    for i = 1:3
        idxs = [1 2 3]+(i-1)*3;
        tmp = out(:,idxs)./repmat(sum(abs(out(:,idxs)),2),1,3);
        tmp = tmp.*repmat(sign(tmp(:,2)),1,3);
        plot(tmp(:,1),tmp(:,2),symbols{i},'MarkerFaceColor',symbols{i}(1),'MarkerSize',3);
    end
    plot([-1 0 1 -1],[0 1 0 0],'k-');
   % title(titles{whichsurf+1});
   xlabel('L-cone weight');
   ylabel('M-cone weight');
   set(gca,'XTick',[-1 0 1],'YTick',[0 1]);
%     if whichsurf == 0
%         legend('Concave (open)','Medium convex','Tight convex');
%     elseif whichsurf == 1
%         legend('Tight concave','Medium concave','Convex');
%     else % whichsurf == 3
%         legend('Weak convex','Medium convex','Tight convex');
%     end
    % How many ellipsoids have axes that we're chucking?
       if (whichsurf == 0)
            sum(any(isnan(out')))
       end
end
h=legend('PC1','PC2','PC3');
set(h,'Position',[6 3.3 1 .6]);

% How many axis comparisons, and of these, how many fall below
% evalratiothresh?
naxiscomparisons = size(eigenvalues,2)*2;
counter = 0;
for i = 1:size(eigenvalues,2)
    ratio12 = max(abs(sqrt(eigenvalues(1,i))./sqrt(eigenvalues(2,i))),sqrt(abs(eigenvalues(2,i))./sqrt(eigenvalues(1,i))));
    ratio23 = max(abs(sqrt(eigenvalues(2,i))./sqrt(eigenvalues(3,i))),sqrt(abs(eigenvalues(3,i))./sqrt(eigenvalues(2,i))));
    if (ratio12 < evalratiothresh)
        counter = counter +1;
    end
    if (ratio23 < evalratiothresh)
        counter = counter +1;
    end
end
counter
%% Section 9.1) Distribution of different types of cells and their principal
% axes.  Including randomization test and data skipping first round.
[fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\Gratings&NT2.txt');
markersize = 3;
%[fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\altinitcoldirs.txt');
%markersize = 3;

AXESWIDTH = 1.7;
evecs1 = []; evecs2 = []; evals = [];
for cellcounter = 1:size(fnames,1)
    NT = {};
    for i = 1:size(fnames{cellcounter},2)
        filename = findfile(char(fnames{cellcounter}(i)));
        paradigmID = getparadigmID(filename);
        if paradigmID == 103
            if (isempty (NT))
                NT = nex2stro(filename)
            else
                tmp = nex2stro(filename);
                NT = strocat(NT, tmp);
            end
        end
    end
    lmsidxs = [find(strcmp(NT.sum.trialFields(1,:),'lcont'))...
        find(strcmp(NT.sum.trialFields(1,:),'mcont'))...
        find(strcmp(NT.sum.trialFields(1,:),'scont'))];
    fundamentals = NT.sum.exptParams.fundamentals;
    fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
    mon_spd = NT.sum.exptParams.mon_spd;
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
    M = fundamentals'*mon_spd;
    load ('T_cones_smj10');
    M10 = T_cones_smj10*mon_spd;
    NT.trial(:,lmsidxs) = ConvertConeContrastBasis(M, M10, NT.sum.exptParams.bkgndrgb, NT.trial(:,lmsidxs));
    
    % CAH asked - how does this change if you equate contrasts for
    % detection threshold?
   % NT.trial(:,find(strcmp(NT.sum.trialFields(1,:),'scont'))) = NT.trial(:,find(strcmp(NT.sum.trialFields(1,:),'scont')))./10;
    % NOTE: Comment out the line above in general!  This is just a test.
    
    for j = 1:2  % 1 = normal, 2 = skipping initial color directions.
        out = NTpreprocess(NT,0,Inf);  % Getting the termination points
        if (i == 2)
            firstroundidx = logical(out(:,6) == 1);
            out(firstroundidx,:) = [];
        end
        scaled = out(:,[2 3 4]).*repmat(out(:,5),1,3);
        Loog = logical(out(:,end));
        
        [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);
        A = [quadparams(1) quadparams(4) quadparams(5);...
            quadparams(4) quadparams(2) quadparams(6);...
            quadparams(5) quadparams(6) quadparams(3)];
        % Transforming quadparams into un-whitened space
        [tmpevecs,tmpevals] = eig(xformmat*A*xformmat');
        [tmpevals,i] = sort(diag(tmpevals),1,'ascend');
        tmpevecs = tmpevecs(:,i);
        if (j == 1)
            evals = cat(2,evals,tmpevals);
            evecs1 = cat(3,evecs1,tmpevecs);
        else
            evecs2 = cat(3,evecs2,tmpevecs);           
        end
    end
    shortfn = NT.sum.fileName(find(NT.sum.fileName == filesep,1,'last')+1:find(NT.sum.fileName == '.',1,'last')-1);
end
surfacetypes = sum(evals<0);

% Setting up the figure
figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
% Nested loops: Which surface, which analysis, which cell
symbols = {'gs','rd','bv'};
titles = {'Ellipsoid','Hyperboloid 1 sheet','Hyperboloid 2 sheets'};
evalratiothresh = 5;

for whichsurf = 0:2; % 0 = ellipsoid, 1 = 1 sheet, 2 = 2 sheets
    for j = [1 2 3]   % Which analysis
        if (j == 1)
            axes('position',[1 1+whichsurf*2 AXESWIDTH AXESWIDTH/sqrt(2)]); hold on;
        elseif (j == 2)
            axes('position',[3.5 1+whichsurf*2 AXESWIDTH AXESWIDTH/sqrt(2)]); hold on;
        else % j == 3
            axes('position',[6 1+whichsurf*2 AXESWIDTH AXESWIDTH/sqrt(2)]); hold on;
        end
        out = [];
        for i = 1:size(evals,2)  % looping over cells.  Creating "out"
            if (any(strcmp(fnames{i}(2), examplefilenames([1 2 3]))))
                highlight = 1;
            else
                highlight = 0;
            end
            if (surfacetypes(i) == whichsurf)
               if (j == 2)
                    switch (whichsurf)
                        case 1
                            whichvecstorotate = [1 2];
                        case 2
                            whichvecstorotate = [1 3];
                        case 0
                            whichvecstorotate = [2 3];
                    end  
                    th = unifrnd(0,2*pi);
                    rotmat = [cos(th) sin(th); -sin(th) cos(th)];
                    tmp = evecs1(:,:,i);
                    tmp(:,whichvecstorotate) = tmp(:,whichvecstorotate)*rotmat;
                    out = [out; reshape(tmp,1,9) highlight];
               end
               if (j == 1)
                   out = [out; reshape(evecs1(:,:,i),1,9) highlight];
               elseif (j == 3)
                   out = [out; reshape(evecs2(:,:,i),1,9) highlight];
               end
               ratio12 = max(abs(sqrt(evals(1,i))./sqrt(evals(2,i))),sqrt(abs(evals(2,i))./sqrt(evals(1,i))));
               ratio23 = max(abs(sqrt(evals(2,i))./sqrt(evals(3,i))),sqrt(abs(evals(3,i))./sqrt(evals(2,i))));
             
               if (ratio12 < evalratiothresh)
                   out(end,[1:6]) = nan;
               end
               if (ratio23 < evalratiothresh)
                   out(end,[4:9]) = nan;
               end
            end
        end
        
        for i = 1:3
            idxs = [1 2 3]+(i-1)*3;
            tmp = out(:,idxs)./repmat(sum(abs(out(:,idxs)),2),1,3);
            tmp = tmp.*repmat(sign(tmp(:,2)),1,3);
            plot(tmp(:,1),tmp(:,2),symbols{i},'MarkerFaceColor',symbols{i}(1),'MarkerSize',markersize);
            L = out(:,10) == 1;
            if (any(L))
               text(tmp(L,1),tmp(L,2),'1');
               plot(tmp(L,1),tmp(L,2),symbols{i},'MarkerEdgeColor','black','MarkerFaceColor',symbols{i}(1),'MarkerSize',markersize);
            end
        end
        plot([-1 0 1 -1],[0 1 0 0],'k-');
        title(titles{whichsurf+1});
        xlabel('L-cone weight');
        ylabel('M-cone weight');
        set(gca,'XTick',[-1 0 1],'YTick',[0 1]);
    end
end
h=legend('PC1','PC2','PC3');
set(h,'Position',[6 3.3 1 .6]);


%% Section 10) Effects of changes in threshold

NCONEWEIGHTS = 3;
NQUADPARAMS = 6;
load ('T_cones_smj10');

% A few examples
filenames = {'K033111002','S031610003','S081010005','S112410005',...
                'S050410006','S041510003','S021810008','S012810004'};

[fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\NTmultiples.txt');
maxnfiles = 0;
for i = 1:size(fnames,1)
    maxnfiles = max(maxnfiles, length(fnames{i}));
end
whichsurface = nan*ones(length(fnames), maxnfiles);
thresholds = nan*ones(length(fnames), maxnfiles);
allconeweights = nan*ones(length(fnames), maxnfiles*NCONEWEIGHTS);
allquadparams = nan*ones(length(fnames), maxnfiles*NQUADPARAMS);
allscaled = cell(size(fnames,1),maxnfiles);
allplaneSSEs = nan*ones(length(fnames),maxnfiles);
allquadSSEs = nan*ones(length(fnames),maxnfiles);
sumnoog = nan*ones(length(fnames),maxnfiles);

for cellcounter = 1:size(fnames,1)
    NT = {}; 
    for i = 1:size(fnames{cellcounter},2)
        filename = findfile(char(fnames{cellcounter}(i)));
        NT{i} = nex2stro(filename);
        lmsidxs = [find(strcmp(NT{i}.sum.trialFields(1,:),'lcont'))...
            find(strcmp(NT{i}.sum.trialFields(1,:),'mcont'))...
            find(strcmp(NT{i}.sum.trialFields(1,:),'scont'))];
        % Converting to 10 degree fundamentals
        fundamentals = NT{i}.sum.exptParams.fundamentals;
        fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
        mon_spd = NT{i}.sum.exptParams.mon_spd;
        mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
        mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
        M = fundamentals'*mon_spd;
        M10 = T_cones_smj10*mon_spd;
        NT{i}.trial(:,lmsidxs) = ConvertConeContrastBasis(M, M10, NT{i}.sum.exptParams.bkgndrgb, NT{i}.trial(:,lmsidxs));
    end
    
    scaled = []; Loog = [];
    for i = 1:length(NT)
        out = NTpreprocess(NT{i},0,Inf);  % .4, 2  or 0, Inf
        scaled = [scaled; out(:,[2:4]).*repmat(out(:,5), 1,3) logical(out(:,7))];
        Loog = [Loog; logical(out(:,7))];
        [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(out(:,[2:4]).*repmat(out(:,5), 1,3), out(:,7));
        A = [quadparams(1) quadparams(4) quadparams(5);...
            quadparams(4) quadparams(2) quadparams(6);...
            quadparams(5) quadparams(6) quadparams(3)];
        B = xformmat*A*xformmat';
        coneweights = planeparams'*xformmat';
        coneweights = coneweights./repmat(sum(abs(coneweights),2),1,3);
        coneweights = coneweights.*repmat(sign(coneweights(1,2)),1,3);
        allconeweights(cellcounter,(i-1)*NCONEWEIGHTS+[1:NCONEWEIGHTS]) = coneweights;
        allquadparams(cellcounter,(i-1)*NQUADPARAMS+[1:NQUADPARAMS]) = [B(1,1) B(2,2) B(3,3) B(1,2) B(1,3) B(2,3)];
        allplaneSSEs(cellcounter,i) = planeSSE;
        allquadSSEs(cellcounter,i) = quadSSE;
        sumnoog(cellcounter,i) = sum(~out(:,7));
        if (sum(~out(:,7)) < 3)
            keyboard
        end
        whichsurface(cellcounter,i) = sum(eig(B)<0);
        thresholds(cellcounter,i) = NT{i}.sum.exptParams.threshold;
        allscaled{cellcounter,i} = scaled;
    end
end

% Setting up the figure (for the EPS part)
figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);

% Cone weights
axes('Position',[1,5,2*sqrt(2),2]); hold on;
% % Old plot - Lennie style cone weights space
% for i = 1:size(allconeweights,1)
%     tmp = reshape(allconeweights(i,:),3,maxnfiles)';
%     plot(tmp(:,1),tmp(:,2),'k-');
%     for j = 1:size(allconeweights,2)/NCONEWEIGHTS
%         h = plot(allconeweights(i,(j-1)*NCONEWEIGHTS+1),allconeweights(i,(j-1)*NCONEWEIGHTS+2),'ko','MarkerFaceColor',[0 0 0],'MarkerSize',4);
%         if (allconeweights(i,3) < 0)
%             set(h,'MarkerFaceColor',[1 1 1]);
%         end
%     end
% end
% plot([-1 1 0 -1],[0 0 1 0],'k-');
% set(gca,'XTick',[-1 -.5 0 .5 1],'YTick',[-1 -.5 0 .5 1]);
% xlabel('L-cone weight','FontSize',14);
% ylabel('M-cone weight','FontSize',14);

% % New plot: histograms of correlation coefficients
r = zeros(size(allconeweights,1),1);
for i = 1:size(allconeweights,1)
    tmp = corrcoef([allconeweights(i,[1 2 3])', allconeweights(i,[4 5 6])']);
    r(i) = tmp(1,2);
end
[n,x] = hist(r,[-1:.1:11]);
bar(x,n,'FaceColor','black')
set(gca,'XTick',[-1 -.5 0 .5 1],'xlim',[-1.1,1.1]);
xlabel('Cone weight correlation','FontSize',14);
ylabel('Count','FontSize',14);

% Find cells that switch from L+M to L-M
v = allconeweights(:,[1 4]);
L = sign(allconeweights(:,1)) ~= sign(allconeweights(:,4));
sum(L)
pfs = [];
for i = find(L)'
    for j = 1:2
        F1 = ((allplaneSSEs(i,j)-allquadSSEs(i,j))/3)/(allquadSSEs(i,j)/(sumnoog(i,j)-6));
        pfs(i == find(L),j) = 1-fcdf(F1,3,sumnoog(i,j)-3);
    end
end

% Flagging the cells whose cone weights change a lot
letter = ['C', 'D', 'E', 'F'];
for i = 1:4 % hard coding four example cells
    for j = 1:size(fnames,1)
        if any(strcmp(filenames(i),fnames{j}))
            [i j]
            for k = 1:size(allconeweights,2)/NCONEWEIGHTS
                h = text(allconeweights(j,(k-1)*NCONEWEIGHTS+1),allconeweights(j,(k-1)*NCONEWEIGHTS+2),letter(i));
                set(h,'HorizontalAlignment','center','FontSize',12,'FontWeight','bold')
            end
        end
    end
end
% Transition matrix    
x = zeros(3,3);
n = zeros(3,1);
for i = 1:size(thresholds,1)
    [a,b] = sort(thresholds(i,:),2,'ascend');
    for j = 1:sum(~isnan(thresholds(i,:)))
        shapes = whichsurface(i,b(~isnan(a)));
        n(shapes(j)+1) = n(shapes(j)+1)+1;
        if (j<sum(~isnan(thresholds(i,:))))
            x(shapes(j)+1,shapes(j+1)+1) = x(shapes(j)+1,shapes(j+1)+1)+1;
        end
    end
end
sum(x(:)) % total number of cells
sum(diag(x)) % number of cells that didn't change shape

axes('Position',[1,1,2*sqrt(2),2*sqrt(2)]); hold on;
transprob = x./repmat(sum(x,2),1,3)  % transition probabilities
%transprob = x./max(x(:)) % transitions - no longer normalized to be probabilities

locs = [0 -1; -1 1; 1 1];  % Order: Ellipsoid, Hyperboloid 1, Hyperboloid 2
text(-1,1.5,'Hyperboloid 1','HorizontalAlignment','center');
text(1,1.5,'Hyperboloid 2','HorizontalAlignment','center');
text(-.5,-0.8,'Ellipsoid','HorizontalAlignment','center');
set(gca,'Xlim',[min(locs(:,1)) max(locs(:,1))]*2,'Ylim',[min(locs(:,2)) max(locs(:,2))]*2)
lineweightfact = 8;
for i = 1:3
    for j = 1:3
        if (transprob(i,j) > 0)
            if (i == j)
                scalefactor = .2;
                tmp = linspace(0,2*pi,50);
                plot(scalefactor*cos(tmp)+locs(i,1)*1.1, scalefactor*sin(tmp)+locs(i,2)*1.1,'k-','LineWidth',lineweightfact*transprob(i,j));
            else
                if (i > j)
                    epsilon = .05;
                else
                    epsilon = -.05;
                end
                plot([locs(i,1) locs(j,1)],[locs(i,2) locs(j,2)]+epsilon,'k-','LineWidth',lineweightfact*transprob(i,j));
            end
        end
    end
end

for i = 1:3
    plot(locs(i,1),locs(i,2),'ko','MarkerSize',22,'MarkerFaceColor','white','Linewidth',2)
end
set(gca,'XTick',[],'YTick',[]);

% Chisquare test on transition probabilities
y = [1 0 0; 1 3 0; 0 2 1]

% Setting up the figure (for the bitmapped part)
figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);

views = [-75 32; -48 24; 33 -20; 42 24; 53 18; 49 22; 53 18; 64 32];
plotcolors = [1 0 1; 0 0 1; 0 1 0];
for whichcell = 1:length(filenames)
    tmp = nan*ones(size(fnames,1),1);
    for i = 1:size(fnames,1)
        tmp(i) = strcmp(fnames{i}(1),filenames(whichcell));
    end
    neuronidx = find(tmp);
    AXISCOL = 1+3*round((whichcell-1)/length(filenames));
    AXISROW = 2+2*mod((whichcell-1),length(filenames)/2);

    axes('Position',[AXISCOL,AXISROW,1.3,1.3]); hold on;
    %subplot(4,2,whichcell); hold on;
    tmpthresh = thresholds(neuronidx,~isnan(thresholds(neuronidx,:)));
    [y,idx] = sort(tmpthresh);

    
    for i = 1:length(idx)
        if (~isempty(allscaled{neuronidx,idx(i)}))
            scaled = allscaled{neuronidx,idx(i)}(:,1:3);
            Loog = allscaled{neuronidx,idx(i)}(:,4);
            plotlim  = max(abs(scaled(~Loog,:)))*1.4;
            [xx yy zz] = meshgrid(linspace(-plotlim(1),plotlim(1),40),...
                linspace(-plotlim(2),plotlim(2),40),...
                linspace(-plotlim(3),plotlim(3),40));
            xformedxyz = [xx(:) yy(:) zz(:)];
            variables = [xformedxyz(:,1).^2 xformedxyz(:,2).^2 xformedxyz(:,3).^2 2*xformedxyz(:,1).*xformedxyz(:,2) 2*xformedxyz(:,1).*xformedxyz(:,3) 2*xformedxyz(:,2).*xformedxyz(:,3)];
            fr = variables*allquadparams(neuronidx,NQUADPARAMS*(i-1)+[1:NQUADPARAMS])';
            surfstruct = isosurface(xx,yy,zz,reshape(fr,size(xx)), 1);
            p = patch(surfstruct);
            col = plotcolors(find(unique(y) == tmpthresh(i)),:);
            plot3(scaled(~Loog,1),scaled(~Loog,2),scaled(~Loog,3),'o','MarkerFaceColor',col/2,'MarkerEdgeColor','none','MarkerSize',3);
            plot3(-scaled(~Loog,1),-scaled(~Loog,2),-scaled(~Loog,3),'o','MarkerFaceColor',col/2,'MarkerEdgeColor','none','MarkerSize',3);
            set(p,'FaceAlpha',1-idx(i)/(length(idx)+1),'FaceColor',col,'Edgealpha',0);
        end
    end
    axis vis3d;
    set(gca,'XTick',[],'YTick',[],'ZTick',[],'View',views(whichcell,:));
    set(gca,'XLim',[-plotlim(1) plotlim(1)],'YLim',[-plotlim(2) plotlim(2)],'ZLim',[-plotlim(3) plotlim(3)]);
    camlight;
    lighting phong;
    axis square;
end
%%
% Section 10.1 Effects of changing threshold version 2
% This cell takes a while to run (because of the constrained model fit).

NCONEWEIGHTS = 3;
NQUADPARAMS = 6;
load ('T_cones_smj10');
[fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\NTmultiples.txt');
[planefnames] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\plane0.01.txt');
maxnfiles = 0;
for i = 1:size(fnames,1)
    maxnfiles = max(maxnfiles, length(fnames{i}));
end
whichsurface = nan*ones(length(fnames), maxnfiles);
thresholds = nan*ones(length(fnames), maxnfiles);
allconeweights = nan*ones(length(fnames), maxnfiles*NCONEWEIGHTS);
allquadparams = nan*ones(length(fnames), maxnfiles*NQUADPARAMS);
allquadSSEs = Inf*ones(length(fnames), maxnfiles);
constrquadSSE = Inf*ones(length(fnames),1);
constquadparams = nan*ones(length(fnames), maxnfiles*NQUADPARAMS);
allscaled = cell(size(fnames,1),maxnfiles);
Lplane = zeros(size(fnames,1),1);
for cellcounter = 1:size(fnames,1)
    cellcounter
    NT = {}; 
    for i = 1:size(fnames{cellcounter},2)
        filename = findfile(char(fnames{cellcounter}(i)));
        NT{i} = nex2stro(filename);
        lmsidxs = [find(strcmp(NT{i}.sum.trialFields(1,:),'lcont'))...
            find(strcmp(NT{i}.sum.trialFields(1,:),'mcont'))...
            find(strcmp(NT{i}.sum.trialFields(1,:),'scont'))];
        % Converting to 10 degree fundamentals
        fundamentals = NT{i}.sum.exptParams.fundamentals;
        fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
        mon_spd = NT{i}.sum.exptParams.mon_spd;
        mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
        mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
        M = fundamentals'*mon_spd;
        M10 = T_cones_smj10*mon_spd;
        NT{i}.trial(:,lmsidxs) = ConvertConeContrastBasis(M, M10, NT{i}.sum.exptParams.bkgndrgb, NT{i}.trial(:,lmsidxs));
    end
    
    % Doing the main analysis (one file at a time)
    for i = 1:length(NT)
        out = NTpreprocess(NT{i},0,Inf);  % .4, 2  or 0, Inf
        [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(out(:,[2:4]).*repmat(out(:,5), 1,3), out(:,7));
        A = [quadparams(1) quadparams(4) quadparams(5);...
            quadparams(4) quadparams(2) quadparams(6);...
            quadparams(5) quadparams(6) quadparams(3)];
        B = xformmat*A*xformmat';
        coneweights = planeparams'*xformmat';
        coneweights = coneweights./repmat(sum(abs(coneweights),2),1,3);
        coneweights = coneweights.*repmat(sign(coneweights(1,2)),1,3);
        allconeweights(cellcounter,(i-1)*NCONEWEIGHTS+[1:NCONEWEIGHTS]) = coneweights;
        allquadparams(cellcounter,(i-1)*NQUADPARAMS+[1:NQUADPARAMS]) = [B(1,1) B(2,2) B(3,3) B(1,2) B(1,3) B(2,3)];
        allquadSSEs(cellcounter,i) = quadSSE;
        thresholds(cellcounter,i) = NT{i}.sum.exptParams.threshold;
        whichsurface(cellcounter,i) = sum(eig(B)<0);
        allscaled{cellcounter,i} = [out(:,[2:4]).*repmat(out(:,5), 1,3) logical(out(:,7))];
    end
    
    % Doing the constrained analysis (both files at once)
    scaled = [allscaled{cellcounter,1};allscaled{cellcounter,2}];
    Loog = scaled(:,4);
    scaled(:,4) = [ones(size(allscaled{cellcounter,1},1),1); 2*ones(size(allscaled{cellcounter,2},1),1)];
    [v,d] = eig(cov([scaled(~Loog,[1 2 3]); -scaled(~Loog,[1 2 3])]));
    d = diag(d);
    if (min(d) < 2*eps)
        disp('Too few data points for whitening');
        whtmat = eye(3);
    else
        whtmat = v*diag(sqrt(1./d));
    end
    newscaled = [scaled(:,[1 2 3])*whtmat scaled(:,4)];
    [th,ph,r] = cart2sph(newscaled(:,1),newscaled(:,2),newscaled(:,3));
    
    errs = zeros(50,50);
    tmp = linspace(0,pi,size(errs,1));
    for i = 1:length(tmp)
        for j = 1:length(tmp)
            tmpresid = 0;
            for k = 1:2
                L = scaled(:,end) == k;
                [tmpa, tmpb, tmpc] = sph2cart(tmp(i),tmp(j),1);
                predr = 1./(tmpa.*cos(ph).*cos(th)+tmpb.*cos(ph).*sin(th)+tmpc.*sin(ph));
                predr = predr.*geomean(r(L&~Loog)./abs(predr(L&~Loog)));
                resid = log(abs(predr(L)))-log(r(L));
                resid(Loog(L)&(abs(predr(L)) > r(L))) = 0;
                tmpresid = tmpresid + sum(resid.^2);
            end
            errs(i,j) = tmpresid;
        end
    end
    [cand_i, cand_j] = ind2sub(size(errs),find(errs(:) < prctile(errs(:),10)));
    for i = 1:size(cand_i,1)
        scalefactor = [nan nan];
        for k = 1:2
            [a, b, c] = sph2cart(tmp(cand_i(i)), tmp(cand_j(i)),1);
            L = scaled(:,end) == k;
            predr = 1./(a.*cos(ph(L)).*cos(th(L))+b.*cos(ph(L)).*sin(th(L))+c.*sin(ph(L)));
            scalefactor(k) = geomean(abs(predr)./r(L));
        end
        
        initplaneguess = [a;b;c].*scalefactor(1);
        initscaleguess = scalefactor(2).^2/scalefactor(1).^2;  % not sure if this is right.
        initquadguess = [initplaneguess(1)^2 initplaneguess(2)^2 initplaneguess(3)^2 initplaneguess(1)*initplaneguess(2) initplaneguess(1)*initplaneguess(3) initplaneguess(2)*initplaneguess(3)]';
        options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-6,'TolX',10^-6,'Display','off');
        [tmpquadparams, tmpquadSSE, exitflag] = fminsearch(@(x) surfacefiterr3(newscaled, x, Loog),[initquadguess' initscaleguess],options);
        if (tmpquadSSE < constrquadSSE(cellcounter))
            constrquadparams(cellcounter,:) = tmpquadparams;
            constrquadSSE(cellcounter) = tmpquadSSE;
        end
    end
    for i = 1:2
        shortfn = NT{i}.sum.fileName(find(NT{i}.sum.fileName == filesep,1,'last')+1:find(NT{i}.sum.fileName == '.',1,'last')-1);
        Lplane(cellcounter) = Lplane(cellcounter) | ismember(shortfn,[planefnames{:}]');
    end
end

% Setting up the figure (for the EPS part)
figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);

%Histograms of correlation coefficients
axes('Position',[1,5,2*sqrt(2),2]); hold on;
r = zeros(size(allconeweights,1),1);
for i = 1:size(allconeweights,1)
    tmp = corrcoef([allconeweights(i,[1 2 3])', allconeweights(i,[4 5 6])']);
    r(i) = tmp(1,2);
end
truemeanr = mean(r);
[n,x] = hist(r,[-1:.1:11]);
bar(x,n,'FaceColor','black')
set(gca,'XTick',[-1 -.5 0 .5 1],'xlim',[-1.1,1.1]);
xlabel('Cone weight correlation','FontSize',14);
ylabel('Count','FontSize',14);

% How many nonplanar cells have correlation coeffiecients > .85
sum(Lplane)
sum(r(logical(Lplane)) > 0.85);

% Permutation test on correlation coefficients
niter = 10;
meanrs = zeros(1,niter);
for i = 1:niter
    x = randperm(size(allconeweights,1));
    for j = 1:size(allconeweights,1)
        tmp = corrcoef([allconeweights(j,[1 2 3])', allconeweights(x(j),[4 5 6])']);
        r(j) = tmp(1,2);
    end
    meanrs(i) = mean(r);
end
p = sum(meanrs >= truemeanr)./niter;
disp(['p-value from permutation test: ',num2str(p)]);

% Trying some F-tests to compare nested (non-linear) models
for i = 1:size(allscaled,1)
    n1 = sum(~allscaled{i,1}(:,4));  % Only considering ~OOGs in df calculation
    n2 = sum(~allscaled{i,2}(:,4));  % Only considering ~OOGs in df calculation
    DF1 = (n1+n2)-12;
    DF2 = (n1+n2)-7;
    SS1(i) = sum(allquadSSEs(i,:));
    SS2(i) = constrquadSSE(i);
    num = (SS2-SS1)/(12-7);
    den = SS1/DF1;
    F(i) = num/den;
    p(i) = 1-fcdf(num/den,12-7,DF1);
end

axes('Position',[1,2,2*sqrt(2),2]); hold on;
firingrateratio = max(thresholds,[],2)./min(thresholds,[],2);
L = p < 0.01;
plot(log10(firingrateratio(L))',SS1(L)./SS2(L),'ro','MarkerFaceColor','red','MarkerSize',5)
plot(log10(firingrateratio(~L))',SS1(~L)./SS2(~L),'ko','MarkerFaceColor','black','MarkerSize',5)
set(gca,'XTick',log10([1 2 3 4]),'XTickLabel',[1 2 3 4]);
set(gca,'Xlim',log10([.9 4.1]));
xlabel('Ratio of target firing rates','FontSize',14);
ylabel('SSE_{full}/SSE_{nested}','FontSize',14);

% Correlation of deviance and the ratio of thresholds
[r,p1] = corr([log10(firingrateratio),(SS1./SS2)'],'type','Spearman')

% Now plotting the inidividual cell examples
filenames = {'K033111002','K072009003','S081010005','S112410005',...
                'S050410006','S021810008','S012910002','K101409003'};
% filenames = {'K033111002','S081010005','S012810004','K101409003',...
%                'S050410006','S021810008','S112410005','S020310002'};
% % F-test on example neurons
% for i = 1:length(filenames)
%    cellnumber = ceil(find(strcmp([fnames{:}],filenames(i)))/2);
%    [p(cellnumber) thresholds(cellnumber,:)]
% end
% % Setting up the figure (for the bitmapped part)
figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);

views = [-75 32; 55 10; -121 24;-142 26; 53 18; 39 10; 53 18; 35 24];
%views = [-75 32; -121 24; -66 -16;-142 26; 53 18; 39 10; 53 18; 35 24];

plotcolors = [0 0 1; 1 0 0];
for whichcell = 1:length(filenames)
    tmp = nan*ones(size(fnames,1),1);
    for i = 1:size(fnames,1)
        tmp(i) = strcmp(fnames{i}(1),filenames(whichcell));
    end
    neuronidx = find(tmp);
    AXISCOL = 1+3*round((whichcell-1)/length(filenames));
    AXISROW = 2+2*mod((whichcell-1),length(filenames)/2);

    axes('Position',[AXISCOL,AXISROW,1.3,1.3]); hold on;
    %subplot(4,2,whichcell); hold on;
    tmpthresh = thresholds(neuronidx,~isnan(thresholds(neuronidx,:)));
    [y,idx] = sort(tmpthresh);
    p = [];
    for i = idx
        if (~isempty(allscaled{neuronidx,i}))
            scaled = allscaled{neuronidx,i}(:,1:3);
            Loog = allscaled{neuronidx,i}(:,4);
            plotlim  = max(abs(scaled(~Loog,:)))*1.4;
            [xx yy zz] = meshgrid(linspace(-plotlim(1),plotlim(1),40),...
                linspace(-plotlim(2),plotlim(2),40),...
                linspace(-plotlim(3),plotlim(3),40));
            xformedxyz = [xx(:) yy(:) zz(:)];
            variables = [xformedxyz(:,1).^2 xformedxyz(:,2).^2 xformedxyz(:,3).^2 2*xformedxyz(:,1).*xformedxyz(:,2) 2*xformedxyz(:,1).*xformedxyz(:,3) 2*xformedxyz(:,2).*xformedxyz(:,3)];
            fr = variables*allquadparams(neuronidx,NQUADPARAMS*(i-1)+[1:NQUADPARAMS])';
            surfstruct = isosurface(xx,yy,zz,reshape(fr,size(xx)), 1);
            p = patch(surfstruct);
            col = plotcolors(i==idx,:);
            plot3(scaled(~Loog,1),scaled(~Loog,2),scaled(~Loog,3),'o','MarkerFaceColor',col,'MarkerEdgeColor','none','MarkerSize',4);
            plot3(-scaled(~Loog,1),-scaled(~Loog,2),-scaled(~Loog,3),'o','MarkerFaceColor',col,'MarkerEdgeColor','none','MarkerSize',4);
            set(p,'FaceAlpha',.2,'FaceColor',max(col,.5),'Edgealpha',0);
        end
    end
    h = title({['\color{blue}',num2str(round(y(1)))],...
               ['\color{red}',num2str(round(y(2)))]},'FontSize',14);
    set(h,'Units','inches','Position',[2 .75])
    axis vis3d;
    set(gca,'XTick',[],'YTick',[],'ZTick',[],'View',views(whichcell,:));
    set(gca,'XLim',[-plotlim(1) plotlim(1)],'YLim',[-plotlim(2) plotlim(2)],'ZLim',[-plotlim(3) plotlim(3)]);
    camlight;
    lighting phong;
    axis square;
end
%print -dtiff -r300 -cmyk junk


% How often do L and M cone weights change opponency?
% M cone weights are positive by convention
% Lconesigns = [sign(allconeweights(:,1)) sign(allconeweights(:,4))];
% L = Lconesigns(:,1).*Lconesigns(:,2) == -1;
% sum(L)
% % which surface types are these?
surfacetypes = [];
for i = 1:size(allquadparams,1)
    A1 = [allquadparams(i,1) allquadparams(i,4) allquadparams(i,5);...
        allquadparams(i,4) allquadparams(i,2) allquadparams(i,6);...
        allquadparams(i,5) allquadparams(i,6) allquadparams(i,3)];              
    A2 = [allquadparams(i,7) allquadparams(i,10) allquadparams(i,11);...
        allquadparams(i,10) allquadparams(i,8) allquadparams(i,12);...
        allquadparams(i,11) allquadparams(i,12) allquadparams(i,9)];              
    surfacetypes = [surfacetypes; sum(eig(A1)>0) sum(eig(A2)>0)];
end        
sum(surfacetypes(:,1) == surfacetypes(:,2))



%%
% Section 10.2 Effects of changing threshold version 2
% Only plotting data from a few example cells.

PLOTSURFS = 1; % 0 = draw points (for eps file) 1 = draw surfaces (for tiff file)

load ('T_cones_smj10');
[fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\NTmultiples.txt');
% These are good
filenames = {'S012910002','S050410006','S012910002','K101309003',...
                'K033111002','S040210002','S021810008','K050211002'};
views = [-33 22;-152 36; 0 0; -128 20;...
         -94 30;-21 14; 61 -72; 24 -18];
     
     
filenames = {'S012910002','S050410006','K082709004','K101309003',...
                'K033111002','S040210002','S021810008','K050211002'};
views = [-33 22;37 -42; 14 80; -128 20;...
         -94 30;-21 14;124 -68; 24 -18];
     
PLOTLIMSCALEFACTOR = .80;
% Now plotting the inidividual cell examples
figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
plotcolors = [0 0.6 0.5;0.8 0.5 0];

for cellcounter = 1:length(filenames)
    cellcounter
    tmp = nan*ones(size(fnames,1),1);
    for i = 1:size(fnames,1)
        tmp(i) = strcmp(fnames{i}(1),filenames(cellcounter)) | strcmp(fnames{i}(2),filenames(cellcounter));
    end
    neuronidx = find(tmp);

    NT = {}; out = {};
    for i = 1:2
        filename = findfile(char(fnames{neuronidx}(i)));
        NT{i} = nex2stro(filename);
        lmsidxs = [find(strcmp(NT{i}.sum.trialFields(1,:),'lcont'))...
            find(strcmp(NT{i}.sum.trialFields(1,:),'mcont'))...
            find(strcmp(NT{i}.sum.trialFields(1,:),'scont'))];
        % Converting to 10 degree fundamentals
        fundamentals = NT{i}.sum.exptParams.fundamentals;
        fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
        mon_spd = NT{i}.sum.exptParams.mon_spd;
        mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
        mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
        M = fundamentals'*mon_spd;
        M10 = T_cones_smj10*mon_spd;
        NT{i}.trial(:,lmsidxs) = ConvertConeContrastBasis(M, M10, NT{i}.sum.exptParams.bkgndrgb, NT{i}.trial(:,lmsidxs));
        out{i} = NTpreprocess(NT{i},0,Inf);
    end
    
    % Changing the order so that the lower threshold is first
    tmpthresh = [NT{1}.sum.exptParams.threshold NT{2}.sum.exptParams.threshold];
    [y,idx] = sort(tmpthresh);
    scaled = [];
    for i = 1:2
        scaled = [scaled; out{idx(i)}(:,[2:4]).*repmat(out{idx(i)}(:,5), 1,3) logical(out{idx(i)}(:,7)) i*ones(size(out{idx(i)},1),1)];
    end
    NT = [NT(idx)];
    % Doing the constrained analysis (both files at once)
    Loog = scaled(:,4);
    [v,d] = eig(cov([scaled(~Loog,[1 2 3]); -scaled(~Loog,[1 2 3])]));
    d = diag(d);
    if (min(d) < 2*eps)
        disp('Too few data points for whitening');
        whtmat = eye(3);
    else
        whtmat = v*diag(sqrt(1./d));
    end
    newscaled = [scaled(:,[1 2 3])*whtmat scaled(:,4)];
    [th,ph,r] = cart2sph(newscaled(:,1),newscaled(:,2),newscaled(:,3));
    
    % Doing the fitting
    errs = zeros(50,50);
    tmp = linspace(0,pi,size(errs,1));
    for i = 1:length(tmp)
        for j = 1:length(tmp)
            tmpresid = 0;
            for k = 1:2
                L = scaled(:,end) == k;
                [tmpa, tmpb, tmpc] = sph2cart(tmp(i),tmp(j),1);
                predr = 1./(tmpa.*cos(ph).*cos(th)+tmpb.*cos(ph).*sin(th)+tmpc.*sin(ph));
                predr = predr.*geomean(r(L&~Loog)./abs(predr(L&~Loog)));
                resid = log(abs(predr(L)))-log(r(L));
                resid(Loog(L)&(abs(predr(L)) > r(L))) = 0;
                tmpresid = tmpresid + sum(resid.^2);
            end
            errs(i,j) = tmpresid;
        end
    end
    [cand_i, cand_j] = ind2sub(size(errs),find(errs(:) < prctile(errs(:),10)));
    quadparams = []; quadSSE = Inf;
    for i = 1:size(cand_i,1)
        scalefactor = [nan nan];
        for k = 1:2
            [a, b, c] = sph2cart(tmp(cand_i(i)), tmp(cand_j(i)),1);
            L = scaled(:,end) == k;
            predr = 1./(a.*cos(ph(L)).*cos(th(L))+b.*cos(ph(L)).*sin(th(L))+c.*sin(ph(L)));
            scalefactor(k) = geomean(abs(predr)./r(L));
        end
        
        initplaneguess = [a;b;c].*scalefactor(1);
        initscaleguess = scalefactor(2).^2/scalefactor(1).^2;  % not sure if this is right.
        initquadguess = [initplaneguess(1)^2 initplaneguess(2)^2 initplaneguess(3)^2 initplaneguess(1)*initplaneguess(2) initplaneguess(1)*initplaneguess(3) initplaneguess(2)*initplaneguess(3)]';
        options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-6,'TolX',10^-6,'Display','off');
        [tmpquadparams, tmpquadSSE, exitflag] = fminsearch(@(x) surfacefiterr3([newscaled(:,[1 2 3]) scaled(:,end)], x, Loog),[initquadguess' initscaleguess],options);
        if (tmpquadSSE < quadSSE)
            quadparams = tmpquadparams;
            quadSSE = tmpquadSSE;
        end
    end
    
    % Doing the plotting
    AXISCOL = .5+2*round((cellcounter-1)/length(filenames));
    AXISROW = .5+2*mod((cellcounter-1),length(filenames)/2);
    axes('Position',[AXISCOL,AXISROW,2,2]); hold on;

    if (sum(Loog) > 3)
        plotlim  = max(abs(newscaled(~Loog,:)))*PLOTLIMSCALEFACTOR;
    else
        plotlim  = max(abs(newscaled(~Loog,:)))*.9;    
    end
    if (PLOTSURFS)
        [xx yy zz] = meshgrid(linspace(-plotlim(1),plotlim(1),50),...
        linspace(-plotlim(2),plotlim(2),50),...
        linspace(-plotlim(3),plotlim(3),50));
        xformedxyz = [xx(:) yy(:) zz(:)];
        variables = [xformedxyz(:,1).^2 xformedxyz(:,2).^2 xformedxyz(:,3).^2 2*xformedxyz(:,1).*xformedxyz(:,2) 2*xformedxyz(:,1).*xformedxyz(:,3) 2*xformedxyz(:,2).*xformedxyz(:,3)];
        fr = variables*quadparams(1:6)';
      %  fr = variables*quadparams(1:6)'*sqrt(quadparams(7));
        
        surfstruct = isosurface(xx,yy,zz,reshape(fr,size(xx)), 1);
        hsurf = patch(surfstruct); 
        set(hsurf,'FaceVertexCData',repmat([.6 1 1],size(getfield(get(hsurf),'Vertices'),1),1))
        set(hsurf,'CDataMapping','direct');
        set(hsurf,'FaceColor','interp','EdgeColor','none');
        set(hsurf,'SpecularColorReflectance',.1);
        set(hsurf,'SpecularExponent',.2,'SpecularStrength',1);
        set(hsurf,'DiffuseStrength',.8,'AmbientStrength',.4);
        set(hsurf,'FaceAlpha',.3);
      %  camlight;
        camlight(45,45);
%        camlight(-45,-45);
        camlight(0, -90)
        lighting phong;
        material shiny
    end
%    if (~PLOTSURFS)
        for i = 1:2
            col = plotcolors(i,:);
            L = scaled(:,5) == i;
            if (i == 1)
                plot3(newscaled(L&~Loog,1),newscaled(L&~Loog,2),newscaled(L&~Loog,3),'.','MarkerEdgeColor',col,'MarkerSize',8);
                plot3(-newscaled(L&~Loog,1),-newscaled(L&~Loog,2),-newscaled(L&~Loog,3),'.','MarkerEdgeColor',col,'MarkerSize',8);
                
                % plot3([-newscaled(L&Loog,1) newscaled(L&Loog,1)]',[-newscaled(L&Loog,2) newscaled(L&Loog,2)]',[-newscaled(L&Loog,3) newscaled(L&Loog,3)]','-','Color',max(col,.8));
                plot3([max(-newscaled(L&Loog,1),-plotlim(1)) min(newscaled(L&Loog,1), plotlim(1))]',[max(-newscaled(L&Loog,2),-plotlim(2)) min(newscaled(L&Loog,2),plotlim(2))]',[max(-newscaled(L&Loog,3),-plotlim(3)) min(newscaled(L&Loog,3),plotlim(3))]','-','Color',min(col+.6,1));
            else
                % The other group of points
                   plot3(newscaled(L&~Loog,1).*sqrt(quadparams(7)),newscaled(L&~Loog,2).*sqrt(quadparams(7)),newscaled(L&~Loog,3).*sqrt(quadparams(7)),'.','MarkerEdgeColor',col,'MarkerSize',8);
                   plot3(-newscaled(L&~Loog,1).*sqrt(quadparams(7)),-newscaled(L&~Loog,2).*sqrt(quadparams(7)),-newscaled(L&~Loog,3).*sqrt(quadparams(7)),'.','MarkerEdgeColor',col,'MarkerSize',8)
                % OOGs
                plot3([max(-newscaled(L&Loog,1),-plotlim(1)) min(newscaled(L&Loog,1), plotlim(1))]',[max(-newscaled(L&Loog,2),-plotlim(2)) min(newscaled(L&Loog,2),plotlim(2))]',[max(-newscaled(L&Loog,3),-plotlim(3)) min(newscaled(L&Loog,3),plotlim(3))]','-','Color',min(col+.6,1));
                % The small (unscaled) points
                %plot3(newscaled(L&~Loog,1),newscaled(L&~Loog,2),newscaled(
                %L&~Loog,3),'.','MarkerEdgeColor',max(col,.5),'MarkerSize',6);
                %plot3(-newscaled(L&~Loog,1),-newscaled(L&~Loog,2),-newscaled(L&~Loog,3),'.','MarkerEdgeColor',max(col,.5),'MarkerSize',6);
                % Lines connecting the unscaled to the scaled points
                %  plot3([newscaled(L&~Loog,1) newscaled(L&~Loog,1).*sqrt(quadparams(7))]',...
                %      [newscaled(L&~Loog,2) newscaled(L&~Loog,2).*sqrt(quadparams(7))]',...
                %      [newscaled(L&~Loog,3) newscaled(L&~Loog,3).*sqrt(quadparams(7))]','-','Color',col);
                %  plot3(-[newscaled(L&~Loog,1) newscaled(L&~Loog,1).*sqrt(quadparams(7))]',...
                %      -[newscaled(L&~Loog,2) newscaled(L&~Loog,2).*sqrt(quadparams(7))]',...
                %      -[newscaled(L&~Loog,3) newscaled(L&~Loog,3).*sqrt(quadparams(7))]','-','Color',col);
                
            end
       % end
        
        drawnow;
        %h = title({['\color{blue}',num2str(round(NT{1}.sum.exptParams.threshold))],...
        %           ['\color{red}',num2str(round(NT{2}.sum.exptParams.threshold))]},'FontSize',14);
        %h = text(0,0,{['\color{blue}',num2str(round(NT{1}.sum.exptParams.threshold))],...
        %    ['\color{red}',num2str(round(NT{2}.sum.exptParams.threshold))]},'FontSize',14);
        
        %set(h,'Units','inches','Position',[2 .75])
    end
    set(gca,'Visible','off')
    set(gca,'XTick',[],'YTick',[],'ZTick',[],'View',views(cellcounter,:));
    set(gca,'XLim',[-plotlim(1) plotlim(1)],'YLim',[-plotlim(2) plotlim(2)],'ZLim',[-plotlim(3) plotlim(3)]);
    axis vis3d;

    lighting phong;
    axis square;
   % keyboard
end
%print -dtiff -r600 -cmyk junk
% Need to make symbols transparent in Illustrator


%% Section 11: Modulation ratios as a function of surface type

data = [];
[planefnames] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\plane0.01.txt');
[fnames, spikeIdx] = fnamesFromTxt2();
for cellcounter = 1:size(fnames,1)
    NT = {}; GT = {};
    for i = 1:size(fnames{cellcounter},2)
        filename = findfile(char(fnames{cellcounter}(i)));
        paradigmID = getparadigmID(filename);
        if paradigmID == 150
            GT = nex2stro(filename);
            GTstruct = getGratingTuning(GT,1);
        else
            NT = nex2stro(filename);
        end
    end
    out = NTpreprocess(NT,0,Inf);
    
    scaled = out(:,[2:4]).*repmat(out(:,5), 1,3);
    Loog = logical(out(:,7));
    [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);

    % What type of quadric?
    A = [quadparams(1) quadparams(4) quadparams(5);...
        quadparams(4) quadparams(2) quadparams(6);...
        quadparams(5) quadparams(6) quadparams(3)];
    % Transforming quadparams into un-whitened space
    [evecs,evals] = eig(A);
    evals = diag(evals);
    if (sum(evals<0) == 1) % 1 sheet
        WHICHSURF = 1;
    elseif (sum(evals<0) == 2) % 2 sheets
        WHICHSURF = 2;
    else %(all(diag(evals)) % ellipsoid
        WHICHSURF = 3;
    end
    shortfn = NT.sum.fileName(find(NT.sum.fileName == filesep,1,'last')+1:find(NT.sum.fileName == '.',1,'last')-1);
    Lplane = ismember(shortfn,[planefnames{:}]');

    data = [data; WHICHSURF GTstruct.modulationratio Lplane];
end

nbins = 14;
bins = linspace(0,2,nbins);
n = zeros(3,nbins);
Lplane = logical(data(:,3));
%Lplane = logical(zeros(size(data,1),1));
X = [];
for whichsurf = 1:3
    L = data(:,1) == whichsurf;
    L = L & ~Lplane;
    n(whichsurf,:) = hist(data(L,2),bins);
    X = [X; data(L,[1 2])];
end
X = [X; 4*ones(sum(Lplane),1) data(Lplane,2)];

Lnotnan = ~isnan(X(:,2));
[p,anovatab,stats] = anova1(X(:,2),X(:,1))
geomean(X(X(:,1) < 3 & Lnotnan,2))
geomean(X(X(:,1) == 1 & Lnotnan,2)) % 1 sheet
geomean(X(X(:,1) == 2 & Lnotnan,2)) % 2 sheet
geomean(X(X(:,1) == 3 & Lnotnan,2)) % ellipsoid
geomean(X(X(:,1) == 4 & Lnotnan,2)) % plane
comparison = multcompare(stats)
n(size(n,1)+1,:) = hist(data(Lplane,2),bins)
% Order is hyp 1, hyp2, ellipsoid, (plane)

% Setting up the figure
figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
colormap([ .75 .75 0; .75 .25 .25; .25 .75 .25; .25 .25 .75;]);
%colormap([.75 .25 .25; .25 .75 .25; .25 .25 .75]);

axes('Position',[2 2 4 4]);
bar(bins,n([3 1 2 4],:)','stacked')
%bar(bins,n([3 1 2],:)','stacked')
h = legend({'Ellipsoid','Hyperboloid 1','Hyperboloid 2','Plane'});
%h = legend({'Ellipsoid','Hyperboloid 1','Hyperboloid 2'});
set(gca,'XLim',[0 2.1]);
ylabel('Count','FontSize',14);
xlabel('Modulation ratio (F1/F0)','FontSize',14);

%%
% Section 12: Finding intertrial intervals and RF locations

[fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\Gratings&NT2.txt');
data = nan*ones(length(fnames),4);
for cellcounter = 1:size(fnames,1)
    NT = {}; GT = {};
    for i = 1:size(fnames{cellcounter},2)
        filename = findfile(char(fnames{cellcounter}(i)));
        paradigmID = getparadigmID(filename);
        if paradigmID == 150
            GT = nex2stro(filename);
            stimon_t = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'stim_on'));
            stimoff_t = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'stim_off'));
            data(cellcounter,1) = median(stimoff_t-stimon_t);
        else
            NT = nex2stro(filename);
            stimon_t = NT.trial(:,strcmp(NT.sum.trialFields(1,:),'stim_on'));
            stimoff_t = NT.trial(:,strcmp(NT.sum.trialFields(1,:),'stim_off'));
            data(cellcounter,2) = median(stimoff_t-stimon_t);
            data(cellcounter,3) = NT.sum.exptParams.rf_x;
            data(cellcounter,4) = NT.sum.exptParams.rf_y;
        end
    end
end

% RF stats
ecc = sqrt((data(:,3)/10).^2+(data(:,4)/10).^2)

%% 
% Section 13: Patrick's simulation
% Getting the parameters for the quadratic model
cd ('C:\Matlab Code\Analysis\Sandbox\Greg\PatrickSim');
WHICHINITDIRS = 3;  % 1 = standard, 2 = alternative, 3 = random
filename = 'K072909005';
filename = 'S031610003';
%filename = 'K090309006';
%Specify how many phases to be completed
phases=5;
%Specify number of grand loops
niter=20;
if (WHICHINITDIRS == 3)
    ncolordirs = 10;
else
    ncolordirs = 1;
end

NT = nex2stro(findfile(filename));
fundamentals = NT.sum.exptParams.fundamentals;
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = NT.sum.exptParams.mon_spd;
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
global M
load ('T_cones_smj10');
M = T_cones_smj10*mon_spd;
global bkgndLMS 
bkgndLMS=M* NT.sum.exptParams.bkgndrgb;

out = NTpreprocess(NT,0,Inf);  % .4, 2  or 0, Inf
scaled = out(:,[2:4]).*repmat(out(:,5), 1,3);
scaled = ConvertConeContrastBasis(fundamentals'*mon_spd, M, NT.sum.exptParams.bkgndrgb, scaled);
Loog = logical(out(:,7));
[planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);
A = [quadparams(1) quadparams(4) quadparams(5);...
    quadparams(4) quadparams(2) quadparams(6);...
    quadparams(5) quadparams(6) quadparams(3)];
B = xformmat*A*xformmat';
originalquadparams = [B(1,1) B(2,2) B(3,3) B(1,2) B(1,3) B(2,3)];
thresh = NT.sum.exptParams.threshold;

symbols = {'o','x','+','*','s','d','v','^','<','>','p','h'};
figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
axes('Position',[1 2 5*sqrt(2) 5]); hold on;
plot([-1 0 1 -1],[0 1 0 0],'k-');
colors = ['g', 'r', 'b'];

for initcolordircounter = 1:ncolordirs
    
    %Specify color directions (in cone contrast units)
    if (WHICHINITDIRS == 1)
        LMSdir(1,:)=[.025 .025 0];
        LMSdir(2,:)=[-.01 .01 0];
        LMSdir(3,:)=[0 0 .1];
    elseif (WHICHINITDIRS == 2)
        LMSdir(1,:)=[.025 .025 .025];
        LMSdir(2,:)=[-.01 .01 .1];
        LMSdir(3,:)=[.01 -.01 .1];
    else
        LMSdir(1,:)=normrnd(0,[.05 .05 .05]);
        LMSdir(2,:)=normrnd(0,[.05 .05 .05]);
        LMSdir(3,:)=normrnd(0,[.05 .05 .05]);
        LMSdir = MakeOrtho(LMSdir)/10;
    end
    %Specify staircasing variables
    numReversals=7;
    stepSize=0.5;
    scale=0.5;
    
    quadparamsmatrix=nan(niter,6);
    
    for q=1:niter
        q
        [trialspec]=trialspecs(LMSdir,numReversals,thresh,stepSize,scale,phases,originalquadparams);
        scaled=cat(1,trialspec.coordinates);
        Loog=[trialspec.oog]';
        [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);
        SSE=(planeSSE/quadSSE);
        A =  [quadparams(1) quadparams(4) quadparams(5);...
            quadparams(4) quadparams(2) quadparams(6);...
            quadparams(5) quadparams(6) quadparams(3)];
        % Transforming quadparams into un-whitened space
        B = xformmat*A*xformmat';
        quadparamsmatrix(q,:) = [B(1,1) B(2,2) B(3,3) B(1,2) B(1,3) B(2,3)];
    end
       
    for q = 1:niter
        quadparams = quadparamsmatrix(q,:);
        A =  [quadparams(1) quadparams(4) quadparams(5);...
            quadparams(4) quadparams(2) quadparams(6);...
            quadparams(5) quadparams(6) quadparams(3)];
        [evecs,evals] = eig(A);
        [evals,i] = sort(diag(evals),1,'ascend');
        evecs = evecs(:,i);
        for j = 1:3
            v = evecs(:,j);
            if (v(2) < 0) v = -v; end
            h = plot(v(1)./sum(abs(v)),v(2)./sum(abs(v)),[colors(j),'s'],'MarkerFaceColor',colors(j),'Marker',char(symbols(initcolordircounter)));
        end
    end
    for j = 1:3  % Plotting initial color directions
        v = LMSdir(j,:);
        if (v(2) < 0); v = -v; end
        h = plot(v(1)./sum(abs(v)),v(2)./sum(abs(v)),'k+','LineWidth',2,'Marker',char(symbols(initcolordircounter)));
    end
    set(gcf,'Name',[filename,': ',num2str(WHICHINITDIRS)]);
    drawnow;
end

% Plotting the true answers
quadparams = originalquadparams;
A =  [quadparams(1) quadparams(4) quadparams(5);...
    quadparams(4) quadparams(2) quadparams(6);...
    quadparams(5) quadparams(6) quadparams(3)];
[evecs,evals] = eig(A);
[evals,i] = sort(diag(evals),1,'ascend');
evecs = evecs(:,i);
for j = 1:3
    v = evecs(:,j);
    if (v(2) < 0) v = -v; end
    h = plot(v(1)./sum(abs(v)),v(2)./sum(abs(v)),[colors(j),'s']);
    set(h,'Color',get(h,'Color')/2);
end

%%
% Section 14
% Schematic model of hierarchical model showing how we can combine linear
% neurons to create nonlinear ones.

x = meshgrid([-1:.005:1]);
y = x';

figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);

cmap = [linspace(57,237,64)'/255,linspace(83,28,64)'/255,linspace(164,36,64)'/255];
colormap(cmap);

% Linear neuron 1
axes('Position',[1 1 2 2]); hold on;
imagesc(2*abs(x)-1);
[c,h] = contour(x,5);
set(h,'Linewidth',4,'Color','black');
axis tight;
set(gca,'XTick',[],'YTick',[]);

% Linear neuron 2
axes('Position',[4 1 2 2]); hold on;
imagesc(2*abs(y)-1);
[c,h] = contour(y,5);
set(h,'Linewidth',4,'Color','black');
axis tight;
set(gca,'XTick',[],'YTick',[]);

% Nonlinear neuron 1
axes('Position',[1 4 2 2]); hold on;
z = sqrt(2*x.^2+y.^2);
imagesc(z);
[c,h] = contour(z,5);
set(h,'Linewidth',4,'Color','black');
axis tight;
set(gca,'XTick',[],'YTick',[]);

% Nonlinear neuron 2
axes('Position',[4 4 2 2]); hold on;
z = sqrt(x.^2./(y.^2+.5));
imagesc(z);
[c,h] = contour(z,5);
set(h,'Linewidth',4,'Color','black');
axis tight;
set(gca,'XTick',[],'YTick',[]);

% scalebar
axes('Position',[6.3 4 .2 2]); hold on;
imagesc(y)
axis tight;
set(gca,'XTick',[],'YTick',[],'Box','off')

%%
% Section 15
% Can I describe the color direction selecting algorithm in 10 deg cone
% contrast space?  In RGB space?
% Bottom line: Averaging lights in RGB intensity space and converting to cone contrasts is
% the same as averaging cone contrasts, but averaging 2 degree cone
% contrasts and converting to 10 degree cone contrasts is *not* the same as
% averaging 10 degree cone contrasts.

filename = 'K082609010.nex';  % pan color
NT = nex2stro(findfile(filename));
bkgndrgb = NT.sum.exptParams.bkgndrgb;

fundamentals = NT.sum.exptParams.fundamentals;
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = NT.sum.exptParams.mon_spd;
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;
load ('T_cones_smj10');
M10 = T_cones_smj10*mon_spd;

bkgndlms2 = M*bkgndrgb;
bkgndlms10 = M10*bkgndrgb;

rgb = unifrnd(0,1,3,3);  % each column is a stimulus.
midpoint_rgb = mean(rgb,2);  % mean in RGB space
% Computing cone contrasts from midpoint_rgb
lms = M*midpoint_rgb;
(lms-bkgndlms2)./bkgndlms2  % midpoint RGB -> 2 deg LMS cc
lms = M10*midpoint_rgb;
(lms-bkgndlms10)./bkgndlms10  % midpoint RGB -> 10 deg LMS cc

% Finding the midpoint in 2 deg cone contrast space
lms = M*rgb;
conecontrasts2 = (lms-repmat(bkgndlms2,[1 3]))./repmat(bkgndlms2,[1 3]);
midpoint_cc2 = mean(conecontrasts2,2)

% Finding the midpoint in 10 deg cone contrast space
lms = M10*rgb;
conecontrasts10 = (lms-repmat(bkgndlms10,[1 3]))./repmat(bkgndlms10,[1 3]);
midpoint_cc10 = mean(conecontrasts10,2)

% Finding the midpoint in 2 deg contrast space and converting to 10 deg
% cone contrast space
conecontrasts10_new= ConvertConeContrastBasis(M, M10,bkgndrgb,midpoint_cc2)

