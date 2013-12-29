% NeuroThresh Population analyses
%
% Section 1: Trying to predict responses to chromatic gratings from the
% NeuroThresh data.  Either using a plane fit or a quadratic fit.
%
% Section 2: Looking at superimposed isoresponse surafaces from a few
% different cells.
%
% Section 3: Bootstrap tests on a bunch of cells to see which ones are
% better described by a curved surface than a plane.  Also classifying
% quadric surfaces on the basis of their eignevalues.  Also using a
% conventional F-test to separate planar from non-planar surfaces - very
% similar to the bootstrap test.
%
% Section 3.1 Computing pseudo-R^2, a generalized goodness of fit statistic
% for nonlinear models.
%
% Section 4: Looking at isoresponse surfaces for a few cells in two 
% color spaces  - one defined by the 2° fundamentals and one defined
% by the 10° fundamentals.
%
% Section 4.1: Looking at coneweights for luminance cells using several
% different versions of the fundamentals.  Also principal axis of
% ellipsoidal cells.
%
% Section 5: Ellipsoidal fits of isoresponse surfaces of "pan-color"
% cells.  Do these ellipsoids have a systematic orientation?
%
% Section 6: Comparing planar and quadratic fits using some kind of 
% F-type statistic.  Are the more planar cells tuned for a particular
% subset of color directions?
%
% Section 7: What are the "preferred color direction" of the pancolor cells
% assessed with grating data? 
%
% Section 8: Looking at the distributions of preferred color directions as
% assessed with gratings and as assessed with NeuroThresh.
%
% Section 8.1: Do Neurothresh and Gratings agree with respect to the L to M
% cone ratio for non-opponent cells?
%
% Section 9: Looking at the F1/F0 modulation ratio for pancolor cells and
% comparing it to non-pancolor cells.
%
% Section 9.1: Looking at the F1/F0 modulation ratio broken down by quadric
% fit type.  Also keeping track of preferred SF in case this matters.
%
% Section 10: Looking at the relationships between direction selectivity,
% orientation selectivity, double-opponency, and threshold as a function
% of isoresonse surface shapes.
%
% Section 11: For some cells, the NeuroThresh data does a poor job at
% explaining the Gratings data.  Is this because the range of
% contrasts/firing rates was very different between these two paradigms?
% Is it because the "optimal" size stimulus was used in the grating color
% tuning?
%
% Section 12: Looking at the curvature of the fitted isoresponse surface in
% the L-M direction and the S direction for luminance-tuned neurons.
%
% Section 13: Looking at the distance between the best fit plane from the
% NeuroThresh data and the isodetection surface as measured by DTNT.
%
% Section 14: Trying to classify isoresponse surface shapes and look for
% consistency in the orientations of these shapes.  Analysis of
% eigenvectors.
%
% Section 15: Looking at the consistency (or lack thereof) of isoresponse
% surfaces across dfferent choices of threshold.
%
% Section 15.1: For the neurons that were tested with multiple thresholds
% we should have a principled approach to selecting one for predicting the
% gratings data.  Maybe the one for which the staircase terminations are
% somehow closest to the gratings.
%
% Section 15.2: Dividing data sets in half and seeing how often quadratic
% isoresponse surface types change.
%
% Section 16: Permutation tests to test the hypothesis that the quadratic
% surface predicts the data better (in the sense of correlation
% coefficient) than the plane does.
%
% Section 17: Trying to fit a model to NTmultiples data that constrains
% isoresponse surfaces to be the same up to a scale factor.
%
% Section 17.1: Comparing constrained fits to unconstrained fits
% quantitatively.  Also randomization tests?  Constraint is that the two
% surfaces have to be scaled versions of each other.
%
% Section 18: How many staircases get chucked because of differences 
% between online and offline spike sorting?  How many stimulus presented
% per trial?
%
% Section 19: Finding the cone contrasts of the gratings in a cone
% contrast space defined by 10 degree cone fundamentals.
%
% Section 20: Are the steps after the first reversal not decremented?
%
% Section 21: Comparing cross-validation errors between linear and
% quadratic predictions. A new way to sort "planar" from "nonplanar"
% neurons.
%
% Section 22: Finding the color of the stimulus used to measure 
% orientation and SF tuning in gratings files.
% -----------------------------
%%
% Section 1.  Predicting responses obtained in grating.d paradigm from data obtained
% in NeuroThresh paradigm.

CORRTYPE = 'Pearson';  % 'Pearson' or 'Spearman'
data = [];
[fnames, spikeIdx] = fnamesFromTxt2();
PLOTSURFS = 0;
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
    out = NTpreprocess(NT,0,Inf);  % .4, 2  or 0, Inf
        
    scaled = out(:,[2:4]).*repmat(out(:,5), 1,3);
    Loog = logical(out(:,7));
    if (sum(~Loog) < 3)
        disp('Fewer than 3 points in gamut');
        continue;
    end
    [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);
    xyz = scaled*xformmat;
    
    colors = GTstruct.color.colors;
    rotcols = colors*xformmat;
    planepredresp = NT.sum.exptParams.threshold*abs(rotcols*planeparams);
    trueresp = GTstruct.color.colresp(:,1);

    % Finding predictions from quadratic model
    [th,ph,r] = cart2sph(rotcols(:,1),rotcols(:,2),rotcols(:,3));
    predr2 = 1./(quadparams(1).*(cos(ph).*cos(th)).^2 +...
    quadparams(2).*(cos(ph).*sin(th)).^2 +...
    quadparams(3).*sin(ph).^2 + ...
    2*quadparams(4).*cos(ph).*cos(th).*cos(ph).*sin(th) +...
    2*quadparams(5).*cos(ph).*cos(th).*sin(ph) +...
    2*quadparams(6).*cos(ph).*sin(th).*sin(ph));
    
    predr2(predr2<0) = Inf;
    quadpredresp = NT.sum.exptParams.threshold.*(r./sqrt(predr2));
    
    % What type of quadric?
    A = [quadparams(1) quadparams(4) quadparams(5);...
        quadparams(4) quadparams(2) quadparams(6);...
        quadparams(5) quadparams(6) quadparams(3)];
    d = eig(A);
    if (sum(d<0) == 1)
        disp('hyperboloid of one sheet');
        plotlim = max(abs(xyz(:)))/1.1;
        WHICHSURF = 1;
    elseif (sum(d<0) == 2)
        disp('hyperboloid of two sheets');
        plotlim = max(abs(xyz(:)))/3;
        WHICHSURF = 2;
    else %(all(d>0))
        disp('Ellipsoid');
        plotlim = max(abs(xyz(:)))*1.1;
        WHICHSURF = 3;
    end
    
    % Computing correlations
    correlations = corr([planepredresp, quadpredresp, trueresp],'Type',CORRTYPE);

    % Computing "MSE"
%    planeMSE = mean((trueresp-planepredresp).^2./trueresp);
%    quadMSE = mean((trueresp-quadpredresp).^2./trueresp);
    planeMSE = mean((trueresp-planepredresp).^2);
    quadMSE = mean((trueresp-quadpredresp).^2);

    if (PLOTSURFS)
        ylabels = {'Linear','Quadratic'};
        %    figure('position',[251 77 325 613]);
        clf
        set(gcf,'Name',NT.sum.fileName(find(NT.sum.fileName == '\',1,'last')+1:end));
        for j = 1:2
            subplot(3,1,j); hold on;
            for i = 1:length(trueresp)
                if (j == 1)
                    h = plot(trueresp(i), planepredresp(i), 'ko','Markersize',12);
                elseif (j == 2)
                    h = plot(trueresp(i), quadpredresp(i), 'ko','Markersize',12);
                else % j == 3
                    h = plot(trueresp(i), nppredresp(i), 'ko','Markersize',12);
                end
                col = GTstruct.color.colors(i,:).*[1 1 .2];
                col = col./(2*max(abs(col))) + .5;
                set(h,'MarkerFaceColor',col,'ButtonDownFcn',['disp(num2str(GTstruct.color.colors(',num2str(i),',:)))']);
            end
            xlabel('Response'); ylabel([ylabels{j},' prediction']);
        end
        subplot(3,1,1);
        title(['r = ',num2str(correlations(1,3))]);
        subplot(3,1,2);
        title(['r = ',num2str(correlations(2,3))]);
        
        % Plotting the fits
        subplot(3,1,3); hold on;
        xyz = scaled(:,[1 2 3])*xformmat;
        plot3(xyz(~Loog,1),xyz(~Loog,2),xyz(~Loog,3),'k.');
        plot3(-xyz(~Loog,1),-xyz(~Loog,2),-xyz(~Loog,3),'k.');
        plot3(xyz(Loog,1),xyz(Loog,2),xyz(Loog,3),'y.');
        plot3(-xyz(Loog,1),-xyz(Loog,2),-xyz(Loog,3),'y.');
        
        % Now plotting surfaces
        [xx yy zz] = meshgrid(linspace(-plotlim,plotlim,20),...
            linspace(-plotlim,plotlim,20),...
            linspace(-plotlim,plotlim,20));
        xformedxyz = [xx(:) yy(:) zz(:)];
        variables = [xformedxyz(:,1).^2 xformedxyz(:,2).^2 xformedxyz(:,3).^2 2*xformedxyz(:,1).*xformedxyz(:,2) 2*xformedxyz(:,1).*xformedxyz(:,3) 2*xformedxyz(:,2).*xformedxyz(:,3)];
        for PLANE = 0:1
            if (PLANE)
                coefficients = [planeparams'.*planeparams' planeparams(1)*planeparams(2) planeparams(1)*planeparams(3) planeparams(2)*planeparams(3)]';
                patchcolor = 'green';
            else
                coefficients = quadparams;
                patchcolor = 'blue';
            end
            
            fr = variables*coefficients;
            surfstruct = isosurface(xx,yy,zz,reshape(fr,size(xx)), 1);
            %  surfstruct.vertices = surfstruct.vertices/xformmat;
            
            p = patch(surfstruct);
            set(p, 'FaceColor', 'blue', 'EdgeColor', 'none');
            set(p,'FaceAlpha',.2,'FaceColor',patchcolor,'Edgealpha',0);
            axis vis3d;
            camlight;
            lighting phong;
            axis square;
        end
    end
    data = [data; correlations(1,3), correlations(2,3), correlations(1,2) (planeparams'*xformmat') WHICHSURF planeMSE quadMSE];
end
% 1) r(plane,true)
% 2) r(quad,true)
% 3) r(plane,quad)

Lcor = data(:,3) <= 1 & ~isnan(data(:,2));  % threshold on correlation coefficient btn lin/quad preds (and no nans)
%Llum = abs(sum(coneweights(:,[1 2]),2)) < .7
%Lcor = Lcor & Llum;

figure; subplot(2,2,1); hold on;
plot(data(Lcor,1),data(Lcor,2),'r.');
plot(data(Lcor,1),data(Lcor,2),'k.');
axis equal;
plot([-1 1],[-1 1],'k-')
xlabel('Plane fit');
ylabel('Quad fit');
[h,p] = ttest(data(Lcor,1)-data(Lcor,2))
p = signrank(data(Lcor,1),data(Lcor,2))
title(['p = ',num2str(p)]);
disp(['median r(plane,true) r(quad, true): ',num2str(median(data(Lcor,[1 2])))]);

% Considering only one shape type at a time
for i = 1:3
    Lshape = data(:,7) == i;
    p = signrank(data(Lcor&Lshape,1),data(Lcor&Lshape,2));
    disp(['median r(plane,true) r(quad, true) p: ',num2str(median(data(Lcor&Lshape,[1 2]))),'  ', num2str(p)]);
end
Lshape = zeros(size(data,1),1);
% 1 = 1 sheet, 2 = 2 sheets, 3 = ellipsoid

% Both hyperboloids together
Lshape = data(:,7) < 3;
p = signrank(data(Lcor&Lshape,1),data(Lcor&Lshape,2))
median(data(Lcor&Lshape,[1 2]))


subplot(2,2,2);
hist(data(Lcor,1)-data(Lcor,2),20);
subplot(2,2,3);
hist(data(Lcor,1),[-1:.1:1]);
set(gca,'XLim',[0 1]);
subplot(2,2,4);
hist(data(Lcor,2),[-1:.1:1]);
set(gca,'XLim',[0 1]);

% Cone weight analysis
figure; axes; hold on;
coneweights = data(:,[4 5 6])./repmat(sum(abs(data(:,[4 5 6])),2),1,3);
colors = {'black','red','blue'};
for i = 1:3
    L = data(:,7) == i;
    h(i) = plot(coneweights(L,1),coneweights(L,2),'ko','MarkerFaceColor',colors{i});
    plot(-coneweights(L,1),-coneweights(L,2),'ko','MarkerFaceColor',colors{i});
end
plot([0 1 0 -1 0],[1 0 -1 0 1],'k-'); axis square;
legend(h,{'Hyper 1','Hyper 2','Ellipsoid'});

% Analysis of MSE
figure; axes; hold on;
plot(data(Lcor&Lshape,8),data(Lcor&Lshape,9),'r.')
plot(data(Lcor&~Lshape,8),data(Lcor&~Lshape,9),'k.')
set(gca,'XScale','log')
set(gca,'YScale','log')
plot([1 10^5],[1 10^5],'k:');
axis equal;
p = signrank(data(Lcor&~Lshape,8)-data(Lcor&~Lshape,9));
title(num2str(p));
xlabel('plane MSE');
ylabel('quad MSE');
disp(['mean MSEplane MSEquad: ',num2str(mean(data(Lcor,[8 9])))]);

% How mny cells (of which type) had linear and quadratic predictions that
% were anticorrelated with the data?
L = data(:,1) < 0 & data(:,2) < 0
data(L,7)
%%
% Section 2: Looking at superimposed isoresponse surfaces from a few
% different cells.  

USE10DEG = 1;
data = [];
hs = [];
[fnames, spikeIdx] = fnamesFromTxt2();
USE10DEG = 1;
if (USE10DEG)
    load ('T_cones_smj10');
end
figure(1); axes; hold on;
figure(2); axes; hold on;
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
    out = NTpreprocess(NT,.4,2);
    
    scaled = out(:,[2:4]).*repmat(out(:,5), 1,3);
    if (USE10DEG)
        fundamentals = NT.sum.exptParams.fundamentals;
        fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
        mon_spd = NT.sum.exptParams.mon_spd;
        mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
        mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
        M = fundamentals'*mon_spd;
        M10 = T_cones_smj10*mon_spd;
        scaled = ConvertConeContrastBasis(M, M10, NT.sum.exptParams.bkgndrgb, scaled);
    end
    Loog = logical(out(:,7));
    
    [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);
    data(cellcounter).planeparams = planeparams;
    data(cellcounter).quadparams = quadparams;
    data(cellcounter).xformM = xformmat;
    data(cellcounter).xyz = scaled;
    data(cellcounter).Loog = Loog;
    
    m = max(max(abs(scaled(~Loog,:))));
    domain = linspace(-m,m,30)/2;
    [x,y,z] = meshgrid(domain,domain,domain);
    xyz = reshape([x(:) y(:) z(:)]*xformmat, [size(x,1),size(x,2),size(x,3),3]);
    x = xyz(:,:,:,1);
    y = xyz(:,:,:,2);
    z = xyz(:,:,:,3);
    variables = [x(:).^2 y(:).^2 z(:).^2 2*x(:).*y(:) 2*x(:).*z(:) 2*y(:).*z(:)];
    fr = variables*quadparams;
    figure(1);
    p = patch(isosurface(x,y,z,reshape(fr,size(x)), 1));
    set(p,'FaceAlpha',.2,'FaceColor',unifrnd(0,1,1,3),'Edgealpha',0);
    axis vis3d;
    camlight;
    lighting phong;
    axis square;
    xlabel('L'); ylabel('M'); zlabel('S');
    hs = [hs; p];
    
    figure(2);
    if (cellcounter == 1)
        plot([0 1 0 -1 0],[1 0 -1 0 1],'k-');
    end
    coneweights = (planeparams'*xformmat')./sum(abs(planeparams'*xformmat'));
    [cellcounter coneweights]

    h = plot(coneweights(1),coneweights(2), 'ko','MarkerSize',10);
    if (coneweights(3) > 0)
        set(h,'MarkerFaceColor','black');  % S going with L/M
    end
    set(h,'ButtonDownfcn',['set(hs(',num2str(cellcounter),'),''Visible'',''on'')']);
    h = plot(-coneweights(1)./sum(abs(coneweights)),-coneweights(2)./sum(abs(coneweights)),'ko','MarkerSize',10);
    if (coneweights(3) < 0)
        set(h,'MarkerFaceColor','black');
    end
    set(h,'ButtonDownfcn',['set(hs(',num2str(cellcounter),'),''Visible'',''on'')']);
end
set(hs,'Visible','off');
figure(1);
plot3(0,0,0,'y*');
uih = uicontrol('Parent', figure(1), 'Style', 'PushButton', 'Callback','set(hs,''Visible'',''off'');','String','Clear');


%%
% Section 3
% Bootstrapping tests.
% Also calculating things for a planarity index, surface types, and "cone
% weights" and a conventional F test.
% 1/1/12 GDLH adding a RESTRICT_S flag. When set to '1' it artificially
% restricts the S-cone axis as per the request of reviewer #4.

RESTRICT_S = 0;

nbootiter = 1000;
ps = [];
planeSSEs = [];
quadSSEs = [];
totSSEs = [];
surfacetypes = [];
coneweights = [];
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
    %out = NTpreprocess(NT,.4,2);
    out = NTpreprocess(NT,0,Inf);
    
    scaled = out(:,[2:4]).*repmat(out(:,5), 1,3);
    Loog = logical(out(:,7));
    
    if (RESTRICT_S)
        L = (abs(scaled(:,3)).*~Loog) > 0.2; % Not going out beyond 20% S-cone contrast
        Loog(L) = 1;
    end
    
    if (sum(~Loog) <= 3)
        disp('Too few points to fit anything');
        ps = [ps; nan nan];
        surfacetypes = [surfacetypes; 4];
       continue
    end
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
    p = sum(nulldist>planeSSE/quadSSE)./nbootiter;
    planeSSEs = [planeSSEs; planeSSE mean(SSEs(:,1)) std(SSEs(:,1))];
    quadSSEs = [quadSSEs; quadSSE mean(SSEs(:,2)) std(SSEs(:,2))];
    totSSEs = [totSSEs; sum(r(~Loog).^2)];
    
    % Conventional F-test
    F = ((planeSSE-quadSSE)/3)/(quadSSE/(sum(~Loog)-6))
    pf = 1-fcdf(F,3,sum(~Loog)-3);
    ps = [ps; p pf]
    surfacetypes = [surfacetypes; WHICHSURF];
    coneweights = [coneweights; (planeparams'*xformmat')./sum(abs(planeparams'*xformmat'));];
end
ps(isnan(ps)) = 1.1;

%planeidx = (planeSSEs(:,1)-quadSSEs(:,1))./totSSEs;
%corr([planeidx ps]) % Not great agreement between plane idx and p-value
alpha = 0.01;
% Which quadric shape is most commonly associated with planes? 
% Hyperboloids of 2 sheets
Lplane = ps(:,2) > alpha;
types = zeros(2,3);
for i = 0:2
    types(1,i+1) = sum(surfacetypes == i & Lplane);
    types(2,i+1) = sum(surfacetypes == i & ~Lplane);
end
figure;
bar(types')
set(gca,'Xticklabel',{'Ellipsoid','1 sheet','2 sheets'});
legend({['p > ',num2str(alpha)],['p < ',num2str(alpha)]})
whichsurface = surfacetypes+1;
whichsurface(Lplane == 1) = 4;
shapenames = {'Ellipsoid','Hyperboloid 1','Hyperboloid 2','Plane'};

% Piechart showing a breakdown of the categories of neurons
piedata = [sum(whichsurface == 1) sum(whichsurface == 2) sum(whichsurface == 3) sum(whichsurface == 4)];
figure;
pie(piedata);
legend(shapenames);

% "Coneweights" for different categories of neuron
% Right now coneweights are calulated from plane fit,
% but really should be on the basis of eigenvectors
% for nonplanar cells.
cmap = colormap;
symbolcolors = cmap(floor(linspace(1,length(cmap),4)),:)
figure;
for j = 1:4
    subplot(2,2,j); hold on;
    plot([-1 1 0 -1],[0 0 1 0],'k-');
    for i = find(whichsurface == j)'
        cw = coneweights(i,:);
        cw = cw*sign(cw(2));  % Positive M convention
        h = plot(cw(1),cw(2),'ko','MarkerEdgeColor',symbolcolors(j,:),'MarkerFaceColor',symbolcolors(j,:));
        if (sign(cw(3)) == -1)
            set(h,'MarkerFaceColor','white');
        end
    end
    title(shapenames{j});
end

% list of planar cells
cellfun(@(x) char(x(2)), fnames(Lplane),'UniformOutput',0)

nsig = sum(ps<0.01);
(nsig(1)-nsig(2))/length(ps)  % difference in the proportion of sigs
% Only 3% of cells go from sig to non sig, and conventional F-test appears
% to be conservative
sum(all(ps<0.01,2))  % all of the significant cells but 2 are significant by both tests
L = surfacetypes == 2;
[sum(L) sum(ps(L,2)<0.01)]

%%
% Section 3.1 
% Pseudo R-squared goodnes-to-fit statistics
% S041310002 - planar example cell from manuscript
PR2_quad = [];
PR2_plane = [];
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
    %out = NTpreprocess(NT,.4,2);
    out = NTpreprocess(NT,0,Inf);
    
    scaled = out(:,[2:4]).*repmat(out(:,5), 1,3);
    Loog = logical(out(:,7));
    [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);
    [th,ph,r] = cart2sph(scaled(:,1),scaled(:,2),scaled(:,3));
    SST = sum((log(r)-mean(log(r))).^2);
    PR2_quad(cellcounter) = 1-quadSSE/SST
    PR2_plane(cellcounter) = 1-planeSSE/SST
end

%%
% Section 4: Comparing plane fits for luminance cells using a color space based on 2°
% fundamentals and a space based on 10° fundamentals.  NTlum is good
% filelist.

load ('T_cones_smj10');
data = [];
hs = [];
[fnames, spikeIdx] = fnamesFromTxt2;
for cellcounter = 1:size(fnames,1)
    NT = {};
    for i = 1:size(fnames{cellcounter},2)
        filename = findfile(char(fnames{cellcounter}(i)));
        paradigmID = getparadigmID(filename);
        if paradigmID == 103
            NT = nex2stro(filename);
        end
    end
    out = NTpreprocess(NT,.4,2);
    fundamentals = NT.sum.exptParams.fundamentals;
    fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
    mon_spd = NT.sum.exptParams.mon_spd;
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
    M = fundamentals'*mon_spd;
    M10 = T_cones_smj10*mon_spd;
    Loog = logical(out(:,7));
    ecc = sqrt((NT.sum.exptParams.rf_x/10).^2 + (NT.sum.exptParams.rf_y/10).^2);
    figure; 
    for i = 1:2
        scaled = out(:,[2:4]).*repmat(out(:,5), 1,3);
        if (i == 2)  % converting scaled to 10 degree fundamentals
            scaled = ConvertConeContrastBasis(M, M10, NT.sum.exptParams.bkgndrgb, scaled);
        end
        [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);
        coneweights = (xformmat*planeparams)';
        if (sign(coneweights(1)) == -1)
            coneweights = -coneweights;
        end
        data(cellcounter,:,i) = coneweights./sum(abs(coneweights));
        subplot(2,1,i); hold on;
        plot3(scaled(~Loog,1), scaled(~Loog,2), scaled(~Loog,3),'k.');
        plot3(-scaled(~Loog,1), -scaled(~Loog,2), -scaled(~Loog,3),'k.');
        plot3([zeros(sum(Loog),1) scaled(Loog,1)]', [zeros(sum(Loog),1) scaled(Loog,2)]', [zeros(sum(Loog),1) scaled(Loog,3)]','y-');
        plot3([zeros(sum(Loog),1) -scaled(Loog,1)]', [zeros(sum(Loog),1) -scaled(Loog,2)]', [zeros(sum(Loog),1) -scaled(Loog,3)]','y-');
        
        signflips = sign(scaled(~Loog,:)*coneweights');
        scaled(~Loog,:) = repmat(signflips,1,3).*scaled(~Loog,:);
        for j = 0:1
            if (j == 1)
                scaled(~Loog,:) = -scaled(~Loog,:);
            end
            if (sum(~Loog) > 3)
                T = DelaunayTri(scaled(~Loog,1),scaled(~Loog,2),scaled(~Loog,3));
                [tri, Xb] = convexHull(T);
                htm = trisurf(tri,scaled(~Loog,1),scaled(~Loog,2),scaled(~Loog,3),'FaceColor', 'cyan', 'FaceAlpha', 0.8, 'EdgeAlpha',.01);
            end
        end
        axis square
        set(gca,'ZLim',[-.7 .7],'XLim',[-.2 .2],'Ylim',[-.2 .2],'View',[-137 0]);
        slashidxs = strfind(filename,'\');
        title(filename(slashidxs(end)+1:end));
    end
    title(num2str(ecc));
end

figure;
for i = 1:2
    subplot(3,2,2*(i-1)+1); hold on;
    plot([-1 0 1],[0 1 0],'k-');
    L = logical(sign(data(:,3,i)) == 1);
    plot(squeeze(data(L,1,i)), squeeze(data(L,2,i)),'ko','MarkerSize',7,'MarkerFaceColor',[0 0 0]);
    plot(squeeze(data(~L,1,i)), squeeze(data(~L,2,i)),'ko','MarkerSize',7,'MarkerFaceColor',[1 1 1]);
    set(gca,'XLim',[0 1],'Ylim',[0 1]);
    axis square;
    subplot(3,2,2*(i-1)+2); hold on;
    hist(data(:,3,i));
    [h,p] = ttest(data(:,3,i));
    title(['p = ',num2str(p)]);
end
subplot(3,2,5); hold on;
plot(data(:,3,1),data(:,3,2),'k.');
plot([min(min(data(:,3,:))) max(max(data(:,3,:)))], [min(min(data(:,3,:))) max(max(data(:,3,:)))],'y-');
set(gca,'XLim',[min(min(data(:,3,:))) max(max(data(:,3,:)))]);
set(gca,'YLim',[min(min(data(:,3,:))) max(max(data(:,3,:)))]);
xlabel('S-cone weight with 2 deg');
ylabel('S-cone weight with 10 deg');

axis square
mean(data)
% With the 2 deg fundamentals, there's significant S-cone input with the same sign as L and M 
% With the 10 deg fundamentals the S-cone input gets smaller.

%% 
% Section 4.1
% How much macular pigment do we have use to get rid of the S-cone
% contribution to luminance cells?  Does this also straighten up the
% ellipsoids?  That would be cool independent confirmation for a good set
% of cone fundamentals for monkeys.
%
% 3/9/12 GDLH Overhaul. Getting rid of outdated code (adjusting
% macular pigment of the 10 deg (or 2 deg) fundamentals to create the other
% set. This was done before I knew that the 2 deg funds were based on the
% 10 deg funds and also included a change in photopigment density.)

load('den_mac_vos.mat'); 
load('den_mac_bone.mat');  den_mac_bone = SplineRaw([390  1  441],den_mac_bone,[380 4 81]); % sucks
load('den_mac_ws.mat'); 
load('T_cones_smj.mat'); T_cones_smj = T_cones_smj';
load('T_cones_smj10.mat'); T_cones_smj10 = T_cones_smj10';
load('T_cones_ss2.mat'); T_cones_ss2 = SplineRaw([390  1  441],T_cones_ss2',[380 4 81]);
load('T_cones_ss10.mat'); T_cones_ss10 = SplineRaw([390  1  441],T_cones_ss10',[380 4 81]);
load('T_cones_synthgh1'); % My new synthetic set of monkey cone fundamentals

% Here's where you can change the fundamentals used
whichfundamentals = T_cones_synthgh1; 

macpig = den_mac_ws;
macpig(4) = 0.0425; %  minor tweak made by SMJ 1993: density of 0.0425 @ 395 nm
macpigpeaks = linspace(0,.5,10); % peak densities
% Computing cone weights for a bunch of cells using a few different
% cone fundamentals.

if all(all(whichfundamentals == T_cones_smj10)) % Removing assumed macular pigment
    macpigtransmittance = 1./(10.^(macpig*.28));
    whichfundamentals = whichfundamentals./repmat(macpigtransmittance,1,3);
    whichfundamentals = whichfundamentals./repmat(max(whichfundamentals),81,1);
end

FILELISTIDX = 2;
filelists = {'NTlum.txt','pancolor2.txt'};
[fnames, spikeIdx] = fnamesFromTxt2(['N:\NexFiles\nexfilelists\Greg\NT\',filelists{FILELISTIDX}]);
data = nan*ones(size(fnames,1),3,length(macpigpeaks));
eccentricities = zeros(size(fnames,1),1);
for cellcounter = 1:size(fnames,1)
    NT = {};
    for i = 1:size(fnames{cellcounter},2)
        filename = findfile(char(fnames{cellcounter}(i)));
        paradigmID = getparadigmID(filename);
        if paradigmID == 103
            NT = nex2stro(filename);
        end
    end
    out = NTpreprocess(NT,0,Inf);
    Loog = logical(out(:,7));
    scaled = out(:,[2:4]).*repmat(out(:,5), 1,3);
    
    fundamentals = NT.sum.exptParams.fundamentals;
    fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
    mon_spd = NT.sum.exptParams.mon_spd;
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
    M = fundamentals'*mon_spd;
    bkgndlms = M*NT.sum.exptParams.bkgndrgb;
    
    for i = 1:length(macpigpeaks)
        i
        macpigtransmittance = 1./(10.^(macpig./max(macpig)*macpigpeaks(i)));
        currentfunds = whichfundamentals.*repmat(macpigtransmittance,1,3);
        currentfunds = currentfunds./repmat(max(currentfunds),81,1);
        currentM = currentfunds'*mon_spd;
        currentscaled = ConvertConeContrastBasis(M, currentM, NT.sum.exptParams.bkgndrgb, scaled);
        
        [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(currentscaled, Loog);
        
        if (strcmp(filelists{FILELISTIDX},'NTlum.txt')) % lum planes
            coneweights = (xformmat*planeparams)';
            if (coneweights(1)+coneweights(2) < 0)
                coneweights = -coneweights;
            end
        else  % ellipsoids
            A = [quadparams(1) quadparams(4) quadparams(5);...
                quadparams(4) quadparams(2) quadparams(6);...
                quadparams(5) quadparams(6) quadparams(3)];
            B = xformmat*A*xformmat';
            [evecs, evals] = eig(B);
            [evals,idxs] = sort(diag(evals));
            evecs = evecs(:,idxs);
            coneweights = evecs(:,1);  % Small eigenvalue = long axis
            if (coneweights(3) < 0)  % Positive S-cone weights
                coneweights = -coneweights;
            end
        end
        data(cellcounter,:,i) = coneweights./sum(abs(coneweights));
        eccentricities(cellcounter) = sqrt(NT.sum.exptParams.rf_x.^2+NT.sum.exptParams.rf_y.^2)/10;
    end
end
if (strcmp(filelists{FILELISTIDX},'pancolor.txt')) 
    data(data(:,3,end) < .8,:,:) = []  % 
end
% Cone weight diagram   
figure;
for i = 1:size(data,3)
    subplot(ceil(sqrt(size(data,3))),ceil(sqrt(size(data,3))),i); hold on;
    for j = 1:size(data,1)
        h = plot(data(j,1,i),data(j,2,i),'ko');
        if sign(data(j,3,i)) == 1
           set(h,'MarkerFaceColor','black'); 
        end
    end
    plot([0 1],[1 0],'k-');
end

% S-cone weight histograms
figure;
bins = linspace(min(min(squeeze(data(:,3,:)))),max(max(squeeze(data(:,3,:)))),10);
p = [];
for i = 1:size(data,3)
    subplot(ceil(sqrt(size(data,3))),ceil(sqrt(size(data,3))),i); hold on;
    hist(data(:,3,i),bins)
    [h,p(i)] = ttest(data(:,3,i));
    title(['mn = ',num2str(mean(data(:,3,i)),2),' p = ',num2str(p(i),2)]);
end

% Line plot
figure; axes; hold on;
for i = 1:size(data,3)
    plot(macpigpeaks(i), mean(data(:,3,i)),'o');
    plot(macpigpeaks(i)*[1 1], mean(data(:,3,i))+std(data(:,3,i))./sqrt(size(data,3)).*[-1 1],'k-')
end
xlabel('Peak macular pigment density');
ylabel('S-cone weight');


% % ------------ OUTDATED STUFF -----------------
% % Finding a macular pigment density (to take away) to make fund2
% % as similar as possible to fund10.  This could be done with linear
% % regression I think.
%
% ADDMACPIG = 1;
% figure; subplot(3,1,1); hold on;
% error = [];
% if (ADDMACPIG) % We add macular pigment to the 10 degree fundamentals
%      %  (as opposed to taking it away from the 2 degree ones).
%     macpigweights = linspace(.55,-.3,20);
%     for macpigweight = macpigweights
%         macpigtransmittance = 1./(10.^(macpig*macpigweight));
%         synthfund = fund10.*repmat(macpigtransmittance,1,3);
%         % Starting with lots macular pigment and taking it away
%         synthfund = synthfund./repmat(max(synthfund),size(synthfund,1),1);
%         plot(synthfund-fund2);
%         error = [error; sum(sum((synthfund(:,3)-fund2(:,3)).^2))];
%         % Starting with lots macular pigment and taking it away
%         synthfunds(:,:,macpigweight==macpigweights) = synthfund./repmat(max(synthfund),size(synthfund,1),1);
%     end
% else % We start with the 2 degree fundamentals and take macular pigment away
%     macpigweights = linspace(-.55,-.15,10);
%     for macpigweight = macpigweights
%         macpigtransmittance = 1./(10.^(macpig*macpigweight));
%         synthfund = fund2.*repmat(macpigtransmittance,1,3);
%         % Starting with lots macular pigment and taking it away
%         synthfund = synthfund./repmat(max(synthfund),size(synthfund,1),1);
%         plot(synthfund-fund10);
%         error = [error; sum(sum((synthfund(:,3)-fund10(:,3)).^2))];
%         % Starting with lots macular pigment and taking it away
%         synthfunds(:,:,macpigweight==macpigweights) = synthfund./repmat(max(synthfund),size(synthfund,1),1);
%     end
% end
% synthfund = synthfunds(:,:,error == min(error));
% subplot(3,1,2);
% plot(macpigweights,error,'b.-');  % -.36 is pretty good for ~ADDMACPIG (with den_mac_ws, T_cones_smj)
% subplot(3,1,3);hold on;
% plot(fund2,'b-')
% plot(fund10,'g-')
% plot(synthfund,'r-') % Red is supposed to be similar to green (or blue) and it is.
% % ------------ OUTDATED STUFF -----------------



%%
% Section 5
% Looking at isoresponse ellipsoids for a bunch of pan-color cells
% and seeing whether there is any consistency in their orientation.
% Option to use the 10 degree fundamentals.

eigenvectors = [];
eigenvalues = [];
USE10DEG = 1;
COLOREDSYMBOLS = 1;

if (USE10DEG)
    load ('T_cones_smj10');
end
[fnames, spikeIdx] = fnamesFromTxt2;
for cellcounter = 1:size(fnames,1)
    NT = {};
    for i = 1:size(fnames{cellcounter},2)
        filename = findfile(char(fnames{cellcounter}(i)));
        paradigmID = getparadigmID(filename);
        if paradigmID == 103
            NT = nex2stro(filename);
        end
    end
    out = NTpreprocess(NT,.2,1);
    scaled = out(:,[2:4]) .*repmat(out(:,5), 1,3);
    if (USE10DEG)
        fundamentals = NT.sum.exptParams.fundamentals;
        fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
        mon_spd = NT.sum.exptParams.mon_spd;
        mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
        mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
        M = fundamentals'*mon_spd;
        M10 = T_cones_smj10*mon_spd;
        scaled = ConvertConeContrastBasis(M, M10, NT.sum.exptParams.bkgndrgb, scaled);
    end
    
    Loog = logical(out(:,end));
    if (any(Loog))
        disp(fnames{cellcounter});
        continue;
    end
    figure; axes; hold on;
    plot3(scaled(:,1),scaled(:,2),scaled(:,3),'r.' );
    plot3(-scaled(:,1),-scaled(:,2),-scaled(:,3),'r.' );
    [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);
    A = [quadparams(1) quadparams(4) quadparams(5);...
         quadparams(4) quadparams(2) quadparams(6);...
         quadparams(5) quadparams(6) quadparams(3)];
    % Transforming quadparams into un-whitened space
    [evecs,evals] = eig(xformmat*A*xformmat')
    [evals,i] = sort(diag(evals),1,'ascend');
    evecs = evecs(:,i);
    radii = sqrt(1./evals);
    [x, y, z] = ellipsoid(0,0,0,radii(1),radii(2), radii(3));
    [newxyz] = [x(:) y(:) z(:)]*evecs';
    h = surf(reshape(newxyz(:,1),size(x)), reshape(newxyz(:,2),size(y)), reshape(newxyz(:,3),size(z)));
    set(h,'FaceAlpha',.2,'FaceColor','green','Edgealpha',0);
    axis vis3d;
    camlight;
    lighting phong;
    
    eigenvectors = cat(3,eigenvectors,evecs);
    eigenvalues = cat(2,eigenvalues,evals);
end

radii = [];
for i = 1:size(eigenvalues,2)
    radii(:,i) = sqrt(1./eigenvalues(:,i));
end

% Don't plot symbols for axes that are too similar to 
% other axes
radiusratiothresh = 2;
Lrad = logical([radii(2,:)./radii(3,:) > radiusratiothresh;...
radii(1,:)./radii(2,:) > radiusratiothresh]);


symbols = {'gs','rd','bv'};
figure; axes; hold on;
plot([-1 0 1],[0 1 0],'k-');
set(gca,'Xlim',[-1 1],'Ylim',[0 1])
for ax = 1:3
    if (ax == 1)
        L = Lrad(1,:);
    elseif (ax == 2)
        L = Lrad(1,:) & Lrad(2,:);        
    else (ax == 3)
        L = Lrad(2,:);
    end
    tmp = squeeze(eigenvectors(:,ax,L));
    tmp = tmp./repmat(sum(abs(tmp)),3,1);
    tmp = tmp.*repmat(sign(tmp(2,:)),3,1);
    if (COLOREDSYMBOLS)
        h = plot(tmp(1,:),tmp(2,:),symbols{ax});
        set(h,'MarkerFaceColor',get(h,'Color'),'MarkerSize',10);
        %L = sign(tmp(:,3)) == 1;
    else
        h = text(tmp(1,:),tmp(2,:),num2str(ax));
        L = sign(tmp(:,3)) == 1;
        set(h(L),'Color','red');
    end
end

%%
% Section 6.  Looking for a good index of planarity 
% (and asking whether planar isoresponse surfaces tend to be oriented in a 
% consistent way).  Bottom line: there are a lot of planar lum cells and a
% few planar L-M cells too.

data = [];
[fnames, spikeIdx] = fnamesFromTxt2();
for cellcounter = 1:size(fnames,1)
    NT = {};
    for i = 1:size(fnames{cellcounter},2)
        filename = findfile(char(fnames{cellcounter}(i)));
        paradigmID = getparadigmID(filename);
        if paradigmID == 103
            NT = nex2stro(filename);
        end
    end
    out = NTpreprocess(NT,0,Inf);
    Loog = logical(out(:,7));
    scaled = out(:,[2:4]) .*repmat(out(:,5), 1,3);
    r = sqrt(sum(scaled.^2,2));
    [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);

    totSSE1 = sum(r(~Loog).^2); % not including OOGs
    totSSE2 = sum(r.^2); % Including OOGs
    data = [data;(planeparams'*xformmat') planeSSE quadSSE totSSE1 totSSE2 sum(Loog) sum(~Loog)];
    data
end

L = data(:,8) > 2;  % At least 2 OOG = not a closed surface
L = L&data(:,8)+data(:,9) > 10 % at least n points survive filtering
coneweights = data(L,[1 2 3])./repmat(sum(abs(data(L,[1 2 3])),2),1,3);
%stat = data(L,5)./data(L,4);  % quad error/plane error
%stat = data(L,4)./data(L,6);   % plane error/total error
%stat = data(L,5)./data(L,6);  % quad error/total error
%stat = data(L,8) + data(L,9); % Total number of points

figure;  % Trying to come up with a reasonable planarity statistic that doesn't depend too much on 'n'
stat = (data(L,4)-data(L,5))./data(L,6);  % normalized quad error/ normalized plane error
plot(data(L,8) ,stat,'k.'); xlabel('N in gamut');ylabel('differences in normalized SSE'); 
set(gca,'YScale','log');
% When this statistic is near zero, that means that quadSSE~planeSSE ie.
% plane and quad fit about the same.  Doesn't seem to depend much on 'n'
fnames{find(stat < 0.05)} % Nothing special about 0.05;

% Any relationship between this statistic and L+M cone weights?
figure; subplot(2,1,1); hold on;
plot([0 1 0 -1 0],[1 0 -1 0 1],'k-');
for i = 1:size(coneweights,1)
    plot(coneweights(i,1),coneweights(i,2),'ko','markersize',[10*sqrt(stat(i)/max(stat))]);
    plot(-coneweights(i,1),-coneweights(i,2),'ko','markersize',[10*sqrt(stat(i)/max(stat))]);
end
axis square;
subplot(2,1,2);
plot(abs(coneweights(:,1)+coneweights(:,2)),stat,'k.')
[rho,p] = corr([abs(sum(coneweights(:,[1 2]),2)),stat],'Type','Spearman');
title(['r = ',num2str(rho(1,2)),'  p = ',num2str(p(1,2))])
xlabel('L+M'); ylabel('plane fit statistic');
set(gca,'YScale','log');

%%
% Section 7.  What are the "preferred color directions" of the pan color
% cells estimated from the gratings data?  Ans: All over the map.
data = [];
[NTfnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\pancolor.txt');
[GTfnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\Gratings&NT.txt');
for cellcounter = 1:size(NTfnames,1)
    for i = 1:size(GTfnames,1)
        if (any(strcmp(NTfnames{cellcounter}, GTfnames{i})))
            GT = {};
            for j = 1:size(GTfnames{cellcounter},2)
                filename = findfile(char(GTfnames{i}(j)));
                paradigmID = getparadigmID(filename);
                if paradigmID == 150
                    GT = nex2stro(filename);
                    GTstruct = getGratingTuning(GT,1);
                    data = [data; GTstruct.color.prefcolor]
                end
            end
        end
    end
end
coneweights = data./repmat(sum(abs(data),2),1,3);
Lsplus = sign(coneweights(:,3)) == 1;
figure; axes; hold on;
plot(coneweights(Lsplus,1), coneweights(Lsplus,2),'ko');
plot(-coneweights(~Lsplus,1), coneweights(~Lsplus,2),'k.');
plot(coneweights(~Lsplus,1), coneweights(~Lsplus,2),'ko');
plot(-coneweights(Lsplus,1), -coneweights(Lsplus,2),'k.');
plot([-1 0 1 0 -1],[0 1 0 -1 0],'k-');


%%
% Section 8
% Comparing estimates of "preferred color direction" from grating data and
% NeuroThresh data.  Are cells with the greatest disagreements the ones
% whose isoresponse surfaces are the least planar?

data = [];
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
    out = NTpreprocess(NT,.4,2);
        
    scaled = out(:,[2:4]).*repmat(out(:,5), 1,3);
    Loog = logical(out(:,7));
 
    [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);
    A = [quadparams(1) quadparams(4) quadparams(5);...
         quadparams(4) quadparams(2) quadparams(6);...
         quadparams(5) quadparams(6) quadparams(3)];
    if (all(eig(A) > 0))
        ELLIPSOID = 1;
    else
        ELLIPSOID = 0;
    end
    GTconeweights = GTstruct.color.prefcolor./sum(abs(GTstruct.color.prefcolor));
    if (GTconeweights(2) < 0)
        GTconeweights = -GTconeweights;
    end
    NTconeweights = (planeparams'*xformmat')./sum(abs(planeparams'*xformmat'));
    if (NTconeweights *  GTconeweights' < 0)
        NTconeweights = -NTconeweights;
    end
    data = [data; GTconeweights NTconeweights planeSSE quadSSE ELLIPSOID]
end

figure; axes; hold on;
L = logical(data(:,end)); % Ellipsoid fit
plot(data(:,1),data(:,2),'k.');
plot(data(~L,4),data(~L,5),'k+');
plot(data(L,1),data(L,2),'ro');
plot(data(L,4),data(L,5),'ro');
plot(data(:,[1 4])',data(:,[2 5])','k-');
plot([0 1 0 -1 0],[1 0 -1 0 1],'k-');

figure; axes; hold on;
subplot(2,1,1); hold on;
dist = sqrt(sum((data(:,[1 2 3])-data(:,[4 5 6])).^2,2)); % distance in normalized cone weights space. yuck
bins = linspace(0,1,20);
[n,x] = hist(dist,bins);
bar(x,n,'k');
[n,x] = hist(dist(L),bins);
bar(x,n,'r');
[h,p] = ttest2(dist(~L),dist(L));
title(['p = ',num2str(p)]);
subplot(2,1,2); hold on;
plot(dist(~L),data(~L,7)./data(~L,8),'k.');
plot(dist(L),data(L,7)./data(L,8),'r.');
set(gca,'Yscale','log');
%[r,p] = corr([dist data(:,7)./data(:,8)],'type','Spearman')
[r,p] = corr([dist(~L) data(~L,7)./data(~L,8)],'type','Spearman')
title(['r = ',num2str(r(1,2)),' p = ',num2str(p(1,2))]);
xlabel('distance in cone space');
ylabel('planeSSE/quadSSE');
% Not surprisingly, pan color cells lead to bigger disparities in cone
% weight estimates than the other cells.
% Weak but correlation between crappiness of plane fit and
% disparity in cone weight estimates (not including pan color cells).

% Plotting "cone weights" (from grating data top and NT data bottom)
% in Lennie diagram. Different symbols for pan color and other cells.

ms = 1.5; % Markersize
figure;
titles = {'All','Ellipsoids in red','Non-ellipsoids'};
for j = 1:3
    for i = 1:2
        subplot(2,3,j+3*(i-1)); hold on;
        lcw = 3*i-2;
        mcw = 3*i-1;
        L = logical(data(:,end)); % Ellipsoid fit
        plot(data(~L,lcw),data(~L,mcw),'ko','MarkerFaceColor','k','MarkerSize',ms);
        plot(-data(~L,lcw),-data(~L,mcw),'ko','MarkerFaceColor','k','MarkerSize',ms);
        if (j == 1)  % Plot OOGs in black
            plot(data(L,lcw),data(L,mcw),'ko','MarkerFaceColor','k','MarkerSize',ms);
            plot(-data(L,lcw),-data(L,mcw),'ko','MarkerFaceColor','k','MarkerSize',ms);
        elseif (j == 2) % plot OOGs in red
            plot(data(L,lcw),data(L,mcw),'ro','MarkerFaceColor','r','MarkerSize',ms);
            plot(-data(L,lcw),-data(L,mcw),'ro','MarkerFaceColor','r','MarkerSize',ms);
        end  % else don't plot OOgs
        plot([0 1 0 -1 0],[1 0 -1 0 1],'k-');
        axis square;
        if (i == 1)
            title(titles{j});
        end
    end
end

% Standard estimates of cone weights
cw = data(:,[1 2 3]);
cw(cw(:,2) < 0,:) = -cw(cw(:,2) < 0,:);
figure; axes; hold on
L = cw(:,3) > 0;
plot(cw(L,1),cw(L,2),'ko','MarkerFaceColor','black','MarkerSize',5);
plot(cw(~L,1),cw(~L,2),'ko','MarkerSize',5);
plot([-1 0 1 -1],[0 1 0 0],'k-');
%%
% Section 8.1
% Is the relative L- to M-cone weight of luminance tuned cells
% consistent among three estimates: 
%    raw L to raw M response (a la Johnson and Shapley)
%    linear model fits (a la Lennie)
%    neurothresh plane orientation


data = [];
[fnames1, spikeIdx1] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\NTlum.txt');
[fnames2, spikeIdx2] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\Gratings&NT.txt');

for cellcounter = 1:size(fnames1,1)
    filename = findfile(char(fnames1{cellcounter}));
    NT = nex2stro(filename);
    % Finding gratings data
    
    GT = {};
    for i = 1:size(fnames2,1)
        for j = 1:size(fnames2{i},2)
            if (strcmp(fnames2{i}(j),char(fnames1{cellcounter})))
                for k = 1:size(fnames2{i},2)
                    paradigmID = getparadigmID(findfile(fnames2{i}(k)));
                    if paradigmID == 150
                        GT = nex2stro(findfile(fnames2{i}(k)));
                    end
                end
                if (isempty(GT))
                    keyboard;
                end
            end
        end
    end
    if (isempty(GT))
        disp(['got here: ',num2str(cellcounter)]);
        continue;
    end
    
    GTstruct = getGratingTuning(GT,1);
    GTconeweights = GTstruct.color.prefcolor./sum(abs(GTstruct.color.prefcolor));
    if (GTconeweights(2) < 0)
        GTconeweights = -GTconeweights;
    end
    L = GTstruct.color.colors(:,1) & ~GTstruct.color.colors(:,2) & ~GTstruct.color.colors(:,3);
    GTrawLresp = GTstruct.color.colresp(L,1);
    L = ~GTstruct.color.colors(:,1) & GTstruct.color.colors(:,2) & ~GTstruct.color.colors(:,3);
    GTrawMresp = GTstruct.color.colresp(L,1);

    out = NTpreprocess(NT,.4,2);        
    scaled = out(:,[2:4]).*repmat(out(:,5), 1,3);
    Loog = logical(out(:,7));
    [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);
    NTconeweights = (planeparams'*xformmat')./sum(abs(planeparams'*xformmat'));
    if (NTconeweights *  GTconeweights' < 0)
        NTconeweights = -NTconeweights;
    end
    data = [data; GTrawLresp GTrawMresp GTconeweights(1) GTconeweights(2) NTconeweights(1) NTconeweights(2)];
end
relLweight = [log(data(:,1)./data(:,2)), log(data(:,3)./data(:,4)), log(data(:,5)./data(:,6))];
[r,p] = corr(relLweight);

figure;
subplot(2,2,1);
plot(relLweight(:,1), relLweight(:,3),'k.');
title(['r=',num2str(r(1,3),2),' p=',num2str(p(1,3),2)])
ylabel('rel. L-cone weight (NT)');
subplot(2,2,2);
plot(relLweight(:,2), relLweight(:,3),'k.');
title(['r=',num2str(r(2,3),2),' p=',num2str(p(2,3),2)])
xlabel('rel. L-cone weight (GT)');
subplot(2,2,3);
plot(relLweight(:,1), relLweight(:,2),'k.')
ylabel('rel. L-cone weight (GT)');
xlabel('rel. L-cone weight (GTraw)');
title(['r=',num2str(r(1,2),2),' p=',num2str(p(1,2),2)])


%%
% Section 9
% Comparing F1/F0 modulation ratios for pan-color and non-pan-color cells.
% Do the pan color cells tend to be complex cells?

data = [];
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
    out = NTpreprocess(NT,.4,2);
        
    scaled = out(:,[2:4]).*repmat(out(:,5), 1,3);
    Loog = logical(out(:,7));
    data = [data; sum(Loog) sum(~Loog) GTstruct.modulationratio]
end
Lpancolor = data(:,1) <=3;  % Pretty lame definition
bins = linspace(0,2,20);
[n1,x] = hist(data(~Lpancolor,3),bins);
[n2,x] = hist(data(Lpancolor,3),bins);
figure; axes; hold on;
bar(x,n1+n2,'FaceColor',[0 0 0]);
bar(x,n2,'FaceColor',[1 0 0]);
[h,p] = ttest2(data(Lpancolor,3),data(~Lpancolor,3))

%%
% Section 9.1
% Looking at F1/F0 modulation ratios on the basis of quadric fit type.
% Will have to deal with ellipsoids that extend outside of the monitor
% gamut and planes.

data = [];
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
  %  out = NTpreprocess(NT,.4,1);
    
    scaled = out(:,[2:4]).*repmat(out(:,5), 1,3);
    Loog = logical(out(:,7));
    [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);
    coneweights = (planeparams'*xformmat')./sum(abs(planeparams'*xformmat'));
    r2 = sum(scaled.^2,2);

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
    
    data = [data; WHICHSURF GTstruct.modulationratio GTstruct.sf.prefSF coneweights]
end
for i = 1:3
   L = logical(data(:,1) == i) & ~isnan(data(:,2));
   figure; 
   hist(data(L,2),linspace(0,2,20));
   geomean(data(L,2))
end
[p,~,a1] = anova1(data(:,2),data(:,1))

% Pooling hyperboloids
L = logical(data(:,1) < 3)
[nanmean(data(L,2)) nanmean(data(~L,2))]

% What's the deal with these ellipsoidal cells with high F1/F0s?
idxs = find(data(:,1) == 3 & data(:,2) > 1)
fnames{idxs}
[fnames{idxs}]'

% Stacked barplot of modulation ratios
USEPLANE = 1;
Lplane = logical(zeros(size(fnames,1),1));
if USEPLANE
[planefnames] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\plane0.01.txt');
for i = 1:size(fnames,1)
    Lplane(i) = ismember(fnames{i}(2),[planefnames{:}]');
end
end
nbins = 14;
bins = linspace(0,2,nbins);
n = zeros(3,nbins);
for whichsurf = 1:3
    L = data(:,1) == whichsurf;
    n(whichsurf,:) = hist(data(L & ~Lplane,2),bins);
end
n(size(n,1)+1,:) = hist(data(Lplane,2),bins);

figure;
bar(bins,n','stacked')
legend({'Hyp 1','Hyp 2','ellipsoid','plane'});
tmpdata = data;
tmpdata(Lplane,1) = 0;
anova1(tmpdata(:,2),tmpdata(:,1))

% Any planar, complex neurons with opponent weights?  Yes.
L = Lplane & data(:,1) == 3 & data(:,2)<1;
data(L, [4,5,6])
%%
% Section 10
% Are direction selective cells unusually planar?
% Yes, but only because pan color cells tend not to be direction
% selective (nor orientation-selective, hmm, this makes them look like
% the product of poor isolation).
data = [];
[fnames, spikeIdx] = fnamesFromTxt2();
planelist = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\plane0.01.txt');
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
    
    % Getting info for direction and orientation tuning indices
    theta = GTstruct.orient.stim(1:end-1);  % getting rid of 2*pi duplication
    resp = GTstruct.orient.resp(1:end-1,:);
    preforient = theta(resp(:,1) == max(resp(:,1)));
    prefresp = resp(resp(:,1) == max(resp(:,1)));
    nullorient = mod(preforient + pi,2*pi);
    nullresp = resp(softeq(nullorient, theta));
    orthorient = mod(preforient + [pi/2 -pi/2],2*pi);
    orthresp = mean([resp(softeq(orthorient(1), theta)), resp(softeq(orthorient(2), theta))]);
    inlineresp = mean([prefresp; nullresp])
    if (prefresp > 10)
        DS = (prefresp-nullresp)/(prefresp+nullresp);
        OS = (inlineresp-orthresp)/(inlineresp+orthresp);
    else
        DS = nan;
        OS = nan;
    end
    
    % Getting things for Johnson and Shapley definition of double-opponency
    if(~isfield(GTstruct.sf,'stim'))
        Lsf = 0;
    elseif any(GTstruct.sf.resp(2:end,1) >  GTstruct.sf.resp(1,1))
        Lsf = 1;
    else
        Lsf = 0;
    end
    if (~isfield(GTstruct,'color'))
        continue;
    end
    if (isnan(GTstruct.color.colors))
        continue;
    end
    dotprods = abs(GTstruct.color.colors*[1 -1 0]');
    rgidx = dotprods == max(dotprods) & GTstruct.color.colors(:,3) == 0;
    rg = GTstruct.color.colresp(rgidx,1);
    lumidx = all(GTstruct.color.colors == repmat(abs(GTstruct.color.colors(rgidx,:)),size(GTstruct.color.colors,1),1),2);
    lum = GTstruct.color.colresp(lumidx,1);
    DO = rg > lum & Lsf &  GTstruct.modulationratio>= 1;
    
    A = [quadparams(1) quadparams(4) quadparams(5);...
        quadparams(4) quadparams(2) quadparams(6);...
        quadparams(5) quadparams(6) quadparams(3)];
    d = eig(A);

    shortfn = NT.sum.fileName(find(NT.sum.fileName == '\',1,'last')+1:find(NT.sum.fileName == '.',1,'last')-1);
    Lplane = ismember(shortfn,[planelist{:}]');
    data = [data; DS OS DO sum(d<0) Lplane GTstruct.modulationratio NT.sum.exptParams.threshold]
end
surfacetype = data(:,4); % 1,2 = hyperboloids
surfacetype(surfacetype == 0) = 3;  % 3 = ellipsoids
surfacetype(logical(data(:,5))) = 0;  % 0 = planes

% Is there a relationship between direction selectivity and surface type?
anovan(data(:,1),surfacetype);
boxplot(data(:,1),surfacetype);
set(gca,'XTick',[1 2 3 4], 'XTickLabel',{'plane','1 sheet','2 sheets','ellipsoid'});
[h,p] = ttest2(data(surfacetype==0,1),data(surfacetype==3,1))  % planes vs ellipsoids
ylabel('Direction selectivity index');

% Is there a relationship between orientation selectivity and surface type?
anovan(data(:,2),surfacetype);
boxplot(data(:,2),surfacetype);
set(gca,'XTick',[1 2 3 4], 'XTickLabel',{'plane','1 sheet','2 sheets','ellipsoid'});
[h,p] = ttest2(data(surfacetype==0,1),data(surfacetype==3,1))  % planes vs ellipsoids
ylabel('Orientation selectivity index');

% Is there a relationship between modulation ratio and surface type?
anovan(data(:,6),surfacetype);
boxplot(data(:,6),surfacetype);
set(gca,'XTick',[1 2 3 4], 'XTickLabel',{'plane','1 sheet','2 sheets','ellipsoid'});
[h,p] = ttest2(data(surfacetype==0,6),data(surfacetype==3,6))  % planes vs ellipsoids
ylabel('F1/F0');

% Proportion of DO cells in each category
figure;
L_DO = logical(data(:,3));
[n_DO,x] = hist(surfacetype(L_DO),[0 1 2 3]);
[n_notDO,x] = hist(surfacetype(~L_DO),[0 1 2 3]);
bar([n_DO; n_notDO]');
set(gca,'XTick',[1 2 3 4], 'XTickLabel',{'plane','1 sheet','2 sheets','ellipsoid'});
% A surprisingly large fraction of 2-sheet hyperboloids are double-opponent
% cells?  Very few double-opponent cells are planes.

% Is there a relationship between threshhold and surface type?
anovan(data(:,7),surfacetype);
boxplot(data(:,7),surfacetype);
set(gca,'XTick',[1 2 3 4], 'XTickLabel',{'plane','1 sheet','2 sheets','ellipsoid'});
[h,p] = ttest2(data(surfacetype==0,7),data(surfacetype==3,7))  % planes vs ellipsoids
ylabel('Threshold');
for i = 1:4
   n(i,:) = hist(data(surfacetype==i-1,7),linspace(0,max(data(:,7)),20)); 
end
figure;
bar(n','stacked')
ylabel('count');
xlabel('threshold');
legend({'plane','1 sheet','2 sheets','ellipsoid'});


%%
% Section 11
% Why are some Gratings color data so poorly fit by the NeuroThresh
% surfaces?
% 1) Different dynamic range in GT and NT?
% 2) Saturating responses in GT (too little variance in firing rates)?
% 3) Surround suppression (size was different in GT and NT)
% 4) The type of surface?
% The two factors than seem to make a difference are making big
% extrapolations and the surface type.  Are these related?

CORRTYPE = 'Pearson';  % 'Pearson' or 'Spearman'
data = [];
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
    % Linear transformations don't change the signs on the eigenvalues,
    % right?
    evals = eig(A);
    
    colors = GTstruct.color.colors;
    rotcols = colors*xformmat;
    planepredresp = NT.sum.exptParams.threshold*abs(rotcols*planeparams);
    trueresp = GTstruct.color.colresp(:,1);

    % Finding predictions from quadratic model
    [th,ph,r] = cart2sph(rotcols(:,1),rotcols(:,2),rotcols(:,3));    
    predr2 = 1./(quadparams(1).*(cos(ph).*cos(th)).^2 +...
                quadparams(2).*(cos(ph).*sin(th)).^2 +...
                quadparams(3).*sin(ph).^2 + ...
                2*quadparams(4).*cos(ph).*cos(th).*cos(ph).*sin(th) +...
                2*quadparams(5).*cos(ph).*cos(th).*sin(ph) +...
                2*quadparams(6).*cos(ph).*sin(th).*sin(ph));
    if (any(predr2<0))
    	predr2(predr2<0) = inf;
    end
    quadpredr = sqrt(predr2);
    quadpredresp = NT.sum.exptParams.threshold.*(r./quadpredr);
    % Computing correlations and other types of error metrics
    correlations = corr([quadpredresp, trueresp],'Type',CORRTYPE);
    Lhitsurface = ~isinf(predr2);
    data = [data; GTstruct.areasummation.prefsize correlations(1,2) mean(r(Lhitsurface)./quadpredr(Lhitsurface)) var(trueresp) sum(evals<0)];
end
Lnan = any(isnan(data),2);
if any(Lnan)
    disp(['Eliminating cell for having a nan: ',num2str(find(Lnan))]);
    data(Lnan,:) = [];
end

figure;
plot(data(:,1),data(:,2),'k.'); lsline;
xlabel('Preferred size'); ylabel('Quality of NT fit');
[r,p] = corrcoef(data(:,[1 2]));
title(['r = ',num2str(r(1,2)),' p = ',num2str(p(1,2))]);

figure;
plot(data(:,3),data(:,2),'k.'); lsline
%set(gca,'Xscale','log');
xlabel('Interpolating --- extrapolating'); ylabel('Quality of NT fit');
[r,p] = corrcoef(data(:,[3 2]));
title(['r = ',num2str(r(1,2)),' p = ',num2str(p(1,2))]);

figure;
plot(sqrt(data(:,4)),data(:,2),'k.');
%set(gca,'Xscale','log');
xlabel('SD of firing rates'); ylabel('Quality of NT fit');
[r,p] = corrcoef(data(:,[4 2]));
title(['r = ',num2str(r(1,2)),' p = ',num2str(p(1,2))]);

figure;
boxplot(data(:,2),data(:,5));
set(gca,'XTick',[1 2 3],'XTickLabel',{'Ellip','1 sheet','2 sheet'});
anovan(data(:,2),data(:,5));
ylabel('Quality of Grating predictions');

%%
% Section 12
% Looking at isoresponse *curves* in the S and L-M directions for luminance
% tuned cells.  Outdated and lame - analysis should be done on fitted
% quadratic surface.  May have bugs. OUTDATED.

[fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\NTlum.txt');
USE10DEG = 1;
load T_cones_smj10;
Sdata = [];
LMdata = [];
p_plane = [];
for cellcounter = 1:size(fnames,1)
    filename = findfile(char(fnames{cellcounter}));
    NT = nex2stro(filename);
    out = NTpreprocess(NT,.4,2);
    scaled = out(:,[2:4]).*repmat(out(:,5), 1,3);
   if (USE10DEG)
        fundamentals = NT.sum.exptParams.fundamentals;
        fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
        mon_spd = NT.sum.exptParams.mon_spd;
        mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
        mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
        M = fundamentals'*mon_spd;
        M10 = T_cones_smj10*mon_spd;
        scaled = ConvertConeContrastBasis(M, M10, NT.sum.exptParams.bkgndrgb, scaled);
    end
    Loog = logical(out(:,7));
    [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);
    
    % Getting the curvature in the S and L-M directions
    npts = 500;
    for i = 1:2
        if (i == 1) % tmp is set of vectors with constant L+M and variable S
            tmp = [ones(npts,2), linspace(-50, 50, npts)'];
        else % tmp is set of vectors with constant L+M and variable L-M
            tmp = [ones(npts,1), ones(npts,1), zeros(npts,1)];
            tmp = tmp+[linspace(-10, 10, npts)', -linspace(-10, 10, npts)', zeros(npts,1)];
        end
        tmp = mkbasis(tmp')';  % Unit vectors in these directions
        xformed = tmp*xformmat;  % vectors (formerly unit) in transformed space
        [th,ph,r] = cart2sph(xformed(:,1),xformed(:,2),xformed(:,3)); % directions in transformed space
        predr2 = 1./(quadparams(1).*(cos(ph).*cos(th)).^2 +...
            quadparams(2).*(cos(ph).*sin(th)).^2 +...
            quadparams(3).*sin(ph).^2 + ...
            2*quadparams(4).*cos(ph).*cos(th).*cos(ph).*sin(th) +...
            2*quadparams(5).*cos(ph).*cos(th).*sin(ph) +...
            2*quadparams(6).*cos(ph).*sin(th).*sin(ph));
        if (any(predr2<0))
            predr2(predr2<0) = inf;
        end
        quadpredr = sqrt(predr2);
        quadpredresp = NT.sum.exptParams.threshold.*(r./quadpredr);
        
        [x,y,z] = sph2cart(th,ph,quadpredr);
        xformed = [x,y,z]*inv(xformmat);
        plot3(xformed(:,1),xformed(:,2),xformed(:,3),'r-');
        plot3(-xformed(:,1),-xformed(:,2),-xformed(:,3),'r-');
        if (i == 1)
            Sdata = cat(3,Sdata,xformed);
        else
            LMdata = cat(3,LMdata,xformed);
        end            
        drawnow;
    end
    
    % Bootstrap test for curvature
    nbootiter = 100;
    
    xyz = scaled*xformmat;    
    [th,ph,r] = cart2sph(xyz(~Loog,1),xyz(~Loog,2),xyz(~Loog,3));
    planer = abs(1./(planeparams(1).*cos(ph).*cos(th)+planeparams(2).*cos(ph).*sin(th)+planeparams(3).*sin(ph)));
    resid = log(r)-log(planer);
    
    nulldist = nan*ones(nbootiter,1);
    wait_h = waitbar(0,'Bootstrapping...');
    disp(['Bootstrapping.  PlaneSSE = ',num2str(planeSSE),' QuadSSE = ',num2str(quadSSE)]);
    options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-6,'TolX',10^-6);
    SSEs = [];
    for j = 1:nbootiter
        waitbar(j/nbootiter, wait_h);
        tmpresid = exp(resid(unidrnd(length(resid),[1 length(resid)])));
        
        tmpx = xyz(:,1);
        tmpy = xyz(:,2);
        tmpz = xyz(:,3);
        
        tmpx(~Loog) = planer.*tmpresid.*cos(ph).*cos(th);
        tmpy(~Loog) = planer.*tmpresid.*cos(ph).*sin(th);
        tmpz(~Loog) = planer.*tmpresid.*sin(ph);
        
        [~, SSE1, exitflag1] = fminsearch(@(x) surfacefiterr2([tmpx tmpy tmpz],x, Loog),planeparams, options);
        [tmp, SSE2, exitflag2] = fminsearch(@(x) surfacefiterr2([tmpx tmpy tmpz],x, Loog),[planeparams;0;0;0], options);
        if (~exitflag1)
            disp('Bonked on a plane fit');
        end
        if (~exitflag2)
            disp('Bonked on a quadratic fit');
        end
        if (~exitflag1 | ~exitflag2)
            j = j-1;
        end
        nulldist(j)= SSE1./SSE2;
        SSEs(j,:) = [SSE1 SSE2];
    end
    close(wait_h);
    p_plane(cellcounter) = sum(nulldist>planeSSE/quadSSE)./nbootiter
    
end
% "data" will have lots of nans and infs that correspond to directions that
% didn't hit the surface

% Looking at curvature in S direction
figure; axes; hold on;
Scurvature = [];
for i = 1:size(Sdata,3)
    L = abs(Sdata(:,3,i)) < .3;  % Don't want to extrapolate too far
    if (sum(L) < 5)
        continue;
    end
    lm = sqrt(Sdata(L,1,i).^2+Sdata(L,2,i).^2);
    plot(Sdata(L,3,i),lm);
    secondderiv = diff(diff(lm));
    Scurvature = [Scurvature; all(secondderiv > 0)];
end
ylabel('L+M contrast'); xlabel('S-cone contrast'); 

[sum(Scurvature > 0) sum(Scurvature == 0)]
p = signtest(Scurvature-.5)

% Looking at curvature in L-M direction
figure; axes; hold on;
LMcurvature = [];
for i = 1:size(LMdata,3)
    rg = LMdata(:,:,i)*[1/sqrt(2) -1/sqrt(2) 0]';
    lum = LMdata(:,:,i)*[1/sqrt(2) +1/sqrt(2) 0]';
    
    L = abs(rg) < .05;  % Don't want to extrapolate too far
    if (sum(L) < 5)
        continue;
    end
    plot(rg(L),lum(L));
    secondderiv = diff(diff(lum(L)));
    LMcurvature = [LMcurvature; all(secondderiv > 0)];
end
ylabel('L+M contrast'); xlabel('L-M contrast');
[sum(LMcurvature > 0) sum(LMcurvature == 0)]
p = signtest(LMcurvature-.5)

Lsig = p_plane < 0.05
Scurvature(Lsig)
LMcurvature(Lsig)

%%
% Section 13
% How far are the best-fit isoresponse planes from the isodetection
% contours (as measured with DTNT)?

DTNTsfs = [0.5008 0.9919 1.9839 3.2238 3.9677];
DTNTparams =[ 1.1817    2.8475    2.8556         0    3.0997;
        0.1203    0.0147   -0.0035         0    0.0012;
       -0.0052   -0.0033    0.0072         0    0.0126;
        0.1445    0.2288    0.1570         0    0.0491;
       -0.0367   -0.5428    0.2218         0   -0.2422;
       -0.3620   -0.2225    0.4693         0    0.2491;
       0.1134    0.0706    0.0446         0    0.0476;
       0.7796    1.2509   -0.7186         0   -0.3305;
       -0.6624   -1.2128    0.5492         0   -0.2676;
      -0.0508   -0.0515    0.0766         0   -0.0105];

[fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\NTlum.txt');
load T_cones_smj10;
data = nan(size(fnames,1),2);
for cellcounter = 1:size(fnames,1)
    filename = findfile(char(fnames{cellcounter}));
    NT = nex2stro(filename);
    sf = NT.trial(:,strcmp(NT.sum.trialFields(1,:),'sf')); % sf wasn't dropped correctly in early files!
    if (sf(1) == 0)
        disp('No sf');
        continue
    end
    out = NTpreprocess(NT,.4,2);
    scaled = out(:,[2:4]).*repmat(out(:,5), 1,3);
    Loog = logical(out(:,7));

    % Conversion to 10 degree fundamentals.
    fundamentals = NT.sum.exptParams.fundamentals;
    fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
    mon_spd = NT.sum.exptParams.mon_spd;
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
    M = fundamentals'*mon_spd;
    M10 = T_cones_smj10*mon_spd;
    scaled = ConvertConeContrastBasis(M, M10, NT.sum.exptParams.bkgndrgb, scaled);
    
    err = abs(DTNTsfs-sf(1));
    whichsf = err == min(err);
    params = DTNTparams(:,whichsf);
    
    [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);
    coneweights = (planeparams\xformmat)';
    planevect =  coneweights'/norm(coneweights).^2*100;  % x100 to adhere to DTNT contrast convention
    [th,phi,r1] = cart2sph(planevect(1),planevect(2),planevect(3));
    options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-10);
    r2 = abs(fminsearch(@(r) colefitpredr(th,phi,r,params),r1,options)); % This is outdated - see Zack's code
    [x,y,z] = sph2cart(th,phi,r2);
    % "normfact" is the scale factor that you have to multiply planeparams by
    % to get to the isodection surface
    normfact = r2./r1
    data(cellcounter,:) = [normfact sf(1)];
end

%%
% Section 14
% Looking at the distribution of isoresponse surface shapes and the
% orientations of these shapes.  

ROTATEEIGS = 0; % Rotate the largest two eigenvectors.  Useful for determining whether the clusters
% in the medium and short axes of the ellipsoids (or other surfaces) could
% happen by chance.  The answer appears to be "no".
evalratiothresh = 0;

data = [];
[fnames, spikeIdx] = fnamesFromTxt2();
[planefnames] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\plane0.01.txt');
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
   % out = NTpreprocess(NT,.4,2);
    out = NTpreprocess(NT,0,Inf);
        
    scaled = out(:,[2:4]).*repmat(out(:,5), 1,3);
    Loog = logical(out(:,7));
    [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);
    r2 = sum(scaled.^2,2);

    % What type of quadric?
    A = [quadparams(1) quadparams(4) quadparams(5);...
        quadparams(4) quadparams(2) quadparams(6);...
        quadparams(5) quadparams(6) quadparams(3)];
    % Transforming quadparams into un-whitened space
    [evecs,evals] = eig(xformmat*A*xformmat');
    [evals,i] = sort(diag(evals),1,'ascend');
    evecs = evecs(:,i);
    if (sum(evals<0) == 1) % 1 sheet
        WHICHSURF = 1;
    elseif (sum(evals<0) == 2) % 2 sheets
        WHICHSURF = 2;
    else %(all(evals>0)) % ellipsoid
        WHICHSURF = 3;
    end
    shortfn = NT.sum.fileName(find(NT.sum.fileName == '\',1,'last')+1:find(NT.sum.fileName == '.',1,'last')-1);

    data(cellcounter).plane = ismember(shortfn,[planefnames{:}]');
    data(cellcounter).r2 = r2;
    data(cellcounter).planeSSE = planeSSE;
    data(cellcounter).quadSSE = quadSSE;
    data(cellcounter).WHICHSURF = WHICHSURF;
    data(cellcounter).evecs = evecs;
    data(cellcounter).evals = evals;
    data(cellcounter).coneweights = xformmat*planeparams;
end
% For plane, nonzero eigenvector points to plane
% For hyperboloid of 2 sheets, positive eigenvector points to surface
% First, let's just see how many of what kinds of cells we've got
out = [];
for i = 1:length(data)
    out(i,:) = [data(i).WHICHSURF (data(i).planeSSE-data(i).quadSSE)/sum(data(i).r2)]; 
end
figure;
subplot(3,1,1);
hist(out(:,1),[1 2 3]); set(gca,'XTick',[1 2 3],'XTicklabel',{'hyp 1','hyp 2','ellip'});
ylabel('Count');
subplot(3,1,2);
plot(out(:,1),out(:,2),'k.');
ylabel('non-planarity index');

% Planar cells tends to get lumped into the two hyperboloid categories

% Verifying that the long axis of the ellipsoids are oriented in the
% S-cone direction.  Note: long axis is small eigenvalue
% Option to set a threshold on the eigenvalue ratios
evalratiothresh = sqrt(5);
for whichsurf = 1:3; % 1 = 1 sheet, 2 = 2 sheets, 3 = ellipsoid
    out = [];
    if (whichsurf == 1)
        whichvecstorotate = [1 2];
    elseif(whichsurf == 2)
        whichvecstorotate = [1 3];
    elseif(whichsurf == 3)
        whichvecstorotate = [2 3];
    end
    
    for i = 1:length(data)
        if (data(i).WHICHSURF == whichsurf)
            if (ROTATEEIGS)
                disp('Rotating eigenvectors');
                th = unifrnd(0,2*pi);
                rotmat = [cos(th) sin(th); -sin(th) cos(th)];
                tmp = data(i).evecs;
                tmp(:,whichvecstorotate) = data(i).evecs(:,whichvecstorotate)*rotmat;
                out = [out; reshape(tmp,1,9)];                
            else
                out = [out; reshape(data(i).evecs,1,9)];
            end
            ratio12 = max(abs(sqrt(data(i).evals(2)./data(i).evals(1))),sqrt(abs(data(i).evals(1)./data(i).evals(2))));
            ratio23 = max(abs(sqrt(data(i).evals(3)./data(i).evals(2))),sqrt(abs(data(i).evals(2)./data(i).evals(3))));
            if (ratio12 < evalratiothresh)
                out(end,[1:6]) = nan;
            end
            if (ratio23 < evalratiothresh)
                out(end,[4:9]) = nan;
            end
        end
    end
    symbols = {'gs','rd','bv'};
    titles = {'Hyperboloid 1 sheet','Hyperboloid 2 sheets','Ellipsoid'};
    figure; axes; hold on;
    for i = 1:3
        idxs = [1 2 3]+(i-1)*3;
        tmp = out(:,idxs)./repmat(sum(abs(out(:,idxs)),2),1,3);
        tmp = tmp.*repmat(sign(tmp(:,2)),1,3);
        plot(tmp(:,1),tmp(:,2),symbols{i},'MarkerFaceColor',symbols{i}(1));
    end
    plot([-1 0 1 -1],[0 1 0 0],'k-');
    title(titles{whichsurf});
    if whichsurf == 1
        legend('Concave (open)','Medium convex','Tight convex');
    elseif whichsurf == 2
        legend('Tight concave','Medium concave','Convex');
    else % whichsurf == 3
        legend('Weak convex','Medium convex','Tight convex');
    end
end

% Comparing eigenvectors to "coneweights" from linear model
whichsurf = 2; % 1 = 1 sheet, 2 = 2 sheets, 3 = ellipsoid
planeidxthreshold = 0;  % Large values are less planar
out = [];
for i = 1:length(data)
    if (data(i).WHICHSURF == whichsurf)
        planeidx = (data(i).planeSSE-data(i).quadSSE)/sum(data(i).r2);
        if (planeidx > planeidxthreshold)
            out = [out; (data(i).coneweights'./norm(abs(data(i).coneweights)))*data(i).evecs];
        end
    end
end
figure; axes; set(gcf,'DefaultAxesColorOrder',[0 1 0; 1 0 0; 0 0 1]);
plot(abs(out)); title(titles{whichsurf});

%For ellipsoidal cells, planes tend to be fitted to the short axis
%For 1 sheet cells, planes tend to be fitted to the axis orthogonal to the
%two axes of curvature.
%For all cells, coneweights tend to be similar to largerest eigenvector,
%which corresponds to the shortest radius convex dimension.

%Looking at the curvature of surface for "luminance cells"
out = [];
whichsurf = 3;
for i = 1:length(data)
    if (data(i).WHICHSURF == whichsurf)
        tmp = data(i).evecs(:,3)'./sum(abs(data(i).evecs(:,3)));
        L = abs(tmp(1)+tmp(2)) > .9; % ad hoc defn of luminance cell
      %  L = abs(tmp(1)-tmp(2)) > .9; % ad hoc defn of L-M cell
       % L = abs(tmp(3)) > .1; % ad hoc defn of an S-dominated cell
        if (L)  
            out = [out; reshape(data(i).evecs,1,9)];
        end
    end
end
figure; axes; hold on;
for i = 1:3
    idxs = [1 2 3]+(i-1)*3;
    tmp = out(:,idxs)./repmat(sum(abs(out(:,idxs)),2),1,3);
    tmp = tmp.*repmat(sign(tmp(:,2)),1,3);
    plot(tmp(:,1),tmp(:,2),symbols{i},'MarkerFaceColor',symbols{i}(1));
end
plot([-1 0 1 -1],[0 1 0 0],'k-');
legend('Tight concave','Medium concave','Convex');
title(titles{whichsurf});

%%
% Section 15
% Looking at how isoresponse surfaces change as a function of the
% threshold.
PLOTFIGS = 0;
[fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\NTmultiples.txt');

% Setting up axes for plotting
trianglefigh = figure; axes; hold on;
plot([-1 1 0 -1],[0 0 1 0],'k-');

maxnfiles = 0;
for i = 1:size(fnames,1)
    maxnfiles = max(maxnfiles, length(fnames{i}));
end
whichsurface = nan*ones(length(fnames), maxnfiles);
thresholds = nan*ones(length(fnames), maxnfiles);


for cellcounter = 1:size(fnames,1)
    NT = {}; 
    for i = 1:size(fnames{cellcounter},2)
        filename = findfile(char(fnames{cellcounter}(i)));
        NT{i} = nex2stro(filename);
    end
    % End of the cell screening part
    allconeweights = []; allquadparams = []; evals = []; scaled = []; Loog = [];
    for i = 1:length(NT)
        out = NTpreprocess(NT{i},0,Inf);  % .4, 2  or 0, Inf
        scaled = [scaled; out(:,[2:4]).*repmat(out(:,5), 1,3) repmat(i,size(out,1),1)];
        Loog = [Loog; logical(out(:,7))];
        [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(out(:,[2:4]).*repmat(out(:,5), 1,3), out(:,7));
        A = [quadparams(1) quadparams(4) quadparams(5);...
            quadparams(4) quadparams(2) quadparams(6);...
            quadparams(5) quadparams(6) quadparams(3)];
        B = xformmat*A*xformmat';
        allconeweights(i,:) = planeparams'*xformmat';
        allquadparams(i,:) = [B(1,1) B(2,2) B(3,3) B(1,2) B(1,3) B(2,3)];
        evals(i,:) = eig(B)';
        whichsurface(cellcounter,i) = sum(eig(B)<0);
        thresholds(cellcounter,i) = NT{i}.sum.exptParams.threshold;
    end
    
    % Now plotting surfaces
    if (PLOTFIGS)
        figure; axes; hold on;
        plotlim = max(abs(scaled(:,[1 2 3])))*1.1;
        [xx yy zz] = meshgrid(linspace(-plotlim(1),plotlim(1),40),...
            linspace(-plotlim(2),plotlim(2),40),...
            linspace(-plotlim(3),plotlim(3),40));
        xformedxyz = [xx(:) yy(:) zz(:)];
        variables = [xformedxyz(:,1).^2 xformedxyz(:,2).^2 xformedxyz(:,3).^2 2*xformedxyz(:,1).*xformedxyz(:,2) 2*xformedxyz(:,1).*xformedxyz(:,3) 2*xformedxyz(:,2).*xformedxyz(:,3)];
        for i = 1:length(NT)
            fr = variables*allquadparams(i,:)';
            surfstruct = isosurface(xx,yy,zz,reshape(fr,size(xx)), 1);
            p = patch(surfstruct);
            col = unifrnd(0,1,3,1);
            L = logical(scaled(:,4) == i);
            plot3(scaled(~Loog&L,1),scaled(~Loog&L,2),scaled(~Loog&L,3),'ko','MarkerFaceColor',col/2,'MarkerEdgeColor','none');
            plot3(-scaled(~Loog&L,1),-scaled(~Loog&L,2),-scaled(~Loog&L,3),'ko','MarkerFaceColor',col/2,'MarkerEdgeColor','none');
            set(p,'FaceAlpha',.5,'FaceColor',col,'Edgealpha',0);
            axis vis3d;
            camlight;
            lighting phong;
            axis square;
        end
        sum(evals<0,2)
        set(gcf,'Name',NT{1}.sum.fileName);
    end
    allconeweights = allconeweights./repmat(sum(abs(allconeweights),2),1,3);
    allconeweights = allconeweights.*repmat(sign(allconeweights(:,2)),1,3);
    figure(trianglefigh);
    h = plot(allconeweights(:,1),allconeweights(:,2),'ko-');
    set(h,'ButtonDownFcn',['disp (''', NT{1}.sum.fileName,''')'])
    set(gca,'Xlim',[-1 1],'YLim',[0 1]);
end

% How do surfaces change as the threshold changes?
x = zeros(3,3);
n = zeros(3,1);
for i = 1:size(thresholds,1)
    [a,b] = sort(thresholds(i,:),2,'ascend')
    for j = 1:sum(~isnan(thresholds(i,:)))
        shapes = whichsurface(i,b(~isnan(a)))
        n(shapes(j)+1) = n(shapes(j)+1)+1;
        if (j<sum(~isnan(thresholds(i,:))))
            x(shapes(j)+1,shapes(j+1)+1) = x(shapes(j)+1,shapes(j+1)+1)+1
        end
    end
end

transprob = x./repmat(sum(x,2),1,3);  % transition probabilities
figure; axes; hold on;
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
                scalefactor = .2
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


%%
% Section 15.1.  Picking which of the multiple thresholds to use for
% comparison with the gratings data.  Either using the NT file that's
% closest in time to teh gratings data or the NT file with the greatest
% number of data points.

[fnames1, spikeIdx1] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\NTmultiples.txt');
[fnames2, spikeIdx2] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\Gratings&NT.txt');
fnamestochoose = []; 
for i = 1:size(fnames1,1)
    i
    for j = 1:size(fnames2,1)
        if any(ismember(fnames1{i}, fnames2{j}))
            % First find gratings data
            for k = 1:length(fnames2{j})
                filename = findfile(char(fnames2{j}(k)));
                paradigmID = getparadigmID(filename);
                if paradigmID ~= 150
                    continue
                else
                    k = size(fnames2{j}); % break out of the loop
                    GT = nex2stro(filename);
                    GTstruct = getGratingTuning(GT,1);
                    GTcols = GTstruct.color.colors;
                    trueresp = GTstruct.color.colresp(:,1);
                end
            end
            % getting the NT data
            data = [];
            for k = 1:length(fnames1{i})
                filename = findfile(char(fnames1{i}(k)));
                NT = nex2stro(filename);
%                data = [data; NT.sum.exptParams.threshold max(trueresp)];
                
                 out = NTpreprocess(NT,0,Inf);  % Getting the termination points
                 scaled = out(:,[2 3 4]).*repmat(out(:,5),1,3);
                 Loog = logical(out(:,end));
                 data = [data; sum(~Loog)];  % basing decision on # of data points
%                 [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);
%                 rotcols = GTcols*xformmat;
%                
%                 % Finding predictions from quadratic model
%                 [th,ph,r] = cart2sph(rotcols(:,1),rotcols(:,2),rotcols(:,3));
%                 predr2 = 1./(quadparams(1).*(cos(ph).*cos(th)).^2 +...
%                     quadparams(2).*(cos(ph).*sin(th)).^2 +...
%                     quadparams(3).*sin(ph).^2 + ...
%                     2*quadparams(4).*cos(ph).*cos(th).*cos(ph).*sin(th) +...
%                     2*quadparams(5).*cos(ph).*cos(th).*sin(ph) +...
%                     2*quadparams(6).*cos(ph).*sin(th).*sin(ph));
%                 
%                 predr2(predr2<0) = Inf;
%                 quadpredresp = NT.sum.exptParams.threshold.*(r./sqrt(predr2));
%                 data = [data; mean(quadpredresp-trueresp)];
%                 
                % Trying some plotting
                figure(k); clf; axes; hold on; set(k,'Name','blank');
                plot3(scaled(~Loog,1),scaled(~Loog,2),scaled(~Loog,3),'k.')
                plot3(-scaled(~Loog,1),-scaled(~Loog,2),-scaled(~Loog,3),'k.')
                plot3(GTcols(:,1),GTcols(:,2),GTcols(:,3),'m*');
                plot3(-GTcols(:,1),-GTcols(:,2),-GTcols(:,3),'m*');
            end
          %  err = abs(data(:,1)-data(:,2));
          %  kidx = find(err == min(err));
            kidx = find(data == max(data),1,'first');
            set(kidx,'Name','This one');
            fnamestochoose = [fnamestochoose; char(fnames1{i}(kidx))];
           % pause
        end
    end
end

%%
% Section 15.2
% Fitting quadratic surfaces to the early and late part of the data files
% separately.  How often do they change?

[fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\Gratings&NT.txt');
data = [];
for cellcounter = 1:size(fnames,1)
    NT = {};
    for i = 1:size(fnames{cellcounter},2)
        filename = findfile(char(fnames{cellcounter}(i)));
        paradigmID = getparadigmID(filename);
        if paradigmID == 103
            NT = nex2stro(filename);
        end
    end
    out = NTpreprocess(NT,0,Inf);  % .4, 2  or 0, Inf
    for i = 1:2
        if (i == 1)
            L = 1:floor(size(out,1)/2);
        else
            L = floor(size(out,1)/2):size(out,1);
        end
        scaled = out(L,[2:4]).*repmat(out(L,5), 1,3);
        Loog = logical(out(L,7));
        if (sum(~Loog) < 3)
            disp('Fewer than 3 points in gamut');
            continue;
        end
        [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);
        % What type of quadric?
        A = [quadparams(1) quadparams(4) quadparams(5);...
            quadparams(4) quadparams(2) quadparams(6);...
            quadparams(5) quadparams(6) quadparams(3)];
        d = eig(A);
        data(cellcounter,i) = sum(d<0);
    end
end

x = zeros(3,3);
n = zeros(3,1);
for i = 0:2
    for j = i:2
        n(i+1) = sum(data(:,1) == i);
        x(i+1,j+1) = sum(data(:,1) == i & data(:,2) == j);
    end
end

x = [8 4 0; 1 5 4; 0 4 2];
y = [4 4 0; 0 10 2; 0 0 3];
z = [1 0 0; 1 3 0; 0 2 1];

expected = sum(x)'*sum(y,2)'./sum([y(:);x(:)]);




%%
% Section 16: Permutation tests on the difference between correlation
% coefficients (quadpred and response vs. plane pred and response).

CORRTYPE = 'Pearson';  % 'Pearson' or 'Spearman'
p = [];
teststat = [];
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
    out = NTpreprocess(NT,0,Inf);  % .4, 2  or 0, Inf
        
    scaled = out(:,[2:4]).*repmat(out(:,5), 1,3);
    Loog = logical(out(:,7));
    if (sum(~Loog) < 3)
        disp('Fewer than 3 points in gamut');
        continue;
    end
    [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);
    xyz = scaled*xformmat;
    
    colors = GTstruct.color.colors;
    rotcols = colors*xformmat;
    planepredresp = NT.sum.exptParams.threshold*abs(rotcols*planeparams);
    trueresp = GTstruct.color.colresp(:,1);

    % Finding predictions from quadratic model
    [th,ph,r] = cart2sph(rotcols(:,1),rotcols(:,2),rotcols(:,3));
    predr2 = 1./(quadparams(1).*(cos(ph).*cos(th)).^2 +...
    quadparams(2).*(cos(ph).*sin(th)).^2 +...
    quadparams(3).*sin(ph).^2 + ...
    2*quadparams(4).*cos(ph).*cos(th).*cos(ph).*sin(th) +...
    2*quadparams(5).*cos(ph).*cos(th).*sin(ph) +...
    2*quadparams(6).*cos(ph).*sin(th).*sin(ph));
    
    predr2(predr2<0) = Inf;
    quadpredresp = NT.sum.exptParams.threshold.*(r./sqrt(predr2));
    
    % Computing correlations
    r = corr([trueresp, planepredresp, quadpredresp],'Type',CORRTYPE);
    teststat(cellcounter) = r(1,3)-r(1,2);
    data = [];
    bothpredresp = [quadpredresp, trueresp];
    ncols = size(trueresp,1);
    for i = 1:2^ncols-1
        screen = zeros(ncols,1)+dec2bin(i,ncols)'-48;
        screen = [screen, ~screen];
        r = corr([trueresp, sum(bothpredresp.*screen,2), sum(bothpredresp.*~screen,2)],'Type',CORRTYPE);
        data(i) = r(1,3)-r(1,2);
    end
    p(cellcounter) = sum(data > teststat(cellcounter))./(2^ncols-1)
end

figure; axes; hold on;
plot(p,teststat,'k.'); xlabel('p'); ylabel('diff in corr coef');
tmp = [];
for i = find(p < 0.05)
    tmp = [tmp; char(fnames{i}(2))];
end
% List of neurons for which quad preds are better than plane preds
tmp

%%
% Section 17
% Trying to fit a constrained model to the data wherein isoresponse
% surfaces are scaled versions of each other.
[fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\NTmultiples.txt');

maxnfiles = 0;
for i = 1:size(fnames,1)
    maxnfiles = max(maxnfiles, length(fnames{i}));
end
thresholds = nan*ones(length(fnames), maxnfiles);
data = [];
for cellcounter = 1:size(fnames,1)
    NT = {}; 
    for i = 1:size(fnames{cellcounter},2)
        filename = findfile(char(fnames{cellcounter}(i)));
        NT{i} = nex2stro(filename);
    end
    
    allplaneparams = []; allquadparams = []; scaled = []; Loog = [];
    for i = 1:length(NT)
        out = NTpreprocess(NT{i},0,Inf);  % .4, 2  or 0, Inf
        scaled = [scaled; out(:,[2:4]).*repmat(out(:,5), 1,3) repmat(i,size(out,1),1)];
        Loog = [Loog; logical(out(:,7))];
        [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(out(:,[2:4]).*repmat(out(:,5), 1,3), out(:,7));
        A = [quadparams(1) quadparams(4) quadparams(5);...
            quadparams(4) quadparams(2) quadparams(6);...
            quadparams(5) quadparams(6) quadparams(3)];
        B = xformmat*A*xformmat';
        allquadparams(i,:) = [B(1,1) B(2,2) B(3,3) B(1,2) B(1,3) B(2,3)];
        allplaneparams(i,:) = (xformmat*planeparams)';
        thresholds(cellcounter,i) = NT{i}.sum.exptParams.threshold;
    end
    
    % Plotting
    colors = [1 0 0; 0 1 0; 0 0 1];
    plotlim = max(abs(scaled(~Loog,[1 2 3])));
    [xx yy zz] = meshgrid(linspace(-plotlim(1),plotlim(1),20),...
        linspace(-plotlim(2),plotlim(2),20),...
        linspace(-plotlim(3),plotlim(3),20));
    figure(cellcounter); axes; hold on;
    for i = 1:length(NT)
        L = scaled(:,end) ==  i;
        h(1) = plot3(scaled(~Loog&L,1),scaled(~Loog&L,2),scaled(~Loog&L,3),'ko');
        h(2) = plot3(-scaled(~Loog&L,1),-scaled(~Loog&L,2),-scaled(~Loog&L,3),'ko');
        plot3([zeros(sum(Loog&L),1) scaled(Loog&L,1)]', [zeros(sum(Loog&L),1) scaled(Loog&L,2)]', [zeros(sum(Loog&L),1) scaled(Loog&L,3)]','-','Color',colors(i,:));
        plot3([zeros(sum(Loog&L),1) -scaled(Loog&L,1)]', [zeros(sum(Loog&L),1) -scaled(Loog&L,2)]', [zeros(sum(Loog&L),1) -scaled(Loog&L,3)]','-','Color',colors(i,:));
        set(h,'MarkerFaceColor',colors(i,:));
        
%         % Now plotting (individually fitted) surfaces
%         xformedxyz = [xx(:) yy(:) zz(:)];
%         variables = [xformedxyz(:,1).^2 xformedxyz(:,2).^2 xformedxyz(:,3).^2 2*xformedxyz(:,1).*xformedxyz(:,2) 2*xformedxyz(:,1).*xformedxyz(:,3) 2*xformedxyz(:,2).*xformedxyz(:,3)];
%         coefficients = allquadparams(i,:);
%         fr = variables*coefficients';
%         surfstruct = isosurface(xx,yy,zz,reshape(fr,size(xx)), 1);
%         
%         p = patch(surfstruct);
%         set(p, 'FaceAlpha',.2,'FaceColor', colors(i,:), 'EdgeColor', 'none');
%         axis vis3d;
%         camlight;
%         lighting phong;
%         axis square;
    end
    
    
    
%     % Trying to fit a constrained model to both data sets
%     % First trying a 2-D grid search with planes
% This has lower quadSSEs than starting with mean of planes,
% but for most neurons the fit is identical. This way seems the best.
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
    quadSSE = Inf;
    quadparams = zeros(1,6);
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
        if (tmpquadSSE < quadSSE)
            quadparams = tmpquadparams;
            quadSSE = tmpquadSSE;
        end
        drawnow;
    end

    % Trying a new way of getting the initial guess - use the individual
    % best-fit planes and combine them by taking the means. #2
%     
%     [v,d] = eig(cov([scaled(~Loog,[1 2 3]); -scaled(~Loog,[1 2 3])]));
%     d = diag(d);
%     if (min(d) < 2*eps)
%         disp('Too few data points for whitening');
%         whtmat = eye(3);
%     else
%         whtmat = v*diag(sqrt(1./d));
%     end
%     newscaled = [scaled(:,[1 2 3])*whtmat scaled(:,4)];
%     whtplaneparams = allplaneparams*inv(whtmat');
%     [th,ph,r] = cart2sph(whtplaneparams(:,1),whtplaneparams(:,2),whtplaneparams(:,3));
%     % Average the angles together and use the same 'r's
%     meanunitvect = mean(whtplaneparams)./sqrt(sum(mean(whtplaneparams).^2));
%     initplaneguess = meanunitvect*r(1);
%     initscaleguess = r(2).^2/r(1).^2;  % not sure if this is right.
%     initquadguess = [initplaneguess(1)^2 initplaneguess(2)^2 initplaneguess(3)^2 initplaneguess(1)*initplaneguess(2) initplaneguess(1)*initplaneguess(3) initplaneguess(2)*initplaneguess(3)]';
%     options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-6,'TolX',10^-6,'Display','off');
%     [quadparams, quadSSE, exitflag] = fminsearch(@(x) surfacefiterr3(newscaled, x, Loog),[initquadguess' initscaleguess],options);
%     
%     % undoing the whitening
%     A = [quadparams(1) quadparams(4) quadparams(5);...
%          quadparams(4) quadparams(2) quadparams(6);...
%          quadparams(5) quadparams(6) quadparams(3)];
%     B = whtmat*A*whtmat';
%    

   % Trying yet another a new way of getting the initial guess - 
   % use the mean of allquadparams.  #3
%     
%     [v,d] = eig(cov([scaled(~Loog,[1 2 3]); -scaled(~Loog,[1 2 3])]));
%     d = diag(d);
%     if (min(d) < 2*eps)
%         disp('Too few data points for whitening');
%         whtmat = eye(3);
%     else
%         whtmat = v*diag(sqrt(1./d));
%     end
%     newscaled = [scaled(:,[1 2 3])*whtmat scaled(:,4)];
% 
%     meanquadparams = mean(allquadparams);
%     A1 = [meanquadparams(1) meanquadparams(4) meanquadparams(5);...
%         meanquadparams(4) meanquadparams(2) meanquadparams(6);...
%         meanquadparams(5) meanquadparams(6) meanquadparams(3)];
%     B = inv(whtmat')*A1*inv(whtmat);
%     
%     scalefactors = [];
%     for i = 1:2
%         A = [allquadparams(i,1) allquadparams(i,4) allquadparams(i,5);...
%             allquadparams(i,4) allquadparams(i,2) allquadparams(i,6);...
%             allquadparams(i,5) allquadparams(i,6) allquadparams(i,3)];
%         d(i) = max(eig(inv(whtmat')*A*inv(whtmat)));
%         scalefactors(i) = max(eig(A))./max(eig(A1));
%     end
%     initquadguess = [B(1,1) B(2,2) B(3,3) B(1,2) B(1,3) B(2,3)]*scalefactors(1);  % This causes max eig = 1. GDLH 6/7/11
%     initscaleguess = scalefactors(2)./scalefactors(1);
%     options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-6,'TolX',10^-6,'Display','off');
%     [quadparams, quadSSE, exitflag] = fminsearch(@(x) surfacefiterr3(newscaled, x, Loog),[initquadguess initscaleguess],options);
%     
%     % undoing the whitening
%     A = [quadparams(1) quadparams(4) quadparams(5);...
%          quadparams(4) quadparams(2) quadparams(6);...
%          quadparams(5) quadparams(6) quadparams(3)];
%     B = whtmat*A*whtmat';
%    


    % Now plotting (jointly fitted) surfaces
    [xx yy zz] = meshgrid(linspace(-plotlim(1),plotlim(1),20),...
        linspace(-plotlim(2),plotlim(2),20),...
        linspace(-plotlim(3),plotlim(3),20));
  
    xformedxyz = [xx(:) yy(:) zz(:)];
    variables = [xformedxyz(:,1).^2 xformedxyz(:,2).^2 xformedxyz(:,3).^2 2*xformedxyz(:,1).*xformedxyz(:,2) 2*xformedxyz(:,1).*xformedxyz(:,3) 2*xformedxyz(:,2).*xformedxyz(:,3)];
    coefficients = [B(1,1) B(2,2) B(3,3) B(1,2) B(1,3) B(2,3)];
    figure(cellcounter);
    for i = 1:2
        if (i ==1)
            fr = variables*coefficients';
        else
            fr = variables*coefficients'.*quadparams(7);           
        end
        surfstruct = isosurface(xx,yy,zz,reshape(fr,size(xx)), 1);
        
        p = patch(surfstruct);
        set(p, 'FaceAlpha',.4,'FaceColor', colors(i,:), 'EdgeColor', 'none');
        axis vis3d;
        camlight;
        lighting phong;
        axis square;
    end
    drawnow;
    data = [data; quadSSE];
end


%%
% Section 17.1
% Trying to fit a constrained model to the data wherein isoresponse
% surfaces are scaled versions of each other.  Pseudo R-squared
[fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\NTmultiples.txt');

maxnfiles = 0;
for i = 1:size(fnames,1)
    maxnfiles = max(maxnfiles, length(fnames{i}));
end
thresholds = nan*ones(length(fnames), maxnfiles);
data = [];
parameters = [];
for cellcounter = 1:size(fnames,1)
    NT = {}; 
    for i = 1:size(fnames{cellcounter},2)
        filename = findfile(char(fnames{cellcounter}(i)));
        NT{i} = nex2stro(filename);
    end
    
    allplaneparams = []; allquadparams = []; scaled = []; Loog = []; allquadSSEs = [];
    for i = 1:length(NT)
        out = NTpreprocess(NT{i},0,Inf);  % .4, 2  or 0, Inf
        scaled = [scaled; out(:,[2:4]).*repmat(out(:,5), 1,3) repmat(i,size(out,1),1)];
        Loog = [Loog; logical(out(:,7))];
    end
    [v,d] = eig(cov([scaled(~Loog,[1 2 3]); -scaled(~Loog,[1 2 3])]));
    d = diag(d);
    if (min(d) < 2*eps)
        disp('Too few data points for whitening');
        whtmat = eye(3);
    else
        whtmat = v*diag(sqrt(1./d));
    end
    newscaled = [scaled(:,[1 2 3])*whtmat scaled(:,4)];
    for i = 1:length(NT)
        L = newscaled(:,4) == i;
        [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled(L,[1 2 3]), Loog(L));
        A = [quadparams(1) quadparams(4) quadparams(5);...
            quadparams(4) quadparams(2) quadparams(6);...
            quadparams(5) quadparams(6) quadparams(3)];
        B = xformmat*A*xformmat';
        allquadparams(i,:) = [B(1,1) B(2,2) B(3,3) B(1,2) B(1,3) B(2,3)];
        allplaneparams(i,:) = (xformmat*planeparams)';
        allquadSSEs(i) = quadSSE;
        thresholds(cellcounter,i) = NT{i}.sum.exptParams.threshold;
    end
 
    
    % Plotting
    colors = [1 0 0; 0 1 0; 0 0 1];
    switch (sum(d<0))
        case 1
        plotlim = max(abs(newscaled))/1.1;
        case 2
        plotlim = max(abs(newscaled))/3;
        otherwise
        plotlim = max(abs(newscaled))*1.1;
    end
%     [xx yy zz] = meshgrid(linspace(-plotlim(1),plotlim(1),20),...
%         linspace(-plotlim(2),plotlim(2),20),...
%         linspace(-plotlim(3),plotlim(3),20));
%     figure(cellcounter); axes; hold on;
%     for i = 1:length(NT)
%         L = scaled(:,end) ==  i;
%         h(1) = plot3(scaled(~Loog&L,1),scaled(~Loog&L,2),scaled(~Loog&L,3),'ko');
%         h(2) = plot3(-scaled(~Loog&L,1),-scaled(~Loog&L,2),-scaled(~Loog&L,3),'ko');
%         plot3([zeros(sum(Loog&L),1) scaled(Loog&L,1)]', [zeros(sum(Loog&L),1) scaled(Loog&L,2)]', [zeros(sum(Loog&L),1) scaled(Loog&L,3)]','-','Color',colors(i,:));
%         plot3([zeros(sum(Loog&L),1) -scaled(Loog&L,1)]', [zeros(sum(Loog&L),1) -scaled(Loog&L,2)]', [zeros(sum(Loog&L),1) -scaled(Loog&L,3)]','-','Color',colors(i,:));
%         set(h,'MarkerFaceColor',colors(i,:));
%     end
%     
    
    % A constrained fit to both sets of data 

  
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
    quadSSE = Inf;
    quadparams = zeros(1,6);
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
    B = [quadparams(1) quadparams(4) quadparams(5);...
        quadparams(4) quadparams(2) quadparams(6);...
        quadparams(5) quadparams(6) quadparams(3)];
    [v,d] = eig(B);
    % For debugging - trying some random start points
    for i = 1:2000
        if (i < 1000)
            tmpd = diag(diag(d)+ normrnd(0,.5*max(abs(diag(d))),3,1));
            tmpB = v*tmpd*v';
            randomguess = [tmpB(1,1) tmpB(2,2) tmpB(3,3) tmpB(1,2) tmpB(1,3) tmpB(2,3) quadparams(7)];
        elseif (i >= 1000)
            newd = diag(diag(d)+ normrnd(0,.5*max(abs(diag(d))),3,1));
            phi = unifrnd(-pi/12, pi/12);
            theta = unifrnd(-pi/12, pi/12);
            psi = unifrnd(-pi/12, pi/12);
            % phi, theta, psi
            rotmat = [cos(theta)*cos(psi) -cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi) sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi);...
                cos(theta)*sin(psi) cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi) -sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi)
                -sin(theta) sin(phi)*cos(theta) cos(phi)*cos(theta)];
            tmpv = rotmat*v;
            tmpB = tmpv*tmpd*tmpv';
            randomguess = [tmpB(1,1) tmpB(2,2) tmpB(3,3) tmpB(1,2) tmpB(1,3) tmpB(2,3) quadparams(7)];
        end
        if (i > 1500)
            randomguess(end) = unifrnd(-500, 500);
        end        
        [tmpquadparams, tmpquadSSE, exitflag] = fminsearch(@(x) surfacefiterr3([newscaled(:,[1 2 3]) scaled(:,end)], x, Loog),randomguess,options);
        if (tmpquadSSE < quadSSE - .01*quadSSE) % more that a 1% difference required
            disp(['Got here on round ',num2str(i)]);
            [tmpquadSSE quadSSE]
            quadparams = tmpquadparams;
            quadSSE = tmpquadSSE;
        end
    end
    B = [quadparams(1) quadparams(4) quadparams(5);...
        quadparams(4) quadparams(2) quadparams(6);...
        quadparams(5) quadparams(6) quadparams(3)];

    % Now plotting (jointly fitted) surfaces
    [xx yy zz] = meshgrid(linspace(-plotlim(1),plotlim(1),50),...
        linspace(-plotlim(2),plotlim(2),50),...
        linspace(-plotlim(3),plotlim(3),50));
    
%     figure(cellcounter);
%     for i = 1:2
%         if (i ==1)
%             fr = variables*coefficients';
%         else
%             fr = variables*coefficients'.*sqrt(quadparams(7));           
%         end
%         surfstruct = isosurface(xx,yy,zz,reshape(fr,size(xx)), 1);
%         
%         p = patch(surfstruct);
%         set(p, 'FaceAlpha',.4,'FaceColor', colors(i,:), 'EdgeColor', 'none');
%         axis vis3d;
%         camlight;
%         lighting phong;
%         axis square;
%     end
%     drawnow;
    [th,ph,r] = cart2sph(newscaled(:,1),newscaled(:,2),newscaled(:,3));
    SST = sum((log(r)-mean(log(r))).^2);
    PR2_quad = 1-quadSSE/SST;
    fr1 = NT{1}.sum.exptParams.threshold;
    fr2 = NT{2}.sum.exptParams.threshold;
    data = [data; sum(allquadSSEs) quadSSE PR2_quad fr1 fr2]
    parameters = [parameters; quadparams];
    
    % Plotting both data sets (one scaled to match the other)
    figure; axes; hold on
    L = scaled(:,4) == 1;
    plot3(newscaled(L&~Loog,1),newscaled(L&~Loog,2),newscaled(L&~Loog,3),'ro','MarkerFaceColor','red');
    plot3(-newscaled(L&~Loog,1),-newscaled(L&~Loog,2),-newscaled(L&~Loog,3),'ro','MarkerFaceColor','red');
    plot3([-newscaled(L&Loog,1) newscaled(L&Loog,1)]',[-newscaled(L&Loog,2) newscaled(L&Loog,2)]',[-newscaled(L&Loog,3) newscaled(L&Loog,3)]','-','Color',[1 .8 .8]);
    % The other group of points
    plot3(newscaled(~L&~Loog,1).*sqrt(quadparams(7)),newscaled(~L&~Loog,2).*sqrt(quadparams(7)),newscaled(~L&~Loog,3).*sqrt(quadparams(7)),'bo','MarkerFaceColor','blue');
    plot3(-newscaled(~L&~Loog,1).*sqrt(quadparams(7)),-newscaled(~L&~Loog,2).*sqrt(quadparams(7)),-newscaled(~L&~Loog,3).*sqrt(quadparams(7)),'bo','MarkerFaceColor','blue');
    h =plot3([-newscaled(~L&Loog,1) newscaled(~L&Loog,1)]',[-newscaled(~L&Loog,2) newscaled(~L&Loog,2)]',[-newscaled(~L&Loog,3) newscaled(~L&Loog,3)]','-','Color',[.8 .8 1]);
    title(NT{1}.sum.fileName);
   % Now plotting (jointly fitted) surfaces
    [xx yy zz] = meshgrid(linspace(-plotlim(1),plotlim(1),50),...
        linspace(-plotlim(2),plotlim(2),50),...
        linspace(-plotlim(3),plotlim(3),50));  
    xformedxyz = [xx(:) yy(:) zz(:)];
    variables = [xformedxyz(:,1).^2 xformedxyz(:,2).^2 xformedxyz(:,3).^2 2*xformedxyz(:,1).*xformedxyz(:,2) 2*xformedxyz(:,1).*xformedxyz(:,3) 2*xformedxyz(:,2).*xformedxyz(:,3)];
    fr = variables*quadparams([1:6])'.*sqrt(quadparams(7));

    surfstruct = isosurface(xx,yy,zz,reshape(fr,size(xx)), 1);
    p = patch(surfstruct);
    set(p, 'FaceAlpha',.3,'FaceColor',[1 .8 1], 'EdgeColor', 'none');
    axis vis3d;
    camlight;
    lighting phong;
    set(p,'SpecularColorReflectance',.4);
    set(p,'SpecularExponent',10,'SpecularStrength',.2);
    set(p,'DiffuseStrength',.7,'AmbientStrength',.4);
    axis square;
    drawnow;
    
end
% Order:  full, restricted, psuedo-R2, target fr 1, target fr 2

% Squared error is log likelihood under Gaussian assumptions.
% x1 = -sum(x_i-mu_i).^2  % Where's the variance?!
% -2*log(exp(x1)/exp(x2)) 
% -2*log(exp(x1))-log(exp(x2))
% -2*(x1-x2)

% Try F-test

D = -2*(data(:,1)-data(:,2));  % This is not right?
p = 1-chi2cdf(D,12-7)
figure; subplot(2,1,1);
hist(p)
subplot(2,1,2);
%plot(abs(thresholds(:,1)-thresholds(:,2)),D,'k.')
% Correlation of deviance and the ratio of thresholds
[r,p1] = corr([max(thresholds,[],2)./min(thresholds,[],2),D],'type','Spearman')
plot(log10(max(thresholds,[],2)./min(thresholds,[],2)),D,'ko','MarkerFaceColor','black','MarkerSize',3)
set(gca,'XTick',log10([1 2 3 4]),'XTickLabel',[1 2 3 4]);
xlabel('Ratio of target firing rates');
ylabel('Deviance');


% As the thresholds get farther apart the D statistic (the deviance)
% increases

% How about if we restrict our attention to the data files where the
% threshold didn't change much at all.
L = max(thresholds,[],2)./min(thresholds,[],2) == 1.01;
% What fraction of the experiments with the small changes in threshold had
% significant changes in surface shape?
sum(p(L) < 0.05)./sum(L) % 1/3, that's a lot
sum(p(~L) < 0.05)./sum(~L) 



% Two sample test on the deviance statistics
[h,p] = ttest2(D(L),D(~L))
[sum(L) sum(~L) mean(D(L)) mean(D(~L)) p]

% One sample test on the deviance stats in expts with small threshold
% changes
[h,p] = ttest(D(L)-5)
% One sample test on the deviance stats in expts with large threshold
% changes
[h,p] = ttest(D(~L)-5)

%%
% Section 18
% How many staircases get culled because of differences between online and
% offline spike sorting?  How many stimuli per trial on average?

[fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\Gratings&NT2.txt');
data = [];
for cellcounter = 1:size(fnames,1)
    NT = {};
    for i = 1:size(fnames{cellcounter},2)
        filename = findfile(char(fnames{cellcounter}(i)));
        paradigmID = getparadigmID(filename);
        if paradigmID == 103
            NT = nex2stro(filename);
        end
    end
    out1 = NTpreprocess(NT,0,Inf);
    out2 = NTpreprocess(NT,0,Inf,Inf);
    data = [data; size(out1,1) size(out2,1) sum(~isnan(NT.trial(:,1))) sum(isnan(NT.trial(:,1)))]
end
sum(data(:,1))./sum(data(:,2)) % fraction of staircases omitted by SORTTHRESH
3*size(data,1)./sum(data(:,1))  % fraction of starcases in the first three directions

a = data(:,3)+data(:,4) % Total number of stimulus presentations
b = data(:,3) % Number of fixations (number of stimulus presentation epochs)
a./b
sum(a)./sum(b)

%%
% Section 19
% Cone contrasts for the gratings in 10 deg space
% Also looking at how step size changes with reversals

[fnames, spikeIdx] = fnamesFromTxt2();
load ('T_cones_smj10');
data = [];
for cellcounter = 1:size(fnames,1)
    GT = {};
    for i = 1:size(fnames{cellcounter},2)
        filename = findfile(char(fnames{cellcounter}(i)));
        paradigmID = getparadigmID(filename);
        if paradigmID == 150
            GT = nex2stro(filename);
            GTstruct = getGratingTuning(GT,1);
        else
            NT = nex2stro(filename);
            bkgndrgb = NT.sum.exptParams.bkgndrgb;  % Ugly hack - shoud drop bkgndrgb in gratings.
        end
    end
    
    fundamentals = GT.sum.exptParams.fundamentals;
    fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
    mon_spd = GT.sum.exptParams.mon_spd;
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
    M = fundamentals'*mon_spd;
    M10 = T_cones_smj10*mon_spd;
%    data(:,:,cellcounter) = GTstruct.color.colors;
    data(:,:,cellcounter) = ConvertConeContrastBasis(M, M10, bkgndrgb, GTstruct.color.colors);
end

figure; axes; hold on;
for i = 1:size(data,1)
    for j = 1:size(data,2)
        plot(squeeze(data(i,j,:)),'k.');
    end
end
% One cell is off the lattice - should pull it

%%
% Section 20: After the first reversal, does the step size decrease?
[fnames, spikeIdx] = fnamesFromTxt2();
for cellcounter = 10:size(fnames,1)
    NT = {};
    for i = 1:size(fnames{cellcounter},2)
        filename = findfile(char(fnames{cellcounter}(i)));
        paradigmID = getparadigmID(filename);
        if paradigmID == 103
            stro = nex2stro(filename);
        end
    end
    % Getting a bunch of important stuff from the stro file.
    spikeidx = strcmp(stro.sum.rasterCells(1,:),getSpikenum(stro,'first'));
    coloridxs = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'coloridx'));
    levels = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'levelidx'));
    reversals = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'reversal'));
    stepsize = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stepsize'));  
    lmsidxs = [find(strcmp(stro.sum.trialFields(1,:),'lcont'))...
        find(strcmp(stro.sum.trialFields(1,:),'mcont'))...
        find(strcmp(stro.sum.trialFields(1,:),'scont'))];
    lms = stro.trial(:,lmsidxs);
    uniquecoloridxs = unique(coloridxs);
    for i = uniquecoloridxs'
        L = logical(coloridxs == i);
        if (any(stepsize(L) == 0))  % Skip color directions with 0 stepsize
            continue;
        end
        rev = reversals(L);
        ss = stepsize(L);
        lev = unique(levels(L));
        [u,s,v] = svd(lms(L,:));
        unitvector = v(:,1);
        contrasts = lms(L,:)*unitvector;
        if(sum(contrasts) < 0)
            unitvector = -unitvector;
            contrasts = -contrasts;
        end
        newc = contrasts(1);
        stepdir = 1;
        for j = 2:length(rev)
            if (rev(j) == 0)
                if (stepdir == 1)
                    newc(j) = newc(j-1)*(1+ss(j-1));
                else
                    newc(j) = newc(j-1)*(1-ss(j-1));
                end
            elseif (rev(j) == 1)
                newc(j) = newc(j-1)*(1+ss(j-1));
                stepdir = -1;
            else % rev(j) == -1
                newc(j) = newc(j-1)*(1-ss(j-1));
                stepdir = 1;
            end
        end
        tmp = [contrasts newc'];
       % if (any(tmp(:,1)-tmp(:,2) > 10.^-4))
            figure; axes; hold on;
            plot(tmp(:,1),'m.-');
            plot(tmp(:,2),'k.-');
            keyboard;
       % end
    end
end

%%
% Section 21
% Using cross-validation to compare the quality of fits between linear and
% nonlinear models. Can I do a Wilcoxon (or something) on the prediction
% errors?

data = [];
[fnames, spikeIdx] = fnamesFromTxt2();
for cellcounter = 1:size(fnames,1)
    NT = {};
    for i = 1:size(fnames{cellcounter},2)
        filename = findfile(char(fnames{cellcounter}(i)));
        paradigmID = getparadigmID(filename);
        if paradigmID == 103
            NT = nex2stro(filename);
        end
    end
    out = NTpreprocess(NT,0,Inf);
    
    scaled = out(:,[2:4]).*repmat(out(:,5), 1,3);
    Loog = logical(out(:,7));
    tmpdata = [];
    for i = find(~Loog)'
        tmpscaled = scaled; tmpscaled(i,:) = [];
        tmpLoog = Loog; tmpLoog(i,:) = [];
        testpoint = scaled(i,:);
        [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(tmpscaled, tmpLoog);
        xformedtestpoint = testpoint*xformmat;
        [th,ph,true_r] = cart2sph(xformedtestpoint(1),xformedtestpoint(2),xformedtestpoint(3));
        planepred_r = true_r./abs(xformedtestpoint*planeparams);
  
       % planepred_r = abs(1./(planeparams(1).*(cos(ph).*cos(th)) +...
       %     planeparams(2).*(cos(ph).*sin(th)) +...
       %     planeparams(3).*sin(ph)))  % Should the same as above,and it is.
        
        % Finding predictions from quadratic model
        predr2 = 1./(quadparams(1).*(cos(ph).*cos(th)).^2 +...
            quadparams(2).*(cos(ph).*sin(th)).^2 +...
            quadparams(3).*sin(ph).^2 + ...
            2*quadparams(4).*cos(ph).*cos(th).*cos(ph).*sin(th) +...
            2*quadparams(5).*cos(ph).*cos(th).*sin(ph) +...
            2*quadparams(6).*cos(ph).*sin(th).*sin(ph));
        
        predr2(predr2<0) = Inf;
        quadpred_r = sqrt(predr2);
        tmpdata = [tmpdata; planepred_r quadpred_r true_r];
    end
    planeerrors = (log10(tmpdata(:,1)) - log10(tmpdata(:,3))).^2;
    quaderrors = (log10(tmpdata(:,2)) - log10(tmpdata(:,3))).^2;
    p = signrank(planeerrors-quaderrors);
    data = [data; p]
end
sum(data <0.01)
% How many of these cells were rejected by the F-test? 
planelist = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\plane0.01.txt');
Lplane = zeros(length(fnames),1);
for cellcounter = 1:size(fnames,1)
    filename = char(fnames{cellcounter}(2));
    Lplane(cellcounter) = ismember(filename,[planelist{:}]');
end 
sum(Lplane & data < 0.01)
% The neurons for which a planar fit is rejected by this new test
% is a strict subset of the ones that are rejected by the F-test.
sum(data < 0.01)  % 39
sum(~Lplane)   % 65


%% 
% Section 22
% Finding the color of the stimulus used to measure orientation and SF
% tuning in the gratings paradigm.

data = [];
[fnames, spikeIdx] = fnamesFromTxt2();
for cellcounter = 1:size(fnames,1)
    GT = {};
    for i = 1:size(fnames{cellcounter},2)
        filename = findfile(char(fnames{cellcounter}(i)));
        paradigmID = getparadigmID(filename);
        if paradigmID == 150
            GT = nex2stro(filename);
        end
    end
    lidx = find(strcmp(GT.sum.trialFields(1,:),'lcont'));
    midx = find(strcmp(GT.sum.trialFields(1,:),'mcont'));
    sidx = find(strcmp(GT.sum.trialFields(1,:),'scont'));
    lms = GT.trial(:,[lidx midx sidx]);
    protocolidx = find(strcmp(GT.sum.trialFields(1,:),'protocol'));
    L = GT.trial(:,protocolidx) < 3;
    data = [data; mode(lms(L,:),1)];
end

Lachrom = (data(:,1) == data(:,2)) & (data(:,1) == data(:,3));
sum(Lachrom)