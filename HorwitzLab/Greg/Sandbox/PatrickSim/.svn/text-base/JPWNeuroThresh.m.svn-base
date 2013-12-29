%close all
%clear all
%%
% Getting the parameters for the quadratic model
filename = 'K082609010';
filename = 'K010610002';
%filename = 'K090309006';

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

%% User-defined variables
WHICHINITDIRS = 2;

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
    LMSdir(1,:)=[.01 0 0];
    LMSdir(2,:)=[0 .01 0];
    LMSdir(3,:)=[0 0 .1];    
end
%Specify staircasing variables
numReversals=7;

stepSize=0.5;
scale=0.5;

%Specify how many phases to be completed
phases=5;

%Specify number of grand loops
grandloops=100;


%% Non-user-defined variables
%load rig1mon
%global M 
%M=monspd*fundamentals;
%global bkgndLMS 
%bkgndLMS=M*bkgndrgb;


%% Start grand loop
quadparamsmatrix=nan(grandloops,6);
onecolumn=[];

for q=1:grandloops
clear trialspec


%% 3+ phases in 1 loop
[trialspec]=trialspecs(LMSdir,numReversals,thresh,stepSize,scale,phases,originalquadparams);


%% Plot these stimuli in color space and fit plane
[xformmat,scaled,Loog,planeparams,quadparams,SSE]=plotincolorspace(trialspec,thresh);
quadparamsmatrix(q,:)= quadparams;


%% End grand loop
end
figure; axes; hold on;
plot([-1 0 1 -1],[0 1 0 0],'k-');
colors = ['g', 'r', 'b'];
for q = 1:grandloops+1
    if (q == grandloops+1)
        quadparams = originalquadparams;        
    else
        quadparams = quadparamsmatrix(q,:);
    end
    A =  [quadparams(1) quadparams(4) quadparams(5);...
        quadparams(4) quadparams(2) quadparams(6);...
        quadparams(5) quadparams(6) quadparams(3)];
    [evecs,evals] = eig(A);
    [evals,i] = sort(diag(evals),1,'ascend');
    evecs = evecs(:,i);
    for j = 1:3
        v = evecs(:,j);
        if (v(2) < 0) v = -v; end
        h = plot(v(1)./sum(abs(v)),v(2)./sum(abs(v)),[colors(j),'s'],'MarkerFaceColor',colors(j));
        if (q == grandloops+1)
            set(h,'MarkerFaceColor',get(h,'MarkerFaceColor')/2)
        end
    end
end
for j = 1:3  % Plotting initial color directions
    v = LMSdir(j,:);
    if (v(2) < 0); v = -v; end
    plot(v(1)./sum(abs(v)),v(2)./sum(abs(v)),'k+','LineWidth',3);
end
set(gcf,'Name',[filename,': ',num2str(WHICHINITDIRS)]);

%% Display results

% Correlation Analysis - THIS IS LOCATED IN FILE NINEPOINTS.M
% correlationmean=mean(correlationMatrix);
% correlationstd=std(correlationMatrix);
% fprintf('The mean of the correlation coefficient distribution composed of %u samples is %f\n', q, correlationmean)
% fprintf('The standard deviation of the correlation coefficient distribution composed of %u samples is %f\n', q, correlationstd)


% Estimated error of plane fitting
% errorAngle=acosd((ConeWeights/norm(ConeWeights))*((planeparams'*xformmat')/norm(planeparams'*xformmat'))');
% errorMatrix=[errorMatrix, errorAngle];
% fprintf('The angle of the error in degrees is %f\n', errorAngle)


% SSE Neurothresh Model
% meanSSE=mean(SSEmatrix);
% stdSSE=std(SSEmatrix);
% figure; hold on;
% title('SSE Neurothresh Model')
% hist(SSEmatrix,100)
% legend(['Mean SSE =',num2str(meanSSE),   'Mean std =',num2str(stdSSE)]); 


% SSE Bootstrapping
% meanSSEmix=mean(mean(SSEmixmatrix)); 
% stdSSEmix=std(meanSSEmix);
% % for i=1:mixloops %Displays results for individual loops
% %     figure; hold on;
% %     title(['SSE Bootstrapping loop',num2str(i)])
% %     hist(SSEmixmatrix(:,i),100)
% %     legend(['Mean SSE =',num2str(meanSSEmix(i)),'Mean std =',num2str(stdSSEmix(i))],1); 
% % end
% figure; hold on;
% title('Combined SSE Bootstrapping loops')
% text(4,400,['Mean SSE =',num2str(meanmeanSSEmix)])
% text(4,350,['Mean std =',num2str(meanstdSSEmix)])
% hist(SSEmixmatrix,100)


% Display the means and standard deviations
% fprintf('The mean SSE composed of %u samples is %f\n', q, meanSSE)
% fprintf('The standard deviation of the SSEs composed of %u samples is %f\n', q, stdSSE)
% fprintf('The mean bootstrapped SSE composed of %u samples (each with %u shuffles)  is %f\n', q, nbootiter, meanmeanSSEmix)
% fprintf('The standard deviation of the bootstrapped SSEs composed of %u samples (each with %u shuffles) is %f\n', q, nbootiter, meanstdSSEmix)

% Display how well Planeparams predicts Cone Weights

%% Make ConeWeightsDist movie without actual ConeWeights

% figure; axes; hold on;
% plot3(ConeWeightsDist(:,1),ConeWeightsDist(:,2),ConeWeightsDist(:,3),'bo','MarkerFaceColor','blue')
% hold on; %grid on;
% %plot3(ConeWeights(1),ConeWeights(2),ConeWeights(3),'oy','MarkerFaceColor','yellow')
% set(gcf,'Color',[0 0 0]);
% set(gca,'XTick',[],'YTick',[],'ZTick',[],'Visible','off');
% %set(gca,'CameraViewAngleMode','manual')
% set(gca,'Position',[.15 .15 .7 .7]);
% %set(gca,'Color','none');
% set(gcf,'Position',[ 10   321   672   604]);
% %set(gcf,'Position',[100   500   400   0]);

% 
% viewangles = [180:3:360]+90;
% viewangles(end) = [];
% 
% clear M;
% for i = 1:length(viewangles)
%     axis vis3d
%     set(gca,'View',[viewangles(i) 22])
%     M(i) = getframe(gcf);
% end
% 
% repeat = 1;     %default = 1
% pSearch = 1;    %default = 0
% bSearch = 1;    %default = 1
% reference = 1;  %default = 0
% pixRange = 10;  %default = 10
% iFrame = 8;     %default = 8
% pFrame = 10;    %default = 10
% bFrame = 25;    %default = 25
% 
% options = [repeat, pSearch, bSearch, reference, pixRange, iFrame, pFrame, bFrame];
% mpgwrite(M, gray, 'ConeWeightsPlot1.mpg', options);


%% Make ConeWeightsDist movie with actual ConeWeights
% figure; axes; hold on;
% plot3(ConeWeightsDist(:,1),ConeWeightsDist(:,2),ConeWeightsDist(:,3),'bo','MarkerFaceColor','blue')
% hold on; %grid on;
% plot3(ConeWeights(1),ConeWeights(2),ConeWeights(3),'oy','MarkerFaceColor','yellow')
% set(gcf,'Color',[0 0 0]);
% set(gca,'XTick',[],'YTick',[],'ZTick',[],'Visible','off');
% %set(gca,'CameraViewAngleMode','manual')
% set(gca,'Position',[.15 .15 .7 .7]);
% %set(gca,'Color','none');
% set(gcf,'Position',[ 10   321   672   604]);
% %set(gcf,'Position',[100   500   400   0]);
% 
% 
% viewangles = [0:3:360]+90;
% viewangles(end) = [];
% 
% clear M;
% for i = 1:length(viewangles)
%     axis vis3d
%     set(gca,'View',[viewangles(i) 22])
%     M(i) = getframe(gcf);
% end
% 
% repeat = 1;     %default = 1
% pSearch = 1;    %default = 0
% bSearch = 1;    %default = 1
% reference = 1;  %default = 0
% pixRange = 10;  %default = 10
% iFrame = 8;     %default = 8
% pFrame = 10;    %default = 10
% bFrame = 25;    %default = 25
% 
% options = [repeat, pSearch, bSearch, reference, pixRange, iFrame, pFrame, bFrame];
% mpgwrite(M, gray, 'ConeWeightsPlot2.mpg', options);
% 
