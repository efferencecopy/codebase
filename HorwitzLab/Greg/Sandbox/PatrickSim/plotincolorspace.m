function [xformmat,scaled,Loog,planeparams,quadparams,SSE]=plotincolorspace(trialspec,~)

global ConeWeights

figure(2),clf, axes, hold on
set(gcf,'Position',[ 10   21   672   604]);
for i=1:length(trialspec)
    if ~isempty(trialspec(i).oog) && trialspec(i).oog==0 && trialspec(i).inqueue==0
        for n=[-1 1]
            if trialspec(i).round==1
                plot3(n*trialspec(i).coordinates(1),n*trialspec(i).coordinates(2),n*trialspec(i).coordinates(3),'ob')
            end
            if trialspec(i).round==2
                plot3(n*trialspec(i).coordinates(1),n*trialspec(i).coordinates(2),n*trialspec(i).coordinates(3),'og')
                plot3(n*[0 trialspec(i).coordinates(1)], n*[0 trialspec(i).coordinates(2)], n*[0 trialspec(i).coordinates(3)],'g:')
            end
            if trialspec(i).round>=3
%                 plot3(n*trialspec(i).parenttriangle([2:3,1],1),n*trialspec(i).parenttriangle([2:3,1],2),n*trialspec(i).parenttriangle([2:3,1],3),'g--')
                plot3(n*trialspec(i).coordinates(1),n*trialspec(i).coordinates(2),n*trialspec(i).coordinates(3),'oc')
                plot3(n*[0 trialspec(i).coordinates(1)], n*[0 trialspec(i).coordinates(2)], n*[0 trialspec(i).coordinates(3)],'c:')
            end
        end
    end
end

%scaledConeWeights=1.3*ConeWeights*(min([trialspec.vectorlength])/norm(ConeWeights));
%plot3([-scaledConeWeights(1) scaledConeWeights(1)],[-scaledConeWeights(2) scaledConeWeights(2)],[-scaledConeWeights(3) scaledConeWeights(3)]) %plot ABC as sanity check

%grid on
%axis normal
%rotate3d

%% Make tetrahedron movie without planes
% 
% set(gcf,'Color',[0 0 0]);
% set(gca,'XTick',[],'YTick',[],'ZTick',[],'Visible','off');
% %set(gca,'CameraViewAngleMode','manual')
% set(gca,'Position',[.3 .3 .4 .4]);
% %set(gca,'Color','none');
% set(gcf,'Position',[ 10   321   672   604]);
% %set(gcf,'Position',[100   500   400   0]);
% 
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
% mpgwrite(M, gray, 'TetrahedronNoPlanes.mpg', options);
% 
% 
% %% Make Tetrahedron movie with planes
% 
scaled=cat(1,trialspec.coordinates);
Loog=[trialspec.oog]';
[planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);
SSE=(planeSSE/quadSSE);
A =  [quadparams(1) quadparams(4) quadparams(5);...
    quadparams(4) quadparams(2) quadparams(6);...
    quadparams(5) quadparams(6) quadparams(3)];
% Transforming quadparams into un-whitened space
B = xformmat*A*xformmat';
quadparams = [B(1,1) B(2,2) B(3,3) B(1,2) B(1,3) B(2,3)]; 
 
% planeparams=planeparams*sign((planeparams'*xformmat')*ConeWeights');
% x=-2:.1:2;
% y=-2:.1:2;
% [X,Y]=meshgrid(x,y);
% Z = (-1-planeparams(1)*X -planeparams(2)*Y)/planeparams(3);
% tmp = permute(cat(3,X,Y,Z),[3 1 2]);
% tmp = reshape(tmp, 3, size(x,2)*size(y,2));
% xformed = tmp'*inv(xformmat);
% h1 = mesh(reshape(xformed(:,1),41,41),reshape(xformed(:,2),41,41),reshape(xformed(:,3),41,41));
% h2 = mesh(reshape(-xformed(:,1),41,41),reshape(-xformed(:,2),41,41),reshape(-xformed(:,3),41,41));
% hidden off
% %figure(2);rotate3d;
% 
% 
% %figure(2); axes; hold on;
% set(gcf,'Color',[0 0 0]);
% set(gca,'XTick',[],'YTick',[],'ZTick',[],'Visible','off');
% %set(gca,'CameraViewAngleMode','manual')
% set(gca,'Position',[.01 .01 .98 .98]);
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
% mpgwrite(M, gray, 'TetrahedronWithPlanes.mpg', options);
