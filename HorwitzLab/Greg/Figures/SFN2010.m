% Figures for the P Bio retreat talk 2009
%
% Section 1) 
% A rotating cube
%
% Section 2) 
% Cube with optional vector
%
% Section 3) 
% A contrast series of gabor images
%
% Section 4) 
% A color series of gabor images
%
% Section 5)
% Cube with colored gradient
%
% Section 6)
% Cube with a colored gradient and an embedded plane
%
% Section 7)
% Cube with a colored gradient at an orientation
%
% Section 8)
% Cube with a colored gradient at an orientation and an embedded 
% isoresponse plane.
%
% Section 9)
% Nonlinear 3-F function.
%
% Section 10)
% Transparent nonlinear 3-D function plus isoresponse surfaces.
%
% Section 11)
% An example stairase
%
% Section 12)
% Threshold estimates in normalized L,M,S space.  With meshs drawn.
%
% Section 13)
% Threshold estimates in normalized L,M,S space with locations of 
% grating stimuli drawn.
%
% Section 14)
% Transparent nonlinear 3-D function plus sample points.
%
% Section 15)
% Transparent LINEAR 3-D function with the same sample points as prev
% section.
%
% Section 16)
% Example monopolar stimuli
%
% Section 17)
% An example monopolar dataset. 
%
% Section 18) 
% A few Gabors that all have the same projection onto some vector


%%
% Section 1
% A rotating cube
ORIGIN = 0;
figure; axes; hold on;
set(gcf,'Color',[1 1 1]);
x = linspace(0,1,100);
plot3([0 1 1 0 0],[0 0 1 1 0],[0 0 0 0 0],'k-'); % face 1
plot3([0 1 1 0 0],[0 0 1 1 0],[1 1 1 1 1],'k-'); % face 2 
plot3([0 0],[0 0],[0 1],'k-'); % edge
plot3([0 0],[0 0],[0 1],'k-'); % edge
plot3([0 0],[1 1],[0 1],'k-'); % edge
plot3([1 1],[1 1],[0 1],'k-'); % edge
plot3([1 1],[0 0],[0 1],'k-'); % edge

if (ORIGIN)
    plot3([.5 .5],[.5 .5],[.5 .5],'ko','MarkerSize',3,'MarkerFaceColor','black'); % origin
end

set(gca,'Xlim',[-.2 1.2],'Ylim',[-.2 1.2],'Zlim',[-.2 1.2])
set(gcf,'Position',[131   537   430   430]);
set(gca,'XTick',[],'YTick',[],'ZTick',[],'Visible','off');
set(gca,'CameraViewAngleMode','manual')
%axis equal
set(gca,'View',[35 22])
viewangles = [0:3:520]+35;
viewangles(end) = [];

clear M;
for i = 1:length(viewangles)
    set(gca,'View',[viewangles(i) 22])
    M(i) = getframe(gcf);
end
repeat = 1;     %default = 1
pSearch = 1;    %default = 0
bSearch = 1;    %default = 1
reference = 1;  %default = 0
pixRange = 10;  %default = 10
iFrame = 8;     %default = 8
pFrame = 10;    %default = 10
bFrame = 25;    %default = 25
options = [repeat, pSearch, bSearch, reference, pixRange, iFrame, pFrame, bFrame];
mpgwrite(M, gray, 'rotcube.mpg', options);

%%
% Section 2
% Cube with optional vector
PLOTVECTOR = 1;
figure; axes; hold on;
set(gcf,'Color',[1 1 1]);
x = linspace(0,1,100);
plot3([0 1 1 0 0],[0 0 1 1 0],[0 0 0 0 0],'k-'); % face 1
plot3([0 1 1 0 0],[0 0 1 1 0],[1 1 1 1 1],'k-'); % face 2 
plot3([0 0],[0 0],[0 1],'k-'); % edge
plot3([0 0],[0 0],[0 1],'k-'); % edge
plot3([0 0],[1 1],[0 1],'k-'); % edge
plot3([1 1],[1 1],[0 1],'k-'); % edge
plot3([1 1],[0 0],[0 1],'k-'); % edge
plot3([.5 .5],[.5 .5],[.5 .5],'ko','MarkerSize',3,'MarkerFaceColor','black'); % origin

set(gca,'Xlim',[-.2 1.2],'Ylim',[-.2 1.2],'Zlim',[-.2 1.2])
set(gcf,'Position',[131   537   430   430]);
set(gca,'XTick',[],'YTick',[],'ZTick',[],'Visible','off');
set(gca,'CameraViewAngleMode','manual')
set(gca,'View',[90+35 22]);
%annotation('arrow',[0.5 .75],[0.5 .78])
if (PLOTVECTOR)
    a = [0 0 .5];
    h = quiver3(.5, .5, .5, a(1), a(2), a(3));
    set(h,'LineWidth',2,'Color','black');
    h = quiver3(.5, .5, .5, -a(1), -a(2), -a(3));
    set(h,'LineWidth',2,'Color','black');
end
%%
% Section 3
% A contrast series

nframes = 4;
bkgndrgb = [.5 .5 .5];
%gaborrgb = [.05 .05 .05];
gaborrgb = [-.05 .05 .05];

figure;
for i = 1:nframes
    subplot(2,nframes,i);
    im = DrawGaborEdge(bkgndrgb, gaborrgb*2*(i-1), [0 0 0], pi/2, 1.5, 1, 1, 0, 0, 0, 0, 0, .999, 45);
    image(im);

    set(gca,'XTick',[],'YTick',[],'visible','off');
    axis square;
end
for i = 1:nframes
    subplot(2,nframes,i+nframes);
    im = DrawGaborEdge(bkgndrgb, -gaborrgb*2*(i-1), [0 0 0], pi/2, 1.5, 1, 1, 0, 0, 0, 0, 0, .999, 45);
    image(im);

    set(gca,'XTick',[],'YTick',[],'visible','off');
    axis square;
end

%%
% Section 4
% A color series

bkgndrgb = [.5 .5 .5];
colordirections = [.09 .09 0;
    .09 -.09 0;
    0 0 0.6364;
    .1273 0 0;
    0 .1273 0;
    0.0636   -0.0636    -0.4500;
    0.0636   -0.0636    0.4500;
    0.0636   0.0636    0.4500;
    0.0636  0.0636    -0.4500;
    ];

% Stimuli that are metameric for K082509007
% colordirections = [[0.4171    0.7739    0.4765]*.3078;
%     [0.2351   -0.7120    0.6617]*0.1859;
%     [0.4080   -0.0336    0.9124]*0.1602;
%     [-0.5162   -0.8308    0.2079]*0.7149];
    

M = [0.0608    0.1219    0.0175;
    0.0220    0.1266    0.0257;
    0.0019    0.0095    0.0976];

bkgndlms = M*bkgndrgb';
for i = 1:2
    colordirections = -colordirections;
    stimlms = (colordirections+1).*repmat(bkgndlms',size(colordirections,1),1);
    gaborrgb = (inv(M)*stimlms')';
    figure;
    for j = 1:size(gaborrgb,1)
        subplot(ceil(sqrt(size(gaborrgb,1))),ceil(sqrt(size(gaborrgb,1))),j);
        im = DrawGaborEdge(bkgndrgb, gaborrgb(j,:)-bkgndrgb, [0 0 0], pi/2, 1.5, 1, 1, 0, 0, 0, 0, 0, .999, 45);
        image(im);
        
        set(gca,'XTick',[],'YTick',[],'visible','off');
        axis square;
    end
end
%%
% Section 5
% A rotating cube with a gradient of color
MONOPOLAR = 0;

if (~MONOPOLAR)
    fun = inline('2*abs(x-.5)');  % fullwave rectification
else
    fun = inline('max((1/.6)*(x-.4),0)');  % halfwave rectification
end

ALPHA = 1;
figure; axes; hold on;
set(gcf,'Color','none','Inverthardcopy','off');
set(gca,'Color','none');
x = linspace(0,1,200);

set(gca,'Xlim',[-.2 1.2],'Ylim',[-.2 1.2],'Zlim',[-.2 1.2])
set(gcf,'Position',[131   537   430   430]);
set(gca,'XTick',[],'YTick',[],'ZTick',[],'Visible','off');
set(gca,'CameraViewAngleMode','manual')
[X,Y] = meshgrid(x,x);

% Top and bottom
h = surf(X,Y,ones(size(X)),fun(X));
set(h,'EdgeColor','none','FaceAlpha',ALPHA);
h = surf(X,Y,zeros(size(X)),fun(X));
set(h,'EdgeColor','none','FaceAlpha',ALPHA);

% Front and back
h = surf(X, zeros(size(X)),Y,fun(X));
set(h,'EdgeColor','none','FaceAlpha',ALPHA);
h = surf(X, ones(size(X)),Y,fun(X));
set(h,'EdgeColor','none','FaceAlpha',ALPHA);

% Sides
if (~MONOPOLAR)
    h = surf(zeros(size(X)),Y,X,ones(size(X)));
else
    h = surf(zeros(size(X)),Y,X,zeros(size(X)));    
end
set(h,'EdgeColor','none','FaceAlpha',ALPHA);
h = surf(ones(size(X)),Y,X,ones(size(X)));
set(h,'EdgeColor','none','FaceAlpha',ALPHA);

set(gca,'View',[35 22]);
viewangles = [0:3:360]+35;
viewangles(end) = [];

clear M;
for i = 1:length(viewangles)
    set(gca,'View',[viewangles(i) 22])
    M(i) = getframe(gcf);
end

repeat = 1;     %default = 1
pSearch = 1;    %default = 0
bSearch = 1;    %default = 1
reference = 1;  %default = 0
pixRange = 10;  %default = 10
iFrame = 8;     %default = 8
pFrame = 10;    %default = 10
bFrame = 25;    %default = 25

options = [repeat, pSearch, bSearch, reference, pixRange, iFrame, pFrame, bFrame];
%mpgwrite(M, gray, 'rotgrad1.mpg', options);

figure;
colorbar;

%%
% Section 6
% A rotating cube with a gradient of color and an couple of embedded
% rectangle indicating isoresponse surfaces
MONOPOLAR = 1;

if (~MONOPOLAR)
    fun = inline('2*abs(x-.5)');  % fullwave rectification
else
    fun = inline('max((1/.6)*(x-.4),0)');  % halfwave rectification
end

ALPHA = .5;
figure; axes; hold on;
set(gcf,'Color','none','InvertHardCopy','off');
set(gca,'Xlim',[-.2 1.2],'Ylim',[-.2 1.2],'Zlim',[-.2 1.2],'color','none')
set(gcf,'Position',[131   537   430   430]);
set(gca,'XTick',[],'YTick',[],'ZTick',[],'Visible','off');
set(gca,'CameraViewAngleMode','manual')

x = linspace(0,1,200);
[X,Y] = meshgrid(x,x);

set(gca,'View',[35 22]);
viewangles = [0:3:720]+35;

clear M;
for i = 1:length(viewangles)
    cla;
    % Top and bottom
    h1 = surf(X,Y,ones(size(X)),fun(X));
    set(h1,'EdgeColor','none','FaceAlpha',.5);

    % Front and back
    if (mod(viewangles(i),360) < 90 | mod(viewangles(i),360) > 270)
        h3 = surf(X, zeros(size(X)),Y,fun(X));
        set(h3,'EdgeColor','none','FaceAlpha',ALPHA);
    else
        h4 = surf(X, ones(size(X)),Y,fun(X));
        set(h4,'EdgeColor','none','FaceAlpha',ALPHA);
    end
    
    % Sides
    if (mod(viewangles(i),360) < 180)
        h5 = surf(ones(size(X)),Y,X,ones(size(X)));
        set(h5,'EdgeColor','none','FaceAlpha',ALPHA);
    else
        if (~MONOPOLAR)
            h6 = surf(zeros(size(X)),Y,X,ones(size(X)));
        else
            h6 = surf(zeros(size(X)),Y,X,zeros(size(X)));
        end
        set(h6,'EdgeColor','none','FaceAlpha',ALPHA);
    end
    
    % The isoresponse planes
    if (~MONOPOLAR)
        p = patch([.25 .25 .25 .25],[0 1 1 0],[1 1 0 0],'w-');
    end
    p = patch([.75 .75 .75 .75],[0 1 1 0],[1 1 0 0],'w-');

    set(gca,'View',[viewangles(i) 22])
    
    M(i) = getframe(gcf);
end

repeat = 1;     %default = 1
pSearch = 1;    %default = 0
bSearch = 1;    %default = 1
reference = 1;  %default = 0
pixRange = 10;  %default = 10
iFrame = 8;     %default = 8
pFrame = 10;    %default = 10
bFrame = 25;    %default = 25
% 
options = [repeat, pSearch, bSearch, reference, pixRange, iFrame, pFrame, bFrame];
mpgwrite(M, gray, 'rotgrad2.mpg', options);

%%
% Section 7
% A rotating cube with a gradient of color at some kind of orientation.
MONOPOLAR = 1;
if (~MONOPOLAR)
    fun = inline('2*abs(x-.5)');  % fullwave rectification
else
    fun = inline('max((1/.6)*(x-.4),0)');  % halfwave rectification
end

ALPHA = 1;
figure; axes; hold on;
set(gcf,'Color','none');
x = linspace(0,1,100);

set(gca,'Xlim',[-.2 1.2],'Ylim',[-.2 1.2],'Zlim',[-.2 1.2])
set(gcf,'Position',[131   537   430   430]);
set(gca,'XTick',[],'YTick',[],'ZTick',[],'Visible','off');
set(gca,'CameraViewAngleMode','manual')
[X,Y] = meshgrid(x,x);

% Top and bottom
h = surf(X,Y,ones(size(X)),fun((X+Y+ones(size(X)))/3));
set(h,'EdgeColor','none','FaceAlpha',ALPHA);
h = surf(X,Y,zeros(size(X)),fun((X+Y+zeros(size(X)))/3));
set(h,'EdgeColor','none','FaceAlpha',ALPHA);

% Front and back
h = surf(X, zeros(size(X)),Y,fun((X+Y+zeros(size(X)))/3));
set(h,'EdgeColor','none','FaceAlpha',ALPHA);
h = surf(X, ones(size(X)),Y,fun((X+Y+ones(size(X)))/3));
set(h,'EdgeColor','none','FaceAlpha',ALPHA);

% Sides
h = surf(zeros(size(X)),Y,X,fun((X+Y+zeros(size(X)))/3));
set(h,'EdgeColor','none','FaceAlpha',ALPHA);
h = surf(ones(size(X)),Y,X,fun((X+Y+ones(size(X)))/3));
set(h,'EdgeColor','none','FaceAlpha',ALPHA);

set(gca,'View',[35 22]);
viewangles = [0:3:720]+35;

clear M;
for i = 1:length(viewangles)
    set(gca,'View',[viewangles(i) 22])
    M(i) = getframe(gcf);
end

repeat = 1;     %default = 1
pSearch = 1;    %default = 0
bSearch = 1;    %default = 1
reference = 1;  %default = 0
pixRange = 10;  %default = 10
iFrame = 8;     %default = 8
pFrame = 10;    %default = 10
bFrame = 25;    %default = 25

options = [repeat, pSearch, bSearch, reference, pixRange, iFrame, pFrame, bFrame];
%mpgwrite(M, gray, 'rotgrad3.mpg', options);

%%
% Section 8
% A rotating cube with a gradient of color at an orientation
% and an couple of embedded rectangles indicating isoresponse surfaces
MONOPOLAR = 1;
if (~MONOPOLAR)
    fun = inline('2*abs(x-.5)');  % fullwave rectification
else
    fun = inline('max((1/.6)*(x-.4),0)');  % halfwave rectification
end

ALPHA = .5;
figure; axes; hold on;
set(gcf,'Color','none');
set(gca,'Xlim',[-.2 1.2],'Ylim',[-.2 1.2],'Zlim',[-.2 1.2])
set(gcf,'Position',[131   537   430   430]);
set(gca,'XTick',[],'YTick',[],'ZTick',[],'Visible','off');
set(gca,'CameraViewAngleMode','manual')

x = linspace(0,1,100);
[X,Y] = meshgrid(x,x);

set(gca,'View',[35 22]);
viewangles = [0:3:720]+35;

clear M;
for i = 1:length(viewangles)
    cla;
    % Top only (bottom is never in front)
    h1 = surf(X,Y,ones(size(X)),fun((X+Y+ones(size(X)))/3));
    set(h1,'EdgeColor','none','FaceAlpha',ALPHA);

    % Front and back
    if (mod(viewangles(i),360) < 90 | mod(viewangles(i),360) > 270)
        h3 = surf(X, zeros(size(X)),Y,fun((X+Y+zeros(size(X)))/3));
        set(h3,'EdgeColor','none','FaceAlpha',ALPHA);
    else
        h4 = surf(X, ones(size(X)),Y,fun((X+Y+ones(size(X)))/3));
        set(h4,'EdgeColor','none','FaceAlpha',ALPHA);
    end
    
    % Sides
    if (mod(viewangles(i),360) < 180)
        h5 = surf(ones(size(X)),Y,X,fun((X+Y+ones(size(X)))/3));
        set(h5,'EdgeColor','none','FaceAlpha',ALPHA);    
    else
        h6 = surf(zeros(size(X)),Y,X,fun((X+Y+zeros(size(X)))/3));
        set(h6,'EdgeColor','none','FaceAlpha',ALPHA);
    end
    
    % The isoresponse planes
    patch([1  .2 1],[.2 1 1],[1 1 .2] ,'w-')
    if (~MONOPOLAR)
        patch([0  .8 0],[.8 0 0],[0 0 .8] ,'w-')
    end
    set(gca,'View',[viewangles(i) 22])
    
    M(i) = getframe(gcf);
end

repeat = 1;     %default = 1
pSearch = 1;    %default = 0
bSearch = 1;    %default = 1
reference = 1;  %default = 0
pixRange = 10;  %default = 10
iFrame = 8;     %default = 8
pFrame = 10;    %default = 10
bFrame = 25;    %default = 25

options = [repeat, pSearch, bSearch, reference, pixRange, iFrame, pFrame, bFrame];
mpgwrite(M, gray, 'rotgrad4.mpg', options);

%%
% Section 9
% A rotating cube with a pattern of color indicating nonlinear 
% combinations of signals.
% Does this even work?

ALPHA = .5;
theta = pi/4; % elev
phi = pi/4; % 4z

figure; axes; hold on;
set(gcf,'Color','none');
x = linspace(-1,1,200);

set(gca,'Xlim',[-1.2 1.2],'Ylim',[-1.2 1.2],'Zlim',[-1.2 1.2])
set(gcf,'Position',[131   537   430   430]);
set(gca,'XTick',[],'YTick',[],'ZTick',[],'Visible','off');
set(gca,'CameraViewAngleMode','manual')
[X,Y,Z] = meshgrid(x,x,x);
yprime = Y.*cos(theta)-Z.*sin(theta);
zprime = Y.*sin(theta)+Z.*cos(theta);
xprime = X.*cos(phi)-yprime.*sin(phi);
yprime = X.*sin(phi)+yprime.*cos(phi);

F = (xprime).^2-(yprime).^2-(zprime).^2;
%if (MONOPOLAR)
%    F(F<0) = 0;
%end


% Top and bottom
h = surf(squeeze(X(:,:,1)),squeeze(Y(:,:,1)),squeeze(Z(:,:,1)),squeeze(F(:,:,1)));
set(h,'EdgeColor','none','FaceAlpha',ALPHA);
h = surf(squeeze(X(:,:,end)),squeeze(Y(:,:,end)),squeeze(Z(:,:,end)),squeeze(F(:,:,end)));
set(h,'EdgeColor','none','FaceAlpha',ALPHA);

% Front and back
h = surf(squeeze(X(:,1,:)),squeeze(Y(:,1,:)),squeeze(Z(:,1,:)),squeeze(F(:,1,:)));
set(h,'EdgeColor','none','FaceAlpha',ALPHA);
h = surf(squeeze(X(:,end,:)),squeeze(Y(:,end,:)),squeeze(Z(:,end,:)),squeeze(F(:,end,:)));
set(h,'EdgeColor','none','FaceAlpha',ALPHA);

% Sides
h = surf(squeeze(X(1,:,:)),squeeze(Y(1,:,:)),squeeze(Z(1,:,:)),squeeze(F(1,:,:)));
set(h,'EdgeColor','none','FaceAlpha',ALPHA);
h = surf(squeeze(X(end,:,:)),squeeze(Y(end,:,:)),squeeze(Z(end,:,:)),squeeze(F(end,:,:)));
set(h,'EdgeColor','none','FaceAlpha',ALPHA);

set(gca,'View',[35 22]);

viewangles = [0:3:720]+35;
clear M;
for i = 1:length(viewangles)
    set(gca,'View',[viewangles(i) 22])
    M(i) = getframe(gcf);
end
% 
repeat = 1;     %default = 1
pSearch = 1;    %default = 0
bSearch = 1;    %default = 1
reference = 1;  %default = 0
pixRange = 10;  %default = 10
iFrame = 8;     %default = 8
pFrame = 10;    %default = 10
bFrame = 25;    %default = 25

options = [repeat, pSearch, bSearch, reference, pixRange, iFrame, pFrame, bFrame];
%mpgwrite(M, gray, 'rotgrad5.mpg', options);

%%
% Section 10
% A rotating cube with a pattern of color indicating nonlinear 
% combinations of signals.
MONOPOLAR = 0;  % Not implemented yet
ISOSURFACE = 1;

ALPHA = .5;
theta = pi/4; % elev
phi = pi/4; % 4z
%theta = 0; phi = 0;
viewangles = [0:3:720]+35;  % hot spot facing viewer
viewangles = [0:3:720]-59;  % viewer can see both hot spots

figure; axes; hold on;
set(gcf,'Color','none','InvertHardCopy','off');
x = linspace(-1,1,100);

set(gca,'Xlim',[-1.2 1.2],'Ylim',[-1.2 1.2],'Zlim',[-1.2 1.2])
set(gcf,'Position',[131   537   430   430]);
set(gca,'XTick',[],'YTick',[],'ZTick',[],'Visible','off');
set(gca,'CameraViewAngleMode','manual')
cmin = -2;
cmax = 2;  % yucky hardcoding
set(gca,'Clim',[cmin cmax]);

[X,Y,Z] = meshgrid(x,x,x);
yprime = Y.*cos(theta)-Z.*sin(theta);
zprime = Y.*sin(theta)+Z.*cos(theta);
xprime = X.*cos(phi)-yprime.*sin(phi);
yprime = X.*sin(phi)+yprime.*cos(phi);

F = (xprime).^2-(yprime).^2-(zprime).^2;
%F = F-min(F(:))/(max(F(:))-min(F(:)))
if (MONOPOLAR)
    F = F;

end

clear M;
for i = 1:length(viewangles)
    cla;
    % Top only (bottom is never in front)
    h = surf(squeeze(X(:,:,end)),squeeze(Y(:,:,end)),squeeze(Z(:,:,end)),squeeze(F(:,:,end)));
    set(h,'EdgeColor','none','FaceAlpha',ALPHA);
    
    % Front and back
    if (mod(viewangles(i),360) < 90 | mod(viewangles(i),360) > 270)
        h = surf(squeeze(X(1,:,:)),squeeze(Y(1,:,:)),squeeze(Z(1,:,:)),squeeze(F(1,:,:)));
        set(h,'EdgeColor','none','FaceAlpha',ALPHA);
    else
        h = surf(squeeze(X(end,:,:)),squeeze(Y(end,:,:)),squeeze(Z(end,:,:)),squeeze(F(end,:,:)));
        set(h,'EdgeColor','none','FaceAlpha',ALPHA);
    end
    
    % Sides
    if (mod(viewangles(i),360) < 180)
        h = surf(squeeze(X(:,end,:)),squeeze(Y(:,end,:)),squeeze(Z(:,end,:)),squeeze(F(:,end,:)));
        set(h,'EdgeColor','none','FaceAlpha',ALPHA);
    else
        h = surf(squeeze(X(:,1,:)),squeeze(Y(:,1,:)),squeeze(Z(:,1,:)),squeeze(F(:,1,:)));
        set(h,'EdgeColor','none','FaceAlpha',ALPHA);
    end
    
    % The isoresponse surfaces
    if (ISOSURFACE)
        iso= isosurface(X,Y,Z,F,.7);
        p = patch(iso);
        set(p,'EdgeColor', 'white','FaceColor','white','FaceAlpha',1);
    end
    set(gca,'View',[viewangles(i) 22])
    
    M(i) = getframe(gcf);
end
%
repeat = 1;     %default = 1
pSearch = 1;    %default = 0
bSearch = 1;    %default = 1
reference = 1;  %default = 0
pixRange = 10;  %default = 10
iFrame = 8;     %default = 8
pFrame = 10;    %default = 10
bFrame = 25;    %default = 25

options = [repeat, pSearch, bSearch, reference, pixRange, iFrame, pFrame, bFrame];
mpgwrite(M, gray, 'rotgrad6.mpg', options);

%%
% Section 11
% An example staircase

%filename = 'K082509007.nex';  % planar L-M+S cell
%whichidx = 39;  % which color index to use
%whichidx = 25;  % which color index to use
filename = 'K082509007.nex';  
whichidx = 4;  % which color index to use
%


stro = nex2stro(findfile(filename));
spikeidx = strcmp(stro.sum.rasterCells(1,:),getSpikenum(stro));
%Lbldone = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'bldone'));
coloridxs = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'coloridx'));
levels = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'levelidx'));
reversals = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'reversal'));
stepsize = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stepsize'));

lmsidxs = [find(strcmp(stro.sum.trialFields(1,:),'lcont'))...
            find(strcmp(stro.sum.trialFields(1,:),'mcont'))...
            find(strcmp(stro.sum.trialFields(1,:),'scont'))];
       
fpacq_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fp_acq'));
stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_on'));
stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_off'));
eot_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'eot'));

lms = stro.trial(:,lmsidxs);
lat = stro.sum.exptParams.latency;

Lbaseline = ~levels & all(lms == zeros(size(lms)),2);
Lprobe = ~levels & ~all(lms == zeros(size(lms)),2);

spikerates = zeros(size(stro.trial,1),1);
for i = 1:size(stro.trial,1)
    spiketimes = stro.ras{i,spikeidx};
    nspikes = sum(spiketimes > stimon_t(i)+stro.sum.exptParams.latency/1000 & spiketimes < stimoff_t(i));
    spikerates(i) = nspikes./(stimoff_t(i)-stimon_t(i)-stro.sum.exptParams.latency/1000);
end

[uniquecoloridxs,tmp1,tmp2] = unique(coloridxs);
uniquecolordirs = mkbasis(lms(tmp1,:)')';

% getting M matrix etc.

fundamentals = stro.sum.exptParams.fundamentals;
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = stro.sum.exptParams.mon_spd;
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;

% Calculating a "data" matrix which is a summary of the data
% to be used in subsequent cells of the script
data = [];
for i = uniquecoloridxs'
    L = coloridxs == i;
    rev = reversals(logical(L));
    if (any(stepsize(L) == 0))
        continue;
    end
    lev = unique(levels(L));
    [u,s,v] = svd(lms(L,:));
    unitvector = v(:,1);
    contrasts = lms(L,:)*unitvector;
    if(sum(contrasts) < 0)
        contrasts = -contrasts;
    end
    ss = stepsize(L);
    GV1 = GamutViolation(contrasts(end)*(1+ss(end)), unitvector, M, stro.sum.exptParams.bkgndrgb);
    GV2 = GamutViolation(contrasts(end)*(1-ss(end)), unitvector, M, stro.sum.exptParams.bkgndrgb);
    thresh = contrasts(end);
    %if (ss(end) > min(stepsize(stepsize > 0)))  % didn't finish for some reason
     if (sum(abs(rev)) < stro.sum.exptParams.nreversals)
        if (GV1)
            outofgamut = 1;
            thresh = thresh*(1+ss(end));
        elseif (GV2)
            outofgamut = 1;
            thresh = thresh*(1-ss(end));
        else
            outofgamut = nan;  % incomplete color direction
        end
    else
        outofgamut = 0;
    end
    if ~isnan(outofgamut)
        data = [data; i unitvector' thresh lev outofgamut];
    end
end
uniquecoloridxs = data(:,1);  % Ugly coding
data(:,1) = [];  % Ugly coding

%
stot = stimoff_t-stimon_t;
figure;
L = coloridxs == whichidx;
[u,s,v] = svd(lms(L,:));
unitvector = v(:,1);
if (sum(sign(unitvector)) < 0)
    unitvector = -unitvector;
end
contrasts = lms(L,:)*unitvector;
if (sum(contrasts < 0))
    contrasts = -contrasts;
end
contrastpeaks = contrasts(reversals(L) == 1);
contrasttroughs = contrasts(reversals(L) == -1);

% Now the actual plotting
subplot(3,2,1); hold on;
plot(contrasts,'ko-','MarkerSize',4,'MarkerFaceColor',uniquecolordirs(whichidx,:)/2+.5);
set(gca,'XLim',[0 length(contrasts)],'XTick',[]);

% Normal orientation rasters
subplot(3,2,2); hold on;
trialcounter = 1;
for j = find(L')
    sp = stro.ras{j,1}-stimon_t(j);
    h = plot([sp sp]',[-.3*ones(length(sp),1) .3*ones(length(sp),1)]'+trialcounter,'k-');
    plot([0 0],[-.5 .5]+trialcounter,'m-','Linewidth',2);
    plot([stro.sum.exptParams.latency stro.sum.exptParams.latency]/1000,[-.5 .5]+trialcounter,'m:','Linewidth',2);
    plot([stot(j) stot(i)],[-.5 .5]+trialcounter,'m-','Linewidth',2);
    trialcounter = trialcounter+1;
end
set(gca,'YLim',[-.5 trialcounter],'XLim',[-.1 max(stimoff_t-stimon_t)+.1]);

% Contrast as a function of trial number
subplot(3,2,3); hold on;
plot(spikerates(L),'ko-','MarkerSize',3,'Color',uniquecolordirs(whichidx,:)/4+.5);
plot([0 length(contrasts)], repmat(stro.sum.exptParams.threshold,1,2),'k:');
set(gca,'XLim',[0 length(contrasts)],'XTick',[]);
Loverthresh =  spikerates(L) > stro.sum.exptParams.threshold;
Lviolation = xor(Loverthresh(1:end-1)', sign(diff(contrasts))'== -1);
if any(Lviolation) % contrast trajectory differs from expected.  Probably due to bad isolation
    plot(find(Lviolation)+1,0,'rx');
end

subplot(3,2,5); hold on;
trialcounter = 1;
for j = find(L')
    sp = stro.ras{j,1}-stimon_t(j);
    h = plot([-.3*ones(length(sp),1) .3*ones(length(sp),1)]'+trialcounter, [sp sp]','k-');
    plot([-.5 .5]+trialcounter,[0 0],'m-','Linewidth',2);
    plot([-.5 .5]+trialcounter, [stro.sum.exptParams.latency stro.sum.exptParams.latency]/1000,'m:','Linewidth',2);
    plot([-.5 .5]+trialcounter,[stot(j) stot(i)],'m-','Linewidth',2);
    trialcounter = trialcounter+1;
end
set(gca,'XLim',[-.5 trialcounter],'YLim',[-.1 max(stimoff_t-stimon_t)+.1]);

subplot(3,2,4); hold on;
hist(contrasts);
if (~isempty(contrastpeaks))
    plot(contrastpeaks,3,'m^');
    plot(contrasttroughs,3,'bv');
end
plot(contrasts(end),3,'g*');

subplot(3,2,6); hold on;
plot(contrasts, spikerates(L),'k*');
plot([contrasts(end), contrasts(end)],[0 max(spikerates(L))],'k:');
plot([min(contrasts) max(contrasts)],[stro.sum.exptParams.threshold stro.sum.exptParams.threshold],'k:');
set(gcf,'Name',['condition ',num2str(whichidx)]);

%%
% Section 12
% 3-D distribution of threshold points.  Mesh optional.  GH: change to 10
% deg is incorrect as coded

PLOTMESH = 0;
PLOTFIRSTTHREE = 0;
MAKEMOVIE = 1;
PLOTPLANE = 0;
SHOWOOG = 1;
WHITEN = 1;
PLOTQUAD = 1;
PLOTELLIPSOID = 0;
USE10DEG = 1;

%filename = 'K082509007.nex';  % planar L-M+S cell
%filename = 'S041310002.nex';  % planar luminance cell
%filename = 'K072009003.nex';  % Luminance funnel cell
filename = 'K070309004.nex';  % Funnel intermediate color
%filename = 'K082609010.nex';  % pan color

if (strcmp(filename,'K082509007.nex'))
    STARTANGLE = 35;
    ENDANGLE = 678;
    ELEVVIEWANGLE = 22;
elseif(strcmp(filename,'K072009003.nex'))
    STARTANGLE = 228;
    ENDANGLE = 588;
    ELEVVIEWANGLE = 22;
elseif(strcmp(filename,'K070309004.nex'))
    STARTANGLE = -126;
    ENDANGLE = 360+STARTANGLE;
    ELEVVIEWANGLE = 16;
elseif(strcmp(filename,'K082609010.nex'))
    STARTANGLE = -30;
    ENDANGLE = 360+STARTANGLE;
    ELEVVIEWANGLE = 12;
end
stro = nex2stro(findfile(filename));
spikeidx = strcmp(stro.sum.rasterCells(1,:),getSpikenum(stro));
%Lbldone = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'bldone'));
coloridxs = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'coloridx'));
levels = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'levelidx'));
reversals = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'reversal'));
stepsize = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stepsize'));

lmsidxs = [find(strcmp(stro.sum.trialFields(1,:),'lcont'))...
            find(strcmp(stro.sum.trialFields(1,:),'mcont'))...
            find(strcmp(stro.sum.trialFields(1,:),'scont'))];
       
fpacq_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fp_acq'));
stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_on'));
stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_off'));
eot_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'eot'));

lms = stro.trial(:,lmsidxs);
lat = stro.sum.exptParams.latency;

Lbaseline = ~levels & all(lms == zeros(size(lms)),2);
Lprobe = ~levels & ~all(lms == zeros(size(lms)),2);

spikerates = zeros(size(stro.trial,1),1);
for i = 1:size(stro.trial,1)
    spiketimes = stro.ras{i,spikeidx};
    nspikes = sum(spiketimes > stimon_t(i)+stro.sum.exptParams.latency/1000 & spiketimes < stimoff_t(i));
    spikerates(i) = nspikes./(stimoff_t(i)-stimon_t(i)-stro.sum.exptParams.latency/1000);
end

[uniquecoloridxs,tmp1,tmp2] = unique(coloridxs);
uniquecolordirs = mkbasis(lms(tmp1,:)')';

% getting M matrix etc.
if (USE10DEG)
    load ('T_cones_smj10');
    fundamentals = T_cones_smj10';
else
    fundamentals = stro.sum.exptParams.fundamentals;
    fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
end
mon_spd = stro.sum.exptParams.mon_spd;
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;
M10 = T_cones_smj10*mon_spd;
M2to10 = M10/M;
    
% Calculating a "data" matrix which is a summary of the data
% to be used in subsequent cells of the script
data = [];
for i = uniquecoloridxs'
    L = coloridxs == i;
    rev = reversals(logical(L));
    if (any(stepsize(L) == 0))
        continue;
    end
    lev = unique(levels(L));
    [u,s,v] = svd(lms(L,:));
    unitvector = v(:,1);
    contrasts = lms(L,:)*unitvector;
    if(sum(contrasts) < 0)
        contrasts = -contrasts;
    end
    ss = stepsize(L);
    GV1 = GamutViolation(contrasts(end)*(1+ss(end)), unitvector, M, stro.sum.exptParams.bkgndrgb);
    GV2 = GamutViolation(contrasts(end)*(1-ss(end)), unitvector, M, stro.sum.exptParams.bkgndrgb);
    thresh = contrasts(end);
    %if (ss(end) > min(stepsize(stepsize > 0)))  % didn't finish for some reason
     if (sum(abs(rev)) < stro.sum.exptParams.nreversals)
        if (GV1)
            outofgamut = 1;
            thresh = thresh*(1+ss(end));
        elseif (GV2)
            outofgamut = 1;
            thresh = thresh*(1-ss(end));
        else
            outofgamut = nan;  % incomplete color direction
        end
    else
        outofgamut = 0;
    end
    if ~isnan(outofgamut)
        data = [data; i unitvector' thresh lev outofgamut];
    end
end
uniquecoloridxs = data(:,1);  % Ugly coding
data(:,1) = [];  % Ugly coding

% ---------------------------
% Filtering bad thresholds
% ---------------------------

RTHRESH = .4;  % threshold on rank contrast-response correlation.
tmpdata = [];
for i = uniquecoloridxs'
    L = coloridxs == i;
        
    [u,s,v] = svd(lms(L,:));
    unitvector = v(:,1);
    if (sum(sign(unitvector)) < 0)
        unitvector = -unitvector;
    end
    contrasts = lms(L,:)*unitvector;
    if (sum(contrasts < 0))
       contrasts = -contrasts; 
    end
    contrastpeaks = contrasts(reversals(L) == 1);
    contrasttroughs = contrasts(reversals(L) == -1);
    
    ss = stepsize(L);
    GV1 = GamutViolation(contrasts(end)*(1+ss(end)), unitvector, M, stro.sum.exptParams.bkgndrgb);
    GV2 = GamutViolation(contrasts(end)*(1-ss(end)), unitvector, M, stro.sum.exptParams.bkgndrgb);
    if (GV1 | GV2)
        outofgamut = 1;
    else
        outofgamut = 0;
    end
    
    Loverthresh =  spikerates(L) > stro.sum.exptParams.threshold;
    skipit = 0;
    if (length(contrasts) > 1)
        Lviolation = xor(Loverthresh(1:end-1)', sign(diff(contrasts))'== -1);
    % Contrast trajectory not in agreement with that calculated from
    % post-hoc sorted spikes
    else
        Lviolation = nan;
    end
    
    if (sum(stepsize(L)) == 0)
        skipit = 1;
    end
    if (ss(end) > min(stepsize(stepsize > 0)) & ~outofgamut)
        skipit = 1;
    end
    if (~skipit)
        r = corr([spikerates(L), contrasts],'type','Spearman');
        tmpdata = [tmpdata;  r(1,2), sum(contrastpeaks < contrasts(end)), sum(contrasttroughs > contrasts(end)), sum(Lviolation), outofgamut];
    end
end
% Columns:
% 1) Correlation between spikerate and contrast
% 2) Number of contrast peaks that exceed final contrast threshold estimate
% 3) Number of contrast troughs smaller than final contrast threshold estimate
% 4) Number of contrasts at which step direction disagrees with the predicted
% 5) Out of gamut

% Modify the Laccept vector which we use as a filter
% No longer destructively modifying the data matrix
Lfilter = logical(tmpdata(:,1) < RTHRESH & tmpdata(:,5) == 0); % Cor. filter if not oog
Lfilter = Lfilter | tmpdata(:,2) > 2;  % getting rid on unstable trajectories
Lfilter = Lfilter | tmpdata(:,4) > 0;  % getting rid of trajectory violations
Laccept = ~Lfilter;

% ---------------------------
% Sorting the data into two groups
% ---------------------------

Loog = data(Laccept,6) == 1;
azs = linspace(-pi,pi,60);
elevs = linspace(-pi/2, pi/2, 60);
oogmat = zeros(60,60);
noogmat = zeros(60,60);
tmp = data(Laccept,[1:3]);
xformed = tmp;
for i = 1:length(azs)
    for j = 1:length(elevs)
        az = azs(i);
        elev = elevs(j);
        rotmat1 = [cos(az) -sin(az) 0; sin(az) cos(az) 0; 0 0 1];
        rotmat2 = [1 0 0; 0 cos(elev) -sin(elev); 0 sin(elev) cos(elev)];
        tmp = rotmat1*[xformed; -xformed]';
        out = (rotmat2*tmp)';
        out = mkbasis(out')';
        tmpoog = out(logical([Loog; Loog]),3);
        tmpnoog = out(logical([~Loog; ~Loog]),3);
        tmpoog(tmpoog < 0) = [];
        tmpnoog(tmpnoog < 0) = [];
        
        tmp = tiedrank([tmpoog; tmpnoog]);
        oogmat(i,j) = sum(tmp(1:length(tmpoog)));
        noogmat(i,j) = sum(tmp(length(tmpoog)+1:end));
    end
end
ratio = (noogmat./oogmat);
[i,j] = ind2sub(size(ratio),find(ratio == max(ratio(:)),1));
bestaz = azs(i);
bestelev = elevs(j);

% Now just giving it a try
rotmat1 = [cos(bestaz) -sin(bestaz) 0; sin(bestaz) cos(bestaz) 0; 0 0 1];
rotmat2 = [1 0 0; 0 cos(bestelev) -sin(bestelev); 0 sin(bestelev) cos(bestelev)];
normalvect = inv(rotmat2*rotmat1)*[0 0 1]'; % could be in PCA space.
normalvect = normalvect./norm(normalvect);
out = rotmat2*rotmat1*xformed';

Lgrp1 = logical(data(:,[1 2 3])*normalvect > 0);  % normalvect is in LMS space
% Destructively modifying data so that all of the points are from 
% a single cluster.  Flip the signs to get the other cluster.
data(Lgrp1,[1 2 3]) = -data(Lgrp1,[1 2 3]);

% ---------------------------
% Plotting thresholds in normalized cone space
% ---------------------------
Loog = logical(data(:,6));
figure; axes; hold on;
set(gcf,'Color',[0 0 0],'InvertHardCopy','off');
set(gcf,'Position',[403   246   430   430]);
set(gca,'XTick',[],'YTick',[],'ZTick',[]);
set(gca,'Color',[0 0 0]);
%set(gca,'XColor','none','YColor','none','ZColor','none');

axis square; set(gca,'View',[35 22]);
set(gca,'XLim',[-1.01 1.01],'Ylim',[-1.01 1.01],'Zlim',[-1.01 1.01]);
set(gca,'CameraViewAngleMode','manual')

%plot3(1.2*[-1 1 1 -1 -1],1.2*[-1 -1 1 1 -1],1.2*[-1 -1 -1 -1 -1],'w-'); % face 1
%plot3(1.2*[-1 1 1 -1 -1],1.2*[-1 -1 1 1 -1],1.2*[1 1 1 1 1],'w-'); % face 2 
%lot3(1.2*[-1 -1],1.2*[-1 -1],1.2*[-1 1],'w-'); % edge
%plot3(1.2*[-1 -1],1.2*[-1 -1],1.2*[-1 1],'w-'); % edge
%plot3(1.2*[-1 -1],1.2*[1 1],1.2*[-1 1],'w-'); % edge
%plot3(1.2*[1 1],1.2*[1 1],1.2*[-1 1],'w-'); % edge
%plot3(1.2*[1 1],1.2*[-1 -1],1.2*[-1 1],'w-'); % edge

plot3([-1 1 1 -1 -1],[-1 -1 1 1 -1],[-1 -1 -1 -1 -1],'w-'); % face 1
plot3([-1 1 1 -1 -1],[-1 -1 1 1 -1],[1 1 1 1 1],'w-'); % face 2 
plot3([-1 -1],[-1 -1],[-1 1],'w-'); % edge
plot3([-1 -1],[-1 -1],[-1 1],'w-'); % edge
plot3([-1 -1],[1 1],[-1 1],'w-'); % edge
plot3([1 1],[1 1],[-1 1],'w-'); % edge
plot3([1 1],[-1 -1],[-1 1],'w-'); % edge


%xlabel('L'); set(get(gca,'XLabel'),'Color',[0 0 0]);
%ylabel('M'); set(get(gca,'YLabel'),'Color',[0 0 0]);
%zlabel('S'); set(get(gca,'ZLabel'),'Color',[0 0 0]);

if (WHITEN)
%    tmp = data(Laccept&~Loog,[1:3]) .*repmat(data(Laccept&~Loog,4), 1,3);
    tmp = data(~Loog,[1:3]) .*repmat(data(~Loog,4), 1,3);  
    [v,d] = eig(cov([tmp; -tmp]));
    scalemat = v*1/sqrt(d)/3;
    scaled = data(:,[1:3]) .*repmat(data(:,4), 1,3);
    scaled = scaled*scalemat;
else
%   scaled = data(Laccept,[1:3]) .*repmat(data(Laccept,4), 1,3);
    scaled = data(:,[1:3]) .*repmat(data(:,4), 1,3);
    scalefactors = max(abs([data(~Loog,[1 2 3]).*repmat(data(~Loog,4),1,3);...
    -data(~Loog,[1 2 3]).*repmat(data(~Loog,4),1,3)]));
    scaled = scaled./repmat(scalefactors, size(scaled,1),1);
end
% Shrinking OOG lines so that they fit in the cube (barely)
%scaled(Loog,:) = scaled(Loog,:)./repmat(max(abs(scaled(Loog,:)),[],2),1,3);
uniquelevs = unique(data(Laccept,5))';
if (PLOTFIRSTTHREE)
    uniquelevs = 1;
end

for i  = uniquelevs
    for j = 0:1  % Out of gamut
        L = logical(Laccept & data(:,5) == i & Loog == j);
        if(any(L))
            if (~j)
                plot3(scaled(L,1), scaled(L,2), scaled(L,3),'y.','MarkerSize',13);
                plot3(-scaled(L,1), -scaled(L,2), -scaled(L,3),'y.','MarkerSize',13);
            else
                if (SHOWOOG)
                    plot3([zeros(sum(L),1) scaled(L,1)]', [zeros(sum(L),1) scaled(L,2)]', [zeros(sum(L),1) scaled(L,3)]','r-');
                    plot3([zeros(sum(L),1) -scaled(L,1)]', [zeros(sum(L),1) -scaled(L,2)]', [zeros(sum(L),1) -scaled(L,3)]','r-');
                end
            end
        end
    end
end

if (PLOTMESH)
%    T = delaunay3(tmp(~Loog,1),tmp(~Loog,2),tmp(~Loog,3));
%    htm = tetramesh(T,[tmp(~Loog,1),tmp(~Loog,2),tmp(~Loog,3)], ones(length(T),1),'FaceAlpha', .1);
    tmpnoog = scaled(Laccept & ~Loog,:);
    K = convhulln(tmpnoog(:,[1 2 3]));
    convhullpts = unique(K(:));
    T = delaunay3(tmpnoog(convhullpts,1),tmpnoog(convhullpts,2),tmpnoog(convhullpts,3));
    htm{1} = tetramesh(T,[tmpnoog(convhullpts,1),tmpnoog(convhullpts,2),tmpnoog(convhullpts,3)], ones(length(T),1),'FaceAlpha', .2);
    set(htm{1},'EdgeAlpha',.01)
    htm{2} = tetramesh(T,[-tmpnoog(convhullpts,1),-tmpnoog(convhullpts,2),-tmpnoog(convhullpts,3)], ones(length(T),1),'FaceAlpha', .2);
    set(htm{2},'EdgeAlpha',.01)
end

if (PLOTPLANE || PLOTQUAD)
    [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);
    A = [quadparams(1) quadparams(4) quadparams(5);...
         quadparams(4) quadparams(2) quadparams(6);...
         quadparams(5) quadparams(6) quadparams(3)];
    coneweights = planeparams'*xformmat';
    B = xformmat*A*xformmat';
end
if (PLOTPLANE)
    plotlims = max(abs(scaled));
    [x,y,z] = meshgrid(linspace(-plotlims(1),plotlims(1),30),linspace(-plotlims(2),plotlims(2),30),linspace(-plotlims(3),plotlims(3),30));
    v = abs(x.*planeparams(1)+y.*planeparams(2)+z.*planeparams(3));
    fv = isosurface(x,y,z,v,1);
    h = patch(fv);
    set(h,'FaceColor','green','EdgeColor','none','FaceAlpha',.8);
end

if (PLOTQUAD)
    plotlims = max(abs(scaled))
    [x,y,z] = meshgrid(linspace(-plotlims(1),plotlims(1),30),linspace(-plotlims(2),plotlims(2),30),linspace(-plotlims(3),plotlims(3),30));
    xformedxyz = [x(:) y(:) z(:)];
    variables = [xformedxyz(:,1).^2 xformedxyz(:,2).^2 xformedxyz(:,3).^2 2*xformedxyz(:,1).*xformedxyz(:,2) 2*xformedxyz(:,1).*xformedxyz(:,3) 2*xformedxyz(:,2).*xformedxyz(:,3)];
    coefficients = [B(1,1) B(2,2) B(3,3) B(1,2) B(1,3) B(2,3)]';
    v = variables*coefficients;
    fv = isosurface(x,y,z,reshape(v,size(x)), 1);
    h = patch(fv);
    set(h,'FaceColor','green','EdgeColor','none','FaceAlpha',.8);
end

if (PLOTELLIPSOID)
   [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog, 'ellipsoid');
   A = [quadparams(1) quadparams(4) quadparams(5);...
       quadparams(4) quadparams(2) quadparams(6);...
       quadparams(5) quadparams(6) quadparams(3)];
   [evecs, evals] = eig(A);
   radii = sqrt(1./diag(evals));
   [x, y, z] = ellipsoid(0,0,0,radii(1),radii(2), radii(3));
   [newxyz] = [x(:) y(:) z(:)]*evecs'*inv(xformmat);
   
   
   h = surf(reshape(newxyz(:,1),size(x)), reshape(newxyz(:,2),size(y)), reshape(newxyz(:,3),size(z)));
   set(h,'FaceAlpha',.5,'FaceColor','green','Edgealpha',0);
   
end

if (MAKEMOVIE)
    clear M;
    if (PLOTFIRSTTHREE)
        viewangles = [[STARTANGLE:3:STARTANGLE+30],[STARTANGLE+30-3:-3:STARTANGLE]];        
    else
        viewangles = [STARTANGLE:3:ENDANGLE];
    end
    for i = 1:length(viewangles)
        set(gca,'View',[viewangles(i) ELEVVIEWANGLE])
        M(i) = getframe(gcf);
    end
    repeat = 1;     %default = 1
    pSearch = 1;    %default = 0
    bSearch = 1;    %default = 1
    reference = 1;  %default = 0
    pixRange = 10;  %default = 10
    iFrame = 8;     %default = 8
    pFrame = 10;    %default = 10
    bFrame = 25;    %default = 25
    options = [repeat, pSearch, bSearch, reference, pixRange, iFrame, pFrame, bFrame];
    mpgwrite(M, gray, 'rotdata.mpg', options);
end


%%
% Section 13
% A circulating movie of grating stimuli in cone space

STARTANGLE = -28;
ENDANGLE = 152;

figure; axes; hold on;
set(gcf,'Color',[0 0 0],'InvertHardCopy','off');
set(gca,'Xlim',[-.2 1.2],'Ylim',[-.2 1.2],'Zlim',[-.2 1.2])
set(gcf,'Position',[131   537   430   430]);
set(gca,'XTick',[],'YTick',[],'ZTick',[]);
set(gca,'Color',[0 0 0]);
axis square; set(gca,'View',[35 22]);
set(gca,'XLim',[-1.21 1.21],'Ylim',[-1.21 1.21],'Zlim',[-1.21 1.21]);
set(gca,'CameraViewAngleMode','manual')

colordirections = [.09 .09 0;
    .09 -.09 0;
    0 0 0.6364;
    .1273 0 0;
    0 .1273 0;
    0.0636   -0.0636    -0.4500;
    0.0636   -0.0636    0.4500;
    0.0636   0.0636    0.4500;
    0.0636  0.0636    -0.4500;
    ];
scalefactors = max(abs(colordirections));

for i = 1:size(colordirections,1)
    plot3(colordirections(i,1)./scalefactors(1),...
        colordirections(i,2)./scalefactors(2),...
        colordirections(i,3)./scalefactors(3),'w*','MarkerFaceColor','white','MarkerSize',3)
    plot3(-colordirections(i,1)./scalefactors(1),...
        -colordirections(i,2)./scalefactors(2),...
        -colordirections(i,3)./scalefactors(3),'w*','MarkerFaceColor','white','MarkerSize',3)
end
% Drawing the frame
plot3(1.2*[-1 1 1 -1 -1],1.2*[-1 -1 1 1 -1],1.2*[-1 -1 -1 -1 -1],'w-'); % face 1
plot3(1.2*[-1 1 1 -1 -1],1.2*[-1 -1 1 1 -1],1.2*[1 1 1 1 1],'w-'); % face 2 
plot3(1.2*[-1 -1],1.2*[-1 -1],1.2*[-1 1],'w-'); % edge
plot3(1.2*[-1 -1],1.2*[-1 -1],1.2*[-1 1],'w-'); % edge
plot3(1.2*[-1 -1],1.2*[1 1],1.2*[-1 1],'w-'); % edge
plot3(1.2*[1 1],1.2*[1 1],1.2*[-1 1],'w-'); % edge
plot3(1.2*[1 1],1.2*[-1 -1],1.2*[-1 1],'w-'); % edge

tmp = data(Laccept,[1:3]) .*repmat(data(Laccept,4), 1,3);
tmp = tmp./repmat(scalefactors, sum(Laccept),1);
[planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(tmp, Loog);
planeparams = planeparams'*xformmat';

tmp1 = [-1 1];
tmp2 = [-1 1];
[xx yy] = meshgrid(tmp1,tmp2);
%planez = -(1+planeparams(1)*xx+planeparams(2)*yy)/planeparams(3); % plane 1
planez = -(1+planeparams(2)*xx+planeparams(3)*yy)/planeparams(1); % plane 1

h = surf(planez,xx,yy);
set(h,'FaceAlpha',.5,'FaceColor','green');
planez =-(-1+planeparams(2)*xx+planeparams(3)*yy)/planeparams(1); % plane 1
h = surf(planez,xx,yy);
set(h,'FaceAlpha',.5,'FaceColor','green');


set(gca,'View',[ENDANGLE 22]);
set(gca,'visible','off')
 
 clear M;
 viewangles = [STARTANGLE:3:ENDANGLE];
 for i = 1:length(viewangles)
     set(gca,'View',[viewangles(i) 22])
     M(i) = getframe(gcf);
 end
 repeat = 1;     %default = 1
 pSearch = 1;    %default = 0
 bSearch = 1;    %default = 1
 reference = 1;  %default = 0
 pixRange = 10;  %default = 10
 iFrame = 8;     %default = 8
 pFrame = 10;    %default = 10
 bFrame = 25;    %default = 25
 options = [repeat, pSearch, bSearch, reference, pixRange, iFrame, pFrame, bFrame];
 mpgwrite(M, gray, 'rotdata.mpg', options);
 
%%
% Section 14
% A rotating cube with a pattern of color indicating nonlinear 
% combinations of signals.  Including a few sample points that 
% are supposed to represent the grating stimuli.  The cube will
% become transparent, but the points will not.

FADE = 1;
SHOWVOLUME = 1;
NROTATIONS = 1;
INITVIEWANGLE = 35 % hot spot facing viewer
INITVIEWANGLE = -59 % viewer can see both hot spots


theta = pi/4; % elev
phi = pi/4; % az
samplepoints = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 -1 0; 1 1 1; 1 1 -1; 1 -1 1; 1 -1 -1];
samplepoints = samplepoints./repmat(sqrt(sum(samplepoints.^2,2)),[1, size(samplepoints,2)]);
samplepoints = samplepoints/1.1;  % scaling down so we're not on the edge of the cube

figure; axes; hold on;
set(gcf,'Color',[0 0 0]);
set(gca,'Xlim',[-1.2 1.2],'Ylim',[-1.2 1.2],'Zlim',[-1.2 1.2])
set(gcf,'Position',[131   537   430   430]);
set(gca,'XTick',[],'YTick',[],'ZTick',[],'Visible','off');
set(gca,'CameraViewAngleMode','manual')
x = linspace(-1,1,200);
[X,Y,Z] = meshgrid(x,x,x);
yprime = Y.*cos(theta)-Z.*sin(theta);
zprime = Y.*sin(theta)+Z.*cos(theta);
xprime = X.*cos(phi)-yprime.*sin(phi);
yprime = X.*sin(phi)+yprime.*cos(phi);
F = (xprime).^2-(yprime).^2-(zprime).^2;
cmin= min(F(:)); cmax = max(F(:));
cmin = -2;
cmax = 2;  % yucky hardcoding
set(gca,'Clim',[cmin cmax]);

%Rotating samplepoints to find response only
yprime = samplepoints(:,2).*cos(theta)-samplepoints(:,3).*sin(theta);
zprime = samplepoints(:,2).*sin(theta)+samplepoints(:,3).*cos(theta);
xprime = samplepoints(:,1).*cos(phi)-yprime.*sin(phi);
yprime = samplepoints(:,1).*sin(phi)+yprime.*cos(phi);
Fsamplepoints = (xprime).^2-(yprime).^2-(zprime).^2;

% Plotting sample points in the appropriate color
cmap = colormap;
cm_length = length(colormap);
for i = 1:size(samplepoints, 1)
     h1 = plot3(samplepoints(i,1),samplepoints(i,2),samplepoints(i,3),'k.','MarkerSize',20);
     h2 = plot3(-samplepoints(i,1),-samplepoints(i,2),-samplepoints(i,3),'k.','MarkerSize',20);
     colormap_index = fix((Fsamplepoints(i)-cmin)/(cmax-cmin)*cm_length)+1
     set(h1,'Color',cmap(colormap_index,:))
     set(h2,'Color',cmap(colormap_index,:))
end

% Top and bottom
h(1) = surf(squeeze(X(:,:,1)),squeeze(Y(:,:,1)),squeeze(Z(:,:,1)),squeeze(F(:,:,1)));
h(2) = surf(squeeze(X(:,:,end)),squeeze(Y(:,:,end)),squeeze(Z(:,:,end)),squeeze(F(:,:,end)));

% Front and back
h(3) = surf(squeeze(X(:,1,:)),squeeze(Y(:,1,:)),squeeze(Z(:,1,:)),squeeze(F(:,1,:)));
h(4) = surf(squeeze(X(:,end,:)),squeeze(Y(:,end,:)),squeeze(Z(:,end,:)),squeeze(F(:,end,:)));

% Sides
h(5) = surf(squeeze(X(1,:,:)),squeeze(Y(1,:,:)),squeeze(Z(1,:,:)),squeeze(F(1,:,:)));
h(6) = surf(squeeze(X(end,:,:)),squeeze(Y(end,:,:)),squeeze(Z(end,:,:)),squeeze(F(end,:,:)));

viewangles = [0:3:360*NROTATIONS]+INITVIEWANGLE;
clear M;
for i = 1:length(viewangles)
    set(h,'EdgeColor','none')
    if (FADE == 1)
        set(h,'FaceAlpha',(i-1)/length(viewangles));
    elseif (FADE == -1)
        set(h,'FaceAlpha',(length(viewangles)+1-i)/length(viewangles));
    else
        set(h,'FaceAlpha',.5);
    end    
    set(gca,'View',[viewangles(i) 22])
    if (~SHOWVOLUME)
        set(h,'Visible','off');
    end
    M(i) = getframe(gcf);
end
% 
repeat = 1;     %default = 1
pSearch = 1;    %default = 0
bSearch = 1;    %default = 1
reference = 1;  %default = 0
pixRange = 10;  %default = 10
iFrame = 8;     %default = 8
pFrame = 10;    %default = 10
bFrame = 25;    %default = 25

options = [repeat, pSearch, bSearch, reference, pixRange, iFrame, pFrame, bFrame];
%mpgwrite(M, gray, 'rotgrad18.mpg', options);
%%
% Section 15
% A rotating cube with a pattern of color indicating linear 
% combination of signals followed by a static nonlinearity.
% Including the same sample points from above section.  The key is to 
% come up with an LN model cell that makes the same predictions as the 
% NL function above.  Set FADE to -1 to fade out the full volume, set
% FADE to 1 to fade in the whole volume.
INITVIEWANGLE = 35 % hot spot facing viewer
INITVIEWANGLE = -59 % viewer can see both hot spots
ISOSURFACE = 0;
FADE = 1;
PLOTSAMPLES = 1;

theta = pi/4; % elev
phi = pi/4; % az
yprime = sin(-phi);  % undoing the rotation of the space above
xprime = cos(-phi);
zprime = yprime*sin(-theta);
yprime = yprime*cos(-theta);
prefdir = [xprime yprime zprime];

samplepoints = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 -1 0; 1 1 1; 1 1 -1; 1 -1 1; 1 -1 -1];
samplepoints = samplepoints./repmat(sqrt(sum(samplepoints.^2,2)),[1, size(samplepoints,2)]);
samplepoints = samplepoints/1.1;  % scaling down so we're not on the edge of the cube
projections = samplepoints*prefdir';

% Getting responses to sample points
yprime = samplepoints(:,2).*cos(theta)-samplepoints(:,3).*sin(theta);
zprime = samplepoints(:,2).*sin(theta)+samplepoints(:,3).*cos(theta);
xprime = samplepoints(:,1).*cos(phi)-yprime.*sin(phi);
yprime = samplepoints(:,1).*sin(phi)+yprime.*cos(phi);
Fsamplepoints = (xprime).^2-(yprime).^2-(zprime).^2;

% Rotating samplepoints to find response only
% Sanity checking
%yprime = samplepoints(:,2).*cos(theta)-samplepoints(:,3).*sin(theta);
%zprime = samplepoints(:,2).*sin(theta)+samplepoints(:,3).*cos(theta);
%xprime = samplepoints(:,1).*cos(phi)-yprime.*sin(phi);
%yprime = samplepoints(:,1).*sin(phi)+yprime.*cos(phi);
%Fsamplepoints = (xprime).^2-(yprime).^2-(zprime).^2;

% Getting coefficients for nonlinearity: y = ax^2+bx+c
b = regress(Fsamplepoints,[projections.^2, projections, ones(length(projections),1)]);

% Setting up figure
figure; axes; hold on;
set(gcf,'Color',[0 0 0]);
x = linspace(-1,1,200);
set(gca,'Xlim',[-1.2 1.2],'Ylim',[-1.2 1.2],'Zlim',[-1.2 1.2])
set(gcf,'Position',[131   537   430   430]);
set(gca,'XTick',[],'YTick',[],'ZTick',[],'Visible','off');
set(gca,'CameraViewAngleMode','manual')
cmin = -2;
cmax = 2;  % yucky hardcoding

set(gca,'Clim',[cmin cmax]);

% Plotting sample points in the appropriate color
if (PLOTSAMPLES)
    cmap = colormap;
    cm_length = length(colormap);
    for i = 1:size(samplepoints, 1)
        h1 = plot3(samplepoints(i,1),samplepoints(i,2),samplepoints(i,3),'k.','MarkerSize',20);
        h2 = plot3(-samplepoints(i,1),-samplepoints(i,2),-samplepoints(i,3),'k.','MarkerSize',20);
        colormap_index = fix((Fsamplepoints(i)-cmin)/(cmax-cmin)*cm_length)+1;
        set(h1,'Color',cmap(colormap_index,:))
        set(h2,'Color',cmap(colormap_index,:))
    end
end
[X,Y,Z] = meshgrid(x,x,x);
linresp = (X*prefdir(1)+Y*prefdir(2)+Z*prefdir(3));
F = b(1)*linresp.^2+b(3);

% Top and bottom
h(1) = surf(squeeze(X(:,:,1)),squeeze(Y(:,:,1)),squeeze(Z(:,:,1)),squeeze(F(:,:,1)));
h(2) = surf(squeeze(X(:,:,end)),squeeze(Y(:,:,end)),squeeze(Z(:,:,end)),squeeze(F(:,:,end)));

% Front and back
h(3) = surf(squeeze(X(:,1,:)),squeeze(Y(:,1,:)),squeeze(Z(:,1,:)),squeeze(F(:,1,:)));
h(4) = surf(squeeze(X(:,end,:)),squeeze(Y(:,end,:)),squeeze(Z(:,end,:)),squeeze(F(:,end,:)));

% Sides
h(5) = surf(squeeze(X(1,:,:)),squeeze(Y(1,:,:)),squeeze(Z(1,:,:)),squeeze(F(1,:,:)));
h(6) = surf(squeeze(X(end,:,:)),squeeze(Y(end,:,:)),squeeze(Z(end,:,:)),squeeze(F(end,:,:)));

set(gca,'View',[35 22]);
set(h,'EdgeColor','none','FaceAlpha',1)


if (ISOSURFACE)
    iso= isosurface(X,Y,Z,F,.7);
    p = patch(iso);
    set(p,'EdgeColor', 'white','FaceColor','white','FaceAlpha',1);
end

viewangles = [0:3:360]+INITVIEWANGLE;

clear M;
for i = 1:length(viewangles)
    set(h,'EdgeColor','none');
    if (FADE == 1)
        set(h,'FaceAlpha',(i-1)/length(viewangles));
    elseif (FADE == -1)
        set(h,'FaceAlpha',(length(viewangles)+1-i)/length(viewangles));
    else
        set(h,'FaceAlpha',.5);
    end    
    set(gca,'View',[viewangles(i) 22])
    M(i) = getframe(gcf);
end
% 
repeat = 1;     %default = 1
pSearch = 1;    %default = 0
bSearch = 1;    %default = 1
reference = 1;  %default = 0
pixRange = 10;  %default = 10
iFrame = 8;     %default = 8
pFrame = 10;    %default = 10
bFrame = 25;    %default = 25

options = [repeat, pSearch, bSearch, reference, pixRange, iFrame, pFrame, bFrame];
%mpgwrite(M, gray, 'rotgrad20.mpg', options);

%%
% Section 16
% Example monopolar stimuli

bkgndrgb = [.5 .5 .5];
colordirections = [.09 .09 0;
    .09 -.09 0;
    0 0 0.6364;
    .1273 0 0;
    0 .1273 0;
    0.0636   -0.0636    -0.4500;
    0.0636   -0.0636    0.4500;
    0.0636   0.0636    0.4500;
    0.0636  0.0636    -0.4500;
    ];

M = [0.0608    0.1219    0.0175;
    0.0220    0.1266    0.0257;
    0.0019    0.0095    0.0976];

bkgndlms = M*bkgndrgb';
for i = 1:2
    colordirections = -colordirections;
    stimlms = (colordirections+1).*repmat(bkgndlms',size(colordirections,1),1);
    gaborrgb = (inv(M)*stimlms')'-repmat(bkgndrgb,size(stimlms,1),1);
    figure;
    for j = 1:size(gaborrgb,1)
        subplot(ceil(sqrt(size(gaborrgb,1))),ceil(sqrt(size(gaborrgb,1))),j);
        theta = pi/2;
        
        stimsizeindeg = norminv([1-.999 .999],0,1)*abs(1)/min([1 1]);
        stimsizeinpix = round(2*stimsizeindeg(2)*45);
        interval = linspace(stimsizeindeg(1), stimsizeindeg(2), stimsizeinpix);
        [X, Y] = meshgrid(interval,interval);
        xprime = X.*cos(-theta)+Y.*sin(-theta);
        yprime = -X.*sin(-theta)+Y.*cos(-theta);
        if (i == 1)
            fittedgabor = exp(-(xprime.^2.*yprime.^2)./(2.*2^2)).*cos(2.*pi.*yprime./1.5);
            fittedgabor(fittedgabor < 0) = 0;
        else
            fittedgabor = exp(-(xprime.^2.*yprime.^2)./(2.*2^2)).*cos(2.*pi.*yprime./1.5-pi);
            fittedgabor(fittedgabor < 0) = 0;
            
        end
        im = zeros(stimsizeinpix,stimsizeinpix,3);
        for plane = 1:3
            im(:,:,plane) = (fittedgabor.*gaborrgb(j,plane))+bkgndrgb(plane);
        end
        image(im);
        set(gca,'XTick',[],'YTick',[],'visible','off');
        axis square;
    end
end

%%
% Section 17
% Example monopolar dataset

PLOTMESH = 1;  % not supported
MAKEMOVIE = 1;
WHITEN = 1;
SHOWOOG = 1;

%filename = 'S090910004.nex';  % plane + extra points
filename = 'S091610002.nex';  % plane

if strcmp(filename, 'S090910004.nex')
    STARTANGLE = 129;
    ENDANGLE = 129+360;
    ELEVVIEWANGLE = 12;
elseif strcmp(filename,'S091610002.nex')
    STARTANGLE = 81;
    ENDANGLE = 81+360;
    ELEVVIEWANGLE = 14;
end
stro = nex2stro(findfile(filename));
out = NTpreprocess(stro,.4,2);
out(:,1) = [];

% getting M matrix etc.
fundamentals = stro.sum.exptParams.fundamentals;
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = stro.sum.exptParams.mon_spd;
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;
if (exist('M2to10'))
    M = M*M2to10;
end

% ---------------------------
% Plotting thresholds in normalized cone space
% ---------------------------
Loog = logical(out(:,6));
figure; axes; hold on;
set(gcf,'Color',[0 0 0],'InvertHardCopy','off');
set(gca,'Xlim',[-.2 1.2],'Ylim',[-.2 1.2],'Zlim',[-.2 1.2])
set(gcf,'Position',[131   537   430   430]);
set(gca,'XTick',[],'YTick',[],'ZTick',[]);
set(gca,'Color',[0 0 0]);
%set(gca,'XColor','none','YColor','none','ZColor','none');

axis square; set(gca,'View',[35 22]);
set(gca,'XLim',[-1.21 1.21],'Ylim',[-1.21 1.21],'Zlim',[-1.21 1.21]);
set(gca,'CameraViewAngleMode','manual')

plot3(1.2*[-1 1 1 -1 -1],1.2*[-1 -1 1 1 -1],1.2*[-1 -1 -1 -1 -1],'w-'); % face 1
plot3(1.2*[-1 1 1 -1 -1],1.2*[-1 -1 1 1 -1],1.2*[1 1 1 1 1],'w-'); % face 2 
plot3(1.2*[-1 -1],1.2*[-1 -1],1.2*[-1 1],'w-'); % edge
plot3(1.2*[-1 -1],1.2*[-1 -1],1.2*[-1 1],'w-'); % edge
plot3(1.2*[-1 -1],1.2*[1 1],1.2*[-1 1],'w-'); % edge
plot3(1.2*[1 1],1.2*[1 1],1.2*[-1 1],'w-'); % edge
plot3(1.2*[1 1],1.2*[-1 -1],1.2*[-1 1],'w-'); % edge

%xlabel('L'); set(get(gca,'XLabel'),'Color',[0 0 0]);
%ylabel('M'); set(get(gca,'YLabel'),'Color',[0 0 0]);
%zlabel('S'); set(get(gca,'ZLabel'),'Color',[0 0 0]);

if (WHITEN)
    tmp = out(~Loog,[1:3]) .*repmat(out(~Loog,4), 1,3);
    [v,d] = eig(cov([tmp; -tmp]));
    scalemat = v*1/sqrt(d)/2;
    scaled = out(:,[1 2 3]) .*repmat(out(:,4), 1,3);
    scaled = scaled*scalemat;
else
    scaled = out(:,[1:3]) .*repmat(out(:,4), 1,3);
    scalefactors = max(abs([out(~Loog,[1 2 3]).*repmat(out(~Loog,4),1,3);...
    -out(~Loog,[1 2 3]).*repmat(out(~Loog,4),1,3)]));
    scaled = scaled./repmat(scalefactors, size(out,1),1);
end
uniquelevs = unique(out(:,5))';
for i  = uniquelevs
    for j = 0:1  % Out of gamut
        L = logical(out(:,5) == i & Loog == j);
        if(any(L))
            if (~j)
                plot3(scaled(L,1), scaled(L,2), scaled(L,3),'y.','MarkerSize',13);
            else
                if (SHOWOOG)
                    plot3([zeros(sum(L),1) scaled(L,1)]', [zeros(sum(L),1) scaled(L,2)]', [zeros(sum(L),1) scaled(L,3)]','y-');
                end
            end
        end
    end
end

if (PLOTMESH)
%    T = delaunay3(tmp(~Loog,1),tmp(~Loog,2),tmp(~Loog,3));
%    htm = tetramesh(T,[tmp(~Loog,1),tmp(~Loog,2),tmp(~Loog,3)], ones(length(T),1),'FaceAlpha', .1);
    tmpnoog = scaled(~Loog,:);
    K = convhulln(tmpnoog(:,[1 2 3]));
    convhullpts = unique(K(:));
    T = delaunay3(tmpnoog(convhullpts,1),tmpnoog(convhullpts,2),tmpnoog(convhullpts,3));
    htm{1} = tetramesh(T,[tmpnoog(convhullpts,1),tmpnoog(convhullpts,2),tmpnoog(convhullpts,3)], ones(length(T),1),'FaceAlpha', .2);
    set(htm{1},'EdgeAlpha',.01)
end


if (MAKEMOVIE)
    clear M;
    viewangles = [STARTANGLE:3:ENDANGLE];
    for i = 1:length(viewangles)
        set(gca,'View',[viewangles(i) ELEVVIEWANGLE])
        M(i) = getframe(gcf);
    end
    repeat = 1;     %default = 1
    pSearch = 1;    %default = 0
    bSearch = 1;    %default = 1
    reference = 1;  %default = 0
    pixRange = 10;  %default = 10
    iFrame = 8;     %default = 8
    pFrame = 10;    %default = 10
    bFrame = 25;    %default = 25
    options = [repeat, pSearch, bSearch, reference, pixRange, iFrame, pFrame, bFrame];
    mpgwrite(M, gray, 'rotdata.mpg', options);
end


%%
% Section 18
% A few Gabors that all have the same projection onto some color vector

bkgndrgb = [.5 .5 .5];
colordirections = [.05*ones(5,1) normrnd(0,.1,5,2)];

M = [0.0608    0.1219    0.0175;
    0.0220    0.1266    0.0257;
    0.0019    0.0095    0.0976];

bkgndlms = M*bkgndrgb';
for i = 1:2
    colordirections = -colordirections;
    stimlms = (colordirections+1).*repmat(bkgndlms',size(colordirections,1),1);
    gaborrgb = (inv(M)*stimlms')'-repmat(bkgndrgb,size(stimlms,1),1);
    figure;
    for j = 1:size(gaborrgb,1)
        subplot(ceil(sqrt(size(gaborrgb,1))),ceil(sqrt(size(gaborrgb,1))),j);
        theta = pi/2;
        
        stimsizeindeg = norminv([1-.999 .999],0,1)*abs(1)/min([1 1]);
        stimsizeinpix = round(2*stimsizeindeg(2)*45);
        interval = linspace(stimsizeindeg(1), stimsizeindeg(2), stimsizeinpix);
        [X, Y] = meshgrid(interval,interval);
        xprime = X.*cos(-theta)+Y.*sin(-theta);
        yprime = -X.*sin(-theta)+Y.*cos(-theta);
        if (i == 1)
            fittedgabor = exp(-(xprime.^2.*yprime.^2)./(2.*2^2)).*cos(2.*pi.*yprime./1.5);
        else
            fittedgabor = exp(-(xprime.^2.*yprime.^2)./(2.*2^2)).*cos(2.*pi.*yprime./1.5-pi);            
        end
        im = zeros(stimsizeinpix,stimsizeinpix,3);
        for plane = 1:3
            im(:,:,plane) = (fittedgabor.*gaborrgb(j,plane))+bkgndrgb(plane);
        end
        image(im);
        set(gca,'XTick',[],'YTick',[],'visible','off');
        axis square;
    end
end

