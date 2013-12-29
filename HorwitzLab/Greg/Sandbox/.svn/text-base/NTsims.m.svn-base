% Looking at isoresponse surfaces
% of various quadratics.
% It is indeed true that 1- and 2-sheet hyperboloids are closely related.
% They differ only in the sign of a scale factor.

[x y z] = meshgrid([-1:.1:1],[-1:.1:1],[-1:.1:1]);
fn = -x.^2-y.^2+z.^2;

figure;
fnvals = [-.1 0 .1];
for i = 1:length(fnvals)
    subplot(3,1,i)
    v = isosurface(x,y,z,fn,fnvals(i));
    patch(v);
    axis square
end

%%
% Plotting two "scaled versions" of the same hyperbolic surface, but with a
% negative scale factor. Can we turn a 1-sheet hyperboloid into a 2-sheet?
% Yes, they are the same, but with a negative scalefactor

[x y z] = meshgrid([-3:.3:3],[-3:.3:3],[-3:.3:3]);
coefs = [1 1 -1];

figure; axes; hold on;
fn = coefs(1)*x.^2+coefs(2)*y.^2+coefs(3)*z.^2;
v = isosurface(x,y,z,fn,1);
h = patch(v);
set(h,'FaceColor',[1 0 0],'FaceAlpha',.5);

scalefactor = -1;
fn = scalefactor*(coefs(1)*x.^2+coefs(2)*y.^2+coefs(3)*z.^2);
v = isosurface(x,y,z,fn,1);
h = patch(v);
set(h,'FaceColor',[0 1 0],'FaceAlpha',.5);




%%
% Simulations for Reviewer #3
% In this simulation, we start with activity in 3 postreceptoral channels.
% We weight signals from these postreceptoral channels randomly, add them
% together, and square them. The resulting isoresponse surfaces are planes
% (because there is only one nonlinear operation which happens at the very
% end of the chain of events) that are randomly oriented.
[l m s] = meshgrid([-1:.1:1],[-1:.1:1],[-1:.1:1]);
ch1 = l+m;
ch2 = l-m;
ch3 = s-.1*(l+m);
figure;
for i = 1:9
    weights = unifrnd(-1,1,1,3);     
    fn = (weights(1).*ch1+weights(2).*ch2+weights(3).*ch3).^2;
    subplot(3,3,i)
    v = isosurface(l,m,s,fn,.5);
    patch(v);
    axis square
end

% Now doing a bunch of simulations and extracting the principal axes by
% brute force.
principalaxes = [];
for i = 1:100
    i
    weights = unifrnd(-1,1,1,3);
    fn = (weights(1).*ch1+weights(2).*ch2+weights(3).*ch3).^2;
    v = isosurface(l,m,s,fn,.5);
    scaled = v.vertices;
    if (isempty(scaled))
        continue;
    end
    scaled = scaled(round(linspace(1,length(scaled),10)),:);
    Loog = zeros(size(scaled,1),1);
    [planeparams, ~, quadparams, ~, xformmat] = NTsurfacefit(scaled, Loog);
    coneweights = (planeparams'*xformmat')./sum(abs(planeparams'*xformmat'));
    if (coneweights(2) < 0)   % convention: positive M-cone weight
        coneweights = -coneweights;
    end
    principalaxes = [principalaxes; coneweights];
end

figure; axes; hold on;
LnegS = principalaxes(:,3) < 0;
plot(principalaxes(~LnegS,1),principalaxes(~LnegS,2),'o','MarkerSize',7,'MarkerFaceColor','black','MarkerEdgeColor','black');
plot(principalaxes(LnegS,1),principalaxes(LnegS,2),'o','MarkerSize',7,'MarkerFaceColor','white','MarkerEdgeColor','black');
plot([-1 1 0 -1],[0 0 1 0],'k-');

%%
% Simulations for Reviewer #3
% In this simulation, we start with activity in 3 postreceptoral channels.
% We weight signals from these postreceptoral channels randomly and add them.
% We then square them, weight them, and add them.
% Unlike simualation #1, the resulting isoresponse surfaces are quadratic 
% (because the squaring happens before the find summing). They are randomly
% oriented because of the initial random weighted sum of the postreceptoral
% signals.

[l m s] = meshgrid([-3:.1:3],[-3:.1:3],[-3:.1:3]);
ch01 = l+m;
ch02 = l-m;
ch03 = s-.1*(l+m);

figure
for i = 1:9
    % First, a linear combination of postreceptoral signals
    weights = unifrnd(-1,1,3,3);
 %   weights = .5*weights./repmat(sqrt(sum(weights.^2,2)),1,3)
    
    ch11 = weights(1,1)*ch01+weights(1,2)*ch02+weights(1,3)*ch03;
    ch12 = weights(2,1)*ch01+weights(2,2)*ch02+weights(2,3)*ch03;
    ch13 = weights(3,1)*ch01+weights(3,2)*ch02+weights(3,3)*ch03;

    % Now squaring and adding
    weights = unifrnd(-1,1,1,3);
%    weights = .5*weights./repmat(sqrt(sum(weights.^2)),1,3);
    fn = weights(1)*ch11.^2+...
         weights(2)*ch12.^2+...
         weights(3)*ch13.^2;

    subplot(3,3,i)
    v = isosurface(l,m,s,fn,1);
    patch(v);
    axis square; set(gca,'XLim',[-3 3],'YLim',[-3 3],'ZLim',[-3 3])
end


% Now doing a bunch of simulations and extracting the principal axes by
% brute force.
principalaxes = []; eigenvalues = []; surfacetypes = [];
for i = 1:200
    i
    % First, a linear combination of postreceptoral signals
    weights = unifrnd(-1,1,3,3);
    
    ch11 = weights(1,1)*ch01+weights(1,2)*ch02+weights(1,3)*ch03;
    ch12 = weights(2,1)*ch01+weights(2,2)*ch02+weights(2,3)*ch03;
    ch13 = weights(3,1)*ch01+weights(3,2)*ch02+weights(3,3)*ch03;

    % Now squaring and adding
    weights = unifrnd(-1,1,1,3);
    fn = weights(1)*ch11.^2+...
         weights(2)*ch12.^2+...
         weights(3)*ch13.^2;
    v = isosurface(l,m,s,fn,.5);
    scaled = v.vertices;
    if (isempty(scaled))
        continue;
    end
    scaled = scaled(round(linspace(1,length(scaled),10)),:);
    Loog = zeros(size(scaled,1),1);
    [planeparams, ~, quadparams, ~, xformmat] = NTsurfacefit(scaled, Loog);
    
    A =  [quadparams(1) quadparams(4) quadparams(5);...
        quadparams(4) quadparams(2) quadparams(6);...
        quadparams(5) quadparams(6) quadparams(3)];
    % Transforming quadparams into un-whitened space
    [evecs,evals] = eig(xformmat*A*xformmat');
    [evals,idx] = sort(diag(evals),1,'ascend');
    evecs = evecs(:,idx);
    
    principalaxes = cat(3,principalaxes,evecs);
    eigenvalues = cat(2,eigenvalues,evals);
    surfacetypes = [surfacetypes; sum(evals < 0)];
end

evalratiothresh = 5;
for whichsurf = 0:2; % 0 = ellipsoid, 1 = 1 sheet, 2 = 2 sheets
    out = [];
    for i = 1:size(eigenvalues,2)
        if (surfacetypes(i) == whichsurf)
            out = [out; reshape(principalaxes(:,:,i),1,9)];
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
    figure; axes; hold on;
    symbols = {'gs','rd','bv'};
    titles = {'Ellipsoid','Hyperboloid 1 sheet','Hyperboloid 2 sheets'};
    for i = 1:3
        idxs = [1 2 3]+(i-1)*3;
        tmp = out(:,idxs)./repmat(sum(abs(out(:,idxs)),2),1,3);
        tmp = tmp.*repmat(sign(tmp(:,2)),1,3);
        plot(tmp(:,1),tmp(:,2),symbols{i},'MarkerFaceColor',symbols{i}(1),'MarkerSize',5);
    end
    plot([-1 1 0 -1],[0 0 1 0],'k-');
end

%%
% An intro-style figure for Reviewer 3. Showing that the same results can be
% obtained from two different "broadly-tuned" neurons that have very
% different isoresponse surfaces.

a = 10;
b = 30;
teststimuli = [cos(linspace(0,2*pi,9)); sin(linspace(0,2*pi,9))]';
teststimuli(end,:) = [];
v1 = [1 2]';
v1 = v1./norm(v1);
v2 = [-2 1]';
v2 = v2./norm(v2);

LM = teststimuli*v1;
LvM = teststimuli*v2;
fr = a.*LM.^2+b*LvM.^2;

figure; axes; hold on;
for i = 1:length(fr)
   plot(teststimuli(i,1),teststimuli(i,2),'ko','MarkerSize',fr(i),'MarkerFaceColor','black');
end
axis square;

figure; axes; hold on;
v = [2 -1]';
plot(abs(teststimuli*v),fr,'ko');
regcoefs = regress(fr,[ones(length(fr),1) abs(teststimuli*v).^2]);

[x,y] = meshgrid([-1.5:.1:1.5],[-1.5:.1:1.5]);
fn1 = regcoefs(1)+regcoefs(2)*abs([x(:) y(:)]*v).^2;
LM = [x(:) y(:)]*v1;
LvM = [x(:) y(:)]*v2;
fn2 = a.*LM.^2+b*LvM.^2;
fnmax = max([fn1(:); fn2(:)]);

figure; subplot(2,1,1); hold on; axis square;
surf(x,y,zeros(size(x)),reshape(fn1,size(x,1),size(x,2))./fnmax);
set(gca,'XLim',[min(x(:)) max(x(:))],'Ylim',[min(y(:)) max(y(:))],'XTick',[],'YTick',[])
caxis([0,1]);

subplot(2,1,2); hold on; axis square; 
surf(x,y,zeros(size(x)),reshape(fn2,size(x,1),size(x,2))./fnmax);
set(gca,'XLim',[min(x(:)) max(x(:))],'Ylim',[min(y(:)) max(y(:))],'XTick',[],'YTick',[])
caxis([0,1]); 
colormap(hsv);

for i = 1:2
    subplot(2,1,i);
    for j = 1:length(fr)
        plot(teststimuli(j,1),teststimuli(j,2),'ko','MarkerSize',.5*fr(j),'MarkerFaceColor','black');
    end
end

%%
% Sums of sinewaves
% What would the isoresponse surface of a linear neuron with a
% spatiochroamtic inseparable RF look like?
% 2-D simulation
w1 = [.1 1];
w2 = [-1 .5];
delay = 3;
stim_template = sin(linspace(0,20,100));
stim_template = stim_template.*linspace(0,1,100);
RF1 = diff(normpdf(linspace(-12,12,101),0,3));
RF2 = diff(normpdf(linspace(-12,12,101)+delay,0,1));
figure;
subplot(2,1,1); hold on;
plot(RF1,'b-');
plot(RF2,'g-');

[x,y] = meshgrid(linspace(-1,1,10), linspace(-1,1,10));
responses = zeros(size(x));

for i = 1:length(x(:));
  %  resp1 = (x(i)*stim_template) .* (w1(1)*RF1) + (y(i)*stim_template) .* (w1(2)*RF1);
  %  resp2 = (x(i)*stim_template) .* (w2(1)*RF2) + (y(i)*stim_template) .* (w2(2)*RF2);
    resp1 = conv(x(i)*stim_template,(w1(1)*RF1))+conv(y(i)*stim_template,(w1(2)*RF1));
    resp2 = conv(x(i)*stim_template,(w2(1)*RF2))+conv(y(i)*stim_template,(w2(2)*RF2));
    responses(i) = sum(resp1+resp2);
end
subplot(2,1,2);
imagesc(responses)

%%
% Space-time inseparable neurons don't have null positions
% for stationary gratings

x = mkgrating(40,40,pi/5,10,pi/6)
envelope = normpdf(linspace(-6,6,40),0,1);
RF = (envelope'*envelope).*x;
RF = zeros(40,40);
RF([3:20],20) = 1;
RF([3:20],21) = 1;
RF([3:20],22) = 1;
RF([10:30],23) = -1;
RF([10:30],24) = -1;
RF([10:30],25) = -1;

figure; 
imagesc(RF);
colormap(gray); 

% testing responses to a stationary grating
phases = linspace(0,2*pi,100);
resp =[];
for i = 1:length(phases)
    tmpresp = zeros(40,1);
    for j = 1:40
        stim = mkgrating(40,40,pi/2,10,phases(i));
        start = j;
        stop = min(40,j+20);
        stim([start:stop],:) = -stim([start:stop],:); % Contrast reversal
        tmpresp(j) = sum(sum(stim.*RF));
    end
    resp(i) = max(tmpresp);
  resp(i) = sum(sum(stim.*RF));
end
%plot(resp);
    
figure;axes; hold on;
L = phases<pi;
plot(phases(L),resp(L));
plot(phases(~L)-pi,resp(~L),'g-')
%%
% A firing rate function can be separable in 'r' and 'theta' and not have
% isoresponse contours that are scalar multiples of each other.

theta = linspace(0,2*pi,80);
r = linspace(0,5,20);
[rs,ths] = meshgrid(r, theta);

x = rs.*cos(ths);
y = rs.*sin(ths);
%plot(x,y,'k.');

fr = max(min(rs,3),0).*sin(2*ths)+5;

figure; axes; hold on;
surf(x,y,fr)
contour(x,y,fr+1,13)