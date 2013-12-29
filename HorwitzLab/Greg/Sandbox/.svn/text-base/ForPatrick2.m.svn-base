% Playing around with radial basis functions

n = 200;
nbasis = n;
x = unifrnd(0,10,n,1);
y = unifrnd(0,10,n,1);
z = poissrnd(x+(y-5).^2);
%z = poissrnd(x+y);
%%
%plot3(x,y,z,'k.');

sigma = 10000;
%X = [];
%for i = 1:nbasis
%    X = [X,normpdf(x,x(i),sigma).*normpdf(y,y(i),sigma)];
%end

X = [];
for i = 1:nbasis
    X = [X,exp(-(x-x(i)).^2/sigma).*exp(-(y-y(i)).^2/sigma)];
end


% least-squares
%b = regress(z,X);
%preds = zeros(n,1);
%for i = 1:nbasis
%    preds = preds+b(i).*normpdf(x,x(i),sigma).*normpdf(y,y(i),sigma);
%end
%(preds - z)
% Not surprisingly, it fits every point perfectly.

% Now trying ridge regression
%b = ridge(z,X,.1);
%preds = zeros(n,1);
%for i = 1:n
%    preds = preds+b(i).*normpdf(x,x(i),sigma).*normpdf(y,y(i),sigma);
%end
% Didn't work

% Let's see if I can use glmfit with these radial basis functions...
[b,dev,stats] = glmfit(X,z,'poisson','link','identity','constant','off');
preds = zeros(n,1);
%for i = 1:nbasis
  %  preds = preds+b(i).*normpdf(x,x(i),sigma).*normpdf(y,y(i),sigma);
%     preds = preds+b(i).*exp(-(x-x(i)).^2/sigma).*exp(-(y-y(i)).^2/sigma);
%end
preds = X*b;
(preds - z)
figure
plot(z,preds-z,'k.');

% well, I guess we can plot the surface...
%tmp = linspace(0,10,30);
%surface = zeros(length(tmp));
%for i = 1:length(tmp)
%    for j = 1:length(tmp)
%        pred = 0;
%        for k = 1:nbasis
%            pred = pred+b(k).*exp(-(tmp(i)-x(k)).^2/sigma).*exp(-(tmp(j)-y(k)).^2/sigma);
%        end
%        surface(j,i) = pred;
%    end
%end

% vectorizing
[tmpx,tmpy] = meshgrid(linspace(0,10,30));
surface = 0;
for k = 1:nbasis
    surface = surface+b(k).*exp(-(tmpx-x(k)).^2/sigma).*exp(-(tmpy-y(k)).^2/sigma);
end

figure; axes; hold on;
surf(tmpx, tmpy, surface);
plot3(x,y,z,'o');


%This is working, but it's way too bumpy a surface
% I could either get rid of basis function or add a ridge term. Can I do
% this with GLM? Or just make sigma realy huge.


%%
% Use CV to find the optimal sigma.
sigmas = logspace(2,12,30);
errs = nan*ones(length(sigmas),n);
for i = 1:length(sigmas)  % sigma loop
    X = [];
    for k = 1:n
%        X = [X,normpdf(x,x(k),sigmas(i)).*normpdf(y,y(k),sigmas(i))];
        X = [X,exp(-(x-x(k)).^2/sigmas(i)).*exp(-(y-y(k)).^2/sigmas(i))];
    end
    for j = 1:n  % which point we're leaving out
        tmpX = X;
        tmpx = x;
        tmpy = y;
        tmpz = z;
        tmpx(j) = [];
        tmpy(j) = [];
        tmpz(j) = [];
        tmpX(:,j) = [];
        tmpX(j,:) = [];
       % [b,dev,stats] = glmfit(tmpX,tmpz,'poisson','link','identity','constant','off');
        b = regress(tmpz,tmpX); % LS to go faster

        pred = 0;
        for k = 1:length(b)  % could use glmval for this
          %  pred = pred+b(k).*normpdf(x(j),tmpx(k),sigmas(i)).*normpdf(y(j),tmpy(k),sigmas(i));
            pred = pred+b(k).*exp(-(x(j)-tmpx(k)).^2/sigmas(i)).*exp(-(y(j)-tmpy(k)).^2/sigmas(i));
        end
       
        errs(i,j) = (pred-z(j)).^2;
    end
end
figure; subplot(2,1,1);
imagesc(log(errs));
subplot(2,1,2);
plot(sigmas,sum(errs'),'.-');
toterr = sum(errs');
prefsigma = sigmas(toterr == min(toterr))
set(gca,'XScale','log','YScale','log');
sigma = prefsigma

% Plotting data and fit with prefsigma

X = [];
for i = 1:nbasis
    X = [X,exp(-(x-x(i)).^2/sigma).*exp(-(y-y(i)).^2/sigma)];
end
[b,dev,stats] = glmfit(X,z,'poisson','link','identity','constant','off');
%b = regress(z,X); % LS to go faster

tmp = linspace(0,10,30);
surface = zeros(length(tmp));
for i = 1:length(tmp)
    for j = 1:length(tmp)
        pred = 0;
        for k = 1:nbasis
            %pred = pred+b(k).*normpdf(tmp(i),x(k),sigma).*normpdf(tmp(j),y(k),sigma);
            pred = pred+b(k).*exp(-(tmp(i)-x(k)).^2/sigma).*exp(-(tmp(j)-y(k)).^2/sigma);
        end
        surface(j,i) = pred;
    end
end
%surface = flipud(surface);
figure; axes; hold on;
surf(tmp, tmp, surface);
plot3(x,y,z,'o');

% How many basis functions do I need to model Patrick's data so far? Where
% should they be placed?

%%
% Trying to make a uniform B-spline
% 'm' in wikipedia is # of knots
% 'n' is the degree of the spline

x = linspace(0,1,100);
nknots = 4;
knots = linspace(x(1),x(end),nknots)
ti = 0;
X = [ones(length(x),1) x',x'.^2,x'.^3];
mat = [-1 3 -3 1; 3 -6 3 0; 3 0 3 0; 1 4 1 0];
P = [-1 -1 0 0]';
greg = X*(1/6)*mat*P
plot(greg)

%%
% trying a zero-order spline
lambda = [5 20 10];
projs = repmat([1 2 3],10,1);
projs = projs(:);
y = [poissrnd(lambda(1),10,1),poissrnd(lambda(2),10,1),poissrnd(lambda(3),10,1)];
y = y(:);
X = [projs == 1, projs == 2, projs == 3];
[b,dev] = glmfit(X,y,'poisson','link','identity','constant','off');
figure;axes; hold on;
plot(projs,y,'.')
plot([1 2 3], b,'g*')
%%
% Trying a first-order B-spline
n = 60;
maxfr = 100;
nknots = 10;
projs = unifrnd(0,maxfr,n,1);
fr= poissrnd(.09*projs+10);
knots = linspace(0,maxfr,nknots);
knots = [knots(1)-(knots(end)-knots(end-1)), knots, knots(end)+(knots(end)-knots(end-1))];
%knots = projs;
X = zeros(n,length(knots)-2);
for i = 1:length(knots)-2
    tmp1 = (projs-knots(i))/(knots(i+1)-knots(i));
    tmp2 = (knots(i+2)-projs)/(knots(i+2)-knots(i+1));
    L1 = projs>=knots(i) & projs < knots(i+1);
    L2 = projs>=knots(i+1) & projs < knots(i+2);
    X(L1,i) = tmp1(L1);
    X(L2,i) = tmp2(L2);
end
[b,dev] = glmfit(X,fr,'poisson','link','identity','constant','off');

figure;axes; hold on;
plot(projs,fr,'.');
x = linspace(0,maxfr);
X = zeros(length(x),length(knots)-2);
for i = 1:length(knots)-2
    tmp1 = (x-knots(i))/(knots(i+1)-knots(i));
    tmp2 = (knots(i+2)-x)/(knots(i+2)-knots(i+1));
    L1 = x>=knots(i) & x < knots(i+1);
    L2 = x>=knots(i+1) & x < knots(i+2);
    X(L1,i) = tmp1(L1);
    X(L2,i) = tmp2(L2);
end
plot(X*b,'k-');
%%
% Trying a quadratic B-spline (upgrading to quadratic and cubic)
% Piecewise linear fit needs three knots per basis function
% Piecewise quadratic fit needs four.
% Will this help glmfits convergence?
n = 300
maxfr = 100;
nknots = 30;
projs = unifrnd(0,maxfr,n,1);
fr= poissrnd(3*projs+1);
%knots = linspace(0,maxfr,nknots);
%knots = [knots(1)-(knots(end)-knots(end-1)), knots, knots(end)+(knots(end)-knots(end-1))];
knots = linspace(-50,150,nknots);  % pretty arbitrary
X2 = zeros(n,length(knots)-3);
X3 = zeros(n,length(knots)-4);
figure; subplot(2,1,1); hold on;
for i = 1:length(knots)-3
    sp = spmak(knots(i:i+3),1);
    X2(:,i)= spval(sp,projs)
    fnplt(sp);
end
subplot(2,1,2); hold on;
for i = 1:length(knots)-4
    sp = spmak(knots(i:i+4),1);
    X3(:,i)= spval(sp,projs)
    fnplt(sp);
end

[b2,dev] = glmfit(X2,fr,'poisson','link','identity','constant','off');
[b3,dev] = glmfit(X3,fr,'poisson','link','identity','constant','off');
% Trying least-squares with the same basis functions (makes only a very
% small difference, but does make a difference).
bls = regress(fr,X3);

figure;axes; hold on;
plot(projs,fr,'.');
x = linspace(0,maxfr);
X2 = zeros(length(x),length(knots)-3);
for i = 1:length(knots)-3
    sp = spmak(knots(i:i+3),1);
    X2(:,i)= spval(sp,x)
end
X3 = zeros(length(x),length(knots)-4);
for i = 1:length(knots)-4
    sp = spmak(knots(i:i+4),1);
    X3(:,i)= spval(sp,x)
end
plot(X2*b2,'m-','linewidth',3);
plot(X3*bls,'c-','linewidth',3);
plot(X3*b3,'y-','linewidth',3);

% weird edge effects. Curve fit goes down too much (because basis functions
% go to 0)
%%
% Plotting a single quadratic B-spline from wikipedia
% Seeing if it's the same as what I get from spmak.
% It is. Good.
knts = [0 1 2 3]
x = linspace(knts(1),knts(end),100);
L1= x > knts(1) & x <= knts(2);
L2= x > knts(2) & x <= knts(3);
L3= x > knts(3) & x <= knts(4);
y = zeros(length(x),1);
y(L1) = .5*(x(L1)-knts(1)).^2;
y(L2) = -(x(L2)-knts(2)).^2+(x(L2)-knts(2))+.5;
y(L3) = .5*(1-(x(L3)-knts(3))).^2;
figure; axes; hold on;
sp = spmak(knts,1);
fnplt(sp);
plot(x,y,'m:');

%%
% Trying to estimate a (single) Poisson spike rate using a Bayesian
% technique. This works.

lambda = 50;
alpha = 1;
maxfr = 200;
beta = maxfr/(2*alpha);
x = linspace(0,maxfr,500);

figure; axes; hold on;
plot(x,gampdf(x,alpha,beta))
plot(lambda,0, 'm*')
% mean is alpha*beta, variance is alpha*beta^2
a = alpha;
b = beta;
for i = 1:10
    pause
    draw = poissrnd(lambda);
    %draw = lambda + round(normrnd(0,sqrt(lambda),1,1));
    a = a+draw;
    b = 1./(1+(1/b));
 
    cla;
    plot(lambda,0,'m*')
    plot(draw,0,'k*')
    plot(x,gampdf(x,a,b))
end

%%
% Trying a 1-D case with Gaussian basis functions 
% Use CV to estimate optimal sigma (smoothing)
% How do I propagate the uncertainty? Linear propagation of a and b?
alpha = 1;
maxfr = 200;
maxcontrast = 10;
beta = maxfr/(2*alpha);
as = alpha*ones(50,1);
bs = beta*ones(50,1);
idxs = linspace(0,maxcontrast,50)';

neuron = inline('100./(1+.2*exp(-(3*x-10)))','x');
%neuron = inline('50*normpdf(x,5,2)','x');

defaultsigma = .3;

figure;
subplot(2,2,1);
plot(idxs,neuron(idxs),'k.');
subplot(2,2,2);
plot(idxs,as);
subplot(2,2,3);
plot(idxs,bs);
subplot(2,2,4);
plot(idxs,as.*bs.^2);

sigma = defaultsigma;
stims = [];
resps = [];
for i = 1:100
    vars = as.*bs.^2;
    contrast = idxs(find(vars == max(vars),1));
    resp = poissrnd(neuron(contrast));
    %resp = neuron(contrast)+normrnd(0,2*0);
    stims = [stims; contrast];
    resps = [resps; resp];
    weighting = exp(-(idxs-contrast).^2/sigma);
    
    % updating the as and bs
    as = as+weighting*resp;
    bs = 1./(weighting+1./bs);
    
    subplot(2,2,1); hold on;
    plot(contrast,resp,'m*')
    subplot(2,2,2);
    plot(idxs,weighting);
    subplot(2,2,3);
    plot(idxs,as.*bs);
    title('mean of lambda')
    subplot(2,2,4);
    plot(idxs,as.*bs.^2);
    title('variance of lambda')
    pause;

     if (i > 8)
%          nbasis = 5;
%          sigmas = linspace(.1,10,10);
%          nCVreps = 5;
%          err = zeros(nCVreps,length(sigmas));
%          basiscenters = linspace(0,maxcontrast,nbasis);
%          for k = 1:length(sigmas)
%              X = [];
%              for j = 1:nbasis % inefficient
%                  X = [X,exp(-(stims-basiscenters(j)).^2/sigmas(k))];
%              end
%              L = logical([0; ones(length(stims)-1,1)]);
%              for j = 1:nCVreps % arbitrary, ten CVs
%                  L = L(randperm(length(L)));
%                  b = X(L,:)\resps(L)
%                  
%                  pred = zeros(sum(~L),1);
%                  for m = 1:nbasis % inefficient
%                      pred = pred + exp(-(stims(~L)-basiscenters(j)).^2/sigmas(k))*b(m);
%                  end
%                  k
%                  err(j,k) = sum((pred-resps(~L)).^2); % Make neg likelihood sometime
%              end
%          end
%          toterr = mean(err);
% %        % plot(toterr);
%          sigma = sigmas(find(toterr == min(toterr)))
% %         sigma = 3;
% 
% 
% % Here's how to do kernel regression
% X = [];
% for i = 1:nbasis
%     X = [X,exp(-(stims-basiscenters(i)).^2/sigma)];
% end
% b = X\resps
% 
% fnest = zeros(length(idxs),1);
% for j = 1:nbasis
%     fnest = fnest+ b(j)*exp(-(idxs-basiscenters(j)).^2/sigma);
% end
%          
%           subplot(2,2,2);
%          plot(idxs,fnest,'b.')
%          title(['sigma = ',num2str(sigma)]);
%          drawnow
%          pause
%    %     keyboard
%     % Need to recalculate a and b here...
%     end


% Instead of incrementally changing as and bs, recalculating them from
% scratch with an arbitrary sigma

%data = [];
%sigmas = linspace(.1,3,20);
%for testsigma = sigmas
 %   as = alpha*ones(50,1);
 %   bs = beta*ones(50,1);
%    for i = 1:length(stims)
%        weighting = exp(-(idxs-stims(i)).^2/testsigma);
        % updating the as and bs
%        as = as+weighting*resps(i);
 %       bs = 1./(weighting+1./bs);
 %   end
 %   mn = as.*bs;
 %   var = as.*bs.^2;
 %   r = corrcoef([mn var]);  % This doesn't make much sense
 %  figure; subplot(2,1,1); plot(as.*bs); subplot(2,1,2); plot(as.*bs.^2);
 %ose all
% title (['r = ',num2str(r(1,2))]);
 %   data = [data; r(1,2)];
%end
%sigma = sigmas(data == max(data))  % Setting it.
sigma = unifrnd(.1, 5);
as = alpha*ones(50,1);
bs = beta*ones(50,1);
for i = 1:length(stims)
    weighting = exp(-(idxs-stims(i)).^2/sigma);
    % updating the as and bs
    as = as+weighting*resps(i);
    bs = 1./(weighting+1./bs);
end
     end
end


%%
% Starting over from above using a spline and using the derivative of the
% spline to propagate alpha and beta. Worked better even thought it's
% ad-hoc. Need a better way of propagating information.
alpha = 1;
maxfr = 200;
maxcontrast = 1;
beta = maxfr/(2*alpha);
as = alpha*ones(50,1);
bs = beta*ones(50,1);
idxs = linspace(0,maxcontrast,50)';

neuron = inline('100./(1+exp(-(50*x-20)))','x');
%neuron = inline('max(0,100.*sin(6*x))','x');
%neuron = inline('100./(1+exp(-(2*x-10)))+100*normpdf(x,6,1)','x');

%neuron = inline('20*normpdf(x,.5,.1)','x');
neuron = inline('50*(x-.5).^2','x');

close all;
defaultsigma = .1;

sigma = defaultsigma;      
stims = [];
resps = [];
for i = 1:50
    vars = as.*bs.^2;
    contrast = idxs(vars == max(vars));
    contrast = contrast(unidrnd(length(contrast),1,1));
    resp = poissrnd(neuron(contrast));
   % resp = neuron(contrast)+normrnd(0,2);
    stims = [stims; contrast];
    resps = [resps; resp];

    subplot(2,2,1); hold on;

    subplot(2,2,3);
    plot(idxs,as.*bs);
    title('mean of lambda')
    subplot(2,2,4);
    plot(idxs,as.*bs.^2);
    title('variance of lambda')

    if (i > 4)
        [pp, p] = csaps(stims,resps);
        dpp = fnder(pp); % derivative
        subplot(2,2,2); cla; hold on;
        fnplt(pp,'g-')
        %plot(idxs,fnval(dpp,idxs),'r-')
        as = alpha*ones(50,1);
        bs = beta*ones(50,1);
        estimvals = fnval(pp,idxs);
        for j = 1:length(stims)
           % currentidx = find(stims(j)==idxs);
           % rhat = estimvals(currentidx);
            
            
           % stimresp = sortrows([stims resps]);
           % L = stims >=stims(i);
           % D = sum((resps-mean(resps)).^2)./mean(resps)
            
            
% Linear ramps - this didn't work terribly well
%             ledge = find(diff(abs(estimvals(currentidx:-1:1)-rhat)<10),1,'first');
%             if isempty(ledge) 
%                 ledge = 1;
%             else
%                 ledge = currentidx +1 - ledge;
%             end
%             redge = find(diff(abs(estimvals(currentidx:end)-rhat)<10),1,'first');
%             if isempty(redge) 
%                 redge = 50;
%             else
%                 redge = currentidx -1 + redge;
%             end
%             [ledge redge]
%             weighting = zeros(size(as));
%             leftidxs = ledge:find(stims(i)==idxs);
%             weighting(leftidxs) = linspace(0,1,length(leftidxs)); % linear ramp from 0
%             
%             rightidxs = find(stims(i)==idxs):redge;
%             weighting(rightidxs) = fliplr(linspace(0,1,length(rightidxs))); % linear ramp to 0

% Inversely related to spline predicted change in fr - also didn't work well
%            weighting = 1-abs(rhat-estimvals)./max(abs(rhat-estimvals));


            sigma = 5./(abs(fnval(dpp,stims(j))));  % numerator is spread parameter
            weighting = exp(-(idxs-stims(j)).^2/sigma.^2); % peaks at 1
            % updating the as and bs
            as = as+weighting*resps(j);
            bs = 1./(weighting+1./bs);
        end
    else
        weighting = exp(-(idxs-contrast).^2/sigma.^2);
        % updating the as and bs
        as = as+weighting*resp;
        bs = 1./(weighting+1./bs);
    end
    subplot(2,2,1); hold on
    plot(idxs,neuron(idxs));
    plot(stims(1:end-1),resps(1:end-1),'m*')
    subplot(2,2,4); hold on; plot(idxs, weighting);
    pause;
    cla;
end
subplot(2,2,4); cla;
[n,x] = hist(stims,50);
bar(x,n);
subplot(2,2,2); cla;
plot(stims,'.-')

%%
% Starting over from above but using a Poisson dispersion test to fix the
% limits of 'a' and 'b' propagation.
% COuld improve this by making the default spread adaptive(?)

% Setting up priors
alpha = 1;
maxfr = 200;
beta = maxfr/(2*alpha);
ndivs = 150; % number of unique contrasts
ntrials = 30; % number of simulated trials to perform
as = alpha*ones(ndivs,1);
bs = beta*ones(ndivs,1);
defaultspread = .2; % Maximum (1-sided) distance over which 'resp' can propagate
ALPHATHRESH = 0.01; % for finding the info propagation neighborhood

% Setting up simulation parameters and variables
maxcontrast = 1;
idxs = linspace(0,maxcontrast,ndivs)'; % the domain
neuron = inline('140./(1+exp(-(40*x-25)))','x');
neuron = inline('20*normpdf(x,.75,.1)','x');
neuron = inline('max(0,100.*sin(6*x))','x');
%neuron = inline('max(0,100.*sin(20*x))','x');


stims = [];
resps = [];
for i = 1:ntrials
    vars = as.*bs.^2;
    contrast = idxs(vars == max(vars));
    contrast = contrast(unidrnd(length(contrast),1,1)); % randomly selecting among the candidates
    resp = poissrnd(neuron(contrast));
    %resp = resp+normrnd(0,10); % extra noise
    %resp = poissrnd(resp); % extra noise
    stims = [stims; contrast];
    resps = [resps; resp];
    
    % a and b propagate as a linear ramp from the stimulus contrast to the
    % nearest stimulus contrast (in either direction) for which a Poisson
    % dispersion test rejects. The maximum distance this is allowed to be
    % away is "defaultspread". If the point is within defaultspread of one
    % of the ends, and the test does not reject, propagation is flat to
    % the edge of the contrast range.
    stimresp = sortrows([stims resps]);
    lowerend = nan*ones(length(stims),1);
    upperend = nan*ones(length(stims),1);
    for j = 1:size(stimresp,1) % Looping over all stimuli ever shown
        lowerend(j) = stimresp(j,1)-defaultspread; % Starting with default
        % finding all the contrasts within the neighborhood that we
        % care about. Contrasts must be greater than the default lower
        % bound and less than the default upper bound
        
        % First going to the left (lower values of contrast)
        Lneighborhood = stimresp(:,1) > lowerend(j) & stimresp(:,1) < stimresp(j,1);
        if (any(Lneighborhood))
            p = nan*ones(size(Lneighborhood,1),1);
            for k = find(Lneighborhood)' % haw many lower contrasts to look at
                D = sum((stimresp(k:j,2)-mean(stimresp(k:j,2))).^2)./mean(stimresp(k:j,2));
                p(k) = 1-chi2cdf(D,j-k);
            end
            k = find(p < ALPHATHRESH,1,'last');
            if (any(k)) 
                disp('Got here1 (lower end)');
                stimresp([j k],:)
                lowerend(j) = stimresp(k,1);
            end
        end
        
        % ---------
        % Now working on the upper end
        % ---------
         upperend(j) = stimresp(j,1)+defaultspread; % temporary
         Lneighborhood = stimresp(:,1) < upperend(j) & stimresp(:,1) > stimresp(j,1);
        if (any(Lneighborhood))
            p = nan*ones(size(Lneighborhood,1),1);
            for k = find(Lneighborhood)' % haw many lower contrasts to look at
                D = sum((stimresp(j:k,2)-mean(stimresp(j:k,2))).^2)./mean(stimresp(j:k,2));
                p(k) = 1-chi2cdf(D,k-j);
            end
            k = find(p < ALPHATHRESH,1,'first');
            if (any(k)) 
                disp('Got here 2 (upper end)');
                stimresp([j k],:)
                upperend(j) = stimresp(k,1);
            end
        end

    end
    
    % Now going to the right (higher values of contrast)
    % Now updating the the as and bs
    % starting from scratch each time
    as = alpha*ones(ndivs,1);
    bs = beta*ones(ndivs,1);
    subplot(2,2,2); cla;
    for j = 1:size(stimresp,1)
        weighting = zeros(length(idxs),1);
        weighting(idxs==stimresp(j,1)) = 1; %learn the most about where we are (avoid round off error)
        % First working on the lower limb of the weighting function
        m = 1./(stimresp(j,1)-lowerend(j));
        b = -m*lowerend(j);
        tmpweight = m*idxs+b;
        weighting(tmpweight > 0 & tmpweight <=1) = tmpweight(tmpweight > 0 & tmpweight <=1);
        
        % Now the upper limb (basically a repeat of the above code)
        m = 1./(stimresp(j,1)-upperend(j));
        b = -m*upperend(j);
        tmpweight = m*idxs+b;
        if upperend(j) > 1 % avoid edge effect
            weighting(tmpweight > 0 & tmpweight <=1) = 1;
        else 
            weighting(tmpweight > 0 & tmpweight <=1) = tmpweight(tmpweight > 0 & tmpweight <=1);
        end
        % updating the as and bs
        as = as+weighting*stimresp(j,2);
        bs = 1./(weighting+1./bs);
        
        subplot(2,2,2); hold on;
        plot(idxs,weighting);
        title('weighting function');
    end
    
    % Doing all the plotting
    subplot(2,2,1); hold on
    plot(idxs,neuron(idxs)); % The real CR function
    plot(stims,resps,'m*')
    subplot(2,2,3);
    plot(idxs,as.*bs);
    title('mean of lambda');

    subplot(2,2,4);
    plot(idxs,as.*bs.^2);
    title('variance of lambda');
    pause

end

subplot(2,2,4); cla;
[n,x] = hist(stims,ndivs);
bar(x,n);
subplot(2,2,2); cla;
plot(stims,'.-')

%%
% Trying a Poisson dispersion test
lambda = 10;
n = 10; % number of draws
niter = 5000;

Ds = [];ps= [];
for i = 1:niter
    x = poissrnd(lambda,1,n);
    %x(end) = x(end)+10;
    
    Ds(i) = sum((x-mean(x)).^2)./mean(x);
    ps(i) = 1-chi2cdf(Ds(i),n-1);
end
figure; hist(ps,40);
%%
% Explaining recursive Bayesian estimation of a Poisson parameter
lambdas = linspace(0,50,1000);
x = [0:100];
fr = [];
for i = 1:length(lambdas)
    fr(i,:) = poisspdf(x,lambdas(i));
end

figure; axes; hold on;
surf(x,lambdas,fr-.1);
shading interp
xlabel('f(x|lambda)');
ylabel('lambda');
axis ij
plot([x(1) x(end)],[20 20],'k-');
plot([10 10],[lambdas(1) lambdas(end)],'k-');
a = get(gca,'YTick');

% looking at one of the Poisson distns
figure; bar(x,fr(300,:));
xlabel('f(x|lambda)');
title(['lambda = ',num2str(lambdas(300))])

% Looking at one of the distns of lambda
figure; plot(lambdas,fr(:,10));
xlabel('f(lambda|x)');
title(['x = ',num2str(x(10))])

%%
% Messing round with models
[x,y] = meshgrid(linspace(-10,10,100), linspace(-10,10,100));

params_x = [100 50 2 2 2 2 2];
f_x = ComputeNakaRushtonJPW(params_x,x);

params_y = [100 50 10 10 2 2 2];
f_y = ComputeNakaRushtonJPW(params_y,y);

params_xy = [100 20 2 3];
a = 1;
b = 2;
f_xy = ComputeNakaRushtonJPW(params_xy,a*x.^2+b*y.^2);

subplot(2,2,1);
surf(x,y,reshape(f_x,size(x))); axis square; set(gca,'View',[0 90]);
subplot(2,2,2);
surf(x,y,reshape(f_y,size(x))); axis square; set(gca,'View',[0 90]);
%subplot(2,2,3);
%surf(x,y,f_x.^2+f_y.^2); axis square; set(gca,'View',[0 90]);
subplot(2,2,4);
surf(x,y,reshape(f_xy,size(x))); axis square; set(gca,'View',[0 90]);



%%
% Looking at V1 retinotopy.
% r in mm, theta in deg, RF_x in deg, RF_y in deg
a = ...
[4	270	-1.1	-2.5;...
6	270	-1.7	-2.4;...
7	270	-2.1	-2.4;...
8	270	-2.1	-2.3;...
3	270	-0.4	-2.2;...
2	270	-0.3	-2.1;...
1	270	0	-2.3;...
0	270    nan   nan;...
5	0	nan	  nan;...
5	315	-0.6	-1.8;...
4	315	-0.4	-1.7;...
2	315	0	-1.6;...
8	330	-0.3	-1.1;...
8	345 nan     nan;...
8	240	-0.8	-3.3;...
5	330	-0.3	-1.4;...
5	300	-0.8	-1.7;...
5	285	-1	-2;...
5	240	-0.8	-3;...
5	225	-0.3	-3.4;...
5	210	0.1	-3.3;...
8	225	-1.3	-4.1;...
8	337.5	-0.1	-1.2];

[ch_x,ch_y] = pol2cart(a(:,2)*pi/180, a(:,1));
rfx = a(:,3); rfy = a(:,4);


% % predicting rf location as a function of chamber location
% b_x = regress(rfx, [ch_x,ch_y, ones(length(ch_x),1)]);
% b_y = regress(rfy, [ch_x,ch_y, ones(length(ch_y),1)]);
% % Looking at the residuals
% pred_x = [ch_x,ch_y, ones(length(ch_x),1)]*b_x;
% pred_y = [ch_x,ch_y, ones(length(ch_y),1)]*b_y;
% figure; axes; hold on;
% plot(rfx,pred_x,'r.');
% plot(rfy,pred_y,'b.');
% 
% % predicting chamber location as a function of RF location
% b_x = regress(ch_x, [rfx,rfy, ones(length(rfx),1)]);
% b_y = regress(ch_y, [rfx,rfy, ones(length(rfy),1)]);
% % Looking at the residuals
% pred_x = [rfx,rfy, ones(length(ch_x),1)]*b_x;
% pred_y = [rfx,rfy, ones(length(ch_y),1)]*b_y;
% figure; axes; hold on;
% plot(ch_x,pred_x,'r.');
% plot(ch_y,pred_y,'b.');

NORFACT = 10; % Normalization factor for pp;
% using color to indicate the mapping
[x,y] = meshgrid (linspace(-4,0,31),linspace(-4,0,31)); % On the screen in deg
b_x = regress(rfx, [ch_x,ch_y, ones(length(ch_x),1)]); % b_x goes from chamber to RF location
b_y = regress(rfy, [ch_x,ch_y, ones(length(ch_y),1)]);
figure; subplot(2,1,1); hold on;
image(cat(3,(x+2)/NORFACT+.5,(y+2)/NORFACT+.5,.5*ones(size(x))));
axis xy; axis image;
set(gca,'Xtick',floor(linspace(1,size(x,1),5)),'XtickLabel',x(1,floor(linspace(1,size(x,1),5))))
set(gca,'Ytick',floor(linspace(1,size(x,1),5)),'YtickLabel',y(floor(linspace(1,size(x,1),5)),1))

[x,y] = meshgrid (linspace(-10,10,31),linspace(-10,10,31)); % In the chamber (in mm)
pred_x = reshape([x(:),y(:), ones(length(x(:)),1)]*b_x, size(x)); % On the screen (in deg)
pred_y = reshape([x(:),y(:), ones(length(y(:)),1)]*b_y, size(y)); % On the screen (in deg)
mask = sqrt(x.^2+y.^2) > 8 | pred_x > 1 | pred_y > 1;
pred_x = pred_x.*~mask;
pred_y = pred_y.*~mask;
subplot(2,1,2); hold on;
image(cat(3,(pred_x+2)/NORFACT+.5,(pred_y+2)/NORFACT+.5,(2/NORFACT)+.5*ones(size(x))));
axis xy; axis image;
set(gca,'Xtick',floor(linspace(1,size(x,1),5)),'XtickLabel',x(1,floor(linspace(1,size(x,1),5))))
set(gca,'Ytick',floor(linspace(1,size(x,1),5)),'YtickLabel',y(floor(linspace(1,size(x,1),5)),1))

% predicting RF location for a given location in the chamber
x = 0;
y = -7;
pred_x = [x y 1]*b_x; % On the screen (in deg)
pred_y = [x y 1]*b_y; % On the screen (in deg)
[pred_x pred_y]



