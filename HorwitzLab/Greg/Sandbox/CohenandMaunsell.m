% Section 1
% Simulations that test the analysis of Cohen and Maunsell (2009)
% Does the d' value vary with correlation?  Seems like it does...
CROSSVALIDATE = 0;
nx = 100; % N trials in cond X
ny = 100; % N trials in cond y
r = .5;  % inter neuronal correlation
d = 1000;  % Number of neurons
niter = 100;
data = [];
for i = 1:niter
    x = normrnd(0,1,d,nx);
    y = normrnd(0,1,d,ny);
    
    if (CROSSVALIDATE)
        v = mean(normrnd(0,1,d,nx),2)-mean(normrnd(0,1,d,ny),2);
    else
        v = mean(x,2)-mean(y,2);
    end
    px = v'*x;
    py = v'*y;
    bins = linspace(min([px';py']),max([px';py']),100);
    %figure; subplot(2,1,1); hold on;
    tmp = hist(px,bins);
    bar(bins, tmp,'r');
    tmp = hist(py,bins);
    bar(bins, tmp,'k');
    dprime = abs(mean(px)-mean(py))/mean([std(px) std(py)]);
    title(num2str(dprime));
    
    % adding intraneuronal correlations
    S2 = repmat(r,d,d)+diag(ones(1,d)-r);
    S = sqrtm(S2);
    x_r = S*x;
    y_r = S*y;
    if (CROSSVALIDATE)
        v_r = mean(S*normrnd(0,1,d,nx),2)-mean(S*normrnd(0,1,d,ny),2);
    else
        v_r = mean(x_r,2)-mean(y_r,2);
    end
    px = v_r'*x_r;
    py = v_r'*y_r;
    bins = linspace(min([px';py']),max([px';py']),100);
    subplot(2,1,2); hold on;
    tmp = hist(px,bins);
    bar(bins, tmp,'r');
    tmp = hist(py,bins);
    bar(bins, tmp,'k');
    dprime_r = abs(mean(px)-mean(py))/mean([std(px) std(py)]);
    title(num2str(dprime_r));
    
    % sanity check
    %figure; axes;
    %imagesc(cov(x_r'));
    
    data = [data; dprime dprime_r];
end
figure
hist(data(:,1)-data(:,2));
%%
% Section 2
% Playing around with the standard estimates of sigma and tau in a random
% effects ANOVA.  Jim DiCarlo uses this in his goodness of fit measure.

mu = 20;
sigma = .5;
tau = 2;
J = 8;
I = 10;
niter = 20000;

    
U1 = normrnd(0,tau,J,niter);  % normally distributed random effects
%U = exprnd(tau,J,niter)-tau; % exponentially distributed random effects
U2 = reshape(permute(repmat(U1,[1, 1, I]),[3 1 2]),[I*J, niter]);
W = normrnd(0,sigma,I*J,niter); % normally distributed errors
%W = exprnd(sigma,I*J,niter)-tau; % exponentially distributed errors

Y = U2+W+mu;
data = []; newdata = [];
for i = 1:niter
    tmp = reshape(Y(:,i),[I,J]);
    mns = mean(tmp);
    modelfit = U1(:,i)'+normrnd(0,.5,1,J); % A goofy model 
    
    
    SSW = sum(sum((tmp-repmat(mns,[I,1])).^2));
    SSB = I*sum((mns-mean(mns)).^2);
    data(i,:) = [SSW/(J*(I-1)) SSB/(I*(J-1))-SSW/(J*I*(I-1))];
    
    % Now using a preset model (instead of within group means) to compute
    % SSW and SSB.  This is supposed to emulate the calculation of
    % sigma_res in the DiCarlo paper.
    
    tmp = tmp - repmat(modelfit,[I,1]);
    mns = mean(tmp);
    SSW = sum(sum((tmp-repmat(mns,[I,1])).^2));
    SSB = I*sum((mns-mean(mns)).^2);
    newdata(i,:) = [SSW/(J*(I-1)) SSB/(I*(J-1))-SSW/(J*I*(I-1))];
 
end
sqrt(mean(data))
%plot(data(:,1),data(:,2),'k.');
[rho,p]= corr(data)

%figure
%subplot(3,1,1);
%hist(data(:,1),30);
%subplot(3,1,2);
%hist(data(:,2),30);
%subplot(3,1,3);
%hist(data(:,2)./data(:,1),30);
%[h,p] = ttest(data(:,2)./data(:,1)-sigma/tau)

% Pretty cool - these equations *do* provide unbiassed estimates of the two
% components of the variance.  Interestingly, the estimates are not
% independendent, which is I guess why the F test is based on SSW and SSB
% instead of SSW and SSB-SSW.
% Also, no assumption of normality is required (for either the block
% effects or the errors).

figure;
subplot(2,1,1);
plot(newdata(:,2),data(:,2),'k.')
subplot(2,1,2);
GOF = 1-(newdata(:,2)./data(:,2));
hist(GOF,20);
%%
% Section 3
% Back to Cohen and Maunsell (2009).  Sampling from two multivariate
% distributions (attend and non-attend).  Calling some of the samples
% "corrects" and some of the samples "incorrect".  Projecting all the
% points onto the axis connecting the means of the corrects.  Are the
% corrects father away from each other along this axis than the incorrects?
d = 100;
n = 80;
pcor = .5;
mu1 = 1;
mu2 = -1;

X = normrnd(mu1,1,n,d);
Y = normrnd(mu2,1,n,d);

Lcorx = logical(binornd(1,pcor,n,1));
Lcory = logical(binornd(1,pcor,n,1));

xbar = mean(X(Lcory,:));
ybar = mean(Y(Lcorx,:));
vect = xbar-ybar;
%vect = vect-mean(vect);
vect = vect./norm(vect);
projX = (X-repmat(mean([xbar; ybar]),n,1))*vect';
projY = (Y-repmat(mean([xbar; ybar]),n,1))*vect';

bins = linspace(-20,20,200);
[xcor,tmp] = hist(projX(Lcory),bins);
[xinc,tmp] = hist(projX(~Lcorx),bins);
[ycor,tmp] = hist(projY(Lcory),bins);
[yinc,tmp] = hist(projY(~Lcory),bins);

figure
subplot(2,1,1); hold on;
bar(bins,xcor)
bar(bins,-xinc,'r')

subplot(2,1,2); hold on;
bar(bins,ycor)
bar(bins,-yinc,'r')

%%
% Marlene's rebuttal
% Simulations starts here\
d = 80;   % number of neurons\
n = 100;     % number of trials\
%d = 1000
%n = 10
pcor = .5;    % Probability correct\
mu1 = 1;    % mean for dist 1\
mu2 = 1.5;    % mean for dist 2\

X = normrnd(mu1,1,n,d);
Y = normrnd(mu2,1,n,d);

Lcorx = logical(binornd(1,pcor,n,1));
Lcory = logical(binornd(1,pcor,n,1));

xbar1 = mean(X(Lcorx,:));  % was Lcory in original code\
xbar2 = mean(Y(Lcory,:));
vect = xbar2-xbar1;

% vect = vect-mean(vect);\
% vect = vect./norm(vect);\
% projX = X*vect';\
% projY = Y*vect';\

% compute projection by normalizing the dot product by the square of the\
% norm of vect, which makes proj(xbar1)=0 and proj(xbar2)=1.  Multiply by 2\
% and subtract 1 to make proj(xbar1)=-1 and proj(xbar2)=1.\
projX=-2*dot((repmat(xbar1,n,1)-X)',repmat(vect,n,1)')./(norm(vect).^2)-1;
projY=-2*dot((repmat(xbar1,n,1)-Y)',repmat(vect,n,1)')./(norm(vect).^2)-1;

bins = linspace(-4,4,20);
[xcor,tmp] = hist(projX(Lcorx),bins); % was Lcory in original code\
[xinc,tmp] = hist(projX(~Lcorx),bins);
[ycor,tmp] = hist(projY(Lcory),bins);
[yinc,tmp] = hist(projY(~Lcory),bins);

figure
subplot(2,1,1); hold on;
bar(bins,xcor)
bar(bins,-xinc,'r')

subplot(2,1,2); hold on;
bar(bins,ycor)
bar(bins,-yinc,'r')

[h,p1] = ttest2(projX(Lcorx), projX(~Lcorx))
[h,p2] = ttest2(projY(Lcory), projY(~Lcory))

%%
% Comparing my way of projecting to Marlene's way

d = 80;   % number of neurons
n = 100;     % number of trials
pcor = .5;    % Probability correct
mu1 = 1;    % mean for dist 1
mu2 = 1.5;    % mean for dist 2
 
X = normrnd(mu1,1,n,d);
Y = normrnd(mu2,1,n,d);

Lcorx = logical(binornd(1,pcor,n,1));
Lcory = logical(binornd(1,pcor,n,1));

xbar1 = mean(X(Lcorx,:));
xbar2 = mean(Y(Lcory,:));
vect = xbar2-xbar1;

% Marlene's way
projX1=-2*dot((repmat(xbar1,n,1)-X)',repmat(vect,n,1)')./(norm(vect).^2)-1;
projY1=-2*dot((repmat(xbar1,n,1)-Y)',repmat(vect,n,1)')./(norm(vect).^2)-1;

% My way
projX2 = (repmat(xbar1,n,1)-X)*(vect./norm(vect))';
projY2 = (repmat(xbar1,n,1)-Y)*(vect./norm(vect))';

% Linearly transforming projections so that the means are -1 and 1
a = [mean(projX2) mean(projY2); 1 1];
b = [-1 1]*inv(a);
projX2 = b(1)*projX2+b(2);
projY2 = b(1)*projY2+b(2);

figure
subplot(2,1,1); hold on;
plot(projX1,projX2,'.');
plot([-.5 -1.5],[-.5 -1.5],'k-');
axis equal;
subplot(2,1,2); hold on;
plot(projY1,projY2,'.');
plot([.5 1.5],[.5 1.5],'k-');
axis equal

bins = linspace(-4,4,100);
[xcor,tmp] = hist(projX1(Lcorx),bins);
[xinc,tmp] = hist(projX1(~Lcorx),bins);
[ycor,tmp] = hist(projY1(Lcory),bins);
[yinc,tmp] = hist(projY1(~Lcory),bins);

figure
subplot(2,1,1); hold on;
bar(bins,xcor)
bar(bins,-xinc,'r')
[h,p1] = ttest2(projX1(Lcorx), projX1(~Lcorx))
title(['p = ',num2str(p1)]);

subplot(2,1,2); hold on;
bar(bins,ycor)
bar(bins,-yinc,'r')
[h,p2] = ttest2(projY1(Lcory), projY1(~Lcory))
title(['p = ',num2str(p2)]);
