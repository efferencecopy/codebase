% Figures for Stats Seminar Lectures
% 
% Section 1:
% Distribution of axons diameters from and conduction velocities
%
% Section 2:
% CDFs
%
% Section 3:
% 3-D rendered 2-D probability distributions.
%
% Section 4:
% Convergence of a proportion to the underlying probability.
% Demonstration of the law of large numbers.
%
% Section 5:
% Normal and chisquared(1) PDFs
%
% Section 6:
% An assortment of 2-D PDFs illustrating the difference between
% independence and uncorrelatedness.
%
% Section 7:
% A few 1-D distributions showing how variance changes with scale.
%
% Section 8:
% How the variance of a bernoulli random variable changes with 'p'.
%
% Section 9:
% A few example binomial PMFs
%
% Section 10:
% Binomial distribution n = 5, p = 0.5.  The null distn for a sign test
% example.  Also n = 5, p = 0.7 (an alternative hypothesis).
%
% Section 11:
% A few binomials coverging on the Poisson distribution.
%
% Section 12:
% A few binomials coverging on the Gaussian distribution.
%
% Section 13:
% N-fold convolution of two PDFs, converging on a Gaussian
%
% Section 14:
% Joint Gaussian and marginals
%
% Section 15:
% Same marginals as section 14, but not a jointly Gaussian (not sure how to
% do this...)
%
% Section 16:
% A few example chi-squared distributions.
%
% Section 17:
% A few example F distributions.
%
% Section 18:
% Comparing the distributions of sample means and variances from a Poisson
% distn (relative efficiency for lambda).
%
% Section 19:
% Comparing the relative efficiency of the sample median and xbar for estimating
% mu from a normal dist (or a mixture of two normals).
%
% Section 20:
% Likelihood function for binomial with X = 12 n = 20.
% Also binomial PMF with n=20, p=0.1 and n=20, p = 0.3
%
% Section 21:
% A pair of Gaussians for the large sample sign test.
%
% Section 22:
% Comparison of tail probabilities from binomial and normal (for some n)
%
% Section 23:
% Progression of overlapping binomial distributions as N increases (showing
% how the power of a test increases with n)
%
% Section 24:
% Stimulus movie (B/W or color)
%
% Section 25:
% Skewed and kurtotic distributions.
%
% Section 26:
% A few made-up stimulus frames.
%
% Section 27:
% An example STA, Zscore, null distribution, map of p-values.
%
% Section 28:
% An example STV, Zscore, null distribution, map of p-values.
% Jusrt like section 28 but for the variance.
%
% Section 29:
% A few spike-triggered stimulus frames.
%
% Section 30:
% Numerator and denomenator distributions for the one-sample t-test example
% (normal and sqrt(chi-squared)
%
% Section 31:
% A few t-distributions
%
% Section 32:
% Computations for the example paired t-test
%
% Section 33:
% Simulations comparing the power and coverage probability of the paired
% t-test and sign test
%
% Section 34:
% Simulations comparing the power and coverage probability of the paired
% and unpaired t-tests.
%
% Section 35:
% Simulations showing that the t-test is robust to departures from normality
% 
% Section 36:
% Looking at the effect of negative correlations on the paired and unpaired
% t-test.
%
% Section 37:
% Simulations shoing the effects of non-independence on the t-test
%
% Section 38:
% An example regression problem (nonlinear with weird errors)
% Just to motivate the idea of regression in general.
%
% Section 39:
% A worked out linear regression problem with t-test/CIs for slope.
%
% Section 40:
% A regression with partitioned sums of squares
%
% Section 41:
% Examples of regressions with various regSS and errorSS
%
% Section 42:
% Residuals from a univariate regression
%
% Section 43:
% Examples of problem residuals
%
% Section 44:
% Regression on transformed data
%
% Section 45:
% Polynomial regressions
%
% Section 46:
% Nearly singular regression (multicollinearity)
%
% Section 47
% uniform iid process
%
% Section 48
% uniform autoregressive(1) process
%
% Section 49
% Curse of dimensionality (average distance as a function of dimensionality)
%
% Section 50
% An example of linear discriminant analysis
%
% Section 51
% Figures for the probelm set: eg. A continuous 2-D non-separable PDF
% 
% Section 52
% A simulation showing that using s instead of sigma results in a
% z-statistic (really a t-statistic) that doesn't have a standard normal
% distn.
%
% Section 53
% The null distribution of the Wilcoxon signed rank statistic
%
% Section 54
% Simulated end points of dart tosses (for permutation test)
%
% Section 55
% Bootstrapping intro
%
% Section 56
% Bootstrap estimates of the standard error of the median of a normal
% distribution
%
% Section 57
% Bootstrap t-tests with various distributions (see section 35)
%
% Section 60
% Distributions of dfferences of rvs.
%%
% Section 1 
% Distribution of axon diameters
% Basser, Peter 3 March 2009, SPIE Newsroom. DOI 10.1117/2.1200902.1515
x = linspace(0,20,1000);
y  = gampdf(x,2,2.2);  % 2, 2.2 looks approximately right
plot(x,y);
ylabel('Probability')
xlabel('Axon Diam (\mum)');

r = 35.4; % Ohms
c = 10^-6;
K = 10470;
cv = sqrt(K*(x/2)*r*c)

figure;
plot(cv,y);
ylabel('Probability')
xlabel('Conduction vel. (m/s)');

a = sqrt(10470*35.4*10^-6)  % equation is vel = a*sqrt(diam) (I think)

%%
% Section 2
% CDFs
x = linspace(0,20,1000);
y  = gampdf(x,2,2.2);  % 2, 2.2 looks approximately right
cdf = cumsum(y)./sum(y);

figure;
plot(x,cdf);



%%
% Section 3
% 3-D rendered 2-D probability distributions

% FIrst a discrete distribution.  Uniform prior on p, binomial in the other
% dimension.
p = linspace(0,.5,10);
data = [];
for i = 1:length(p)
    x = binopdf(round(10*p),10*ones(1,length(p)),p(i));
    data = [data; x];
end
data = data./sum(data(:));

figure;
bar3(data);
set(gca,'color',[.5 .5 .5]);
set(gca,'XTick',[],'YTick',[]);
shading faceted;
view(-55,20);
% Continuous  (2-D Gaussian)

[x,y] = meshgrid(linspace(-4,4,100),linspace(-4,4,100));
rho = .4;
normfact = 1/(2*pi*sqrt(1-rho^2));
z = normfact*exp(-1/(2*(1-rho^2)).*(x.^2+y.^2-(2.*rho.*x.*y)))

figure;
h = surf(x,y,z);
set(h,'EdgeAlpha',0);
colormap(jet);
set(gca,'color',[.5 .5 .5]);
shading interp;
set(gca,'XTick',[],'YTick',[]);
view(68,46)

%%
% Section 4
% Demonstration of the law of large numbers with simulated coin flips
p = .5;
n = 6000;
x = binornd(1,p,n,1);
plot([1:n],cumsum(x)./[1:n]','k-','linewidth',2);
set(gca,'Ylim',[.2 .8]);

%%
% Section 5
% Normal and chisquared(1) PDFs
x1 = linspace(-3, 3, 400);
x2 = linspace(.1, 9, 400);
y = normpdf(x1,0,1);
z = chi2pdf(x2,1);
y = y./sum(y);
z = z./sum(z);

figure;
subplot(2,1,1);
plot(x1,y)
subplot(2,1,2);
plot(x2,z)

%%
% Section 6
% PDFs of independent and uncorrelated random variables

% Indep, uncorrelated Gaussian (cov = 0)
[x,y] = meshgrid(linspace(-4,4,100),linspace(-4,4,100));
z1 = normpdf(x,0,1.5).*normpdf(y,0,.5);
z1 = z1./sum(z1(:));

% Correlated, nonindep Gaussian (cov > 0)
rho = .6;
normfact = 1/(2*pi*sqrt(1-rho^2));
z2 = normfact*exp(-1/(2*(1-rho^2)).*(x.^2+y.^2-(2.*rho.*x.*y)));
z2 = z2./sum(z2(:));

% Indep, uniform square (cov = 0)
z3 = (abs(x)<2&abs(y)<2);
z3 = z3./sum(z3(:));

% Nonindep, uniform rectangle (cov > 0)
z4 = (abs(x-y)<1&abs(x+y)<5);
z4 = z4./sum(z4(:));

% Nonindep, uniform square (cov = 0)
z5 = (abs(x-y)<3&abs(x+y)<3);
z5 = z5./sum(z5(:));

% Nonindep Gaussian with hole in the middle (cov = 0)
z6 = (normpdf(x,0,1).*normpdf(y,0,1))-.55*(normpdf(x,0,.75).*normpdf(y,0,.75));
z6 = z6./sum(z6(:));

figure;
subplot(2,3,1);
imagesc(z1);
axis square;
set(gca,'XTick',[],'YTick',[]);

subplot(2,3,2);
imagesc(z6);
axis square;
set(gca,'XTick',[],'YTick',[]);

subplot(2,3,3);
imagesc(z2);
axis square;
set(gca,'XTick',[],'YTick',[]);

subplot(2,3,4);
imagesc(z3);
axis square;
set(gca,'XTick',[],'YTick',[]);

subplot(2,3,5);
imagesc(z5);
axis square;
set(gca,'XTick',[],'YTick',[]);

subplot(2,3,6);
imagesc(z4);
axis square;
set(gca,'XTick',[],'YTick',[]);

colormap(gray);

% Checking covariances

sum(sum((x.*y).*z1))
sum(sum((x.*y).*z2))
sum(sum((x.*y).*z3))
sum(sum((x.*y).*z4))
sum(sum((x.*y).*z5))
sum(sum((x.*y).*z6))

%%
% Section 7
% A bunch of PDFs to give a feeling for what "variance" means

figure;
x = linspace(-10,10,500);

% Normals
subplot(3,3,1)
y = normpdf(x,0,1);
plot(x,y);
set(gca,'Ylim',[0 .4]);

subplot(3,3,2)
y = normpdf(x,0,2);
plot(x,y);
set(gca,'Ylim',[0 .4]);

subplot(3,3,3)
y = normpdf(x,0,3);
plot(x,y);
set(gca,'Ylim',[0 .4]);

% exponentials
subplot(3,3,4)
y = exppdf(x,1);
plot(x,y);
set(gca,'Ylim',[0 1]);
[m,v] = expstat(1)

subplot(3,3,5)
y = exppdf(x,1);
plot(-x,y);
set(gca,'Ylim',[0 1]);
[m,v] = expstat(1)

% uniform
subplot(3,3,6)
y = unifpdf(x,-1.732,1.732);
plot(x,y);
set(gca,'Ylim',[0 .4]);
[m,v] = unifstat(-1.732, 1.732)

% Now a few discrete distributions
x = linspace(-1,1,3);
p = ones(length(x),1)./length(x);
var = x.^2*p

x = linspace(-1.732,1.732,5);
p = ones(length(x),1)./length(x);
var = x.^2*p

x = linspace(-1.732,1.732,10);
p = ones(length(x),1)./length(x);
var = x.^2*p

%%
% Section 8
% How the variance of a bernoulli random variable changes with 'p'.
p = linspace(0,1,100);
plot(p,p.*(1-p));

%%
% Section 9
% A few example binomials 
figure;
x = [0:5];
y = binopdf(x,5,0.5);
subplot(1,3,1);
h = bar(x,y,.2);
set(h,'FaceColor','yellow','EdgeColor','yellow');
set(gca,'XLim',[x(1)-.5 x(end)+.5],'YLim',[0 0.35]);

x = [0:20];
y = binopdf(x,20,0.5);
subplot(1,3,2);
h = bar(x,y,.2);
set(h,'FaceColor','yellow','EdgeColor','yellow');
set(gca,'XLim',[x(1)-.5 x(end)+.5],'YLim',[0 0.35]);

x = [0:20];
y = binopdf(x,20,0.9);
subplot(1,3,3);
h = bar(x,y,.2);
set(h,'FaceColor','yellow','EdgeColor','yellow');
set(gca,'XLim',[x(1)-.5 x(end)+.5],'YLim',[0 0.35]);
%%
% Section 10:
% Binomial distribution n = 5, p = 0.5.  The null distn for a sign test
% example.
n = 5;
p = .6;
figure; axes; hold on;

x = [0:n];
y = binopdf(x,n,p);
h = bar(x,y,.2);
set(h,'FaceColor','yellow','EdgeColor','yellow');

set(gca,'XLim',[x(1)-.5 x(end)+.5],'YLim',[0 0.40]);
p = sum(y(L))

%%
% Section 11
% Binomial distributions converging to a Poisson
figure;
subplot(2,2,1);
n = 10;
x = [0:n];
y = binopdf(x,n,0.5);
h = bar(x,y,.2);
set(h,'FaceColor','yellow','EdgeColor','yellow');
set(gca,'XLim',[x(1)-.5 15+.5],'XTick',[0 5 10 15],'YLim',[0 .25]);

subplot(2,2,2);
n = 20;
x = [0:n];
y = binopdf(x,n,0.25);
h = bar(x,y,.2);
set(h,'FaceColor','yellow','EdgeColor','yellow');
set(gca,'XLim',[x(1)-.5 15+.5],'XTick',[0 5 10 15],'YLim',[0 .25]);

subplot(2,2,3);
n = 100;
x = [0:n];
y = binopdf(x,n,0.05);
h = bar(x,y,.2);
set(h,'FaceColor','yellow','EdgeColor','yellow');
set(gca,'XLim',[x(1)-.5 15+.5],'XTick',[0 5 10 15],'YLim',[0 .25]);


subplot(2,2,4);
n = 1000;
x = [0:n];
y = poisspdf(x,5);
h = bar(x,y,.2);
set(h,'FaceColor','yellow','EdgeColor','yellow');
set(gca,'XLim',[x(1)-.5 15+.5],'XTick',[0 5 10 15],'YLim',[0 .25]);

%%
% Section 12:
% A few binomials coverging on the Gaussian distribution.
figure;
subplot(2,2,1);
n = 10;
x = [0:n];
y = binopdf(x,n,0.5);
h = bar(x,y,.2);
set(h,'FaceColor','yellow','EdgeColor','yellow');
sd = sqrt(n*.25);
set(gca,'XLim',[(0.5*n)-(4*sd) (0.5*n)+(4*sd)],'XTick',[0 5 10 15],'YLim',[0 .25]);

subplot(2,2,2);
n = 20;
x = [0:n];
y = binopdf(x,n,0.5);
h = bar(x,y,.2);
sd = sqrt(n*.25);
set(h,'FaceColor','yellow','EdgeColor','yellow');
set(gca,'XLim',[(0.5*n)-(4*sd) (0.5*n)+(4*sd)],'XTick',[0 5 10 15],'YLim',[0 .2]);

subplot(2,2,3);
n = 100;
x = [0:n];
y = binopdf(x,n,0.5);
h = bar(x,y,.2);
sd = sqrt(n*.25);
set(h,'FaceColor','yellow','EdgeColor','yellow');
set(gca,'XLim',[(0.5*n)-(4*sd) (0.5*n)+(4*sd)],'YLim',[0 .1]);

subplot(2,2,4);
n = 100;
x = [0:.1:100]
sd = sqrt(n*.25);
y = normpdf(x,n*0.5,sd);
y = y./sum(y)*10;
h = plot(x,y,'y-','LineWidth',2);
set(gca,'XLim',[(0.5*n)-(4*sd) (0.5*n)+(4*sd)],'YLim',[0 .1]);

%%
% Section 13:
% N-fold convolution of two PDFs, converging on a Gaussian
% 5-fold and then 30-fold

% Uniform
pdf = [ones(17,1);zeros(10,1)];
pdf = pdf./sum(pdf);
data = pdf;
figure;
subplot(3,2,1);
plot(pdf);
for i = 1:5
    data = conv(data,pdf);
end
subplot(3,2,3)
plot(data);
set(gca,'XLim',[0 100]);
for i = 1:25
    data = conv(data,pdf);
end
subplot(3,2,5)
plot(data)
set(gca,'XLim',[100 400]);

% Exponentials
pdf = [zeros(4,1);exppdf([0:25],4.6)'];
pdf = pdf./sum(pdf);
data = pdf;
subplot(3,2,2);
plot(pdf);
for i = 1:5
    data = conv(data,pdf);
end
subplot(3,2,4)
plot(data)
for i = 1:25
    data = conv(data,pdf);
end
set(gca,'XLim',[0 100]);
subplot(3,2,6)
plot(data)
set(gca,'XLim',[100 400]);
set(gca,'YLim',[0 0.018]);

%%
% Section 14:
% Jointly Gaussian PDF and marginals

[x,y] = meshgrid(linspace(-4,4,100),linspace(-4,4,100));
sx = 1;
sy = 1;
rho = 0;
normfact = 1/(2*pi*sx*sy*sqrt(1-rho^2));
z = normfact*exp(-1/(2*(1-rho^2)).*(x.^2/sx^2+y.^2/sy^2-(2.*rho.*x.*y)/(sx*sy)))
subplot(2,2,1);
imagesc(z); set(gca,'XTick',[],'YTick',[]);
colormap(gray)
axis square;
pdf1 = sum(z);
pdf2 = sum(z');
subplot(2,2,2);
plot(pdf1);
subplot(2,2,3);
plot(pdf2);

%%
% Section 15:
% Marginals are Gaussian, but joint density is not

[x,y] = meshgrid(linspace(-4,4,100),linspace(-4,4,100));
sx = 1;
sy = .5;
rho = 0;
normfact = 1/(2*pi*sx*sy*sqrt(1-rho^2));
z = normfact*exp(-1/(2*(1-rho^2)).*(x.^2/sx^2+y.^2/sy^2-(2.*rho.*x.*y)/(sx*sy)))
%win = x.*y>0;
%z = z+(z.*win);
%min(z(:))
z = z./sum(z(:));

figure;
subplot(2,2,1);
imagesc(z); set(gca,'XTick',[],'YTick',[]);
axis square;
pdf1 = sum(z);
pdf2 = sum(z');
subplot(2,2,2);
plot(pdf1);
subplot(2,2,3);
plot(pdf2);
%%
% Section 16:
% a few chisquared distributions
figure;
subplot(1,3,1);
x = linspace(0,10,100);
plot(x,chi2pdf(x,1));
set(gca,'XLim',[0 x(end)]);

subplot(1,3,2);
x = linspace(0,25,100);
plot(x,chi2pdf(x,5));
set(gca,'XLim',[0 x(end)]);

subplot(1,3,3);
x = linspace(0,70,100);
plot(x,chi2pdf(x,30));
set(gca,'XLim',[0 x(end)]);


%%
% Section 17:
% a few F distributions
figure;
subplot(2,2,1);
x = linspace(0,10,200);
plot(x,fpdf(x,2,2));
set(gca,'XLim',[0 x(end)]);

subplot(2,2,2);
x = linspace(0,10,200);
plot(x,fpdf(x,2,20));
set(gca,'XLim',[0 x(end)]);

subplot(2,2,3);
x = linspace(0,10,200);
plot(x,fpdf(x,20,2));
set(gca,'XLim',[0 x(end)]);

subplot(2,2,4);
x = linspace(0,10,200);
plot(x,fpdf(x,20,20));
set(gca,'XLim',[0 x(end)]);
%%
% Section 18
% Looking at the relative efficiency of xbar and s^2 for estimating lambda

a = poissrnd(5,10,10000);
m = mean(a);
v = var(a);
bins = linspace(0,max([m';v']),50);
figure;
subplot(2,2,1);
h = bar([0:25],poisspdf([0:25],5));
set(h,'FaceColor','black','EdgeColor','black');

subplot(2,2,3);
[n,x] = hist(m,bins);
[mean(m) var(m)]
h = bar(x,n);
set(h,'FaceColor','yellow','EdgeColor','yellow');
set(gca,'Xlim',[0 20]);
subplot(2,2,4)
[n,x] = hist(v,bins);
[mean(v) var(v)]
h = bar(x,n);
set(h,'FaceColor','yellow','EdgeColor','yellow');
set(gca,'Xlim',[0 20]);


%%
% Section 19
% Comparing the relative efficiency of the sample mean and median from a
% normal distribution or a mixture distribution.

mu1=5;
mu2=10;
thresh = .90;  % 1/10 of measurements are bad

a = normrnd(mu1,1,10,10000);
m1= median(a);
m2 = mean(a);
bins = linspace(min([m1';m2']),max([m1';m2']),50);
figure;
subplot(2,1,1);
hist(m1,bins);
text(5.5,700,'median');
text(5.5,600,['mean: ',num2str(mean(m1))]);
text(5.5,500,['var: ',num2str(var(m1))]);

subplot(2,1,2)
hist(m2,bins);
text(5.5,700,'mean');
text(5.5,600,['mean: ',num2str(mean(m2))]);
text(5.5,500,['var: ',num2str(var(m2))]);

% Now from a mixture of normals
a = normrnd(mu1,1,10,10000);
b = normrnd(mu2,1,10,10000);
mask = unifrnd(0,1,10,10000) > thresh;
a(mask) = b(mask);
m1= median(a);
m2 = mean(a);
bins = linspace(min([m1';m2']),max([m1';m2']),50);
figure;
subplot(2,1,1);
hist(m1,bins);
text(7,700,'median');
text(7,600,['mean: ',num2str(mean(m1))]);
text(7,500,['var: ',num2str(var(m1))]);

subplot(2,1,2)
hist(m2,bins);
text(7,700,'mean');
text(7,600,['mean: ',num2str(mean(m2))]);
text(7,500,['var: ',num2str(var(m2))]);

% Here's the mixture distribution we're pulling from
figure;
x = linspace(0,14,1000);
y1 = normpdf(x,mu1,1);
y2 = normpdf(x,mu2,1);
y = thresh*y1+(1-thresh)*y2
plot(x,y,'linewidth',2);
%%
% Section 20
% Likelihood function for binomial
p = [0:.001:1];
y = binopdf(12,20,p);
figure
plot(p,y);

x = [0:20];
figure;
subplot(2,1,1);
bar(x,binopdf(x,20,.1));
subplot(2,1,2);
bar(x,binopdf(x,20,.4));

%%
% Section 21
% A pair of Gaussians for the large sample sign test
% First the distribution under Ho
x = [200:.1:350]
n = 500;
p = 0.5;
y = normpdf(x,n*p,sqrt(n*p*(1-p)));
figure; subplot(2,1,1);
plot(x,y,'LineWidth',2);
% Now the distribution under Ha
subplot(2,1,2);
p = 0.6;
y = normpdf(x,n*p,sqrt(n*p*(1-p)));
plot(x,y,'LineWidth',2);

%%
% Section 22
% Comparison of tail probabilities from binomial and normal (for some n)
n = 500;
p = .5;
x = logspace(-1, -4, 1000);

y = norminv(x,n*p,sqrt(n*p*(1-p)));
z = binocdf(y,n,p);
figure;
plot(x,x./z);
set(gca,'XScale','log');
set(gca,'XDir','reverse')
%%
% Section 23
% Showing how the power of a test increases with n using a pair of binomial
% distributions

n = [50,100,500];
p = [.5 .6];
figure;
for i = 1:length(n)
    x = [0:n(i)];
    y1 = binopdf(x,n(i),p(1));
    y2 = binopdf(x,n(i),p(2));
    tmp = y1+y2;    
    subplot(2,length(n),i);
    if (n(i) > 70);
        plot(x,y1,'Color','yellow','LineWidth',3);
        L = find(tmp > max(tmp)/1000);
        lims = [L(1) L(end)];
    else
        h = bar(x,y1,.3);
        set(h, 'FaceColor','yellow','EdgeColor','yellow');
        lims = [0 n(i)];
    end
    set(gca,'XLim',[lims(1)-.5 lims(2)+.5],'YLim',[0 max(tmp)*1.1]);
    
    subplot(2,length(n),i+length(n)); 
    if (n(i) > 70);
        plot(x,y2,'Color','yellow','LineWidth',3);
    else
        h = bar(x,y2,.3);
        set(h, 'FaceColor','yellow','EdgeColor','yellow');
    end
    set(gca,'XLim',[lims(1)-.5 lims(2)+.5],'YLim',[0 max(tmp)*1.1]);
end

%% 
% Section 24
% Stimulus movie
nframes = 1000;
BandW = 0;
figure;
axes;
clear M;
for i = 1:nframes
    a = uint8(normrnd(128,30,10,10,3));
    if (BandW == 1)
        a(:,:,2) = a(:,:,1); 
        a(:,:,3) = a(:,:,1); 
    end
    image(a);
    colormap(gray)
    axis square;

    set(gca,'XTick',[],'YTick',[],'Visible','off');
    axis image;
    %M(i) = im2frame(a);
    M(i) = getframe(gca);
end

repeat = 1;     %default = 1
pSearch = 0;    %default = 0
bSearch = 0;    %default = 1
reference = 1;  %default = 0
pixRange = 1;  %default = 10
iFrame = 1;     %default = 8
pFrame = 5;    %default = 10
bFrame = 5;    %default = 25

options = [repeat, pSearch, bSearch, reference, pixRange, iFrame, pFrame, bFrame];
mpgwrite(M, gray, 'testmpg.mpg', options)

%%
% Section 25 
% Skewed and kurtotic distributions
x = linspace(0,10,10000);
% Skewness
figure;
plot(x,gampdf(x,3,1),'linewidth',2);

% Kurtosis
figure;
df = 4;
x = linspace(-6,6,10000);
tpdf1 = tpdf(x,df-1);
tpdf1 = tpdf1./sum(tpdf1);
tpdf2 = normpdf(x,0,sqrt(df/(df-2)));
tpdf2 = tpdf2./sum(tpdf2);

subplot(2,2,3);
plot(x,tpdf1,'linewidth',2);
set(gca,'XLim',[x(1) x(end)],'YLim',[0 .0005]);

subplot(2,2,4);
plot(x,tpdf2,'linewidth',2);
set(gca,'XLim',[x(1) x(end)],'YLim',[0 .0005]);

% See if they're close
tpdf1*(x.^2)'  % "Variance"
tpdf2*(x.^2)'

%%
% Section 26
% A few stimulus frames
nframes = 20;
BandW = 1;
figure;
for i = 1:nframes
    a = uint8(normrnd(128,30,10,10,3));
    if (BandW == 1)
        a(:,:,2) = a(:,:,1); 
        a(:,:,3) = a(:,:,1); 
    end
    subplot(ceil(sqrt(10)),ceil(sqrt(10)),i);
    image(a);
    colormap(gray)
    axis square;

    set(gca,'XTick',[],'YTick',[],'Visible','off');
    axis image;
end

%%
% Section 27
% An example STA
% Pretending that the mean luminance is at 50 cd/m^2
% and that the maximum range of the monitor is 0-100 cd/m^2
% This probably isn't that far off.

filename = 'K051509005';  % Simple cell
%filename = 'K021709003' % Complex cell
BW = 1;  % if BW = 1 only look at the green channel

WN=nex2stro(findfile(filename));
nstixperside = WN.sum.exptParams.nstixperside;
ntrials = length(WN.sum.absTrialNum);
noisetypeidx = find(strcmp(WN.sum.trialFields(1,:),'noise_type'));
sigmaidxs = strmatch('sigma',WN.sum.trialFields(1,:));
sigmas = WN.trial(:,sigmaidxs)/1000;

spikename = getSpikenum(WN);
spikeidx = find(strcmp(WN.sum.rasterCells(1,:),spikename));
maxT = 9;
Lconenoise = WN.trial(:,noisetypeidx) == 2;

tmpstro = WN;
if (any(Lconenoise))  % Removing cone noise
    tmpstro.ras(Lconenoise,:) = [];
    tmpstro.trial(Lconenoise,:) = [];
end
out = getWhtnsStats(tmpstro,maxT,'STCOVmex', {nstixperside^2, 3, maxT}, spikename);
tmpstro = [];
STAs = out{1};
STCs = out{2};
nspikes = out{3};
%
% Plotting the STA and PCs
%
% Normalizing images
nstixperside = WN.sum.exptParams.nstixperside;
template = reshape([1:nstixperside^2],nstixperside,nstixperside);
edgepixels = [template(:,1); template(1,[2:end-1])'; template(end,[2:end-1])'; template(:,end)];

rowidxs = reshape([1:3*nstixperside^2],[nstixperside^2 3]);
maxes = []; mins = [];
for i = 1:3
    maxes = [maxes; max(max(STAs(rowidxs(:,i),:)))];
    mins = [mins; min(min(STAs(rowidxs(:,i),:)))];
end
potentialnormfactors = [(1-[.5; .5; .5]-eps)./(maxes-[.5; .5; .5]); (-[.5; .5; .5]+eps)./(mins-[.5; .5; .5])];
% 'eps' in above line is a kludge that is required for avoiding
% out of bounds errors.
potentialnormfactors(potentialnormfactors < 0) = []; % if min > mu or max < mu
normfactor = min(potentialnormfactors);
muvect = reshape(repmat([.5 .5 .5],nstixperside^2,1),nstixperside^2*3,1);

% gausslims = [WN.sum.exptParams.gauss_locut WN.sum.exptParams.gauss_hicut]/1000;
% NPOINTS = 65536;
% x = linspace(gausslims(1),gausslims(2),NPOINTS);
% Fx = norminv(x)*sigmas(1,2);
% sigmacorrectionfactor = std(Fx)./sigmas(1,2);
% muvar = (sigmas(1,2)*sigmacorrectionfactor)^2;

% Plotting
figure;
for i = 1:size(STAs,2)
    STA = normfactor*(STAs(:,i)-muvect)+muvect;
    STA = reshape(STA,[nstixperside nstixperside 3]);
    subplot(2,size(STAs,2),i);
    if (BW)
       STA(:,:,3) = STA(:,:,2);
       STA(:,:,1) = STA(:,:,2);
    end
    image(STA);
    set(gca,'XTick',[],'YTick',[]); axis square;
    
    % Significance tests
    STA = reshape(STAs(:,i),[nstixperside nstixperside 3]);
    STA = squeeze(STA(:,:,2));
    zmat = STA/(sigmas(1,2).*sqrt(nspikes));  % STA is really STS
    pmat = 2*(1-normcdf(abs(zmat),0,1));
    subplot(2,size(STAs,2),i+size(STAs,2));
    image((pmat<0.05)+1);
    c = colormap;
    c(1,:) = [0 0 0];
    c(2,:) = [1 1 0];
    colormap(c);
    set(gca,'XTick',[],'YTick',[]); axis square;    
end

% Plotting the null distribution

figure;
x = linspace(-6,6,20000);
y = normpdf(x);
plot(x,y);

% Getting the z-scores
whichframe = 3;
STA = reshape(STAs(:,whichframe),[nstixperside nstixperside 3]);
STA = squeeze(STA(:,:,2));
zmat = STA/(sigmas(1,2).*sqrt(nspikes));
figure;
imagesc(zmat);
hold on;
for i = 1:nstixperside
    for j = 1:nstixperside
        h=text(i,j,num2str(zmat(j,i),3));
        set(h,'HorizontalAlignment','center','FontSize',8);
    end
end
%%
% Section 28
% An example STV for a complex cell

filename = 'K021709003'; % Complex cell
BW = 1;  % if BW = 1 only look at the green channel

WN=nex2stro(findfile(filename));
nstixperside = WN.sum.exptParams.nstixperside;
ntrials = length(WN.sum.absTrialNum);
noisetypeidx = find(strcmp(WN.sum.trialFields(1,:),'noise_type'));
sigmaidxs = strmatch('sigma',WN.sum.trialFields(1,:));
sigmas = WN.trial(:,sigmaidxs)/1000;

spikename = getSpikenum(WN);
spikeidx = find(strcmp(WN.sum.rasterCells(1,:),spikename));
maxT = 9;
Lconenoise = WN.trial(:,noisetypeidx) == 2;

tmpstro = WN;
if (any(Lconenoise))  % Removing cone noise
    tmpstro.ras(Lconenoise,:) = [];
    tmpstro.trial(Lconenoise,:) = [];
end
out = getWhtnsStats(tmpstro,maxT,'STCOVmex', {nstixperside^2, 3, maxT}, spikename);
tmpstro = [];
STAs = out{1};
STCs = out{2};
nspikes = out{3};
STS2s = [];
for i = 1:size(STAs,2)
    STC = reshape(STCs(:,i),[sqrt(length(STCs(:,i))),sqrt(length(STCs(:,i)))]);
    STS2s(:,i) = diag(STC);
end

gausslims = [WN.sum.exptParams.gauss_locut WN.sum.exptParams.gauss_hicut]/1000;
NPOINTS = 65536;
x = linspace(gausslims(1),gausslims(2),NPOINTS);
Fx = norminv(x)*sigmas(1,2);
sigmacorrectionfactor = std(Fx)./sigmas(1,2);
muvar = (sigmas(1,2)*sigmacorrectionfactor)^2;
sigmavar = muvar*sqrt(2/nspikes);

% Plotting
figure;
for i = 1:size(STAs,2)
    STV = STS2s(:,i)./(sigmas(1,2)*sigmacorrectionfactor)^2;
    STV = reshape(STV,[nstixperside, nstixperside, 3]);
    subplot(2,size(STVs,2),i);
    if (BW)
       STV(:,:,3) = STV(:,:,2);
       STV(:,:,1) = STV(:,:,2);
    end
    image((STV-mean(STV(:)))./(2*range(STV(:)))+.5);
    set(gca,'XTick',[],'YTick',[]); axis square;
    
    % Significance tests
    
    pmat = 2*(1-chi2cdf(STV(:,:,2),nspikes));
    subplot(2,size(STS2s,2),i+size(STS2s,2));
    image((pmat<0.05)+1);
    c = colormap;
    c(1,:) = [0 0 0];
    c(2,:) = [1 1 0];
    colormap(c);
    set(gca,'XTick',[],'YTick',[]); axis square;    
end

% Plotting the null distribution

figure;
x = linspace(-nspikes/40,nspikes/40,20000)+nspikes;
y = chi2pdf(x,nspikes);
plot(x,y);
set(gca,'XLim',[x(1), x(end)]);

% Getting the z-scores
whichframe = 3;
STV = STS2s(:,whichframe)./(sigmas(1,2)*sigmacorrectionfactor)^2;
STV = reshape(STV,[nstixperside, nstixperside, 3]);
STV = STV(:,:,2);
figure;
imagesc(STV);
hold on;
for i = 1:nstixperside
    for j = 1:nstixperside
        h=text(i,j,num2str(STV(j,i),5));
        set(h,'HorizontalAlignment','center','FontSize',8);
    end
end

%% 
% Section 29
% A few spike-triggered stimulus frames

filename = 'K021709003'; % Complex cell
filename = 'K051509005';  % Simple cell
BW = 1;  % if BW = 1 only look at the green channel
nframesback = 4;   % Number of frames back in time
nSTSs = 20;

stro=nex2stro(findfile(filename));
nstixperside = stro.sum.exptParams.nstixperside;
gammaTable = stro.sum.exptParams.gamma_table;
gammaTable = reshape(gammaTable,length(gammaTable)/3,3);
invgammaTable = InvertGamma(gammaTable,1);
ngammasteps = size(invgammaTable,1);
seedidx = strcmp(stro.sum.trialFields(1,:),'seed');
noisetypeidx = strcmp(stro.sum.trialFields(1,:),'noise_type');
nframesidx = strcmp(stro.sum.trialFields(1,:),'num_frames');
stimonidx = strcmp(stro.sum.trialFields(1,:),'stim_on');
muidxs = [find(strcmp(stro.sum.trialFields(1,:),'mu1')),...
         find(strcmp(stro.sum.trialFields(1,:),'mu2')),...
         find(strcmp(stro.sum.trialFields(1,:),'mu3'))];
sigmaidxs = [find(strcmp(stro.sum.trialFields(1,:),'sigma1')),...
            find(strcmp(stro.sum.trialFields(1,:),'sigma2')),...
            find(strcmp(stro.sum.trialFields(1,:),'sigma3'))];
spikeidx = find(strcmp(stro.sum.rasterCells(1,:),getSpikenum(stro)));

msperframe = 1000/stro.sum.exptParams.framerate;
ntrials = size(stro.trial,1);

data = [];
i = 1;
while (size(data,2) < nSTSs)  % Getting 20 STSs
    seed = stro.trial(i,seedidx);
    noisetype = stro.trial(i,noisetypeidx);
    nframes = stro.trial(i,nframesidx);
    mu = stro.trial(i,muidxs)/1000;
    sigma = stro.trial(i,sigmaidxs)/1000;
    
    % Assuming gun noise
    x = linspace(stro.sum.exptParams.gauss_locut/1000, stro.sum.exptParams.gauss_hicut/1000, ngammasteps);   
    for gun = 1:3
        invnormcdf(:,gun) = norminv(x)*sigma(gun)+mu(gun);
    end
    randnums = getEJrandnums(3*nstixperside^2*nframes, seed);
    randnums = reshape(randnums, [nstixperside^2*3, nframes]);
    for j = 1:3
        idxs = [1:nstixperside^2]+nstixperside^2*(j-1);
        randnums(idxs,:) = reshape(invnormcdf(randnums(idxs,:)+1,j),[length(idxs),nframes]);
    end
    
    t_stimon = stro.trial(i, stimonidx);
    spiketimes = (stro.ras{i,spikeidx}-t_stimon)*1000;  % converting to ms
    frametimes = linspace(0, nframes*msperframe, nframes)+(msperframe/2)';
    spiketimes(spiketimes < nframesback*msperframe) = [];
    spiketimes(spiketimes > frametimes(end)) = [];
    [n,x] = hist(spiketimes, frametimes);
    idxs = find(n>0,nSTSs)-nframesback;
    data = [data, randnums(:,idxs)];
    i = i +1;
end
figure;
for i = 1:nSTSs
    a = reshape(data(:,i),[nstixperside, nstixperside, 3]);
    if (BW)
       a(:,:,1) = a(:,:,2);
       a(:,:,3) = a(:,:,2);
    end
    a = a+.5;
    subplot(ceil(sqrt(nSTSs)),ceil(sqrt(nSTSs)),i);
    image(a);
    set(gca,'Visible','off');
    axis square;
end

%% 
% Section 30
% Numerator and denomenator of the t-statistic

n = 3;
sigma = 1;
x = linspace(-3,3,1000);
y = normpdf(x,0,n*sigma);
figure;
plot(x,y,'y-','LineWidth',2);

figure;
x = linspace(0,15,1000);
y = chi2pdf(x,n-1);
plot(sigma*sqrt(x)/sqrt(n-1),y,'y-','LineWidth',2);

%% 
% Section 31
% A few t-distributions

x = linspace(-6,6,1000);

y1 = tpdf(x,2);
y2 = tpdf(x,10);
y3 = normpdf(x);

figure; axes; hold on;
plot(x,y1,'b-','LineWidth',2);
plot(x,y2,'r-','LineWidth',2);
plot(x,y3,'k-','LineWidth',2);

figure; axes; hold on;
t = sqrt(3)*(1.32-1)/0.25;
plot(x,y1,'b-','LineWidth',2);
plot(x,y3,'k-','LineWidth',2);
plot(t,0,'m*');
plot(-t,0,'m*');

% Here are the p-values for the t- and z-test
pz = 2*(1-normcdf(t))
pt = 2*(1-tcdf(t,n-1))

%% 
% Section 32
% Example of a paired t-test
x = [1.5 0.7 -.2 1.5 2.3];
mean(x);
std(x);
n = length(x);
t = sqrt(n)*mean(x)/std(x)
[h,p] = ttest(x);

x1 = linspace(-6,6,1000);
y1 = tpdf(x1,4);
figure; axes; hold on;
plot(x1,y1,'k-','LineWidth',2);

%% 
% Section 33
% Simulations comparing the coverage probability and power of the sign test
% and paired t-test
niter = 20000;
n = 10;
mu1 = 0;
mu2 = 1;
sigma = 1;
%2*(1-binocdf(8,10,.5))  % what critical value to use

data = zeros(niter,2);
for i = 1:niter
    x = normrnd(mu1,sigma,n,1);
    y = normrnd(mu2,sigma,n,1);
    [h,p1] = ttest(x-y);
   % k = sum(x-y > 0);
   % p2 = (2*binocdf(min(k, n-k),n,.5))-binopdf(k,n,.5);
   p2 = signtest(x-y);
   data(i,:) = [p1 p2];
end

ps = unique(data(:,2));
ns = zeros(length(ps),1);
for i = 1:length(ps)
    ns(i) = sum(data(:,2) == ps(i));
end

binswzero = [0;ps]
binwidths = diff(binswzero);
prob = ns./binwidths;
prob = prob./sum(prob);
count = niter.*prob;

crit = max(ps(ps<0.1));
% These should be close
sum(data(:,2) <= crit)/niter
crit

figure;
subplot(2,1,1);
for i = 1:length(binswzero)-1
    l = binswzero(i)+.005;
    r = binswzero(i+1)-.005;
    if (l< r)
        h = patch([l r r l l],[0 0 count(i) count(i) 0],'black');
    end
end

subplot(2,1,2);
[n,tmp] = hist(data(:,1));
bar(tmp,n,.9,'k');


sum(data(:,1) <= crit)/niter

%% 
% Section 34
% Simulations comparing the paired and unparied t-tests
% rho = .2, mu2 = 1, n = 5 gives more power to unpaired!
s2 = 1;
rho = .2;
mu1 = 0;
mu2 = 0;
n = 10;
niter = 40000;

data = zeros(niter,2);
A = normrnd(0,rho*s2,n,niter);
B = normrnd(mu1,(1-rho)*s2,n,niter);
C = normrnd(mu2,(1-rho)*s2,n,niter);
X=A+B;
Y=A+C;

[h,p1] = ttest(X-Y);
[h,p2] = ttest2(X,Y);
data = [p1' p2'];

figure;
subplot(2,1,1);
hist(data(:,1));
sum(data(:,1)<0.05)/niter

subplot(2,1,2);
hist(data(:,2));
sum(data(:,2)<0.05)/niter

%%
% Section 35
% Simulations showing that the t-test is robust to departures from
% normality
niter = 300000;
ns = [3:1:20];
data05 = [];
data01 = [];
for n = ns
    % First exponential
    randnums = exprnd(1,n,niter)-1; % mu = 0, v = 1
    [h,p1] = ttest(randnums);
    randnums = normrnd(0,1,n,niter); % mu = 0, v = 1
    [h,p2] = ttest(randnums);
    randnums = unifrnd(-sqrt(3),sqrt(3),n,niter); % mu = 0, v = 1
    [h,p3] = ttest(randnums);
    data05 = [data05; sum(p1<0.05) sum(p2<0.05) sum(p3<0.05)];
    data01 = [data01; sum(p1<0.01) sum(p2<0.01) sum(p3<0.01)];
end

data05./niter
data01./niter

figure;
plot(ns, data05./niter,'.-');
legend('exp','norm','unif');
title('p<0.05');
figure;
plot(ns, data01./niter,'.-');
legend('exp','norm','unif');
title('p<0.01');

% Now have to aplot the actual PDFs
figure;
x = [-4:.01:4];
subplot(3,1,1);
plot(x,exppdf(x+1,1),'LineWidth',2);
subplot(3,1,2);
plot(x,normpdf(x,0,1),'LineWidth',2);
subplot(3,1,3);
plot(x,unifpdf(x,-sqrt(3),sqrt(3)),'LineWidth',2);

%%
% Section 36
% Looking at the paired and unpaired t-tests in the case of anticorrelated
% X and Y
% ) Too many type I errors in unpaired t-test (- cov)
% ) paired t-test keeps a constant type I error

s2 = 1;
rho = -.8;
mu1 = 0;
mu2 = 0;
n = 10;
niter = 5000;

covmat = [s2, rho*s2;  rho*s2, s2];
[v,d] = eig(covmat);
M = v*sqrt(d);

data = zeros(niter,2);
for i = 1:niter
    X = normrnd(mu1,1,2,n);
    Xprime = M*X;
   
    Xprime = Xprime+repmat([mu1; mu2], 1, n);
    %plot(Xprime(1,:),Xprime(2,:),'k.');
    [h,p1] = ttest(Xprime(1,:)-Xprime(2,:));
    [h,p2] = ttest2(Xprime(1,:),Xprime(2,:));
    data(i,:) = [p1 p2];
end
figure;
subplot(2,1,1);
hist(data(:,1),20);
subplot(2,1,2);
hist(data(:,2),20);
%%
% Section 37
% Simulations showing the effects of correlations among variables on the
% outcome of the t-test.

s2 = 1;
rho = 0.1;
mu1 = 0;
mu2 = 0;
ns = 10;
niter = 300000;
ns = [3:1:20];
data05 = [];
data01 = [];

for n = ns
    covmat = repmat(rho,n,n)+eye(n)-rho*eye(n);
    [v,d] = eig(covmat);
    M = v*sqrt(d);
    
    X = normrnd(mu1,1,n,niter);
    Xprime = M*X;
    [h,p1] =ttest(Xprime);
    [h,p2] =ttest(X);
    %hist(p,20);
   % corrcoef(Xprime')
    data05 = [data05; sum(p1<0.05)./niter sum(p2<0.05)./niter];
    data01 = [data01; sum(p1<0.01)./niter sum(p2<0.01)./niter];
end

plot(ns, data05);

%%
% Section 38
% An example nonlinear regression
b = [30 8 2 -.13]/10;
x = linspace(0,18,100);
X = [ones(1,size(x,2)); x; x.^2; x.^3];
eyx = b*X;

xsamples = unifrnd(0,x(end),100,1);
X = [ones(size(xsamples,1),1), xsamples, xsamples.^2, xsamples.^3];
linpart = X*b'
ysamples = exprnd(linpart/2)+linpart/2

figure;
plot(xsamples,ysamples,'k.','MarkerSize',15);
hold on;
plot(x,eyx)

%%
% Section 39
% A univariate regression problem with a t-test/CI on the slope.

n = 10;
beta = [10 5];
%sigma = 10;
%x = unifrnd(0, 10, n,1);

x = [6.2790;...
    7.7198;...
    9.3285;...
    9.7274;...
    1.9203;...
    1.3887;...
    6.9627;...
    0.9382;...
    5.2540;...
    5.3034];

y = [63.6665;...
   47.9069;...
   51.5694;...
   60.9951;...
   22.0595;...
   17.6442;...
   38.7275;...
    2.4651;...
   39.4352;...
   23.0885]';

X = [ones(n,1), x];
%y =beta*X'+normrnd(0,sigma,1,n);

b = inv(X'*X)*X'*y';
yhat = X*b;
s = sqrt(sum((y'-yhat).^2)/(n-2));
seb1 = s/sqrt(sum((x-mean(x)).^2));

tcrit = tinv(.975, n-2);
CI = [b(2)-seb1*tcrit, b(2)+seb1*tcrit];

% Sanity check
[bnew,bint] = regress(y',X);

figure;
plot(x,y,'k.','MarkerSize',15);
hold on;
plot([0 10],[b(1) b(1)+10*b(2)],'k-');
title(['slope = ',num2str(b(2)),' 95% CI = ',num2str(CI)])

figure; hold on;
plot(linspace(-8*seb1,8*seb1,1000), tpdf(linspace(-8*seb1,8*seb1,1000), n-2));
plot(bnew(2),0,'y*');
plot(tcrit*seb1*[-1 1],[0 0],'m*');

p = 2*(1-tcdf(b(2)/seb1, n-2));
title(['p = ',num2str(p)]);

%%
% Section 40
% Partitioned sums of squares in a regression 

n = 10;
beta = [10 5];
sigma = 10;
%x = unifrnd(0, 10, n,1);

x = [6.2790;...
    7.7198;...
    9.3285;...
    9.7274;...
    1.9203;...
    1.3887;...
    6.9627;...
    0.9382;...
    5.2540;...
    5.3034];

y = [63.6665;...
   47.9069;...
   51.5694;...
   60.9951;...
   22.0595;...
   17.6442;...
   38.7275;...
    2.4651;...
   39.4352;...
   23.0885]';

X = [ones(n,1), x];
%y =beta*X'+normrnd(0,sigma,1,n);


b = inv(X'*X)*X'*y';
yhat = X*b;

figure;
plot(x,y,'k.','MarkerSize',15);
hold on;
plot([0 10],[b(1) b(1)+10*b(2)],'k-');
plot([0 10],[mean(y) mean(y)],'m-');

% F test stuff
num = sum((yhat-repmat(mean(y),n,1)).^2);
denom = sum((y'-yhat).^2)/(n-2);
Fstat = num./denom;
% Sanity check
[tmp1, tmp2, tmp3, tmp4, stats] = regress(y',X);

figure; axes; hold on;
plot(linspace(0,50,1000),fpdf(linspace(0,50,1000),1,n-2),'k-');
plot(Fstat,0,'y*');
title(num2str(1-fcdf(Fstat,1,n-2)));
%%
% Section 41
% Regressions with various regression sums of squares and error sums of
% squares.

n = 20;
x = unifrnd(0, 10, n,1);

slopes = [5 1 5 1];
intercepts = [10 30 10 30];
sigmas = [10 3 3 10];
for i = 1:4
    subplot(2,2,i); hold on;
    beta = [intercepts(i) slopes(i)];
    sigma = sigmas(i);
    X = [ones(n,1), x];
    y = beta*X'+normrnd(0,sigma,1,n);
    plot(x,y,'k.','MarkerSize',17);
    b = regress(y',X);
    plot([0 10],[b(1) b(1)+10*b(2)],'k-');
    plot([0 10],[mean(y) mean(y)],'k:');
    set(gca,'YLim',[0 80]);
end

%%
% Section 42
% Residuals from the canonical univariate regression

x = [6.2790;...
    7.7198;...
    9.3285;...
    9.7274;...
    1.9203;...
    1.3887;...
    6.9627;...
    0.9382;...
    5.2540;...
    5.3034];

y = [63.6665;...
   47.9069;...
   51.5694;...
   60.9951;...
   22.0595;...
   17.6442;...
   38.7275;...
    2.4651;...
   39.4352;...
   23.0885]';

X = [ones(length(x),1), x];

b = inv(X'*X)*X'*y';
yhat = X*b;
residual = y-yhat';
figure; axes; hold on;
plot(x,residual,'k.','MarkerSize',15);
plot([0 10],[0 0],'k:');

%%
% Section 43
% Problematic patterns of residuals

n = 30;
p = .2;
sigma1 = 5;
sigma2 = 40;

figure;
% Mixture of Gaussians
subplot(3,1,1); hold on;
x = unifrnd(0,10,n,1);
res = normrnd(0,sigma1,n,1);
L = logical(binornd(1,p,n,1));
res(L) = normrnd(0,sigma2,sum(L),1);
res = res - mean(res);
plot(x,res,'k.','MarkerSize',15);
plot([0 10],[0 0],'k:');

% Heteroskedasticity
subplot(3,1,2); hold on;
res = normrnd(0,sigma1*x/2,n,1);
res = res - mean(res);
plot(x,res,'k.','MarkerSize',15);
plot([0 10],[0 0],'k:');

% Nonlinear relationship
subplot(3,1,3); hold on;
res = normrnd(4*(x-5).^2,2*sigma1,n,1);
res = res - mean(res);
plot(x,res,'k.','MarkerSize',15);
plot([0 10],[0 0],'k:');

maxlim = 0;
for i = 1:3
    subplot(3,1,i);
    a = get(gca,'YLim');
    maxlim = max(maxlim, max(abs(a)));
end;
for i = 1:3
    subplot(3,1,i);
    set(gca,'YLim',[-maxlim maxlim]);
end

%% 
% Section 44
% Regression on transformed variables
n = 30;
b = [2 .5];
sigma = .5;
x = unifrnd(0,10,n,1);
X = [ones(length(x),1), x];
y = (b*X')+normrnd(0,sigma,1,n);
%y = (2*x.^1.8)+normrnd(0,sigma,n,1);
figure(1); axes;
hold on;
plot(x,exp(y),'k.')
figure(2); axes;
hold on;
plot(x,y,'k.')

b = regress(y',X);
plot([0 10],[b(1) b(1)+10*b(2)],'k-');

figure(1);
domain = linspace(0,10,1000);
plot(domain, exp(b(2)*domain+b(1)),'k-');

%%
% Section 45
% Polynomial regressions
n = 10;
b = [2 5 -4 .35];
sigma = 6;
x = unifrnd(0,10,n,1);
x = linspace(0,10,n)'
X = [ones(length(x),1), x, x.^2, x.^3];
y = (b*X')+normrnd(0,sigma,1,n);
X = [ones(length(x),1), x, x.^2, x.^3, x.^4, x.^5, x.^6, x.^7, x.^8, x.^9, x.^10];

domain = linspace(0,10,1000);
plotX = [ones(length(domain),1), domain', domain'.^2, domain'.^3, domain'.^4, domain'.^5, domain'.^6, domain'.^7, domain'.^8, domain'.^9, domain'.^10];

figure;
axes; hold on;
plot(x,y,'k.','MarkerSize',15);
bhat = regress(y',X(:,[1:2]));
plot(domain,plotX(:,[1:length(bhat)])*bhat,'m-')
bhat = regress(y',X(:,[1:3]));
plot(domain,plotX(:,[1:length(bhat)])*bhat,'g-')
bhat = regress(y',X(:,[1:4]));
plot(domain,plotX(:,[1:length(bhat)])*bhat,'b-')
bhat = regress(y',X);
plot(domain,plotX(:,[1:length(bhat)])*bhat,'k-')
set(gca,'Ylim',1.3*[min(y) max(y)]);

%%
% Section 46
% Multicollinearity nearly singular regression
sigma1 = .5;
sigma2 = .001;
meanweightinkg = 10;
n = 50;

x1 = normrnd(meanweightinkg,sigma1,n,1)+normrnd(0,sigma2,n,1);  % weight in kg
x2 = (x1*2.2)+normrnd(0,sigma2,n,1);
y = x1+normrnd(0,sigma1,n,1);
[b1,bint1] = regress(y,[ones(n,1) x1]);
[b2,bint2] = regress(y,[ones(n,1) x2]);
[b3,bint3] = regress(y,[ones(n,1) x1 x2]);

% Now a little simulation
niter = 10000;
data = zeros(niter,3);
for i = 1:niter
    x1 = normrnd(meanweightinkg,sigma1,n,1)+normrnd(0,sigma2,n,1);  % weight in kg
    x2 = (x1*2.2)+normrnd(0,sigma2,n,1);
    y = x1+normrnd(0,sigma1,n,1);
    [b1,bint1] = regress(y,[ones(n,1) x1]);
    [b2,bint2] = regress(y,[ones(n,1) x2]);
    [b3,bint3] = regress(y,[ones(n,1) x1 x2]); 
    data(i,:) =[b1(2) b3(2) b3(3)];
end

figure;
hist(data(:,1),[.4:.05:1.6])

figure;
hist(data(:,2),linspace(-800, 800,30))
set(gca,'XLim',[-800 800]);

figure;
hist(data(:,3),linspace(-300, 300,30))
set(gca,'XLim',[-300 300]);

figure;
plot(data(:,2),data(:,3),'k.','MarkerSize',15);
axis equal;

%%
% Section 47
% iid uniform process
x = unifrnd(-1,1,30,1)
subplot(2,1,1);
plot(x,'k.','MarkerSize',15);
set(gca,'YLim',[-1 1]);
%%
% Section 48
% autoregressive uniform process
n = 30;
phi = .3;
x = unifrnd(-1,1,n,1)
y = zeros(size(x));
y(1) = x(1);
for i = 2:n    
    y(i) = phi*y(i-1)+x(i);
end
figure
subplot(2,1,1);
plot(y,'k.','MarkerSize',15);

%%
% Section 49
% Curse of dimensionality
N = 100;
n = 1000;
x = normrnd(0,1,N,n);
meandist = zeros(N,n)
for i = 1:N
   for j = 1:n
      tmp = repmat(x(1:i,j),[1,size(x,2)])-x(1:i,:);
      tmp(:,j) = [];
      if (i == 1)
          dists = abs(tmp);
      else
        dists = sqrt(sum(tmp.^2));
      end
      meandist(i,j) = mean(dists);
   end
end
figure;
subplot(2,1,1);
plot(mean(meandist,2),'k-');
set(gca,'ylim',[0 15]);

subplot(2,1,2);
plot(1:N,sqrt(1:N),'k-');
set(gca,'ylim',[0 15]);

%%
% Section 50
% An example of linear discriminant analysis

n = 100;
sigma = .7;
Sigma = [1 sigma; sigma 1];
mu1 = [-1 0]';
mu2 = [1 -1]';
x = normrnd(0,1,2,n);
y = normrnd(0,1,2,n);
x = sqrtm(Sigma)*x;
x = x+repmat(mu1,1,n);

y = sqrtm(Sigma)*y;
y = y+repmat(mu2,1,n);

figure; axes; hold on;
plot(x(1,:),x(2,:),'r.','MarkerSize',20);
plot(y(1,:),y(2,:),'g.','MarkerSize',20);
axis square
set(gca,'XLim',[-4 4],'YLim',[-4 4]);
w = inv(Sigma)*(mu1-mu2);
centroid = (mu1+mu2)/2;

plot([-3 3]*w(2)+centroid(1),[-3 3]*-w(1)+centroid(2))

%%
% Section 51
% Figures for the 1st (only?) problem set
% a) A 2-D inseparable PDF
% b) a 1-D PMF
% c) a 1-D PDF (expontential)
% d) 
theta = -pi/4;
xaxistick = linspace(-10,30,500);
yaxistick = linspace(-10,30,500);

[x,y] = meshgrid(xaxistick,yaxistick);
xprime = x.*cos(-theta)+y.*sin(-theta);
%yprime = -x.*sin(-theta)+y.*cos(-theta);
jointpdf = gampdf(xprime,3,2).*normpdf(x,1,2);
imagesc(jointpdf)
axis square;
axis image;
axis xy;
xlabel('X'); ylabel('Y');

figure;
plot(sum(jointpdf)); % x marginal
figure;
plot(jointpdf(:,150)); % dist Y|x=150

% 1-D PMF
x = [1 2 3 4 5];
pmf = [1 2 3 2 1];
pmf = pmf./sum(pmf);
plot(x,pmf);
mn = x*pmf' % mean
(x-mn)*pmf' % mean absolute deviation from mean
(x-mn).^2*pmf' % mean absolute deviation from mean
sqrt(x)*pmf'  % E(sqrt(X));
%
% 1-D PDF (exponential)
x = linspace(0,30,1000);
figure;
plot(x,expdf(x,3));

%%
% Section 52
% Simulations showing that estimating the standard deviation in the z
% statistic creates a t distribution, not a standard normal.
niter = 20000;
n = 3;
x = normrnd(0,1,niter,n)
z = sqrt(n)*mean(x,2);
t = sqrt(n)*mean(x,2)./std(x,0,2)
figure; axes; hold on;
[n,x] = hist(z,[-6:.2:6]);
h = bar(x,n./max(n),'y');
set(h,'FaceCOlor','yellow');
tmpx = [-6:.01:6];
tmpy = normpdf(tmpx,0,1);
plot(tmpx,tmpy./max(tmpy),'k-','linewidth',3);

figure; axes; hold on;
[n,x] = hist(t,[-6:.2:6]);
n(1) = 0; n(end) = 0;
h = bar(x,n./max(n),'y');
set(h,'FaceCOlor','yellow');
tmpx = [-6:.01:6];
tmpy = normpdf(tmpx,0,1);
plot(tmpx,tmpy./max(tmpy),'k-','linewidth',3);

%%
% Section 53
% Null distribution for the Wilcoxon signed-rank statistic
n = 5
allposs = (ff2n(n))';
idx = (1:n)';
idx = idx(:,ones(2.^n,1));
pranks = sum(allposs.*idx,1);

%%
% Section 54
% Simulated end points of dart tosses/saccade endpoints (for permutation test)
n1 = 20;
n2 = 10;
z1 = unidrnd(2,n1,2)*2-3;
x1 = exprnd(1,n1,1).*z1(:,1);
y1 = exprnd(1,n1,1).*z1(:,2);

z2 = unidrnd(2,n2,2)*2-3;
x2 = exprnd(1,n2,1).*z2(:,1)+.7;
y2 = exprnd(1,n2,1).*z2(:,2)-.5;
figure; axes; hold on;
plot(x1,y1,'ko','MarkerSize',10,'MarkerFaceColor','black','MarkerEdgeColor','black');
plot(mean(x1),mean(y1),'kx','MarkerSize',10);
plot(x2,y2,'ko','MarkerSize',10,'MarkerFaceColor','yellow','MarkerEdgeColor','yellow');
plot(mean(x2),mean(y2),'yx','MarkerSize',10);
truedist = sqrt((mean(x1)-mean(x2)).^2+(mean(y1)-mean(y2)).^2);
axis equal;
title('unshuffled');

whichgroup = [zeros(n1,1); ones(n2,1)];
niter = 2000;
dists = zeros(niter,1);
allx = [x1; x2];
ally = [y1; y2];
print -dpsc junk

for i = 1:niter 
    tmp = logical(whichgroup(randperm(length(whichgroup))));
    if (i < 5)
        figure; axes;hold on;
        plot(allx(tmp),ally(tmp),'ko','MarkerSize',10,'MarkerFaceColor','black','MarkerEdgeColor','black');
        plot(mean(allx(tmp)),mean(ally(tmp)),'kx','MarkerSize',10);
        plot(allx(~tmp),ally(~tmp),'ko','MarkerSize',10,'MarkerFaceColor','yellow','MarkerEdgeColor','yellow');
        plot(mean(allx(~tmp)),mean(ally(~tmp)),'yx','MarkerSize',10);
        axis equal;
        eval(['print -dpsc junk',num2str(i)]);
    end
    dists(i) = sqrt((mean(allx(tmp))-mean(allx(~tmp))).^2+(mean(ally(tmp))-mean(ally(~tmp))).^2);
end
figure; axes; hold on;
[n,x] = hist(dists,20);
bar(x,n,'yellow');
plot(truedist,0,'b^')
title(['p = ',num2str(sum(dists>truedist)/niter)]);

%%
% Section 55
% Bootstrapping intro
n = 10;
x = normrnd(0,1,n,1);
bins = [-4:.5:4];
[counts,~] = hist(x,bins);
figure; axes;
bar(bins,counts,'yellow');
title('real')
print -dpsc junk;

for i = 1:5
    bootdata = x(unidrnd(n,n,1));
    [counts,~] = hist(bootdata,bins);
    figure; axes;
    bar(bins,counts,'yellow');
    eval(['print -dpsc junk',num2str(i)]);
end

% More bootstrapping
n = 20;
bins = linspace(0,4,20);
x = linspace(0,4,200);
p = wblpdf(x,2,4);
cdf = wblcdf(linspace(0,4,200),2,4);
medianidx = find((cdf-.5).^2 == min((cdf-.5).^2));
median = x(medianidx);
figure;axes; hold on;
plot(x,p,'k-','color','yellow','linewidth',3);
plot(median,0,'m^');
title(['median = ',num2str(median)]);
print -dpsc junk0

% the "real" data
figure; axes; hold on;
dat = wblrnd(2,4,n,1);
[counts,~] = hist(dat,bins);
bar(bins,counts,'yellow');
title(['median (real sample) = ',num2str(prctile(dat,50))]);
plot(prctile(dat,50),0,'m^');
print -dpsc junk00

% Bootstrapping
niter = 5000;
bootmeds = zeros(niter,1);
for i = 1:niter
    bootdat = dat(unidrnd(n,n,1));
    if (i < 5)
        [counts,~] = hist(bootdat,bins);
        figure; axes; hold on;
        bar(bins,counts,'yellow');
        plot(prctile(bootdat,50),0,'m^');
        title(['median = ',num2str(prctile(bootdat,50))]);
        eval(['print -dpsc junk',num2str(i)]);
    end
    bootmeds(i) = prctile(bootdat,50);
end

figure; axes; hold on;
[counts,~] = hist(bootmeds,bins);
bar(bins,counts,'yellow');
plot(prctile(bootdat,50),0,'m^');
title(['STD = ',num2str(std(bootmeds)),' bias = ',num2str(mean(bootmeds)-mean(dat))])
print -dpsc junk10

%%
% Section 60
% Distributions of difference of rvs.
npts = 200;
for i = 1:5
    t = linspace(-6,6,npts);
    if (i==1)
        x = normpdf(t,0,1);
        y = normpdf(t,0,1);
    elseif (i == 2)
        x = exppdf(t+1,1);
        y = exppdf(-t+1,1);  
    elseif (i == 3)
        x = chi2pdf(2*t+4,3);
        y = chi2pdf(-2*t+4,3);
    elseif (i == 4)
        x = betapdf(t/6+.5,.5,.5);
        y = betapdf(-t/6+.5,.5,.5);
    elseif (i == 5)
        x = unifpdf(t,-5,5);
        y = unifpdf(-t,-5,5);
    end
    figure;
    subplot(2,1,1);
    plot(t,x,'LineWidth',3);
    set(gca,'XLim',[-10 10]);
    subplot(2,1,2);
    z = conv(x,y);
    plot(linspace(2*t(1), 2*t(end),length(z)),z,'Linewidth',3)
    set(gca,'XLim',[-10 10]);
end    

%%
% Section 61
niter = 1000;
n = 20;
p = zeros(niter,2);
for i =1:niter
    p(i,1)= signrank(exprnd(1,n,1)-log(2));
    p(i,2)= signrank(normrnd(0,1,n,1));
end

%%
% Section 56
% Bootstrapping the SE of the median of a normal distribution
% The original data set:
n = 30;
x = normrnd(0,1,n,1);
bounds = max(abs(x));
bins = linspace(-bounds,bounds,10);
[n,~] = hist(x,bins);
figure; axes; hold on;
h = bar(bins,n);
set(h,'FaceColor','yellow');
plot(median(x),0,'k*'); title(['full data ',num2str(median(x))]);
% Now bootstrapping
nboot = 200;
data = [];
for i = 1:nboot
    y = x(unidrnd(length(x),length(x),1));
    
    if (i<3)
        figure; axes; hold on;
        [n,~] = hist(y,bins);
        h = bar(bins,n);
        set(h,'FaceColor','yellow');
        plot(median(y),0,'k*'); title(num2str(median(y)));
    end
    data(i) = median(y);
end
figure; axes; hold on;
[n,~] = hist(data,linspace(bins(1),bins(end),30));
h = bar(linspace(bins(1),bins(end),30),n);
set(h,'FaceColor','yellow','LineWidth',1);
title(num2str(std(data)));

% Getting the true value (as close as I can figure)
x = normrnd(0,1,n,1000000);
std(median(x))

%%
% Section 57
% Bootstrap t-tests
nbigloops = 30;
niter = 10000;
nboot = 5000;
ns = [3:1:20];
data = zeros(length(ns),3);
for i = 1:nbigloops
    i
    for j = 1:length(ns)
        n = ns(j)
        % First exponential
        randnums = exprnd(1,n,niter)-1; % mu = 0, v = 1
        ci = bootci(nboot,@(x)(sqrt(n)*mean(x)./std(x)),randnums);
        h1 = sum(sign(ci(1,:)) == sign(ci(2,:)));
        randnums = normrnd(0,1,n,niter); % mu = 0, v = 1
        ci = bootci(nboot,@(x)(mean(x)),randnums);
        h2 = sum(sign(ci(1,:)) == sign(ci(2,:)));
        randnums = unifrnd(-sqrt(3),sqrt(3),n,niter); % mu = 0, v = 1
        ci = bootci(nboot,@(x)(mean(x)),randnums);
        h3 = sum(sign(ci(1,:)) == sign(ci(2,:)));
        data(j,:) = data(j,:)+[sum(h1) sum(h2) sum(h3)];
    end
end
figure;
plot(ns, data./(niter*nbigloops),'.-');
legend('exp','norm','unif');
title('p<0.05');
set(gca,'Ylim',[0.045 .125])

% Now have to aplot the actual PDFs
figure;
x = [-4:.01:4];
subplot(3,1,1);
plot(x,exppdf(x+1,1),'LineWidth',2);
subplot(3,1,2);
plot(x,normpdf(x,0,1),'LineWidth',2);
subplot(3,1,3);
plot(x,unifpdf(x,-sqrt(3),sqrt(3)),'LineWidth',2);
