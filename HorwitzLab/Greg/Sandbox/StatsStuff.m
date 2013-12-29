% Section 1
% Code to compute the count distribution for a renewal process with
% Gamma-distributed intervals.
%
% Section 2
% A test of the formula to compute bounds on r23 given r12 and r23.
%
% Section 3
% Playing around with Cauchy random variables and confidence intervals for
% the ratio of two independent normals.

%%
% Section 1
%
% Code to compute the count distribution for a renewal process with
% Gamma-distributed intervals.

% User-defined parmeters:
t = 50;  % Counting window length
a = 5; % First parameter of Gamma distn. (mean of exponentials).
b = 5; % Second parameter of Gamma (number of exponentials convolved).
maxspikes = 5*t/(a*b); 
% Largest number of spikes that we're going to consider as possible in a
% window of length 't'.  No just a big number - no principled reason for
% this choice.

%----------------------------
% Computing 1-CDF of count distribution using fact that:
% P(N(t)>=k) = P(Wk<=t)
% where N(t) is count in a window of length t
% and Wk is waiting time to kth renewal (in this case gamma with parameter
% k*b.
x = zeros(maxspikes,1);
for k = 1:maxspikes
    x(k) = gamcdf(t,k*b,a);
end
x = [1;x];  % P(N(t) >= 0) = 1
% x(k) is now P(N(t) > k)
countcdf = 1-x;
countpdf = diff(countcdf);
bar(countpdf);  % here it is

% As a sanity check we comparing this
% distribution against the Poisson (set b = 1 for Poisson process)
if (b == 1)
    pdf = poisspdf([0:maxspikes],t/a);
    hold on;
    plot(pdf,'y*');
end

%%
% Section 2

% Assume three random variables.  If we fix the correlation between
% variables 1 and 2 and we fix the correlation between variables 1 and 3,
% what can we say about the allowable correlations between variables 2 and
% 3?

r12 = -.9;
r13 = 0;
% Here are the bounds on r23 (the correlation between rvs 2 and 3)
bnds = [cos(acos(r12)-acos(r13)) cos(acos(r12)+acos(r13))];

% Now let's test it to make sure i haven't gotten it wrong.
% Building a series of correlation matrices with the values of r12 and r13
% given above.  Then, trying many different values for r23 and seeing if we
% get a legitimate correlation matrix (postive semi-definite).
cormat = eye(3);
cormat(1,2) = r12; cormat(2,1) = r12;
cormat(1,3) = r13; cormat(3,1) = r13;

testrs = linspace(-1,1,100);
success = nan*ones(1,length(testrs));
for i = 1:length(testrs)
    cormat(2,3) = testrs(i); cormat(3,2) = testrs(i);
    success(i) = det(cormat)>=0;
end

figure; axes; hold on;
plot(testrs,success);
plot(bnds,[1 1],'m*');

%%
% Section 3
% Playing around with Cauchy random variables and confidence intervals for
% the ratio of two independent normals.

muX = 10;
muY = 40;
sigmaX = 10;
sigmaY = 4;

% First, a Monte Carlo simulation
niter = 10000;
X = normrnd(muX, sigmaX, niter, 1);
Y = normrnd(muY, sigmaY, niter, 1);
Z = X./Y;
figure;
[n,x] = hist(Z,100);
bar(x,n./max(n));
hold on;

% Now using the formulas from Kamerund 1978
z = linspace(min(Z),max(Z),1000);
w = (sigmaY/sigmaX)*z;  % w is a scaled version of z
s = 1./sqrt(w.^2+1);
k = (muX/sigmaX.*w+muY/sigmaY).*s.^2;
M = -.5*(muY/sigmaY.*w-muX/sigmaX).^2.*s.^2;
Q = k.*s.*sqrt(2*pi).*(1-2.*normcdf(-k./s))+(2*s.^2.*exp(-k.^2./(2*s.^2)));

fg = 1/(2*pi).*Q.*exp(M);
fz = (sigmaY/sigmaX)*fg;
fz = fz./sum(fz);
plot(z,fz./max(fz),'m-');

% CDF
Fz = cumsum(fz);
CI = interp1(Fz,z,[.025 .975])
