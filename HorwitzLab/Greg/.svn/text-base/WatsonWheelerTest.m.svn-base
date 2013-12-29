function [h,p] = WatsonWheelerTest(in1,in2, alpha)
% [h,p] = function WatsonWheelerTest (in1,in2, <alpha>)
%   Test of mu1 = mu2 for circular distributions.
%   See Zar, p.633 This test is probably going to be
%   superior to the WatsonWilliams test, which appears to 
%   be applicable only when the underlying distributions 
%   have very low dispersion.  
%
%   Note: Doesn't work with ties or with n<10  
%
%   GDLH 3/23/09

if (nargin < 3)
    alpha = 0.05;
end
if (size(in1,1) < size(in1,2))
    in1 = in1';
end
if (size(in2,1) < size(in2,2))
    in2 = in2';
end
n1 = size(in1,1);
n2 = size(in2,1);
n = n1+n2;
if (n < 10)
    disp('Warning n > 10.  Do not trust p-value.');
end

alldata = [in1; in2];
idxs = randperm(length(alldata));  % for breaking up ties.
%idxs = 1:length(alldata);
% Otherwise all the tied values from 'in1' are considered to be
% smaller than the identical values from 'in2'.
[y,i] = sort(alldata(idxs));
L = logical([ones(n1,1); zeros(n2,1)]);
d = 2*pi*idxs(i)./n;
c1 = sum(cos(d(L)));
s1 = sum(sin(d(L)));
c2 = sum(cos(d(~L)));
s2 = sum(sin(d(~L)));

W = 2*((c1^2+s1^2)/n1 + (c2^2+s2^2)/n2);
p = 1-chi2cdf(W,2);
h = p<alpha;
end
% 
% % Some stuff for testing
% ps = []; rs = [];
% for iter = 1:10000
%     n1 = round(normrnd(100,10));
%     n2 = round(normrnd(100,10));
% 
%     x = unifrnd(0,2*pi,n1,1);
%     y = unifrnd(0,2*pi,n2,1);
%     [h,p] = WatsonWheelerTest(x,y);
%     ps = [ps; p];
%     rs = [rs; h];
% end