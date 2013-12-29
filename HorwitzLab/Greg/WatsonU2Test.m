function [h,p] = WatsonU2Test(in1,in2, alpha)
% [h,p] = function WatsonU2Test (in1,in2, <alpha>)
% Nonparametric test on the median of two circular distributions.
% Recommended when at least one of the populations
% is not unimodal and when there are ties.
% See Zar, p.630.  Assuming both sample sizes > 100 individually.
%
% This implementation returns p = 1 for any p > 0.5 because I don't 
% have critical values for p > 0.5 (see Zar Table B.38 p. App197)
%
% GDLH 4/7/09

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

alldata = mod([in1; in2],2*pi);
L = logical([ones(n1,1); zeros(n2,1)]);
[b,i] = sort(alldata);
cdf1 = cumsum(L(i))./sum(L);
cdf2 = cumsum(~L(i))./sum(~L);

d = sum(cdf1-cdf2);
d2 = sum((cdf1-cdf2).^2);
U2 = ((n1*n2)/n^2)*(d2-d^2/n);

ptable = [0:.001:1];
crittable= -(1/(2*pi^2)).*(log(ptable./2)-log(1+(ptable./2).^3));
p = interp1(crittable, ptable, U2);
h = p<alpha;


% 
 % Some stuff for testing
%  ps = [];
%  for iter = 1:10000
%      n1 = round(normrnd(200,10));
%      n2 = round(normrnd(150,10));
%  
%      x = unifrnd(0,2*pi,n1,1);
%      y = unifrnd(0,2*pi,n2,1);
%      [h,p] = WatsonU2Test(x,y);
%      ps = [ps; p];
%  end