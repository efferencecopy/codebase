function [h,p] = WatsonWilliamsTest(in1,in2, alpha)
% [h,p] = function WatsonWilliamsTest (in1,in2, <alpha>)
%   Test of mu1 = mu2 for circular distributions.
%   See Zar, p.625
%
%   GDLH 3/16/09

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

x1 = sum(cos(in1));
y1 = sum(sin(in1));
x2 = sum(cos(in2));
y2 = sum(sin(in2));

R = (n1+n2)*norm([(x1+x2)./n (y1+y2)./n]);
R1 = n1*norm([x1./n1 y1./n1]);
R2 = n2*norm([x2./n2 y2./n2]);

stat = ((n1+n2-2)*(R1+R2-R))/(n1+n2-R1-R2);
rw = (R1+R2)/(n1+n2);
if  rw < .7 || n/2 < 10
    disp('Test not applicable. Number of samples to low or average resultant vector length to low.')
    p = nan;
    h = nan;
    return;
end

% There's a lookup table (B.37) on page App193 of Zar that I'm trying to implement
% here.
rs = [0.001:0.001:0.049];
ks = [188.4989 94.7472 63.5015 47.8749 38.4992 32.2498 27.7851 24.4367 21.8325...
    19.7489 18.0444 16.6239 15.4219 14.3916 13.4986 12.7173 12.0278 11.4150 10.8667...
    10.3731 9.9266 9.5206 9.1500 8.8103 8.4976 8.2091 7.9419 7.6938 7.4628...
    7.2472 7.0455 6.8564 6.6787 6.5115 6.3539 6.2050 6.0641 5.9306 5.8040...
    5.6837 5.5693 5.4603 5.3564 5.2572 5.1625 5.0718 4.9850 4.9017 4.8219];
rs = [rs, [0.05:.05:.99]];
ks = [ks, 4.7453 2.8656 2.2358 1.9185 1.7261 1.596 1.5014 1.429 1.3712 1.3235 1.2829 ...
    1.2474 1.2156 1.1862 1.1583 1.1306 1.1019 1.0707 1.0365];
K = interp1(rs, ks, rw, 'spline');
F = K*stat;
p = 1-fcdf(F, 1, n-2);
h = p<alpha;
h = rw;
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
%     [h,p] = WatsonWilliamsTest(x,y);
%     ps = [ps; p];
%     rs = [rs; h];
% end