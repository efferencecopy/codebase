function out = circleVar(resp, theta)
%
% EXAMPLE: CIRVAR = circVar(rates, theta);
%
% Computes the circular variance according to equations in Johnson &
% Shapley 2001, and their reference in their 2008 paper. CircVar is just
% 1-R where R is the vector sum of all the data (when represented as
% vectors)
%
% CAH 10/11


R = sum(resp .* exp(sqrt(-1)*2.*theta))./sum(resp); %a complex number
out = 1-abs(R); %from ringach paper cited by Johnson 2008