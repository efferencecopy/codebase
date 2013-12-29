function [low, hi] = rhoCI(x,y, type);

%accepts bi-variate data and calculates 95% confidence intervals for rho

if ~exist('type', 'var')
    type = 'spearman';
end

if numel(x) ~= numel(y);
    error('X and Y have unequal sizes')
end

N = numel(x);
iters = 100;
fakeX = nan(N, iters);
shuffIndX = unidrnd(N, size(fakeX));
fakeY = nan(N, iters);
shuffIndY = unidrnd(N, size(fakeY));

keyboard