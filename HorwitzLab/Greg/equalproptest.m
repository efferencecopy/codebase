function [h,p] = equalproptest(counts,ns,alpha)
%   This function will perform a test
% to determine whether 'n' proportions
% are the same or different.
%
% function [h,p] = EqualPropTest(counts,ns,alpha)
%
% (See Snedecor and Cochran for details)

if (length(ns) ~= length(counts))
  error('The list of counts and ns have to have the same length');
end

props = counts./ns;

% From Snedecor and Cochran p.203
pbar = sum(counts)/sum(ns);
X2 = (sum(props.*counts) - pbar*sum(counts))/(pbar*(1-pbar));
p = 1-chi2cdf(X2,length(counts)-1);
h = p<alpha;

% Testing
%ngroups = 5;
%n = 40;
% for i = 1:2000
%   counts = binornd(n,.2,ngroups,1);
%   [h,p] = EqualPropTest(counts,repmat(n,ngroups,1),.05);
%   ps = [ps p];
% end
% This is a meaningless change to the code to test out SVN
% And another