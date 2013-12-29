function [N,X] = hist2(y,x)
%  HIST2   2-D histogram
%
%      This function will be something like the Matlab function,
% HIST, but will do things in 2-D (or maybe someday N-D).
%
%      [N,X] = HIST2(Y,M), where M is a scalar, uses M bins per dimension.
%
%      [N,X] = HIST2(Y,X), where X is a 2-vector, uses X(1) bins for 
% dimension 1 and X(2) bins for dimension 2.
%
%      [N,X] = HIST2(Y,X), where X is a 2 x 2 matrix, uses X(1,1) bins for 
% dimension 1 and X(1,2) bins for dimension 2 and ensures that one of the
% bins has its center at [X(2,1) X(2,2)].
%
%      [N,X] = HIST2(Y,X), where X is an 3 x 2 matrix, uses X(1,1) bins for 
% dimension 1 and X(1,2) bins for dimension 2.  The center of the lowest
% bin is X(2,:) and the center of the highest bin in X(3,:).
%
%      By default this function positions bins with centers at the 
% maximum and minimum of both dimensions.  This biases the counts in these
% bins downward and imposes the contraint of at least one count on each
% edge.
%
%      OUTPUTS:
%          N is a matrix of counts
%          X is a cell array containing the bin centers

if (size(y,2) ~= 2)
   error('Y argument must be N x 2');
end

t0 = [];
miny = [];
maxy = [];
if (nargin == 1)
   nbins = [10 10];
elseif (nargin == 2)
   if (size(x,1) == 2 & size(x,2) == 2)
      nbins = x(1,:);
      t0 = x(2,:);
   elseif (size(x,1) == 3 & size(x,2) == 2)
      nbins = x(1,:);
      miny = x(2,:);
      maxy = x(3,:);
   elseif (length(x) == 2)
      nbins = x;
   else  % if (length(x) == 1)
      nbins = [x x];
   end
end

if (isempty(miny) & isempty(maxy))
   miny = min(y);
   maxy = max(y);
end
binwidth = (maxy-miny)./(nbins-1);
leftedges = miny - binwidth/2;
rightedges = maxy + binwidth/2;
L = y(:,1)>rightedges(1) |...
    y(:,1)<leftedges(1) |...
    y(:,2)>rightedges(2) |...
    y(:,2)<leftedges(2);
y(L,:) = [];

n = size(y,1);

if (n > 0)
   scaleddata = y-repmat(leftedges,n,1);
   scalefactors = nbins./(rightedges-leftedges);  % to get to [0:nbins-1]
   scaleddata = scaleddata.*repmat(scalefactors,n,1);
   % scaleddata maps leftedges to 0 and rightedges to nbins-1
   if (~isempty(t0))
      % Applying the same transformation to t0
      scaledt0 = t0-leftedges;
      scaledt0 = scaledt0.*scalefactors;
      % shifting scaleddata so that scaledt0 is an integer vector
      % (so it's the center of a bin)
      additionalshift = round(scaledt0)-scaledt0;
      scaleddata = scaleddata + repmat(additionalshift,n,1);
   else
      additionalshift = [0 0];
   end
%      indices now goes from 0 to nbins-1 
%      (+ e if we've shifted to keep t0 as a bin center) 
   indices = ceil(scaleddata);
%      min(indices) goes to 1 and max(indices) goes to nbins
   subs = sub2ind(nbins,indices(:,1),indices(:,2));
   [sortedsubs, i]=sort(subs);
   s=[1; diff(sortedsubs)];
   uniques = subs(i(s~=0),:);
   nreps=diff(find([s ;1]~=0));

   N = zeros(nbins); 
   for i = 1:max(nreps)
      N(uniques(nreps == i)) = i;
   end

% At this point, the first dimension is on the rows of N
% (represented as variations along the y-dimension) with 
% small numbers (low numbered rows) at the top of the matrix.
% The second dimension is the columns or x-dimension with
% low numbers (low numbered columns) at the right edge.
%
% Before we plot anything or hand anything back, we have to transpose
% N to get the first dimension on the x-axis and the second
% on the y-axis and then flipud to get the y-dimension set up
% so that the small numbers are at the bottom and the large numbers
% are at the top.

   N = flipud(N');

   X{1} = linspace(miny(1),maxy(1),nbins(1))-(additionalshift(1)/scalefactors(1));
   X{2} = linspace(miny(2),maxy(2),nbins(2))-(additionalshift(2)/scalefactors(2));

else   % There is no data in the appropriate range
   N = zeros(nbins); X = [];
end

if (nargout == 0)
   figure;
   imagesc(N);
   axis('square');
   colormap(gray);
   set(gca,'XTick',[1,nbins(1)]);
   set(gca,'YTick',[1,nbins(2)]);
   set(gca,'XTickLabel',[X{1}(1) X{1}(end)]);
   set(gca,'YTickLabel',[X{2}(end) X{2}(1)]);
   clear N;
end

