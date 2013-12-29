function [out] = mkbasis(in,opt)
%
%     Called by EeroStuff.m  This function will create a 
% basis set for use in visualizing and analysizing the 
% response of the neuron as a function of projections 
% onto 2-D linear subspaces.  
%
%     INPUT
%     in: A column matrix of vectors that span the subspace
% of interest.
%     opt: Option flag that specifies operations to be 
% performed on the input vectors
%          'RotOrtho': Rotate columns [2:end] so that 
%     they are orthogonal to column 1
%          'ProjOrtho': If called with a three column
%     input matrix - returns the linear combination 
%     of the 2nd and 3rd input columns that is orthogonal 
%     to the first input column.
% 
%          If the 'opt' field is left blank, we just scale
%     the columns of 'in' to unit norm.
%
ndim = size(in,1);
nbasisvects = size(in,2);
if (nargin < 2)
   opt = 'Nothing';
end

% First, rescalling vectors to unit length
normfactors = sqrt(sum(in.^2));
in = in./repmat(normfactors, ndim,1);
in(:,normfactors == 0) = 0;

if (strcmp(opt,'ProjOrtho'))
   if (nbasisvects == 3)
      proj = [in(:,1)'*in(:,2); in(:,1)'*in(:,3)];
      proj = rotateXY(proj'./norm(proj),pi/2);
      in = proj(1)*in(:,2)+proj(2)*in(:,3);
   else
      error('Only three columns allowed for ProjOrtho');
   end
end

if (strcmp(opt,'RotOrtho'))
   projs = in(:,1)'*in(:,[2:end]);
   in(:,[2:end]) = in(:,[2:end])-repmat(projs,ndim,1).*...
                   repmat(in(:,1),1,nbasisvects-1);
end

% Rescalling vectors to unit length again
normfactors = sqrt(sum(in.^2));
out = in./repmat(normfactors, ndim,1);
out(:,normfactors == 0) = 0;

