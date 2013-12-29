function t_out = coleThresh(mech, beta, colors)
%
%   t_out = coleThresh(mech, beta, colors)
%
% colors:  an Nx3 matrix of LMS cc colors. Each row is a color
% beta:    an Mx1 vector of betas
% mech:    a  3x3 matrix of mechanism directions. Each COLUMN is a mechanism

% make sure the colors are unit vectors
colors = bsxfun(@rdivide, colors, sqrt(sum(colors.^2, 2)));

% calculate the threshold according to the cole 1993 model
t_out = (1./sum(abs(colors * mech).^beta, 2)).^(1./beta);
