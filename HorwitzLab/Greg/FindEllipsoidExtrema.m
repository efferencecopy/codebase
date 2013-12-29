function [bounds] = FindEllipsoidExtrema(M,crit)

% function [bounds] = FindEllipsoidExtrema(M,crit)
%
% This function takes a covariance matrix, M, and a critical
% value (the number of standard deviations to go out in every
% direction), crit, and finds the extreme points on the ellipsoid
% in all three cardinal directions.  Because we're only considering
% the covariance, the mean of the ellipsoid is constrained to be
% the origin of the space.  Thus the output argument, bounds,
% can be interpreted as the maximum of the centered ellipsoid and
% -bounds can be interpreted as the minimum.
%
% Solving this problem required doing Lagrange multipliers.
% let inv(M) = A D E
%              D B F
%              E F C
%
% Then the ellipsoid can be written:
%    Ax^2 + By^2 + Cz^2 + 2Dxy + 2Exz + 2Fyz = crit
%
% And using Lagrange multipliers we can get the following two equations
% expressing y and z as a function of x.
%
% y = (EF-CD)*x/(BC-F^2)
% z = (EF-CD)*x/(BC-F^2)
%
% GDLH 5/17/04

invM = inv(M);

orders = [1 5 9 2 3 6; 5 1 9 2 6 3; 9 5 1 6 3 2];
bounds = [];
for i = 1:3      
   a = invM(orders(i,1));
   b = invM(orders(i,2));
   c = invM(orders(i,3));
   d = 2*invM(orders(i,4));
   e = 2*invM(orders(i,5));
   f = 2*invM(orders(i,6));

   numerator = 4*b*c-f^2;
   denom = 4*a*b*c - c*d^2 - b*e^2 - a*f^2 + d*e*f;
   bounds(i) = sqrt(crit*numerator/denom);
end