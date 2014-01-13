function [c,ceq] = mynlcon(x)
c = [];
ceq = norm(x(1:3)) - 1;