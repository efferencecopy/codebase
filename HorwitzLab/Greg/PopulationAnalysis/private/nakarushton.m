function out = nakarushton(beta, x)
    rmax = beta(1);
    n = beta(2);
    c50 = beta(3);
    
    out = rmax.*(x.^n./(x.^n+c50.^n));
end