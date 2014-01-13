function out = modelfun(x, lms, fr, r, predr2, opt, err)
x = x(:); % force to column vector
switch opt
    case 'abs'
        pred = abs(lms*x);
    case 'exponent'
        pred = abs(lms*x(1:3)).^x(4);
    case 'baseline'
        pred = abs(lms*x(1:3)) + x(4);
    case 'both'
        pred = abs(lms*x(1:3)).^x(4) + x(5);
    case 'linplusquad'
        pred = x(4)*(abs(lms*x(1:3))) + x(5)*(lms*x(1:3)).^2 + x(6);
    case 'newguy'
        pred = x(1)*(r./sqrt(predr2)) + x(2)*(r./sqrt(predr2)) + x(3);
end
switch err
    case 'fish'
        out = -sum(-pred + fr.*log(pred));
    case 'ss'
        out = sum((fr - pred).^2);
end
if isnan(out) || isinf(out) || ~isreal(out)
    out = 1e10;
end