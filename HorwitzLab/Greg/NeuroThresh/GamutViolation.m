function out = GamutViolation(norm, colordir, M, bkgndrgb)

if (size(colordir,1) < size(colordir,2))
    colordir = colordir';
end

bkgndlms = M*bkgndrgb;
rgbs = [M\(bkgndlms.*(1+norm*colordir)); % zack changed this to use '\'
    M\(bkgndlms.*(1-norm*colordir))];
if (any(rgbs > 1) || any(rgbs < 0))
    out = 1;
else
    out = 0;
end
end