function [prefTheta, expnt, gain] = raisedSin(thetas, resp)

options.MaxIter = 100000;
options.MaxFunEvals = 100000;
[gain, maxIdx] = max(resp);
prefThetaGuess = thetas(maxIdx);
expntGuess = 2;


[out, ~, flag, verbage] = fminsearch(@sinErr, [prefThetaGuess, expntGuess], options);
prefTheta = out(1);
expnt = out(2);
if ~flag
    keyboard
end



    function SSE = sinErr(input)
        pref = input(1);
        ex = input(2);
        
        testAngles = thetas+pref;
        modResp = (abs(sin(testAngles))).^ex;  %exponentiate the absolute val of the sin(theta)
        modResp = gain.*modResp;
        err = resp-modResp;
        
        if ex<0;
            SSE = inf;
        else
            SSE = sum(err.^2);
        end
    end
end