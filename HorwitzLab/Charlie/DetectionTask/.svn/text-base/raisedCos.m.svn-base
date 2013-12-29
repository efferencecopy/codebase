function [prefTheta, expnt, gain] = raisedCos(thetas, resp)

options.MaxIter = 1e6;
options.MaxFunEvals = 1e6;
[gainGuess, maxIdx] = max(resp);
prefThetaGuess = thetas(maxIdx);
expntGuess = 2;


[out, ~, flag] = fminsearch(@cosSSE, [prefThetaGuess, expntGuess, gainGuess], options);
if flag>0
    prefTheta = out(1);
    expnt = out(2);
    gain = out(3);
else
    [prefTheta, expnt, gain] = deal(nan);
end



    function SSE = cosSSE(input)
        pref = input(1);
        ex = input(2);
        g = input(3);
        
        testAngles = thetas-pref;
        modResp = abs((cos(testAngles)).^ex);  %exponentiate the absolute val of the sin(theta)
        modResp = g.*modResp;
        err = resp(:)-modResp(:);
        
        if ex<0;
            SSE = inf;
        else
            SSE = sum(err.^2);
        end
    end
end

% ********************************** %
% testing code for raisedCos:
% cut/copy/paste the guts of this code
% to command window...
% *********************************** %
function testRaisedCos
    fin
    nparams = 10; %the number of different "true" exponents and thetas to test
    nColors = 4;  % mimics the number of different color dirs that GT tests
    trueTheta = linspace(0, pi, nparams);
    trueTheta = repmat(trueTheta,nparams,1);
    trueExp = linspace(0.1, 20, nparams);
    trueExp = repmat(trueExp, nparams, 1)'; %the transpose is critical to do all pairwise comparisons
    thetaEstimate = nan(size(trueTheta));
    expEstimate = nan(size(trueTheta));
    
    for t = 1:nparams;
        for e = 1:nparams;
            colorDirs = linspace(0,pi-(pi/nColors), nColors);
            synthR = 1 * abs(cos(colorDirs-trueTheta(t,e)).^trueExp(t,e));
            [thetaEstimate(t,e), expEstimate(t,e)] = raisedCos(colorDirs, synthR);
            %if expEstimate(t,e) < (trueExp(t,e)/2);disp('exp wrong'); keyboard; end
        end
    end
    figure
    plot(trueTheta', thetaEstimate', '-o');
    figure
    plot(trueExp, expEstimate, '-o');

end