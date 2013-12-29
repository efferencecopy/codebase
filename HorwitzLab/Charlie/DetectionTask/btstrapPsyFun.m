function errs = btstrapPsyFun(alpha, beta, gamma, contrasts, nTrials, nStraps, timeout)

%
%   EXAMPLE errs = btstrapPsyFun(alpha, beta, gamma, contrasts, nTrials, nStraps, [timeout])
%
%


if ~exist('timeout', 'var')
    timeout = inf;
end


a = 1;
tic %start the timer
p = gamma - (gamma-0.5).*exp(-(contrasts./alpha).^beta);
while (a <= nStraps);
    btStrap = binornd(nTrials, p);
    percentCorrect = btStrap ./ nTrials;
    deviates = abs(percentCorrect - 0.82);
    initAlphaIdx = find(deviates == min(deviates), 1, 'last');
    initAlpha = contrasts(initAlphaIdx);
    [aSSE, ~, ~, success(1)] = weibullFit(contrasts, percentCorrect, 'sse' ,[initAlpha 1]);
    [errs.alpha(a), errs.beta(a), errs.gamma(a), success(2)] = weibullFit(contrasts, [btStrap', (nTrials-btStrap)'], 'mle', [aSSE, 1]);
    
    if (sum(success) == 2)
        a = a+1; %only include successfull fits
    end
    if toc > timeout
        errs = nan;
        break
    end
end


%  NOTE AND OBSERVATIONS
%
% establishing a best guess prior to SSE fitting is important. I've just
% used the contrasts closest to 82% correct. It doesn't seem to matter if
% the psychometric fxn is non-monotonic I just have to get the order of
% magnitude correct on the initial alpha guess. 
%
% setting the beta estimate for SSE to 1 seems to minimize the number of
% times fminsearch converges on a hugely to large beta estimate. 
%
% Occasionally SSE returns a plausible alpha estimate but a slope estimate
% that's way to high. Feeding the MLE algorithm the best guess of threshold
% (from SSE) but DISREGARDING the beta estimate (bSSE) from SSE seems to
% help a lot. In these situations MLE is able to converge on plausible
% estimates for both threshold and slope.
%
% for an N of 1, it didn't seem to make a difference if I seeded the MLE
% with the aSSE estimate. If I just seed MLE with a reasonable guess things
% comeout quite similarly.





