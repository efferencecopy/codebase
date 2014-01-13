function [baseline, thresh, slope, success] = contrastrespfit(contrasts, spikecounts, type)
    
    if (nargin < 3)
        type = 'quadratic';
    end
    options.MaxIter = 100000;
    options.MaxFunEvals = 100000;
    
    guesses(1) = mean(spikecounts(contrasts == min(contrasts)));
    guesses(2) = min(contrasts(spikecounts > guesses(1)));
    %guesses(2) = mean(contrasts);
    guesses(3) = regress(spikecounts, contrasts);
    
    [model, fval, success] = fminsearch(@mlerr, guesses, options);
    baseline = model(1);
    thresh = model(2);
    slope = model(3);
    
    
    %
    %       NESTED SUBFUNCTIONS:
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % likelihood error function
    function lik = mlerr(input)        
         a = input(1);
         b = input(2);
         c = input(3);
          if (a <= 0)
             lik = inf;
             return 
          end
         
        L = logical(contrasts <= b);
        lambda = a;
        lik1 = spikecounts(L).*log(lambda)-lambda-log(factorial(spikecounts(L)));
        if (strcmp(type, 'linear'))
            lambda = c*(contrasts(~L)-b)+a;
        elseif (strcmp(type, 'quadratic'))
            lambda = c*(contrasts(~L)-b).^2+a;
        else
            lambda = 0;
        end
        lik2 = spikecounts(~L).*log(lambda)-lambda-log(factorial(spikecounts(~L)));
        lik = -sum([lik1;lik2]);
    end
end