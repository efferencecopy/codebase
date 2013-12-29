function [alpha, beta, gamma, success, modErrs] = weibullFit(contrasts, data, errFun, guesses)
    
    options.MaxIter = 100000;
    options.MaxFunEvals = 100000;
    
    %determine guesses for alpha and beta
    useModGamma = 1;
    if nargin<4
        alphaGuess = .1; 
        betaGuess = 1; %the apparent slope if fit with line
        gammaGuess = 1; %where performace asymtopes
        guesses = [alphaGuess, betaGuess, gammaGuess];
    elseif nargin == 4
        if length(guesses) < 3;
            useModGamma = 0;
        end
    end
    
    %clean up the inputs. We won't consider data that comes from zero
    %contrast trials, so chuck these data.
    if numel(data) == length(data)
        data = data(:); %make SSE data a column so that the following deletions work for SSE & MLE
    end
    l_zeros = softEq(0, contrasts, 5); %anything less than 0.00001 gets considered a zero
    contrasts(l_zeros) = [];
    data(l_zeros,:) = [];
    
    
    
    %call the appropriate fitting routine and perform the fit
    if strcmpi(errFun, 'sse')
        data = data(:);
        contrasts = contrasts(:);
        [model, ~, success] = fminsearch(@sserr, guesses, options);
    elseif strcmpi(errFun, 'mle')
       [model, ~, success] = fminsearch(@mlerr, guesses, options);
       if nargout > 4 %if the SE of the model estimates are desired
          [hess, ~] = hessian(@mlerr, model);
          if rcond(hess) > (2*eps) %i.e., the matrix is poorly invertible
              invHess = inv(hess);
              modErrs = sqrt(abs(diag(invHess)));
          else
              modErrs = nan(2,1);
          end
       end
    end
    
    
    %determine the outputs of the function
    alpha = model(1);
    beta = model(2);
    if useModGamma
        gamma = model(3);
    else
        gamma = 1;
    end
    
    
    
    
    %
    %       NESTED SUBFUNCTIONS:
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    % sum of squares error function
    function sse = sserr(input)
        a = input(1);
        b = input(2);
        if (length(input) < 3)
            g = 1;
        else
            g = input(3);
        end
        fit = g +(0.5-g).*exp(-((contrasts./a).^b));
        sse = sum((fit-data).^2);
        if (g > 1); sse = Inf; end
        if (a < 0); sse = Inf; end
    end

    % likelihood error function
    function lik = mlerr(input)        
        a = input(1);
        b = input(2);
        if (length(input) < 3)
            g = 1;
        else
            g = input(3);
        end
        fit = g +(0.5-g).*exp(-((contrasts./a).^b));
        if (any(fit > 1) || g > 1)
           lik = inf; 
           return
        end
        if a < 0;
           lik = inf; 
           return
        end
       fit(fit>1-eps) = 1-eps;  % Ugly hack to avoid log10(0) = -Inf errors.
       logP = log10([fit(:), 1-fit(:)]);
       lik = sum(sum(data.*logP));
       lik = -lik;
    end
end