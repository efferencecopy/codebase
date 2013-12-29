function [fittedparams, success, modErrs] = weibullFit(contrasts, data, errFun, guesses)
%
% Contrasts should be a column vector of contrasts or optionally an N x 2
% matrix of contrasts (in the first column) and some boolean grouping
% vector (in the second column).  If contrasts is one column, we fit this
% model:
% y = exp(-contrast/alpha).^beta
%
% If contrasts is two columns we fit this model:
% y = exp(-contrast(:,1)/(alpha+delta*contrast(:,2)).^beta
%
% 'Guesses' is a row vector of initial guesses for the fitted parameters.
% Order of parameters: alpha (threshold), beta (slope), gamma (lapse rate)
% delta (change in threshold due to some condition).  Set guesses to 'nan'
% if you don't want to estimate those parameters (important for gamma).
% 

    options.MaxIter = 100000;
    options.MaxFunEvals = 100000;
    
    %determine guesses for alpha and beta
    if nargin<4
        alphaGuess = .1; 
        betaGuess = 1; %the apparent slope if fit with line
        gammaGuess = 1; %where performace asymptotes
        deltaGuess = 0;
        guesses = [alphaGuess, betaGuess, gammaGuess, deltaGuess];
    elseif nargin == 4
        if length(guesses) < 2;  % only alpha guess supplied
            guesses = [guesses, 1 nan nan];
        elseif length(guesses) < 3;  % alpha and beta guesses supplied
            guesses = [guesses, nan nan];
        elseif length(guesses) < 4; % alpha, beta, and gamma guesses supplied
            guesses = [guesses, nan];
        end
    end
    if (size(contrasts,2) == 1)  % No need to estimate delta if no second column in "contrasts"
        guesses(4) = nan;
    end
    if (size(contrasts,2) == 2 & isnan(guesses(4)))  % Do estimate delta if there is second column in "contrasts"
        guesses(4) = 0;
    end
    
    %clean up the inputs. We won't consider data that comes from zero
    %contrast trials, so chuck these data.
    if numel(data) == length(data)
        data = data(:); %make SSE data a column so that the following deletions work for SSE & MLE
    end
    l_zeros = softEq(0, contrasts(:,1), 5); %anything less than 0.00001 gets considered a zero
    contrasts(l_zeros,:) = [];
    data(l_zeros,:) = [];
    
    %call the appropriate fitting routine and perform the fit
    if strcmpi(errFun, 'sse')
        data = data(:);
        contrasts = contrasts(:);
        [model, fval, success] = fminsearch(@sserr, guesses, options);
    elseif strcmpi(errFun, 'mle')
        [model, fval, success] = fminsearch(@mlerr, guesses, options);
        if nargout > 2 %if the SE of the model estimates are desired
            [hess, hessErr] = hessian(@mlerr, model);
            hess = hess(~isnan(model), ~isnan(model));
            invHess = inv(hess);
            modErrs = sqrt(abs(diag(invHess)));
        end
    end
    
    
    %determine the outputs of the function
    fittedparams = model;
    
    
    
    
    
    %
    %       NESTED SUBFUNCTIONS:
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    % sum of squares error function
    % doesn't take optional 4th argument
    function sse = sserr(input)
        
        a = input(1);
        b = input(2);
        if (isnan(input(3)))
            g = 1;
        else
            g = input(3);
        end
        if (isnan(input(4)))
            d = 0;
        else
            d = input(4);
        end
        
        if (size(contrasts,2) == 1)
            fit = g +(0.5-g).*exp(-((contrasts./a).^b));
        elseif (size(contrasts,2) == 2)
            fit = g +(0.5-g).*exp(-((contrasts(:,1)./(a+d*contrasts(:,2))).^b));    
        else
            error('''contrasts'' is peculiarly sized');
        end
        sse = sum((fit-data).^2);
        if (g > 1) sse = Inf; end
        if (a < 0) sse = Inf; end
    end

    % likelyhood error function
    function lik = mlerr(input)   
        a = input(1);
        b = input(2);
        if (isnan(input(3)))
            g = 1;
        else
            g = input(3);
        end
        if (isnan(input(4)))
            d = 0;
        else
            d = input(4);
        end
        
        if (size(contrasts,2) == 1)
            fit = g +(0.5-g).*exp(-((contrasts./a).^b));
        elseif (size(contrasts,2) == 2)
            fit = g +(0.5-g).*exp(-((contrasts(:,1)./(a+d*contrasts(:,2))).^b));    
        else
            error('''contrasts'' is peculiarly sized');
        end
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
       lik = sum(sum([data ~data].*logP));
       lik = -lik;
    end
end