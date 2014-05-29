function [Kc, SIGMAc, Ks, SIGMAs, fittype, exitflag] = DoGfit(R0, SF, rates)
    
    

    
    % compiling the mean responses to assist in finding good initial
    % guesses
    uniqueSF = unique(SF);
    for a = 1:numel(uniqueSF)
        l = SF == uniqueSF(a);
        meanRates(a) = mean(rates(l));
    end
    
    % amplitude guesses
    Kc_guess = max(meanRates);
    Ks_guess = Kc_guess - meanRates(1);
    
    
    %sigma guesses
    if all(meanRates(end)>meanRates(1:end-1)) % high pass
        fittype = 1;
        sigma_c_guess = 100;
        sigma_s_guess = sigma_c_guess/10;
    elseif all(meanRates(1)>meanRates(2:end)) % low pass?
        fittype = 2;
        Ks_guess = 0;
        Kc_guess = Kc_guess.*10;        
        sigma_s_guess = 0;
        sigma_c_guess = 20;
    else
        fittype = 3;
        sigma_c_guess = 6;
        sigma_s_guess = 1;
    end
    
    
    options.MaxIter = 100e3;
    options.MaxFunEvals = 100e3;
    options.Algorithm = 'interior-point';
    
    problem.objective = @MSE;
    problem.x0 = [Kc_guess, sigma_c_guess, Ks_guess, sigma_s_guess];
    problem.solver = 'fminsearch';
    problem.options = options;
    [out, ~, exitflag] = fminsearch(problem);
    
    
    
    Kc = out(1);
    SIGMAc = out(2);
    Ks = out(3);
    SIGMAs = out(4);
    exitflag = exitflag == 1;
    
    
    function err = MSE(in)
        k_c = in(1);
        mu_c = 0;
        sigma_c = in(2);
        k_s = in(3);
        mu_s = 0;
        sigma_s = in(4);
        
        Rpred = R0 + (k_c .* exp(-((SF-mu_c)./(2.*sigma_c)).^2)) - (k_s .* exp(-((SF-mu_s)./(2.*sigma_s)).^2));
        err = mean((rates-Rpred).^2);
        
        if any([k_c, k_s, sigma_c, sigma_s]<0)
            err = inf;
        end
        
        % set a stricter criterion for very weakly band pass cells
        %if (fittype == 3) && (meanRates(1)/meanRates(2) > 0.8)
            Rpred_0 = R0 + (k_c .* exp(-((0-mu_c)./(2.*sigma_c)).^2)) - (k_s .* exp(-((0-mu_s)./(2.*sigma_s)).^2));
            if Rpred_0 < 0;
                err = inf;
            end
        %end
        
    end
end