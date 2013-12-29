function [Kc, SIGMAc, Ks, SIGMAs] = DoGfit(R0, SF, rates)
    
    
    options.MaxIter = 100000;
    options.MaxFunEvals = 100000;
    
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
    if all(meanRates(end)>meanRates(1:end-1))
        sigma_c_guess = 300;
        sigma_s_guess = sigma_c_guess/6;
    elseif all(meanRates(2)>meanRates(2:end))
        sigma_s_guess = 300;
        sigma_c_guess = sigma_s_guess/6;
    else
        Kc_guess = max(meanRates);
        Ks_guess = Ks_guess .* 1.2;
        sigma_c_guess = 1;
        sigma_s_guess = 1;
    end
    
    
    out = fminsearch(@MSE, [Kc_guess, sigma_c_guess, Ks_guess, sigma_s_guess], options);
    Kc = out(1);
    SIGMAc = out(2);
    Ks = out(3);
    SIGMAs = out(4);
    
    
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
    end
end