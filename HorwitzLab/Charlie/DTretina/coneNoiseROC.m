function [idlob, cones] = coneNoiseROC(params, idlob, cones, gab)

nContrasts = size(idlob.resp,2);
nColors = size(idlob.resp,1);
idlob.roc = {};
for clr = 1:nColors;
    idlob.roc{clr, 1} = nan(1,nContrasts);
    for cnt = 1:nContrasts
                
        % first, based on the repeated trials (if present)
        if size(idlob.resp,3) > 0;
            switch params.obsMethod
                case {'obsMethod_noClrEqSpace', 'obsMethod_phaseInvariant', 'obsMethod_filteredWtFxn'}

                    
                    % this is a hack to make the 6 wt fxn analysis work
                    % with the existing code from the 3 wt fxn analysis.
                    dims = size(idlob.resp{clr,cnt,1});
                    catdim = find(dims == 1);
                    needTranspose = catdim == 2;
                    
                    % Do LDA on the 3D clouds of data points
                    noise = cat(catdim, idlob.resp{clr,1,:}); % <nTrials x nCones (3or6)>
                    if needTranspose; noise = noise'; end
                    mu_noise = mean(noise,1);
                    cov_noise = cov(noise);
                    
                    sig = cat(catdim, idlob.resp{clr,cnt,:});
                    if needTranspose; sig = sig'; end
                    mu_sig = mean(sig,1);
                    cov_sig = cov(sig);
                    
                    % the "within groups" and "between groups" scatter
                    S_within = (cov_noise) + (cov_sig);
                    S_between = (mu_noise - mu_sig)' * (mu_noise - mu_sig); % the outer product
                    
                    % calculate the discriminant vector
                    [vec, val] = eig(S_within \ S_between); % inv(S_within) * S_between
                    [~, maxEigVal] = max(diag(val));
                    W_mc = vec(:,maxEigVal);
                    if W_mc(1)< 0;
                        W_mc = -W_mc; % standardize the representation of eigenvectors
                    end
                    newSig_mc = sig * W_mc;
                    newNoise_mc = noise * W_mc;
                    
                    
                    auc_mc = roc(newNoise_mc, newSig_mc);
                    if mean(newSig_mc) < mean(newNoise_mc);  % LDA eigenvectors can point in one of two directions. Correct for that here.
                        auc_mc = 1-auc_mc;
                    end
                    idlob.roc{clr}(cnt) = auc_mc; %store the ROC value
                    
                    
                case 'obsMethod_all'
                    % the signal and noise distributions are already
                    % univariate, so no LDA is needed
                    idlob.roc{clr}(cnt) = roc(idlob.resp(clr, 1, :), idlob.resp(clr, cnt,:));
            end
            
        end
        
        % Second, based on the analytic solutions to the mean and variance.
        % This will depend on how the ideal observer worked.
        switch params.obsMethod
            case {'obsMethod_noClrEqSpace', 'obsMethod_phaseInvariant', 'obsMethod_filteredWtFxn'}
                
                % need to do LDA. Start by setting things up. Some
                % simulations disregard the S-cones, so determine if we
                % should pay attention to the S-cones.
                if params.enableScones
                    
                    switch params.obsMethod
                        case {'obsMethod_noClrEqSpace', 'obsMethod_filteredWtFxn'}
                        numcones = 3;
                        coneidx = logical([1 1 1]);
                        
                        case 'obsMethod_phaseInvariant'
                        numcones = 6;
                        coneidx = logical([1 1 1 1 1 1]);
                    end
                    
                elseif ~params.enableScones
                    
                    switch params.obsMethod
                        case {'obsMethod_noClrEqSpace', 'obsMethod_filteredWtFxn'}
                            numcones = 2;
                            coneidx = logical([1 1 0]);
                            
                        case 'obsMethod_phaseInvariant'
                            numcones = 4;
                            coneidx = logical([1 1 0 1 1 0]);
                    end
                    
                    % trow a warning if the analysis ignores the S-cone mosaic
                    if clr == 1; disp('Ignoring the S-cone mosiac'); end
                    
                else
                    error('unclear what to do with the S-cones')
                end
                
                % assemble the data. enforce  a particular dimensionality
                % (row, or column vectors)...
                mu_noise = idlob.analyticMean{clr,1}(coneidx);
                mu_noise = mu_noise(:)'; % needs to be a row vector
                tmp = idlob.analyticVar{clr,1}(coneidx); % turned into a column vector in the line below
                cov_noise = eye(numcones) .* repmat(tmp(:),1,numcones);
                
                mu_sig = idlob.analyticMean{clr,cnt}(coneidx);
                mu_sig = mu_sig(:)'; % needs to be a row vector
                tmp = idlob.analyticVar{clr,cnt}(coneidx); % turned into a column vector in the line below
                cov_sig = eye(numcones) .* repmat(tmp(:),1,numcones);
                
                % the "scatter" matricies
                S_within = (cov_noise) + (cov_sig);
                S_between = (mu_noise - mu_sig)' * (mu_noise - mu_sig); % the outer product
                
                % finding the discriminant vector
                [vec, val] = eig(S_within \ S_between); % inv(S_within) * S_between
                [~, maxEigVal] = max(diag(val));
                W_anly = vec(:,maxEigVal); % not assuming the first column corresponds to the largest eigenvalue
                if W_anly(1)< 0;
                    W_anly = -W_anly; % standardize the representation of eigenvectors
                end
                
                % calculate the necessary elements for ROC analysis
                mu_signal = mu_sig * W_anly;
                sigma_signal = sqrt(diag(cov_sig)' * W_anly.^2);
                mu_noise = mu_noise * W_anly;
                sigma_noise = sqrt(diag(cov_noise)' * W_anly.^2);
                
            case 'obsMethod_all'
                
                % no calculations necessary, just pull out the appropriate
                % items from the data
                mu_noise = idlob.analyticMean(clr, 1);
                sigma_noise = sqrt(idlob.analyticVar(clr, 1));
                mu_signal = idlob.analyticMean(clr, cnt);
                sigma_signal = sqrt(idlob.analyticVar(clr, cnt));
                
            case 'obsMethod_absThresh'
                % no calculations necessary, just pull out the appropriate
                % items from the data
                mu_noise = idlob.analyticMean(clr, 1);
                sigma_noise = sqrt(idlob.analyticVar(clr, 1));
                mu_signal = idlob.analyticMean(clr, cnt);
                sigma_signal = sqrt(idlob.analyticVar(clr, cnt));
        end
        
        
        % determine the edges of each distribution
        noise_icdf = norminv([0.001, 0.999], mu_noise, sigma_noise); % the upper and lower portions of the noise distribution
        signal_icdf = norminv([0.001, 0.999], mu_signal, sigma_signal); % the upper and lower portions of the signal distribution
        lowVal = min([noise_icdf(1), signal_icdf(1)]);
        highVal = max([noise_icdf(2), signal_icdf(2)]);
        xx = linspace(lowVal, highVal, 500); % a domain that covers both distributions
        xx = [-inf, xx, inf]; %manually add the p=0 and p=1 condition
        
        % calculate the analytic solution to the ROC
        pFA = 1-normcdf(xx, mu_noise, sigma_noise);
        pHit = 1 - normcdf(xx, mu_signal, sigma_signal);
        auc_anly = -trapz(pFA, pHit);
        if mu_signal < mu_noise; % LDA eigenvectors can point in one of two directions. Correct for that here.
            auc_anly = 1-auc_anly;
        end
        idlob.roc_analytic{clr}(cnt) = auc_anly;
    end
    
end


% fit the ROC data with a cumulative Weibull
[cones.alpha, cones.beta, cones.alpha_analytic, cones.beta_analytic] = deal(nan(nColors,1));
for clr = 1:nColors;
    
    % first, for the repeated trials.
    if size(idlob.resp,3)>0
        tmp_norms = gab.contrasts{clr};
        tmp_area = idlob.roc{clr};
        if ~any(tmp_area>0.90) || ~any(tmp_area<0.55)
            [cones.alpha_analytic(clr), cones.beta_analytic(clr)] = deal(NaN);
        else
            [~, idx] = min(abs(tmp_area - 0.816));
            alphaGuess = tmp_norms(idx);
            [aSSE, bSSE] = weibullFit(tmp_norms, tmp_area, 'sse', [alphaGuess 1]);
            trialsByContrast = ones(1, nContrasts) .* gab.nTrials;
            correctByContrast = (tmp_area.*trialsByContrast);
            wrongByContrast = (trialsByContrast - correctByContrast);
            [cones.alpha(clr), cones.beta(clr)] = weibullFit(tmp_norms, [correctByContrast(:), wrongByContrast(:)], 'mle', [aSSE, bSSE]);
        end
    end
    
    % second, for the analytic solutions
    tmp_norms = gab.contrasts{clr};
    tmp_area = idlob.roc_analytic{clr};
    if ~any(tmp_area>0.90) || ~any(tmp_area<0.55)
        [cones.alpha_analytic(clr), cones.beta_analytic(clr)] = deal(NaN);
    else
        [~, idx] = min(abs(tmp_area - 0.816));
        alphaGuess = tmp_norms(idx);
        [cones.alpha_analytic(clr), cones.beta_analytic(clr)] = weibullFit(tmp_norms, tmp_area, 'sse', [alphaGuess 1]);
    end
end


