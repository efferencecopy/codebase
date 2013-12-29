function fpars_bt = bootstrapDetectSurf(colors, alphas, fpar_emp, nstraps)

% assign some values
coleThresh = @(mech, beta, colors) (1./sum(abs(colors * mech).^beta, 2)).^(1./beta);

% find the residuals
colors = bsxfun(@rdivide, colors, sqrt(sum(colors.^2, 2))); % make sure they're unit vecs
mod_thresh_emp = coleThresh(reshape(fpar_emp(2:end),3,3), fpar_emp(1), colors);
log_residual = log(alphas) - log(mod_thresh_emp);
N = numel(alphas);

% reassign the residuals to the model thresholds based on empirical data
% and then re-run the fit. save the fparams as the bootstrap versions.
starttime = GetSecs;
strapNum = 1;
while (strapNum <= nstraps)
    
    idx = unidrnd(N, N, 1); % resampling with replacement! 
    alphas_bt = mod_thresh_emp .* exp(log_residual(idx)); % bootstraped thresholds
    
    % run the fit using the three typical methods
    [fpar_default, fval_default] = fitDetectionSurface(colors, alphas_bt, 'ellipsoid');
    [fpar_initParam, fval_initParam] = fitDetectionSurface(colors, alphas_bt, fpar_emp);
    
    
    % find the best method, store it as a boot straped version.
    fvals_tmp = [fval_default; fval_initParam]
    fpars_tmp = [fpar_default'; fpar_initParam];
    [~, idx] = min(fvals_tmp)
    bestPars = fpars_tmp(idx,:);    
    
    if bestPars(1)<50;
        fpars_bt{strapNum} = bestPars;
        strapNum = strapNum+1;
    end
    
    % alert the user to how much time things are taking.
    if rem(strapNum,100)==0
        clc
        fprintf('bootstrap number: %d in %.3f seconds\n', strapNum, GetSecs-starttime)
    end
end

