function [fit, liklihood, CI] = halfSquareFit(designMtx, response, b0, fitType)

%
% Fits a half-squaring contrast response function (with threshold) to data
% asuming Poisson error. "fit" provides the fitted parameters, "liklihood"
% is the negative of the liklihood given the fitted parameters. "CI" is the
% 95% confidence for fit(4) (i.e., the scaleFactor) when the fitType is
% 'yoked'. There are no CIs for 'sharedBaseline' or 'simple'.
%
% One of three different models is fit as specified by "fitType"
% 
%
% fitType == 'simple'
%
%   fit = b(1)+b(3)*(max(x-b(2),0).^2);
%
%   "designMtx" should be a vector of contrasts
%   "response" should be a vector of spike counts
%   "b0" is a vector of three elements
%
%
% fitType == 'yoked'
%
%	for I = 0: pred = b(1)+b(3)*(max(x-b(2),0).^2);
%	for I = 1: pred = b(1)+b(3)*(max(x.*b(4)-b(2),0).^2);
%
%   "designMtx" should be an nTrialsx2 mtx: [contrasts, Indicator]
%   The Indicator should be a 1 or 0, which codes for different color dirs
%   "response" is a vector of spike counts
%   "b0" is a vector of four elements
%
%
% fitType == 'sharedBaseline'
%
%	for I = 0: pred = b(1)+b(3)*(max(x-b(2),0).^2);
%	for I = 1: pred = b(1)+b(5)*(max(x-b(4),0).^2);
%
%   "designMtx" should be an nTrialsx2 mtx: [contrasts, Indicator]
%   The Indicator should be a 1 or 0, which codes for different color dirs
%   "response" is a vector of spike counts
%   "b0" is a vector of five elements
%
% CAH & GDLH 2/2012

% test this function by setting the first input argument to 'testing'
if nargin==1 && strcmpi(designMtx, 'testing')
    test_halfSquareFit()
    return
end


% change the default values for fmincon
options = optimset('Diagnostics','off','Display','off', 'Algorithm', 'active-set', 'TolFun', 1e-15, 'TolX', 1e-15);
vlb = ones(size(b0)) .* 0;
vub = ones(size(b0)) .* Inf;

if strcmpi(fitType, 'simple')
    
    %check the inputs
    if (size(designMtx,2) ~= 1) || (numel(b0)~=3); error('too many predictors for simple fit');end
    
    %specify the upper and lower limits
    x = designMtx(:);
    response = response(:); %make this a column
    
    %constrain the threshold guess to be in the range of the data;
    vlb(2) = min(x);
    vub(2) = max(x);
    
    
    %do the fitting
    [fit, liklihood] = fmincon(@simple,b0,[],[],[],[],vlb,vub,[],options);
    
    %check the output
    if any(fit<vlb) || any(fit>vub); error('Fitted parameters out of bounds'); end %values out of bounds
    
    
elseif strcmpi(fitType, 'yoked')
    
    %check the inputs
    if (size(designMtx,2) ~=2) || (numel(b0)~=4); error('expecting 2 columns for design matrix');end
    
    %identify color directions
    L = logical(designMtx(:,2));
    x = designMtx(:,1);
    response = response(:); %make this a column
    
    %constrain the threshold guess to be in the range of the data;
    vlb(2) = min(x);
    vub(2) = max(x);
    
    %do the fitting
    [fit, liklihood, flag] = fmincon(@yoked,b0,[],[],[],[],vlb,vub,[],options);
    
    %check the output
    if any(fit<vlb) || any(fit>vub); error('Fitted parameters out of bounds'); end %values out of bounds
    
    %compute the confidence interval
    [hess, ~] = hessian(@yoked, fit);
    if rcond(hess) < (2*eps) %i.e. hess is close to singular
        CI = [nan, nan];
    else
        invHess = inv(hess);
        SEM = sqrt(abs(diag(invHess)));
        bounds = [-1.96, 1.96] .* SEM(4);
        CI = fit(4) + bounds;
    end
    
elseif strcmpi(fitType, 'sharedBaseline')
    
    %check the inputs
    if (size(designMtx,2) ~= 2) || (numel(b0) ~= 5); error('size of inputs do not match model');end
    
    %identify color directions
    L = logical(designMtx(:,2));
    x = designMtx(:,1);
    response = response(:); %make this a column
    
    %constrain the threshold guess to be in the range of the data;
    vlb([2,4]) = min(x);
    vub([2,4]) = max(x);
    
    %do the fitting
    [fit, liklihood] = fmincon(@sharedBaseline,b0,[],[],[],[],vlb,vub,[],options);
    
    %check the output
    if any(fit<vlb) || any(fit>vub); error('Fitted parameters out of bounds'); end %values out of bounds
    
end

%check the fits to make sure they conform to the upper and lower bounds
if any(fit<vlb) || any(fit>vub)
    error(['The fitted parameters [', num2str(fit), '] are outside the allowable range'])
end



%
% NESTED SUBFUNCTIONS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function lik = simple(b)
        if any(isnan(b))
            lik = Inf;
        else
            prediction = (b(1)+b(3)*(max(x-b(2),0).^2)) + eps;
            lik = sum(prediction)-(log(prediction')*response);  % -1 * log-likelihood (Poisson)
            
            %Filter these out
            if ~isreal(lik); lik = Inf;end %liklihood is complex when prediction<0.
        end
        
    end

    function lik = yoked(b)
        if any(isnan(b))
            lik = Inf;
        else
            pred1 = (b(1)+b(3)*(max(x(~L)-b(2),0).^2)) + eps;
            pred2 = (b(1)+b(3)*(max(x(L).*b(4)-b(2),0).^2)) + eps;
            lik1 = sum(pred1)-(log(pred1')*response(~L));
            lik2 = sum(pred2)-(log(pred2')*response(L));
            lik = lik1+lik2;
            
            %Filter these out
            if ~isreal(lik); lik = Inf;end %liklihood is complex when prediction<0.
        end
    end

    function lik = sharedBaseline(b)
        if any(isnan(b))
            lik = Inf;
        else
            pred1 = (b(1)+b(3)*(max(x(~L)-b(2),0).^2)) + eps;
            pred2 = (b(1)+b(5)*(max(x(L)-b(4),0).^2)) + eps;
            lik1 = sum(pred1)-(log(pred1')*response(~L));
            lik2 = sum(pred2)-(log(pred2')*response(L));
            lik = lik1+lik2;
            
            %Filter these out
            if ~isreal(lik); lik = Inf;end %liklihood is complex when prediction<0. 
        end
    end
end

%
% TESTING ROUTINES
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function test_halfSquareFit()
    test_pvalues();
    test_confidenceIntervals()
end

function test_pvalues()
    MAKEPLOT = 1;
    nIters = 20;
    ntrials = 8;
    [devstat, devstat_shared, p, p_shared] = deal(nan(nIters,1));

    %alert the user to what's about to happen
    fprintf('\n\n\n    *************** TESTING **************\n')
    fprintf('  About to test halfSquareFit <%d> times\n', nIters)

    try
        for a = 1:nIters
            %assign the contrast scaling
            contrast_scaling = unifrnd(1,10);  % bounded between .1 and 10

            %fake data for "card case"
            b = [unidrnd(5) unifrnd(0,2) unifrnd(2,10)];
            contrasts = linspace(0,5,6);
            lambda = b(1)+b(3)*(max(contrasts-b(2),0)).^2;
            nspikes = poissrnd(repmat(lambda, ntrials, 1));

            %fake data for "int case"
            b_int = [b(1) b(2) b(3)]; %same as card for now (i.e., assume the null hypothesis)
            contrasts_int = contrasts;
            lambda_int = b_int(1)+b_int(3)*(max(contrasts_int.*contrast_scaling-b_int(2),0).^2);
            nspikes_int = poissrnd(repmat(lambda_int, ntrials, 1));

            % --------------------
            % Fitting them individually
            spkThresh = mean(nspikes(:,1)) + (std(nspikes(:,1))./sqrt(size(nspikes,1)));
            threshGuess = contrasts(max((find(mean(nspikes) > spkThresh, 1, 'first'))-1, 0));
            if isempty(threshGuess); threshGuess = 0; end
            cnt_card = repmat(contrasts,ntrials,1); % making a design matrix
            cnt_card = cnt_card(:);
            spikes_card = nspikes(:);
            l_ols = cnt_card >= threshGuess;
            designMtx = [ones(sum(l_ols),1), cnt_card(l_ols).^2];
            beta1 = designMtx \ spikes_card(l_ols);
            params0 = [mean(nspikes(:,1)) threshGuess beta1(2)];
            [fit_card, lik_card] = halfSquareFit(cnt_card, spikes_card, params0, 'simple');


            spkThresh = mean(nspikes_int(:,1)) + (std(nspikes_int(:,1))./sqrt(size(nspikes_int,1)));
            threshGuess = contrasts_int(max((find(mean(nspikes_int) > spkThresh, 1, 'first'))-1, 0));
            if isempty(threshGuess); threshGuess = 0; end
            cnt_int = repmat(contrasts_int,ntrials,1); % making a design matrix
            cnt_int = cnt_int(:);
            spikes_int = nspikes_int(:);
            l_ols = cnt_int >= threshGuess;
            designMtx = [ones(sum(l_ols),1), cnt_int(l_ols).^2];
            beta2 = designMtx \ spikes_int(l_ols);
            params0 = [mean(nspikes_int(:,1)) threshGuess beta2(2)];
            [fit_int, lik_int] = halfSquareFit(cnt_int, spikes_int, params0, 'simple');


            % Trying to fit both yoked
            yoke_b0 = mean([fit_card(1), fit_int(1)]); %guess for baseline counts
            yoke_b1 = fit_card(2); %guess for spike threshold
            yoke_b2 = fit_card(3); %guess for gain
            yoke_b3 = sqrt(fit_int(3)/fit_card(3)); %guess for contrast scaling
            params1 = [yoke_b0, yoke_b1, yoke_b2, yoke_b3];
            X = [cnt_card ; cnt_int];
            X(:,2) = [zeros(numel(cnt_card),1); ones(numel(cnt_int),1)];
            y = [spikes_card; spikes_int];
            [fit_yoke, lik_yoke] = halfSquareFit(X,y,params1,'yoked');

            if any([fit_yoke, fit_card, fit_int]<0)
                disp('fits exceed limits')
            end


            % Re-fit the CRFs individually using the yoked parameters as
            % estimates.
            params0 = fit_yoke(1:3);
            [fit_card2, lik_card2] = halfSquareFit(cnt_card, spikes_card, params0, 'simple');
            params0 = [fit_yoke(1), fit_yoke(2)./fit_yoke(3), fit_yoke(3).*fit_yoke(4)^2];
            [fit_int2, lik_int2] = halfSquareFit(cnt_int, spikes_int, params0, 'simple');

            %select the best separate fits based off their liklihoods
            if lik_int2 < lik_int
                lik_int = lik_int2;
                fit_int = fit_int2;
            end
            if lik_card2< lik_card
                lik_card = lik_card2;
                fit_card = fit_card2;
            end


            %fit both sets with a model that has a shared baseline
            params0 = [mean([fit_card(1), fit_int(1)]), fit_card(2:3), fit_int(2:3)];
            [fit_sharedBaseline, lik_sharedBaseline] = halfSquareFit(X, y, params0, 'sharedBaseline');

            %compute the deviances and p-values
            devstat(a) = 2*(lik_yoke - (lik_card + lik_int));  % "full model" goes second. dev1 and dev2 are -1*llik
            p(a) = 1- chi2cdf(devstat(a),2);
            devstat_shared(a) = 2*(lik_yoke - lik_sharedBaseline);
            p_shared(a) = 1- chi2cdf(devstat_shared(a),1);

            if MAKEPLOT
                figure
                set(gcf, 'position', [26   281   977   407])
                subplot(1,2,1), cla
                hold on
                plot(contrasts,nspikes','k.');
                plot(contrasts_int,nspikes_int','m.')
                x = linspace(0,max([contrasts, contrasts_int]),100);
                plot(x,fit_card(1)+fit_card(3)*(max(x-fit_card(2),0).^2),'-k');
                plot(x,fit_int(1)+fit_int(3)*(max(x-fit_int(2),0).^2),'-m');
                set(gca,'Ylim',[0 max([nspikes_int(:);nspikes(:)])]);
                hold off
                title('Fit Separately')

                subplot(1,2,2), cla
                hold on,
                plot(contrasts,nspikes','k.');
                plot(contrasts_int,nspikes_int','m.');
                x = linspace(0,max([contrasts, contrasts_int]),30);
                plot(x,fit_yoke(1)+fit_yoke(3)*(max(x-fit_yoke(2),0).^2),'b--'); % fmincon fit
                plot(x,fit_yoke(1)+fit_yoke(3)*(max(x*fit_yoke(4)-fit_yoke(2),0).^2),'b--'); % fmincon fit
                plot(x,fit_sharedBaseline(1)+fit_sharedBaseline(3)*(max(x-fit_sharedBaseline(2),0).^2),'c'); % shared baseline for "card"
                plot(x,fit_sharedBaseline(1)+fit_sharedBaseline(5)*(max(x-fit_sharedBaseline(4),0).^2),'c'); % shared baseline for "card"
                set(gca,'Ylim',[0 max([nspikes(:);nspikes_int(:)])]);
                title('Yoked = Blue, Shared Baseline = Cyan');
                hold off
            end
            if rem(a,50) == 0
                fprintf('iteration %d of %d\n', a, nIters);
            end
        end
    catch
        keyboard
    end


    %
    % for yoked vs two single fits
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %plot the results
    figure
    subplot(1,2,1)
    hist(p,15)
    title('H_{o}: Yoked vs. Full')
    xlabel('p value')

    %plot the distribution of deviances
    subplot(1,2,2), hold on,
    edges = linspace(min(devstat), max(devstat), 40);
    counts = histc(devstat, edges);
    normCounts = counts./sum(counts);
    binSize = edges(2)-edges(1);
    edges = edges+binSize/2;
    chi2pred = chi2pdf(edges, 2) ./ sum(chi2pdf(edges, 2));
    bar(edges, normCounts, 'type', 'histc')
    plot(edges, chi2pred, 'r', 'linewidth', 2)
    title('H_{o}: Yoked vs. Full')
    xlabel('deviance')


    %
    % for yoked vs shared baseline fits
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %plot the results
    figure
    subplot(1,2,1)
    hist(p_shared,15)
    title('H_{o}: Yoked vs. Shared Baseline')
    xlabel('p value')

    %plot the distribution of deviances
    subplot(1,2,2), hold on,
    edges = linspace(min(devstat_shared), max(devstat_shared), 40);
    counts = histc(devstat_shared, edges);
    normCounts = counts./sum(counts);
    binSize = edges(2)-edges(1);
    edges = edges+binSize/2;
    chi2pred = chi2pdf(edges, 1) ./ sum(chi2pdf(edges, 1));
    bar(edges, normCounts, 'type', 'histc')
    plot(edges, chi2pred, 'r', 'linewidth', 2)
    title('H_{o}: Yoked vs. Shared Baseline')
    xlabel('deviance')
end

function test_confidenceIntervals()
    MAKEPLOT = 0;
    nIters = 4;
    nTestsPerIter = 200;
    ntrials = 8;
    
    percentTrueValueInCI = nan(nIters,1);
    for a = 1:nIters
        [estScaleFactor] = nan(nTestsPerIter,1);
        CI = nan(nTestsPerIter,2);
        
        %assign the contrast scaling
        contrast_scaling = unifrnd(1,10);  % bounded between .1 and 10
        b = [unidrnd(5) unifrnd(0.2,2) unifrnd(2,10)];
        b_int = [b(1) b(2) b(3)]; %same as card for now (i.e., assume the null hypothesis)
        
        
        for j = 1:nTestsPerIter
            
            %fake data for "card case"
            contrasts = linspace(0,5,6);
            lambda = b(1)+b(3)*(max(contrasts-b(2),0)).^2;
            nspikes = poissrnd(repmat(lambda, ntrials, 1));
            
            %fake data for "int case"
            contrasts_int = contrasts/contrast_scaling;
            lambda_int = b_int(1)+b_int(3)*(max(contrasts_int.*contrast_scaling-b_int(2),0).^2);
            nspikes_int = poissrnd(repmat(lambda_int, ntrials, 1));
            
            % --------------------
            % Fitting them individually
            spkThresh = mean(nspikes(:,1)) + (std(nspikes(:,1))./sqrt(size(nspikes,1)));
            threshGuess = contrasts(max((find(mean(nspikes) > spkThresh, 1, 'first'))-1, 0));
            if isempty(threshGuess); threshGuess = 0; end
            cnt_card = repmat(contrasts,ntrials,1); % making a design matrix
            cnt_card = cnt_card(:);
            spikes_card = nspikes(:);
            l_ols = cnt_card >= threshGuess;
            designMtx = [ones(sum(l_ols),1), cnt_card(l_ols).^2];
            beta1 = designMtx \ spikes_card(l_ols);
            params0 = [mean(nspikes(:,1)) threshGuess beta1(2)];
            [fit_card, ~] = halfSquareFit(cnt_card, spikes_card, params0, 'simple');
            
            
            spkThresh = mean(nspikes_int(:,1)) + (std(nspikes_int(:,1))./sqrt(size(nspikes_int,1)));
            threshGuess = contrasts_int(max((find(mean(nspikes_int) > spkThresh, 1, 'first'))-1, 0));
            if isempty(threshGuess); threshGuess = 0; end
            cnt_int = repmat(contrasts_int,ntrials,1); % making a design matrix
            cnt_int = cnt_int(:);
            spikes_int = nspikes_int(:);
            l_ols = cnt_int >= threshGuess;
            designMtx = [ones(sum(l_ols),1), cnt_int(l_ols).^2];
            beta2 = designMtx \ spikes_int(l_ols);
            params0 = [mean(nspikes_int(:,1)) threshGuess beta2(2)];
            [fit_int, ~] = halfSquareFit(cnt_int, spikes_int, params0, 'simple');
            
            
            % Trying to fit both yoked
            yoke_b0 = mean([fit_card(1), fit_int(1)]); %guess for baseline counts
            yoke_b1 = fit_card(2); %guess for spike threshold
            yoke_b2 = fit_card(3); %guess for gain
            yoke_b3 = sqrt(fit_int(3)/fit_card(3)); %guess for contrast scaling
            params1 = [yoke_b0, yoke_b1, yoke_b2, yoke_b3];
            X = [cnt_card ; cnt_int];
            X(:,2) = [zeros(numel(cnt_card),1); ones(numel(cnt_int),1)];
            y = [spikes_card; spikes_int];
            [fit_yoke, ~, CI(j,:)] = halfSquareFit(X,y,params1,'yoked');
            estScaleFactor(j) = fit_yoke(4);
            
            if any([fit_yoke, fit_card, fit_int]<0)
                disp('fits exceed limits')
            end
        end %nTestsPerIter
        
        %determine how often the CI includes the real value
        N = sum((contrast_scaling >= CI(:,1))& (contrast_scaling <= CI(:,2)));
        percentTrueValueInCI(a) = N./nTestsPerIter;
        
        %plot if need be
        if MAKEPLOT
            figure, hold on,
            plot(CI', repmat([1:size(CI,1)], 2, 1), 'b')
            plot(estScaleFactor, [1:size(CI,1)], 'b.');
            plot([contrast_scaling, contrast_scaling], [0, size(CI,1)], 'm:')
            ylim([-2, size(CI,1)+2])
            title(sprintf('Percent of CIs that contain true value: %g', percentTrueValueInCI(a)))
            ylabel('Cell Number')
            xlabel('95% CI for Scale Factor')
            hold off
        end
    end %nIters
    
    %report a histogram of %within CI
    figure
    plot([1:nIters], percentTrueValueInCI', 'b-o')
    xlim([0-1, nIters+2])
    ylim([0.85, 1])
end %function