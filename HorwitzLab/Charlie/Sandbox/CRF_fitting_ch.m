
fin
MAKEPLOT = 0;
FITOLD = 0;
nIters = 2000;
ntrials = 8;
options = optimset('Diagnostics','off','Display','off', 'Algorithm', 'active-set', 'TolFun', 1e-15, 'TolX', 1e-15);
[devstat, devstat_shared, p, p_shared] = deal(nan(nIters,1));

for a = 1:nIters
    %assign the contrast scaling
    contrast_scaling = unifrnd(1,10);  % bounded between .1 and 10
    
    %fake data for "card case"
    b = [5 1 5];
    b = [unidrnd(10) unifrnd(0,3) unifrnd(5,20)];
    contrasts = linspace(0,4,6);
    lambda = b(1)+b(3)*(max(contrasts-b(2),0)).^2;
    nspikes = poissrnd(repmat(lambda, ntrials, 1));

    %fake data for "int case"
    b_int = [b(1) b(2) b(3)]; %same as card for now
    contrasts_int = contrasts/contrast_scaling;
    lambda_int = b_int(1)+b_int(3)*(max(contrasts_int.*contrast_scaling-b_int(2),0).^2);
    nspikes_int = poissrnd(repmat(lambda_int, ntrials, 1));
    
    % --------------------
    % Fitting them individually
    
    spkThresh = mean(nspikes(:,1)) + (std(nspikes(:,1))./sqrt(size(nspikes,1)));
    threshGuess = contrasts(max((find(mean(nspikes) > spkThresh, 1, 'first'))-1, 0));
    if isempty(threshGuess); threshGuess = 0; end
    x = repmat(contrasts,ntrials,1); % making a design matrix
    x = x(:);
    y = nspikes(:);
    l_ols = x >= threshGuess;
    beta1 = regress(y(l_ols),[ones(length(x(l_ols)),1) x(l_ols).^2]);
    params0 = [mean(nspikes(:,1)) threshGuess beta1(2)];
    c = repmat(contrasts,ntrials,1);
    [f1_ch, dev1_ch] = halfSquareFit(c(:), nspikes(:), params0, 'simple');
    
    if FITOLD
        vlb = [0 0 0];
        vub = [50 50 50];
        f1 = fmincon('HalfSquareFit_old',params0,[],[],[],[],vlb,vub,[],options,c(:),nspikes(:));
        dev1 = HalfSquareFit_old(f1,c(:),nspikes(:));
    end
    
    spkThresh = mean(nspikes_int(:,1)) + (std(nspikes_int(:,1))./sqrt(size(nspikes_int,1)));
    threshGuess = contrasts_int(max((find(mean(nspikes_int) > spkThresh, 1, 'first'))-1, 0));
    if isempty(threshGuess); threshGuess = 0; end
    x = repmat(contrasts_int,ntrials,1); % making a design matrix
    x = x(:);
    y = nspikes_int(:);
    l_ols = x >= threshGuess;
    beta2 = regress(y(l_ols),[ones(length(x(l_ols)),1) x(l_ols).^2]);
    params0 = [mean(nspikes_int(:,1)) threshGuess beta2(2)];
    c = repmat(contrasts_int,ntrials,1);
    [f2_ch, dev1_ch(2)] = halfSquareFit(c(:), nspikes_int(:), params0, 'simple');
    
    if FITOLD
        f2 = fmincon('HalfSquareFit_old',params0,[],[],[],[],vlb,vub,[],options,c(:),nspikes_int(:));
        dev1(2)  = HalfSquareFit_old(f2,c(:),nspikes_int(:));
        
        difs = [dev1, f1, f2] - [dev1_ch, f1_ch, f2_ch];
        if any(difs > 0)
            fprintf('Difference:\t')
            disp(difs)
        end
    end
    

    % Trying to fit both yoked
    %params1 = [mean([f1_ch(1); f2_ch(1)]) f1_ch(2) f1_ch(3) sqrt(f2_ch(3)./f1_ch(3))];
    params1 = [mean([f1_ch(1); f2_ch(1)]) f1_ch(2) f1_ch(3) contrast_scaling];
    c1 = repmat(contrasts,ntrials,1);
    c2 = repmat(contrasts_int,ntrials,1);
    X = [c1(:);c2(:)];
    X(:,2) = [zeros(numel(c1),1); ones(numel(c2),1)];
    y = [nspikes(:); nspikes_int(:)];
    [f_yoke, dev_yoke_ch] = halfSquareFit(X,y,params1,'yoked');
    
    if any([f_yoke, f1_ch, f2_ch]<0)
        disp('fits exceed limits')
    end
    
    if FITOLD
        vlb = [0 0 0 .1];
        vub = [50 50 50 10];
        f = fmincon('HalfSquareFit_old',params1,[],[],[],[],vlb,vub,[],options,X,y);
        dev2 = HalfSquareFit_old(f,X,y);
    end
    
    devstat(a) = 2*(dev_yoke_ch-sum(dev1_ch));  % "full model" goes second. dev1 and dev2 are -1*llik
    p(a) = 1- chi2cdf(devstat(a),2);
    
    %fit both sets with a model that has a shared baseline
    params0 = [mean([f1_ch(1), f2_ch(1)]), f1_ch(2:3), f2_ch(2:3)]; 
    [fit_sharedBaseline, dev_sharedBaseline] = halfSquareFit(X, y, params0, 'sharedBaseline');
    devstat_shared(a) = 2*(dev_yoke_ch - dev_sharedBaseline);
    p_shared(a) = 1- chi2cdf(devstat_shared(a),1);
    
    if MAKEPLOT && (any([devstat(a), devstat_shared(a)]>50) || any([devstat(a), devstat_shared(a)]<0))
        figure(1)
        subplot(1,2,1), cla
        hold on
        plot(contrasts,nspikes','k.');
        plot(contrasts_int,nspikes_int','m.')
        x = linspace(0,max([contrasts, contrasts_int]),100);
        plot(x,f1_ch(1)+f1_ch(3)*(max(x-f1_ch(2),0).^2),'k');
        plot(x,f2_ch(1)+f2_ch(3)*(max(x-f2_ch(2),0).^2),'m');
        if FITOLD
            plot(x,f1(1)+f1(3)*(max(x-f1(2),0).^2),'k--');
            plot(x,f2(1)+f2(3)*(max(x-f2(2),0).^2),'m--');
        end
        set(gca,'Ylim',[0 max([nspikes_int(:);nspikes(:)])]);
        hold off
        title('Fit Separately')
        
        subplot(1,2,2), cla
        hold on,
        plot(contrasts,nspikes','k.');
        plot(contrasts_int,nspikes_int','m.');
        x = linspace(0,max([contrasts, contrasts_int]),30);
        plot(x,f_yoke(1)+f_yoke(3)*(max(x-f_yoke(2),0).^2),'b.'); % fmincon fit
        plot(x,f_yoke(1)+f_yoke(3)*(max(x*f_yoke(4)-f_yoke(2),0).^2),'b.'); % fmincon fit
        plot(x,fit_sharedBaseline(1)+fit_sharedBaseline(3)*(max(x-fit_sharedBaseline(2),0).^2),'c'); % shared baseline for "card"
        plot(x,fit_sharedBaseline(1)+fit_sharedBaseline(5)*(max(x-fit_sharedBaseline(4),0).^2),'c'); % shared baseline for "card"
        if FITOLD
            plot(x,f(1)+f(3)*(max(x-f(2),0).^2),'k--'); % fmincon fit
            plot(x,f(1)+f(3)*(max(x*f(4)-f(2),0).^2),'m--'); % fmincon fit
            title(['p = ',num2str(p(a))]);
        end
        set(gca,'Ylim',[0 max([nspikes(:);nspikes_int(:)])]);
        title(sprintf('yoked: %.3f shared: %.3f', devstat(a), devstat_shared(a)));
        hold off
        keyboard
    end
    if rem(a,20) == 0
        disp(a)
    end
end


%
% for yoked vs two single fits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot the results
figure
subplot(1,2,1)
hist(p,15)
title('p vaules under the null hypothesis')
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
title('dev stats under null hypothesis')
xlabel('deviance')


%
% for yoked vs shared baseline fits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot the results
figure
subplot(1,2,1)
hist(p_shared,15)
title('p vaules under the null hypothesis')
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
title('dev stats under null hypothesis')
xlabel('deviance')