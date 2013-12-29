%% making myself a tutorial about confidence intervals
fin

%set up the true value of the population
mu = 10;
sigma2 = 10;
sigma = sqrt(sigma2);
N = 30;

%run a simulation to see what the CI is for different draws
iters = 1e2;
CI = nan(iters, 2); %lower and upper values;
for a = 1:iters
    sample = normrnd(mu, sigma, 1, N);
    xbar = mean(sample);
    sd = std(sample);
    sem = sd./sqrt(N);
    z = 1.96 * sem; %assuming 95% CI
    CI(a,:) = [xbar-z, xbar+z];    
end
p = sum((CI(:,1) <= mu) & (CI(:,2) >= mu)) ./ iters

figure, hold on,
plot([1:iters ; 1:iters], CI', 'b-');
plot([1, iters], [mu, mu], 'k--')
title(sprintf('p = %.4f', p))


%% Trying to understand spearman's rho

fin

% make some fake data where the y variable is a monotonic function of x
N = 50;
xx = linspace(0,3, N);
fakeData = [xx(:), xx(:).^1.2] + normrnd(0, 3, N, 2);
[rho_observed, p_observed] = corr(fakeData(:,1), fakeData(:,2), 'type', 'spearman');

% The function that I'll use to caluculate p values
p_from_t = @(n,r) (2.*(1-tcdf(abs(r.*sqrt((n-2)./(1-r.^2))), n-2)));


% plot the two RVs against one another.
figure
set(gcf, 'position', [10,71,1431,385])
subplot(1,3,1)
plot(fakeData(:,1), fakeData(:,2), 'bo', 'markerfacecolor', 'b')
title(sprintf('Spearman''s r = %.3f', rho_observed));
xlabel('RV one')
ylabel('RV two');


% verify that the entry on Wikipedia gives a rho that's consistent with
% MatLab
t_pvalue = p_from_t(N, rho_observed)
p_observed

% now trying to understand power analysis. First, trying to find a value of
% rho that would be significant for a given N.
rho_synth = linspace(0, 0.5, 500); %only considering positive rho
p_synth = p_from_t(N, rho_synth);

subplot(1,3,2), hold on
plot(rho_synth, p_synth, 'k.-')
plot(rho_synth, ones(size(rho_synth)).*0.05, 'r--')
xlabel('synthetic rho')
ylabel('synthetic p with actual N')
title('How big does rho need to be with actual N?')


% now trying to find an N that would give a p<0.05 for the observed value
% of rho
N_synth = 1:500;
p_synth = p_from_t(N_synth, rho_observed);

subplot(1,3,3), hold on,
plot(N_synth, p_synth, 'k.-')
plot(N_synth, ones(size(N_synth))*0.05, 'r--')
xlabel('synthetic N')
ylabel('synthetic p with actual rho')
title('how big does N need to be with actual rho?')





%% a module to find the critical values for rho and N that yield significance

% The function that I'll use to caluculate p values
fin
p_from_t = @(n,r) (2.*(1-tcdf(abs(r.*sqrt((n-2)./(1-r.^2))), n-2)));
rho_observed = -0.076;
N_observed = 86;



% **********************************
% now trying to find an N that would
% give a p<0.05 for the observed value
% of rho
% **********************************
N_synth = 1:25000;
p_synth = p_from_t(N_synth, rho_observed);
ind = find(p_synth<0.05, 1, 'first');
fprintf('Ncrit = %d \n\n', N_synth(ind));
subplot(1,2,1), hold on,
plot(N_synth, p_synth, 'k.-')
plot(N_synth, ones(size(N_synth))*0.05, 'r--')
plot(N_synth(ind), p_synth(ind), 'm*')
xlabel('synthetic N')
ylabel('synthetic p with actual rho')
title('how big does N need to be with actual rho?')
hold off




% **********************************
% now trying to find an rho that would
% give a p<0.05 for the observed value
% of N
% **********************************
rho_synth = linspace(0, 0.5, 500);
p_synth = p_from_t(N_observed, rho_synth);
ind = find(p_synth<0.05, 1, 'first');
fprintf('Rho_crit = %.3f \n\n', rho_synth(ind));
subplot(1,2,2), hold on
plot(rho_synth, p_synth, 'k.-')
plot(rho_synth, ones(size(rho_synth)).*0.05, 'r--')
plot(rho_synth(ind), p_synth(ind), 'm*')
xlabel('synthetic rho')
ylabel('synthetic p with actual N')
title('How big does rho need to be with actual N?')




