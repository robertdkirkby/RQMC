%% Monte Carlo Demo: estimating E[U], U~Uniform(0,1)
% The true mean is 1/2. We estimate it with a sample mean of N i.i.d. draws.
% Repeating the experiment M times gives us the sampling distribution of the
% estimator, so we can look at its mean, variance, and histogram.
%
% This is the baseline (plain Monte Carlo). Later we compare against QMC and RQMC.

fprintf('\nMC results:\n')

%% N = 20 draws
% N1 = 20;
estimates_N1 = zeros(M,1);
for m = 1:M
    U = rand(N1,1);            % N i.i.d. Uniform(0,1) draws
    estimates_N1(m) = mean(U); % the estimator for this repetition
end

mean_N1 = mean(estimates_N1); % should be close to the true value 1/2
var_N1  = var(estimates_N1);  % sampling variance of the estimator

fprintf('N=20:  mean of estimator = %.6f, variance of estimator = %.8f\n', mean_N1, var_N1)
fprintf('       (theoretical: mean = 0.5, variance = Var(U)/N = (1/12)/%d = %.8f)\n', N1, (1/12)/N1)

%% N = 100 draws
% N2 = 100;
estimates_N2 = zeros(M,1);
for m = 1:M
    U = rand(N2,1);
    estimates_N2(m) = mean(U);
end

mean_N2 = mean(estimates_N2);
var_N2  = var(estimates_N2);

fprintf('N=100: mean of estimator = %.6f, variance of estimator = %.8f\n', mean_N2, var_N2)
fprintf('       (theoretical: mean = 0.5, variance = Var(U)/N = (1/12)/%d = %.8f)\n', N2, (1/12)/N2)

%% Histograms of the sampling distribution
figure(1)
subplot(2,1,1)
histogram(estimates_N1, 30)
xline(0.5, 'r--', 'LineWidth', 1.5)
title(sprintf('MC estimator of E[U],  N=20  (var = %.2e)', var_N1))
xlabel('estimate'); ylabel('count')
xlim([0 1])

subplot(2,1,2)
histogram(estimates_N2, 30)
xline(0.5, 'r--', 'LineWidth', 1.5)
title(sprintf('MC estimator of E[U],  N=100  (var = %.2e)', var_N2))
xlabel('estimate'); ylabel('count')
xlim([0 1])
