%% RQMC_CheckNormal: is the RQMC estimator normally distributed?
% Plain Monte Carlo estimators are approximately normal by the central limit
% theorem. RQMC estimators need not be: the low-discrepancy structure that
% shrinks their variance also distorts the shape of the sampling distribution.
% Here we check the Sobol + Owen scramble estimator, using both a formal test
% and a visual comparison against a fitted normal.
%
% Requires est_owen and Nvals from RandomizedQuasiMonteCarlo_Demo.

%% Lilliefors normality test
% (Lilliefors = Kolmogorov-Smirnov with the mean and variance estimated from
% the data, which is the right variant here since we do not know them a priori.)
% A small p-value means we reject the hypothesis that the estimates are normal.
% A warning that "P is less than the smallest tabulated value" means the true
% p is below 0.001, i.e. normality is rejected even more strongly.
fprintf('\nNormality of the Sobol + Owen RQMC estimator (Lilliefors test):\n')
for jN = 1:2
    [h, pL] = lillietest(est_owen(:,jN));
    fprintf('   N=%d: h = %d, p = %.4f  (h=1 => reject normality at 5%%)\n', ...
            Nvals(jN), h, pL)
    if h == 1
        fprintf('   ==> Lilliefors test rejects normal distribution\n')
    else
        fprintf('   ==> Lilliefors test cannot reject normal distribution\n')
    end
end

%% Figure 4: histograms with a fitted normal distribution
% histfit overlays the maximum-likelihood normal fit on the histogram, so we
% can visually judge how close the Owen-scrambled estimator is to normal.
figure(4)
subplot(2,1,1)
histfit(est_owen(:,1), 30, 'normal')
xline(0.5, 'r--', 'LineWidth', 1.5)
title(sprintf('Sobol + Owen scramble, N=20  (var = %.2e)', var(est_owen(:,1))))
xlabel('estimate'); ylabel('count')

subplot(2,1,2)
histfit(est_owen(:,2), 30, 'normal')
xline(0.5, 'r--', 'LineWidth', 1.5)
title(sprintf('Sobol + Owen scramble, N=100  (var = %.2e)', var(est_owen(:,2))))
xlabel('estimate'); ylabel('count')
