%% MarkovRQMC_CheckNormal: is the RQMC estimator of E[S] normally distributed?
% For plain Monte Carlo the estimator of E[S] from a path of a Markov chain is
% approximately normal (central limit theorem for Markov chains). For the RQMC
% methods we check, using both a formal test and a visual comparison against a
% fitted normal.
%
% Requires est_ES, methods, nmethods and Tvals from MarkovRandomizedQuasiMonteCarlo_Demo.

%% Lilliefors normality test
% (Lilliefors = Kolmogorov-Smirnov with the mean and variance estimated from
% the data, which is the right variant here since we do not know them a priori.)
% A small p-value means we reject the hypothesis that the estimates are normal.
% A warning that "P is less than the smallest tabulated value" means the true
% p is below 0.001, i.e. normality is rejected even more strongly. A warning
% that "P is greater than the largest tabulated value" means the true p is above 0.5.
fprintf('\nNormality of the RQMC estimators of E[S] (Lilliefors test):\n')
for k = 1:nmethods
    fprintf('%s\n', methods{k})
    for jT = 1:2
        [h, pL] = lillietest(est_ES(:,jT,k));
        fprintf('   T=%d: h = %d, p = %.4f  (h=1 => reject normality at 5%%)\n', ...
                Tvals(jT), h, pL)
        if h == 1
            fprintf('   ==> Lilliefors test rejects normal distribution\n')
        else
            fprintf('   ==> Lilliefors test cannot reject normal distribution\n')
        end
    end
end

%% Figure 4: histograms with a fitted normal distribution, T=T2
% histfit overlays the maximum-likelihood normal fit on the histogram, so we
% can visually judge how close each estimator is to normal.
figure(4)
for k = 1:nmethods
    subplot(nmethods, 1, k)
    histfit(est_ES(:,2,k), 30, 'normal')
    xline(trueES, 'r--', 'LineWidth', 1.5)
    title(sprintf('%s, T=%d  (var = %.2e)', methods{k}, Tvals(2), var(est_ES(:,2,k))))
    xlabel('estimate'); ylabel('count')
end
