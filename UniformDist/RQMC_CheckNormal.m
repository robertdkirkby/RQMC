%% RQMC_CheckNormal: is the RQMC estimator normally distributed?
% Plain Monte Carlo estimators are approximately normal by the central limit
% theorem. RQMC estimators need not be: the low-discrepancy structure that
% shrinks their variance can also distort the shape of the sampling distribution.
% Theory says Sobol + Owen (nested uniform) scrambling does satisfy a CLT
% [https://doi.org/10.1214/aos/1059655914]. Here we check this, and contrast it
% with Matlab's 'MatousekAffineOwen' scramble, which has the same variance but is
% not Owen scrambling. We use both a formal test and a visual comparison against
% a fitted normal.
%
% Requires est_owen, est_affine and Nvals from RandomizedQuasiMonteCarlo_Demo.

%% Lilliefors normality test
% (Lilliefors = Kolmogorov-Smirnov with the mean and variance estimated from
% the data, which is the right variant here since we do not know them a priori.)
% A small p-value means we reject the hypothesis that the estimates are normal.
% A warning that "P is less than the smallest tabulated value" means the true
% p is below 0.001, i.e. normality is rejected even more strongly. A warning
% that "P is greater than the largest tabulated value" means the true p is above 0.5.
normmethods = {'(iii) Sobol + Owen scramble', '(iv)  Sobol + MatousekAffineOwen'};
normest     = {est_owen, est_affine};

fprintf('\nNormality of the RQMC estimators (Lilliefors test):\n')
for k = 1:2
    E = normest{k};
    fprintf('%s\n', normmethods{k})
    for jN = 1:2
        [h, pL] = lillietest(E(:,jN));
        fprintf('   N=%d: h = %d, p = %.4f  (h=1 => reject normality at 5%%)\n', ...
                Nvals(jN), h, pL)
        if h == 1
            fprintf('   ==> Lilliefors test rejects normal distribution\n')
        else
            fprintf('   ==> Lilliefors test cannot reject normal distribution\n')
        end
    end
end

%% Figure 4: histograms with a fitted normal distribution
% histfit overlays the maximum-likelihood normal fit on the histogram, so we
% can visually judge how close each estimator is to normal.
% Columns are (iii) Owen and (iv) MatousekAffineOwen, rows are N=N1 and N=N2.
figure(4)
for k = 1:2
    E = normest{k};
    for jN = 1:2
        subplot(2, 2, (jN-1)*2 + k)
        histfit(E(:,jN), 30, 'normal')
        xline(0.5, 'r--', 'LineWidth', 1.5)
        title(sprintf('%s, N=%d  (var = %.2e)', normmethods{k}, Nvals(jN), var(E(:,jN))))
        xlabel('estimate'); ylabel('count')
    end
end
