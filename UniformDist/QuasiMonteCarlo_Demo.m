%% Quasi-Monte Carlo Demo: estimating E[U], U~Uniform(0,1)
% QMC replaces the random draws of Monte Carlo with a deterministic set of
% points that is deliberately spread evenly across [0,1] (low discrepancy).
% Because the point set is deterministic there is no randomness to repeat: a
% single "simulation" of N points gives a single estimate, so M is not used
% here. The pay-off is that the estimation error typically shrinks faster in N
% than the O(1/sqrt(N)) of plain Monte Carlo.
%
% We show two ways of constructing the points:
%   (1) Stratification - split [0,1] into N equal strata and take the midpoint
%                        of each. Points are (i-0.5)/N, i=1..N.
%   (2) Sobol sequence - a classic low-discrepancy sequence.
% The true mean is 1/2.

fprintf('\nQMC results:\n')

%% (1) Stratification
strat_N1 = ((1:N1)'-0.5)/N1;   % N1 equally spaced midpoints on [0,1]
strat_N2 = ((1:N2)'-0.5)/N2;   % N2 equally spaced midpoints on [0,1]

est_strat_N1 = mean(strat_N1);
est_strat_N2 = mean(strat_N2);

fprintf('QMC - Stratification (midpoints):\n')
fprintf('  N=20:  estimate = %.6f, error = %+.2e\n', est_strat_N1, est_strat_N1-0.5)
fprintf('  N=100: estimate = %.6f, error = %+.2e\n', est_strat_N2, est_strat_N2-0.5)

%% (2) Sobol sequence
p = sobolset(1);               % 1-dimensional Sobol point set
sobol_N1 = net(p, N1);         % first N1 points of the sequence
sobol_N2 = net(p, N2);         % first N2 points of the sequence

est_sobol_N1 = mean(sobol_N1);
est_sobol_N2 = mean(sobol_N2);

fprintf('QMC - Sobol sequence:\n')
fprintf('  N=20:  estimate = %.6f, error = %+.2e\n', est_sobol_N1, est_sobol_N1-0.5)
fprintf('  N=100: estimate = %.6f, error = %+.2e\n', est_sobol_N2, est_sobol_N2-0.5)

%% Visualize the point sets (N=20) to see how evenly they cover [0,1]
% There is no sampling distribution to histogram here (QMC gives a single
% deterministic estimate), so instead we plot where the N points land.
figure(2)
subplot(2,1,1)
plot(strat_N1, zeros(N1,1), 'o', 'MarkerFaceColor', 'b')
xline(0.5, 'r--', 'LineWidth', 1.5)
title(sprintf('QMC Stratification points, N=20  (estimate = %.4f)', est_strat_N1))
xlabel('point location on [0,1]'); xlim([0 1]); yticks([])

subplot(2,1,2)
plot(sobol_N1, zeros(N1,1), 'o', 'MarkerFaceColor', 'b')
xline(0.5, 'r--', 'LineWidth', 1.5)
title(sprintf('QMC Sobol points, N=20  (estimate = %.4f)', est_sobol_N1))
xlabel('point location on [0,1]'); xlim([0 1]); yticks([])
