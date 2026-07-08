%% Randomized Quasi-Monte Carlo Demo: estimating E[U], U~Uniform(0,1)
% RQMC takes the evenly-spread QMC points and applies a *random* transformation
% that (a) preserves the low-discrepancy structure, so the estimator keeps QMC's
% low variance, but (b) makes each point marginally Uniform(0,1), so the
% estimator is unbiased and we can measure its variance by repeating M times.
% Hence -- unlike deterministic QMC -- we again get a sampling distribution and
% can draw a histogram (compare figure 1 for plain MC), but with much smaller
% variance than MC at the same N.
%
% Three examples:
%   (i)   Stratification with a random shift (Cranley-Patterson rotation):
%         take the equally-spaced grid and shift all points by a common
%         U~Uniform(0,1), wrapping around mod 1.
%   (ii)  Sobol sequence with a random digital shift in base 2: XOR the binary
%         digits of every Sobol point with a common random bit string.
%   (iii) Sobol sequence with Owen (nested) scrambling.
% The true mean is 1/2.

Nvals = [N1 N2];

est_strat  = zeros(M, 2);   % (i)   stratification + random shift
est_dshift = zeros(M, 2);   % (ii)  Sobol + digital shift in base 2
est_owen   = zeros(M, 2);   % (iii) Sobol + Owen scramble

p     = sobolset(1);        % 1-D Sobol point set (reused for all randomizations)
nbits = 32;                 % number of bits used for the base-2 digital shift
scale = 2^nbits;

for jN = 1:2
    N = Nvals(jN);

    % Deterministic QMC building blocks for this N
    grid    = ((1:N)'-0.5)/N;                % equally-spaced stratification grid
    sob     = net(p, N);                     % first N Sobol points
    sob_int = uint64(floor(sob*scale));      % integer form, for the digital shift

    for m = 1:M
        % (i) stratification with a common random shift (Cranley-Patterson)
        U = rand;
        pts = mod(grid + U, 1);
        est_strat(m, jN) = mean(pts);

        % (ii) Sobol with a random digital shift in base 2 (common XOR mask)
        shift = uint64(randi([0 scale-1]));
        pts = double(bitxor(sob_int, shift))/scale;
        est_dshift(m, jN) = mean(pts);

        % (iii) Sobol with Owen scrambling (a fresh scramble each repetition)
        ps = scramble(p, 'MatousekAffineOwen');
        pts = net(ps, N);
        est_owen(m, jN) = mean(pts);
    end
end

%% Report mean and variance of each estimator
methods = {'(i)   Stratification + random shift', ...
           '(ii)  Sobol + digital shift (base 2)', ...
           '(iii) Sobol + Owen scramble'};
allest  = {est_strat, est_dshift, est_owen};

fprintf('\nRQMC results:\n')
for k = 1:3
    E = allest{k};
    fprintf('%s\n', methods{k})
    fprintf('   N=20:  mean = %.6f, variance = %.3e\n', mean(E(:,1)), var(E(:,1)))
    fprintf('   N=100: mean = %.6f, variance = %.3e\n', mean(E(:,2)), var(E(:,2)))
end

%% Histograms of the M estimates (figure 3)
% Rows are the three randomization methods, columns are N=20 and N=100.
% Note the variance is far smaller than the MC histograms in figure 1.
figure(3)
for k = 1:3
    E = allest{k};
    for jN = 1:2
        subplot(3, 2, (k-1)*2 + jN)
        histogram(E(:,jN), 30)
        xline(0.5, 'r--', 'LineWidth', 1.5)
        title(sprintf('%s, N=%d  (var = %.2e)', methods{k}, Nvals(jN), var(E(:,jN))))
        xlabel('estimate'); ylabel('count')
    end
end
