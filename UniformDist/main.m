% MC, QMC and RQMC of a Uniform(0,1) i.i.d. random variable
% This is intended as an illustration of the concepts at
% http://discourse.vfitoolkit.com/t/simulating-a-markov-chain-mc-qmc-rqmc/687/2

% All three examples will use a Uniform(0,1) distribution.
% They use N=16 points to estimate the mean (the estimator is the sample mean)
% The results are based on M=1000 repetitions of the experiment, and report
%  - the mean (of the estimator of the mean), which can be evaluated as unbiased.
%  - the variance of the estimator
%  - a histogram of the estimates (with vertical line showing the mean)
% This is then repeated for N=128
% N is a power of 2 because the Sobol sequence is only evenly balanced when N=2^k
% (the first 2^k points place exactly one point in each interval [j/2^k,(j+1)/2^k)).

rng(1) % for reproducibility

N1=16;  % 2^4
N2=128; % 2^7

M = 1000; % number of independent repetitions of the experiment

fprintf('\nNote, the true mean = 0.5:\n')

%% First, we do Monte-Carlo simulation
% The code uses a random number generator to create N=16 and N=128 points.

MonteCarlo_Demo
% Figure 1 is the MC: plots the histogram over the M estimates

%% Second, we do Quasi-Monte Carlo simulation
% Two examples: stratification and Sobol sequence
%
% Note that M is not relevant here, there is a single simulation
% This is because QMC is deterministic, not stochastic, so repeating it would just give the identical result.

QuasiMonteCarlo_Demo
% Figure 2 is the QMC: plots the points used with N=16 by the QMC

%% Third, we do Randomized Quasi-Monte Carlo simulation
% Four examples: (i) stratification with random shift, (ii) Sobol sequence
% with 'random diagonal shift in base 2', (iii) Sobol sequence with Owen
% (nested uniform) scrambling, (iv) Sobol sequence with Matlab's built-in
% 'MatousekAffineOwen' scramble, which despite the name is not Owen scrambling
% (it is a random affine matrix scramble plus a random digital shift).

RandomizedQuasiMonteCarlo_Demo
% Figure 3 is the RQMC: plots the histogram over the M estimates

%% Fourth, asymptotics
% RQMC (iii), the Sobol sequence with Owen scrambling
% Theory says this is asymptotically normal. We contrast it with RQMC (iv).
RQMC_CheckNormal
% Figure 4 plots a fitted normal distribution alongside a repeat of the
% histograms for N=16 and N=128, for RQMC (iii) and RQMC (iv).

% The codes also run a 'Lilliefors test' for whether the RQMC (iii) and (iv)
% histograms are normally distributed. The test cannot reject normality for
% (iii), as the theory says. It strongly rejects normality for (iv): as seen in
% the figures, (iv) has a massive 'additional' mass at the exact solution
% (population value). (iii) and (iv) have about the same variance, so the
% difference is in the shape of the distribution, not its spread.
% In one dimension, Owen scrambling the first N=2^k Sobol points is exactly
% stratified sampling (one independent Uniform draw in each interval
% [j/N,(j+1)/N)), so the estimator is a mean of N independent terms and the
% usual CLT applies. The affine scramble makes the deeper digits of the points
% linear functions of their leading digits, so they are not independent.

%% Some things to note
% All these estimates are pretty accurate, because the exercise is pretty simple.
% QMC is very accurate, but also is deterministic.
% RQMC has much lower variances than MC: MC is around 5e-3 (N=16) and 7e-4 (N=128),
% RQMC (i) and (ii) are around 3e-4 and 5e-6 [rate 1/N^2], and RQMC (iii) and (iv)
% are around 2e-5 and 4e-8 to 6e-8 [rate 1/N^3].
% Unsurprisingly the N=128 estimates have lower variance than the N=16 estimates.
