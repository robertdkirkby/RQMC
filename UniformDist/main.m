% MC, QMC and RQMC of a Uniform(0,1) i.i.d. random variable
% This is intended as an illustration of the concepts at
% http://discourse.vfitoolkit.com/t/simulating-a-markov-chain-mc-qmc-rqmc/687/2

% All three examples will use a Uniform(0,1) distribution.
% They use N=20 points to estimate the mean (the estimator is the sample mean)
% The results are based on M=1000 
%  - the mean (of the estimator of the mean), which can be evaluated as unbiased.
%  - the variance of the estimator
%  - a histogram of the estimates (with vertical line showing the mean)
% This is then repeated for N=100

rng(1) % for reproducibility

N1=20;
N2=100;

M = 1000; % number of independent repetitions of the experiment

fprintf('\nNote, the true mean = 0.5:\n')

%% First, we do Monte-Carlo simulation
% The code uses a random number generator to create N=20 and N=100 points.

MonteCarlo_Demo
% Figure 1 is the MC: plots the histogram over the M estimates

%% Second, we do Quasi-Monte Carlo simulation
% Two examples: stratification and Sobol sequence
%
% Note that M is not relevant here, there is a single simluation
% This is because QMC is deterministic, not stochastic, so repeating it would just give the identical result.

QuasiMonteCarlo_Demo
% Figure 2 is the QMC: plots the points used with N=20 by the QMC

%% Third, we do Randomized Quasi-Monte Carlo simulation
% Three examples: (i) stratification with random shift, (ii) Sobol sequence
% with 'random diagonal shift in base 2', (iii) Sobol sequence with Owen
% scrambling.

RandomizedQuasiMonteCarlo_Demo
% Figure 3 is the RQMC: plots the histogram over the M estimates

%% Fourth, asymptotics
% RQMC (iii), the Sobol sequence with Owen scrambling
% Asymptotics suggest that this will be asymptotically normal
RQMC_CheckNormal
% Figure 4 plots a fitted normal distribution alongside a
% repeat of the histograms for N=20 and N=100 with RQMC (iii).

% The codes also run a 'Lilliefors test' for whether the RQMC (iii) histogram
% is normally distributed. This test strongly rejects normal distribution.
% As seen in the figures the issue is a massive 'additional' mass around the
% exact solution (population value), and that the rest is close to normal.

%% Some things to note
% All these estimates are pretty accurate, because the exercise is pretty simple.
% QMC is very accurate, but also is determinisitic.
% RQMC has much lower variances (around 1e-4 or 1e-6) than MC (around 1e-3)
% Unsurprsingly the N=100 estimates have lower variance than the N=20 estimates.
