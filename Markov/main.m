% MC, QMC and RQMC of a Markov chain
% This is intended as an illustration of the concepts at
% http://discourse.vfitoolkit.com/t/simulating-a-markov-chain-mc-qmc-rqmc/687/3
% It follows the same structure as the codes for the Uniform(0,1) i.i.d. case (see ../UniformDist/)
%
% The five-state chain (whichchain=2) uses the Farmer-Toda discretization command from VFI Toolkit,
% so VFI Toolkit must be on the Matlab path. You can find the command on github at
% https://github.com/vfitoolkit/VFIToolkit-matlab/blob/master/DiscretizationMethods/discretizeAR1_FarmerToda/discretizeAR1_FarmerToda.m

% We simulate a single path of T periods of a Markov chain, S_1,...,S_T, as is
% needed by the Matched-Expectations Path (MEP) algorithm. Every method uses
% the same 'inverse-cdf' step: given current state s and a number u in [0,1),
% next period state is the first s' with cumsum(pi_S(s,:))>u. The methods
% differ only in where the u come from.
% The results are based on M=1000 repetitions of the experiment, and report
%  - estimate of E[S] (the mean of the path), which can be evaluated as unbiased.
%  - the variance of the estimator of E[S]
%  - 'transition error': max abs difference between the empirical transition matrix of the path and pi_S
%  - 'stationary error': max abs difference between the empirical distribution of the path and the stationary distribution
%  - 'count error': max abs difference between the number of each (s,s') transition and its expected number
%  - 'first-half transition error': transition error using only the first half of the path
%  - 'duration error': max abs difference between the cdf of the lengths of spells
%    in a state, and the geometric distribution they should follow given pi_S.
%    A path can have exactly the right transition counts, but with every spell
%    the same length, so we need to check this separately.
%  - a histogram of the estimates of E[S] (with vertical line showing the true value)
% This is done for T=256 and T=2048

rng(1) % for reproducibility

T1=256;  % 2^8
T2=2048; % 2^11

M = 1000; % number of independent repetitions of the experiment

%% The Markov chain
% whichchain=1: two-state chain, like the aggregate shock in Krusell & Smith (1998), expected duration of 8 periods
% whichchain=2: five-state Farmer-Toda discretization of an AR(1) with autocorrelation 0.95
whichchain=1;

if whichchain==1
    S_grid=[0.99; 1.01];
    pi_S=[0.875, 0.125; 0.125, 0.875];
elseif whichchain==2
    rho=0.95;
    sigma_epsilon=0.007;
    n_S=5;
    % Farmer-Toda discretization (uses the VFI Toolkit, so it must be on the path)
    farmertodaoptions.parallel=0; % cpu, nothing in these codes uses the gpu
    [S_grid,pi_S]=discretizeAR1_FarmerToda(0,rho,sigma_epsilon,n_S,farmertodaoptions);
    S_grid=gather(S_grid);
    pi_S=gather(pi_S);
end
N_S=length(S_grid);
cumsum_pi_S=cumsum(pi_S,2);
cumsum_pi_S(:,end)=1; % guard against rounding error leaving the cumsum slightly below 1

% Stationary distribution of the chain (eigenvector of pi_S' with eigenvalue 1)
[V,D]=eig(pi_S');
[~,unitroot]=min(abs(diag(D)-1));
statdist=V(:,unitroot)/sum(V(:,unitroot));
cumsum_statdist=cumsum(statdist);
cumsum_statdist(end)=1;

trueES=statdist'*S_grid;
fprintf('\nNote, the true mean E[S] = %.6f:\n', trueES)

% In every method the initial state S_1 is drawn from the stationary
% distribution (using the inverse-cdf of statdist), so there is no burn-in.

%% First, we do Monte-Carlo simulation
% The code uses a random number generator to create the T-1 u's that drive the path
% (this is what the VFI Toolkit currently does for MEP).

MarkovMonteCarlo_Demo
% Figure 1 is the MC: plots the histogram over the M estimates of E[S]

%% Second, we do Quasi-Monte Carlo simulation
% Four examples:
%   (a) a single Sobol sequence drives the whole path
%   (b) one Sobol sequence per state: the k-th time the path is in state s it uses the k-th point of sequence s
%   (c) a CUD sequence: a full-period linear congruential generator, run over its whole period
%   (d) a rotor-router: each state cycles through its next states in a fixed, balanced, order
%
% Note that M is not relevant here, there is a single simulation
% This is because QMC is deterministic, not stochastic, so repeating it would just give the identical result.

MarkovQuasiMonteCarlo_Demo
% Figure 2 is the QMC: plots the first 100 periods of the path from each

%% Third, we do Randomized Quasi-Monte Carlo simulation
% Seven examples:
%   (i)   CUD sequence with a random shift
%   (ii)  one Sobol sequence per state, each with a random digital shift in base 2
%   (iii) one Sobol sequence per state, each with Owen (nested uniform) scrambling
%   (iv)  randomized rotor-router: each state starts at a random point of its cycle
%   (v)   fixed counts: a uniformly random path among all paths with exactly the
%         same transition counts as the path from (iv). [We take the counts from
%         the path from (iv), rather than rounding (T-1)*pi_i*p_ij, as the counts of
%         an actual path automatically have every state entered as often as it is
%         left.] This is not RQMC, it is
%         'conditional Monte Carlo'. (iv) is the RQMC version of (v): same counts,
%         but each state's next states come in a balanced order, rather than a random order.
%   (vi)  fixed two-step counts: as (v), but a uniformly random path among all
%         paths with exactly the same counts of each (s,s',s'') as a
%         'second-order' randomized rotor-router path (one rotor per pair of
%         previous and current state). Again this is 'conditional Monte Carlo'.
%   (vii) fixed counts, stratified spell lengths: the same transition counts as
%         (iv) and (v), with the order of the states of the spells drawn as in
%         (v), but the lengths of the spells in each state drawn by stratified
%         sampling from a geometric distribution (rather than at random).

MarkovRandomizedQuasiMonteCarlo_Demo
% Figure 3 is the RQMC: plots the histogram over the M estimates of E[S]

%% Fourth, asymptotics
% Is the RQMC estimator of E[S] normally distributed?
MarkovRQMC_CheckNormal
% Figure 4 plots a fitted normal distribution alongside a repeat of the
% histograms for T=2048, for each RQMC method.

%% Some things to note
% Results from the two-state chain (whichchain=1), T=2048, with whichchain=2 in square brackets.
% Variance of the estimator of E[S]:
%   MC 3e-7 [9e-6], (i) CUD 2e-8 [2e-6], (ii)-(v) about 1e-9 [about 4e-7]
% Count error (worst over the M repetitions):
%   MC 192 [352], (i) CUD 43 [172], (ii)-(v) 7 [39-49]
% Duration error (mean):
%   MC 0.05 [0.06], (i) CUD 0.03 [0.05], (ii)-(iv) 0.60 [0.28-0.39], (v) 0.04 [0.04]
%
% QMC (a), a single Sobol sequence, fails badly: the path spends far too much
% time in one state. You cannot just swap rand for a low-discrepancy sequence.
%
% Balancing the departures from each state, which is what (ii), (iii) and
% (iv) do, gets the transition counts almost exactly right, but it also makes
% the path far too regular: in the two-state chain every spell lasts exactly 8
% periods, when they should be geometric with mean 8. Owen scrambling does not
% fix this, as every block of 8 consecutive points in a state's Sobol sequence
% gets the same scramble of its leading digits, so the departure is in the same
% place in every block.
%
% (v) has exactly the same transition counts as (iv), but because the order is
% random the spell lengths are as good as MC. This is because every path with
% the same transition counts (and first state) has the same probability under
% the Markov chain, so (v) is exactly MC conditional on the transition counts.
%
% So the RQMC version of (v), which is (iv), is better on the first half of the
% path (first-half transition error) but much worse on the spell lengths.
%
% (vi), fixing the two-step counts rather than the one-step counts, is not an
% improvement on (v). The duration error is about the same (0.034 vs 0.036
% [0.044 vs 0.044]), while the count error is a bit worse (worst 8 vs 7 [65 vs 39])
% and the variance of the estimator of E[S] is higher (1.3e-9 vs 0.9e-9 [6.3e-7 vs 3.8e-7]).
% The second-order rotor has a separate rotor for each pair of states, each
% visited less often, so the one-step counts are balanced less well. And (v)
% already has the right distribution of spell lengths (it is MC conditional on
% the counts), so what is left of its duration error is sampling noise in the
% longer spells, which fixing the two-step counts does not remove (it pins down
% only the number of spells of length one).
%
% (vii) is the improved version of (v). It has exactly the same transition counts
% as (iv) and (v), and the same order of the states of the spells as (v) would
% (random), but the spell lengths in each state are a stratified sample of a
% geometric distribution, rather than a random one. Duration error (mean, worst):
%   MC 0.045, 0.145 [0.056, 0.163]; (v) 0.037, 0.092 [0.044, 0.092]; (vii) 0.011, 0.036 [0.022, 0.070]
% So (vii) gets the transition counts as well as the rotor-router does, and the
% spell lengths much better than MC does. Its first-half transition error is
% the same as (v), as the order is still random.
%
% (i) CUD improves on MC for everything, but by much less than (ii)-(v) do
% on the transition counts.
%
% Unlike the Uniform(0,1) i.i.d. case, Lilliefors rejects normality for every
% method in the two-state chain (the estimator only takes a few distinct
% values). In the five-state chain it cannot reject normality for any method
% when T=2048 (when T=256 it rejects for (iv), (v) and (vii)).
