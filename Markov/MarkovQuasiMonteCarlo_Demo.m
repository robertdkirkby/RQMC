%% Quasi-Monte Carlo Demo: simulating a T-period path of a Markov chain
% QMC replaces the random u's of Monte Carlo with deterministic ones. Because
% the u's are deterministic there is no randomness to repeat: a single
% "simulation" gives a single path, so M is not used here.
%
% We show four ways of choosing the u's (in all four, S_1 is the state at u=0.5
% of the stationary distribution, except (a) and (c) which use their first point):
%   (a) Single Sobol sequence - the t-th point of the Sobol sequence drives
%       period t. This deliberately shows that you cannot just swap rand for
%       Sobol: consecutive 1-D Sobol points are far from independent (e.g.
%       they alternate between [0,0.5) and [0.5,1)), and a Markov chain uses
%       consecutive u's in combination.
%   (b) One Sobol sequence per state - the k-th time the path is in state s, it
%       uses the k-th point of sequence s. So the departures from each state
%       are balanced, whatever the order the states are visited in.
%   (c) CUD sequence - a sequence is 'completely uniformly distributed' (CUD) if
%       every block of d consecutive points is evenly spread in [0,1)^d, for all
%       d. Small full-period linear congruential generators (LCG), run over their
%       whole period, behave like this (Owen & Tribble, 2005). We use
%       x_{t+1}=mod(a*x_t+1,T), u_t=x_t/T, which has period exactly T (as T is a
%       power of 2 and a=1 mod 4). The multiplier a is chosen by the spectral test
%       (a=157 for T=256, a=977 for T=2048).
%   (d) Rotor-router - each state has a fixed cycle listing its next states, and
%       the k-th time the path is in state s it goes to the k-th entry of the
%       cycle of s. The cycle is built so that every stretch of it is
%       balanced: each entry goes to the next state that is most behind its quota
%       (k*pi_S(s,s') minus the number of times s' has been used so far).
%       This needs pi_S to be rational, so we round each row of pi_S to
%       multiples of 1/L with L=2^16.

fprintf('\nQMC results:\n')

Tvals = [T1 T2];
lcg_a = [157 977];  % LCG multipliers for T1 and T2

methods_qmc = {'(a) Single Sobol sequence', ...
               '(b) One Sobol sequence per state', ...
               '(c) CUD (full-period LCG)', ...
               '(d) Rotor-router'};

% Rotor-router cycles (used by (d), and also by RQMC (iv))
L = 2^16;
rotorcycle = zeros(N_S, L);
for ii = 1:N_S
    % Round L*pi_S(ii,:) to integers that sum to L (largest remainder)
    nL = floor(L*pi_S(ii,:));
    [~, ord] = sort(L*pi_S(ii,:)-nL, 'descend');
    nshort = L-sum(nL);
    nL(ord(1:nshort)) = nL(ord(1:nshort))+1;
    % Entry k goes to the next state that is most behind its quota
    used = zeros(1, N_S);
    for k = 1:L
        [~, jj] = max(k*nL/L-used);
        used(jj) = used(jj)+1;
        rotorcycle(ii,k) = jj;
    end
end

p = sobolset(1);  % 1-D Sobol point set

est_ES_qmc   = zeros(4, 2); % estimate of E[S]
transerr_qmc = zeros(4, 2); % transition error
staterr_qmc  = zeros(4, 2); % stationary error
counterr_qmc = zeros(4, 2); % count error
halferr_qmc  = zeros(4, 2); % first-half transition error
durerr_qmc   = zeros(4, 2); % duration error
Spaths_qmc   = cell(2, 1);  % the paths, for figure 2

for jT = 1:2
    T = Tvals(jT);
    Spaths = zeros(T, 4);

    % (a) Single Sobol sequence
    U = net(p, T);
    S = zeros(T,1);
    [~,S(1)] = max(cumsum_statdist'>U(1));
    for t = 2:T
        [~,S(t)] = max(cumsum_pi_S(S(t-1),:)>U(t));
    end
    Spaths(:,1) = S;

    % (b) One Sobol sequence per state
    Ustate = net(p, T);        % the same deterministic sequence for every state
    visits = zeros(N_S,1);     % number of times we have left each state so far
    S = zeros(T,1);
    [~,S(1)] = max(cumsum_statdist'>0.5);
    for t = 2:T
        visits(S(t-1)) = visits(S(t-1))+1;
        [~,S(t)] = max(cumsum_pi_S(S(t-1),:)>Ustate(visits(S(t-1))));
    end
    Spaths(:,2) = S;

    % (c) CUD: full-period LCG
    x = zeros(T,1);
    for t = 2:T
        x(t) = mod(lcg_a(jT)*x(t-1)+1, T);
    end
    U = x/T;
    S = zeros(T,1);
    [~,S(1)] = max(cumsum_statdist'>U(1));
    for t = 2:T
        [~,S(t)] = max(cumsum_pi_S(S(t-1),:)>U(t));
    end
    Spaths(:,3) = S;

    % (d) Rotor-router, every state starts at the beginning of its cycle
    visits = zeros(N_S,1);
    S = zeros(T,1);
    [~,S(1)] = max(cumsum_statdist'>0.5);
    for t = 2:T
        visits(S(t-1)) = visits(S(t-1))+1;
        S(t) = rotorcycle(S(t-1), mod(visits(S(t-1))-1, L)+1);
    end
    Spaths(:,4) = S;

    % Measures of how representative each path is
    for k = 1:4
        S = Spaths(:,k);
        counts = accumarray([S(1:T-1),S(2:T)], 1, [N_S,N_S]);
        Phat = counts./max(sum(counts,2),1);
        transerr_qmc(k,jT) = max(max(abs(Phat-pi_S)));
        statdisthat = accumarray(S, 1, [N_S,1])/T;
        staterr_qmc(k,jT) = max(abs(statdisthat-statdist));
        counterr_qmc(k,jT) = max(max(abs(counts-(T-1)*statdist.*pi_S)));
        counts_half = accumarray([S(1:T/2-1),S(2:T/2)], 1, [N_S,N_S]);
        Phat_half = counts_half./max(sum(counts_half,2),1);
        halferr_qmc(k,jT) = max(max(abs(Phat_half-pi_S)));
        % Duration error: spell lengths vs the geometric distribution implied by pi_S (drop the first and last spell, as they are cut off)
        spellstart = [1; find(diff(S)~=0)+1];
        spelllength = diff([spellstart; T+1]);
        spellstate = S(spellstart);
        spelllength = spelllength(2:end-1);
        spellstate = spellstate(2:end-1);
        if isempty(spelllength)
            durerr_qmc(k,jT) = 1; % path (almost) never changes state
        else
            pss = pi_S(sub2ind([N_S,N_S], spellstate, spellstate));     % probability of staying, for each spell
            ells = 1:max(spelllength);
            durerr_qmc(k,jT) = max(abs(mean(spelllength<=ells,1)-mean(1-pss.^ells,1)));
        end
        est_ES_qmc(k,jT) = mean(S_grid(S));
    end
    Spaths_qmc{jT} = Spaths;
end

for k = 1:4
    fprintf('%s\n', methods_qmc{k})
    for jT = 1:2
        fprintf('   T=%d: estimate of E[S] = %.6f, error = %+.2e, transition error = %.4f, stationary error = %.4f, count error = %.2f, first-half transition error = %.4f, duration error = %.4f\n', ...
                Tvals(jT), est_ES_qmc(k,jT), est_ES_qmc(k,jT)-trueES, transerr_qmc(k,jT), staterr_qmc(k,jT), counterr_qmc(k,jT), halferr_qmc(k,jT), durerr_qmc(k,jT))
    end
end

%% Plot the first 100 periods of each path (T=2048)
% There is no sampling distribution to histogram here (QMC gives a single
% deterministic path), so instead we plot the paths.
figure(2)
for k = 1:4
    subplot(4,1,k)
    stairs(1:100, S_grid(Spaths_qmc{2}(1:100,k)), 'b')
    title(sprintf('%s, T=%d, first 100 periods  (estimate of E[S] = %.4f)', methods_qmc{k}, T2, est_ES_qmc(k,2)))
    xlabel('period'); ylabel('S')
    ylim([min(S_grid) max(S_grid)]+[-0.1 0.1]*(max(S_grid)-min(S_grid)))
end
