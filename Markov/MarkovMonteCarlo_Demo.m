%% Monte Carlo Demo: simulating a T-period path of a Markov chain
% The path is S_1,...,S_T. S_1 is drawn from the stationary distribution, then
% each S_t is drawn from row S_{t-1} of pi_S, using T i.i.d. Uniform(0,1) draws.
% Repeating the experiment M times gives us the sampling distribution of the
% estimator of E[S] and of the other measures of how representative the path is.
%
% This is the baseline (plain Monte Carlo). Later we compare against QMC and RQMC.

fprintf('\nMC results:\n')

Tvals = [T1 T2];

est_ES_mc    = zeros(M, 2); % estimate of E[S]
transerr_mc  = zeros(M, 2); % transition error
staterr_mc   = zeros(M, 2); % stationary error
counterr_mc  = zeros(M, 2); % count error
halferr_mc   = zeros(M, 2); % first-half transition error
durerr_mc    = zeros(M, 2); % duration error

for jT = 1:2
    T = Tvals(jT);
    for m = 1:M
        U = rand(T,1);             % T i.i.d. Uniform(0,1) draws
        S = zeros(T,1);
        [~,S(1)] = max(cumsum_statdist'>U(1));
        for t = 2:T
            [~,S(t)] = max(cumsum_pi_S(S(t-1),:)>U(t));
        end

        % Measures of how representative the path is
        counts = accumarray([S(1:T-1),S(2:T)], 1, [N_S,N_S]);      % number of each (s,s') transition
        Phat = counts./max(sum(counts,2),1);                         % empirical transition matrix
        transerr_mc(m,jT) = max(max(abs(Phat-pi_S)));
        statdisthat = accumarray(S, 1, [N_S,1])/T;                   % empirical distribution
        staterr_mc(m,jT) = max(abs(statdisthat-statdist));
        counterr_mc(m,jT) = max(max(abs(counts-(T-1)*statdist.*pi_S)));
        counts_half = accumarray([S(1:T/2-1),S(2:T/2)], 1, [N_S,N_S]);
        Phat_half = counts_half./max(sum(counts_half,2),1);
        halferr_mc(m,jT) = max(max(abs(Phat_half-pi_S)));
        % Duration error: spell lengths vs the geometric distribution implied by pi_S (drop the first and last spell, as they are cut off)
        spellstart = [1; find(diff(S)~=0)+1];
        spelllength = diff([spellstart; T+1]);
        spellstate = S(spellstart);
        spelllength = spelllength(2:end-1);
        spellstate = spellstate(2:end-1);
        if isempty(spelllength)
            durerr_mc(m,jT) = 1; % path (almost) never changes state
        else
            pss = pi_S(sub2ind([N_S,N_S], spellstate, spellstate));     % probability of staying, for each spell
            ells = 1:max(spelllength);
            durerr_mc(m,jT) = max(abs(mean(spelllength<=ells,1)-mean(1-pss.^ells,1)));
        end
        est_ES_mc(m,jT) = mean(S_grid(S));                           % the estimator of E[S] for this repetition
    end

    fprintf('T=%d: mean of estimator of E[S] = %.6f, variance of estimator = %.3e\n', T, mean(est_ES_mc(:,jT)), var(est_ES_mc(:,jT)))
    fprintf('       transition error: mean = %.4f, worst = %.4f\n', mean(transerr_mc(:,jT)), max(transerr_mc(:,jT)))
    fprintf('       stationary error: mean = %.4f, worst = %.4f\n', mean(staterr_mc(:,jT)), max(staterr_mc(:,jT)))
    fprintf('       count error:      mean = %.2f, worst = %.2f\n', mean(counterr_mc(:,jT)), max(counterr_mc(:,jT)))
    fprintf('       first-half transition error: mean = %.4f, worst = %.4f\n', mean(halferr_mc(:,jT)), max(halferr_mc(:,jT)))
    fprintf('       duration error:   mean = %.4f, worst = %.4f\n', mean(durerr_mc(:,jT)), max(durerr_mc(:,jT)))
end

%% Histograms of the sampling distribution of the estimator of E[S]
figure(1)
for jT = 1:2
    subplot(2,1,jT)
    histogram(est_ES_mc(:,jT), 30)
    xline(trueES, 'r--', 'LineWidth', 1.5)
    title(sprintf('MC estimator of E[S],  T=%d  (var = %.2e)', Tvals(jT), var(est_ES_mc(:,jT))))
    xlabel('estimate'); ylabel('count')
end
