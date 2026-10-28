%% Randomized Quasi-Monte Carlo Demo: simulating a T-period path of a Markov chain
% RQMC takes the deterministic QMC constructions and applies a *random*
% transformation that keeps their balance, but makes each u marginally
% Uniform(0,1). Hence -- unlike deterministic QMC -- we again get a sampling
% distribution and can draw a histogram (compare figure 1 for plain MC).
%
% Seven examples:
%   (i)   CUD (full-period LCG) with a random shift: add a common U~Uniform(0,1)
%         to every point, wrapping around mod 1, and start the LCG at a random point.
%   (ii)  One Sobol sequence per state, each with its own random digital shift in base 2.
%   (iii) One Sobol sequence per state, each with its own Owen (nested uniform)
%         scrambling (same code as in ../UniformDist/RandomizedQuasiMonteCarlo_Demo.m).
%   (iv)  Randomized rotor-router: each state starts at a random point of its cycle.
%   (v)   Fixed counts: take the transition counts of the path from (iv), and
%         draw a path uniformly at random from all paths with exactly those
%         transition counts (and the same first and last state). This is
%         'conditional Monte Carlo' rather than RQMC: the counts are fixed but
%         the order is random. To draw the path: for every state other than the
%         last state, pick which next state it goes to on its final departure,
%         by drawing a random spanning tree that points towards the last state
%         (Wilson's algorithm); put the other departures from each state in a
%         random order, with the final departure last; then walk the path. This
%         gives every path with these counts the same probability (BEST theorem).
%         Note: rather than rounding (T-1)*pi_i*p_ij to integers, we take the
%         counts from the path from (iv). Rounding so that every state is entered
%         as often as it is left (allowing for the first and last states) is
%         awkward, while the counts of an actual path satisfy this automatically,
%         and the rotor-router counts are within a few transitions of (T-1)*pi_i*p_ij.
%         (iv) is the RQMC version of (v): same counts, but each state's
%         departures come in a balanced order rather than a random order.
%         Because every path with the same transition counts (and first state)
%         has the same probability under the Markov chain, (v) is exactly MC
%         conditional on the transition counts.
%   (vi)  Fixed two-step counts: as (v), but fixing the number of times each
%         (s,s',s'') occurs, rather than each (s,s'). Every path with the same
%         two-step counts (and the same first two states and last two states)
%         has the same probability under the Markov chain, so (vi) is exactly
%         MC conditional on the two-step counts. The counts come from a
%         'second-order' randomized rotor-router path, which has a separate
%         rotor for each pair (previous state, current state), each using the
%         cycle of the current state and starting at a random point of it.
%         The path is then drawn exactly as in (v), but on the graph whose
%         nodes are the pairs (s,s') and whose edges are the (s,s',s'').
%         Fixing the two-step counts fixes more of the spell lengths: e.g. the
%         number of spells of length one in state s' is the number of
%         (s,s',s'') with s~=s' and s''~=s'.
%   (vii) Fixed counts, stratified spell lengths: think of the path as a
%         sequence of spells, each spell being some periods in the same state.
%         Take the transition counts of the path from (iv). Draw the order of
%         the states of the spells uniformly at random among all orders with the
%         same counts of jumps (s to s'~=s), exactly as in (v) but without the
%         stays. Then for each state s, with K_s spells and a total of tau_s
%         periods, draw the K_s spell lengths by stratified sampling (one
%         uniform draw in each of K_s equal intervals) and the inverse-cdf of the
%         geometric distribution with stay probability (tau_s-K_s)/tau_s, add or
%         subtract one period from randomly chosen spells until they add up to
%         tau_s, and give them to the spells of s in a random order.
%         This has exactly the same transition counts as (iv) and (v). Under MC,
%         conditional on the counts, the spell lengths of s would be a uniformly
%         random way of splitting tau_s periods into K_s spells; (vii) replaces this
%         with a stratified version, so it is the RQMC version of (v) for the
%         spell lengths.
% Unlike the Uniform(0,1) i.i.d. case, the u's used by (i)-(iii) are not
% independent across periods, so the path is not exactly a draw of the Markov
% chain and the estimator of E[S] is not guaranteed to be unbiased for a given T.
%
% Requires rotorcycle and L from MarkovQuasiMonteCarlo_Demo.

Tvals = [T1 T2];
lcg_a = [157 977];  % LCG multipliers for T1 and T2 (see MarkovQuasiMonteCarlo_Demo)

methods = {'(i)   CUD + random shift', ...
           '(ii)  Sobol per state + digital shift (base 2)', ...
           '(iii) Sobol per state + Owen scramble', ...
           '(iv)  Randomized rotor-router', ...
           '(v)   Fixed counts, random order', ...
           '(vi)  Fixed two-step counts, random order', ...
           '(vii) Fixed counts, stratified spell lengths'};
nmethods = length(methods);

est_ES   = zeros(M, 2, nmethods); % estimate of E[S]
transerr = zeros(M, 2, nmethods); % transition error
staterr  = zeros(M, 2, nmethods); % stationary error
counterr = zeros(M, 2, nmethods); % count error
halferr  = zeros(M, 2, nmethods); % first-half transition error
durerr   = zeros(M, 2, nmethods); % duration error

p     = sobolset(1);        % 1-D Sobol point set (reused for all randomizations)
nbits = 32;                 % number of bits used for the base-2 digital shift and the Owen scramble
scale = 2^nbits;

for jT = 1:2
    T = Tvals(jT);

    % Deterministic QMC building blocks for this T
    sob     = net(p, T);                     % first T Sobol points
    sob_int = uint64(floor(sob*scale));      % integer form, for the digital shift and Owen scramble

    for m = 1:M
        Spaths = zeros(T, nmethods);

        % (i) CUD with a random shift (and random starting point of the LCG)
        x = zeros(T,1);
        x(1) = randi([0 T-1]);
        for t = 2:T
            x(t) = mod(lcg_a(jT)*x(t-1)+1, T);
        end
        U = mod(x/T + rand, 1);
        S = zeros(T,1);
        [~,S(1)] = max(cumsum_statdist'>U(1));
        for t = 2:T
            [~,S(t)] = max(cumsum_pi_S(S(t-1),:)>U(t));
        end
        Spaths(:,1) = S;

        % (ii) Sobol per state, each with its own random digital shift in base 2
        Ustate = zeros(T, N_S);
        for ii = 1:N_S
            shift = uint64(randi([0 scale-1]));
            Ustate(:,ii) = double(bitxor(sob_int, shift))/scale;
        end
        visits = zeros(N_S,1);
        S = zeros(T,1);
        [~,S(1)] = max(cumsum_statdist'>rand);
        for t = 2:T
            visits(S(t-1)) = visits(S(t-1))+1;
            [~,S(t)] = max(cumsum_pi_S(S(t-1),:)>Ustate(visits(S(t-1)),S(t-1)));
        end
        Spaths(:,2) = S;

        % (iii) Sobol per state, each with its own Owen (nested uniform) scrambling
        % Loop over digits j=1 (most significant) to nbits. All points that share the
        % same first j-1 digits get the same random flip of digit j; points with
        % different first j-1 digits get independent flips.
        for ii = 1:N_S
            pts_int = sob_int;
            for jj = 1:nbits
                bitpos = nbits-jj;                               % digit jj sits at bit position nbits-jj (0=least significant)
                prefix = bitshift(sob_int, -(bitpos+1));         % first jj-1 digits of each (unscrambled) point
                [~, ~, prefixgroup] = unique(prefix);            % which distinct prefix each point has
                flip = uint64(rand(max(prefixgroup),1) < 0.5);   % one random flip per distinct prefix
                pts_int = bitxor(pts_int, bitshift(flip(prefixgroup), bitpos));
            end
            Ustate(:,ii) = double(pts_int)/scale;
        end
        visits = zeros(N_S,1);
        S = zeros(T,1);
        [~,S(1)] = max(cumsum_statdist'>rand);
        for t = 2:T
            visits(S(t-1)) = visits(S(t-1))+1;
            [~,S(t)] = max(cumsum_pi_S(S(t-1),:)>Ustate(visits(S(t-1)),S(t-1)));
        end
        Spaths(:,3) = S;

        % (iv) Randomized rotor-router: each state starts at a random point of its cycle
        startpos = randi([0 L-1], N_S, 1);
        visits = zeros(N_S,1);
        S = zeros(T,1);
        [~,S(1)] = max(cumsum_statdist'>rand);
        for t = 2:T
            visits(S(t-1)) = visits(S(t-1))+1;
            S(t) = rotorcycle(S(t-1), mod(startpos(S(t-1))+visits(S(t-1))-1, L)+1);
        end
        Spaths(:,4) = S;

        % (v) Fixed counts (those of the path from (iv)), random order
        counts = accumarray([S(1:T-1),S(2:T)], 1, [N_S,N_S]);
        firststate = S(1);
        laststate = S(T);
        ndepart = sum(counts,2);
        % Final departure from each state (other than laststate): random spanning tree pointing to laststate (Wilson's algorithm)
        intree = false(N_S,1);
        intree(laststate) = true;
        finalnext = zeros(N_S,1);
        for ii = 1:N_S
            if ndepart(ii)>0
                ss = ii;
                while ~intree(ss) % random walk, using the transition counts as weights, until we hit the tree
                    [~,finalnext(ss)] = max(cumsum(counts(ss,:))/ndepart(ss)>rand);
                    ss = finalnext(ss);
                end
                ss = ii;
                while ~intree(ss) % add the (loop-erased) walk to the tree
                    intree(ss) = true;
                    ss = finalnext(ss);
                end
            end
        end
        % Departures from each state: the others in random order, then the final one
        departures = cell(N_S,1);
        for ii = 1:N_S
            deps = repelem(1:N_S, counts(ii,:));
            if ii~=laststate && ndepart(ii)>0
                deps(find(deps==finalnext(ii),1)) = [];
                deps = [deps(randperm(length(deps))), finalnext(ii)];
            else
                deps = deps(randperm(length(deps)));
            end
            departures{ii} = deps;
        end
        visits = zeros(N_S,1);
        S = zeros(T,1);
        S(1) = firststate;
        for t = 2:T
            visits(S(t-1)) = visits(S(t-1))+1;
            S(t) = departures{S(t-1)}(visits(S(t-1)));
        end
        Spaths(:,5) = S;

        % (vi) Fixed two-step counts, random order
        % A pair (s,s') is indexed as (s-1)*N_S+s'. Pair pp=(s,s') followed by s'' moves to pair (s',s'').
        % First, the second-order randomized rotor-router path, which gives the two-step counts
        startpos2 = randi([0 L-1], N_S^2, 1);
        visits2 = zeros(N_S^2,1);
        S = zeros(T,1);
        [~,S(1)] = max(cumsum_statdist'>rand);
        [~,S(2)] = max(cumsum_pi_S(S(1),:)>rand);
        for t = 3:T
            pp = (S(t-2)-1)*N_S+S(t-1);
            visits2(pp) = visits2(pp)+1;
            S(t) = rotorcycle(S(t-1), mod(startpos2(pp)+visits2(pp)-1, L)+1);
        end
        pairs = (S(1:T-1)-1)*N_S+S(2:T);                              % the pair at t=1,...,T-1
        counts2 = accumarray([pairs(1:T-2),S(3:T)], 1, [N_S^2,N_S]); % number of times pair pp is followed by s''
        firstpair = pairs(1);
        lastpair = pairs(T-1);
        ndepart2 = sum(counts2,2);
        % Final departure from each pair (other than lastpair): random spanning tree pointing to lastpair (Wilson's algorithm)
        intree = false(N_S^2,1);
        intree(lastpair) = true;
        finalnext = zeros(N_S^2,1); % the s'' of the final departure from each pair
        for pp = 1:N_S^2
            if ndepart2(pp)>0
                qq = pp;
                while ~intree(qq) % random walk, using the two-step counts as weights, until we hit the tree
                    [~,finalnext(qq)] = max(cumsum(counts2(qq,:))/ndepart2(qq)>rand);
                    qq = (mod(qq-1,N_S))*N_S+finalnext(qq); % pair (s,s') followed by s'' is pair (s',s'')
                end
                qq = pp;
                while ~intree(qq) % add the (loop-erased) walk to the tree
                    intree(qq) = true;
                    qq = (mod(qq-1,N_S))*N_S+finalnext(qq);
                end
            end
        end
        % Departures from each pair: the others in random order, then the final one
        departures2 = cell(N_S^2,1);
        for pp = 1:N_S^2
            deps = repelem(1:N_S, counts2(pp,:));
            if pp~=lastpair && ndepart2(pp)>0
                deps(find(deps==finalnext(pp),1)) = [];
                deps = [deps(randperm(length(deps))), finalnext(pp)];
            else
                deps = deps(randperm(length(deps)));
            end
            departures2{pp} = deps;
        end
        visits2 = zeros(N_S^2,1);
        S = zeros(T,1);
        S(1) = ceil(firstpair/N_S);
        S(2) = mod(firstpair-1,N_S)+1;
        for t = 3:T
            pp = (S(t-2)-1)*N_S+S(t-1);
            visits2(pp) = visits2(pp)+1;
            S(t) = departures2{pp}(visits2(pp));
        end
        Spaths(:,6) = S;

        % (vii) Fixed counts (those of the path from (iv)), random order of the jumps, stratified spell lengths
        S = Spaths(:,4);
        counts = accumarray([S(1:T-1),S(2:T)], 1, [N_S,N_S]);
        firststate = S(1);
        laststate = S(T);
        timein = accumarray(S, 1, [N_S,1]);   % number of periods in each state
        jumps = counts-diag(diag(counts));    % number of jumps from s to s'~=s
        njump = sum(jumps,2);
        % Final jump from each state (other than laststate): random spanning tree pointing to laststate (Wilson's algorithm)
        intree = false(N_S,1);
        intree(laststate) = true;
        finalnext = zeros(N_S,1);
        for ii = 1:N_S
            if njump(ii)>0
                ss = ii;
                while ~intree(ss) % random walk, using the jump counts as weights, until we hit the tree
                    [~,finalnext(ss)] = max(cumsum(jumps(ss,:))/njump(ss)>rand);
                    ss = finalnext(ss);
                end
                ss = ii;
                while ~intree(ss) % add the (loop-erased) walk to the tree
                    intree(ss) = true;
                    ss = finalnext(ss);
                end
            end
        end
        % Jumps from each state: the others in random order, then the final one
        departures = cell(N_S,1);
        for ii = 1:N_S
            deps = repelem(1:N_S, jumps(ii,:));
            if ii~=laststate && njump(ii)>0
                deps(find(deps==finalnext(ii),1)) = [];
                deps = [deps(randperm(length(deps))), finalnext(ii)];
            else
                deps = deps(randperm(length(deps)));
            end
            departures{ii} = deps;
        end
        % The state of each spell
        spellstate = zeros(sum(njump)+1,1);
        spellstate(1) = firststate;
        visits = zeros(N_S,1);
        for kk = 2:length(spellstate)
            visits(spellstate(kk-1)) = visits(spellstate(kk-1))+1;
            spellstate(kk) = departures{spellstate(kk-1)}(visits(spellstate(kk-1)));
        end
        % The length of each spell: stratified draws from a geometric distribution, adjusted to add up to the time in the state
        spelllength = zeros(length(spellstate),1);
        for ii = 1:N_S
            K_s = sum(spellstate==ii);
            if K_s>0
                qstay = (timein(ii)-K_s)/timein(ii);    % probability of staying implied by the counts
                ustrat = ((1:K_s)'-1+rand(K_s,1))/K_s;  % stratified: one uniform draw in each of K_s equal intervals
                if qstay>0
                    ell = max(1, ceil(log(1-ustrat)/log(qstay))); % inverse-cdf of the geometric distribution
                else
                    ell = ones(K_s,1);
                end
                while sum(ell)<timein(ii) % add one period to a randomly chosen spell
                    kk = randi(K_s);
                    ell(kk) = ell(kk)+1;
                end
                while sum(ell)>timein(ii) % remove one period from a randomly chosen spell (that is longer than one period)
                    longer = find(ell>1);
                    kk = longer(randi(length(longer)));
                    ell(kk) = ell(kk)-1;
                end
                spelllength(spellstate==ii) = ell(randperm(K_s));
            end
        end
        S = repelem(spellstate, spelllength);
        Spaths(:,7) = S;

        % Measures of how representative each path is
        for k = 1:nmethods
            S = Spaths(:,k);
            counts = accumarray([S(1:T-1),S(2:T)], 1, [N_S,N_S]);
            Phat = counts./max(sum(counts,2),1);
            transerr(m,jT,k) = max(max(abs(Phat-pi_S)));
            statdisthat = accumarray(S, 1, [N_S,1])/T;
            staterr(m,jT,k) = max(abs(statdisthat-statdist));
            counterr(m,jT,k) = max(max(abs(counts-(T-1)*statdist.*pi_S)));
            counts_half = accumarray([S(1:T/2-1),S(2:T/2)], 1, [N_S,N_S]);
            Phat_half = counts_half./max(sum(counts_half,2),1);
            halferr(m,jT,k) = max(max(abs(Phat_half-pi_S)));
            % Duration error: spell lengths vs the geometric distribution implied by pi_S (drop the first and last spell, as they are cut off)
            spellstart = [1; find(diff(S)~=0)+1];
            spelllength = diff([spellstart; T+1]);
            spellstate = S(spellstart);
            spelllength = spelllength(2:end-1);
            spellstate = spellstate(2:end-1);
            if isempty(spelllength)
                durerr(m,jT,k) = 1; % path (almost) never changes state
            else
                pss = pi_S(sub2ind([N_S,N_S], spellstate, spellstate));     % probability of staying, for each spell
                ells = 1:max(spelllength);
                durerr(m,jT,k) = max(abs(mean(spelllength<=ells,1)-mean(1-pss.^ells,1)));
            end
            est_ES(m,jT,k) = mean(S_grid(S));
        end
    end
end

%% Report mean and variance of each estimator, and the mean and worst of each measure
fprintf('\nRQMC results:\n')
for k = 1:nmethods
    fprintf('%s\n', methods{k})
    for jT = 1:2
        fprintf('   T=%d: mean of estimator of E[S] = %.6f, variance of estimator = %.3e\n', Tvals(jT), mean(est_ES(:,jT,k)), var(est_ES(:,jT,k)))
        fprintf('          transition error: mean = %.4f, worst = %.4f\n', mean(transerr(:,jT,k)), max(transerr(:,jT,k)))
        fprintf('          stationary error: mean = %.4f, worst = %.4f\n', mean(staterr(:,jT,k)), max(staterr(:,jT,k)))
        fprintf('          count error:      mean = %.2f, worst = %.2f\n', mean(counterr(:,jT,k)), max(counterr(:,jT,k)))
        fprintf('          first-half transition error: mean = %.4f, worst = %.4f\n', mean(halferr(:,jT,k)), max(halferr(:,jT,k)))
        fprintf('          duration error:   mean = %.4f, worst = %.4f\n', mean(durerr(:,jT,k)), max(durerr(:,jT,k)))
    end
end

%% Histograms of the M estimates of E[S] (figure 3)
% Rows are the seven methods, columns are T=T1 and T=T2.
figure(3)
for k = 1:nmethods
    for jT = 1:2
        subplot(nmethods, 2, (k-1)*2 + jT)
        histogram(est_ES(:,jT,k), 30)
        xline(trueES, 'r--', 'LineWidth', 1.5)
        title(sprintf('%s, T=%d  (var = %.2e)', methods{k}, Tvals(jT), var(est_ES(:,jT,k))))
        xlabel('estimate'); ylabel('count')
    end
end
