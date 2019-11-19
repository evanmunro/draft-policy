using Random, Combinatorics, Distributions, StatsBase

function game_prob(a,sigma=1)
    d = Normal(0,sigma)
    N = length(a)
    p = zeros(N,N)
    for i in 1:N
        for j in 1:N
            p[i,j] = 1-cdf(d,(a[j]-a[i])/(2*sigma))
        end
    end
    return p
end

function win_prob_ex(existing_wins,teams,matches,win_prob)
    all_outcomes = collect(Iterators.product(matches...))
    victory_prob = zeros(length(teams))
    for out =1:length(all_outcomes)
        prob = 1
        outcome = all_outcomes[out]
        wins=[count(x->x==i,outcome) for i in 1:length(teams)] + existing_wins
        for m =1:length(outcome)
            match = matches[m]
            winner = outcome[m]
            loser = match[findall(x->x!=winner,match)[1]]
            prob = prob*win_prob[winner,loser]
        end
        win_record = maximum(wins)
        winners = findall(x->x==win_record,wins)
        for j=1:length(winners)
            victory_prob[winners[j]] += prob/length(winners)
        end
    end
    return victory_prob
end

function win_prob_sim(existing_wins,teams,matches,win_prob)
    return 0
end

function add_win(wins,team)
    new_wins = copy(wins)
    new_wins[team] +=1
    return new_wins
end

#outcome_probs = get_prob_prize([0,0,0],teams,matches,win_prob)

function simulate_draft_rule(N,G,abilities,multiple=1,seed=0)
    #Random.seed!(seed)
    teams = Array(1:N)
    matches = repeat(collect(combinations(teams,2)),G)
    M = length(matches)
    win_prob = game_prob(abilities)
    ys = zeros(N,M+1) .+ 1/N
    cwins = zeros(N)
    wins = zeros(N,M)

    for m in 1:M
        match = matches[m]
        t1 = match[1]
        t2 = match[2]
        remaining_matches = matches[(m+1):length(matches)]
        probw = win_prob_ex(add_win(cwins,t1),teams,remaining_matches,win_prob)
        probl = win_prob_ex(add_win(cwins,t2),teams,remaining_matches,win_prob)
        diff1 = multiple*(probw[t1]-probl[t1])
        diff2 = multiple*(probl[t2]-probw[t2])

        #added ys[t1],ys[t2] in earlier version
        min_diff = min(diff1,diff2,2*ys[t1,m],2*ys[t2,m])
        boundl1 = min_diff/2+ys[t1,m]
        y1L = min(boundl1,1)
        y2L = min(min_diff/2 + ys[t2,m],1)
        y1W = max(2*ys[t1,m] - y1L,0)
        y2W = max(2*ys[t2,m] -y2L,0)
        weights = ProbabilityWeights([win_prob[match...], 1-win_prob[match...]])
        g_winner = sample(match,weights)
        cwins[g_winner] += 1
        wins[g_winner,m] = 1
        if g_winner==t1
            ys[t1,(m+1):(M+1)] .= y1W
            ys[t2,(m+1):(M+1)] .= y2L
        else
            ys[t1,(m+1):(M+1)] .= y1L
            ys[t2,(m+1):(M+1)] .= y2W
        end

    end
    println("Final standings in the season: ")
    println(cwins)
    println("Final Draft Probs: ")
    println(ys[:,M+1])
    return ys, wins, cwins
end

function simulate_draft_rule_N2(G,multiple=1,seed=0)
    Random.seed!(seed)
    y1s = zeros(G+1) .+ 1/2
    wins = zeros(G)
    k = (G+1)/2
    for g in 1:G
        d_wins = Binomial(Int(G-g),0.5)
        v1 = sum(wins)
        v2 = g-v1-1
        diff1 = multiple*(cdf(d_wins,k-v1-1) - cdf(d_wins,k-v1-2))
        diff2 = multiple*(cdf(d_wins,k-v2-1) - cdf(d_wins,k-v2-2))
        min_diff = min(diff1,diff2)
        boundl1 = min_diff/2+y1s[g]
        y1L = min(boundl1,1)
        y1W = max(2*y1s[g] - y1L,0)
        wins[g] = sample([0,1],ProbabilityWeights([0.5,0.5]))
        if wins[g]==1
            y1s[g+1] = y1W
        else
            y1s[g+1] = y1L
        end
    end
    return wins, cumsum(wins)./Array(1:G), y1s

end
