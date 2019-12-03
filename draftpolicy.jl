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

#to make parallel
#export JULIA_NUM_THREADS=4

function win_prob_ex(existing_wins,teams,matches,win_prob,top_k=1)
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

function win_prob_sim(existing_wins,teams,matches,p,top_k,sims=1000)
    win_record = zeros(length(teams))
    bottom_record = zeros(length(teams))
    M = size(matches,1)
    for s in 1:sims
        wins = zeros(length(teams))
        wins += existing_wins
        for m in 1:M
            match =matches[m,:]
            t1,t2 = match
            weights = ProbabilityWeights([p[t1,t2], 1-p[t1,t2]])
            gwin = sample(match,weights)
            wins[gwin]+=1
        end
        cutoff = sort(wins,rev=true)[top_k]
        bottom = sort(wins)[1]
        losers = (wins .<= bottom)
        winners = (wins .>= cutoff)
        win_record = win_record .+ winners/(sum(winners)/top_k)
        bottom_record = bottom_record .+ losers/(sum(losers))
    end
    win_record = win_record./sims
    loss_record = bottom_record ./sims
    #re-normalize in case of ties
    return win_record, loss_record
end

function add_win(wins,team)
    new_wins = copy(wins)
    new_wins[team] +=1
    return new_wins
end

function calculate_draft_rule(teams,matches,result,abilities,top_k=1,multiple=1,seed=0)
    #Random.seed!(seed)
    M = size(matches,1)
    N = length(teams)
    win_prob = game_prob(abilities)
    ys = zeros(N,M+1) .+ 1/N
    cwins = zeros(N)
    wins = zeros(N,M)
    losses = zeros(N,M)
    adjust = true
    for m in 1:M
        match = matches[m,:]
        t1 = match[1]
        t2 = match[2]
        remaining_matches = matches[(m+1):M,:]
        probw,problw = win_prob_sim(add_win(cwins,t1),teams,remaining_matches,win_prob,top_k)
        probl,probll = win_prob_sim(add_win(cwins,t2),teams,remaining_matches,win_prob,top_k)
        diff1 = max(multiple*(probw[t1]-probl[t1]),0)
        diff2 = max(multiple*(probl[t2]-probw[t2]),0)
        cwins[result[m]] += 1
        wins[result[m],m] = 1
        if adjust && min(diff1,diff2) >= max(probll[t2]-problw[t2],probll[t1]-problw[t1])
            if result[m]==t1
                ys[:,(m+1):(M+1)] .= problw
                losses[t2,m] = 1
            else
                ys[:,(m+1):(M+1)] .= probll
                losses[t1,m] = 1
            end

        else
            println("stopped adjusting: ",m)
            adjust = false
        end
    end
    println("Final standings in the season: ")
    println(cwins)
    return ys, wins, cwins, losses
end


function simulate_draft_rule(N,G,abilities,multiple=1,top_k=1,seed=0)
    #Random.seed!(seed)
    teams = Array(1:N)
    team_combos = collect(Iterators.flatten(combinations(teams,2)))
    team_combos = transpose(reshape(team_combos,(2,Integer(length(team_combos)/2))))
    matches = repeat(team_combos,G)
    M = size(matches,1)
    println(M)
    win_prob = game_prob(abilities)
    wp0, lp0 = win_prob_sim(zeros(N),teams,matches,win_prob,top_k)
    ys = zeros(N,M+1) .+ lp0
    cwins = zeros(N)
    wins = zeros(N,M)
    losses = zeros(N,M)
    adjust = true

    for m in 1:M
        match = matches[m,:]
        t1 = match[1]
        t2 = match[2]
        remaining_matches = matches[(m+1):M,:]
        probw,problw = win_prob_sim(add_win(cwins,t1),teams,remaining_matches,win_prob,top_k)
        probl,probll = win_prob_sim(add_win(cwins,t2),teams,remaining_matches,win_prob,top_k)
        diff1 = max(multiple*(probw[t1]-probl[t1]),0)
        diff2 = max(multiple*(probl[t2]-probw[t2]),0)
        weights = ProbabilityWeights([win_prob[match...], 1-win_prob[match...]])
        gwinner = sample(match,weights)
        cwins[gwinner] += 1
        wins[gwinner,m] = 1
        if adjust && min(diff1,diff2) >= max(probll[t2]-problw[t2],probll[t1]-problw[t1])
            if gwinner==t1
                ys[:,(m+1):(M+1)] .= problw
                losses[t2,m] = 1
            else
                ys[:,(m+1):(M+1)] .= probll
                losses[t1,m] = 1
            end

        else
            println("stopped adjusting: ",m)
            adjust= false
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
