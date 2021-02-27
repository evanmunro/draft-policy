###############################################################################
# Functions for calculating draft rules and estimating abilities
# for simulated and actual sports seasons
# Author: Evan Munro
# Date: January, 2020
###############################################################################

using Random, Combinatorics, Distributions, StatsBase, JuMP, Ipopt, MosekTools, SpecialFunctions, Suppressor


############################################################################
# Functions for calculating abilities
#############################################################################

#calculate updated abilities depending on if you win or lose the match
function contingent_abilities(matches,winners,losers,N,lambda)
    m = length(winners)
    t1,t2 = matches[m,:]
    winners[m] = t1
    losers[m] = t2
    ahatw = estimate_alpha(matches,winners,losers,N,lambda)
    winners[m] = t2
    losers[m] = t1
    ahatl = estimate_alpha(matches,winners,losers,N,lambda)
    winners[m] = 0
    losers[m] = 0
    return ahatw, ahatl
end

#MLE estimation of abilities based on records
function estimate_alpha(matches,wins,loses,N,lambda=10)
    model=Model(Ipopt.Optimizer)
    @variable(model,a[1:N],start=1/N)
    ncdf(x) = (1+erf(x/sqrt(2)))/2
    register(model,:ncdf,1,ncdf,autodiff=true)
    @constraint(model, sum(a)==1)
    @constraint(model, a .>=0)
    @NLobjective(model,Max,sum(log(1-ncdf((a[j]-a[i])/sqrt(2))) for (i,j) in zip(wins,loses)) -lambda*sum((a[i] - 1/N)^2 for i in 1:N))
    @suppress_out begin
        status = optimize!(model);
    end
    return value.(a)
end

#to make parallel
#export JULIA_NUM_THREADS=4


############################################################################
# Functions for calculating expected values of league objective
#############################################################################


function expected_objective_sim(s,abilities,calc_ability=false,mult=10,sims=1000)
    exp_obj = 0.0
    exp_stop = 0.0
    p = game_prob(abilities)
    M = size(s.matches,1)
    for i in 1:sims
        winners = zeros(Int64,M)
        for m in 1:M
            match = s.matches[m,:]
            t1,t2 = match
            weights = ProbabilityWeights([p[t1,t2], 1-p[t1,t2]])
            gwin = sample(match,weights)
            winners[m]=gwin
        end
        record = calculate_draft_rule(s, winners, calc_ability, abilities,mult)
        exp_obj += objective(record)
        exp_stop += record.stop
    end
    return exp_obj/sims, exp_stop/sims
end
function expected_objective(s,mult=10)

    M = size(s.matches,1)
    match_list = fill(Int[], M)
    for m in 1:M
        match_list[m] = s.matches[m,:]
    end
    all_outcomes = collect(Iterators.product(match_list...))
    exp_obj = 0
    exp_stop = 0
    p_total = 0
    for out=1:length(all_outcomes)
        winners = [all_outcomes[out]...]
        p_outcome = prob_winners(s, winners, s.abilities)
        record = calculate_draft_rule(s, winners, mult)
        exp_obj += p_outcome* objective(record)
        exp_stop += p_outcome*record.stop
        p_total += p_outcome
    end
    return exp_obj, exp_stop#, p_total
end

############################################################################
# Functions for calculating win and loss probabilities in single games
# and in seasons
#############################################################################

#probability of a certain history
function prob_winners(s, winners, abilities)
    prob = 1
    pwins = game_prob(abilities)
    for m in 1:size(s.matches,1)
        match = s.matches[m,:]
        w = winners[m]
        l = match[match.!=w][1]
        prob = prob*pwins[w,l]
    end
    return prob
end

#exact calculation of win probability based on enumerating all possible
#results from the season
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

#approximate calculation of win probability based on simulation
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

function game_prob(a,sigma=1)
    d = Normal(0,1)
    N = length(a)
    p = zeros(N,N)
    for i in 1:N
        for j in 1:N
            p[i,j] = 1-cdf(d,(a[j]-a[i])/(sqrt(2)*sigma))
        end
    end
    return p
end

#############################################################################
# Structures to hold details and results for a season
#############################################################################

struct Season
    teams::Array{Int}
    matches::Array{Int,2}
    abilities::Array{Float64}
    pwin::Array{Float64,2}
end

function Season(N::Int, G::Int, abilities::Array{Float64})
    teams = Array(1:N)
    team_combos = collect(Iterators.flatten(combinations(teams, 2)))
    team_combos = transpose(reshape(team_combos, (2, Integer(length(team_combos)/2))))
    matches = repeat(team_combos, G)
    pwin = game_prob(abilities)
    return Season(teams, matches, abilities, pwin)
end

struct Record
    draft_prob::Array{Float64}
    stop::Int
    cwins::Array
    abilities::Array{Float64}
end

function summarize(r::Record)
    println("Final Draft probabilities: ")
    println(r.draft_prob[:,size(r.draft_prob,2)])
    println("Stopped at match ", r.stop)
    println("Final record: ")
    println(r.cwins[:,size(r.cwins,2)])
    println("Ability estimates: ")
    println(r.abilities[:,size(r.abilities,2)])
end

function objective(r::Record)
    d = draft(r)
    #adjust = d[r.cwins[:,size(r.cwins,2)].==maximum(r.cwins[:,size(r.cwins,2)])][1]
    obj = d[r.cwins[:,size(r.cwins,2)].==minimum(r.cwins[:,size(r.cwins,2)])][1]
    return obj #+ adjust/2
end

function draft(r::Record)
    return r.draft_prob[:,size(r.draft_prob,2)]
end

#function ability(r::Record)

function add_win(wins,team::Int)::Array{Int}
    new_wins = copy(wins)
    new_wins[team] +=1
    return new_wins
end

function get_losers(winners::Array, matches::Array{Int,2})::Array{Int}
    losers = zeros(Int,length(winners))
    for i in 1:length(winners)
        match = matches[i,:]
        w = winners[i]
        l = match[match.!=w][1]
    end
    return losers
end

#############################################################################
# Functions for calculating evolution of rule
##############################################################################

function calculate_draft_rule(s::Season, winners=nothing,
                              multiple=10, top_k = 1, seed=0)
    #Random.seed!(seed)
    #default is to not calculate ability updates
    adjust = true
    matches = s.matches
    teams = s.teams
    abilities = s.abilities
    M = size(matches,1)
    N = length(teams)

    if winners == nothing
        winners = zeros(Int64,M)
        losers = zeros(Int64,M)
    else
        losers = get_losers(winners,matches)
    end

    draft_prob = zeros(N,M+1) .+ 1/N
    cwins = zeros(Int64,N,M)

    stop_time = M
    for m in 1:M
        if m>1
            lag= m-1
        else
            lag = 1
        end
        match = matches[m,:]
        t1, t2 = match

        #if simulating season, then sample result
        if winners[m] == 0
            win_prob = game_prob(abilities)
            weights = ProbabilityWeights([win_prob[match...], 1-win_prob[match...]])
            w = sample(match,weights)
            winners[m] = w
        end

        if adjust
            probw1,probl1 = win_prob_sim(add_win(cwins[:,lag],t1),
                                        teams, matches[(m+1):M,:],
                                        s.pwin,top_k)
            probw2,probl2 = win_prob_sim(add_win(cwins[:,lag],t2),
                                        teams, matches[(m+1):M,:],
                                        s.pwin,top_k)

            diff1 = max(multiple*(probw1[t1]-probw2[t1]),0)
            diff2 = max(multiple*(probw2[t2]-probw1[t2]),0)

            if  (diff1 <  (probl2[t1]-probl1[t1])) || (diff2 < (probl1[t2]-probl2[t2]))
                adjust = false
                stop_time = m
            else
                if winners[m]==t1
                    draft_prob[:,(m+1):(M+1)] .= probl1
                else
                    draft_prob[:,(m+1):(M+1)] .= probl2
                end
            end
        end

        cwins[:,m] = cwins[:,lag]
        cwins[winners[m],m] += 1
        losers[m] = match[match.!=winners[m]][1]
    end
    return Record(draft_prob, stop_time, cwins, abilities)
end


##OLDWITHH ABILITIES
#estimate evolution of draft rule based on real data
function calculate_draft_rule(s::Season, winners, multiple,
                                top_k, calc_ability, abilities,
                                lambda=10, seed=0)
    #Random.seed!(seed)
    #default is to not calculate ability updates
    println(multiple)
    println(top_k)
    adjust = true
    matches = s.matches
    teams = s.teams
    M = size(matches,1)
    N = length(teams)
    ahat = zeros(N,M) .+ 1/N

    if abilities==nothing && winners == nothing
        error("To simulate a season, must provide abilities to simulate from")
    end

    if abilities == nothing
        abilities = zeros(N) .+ 1/N
    end

    if calc_ability==false
        ahat =  zeros(N,M) .+ abilities
        println("not calculating")
    end
    if winners == nothing
        winners = zeros(Int64,M)
        losers = zeros(Int64,M)
    else
        losers = get_losers(winners,matches)
    end

    draft_prob = zeros(N,M+1) .+ 1/N
    cwins = zeros(Int64,N,M)

    stop_time = M
    for m in 1:M
        if m>1
            lag= m-1
        else
            lag = 1
        end
        match = matches[m,:]
        t1, t2 = match
        if calc_ability
            ahat1,ahat2 = contingent_abilities(matches[1:m,:],
                                              copy(winners[1:m]),
                                              copy(losers[1:m]), N,lambda)
        else
            ahat1 = abilities
            ahat2 = abilities
        end
        #if simulating season, then sample result
        if winners[m] == 0
            win_prob = game_prob(abilities)
            weights = ProbabilityWeights([win_prob[match...], 1-win_prob[match...]])
            w = sample(match,weights)
            winners[m] = w
        end

        if adjust
            probw1,probl1 = win_prob_sim(add_win(cwins[:,lag],t1),
                                        teams, matches[(m+1):M,:],
                                        game_prob(ahat1),top_k)
            probw2,probl2 = win_prob_sim(add_win(cwins[:,lag],t2),
                                        teams, matches[(m+1):M,:],
                                        game_prob(ahat2),top_k)

            diff1 = max(multiple*(probw1[t1]-probw2[t1]),0)
            diff2 = max(multiple*(probw2[t2]-probw1[t2]),0)
            if t1==1
                incentives[m] = diff1/multiple
            elseif t2==1
                incentives[m] = diff2/multiple
            else
                incentives[m] = incentives[lag]
            end
            if  (diff1 < (probl2[t1]-probl1[t1])) || (diff2 <(probl1[t2]-probl2[t2]))
                adjust = false
                stop_time = m
            else
                if winners[m]==t1
                    draft_prob[:,(m+1):(M+1)] .= probl1
                else
                    draft_prob[:,(m+1):(M+1)] .= probl2
                end
            end
        end

        if calc_ability
            if winners[m]==t1
                ahat[:,m] = ahat1
            else
                ahat[:,m] = ahat2
            end
        end
        cwins[:,m] = cwins[:,lag]
        cwins[winners[m],m] += 1
        losers[m] = match[match.!=winners[m]][1]
    end
    return Record(draft_prob, stop_time, cwins, ahat)
end


#############################################################################
# Functions for calculating globally optimal rule
##############################################################################

function calc_objective(s, mult)
    println(mult)
    topk = 1
    N = length(s.teams)
    y = zeros(N) .+ 1/N
    M = size(s.matches, 1)
    ntd_val = zeros(2^M)
    ntd_wts = zeros(M)
    losers = ones(Int,2^M)
    adjust = true
    model = Model(Mosek.Optimizer)

    #draft rule for all partial histories
    @variable(model, 0<=x[i=0:M, j=1:2^i, k=1:N]<=1)

    #x = model[:x]
    #FAIR constraint
    for k in 1:N
        @constraint(model,x[0,1,k]==1/N)
    end

    #PROB constraint
    for i in 1:M
        for j in 1:2^i
            @constraint(model, sum(x[i,j,k] for k in 1:N)==1)
        end
    end

    #DC Constraint under equal ability
    for i in 0:(M-1)
        for j in 1:2^i
            for k in 1:N
                @constraint(model, x[i,j,k] == 1/2*x[i+1,2*j,k]+ 1/2*x[i+1,2*j-1,k])
            end
        end
    end

    history = zeros(N)
    step_forward(model, history, y, adjust, ntd_val, ntd_wts, losers,
                  s, 1, 1, mult, topk)

    @objective(model, Min, sum((1-x[M, j, losers[j]])^2 - (0-x[M, j, losers[j]])^2 +
                             sum((0-x[M, j, k])^2 for k in 1:N)  for j in 1:2^M)/2^M)
    #@NLobjective(model, Min, 1/2*sum(abs(1-x[M, j, losers[j]]) - abs(0-x[M, j, losers[j]]) +
    #                         sum(abs(0-x[M, j, k]) for k in 1:N)  for j in 1:2^M)/2^M)

    optimize!(model)
    println("here")
    opt_wts = zeros(M)
    for m in 1:M
        wt = 0
        for j in 1:2^(m-1)
            for n in 1:N
                wt += abs(value(x[m,2*j-1,n]) - value(x[m-1,j,n]))
                wt += abs(value(x[m,2*j,n]) -  value(x[m-1,j,n]))
            end
        end
        opt_wts[m] = wt/(N*2^m)
    end
    println("OBJECTIVE VALUE: ", objective_value(model))
    return objective_value(model),
            mean(ntd_val), opt_wts, ntd_wts

end

#adds NTD constraint and calculates draft rule as step forward through
#tree of all possible season outcomes
#only for equal abilities, for now
function step_forward(model, history, y, adjust, ntd_val, ntd_wts, losers, s, m, j, pi, topk)
    M = size(s.matches,1)
    t1,t2 = s.matches[m,:]
    probw1, probl1 = win_prob_sim(add_win(history,t1), s.teams,
                    s.matches[(m+1):M,:], s.pwin, topk)
    probw2, probl2 = win_prob_sim(add_win(history,t2), s.teams,
                    s.matches[(m+1):M,:], s.pwin, topk)

    diff1 = pi*(probw1[t1]-probw2[t1])
    diff2 = pi*(probw2[t2]-probw1[t2])

    x = model[:x]
    @constraint(model, x[m, 2*j, t1] - x[m, 2*j-1, t1] <= diff1)
    @constraint(model, x[m, 2*j-1, t2] - x[m, 2*j, t2] <= diff2)

    if !adjust || diff1 < (probl2[t1]-probl1[t1]) || diff2 < (probl1[t2]-probl2[t2])
        adjust = false
        y1 = copy(y)
        y2 = copy(y)
    else
        y1 = probl1
        y2 = probl2
    end

    ntd_wts[m] += (mean(abs.(y1.-y)) + mean(abs.(y2.-y)))/2^m

    if m < M
        step_forward(model, add_win(history, t1), y1, adjust,
                     ntd_val, ntd_wts, losers, s, m+1, 2*j-1, pi, topk)
        step_forward(model, add_win(history, t2), y2, adjust,
                     ntd_val, ntd_wts, losers, s, m+1, 2*j, pi, topk)
    else
        #which loser for the given history
        N = length(s.teams)
        l1 = argmin(add_win(history,t1))
        l2 = argmin(add_win(history,t2))
        losers[2*j-1] = l1
        losers[2*j] = l2
        ntd_val[2*j-1] = ((1-y1[l1])^2 - (0-y1[l1])^2 + sum( (0-y1[k])^2 for k in 1:N))
        ntd_val[2*j] = ((1-y2[l2])^2  - (0-y2[l2])^2 + sum((0-y2[k])^2 for k in 1:N))
    end
end
