###############################################################################
# Functions for calculating draft rules and estimating abilities
# for simulated and actual sports seasons
# Author: Evan Munro
# Date: January, 2020
###############################################################################

using Random, Combinatorics, Distributions, StatsBase, JuMP, Ipopt, MosekTools, SpecialFunctions, Suppressor


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


#approximate calculation of win probability based on simulation
function win_prob_sim(existing_wins, teams, matches, p, top_k, sims=1000)
    win_record = zeros(length(teams))
    bottom_record = zeros(length(teams))
    M = size(matches,1)
    for s in 1:sims
        wins = zeros(length(teams))
        wins += existing_wins
        for m in 1:M
            match =matches[m,:]
            t1, t2 = match
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

function game_prob(a, sigma=1)
    N = length(a)
    p = zeros(N,N)
    for i in 1:N
        for j in 1:N
            p[i, j] = a[i]/(a[i] + a[j])
        end
    end
    return p
end

#############################################################################
# Structures to hold details and results for a season
#############################################################################

struct Season
    teams::Array{Int}
    matches::Array{Int, 2}
    pwin::Array{Float64,2}
end

struct Record
    draft_prob::Array{Float64}
    stop::Int
    cwins::Array
end

function summarize(r::Record)
    println("Final Draft probabilities: ")
    println(r.draft_prob[:,size(r.draft_prob,2)])
    println("Stopped at match ", r.stop)
    println("Final record: ")
    println(r.cwins[:,size(r.cwins,2)])
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
                              multiple=10, top_k = 1, seed=100)
    #Random.seed!(seed)
    adjust = true
    matches = s.matches
    teams = s.teams
    M = size(matches,1)
    N = length(teams)

    losers = get_losers(winners,matches)

    pw, pl = win_prob_sim(zeros(N),
                        teams, matches,
                        s.pwin, top_k)
    println(pl)
    draft_prob = zeros(N,M+1) .+ pl
    cwins = zeros(Int64, N, M)

    stop_time = M
    stop_prob = 100
    for m in 1:M
        if m>1
            lag= m-1
        else
            lag = 1
        end
        match = matches[m,:]
        t1, t2 = match

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
                println(stop_time)
                stopwins, _ = win_prob_sim(cwins[:, lag],
                                            teams, matches[m:M,:],
                                            s.pwin, top_k)
                stop_prob = minimum(stopwins)
                println(stopwins)
                println(stop_prob)
                println(draft_prob[11, m+1])

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

        #stopwins, _ = win_prob_sim(cwins[:, lag],
        #                            teams, matches[m:M,:],
        #                            s.pwin, top_k)
        #stop_prob = minimum(stopwins)
        #if stop_prob < 0.0001
        ##    println(stop_prob)
        #    println(m)
        #end
    end
    return Record(draft_prob, stop_time, cwins)
end
