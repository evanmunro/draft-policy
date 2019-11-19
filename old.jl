#g is how many games have happened already
#N is how many games there are
#v1 is how many games won in the g games by Team 1
#v2 is how many games won in the g games by Team 2
function pi_d_bound(N,g,p_d,pV=10000,start=0)
    k = (N+1)/2
    v1s = range(start,min(g,k-1)-start,step=1)
    bound = zeros(length(v1s))
    for i in 1:length(v1s)
        v1 = v1s[i]
        v2 = g-v1
        d_wins = Binomial(N-g,0.5)
        p_diff = cdf(d_wins,k-v1) - cdf(d_wins,k-v1-1)
        bound[i] = pV - 1/(p_d*p_diff)
    end
    return v1s, bound
end

N=81
p_d = pdf(Normal(),1)
g=10
v1s, bound = pi_d_bound(N,g,p_d)
plot(v1s,bound,xlabel="v1",label=["g=10","g=1","g=20"])
v1s, bound = pi_d_bound(N,1,p_d)
plot!(v1s,bound,label=["g=10","g=1","g=20"])

v1s, bound = pi_d_bound(N,20,p_d)
plot!(v1s,bound,label=["g=10","g=1","g=20"])


savefig("bound_early.pdf")

v1s, bound = pi_d_bound(N,35,p_d,10000,5)
plot(v1s,bound,label=["g=35"])
savefig("bound_later.pdf")

gs,plimit, p_diff = prob_spread_bound(81,10)
plot(gs,plimit,xlabel="Game g",
    ylabel="Probability Change",label=["T1 Draft Change Prob Limit","T1 Change in Win Prob if Win"],title="Team 1 Wins First 10 Games, Loses the Rest")
plot!(gs,p_diff,label=["T1 Draft Change Prob Limit","T1 Change in Win Prob if Win "])
savefig("prob_bound.pdf")

win_effect = simulate_season(teams,matches,win_prob)
condense_effect= zeros(length(teams),G)
for team in teams
    g =1
    for m in 1:length(matches)
        match = matches[m]
        if !isnothing(findfirst(x->x==team,match))
            condense_effect[team,g] = win_effect[team,m]
            g=g+1
        end
    end
end
condense_effect

plot(1:size(condense_effect,2),transpose(condense_effect),xlabel="Game Number",ylabel="Change in Prize Probability if Win",
    label=["Weak","Strong","Strong"])

    
