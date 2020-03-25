#using Revise
using Plots, CSV, DataFrames
include("draftpolicy.jl")
path = "/scratch/users/munro/draft-policy/"

N= 3
G= 20
s = Season(N,G)
abilities = zeros(N) .+ 1/N
record,incentives = calculate_draft_rule(s,nothing,10,1,false,abilities)

plot(1:length(incentives),incentives,label=[""],xlabel="Match",ylabel="Bound on Draft Probability Change")
savefig("sim3_inc.pdf")

plot(1:(size(record.draft_prob,2)),transpose(record.draft_prob),xlabel="Match",ylabel="Percentage",
    color=[:blue :red :green], legend=:topright,label=["T1","T2","T3"])
savefig("sim3_prod.pdf")

plot(1:(size(record.cwins,2)),transpose(record.cwins),xlabel="Match",ylabel="Wins",
    color=[:blue :red :green], legend=:top,label=["T1","T2","T3"])
savefig("sim3_wins.pdf")

plot(1:(size(record.abilities,2)),transpose(record.abilities),xlabel="Match",ylabel="Abilities",
    color=[:blue :red :green], legend=:top,label=["T1","T2","T3"])
savefig("sim3_abilities.pdf")


function simulate_season(year,mult=10)
    team_id = CSV.read(string(path,"data/cleaned/teams_",
                            string(year),".csv"))
    season = CSV.read(string(path,"data/cleaned/season_",year,".csv"))
    matches = convert(Matrix, season[[:hID,:aID]])
    wins = convert(Array, season[:win])
    teams = convert(Array,team_id[:id])
    s = Season(teams,matches)
    record,_=calculate_draft_rule(s,wins,mult,16)
    stoptime=record.stop
    end_wins = record.cwins[:,size(record.cwins,2)]
    end_draft = record.draft_prob[:,size(record.draft_prob,2)]
    print("got here")
    df = DataFrame(id = copy(team_id[:id]),team = copy(team_id[:team]),wins = copy(end_wins),score = copy(end_draft))
    sort!(df, (:wins), rev=(false))
    bottom_id = df[:id][1:3]
    bottom_name = df[:team][1:3]
    prob = record.draft_prob[bottom_id,:]
    wins = record.cwins[bottom_id,:]
    plot(transpose(prob[:,1:stoptime]),label=bottom_name,
        size=(800,250),
        linewidth=4,legend=:topright,xlabel="Game",title="Evolution of Draft Probability")
    savefig(string("prob",string(year),".pdf"))

    plot(transpose(wins[:,1:stoptime]),label=vcat(bottom_name,["Stopping Time"]),
        size=(800,250),
        linewidth=4,legend=:topleft,title="Cumulative Wins")
    vline!([stoptime],color="black",label=vcat(bottom_name,["Stopping Time"]))
    savefig(string("wins",string(year),".pdf"))

    CSV.write(string("record",string(year),".csv"),df)
    df
end

for year in [1985,1986,1987,1988,1989]
    simulate_season(year)
end
