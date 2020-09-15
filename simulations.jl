using Revise
using Plots, CSV, DataFrames
includet("draftpolicy.jl")

16+21+23

N=3
G= 20
s = Season(N,G)
size(s.matches)

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
#plot(1:size(record.cwins,1),record.cwins,linestyle=:dash,color=[:blue :red :green], legend=:top,
  #  title="Win Records", 
   # label=["T1","T2","T3"],xlabel="Match",ylabel="Cumulative Wins")
savefig("sim3_abilities.pdf")




#



N=3
G=3
abilities = zeros(N).+ 1/N 
s = Season(N,G,abilities)

stops = zeros(20)
obj = zeros(20)

for m in 1:20
    o,st = expected_objective(s,m)
    obj[m] = o
    stops[m] = st 
    println("done ", m)
end 

obj[5]

serialize("results/draft_obj.jb",obj)

plot(1:20,obj,color=[:blue],title="Expected Gain over Uniform Lottery",xlabel="B",ylabel="Probability",label=[""])
#savefig("sim_obj.pdf")
##plot(1:20,stops/size(s.matches,1),color=[:red],title="Games Taken Into Account",label=[""],xlabel="B",ylabel="Percent")
#savefig("sim_stops.pdf")

using Serialization
record = deserialize("results/run1/object1985.jb")

year = 1985
stoptime = record.stop
team_id = CSV.read(string("/Users/evanmunro/Documents/Github/draft-policy/data/cleaned/teams_",
                            string(year),".csv"))
end_wins = record.cwins[:,size(record.cwins,2)]
end_draft = record.draft_prob[:,size(record.draft_prob,2)]
df = DataFrame(id = copy(team_id[:id]),team = copy(team_id[:team]),wins = copy(end_wins),score = copy(end_draft))
sort!(df, (:wins), rev=(false))
bottom_id = df[:id][1:3]
bottom_name = df[:team][1:3]
prob = record.draft_prob[bottom_id,:]
wins = record.cwins[bottom_id,:]
plot(transpose(prob[:,1:stoptime]),label=bottom_name,
    size=(800,250),
        linewidth=4,legend=:topright,xlabel="Game",title="Evolution of Draft Probability")

    plot(transpose(wins[:,1:stoptime]),label=vcat(bottom_name,["Stopping Time"]),
        size=(800,250),
        linewidth=4,legend=:topleft,title="Cumulative Wins")
    vline!([stoptime],color="black",label=vcat(bottom_name,["Stopping Time"]))

function simulate_season(year,mult=10)
    team_id = CSV.read(string("/Users/evanmunro/Documents/Github/draft-policy/data/cleaned/teams_",
                            string(year),".csv"))
    season = CSV.read(string("/Users/evanmunro/Documents/Github/draft-policy/data/cleaned/season_",year,".csv")) 
    matches = convert(Matrix, season[[:hID,:aID]])
    wins = convert(Array, season[:win])
    teams = convert(Array,team_id[:id])
    s = Season(teams,matches)
    record,_=calculate_draft_rule(s,wins,mult,16)
    stoptime=record.stop
    println(record.stop)
    end_wins = record.cwins[:,size(record.cwins,2)]
    end_draft = record.draft_prob[:,size(record.draft_prob,2)] 

    df = DataFrame(team = copy(team_id[:team]),wins = copy(end_wins),score = copy(end_draft))
    sort!(df, (:wins), rev=(true))

    added_amt = sum(df[1:16,:score])/(length(teams)-16)
    df[17:length(teams),:score] = df[17:length(teams),:score] .+ added_amt 
    df[1:16,:score] =  0 
    df[17:length(teams),:]
end

year = 2015 
simulate_season(year)

#1985,1987,1988,1989
year = 1985 
team_id = CSV.read(string("/Users/evanmunro/Documents/Github/draft-policy/data/cleaned/teams_",
                            string(year),".csv"))
season = CSV.read(string("/Users/evanmunro/Documents/Github/draft-policy/data/cleaned/season_",year,".csv")) 
team_id

sort!(team_id,:id,rev=true)

function simulate_season(year,mult=10)
    team_id = CSV.read(string("/Users/evanmunro/Documents/Github/draft-policy/data/cleaned/teams_",
                            string(year),".csv"))
    season = CSV.read(string("/Users/evanmunro/Documents/Github/draft-policy/data/cleaned/season_",year,".csv"))
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

simulate_season(1986)

matches = convert(Matrix, season[[:hID,:aID]])
wins = convert(Array, season[:win])
teams = convert(Array,team_id[:id])
s = Season(teams,matches)
#0.5 had good results 
#wins
record=calculate_draft_rule(s,wins,1.05,16)
stoptime=record.stop
println(record.stop)
end_wins = record.cwins[:,size(record.cwins,2)]
end_draft = record.draft_prob[:,size(record.draft_prob,2)] 

df = DataFrame(team = copy(team_id[:team]),wins = copy(end_wins),score = copy(end_draft))
sort!(df, (:wins), rev=(true))

added_amt = sum(df[1:16,:score])/(length(teams)-16)
df[17:length(teams),:score] = df[17:length(teams),:score] .+ added_amt 
df[1:16,:score] =  0 
df[17:length(teams),:]
# reduce at 69 when multiple=5

bottom3 = [16,8,10]
prob = record.draft_prob[bottom3,:]
wins = record.cwins[bottom3,:]
plot(transpose(prob[:,1:250]).+added_amt,label=["Knicks","Warriors","Spurs"],
    size=(800,250),
    #ylims=(0.05,0.3),
    linewidth=4,legend=:topright,xlabel="Game",title="Evolution of Draft Probability")
savefig("prob85.pdf")

plot(transpose(wins[:,1:250]),label=["Knicks","Warriors","Spurs","Stopping Time"],
    size=(800,250),
    linewidth=4,legend=:topleft,title="Cumulative Wins")
vline!([stoptime],color="black",label=["Knicks","Warriors","Spurs","Stopping Time"])
savefig("wins85.pdf")

year = 1989
team_id = CSV.read(string("/Users/evanmunro/Documents/Github/draft-policy/data/cleaned/teams_",
                            string(year),".csv"))
season = CSV.read(string("/Users/evanmunro/Documents/Github/draft-policy/data/cleaned/season_",year,".csv")) 
team_id

matches = convert(Matrix, season[[:hID,:aID]])
wins = convert(Array, season[:win])
teams = convert(Array,team_id[:id])
s = Season(teams,matches)
#0.5 had good results 
#wins
record=calculate_draft_rule(s,wins,15,16)
println(record.stop)
end_wins = record.cwins[:,size(record.cwins,2)]
end_draft = record.draft_prob[:,size(record.draft_prob,2)] 
df = DataFrame(team = copy(team_id[:team]),wins = copy(end_wins),score = copy(end_draft))
sort!(df, (:wins), rev=(true))

added_amt = sum(df[1:16,:score])/(length(teams)-16)
df[17:length(teams),:score] = df[17:length(teams),:score] .+ added_amt 
df[1:16,:score] =  0 
df[17:length(teams),:]

bottom3 = [14,3,22]
prob = record.draft_prob[bottom3,:]
wins = record.cwins[bottom3,:]
plot(transpose(prob[:,1:500]).+added_amt,label=["Heat","Hornets","Spurs"],
    size=(800,250), linewidth=4,
    legend=:topleft,xlabel="Game",title="Evolution of Draft Probability")
savefig("prob89.pdf")

plot(transpose(wins[:,1:500]),label=["Heat","Hornets","Spurs","Stopping Time"],
    linewidth=4,size=(800,250),
    legend=:topleft,title="Cumulative Wins")
vline!([333],color="black",label=["Heat","Hornets","Spurs","Stopping Time"])
savefig("wins89.pdf")

year = 1987
team_id = CSV.read(string("/Users/evanmunro/Documents/Github/draft-policy/data/cleaned/teams_",
                            string(year),".csv"))
season = CSV.read(string("/Users/evanmunro/Documents/Github/draft-policy/data/cleaned/season_",year,".csv")) 
team_id

matches = convert(Matrix, season[[:hID,:aID]])
wins = convert(Array, season[:win])
teams = convert(Array,team_id[:id])
s = Season(teams,matches)
#0.5 had good results 
#wins
record=calculate_draft_rule(s,wins,5,16)
println(record.stop)
end_wins = record.cwins[:,size(record.cwins,2)]
end_draft = record.draft_prob[:,size(record.draft_prob,2)] 
df = DataFrame(team = copy(team_id[:team]),wins = copy(end_wins),score = copy(end_draft))
sort!(df, (:wins), rev=(true))

added_amt = sum(df[1:16,:score])/(length(teams)-16)
df[17:length(teams),:score] = df[17:length(teams),:score] .+ added_amt 
df[1:16,:score] =  0 
df[17:length(teams),:]

bottom3 = [14,3,22]
prob = record.draft_prob[bottom3,:]
wins = record.cwins[bottom3,:]
plot(transpose(prob[:,1:500]).+added_amt,label=["Heat","Hornets","Spurs"],
    size=(800,250), linewidth=4,
    legend=:topleft,xlabel="Game",title="Evolution of Draft Probability")
savefig("prob87.pdf")

plot(transpose(wins[:,1:500]),label=["Heat","Hornets","Spurs","Stopping Time"],
    linewidth=4,size=(800,250),
    legend=:topleft,title="Cumulative Wins")
vline!([333],color="black",label=["Heat","Hornets","Spurs","Stopping Time"])
savefig("wins87.pdf")
