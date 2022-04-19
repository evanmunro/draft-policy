#using Revise
using Plots, CSV, DataFrames, Serialization
include("draftpolicy.jl")
#path = "/scratch/users/munro/draft-policy/"
path  = ""

function simulate_season(year, mult=10)
    team_id = DataFrame(CSV.File(string(path,"data/cleaned/teams_",
                            string(year),".csv")))
    #ID 11 is Clippers
    season = DataFrame(CSV.File(string(path,"data/cleaned/season_",year,".csv")))
    matches = convert(Matrix, season[[:hID,:aID]])
    wins = convert(Array, season[:win])
    teams = convert(Array,team_id[:id])
    abilities = ones(length(teams)) .+ 1
    abilities[11] = 1.8
    abilities[14] = 1.8
    abilities[15] = 1.8

    s = Season(teams, matches, game_prob(abilities))
    record =calculate_draft_rule(s, wins, mult, 16)
    stoptime=record.stop
    println(year)
    println(stoptime)
    println(record.stop_prob)
    end_wins = record.cwins[:,size(record.cwins,2)]
    end_draft = record.draft_prob[:,size(record.draft_prob,2)]
    print("got here")
    df = DataFrame(id = copy(team_id[:id]),team = copy(team_id[:team]),wins = copy(end_wins),score = copy(end_draft))
    sort!(df, (:wins), rev=(false))
    bottom_id = df[:id][1:3]
    bottom_name = df[:team][1:3]
    prob = record.draft_prob[bottom_id,:]
    wins = record.cwins[bottom_id,:]
    serialize(string(path,"object",string(year),".jb"),record)
    CSV.write(string(path,"record",string(year),".csv"),df)
    df
end

#for year in [1985,1986,1987,1988,1989]
for year in [1987]
    simulate_season(year)
end
