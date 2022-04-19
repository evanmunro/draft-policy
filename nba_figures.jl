using Plots, CSV, DataFrames, Serialization
include("draftpolicy.jl")

CSV.read(string("results/run2/record",1987,".csv"))

record = deserialize(string("results/run2/object",1987,".jb"))
bottom3 = [11,14,15]
prob = record.draft_prob[bottom3,:]
wins = record.cwins[bottom3,:]
stoptime=record.stop
plot(transpose(wins[:,1:(stoptime+50)]),label=["Clippers","Nets","Knicks","Stopping Time"],
    size=(800,250),
    linewidth=4,legend=:topleft,ylabel="Cumulative Wins")
vline!([stoptime],color="black",label=["Clippers","Nets","Knicks","Stopping Time"])
savefig("wins87.pdf")

plot(transpose(prob[:,1:(stoptime+50)]),label=["Clippers","Nets","Knicks","Stopping Time"],
    size=(800,250), ylabel="Draft Probability", 
    #ylims=(0.05,0.3),
    linewidth=4,legend=:topleft,xlabel="Game")
vline!([stoptime],color="black",label=["Clippers","Nets","Knicks","Stopping Time"])
savefig("prob87.pdf")

years = [1985,1986,1987,1988,1989]

results = DataFrame(ranking=1:8)

results
for i in 1:length(years)
    yr = years[i] 
    record = deserialize(string("results/run2/object",yr,".jb"))
    standings = CSV.read(string("results/run2/record",yr,".csv"))
    results[Symbol(string("score_",yr))] = standings[:score][1:8].*100
    results[Symbol(string("score_",yr))][8] = sum(standings[:score][8:length(standings[:score])])*100
    println(size(record.draft_prob))
    println(record.stop)
end



((527/944)+435/944+371/944+465/944+353/1026)/5

58.96+5.1+6.5+17.4+7.4+0.92+0.695

results

sum(results[1,:][2:length(results[1,:])])/5


