using Plots, Suppressor
include("draftpolicy.jl")
ENV["GKSwstype"]=100
#N=3 and G=5 is feasible ~10 minutes

N = 3
G = 4
abilities = zeros(N).+1/N
mult = 10
obj_opt = zeros(G)
obj_ntd = zeros(G)

#a,b, opt_wts, ntd_wts = calc_objective(s,5)
for g in 1:G
    s = Season(N, G, abilities)
    @suppress_out begin
        obj_opt[g], obj_ntd[g], a, b = calc_objective(s, mult)
    end
end
println(obj_opt .- obj_ntd)
obj_unif= zeros(G) .+ (0.5*abs(1- 1/(N-1)) + 0.5*(N-2)*abs(0-1/(N-1)))
obj_ante= zeros(G) .+ (0.5*abs(1- 1/(N)) + 0.5*(N-1)*abs(0-1/(N)))
plot(1:G, [ obj_ante obj_unif obj_opt obj_ntd], ylims=[0,1.1], xlabel="G", ylabel="Objective",label=["Ex-Ante Uniform" "Ex-Post Uniform" "Optimal" "T-IC"])
savefig("figures/gain_g.pdf")
