using Plots, Suppressor
includet("draftpolicy.jl")
ENV["GKSwstype"]=100
#N=3 and G=5 is feasible ~10 minutes

N = 3
G = 4
abilities = zeros(N).+1/N
s = Season(N,G,abilities)
pi_t = 10
pi_b = 1
nsims = pi_t - pi_b+1
obj_opt = zeros(nsims)
obj_ntd = zeros(nsims)

#a,b, opt_wts, ntd_wts = calc_objective(s,5)
for mult in pi_b:pi_t
    opt_wts = zeros(length(s.matches))
    ntd_wts= zeros(length(s.matches))
    i = mult-pi_b+1
    obj_opt[i], obj_ntd[i], opt_wts, ntd_wts = calc_objective(s, mult)
    if mult ==10
        println(obj_opt)
        println(obj_ntd)
        plot([ opt_wts ntd_wts],xlabel="Match",ylabel="Weight",ylims=[0,0.2],labels=["Optimal Weight" "T-IC Weight"])
        savefig("figures/weights10.pdf")
    end
end

obj_unif= zeros(nsims) .+ (1- 1/(N-1))^2 .+ (N-2)*(0-1/(N-1))^2
obj_ante= zeros(nsims) .+ (1- 1/(N))^2 .+ (N-1)*(0-1/(N))^2
plot(pi_b:pi_t,[ obj_ante obj_unif obj_opt obj_ntd], ylims=[0,1.1], xlabel="B", ylabel="Objective",label=["Ex-Ante Uniform" "Ex-Post Uniform" "Optimal" "T-IC"])
savefig("figures/simulated_gain.pdf")



obj_ante
