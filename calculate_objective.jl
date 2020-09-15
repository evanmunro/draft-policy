using Plots, Suppressor
include("draftpolicy.jl")

#N=3 and G=5 is feasible ~10 minutes

N = 3
G = 3
abilities = zeros(N).+1/N
s = Season(N,G,abilities)
pi_t = 10
pi_b = 0
nsims = pi_t - pi_b+1
obj_opt = zeros(nsims)
obj_ntd = zeros(nsims)
opt_wts = zeros(length(s.matches))
ntd_wts= zeros(length(s.matches))
#a,b, opt_wts, ntd_wts = calc_objective(s,5)
for mult in pi_b:pi_t
    @suppress_out begin
        i = mult-pi_b+1
        obj_opt[i], obj_ntd[i], opt_wts, ntd_wts = calc_objective(s, mult)
    end
    if mult ==10
        plot([ opt_wts ntd_wts],xlabel="Match",ylabel="Weight",ylims=[0,0.2],labels=["optimal weights","draft weights"])
        savefig("figures/weights10.pdf")
    end
end

obj_unif= zeros(nsims) .+ 1/(N-1)
obj_ante= zeros(nsims) .+ 1/(N)
plot(pi_b:pi_t,[ obj_ante obj_unif obj_opt obj_ntd], ylims=[0,1.1], xlabel="B", ylabel="Expected Objective",label=["Ex-Ante Uniform", "Ex-Post Uniform","Optimal","R-NTD"])
savefig("figures/simulated_gain.pdf")



obj_ante
