includet("aidpolicy.jl")
using Plots
#6 (or maybe 7) is the limit on my laptop 
N=3
G=6
pi_t = 10
pi_b = 1
nsims = pi_t - pi_b+1
obj_opt = zeros(nsims)
obj_tic = zeros(nsims)

for (pi, i) in zip(pi_b:pi_t, 1:nsims)
    tree = initializeGame(N, G)
    calculate_TIC(tree, pi)
    obj_tic[i] = get_objective(tree)
    obj_opt[i] = optimalRule(tree, pi)
end

obj_unif = zeros(nsims) .+ (1- 1/(N-1))^2 .+ (N-2)*(0-1/(N-1))^2
obj_ante = zeros(nsims) .+ (1- 1/(N))^2 .+ (N-1)*(0-1/(N))^2
plot(pi_b:pi_t,[ obj_ante obj_unif obj_opt obj_tic], ylims=[0,1.1], xlabel="B", ylabel="Objective",label=["Ex-Ante Uniform" "Ex-Post Uniform" "Optimal" "T-IC"])
