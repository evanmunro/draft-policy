using Random, Distributions, MosekTools, JuMP

mutable struct GameNode
    remain::Int
    pwin::Array{Float64}
    plose::Array{Float64}
    score::Array{Int}
    y::Array{Float64}
    children
end

function addNode(score, remain)
    N = length(score)
    node = GameNode(remain, zeros(Float64, N), zeros(Float64, N), score, zeros(Float64, N), nothing)
    addChildren(node)
    return node
end

function initializeGame(N, G)
    return addNode(zeros(Int, N), G)
end

function addChildren(parent)
    remain = parent.remain
    score = parent.score
    N = length(score)
    outcomes = collect(Iterators.product(repeat([[0,1]], N)...))
    pwin = zeros(N)
    plose = zeros(N)
    if remain == 0
        children = nothing
        losers = findall(==(minimum(score)), score)
        winners = findall(==(maximum(score)), score)
        pwin = [a in winners for a in 1:N ]./length(winners)
        plose = [ a in losers for a in 1:N ]./length(losers)
    else
        children = []
        for i in 1:length(outcomes)
            node = addNode(score .+ outcomes[i], remain-1)
            push!(children, node)
            pwin .+= 1/length(outcomes)*node.pwin
            plose .+= 1/length(outcomes)*node.plose
        end
    end
    parent.children = children
    parent.pwin = pwin
    parent.plose = plose
end

function get_objective(node)
    N = length(node.score)
    obj = 0
    if node.remain == 0
        losers = findall(==(minimum(node.score)), node.score)
        obj += sum((node.plose[l] - node.y[l])^2 - (node.y[l] - 0)^2 for l in losers)
        obj += sum((node.y[i] - 0)^2 for i in 1:N)
        return obj
    else
        for child in node.children
            obj += get_objective(child)/length(node.children)
        end
        return obj
    end
end

#recursive TIC calculation
function calculate_TIC(node::GameNode, pi, stopped = false, yfixed = nothing)
    if stopped
        node.y = copy(yfixed)
    else
        node.y = copy(node.plose)
    end
    # stop recursion
    if node.remain==0
        return 0
    end
    if stopped == false
        winIC, loseIC = calc_incentive(node)
        if minimum(winIC*pi .- loseIC) < 0
            stopped = true
            yfixed = copy(node.plose)
        end
    end
    for child in node.children
        calculate_TIC(child, pi, stopped, yfixed)
    end
end

#check conditional probabilities of winning market prize and getting aid
# when winning vs losing this game
function calc_incentive(node::GameNode)
    N = length(node.score)
    pwin_if_win = zeros(N); pwin_if_lose = zeros(N)
    plose_if_lose = zeros(N); plose_if_win = zeros(N)
    outcomes = collect(Iterators.product(repeat([[0,1]], N)...))
    for (out, child) in zip(outcomes, node.children)
        winners = collect(out .== 1)
        losers = collect(out .== 0)
        pwin_if_win[winners] .+=  child.pwin[winners]/(length(outcomes)/2)
        pwin_if_lose[losers] .+= child.pwin[losers]/(length(outcomes)/2)
        plose_if_lose[losers] .+= child.plose[losers]/(length(outcomes)/2)
        plose_if_win[winners] .+= child.plose[winners]/(length(outcomes)/2)
    end
    win_incentive = pwin_if_win - pwin_if_lose
    lose_incentive = plose_if_lose - plose_if_win
    return win_incentive, lose_incentive
end

function add_IC(node::GameNode, model, g, j, pi)
    N = length(node.score)
    O = 2^N
    outcomes = collect(Iterators.product(repeat([[0,1]], N)...))
    if node.remain == 0
        return 0
    else
        winIC, loseIC = calc_incentive(node)
        for i in 1:N
            record = collect(Iterators.flatten([out[i] for out in outcomes]))
            signs = ones(length(outcomes))
            signs[record.== 1] .= -1
            x = model[:x]
            @constraint(model, sum(1/O*x[g+1, k+O*(j-1), i]*signs[k] for k in 1:O) <= winIC[i]*pi)
        end

        for (child, jnew) in zip(node.children, (O*(j-1)+1):O*j)
            add_IC(child, model, g+1, jnew, pi)
        end
    end
end

function optimalRule(tree::GameNode, pi)
    N = length(tree.score)
    G = tree.remain
    losers = findLosers(tree)
    model = Model(Mosek.Optimizer)
    O = 2^N
    #draft rule for all partial histories
    @variable(model, 0<=x[i=0:G, j=1:O^i, k=1:N]<=1)

    #x = model[:x]
    #FAIR constraint
    for k in 1:N
        @constraint(model, x[0, 1, k]==1/N)
    end

    #PROB constraint
    for i in 1:G
        for j in 1:O^i
            @constraint(model, sum(x[i, j, k] for k in 1:N)==1)
        end
    end

    #DC Constraint under equal ability
    for i in 0:(G-1)
        for j in 1:O^i
            for k in 1:N
                @constraint(model, x[i, j, k] == 1/O*sum(x[i+1, s, k] for s in (O*(j-1)+1):O*j))
            end
        end
    end

    add_IC(tree, model, 0, 1, pi)

    @objective(model, Min, sum( sum((1/length(losers[j])-x[G, j, l])^2 - (0-x[G, j, l])^2 for l in losers[j]) +
                                sum( (0-x[G, j, k])^2 for k in 1:N)
                                for j in 1:O^G)/O^G)
    optimize!(model)

    return objective_value(model)
end

function findLosers(node)
    loser_record = []
    loserPass(node, loser_record)
    return loser_record
end

function loserPass(node, loser_record)
    if node.remain == 0
        losers = findall(==(minimum(node.score)), node.score)
        push!(loser_record, losers)
    else
        for child in node.children
            loserPass(child, loser_record)
        end
    end
end
