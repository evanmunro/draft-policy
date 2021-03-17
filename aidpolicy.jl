using Random, Distributions, Mosek, JuMP

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


function optimalRule(tree::GameNode)
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
        @constraint(model,x[0, 1, k]==1/N)
    end

    #PROB constraint
    for i in 1:M
        for j in 1:O^i
            @constraint(model, sum(x[i, j, k] for k in 1:N)==1)
        end
    end

    #DC Constraint under equal ability
    for i in 0:(M-1)
        for j in 1:O^i
            for k in 1:N
                @constraint(model, x[i, j, k] == 1/O*sum(x[i+1, s, k] for s in (O*(j-1)+1):O*j))
            end
        end
    end

    @objective(model, Min, sum( sum((1-x[G, j, l])^2 - (0-x[G, j, l])^2 for l in losers[j]) +
                                sum( (0-x[G, j, k])^2 for k in 1:N)
                                for j in 1:O^G)/O^G)
    optimize!(model)
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
