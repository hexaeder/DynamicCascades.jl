using Random
using Distributions
using NetworkLayout

# NOTE
#=Inconsistency of RNG across the different Julia versions for WS-graphs (and probably ER-graphs (not tested)).
Power injection are consistent between Julia v1.8.4 and v1.11.0 but but are saved to be on the save side.=#

"""
    import_system(:erdosrenyi)

Creates Erdős–Rényi-graph.
"""
function import_system(::Val{:erdosrenyi}; N=20, p=0.25, graph_seed=123, μ=0.0, σ=1.0, distr_seed=1234, K=5.0, α=0.7, M=1u"s^2", γ=0.1u"s", τ=0.01u"s")
    @info "Import system Erdős–Rényi"
    if graph_seed !== nothing
        g = MetaGraph(erdos_renyi(N, p, seed=graph_seed))
    else
        g = MetaGraph(erdos_renyi(N, p))
    end

    set_params_general_neworks!(g, distr_seed, N, M, γ, τ, K, α, μ, σ)
    return g
end

"""
    import_system(:wattsstrogatz)

Creates Watts-Strogatz-graph.
"""
function import_system(::Val{:wattsstrogatz}; N=20, k=4, β=0.5, graph_seed=123, μ=0.0, σ=1.0, distr_seed=1234, K=5.0, α=0.7, M=1u"s^2",  γ=0.1u"s", τ=0.01u"s")
    @info "Import system Watts-Strogatz"
    if graph_seed !== nothing
        g = MetaGraph(watts_strogatz(N, k, β, seed=graph_seed))
    else
        g = MetaGraph(watts_strogatz(N, k, β))
    end
    set_prop!(g, :graph_seed, graph_seed)

    set_params_general_neworks!(g, distr_seed, N, M, γ, τ, K, α, μ, σ)
    return g
end

"""
Generates parameters for Erdős–Rényi and Watts-Strogatz-graph, includes distribution for power injections.
    NOTE: `positions` might not fit to Erdős–Rényi-graph
"""
function set_params_general_neworks!(g, distr_seed, N, M, γ, τ, K, α, μ, σ)
    set_prop!(g, :NodeProps, [:n, :Pmech, :Pload, :_M])
    set_prop!(g, :EdgeProps, [:src, :dst, :X, :rating])

    set_prop!(g, 1:nv(g), :Vm, 1.0)

    set_prop!(g, 1:nv(g), :type, :gen)
    set_prop!(g, 1:nv(g), :_M, M)

    set_prop!(g, 1:nv(g), :damping, γ)
    set_prop!(g, 1:nv(g), :timeconst, τ)

    # NOTE
    #= Would have been easier to set coupling directly here: `set_prop!(network, e, :_K, K[i])`
    However this would be problematic: For `calculate_apparent_power!` we need `:_Y`, which is only assigned
    a value in `set_coupling!` if `:_K` is not assigned.=#
    set_prop!(g, edges(g), :R, 0.0)
    set_prop!(g, edges(g), :K, K)  
    set_prop!(g, edges(g), :α, α)
    set_prop!(g, edges(g), :X, 1/K)
    set_prop!(g, edges(g), :rating, α*K)

    set_prop!(g, :distr_seed, distr_seed)
    set_prop!(g, :μ, μ)
    set_prop!(g, :σ, σ)
    d = Distributions.Normal(μ, σ)
    if distr_seed !== nothing
        Random.seed!(distr_seed)
    end
    x = rand(d, 2N)
    Pmech = x[1:N]
    Pload = x[N+1:2N]

    set_prop!(g, 1:nv(g), :Pmech, Pmech)
    set_prop!(g, 1:nv(g), :Pload, Pload)

    balance_power!(g)
    # set_prop!(g, 1:nv(g), :pos, spring(g; C=5.0))
    angles = range(0, stop=2π, length=nv(g)+1)[1:end-1]
    posx = cos.(angles)
    posy = sin.(angles)
    positions = Point2f.([(i, j) for (i, j) in zip(posx, posy)])
    set_prop!(g, 1:nv(g), :pos, positions)
    return g
end


# angles = range(0, stop=2π, length=nv(g)+1)[1:end-1]
# x = cos.(angles)
# y = sin.(angles)
# set_prop!(g, 1:nv(g), :pos, [x y])
