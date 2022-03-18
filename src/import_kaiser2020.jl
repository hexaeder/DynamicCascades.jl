using NPZ
using Random
using LinearAlgebra
using NetworkLayout

"""
    import_system(:kaiser2020)

Imports the Laplacian from the Kaiser2020 paper.
"""
function import_system(::Val{:kaiser2020}; gen_γ, load_τ, seed=1)
    rng = MersenneTwister(seed)
    file = joinpath(@__DIR__, "..", "data", "Kaiser2020", "Laplacian.npy")
    # file = joinpath(DATA_DIR, "Kaiser2020", "Laplacian.npy")
    L = npzread(file)
    A = - L + Diagonal(L)
    @assert unique(A[:]) == [0, 1] "There are weights?"

    g = MetaGraph(A)

    N = nv(g)
    Ngen = 15
    Nload = 25
    @assert N == Ngen + Nload

    Pgen = 10.0/Ngen
    Pload = Pgen * Ngen / Nload
    # Pgen = rand_with_sum(rng, Ngen, 20)
    # Pload = rand_with_sum(rng, Nload, 20)

    # lets split up the vertices in loads and generators
    mask = shuffle(rng, vcat(ones(Bool, Ngen), zeros(Bool, Nload)))
    gen_idxs = (1:N)[mask]
    load_idxs = (1:N)[map(!, mask)]

    # baseP = 1
    # set_prop!(g, :Pbase, baseP)
    set_prop!(g, :seed, seed)
    set_prop!(g, :NodeProps, [:n, :type, :P, :Q, :inertia, :Vm])
    set_prop!(g, :EdgeProps, [:src, :dst, :R, :X, :rating])

    set_prop!(g, 1:N, :Vm, 1.0)

    set_prop!(g, gen_idxs, :type, :gen)
    set_prop!(g, gen_idxs, :P, Pgen)
    set_prop!(g, gen_idxs, :Q, 0.0)
    set_prop!(g, gen_idxs, :inertia, 1.0)
    set_prop!(g, gen_idxs, :damping, gen_γ)

    set_prop!(g, load_idxs, :type, :load)
    set_prop!(g, load_idxs, :P, -Pload)
    set_prop!(g, load_idxs, :Q, 0.0)
    set_prop!(g, load_idxs, :inertia, 1.0)
    set_prop!(g, load_idxs, :timeconst, load_τ)

    set_prop!(g, edges(g), :R, 0.0)
    set_prop!(g, edges(g), :X, 0.01)
    set_prop!(g, edges(g), :rating, 2.0)

    # set the regions
    set_prop!(g, 1:20, :region, 1)
    set_prop!(g, 21:40, :region, 2)
    for e in edges(g)
        r1 = get_prop(g, e.src, :region)
        r2 = get_prop(g, e.dst, :region)
        if r1 == r2
            set_prop!(g, e, :region, r1)
        else
            set_prop!(g, e, :region, 0)
        end
    end

    set_prop!(g, 1:nv(g), :pos, spring(g))

    return g
end

function rand_with_sum(rng, N, total)
    a = [2*rand(rng)*total/N for i in 1:N]
    last = total - sum(a[1:N-1])
    @assert last > 0
    a[end] = last
    return a
end
