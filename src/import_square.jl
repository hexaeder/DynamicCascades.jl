using NPZ
using Random
using LinearAlgebra
using NetworkLayout

"""
    import_system(:square)
"""
function import_system(::Val{:square}; seed=1, M=1u"s^2", D=0.1u"s", K=1, size=10, x=size, y=size,
                       add_e=Tuple{Int,Int}[], genbase=1u"pu", loadbase=1u"pu")
    rng = MersenneTwister(seed)

    N = x*y
    sg = Grid((x, y))

    for (src, dst) in add_e
        add_edge!(sg, src, dst)
    end

    g = MetaGraph(sg)

    Ngen = ceil(Int, N/2)
    Nload = N - Ngen
    @assert N == Ngen + Nload

    Pgen  =   N/(2*Ngen) *genbase
    Pload = - N/(2*Nload)*loadbase

    # lets split up the vertices in loads and generators
    mask = shuffle(rng, vcat(ones(Bool, Ngen), zeros(Bool, Nload)))
    gen_idxs = (1:N)[mask]
    load_idxs = (1:N)[map(!, mask)]

    # baseP = 1
    # set_prop!(g, :Pbase, baseP)
    set_prop!(g, :seed, seed)
    set_prop!(g, :NodeProps, [:n, :type, :P, :Q, :_M, :Vm])
    set_prop!(g, :EdgeProps, [:src, :dst, :X, :rating])

    set_prop!(g, 1:N, :Vm, 1.0)

    set_prop!(g, gen_idxs, :type, :gen)
    set_prop!(g, gen_idxs, :P, Pgen)
    set_prop!(g, gen_idxs, :Q, 0.0u"pu")
    set_prop!(g, gen_idxs, :_M, M)
    set_prop!(g, gen_idxs, :damping, D)

    set_prop!(g, load_idxs, :type, :gen)
    set_prop!(g, load_idxs, :P, Pload)
    set_prop!(g, load_idxs, :Q, 0.0u"pu")
    set_prop!(g, load_idxs, :_M, M)
    set_prop!(g, load_idxs, :damping, D)

    X = -1/K
    set_prop!(g, edges(g), :X, X)
    set_prop!(g, edges(g), :R, 0.0)
    set_prop!(g, edges(g), :rating, 10.0)

    set_prop!(g, 1:nv(g), :pos, squaregrid(g, cols=x))

    return g
end
