"""
    import_system(:nadir_sim)

Creates a toymodel in which two nodes are connected by a line. One node is a
swing equation node, the other node is a slack node.
"""
# `tconst` only needed for syntax, tconst is not used in this model
function import_system(::Val{:nadir_sim}; M=1u"s^2", γ=0.5u"s", tconst=0.0u"s")
    g = MetaGraph(2)
    add_edge!(g, 1,2)

    # set_prop!(g, :NodeProps, [:n, :P, :inertia])
    set_prop!(g, :NodeProps, [:n, :P, :_M])
    set_prop!(g, :EdgeProps, [:src, :dst, :X, :rating])

    set_prop!(g, 1:2, :Vm, 1.0)

    gen_idxs = [1]
    load_idxs = [2]

    set_prop!(g, gen_idxs, :P, 0.0)
    set_prop!(g, gen_idxs, :type, :gen)
    set_prop!(g, gen_idxs, :_M, M)
    set_prop!(g, load_idxs, :P, 0.0)
    set_prop!(g, load_idxs, :type, :load)
    set_prop!(g, load_idxs, :_M, M)
    set_prop!(g, 1:2, :Q, 0.0)
    # set_prop!(g, 1:5, :inertia, 1.0)
    set_prop!(g, 1:2, :damping, γ)
    # set_prop!(g, gen_idxs, :damping, 1.0u"s")
    # set_prop!(g, load_idxs, :damping, γ)
    set_prop!(g, 1:2, :timeconst, tconst)

    K = 1.0
    set_prop!(g, edges(g), :R, 0.0)
    set_prop!(g, edges(g), :X, 1/K)
    set_prop!(g, edges(g), :rating, 1.0)

    positions = Point2f.([(0,0),
                           (1,0)])
    set_prop!(g, 1:2, :pos, positions)
    return g
end
