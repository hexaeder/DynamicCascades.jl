"""
    import_system(:toymodel)

Creates simple toymodel.
"""
# `tconst` only needed for syntax, tconst is not used in this model
function import_system(::Val{:toymodel_square}; M=1u"s^2", γ=0.1u"s", tconst=0.01u"s")
    g = MetaGraph(4)
    add_edge!(g, 1,2)
    add_edge!(g, 1,3)
    add_edge!(g, 2,4)
    add_edge!(g, 3,4)

    # set_prop!(g, :NodeProps, [:n, :P, :inertia])
    set_prop!(g, :NodeProps, [:n, :P, :_M])
    set_prop!(g, :EdgeProps, [:src, :dst, :X, :rating])

    set_prop!(g, 1:4, :Vm, 1.0)

    gen_idxs = [1, 4]
    load_idxs = [2, 3]

    set_prop!(g, gen_idxs, :P, 1.0)
    set_prop!(g, gen_idxs, :type, :gen)
    set_prop!(g, gen_idxs, :_M, M)
    set_prop!(g, load_idxs, :P, -1.0)
    set_prop!(g, load_idxs, :type, :load)
    set_prop!(g, load_idxs, :_M, M)
    set_prop!(g, 1:4, :Q, 0.0)
    # set_prop!(g, 1:5, :inertia, 1.0)
    set_prop!(g, 1:4, :damping, γ)
    # set_prop!(g, gen_idxs, :damping, 1.0u"s")
    # set_prop!(g, load_idxs, :damping, γ)
    set_prop!(g, 1:4, :timeconst, tconst)

    K = 1.0
    set_prop!(g, edges(g), :R, 0.0)
    set_prop!(g, edges(g), :X, 1/K)
    set_prop!(g, edges(g), :rating, 1.0)

    positions = Point2f.([(0, 1),
                           (1,1),
                           (0,0),
                           (1,0)])
    set_prop!(g, 1:4, :pos, positions)
    return g
end

function import_system(::Val{:toymodel_house}; M=1u"s^2", γ=0.1u"s", tconst=0.01u"s")
    g = MetaGraph(5)
    add_edge!(g, 1,2)
    add_edge!(g, 1,3)
    add_edge!(g, 2,4)
    add_edge!(g, 3,4)
    add_edge!(g, 1,5)
    add_edge!(g, 2,5)

    # set_prop!(g, :NodeProps, [:n, :P, :inertia])
    set_prop!(g, :NodeProps, [:n, :P, :_M])
    set_prop!(g, :EdgeProps, [:src, :dst, :X, :rating])

    set_prop!(g, 1:5, :Vm, 1.0)

    gen_idxs = [1, 4]
    load_idxs = [2, 3, 5]

    set_prop!(g, gen_idxs, :P, 1.5)
    set_prop!(g, gen_idxs, :type, :gen)
    set_prop!(g, gen_idxs, :_M, M)
    set_prop!(g, load_idxs, :P, -1.0)
    set_prop!(g, load_idxs, :type, :gen)
    set_prop!(g, load_idxs, :_M, M)
    set_prop!(g, 1:5, :Q, 0.0)
    # set_prop!(g, 1:5, :inertia, 1.0)
    set_prop!(g, 1:5, :damping, γ)
    set_prop!(g, 1:5, :timeconst, tconst)

    K = 1.0
    set_prop!(g, edges(g), :R, 0.0)
    set_prop!(g, edges(g), :X, 1/K)
    set_prop!(g, edges(g), :rating, 1.0)

    positions = Point2f.([(0, 1),
                           (1,1),
                           (0,0),
                           (1,0),
                           (0.5,1.5)])
    set_prop!(g, 1:5, :pos, positions)
    return g
end
