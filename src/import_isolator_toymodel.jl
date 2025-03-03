"""
    import_system(:isolator_toymodel)

Imports the a toy model that consists of two identical toymodels from the
Schaefer 2018 paper that are connected by a network isolator introduced by
Kaiser et al. Network isolators inhibit failure spreading in complex networks.
"""
# `tconst` only needed for syntax, tconst is not used in this model
function import_system(::Val{:isolator_toymodel}; γ=0.1u"s", M=1u"s^2", tconst=0.01u"s")
    g = MetaGraph(10)
    # orginal toymodel
    add_edge!(g, 1,2)
    add_edge!(g, 1,3)
    add_edge!(g, 1,5)
    add_edge!(g, 2,3)
    add_edge!(g, 2,4)
    add_edge!(g, 3,4)
    add_edge!(g, 4,5)

    #toymodel duplicate
    # `nae` stands for number of additional edges
    nae = 5
    add_edge!(g, 1 + nae, 2 + nae)
    add_edge!(g, 1 + nae, 3 + nae)
    add_edge!(g, 1 + nae, 5 + nae)
    add_edge!(g, 2 + nae, 3 + nae)
    add_edge!(g, 2 + nae, 4 + nae)
    add_edge!(g, 3 + nae, 4 + nae)
    add_edge!(g, 4 + nae, 5 + nae)

    # isolator graph
    add_edge!(g, 3,8)
    add_edge!(g, 3,7)
    add_edge!(g, 2,8)
    add_edge!(g, 2,7)

    # set_prop!(g, :NodeProps, [:n, :P, :inertia])
    set_prop!(g, :NodeProps, [:n, :P, :_M])
    set_prop!(g, :EdgeProps, [:src, :dst, :X, :rating])

    set_prop!(g, 1:10, :Vm, 1.0)

    gen_idxs = [2, 5, 2 + nae, 5 + nae]
    load_idxs = [1, 3, 4, 1 + nae, 3 + nae, 4 + nae]

    set_prop!(g, gen_idxs, :P, 1.5)
    set_prop!(g, gen_idxs, :type, :gen)
    set_prop!(g, gen_idxs, :_M, M)
    set_prop!(g, load_idxs, :P, -1.0)
    set_prop!(g, load_idxs, :type, :gen) # NOTE due to restrutured `nd_model`
    set_prop!(g, load_idxs, :_M, M)
    set_prop!(g, 1:10, :model, :swing)
    set_prop!(g, 1:10, :Q, 0.0)
    # set_prop!(g, 1:5, :inertia, 1.0)
    set_prop!(g, 1:10, :damping, γ)
    set_prop!(g, 1:10, :timeconst, tconst)

    K = 1.63
    set_prop!(g, edges(g), :R, 0.0)
    set_prop!(g, edges(g), :X, 1/K)
    set_prop!(g, edges(g), :rating, 0.535*K)

    positions = Point2f.([(0, 1.4),
                           (1,.1),
                           (1,1.3),
                           (0,0),
                           (-0.9,.7),
                           # vertical symmetry axis at x = 1.5
                           (3, 1.4),
                           (2,.1),
                           (2,1.3),
                           (3,0),
                           (3.9,.7)])
    set_prop!(g, 1:10, :pos, positions)
    return g
end
