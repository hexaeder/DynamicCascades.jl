"""
    import_system(:schaefer2018)

Imports the toy model from the Schaefer 2018 paper.
"""
function import_system(::Val{:schaefer2018}; γ=0.1)
    g = MetaGraph(5)
    add_edge!(g, 1,2)
    add_edge!(g, 1,3)
    add_edge!(g, 1,5)
    add_edge!(g, 2,3)
    add_edge!(g, 2,4)
    add_edge!(g, 3,4)
    add_edge!(g, 4,5)

    set_prop!(g, :NodeProps, [:n, :P, :inertia])
    set_prop!(g, :EdgeProps, [:src, :dst, :X, :rating])

    set_prop!(g, 1:5, :Vm, 1.0)

    gen_idxs = [2, 5]
    load_idxs = [1, 3, 4]

    set_prop!(g, gen_idxs, :P, 1.5)
    set_prop!(g, gen_idxs, :type, :gen)
    set_prop!(g, load_idxs, :P, -1.0)
    set_prop!(g, load_idxs, :type, :load)
    set_prop!(g, 1:5, :model, :swing)
    set_prop!(g, 1:5, :Q, 0.0)
    set_prop!(g, 1:5, :inertia, 1.0)
    set_prop!(g, 1:5, :damping, γ)

    K = 1.63
    set_prop!(g, edges(g), :R, 0.0)
    set_prop!(g, edges(g), :X, 1/K)
    set_prop!(g, edges(g), :rating, 0.6*K)

    positions = Point2f0.([(0, 1.4),
                           (1,.1),
                           (1,1.3),
                           (0,0),
                           (-0.9,.7)])
    set_prop!(g, 1:5, :pos, positions)
    return g
end
