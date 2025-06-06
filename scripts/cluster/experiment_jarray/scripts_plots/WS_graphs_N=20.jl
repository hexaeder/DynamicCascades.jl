"""
Watts-Strogatz example networks with N=10 for showing the effect of different
β-values
"""

using Graphs
using CairoMakie
using GraphMakie
using GraphMakie.NetworkLayout
# It is `watts_strogatz(n, k, β)`

watts_strogatz(20, 4, 0.1, seed=1) # {100, 200} undirected simple Int64 graph
watts_strogatz(20, 4, 0.9, seed=1) # {100, 200} undirected simple Int64 graph
# As expected the graphs have the same number of edges and nodes as n and k are not changed

# Now plot the graphs
function simple_WS_plot(g, β)
    f, ax, p = graphplot(g, layout=Shell())

    hidedecorations!(ax); hidespines!(ax)
    ax.aspect = DataAspect()

    # # Getting the positions on a circle
    # angles = range(0, stop=2π, length=nv(g)+1)[1:end-1]
    # posx = cos.(angles)
    # posy = sin.(angles)
    # positions = Point2f.([(i, j) for (i, j) in zip(posx, posy)])
    #
    # fixed_layout(_) = positions
    # # set new layout
    # p.layout = fixed_layout; autolimits!(ax)
    p.edge_width = 3.0
    p.node_size=20
    ax.title = "β=$β"
    ax.titlesize = 40

    return f
end

N=10
β = [0.1, 0.25, 0.5]
for i in β
    CairoMakie.save(joinpath(MA_DIR, "WS", "WS_β=$i.pdf"),simple_WS_plot(watts_strogatz(N, 4, i, seed=1),i))
    CairoMakie.save(joinpath(MA_DIR, "WS", "WS_β=$i.pdf"),simple_WS_plot(watts_strogatz(N, 4, i, seed=1),i))
    CairoMakie.save(joinpath(MA_DIR, "WS", "WS_β=$i.pdf"),simple_WS_plot(watts_strogatz(N, 4, i, seed=1),i))
end
# Explanation HW
f, ax, p = graphplot(watts_strogatz(100, 4, 0.1, seed=1); layout=Shell())
length(get_edge_plot(p).paths[])
# 200

f, ax, p = graphplot(watts_strogatz(100, 4, 0.9, seed=1); layout=Shell())
length(get_edge_plot(p).paths[])
# 200

#= Both graphs actually have the same number of lines. It just looks so different
because in the first example you connect mostly to your direct neighbours in Shell
layout, so you don't really see all of the edges.
You could make them visible permuting the nodes on the circle:=#

f, ax, p = graphplot(watts_strogatz(100, 4, 0.1, seed=1); layout=Shell())
p.node_pos[] = shuffle(p.node_pos[])

# NB One simply does not see the edges but they are there.
