using DynamicCascades
using Graphs
using MetaGraphs
using Statistics
using GraphMakie
using Colors
using ColorSchemes
using DynamicCascades: PLOT_DIR
using Unitful
using NetworkDynamics
using DataFrames

# using CairoMakie

using GLMakie
GLMakie.activate!()
set_theme!(theme_minimal(), fontsize = 20)
# set_theme!(theme_light())

# 20 => 370
# 30 => 855

x = 20
y = 3
add_e = [(10, 31),
         (10, 51),
         (30, 11),
         (30, 51),
         (50, 11),
         (50, 31)]
network = import_system(:square; x, y, M=0.5u"s^2", D=0.1u"s", K=-10u"pu", seed=2, add_e)
region1 = vcat(1:10, 21:30, 41:50);
region2 = vcat(11:20, 31:40, 51:60);

edgetype = map(edges(network)) do e
    if e.src ∈ region1 && e.dst ∈ region1
        1
    elseif e.src ∈ region2 && e.dst ∈ region2
        2
    else
        0
    end
end

function max_flowdiff(sol::SolutionContainer)
    init_region = edgetype[sol.initial_fail[]]
    other_region = init_region == 1 ? 2 : 1

    rege = [findall(isequal(i), edgetype) for i in 1:2]

    static = sol.load_P.saveval[begin]
    diff = [(abs.(val)-abs.(static))/static[sol.initial_fail[]] for val in sol.load_P.saveval]

    t = sol.load_P.t
    this  = [maximum(d[rege[init_region]]) for d in diff]
    other = [maximum(d[rege[other_region]]) for d in diff]

    return (t, this, other)
end

df = DataFrame(; edge = Int[], direction=Symbol[], thismax=Float64[], othermax=Float64[],
               endthis=Float64[], endother=Float64[])
for (i, edge) in enumerate(edges(network))
    @info "Simulate for edge $i"
    edgetype[i] == 0 && continue
    direction = edge.dst == edge.src+1 ? :h : :v
    sol = simulate(network;
                   initial_fail = [i],
                   failtime = 0.1,
                   trip_lines = :none,
                   tspan = (0.0, 100.0),);
    (t, this, other) = max_flowdiff(sol)
    push!(df, (; edge=i, direction, thismax=maximum(this), othermax=maximum(other), endthis=this[end], endother=other[end]))
end

sort(df, :endother)

sort(df[df.direction.==:v, :], :othermax)
sort(df[df.direction.==:h, :], :othermax)

extrema(df[df.direction.==:v, :othermax])
extrema(df[df.direction.==:v, :thismax])

extrema(df[df.direction.==:h, :othermax])
extrema(df[df.direction.==:h, :thismax])


sol1 = simulate(network;
               initial_fail = [60],
               failtime = 0.1,
               trip_lines = :none,
               tspan = (0.0, 50.0),
               solverargs = (; dtmax = 0.01));
inspect_solution(sol1)

sol2 = simulate(network;
               initial_fail = [18],
               failtime = 0.1,
               trip_lines = :none,
               tspan = (0.0, 50.0),
               solverargs = (; dtmax = 0.01));
inspect_solution(sol2)

fig = Figure(resolution=(1500, 600))
ts = [Observable(0.25),
      Observable(0.5),
      Observable(1.5),
      Observable(50.0)]
function graphplot_axis_edges(fig, sol, t::Observable, scaling)
    ax = Axis(fig)
    gpargs = gparguments(sol, t;
                         ecolortype = Observable(:abssteadyboth),
                         activeP=true,
                         ecolorscaling=scaling,
                         )
    p = graphplot!(ax, network; gpargs...,
                   node_size=0,
                   edge_width=5)

    ax.aspect = DataAspect()
    hidedecorations!(ax); hidespines!(ax)
    return ax
end
scaling = Observable(0.25)
for (i,t) in pairs(ts)
    fig[i,1] = Label(fig, "t = $(t[])", tellheight=false)
    fig[i,2] = graphplot_axis_edges(fig, sol1, t, scaling)
    fig[i,3] = graphplot_axis_edges(fig, sol2, t, scaling)
end
fig[0, 2] = Label(fig, "horizontal failure", tellwidth=false)
fig[1, 3] = Label(fig, "vertical failure", tellwidth=false)

# CairoMakie.activate!()
save(joinpath(PLOT_DIR, "insulator.pdf"), fig)
