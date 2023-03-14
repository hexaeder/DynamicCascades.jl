#=
# Dynamic disturbances and topological insulators
![animation](insulator.mp4)
=#
# load used packages
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
using CairoMakie

using GLMakie #jl
GLMakie.activate!() #jl
set_theme!(theme_minimal(), fontsize = 20)

# set up the grid
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
nothing #hide

#=
Simulate failure of all lines, analyise solution and save maximum pertubation on lines
=#
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
    @info "Simulate for edge $i" #jl
    edgetype[i] == 0 && continue
    direction = edge.dst == edge.src+1 ? :h : :v
    sol = simulate(network;
                   initial_fail = [i],
                   failtime = 0.1,
                   trip_lines = :none,
                   tspan = (0.0, 100.0),warn=false,verbose=false);
    (t, this, other) = max_flowdiff(sol)
    push!(df, (; edge=i, direction, thismax=maximum(this), othermax=maximum(other), endthis=this[end], endother=other[end]))
end

sort(df, :endother) #jl
sort(df[df.direction.==:v, :], :othermax) #jl
sort(df[df.direction.==:h, :], :othermax) #jl
nothing #hide

# maximum disturbane in *other* region after vertical fail
extrema(df[df.direction.==:v, :othermax])
# maximum disturbane in *same* region after vertical fail
extrema(df[df.direction.==:v, :thismax])

# maximum disturbane in *other* region after horizontal fail
extrema(df[df.direction.==:h, :othermax])
# maximum disturbane in *same* region after horizontal fail
extrema(df[df.direction.==:h, :thismax])

#=
Create plot for thesis
=#
sol1 = simulate(network;
               initial_fail = [60],
               failtime = 0.1,
               trip_lines = :none,
               tspan = (0.0, 50.0),
               solverargs = (; dtmax = 0.01));
inspect_solution(sol1) #jl

sol2 = simulate(network;
               initial_fail = [18],
               failtime = 0.1,
               trip_lines = :none,
               tspan = (0.0, 50.0),
               solverargs = (; dtmax = 0.01));
inspect_solution(sol2) #jl

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

save(joinpath(PLOT_DIR, "insulator.pdf"), fig) #jl
fig

# create the video
fig = Figure(resolution=(1000, 400))
scaling = Observable(0.5)
t = Observable(0.0)
fig[1,1] = Label(fig, @lift("t = "*repr(round($t,digits=2))*" s"), tellwidth=false)
fig[2,1] = graphplot_axis_edges(fig, sol1, t, scaling)
fig[3,1] = graphplot_axis_edges(fig, sol2, t, scaling)
fig[2, 0] = Label(fig, "horizontal failure", tellheight=false)
fig[3, 1] = Label(fig, "vertical failure", tellheight=false)


T = 10
tmax = 3
fps = 30
trange = range(0.0, tmax, length=Int(T * fps))

record(fig, "insulator.mp4", trange; framerate=30) do time
    t[] = time
end

#=
Create plot for Defense
=#
using GLMakie; GLMakie.activate!() #jl
x = 30
y = 3
add_e = [(15, 46),
         (15, 76),
         (45, 16),
         (45, 76),
         (75, 16),
         (75, 46)]
network1 = import_system(:square; x, y, M=0.5u"s^2", D=0.1u"s", K=-10u"pu", seed=1)
network2 = import_system(:square; x, y, M=0.5u"s^2", D=0.1u"s", K=-10u"pu", seed=1, add_e)
@assert describe_nodes(network1).P[15] != describe_nodes(network1).P[16]

sol1 = simulate(network1;
               initial_fail = [86],
               failtime = 0.1,
               trip_lines = :none,
               tspan = (0.0, 150.0),
               solverargs = (; dtmax = 0.01));
# inspect_solution(sol1) #jl

sol2 = simulate(network2;
               initial_fail = [90],
               failtime = 0.1,
               trip_lines = :none,
               tspan = (0.0, 150.0),
               solverargs = (; dtmax = 0.01));
# inspect_solution(sol2) #jl

right1 = collect(88:2:114)
# left1  = reverse(collect(62:2:84))

right2 = vcat(92,95, 98:2:120)
# left2  = reverse(collect(66:2:88))

fig = Figure(resolution=(1800, 600))
ax1a = fig[1,1] = Axis(fig, xlabel="Time in s", ylabel="Flow Disturbance", title="w/o insulator")
max_dist_1 = Point2[]
maxt_dist_1 = Point2[]
steady_dist_1 = Point2[]
for (dist, idx) in enumerate(right1)
    t, P = seriesforidx(sol1.load_P, idx, f = P ->(P - sol1.load_P.saveval[1][idx]))
    if dist in [2,4,6]
        lines!(ax1a, t, P, label="dist = $dist")
    end
    P = abs.(P)
    pmax, i = firstmax(P)
    tmax = t[i]
    push!(max_dist_1, Point(dist, pmax))
    push!(maxt_dist_1, Point(dist, tmax))
    push!(steady_dist_1, Point(dist, P[end]))
end

ax1b = fig[1,2]=Axis(fig, xlabel="Time in s", ylabel="Flow Disturbance", title="with insulator")
max_dist_2 = Point2[]
maxt_dist_2 = Point2[]
steady_dist_2 = Point2[]
for (dist, idx) in enumerate(right2)
    t, P = seriesforidx(sol2.load_P, idx, f = P -> (P - sol2.load_P.saveval[1][idx]))
    if dist in [2,4,6]
        lines!(ax1b, t, P, label="dist = $dist")
    end
    P = abs.(P)
    pmax, i = firstmax(P)
    tmax = t[i]
    push!(max_dist_2, Point(dist, pmax))
    push!(maxt_dist_2, Point(dist, tmax))
    push!(steady_dist_2, Point(dist, P[end]))
end
axislegend(ax1a)
axislegend(ax1b)
xlims!(ax1a, 0, 3)
ylims!(ax1a, -0.15, 0.3)
xlims!(ax1b, 0, 3)
ylims!(ax1b, -0.15, 0.3)

save("insulator_transients.png", fig)

fig = Figure(resolution=(1800, 600))
ax1 = fig[1,1] = Axis(fig, xlabel="Distance from incidence", ylabel="Maximum disturbance in transient")
scatterlines!(ax1, max_dist_1, label="w/o insulator")
scatterlines!(ax1, max_dist_2, label="with insulator")
axislegend(ax1)

ax2 = fig[1,2] = Axis(fig, xlabel="Distance from incidence", ylabel="steady state")
scatterlines!(ax2, steady_dist_1, label="w/o insulator")
scatterlines!(ax2, steady_dist_2, label="with insulator")
axislegend(ax2)
save("insulator_pmax.png", fig)

# video
t = Observable(0.0)
fig = Figure(resolution=(1200, 400))
ax1 = fig[1,1] = Axis(fig)
ax2 = fig[2,1] = Axis(fig)
hidedecorations!(ax1), hidespines!(ax1)
hidedecorations!(ax2), hidespines!(ax2)
ax1.aspect = DataAspect()
ax2.aspect = DataAspect()
gpargs1 = gparguments(sol1, t;
                     ecolortype=Observable(:abssteady),
                     Δω=Observable(3.0),
                     offlinewidth=5,
                     offlinecolor=colorant"black",
                     ecolorscaling=Observable(0.5),
                     node_size=10)
gpargs2 = gparguments(sol2, t;
                     ecolortype=Observable(:abssteady),
                     Δω=Observable(3.0),
                     offlinewidth=5,
                     offlinecolor=colorant"black",
                     ecolorscaling=Observable(0.5),
                     node_size=10)
p1 = graphplot!(ax1, network1; gpargs1..., edge_width=5)
p2 = graphplot!(ax2, network2; gpargs2..., edge_width=5)

T = 20
fps = 30
tmax = 5
trange = range(0, tmax, length=round(Int, T * fps))
record(fig, "insulator_video.mp4", trange; framerate=30) do time
    t[] = time
end

t[] = 100.0
save("insulator_final_state.png", fig)
