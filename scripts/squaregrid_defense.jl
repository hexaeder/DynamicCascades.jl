
using DynamicCascades
using Graphs
using MetaGraphs
using Statistics
using GraphMakie
using Colors
using ColorSchemes
using DynamicCascades: PLOT_DIR, simulate_pdis
using Unitful
using NetworkDynamics
using CairoMakie

using GLMakie #jl
GLMakie.activate!() #jl
set_theme!(theme_minimal(), fontsize = 20)

x = 50
y = 3
network = import_system(:square; x, y, M=0.5u"s^2", D=0.1u"s", K=-10u"pu", seed=2,
                        genbase=1.0u"pu", loadbase=1.0u"pu")

sol = simulate(network;
               initial_fail = [150],
               failtime = 0.1,
               trip_lines = :none,
               tspan = (0.0, 5.0),
               solverargs = (; dtmax = 0.01));
inspect_solution(sol)

fig = Figure(resolution=(1000, 600))
ax = fig[1,1]=Axis(fig, xlabel="Time in s", ylabel="Flow Disturbance")
for idx in [164, 168, 172]
    t, P = seriesforidx(sol.load_P, idx, f = P -> (P - sol.load_P.saveval[1][idx]))
    lines!(ax, t, P)
end


fig = Figure(resolution=(1000, 600))
ax = fig[1,1]=Axis(fig, xlabel="Time in s", ylabel="Flow Disturbance")
for Kraw in 10:50:200
    x = 50; y = 3
    network = import_system(:square; x, y, M=0.5u"s^2", D=0.1u"s", K=-1*Kraw*u"pu", seed=2,
                            genbase=1.0u"pu", loadbase=1.0u"pu")
    sol = simulate(network;
                   initial_fail = [150],
                   failtime = 0.1,
                   trip_lines = :none,
                   tspan = (0.0, 5.0),
                   solverargs = (; dtmax = 0.01));
    idx = 164
    t, P = seriesforidx(sol.load_P, idx, f = P -> (P - sol.load_P.saveval[1][idx]))
    lines!(ax, t, P)
end


idx = 168
t, P = seriesforidx(sol.load_P, idx, f = P -> abs(P - sol.load_P.saveval[1][idx]))
plot!(t, P)

# x = 50
# y = 1
# network = import_system(:square; x, y, M=10.5u"s^2", D=0.1u"s", K=-100u"pu", seed=2,
#                         genbase=0.0u"pu", loadbase=0.0u"pu")
# sol = simulate_pdis(network; node=1, disturbance=0.1u"pu",
#                solverargs = (; dtmax = 0.01), tspan=(0.0, 4.0));
# inspect_solution(sol)

set_theme!(theme_minimal(), fontsize = 20, linewidth=4)

function firstmax(P)
    last = P[begin]
    for i in eachindex(P)
        new = P[i]
        if new < last
            return last, i-1
        end
        last = new
    end
    return P[end], lastindex(P)
end

using GLMakie; GLMakie.activate!()
x = 30
y = 30
Kbase = 10
network = import_system(:square; x, y, M=0.5u"s^2", D=0.1u"s", K=-1*Kbase*u"pu", seed=2,
                        genbase=1.0u"pu", loadbase=1.0u"pu")

sol = simulate(network;
               initial_fail = [855],
               failtime = 0.1,
               trip_lines = :none,
               tspan = (0.0, 100.0),
               solverargs = (; dtmax = 0.01));
inspect_solution(sol)

vertical = reverse(collect(29:59:855))[2:end-1]
horizontal = collect(855:2:883)[2:end-1]

fig = Figure(resolution=(1800, 600))
ax1a = fig[1,1] = Axis(fig, xlabel="Time in s", ylabel="Flow Disturbance", title="vertical")
xlims!(ax1a, 0, 3)
ylims!(ax1a, -.1, 0.8)
max_dist_vert = Point2[]
maxt_dist_vert = Point2[]
steady_dist_vert = Point2[]
for (dist, idx) in enumerate(vertical)
    t, P = seriesforidx(sol.load_P, idx, f = P ->(P - sol.load_P.saveval[1][idx]))
    if dist in [1,2,3]
        lines!(ax1a, t, -1 .* P, label="dist = $dist")
    end
    P = abs.(P)
    pmax, i = firstmax(P)
    tmax = t[i]
    push!(max_dist_vert, Point(dist, pmax))
    push!(maxt_dist_vert, Point(dist, tmax))
    push!(steady_dist_vert, Point(dist, P[end]))
end
axislegend(ax1a)

ax1b = fig[1,2]=Axis(fig, xlabel="Time in s", ylabel="Flow Disturbance", title="horizontal")
xlims!(ax1b, 0, 3)
ylims!(ax1a, -.1, 0.8)
max_dist_hori = Point2[]
maxt_dist_hori = Point2[]
steady_dist_hori = Point2[]
for (dist, idx) in enumerate(horizontal)
    t, P = seriesforidx(sol.load_P, idx, f = P -> (P - sol.load_P.saveval[1][idx]))
    if dist in [1,2,3]
        lines!(ax1b, t, P, label="dist = $dist")
    end
    P = abs.(P)
    pmax, i = firstmax(P)
    tmax = t[i]
    push!(max_dist_hori, Point(dist, pmax))
    push!(maxt_dist_hori, Point(dist, tmax))
    push!(steady_dist_hori, Point(dist, P[end]))
end
axislegend(ax1b)
save("grid_transients.png", fig)


fig = Figure(resolution=(1800, 600))
ax1 = fig[1,1] = Axis(fig, xlabel="Distance from incidence", ylabel="Maximum disturbance in transient")
scatterlines!(ax1, max_dist_vert, label="vertical")
scatterlines!(ax1, max_dist_hori, label="horizontal")
axislegend(ax1)

ax2 = fig[1,2] = Axis(fig, xlabel="Distance from incidence", ylabel="steady state")
scatterlines!(ax2, steady_dist_vert, label="vertical")
scatterlines!(ax2, steady_dist_hori, label="horizontal")
axislegend(ax2)

save("grid_max.png", fig)

####
#### Video
####
t = Observable(0.0)

fig = Figure(resolution=(1000,1000))
ax = fig[1,1] = Axis(fig)
hidedecorations!(ax), hidespines!(ax)
gpargs = gparguments(sol, t;
                     ecolortype=Observable(:abssteady),
                     Δω=Observable(3.0),
                     offlinewidth=5,
                     offlinecolor=colorant"black",
                     ecolorscaling=Observable(0.2),
                     node_size=10)
p = graphplot!(ax, network; gpargs..., edge_width=10)

T = 20
fps = 30
tmax = 20
trange = range(0, tmax, length=round(Int, T * fps))
record(fig, "grid_video.mp4", trange; framerate=30) do time
    t[] = time
end

t[] = 3.0
save("grid_overview.png", fig)

t = Observable(0.0)

fig = Figure(resolution=(1000,1000))
ax = fig[1,1] = Axis(fig)
hidedecorations!(ax), hidespines!(ax)
gpargs = gparguments(sol, t;
                     ecolortype=Observable(:abssteady),
                     Δω=Observable(0.05),
                     offlinewidth=5,
                     offlinecolor=colorant"black",
                     ecolorscaling=Observable(1000),
                     node_size=30)
p = graphplot!(ax, network; gpargs..., edge_width=10)

T = 10
fps = 30
tmax = 20
trange = range(0, tmax, length=round(Int, T * fps))
record(fig, "grid_video_nodes.mp4", trange; framerate=30) do time
    t[] = time
end

t[] = 0.6
save("grid_overview_nodes.png", fig)

inspect_solution(sol)

nidxs = [436, 437, 439]

nd = nd_model(sol.network)[1]

states = [findfirst(isequal(Symbol("ω_$i")), nd.syms) for i in nidxs]

fig = Figure()
ax = fig[1,1] = Axis(fig; xlabel="time in s", ylabel="node frequency in hz")
for (i,n) in enumerate(states)
    t, ω = seriesforidx(sol.sol, n)
    lines!(ax, t, ω; label="Node $i")
end
xlims!(0, 3)
axislegend(ax)
save("grid_node_frequencies.png", fig)
