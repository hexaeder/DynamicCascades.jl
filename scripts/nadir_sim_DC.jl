using Revise
using DynamicCascades
using NetworkDynamics
using Graphs
using MetaGraphs
using Unitful
using Statistics
using GraphMakie
using Colors
using DynamicCascades: PLOT_DIR

using CairoMakie

using GLMakie
GLMakie.activate!()

I = 1
D = 3
K = 1
Z = 4*I*K - D^2

network = import_system(:nadir_sim; M=1.0u"s^2", γ=0.5u"s", tconst=0.0u"s")
sol = simulate(network;
    initial_fail = [1],
    failtime = 1.0,
    init_pert = :power_perturbation,
    P_perturb = 0.5,
    trip_lines = :none,
    trip_nodes = :none,
    trip_load_nodes = :none,
    tspan = (0.0, 50.0),
    solverargs = (; dtmax = 0.01));

display(inspect_solution(sol))

# plot solution
fig = Figure(resolution=(1800,1000))

# plot power flow
fig[1,1] = ax = Axis(fig; xlabel="time t in s", ylabel="apparent power flow in p.u.", title="power flow")
for i in 1:length(sol.load_S.saveval[1])
    t = sol.load_S.t
    y = [sol.load_S.saveval[t][i] for t in 1:length(sol.load_S.t)]
    lines!(ax, t, y; label="power flow on line ($i)", linewidth=3)
end
vlines!(ax, tobs; color=:black, linewidth=1)


# plot frequency of gen node
fig[1,2] = ax = Axis(fig; xlabel="time t in s", ylabel="angular frequency ω", title="frequency transients of all generator nodes")
(nd, p, overload_cb) = nd_model(network)
state_idx = idx_containing(nd, "ω") # array: indices of ω-states
for i in state_idx
    t = sol.sol.t
    y = [sol.sol.u[t][i] for t in 1:length(sol.sol.t)]
    lines!(ax, t, y; label="frequency on node ($i)", linewidth=3)
end
vlines!(ax, tobs; color=:black, linewidth=1)

fig[2,1] = ax = Axis(fig; xlabel="time t in s", ylabel="phase angle θ", title="frequency transients of all generator nodes")
(nd, p, overload_cb) = nd_model(network)
state_idx = idx_containing(nd, "θ") # array: indices of ω-states
for i in state_idx
    t = sol.sol.t
    y = [sol.sol.u[t][i] for t in 1:length(sol.sol.t)]
    lines!(ax, t, y; label="frequency on node ($i)", linewidth=3)
end
vlines!(ax, tobs; color=:black, linewidth=1)

# plot frequency of slack
fig[1,3] = ax = Axis(fig; xlabel="time t in s", ylabel="angular frequency ω", title="frequency of slack")
for i in 1:length(sol.frequencies_load_nodes.saveval[1])
    t = sol.frequencies_load_nodes.t
    # t = sol.frequencies_load_nodes.t[1:20]
    y = [sol.frequencies_load_nodes.saveval[t][i] for t in 1:length(sol.frequencies_load_nodes.t)]
    # y = [sol.frequencies_load_nodes.saveval[t][i] for t in 1:20]
    lines!(ax, t, y; label="frequency on load node ($i)", linewidth=3)
    # scatter!(ax, t, y; label="frequency on load node ($i)", markersize=5)
end
vlines!(ax, tobs; color=:black, linewidth=1)

fig

# save(joinpath(PLOT_DIR, "nadir_sim.pdf"), fig)
