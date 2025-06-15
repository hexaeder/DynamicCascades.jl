"""
See scripts/archive_nb/porting_to_new_ND/testing_new_ND.jl
"""



# Save steady state (only for θ) for comaprison with new ND
zeroidx=1
(nd, p) = nd_model(network)
x0 = zeros(length(nd.syms));
x_static = solve(SteadyStateProblem(nd, x0, p), NLSolveJL())

θidx = idx_containing(nd, "θ")
offset = x_static[θidx[zeroidx]]
x_static[θidx] .= x_static[θidx] .- offset
@assert iszero(x_static[θidx[zeroidx]])

steady_state_dict = Dict(:SteadyState => x_static[θidx])
CSV.write("/home/brandner/.julia/dev/for_testing_delete_soon/RTS_test_graph_old_ND_steady_state.csv", steady_state_dict)


#
using SciMLNLSolve
damping = 0.1u"s"
scale_inertia = 1.1 # NOTE this parameter has changed
tconst = 0.01u"s" # NOTE this parameter has changed
network = import_system(:rtsgmlc; damping, scale_inertia, tconst = tconst)
zeroidx=1
(nd, p) = nd_model(network)
x0 = zeros(length(nd.syms));
x_static = solve(SteadyStateProblem(nd, x0, p), NLSolveJL())

θidx = idx_containing(nd, "θ")
offset = x_static[θidx[zeroidx]]
x_static[θidx] .= x_static[θidx] .- offset
@assert iszero(x_static[θidx[zeroidx]])

steady_state_dict = Dict(:SteadyState => x_static[θidx])
CSV.write("/home/brandner/.julia/dev/for_testing_delete_soon/RTS_new_stready_state_func_old_ND.csv", steady_state_dict)


###
### Plot trajectories
###

using CairoMakie
fig = Figure(resolution=(1800,1000))

# plot network

# plot all power flows
fig[1,1] = ax = Axis(fig; xlabel="time t in s", ylabel="apparent power flow in p.u.", title="power flow of all lines")
for i in 1:length(sol.load_S.saveval[1])
    t = sol.load_S.t
    y = [sol.load_S.saveval[t][i] for t in 1:length(sol.load_S.t)]
    lines!(ax, t, y; label="frequency on load node ($i)", linewidth=3)
end

fig[2,1] = ax = Axis(fig; xlabel="time t in s", ylabel="apparent power flow in p.u.", title="flow transients of failing lines")
for (i, l) in pairs(sol.failures.saveval)
    t, s = seriesforidx(sol.load_S, l)
    lines!(ax, t, s; label="flow on edge ($i)", linewidth=3)
    scatter!(ax, (t[end], s[end]); marker=:star5, markersize=25)
end
fig

# plot frequencies of all gen nodes
fig[3,1] = ax = Axis(fig; xlabel="time t in s", ylabel="angular frequency ω", title="frequency transients of all generator nodes")
(nd, p, overload_cb) = nd_model(network)
state_idx = idx_containing(nd, "ω") # array: indices of ω-states
for i in state_idx
    t = sol.sol.t
    y = [sol.sol.u[t][i] for t in 1:length(sol.sol.t)]
    lines!(ax, t, y; label="frequency on node ($i)", linewidth=3)
end
fig

# plot frequencies of failed gen nodes
fig[2,2] = ax = Axis(fig; xlabel="time t in s", ylabel="angular frequency ω", title="frequency transients of failing generator nodes")
(nd, p, overload_cb) = nd_model(network)
state_idx = idx_containing(nd, "ω") # array: indices of ω-states
node_idx = map(s -> parse(Int, String(s)[4:end]), nd.syms[state_idx]) # array: indices of vertices that are generators
for (i, l) in pairs([state_idx[findfirst(x -> x == i, node_idx)] for i in sol.failures_nodes.saveval])
    t, s = seriesforidx(sol.sol, l)
    lines!(ax, t, s; label="frequency on node ($i)", linewidth=3)
    scatter!(ax, (t[end], s[end]); marker=:star5, markersize=25)
end
vlines!(ax, tobs; color=:black, linewidth=1)

# plot frequencies of all load nodes
fig[1,3] = ax = Axis(fig; xlabel="time t in s", ylabel="angular frequency ω", title="frequency transients of all load nodes")
for i in 1:length(sol.frequencies_load_nodes.saveval[1])
    t = sol.frequencies_load_nodes.t
    # t = sol.frequencies_load_nodes.t[1:20]
    y = [sol.frequencies_load_nodes.saveval[t][i] for t in 1:length(sol.frequencies_load_nodes.t)]
    # y = [sol.frequencies_load_nodes.saveval[t][i] for t in 1:20]
    lines!(ax, t, y; label="frequency on load node ($i)", linewidth=3)
    # scatter!(ax, t, y; label="frequency on load node ($i)", markersize=5)
end
vlines!(ax, tobs; color=:black, linewidth=1)

# plot frequencies of failed load nodes
load_node_idxs = findall(x -> x==:load, get_prop(network, 1:nv(network), :type))
fig[2,3] = ax = Axis(fig; xlabel="time t in s", ylabel="angular frequency ω", title="frequency transients of failing load nodes")
for (i, l) in pairs([findfirst(x -> x == i, load_node_idxs) for i in sol.failures_load_nodes.saveval])
    t, s = seriesforidx(sol.frequencies_load_nodes, l)
    # scatter!(ax, t, s; label="frequency on load node ($i)", linewidth=3)
    lines!(ax, t, s; label="frequency on load node ($i)", linewidth=3)
    scatter!(ax, (t[end], s[end]); marker=:star5, markersize=25)
end
vlines!(ax, tobs; color=:black, linewidth=1)
fig

