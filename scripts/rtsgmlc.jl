using DynamicCascades
using Graphs
using MetaGraphs
using Statistics
using GraphMakie
using Colors
using DynamicCascades: PLOT_DIR
using NetworkDynamics
using Unitful
using Unitful: @u_str

# using CairoMakie

using GLMakie
GLMakie.activate!()

####
#### do some checks on the steady state
####
network = import_system(:rtsgmlc; damping = 0.1u"s", tconst = 0.01u"s")
n = describe_nodes(network)
u0 = steadystate(network; zeroidx=13)

nd, = nd_model(network)
@assert all(isapprox.(u0[idx_containing(nd, "ω")], 0, atol = 1e-10))
diff = u0[idx_containing(nd, "θ")] .- ustrip.(u"rad", n.Va)
reldiff = diff./ustrip.(u"rad", n.Va)
reldiff[13] = 0
Makie.hist(diff)
Makie.hist(reldiff)

# the angle ranges are in the same ballpark at least
extrema(u0[idx_containing(nd, "θ")])
extrema(ustrip.(u"rad", n.Va))

####
#### Actual simulation
####
network = import_system(:rtsgmlc; damping = 0.1u"s", tconst = 0.01u"s")
sol = simulate(network;
    initial_fail = Int[27],
    failtime = 0.1,
    trip_lines = true,
    tspan = (0.0, 40.0),
    solverargs = (; dtmax = 0.01));

inspect_solution(sol)

plot_failing_lines(sol)

c = sol
maxt = maximum(c.failures.t)
ln = c.failures.saveval


vlines!(fax, lsgrid.sliders[1].value, color = :black)




gparguments(sol, Observable(0.0))

sol.load_S.saveval
sol.load_P.saveval

sol.sol
