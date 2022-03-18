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
using DataFrames
using CSV

# using CairoMakie

using GLMakie
GLMakie.activate!()

####
#### do some checks on the steady state
####
network = import_system(:rtsgmlc; damping = 0.1u"s", tconst = 0.01u"s")
n = describe_nodes(network)
u0 = steadystate(network; zeroidx = 13)

nd, = nd_model(network)
@assert all(isapprox.(u0[idx_containing(nd, "ω")], 0, atol = 1e-10))
diff = u0[idx_containing(nd, "θ")] .- ustrip.(u"rad", n.Va)
reldiff = diff ./ ustrip.(u"rad", n.Va)
reldiff[13] = 0
Makie.hist(diff)
Makie.hist(reldiff)

# the angle ranges are in the same ballpark at least
extrema(u0[idx_containing(nd, "θ")])
extrema(ustrip.(u"rad", n.Va))

network = import_system(:rtsgmlc; damping = 0.1u"s", tconst = 0.01u"s")
sol = simulate(network;
    initial_fail = Int[11],
    failtime = 0.1,
    trip_lines = :dynamic,
    tspan = (0.0, 100.0),
    solverargs = (;), verbose = true);
inspect_solution(sol)

nd, = nd_model(network)
nd.syms[2]


####
#### Actual simulation
####
function trip_all_lines()
    network = import_system(:rtsgmlc; damping = 0.1u"s", tconst = 0.01u"s")
    df = DataFrame(; initial = Int[], t = Vector{Float64}[], line = Vector{Int}[], type=Symbol[])
    for i in 1:ne(network)
        println("Simulate for line $i (dyn)")
        sol = simulate(network;
            initial_fail = Int[i],
            trip_lines = :dynamic,
            tspan = (0.0, 100.0),
            verbose = false)
        terminated(sol) || @warn("i = $i did not termiante early")
        push!(df, (; initial = i, t = sol.failures.t, line = sol.failures.saveval, type=:dynamic))

        println("Simulate for line $i (static)")
        sol = simulate(network;
            initial_fail = Int[i],
            trip_lines = :static,
            tspan = (0.0, 1000.0),
            verbose = false)
        terminated(sol) || @warn("i = $i did not termiante early")
        push!(df, (; initial = i, t = sol.failures.t, line = sol.failures.saveval, type=:static))
    end
    return df
end
res = trip_all_lines()
f = joinpath(RESULTS_DIR, "trip_all_lines_results_dyn_and_static.csv")
CSV.write(f, res)
# res = CSV.read(f, DataFrame)

# problems with early termination
probidx_dynamic = [11, 12, 49, 50, 63, 84, 85]
probidx_static = [11, 12, 49, 84]

#
res.ntriggered = length.(res.t) .- 1
res = res[res.ntriggered .> 0, :]
res.initial

init = [11, 12, 27, 29, 50, 63, 85, 98, 100]

#=
Failure on 27 is super! Keine statisch cascade, aber dnamische
mit sehr schönem beispiel wie die nachbarline dynamisch trippt sonst aber nicht
=#
sol = simulate(network;
    initial_fail = Int[27],
    trip_lines = :dynamic,
    solverargs = (;dtmax=0.01), verbose = true);
inspect_solution(sol)

sol = simulate(network;
    initial_fail = Int[27],
    trip_lines = :static,
    solverargs = (;), verbose = true);
inspect_solution(sol)
