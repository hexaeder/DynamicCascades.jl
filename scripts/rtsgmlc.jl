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
using LinearAlgebra
using GairoMakie

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
function trip_all_lines(;damping = 0.1u"s")
    network = import_system(:rtsgmlc; damping, tconst = 0.01u"s")
    df = DataFrame(; initial = Int[], t = Vector{Float64}[], line = Vector{Int}[], type=Symbol[], D=Float64[], hittlimit=Bool[])
    for i in 1:ne(network)
        println("Simulate for line $i (dyn)")
        sol = simulate(network;
            initial_fail = Int[i],
            trip_lines = :dynamic,
            tspan = (0.0, 100.0),
            verbose = false)
        terminated(sol) || @warn("i = $i did not termiante early")
        result = (; initial = i, t = sol.failures.t, line = sol.failures.saveval,
                  type=:dynamic, D=ustrip(damping), hittlimit=!terminated(sol))
        println(result)
        push!(df, result)

        println("Simulate for line $i (static)")
        sol = simulate(network;
            initial_fail = Int[i],
            trip_lines = :static,
            tspan = (0.0, 100.0),
            verbose = false)
        if !terminated(sol) && !all(sol.load_S.saveval[end] .< describe_edges(network).rating)
            @warn("i = $i did not termiante early and there where overloads")
            hittlimit = true
        else
            hittlimit = false
        end
        result = (; initial = i, t = sol.failures.t, line = sol.failures.saveval,
                  type=:static, D=ustrip(damping), hittlimit)
        println(result)
        push!(df, result)
    end

    f = joinpath(RESULTS_DIR, "trip_all_lines_d=$(ustrip(damping)).csv")
    CSV.write(f, df)
    return df
end
trip_all_lines(damping=0.01u"s")
trip_all_lines(damping=0.03u"s")
trip_all_lines(damping=0.1u"s")
trip_all_lines(damping=0.3u"s")
trip_all_lines(damping=1.0u"s")
exit()

####
#### Auswertung
####
# type piracy :(
Base.tryparse(::Type{Vector{T}}, str) where T = map(x -> parse(T, x), split(str[2:end-1], ','))

df = DataFrame()
for f in filter(contains("trip_all_lines"), readdir(RESULTS_DIR, join=true))
    add = CSV.read(f, DataFrame; types=Dict(:t=>Vector{Float64},
                                            :line=>Vector{Int},
                                            :type=>Symbol))
    df = vcat(df, add)
end
# set the triggered argument
df.ntriggered = length.(df.t) .- 1;

df = df[df.ntriggered .> 0, :]

#=
Check, if the static cascades are independent from D
=#
@testset "static independent from D" begin
    Ds = unique(df.D)
    statics = df[df.type .== :static, :]
    statics.waves = map(statics.t) do ts
        x = Vector{Int}(undef, length(ts))
        x[1] = 1
        for i in 2:length(ts)
            x[i] = ts[i-1]==ts[i] ? x[i-1] : x[i-1]+1
        end
    end

    idxs = unique(statics.initial)
    for i in idxs
        @test size(unique!(statics[statics.initial .== i, [:waves, :line]])) == (1,2)
    end
end

#=
Check out the dependency of the dynamic cascades to D
=#
dyns = sort!(df[df.type .== :dynamic, Not(:type)], :initial)

hist(dyns[dyns.D .== 0.01, :ntriggered])
hist!(dyns[dyns.D .== 0.03, :ntriggered])
hist!(dyns[dyns.D .== 0.1, :ntriggered])
hist!(dyns[dyns.D .== 0.3, :ntriggered])

hist!(dyns[dyns.D .== 0.3, :ntriggered])


A = res[res.type.==:dynamic, :].ntriggered;
B = res[res.type.==:static, :].ntriggered;
pie_dynamic = [(i, count(==(i), A)) for i in unique(A)]
pie_static = [(i, count(==(i), B)) for i in unique(B)]

res[res.ntriggered .> 0, :]

#=
Failure on 27 is super! Keine statisch cascade, aber dnamische
mit sehr schönem beispiel wie die nachbarline dynamisch trippt sonst aber nicht
=#
init = 27
damping = 0.1u"s"
network = import_system(:rtsgmlc; damping, tconst = 0.01u"s")
sol = simulate(network;
    initial_fail = Int[init],
    trip_lines = :dynamic,
    solverargs = (;dtmax=0.01), verbose = true);
GLMakie.activate!()
inspect_solution(sol)

times = [0.0, sol.failures.t..., 6.0]

function plotnetwork(fig, sol, t; line=nothing, offset=nothing, text=nothing, apos=.5)
    ax = Axis(fig)
    hidedecorations!(ax), hidespines!(ax)
    xlims!(ax, -10, 7)
    ylims!(ax, -5, 7)
    gpargs = gparguments(sol, Observable(t);
                         colortype=:abssteady,
                         Δω=Observable(3.0),
                         offlinewidth=5,
                         offlinecolor=colorant"lightgray",
                         ecolorscaling = Observable(1.0),
                         node_size=15)
    p = graphplot!(ax, network; gpargs...)

    if line !==nothing
        pos = GraphMakie.interpolate(get_edge_plot(p).paths[][line], apos)
        offset = Point2(offset)
        gapA = 0.0
        gapB = 0.6
        posA = pos + gapA*normalize(offset)
        posT = pos + offset
        posB = posT - gapB*normalize(offset)
        lines!(ax, [posA, posB], color=:black, linewidth=2)
        text!(ax, text; position=posT, align=(:center, :center))
    end
    return ax
end

CairoMakie.activate!()
set_theme!(theme_minimal(), fontsize = 20)
fig = Figure(resolution=(1500,1500))

fig[1,1] = nwax = plotnetwork(fig, sol, times[1]; line=11, offset=(-4,-.5), text="(3)")
fig[1,1] = Label(fig, "(a) t = $(round(times[1],digits=2)) s", tellwidth=false, tellheight=false, halign=:left, valign=:top)

fig[1,2] = plotnetwork(fig, sol, times[2]; line=27, offset=(0.75,2.5), text="(1)")
fig[1,2] = Label(fig, "(b) t = $(round(times[2],digits=2)) s", tellwidth=false, tellheight=false, halign=:left, valign=:top)

fig[1,3] = plotnetwork(fig, sol, times[3]; line=29, offset=(-2,2.5), text="(2)",apos=.4)
fig[1,3] = Label(fig, "(c) t = $(round(times[3],digits=2)) s", tellwidth=false, tellheight=false, halign=:left, valign=:top)

fig[2,1] = plotnetwork(fig, sol, times[4]; line=11, offset=(-4,-.5), text="(3)")
fig[2,1] = Label(fig, "(d) t = $(round(times[4],digits=2)) s", tellwidth=false, tellheight=false, halign=:left, valign=:top)

fig[2,2] = plotnetwork(fig, sol, times[5]; line=12, offset=(-4,0), text="(4)", apos=.25)
fig[2,2] = Label(fig, "(e) t = $(round(times[5],digits=2)) s", tellwidth=false, tellheight=false, halign=:left, valign=:top)

fig[2,3] = plotnetwork(fig, sol, times[6]; line=37, offset=(1.75,1.75), text="(5)")
fig[2,3] = Label(fig, "(f) t = $(round(times[6],digits=2)) s", tellwidth=false, tellheight=false, halign=:left, valign=:top)

fig[3,1] = plotnetwork(fig, sol, times[7])
fig[3,1] = Label(fig, "(g) t → ∞ (steady state)", tellwidth=false, tellheight=false, halign=:left, valign=:top)

fig[3,2:3] = ax = Axis(fig; xlabel="time t in s", ylabel="apparent power flow in p.u.", title="flow transients of failing lines")
for (i, l) in pairs(sol.failures.saveval)
    t, s = seriesforidx(sol.load_S, l)
    lines!(ax, t, s; label="flow on edge ($i)", linewidth=5)
end
for (i, l) in pairs(sol.failures.saveval)
    t, s = seriesforidx(sol.load_S, l)
    scatter!(ax, (t[end], s[end]); marker=:star5, markersize=35)
end
axislegend(ax; position=:lt)

p = nwax.scene.plots[2]
@assert p isa GraphPlot

fig[0,:] = Colorbar(fig, get_node_plot(p), height=25, vertical=false, label="node frequency deviation in hz")

fig[0, :] = Colorbar(fig, height=25, vertical=false,
                     colormap=Makie.ColorScheme([colorant"yellow", colorant"red"]), label="line load relative to rating")

save(joinpath(PLOT_DIR, "rts_cascade.pdf"), fig)
