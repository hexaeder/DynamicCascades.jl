using DynamicCascades
using Graphs
using NetworkDynamics
using DynamicCascades: swing_equation!, dynamic_load!, powerflow!, algebraic_load!
using OrdinaryDiffEq
using DiffEqCallbacks
using SteadyStateDiffEq
using GLMakie
using GLMakie.GeometryBasics
using Base: rad2deg, deg2rad
using Unitful
using MetaGraphs

# using CairoMakie
# CairoMakie.activate!()

orange = Makie.RGB([227, 114, 34]./255...)
gray = Makie.RGB([142, 144, 143]./255...)
cyan = Makie.RGB([0, 159, 218]./255...)
green = Makie.RGB([105, 146, 58]./255...)
blueish = Makie.RGB([124, 174, 175]./255...)

####
#### analysis part
####

aperiodicD(H, K) = sqrt(- 2*H * K / (2π*50u"1/s"))
int_ω(H, K) = sqrt(-2*2π*50u"1/s" * K / H)
int_period(H, K) = 2pi/int_ω(H, K)
calc_τ(H, K; r=0.05, f=10) = π*K/(2*f*log(r)*int_ω(H, K))

network = import_system(:rtsgmlc; damping = 0.01u"s", tconst = 0.1u"s")
DynamicCascades.set_coupling_sum!(network)
nodes = describe_nodes(network)

####
#### find the right τ
####
Hs = ustrip.(skipmissing(nodes.H))
Ks = nodes.Ksum[.!ismissing.(nodes.H)]
apDs = ustrip.(u"s", skipmissing(aperiodicD.(nodes.H, nodes.Ksum)))
τs = ustrip.(calc_τ.(Hs, Ks; r=0.01, f=10))
Makie.plot(Hs, apDs)
Makie.plot(Hs, Ks)
Makie.plot(Ks, apDs)
extrema(τ)

####
#### Plotting part
####
g = SimpleGraph(3)
add_edge!(g, 1, 2)
add_edge!(g, 1, 3)

function static_node!(dv, _, _, _, _)
    dv[1] = 0.0
end

fixed = ODEVertex(f=static_node!, dim=1, sym = [:θ])

edge = StaticEdge(f=powerflow!, dim=1, coupling=:antisymmetric)
dynload = ODEVertex(f=dynamic_load!, dim=1, sym=[:θ])
swing = ODEVertex(f=swing_equation!, dim=2, sym=[:θ, :ω])
nd = network_dynamics([fixed, dynload, swing], edge, g)

function simulate(; P_load, τ, H, D, K, step, tmax=1.0)
    P_load = ustrip(u"pu", P_load)
    τ = ustrip(u"s", τ)
    H = ustrip(u"MJ/MW", H)
    D = ustrip(u"s", D)
    K = ustrip(u"pu", K)

    node_p = [(NaN, NaN, NaN),  # static p
              (P_load, τ, NaN), # dyn laod p (power, timeconstant)
              (P_load, H, D)] # swing p (power, inertia, damping)
    edge_p = K
    p = (node_p, edge_p)

    theta0 = -asin(P_load/K)
    # @show rad2deg(theta0)
    @assert -π/4 < theta0 < π/4
    x0 = [0.0, theta0, theta0, 0.0]

    dx = similar(x0)
    nd(dx, x0, p, 0.0)
    isapprox(dx, zeros(4); atol=0.01) || error("dx not zero $dx")

    tstep = 0.1
    prob = ODEProblem(nd, x0, (0, tmax), p)
    cb = PresetTimeCallback(tstep, (integrator) -> integrator.u[1] = step)
    sol = solve(prob, Tsit5(); callback = cb, dtmax = 0.0001);

    # find the first crossing of swing
    t, setpoint = seriesforidx(sol, 1)
    tidx = findfirst(t->t>tstep, sol.t)
    Δdyn   = (seriesforidx(sol, 2)[2] .- setpoint .- theta0)[tidx:end]
    Δswing = (seriesforidx(sol, 3)[2] .- setpoint .- theta0)[tidx:end]

    tdyn    = Dict{Int, Float64}()
    tswing  = Dict{Int, Float64}()
    tfactor = Dict{Int, Float64}()
    # println("       dynload  swing   factor")
    for percent in 90:99
        factor = (100-percent)/100
        dynidx   = findfirst(x -> x > -factor*step, Δdyn)
        swingidx = findfirst(x -> x > -factor*step, Δswing)
        tdyn[percent] = isnothing(dynidx) ? Inf : t[tidx+dynidx-1] - tstep
        tswing[percent] = isnothing(swingidx) ? Inf : t[tidx+swingidx-1] - tstep
        tfactor[percent] = tdyn[percent]/tswing[percent]
        # println(" $percent%   ", round(tdyn[percent];digits=5),
        #         "  ", round(tswing[percent];digits=5),
        #         "  ", round(tfactor[percent];digits=5))
    end
    tupcross = Float64[]
    for i in eachindex(Δswing)[1:end-1]
        if Δswing[i] < 0 && Δswing[i+1] > 0
            tup = t[tidx - 1 + i]
            push!(tupcross, tup)
        end
    end

    return (;sol, step, tstep, tdyn, tswing, tfactor, tupcross)
end

# plot
function plot(ntup)
    sol = ntup.sol
    step = ntup.step
    tstep = ntup.tstep

    set_theme!(theme_minimal(), fontsize = 20)
    fig = Figure(resolution = (1000, 400))
    fig[1, 1] = ax = Axis(fig, xlabel = "time t in s", ylabel = "voltage angel in degree")
    t = [sol.t[begin], tstep, tstep, sol.t[end]]
    x = [0, 0, step, step] .+ sol[begin][2]
    alg = lines!(ax, t, rad2deg.(x), color = gray, linewidth = 5, label = "algebraic load")
    dyn = lines!(ax, seriesforidx(sol, 2; f=rad2deg)...; color = orange, linewidth = 5, label = "dynamic load")
    swing = lines!(ax, seriesforidx(sol, 3; f=rad2deg)...; color = blueish, linewidth = 5, label = "swing load")
    axislegend(ax; position = :rb, framevisible = false)

    # fig[1, 2] = ax2 = Axis(fig, xlabel = "time", ylabel = "angel difference in deg")
    # t, angle = seriesforidx(sol, 1)
    # lines!(ax2,t, rad2deg.(seriesforidx(sol, 3)[2].-angle); color = blueish, linewidth = 5, label = "swing load")
    # lines!(ax2, t, rad2deg.(seriesforidx(sol, 2)[2].-angle); color = orange, linewidth = 5, label = "dynamic load")

    fig
end
plot(;kwargs...) = plot(simulate(;kwargs...))

mean(abs.(nodes.P)) # so take P_load = -1.5?
mean(skipmissing(nodes.H)) # so take H = 10 MJ/MW?
mean(nodes.Ksum) # so take K = -70 pu?

using CairoMakie
fig = plot(P_load=-1, τ=0.1u"s", H=10u"MJ/MW", D=0.1u"s", K=-20, step=Base.deg2rad(1))
save(joinpath(PLOT_DIR, "loadmodels.pdf"), fig)
