using DynamicCascades
using Graphs
using NetworkDynamics
using DynamicCascades: swing_equation!, dynamic_load!, powerflow!, algebraic_load!
using OrdinaryDiffEq
using DiffEqCallbacks
using SteadyStateDiffEq
using GLMakie
using GLMakie.GeometryBasics

using CairoMakie
CairoMakie.activate!()

DIR = "/Users/hw/MA/Forschungsbeleg/figures/"

orange = Makie.RGB([227, 114, 34]./255...)
gray = Makie.RGB([142, 144, 143]./255...)
cyan = Makie.RGB([0, 159, 218]./255...)
green = Makie.RGB([105, 146, 58]./255...)
blueish = Makie.RGB([124, 174, 175]./255...)

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
node_p = [(NaN, NaN, NaN),  # static p
          (-1.0, NaN, 0.1), # dyn laod p (power, _, timeconstant)
          (-1.0, 1.0, 0.3)] # swing p (power, inertia, damping)
edge_p = -sqrt(2)
p = (node_p, edge_p)
theta0 = -π/4
x0 = [0.0,        # static
      theta0,       # dynload
      theta0, 0.0]  # swing


dx = similar(x0)
nd(dx, x0, p, 0.0)
@assert dx ≈ zeros(4)

step = 1
tstart = 3.0
prob = ODEProblem(nd, x0, (0, 25.), p)
cb = PresetTimeCallback(tstart, (integrator) -> integrator.u[1] = step)
sol = solve(prob, Tsit5(); callback=cb, dtmax=0.05)

# plot
set_theme!(theme_minimal(), fontsize=20)
fig = Figure(resolution=(1200,600))
fig[1,1] = ax1 = Axis(fig, aspect=DataAspect())
fig[1,2] = ax2 = Axis(fig, xlabel="time", ylabel="voltage angel")

θidx = [1,2,1,3]
tnode = Observable(sol.t[begin])

s1 = @lift(thetashape(sol($tnode)[1]; loc=(0,2)))
s2 = @lift(thetashape(sol($tnode)[2]; loc=(4,2)))
s3 = @lift(thetashape(sol($tnode)[1]; loc=(0,0)))
s4 = @lift(thetashape(sol($tnode)[3]; loc=(4,0)))

poly!(ax1, s1, color=gray)
poly!(ax1, s2, color=orange)
poly!(ax1, s3, color=gray)
poly!(ax1, s4, color=blueish)
lines!(ax1, [(0.6,2),(3.4,2)], color=gray, linewidth=5)
lines!(ax1, [(0.6,0),(3.4,0)], color=gray, linewidth=5)
hidedecorations!(ax1); hidespines!(ax1)

text!(ax1, ["generator", "dynamic load", "swing load"]; align=(:center,:center), position=[(0,1), (4,2.7), (4,-0.7)])

t = [0.0, tstart, tstart, 25.0]
x = [theta0, theta0, theta0+step, theta0+step]
alg = lines!(ax2,t, x, color=gray, linewidth=5, label="algebraic load")
dyn = lines!(ax2,seriesforidx(sol, 2)...; color=orange, linewidth=5, label="dynamic load")
swing = lines!(ax2,seriesforidx(sol, 3)...; color=blueish, linewidth=5, label="swing load")
vlines!(ax2, tnode, color=:gray, linewidth=3, visible=@lift($tnode>0))
axislegend(ax2; position=:rb, framevisible=false)

time, fps = 10, 30
tspan = range(sol.t[begin], sol.t[end], length=fps*time)
save(joinpath(DIR, "../videos", "swing_and_dyn.png"), fig)
record(fig, joinpath(DIR, "../videos", "swing_and_dyn.mp4"), tspan) do time
    tnode[] = time
end
