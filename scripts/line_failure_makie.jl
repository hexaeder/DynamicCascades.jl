#!/usr/bin/env julia

using OrdinaryDiffEq
using DiffEqCallbacks
using NetworkDynamics
using IEEE_RTS_96
using GLMakie
using LightGraphs
using MakieGraphs
using Random
using Colors
using ColorSchemes
import GraphPlot

@info "Find steady state..."
gen_τ   = 1.0
slack_τ = 1.0
load_τ  = 0.1

rts96 = load_rts96()
(nd, p) = nd_model(rts96; gen_τ, slack_τ, load_τ);
x0 = zeros(rts96.num_of_buses + rts96.num_of_generators);
tspan = (0., 2000.);
prob = ODEProblem(nd, x0, tspan, p);
solinit = solve(prob, Rosenbrock23(), callback=TerminateSteadyState(), save_everystep=false);
if solinit[end] == tspan[2]
    @warn "Simulation time ended without reaching steady state!"
else
    @info "Found steadystate at t = $(solinit.t[end])..."
end
x_static = copy(solinit[end]);

initial_fail = [27]
failtime = 1.0

(nd, p, overload_cb) = nd_model(rts96; gen_τ, slack_τ, load_τ);
tspan = (0., 100.)
prob = ODEProblem(nd, copy(x_static), tspan, p);

line_failures = SavedValues(Float64, Int)
S_values = SavedValues(Float64, Vector{Float64})
P_values = SavedValues(Float64, Vector{Float64})
cbs = CallbackSet(initial_fail_cb(initial_fail, failtime),
                  overload_cb(load_S=S_values, load_P=P_values, failures=line_failures))

sol = solve(prob, Rosenbrock23(), callback=cbs);

##################
#     Plot       #
##################

fig = Figure(resolution = (2000, 1500))

state_idx = idx_containing(nd, "ω");
node_idx = map(s -> parse(Int, String(s)[4:end]), nd.syms[state_idx]);
Δωmax = maximum(map(x -> max(abs.(extrema(x[state_idx]))...), sol))

lsgrid = labelslidergrid!(
    fig,
    ["Time", "Δω scaling"],
    [S_values.t, range(0, stop=Δωmax, length=100)];
    formats = [x -> "$(round(x, digits = 3)) $s" for s in ["s", "hz"]])
fig[2, 1:2] = lsgrid.layout
set_close_to!(lsgrid.sliders[1], 0.0)
set_close_to!(lsgrid.sliders[2], 1.0)

t = Node{Float64}(0.0)
Δω = Node{Float64}(1.0)
connect!(t, lsgrid.sliders[1].value)
connect!(Δω, lsgrid.sliders[2].value)

nw_sublayout = fig[1,1] = GridLayout()

tg_ωload = Toggle(fig, active=true)
tg_absolute = Toggle(fig)
tg_activeP = Toggle(fig)
toggles = [tg_ωload Label(fig, "calculate ω load") tg_absolute Label(fig, "absolute flow") tg_activeP Label(fig, "Active Power")]

nw_sublayout[1, 1] = grid!(toggles, tellwidth = false)

state = @lift(sol($t))
load_S = @lift begin
    tidx = findfirst(x->x > $t, S_values.t)
    tidx!==nothing ? S_values.saveval[tidx-1] : S_values.saveval[end]
end
load_P = @lift begin
    tidx = findfirst(x->x > $t, P_values.t)
    tidx!==nothing ? P_values.saveval[tidx-1] : P_values.saveval[end]
end

# load, generator, slack
markers = [:rect, :utriangle, :circle];
nodemarker = Node(markers[rts96.bustype]);

sol(sol.t[end]+0.1)

nodescheme = ColorSchemes.delta
nodecolor = @lift begin
    # hacky way to caclulate ω from θ as finite differnece
    if $(tg_ωload.active)
        idx = idx_containing(nd, "θ")
        ω = (sol($t+0.1)[idx] - sol($t-0.1)[idx]) ./ 0.2
        relative = ω./(2 * $Δω ) .+ .5
        colors = [nodescheme[r] for r in relative]
    else
        colors = [ColorSchemes.RGB{Float64}(1,0,0) for i in 1:rts96.num_of_buses]
    end
    # overwrite with ω where there are "true" values in state
    for (n, s) in zip(node_idx, state_idx)
        relative = ($state)[s]/(2 * $Δω ) + .5
        colors[n] = nodescheme[relative]
    end
    colors
end

alltime_max_S = maximum(map(load -> maximum(abs.(load)), S_values.saveval))
alltime_max_P = maximum(map(load -> maximum(abs.(load)), P_values.saveval))

edgescheme = reverse(ColorSchemes.autumn1)
edgecolor = @lift begin
    load = $(tg_activeP.active) ? $load_P : $load_S
    if $(tg_absolute.active)
        max = $(tg_activeP.active) ? alltime_max_P : alltime_max_S
        colors = get(edgescheme, abs.(load) ./ max)
    else
        colors = get(edgescheme, abs.(load) ./ rts96.emergency_rating)
    end
    offline = @view colors[findall(iszero, $load_S)]
    offline .= ColorSchemes.RGB{Float64}(0,0,0)
    colors
end

# observables which hold the selected nodes and edges
sel_nodes = Node(Set{Int}())
sel_edges = Node(Set{Int}())

width_src, width_dst = 9.0, 3.0
widthdata = zeros(2*rts96.num_of_lines)
edgewidth = @lift begin
    # line segements need two width arguments for each line segement
    for idx in 1:rts96.num_of_lines
        if $load_P[idx] == 0
            w1 = w2 = (width_src + width_dst) / 2
        elseif $load_P[idx] > 0
            w1, w2 = width_src, width_dst
        else
            w1, w2 = width_dst, width_src
        end
        magnify = idx ∈ $sel_edges ? 2.5 : 1.0

        widthdata[2*idx-1] = w1 * magnify
        widthdata[2*idx]   = w2 * magnify
    end
    widthdata
end
(node_s, node_b) = (20.0, 40.0)
nodesize = @lift begin
    size = [node_s for i in 1:rts96.num_of_buses]
    for idx in $sel_nodes
        size[idx] = node_b
    end
    size
end


nwax = nw_sublayout[2,1] = Axis(fig)
nwax.aspect = AxisAspect(1)
hidedecorations!(nwax)
MakieGraphs.networkplot!(nwax, rts96.graph,
                         marker=nodemarker,
                         nodecolor=nodecolor,
                         nodesize=nodesize,
                         edgecolor=edgecolor,
                         edgewidth=edgewidth)

limits = @lift( (-$Δω, $Δω) )
nw_sublayout[3, 1] = Colorbar(fig, height=25, vertical=false, limits=limits, colormap=nodescheme, label="Node frequency")
nw_sublayout[4, 1] = Colorbar(fig, height=25, vertical=false, colormap=edgescheme, label="Line load_S")

# delete other interactions on scene
# empty!(nwax.interactions)
deregister_interaction!(nwax, :rectanglezoom)
# deregister_interaction!(nwax, :dragpan)
# deregister_interaction!(nwax, :limitreset)
# deregister_interaction!(nwax, :scrollzoom)

register_interaction!(nwax, :toggle_elements) do event::MouseEvent, axis
    if event.type === MouseEventTypes.leftclick
       (element, idx) = mouse_selection(axis.scene)
       if element isa LineSegments
           idx = Int(idx/2)
           action = idx ∈ sel_edges[] ? delete! : push!
           action(sel_edges[], idx)
           sel_edges[] = sel_edges[]
           println("Edge $idx toggled.")
           println(" ├ conductance = ", rts96.conductance[idx])
           println(" ├ susceptance = ", rts96.susceptance[idx])
           println(" ├ ap coupling = ", rts96.active_power_coupling[idx])
           println(" ├ rating      = ", rts96.emergency_rating[idx])
           println(" ├ load P      = ", load_P[][idx])
           println(" └ load S      = ", load_S[][idx])
       elseif element isa Scatter
           action = idx ∈ sel_nodes[] ? delete! : push!
           action(sel_nodes[], idx)
           sel_nodes[] = sel_nodes[]
           println("Node $idx toggled:")
           println(" ├ type = ", rts96.bustype[idx], " (1=load, 2=generator, 3=slack)")
           println(" ├ voltage  = ", rts96.voltage[idx])
           println(" ├ interita  = ", rts96.inertia[idx])
           println(" └ active p = ", rts96.active_power[idx])
       end
    end
end

sublayout = fig[1,2] = GridLayout()
θax = sublayout[1,1] = Axis(fig)
θax.title = "Phase angle θ"
θax.xzoomlock = true
ωax = sublayout[2,1] = Axis(fig)
ωax.title = "Frequency ω"
ωax.xzoomlock = true
fax = sublayout[3,1] = Axis(fig)
fax.title = "Powerflow"
fax.xzoomlock = true

islider = sublayout[4,1] = IntervalSlider(fig, range=range(S_values.t[begin], stop=S_values.t[end], length=500))

function set_x_lims(interval)
    xlims!(θax, interval)
    xlims!(ωax, interval)
    xlims!(fax, interval)
end
set_x_lims(islider.interval[]) # call first time
on(set_x_lims, islider.interval)

function plot_nodes(selected, ωload)
    empty!(θax); empty!(ωax)
    colors = copy(fig.scene.palette[:color][])
    for idx in selected
        c = pop!(colors)
        θidx = findfirst(sym -> occursin(Regex("θ_$idx\$"), String(sym)), nd.syms)
        if θidx !== nothing
            lines!(θax, sol.t, map(u->u[θidx], sol), label="θ_$idx", color=c, linewidth=3)
        end
        ωidx = findfirst(sym -> occursin(Regex("ω_$idx\$"), String(sym)), nd.syms)
        if ωidx !== nothing
            lines!(ωax, sol.t, map(u->u[ωidx], sol), label="ω_$idx", color=c, linewidth=3)
        elseif θidx !== nothing && ωload
            θs = map(u->u[θidx], sol)
            lines!(ωax, sol.t[begin:end-1], diff(θs) ./ diff(sol.t), label="ω_$idx", color=c, linewidth=3)
        end
    end
    vlines!(θax, t)
    vlines!(ωax, t)
end
plot_nodes([], tg_ωload.active[]) #plot first time
onany(plot_nodes, sel_nodes, tg_ωload.active)

function plot_edges(selected)
    empty!(fax)
    colors = copy(fig.scene.palette[:color][])
    for idx in selected
        c = pop!(colors)
        lines!(fax, S_values.t, map(l->abs(l[idx]), S_values.saveval), color=c, linewidth=3)
        lines!(fax, P_values.t, map(l->abs(l[idx]), P_values.saveval), color=c, linewidth=3, linestyle=:dash)
        rating = rts96.emergency_rating[idx]
        hlines!(fax, rating, color=c, linewidth=3)
    end
    vlines!(fax, t)
end
plot_edges([]) #plot first time
on(plot_edges, sel_edges)

tslider = lsgrid.sliders[1]

function set_time_interaction(event::MouseEvent, axis)
    if event.type === MouseEventTypes.leftclick
        t = mouseposition(axis.scene)[1]
        set_close_to!(tslider, t)
    end
end

register_interaction!(set_time_interaction, θax, :set_time)
register_interaction!(set_time_interaction, ωax, :set_time)
register_interaction!(set_time_interaction, fax, :set_time)

function move_left()
    newpos = tslider.selected_index[] - 1
    if newpos ≥ 1
        tslider.selected_index[] = newpos
    end
end
function move_right()
    newpos = tslider.selected_index[] + 1
    if newpos ≤ length(tslider.range[])
        tslider.selected_index[] = newpos
    end
end

on(fig.scene.events.keyboardbuttons) do button
    ispressed(button, Keyboard.left) && move_left()
    ispressed(button, Keyboard.right) && move_right()
end
