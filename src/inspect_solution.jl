using GLMakie
using LightGraphs
using MetaGraphs
using GraphMakie
using Colors
using ColorSchemes
using NetworkDynamics
using NetworkLayout: Spring

export inspect_solution

function inspect_solution(network, nd, sol, S_values, P_values)
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
    # t = lsgrid.sliders[1].value
    # Δω = lsgrid.sliders[2].value

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
    markers = Dict(:load => :rect, :gen => :utriangle, :slack => :circle);
    node_marker = [markers[k] for k in get_prop(network, vertices(network), :type)];

    # nodescheme = ColorSchemes.delta
    node_color = @lift begin
        ω = [0.0 for i in 1:nv(network)]
        # hacky way to caclulate ω from θ as finite differnece
        if $(tg_ωload.active)
            idx = idx_containing(nd, "θ")
            ω = (sol($t+0.1)[idx] - sol($t-0.1)[idx]) ./ 0.2
        end
        # overwrite with ω where there are "true" values in state
        for (n, s) in zip(node_idx, state_idx)
            ω[n] = ($state)[s]
        end
        ω
    end
    node_colorscheme = ColorSchemes.delta
    ω_range = @lift (-$Δω, $Δω)

    alltime_max_S = maximum(map(load -> maximum(abs.(load)), S_values.saveval))
    alltime_max_P = maximum(map(load -> maximum(abs.(load)), P_values.saveval))
    edgescheme = reverse(ColorSchemes.autumn1)
    emergency_rating = get_prop(network, edges(network), :rating)
    edge_color = @lift begin
        load = $(tg_activeP.active) ? $load_P : $load_S
        if $(tg_absolute.active)
            max = $(tg_activeP.active) ? alltime_max_P : alltime_max_S
            colors = get(edgescheme, abs.(load) ./ max)
        else
            colors = get(edgescheme, abs.(load) ./ emergency_rating)
        end
        offline = @view colors[findall(iszero, $load_S)]
        offline .= ColorSchemes.RGB{Float64}(0,0,0)
        colors
    end

    # observables which hold the selected nodes and edges
    sel_nodes = Node(Set{Int}())
    sel_edges = Node(Set{Int}())

    width_src, width_dst = 9.0, 3.0
    widthdata = zeros(2*ne(network))
    edge_width = @lift begin
        # line segements need two width arguments for each line segement
        for idx in 1:ne(network)
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
    (node_small, node_big) = (20.0, 40.0)
    node_size = @lift begin
        size = [node_small for i in 1:nv(network)]
        for idx in $sel_nodes
            size[idx] = node_big
        end
        size
    end

    nwax = nw_sublayout[2,1] = Axis(fig)
    nwax.aspect = AxisAspect(1)
    hidedecorations!(nwax)
    nwplot = graphplot!(nwax, network;
                        layout=Spring(),
                        node_marker,
                        node_color,
                        node_size,
                        node_attr=(colorrange=ω_range, colormap=node_colorscheme),
                        edge_color,
                        edge_width)

    nw_sublayout[3, 1] = Colorbar(fig, nwplot.plots[3], height=25, vertical=false, label="Node frequency")
    nw_sublayout[4, 1] = Colorbar(fig, height=25, vertical=false, colormap=edgescheme, label="Line load_S")

    # delete other interactions on scene
    deregister_interaction!(nwax, :rectanglezoom)

    edgeclick = EdgeClickHandler() do idx, even, axis
        sel_edges[] = idx ∈ sel_edges[] ? delete!(sel_edges[], idx) : push!(sel_edges[], idx)
        println("Edge $idx toggled.")
        # println(" ├ conductance = ", network.conductance[idx])
        # println(" ├ susceptance = ", network.susceptance[idx])
        # println(" ├ ap coupling = ", network.active_power_coupling[idx])
        println(" ├ rating      = ", emergency_rating[idx])
        println(" ├ load P      = ", load_P[][idx])
        println(" └ load S      = ", load_S[][idx])
    end
    register_interaction!(nwax, :eclick, edgeclick)

    nodeclick = NodeClickHandler() do idx, even, axis
        sel_nodes[] = idx ∈ sel_nodes[] ? delete!(sel_nodes[], idx) : push!(sel_nodes[], idx)
        println("Node $idx toggled:")
        println(" ├ type = ", get_prop(network, idx, :type), " (1=load, 2=generator, 3=slack)")
        println(" ├ voltage  = ", get_prop(network, idx, :Vm))
        println(" ├ interita  = ", has_prop(network, idx, :inertia) ? get_prop(network, idx, :inertia) : "missing")
        println(" ├ active P = ", get_prop(network, idx, :P))
        println(" └ reactive Q = ", get_prop(network, idx, :Q))
    end
    register_interaction!(nwax, :nclick, nodeclick)

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

    on(islider.interval) do interval
        xlims!(θax, interval)
        xlims!(ωax, interval)
        xlims!(fax, interval)
    end

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
        vlines!(θax, t, color=:black)
        vlines!(ωax, t, color=:black)
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
            rating = emergency_rating[idx]
            hlines!(fax, rating, color=c, linewidth=3)
        end
        vlines!(fax, t, color=:black)
    end
    plot_edges([]) #plot first time
    on(plot_edges, sel_edges)

    tslider = lsgrid.sliders[1]

    function set_time_interaction(event::MouseEvent, axis)
        if event.type === MouseEventTypes.leftclick
            t = mouseposition(axis.scene)[1]
            set_close_to!(tslider, t)
            return true
        end
        return false
    end

    register_interaction!(set_time_interaction, θax, :set_time)
    register_interaction!(set_time_interaction, ωax, :set_time)
    register_interaction!(set_time_interaction, fax, :set_time)

    on(fig.scene.events.keyboardbutton) do e
        if e.action == Keyboard.press || e.action == Keyboard.repeat
            if e.key == Keyboard.left
                dec_slider(tslider)
                return true
            elseif e.key == Keyboard.right
                inc_slider(tslider)
                return true
            end
        end
        return false
    end

    return fig
end

function inc_slider(slider)
    newpos = slider.selected_index[] + 1
    if newpos ≤ length(slider.range[])
        slider.selected_index[] = newpos
    end
end

function dec_slider(slider)
    newpos = slider.selected_index[] - 1
    if newpos ≥ 1
        slider.selected_index[] = newpos
    end
end
