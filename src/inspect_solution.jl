using GLMakie
using LightGraphs
using MetaGraphs
using GraphMakie
using Colors
using ColorSchemes
using NetworkDynamics
using NetworkLayout: spring

export inspect_solution

function inspect_solution(c::SolutionContainer)
    (nd,) = nd_model(c.network)
    inspect_solution(c.network, nd, c.sol, c.load_S, c.load_P)
end

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
    set_close_to!(lsgrid.sliders[2], 0.5)

    t = lsgrid.sliders[1].value
    Δω = lsgrid.sliders[2].value

    nw_sublayout = fig[1,1] = GridLayout()

    tg_ωload = Toggle(fig, active=true)
    tg_activeP = Toggle(fig)
    edgemenu = Menu(fig, options = [("rel to rating", :relrating),
                                    ("diff to steady", :abssteady),
                                    ("absolut", :abs)],
                    i_selected=1,
                    selection=:relrating)
    toggles = reshape([Label(fig, "calculate ω load:"), tg_ωload,
                       Label(fig, "Active Power:"), tg_activeP,
                       Label(fig, "Edge colors:"), edgemenu], 1, :)

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
    node_colorscheme = ColorSchemes.diverging_bkr_55_10_c35_n256
    ω_range = @lift (-$Δω, $Δω)

    alltime_max_S = maximum(map(load -> maximum(abs.(load)), S_values.saveval))
    alltime_max_P = maximum(map(load -> maximum(abs.(load)), P_values.saveval))

    S_steady = S_values.saveval[begin]
    P_steady = P_values.saveval[begin]
    alltime_max_S_diff = maximum(map(load -> maximum(abs.(extrema(load-S_steady))), S_values.saveval))
    alltime_max_P_diff = maximum(map(load -> maximum(abs.(extrema(load-P_steady))), P_values.saveval))

    edgescheme = @lift begin
        if $(edgemenu.selection) === :abssteady
            # ColorSchemes.linear_tritanopic_krjcw_5_95_c24_n256
            ColorScheme([colorant"yellow", colorant"gray", colorant"red"])
        else
            # ColorSchemes.linear_wyor_100_45_c55_n256
            ColorScheme([colorant"yellow", colorant"red"])
            # reverse(ColorSchemes.autumn1)
        end
    end
    emergency_rating = get_prop(network, edges(network), :rating)
    edge_color = @lift begin
        load = $(tg_activeP.active) ? $load_P : $load_S

        mode = $(edgemenu.selection)
        if mode === :relrating
            cvals = abs.(load) ./ emergency_rating
        elseif mode === :abs
            max = $(tg_activeP.active) ? alltime_max_P : alltime_max_S
            cvals = abs.(load) ./ max
        elseif mode === :abssteady
            max = $(tg_activeP.active) ? alltime_max_P_diff : alltime_max_S_diff
            steady = $(tg_activeP.active) ? P_steady : S_steady
            cvals = (load .- steady)./(2*max) .+ 0.5
        else # nothing
            cvals = zeros(length(load))
        end

        colors = get(edgescheme[], cvals)

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
    hidedecorations!(nwax); hidespines!(nwax)
    nwplot = graphplot!(nwax, network;
                        layout=read_pos_or_spring,
                        node_marker,
                        node_color,
                        node_size,
                        node_attr=(colorrange=ω_range, colormap=node_colorscheme),
                        edge_color,
                        edge_width)

    HOVER_DEFAULT = "Hover node/edge to see info!"
    description_text = Node{String}(HOVER_DEFAULT)
    nw_sublayout[3, 1] = Colorbar(fig, nwplot.plots[3], height=25, vertical=false, label="Node frequency")
    nw_sublayout[4, 1] = Colorbar(fig, height=25, vertical=false, colormap=edgescheme, label="Line load_S")
    nw_sublayout[2, 1] = Label(fig, description_text, tellwidth=false, tellheight=false, halign=:left, valign=:top)

    # delete other interactions on scene
    deregister_interaction!(nwax, :rectanglezoom)

    edgeclick = EdgeClickHandler() do idx, event, axis
        sel_edges[] = idx ∈ sel_edges[] ? delete!(sel_edges[], idx) : push!(sel_edges[], idx)
        # println("Edge $idx toggled.")
        # # println(" ├ conductance = ", network.conductance[idx])
        # # println(" ├ susceptance = ", network.susceptance[idx])
        # # println(" ├ ap coupling = ", network.active_power_coupling[idx])
        # println(" ├ rating      = ", emergency_rating[idx])
        # println(" ├ load P      = ", load_P[][idx])
        # println(" └ load S      = ", load_S[][idx])
    end
    register_interaction!(nwax, :eclick, edgeclick)

    nodeclick = NodeClickHandler() do idx, event, axis
        sel_nodes[] = idx ∈ sel_nodes[] ? delete!(sel_nodes[], idx) : push!(sel_nodes[], idx)
        # println("Node $idx toggled:")
        # println(" ├ type = ", get_prop(network, idx, :type), " (1=load, 2=generator, 3=slack)")
        # println(" ├ voltage  = ", get_prop(network, idx, :Vm))
        # println(" ├ interita  = ", has_prop(network, idx, :inertia) ? get_prop(network, idx, :inertia) : "missing")
        # println(" ├ active P = ", get_prop(network, idx, :P))
        # println(" └ reactive Q = ", get_prop(network, idx, :Q))
    end
    register_interaction!(nwax, :nclick, nodeclick)

    nodehover = NodeHoverHandler() do state, idx, event, axis
        if state
            string = """Node $idx
                        ├ type = $(get_prop(network, idx, :type))
                        ├ voltage  = $(get_prop(network, idx, :Vm))
                        ├ interita  = $(has_prop(network, idx, :inertia) ? get_prop(network, idx, :inertia) : "missing")
                        ├ active P = $(get_prop(network, idx, :P))
                        └ reactive Q = $(get_prop(network, idx, :Q))
                     """
        else
            string = HOVER_DEFAULT
        end
        description_text[] = string
    end
    register_interaction!(nwax, :nhover, nodehover)

    edgehover = EdgeHoverHandler() do state, idx, event, axis
        if state
            string = """Edge $idx
                        ├ coupling = $(get_prop(network, collect(edges(network))[idx], :_K))
                        ├ rating  = $(emergency_rating[idx])
                        ├ load P = $(load_P[][idx])
                        └ load S = $(load_S[][idx])
                     """
        else
            string = HOVER_DEFAULT
        end
        description_text[] = string
    end
    register_interaction!(nwax, :ehover, edgehover)

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

    plot_nodes = (selected, ωload) -> begin
        empty!(θax); empty!(ωax)
        colors = reverse(fig.scene.palette[:color][])
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
        vlines!(θax, lsgrid.sliders[1].value, color=:black)
        vlines!(ωax, lsgrid.sliders[1].value, color=:black)
    end

    plot_nodes([], tg_ωload.active[]) #plot first time
    onany(plot_nodes, sel_nodes, tg_ωload.active)

    plot_edges = selected -> begin
        empty!(fax)
        colors = reverse(fig.scene.palette[:color][])
        for idx in selected
            c = pop!(colors)
            lines!(fax, S_values.t, map(l->abs(l[idx]), S_values.saveval), color=c, linewidth=3)
            lines!(fax, P_values.t, map(l->abs(l[idx]), P_values.saveval), color=c, linewidth=3, linestyle=:dash)
            rating = emergency_rating[idx]
            hlines!(fax, rating, color=c, linewidth=3)
        end
        vlines!(fax, lsgrid.sliders[1].value, color=:black)
    end
    plot_edges([]) #plot first time
    on(plot_edges, sel_edges)

    tslider = lsgrid.sliders[1]

    set_time_interaction = (event::MouseEvent, axis) -> begin
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
    idx = slider.selected_index[]
    val = slider.range[][idx]
    while idx < length(slider.range[]) && slider.range[][idx] <= val
        idx += 1
    end
    set_close_to!(slider, slider.range[][idx])
end

function dec_slider(slider)
    idx = slider.selected_index[]
    val = slider.range[][idx]
    while idx > 1 && slider.range[][idx] >= val
        idx -= 1
    end
    set_close_to!(slider, slider.range[][idx])
end

function read_pos_or_spring(g::MetaGraph)
    pos = get_prop(g, 1:nv(g), :pos)
    if any(ismissing.(pos))
        return spring(g)
    else
        return pos
    end
end
