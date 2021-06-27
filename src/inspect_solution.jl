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
    network = c.network
    sol = c.sol
    S_values = c.load_S
    P_values = c.load_P

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

    load_S = @lift S_values($t)
    load_P = @lift P_values($t)

    emergency_rating = get_prop(network, edges(network), :rating)

    # observables which hold the selected nodes and edges
    sel_nodes = Node(Set{Int}())
    sel_edges = Node(Set{Int}())

    nwax = nw_sublayout[2,1] = Axis(fig)
    nwax.aspect = AxisAspect(1)
    hidedecorations!(nwax); hidespines!(nwax)
    gpargs = gparguments(c, t;
                         ωload=tg_ωload.active,
                         Δω,
                         activeP=tg_activeP.active,
                         colortype=edgemenu.selection,
                         load_S,
                         load_P,
                         sel_nodes,
                         sel_edges,
                         nd,
                         state_idx,
                         node_idx)
    nwplot = graphplot!(nwax, network; gpargs...)

    HOVER_DEFAULT = "Hover node/edge to see info!"
    description_text = Node{String}(HOVER_DEFAULT)
    nw_sublayout[3, 1] = Colorbar(fig, nwplot.plots[3], height=25, vertical=false, label="Node frequency")
    nw_sublayout[4, 1] = Colorbar(fig, height=25, vertical=false,
                                  colormap=@lift(edge_colorsheme($(edgemenu.selection))), label="Line load_S")
    nw_sublayout[2, 1] = Label(fig, description_text, tellwidth=false, tellheight=false, halign=:left, valign=:top)

    # delete other interactions on scene
    deregister_interaction!(nwax, :rectanglezoom)

    edgeclick = EdgeClickHandler() do idx, event, axis
        sel_edges[] = idx ∈ sel_edges[] ? delete!(sel_edges[], idx) : push!(sel_edges[], idx)
    end
    register_interaction!(nwax, :eclick, edgeclick)

    nodeclick = NodeClickHandler() do idx, event, axis
        sel_nodes[] = idx ∈ sel_nodes[] ? delete!(sel_nodes[], idx) : push!(sel_nodes[], idx)
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

    ####
    #### Plot graphs on the right side
    ####
    sublayout = fig[1,2] = GridLayout()

    tg_showP = Toggle(fig, active=false)
    tg_showRating = Toggle(fig, active=true)
    righttoggles = reshape([Label(fig, "show P: "), tg_showP,
                            Label(fig, "show rating: "), tg_showRating,
                            ], 1, :)
    sublayout[1,1] = grid!(righttoggles, tellwidth = false)

    θax = sublayout[2,1] = Axis(fig)
    θax.title = "Phase angle θ"
    θax.xzoomlock = true
    ωax = sublayout[3,1] = Axis(fig)
    ωax.title = "Frequency ω"
    ωax.xzoomlock = true
    fax = sublayout[4,1] = Axis(fig)
    fax.title = "Powerflow"
    fax.xzoomlock = true

    islider = sublayout[5,1] = IntervalSlider(fig, range=range(S_values.t[begin], stop=S_values.t[end], length=500))

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
            trange = range(sol.t[begin], sol.t[end], length=1000)
            if θidx !== nothing
                lines!(θax, trange, Float32[sol(t)[θidx] for t in trange], label="θ_$idx", color=c, linewidth=3)
            end
            ωidx = findfirst(sym -> occursin(Regex("ω_$idx\$"), String(sym)), nd.syms)
            if ωidx !== nothing
                lines!(ωax, trange, Float32[sol(t)[ωidx] for t in trange], label="ω_$idx", color=c, linewidth=3)
            elseif θidx !== nothing && ωload
                # TODO: prob use https://diffeq.sciml.ai/stable/basics/solution/
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

            (tS, S) = seriesforidx(S_values, idx)
            (tP, P) = seriesforidx(P_values, idx)
            lines!(fax, tS, S, color=c, linewidth=3)
            lines!(fax, tP, P, color=c, linewidth=3, linestyle=:dash, visible=tg_showP.active)

            rating = emergency_rating[idx]
            hlines!(fax, rating, color=c, linewidth=3, visible=tg_showRating.active)
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

export gparguments

gparguments(c::SolutionContainer, t::Number; kwargs...) = gparguments(c, Observable(t); kwargs...)

function gparguments(c::SolutionContainer, t::Observable;
                     ωload = Observable(true),
                     Δω = Observable(1.0),
                     activeP = Observable(false),
                     colortype = Observable(:relrating),
                     load_S = @lift(c.load_S($t)),
                     load_P = @lift(c.load_P($t)),
                     sel_nodes = Node(Set{Int}()),
                     sel_edges = Node(Set{Int}()),
                     nd=nd_model(c.network)[1],
                     state_idx = idx_containing(nd, "ω"),
                     node_idx = map(s -> parse(Int, String(s)[4:end]), nd.syms[state_idx]),
                     kwargs...)
    sol = c.sol
    network = c.network
    S_values = c.load_S
    P_values = c.load_P

    if colortype isa Symbol
        colortype = Observable(colortype)
    end

    state = @lift(sol($t))

    # load, generator, slack
    markers = Dict(:load => :rect, :gen => :utriangle, :slack => :circle);
    node_marker = [markers[k] for k in get_prop(network, vertices(network), :type)];

    node_color = @lift begin
        ω = [0.0 for i in 1:nv(network)]
        # hacky way to caclulate ω from θ as finite differnece
        if $ωload
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

    edgescheme = @lift edge_colorsheme($colortype)
    emergency_rating = get_prop(network, edges(network), :rating)
    edge_color = @lift begin
        load = $(activeP) ? $load_P : $load_S

        mode = $(colortype)
        if mode === :relrating
            cvals = abs.(load) ./ emergency_rating
        elseif mode === :abs
            max = $(activeP) ? alltime_max_P : alltime_max_S
            cvals = abs.(load) ./ max
        elseif mode === :abssteady
            max = $(activeP) ? alltime_max_P_diff : alltime_max_S_diff
            steady = $(activeP) ? P_steady : S_steady
            cvals = (load .- steady)./(2*max) .+ 0.5
        else # nothing
            cvals = zeros(length(load))
        end

        colors = get(edgescheme[], cvals)

        offline = @view colors[findall(iszero, $load_S)]
        offline .= ColorSchemes.RGB{Float64}(0,0,0)
        colors
    end

    width_src, width_dst = 9.0, 3.0
    widthdata = zeros(2*ne(network))

    iscairo = repr(typeof(Makie.current_backend[])) == "CairoMakie.CairoBackend"

    edge_width = @lift begin
        # line segements need two width arguments for each line segement
        for idx in 1:ne(network)
            if iscairo || $load_P[idx] == 0
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
    return (;layout=read_pos_or_spring,
            node_marker,
            node_color,
            node_size,
            node_attr=(colorrange=ω_range, colormap=node_colorscheme),
            edge_color,
            edge_width, kwargs...)
end

function edge_colorsheme(type)
    if type === :abssteady
        ColorScheme([colorant"yellow", colorant"gray", colorant"red"])
    else
        ColorScheme([colorant"yellow", colorant"red"])
    end
end

function GraphMakie.graphplot(c::SolutionContainer, t; kwargs...)
    graphplot(c.network; gparguments(c, t; kwargs...)...)
end

function GraphMakie.graphplot!(ax, c::SolutionContainer, t; kwargs...)
    graphplot!(ax, c.network; gparguments(c, t; kwargs...)...)
end
