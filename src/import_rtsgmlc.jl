"""
    import_system(:rtsgmlc)

Import the RTS-GMLC system as a MetaGraph.
"""
function import_system(::Val{:rtsgmlc}; losses=false, scale_inertia=1.0, kwargs...)
    @info "Import system RTS-GMLC"
    data = joinpath(DATA_DIR, "RTS-GMLC")
    bus_data    = CSV.read(joinpath(data,"bus.csv"), DataFrame)
    gen_data    = CSV.read(joinpath(data,"gen.csv"), DataFrame)
    branch_data = CSV.read(joinpath(data,"branch.csv"), DataFrame)

    sort!(bus_data, "Bus ID")
    N = nrow(bus_data)
    g = MetaGraph(N)

    baseP = 100u"MW"
    set_prop!(g, :Pbase, baseP)
    set_prop!(g, :NodeProps, [:n, :id, :type, :P, :Q, :H, :Vm, :Vbase,
                              :Va, :P_load, :P_inj, :Q_load, :Q_inj, :x, :y])
    set_prop!(g, :EdgeProps, [:src, :dst, :R, :X, :rating, :id])

    set_prop!(g, 1:N, :id, bus_data."Bus ID")
    set_prop!(g, 1:N, :Vm, bus_data."V Mag"u"pu")
    set_prop!(g, 1:N, :Va, bus_data."V Angle"u"°")
    set_prop!(g, 1:N, :Vbase, bus_data."BaseKV"u"kV")
    set_prop!(g, 1:N, :x, bus_data."lng")
    set_prop!(g, 1:N, :y, bus_data."lat")
    # bustypes = bustype.(bus_data."Bus Type")
    # set_prop!(g, 1:N, :type, bustypes)

    set_prop!(g, 1:N, :P_load, bus_data."MW Load"u"MW")
    set_prop!(g, 1:N, :Q_load, bus_data."MVAR Load"u"MW")

    for (n, busid) in enumerate(bus_data."Bus ID")
        P_load = get_prop(g, n, :P_load)
        Q_load = get_prop(g, n, :Q_load)
        # select attached generators
        controltype = bus_data."Bus Type"[n]
        generators = gen_data[gen_data."Bus ID".==busid, ["Category", "MW Inj", "MVAR Inj", "Inertia MJ/MW"]]
        if isempty(generators)
            P_inj = 0.0u"MW"
            Q_inj = 0.0u"MW"
        else # has generators
            P_inj = sum(generators."MW Inj")u"MW"
            Q_inj = sum(generators."MVAR Inj")u"MW"
            set_prop!(g, n, :P_inj, P_inj)
            set_prop!(g, n, :Q_inj, Q_inj)
        end

        if controltype ∈ ("PV", "Ref")
            if bus_data."Bus ID"[n] ∈ [114, 214, 314]
                type = :syncon
                set_prop!(g, n, :H, scale_inertia * 5u"MJ/MW")
            else
                type = :gen
                inertia = scale_inertia * sum(generators."Inertia MJ/MW")u"MJ/MW"
                @assert !iszero(P_inj) "Generator $n does not inject power?"
                @assert !iszero(inertia) "Generator $n does not have inertia?"
                set_prop!(g, n, :H, inertia)
            end
        elseif controltype == "PQ"
            type = :load
        else
            error("Found unkonown control type $controltype")
        end

        set_prop!(g, n, :type, type)
        set_prop!(g, n, :P, (P_inj - P_load) / baseP * u"pu")
        set_prop!(g, n, :Q, (Q_inj - Q_load) / baseP * u"pu")
    end

    # set the inertia for sync condenser
    # scidxs = findall(s -> !ismissing(s) && iszero(s), describe_nodes(g).H)
    # @assert bus_data."Bus ID"[scidxs] == [114, 214, 314]
    # set_prop!(g, scidxs, :H, 5u"MJ/MW")

    for row in eachrow(branch_data)
        src = findfirst(x->x==row."From Bus", bus_data."Bus ID")
        dst = findfirst(x->x==row."To Bus", bus_data."Bus ID")
        propertys = Dict(:R => row."R"u"pu",
                         :X => row."X"u"pu",
                         # :rating => 0.5*(row."STE Rating"u"MW"/baseP * u"pu"),
                         :rating => row."STE Rating"u"MW"/baseP * u"pu",
                         :id => row."UID")
        add_edge!(g, src, dst, propertys)
    end

    set_gmlc_pos_relaxed!(g)

    losses || lossless!(g)
    set_missing_props!(g; kwargs...)

    return g
end

export balance_power!, lossless!
function balance_power!(network)
    nodes = describe_nodes(network)
    imbalance = sum(nodes.P)
    if imbalance ≈ 0
        println("already balanced!")
        return network
    end
    genidx = findall(!ismissing, nodes.P_inj)

    relative_inj = nodes.P_inj[genidx] ./ sum(nodes.P_inj[genidx])

    newp = copy(nodes.P)
    newp[genidx] .-= relative_inj .* imbalance

    @assert isapprox(sum(newp), 0, atol=1e-8) "Could not balance power! Sum is $(sum(newp))"

    set_prop!(network, 1:nv(network), :P, newp)
end

function lossless!(network)
    balance_power!(network)
    set_prop!(network, edges(network), :R, zeros(ne(network)) * u"pu")
end

function set_missing_props!(network; damping = nothing, tconst = nothing)
    nodes = describe_nodes(network)
    if !isnothing(damping)
        idxs = findall(x -> x === :gen || x === :syncon, bustype.(nodes.type))
        set_prop!(network, idxs, :damping, damping)
    end
    if !isnothing(tconst)
        # idxs = findall(x -> x === :load, bustype.(nodes.type))
        # set_prop!(network, idxs, :timeconst, tconst)
        set_prop!(network, 1:nv(network), :timeconst, tconst)
    end
end

function set_gmlc_pos_relaxed!(g)
    # set location property based on rts gmlc data
    data = joinpath(DATA_DIR, "RTS-GMLC")
    bus_data = CSV.read(joinpath(data,"bus.csv"), DataFrame)
    x = bus_data."lng"
    y = bus_data."lat"
    xn = 10*(x .- minimum(x))./(maximum(x)-minimum(x)).-5
    yn = 10*(y .- minimum(y))./(maximum(y)-minimum(y)).-5
    pos2 = spring(g, initialpos=Point2f.(zip(xn, yn)))

    pos2[21] += Point2(-.3,.25)
    pos2[15] += Point2(-.5,.0)
    pos2[17] += Point2(-.3,-.3)
    pos2[24] += Point2(-.3,-.3)
    pos2[18] += Point2(-.2,-.0)

    set_prop!(g, 1:nv(g), :pos, pos2)
end
