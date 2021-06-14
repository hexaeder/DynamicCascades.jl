"""
    import_system(:rtsgmlc)

Import the RTS-GMLC system as a MetaGraph.
"""
function import_system(::Val{:rtsgmlc})
    @info "Import system RTS-GMLC"
    data = joinpath(DATA_DIR, "RTS-GMLC")
    bus_data    = CSV.read(joinpath(data,"bus.csv"), DataFrame)
    gen_data    = CSV.read(joinpath(data,"gen.csv"), DataFrame)
    branch_data = CSV.read(joinpath(data,"branch.csv"), DataFrame)

    sort!(bus_data, "Bus ID")
    N = nrow(bus_data)
    g = MetaGraph(N)

    baseP = 100
    set_prop!(g, :Pbase, baseP)
    set_prop!(g, :NodeProps, [:n, :id, :type, :P, :Q, :inertia, :Vm, :Vbase,
                              :Va, :P_load, :P_inj, :Q_load, :Q_inj, :x, :y])
    set_prop!(g, :EdgeProps, [:src, :dst, :R, :X, :rating, :id])

    set_prop!(g, 1:N, :id, bus_data."Bus ID")
    set_prop!(g, 1:N, :Vm, bus_data."V Mag")
    set_prop!(g, 1:N, :Va, bus_data."V Angle")
    set_prop!(g, 1:N, :Vbase, bus_data."BaseKV")
    set_prop!(g, 1:N, :x, bus_data."lng")
    set_prop!(g, 1:N, :y, bus_data."lat")
    set_prop!(g, 1:N, :type, bustype.(bus_data."Bus Type"))

    set_prop!(g, 1:N, :P_load, bus_data."MW Load"/baseP)
    set_prop!(g, 1:N, :Q_load, bus_data."MVAR Load"/baseP)

    for (n, busid) in enumerate(bus_data."Bus ID")
        # select attached generators
        generators = gen_data[gen_data."Bus ID".== busid, ["Category", "MW Inj", "MVAR Inj", "Inertia MJ/MW"]]
        P_inj   = isempty(generators) ? 0.0 : sum(generators."MW Inj") / baseP
        Q_inj   = isempty(generators) ? 0.0 : sum(generators."MVAR Inj") / baseP
        inertia = isempty(generators) ? 0.0 : sum(generators."Inertia MJ/MW")
        if P_inj!=0.0 || Q_inj!=0.0 || inertia!=0.0
            set_prop!(g, n, :P_inj, P_inj)
            set_prop!(g, n, :Q_inj, Q_inj)
            set_prop!(g, n, :inertia, inertia)
        end
        P_load = get_prop(g, n, :P_load)
        Q_load = get_prop(g, n, :Q_load)
        set_prop!(g, n, :P, P_inj - P_load)
        set_prop!(g, n, :Q, Q_inj - Q_load)
    end

    for row in eachrow(branch_data)
        src = findfirst(x->x==row."From Bus", bus_data."Bus ID")
        dst = findfirst(x->x==row."To Bus", bus_data."Bus ID")
        propertys = Dict(:R => row."R", :X => row."X",
                         :rating => row."STE Rating"/baseP,
                         :id => row."UID")
        add_edge!(g, src, dst, propertys)
    end

    return g
end
