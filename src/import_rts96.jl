using NetworkLayout: spring
"""
    import_system(:rts96raw)

Import the RTS96 system from raw data as a MetaGraph.
"""
function import_system(::Val{:rts96raw}; losses=false, kwargs...)
    @info "Import system RTS-96 (from raw data)"
    # read data from csv files
    data = joinpath(DATA_DIR, "RTS-96")

    BusData = CSV.read(joinpath(data, "Table-01_BusData.txt"), DataFrame,
                       skipto=7, delim=' ', ignorerepeated=true, footerskip=9,
                       header=[:id, :name, :type, :MW_load, :MVar_load, :GL, :BL, :sub_area, :base_kV, :zone])

    BusLoadData = CSV.read(joinpath(data, "Table-05_BusLoadData.txt"), DataFrame,
                       skipto=6, delim=' ', ignorerepeated=true, footerskip=2,
                       header=[:bid_1, :bid_2, :bid_3, :sys_load_percentage, :MW, :MVar, :MW_peak, :Mvar_peak])

    GeneratorData = CSV.read(joinpath(data, "Table-07_GeneratorData.txt"), DataFrame,
                       skipto=6, delim=' ', ignorerepeated=true, footerskip=4,
                       header=[:id, :unit, :unit_id, :MW_inj, :MVar_inj, :MVar_max, :MVar_min, :V_pu])

    BranchData = CSV.read(joinpath(data, "Table-12_BranchData.txt"), DataFrame,
                       skipto=18, delim=' ', ignorerepeated=true, footerskip=2,
                       header=[:id, :from, :to, :L, :Lamp, :Dur, :Lamt, :R, :X, :B, :Con, :LTE, :STE, :Tr])

    SystemDynamicsData = CSV.read(joinpath(data, "Table-15_SystemDynamicsData.txt"), DataFrame,
                       skipto=7, delim=' ', ignorerepeated=true, footerskip=6,
                       header=[:group, :size, :type, :MVAbase, :reactance_unit, :reactance_transformer, :inertia, :damping])

    # number of components
    N = nrow(BusData)
    g = MetaGraph(N)

    baseP = 100u"MW"
    set_prop!(g, :Pbase, baseP)
    set_prop!(g, :NodeProps, [:n, :id, :P, :type, :Q, :H, :Vm, :Vbase,
                              :P_load, :P_inj, :Q_load, :Q_inj])
    set_prop!(g, :EdgeProps, [:src, :dst, :R, :X, :rating, :id])

    set_prop!(g, 1:N, :id, BusData.id)
    set_prop!(g, 1:N, :type, bustype.(BusData.type))
    set_prop!(g, 1:N, :Vbase, BusData.base_kV*u"kV")

    for (n, busid) in enumerate(BusData.id)
        # collect all generators which are attached to bus
        generators = GeneratorData[GeneratorData.:id .== busid, :]
        generators = leftjoin(generators, SystemDynamicsData, on = :unit => :group)
        generators = generators[:, [:MW_inj, :MVar_inj, :inertia, :V_pu]]

        # extract voltage setpoints in PU, make sure all generators at bus have same setpoint
        vpus = unique(generators.V_pu)
        if length(vpus) == 1
            set_prop!(g, n, :Vm, vpus[1]*u"pu")
        elseif length(vpus) != 0
            error("There are multiple Vpu at single bus!")
        end

        # sum up inertia and P/Q of all attached generators
        P_inj = sum(generators.MW_inj) * u"MW"
        Q_inj = sum(generators.MVar_inj) * u"MW"
        inertia = sum(generators.inertia) * u"MJ/MW"
        if !iszero(P_inj) || !iszero(Q_inj) || !iszero(inertia)
            set_prop!(g, n, :P_inj, P_inj)
            set_prop!(g, n, :Q_inj, Q_inj)
            set_prop!(g, n, :H, inertia)
        end

        # find right bus in busloaddata
        busload = BusLoadData[busid .∈ zip(BusLoadData.bid_1, BusLoadData.bid_2, BusLoadData.bid_3), :]
        if nrow(busload) == 0
            P_load = 0.0u"MW"
            Q_load = 0.0u"MW"
        elseif nrow(busload) == 1
            P_load = busload.MW[1] * u"MW"
            Q_load = busload.MVar[1] * u"MW"
            set_prop!(g, n, :P_load, P_load)
            set_prop!(g, n, :Q_load, Q_load)
        else
            error("nrow(busload) = $(nrow(busload)) should be 0 or 1")
        end

        set_prop!(g, n, :P, (P_inj - P_load)/baseP * u"pu")
        set_prop!(g, n, :Q, (Q_inj - Q_load)/baseP * u"pu")
    end

    for row in eachrow(BranchData)
        src = findfirst(x->x==row.from, BusData.id)
        dst = findfirst(x->x==row.to, BusData.id)
        propertys = Dict(:R => row."R"u"pu",
                         :X => row."X"u"pu",
                         :rating => row.STE*u"MW"/baseP * u"pu",
                         :id => row.id)
        add_edge!(g, src, dst, propertys)
    end

    losses || lossless!(g)
    set_missing_props!(g; kwargs...)

    return g
end

"""
    import_system(:rts96prepared)

Import the RTS96 system from the prepared data as a MetaGraph.
"""
function import_system(::Val{:rts96prepared}; losses=false, kwargs...)
    @info "Import system RTS-96 (from prepared csv files)"
    # read data from csv files
    data = joinpath(DATA_DIR, "RTS-96")
    BusData = CSV.read(joinpath(data,"Bus.csv"), DataFrame)
    LineData = CSV.read(joinpath(data,"Line.csv"), DataFrame)
    GeneratorData = CSV.read(joinpath(data,"Generator.csv"), DataFrame)
    LoadData = CSV.read(joinpath(data,"Load.csv"), DataFrame)
    FlowData = CSV.read(joinpath(data,"Flow.csv"), DataFrame)

    # number of components
    N = nrow(BusData)
    L = nrow(LineData)
    M = nrow(GeneratorData)

    # graph structure
    g = MetaGraph(N)
    # set base parameters
    baseP = 100
    set_prop!(g, :Pbase, baseP)
    set_prop!(g, :NodeProps, [:n, :id, :P, :type, :Q, :H, :Vm, :Vbase,
                              :P_load, :P_inj, :Q_load, :Q_inj])
    set_prop!(g, :EdgeProps, [:src, :dst, :R, :X, :rating])

    # join all data on bus id
    joined = leftjoin(BusData, GeneratorData, on=:ID)
    joined = leftjoin(joined, LoadData, on=:ID)
    joined = leftjoin(joined, FlowData, on=:ID, makeunique=true)
    sort!(joined, :ID)

    set_prop!(g, 1:N, :id, joined.ID)
    set_prop!(g, 1:N, :Vbase, joined.Base_V*u"kV")
    set_prop!(g, 1:N, :type, bustype.(joined.Type))
    set_prop!(g, 1:N, :Vm, joined.Vm*u"pu")
    set_prop!(g, 1:N, :Va, joined.Va*u"°")
    set_prop!(g, 1:N, :H, joined.Inertia*u"MJ/MW")
    set_prop!(g, 1:N, :P_load, joined.P_Load*u"pu")
    set_prop!(g, 1:N, :Q_load, joined.Q_Load*u"pu")

    P_inj = joined.P_Gen * u"pu"
    Q_inj = zeros(length(P_inj))
    set_prop!(g, 1:N, :P_inj, P_inj)
    set_prop!(g, 1:N, :Q_inj, Q_inj)

    # calculate total P and Q
    P = (Missings.replace(P_inj, 0.0) .- joined.P_Load) * u"pu"
    Q = (Missings.replace(Q_inj, 0.0) .- joined.Q_Load) * u"pu"
    set_prop!(g, 1:N, :P, P)
    set_prop!(g, 1:N, :Q, Q)

    # add the edge data
    for line in eachrow(LineData)
        add_edge!(g, line.source, line.dest)
    end

    # resistance, reactance & emergency rating
    R = losses ? LineData.r : zeros(length(LineData.r))
    set_prop!(g, edges(g), :R, R * u"pu")
    set_prop!(g, edges(g), :X, LineData.x * u"pu")
    set_prop!(g, edges(g), :rating, LineData.u * u"pu")

    # check that the orders if edges in graph equals the order of edges in csv
    lineorder1 = [(e.source, e.dest) for e in eachrow(LineData)]
    lineorder2 = [(e.src, e.dst) for e in edges(g)]
    @assert lineorder1 == lineorder2 "Different order of edges in graph and csv!"

    # set location property based on rts gmlc data
    data = joinpath(DATA_DIR, "RTS-GMLC")
    bus_data    = CSV.read(joinpath(data,"bus.csv"), DataFrame)
    x = bus_data."lng"
    y = bus_data."lat"
    xn = 10*(x .- minimum(x))./(maximum(x)-minimum(x)).-5
    yn = 10*(y .- minimum(y))./(maximum(y)-minimum(y)).-5
    pos2 = spring(g, initialpos=Point2f0.(zip(xn, yn)))

    # manual tweak som positions
    # pos2[70] = (-3,3.7)
    # pos2[53] = (-10,1)
    # pos2[54] = (-9,1)
    # pos2[58] = (-9,0.5)
    # pos2[56] = (-8.5,0.7)
    # pos2[55] = (-8.6,0.8)
    # pos2[42] = (4,0.5)
    # pos2[36] = (2.5,5)
    # pos2[32] = (4.5,5)
    # pos2[31] = (4.6,4.9)
    # pos2[29] = (6,5.5)
    # pos2[28] = (5.0,4.2)
    # pos2[30] = (5.0,5)
    # pos2[26] = (5.5,5)
    # pos3 = spring(g, initialpos=pos2, initialtemp=0.1)

    set_prop!(g, 1:N, :pos, pos2)

    losses || lossless!(g)
    set_missing_props!(g; kwargs...)

    return g
end
