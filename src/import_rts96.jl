"""
    import_system(:rts96raw)

Import the RTS96 system from raw data as a MetaGraph.
"""
function import_system(::Val{:rts96raw}; losses=false)
    @info "Import system RTS-96 (from raw data)"
    # read data from csv files
    data = joinpath(DATA_DIR, "RTS-96")

    BusData = CSV.read(joinpath(data, "Table-01_BusData.txt"), DataFrame,
                       datarow=7, delim=' ', ignorerepeated=true, footerskip=9,
                       header=[:id, :name, :type, :MW_load, :MVar_load, :GL, :BL, :sub_area, :base_kV, :zone])

    BusLoadData = CSV.read(joinpath(data, "Table-05_BusLoadData.txt"), DataFrame,
                       datarow=6, delim=' ', ignorerepeated=true, footerskip=2,
                       header=[:bid_1, :bid_2, :bid_3, :sys_load_percentage, :MW, :MVar, :MW_peak, :Mvar_peak])

    GeneratorData = CSV.read(joinpath(data, "Table-07_GeneratorData.txt"), DataFrame,
                       datarow=6, delim=' ', ignorerepeated=true, footerskip=4,
                       header=[:id, :unit, :unit_id, :MW_inj, :MVar_inj, :MVar_max, :MVar_min, :V_pu])

    BranchData = CSV.read(joinpath(data, "Table-12_BranchData.txt"), DataFrame,
                       datarow=18, delim=' ', ignorerepeated=true, footerskip=2,
                       header=[:id, :from, :to, :L, :Lamp, :Dur, :Lamt, :R, :X, :B, :Con, :LTE, :STE, :Tr])

    SystemDynamicsData = CSV.read(joinpath(data, "Table-15_SystemDynamicsData.txt"), DataFrame,
                       datarow=7, delim=' ', ignorerepeated=true, footerskip=6,
                       header=[:group, :size, :type, :MVAbase, :reactance_unit, :reactance_transformer, :inertia, :damping])

    # number of components
    N = nrow(BusData)
    g = MetaGraph(N)

    baseP = 100
    set_prop!(g, :Pbase, baseP)
    set_prop!(g, :NodeProps, [:n, :id, :P, :type, :Q, :inertia, :Vm, :Vbase,
                              :P_load, :P_inj, :Q_load, :Q_inj])
    set_prop!(g, :EdgeProps, [:src, :dst, :R, :X, :rating, :id])

    set_prop!(g, 1:N, :id, BusData.id)
    set_prop!(g, 1:N, :type, bustype.(BusData.type))
    set_prop!(g, 1:N, :Vbase, BusData.base_kV)

    for (n, busid) in enumerate(BusData.id)
        # collect all generators which are attached to bus
        generators = GeneratorData[GeneratorData.:id .== busid, :]
        generators = leftjoin(generators, SystemDynamicsData, on = :unit => :group)
        generators = generators[:, [:MW_inj, :MVar_inj, :inertia, :V_pu]]

        # extract voltage setpoints in PU, make sure alle generators at bus have same setpoint
        vpus = unique(generators.V_pu)
        if length(vpus) == 1
            set_prop!(g, n, :Vm, vpus[1])
        elseif length(vpus) != 0
            error("There are multiple Vpu at single bus!")
        end

        # sum up inertia and P/Q of all attached generators
        P_inj   = sum(generators.MW_inj) / baseP
        Q_inj = sum(generators.MVar_inj) / baseP
        inertia  = sum(generators.inertia)
        if P_inj!=0.0 || Q_inj!=0.0 || inertia != 0.0 # set only if there are generators
            set_prop!(g, n, :P_inj, P_inj)
            set_prop!(g, n, :Q_inj, Q_inj)
            set_prop!(g, n, :inertia, inertia)
        end

        # find right bus in busloaddata
        busload = BusLoadData[busid .∈ zip(BusLoadData.bid_1, BusLoadData.bid_2, BusLoadData.bid_3), :]
        if nrow(busload) == 0
            P_load = 0.0
            Q_load = 0.0
        elseif nrow(busload) == 1
            P_load = busload.MW[1] / baseP
            Q_load = busload.MVar[1] / baseP
            set_prop!(g, n, :P_load, P_load)
            set_prop!(g, n, :Q_load, Q_load)
        else
            error("nrow(busload) = $(nrow(busload)) should be 0 or 1")
        end

        P = P_inj - P_load
        set_prop!(g, n, :P, P)
        Q = Q_inj - Q_load
        set_prop!(g, n, :Q, Q)
    end

    for row in eachrow(BranchData)
        src = findfirst(x->x==row.from, BusData.id)
        dst = findfirst(x->x==row.to, BusData.id)
        propertys = Dict(:R => row.R, :X => row.X,
                         :rating => row.STE/100,
                         :id => row.id)
        add_edge!(g, src, dst, propertys)
    end
    return g
end

"""
    import_system(:rts96prepared)

Import the RTS96 system from the prepared data as a MetaGraph.
"""
function import_system(::Val{:rts96prepared}; gen_γ, slack_γ, load_τ, losses=false)
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
    set_prop!(g, :NodeProps, [:n, :id, :P, :type, :Q, :inertia, :Vm, :Vbase,
                              :P_load, :P_inj, :Q_load, :Q_inj])
    set_prop!(g, :EdgeProps, [:src, :dst, :R, :X, :rating])

    # join all data on bus id
    joined = leftjoin(BusData, GeneratorData, on=:ID)
    joined = leftjoin(joined, LoadData, on=:ID)
    joined = leftjoin(joined, FlowData, on=:ID, makeunique=true)
    sort!(joined, :ID)

    set_prop!(g, 1:N, :id, joined.ID)
    set_prop!(g, 1:N, :Vbase, joined.Base_V)
    set_prop!(g, 1:N, :type, bustype.(joined.Type))
    set_prop!(g, 1:N, :Vm, joined.Vm)
    set_prop!(g, 1:N, :Va, joined.Va)
    set_prop!(g, 1:N, :inertia, joined.Inertia)
    set_prop!(g, 1:N, :P_load, joined.P_Load)
    set_prop!(g, 1:N, :Q_load, joined.Q_Load)

    # TODO what ist P_corr?
    P_inj = losses ? joined.P_Gen : joined.P_Corr
    @warn "Set Q_inj to zero because there is no data!"
    Q_inj = zeros(length(P_inj))
    set_prop!(g, 1:N, :P_inj, P_inj)
    set_prop!(g, 1:N, :Q_inj, Q_inj)

    # calculate total P and Q
    P = Missings.replace(P_inj, 0.0) .- joined.P_Load
    Q = Missings.replace(Q_inj, 0.0) .- joined.Q_Load
    set_prop!(g, 1:N, :P, P)
    set_prop!(g, 1:N, :Q, Q)

    for v in 1:nv(g)
        type = get_prop(g, v, :type)
        if type === :gen
            set_prop!(g, v, :damping, gen_γ)
        elseif type === :slack
            set_prop!(g, v, :damping, slack_γ)
        elseif type === :load
            set_prop!(g, v, :timescale, load_τ)
        end
    end

    # add the edge data
    for line in eachrow(LineData)
        add_edge!(g, line.source, line.dest)
    end

    # resistance, reactance & emergency rating
    R = losses ? LineData.r : zeros(length(LineData.r))
    set_prop!(g, edges(g), :R, R)
    set_prop!(g, edges(g), :X, LineData.x)
    set_prop!(g, edges(g), :rating, LineData.u)

    # check that the orders if edges in graph equals the order of edges in csv
    lineorder1 = [(e.source, e.dest) for e in eachrow(LineData)]
    lineorder2 = [(e.src, e.dst) for e in edges(g)]
    @assert lineorder1 == lineorder2 "Different order of edges in graph and csv!"

    return g
end
