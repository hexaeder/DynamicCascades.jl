module DynamicCascades

using CSV
using DataFrames
using LightGraphs
using LightGraphs: SimpleGraphs.SimpleEdgeIter
using MetaGraphs
using Missings

export DATA_DIR
const DATA_DIR = abspath(@__DIR__, "..", "Data")

export import_system, describe_nodes, describe_edges, bustype

"""
    bustype(i::Int)
    bustype(s::String)
    bustype(s::Symbol)

Helper to convert the different representations of bus types into eachother.
Base representation is Symbol:

    :load  == bustype(1) == bustype("PQ") == bustype(:PQ)
    :gen   == bustype(2) == bustype("PV") == bustype(:PV)
    :slack == bustype(3) == bustype("Ref") == bustype(:Ref)
"""
bustype(i::Int) = bustype(Val(i))
bustype(s::Union{String, Symbol}) = bustype(Val(Symbol(s)))
bustype(::Union{Val{:PQ},  Val{1}, Val{:load}})  = :load
bustype(::Union{Val{:PV},  Val{2}, Val{:gen}})   = :gen
bustype(::Union{Val{:Ref}, Val{3}, Val{:slack}}) = :slack

"""
    import_system(sym::Symbol; kwargs...)::MetaGraph

Main entry point to load the systems. New systems should overload this function. Known implementations
- `import_system(:rtsgmlc)`: loads the GMLC update for the rts96
- `import_system(:rts96raw)`: loads the RTS96 system based on the raw files
- `import_system(:rts96prepared)`: loads the RTS96 system based on Antons prepared files

Those functions return a `MetaGraph` which has properties attached to the Nodes/Edges.

Graph properties:
- `Pbase` : Base power for PU
optional:
- `NodeProps` : tuple of node property names which show appear first
- `EdgeProps` : tuple of edge property names which show appear first

Node properties:
- `P` : active power in PU
- `Q` : reactive power in PU
- `Vbase` : base voltage in kV
- `Vm` : voltage magnitude in PU
optional:
- `inertia` : inertia in MJ/MW (see H in Wikipedia page)
- `x`, `y` : position of Bus for plotting purposes
- `Va` : voltage angle in rad

Edge properties:
- `rating` : short term emergency rating in PU
- `R` : resistance in PU
- `X` : reactance in PU
"""
import_system(sym::Symbol; kwargs...) = import_system(Val(sym); kwargs...)

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
        P_inj   = sum(generators."MW Inj") / baseP
        Q_inj = sum(generators."MVAR Inj") / baseP
        inertia  = sum(generators."Inertia MJ/MW")
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
function import_system(::Val{:rts96prepared}; losses=false)
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

    # add the edge data
    for line in eachrow(LineData)
        add_edge!(g, line.source, line.dest)
    end

    # resistance, reactance & emergency rating
    set_prop!(g, edges(g), :R, LineData.r)
    set_prop!(g, edges(g), :X, LineData.x)
    set_prop!(g, edges(g), :rating, LineData.u)

    # check that the orders if edges in graph equals the order of edges in csv
    lineorder1 = [(e.source, e.dest) for e in eachrow(LineData)]
    lineorder2 = [(e.src, e.dst) for e in edges(g)]
    @assert lineorder1 == lineorder2 "Different order of edges in graph and csv!"

    return g
end

"""
    MetaGraphs.set_prop!(g, vs::Iterable, prop::Symbol, vals::Iterable)

Set same property `prop` with different values `vals` for differet vertices `vs`.
"""
function MetaGraphs.set_prop!(g, vs::Union{AbstractUnitRange, AbstractArray},
                              prop::Symbol, vals::Union{AbstractUnitRange, AbstractArray})
    length(vs) == length(vals) || throw(ArgumentError("vs and vals needs to be of same length!"))
    for (v, val) in zip(vs, vals)
        if !ismissing(val)
            set_prop!(g, v, prop, val)
        end
    end
end

"""
    MetaGraphs.set_prop!(g, edges::SImpleEdgeIter, prop::Symbol, vals::Iterable)

Set same property `prop` with different values `vals` for different edges in EdgeIter `edges`.
"""
function MetaGraphs.set_prop!(g, edges::SimpleEdgeIter, prop::Symbol, vals)
    length(edges) == length(vals) || throw(ArgumentError("EdgeIter and vals needs to be of same length!"))
    for (e, val) in zip(edges, vals)
        if !ismissing(val)
            set_prop!(g, e, prop, val)
        end
    end
end

"""
    describe_nodes(g::MetaGraph; firstcols=Vector{String}())

Returns DataFrame with all the node meta data.
"""
function describe_nodes(g::MetaGraph; firstcols=Vector{String}())
    df = DataFrame(; n=Int[])
    for n in 1:nv(g)
        row = push!(props(g, n), :n=>n)
        push!(df, row, cols=:union)
    end
    firstcols = String.(firstcols) # convert to string if given as Symbol
    has_prop(g, :NodeProps) &&
        append!(firstcols, String.(get_prop(g, :NodeProps)))
    append!(firstcols, names(df)) |> unique!
    select!(df, firstcols)
end

"""
    describe_nodes(g::MetaGraph; firstcols=Vector{String}())

Returns DataFrame with all the edge meta data.
"""
function describe_edges(g::MetaGraph; firstcols=Vector{String}())
    df = DataFrame(; src=Int[], dst=Int[])
    for e in edges(g)
        row = push!(props(g, e), :src=>e.src, :dst=>e.dst)
        push!(df, row, cols=:union)
    end
    firstcols = String.(firstcols) # convert to string if given as Symbol
    has_prop(g, :EdgeProps) &&
        append!(firstcols, String.(get_prop(g, :EdgeProps)))
    append!(firstcols, names(df)) |> unique!
    select!(df, firstcols)
end


include("ND_model.jl")

end
