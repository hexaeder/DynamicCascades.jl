module DynamicCascades

using CSV
using DataFrames
using LightGraphs
using MetaGraphs
using Colors
import ColorSchemes
using GraphPlot
using Random: MersenneTwister
using Printf

export DATA_DIR
const DATA_DIR = abspath(@__DIR__, "..", "Data")

export IEEE_Network, load_rts96, import_system, describe_nodes, describe_edges, BusType

struct IEEE_Network
    graph::SimpleGraph{Int}
    num_of_buses::Int
    num_of_lines::Int
    num_of_generators::Int
    busdict::Dict{Int,Int}
    bustype::Vector{Int} # 1::Load, 2::Generator, 3::Slack
    voltage::Vector{Float64}
    active_power::Vector{Float64}
    inertia::Vector{Float64}
    conductance::Vector{Float64}
    susceptance::Vector{Float64}
    active_power_coupling::Vector{Float64}
    emergency_rating::Vector{Float64}
end

@enum BusType Load=1 Generator=2 Slack=3

function MetaGraphs.set_prop!(g, vs::Union{AbstractUnitRange, AbstractArray},
                              prop::Symbol, vals::Union{AbstractUnitRange, AbstractArray})
    length(vs) == length(vals) || throw(ArgumentError("vs and vals needs to be of same length!"))
    for (v, val) in zip(vs, vals)
        set_prop!(g, v, prop, val)
    end
end

"""
    import_system(sym::Symbol; kwargs...)::MetaGraph

Main entry point to load the systems. New systems should overload this function. Known implementations
- `import_system(:rtsgmlc)`: loads the GMLC update for the rts96
- `import_system(:rts96raw)`: loads the RTS96 system based on the raw files
- `import_system(:rts96prepared)`: loads the RTS96 system based on Antons prepared files

Those functions return a `MetaGraph` which has properties attached to the Nodes/Edges.

Graph properties:
- `BaseP` : Base power for PU
optional:
- `NodeP` : tuple of node property names which show appear first
- `EdgeP` : tuple of edge property names which show appear first

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

function describe_nodes(g::MetaGraph; kw=(id=Int[], x=Int[], y=Int[], type=String[]))
    df = DataFrame(;n=Int[], kw...)
    for n in 1:nv(g)
        row = push!(props(g, n), :n=>n)
        push!(df, row, cols=:union)
    end
    df
end

function describe_edges(g::MetaGraph)
    df = DataFrame(src=Int[], dst=Int[])
    for e in edges(g)
        row = push!(props(g, e), :src=>e.src, :dst=>e.dst)
        push!(df, row, cols=:union)
    end
    df
end

function import_system(::Val{:rtsgmlc})
    @info "Import system RTS-GMLC"
    data = joinpath(DATA_DIR, "RTS-GMLC")
    bus_data    = CSV.read(joinpath(data,"bus.csv"), DataFrame)
    gen_data    = CSV.read(joinpath(data,"gen.csv"), DataFrame)
    branch_data = CSV.read(joinpath(data,"branch.csv"), DataFrame)

    sort!(bus_data, "Bus ID")
    N = nrow(bus_data)
    g = MetaGraph(N)

    set_prop!(g, 1:N, :id, bus_data."Bus ID")
    set_prop!(g, 1:N, :V_mag, bus_data."V Mag")
    set_prop!(g, 1:N, :V_angle, bus_data."V Angle")
    set_prop!(g, 1:N, :x, bus_data."lng")
    set_prop!(g, 1:N, :y, bus_data."lat")
    set_prop!(g, 1:N, :type, bus_data."Bus Type")

    set_prop!(g, 1:N, :MW_load, bus_data."MW Load")
    set_prop!(g, 1:N, :MVAR_load, bus_data."MVAR Load")

    for (n, busid) in enumerate(bus_data."Bus ID")
        generators = gen_data[gen_data."Bus ID".== busid, ["Category", "MW Inj", "MVAR Inj", "Inertia MJ/MW"]]
        mw_inj   = sum(generators."MW Inj")
        mvar_inj = sum(generators."MVAR Inj")
        inertia  = sum(generators."Inertia MJ/MW")
        if mw_inj!=0.0 || mvar_inj!=0.0 || inertia != 0.0
            set_prop!(g, n, :MW_inj, mw_inj)
            set_prop!(g, n, :MVAR_inj, mvar_inj)
            set_prop!(g, n, :inertia, inertia)
        end
        mw_load = get_prop(g, n, :MW_load)
        mvar_load = get_prop(g, n, :MVAR_load)
        set_prop!(g, n, :MW, mw_inj - mw_load)
        set_prop!(g, n, :MVAR, mvar_inj - mvar_load)
    end

    for row in eachrow(branch_data)
        src = findfirst(x->x==row."From Bus", bus_data."Bus ID")
        dst = findfirst(x->x==row."To Bus", bus_data."Bus ID")
        # TODO B does not matches the calculated value at all!
        propertys = Dict(:R => row."R", :X => row."X",
                         :rating => row."STE Rating",
                         :id => row."UID")
        add_edge!(g, src, dst, propertys)
    end

    return g
end

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

    set_prop!(g, 1:N, :id, BusData.id)
    set_prop!(g, 1:N, :type, BusData.type)

    for (n, busid) in enumerate(BusData.id)
        generators = GeneratorData[GeneratorData.:id .== busid, :]
        generators = leftjoin(generators, SystemDynamicsData, on = :unit => :group)
        generators = generators[:, [:MW_inj, :MVar_inj, :inertia]]

        mw_inj   = sum(generators.:MW_inj)
        mvar_inj = sum(generators.:MVar_inj)
        inertia  = sum(generators.:inertia)

        if mw_inj!=0.0 || mvar_inj!=0.0 || inertia != 0.0
            # set_prop!(g, n, :MW_inj, mw_inj)
            # set_prop!(g, n, :MVAR_inj, mvar_inj)
            set_prop!(g, n, :inertia, inertia)
        end

        busload = BusLoadData[busid .∈ zip(BusLoadData.bid_1, BusLoadData.bid_2, BusLoadData.bid_3), :]
        if nrow(busload) == 0
            mw_load = 0.0
            mvar_load = 0.0
        elseif nrow(busload) == 1
            mw_load = busload.MW[1]
            mvar_load = busload.MVar[1]
        end

        P = (mw_inj - mw_load) / 100
        set_prop!(g, n, :P, P)
        Q = (mvar_inj - mvar_load) / 100
        set_prop!(g, n, :Q, Q)
    end

    # set_prop!(g, 1:N, :V_mag, voltage)

    for row in eachrow(BranchData)
        src = findfirst(x->x==row.from, BusData.id)
        dst = findfirst(x->x==row.to, BusData.id)
        # TODO B does not matches the calculated value at all!
        propertys = Dict(:R => row.R, :X => row.X,
                         :B => row.B, :rating => row.STE/100,
                         :id => row.id)
        add_edge!(g, src, dst, propertys)
    end
    return g
end

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
    for line in eachrow(LineData)
        add_edge!(g, line.source, line.dest)
    end

    # busdict connects bus ID form rts96 with internal bus number
    busdict = Dict{Int64,Int64}()
    for bus in eachrow(BusData)
        busdict[bus.ID] = bus.Number
    end

    # set busid
    set_prop!(g, 1:N, :id, BusData.ID)

    # bus types & dictionary for bus IDs
    bustype = convert(Array{Int64}, BusData.Type)
    set_prop!(g, 1:N, :type, bustype)

    # voltage magnitude
    voltage = FlowData.Vm
    set_prop!(g, 1:N, :V_mag, voltage)

    # active power
    P = - LoadData.P_Load
    P_Gen = losses ? GeneratorData.P_Gen : GeneratorData.P_Corr
    for i in 1:M
        idx = busdict[GeneratorData[i,:ID]]
        P[idx] += P_Gen[i]
    end
    set_prop!(g, 1:N, :P, P)

    # inertia
    I_Gen = GeneratorData.Inertia
    inertia = zeros(N)
    for i in 1:M
        idx = busdict[GeneratorData[i,:ID]]
        inertia[idx] += I_Gen[i]
    end
    set_prop!(g, 1:N, :inertia, inertia)

    # resistance & reactance
    X = LineData.x
    R = losses ? LineData.r : zeros(L)
    # Y = G + im B, Y=1/Z
    # TODO added -1 missing in B, is it correct now? seems so!
    G = R./(R.^2+X.^2)
    B = -X./(R.^2+X.^2) # calculating the absolute value here

    # caclulate active_power_coupling for each line
    K = [voltage[e.source] * voltage[e.dest] for e in eachrow(LineData)] .* -1 .* B

    # emergency rating
    rating = LineData.u

    for (e, r, x, b, rat) in zip(edges(g), R, X, B, rating)
        set_props!(g, e, Dict(:R=>r, :X=>x, :B=>b, :rating=>rat))
    end

    # check that the orders if edges in graph equals the order of edges in csv
    lineorder1 = [(e.source, e.dest) for e in eachrow(LineData)]
    lineorder2 = [(e.src, e.dst) for e in edges(g)]
    @assert lineorder1 == lineorder2 "Different order of edges in graph and csv!"

    return g
end

function load_rts96(;losses=false)
    # read data from csv files
    BusData = CSV.read(joinpath(DATA_DIR,"Bus.csv"), DataFrame)
    LineData = CSV.read(joinpath(DATA_DIR,"Line.csv"), DataFrame)
    GeneratorData = CSV.read(joinpath(DATA_DIR,"Generator.csv"), DataFrame)
    LoadData = CSV.read(joinpath(DATA_DIR,"Load.csv"), DataFrame)
    FlowData = CSV.read(joinpath(DATA_DIR,"Flow.csv"), DataFrame)

    # number of components
    N = nrow(BusData)
    L = nrow(LineData)
    M = nrow(GeneratorData)

    # graph structure
    g = SimpleGraph(N)
    for line in eachrow(LineData)
        add_edge!(g, line.source, line.dest)
    end

    # busdict connects bus ID form rts96 with internal bus number
    busdict = Dict{Int64,Int64}()
    for bus in eachrow(BusData)
        busdict[bus.ID] = bus.Number
    end

    # bus types & dictionary for bus IDs
    bustype = convert(Array{Int64}, BusData.Type)

    # voltage magnitude
    voltage = FlowData.Vm

    # active power
    P = - LoadData.P_Load
    GeneratorData
    P_Gen = losses ? GeneratorData.P_Gen : GeneratorData.P_Corr
    for i in 1:M
        idx = busdict[GeneratorData[i,:ID]]
        P[idx] += P_Gen[i]
    end

    # inertia
    I_Gen = GeneratorData.Inertia
    inertia = zeros(N)
    for i in 1:M
        idx = busdict[GeneratorData[i,:ID]]
        inertia[idx] += I_Gen[i]
    end

    # resistance & reactance
    X = LineData.x
    R = losses ? LineData.r : zeros(L)
    # Y = G + im B, Y=1/Z
    # TODO added -1 missing in B, is it correct now?
    G = R./(R.^2+X.^2)
    B = -X./(R.^2+X.^2) # calculating the absolute value here

    # caclulate active_power_coupling for each line
    K = [voltage[e.source] * voltage[e.dest] for e in eachrow(LineData)] .* -1 .* B

    # emergency rating
    rating = LineData.u

    # check that the orders if edges in graph equals the order of edges in csv
    lineorder1 = [(e.source, e.dest) for e in eachrow(LineData)]
    lineorder2 = [(e.src, e.dst) for e in edges(g)]
    @assert lineorder1 == lineorder2 "Different order of edges in graph and csv!"

    # constructor
    IEEE_Network(g, N, L, M, busdict, bustype, voltage, P, inertia, G, B, K, rating)
end

include("ND_model.jl")

end
