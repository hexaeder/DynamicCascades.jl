module DynamicCascades

using CSV
using DataFrames
using Graphs
using MetaGraphs
using Missings

export DATA_DIR, RAWRESULTS_DIR, RESULTS_DIR
const DATA_DIR = abspath(@__DIR__, "..", "data")
const RAWRESULTS_DIR = abspath("/Users/hw/MAScratch")
const RESULTS_DIR = abspath(@__DIR__, "..", "result")
const PLOT_DIR = abspath(@__DIR__, "..", "..", "thesis", "figures")

export import_system, describe_nodes, describe_edges, bustype, is_static_state

include("utils.jl")

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
bustype(s) = bustype(Val(Symbol(s)))
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
- `NodeProps` : tuple of node property names which should appear first in describe functions
- `EdgeProps` : tuple of edge property names which should appear first in describe functions

Node properties:
- `P` : active power in PU
- `Q` : reactive power in PU
- `Vbase` : base voltage in kV
- `Vm` : voltage magnitude in PU
optional:
- `inertia` : inertia in MJ/MW (see H in Wikipedia page)
- `damping` : damping factor γ for swing equations
- `timescale` : time scale factor `D` for dynamic loads
- `x`, `y` : position of Bus for plotting purposes
- `Va` : voltage angle in rad

Edge properties:
- `rating` : short term emergency rating in PU
- `R` : resistance in PU
- `X` : reactance in PU
"""
import_system(sym::Symbol; kwargs...) = import_system(Val(sym); kwargs...)

include("import_rtsgmlc.jl")
include("import_rts96.jl")
include("import_kaiser2020.jl")
include("import_schaefer2018.jl")

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
include("analyse_solution.jl")
include("inspect_solution.jl")

end
