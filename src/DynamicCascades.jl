module DynamicCascades

using CSV
using DataFrames
using LightGraphs
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

include("import_rtsgmlc.jl")
include("import_rts96.jl")

import MetaGraphs: set_prop!, get_prop

"""
    set_prop!(g, vertices::Iterable, prop::Symbol, vals::Iterable)

Set same property `prop` with different values `vals` for differet vertices `vertices`.
"""
function set_prop!(g, vertices, prop::Symbol, vals)
    length(vertices) == length(vals) || throw(ArgumentError("vertices and vals needs to be of same length!"))
    for (v, val) in zip(vertices, vals)
        if !ismissing(val)
            set_prop!(g, v, prop, val)
        end
    end
end

"""
    MetaGraphs.set_prop!(g, edges::AbstractEdgeIter, prop::Symbol, vals::Iterable)

Set same property `prop` with different values `vals` for different edges in EdgeIter `edges`.
"""
function set_prop!(g, edges::AbstractEdgeIter, prop::Symbol, vals)
    length(edges) == length(vals) || throw(ArgumentError("EdgeIter and vals needs to be of same length!"))
    for (e, val) in zip(edges, vals)
        if !ismissing(val)
            set_prop!(g, e, prop, val)
        end
    end
end

"""
    MetaGraphs.get_prop(g, vertices::Iterable, prop::Symbol)

Get same property `prop` with different values `vals` for differet vertices `vertices`.
"""
get_prop(g, vertices, prop::Symbol) = [get_prop(g, v, prop) for v in vertices]

"""
    MetaGraphs.get_prop(g, edges::AbstractEdgeIterable, prop::Symbol)

Get same property `prop` with different values `vals` for differet vertices `vertices`.
"""
get_prop(g, edges::AbstractEdgeIter, prop::Symbol) = [get_prop(g, e, prop) for e in edges]

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
