module DynamicCascades

using CSV
using DataFrames
using Graphs
using MetaGraphs
using Missings
using Unitful

export DATA_DIR, MA_DATA, RESULTS_DIR, PLOT_DIR, MA_DIR, RES_GEN_NET # , F_BELEG_DIR
export ON_YOGA, ON_PIK_HPC, ON_POOL

# @assert VERSION == v"1.8.4"
const ON_YOGA = occursin("L7440", gethostname())
const ON_PIK_HPC = occursin("cs", gethostname())
const ON_POOL = occursin("pool", gethostname())

if ON_YOGA
    const MA_DATA = "/home/brandner/nb_data/HU_Master/2122WS/MA/MA_data/"
    const RESULTS_DIR = "/home/brandner/nb_data/HU_Master/2122WS/MA/MA_data/results_NB/" # generated data
    const PLOT_DIR = "/home/brandner/nb_data/HU_Master/2122WS/MA/MA_data/figures/"
    const RES_GEN_NET = "/home/brandner/nb_data/HU_Master/2122WS/MA/MA_data/res_gen_net" # results general networks
elseif ON_PIK_HPC
    const MA_DATA = abspath(@__DIR__, "..", "..", "MA_data")
    const RESULTS_DIR = abspath(@__DIR__, "..", "..", "MA_data", "results_NB") # generated data
    const PLOT_DIR = abspath(@__DIR__, "..", "..", "MA_data", "figures")
    const RES_GEN_NET = abspath(@__DIR__, "..", "..", "MA_data", "res_gen_net") # results general networks
elseif ON_POOL
    # TODO
end

const DATA_DIR = abspath(@__DIR__, "..", "data") # data used for simulations
const MA_DIR = "/home/brandner/nb_data/repos/NLCp/paper_figs/"
# const F_BELEG_DIR = abspath(@__DIR__, "..", "..", "..", "..", "repos", "/Private_MA/F_BELEG/F_BELEG_figs")

export import_system, describe_nodes, describe_edges, bustype, is_static_state

# define perunit as unit
@unit pu "p.u." PerUnit 1 false
Unitful.register(DynamicCascades)
function __init__()
    Unitful.register(DynamicCascades)
end

include("helpers.jl")

"""
    bustype(i::Int)
    bustype(s::String)
    bustype(s::Symbol)

Helper to convert the different representations of bus types into eachother.
syncon stands for synchronous condenser.
Base representation is Symbol:

    :load   == bustype(1) == bustype("PQ") == bustype(:PQ)
    :gen    == bustype(2) == bustype("PV") == bustype(:PV)
    :syncon == bustype(3) == bustype("Ref") == bustype(:Ref)
"""
bustype(i::Int) = bustype(Val(i))
bustype(s) = bustype(Val(Symbol(s)))
bustype(::Union{Val{:PQ},  Val{1}, Val{:load}})   = :load
bustype(::Union{Val{:PV},  Val{2}, Val{:gen}})    = :gen
bustype(::Union{Val{:Ref}, Val{3}, Val{:syncon}}) = :syncon

"""
    import_system(sym::Symbol; kwargs...)::MetaGraph

Main entry point to load the systems. New systems should overload (add a method to)
this function. Known implementations
- `import_system(:rtsgmlc)`: loads the GMLC update for the rts96
- `import_system(:rts96raw)`: loads the RTS96 system based on the raw files
- `import_system(:rts96prepared)`: loads the RTS96 system based on Antons prepared files

Those functions return a `MetaGraph` which has properties attached to the Nodes/Edges.

Graph properties:
- `Pbase` : Base power in PU
optional:
- `NodeProps` : tuple of node property names which should appear first in describe functions
- `EdgeProps` : tuple of edge property names which should appear first in describe functions

Node properties:
- `n` : node index
- `P` : active power in PU
- `Q` : reactive power in PU
- `Vbase` : base voltage in kV
- `Vm` : voltage magnitude in PU
optional:
- `inertia` : inertia in MJ/MW (see H in Wikipedia page)
- `damping` : damping factor γ for swing equations
- `timeconstant` : time cosntant `τ` for dynamic loads
- `x`, `y` : position of Bus for plotting purposes
- `Va` : voltage angle in rad

Edge properties:
- `rating` : short term emergency rating (= maximum capacity) in PU
- `R` : resistance in PU
- `X` : reactance in PU
"""
import_system(sym::Symbol; kwargs...) = import_system(Val(sym); kwargs...)

include("import_rtsgmlc.jl")
include("import_general_networks.jl")

include("ND_model.jl")
include("visualization.jl")

end
