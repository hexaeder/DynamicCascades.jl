using DynamicCascades
using Graphs
using MetaGraphs
using Statistics
using GraphMakie
using Colors
using DynamicCascades: PLOT_DIR
using NetworkDynamics
using Unitful
using Unitful: @u_str
using DataFrames
using CSV

# using CairoMakie

using GLMakie
GLMakie.activate!()

network = import_system(:rtsgmlc; damping = 0.1u"s", tconst = 0.01u"s")
n = describe_nodes(network)

layout = _ -> n.pos

node_color = map(n.id) do id
    if 100 < id < 200
        :blue
    elseif 200 < id < 300
        :red
    elseif 300 < id < 400
        :green
    end
end

graphplot(network; layout, node_color)
