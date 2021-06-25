using DynamicCascades
using LightGraphs
using MetaGraphs
using Statistics
using GLMakie
using GraphMakie

network = import_system(:schaefer2018)

x = steadystate(network)
# x0 = [-0.12692637482862684, -1.3649456633810975e-6, 0.14641121510104085, 4.2191082676726005e-7, -0.24376507587890778, 1.567589744768255e-6, -0.12692637482862684, -1.3649456633810975e-6, 0.35120661043511864, 7.403907552948938e-7]

# issteadystate(network, x0)
# issteadystate(network, x)

sol = simulate(network;
               initial_fail=[5],
               trip_lines=true,
               solverargs=(;dtmax=0.1));
inspect_solution(sol;)

get_prop(network, 1:5, :model)

nd_model(network)

describe_nodes(network)
describe_edges(network)
