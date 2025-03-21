"""
Debugging new ND model using MWE consisting of a 4 node WS network for the old ND model: Checking the initial steady state. 
"""

using Revise
include(abspath(@__DIR__, "..", "helpers_jarray.jl"))
using DynamicCascades
using MetaGraphs

###
### assign WS network and parameters
###
K = 3
α = 0.7
M = 1
γ = 1
τ = 1
β = 0.5
freq_bound = 0.3
network = import_system(:wattsstrogatz; N=4, k=2, β=β, graph_seed=1, distr_seed=1, K=K, α=α, M=M*u"s^2",  γ=γ*u"s", τ=τ*u"s")


using SteadyStateDiffEq
(nd, p) = nd_model(network)
x0 = zeros(length(nd.syms));
x_static = solve(SteadyStateProblem(nd, x0, p), SSRootfind())
# x_static = solve(NonlinearProblem(nd, x0, p), NLSolveJL())
θidx = idx_containing(nd, "θ")
zeroidx=1
if zeroidx !== nothing
    offset = x_static[θidx[zeroidx]]
    x_static[θidx] .= x_static[θidx] .- offset
    @assert iszero(x_static[θidx[zeroidx]])
end
x_static[1:2:end]