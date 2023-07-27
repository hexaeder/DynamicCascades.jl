using Pkg
Pkg.activate("/home/brandner/DynamicCascades.jl")
# Pkg.instantiate()

using LinearAlgebra
print("Number of threads before setting"); print(LinearAlgebra.BLAS.get_num_threads()); print("\n")
BLAS.set_num_threads(1)
print("Number of threads after setting"); print(LinearAlgebra.BLAS.get_num_threads()); print("\n")


using NetworkDynamics
using Graphs
using OrdinaryDiffEq
using SciMLNLSolve
using Random
using NetworkLayout

function swing_equation(dv, v, edges, p,t)
    P, I, γ = p
    dv[1] = v[2]
    dv[2] = P - γ * v[2] + flow_sum(edges)
    dv[2] = dv[2] / I
    nothing
end

function flow_sum(edges)
    sum = 0.0
    for e in edges
        sum -= e[1]
    end
    return sum
end

function simple_edge(e, v_s, v_d, K, t)
    e[1] = K * sin(v_d[1] - v_s[1])
end

N=50 # choose even integer
k=4; β=0.7
network = watts_strogatz(N, k, β; seed=124)

# Set parameter
I = 1.0
γ = 1.0
K = 10

node_p = Vector{NTuple{3,Float64}}(undef, nv(network)) # (power, inertia, coupling)
node_p = [(isodd(v) ? -1.0 : 1.0, I, γ) for v in 1:nv(network)]

edge_p = K

p = (node_p, edge_p)

# Define nodes/edges and network
odevertex  = ODEVertex(; f=swing_equation, dim=2, sym=[:δ, :ω])
staticedge = StaticEdge(; f=simple_edge, dim=1, sym=[:P], coupling=:antisymmetric)
nd = network_dynamics(odevertex, staticedge, network)

x0 = zeros(length(nd.syms));
# x_static = solve(SteadyStateProblem(nd, x0, p), SSRootfind())
x_static = @time solve(NonlinearProblem(nd, x0, p), NLSolveJL())
println("x_static $x_static \n")

# tol=1e-7
# dx = similar(x_static)
# nd(dx, x_static, p, 0.0)
# residuum = maximum(abs.(dx))
# @assert residuum < tol "No steady state found $residuum"
