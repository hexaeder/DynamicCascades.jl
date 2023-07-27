# using Pkg
# Pkg.activate("/home/brandner/DynamicCascades.jl")
# # Pkg.instantiate()

# using LinearAlgebra
# print("Number of threads before setting"); print(LinearAlgebra.BLAS.get_num_threads()); print("\n")
# BLAS.set_num_threads(1)
# print("Number of threads after setting"); print(LinearAlgebra.BLAS.get_num_threads()); print("\n")

using DynamicCascades
using Graphs
using MetaGraphs
using Unitful
using Statistics
using GraphMakie
using Colors
using DynamicCascades: PLOT_DIR
using CairoMakie
using Dates
using DataFrames
using CSV
using DelimitedFiles

damping = 0.1u"s"
scale_inertia = 0.2

network = import_system(:rtsgmlc; damping, scale_inertia, tconst = 0.01u"s")



(nd, p) = nd_model(network)


x0 = zeros(length(nd.syms));
# x_static = solve(SteadyStateProblem(nd, x0, p), SSRootfind())
x_static = solve(NonlinearProblem(nd, x0, p), NLSolveJL())

tol=1e-7
dx = similar(x_static)
nd(dx, x_static, p, 0.0)
residuum = maximum(abs.(dx))
@assert residuum < tol "No steady state found $residuum"

print("x_static \n"); print(x_static); print("\n")
