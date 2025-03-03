using Pkg
Pkg.activate("/home/brandner/DynamicCascades.jl")
# Pkg.instantiate()

using LinearAlgebra
print("Number of threads before setting"); print(LinearAlgebra.BLAS.get_num_threads()); print("\n")
BLAS.set_num_threads(1)
print("Number of threads after setting"); print(LinearAlgebra.BLAS.get_num_threads()); print("\n")


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

network = import_system(:rtsgmlc; damping = 0.1u"s", scale_inertia = 1, tconst = 0.01u"s")
x_static = @time steadystate(network)
print("x_static \n"); print(x_static); print("\n")

N = 100
network = import_system(:wattsstrogatz; N=N, β=0.7, graph_seed=124, distr_seed=1230, K=10, γ=1u"s", τ=1u"s", σ=1.0)
x_static = @time steadystate(network)
print("x_static \n"); print(x_static); print("\n")
