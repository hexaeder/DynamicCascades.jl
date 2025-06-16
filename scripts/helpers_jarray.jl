"""
utilities for running the same experiments on different machines.
"""

@assert VERSION == v"1.11.0"
const ON_YOGA = occursin("L7440", gethostname())
const ON_PIK_HPC = occursin("cs", gethostname())
const ON_POOL = occursin("pool", gethostname())

@info "Initialize environment"
if ON_YOGA
    PKG_DIR = "/home/brandner/.julia/dev/DynamicCascades"
    server_string = "L7440_"
elseif ON_PIK_HPC
    PKG_DIR = "/home/brandner/DynamicCascades.jl"
    server_string = "PIK_HPC_"
elseif ON_POOL
    PKG_DIR = "/users/stud/brandner/MA/repos/DynamicCascades.jl"
    server_string = "POOL_"
end


using Pkg
@assert typeof(PKG_DIR) == String "Please check `PKG_DIR`"
Pkg.activate(PKG_DIR) 

#= assuring that the nunber of BLAS threads is the same on different machines
(different number of threads may lead to different results) =#
using LinearAlgebra
print("Number of threads before setting"); print(LinearAlgebra.BLAS.get_num_threads()); print("\n")
BLAS.set_num_threads(1)
print("Number of threads after setting"); print(LinearAlgebra.BLAS.get_num_threads()); print("\n")

using DynamicCascades
using NetworkDynamics
using Graphs
using MetaGraphs
using Unitful
using Statistics
using Dates
using DataFrames
using CSV
using Serialization