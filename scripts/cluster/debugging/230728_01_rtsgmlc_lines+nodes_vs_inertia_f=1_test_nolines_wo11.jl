const ON_CLUSTER = occursin("login", gethostname())

using Pkg
PKG_DIR = ON_CLUSTER ? "/home/brandner/DynamicCascades.jl" : joinpath(@__DIR__, "..", "..", "..")
Pkg.activate(PKG_DIR)
# Pkg.instantiate()

using LinearAlgebra
print("Number of threads before setting"); print(LinearAlgebra.BLAS.get_num_threads()); print("\n")
BLAS.set_num_threads(1)
print("Number of threads after setting"); print(LinearAlgebra.BLAS.get_num_threads()); print("\n")

using Revise
using DynamicCascades
using NetworkDynamics
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

freq_bound = round(1.0/(2*π), digits=2)
# create folder
t=now()
datetime = Dates.format(t, "yyyymmdd_HHMMSS.s") # https://riptutorial.com/julia-lang/example/20476/current-time
folder = string("/",datetime,"inertia_vs_line+node_failures_f_bound=$freq_bound")
directory = string(RESULTS_DIR,folder)
mkpath(directory)

damping = 0.1u"s"
network0 = import_system(:rtsgmlc; damping, scale_inertia = 0.2, tconst = 0.01u"s")
nd, = nd_model(network0)
ω_state_idxs = idx_containing(nd, "ω")
gen_node_idxs = map(s -> parse(Int, String(s)[4:end]), nd.syms[ω_state_idxs])

scale_inertia_values = [0.2, 0.5, 0.8] # varying parameter
df_all_failures = DataFrame()
df_all_failures_nodes = DataFrame()
@time for scale_inertia in scale_inertia_values
    network = import_system(:rtsgmlc; damping, scale_inertia, tconst = 0.01u"s")
    x_static = steadystate(network)
    println("Scaling of inertia $scale_inertia \n steady state \n $x_static \n ")
    number_failures = Float64[]
    number_failures_nodes = Float64[]
    for i in 1:ne(network)
        if i==3
            continue
        end
        sol = simulate(network;
                       initial_fail = Int[i],
                       init_pert = :line,
                       tspan = (0, 500),
                       trip_lines = :none,
                       trip_nodes = :dynamic,
                       trip_load_nodes = :none,
                       f_min = -freq_bound,
                       f_max = freq_bound,
                       solverargs = (;dtmax=0.01),
                       verbose = true);
        push!(number_failures, length(sol.failures.saveval)-1) # `-1` as we don't want to count the initial failure
        push!(number_failures_nodes, length(sol.failures_nodes.saveval))
    end
    df_all_failures[!, string(scale_inertia)] = number_failures
    df_all_failures_nodes[!, string(scale_inertia)] = number_failures_nodes
end


# write failures for each node/line to .csv
CSV.write(string(directory,"/all_failures.csv"), df_all_failures)
CSV.write(string(directory,"/all_failures_nodes.csv"), df_all_failures_nodes)

# calculate relative number of failures
#lines
df_all_failures = DataFrame(CSV.File(string(directory,"/all_failures.csv")))
rel_number_failures = Float64[]
network = import_system(:rtsgmlc; damping, tconst = 0.01u"s")
for scale_inertia in scale_inertia_values
    number_failures = df_all_failures[!, string(scale_inertia)]
    push!(rel_number_failures, mean(number_failures)/(ne(network)-1))
end
df_inertia_vs_failures = DataFrame("scale_inertia_values" => scale_inertia_values, "rel_failures" => rel_number_failures)
CSV.write(string(directory,"/inertia_vs_failures.csv"), df_inertia_vs_failures)

# nodes
df_all_failures_nodes = DataFrame(CSV.File(string(directory,"/all_failures_nodes.csv")))
rel_number_failures = Float64[]
for scale_inertia in scale_inertia_values
    number_failures_nodes = df_all_failures_nodes[!, string(scale_inertia)]
    push!(rel_number_failures, mean(number_failures_nodes)/length(gen_node_idxs))
end
df_inertia_vs_failures_nodes = DataFrame("scale_inertia_values" => scale_inertia_values, "rel_failures" => rel_number_failures)
CSV.write(string(directory,"/inertia_vs_failures_nodes.csv"), df_inertia_vs_failures_nodes)

# load data
df_inertia_vs_failures = DataFrame(CSV.File(string(directory,"/inertia_vs_failures.csv")))
df_inertia_vs_failures_nodes = DataFrame(CSV.File(string(directory,"/inertia_vs_failures_nodes.csv")))

# plot data
fig = Figure(fontsize = 30)
Axis(fig[1, 1],
    # title = L"Decreasing $G_{av}$ for one sample grid",
    # titlesize = 30,
    xlabel = "inertia M [s^2]",
    # xlabelsize = 30,
    ylabel = "normalized average of failures",
    # ylabelsize = 30
)

x = df_inertia_vs_failures.scale_inertia_values
y_lines = df_inertia_vs_failures.rel_failures
y_nodes = df_inertia_vs_failures_nodes.rel_failures


scatter!(x, y_lines, label = "lines")
scatter!(x, y_nodes, label = "nodes")
axislegend(position = :rb)

CairoMakie.save(string(directory,"/inertia_vs_number_line+node_failures.pdf"),fig)
