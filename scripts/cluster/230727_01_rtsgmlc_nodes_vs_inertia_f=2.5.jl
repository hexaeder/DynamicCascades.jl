using Pkg
Pkg.activate("/home/brandner/DynamicCascades.jl")
Pkg.instantiate()

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

freq_bound = round(2.5/(2*π), digits=2)
# create folder
t=now()
datetime = Dates.format(t, "yyyymmdd_HHMMSS.s") # https://riptutorial.com/julia-lang/example/20476/current-time
folder = string("/",datetime,"inertia_vs_node_failures_f_bound=$freq_bound")
directory = string(RESULTS_DIR,folder)
mkpath(directory)

nd, = nd_model(network)
ω_state_idxs = idx_containing(nd, "ω")
gen_node_idxs = map(s -> parse(Int, String(s)[4:end]), nd.syms[ω_state_idxs])


damping = 0.1u"s"
scale_inertia_values = [0.2, 0.5, 0.8, 1.1, 1.4, 1.7, 2.0, 3.0, 4.0, 5.1, 6.1, 7.1, 8.0, 9.0, 10.0, 15.0, 21.0] # varying parameter
df_all_failures = DataFrame()
@time for scale_inertia in scale_inertia_values
    network = import_system(:rtsgmlc; damping, scale_inertia, tconst = 0.01u"s")
    x_static = steadystate(network)
    println("Scaling of inertia $scale_inertia \n steady state \n $x_static \n ")
    number_failures = Float64[]
    for i in gen_node_idxs
        sol = simulate(network;
                       x_static = x_static,
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
        push!(number_failures, length(sol.failures_nodes.saveval))
    end
    df_all_failures[!, string(scale_inertia)] = number_failures
end
# write failures for each node/line to .csv
CSV.write(string(directory,"/all_failures.csv"), df_all_failures)

# calculate relative number of failures
df_all_failures = DataFrame(CSV.File(string(directory,"/all_failures.csv")))
rel_number_failures = Float64[]
for scale_inertia in scale_inertia_values
    number_failures = df_all_failures[!, string(scale_inertia)]
    push!(rel_number_failures, mean(number_failures)/length(gen_node_idxs))
end


df_inertia_vs_failures = DataFrame("scale_inertia_values" => scale_inertia_values, "rel_failures" => rel_number_failures)
CSV.write(string(directory,"/inertia_vs_failures.csv"), df_inertia_vs_failures)

# load data
df_inertia_vs_failures = DataFrame(CSV.File(string(directory,"/inertia_vs_failures.csv")))

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
y = df_inertia_vs_failures.rel_failures

# scatter!(x, y, color = :blue,label = "Test")
scatter!(x, y, color = :blue)
# axislegend()

CairoMakie.save(string(directory,"/inertia_vs_number_node_failures.pdf"),fig)
