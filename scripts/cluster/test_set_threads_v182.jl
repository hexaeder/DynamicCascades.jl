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


# damping = 0.1u"s"
# network = import_system(:rtsgmlc; damping, scale_inertia=1.0, tconst = 0.01u"s")
#

# create folder
t=now()
datetime = Dates.format(t, "yyyymmdd_HHMMSS.s") # https://riptutorial.com/julia-lang/example/20476/current-time
folder = string("/",datetime,"inertia_vs_line_failures")
directory = string(RESULTS_DIR,folder)
mkpath(directory)
damping = 0.1u"s"
# scale_inertia_values = [0.2, 0.5, 1, 1.5, 2, 7, 10, 20] # varying parameter
scale_inertia_values = [0.2] # varying parameter
df_all_failures = DataFrame()
@time for scale_inertia in scale_inertia_values
    network = import_system(:rtsgmlc; damping, scale_inertia, tconst = 0.01u"s")
    # ###############
    sol = simulate(network;
                   tspan = (0, 2000),
                   solverargs = (;dtmax=0.01),
                   terminate_steady_state=true,
                   verbose = true);
    x_static = sol.sol[end]
    print("Steady state \n"); print(x_static); print("\n")
    # ###############
    number_failures = Float64[]
    # for i in 1:ne(network)
    for i in 1:4
        sol = simulate(network;
                       x_static = x_static,
                       initial_fail = Int[i],
                       init_pert = :line,
                       tspan = (0, 50),
                       trip_lines = :dynamic,
                       trip_nodes = :none,
                       trip_load_nodes = :none,
                       f_min = -2.5,
                       f_max = 1.5,
                       solverargs = (;dtmax=0.01),
                       verbose = true);
        push!(number_failures, length(sol.failures.saveval)-1) # `-1` as we don't want to count the initial failure
    end
    df_all_failures[!, string(scale_inertia)] = number_failures
end
# write failures for each node/line to .csv
CSV.write(string(directory,"/all_failures.csv"), df_all_failures)

# calculate relative number of failures
df_all_failures = DataFrame(CSV.File(string(directory,"/all_failures.csv")))
rel_number_failures = Float64[]
network = import_system(:rtsgmlc; damping, tconst = 0.01u"s")
for scale_inertia in scale_inertia_values
    number_failures = df_all_failures[!, string(scale_inertia)]
    push!(rel_number_failures, mean(number_failures)/(ne(network)-1))
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

CairoMakie.save(string(directory,"/inertia_vs_number_line_failures.pdf"),fig)
