using DataFrames
using Dates

# PARAMETERS ###################################################################
k = [4, 10]
# beta = [0.1, 0.5, 0.9]
beta = [0.1, 0.5]

#= frequency bounds [narrow bounds, wide bounds] bounds.
The value in numerator of round(0.1/(2*π) is the angular frequency =#
# freq_bounds = [round(0.1/(2*π), digits=2), round(0.5/(2*π), digits=2)]
freq_bounds = [round(0.1/(2*π), digits=2)]

# failure_modes = [trip_lines, trip_nodes]
# failure_modes = [[:dynamic, :dynamic], [:dynamic, :none], [:none, :dynamic]]
failure_modes = [[:dynamic, :dynamic]]

# inertia_values = [0.2, 0.5, 0.8, 1.1, 1.4, 1.7, 2.0, 3.0, 4.0, 5.1, 6.1, 7.1, 8.0, 9.0, 10.0, 15.0, 21.0]
# inertia_values = [0.2, 1.0, 2.1, 3.0, 4.0, 5.0, 7.2, 11.0, 15.0, 21.0]
# inertia_values = [0.2, 1.0, 5.0, 10.0, 15.0, 20.0]
inertia_values = [0.2, 5.0, 7.3]

# constant parameters
N_ensemble_size = 2 # 100
N_nodes = 100

K = 2 # coupling K
gamma = 1 # damping swing equation nodes γ
tau = 1 # time constant τ

alpha = 0.7 # tuning parameter α, :rating = α*K
init_pert = [:line] # initial perturbation set constant to an initial line failure

sigma = 1 # standard deviation σ
mu = 0 # mean μ
################################################################################
#= create hyperparameter: the order is chosen such, that with an increasing number
of finished jobs of the job array the size of the ensemble increases equally,
e.g. when half of all jobs are finished, for each ensemble half of the grids
shall be simulated=#
hyperparam = collect(Iterators.product(inertia_values, beta, freq_bounds, failure_modes, k,
    N_nodes, K, gamma, tau, alpha, init_pert, sigma, mu))[:]

# Repeat hyperparam N_ensemble_size times
hyperparam_ensemble = repeat(hyperparam, N_ensemble_size)
# println(hyperparameter[1])

# TODO check
# TODO maybe separate the two Symbols: [:dynamic, :dynamic]
# create dataframe (hpe stands for hyperparam_ensemble)
#= NOTE `:inertia_values` needs to be in the second column, otherwise the big while
loop below is corrupted =#
df_hpe = DataFrame(map(idx -> getindex.(hyperparam_ensemble, idx), eachindex(first(hyperparam_ensemble))),
    [:inertia_values, :beta, :freq_bounds, :failure_modes, :k, :N_nodes, :K, :gamma, :tau, :alpha, :init_pert, :sigma, :mu])

# add "ArrayTaskID" as first column of df
df_hpe = hcat(DataFrame(ArrayTaskID = 1:length(hyperparam_ensemble)), df_hpe)

# add columns `graph_seed`, `distr_seed` and `filepath` filled with `nan` to df
df_hpe[!, :graph_seed] .= 1; df_hpe[!, :distr_seed] .= 1; df_hpe[!, :filepath] .= "filepath"


function get_network_args(df_hpe::DataFrame, task_id::Int)
    N=df_hpe[task_id,:N_nodes]
    k=df_hpe[task_id,:k]
    β=df_hpe[task_id,:beta]
    graph_seed=df_hpe[task_id,:graph_seed]
    μ=df_hpe[task_id,:mu]
    σ=df_hpe[task_id,:sigma]
    distr_seed=df_hpe[task_id,:distr_seed]
    K=df_hpe[task_id,:K]
    α=df_hpe[task_id,:alpha]
    M=df_hpe[task_id,:inertia_values]*1u"s^2"
    γ=df_hpe[task_id,:gamma]*1u"s"
    τ=df_hpe[task_id,:tau]*1u"s"
    freq_bound=df_hpe[task_id,:freq_bounds]
    failure_mode=df_hpe[task_id,:failure_modes]
    init_pert=df_hpe[task_id,:init_pert]

    return N,k,β,graph_seed,μ,σ,distr_seed,K,α,M,γ,τ,freq_bound,failure_mode,init_pert
end

#= NOTE Could be written as `import_system_wrapper(df_hpe::DataFrame, task_id::Int)`
using multiple dispatch, however not useful here in a script.=#
# TODO Change function name?
function import_system_wrapper(df_hpe::DataFrame, task_id::Int)
    N,k,β,graph_seed,μ,σ,distr_seed,K,α,M,γ,τ,_,_,_ = get_network_args(df_hpe, task_id)
    return import_system(:wattsstrogatz; N=N, k=k, β=β, graph_seed=graph_seed,
        μ=μ, σ=σ, distr_seed=distr_seed, K=K, α=α, M=M, γ=γ, τ=τ)
end

function string_network_args(df_hpe::DataFrame, task_id::Int)
    N,k,β,graph_seed,μ,σ,distr_seed,K,α,M,γ,τ,freq_bound,failure_mode,init_pert = get_network_args(df_hpe, task_id)
    # return "N=$N, k=$k, β=$β, graph_seed=$graph_seed, μ=$μ, σ=$σ, distr_seed=$distr_seed, K=$K, α=$α, M=$M, γ=$γ, τ=$τ"
    return "failure_mode=$failure_mode,freq_bound=$freq_bound,N=$N,k=$k,β=$β,graph_seed=$graph_seed,μ=$μ,σ=$σ,distr_seed=$distr_seed,K=$K,α=$α,M=$M,γ=$γ,τ=$τ,init_pert=$init_pert"
end

string_network_args(df_hpe, task_id)

""" Generation of networks """
task_id = 1
# Loop over each ArrayTaskID:
while task_id <= length(df_hpe.ArrayTaskID)
    #= For every new parameter configuration (excluding the variation of the inertia
    parameter), generate a new network:=#
    if ((task_id-1) % length(inertia_values)) == 0
        string = string_network_args(df_hpe, task_id)
        println("Generate new network for ArrayTaskID=$task_id for the parameters $string"
        network = import_system_wrapper(df_hpe, task_id)
    end
    #= Check for every inertia value of a single parameter configuration if steady
    state within tolerance exists. If for one inertia value no steady state exists,
    go back and generate new network.=#
    try
        steadystate(network) # test if steady state exists

        # save/overwrite filepath, graph_seed, distr_seed in df
        df_hpe[task_id,:graph_seed] = graph_seed; df_hpe[task_id,:distr_seed] = distr_seed
        N,k,β,graph_seed,μ,σ,distr_seed,K,α,M,γ,τ,freq_bound,failure_mode,init_pert = get_network_args(df_hpe, task_id)
        df_hpe[task_id,:filepath] = folder/"failure_mode=$failure_mode,freq_bound=$freq_bound,k=$k,β=$β" + seeds

        # save grid with parameter configuration if last steady state for last inertia value exists
        if ((task_id-1) % length(inertia_values)) == (length(inertia_values) - 1)

            println("Save network for jobs with ArrayTaskIDs $(task_id-length(inertia_values)+1) to $task_id.")


            # TODO create folder if not already existing

            t=now()
            datetime = Dates.format(t, "yyyymmdd_HHMMSS.s")

            directory_params = "failure_mode=$failure_mode,freq_bound=$freq_bound,k=$k,β=$β"

            # create folder
            folder = string(directory_params, "/",datetime,"WS_inertia_vs_line+node_failures_f_bound=$freq_bound")
            directory = string(RESULTS_DIR,folder)
            mkpath(directory)


            isdir(dir) || mkdir(dir)

            # TODO save network
            folder/"failure_mode=$failure_mode,freq_bound=$freq_bound,k=$k,β=$β" + seeds
            graph_seed += 1; distr_seed += 1
        end
        task_id += 1
    # in case `steadystate(network)` throws exeption
    # TODO test this
    catch
        println("No static solution found within tolerance.") # TODO add params for which no solution was found
        # generate new grid by changing seeds
        graph_seed += 1; distr_seed += 1
        # start again with first inertia value
        task_id -= ((task_id-1) % length(inertia_values)) # TODO dann wieder zurück zu Zeile 65
    end
end

##############
sdir(dir) || mkdir(dir)

for i in [0.2, 1.0, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 15.0, 20.0]
    println("i $i \n")
    network = import_system(:wattsstrogatz; N=100, k=4, β=0.7, M=(i * 1u"s^2"), graph_seed=124, distr_seed=1230, K=K, γ=(gamma * 1u"s"), τ=(tau * 1u"s"), μ=mu, σ=sigma)
    savegraph("/home/vollmich/Desktop/delete_soon/test_graphs/test_graph$i.lg", network)
    # x_static = steadystate(network)
end
################

""" loop over parameter configurations """
# loop over all parameter configurations of β,freq_bounds,failure_modes,k
    # loop over `N_ensemble_size` grids and generate grid
        # for this grid loop over all `inertia_values`
            # if steady state for all inertia value (maybe check by relaxation)
                # in newly added columns change `nan` to
                    # graph_seed
                    # distr_seed
                    # grid file path (don't save the same grids with different inertia_values)
            # else
                # generate new grid by changing seeds
                # graph_seed += 1
                # distr_seed += 1


# Design choice: looping over ArrayTaskID vs over parameter configurations
#  - To read the data out I will be looping over ArrayTastID as well.
#  - => I can use the same code.
#  - Looping over ArrayTastID seems to be more natural!


# Save to CSV
CSV.write(joinpath(@__DIR__, "config.csv"), df_hpe)
