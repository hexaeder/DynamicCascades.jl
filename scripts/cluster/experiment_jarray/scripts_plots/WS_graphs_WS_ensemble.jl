"""
Plotting graphs and initial power flow of ensemble
"""
using Revise
using DynamicCascades
using Graphs
using MetaGraphs
using Unitful
using Statistics
using GraphMakie
using Colors
using DynamicCascades: PLOT_DIR

using CairoMakie

function plotnetwork(fig, sol, t)
    ax = Axis(fig)
    hidedecorations!(ax), hidespines!(ax)

    gpargs = gparguments(sol, t;
                         colortype=:abssteady,
                         Δω=Observable(2.5),
                         offlinewidth=3,
                         offlinecolor=colorant"lightgray",
                         ecolorscaling = Observable(1.0),
                         node_size=15,
                         show_labels=false)
    p = graphplot!(ax, network; gpargs...)
    return ax, p
end

exp_name_date = "WS_k=4_exp03_I_over_Dsq_test_PIK_HPC_K_=3,N_G=32_20250121_084915.253"
exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
df_config = DataFrame(CSV.File(joinpath(exp_data_dir, "config.csv")))
exp_params_dict = Serialization.deserialize(joinpath(exp_data_dir, "exp.params"))

inertia_values = exp_params_dict[:inertia_values]
γ_vals = exp_params_dict[:γ]
τ_vals = exp_params_dict[:τ]
K_vals = exp_params_dict[:K]
α_vals = exp_params_dict[:α]
freq_bounds = exp_params_dict[:freq_bounds]
failure_modes = exp_params_dict[:failure_modes]
init_pert = exp_params_dict[:init_pert]
number_of_task_ids_between_graphs = length(inertia_values)*length(γ_vals)*length(τ_vals)*length(K_vals)*length(α_vals)*length(freq_bounds)*length(failure_modes)*length(init_pert)

for task_id in df_config.ArrayTaskID
    if ((task_id-1) % number_of_task_ids_between_graphs) == 0
        add_task_id = 10
        N,k,β,graph_seed,μ,σ,distr_seed,K,α,M,γ,τ,freq_bound,trip_lines,trip_nodes,init_pert,ensemble_element = get_network_args_stripped(df_config, task_id + add_task_id)
        monitored_power_flow = exp_params_dict[:monitored_power_flow]
        steadystate_choice = exp_params_dict[:steadystate_choice]

        network = import_system_wrapper(df_config, task_id + add_task_id)
        steady_state_dict  = CSV.File(df_config[task_id + add_task_id,:filepath_steady_state])
        x_static = steady_state_dict[:SteadyState]
        sol = simulate(network;
                       x_static=x_static,
                       initial_fail = [1],
                       init_pert = init_pert,
                       tspan = (0, 0.11),
                       trip_lines = trip_lines,
                       trip_nodes = trip_nodes,
                       trip_load_nodes = :none,
                       monitored_power_flow = monitored_power_flow,
                       f_min = -freq_bound,
                       f_max = freq_bound,
                       solverargs = (;dtmax=0.01),
                       verbose = true);

        # plot network
        tobs = Observable(0.05)
        fig = Figure(resolution=(1800,1800))
        fig[1,1], p = plotnetwork(fig, sol, tobs)
        graph_folder_path = joinpath(exp_data_dir, "k=$k,β=$β", "graphs")
        ispath(graph_folder_path) || mkdir(graph_folder_path)
        CairoMakie.save(joinpath(graph_folder_path, "graph_seed=$graph_seed,distr_seed=$distr_seed,k=$k,β=$β,ensemble_element=$ensemble_element+10.pdf"),fig)
    end
end
