"""
Plotting graphs and initial power flow of ensemble
"""

include(abspath(@__DIR__, "..", "helpers_jarray.jl"))

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
CairoMakie.activate!()
using CSV
using DataFrames
using Serialization

function plotnetwork(network, fig, sol, t)
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


exp_name_date_dirs = ["WS_k=4_exp03_1_vary_I_only_lines_PIK_HPC_K_=3,N_G=32_20250326_231646.227",
    "WS_k=4_exp03_2_vary_I_only_nodes_PIK_HPC_K_=3,N_G=32_20250326_231746.273",
    "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250321_171511.976",
    "WS_k=4_exp05_1_I_over_Dsq_lines_PIK_HPC_K_=3,N_G=32_20250326_230911.466",
    "WS_k=4_exp05_2_I_over_Dsq_nodes_PIK_HPC_K_=3,N_G=32_20250326_230941.02",
    "WS_k=4_exp05_3_I_over_Dsq_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250326_230741.529",
    "WS_k=4_exp06_1_I_over_D_lines_PIK_HPC_K_=3,N_G=32_20250326_230610.766",
    "WS_k=4_exp06_2_I_over_Dsq_nodes_PIK_HPC_K_=3,N_G=32_20250326_230642.483",
    "WS_k=4_exp06_3_I_over_D_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250326_230522.458",
    "WS_k=4_exp07_1_vary_D_only_lines_PIK_HPC_K_=3,N_G=32_20250326_231429.372",
    "WS_k=4_exp07_2_vary_D_only_nodes_PIK_HPC_K_=3,N_G=32_20250326_231547.118",
    "WS_k=4_exp07_3_vary_D_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250326_231211.092",
    "WS_k=4_exp08_vary_alpha_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250326_231041.385"
    ]

for exp_name_date in exp_name_date_dirs
    exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
    df_config = DataFrame(CSV.File(joinpath(exp_data_dir, "config.csv")))

    # adjust filepaths 
    df_config[!, :filepath_graph] = replace.(df_config[!, :filepath_graph],"/home/brandner" => "/home/brandner/nb_data/HU_Master/2122WS/MA")
    df_config[!, :filepath_steady_state] = replace.(df_config[!, :filepath_steady_state],"/home/brandner" => "/home/brandner/nb_data/HU_Master/2122WS/MA")

    exp_params_dict = Serialization.deserialize(joinpath(exp_data_dir, "exp.params"))


    # TODO better use:
    # number_of_task_ids_between_graphs = length(df_config.ArrayTaskID)/N_ensemble_size
    # for task_id in 1:number_of_task_ids_between_graphs:length(df_config.ArrayTaskID)
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
            N,k,β,graph_seed,μ,σ,distr_seed,K,α,M,γ,τ,freq_bound,trip_lines,trip_nodes,init_pert,ensemble_element = get_network_args_stripped(df_config, task_id)
            monitored_power_flow = exp_params_dict[:monitored_power_flow]

            network = import_system_wrapper(df_config, task_id)
            # NOTE 
            #= `simulate_new_ND_single_model_port` might use a steadystate that is
            globally phase shifted. This shouldn't make a difference for the flows and the
            simulations however. 
            =#
            steady_state_dict  = CSV.File(df_config[task_id,:filepath_steady_state])
            x_static = steady_state_dict[:SteadyState]
            sol = simulate(network;
                        x_static=x_static,
                        initial_fail = [1],
                        init_pert = init_pert,
                        tspan = (0, 0.100001),
                        trip_lines = trip_lines,
                        trip_nodes = trip_nodes,
                        trip_load_nodes = :none,
                        monitored_power_flow = monitored_power_flow,
                        f_min = -freq_bound,
                        f_max = freq_bound,
                        solverargs = (;dtmax=0.01),
                        verbose = true);

            ###
            ### plot network
            ###
            tobs = Observable(0.05)
            fig = Figure(resolution=(1800,1800))
            fig[1,1], p = plotnetwork(network, fig, sol, tobs)
            dir = joinpath(exp_data_dir, "k=$k,β=$β", "graphplots")
            ispath(dir) || mkdir(dir)
            CairoMakie.save(joinpath(dir, "graph_seed=$graph_seed,distr_seed=$distr_seed,k=$k,β=$β,ensemble_element=$ensemble_element.pdf"),fig)
        end
    end
end