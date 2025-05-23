"""
Plotting Pmech vs Pflow of ensemble
"""

include(abspath(@__DIR__, "..", "helpers_jarray.jl"))
include(abspath(@__DIR__, "WS_trajectories_new_ND_single_model_port.jl"))

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

exp_name_date_dirs = ["WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250321_171511.976"]

for exp_name_date in exp_name_date_dirs
    exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
    df_config = DataFrame(CSV.File(joinpath(exp_data_dir, "config.csv")))

    # adjust filepaths 
    df_config[!, :filepath_graph] = replace.(df_config[!, :filepath_graph],"/home/brandner" => "/home/brandner/nb_data/HU_Master/2122WS/MA")
    df_config[!, :filepath_steady_state] = replace.(df_config[!, :filepath_steady_state],"/home/brandner" => "/home/brandner/nb_data/HU_Master/2122WS/MA")

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
            N,k,β,graph_seed,μ,σ,distr_seed,K,α,M,γ,τ,freq_bound,trip_lines,trip_nodes,init_pert,ensemble_element = get_network_args_stripped(df_config, task_id)

            # TODO Can this be done cleaner without starting a simulation each time by only looking at steady state?
            sol = simulate_new_ND(exp_data_dir, task_id, 1;
                    verbose=true,
                    failtime=0.1,
                    tspan = (0, 0.100001),
                    trip_lines=:dynamic,
                    trip_nodes=:dynamic,
                    terminate_steady_state=true,
                    solverargs=(;),
                    warn=true)

            ###
            ### plot Pmech vs Pflow at src/dst
            ###
            fig = Figure()
            fig[1,1] = ax = Axis(fig; xlabel="Pmech", ylabel="Pflow", title="Correlation between Pflow and Pmech for both ends of each edge")


            # get graph
            nw = NetworkDynamics.extract_nw(sol)
            graph = nw.im.g

            # # old way of retrieving graph
            # df_config = DataFrame(CSV.File(joinpath(exp_data_dir, "config.csv")))
            # df_config[!, :filepath_graph] = replace.(df_config[!, :filepath_graph],"/home/brandner" => "/home/brandner/nb_data/HU_Master/2122WS/MA")
            # df_config[!, :filepath_steady_state] = replace.(df_config[!, :filepath_steady_state],"/home/brandner" => "/home/brandner/nb_data/HU_Master/2122WS/MA")
            # graph = loadgraph(df_config[task_id,:filepath_graph])

            for e_index in 1:ne(graph)
                # get src and dst vertices
                edge = collect(edges(graph))[e_index]
                src_idx = src(edge)
                dst_idx = dst(edge)

                # get Plow at source and dst of edge
                Pflow_src = sol(0, idxs=eidxs(e_index, :P))[1]
                Pflow_dst = - Pflow_src

                # get Pmech at source and dst vertices
                Pmech_src = sol(0, idxs=vidxs(src_idx, :Pmech))[1]
                Pmech_dst = sol(0, idxs=vidxs(dst_idx, :Pmech))[1]

                scatter!(ax, (Pmech_src, Pflow_src); color=:blue, markersize=5)
                scatter!(ax, (Pmech_dst, Pflow_dst); color=:blue, markersize=5)
            end
            dir = joinpath(exp_data_dir, "k=$k,β=$β", "Pmech_vs_Pflow")
            ispath(dir) || mkdir(dir)
            CairoMakie.save(joinpath(dir, "graph_seed=$graph_seed,distr_seed=$distr_seed,k=$k,β=$β,ensemble_element=$ensemble_element.pdf"),fig)
        end
    end
end


###
exp_name_date = "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250321_171511.976"


exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
df_config = DataFrame(CSV.File(joinpath(exp_data_dir, "config.csv")))

sol = simulate_new_ND(exp_data_dir, 1, 1;
    verbose=true,
    failtime=0.1,
    tspan = (0, 0.100001),
    trip_lines=:dynamic,
    trip_nodes=:dynamic,
    terminate_steady_state=true,
    solverargs=(;),
    warn=true)

nw = NetworkDynamics.extract_nw(sol)
graph = nw.im.g

edge = collect(edges(graph))[2]
src_idx = src(edge)
dst_idx = dst(edge)

# get Plow at source and dst of edge
e_index = 2
Pflow_src = sol(0, idxs=eidxs(e_index, :P))
Pflow_dst = - Pflow_src

# get Pmech at source and dst vertices
Pmech_src = sol(0, idxs=vidxs(src_idx, :Pmech))[1]
Pmech_dst = sol(0, idxs=vidxs(dst_idx, :Pmech))[1]