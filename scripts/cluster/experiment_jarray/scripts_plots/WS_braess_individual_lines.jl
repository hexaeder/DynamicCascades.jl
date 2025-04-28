"""
Scatterplot: Braessness of individual lines.
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

"""
Returns for each initially triggered edge the number of line and node failures separately.
"""
function get_line_failures(df_config, task_id)
    N,k,β,graph_seed,μ,σ,distr_seed,K,α,M,γ,τ,freq_bound,trip_lines,trip_nodes,init_pert,ensemble_element = get_network_args_stripped(df_config, task_id)
    graph_combinations_path = joinpath(exp_data_dir, "k=$k,β=$β")
    failure_mode_string = joinpath(graph_combinations_path, "trip_lines=$trip_lines,trip_nodes=$trip_nodes")
    failure_mode_frequ_bound = joinpath(failure_mode_string, "trip_lines=$trip_lines,trip_nodes=$trip_nodes,freq_bound=$freq_bound")
    filename = "trip_lines=$trip_lines,trip_nodes=$trip_nodes,freq_bound=$freq_bound,N=$N,k=$k,β=$β,graph_seed=$graph_seed,μ=$μ,σ=$σ,distr_seed=$distr_seed,K=$K,α=$α,M=$M,γ=$γ,τ=$τ,init_pert=$init_pert,ensemble_element=$ensemble_element.csv"
    df_result = DataFrame(CSV.File(joinpath(graph_combinations_path,failure_mode_string,failure_mode_frequ_bound,filename)))
    number_failures_lines = df_result[!, :number_failures_lines]
    number_failures_nodes = df_result[!, :number_failures_nodes]
    return number_failures_lines, number_failures_nodes
end

"""
Calculates for each edge of the network the Braessness for lines and nodes separately.
"""
function get_braessness(df_config, f_b_narrow, f_b_wide, M, ensemble_element)
    task_id_wide = 0
    line_failures_wide = Float64[]; node_failures_wide = Float64[]
    line_failures_narrow = Float64[]; node_failures_narrow = Float64[]
    for task_id in df_config.ArrayTaskID
        N,k,β,graph_seed,μ,σ,distr_seed,K,α,M_,γ,τ,freq_bound_,trip_lines,trip_nodes,init_pert,ensemble_element_ = get_network_args_stripped(df_config, task_id)

        ## get failures for each line
        # wide bound
        if M==M_&& f_b_wide==freq_bound_ && ensemble_element==ensemble_element_
            line_failures_wide, node_failures_wide = get_line_failures(df_config, task_id)
            println("wide task_id=$task_id")
            #= We need the `task_id` of either wide or narrow for the network metrics `ρ_Pmech_Pflow` and `dist`.
            It doesn't matter if we choose the `task_id` of the wide or the narrow bound as these task_ids refer
            to the same networks.=#
            task_id_wide = task_id 
        end
        # narrow bound
        if M==M_&& f_b_narrow==freq_bound_  && ensemble_element==ensemble_element_
            line_failures_narrow, node_failures_narrow = get_line_failures(df_config, task_id)
            println("narrow task_id=$task_id")
        end
    end
    # braessness_sum_lines_and_nodes = (line_failures_wide .+ node_failures_wide) - (line_failures_narrow .+ node_failures_narrow)
    braessness_lines = line_failures_wide .- line_failures_narrow
    braessness_nodes = node_failures_wide .- node_failures_narrow
    return task_id_wide, braessness_lines, braessness_nodes
end

"""
Calculates for each edge of the network the correlation 
    Pflow_src * (Pmech_src - Pmech_dst)
"""
function ρ_Pmech_Pflow(sol)
    nw = NetworkDynamics.extract_nw(sol)
    graph = nw.im.g
    Pmech_Pflow = []
    for e_index in 1:ne(graph)
        # get src and dst vertices
        edge = collect(edges(graph))[e_index]
        src_idx = src(edge)
        dst_idx = dst(edge)

        # get Plow at source and dst of edge
        Pflow_src = sol(0, idxs=eidxs(e_index, :P))[1] # NOTE `sol(0, idxs=eidxs(e_index, :P))` returns 1-Element vector

        # get Pmech at source and dst vertices
        Pmech_src = sol(0, idxs=vidxs(src_idx, :Pmech))[1]
        Pmech_dst = sol(0, idxs=vidxs(dst_idx, :Pmech))[1]

        # Pmech_Pflow_e_index = Pflow_src * (Pmech_src - Pmech_dst) / abs(Pflow_src)^2
        Pmech_Pflow_e_index = Pflow_src * (Pmech_src - Pmech_dst)
        push!(Pmech_Pflow, Pmech_Pflow_e_index)
    end
    Pmech_Pflow
end

"""
Calculates for each edge of the network the distance between its src and its 
    dst vertex after its removal.
"""
function dist_vertices_after_edge_removal(graph)
    dist = []
    for e_index in 1:ne(graph)
        g = deepcopy(graph)

        # get src and dst vertices
        edge = collect(edges(g))[e_index]

        src_idx = src(edge)
        dst_idx = dst(edge)
        rem_edge!(g, src_idx, dst_idx)

        # compute shortest path distances from vertex src_idx
        state = dijkstra_shortest_paths(g, src_idx)
        # retrieve the distance from vertex src_idx to vertex dst_idx
        dist_e_index = state.dists[dst_idx]
        push!(dist, dist_e_index)
    end
    dist
end





exp_name_date = "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250321_171511.976"
exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
df_config = DataFrame(CSV.File(joinpath(exp_data_dir, "config.csv")))

# adjust filepaths 
df_config[!, :filepath_graph] = replace.(df_config[!, :filepath_graph],"/home/brandner" => "/home/brandner/nb_data/HU_Master/2122WS/MA")
df_config[!, :filepath_steady_state] = replace.(df_config[!, :filepath_steady_state],"/home/brandner" => "/home/brandner/nb_data/HU_Master/2122WS/MA")

exp_params_dict = Serialization.deserialize(joinpath(exp_data_dir, "exp.params"))
inertia_values = exp_params_dict[:inertia_values]
freq_bounds = exp_params_dict[:freq_bounds]
N_ensemble_size = exp_params_dict[:N_ensemble_size]
inertia_values
# choose inertia
M = 7.5
# choose `f_b_narrow` and `f_b_wide`
f_b_narrow = 0.03
f_b_wide = 0.15

################################################################################
fontsize = 25
titlesize = (fontsize-5)
markersize = 8
model = "BH+Pmech"
# model = exp_params_dict[:node_failure_model]

fig = Figure(size=(2100,1500), fontsize= fontsize)
fig[1,1] = ax11 = Axis(fig; ylabel="N_wide-N_narrow: node & line failures", title="Braessness vs. ρ_Pmech_Pflow: $model,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$M", titlealign = :left, titlesize = titlesize)
fig[2,1] = ax21 = Axis(fig; ylabel="N_wide-N_narrow: line failures")
fig[3,1] = ax31 = Axis(fig; xlabel="ρ_Pmech_Pflow", ylabel="N_wide-N_narrow: node failures")
fig[1,2] = ax12 = Axis(fig; ylabel="N_wide-N_narrow: node & line failures", title="Braessness vs. dist(src,dst): $model,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$M", titlealign = :left, titlesize = titlesize)
fig[2,2] = ax22 = Axis(fig; ylabel="N_wide-N_narrow: line failures")
fig[3,2] = ax32 = Axis(fig; xlabel="dist(src,dst)", ylabel="N_wide-N_narrow: node failures")

# for ensemble_element in 1:1
for ensemble_element in 1:N_ensemble_size
    ## y-Axis
    task_id_wide, braessness_lines, braessness_nodes = get_braessness(df_config, f_b_narrow, f_b_wide, M, ensemble_element)
    # sum: `braessness_lines .+ braessness_nodes``

    ## x-Axis
    sol = simulate_new_ND_single_model_port(exp_data_dir, task_id_wide, 1;
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
    x_ρ = ρ_Pmech_Pflow(sol)
    x_dist = dist_vertices_after_edge_removal(graph)

    scatter!(ax11, x_ρ, braessness_lines .+ braessness_nodes; color=:blue, markersize=markersize)
    scatter!(ax21, x_ρ, braessness_lines; color=:orange, markersize=markersize)
    scatter!(ax31, x_ρ, braessness_nodes; color=:red, markersize=markersize)
    scatter!(ax12, x_dist, braessness_lines .+ braessness_nodes; color= (:blue, 0.1), markersize=markersize+10)
    scatter!(ax22, x_dist, braessness_lines; color= (:orange, 0.1), markersize=markersize+10)
    scatter!(ax32, x_dist, braessness_nodes; color = (:red, 0.1), markersize=markersize+10)

end
fig

CairoMakie.save(joinpath(exp_data_dir, "braessness_lines", "braessness_$model,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$M,normalized.pdf"),fig)
CairoMakie.save(joinpath(exp_data_dir, "braessness_lines", "braessness_$model,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$M,normalized.png"),fig)

