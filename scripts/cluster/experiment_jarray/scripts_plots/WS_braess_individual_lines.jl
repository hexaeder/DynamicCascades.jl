"""
Plotting scripts: Braessness of individual lines.
"""

include(abspath(@__DIR__, "..", "helpers_jarray.jl"))
include("/home/brandner/.julia/dev/DynamicCascades/src/ND_model_new_ND.jl")

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

################################################################################
############################## Data extraction #################################
################################################################################

"""
Returns for each initially triggered edge the number of line and node failures separately.
"""
function get_line_failures(df_config, task_id)
    N,k,β,graph_seed,μ,σ,distr_seed,K,α,I,γ,τ,freq_bound,trip_lines,trip_nodes,init_pert,ensemble_element = get_network_args_stripped(df_config, task_id)
    graph_combinations_path = joinpath(exp_data_dir, "k=$k,β=$β")
    failure_mode_string = joinpath(graph_combinations_path, "trip_lines=$trip_lines,trip_nodes=$trip_nodes")
    failure_mode_frequ_bound = joinpath(failure_mode_string, "trip_lines=$trip_lines,trip_nodes=$trip_nodes,freq_bound=$freq_bound")
    filename = "trip_lines=$trip_lines,trip_nodes=$trip_nodes,freq_bound=$freq_bound,N=$N,k=$k,β=$β,graph_seed=$graph_seed,μ=$μ,σ=$σ,distr_seed=$distr_seed,K=$K,α=$α,M=$I,γ=$γ,τ=$τ,init_pert=$init_pert,ensemble_element=$ensemble_element.csv"
    df_result = DataFrame(CSV.File(joinpath(graph_combinations_path,failure_mode_string,failure_mode_frequ_bound,filename)))
    number_failures_lines = df_result[!, :number_failures_lines]
    number_failures_nodes = df_result[!, :number_failures_nodes]
    return number_failures_lines, number_failures_nodes
end

"""
Get taks_id for given parameters.
"""
function get_task_id(df_config, f_b, I, ensemble_element)
    for task_id in df_config.ArrayTaskID
        _,_,_,_,_,_,_,_,_,I_,_,_,freq_bound_,_,_,_,ensemble_element_ = get_network_args_stripped(df_config, task_id)
        if I==I_&& f_b==freq_bound_ && ensemble_element==ensemble_element_
            println("task_id=$task_id")
            return task_id
        end
    end
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

"""
Saves a df for each frequency bound f_b  in `freq_bounds` (for given inertia I)
Output-columns in df: ArrayTaskID,ensemble_element,initially_failed_line,number_failures_lines,number_failures_nodes
(This function only preprocesses the result data)
"""
function write_failures_f_b_to_df(exp_name_date, freq_bounds, I)
    exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
    df_config = DataFrame(CSV.File(joinpath(exp_data_dir, "config.csv")))

    # adjust filepaths 
    df_config[!, :filepath_graph] = replace.(df_config[!, :filepath_graph],"/home/brandner" => "/home/brandner/nb_data/HU_Master/2122WS/MA")
    df_config[!, :filepath_steady_state] = replace.(df_config[!, :filepath_steady_state],"/home/brandner" => "/home/brandner/nb_data/HU_Master/2122WS/MA")

    exp_params_dict = Serialization.deserialize(joinpath(exp_data_dir, "exp.params"))
    N_ensemble_size = exp_params_dict[:N_ensemble_size]

    dir = joinpath(exp_data_dir, "braessness_lines", "braessness_data")
    ispath(dir) || mkdir(dir)

    for f_b in freq_bounds
        # create df
        df = DataFrame(
            ArrayTaskID = Int[],
            ensemble_element = Int[],
            initially_failed_line = Int[],
            number_failures_lines = Int[],  
            number_failures_nodes = Int[]  
        )

        # stack ensemble elements in df
        for ensemble_element in 1:N_ensemble_size
            ## y-Axis
            task_id = get_task_id(df_config, f_b, I, ensemble_element)
            line_failures, node_failures = get_line_failures(df_config, task_id)

            # number of edges
            ne_graph = length(line_failures)
            lines = collect(1:ne_graph)

            # build a small DataFrame for this ensemble member
            df_iter = DataFrame(
                ArrayTaskID = fill(task_id, ne_graph),
                ensemble_element = fill(ensemble_element, ne_graph),
                initially_failed_line = lines,
                number_failures_lines = line_failures,
                number_failures_nodes = node_failures
            )

            # append to master dataframe
            append!(df, df_iter)
        end
        CSV.write(joinpath(exp_data_dir,"braessness_lines", "braessness_data", "I=$I,f_b=$f_b.csv"), df)
    end
end

"""
Saves a df with the network measures for the whole ensemble.
Output-columns in df: ensemble_element,initially_failed_line,rho,dist
"""
function write_network_measures_to_df(exp_name_date)
    exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
    df_config = DataFrame(CSV.File(joinpath(exp_data_dir, "config.csv")))
    
    # adjust filepaths 
    df_config[!, :filepath_graph] = replace.(df_config[!, :filepath_graph],"/home/brandner" => "/home/brandner/nb_data/HU_Master/2122WS/MA")
    df_config[!, :filepath_steady_state] = replace.(df_config[!, :filepath_steady_state],"/home/brandner" => "/home/brandner/nb_data/HU_Master/2122WS/MA")
    
    exp_params_dict = Serialization.deserialize(joinpath(exp_data_dir, "exp.params"))
    N_ensemble_size = exp_params_dict[:N_ensemble_size]

    df = DataFrame(
        ensemble_element = Int[],
        initially_failed_line = Int[],
        rho = Float64[],  
        dist = Int[]  
    )

    number_of_task_ids_between_graphs = length(df_config.ArrayTaskID)/N_ensemble_size
    task_ids = 1:number_of_task_ids_between_graphs:length(df_config.ArrayTaskID)
    for (ensemble_element, task_id) in enumerate(task_ids)
        #= #HACK (that saves effort) `sol` passed to `ρ_Pmech_Pflow`. For `ρ_Pmech_Pflow`
        one could also only generate `nw` and then use `s0 = NWState(nw, x_static, pflat(NWParameter(nw)))`.
        =#
        sol = simulate_new_ND(exp_name_date, Int(task_id), 1;
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

        # number of edges
        ne_graph = ne(graph)
        lines = collect(1:ne_graph)

        # build a small DataFrame for this ensemble member
        df_iter = DataFrame(
            ensemble_element = fill(ensemble_element, ne_graph),
            initially_failed_line = lines,
            rho = x_ρ,  
            dist = x_dist 
        )

        # append to master dataframe
        append!(df, df_iter)
    end
    CSV.write(joinpath(exp_data_dir,"braessness_lines", "braessness_data", "measures_ensemble.csv"), df)
end

"""
loads dataframes with respective frequencies and calculates braessness
"""
function get_braessness(exp_name_date, f_b_narrow, f_b_wide, I)
    exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
    df_narrow = DataFrame(CSV.File(joinpath(exp_data_dir,"braessness_lines", "braessness_data", "I=$I,f_b=$f_b_narrow.csv")))
    df_wide = DataFrame(CSV.File(joinpath(exp_data_dir,"braessness_lines", "braessness_data", "I=$I,f_b=$f_b_wide.csv")))
    braessness_lines = df_wide.number_failures_lines .- df_narrow.number_failures_lines
    braessness_nodes = df_wide.number_failures_nodes .- df_narrow.number_failures_nodes

    return braessness_lines, braessness_nodes
end

"""
loads dataframe with network measures
"""
function get_network_measures(exp_name_date)
    exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
    df_measures = DataFrame(CSV.File(joinpath(exp_data_dir,"braessness_lines", "braessness_data", "measures_ensemble.csv")))
    x_ρ = df_measures.rho
    x_dist = df_measures.dist
    return x_ρ, x_dist
end


################################################################################
############################## Plotting functions###############################
################################################################################
"""
freqs maps each point (x,y) (stored as a Tuple{Int,Int}) to how many times it occurs.
`freqs[(x,y)] = get(freqs, (x,y), 0) + 1`
`get(collection, key, default)` Return the value stored for the given key, or the given default value if no mapping for the key is present.
If the key is already in freqs, you get its stored count.
If it’s not yet in freqs, you get the default 0.
You then add 1 to that number, producing the updated count.
freqs[(x,y)] = … assigns it back into the dictionary.
"""
function get_log_count(xs::AbstractVector{Int}, ys::AbstractVector{Int})
    @assert length(xs) == length(ys)
    freqs = Dict{Tuple{Int,Int}, Int}()
    for (x,y) in zip(xs, ys)
        freqs[(x,y)] = get(freqs, (x,y), 0) + 1
    end

    # build a vector of log10(count) for coloring
    counts = [freqs[(x,y)] for (x,y) in zip(xs, ys)]
    return log10.(counts)
end
# 3D version
function get_log_count(xs::AbstractVector{Int}, ys::AbstractVector{Int}, zs::AbstractVector{Int})
    @assert length(xs) == length(ys) == length(zs)
    freqs = Dict{NTuple{3,Int}, Int}()
    for (x,y,z) in zip(xs, ys, zs)
        freqs[(x,y,z)] = get(freqs, (x,y,z), 0) + 1
    end
    counts = [freqs[(x,y,z)] for (x,y,z) in zip(xs, ys, zs) ]
    return log10.(counts)
end


function colorswitcher(z; fancy_colors=true)
    lo, hi = minimum(z), maximum(z) # compute bounds
    colorrange = (lo, hi)
    color = z
    mid = 0  # the point where colors flip
    # Build a diverging colormap that places the “neutral” color at z = 0 see https://github.com/MakieOrg/Makie.jl/issues/2485

    if fancy_colors == true
        colormap = Makie.diverging_palette(230, 0; mid = (mid - lo)/(hi - lo), w=0.8, d2=0.7, s=0.9)
        strokewidth = 1.0
    else
        colormap = :viridis
        strokewidth = 0.0
    end
    return color, colormap, colorrange, strokewidth
end


################################################################################
############################## Plotting scripts ################################
################################################################################
fontsize = 25
titlesize = (fontsize-5)
markersize = 8

###
### preprocess and load data
###

######
###### WS
######
exp_name_date = "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250321_171511.976"
exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)

# choose inertia
I = 7.5
freq_bounds = [0.03, 0.15]

# write_failures_f_b_to_df(exp_name_date, freq_bounds, I)
# write_network_measures_to_df(exp_name_date) # NOTE For WS set `res_tol=1e-5,` in `steadystate_new_ND`

# choose `f_b_narrow` and `f_b_wide`
f_b_narrow = 0.03
f_b_wide = 0.15
braessness_lines, braessness_nodes = get_braessness(exp_name_date, f_b_narrow, f_b_wide, I);
x_ρ, x_dist = get_network_measures(exp_name_date);

# model = Serialization.deserialize(joinpath(RESULTS_DIR, exp_name_date, "exp.params"))[:node_failure_model]
model = "full_failure"
N_nodes = Serialization.deserialize(joinpath(RESULTS_DIR, exp_name_date, "exp.params"))[:N_nodes]


###
### Braessness histograms
###
fig = Figure(size = (2100,1500), fontsize = fontsize)
xscale = Makie.pseudolog10  # maps 0→0, then log10(count+1)
# Top‐left: log‐count histogram of x_ρ
fig[1, 1] = ax11 = Axis(fig; yscale=xscale, ylabel = "count", titlealign = :left, titlesize  = titlesize,
    title  = "f_b_n=$f_b_narrow, f_b_w=$f_b_wide, I=$I")
# Bottom‐left: scatter ρ vs Braessness
fig[2, 1] = ax21 = Axis(fig; xlabel = "ρ_Pmech_Pflow", ylabel = "Braessness")
# Bottom‐right: horizontal log‐count histogram of Braessness
fig[2, 2] = ax22 = Axis(fig; xscale=xscale, xlabel = "count")
fig[2, 3] = ax23 = Axis(fig; xscale=xscale, xlabel = "count")
fig[2, 4] = ax24 = Axis(fig; xscale=xscale, xlabel = "count")

# link the marginal histograms to the scatter
linkxaxes!(ax21, ax11); linkyaxes!(ax21, ax22, ax23, ax24); linkxaxes!(ax22, ax23, ax24)

# center scatter
labels = ["line & node failures", "line failures", "node failures"]
ys     = [braessness_lines .+ braessness_nodes, braessness_lines, braessness_nodes]
cols   = (:blue, :orange, :red)
for (lbl, y, col) in zip(labels, ys, cols)
    scatter!(ax21, x_ρ, y; label=lbl, color=col, markersize = markersize)
end

# top histogram: counts of x_ρ
strokewidth = 0
μ11= round(mean(x_ρ), digits=3)
hist!(ax11, x_ρ; bins=200, label="<ρ_Pmech_Pflow>=$μ11", color=:black, strokewidth=strokewidth)

# bottom right histogram: counts of Braessness
bins22 = maximum(braessness_lines .+ braessness_nodes) - minimum(braessness_lines .+ braessness_nodes)
bins23 = maximum(braessness_lines) - minimum(braessness_lines)
bins24 = maximum(braessness_nodes) - minimum(braessness_nodes)

μ22 = round(mean(braessness_lines .+ braessness_nodes), digits=3)
μ23 = round(mean(braessness_lines), digits=3)
μ24 = round(mean(braessness_nodes), digits=3)

hist!(ax22, braessness_lines .+ braessness_nodes; label="<Braessness>=$μ22", color=:blue , bins=bins22, direction=:x, strokewidth=strokewidth)
hist!(ax23, braessness_lines; label="<Braessness>=$μ23", color=:orange, bins=bins23, direction=:x, strokewidth=strokewidth)
hist!(ax24, braessness_nodes; label="<Braessness>=$μ24", color=:red , bins=bins24, direction=:x, strokewidth=strokewidth)

hidexdecorations!(ax11)
hideydecorations!(ax22); hideydecorations!(ax23); hideydecorations!(ax24)
axislegend(ax11; position = :rt, labelsize=fontsize)
axislegend(ax21; position = :rt, labelsize=fontsize-8)
axislegend(ax22; position = :rt, labelsize=fontsize)
axislegend(ax23; position = :rt, labelsize=fontsize)
axislegend(ax24; position = :rt, labelsize=fontsize)
ax11.yticks = [0, 1, 10, 100, 600]
ax22.xticks = ax23.xticks = ax24.xticks = [0, 1, 10, 100, 1000, 5000]
# CairoMakie.save(joinpath(exp_data_dir, "braessness_lines", "scatter_histograms_and_rho_braessness_$model,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.pdf"),fig)
# CairoMakie.save(joinpath(exp_data_dir, "braessness_lines", "scatter_histograms_and_rho_braessness_$model,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.png"),fig)
fig


###
### histogramms for x_dist
###
fig = Figure(size = (2100, 1700), fontsize = fontsize)
rows = [1,2,3]
data_arr = [(braessness_lines .+ braessness_nodes), braessness_lines, braessness_nodes]
ylabels = ["Braessness lines and nodes", "Braessness lines", "Braessness nodes"]
colors = [:blue, :orange, :red]

for (row, data, ylabel, color) in zip(rows, data_arr, ylabels, colors)
    # get categories (here the different distances)
    cats = sort(unique(x_dist))

    # Create a 1×length(cats) grid of Axes for all n categories
    axs = [Axis(fig[row, c];
                xscale=Makie.pseudolog10,
                xlabel = c == cats[length(cats)÷2] ? "Count" : "",
                ylabel = c == minimum(cats) ? ylabel : "",
                title  = row == 1 ? "d=$c" : ""
                )
            for c in cats]

    # Optionally link all x-axes so panning one pans the rest
    linkyaxes!(axs...); linkxaxes!(axs...)
    row < maximum(rows) ? hidexdecorations!.(axs[1:end], grid=false) : nothing
    hideydecorations!.(axs[2:end])

    # Draw each histogram in its own Axis
    for (c, ax) in zip(cats, axs)
        # Braessness of the lines with dist = c
        ys = data[findall(x -> x == c, x_dist)]
        μ = round(mean(ys), digits=3)
        min_ys = 
        hist!(ax, ys;
            bins      = maximum(ys) - minimum(ys)+1,
            direction = :x,
            label     = "<Braessness>=$μ",
            color     = color,
            strokewidth = 0)
        ax.xticks = [0, 1, 10, 100, 1000]
    end

    [axislegend(i; position = :rt, labelsize=fontsize-6) for i in axs]
end
# CairoMakie.save(joinpath(exp_data_dir, "braessness_lines", "braessness_dist_braessness_histogram_$model,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.pdf"),fig)
# CairoMakie.save(joinpath(exp_data_dir, "braessness_lines", "braessness_dist_braessness_histogram_$model,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.png"),fig)
fig


### 
### Braessness vs Braessness
###

## BH+Pmech
exp_name_date = "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250321_171511.976"
exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
braessness_lines, braessness_nodes = get_braessness(exp_name_date, f_b_narrow, f_b_wide, I);
x_ρ, x_dist = get_network_measures(exp_name_date);
exp_nr_full = exp_name_date[11:12]

## BH only
exp_name_date = "WS_k=4_exp11_vary_I_only_lines_and_nodes_change_to_BH_complement_f_bPIK_HPC_K_=3,N_G=32_20250420_145304.04"
exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
braessness_lines_BH, braessness_nodes_BH = get_braessness(exp_name_date, f_b_narrow, f_b_wide, I);
x_ρ_BH, x_dist_BH = get_network_measures(exp_name_date);
exp_nr_BH = exp_name_date[11:12]

## Pmech only
exp_name_date = "WS_k=4_exp12_vary_I_only_lines_and_nodes_change_Pmech_complement_f_bPIK_HPC_K_=3,N_G=32_20250420_142151.671"
exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
braessness_lines_Pmech, braessness_nodes_Pmech = get_braessness(exp_name_date, f_b_narrow, f_b_wide, I);
x_ρ_Pmech, x_dist_Pmech = get_network_measures(exp_name_date);
exp_nr_Pmech = exp_name_date[11:12]

prefix = exp_name_date[1:10]
exp_nrs = "$prefix,$exp_nr_full,$exp_nr_BH,$exp_nr_Pmech"

###
### 2D-Plots with colorscaling for z-Axis leaving out scatterpoints
###


# data
xs = braessness_lines_Pmech .+ braessness_nodes_Pmech
ys = braessness_lines_BH .+ braessness_nodes_BH
zs = braessness_lines .+ braessness_nodes
failures = "line+node"

# xs = braessness_lines_Pmech
# ys = braessness_lines_BH
# zs = braessness_lines
# failures = "lines"

# xs = braessness_nodes_Pmech
# ys = braessness_nodes_BH
# zs = braessness_nodes
# failures = "nodes"

dimensions = 3 # 2D vs 3D selection of scatterpoints
plot_center_only = true; lim = 10.5 # axis limits
count = 2 # plot only scatterpoints that occur less often than `count`
r = 2 # plot only scatterpoints with (|power failure|,|inertia failure|,|full failure|) > r
fancy_colors = true # two different color codes for positive and negative Braessness
markersize = 15

fig = Figure(size=(1000,1500),fontsize=fontsize)
# No scatterpoints left out
fig[1,1] = ax11 = Axis(fig; title="Braessness $failures failures: f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I (all datapoints)", titlealign = :left, titlesize = titlesize)
color, colormap, colorrange, strokewidth = colorswitcher(zs; fancy_colors)
sc11 = scatter!(ax11, xs, ys; markersize=markersize, color=color, colormap=colormap, colorrange=colorrange, strokecolor=:black, strokewidth=strokewidth)
Colorbar(fig[1,2], sc11; label = "full failure", width = 30)

# Only scatterpoints with 2D/3D count
dimensions == 3 ? note="(there might be points with count > 0)" : note=""
fig[2,1] = ax21 = Axis(fig; ylabel="inertia failure", title="Only scatterpoints with count < $count in $dimensions-D (outliers) $note", titlealign = :left, titlesize = titlesize)
freqs2D = Dict{Tuple{Int,Int}, Int}()
for (x,y) in zip(xs, ys)
    freqs2D[(x,y)] = get(freqs2D, (x,y), 0) + 1
end
freqs3D = Dict{NTuple{3,Int}, Int}()
for (x,y,z) in zip(xs, ys, zs)
    freqs3D[(x,y,z)] = get(freqs3D, (x,y,z), 0) + 1
end
xs_outlier = Int64[]; ys_outlier = Int64[]; zs_outlier = Int64[]
for (x,y,z) in zip(xs, ys, zs)
    if dimensions==2 && freqs2D[(x,y)] < count
        push!(xs_outlier, x); push!(ys_outlier, y); push!(zs_outlier, z)
    end
    if dimensions==3 && freqs3D[(x,y,z)] < count
        push!(xs_outlier, x); push!(ys_outlier, y); push!(zs_outlier, z)
    end
end
color, colormap, colorrange, strokewidth = colorswitcher(zs_outlier; fancy_colors)
sc21 = scatter!(ax21, xs_outlier, ys_outlier; markersize=markersize, color=color, colormap=colormap, colorrange=colorrange, strokecolor=:black, strokewidth=strokewidth)
Colorbar(fig[2,2], sc21; label = "full failure", width = 30)

# Only scatterpoints with (|power failure|,|inertia failure|,|full failure|) > r=0,1,2
dimensions == 3 ? dim_r="(|power failure| or |inertia failure| or |full failure|)" : dim_r="(|power failure| or |inertia failure|)"
fig[3,1] = ax31 = Axis(fig; xlabel="power failure", title="Only scatterpoints with $dim_r > $r (there might be points with count > 1)", titlealign = :left, titlesize = titlesize-5)
xs_wo_center = Int64[]; ys_wo_center = Int64[]; zs_wo_center = Int64[]
for (x,y,z) in zip(xs, ys, zs)
    if dimensions==2 && (abs(x)>r || abs(y)>r)
        push!(xs_wo_center, x); push!(ys_wo_center, y); push!(zs_wo_center, z)
    end
    if dimensions==3 && (abs(x)>r || abs(y)>r || abs(z)>r)
        push!(xs_wo_center, x); push!(ys_wo_center, y); push!(zs_wo_center, z)
    end
end
color, colormap, colorrange, strokewidth = colorswitcher(zs_wo_center; fancy_colors)
sc31 = scatter!(ax31, xs_wo_center, ys_wo_center; markersize=markersize, color=color, colormap=colormap, colorrange=colorrange, strokecolor=:black, strokewidth=strokewidth)
Colorbar(fig[3,2], sc11; label = "full failure", width = 30)

if plot_center_only == true
    xlims!(ax11, -lim, lim); ylims!(ax11, -lim, lim)
    xlims!(ax21, -lim, lim); ylims!(ax21, -lim, lim)
    xlims!(ax31, -lim, lim); ylims!(ax31, -lim, lim)
end
# CairoMakie.save(joinpath(MA_DIR, "braessness", "exp=$exp_nrs,braessness_vs_braessness_failures=$failures,logcount_colorscale,,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I,dim=$dimensions,fancy=$fancy_colors,center=$plot_center_only.pdf"),fig)
# CairoMakie.save(joinpath(MA_DIR, "braessness", "exp=$exp_nrs,braessness_vs_braessness_failures=$failures,logcount_colorscale,,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I,dim=$dimensions,fancy=$fancy_colors,center=$plot_center_only.png"),fig)
fig


###
### 3D scatter plot counts
###
using GLMakie # for interactivity
GLMakie.activate!()
fontsize = 25
titlesize = (fontsize-5)
markersize = 12

# lines and nodes
xs = braessness_lines_BH .+ braessness_nodes_BH
ys = braessness_lines_Pmech .+ braessness_nodes_Pmech
zs = braessness_lines .+ braessness_nodes

fig = Figure(size=(1470,1050), fontsize=fontsize)
fig[1,1] = ax11 = Axis3(fig; xlabel="inertia failure", ylabel="power failure", zlabel="full failure", title="Braessness node & line failures: f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I", titlealign = :left, titlesize = titlesize)
logc = get_log_count(xs, ys, zs)
sc11 = scatter!(ax11, xs, ys, zs; markersize=markersize, color=logc, colormap=:viridis, colorrange=(minimum(logc), maximum(logc)))
Colorbar(fig[1,2], sc11; label = "log₁₀(counts)", width = 30)
# CairoMakie.save(joinpath(MA_DIR, "braessness", "braessness_vs_braessness_lines_and_nodes_3D_logcount_colorscale,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.pdf"),fig)
# CairoMakie.save(joinpath(MA_DIR, "braessness", "braessness_vs_braessness_lines_and_nodes_3D_logcount_colorscale,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.png"),fig)
fig

###
### 3D scatter plot only scatterpoints with count=1
###
freqs3D = Dict{NTuple{3,Int}, Int}()
for (x,y,z) in zip(xs, ys, zs)
    freqs3D[(x,y,z)] = get(freqs3D, (x,y,z), 0) + 1
end

fig = Figure(size=(1470,1050), fontsize=fontsize)
fig[1,1] = ax11 = Axis3(fig; xlabel="inertia failure", ylabel="power failure", zlabel="full failure", title="Braessness node & line failures: f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I", titlealign = :left, titlesize = titlesize)

xs_outlier = Int64[]; ys_outlier = Int64[]; zs_outlier = Int64[]
for (x,y,z) in zip(xs, ys, zs)
    if freqs3D[(x,y,z)] < 2
        push!(xs_outlier, x); push!(ys_outlier, y); push!(zs_outlier, z)
    end
end
color, colormap, colorrange, strokewidth = colorswitcher(zs_outlier; fancy_colors)
sc11 = scatter!(ax11, xs_outlier, ys_outlier, zs_outlier; markersize=markersize, color=color, colormap=colormap, colorrange=colorrange, strokewidth=strokewidth)
Colorbar(fig[1,2], sc11; label = "full failure", width = 30)
# CairoMakie.save(joinpath(MA_DIR, "braessness", "braessness_vs_braessness_lines_and_nodes_3D_logcount_colorscale,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.pdf"),fig)
# CairoMakie.save(joinpath(MA_DIR, "braessness", "braessness_vs_braessness_lines_and_nodes_3D_logcount_colorscale,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.png"),fig)
fig

################################################################################
########################## Old plotting scripts ################################
################################################################################

###
### first scatter plot
###
fontsize = 25
titlesize = (fontsize-5)
markersize = 8
model = "full failure"
# model = exp_params_dict[:node_failure_model]

fig = Figure(size=(2100,1500), fontsize= fontsize)
fig[1,1] = ax11 = Axis(fig; ylabel="N_wide-N_narrow: node & line failures", title="Braessness vs. ρ_Pmech_Pflow: $model,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I", titlealign = :left, titlesize = titlesize)
fig[2,1] = ax21 = Axis(fig; ylabel="N_wide-N_narrow: line failures")
fig[3,1] = ax31 = Axis(fig; xlabel="ρ_Pmech_Pflow", ylabel="N_wide-N_narrow: node failures")
fig[1,2] = ax12 = Axis(fig; ylabel="N_wide-N_narrow: node & line failures", title="Braessness vs. dist(src,dst): $model,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I", titlealign = :left, titlesize = titlesize)
fig[2,2] = ax22 = Axis(fig; ylabel="N_wide-N_narrow: line failures")
fig[3,2] = ax32 = Axis(fig; xlabel="dist(src,dst)", ylabel="N_wide-N_narrow: node failures")

scatter!(ax11, x_ρ, braessness_lines .+ braessness_nodes; color=:blue, markersize=markersize)
scatter!(ax21, x_ρ, braessness_lines; color=:orange, markersize=markersize)
scatter!(ax31, x_ρ, braessness_nodes; color=:red, markersize=markersize)
scatter!(ax12, x_dist, braessness_lines .+ braessness_nodes; color= (:blue, 0.1), markersize=markersize+10)
scatter!(ax22, x_dist, braessness_lines; color= (:orange, 0.1), markersize=markersize+10)
scatter!(ax32, x_dist, braessness_nodes; color = (:red, 0.1), markersize=markersize+10)
# CairoMakie.save(joinpath(exp_data_dir, "braessness_lines", "test_braessness_$model,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.pdf"),fig)
# CairoMakie.save(joinpath(exp_data_dir, "braessness_lines", "test_braessness_$model,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.png"),fig)
fig

###
### density plots Braessness vs ρ_Pmech_Pflow
###
fig = Figure(size=(2100,1500), fontsize= fontsize)
fig[1,1] = ax11 = Axis(fig; ylabel="Density", title="Braessness line failures: f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I", titlealign = :left, titlesize = titlesize)
fig[2,1] = ax21 = Axis(fig; xlabel="ρ_Pmech_Pflow", ylabel="Braessness")
fig[2,2] = ax22 = Axis(fig; xlabel="Density")

# link their axes
linkxaxes!(ax21, ax11)
linkyaxes!(ax21, ax22)

labels = ["line & node failures", "line failures", "node failures"]
ys     = [braessness_lines .+ braessness_nodes, braessness_lines, braessness_nodes]
cols   = (:blue, :orange, :red)

for (lbl, y, col) in zip(labels, ys, cols)
    density!(ax11, x_ρ; color=col)
    scatter!(ax21, x_ρ, y; label=lbl, color=col, markersize=markersize)
    density!(ax22, y; direction=:y, color=col)
end

hidexdecorations!(ax11)
hideydecorations!(ax22)
axislegend(ax21)
fig

###
### violin plot
###
scale = :count
gap = -0.2
ymin = -3
ymax = 3
# model = exp_params_dict[:node_failure_model]

fig = Figure(size=(2100,1500), fontsize= fontsize)
fig[1,1] = ax11 = Axis(fig; ylabel="N_wide-N_narrow: node & line failures", title="Braessness vs. dist(src,dst): $model,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I,scale=$scale", titlealign = :left, titlesize = titlesize)
fig[2,1] = ax21 = Axis(fig; ylabel="N_wide-N_narrow: line failures")
fig[3,1] = ax31 = Axis(fig; xlabel="dist(src,dst)", ylabel="N_wide-N_narrow: node failures")

scale = :count
violin!(ax11, x_dist, braessness_lines .+ braessness_nodes; scale=scale, gap=gap, color=:blue)
violin!(ax21, x_dist, braessness_lines; scale=scale, gap=gap, color=:orange)
violin!(ax31, x_dist, braessness_nodes; scale=scale, gap=gap, color=:red)
# ylims!(ax11, ymin, ymax); ylims!(ax21, ymin, ymax); ylims!(ax31, ymin, ymax)
# CairoMakie.save(joinpath(exp_data_dir, "braessness_lines", "braessness_dist_braessness_$model,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.pdf"),fig)
# CairoMakie.save(joinpath(exp_data_dir, "braessness_lines", "braessness_dist_braessness_$model,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.png"),fig)
fig


### 
### Braessness vs Braessness: Scatter without colorscaling (old)
###
fig = Figure(size=(2100,1500), fontsize= fontsize)

# lines and nodes
fig = Figure(size=(2100,1500), fontsize= fontsize)
fig[1,1] = ax11 = Axis(fig; xlabel="inertia failure",ylabel="full failure", title="Braessness node & line failures: f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I", titlealign = :left, titlesize = titlesize)
scatter!(ax11, braessness_lines_BH .+ braessness_nodes_BH, braessness_lines .+ braessness_nodes; color=:blue, markersize=markersize)
fig[2,1] = ax21 = Axis(fig; xlabel="power failure",ylabel="full failure",)
scatter!(ax21, braessness_lines_Pmech .+ braessness_nodes_Pmech, braessness_lines .+ braessness_nodes; color=:blue, markersize=markersize)
fig[1,2] = ax12 = Axis(fig; xlabel="power failure",ylabel="inertia failure",)
scatter!(ax12, braessness_lines_Pmech .+ braessness_nodes_Pmech, braessness_lines_BH .+ braessness_nodes_BH; color=:blue, markersize=markersize)
# CairoMakie.save(joinpath(exp_data_dir, "braessness_lines", "braessness_vs_braessness_lines_and_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.pdf"),fig)
# CairoMakie.save(joinpath(exp_data_dir, "braessness_lines", "braessness_vs_braessness_lines_and_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.png"),fig)
fig

# lines 
fig = Figure(size=(2100,1500), fontsize= fontsize)
fig[1,1] = ax11 = Axis(fig; xlabel="inertia failure",ylabel="full failure", title="Braessness line failures: f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I", titlealign = :left, titlesize = titlesize)
scatter!(ax11, braessness_lines_BH, braessness_lines; color=:orange, markersize=markersize)
fig[2,1] = ax21 = Axis(fig; xlabel="power failure",ylabel="full failure",)
scatter!(ax21, braessness_lines_Pmech, braessness_lines; color=:orange, markersize=markersize)
fig[1,2] = ax12 = Axis(fig; xlabel="power failure",ylabel="inertia failure",)
scatter!(ax12, braessness_lines_Pmech, braessness_lines_BH; color=:orange, markersize=markersize)
# CairoMakie.save(joinpath(exp_data_dir, "braessness_lines", "braessness_vs_braessness_lines,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.pdf"),fig)
# CairoMakie.save(joinpath(exp_data_dir, "braessness_lines", "braessness_vs_braessness_lines,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.png"),fig)
fig

# nodes
fig = Figure(size=(2100,1500), fontsize= fontsize)
fig[1,1] = ax11 = Axis(fig; xlabel="inertia failure",ylabel="full failure", title="Braessness node failures: f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I", titlealign = :left, titlesize = titlesize)
scatter!(ax11, braessness_nodes_BH, braessness_nodes; color=:red, markersize=markersize)
fig[2,1] = ax21 = Axis(fig; xlabel="power failure",ylabel="full failure",)
scatter!(ax21, braessness_nodes_Pmech, braessness_nodes; color=:red, markersize=markersize)
fig[1,2] = ax12 = Axis(fig; xlabel="power failure",ylabel="inertia failure",)
scatter!(ax12, braessness_nodes_Pmech, braessness_nodes_BH; color=:red, markersize=markersize)
# CairoMakie.save(joinpath(exp_data_dir, "braessness_lines", "braessness_vs_braessness_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.pdf"),fig)
# CairoMakie.save(joinpath(exp_data_dir, "braessness_lines", "braessness_vs_braessness_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.png"),fig)
fig

### 
### Braessness vs Braessness: Scatter using colorscaling for counts
###
fontsize = 25
titlesize = (fontsize-5)
markersize = 12

# lines and nodes
xs = braessness_lines_BH .+ braessness_nodes_BH
ys = braessness_lines .+ braessness_nodes
fig = Figure(size=(1000,1500),fontsize=fontsize)
fig[1,1] = ax11 = Axis(fig; xlabel="inertia failure",ylabel="full failure", title="Braessness line & node failures: f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I", titlealign = :left, titlesize = titlesize)
logc = get_log_count(xs, ys)
sc11 = scatter!(ax11, xs, ys; markersize=markersize, color=logc, colormap=:viridis, colorrange=(minimum(logc), maximum(logc)))
Colorbar(fig[1,2], sc11; label = "log₁₀(counts)", width = 30)

xs = braessness_lines_Pmech .+ braessness_nodes_Pmech
ys = braessness_lines .+ braessness_nodes
fig[2,1] = ax21 = Axis(fig; xlabel="power failure",ylabel="full failure",)
logc = get_log_count(xs, ys)
sc21 = scatter!(ax21, xs, ys; markersize=markersize, color=logc, colormap=:viridis, colorrange=(minimum(logc), maximum(logc)))
Colorbar(fig[2,2], sc21; label = "log₁₀(counts)", width = 30)

xs = braessness_lines_Pmech .+ braessness_nodes_Pmech
ys = braessness_lines_BH .+ braessness_nodes_BH
fig[3,1] = ax31 = Axis(fig; xlabel="power failure",ylabel="inertia failure",)
logc = get_log_count(xs, ys)
sc31 = scatter!(ax31, xs, ys; markersize=markersize, color=logc, colormap=:viridis, colorrange=(minimum(logc), maximum(logc)))
Colorbar(fig[3,2], sc31; label = "log₁₀(counts)", width = 30)
CairoMakie.save(joinpath(MA_DIR, "braessness", "exp=$exp_nrs,braessness_vs_braessness_lines_and_nodes_logcount_colorscale,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.pdf"),fig)
CairoMakie.save(joinpath(MA_DIR, "braessness", "exp=$exp_nrs,braessness_vs_braessness_lines_and_nodes_logcount_colorscale,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.png"),fig)
fig

# lines 
xs = braessness_lines_BH
ys = braessness_lines
fig = Figure(size=(1000,1500),fontsize=fontsize)
fig[1,1] = ax11 = Axis(fig; xlabel="inertia failure",ylabel="full failure", title="Braessness line failures: f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I", titlealign = :left, titlesize = titlesize)
logc = get_log_count(xs, ys)
sc11 = scatter!(ax11, xs, ys; markersize=markersize, color=logc, colormap=:viridis, colorrange=(minimum(logc), maximum(logc)))
Colorbar(fig[1,2], sc11; label = "log₁₀(counts)", width = 30)

xs = braessness_lines_Pmech
ys = braessness_lines
fig[2,1] = ax21 = Axis(fig; xlabel="power failure",ylabel="full failure",)
logc = get_log_count(xs, ys)
sc21 = scatter!(ax21, xs, ys; markersize=markersize, color=logc, colormap=:viridis, colorrange=(minimum(logc), maximum(logc)))
Colorbar(fig[2,2], sc21; label = "log₁₀(counts)", width = 30)

xs = braessness_lines_Pmech
ys = braessness_lines_BH
fig[3,1] = ax31 = Axis(fig; xlabel="power failure",ylabel="inertia failure",)
logc = get_log_count(xs, ys)
sc31 = scatter!(ax31, xs, ys; markersize=markersize, color=logc, colormap=:viridis, colorrange=(minimum(logc), maximum(logc)))
Colorbar(fig[3,2], sc31; label = "log₁₀(counts)", width = 30)
CairoMakie.save(joinpath(MA_DIR, "braessness", "exp=$exp_nrs,braessness_vs_braessness_lines_logcount_colorscale,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.pdf"),fig)
CairoMakie.save(joinpath(MA_DIR, "braessness", "exp=$exp_nrs,braessness_vs_braessness_lines_logcount_colorscale,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.png"),fig)
fig

# nodes
xs = braessness_nodes_BH
ys = braessness_nodes
fig = Figure(size=(1000,1500),fontsize=fontsize)
fig[1,1] = ax11 = Axis(fig; xlabel="inertia failure",ylabel="full failure", title="Braessness node failures: f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I", titlealign = :left, titlesize = titlesize)
logc = get_log_count(xs, ys)
sc11 = scatter!(ax11, xs, ys; markersize=markersize, color=logc, colormap=:viridis, colorrange=(minimum(logc), maximum(logc)))
Colorbar(fig[1,2], sc11; label = "log₁₀(counts)", width = 30)

xs = braessness_nodes_Pmech
ys = braessness_nodes
fig[2,1] = ax21 = Axis(fig; xlabel="power failure",ylabel="full failure",)
logc = get_log_count(xs, ys)
sc21 = scatter!(ax21, xs, ys; markersize=markersize, color=logc, colormap=:viridis, colorrange=(minimum(logc), maximum(logc)))
Colorbar(fig[2,2], sc21; label = "log₁₀(counts)", width = 30)

xs = braessness_nodes_Pmech
ys = braessness_nodes_BH
fig[3,1] = ax31 = Axis(fig; xlabel="power failure",ylabel="inertia failure",)
logc = get_log_count(xs, ys)
sc31 = scatter!(ax31, xs, ys; markersize=markersize, color=logc, colormap=:viridis, colorrange=(minimum(logc), maximum(logc)))
Colorbar(fig[3,2], sc31; label = "log₁₀(counts)", width = 30)
CairoMakie.save(joinpath(MA_DIR, "braessness", "exp=$exp_nrs,braessness_vs_braessness_nodes_logcount_colorscale,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.pdf"),fig)
CairoMakie.save(joinpath(MA_DIR, "braessness", "exp=$exp_nrs,braessness_vs_braessness_nodes_logcount_colorscale,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.png"),fig)
fig

###
### (Pmech=-1, BH=-1) has high count
### 
xs = braessness_lines_Pmech .+ braessness_nodes_Pmech
ys = braessness_lines_BH .+ braessness_nodes_BH
zs = braessness_lines.+ braessness_nodes

freqs = Dict{Tuple{Int,Int}, Int}()
for (x,y) in zip(xs, ys)
    freqs[(x,y)] = get(freqs, (x,y), 0) + 1
end

freqs[(0,0)] # 5747
freqs[(1,0)] # 21
freqs[(0,1)] # 8
freqs[(1,1)] # 21
freqs[(1,-1)] # 22
freqs[(-1,-1)] # 156

freqs = Dict{NTuple{3,Int}, Int}()
for (x,y,z) in zip(xs, ys, zs)
    freqs[(x,y,z)] = get(freqs, (x,y,z), 0) + 1
end

freqs[(0,0,0)] # 5731
freqs[(1,0,0)] # 0
freqs[(0,1,0)] # 1
freqs[(1,1,0)] # 0
freqs[(1,-1,0)] # 0
freqs[(-1,-1,0)] # 5

xs_outlier = Int64[]; ys_outlier = Int64[]; zs_outlier = Int64[]
for (x,y,z) in zip(xs, ys, zs)
    if freqs[(x,y,z)] < count
        push!(xs_outlier, x); push!(ys_outlier, y); push!(zs_outlier, z)
    end
end

###
### Looking at the center (small Braessness) 
###
markersize = 20

xs = braessness_lines_Pmech .+ braessness_nodes_Pmech
ys = braessness_lines_BH .+ braessness_nodes_BH
zs = braessness_lines.+ braessness_nodes


xs_wo_center = Int64[]; ys_wo_center = Int64[]; zs_wo_center = Int64[]
for (x,y,z) in zip(xs, ys, zs)
    if x != 0 || y != 0
        push!(xs_wo_center, x); push!(ys_wo_center, y); push!(zs_wo_center, z)
    end
end

xs = xs_wo_center
ys = ys_wo_center
zs = zs_wo_center

## count
fig = Figure()
fig[1,1] = ax11 = Axis(fig; xlabel="inertia failure",ylabel="power failure", title="Braessness line & node failures: f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I", titlealign = :left)
logc = get_log_count(xs, ys)
sc11 = scatter!(ax11, xs, ys; markersize=markersize, color=logc, colormap=:viridis, colorrange=(minimum(logc), maximum(logc)))
Colorbar(fig[1,2], sc11; label = "log₁₀(counts)", width = 30)
lim = 10.5
xlims!(ax11, -lim, lim); ylims!(ax11, -lim, lim)
# CairoMakie.save(joinpath(exp_data_dir, "braessness_lines", "braessness_vs_braessness_lines_and_nodes_center_logcount_colorscale,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.pdf"),fig)
# CairoMakie.save(joinpath(exp_data_dir, "braessness_lines", "braessness_vs_braessness_lines_and_nodes_center_logcount_colorscale,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.png"),fig)
fig