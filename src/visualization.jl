using Colors

export plot_simulation

function plot_simulation(sol)
    nw = NetworkDynamics.extract_nw(sol)

    # calculate indices of failing lines and nodes
    idxs_init_swing = map(idx -> idx.compidx, vidxs(nw, :, "ω")) # indices that are initially swing nodes
    all_failing_nodes_idxs = [i for i in idxs_init_swing if sol(sol.t[end], idxs=vidxs(i, :node_swing_stat))[1] != 1]
    all_failing_lines_idxs = [i for i in 1:ne(nw) if sol(sol.t[end], idxs=eidxs(i, :line_stat))[1] == 0]
    node_colors_failing = distinguishable_colors(length(all_failing_nodes_idxs))
    line_colors_failing = distinguishable_colors(length(all_failing_lines_idxs))

    # calculate indices of all lines and nodes
    all_nodes_idxs = [i for i in idxs_init_swing]
    all_lines_idxs = [i for i in 1:ne(nw)]
    node_colors = distinguishable_colors(length(all_nodes_idxs))
    line_colors = distinguishable_colors(length(all_lines_idxs))

    ################################################################################
    ############################ Line and nodes ####################################
    ################################################################################
    fontsize = 35
    titlesize = (fontsize+5)
    linewidth = 3.5
    fig = Figure(size=(3100,1500), fontsize=fontsize)
    # Add a global title in the first row spanning all columns.
    # xlim = sol.t[end]/3
    xlim = sol.t[end]

    # CHECK probaly move to plot_simulation(simulation_dir, plot_dir, task_id, initial_fail) using if branch for "WS"
    # and adapt title below in fig[1,1] with ternary operator or so
    # N,k,β,graph_seed,μ,σ,distr_seed,K,α,M,γ,τ,freq_bound,trip_lines,trip_nodes,init_pert,ensemble_element = get_network_args_stripped(df_config, task_id)
    # Maybe even get parameters from `sol`

    # FREQUENCIES ########################################################################
    # fig[1,1] = ax = Axis(fig; ylabel="Frequency [Hz]", title="Lines & nodes that fail (left column), all lines & nodes (right column):  I=$M,D=$γ,τ=$τ,f_b=$freq_bound,α=$α,K=$K,N=$N,k=$k,β=$β,ensemble element=$ensemble_element", titlealign = :left, titlesize = titlesize)
    fig[1,1] = ax = Axis(fig; ylabel="Frequency [Hz]", title="Lines & nodes that fail (left column), all lines & nodes (right column)", titlealign = :left, titlesize = titlesize)
    # NOTE #TODO This works now: `sol(sol.t, idxs=vidxs(1:4, :ω))`
    for i in all_failing_nodes_idxs
        x, y = remove_zero_tail!(deepcopy(sol.t), deepcopy(map(first, sol(sol.t, idxs=vidxs(i, :f)))))
        color = node_colors_failing[findfirst(x -> x == i, all_failing_nodes_idxs)]
        lines!(ax, x, y; label="Node $i", color=color, linewidth=linewidth) # consider plotdensity=Int(1e5)
        scatter!(ax, (x[end], y[end]); color=color, marker=:star5, markersize=25)
    end

    fig[1,2] = ax = Axis(fig; titlealign = :left, titlesize = titlesize)
    for i in all_nodes_idxs
        x, y = remove_zero_tail!(deepcopy(sol.t), deepcopy(map(first, sol(sol.t, idxs=vidxs(i, :f)))))
        color = node_colors[findfirst(x -> x == i, all_nodes_idxs)]
        lines!(ax, x, y; label="Node $i", color=color, linewidth=linewidth) # consider plotdensity=Int(1e5)
        scatter!(ax, (x[end], y[end]); color=color, marker=:star5, markersize=25)
    end
    xlims!(ax, 0, xlim)

    # FLOWS ########################################################################
    fig[2,1] = ax = Axis(fig; xlabel="Time [s]", ylabel="App. power flow [p.u.]")
    for i in all_failing_lines_idxs
        x, y = remove_zero_tail!(deepcopy(sol.t), deepcopy(map(first, sol(sol.t, idxs=eidxs(i, :S)))))
        color = line_colors_failing[findfirst(x -> x == i, all_failing_lines_idxs)]
        lines!(ax, x, y; label="Node $i", color=color, linewidth=linewidth) # consider plotdensity=Int(1e5)
        scatter!(ax, (x[end], y[end]); color=color, marker=:star5, markersize=25)
    end

    fig[2,2] = ax = Axis(fig; xlabel="Time [s]")
    for i in all_lines_idxs
        x, y = remove_zero_tail!(deepcopy(sol.t), deepcopy(map(first, sol(sol.t, idxs=eidxs(i, :S)))))
        color = line_colors[findfirst(x -> x == i, all_lines_idxs)]
        lines!(ax, x, y; label="Node $i", color=color, linewidth=linewidth) # consider plotdensity=Int(1e5)
        scatter!(ax, (x[end], y[end]); color=color, marker=:star5, markersize=25)
    end
    xlims!(ax, 0, xlim)

    return fig
end

function plot_simulation(simulation_dir, plot_dir, task_id, initial_fail)
    sol = Serialization.deserialize(joinpath(simulation_dir, "task_id=$task_id,initial_fail=$initial_fail.sol"));
    fig = plot_simulation(sol)
    # BUG read in parameters here or get them from `sol` 
    CairoMakie.save(joinpath(plot_dir, "ensemble_element=$ensemble_element,I=$M,D=$γ,f_b=$freq_bound,task_id=$task_id,initial_fail=$initial_fail.pdf"),fig)
    CairoMakie.save(joinpath(plot_dir, "ensemble_element=$ensemble_element,I=$M,D=$γ,f_b=$freq_bound,task_id=$task_id,initial_fail=$initial_fail.png"),fig)
end


export remove_zero_tail!
function remove_zero_tail!(x, y)
    @assert length(x) == length(y)
    tail = lastindex(y) + 1
    while y[tail - 1] ≈ 0
        tail -= 1
    end
    if tail <= length(y)
        idxs = tail:lastindex(y)
        deleteat!(x, idxs)
        deleteat!(y, idxs)
    end
    return (x, y)
end

###
### Braessness plots
###

using Graphs
using MetaGraphs
using Unitful
using GraphMakie
using CairoMakie
using DataFrames
using CSV

export write_failures_f_b_to_df, write_network_measures_to_df
export get_braessness, get_network_measures, get_log_count
export plot_braessness_vs_rho_scatter_and_histograms, plot_braessness_vs_dist_histograms, plot_braessness_power_vs_inertia_vs_full_2D
export plot_braessness_vs_braessness_2D_count, plot_braessness_power_vs_inertia_vs_full_3D

## data extraction

"""
Returns for each initially triggered edge the number of line and node failures separately.
"""
function get_line_failures(exp_data_dir, df_config, task_id)
    if "graph_seed" in names(df_config)
        N,k,β,graph_seed,μ,σ,distr_seed,K,α,I,γ,τ,freq_bound,trip_lines,trip_nodes,init_pert,ensemble_element = get_network_args_stripped(df_config, task_id)
        graph_combinations_path = joinpath(exp_data_dir, "k=$k,β=$β")
        failure_mode_string = joinpath(graph_combinations_path, "trip_lines=$trip_lines,trip_nodes=$trip_nodes")
        failure_mode_frequ_bound = joinpath(failure_mode_string, "trip_lines=$trip_lines,trip_nodes=$trip_nodes,freq_bound=$freq_bound")
        filename = "trip_lines=$trip_lines,trip_nodes=$trip_nodes,freq_bound=$freq_bound,N=$N,k=$k,β=$β,graph_seed=$graph_seed,μ=$μ,σ=$σ,distr_seed=$distr_seed,K=$K,α=$α,M=$I,γ=$γ,τ=$τ,init_pert=$init_pert,ensemble_element=$ensemble_element.csv"
        df_result = DataFrame(CSV.File(joinpath(graph_combinations_path,failure_mode_string,failure_mode_frequ_bound,filename)))
    else
        I,γ,τ,freq_bound,trip_lines,trip_nodes,init_pert,ensemble_element = RTS_get_network_args_stripped(df_config, task_id)
        failure_mode_string = joinpath(exp_data_dir, "trip_lines=$trip_lines,trip_nodes=$trip_nodes")
        failure_mode_frequ_bound = joinpath(failure_mode_string, "trip_lines=$trip_lines,trip_nodes=$trip_nodes,freq_bound=$freq_bound")
        filename = "trip_lines=$trip_lines,trip_nodes=$trip_nodes,freq_bound=$freq_bound,M=$I,γ=$γ,τ=$τ,init_pert=$init_pert,ensemble_element=$ensemble_element.csv"
        df_result = DataFrame(CSV.File(joinpath(failure_mode_string,failure_mode_frequ_bound,filename)))
    end
    number_failures_lines = df_result[!, :number_failures_lines]
    number_failures_nodes = df_result[!, :number_failures_nodes]

    return number_failures_lines, number_failures_nodes
end

"""
Get taks_id for given parameters.
"""
function get_task_id(df, f_b, I, ensemble_element)
    for task_id in df.ArrayTaskID
        if I==df[task_id,:inertia_values] && f_b==df[task_id,:freq_bounds] && ensemble_element==df[task_id,:ensemble_element]
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

        # get Pmech at source and dst vertices, `catch` for load nodes with RTS network that don't have `Pmech`
        Pmech_src = try sol(0, idxs=vidxs(src_idx, :Pmech))[1] catch; 0. end
        Pmech_dst = try sol(0, idxs=vidxs(dst_idx, :Pmech))[1] catch; 0. end

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

        # set distance to zero for dead ends
        if dist_e_index > ne(graph)
            dist_e_index = 0
        end

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

    exp_params_dict = Serialization.deserialize(joinpath(exp_data_dir, "exp.params"))
    N_ensemble_size = exp_params_dict[:N_ensemble_size]

    dir = joinpath(exp_data_dir, "braessness_lines", "braessness_data")
    ispath(dir) || mkpath(dir)

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
            line_failures, node_failures = get_line_failures(exp_data_dir, df_config, task_id)

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
        sol = simulate(exp_name_date, Int(task_id), 1;
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


## plotting functions

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

function plot_braessness_vs_rho_scatter_and_histograms(exp_name_date, f_b_narrow, f_b_wide, I;
    fontsize = 25,
    titlesize = (fontsize-5),
    markersize = 8,
    )

    braessness_lines, braessness_nodes = get_braessness(exp_name_date, f_b_narrow, f_b_wide, I)
    x_ρ,_ = get_network_measures(exp_name_date)

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

    # bottom right histograms: counts of Braessness
    μ22 = round(mean(braessness_lines .+ braessness_nodes), digits=3)
    μ23 = round(mean(braessness_lines), digits=3)
    μ24 = round(mean(braessness_nodes), digits=3)
    
    # bins
    all_breaessnesses = vcat(braessness_lines .+ braessness_nodes, braessness_lines, braessness_nodes)
    max_braessness = maximum(all_breaessnesses)
    min_braessness = minimum(all_breaessnesses)
    bins= (min_braessness-0.5):1:(max_braessness+0.5)
    
    hist!(ax22, braessness_lines .+ braessness_nodes; label="<Braessness>=$μ22", color=:blue , bins=bins, direction=:x, strokewidth=strokewidth)
    hist!(ax23, braessness_lines; label="<Braessness>=$μ23", color=:orange, bins=bins, direction=:x, strokewidth=strokewidth)
    hist!(ax24, braessness_nodes; label="<Braessness>=$μ24", color=:red , bins=bins, direction=:x, strokewidth=strokewidth)

    hidexdecorations!(ax11)
    hideydecorations!(ax22); hideydecorations!(ax23); hideydecorations!(ax24)
    axislegend(ax11; position = :rt, labelsize=fontsize)
    axislegend(ax21; position = :rt, labelsize=fontsize-8)
    axislegend(ax22; position = :rt, labelsize=fontsize)
    axislegend(ax23; position = :rt, labelsize=fontsize)
    axislegend(ax24; position = :rt, labelsize=fontsize)
    ax11.yticks = [0, 1, 10, 100, 600]
    ax22.xticks = ax23.xticks = ax24.xticks = [0, 1, 10, 100, 1000, 5000]
    
    return fig
end

function plot_braessness_vs_dist_histograms(exp_name_date, f_b_narrow, f_b_wide, I;
    fontsize = 25
    )

    braessness_lines, braessness_nodes = get_braessness(exp_name_date, f_b_narrow, f_b_wide, I)
    _, x_dist = get_network_measures(exp_name_date)

    fig = Figure(size = (2100, 1700), fontsize = fontsize)
    rows = [1,2,3]
    data_arr = [(braessness_lines .+ braessness_nodes), braessness_lines, braessness_nodes]
    ylabels = ["Braessness lines and nodes", "Braessness lines", "Braessness nodes"]
    colors = [:blue, :orange, :red]

    for (row, data, ylabel, color) in zip(rows, data_arr, ylabels, colors)
        # get categories (here the different distances), exclude 0
        cats = filter(c -> c != 0, sort(unique(x_dist)))

        # Create a 1×length(cats) grid of Axes for all n categories
        axs = [Axis(fig[row, i];
                    xscale = Makie.pseudolog10,
                    xlabel = i == ceil(Int, length(cats)/2) && row == maximum(rows) ? "Count" : "",
                    ylabel = i == 1 ? ylabel : "",
                    title  = row == 1 ? "d=$(cats[i])" : ""
                )
                for (i, c) in enumerate(cats)]

        linkyaxes!(axs...); linkxaxes!(axs...)
        row < maximum(rows) ? hidexdecorations!.(axs[1:end], grid=false) : nothing
        hideydecorations!.(axs[2:end])

        #= compute bin count for each row so all d-panels share the same bin width.
        `0.5` makes sure that bins centrally placed rigt to the y-ticks =#
        const_nbins = (minimum(data)-0.5):1:(maximum(data)+0.5)

        # Draw each histogram in its own Axis
        for (c, ax) in zip(cats, axs)
            # Braessness of the lines with dist = c
            ys = data[findall(x -> x == c, x_dist)]
            μ = round(mean(ys), digits=3)
            # println("c=$c, ys=$ys, bins=$const_nbins") # NOTE use this for checking
            hist!(ax, ys;
                bins      = const_nbins,
                direction = :x,
                label     = "<Braessness>=$μ",
                color     = color,
                strokewidth = 0)
            ax.xticks = [0, 1, 10, 100, 1000]
        end

        [axislegend(i; position = :rt, labelsize=fontsize-6) for i in axs]
    end
    return fig
end

function plot_braessness_power_vs_inertia_vs_full_2D(xs, ys, zs, failures, f_b_narrow, f_b_wide, I;
    fontsize = 25,
    titlesize = (fontsize-5),
    dimensions = 3, # 2D vs 3D selection of scatterpoints
    plot_center_only = false,
    lim = 10.5, # axis limits (only applies for `plot_center_only == true`)
    count = 2, # plot only scatterpoints that occur less often than `count`
    r = 2, # plot only scatterpoints with (|power failure|,|inertia failure|,|full failure|) > r
    fancy_colors = true, # two different color codes for positive and negative Braessness
    markersize = 15
    )   

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

    return fig
end

function plot_braessness_vs_braessness_2D_count(xs, ys, zs, failures, f_b_narrow, f_b_wide, I;
    fontsize = 25,
    titlesize = (fontsize-5),
    markersize = 12,
    plot_coordinates = false
    )

    fig = Figure(size=(1000,1500),fontsize=fontsize)
    fig[1,1] = ax11 = Axis(fig; xlabel="inertia failure",ylabel="full failure", 
        title = plot_coordinates ? "Braessness $failures failures: f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I; (x-axis, y-axis, z-axis)" : "Braessness $failures failures: f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I", 
        titlealign = :left, titlesize = titlesize)
    logc = get_log_count(ys, zs)
    sc11 = scatter!(ax11, ys, zs; markersize=markersize, color=logc, colormap=:viridis, colorrange=(minimum(logc), maximum(logc)))
    plot_coordinates ? add_coordinates(ax11, ys,zs,xs) : nothing
    Colorbar(fig[1,2], sc11; label = "log₁₀(counts)", width = 30)

    fig[2,1] = ax21 = Axis(fig; xlabel="power failure",ylabel="full failure",)
    logc = get_log_count(xs, zs)
    sc21 = scatter!(ax21, xs, zs; markersize=markersize, color=logc, colormap=:viridis, colorrange=(minimum(logc), maximum(logc)))
    plot_coordinates ? add_coordinates(ax21, xs,zs,ys) : nothing
    Colorbar(fig[2,2], sc21; label = "log₁₀(counts)", width = 30)

    fig[3,1] = ax31 = Axis(fig; xlabel="power failure",ylabel="inertia failure",)
    logc = get_log_count(xs, ys)
    sc31 = scatter!(ax31, xs, ys; markersize=markersize, color=logc, colormap=:viridis, colorrange=(minimum(logc), maximum(logc)))
    plot_coordinates ? add_coordinates(ax31, xs,ys,zs) : nothing
    Colorbar(fig[3,2], sc31; label = "log₁₀(counts)", width = 30)
    
    return fig
end

"""
adds coordinates to fig
 - as,bs: determine position on plot
"""
function add_coordinates(ax, as, bs, cs)
    for (a,b,c) in zip(as, bs, cs)
        text!(ax, a, b+1.5, text="($a,$b,$c)";  align = (:center, :top), fontsize=5)
    end
end

function plot_braessness_power_vs_inertia_vs_full_3D(xs, ys, zs, f_b_narrow, f_b_wide, I;
    fontsize = 25,
    titlesize = (fontsize-5),
    markersize = 12
    )

    fig = Figure(size=(1470,1050), fontsize=fontsize)
    fig[1,1] = ax11 = Axis3(fig; xlabel="inertia failure", ylabel="power failure", zlabel="full failure", title="Braessness node & line failures: f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I", titlealign = :left, titlesize = titlesize)
    logc = get_log_count(xs, ys, zs)
    sc11 = scatter!(ax11, xs, ys, zs; markersize=markersize, color=logc, colormap=:viridis, colorrange=(minimum(logc), maximum(logc)))
    Colorbar(fig[1,2], sc11; label = "log₁₀(counts)", width = 30)
    
    return fig
end

export get_task_ids_and_failed_line_from_coordinates
function get_task_ids_and_failed_line_from_coordinates(exp_name_date_x, exp_name_date_y, exp_name_date_z, braessness_x, braessness_y, f_b_narrow, f_b_wide, I;
    run_save_simulation = false,
    generate_plot = false
    )

    braessness_lines_x, braessness_nodes_x = get_braessness(exp_name_date_x, f_b_narrow, f_b_wide, I)
    braessness_lines_y, braessness_nodes_y = get_braessness(exp_name_date_y, f_b_narrow, f_b_wide, I)
    xs = braessness_lines_x .+ braessness_nodes_x
    ys = braessness_lines_y .+ braessness_nodes_y

    idx = only([i for i in eachindex(ys) if xs[i] == braessness_x && ys[i] == braessness_y])
    
    exp_data_dir_x = joinpath(RESULTS_DIR, exp_name_date_x)
    exp_data_dir_y = joinpath(RESULTS_DIR, exp_name_date_y)
    exp_data_dir_z = joinpath(RESULTS_DIR, exp_name_date_z)

    df_narrow_x = DataFrame(CSV.File(joinpath(exp_data_dir_x,"braessness_lines", "braessness_data", "I=$I,f_b=$f_b_narrow.csv")))
    df_wide_x = DataFrame(CSV.File(joinpath(exp_data_dir_x,"braessness_lines", "braessness_data", "I=$I,f_b=$f_b_wide.csv")))
    df_narrow_y = DataFrame(CSV.File(joinpath(exp_data_dir_y,"braessness_lines", "braessness_data", "I=$I,f_b=$f_b_narrow.csv")))
    df_wide_y = DataFrame(CSV.File(joinpath(exp_data_dir_y,"braessness_lines", "braessness_data", "I=$I,f_b=$f_b_wide.csv")))
    df_narrow_z = DataFrame(CSV.File(joinpath(exp_data_dir_z,"braessness_lines", "braessness_data", "I=$I,f_b=$f_b_narrow.csv")))
    df_wide_z = DataFrame(CSV.File(joinpath(exp_data_dir_z,"braessness_lines", "braessness_data", "I=$I,f_b=$f_b_wide.csv")))

    # check
    braessness_line_idx_x = df_wide_x[idx, :].number_failures_lines .- df_narrow_x[idx, :].number_failures_lines
    braessness_node_idx_x = df_wide_x[idx, :].number_failures_nodes .- df_narrow_x[idx, :].number_failures_nodes
    @assert braessness_line_idx_x + braessness_node_idx_x == braessness_x
    braessness_line_idx_y = df_wide_y[idx, :].number_failures_lines .- df_narrow_y[idx, :].number_failures_lines
    braessness_node_idx_y = df_wide_y[idx, :].number_failures_nodes .- df_narrow_y[idx, :].number_failures_nodes
    @assert braessness_line_idx_y + braessness_node_idx_y == braessness_y
    braessness_line_idx_z = df_wide_z[idx, :].number_failures_lines .- df_narrow_z[idx, :].number_failures_lines
    braessness_node_idx_z = df_wide_z[idx, :].number_failures_nodes .- df_narrow_z[idx, :].number_failures_nodes

    exp_name_date_arr = []
    task_ids = []
    initial_fail = []

    push!(exp_name_date_arr, exp_name_date_x)
    push!(task_ids, df_wide_x[idx, :].ArrayTaskID)
    push!(initial_fail, df_wide_x[idx, :].initially_failed_line)
    push!(exp_name_date_arr, exp_name_date_x)
    push!(task_ids, df_narrow_x[idx, :].ArrayTaskID)
    push!(initial_fail, df_narrow_x[idx, :].initially_failed_line)

    push!(exp_name_date_arr, exp_name_date_y)
    push!(task_ids, df_wide_y[idx, :].ArrayTaskID)
    push!(initial_fail, df_wide_y[idx, :].initially_failed_line)
    push!(exp_name_date_arr, exp_name_date_y)
    push!(task_ids, df_narrow_y[idx, :].ArrayTaskID)
    push!(initial_fail, df_narrow_y[idx, :].initially_failed_line)

    push!(exp_name_date_arr, exp_name_date_z)
    push!(task_ids, df_wide_z[idx, :].ArrayTaskID)
    push!(initial_fail, df_wide_z[idx, :].initially_failed_line)
    push!(exp_name_date_arr, exp_name_date_z)
    push!(task_ids, df_narrow_z[idx, :].ArrayTaskID)
    push!(initial_fail, df_narrow_z[idx, :].initially_failed_line)

    println("x-Axis: narrow:  ArrayTaskID=$(df_narrow_x[idx, :].ArrayTaskID), initially_failed_line=$(df_narrow_x[idx, :].initially_failed_line)
         wide:   ArrayTaskID=$(df_wide_x[idx, :].ArrayTaskID), initially_failed_line=$(df_wide_x[idx, :].initially_failed_line)")
    println("y-Axis: narrow:  ArrayTaskID=$(df_narrow_y[idx, :].ArrayTaskID), initially_failed_line=$(df_narrow_y[idx, :].initially_failed_line)
         wide:   ArrayTaskID=$(df_wide_y[idx, :].ArrayTaskID), initially_failed_line=$(df_wide_y[idx, :].initially_failed_line)")
    println("z-Axis: narrow:  ArrayTaskID=$(df_narrow_z[idx, :].ArrayTaskID), initially_failed_line=$(df_narrow_z[idx, :].initially_failed_line)
         wide:   ArrayTaskID=$(df_wide_z[idx, :].ArrayTaskID), initially_failed_line=$(df_wide_z[idx, :].initially_failed_line), Braessness=$(braessness_line_idx_z + braessness_node_idx_z)")

    if run_save_simulation == true
        for (exp_name_date, task_id, initial_fail) in zip(exp_name_date_arr, task_ids, initial_fail)
            println("dir=$exp_name_date,task_id=$task_id, initial_fail=$initial_fail")

            #= Assign node failure model
            TODO For future simulations remove this block. (In WS simulations I did not assign `:gen_model`,
            where initially only one node failure model was simulated).=#
            if exp_name_date[11:12] == "04"
                gen_model = SwingDynLoadModel
            elseif exp_name_date[11:12] == "11"
                gen_model = SwingDynLoadModel_change_to_BH_only
            elseif exp_name_date[11:12] == "12"
                gen_model = SwingDynLoadModel_change_Pmech_only
            end
            # TODO For new WS simulations replace by `gen_model = Serialization.deserialize(joinpath(RESULTS_DIR, exp_name_date, "exp.params"))[:gen_model]`
            gen_model = exp_name_date[1:7] == "RTS_exp" ? Serialization.deserialize(joinpath(RESULTS_DIR, exp_name_date, "exp.params"))[:gen_model] : gen_model

            sol = simulate(exp_name_date, task_id, initial_fail;
                gen_model = gen_model,
                initial_fail = initial_fail,
                verbose = true);
            
            exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
            ispath(joinpath(exp_data_dir, "example_trajectories", "sims")) || mkpath(joinpath(exp_data_dir, "example_trajectories", "sims"))
            ispath(joinpath(exp_data_dir, "example_trajectories", "plots")) || mkpath(joinpath(exp_data_dir, "example_trajectories", "plots"))
            Serialization.serialize(joinpath(exp_data_dir, "example_trajectories", "sims", "gen_model=$gen_model,task_id=$task_id,initial_fail=$initial_fail.sol"), sol)
            if generate_plot == true
                CairoMakie.activate!()
                save(joinpath(exp_data_dir, "example_trajectories", "plots", "gen_model=$gen_model,task_id=$task_id,initial_fail=$initial_fail.pdf"), plot_simulation(sol))
            end

            # generate_plot ? save(joinpath(exp_data_dir, "example_trajectories", "plots", "gen_model=$gen_model,task_id=$task_id,initial_fail=$initial_fail.pdf"), plot_simulation(sol)) : nothing
        end
    end

end


###
### wrappers for `inspect()`
###
using NetworkDynamicsInspector
export inspect_wrapper

"""
`which_trajectories`: `:all`, `:failing`
"""
function inspect_wrapper(sol;
    which_trajectories = :failing,
    tmax = Float64[]
    )

    nw = NetworkDynamics.extract_nw(sol)

    idxs_init_swing = map(idx -> idx.compidx, vidxs(nw, :, "ω")) # indices that are initially swing nodes
    if which_trajectories == :failing
        vindices = [i for i in idxs_init_swing if sol(sol.t[end], idxs=vidxs(i, :node_swing_stat))[1] != 1]
        eindices = [i for i in 1:ne(nw) if sol(sol.t[end], idxs=eidxs(i, :line_stat))[1] == 0]
        # set marker for failed vertices
        for i in vindices
            set_marker!(nw[VIndex(i)], :xcross)
        end
    elseif which_trajectories == :all
        vindices = [i for i in idxs_init_swing]
        eindices = [i for i in 1:ne(nw)]
    end

    inspect(sol; restart=true, reset=true)
    set_graphplot!(; nstate=[:ω], estate=[:S], nstate_rel=false, estate_rel=false)
    if !isempty(tmax)
        set_state!(; t=0.0, tmin=0.0, tmax=tmax)
    end
    define_timeseries!([
        (; selcomp=[VIndex(i) for i in vindices], states=[:ω, :ωmax], rel=false),
        (; selcomp=[EIndex(i) for i in eindices], states=[:S, :rating], rel=false),
    ])
end

function inspect_wrapper(dir, task_id, initial_fail; kwargs...)
    sol = Serialization.deserialize(joinpath(dir, "task_id=$task_id,initial_fail=$initial_fail.sol"));
    inspect_wrapper(sol; kwargs...)
end


