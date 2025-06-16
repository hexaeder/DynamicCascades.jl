using Colors

export plot_simulation

function plot_simulation(sol)
    nw = NetworkDynamics.extract_nw(sol)

    # calculate indices of failing lines and nodes
    idxs_init_swing = map(idx -> idx.compidx, vidxs(nw, :, "ω")) # indices that are initially swing nodes
    all_failing_nodes_idxs = [i for i in idxs_init_swing if sol(sol.t[end], idxs=vidxs(i, :node_swing_stat))[1] == 0]
    all_failing_lines_idxs = [i for i in 1:ne(nw) if sol(sol.t[end], idxs=eidxs(i, :line_stat))[1] == 0]
    node_colors_failing = distinguishable_colors(length(all_failing_nodes_idxs))
    line_colors_failing = distinguishable_colors(length(all_failing_lines_idxs))

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