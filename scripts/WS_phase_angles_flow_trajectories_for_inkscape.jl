"""
Snapshots of illustrative example network with narrow and wide frequency bounds. Simplified version.
"""

using Revise

include(abspath(@__DIR__, "cluster/experiment_jarray", "helpers_jarray.jl"))

using DynamicCascades
using Graphs
using MetaGraphs
using Statistics
using GraphMakie
using Colors
using DynamicCascades: PLOT_DIR
using NetworkDynamics
using Unitful
using Unitful: @u_str
using DataFrames
using CSV
using LinearAlgebra
using CairoMakie
CairoMakie.activate!()

using Serialization
using ColorSchemes
using NetworkLayout: spring


function simulate_wrapper(exp_name_date, task_id, initial_fail;
                            save_simulation = false,
                            tweak_power_injections = false)

    exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
    df_config = DataFrame(CSV.File(joinpath(exp_data_dir, "config.csv")))
    # adjust filepaths 
    df_config[!, :filepath_graph] = replace.(df_config[!, :filepath_graph],"/home/brandner" => "/home/brandner/nb_data/HU_Master/2122WS/MA")
    df_config[!, :filepath_steady_state] = replace.(df_config[!, :filepath_steady_state],"/home/brandner" => "/home/brandner/nb_data/HU_Master/2122WS/MA")
    exp_params_dict = Serialization.deserialize(joinpath(exp_data_dir, "exp.params"))

    N,k,β,graph_seed,μ,σ,distr_seed,K,α,M,γ,τ,freq_bound,trip_lines,trip_nodes,init_pert,ensemble_element = get_network_args_stripped(df_config, task_id)

    #= Generate sol-Objects, the data for the plots. Only use this for recreating
    the solution objects =#
    network = import_system_wrapper(df_config, task_id)

    # read in steady state
    steady_state_dict  = CSV.File(df_config[task_id,:filepath_steady_state])
    x_static = steady_state_dict[:SteadyState]

    ## Tweak power injections
    if tweak_power_injections == true
        # Adapt v5
        P_v5 = get_prop(network, 5, :P)
        P_load_v5 = get_prop(network, 5, :P_load)
        P_inj_v5 =  P_v5 + P_load_v5
        P_inj_v5_new = P_inj_v5 - 2*P_inj_v5
        P_load_v5_new = P_load_v5 - P_inj_v5
        set_prop!(network, 5, :P_load, P_load_v5_new)
        set_prop!(network, 5, :P, P_inj_v5_new - P_load_v5_new)

        # Restore global power balance
        v_restore_id = 35
        P_v_restore = get_prop(network, v_restore_id, :P)
        P_load_v_restore = get_prop(network, v_restore_id, :P_load)
        P_inj_v_restore = P_v_restore + P_load_v_restore
        P_inj_v_restore_new = P_inj_v_restore + 2*P_inj_v5
        P_load_v_restore_new = P_load_v_restore + P_inj_v5
        set_prop!(network, v_restore_id, :P_load, P_load_v_restore_new)
        set_prop!(network, v_restore_id, :P, P_inj_v_restore_new - P_load_v_restore_new)

        x_static = steadystate(network; tol=1e-7, zeroidx=1) # CHECK 
    end

    # network is balanced
    sum([get_prop(network, i, :P) for i in 1:nv(network)]) # 3.3306690738754696e-15

    if exp_name_date[11:12] == "04"
        node_failure_model = :change_to_BH_and_change_Pmech
    elseif exp_name_date[11:12] == "11"
        node_failure_model = :change_to_BH
    elseif exp_name_date[11:12] == "12"
        node_failure_model = :change_Pmech
    end

    sol = simulate(network;
                x_static=x_static,
                initial_fail = [initial_fail],
                node_failure_model = node_failure_model,
                f_min = -freq_bound,
                f_max = freq_bound,
                solverargs = (;reltol=1e-8, abstol=1e-6),
                verbose = true);
    dir = joinpath(exp_data_dir, "trajectories")
    ispath(dir) || mkpath(dir)
    save_simulation ? Serialization.serialize(joinpath(exp_data_dir, "trajectories", "task_id=$task_id,initial_fail=$initial_fail,tweak_power_injections=$tweak_power_injections.sol"), sol) : nothing

    return sol
end

###
### collect data
###
exp_name_date_full = "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250321_171511.976"

initial_fail = 15
task_id_full_n = 7584
task_id_full_w = 7755

## string for saving plot
exp_nr_full = exp_name_date_full[11:12]
prefix = exp_name_date_full[1:10]

df_config = DataFrame(CSV.File(joinpath(RESULTS_DIR, exp_name_date_full, "config.csv")));


tweak_power_injections = true
# Network from ensemble
# sol1 = simulate_wrapper(exp_name_date_full, task_id_full_n, initial_fail; save_simulation = true, tweak_power_injections=tweak_power_injections); # narrow
# sol4 = simulate_wrapper(exp_name_date_full, task_id_full_w, initial_fail; save_simulation = true, tweak_power_injections=tweak_power_injections); # wide

# load .sol-files
sol1 = Serialization.deserialize(joinpath(RESULTS_DIR, exp_name_date_full, "trajectories", "task_id=$task_id_full_n,initial_fail=$initial_fail,tweak_power_injections=$tweak_power_injections.sol"));
sol4 = Serialization.deserialize(joinpath(RESULTS_DIR, exp_name_date_full, "trajectories", "task_id=$task_id_full_w,initial_fail=$initial_fail,tweak_power_injections=$tweak_power_injections.sol"));



###
### global plotting parameters
###
x_max = 14.1

cols = distinguishable_colors(length(selected_lines) + 3)
line_colors = cols[[8, 3, 4, 5, 7]]   # chosen directly, no deleteat!
# hex.(line_colors)
#  "000000"
#  "FF9BFF"
#  "00D3FF"
#  "E2630D"
#  "007E00"
#  "0050E6"

# new
#  "AC0047"
#  "FF9BFF"
#  "00D3FF"
#  "E2630D"
#  "0050E6"





###
### Trajectories phase angles
###
selected_nodes = [4, 5, 6, 82]

# -------------------- GLOBAL (shared) axis style --------------------
const AXSTYLE = (;
    xlabel        = "Time (s)",
    ylabel        = "Phase θ (rad)",
    xticklabelsize = 45,
    yticklabelsize = 45,
    xlabelsize     = 45,
    ylabelsize     = 45,
    titlesize      = 25,
)


# -------------------- One figure per node --------------------
function node_theta_figure(sol_narrow::SolutionContainer, sol_wide::SolutionContainer, nodeid::Int;
        resolution = (900, 600),
        # individually adjustable:
        xlim = (-0.1, 20.0),
        ylim = (-1.0, 2.0),
        yticks = ylim[1]:0.2:ylim[2], 
        vlines_narrow = sol_narrow.failures_nodes.t[1],
        vlines_wide   = sol_wide.failures_nodes.t,
        show_vlines_narrow = true,
        show_vlines_wide   = false,
        lw_curve = 3,
        lw_vline = 2,
        save_fig = true,
    )

    fig = Figure(resolution = resolution, figure_padding = 30)
    ax  = Axis(fig[1,1]; AXSTYLE...)
    # Box(fig[1,1]; color=:transparent, strokecolor=:black, strokewidth=4)

    # place title inside plot area
    node_label_map = Dict(82 => 1, 5 => 2, 4 => 3, 6 => 4)
    node_label_id = node_label_map[nodeid]
    title = "Node $node_label_id"
    text!(ax, title;
    position = (0.04, 0.97),      # ← vertical position inside box
    align = (:left, :top),
    space = :relative,
    textsize = 40,
    color = :black)


    (nd, _,_) = nd_model(sol_narrow.network)
    state_idx = idx_containing(nd, "θ")
    node_idx_all = map(s -> parse(Int, String(s)[4:end]), nd.syms[state_idx])

    k = findfirst(x -> x == nodeid, node_idx_all)
    θidx = state_idx[k]

    # narrow (solid)
    t, s = seriesforidx(sol_narrow.sol, θidx)
    lines!(ax, t, s; color=:black, linewidth=lw_curve)

    # wide (dashed)
    t, s = seriesforidx(sol_wide.sol, θidx)
    lines!(ax, t, s; color=:black, linestyle=:dash, linewidth=lw_curve)

    # star for line failures
    if (nodeid == 4 || nodeid == 6) # line 13 connects node 4 and 6
        t_fail_e13 = sol_wide.failures.t[3]
        θidx_node_adjacent_to_e13 = state_idx[findfirst(x -> x == nodeid, node_idx_all)]
        scatter!(ax, (t_fail_e13, sol_wide.sol(t_fail_e13)[θidx_node_adjacent_to_e13]); color=line_colors[3], marker=:star5, markersize=35)
    end
    if nodeid == 5 # line 13 connects node 4 and 6
        t_fail_e9 = sol_narrow.failures.t[2]
        θidx_node_adjacent_to_e9 = state_idx[findfirst(x -> x == nodeid, node_idx_all)]
        scatter!(ax, (t_fail_e9, sol_wide.sol(t_fail_e9)[θidx_node_adjacent_to_e9]); color=line_colors[1], marker=:star5, markersize=35)
    end


    # vertical lines for node failures
    show_vlines_narrow && vlines!(ax, vlines_narrow; color=:black, linewidth=lw_vline, linestyle=:dash)
    show_vlines_wide   && vlines!(ax, vlines_wide;   color=:black, linewidth=lw_vline, linestyle=:dash)
    if nodeid == 5
        text!(ax, sol1.failures_nodes.t[1]+0.3, 0.6; text="← node $nodeid fails \n     narrow bounds", textsize=35)
    end


    xlims!(ax, xlim...)
    ylims!(ax, ylim...)
    # ax.xticks = xticks
    ax.yticks = yticks

    if save_fig == true 
        save(joinpath(MA_DIR, string("WS/illustrative_plot/node_$nodeid.svg")), fig)
        save(joinpath(MA_DIR, string("WS/illustrative_plot/node_$nodeid.pdf")), fig)
        save(joinpath(MA_DIR, string("WS/illustrative_plot/node_$nodeid.png")), fig)
    end
    
    return fig
end

y_axis_length = 1.0
x_lowlim4 = 0.7; x_lowlim5 = 0.4; x_lowlim6 = 0.0; x_lowlim82 = 0.3;
fig4  = node_theta_figure(sol1, sol4, 4;  xlim=(0.1,x_max), ylim=(x_lowlim4, x_lowlim4 + y_axis_length))
fig5  = node_theta_figure(sol1, sol4, 5;  xlim=(0.1,x_max), ylim=(x_lowlim5, x_lowlim5 + y_axis_length))
fig6  = node_theta_figure(sol1, sol4, 6;  xlim=(0.1,x_max), ylim=(x_lowlim6, x_lowlim6 + y_axis_length))
fig82 = node_theta_figure(sol1, sol4, 82; xlim=(0.1,x_max), ylim=(x_lowlim82, x_lowlim82 + y_axis_length))



###
### Trajectories power flow
###
selected_lines = [9, 12, 13, 14, 16]

fig = Figure(resolution = (1800, 1200), figure_padding = 30)
ax = Axis(fig[1,1];
    xlabel = "Time (s)",
    ylabel = "Apparent power flow (p.u.)",
    title = "Full failure model: narrow frequency bounds (solid) wide frequency bounds (dashed)",
    titlealign = :left,
    xticklabelsize = 45,
    yticklabelsize = 45,
    xlabelsize = 45,
    ylabelsize = 45,
    titlesize = 25,
)

# full failures narrow
for (i, l) in pairs(selected_lines)
    t, s = seriesforidx(sol1.load_S, l)
    lines!(ax, t, s; label="line $l",color=line_colors[i], linewidth=5)
end
for (i, l) in pairs(selected_lines)
    t, s = seriesforidx(sol1.load_S, l)
    scatter!(ax, (t[end], s[end]); marker=:star5,color=line_colors[i],  markersize=35)
end

# full failure wide
for (i, l) in pairs(selected_lines)
    t, s = seriesforidx(sol4.load_S, l)
    lines!(ax, t, s; linestyle=:dash, color=line_colors[i], linewidth=5)
end
for (i, l) in pairs(selected_lines)
    t, s = seriesforidx(sol4.load_S, l)
    scatter!(ax, (t[end], s[end]); marker=:star5,color=line_colors[i], markersize=35)
end

rating = get_prop(sol1.network, Graphs.edges(sol1.network), :rating)[1]
hlines!(ax, rating; color=:black, linestyle=:dot, linewidth=5, label="line rating")
# vlines!(ax, times; color=:black, linestyle=:dot, linewidth=1)

# ax.xticks = xticks_ax_gc21
ax.xticksvisible = true

# xlims!(ax, -0.1, x_max)
xlims!(ax, 0.1, x_max)
# Legend(fig[1,2], ax, labelsize=35) # plots show the same lines in all plots

fig
save(joinpath(MA_DIR, string("WS/illustrative_plot/apparent_power_flow.svg")), fig)
save(joinpath(MA_DIR, string("WS/illustrative_plot/apparent_power_flow.pdf")), fig)
save(joinpath(MA_DIR, string("WS/illustrative_plot/apparent_power_flow.png")), fig)



###
### Plot legend separtely
###
fig = Figure(resolution=(275, 270), figure_padding=30)

# put axis in a tiny cell so nothing meaningful is visible
ax = Axis(fig[1,1], width=0, height=0)
hidedecorations!(ax); hidespines!(ax)

# lines in network 16, 9, 12, 14, 13
labels = ["line 1↔2", "line 2↔network", "line 2↔3", "line 2↔4", "line 3↔4"]
# switch order of line colors
line_colors = [line_colors[5], line_colors[1], line_colors[2], line_colors[4], line_colors[3]]
for (i, l) in pairs(selected_lines)
    lines!(ax, [0, 1], [i, i]; color=line_colors[i], linewidth=5, label=labels[i])
end
lines!(ax, [0, 1], [0, 0]; color=:black, linewidth=5, linestyle=:dot, label="line rating")

# legend as the ONLY visible thing
Legend(fig[1,1], ax; labelsize=35, framevisible=false)

save(joinpath(MA_DIR, string("WS/illustrative_plot/legend.svg")), fig)
save(joinpath(MA_DIR, string("WS/illustrative_plot/legend.pdf")), fig)
save(joinpath(MA_DIR, string("WS/illustrative_plot/legend.png")), fig)


###
### save intermediate files
###
datetime = Dates.format(now(), "_yyyymmdd_HHMMSS.s")
# save different versions of script
script_dest = joinpath(MA_DIR, string("WS/illustrative_plot/old_versions/WS_phase_angles_flow_trajectories_for_inkscape", datetime, ".jl"))
cp(@__FILE__, script_dest; force=true);