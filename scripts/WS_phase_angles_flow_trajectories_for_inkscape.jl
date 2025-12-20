"""
Illustrative example narrow vs wide failures frequency bounds. This script generates phase plots for 4 node and another
with the flow trajectories of adjacent lines. The plots are then used for a figure in inkscape that adds a cut out of the network.
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

        x_static = steadystate(network; tol=1e-7, zeroidx=1)
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
                solverargs = (;reltol=1e-8, abstol=1e-6, dtmax=0.01),
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

## fig size
# DIN A4 size in mm: 210mm x 297mm
# Natur Comms: 210mm x 276 mm
mm_to_pt(mm) = mm * 72 / 25.4
x_max = 14.1
fontsize = 9
lw_vline = 0.5
lw_flow_trajectories = 2.5
lw_phase = 1
starsize = 12
selected_lines = [9, 12, 13, 14, 16]
backgroundcolor = RGBf(242/255, 245/255, 250/255) # figure background #F2F5FA
# backgroundcolor = RGBf(235/255, 239/255, 247/255) # slightly darker
# backgroundcolor = RGBf(222/255, 228/255, 239/255) # darker
# hex(backgroundcolor) # DEE4EF

cols = distinguishable_colors(length(selected_lines) + 3)
line_colors = cols[[8, 3, 4, 5, 7]]
# hex.(line_colors)
#  "AC0047"
#  "FF9BFF"
#  "00D3FF"
#  "E2630D"
#  "0050E6"

function adjust_axis(ax)
    ax.xticksize  = ax.yticksize = 2
    ax.xtickwidth = ax.ytickwidth = 0.5

    # ticks point into the axis
    ax.xtickalign = 1
    ax.ytickalign = 1

    ax.spinewidth = 0.4
    ax.xgridvisible = true
    ax.ygridvisible = true
    ax.xgridwidth   = 0.3
    ax.ygridwidth   = 0.3
    # ax.xgridcolor   = RGBAf(0,0,0,0.08) # TODO also in Blue-grey
    # ax.ygridcolor   = RGBAf(0,0,0,0.08)

    # hide grid
    ax.xgridvisible = true
    ax.ygridvisible = true

    return ax    
end

###
### Trajectories phase angles
###
function node_theta_figure(sol_narrow::SolutionContainer, sol_wide::SolutionContainer, nodeid::Int;
        xlim = (-0.1, 20.0),
        ylim = (-1.0, 2.0),
        save_fig = true)

    # map node ids to labels 1-4
    node_label_map = Dict(82 => 1, 5 => 2, 4 => 3, 6 => 4)
    node_label_id = node_label_map[nodeid]

    # figure size
    width_pt  = mm_to_pt(62.5*0.9)
    height_pt = mm_to_pt(41.5*0.9)


    fig = Figure(resolution = (width_pt, height_pt), fontsize=fontsize, figure_padding = 8, backgroundcolor = backgroundcolor)
    ax  = Axis(fig[1,1]; xlabel="Time (s)", ylabel="Phase θ (rad)", backgroundcolor = backgroundcolor)


    xlims!(ax, xlim...)
    ylims!(ax, ylim...)
    adjust_axis(ax)

    # show only x- and y-axis (left & bottom)
    hidespines!(ax, :t, :r)

    # move tick labels inside
    ax.xticklabelpad = -14 # digits x-ticks
    ax.ylabelpadding = -14  # move y-label a bit further right
    ax.xlabelvisible = false
    if node_label_id == 2
        text!(ax, "Time (s)"; position = (0.52, 0.155), align = (:left, :top), space = :relative, textsize = fontsize, color = :black)
    end


    # x-ticks: show only first and last ticklabel
    ax.xticks = ([2.5, 5.0, 7.5, 10.0, 12.5], ["", "", "", "", "12.5"])

    # y-ticks: show only first and last ticklabel
    ys = collect(ylim[1]:0.2:ylim[2])
    labels = fill("", length(ys))
    labels[1] = string(round(ys[1], digits=2))
    labels[end] = string(round(ys[end], digits=2))
    ax.yticks = (ys, labels)

    # hide decorations for nodes 2/3/4
    if node_label_id != 2
        ax.xlabelvisible = false
        ax.ylabelvisible = false
        ax.xticklabelsvisible = false
    end


    # place title inside plot area
    text!(ax, "Node $node_label_id"; position = (0.13, 1.00), align = (:left, :top), space = :relative, textsize = fontsize, color = :black)

    (nd, _,_) = nd_model(sol_narrow.network)
    state_idx = idx_containing(nd, "θ")
    node_idx_all = map(s -> parse(Int, String(s)[4:end]), nd.syms[state_idx])

    k = findfirst(x -> x == nodeid, node_idx_all)
    θidx = state_idx[k]

    # narrow (solid)
    t, s = seriesforidx(sol_narrow.sol, θidx)
    lines!(ax, t, s; color=:black, linewidth=lw_phase)

    # wide (dashed)
    t, s = seriesforidx(sol_wide.sol, θidx)
    lines!(ax, t, s; color=:black, linestyle=:dash, linewidth=lw_phase)

    # star for line failures
    if (nodeid == 4 || nodeid == 6) # line 13 connects node 4 and 6
        t_fail_e13 = sol_wide.failures.t[3]
        θidx_node_adjacent_to_e13 = state_idx[findfirst(x -> x == nodeid, node_idx_all)]
        scatter!(ax, (t_fail_e13, sol_wide.sol(t_fail_e13)[θidx_node_adjacent_to_e13]); color=line_colors[3], marker=:star5, markersize=starsize)
    end
    if nodeid == 5 # line 13 connects node 4 and 6
        t_fail_e9 = sol_narrow.failures.t[2]
        θidx_node_adjacent_to_e9 = state_idx[findfirst(x -> x == nodeid, node_idx_all)]
        scatter!(ax, (t_fail_e9, sol_wide.sol(t_fail_e9)[θidx_node_adjacent_to_e9]); color=line_colors[1], marker=:star5, markersize=starsize)
        text!(ax, sol1.failures_nodes.t[1]+0.35, 0.74; text="← node $node_label_id fails for\n     narrow bounds", textsize=fontsize-1)
        # text!(ax, 8.0, 1.28; text="wide ↓", textsize=fontsize)
        # text!(ax, 9.7, 1.07; text="↑ narrow", textsize=fontsize)
    end

    # vertical lines for node failures
    vlines!(ax, sol_narrow.failures_nodes.t[1]; color=:black, linewidth=lw_vline, linestyle=:dash)

    if save_fig == true
        save(joinpath(MA_DIR, "WS/illustrative_plot/node_$nodeid.png"), fig)
        save(joinpath(MA_DIR, "WS/illustrative_plot/node_$nodeid.svg"), fig; pt_per_unit=1)
        save(joinpath(MA_DIR, "WS/illustrative_plot/node_$nodeid.pdf"), fig; pt_per_unit=1)
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

width_mm  = 72.5*1.8 # 62.5*2
height_mm = 41.5*1.8
fig = Figure(resolution = (mm_to_pt(width_mm), mm_to_pt(height_mm)), fontsize=fontsize, figure_padding = 0)

ax = Axis(fig[1,1]; xlabel = "Time (s)", ylabel = "Apparent power flow (p.u.)")
    # title = "Full failure model: narrow frequency bounds (solid) wide frequency bounds (dashed)",

# full failures narrow
for (i, l) in pairs(selected_lines)
    t, s = seriesforidx(sol1.load_S, l)
    lines!(ax, t, s; label="line $l",color=line_colors[i], linewidth=lw_flow_trajectories)
end
for (i, l) in pairs(selected_lines)
    t, s = seriesforidx(sol1.load_S, l)
    scatter!(ax, (t[end], s[end]); marker=:star5,color=line_colors[i],  markersize=starsize)
end

# full failure wide
for (i, l) in pairs(selected_lines)
    t, s = seriesforidx(sol4.load_S, l)
    lines!(ax, t, s; linestyle=:dash, color=line_colors[i], linewidth=lw_flow_trajectories)
end
for (i, l) in pairs(selected_lines)
    t, s = seriesforidx(sol4.load_S, l)
    scatter!(ax, (t[end], s[end]); marker=:star5,color=line_colors[i], markersize=starsize)
end

rating = get_prop(sol1.network, Graphs.edges(sol1.network), :rating)[1]
hlines!(ax, rating; color=:black, linestyle=:dot, linewidth=lw_flow_trajectories-1, label="line rating")
vlines!(ax, sol1.failures_nodes.t[1]; color=:black, linewidth=lw_vline, linestyle=:dash)
# text!(ax, 9.4, 0.68; text="wide ↓", textsize=fontsize)
# text!(ax, 10.1, 0.25; text="↑ narrow", textsize=fontsize)

xlims!(ax, 0.1, x_max)
adjust_axis(ax)

save(joinpath(MA_DIR, "WS/illustrative_plot/apparent_power_flow.png"), fig)
save(joinpath(MA_DIR, "WS/illustrative_plot/apparent_power_flow.pdf"), fig; pt_per_unit=1)
save(joinpath(MA_DIR, "WS/illustrative_plot/apparent_power_flow.svg"), fig; pt_per_unit=1)


###
### Plot legends separately
###
patchsize = (24,9)

##
## Legend flow trajectories
##
width_mm  = 35
height_mm = 30
fig = Figure(resolution = (mm_to_pt(width_mm), mm_to_pt(height_mm)), fontsize=fontsize)

# put axis in a tiny cell so nothing meaningful is visible
ax = Axis(fig[1,1], width=0, height=0)
hidedecorations!(ax); hidespines!(ax)

# lines in network 16, 9, 12, 14, 13
labels = ["line 1↔2", "line 2↔network", "line 2↔3", "line 2↔4", "line 3↔4"]
# switch order of line colors
line_colors = [line_colors[5], line_colors[1], line_colors[2], line_colors[4], line_colors[3]]
# lines!(ax, [0, 1], [0, 0]; color=:black, linewidth=lw_flow_trajectories, label="narrow bounds")
# lines!(ax, [0, 1], [0, 0]; color=:black, linewidth=lw_flow_trajectories, linestyle=:dash, label="wide bounds")
for (i, l) in pairs(selected_lines)
    lines!(ax, [0, 1], [i, i]; color=line_colors[i], linewidth=lw_flow_trajectories, label=labels[i])
end
lines!(ax, [0, 1], [0, 0]; color=:black, linewidth=lw_flow_trajectories-1, linestyle=:dot, label="line rating")

# legend as the ONLY visible thing
Legend(fig[1,1], ax; labelsize=fontsize, framevisible=false, patchsize=patchsize, orientation = :vertical)

save(joinpath(MA_DIR, "WS/illustrative_plot/legend.png"), fig)
save(joinpath(MA_DIR, "WS/illustrative_plot/legend.pdf"), fig; pt_per_unit=1)
save(joinpath(MA_DIR, "WS/illustrative_plot/legend.svg"), fig; pt_per_unit=1)


##
## Legend narrow/wide
##
width_mm  = 35
height_mm = 15
fig = Figure(resolution = (mm_to_pt(width_mm), mm_to_pt(height_mm)), fontsize=fontsize)

# put axis in a tiny cell so nothing meaningful is visible
ax = Axis(fig[1,1], width=0, height=0)
hidedecorations!(ax); hidespines!(ax)

lines!(ax, [0, 1], [0, 0]; color=:black, linewidth=lw_flow_trajectories, linestyle=:dash, label="wide bounds")
lines!(ax, [0, 1], [0, 0]; color=:black, linewidth=lw_flow_trajectories, label="narrow bounds")
scatter!(ax, (NaN, NaN); marker=:star5,color=:black, markersize=starsize, label="line failure")

# legend as the ONLY visible thing
Legend(fig[1,1], ax; labelsize=fontsize, framevisible=false, patchsize=patchsize, orientation = :vertical)

save(joinpath(MA_DIR, "WS/illustrative_plot/legend_narrow_wide.png"), fig)
save(joinpath(MA_DIR, "WS/illustrative_plot/legend_narrow_wide.pdf"), fig; pt_per_unit=1)
save(joinpath(MA_DIR, "WS/illustrative_plot/legend_narrow_wide.svg"), fig; pt_per_unit=1)


###
### save intermediate files
###
datetime = Dates.format(now(), "_yyyymmdd_HHMMSS.s")
# save different versions of script
script_dest = joinpath(MA_DIR, string("WS/illustrative_plot/old_versions/WS_phase_angles_flow_trajectories_for_inkscape", datetime, ".jl"))
cp(@__FILE__, script_dest; force=true);
