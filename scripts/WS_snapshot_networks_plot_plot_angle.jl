"""
Snapshots of illustrative example network. Phase angles as disks on nodes.
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


function plotnetwork(fig, sol, t; line=nothing, offset=nothing, text=nothing, apos=.5, Δω=Observable(Δωmax), ecolortype=:relrating,
                     ecolorscaling=Observable(1.0), offlinecolor=colorant"lightgray", node_size=40, edge_width=5, 
                     show_angles=true, show_labels=false, inset=nothing)
    # ax = Axis(fig)
    ax = Axis(fig, aspect = DataAspect())   # 1:1 scaling in x and y

    hidedecorations!(ax), hidespines!(ax)
    if inset !== nothing
        xlims!(ax, inset_lims[1], inset_lims[2])
        ylims!(ax, inset_lims[3], inset_lims[4])
    end

    gpargs = gparguments(sol, t isa Observable ? t : Observable(t);
                         Δω=Δω,
                         ecolortype=Observable(ecolortype),
                         offlinewidth=5,
                         offlinecolor=offlinecolor,
                         ecolorscaling=ecolorscaling,
                         show_labels=show_labels,
                         node_size=node_size,
                         edge_width=edge_width,
                         layout=spring # NOTE for adapting C use `set_prop!(sol.network, 1:nv(sol.network), :pos, spring(sol.network; C=1.0))`
                         )
    p = graphplot!(ax, sol.network; gpargs...)

    # overlay phase-angle disks
    if show_angles
        plot_angle_disks!(ax, sol, t)
    end

    if line !==nothing
        pos = GraphMakie.interpolate(get_edge_plot(p).paths[][line], apos)
        offset = Point2(offset)
        gapA = 0.0
        gapB = 0.6
        posA = pos + gapA*normalize(offset)
        posT = pos + offset
        posB = posT - gapB*normalize(offset)
        lines!(ax, [posA, posB], color=:black, linewidth=2)
        text!(ax, text; position=posT, align=(:center, :center))
    end
    return ax
end





## collect data
exp_name_date_full = "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250321_171511.976"
exp_name_date_inertia = "WS_k=4_exp11_vary_I_only_lines_and_nodes_change_to_BH_complement_f_bPIK_HPC_K_=3,N_G=32_20250420_145304.04"
exp_name_date_power = "WS_k=4_exp12_vary_I_only_lines_and_nodes_change_Pmech_complement_f_bPIK_HPC_K_=3,N_G=32_20250420_142151.671"

initial_fail = 15
task_id_full_n = 7584
task_id_inertia_n = 2229
task_id_power_n = 2229
task_id_full_w = 7755

## string for saving plot
exp_nr_full = exp_name_date_full[11:12]
exp_nr_BH = exp_name_date_inertia[11:12]
exp_nr_Pmech = exp_name_date_power[11:12]
prefix = exp_name_date_full[1:10]

df_config = DataFrame(CSV.File(joinpath(RESULTS_DIR, exp_name_date_full, "config.csv")));
f_b_narrow=df_config[task_id_full_n,:freq_bounds]
f_b_wide=df_config[task_id_full_w,:freq_bounds]
I=df_config[task_id_full_n,:inertia_values]

tweak_power_injections = true
# Network from ensemble
# sol1 = simulate_wrapper(exp_name_date_full, task_id_full_n, initial_fail; save_simulation = true, tweak_power_injections=tweak_power_injections); # narrow
# sol2 = simulate_wrapper(exp_name_date_inertia, task_id_inertia_n, initial_fail; save_simulation = true, tweak_power_injections=tweak_power_injections); # narrow
# sol3 = simulate_wrapper(exp_name_date_power, task_id_power_n, initial_fail; save_simulation = true, tweak_power_injections=tweak_power_injections); # narrow
# sol4 = simulate_wrapper(exp_name_date_full, task_id_full_w, initial_fail; save_simulation = true, tweak_power_injections=tweak_power_injections); # wide

# load .sol-files
sol1 = Serialization.deserialize(joinpath(RESULTS_DIR, exp_name_date_full, "trajectories", "task_id=$task_id_full_n,initial_fail=$initial_fail,tweak_power_injections=$tweak_power_injections.sol"));
sol2 = Serialization.deserialize(joinpath(RESULTS_DIR, exp_name_date_inertia, "trajectories", "task_id=$task_id_inertia_n,initial_fail=$initial_fail,tweak_power_injections=$tweak_power_injections.sol"));
sol3 = Serialization.deserialize(joinpath(RESULTS_DIR, exp_name_date_power, "trajectories", "task_id=$task_id_power_n,initial_fail=$initial_fail,tweak_power_injections=$tweak_power_injections.sol"));
sol4 = Serialization.deserialize(joinpath(RESULTS_DIR, exp_name_date_full, "trajectories", "task_id=$task_id_full_w,initial_fail=$initial_fail,tweak_power_injections=$tweak_power_injections.sol"));

# # find max ω (old way)
# (nd,) = nd_model(sol1.network);
# state_idx = idx_containing(nd, "ω");
# # CHECK exp_name_date_power: trajectories might exceed frequency bounds!
# Δωmax1 = round(maximum(map(x -> max(abs.(extrema(x[state_idx]))...), sol1.sol)) / (2π), digits=2)
# Δωmax2 = round(maximum(map(x -> max(abs.(extrema(x[state_idx]))...), sol2.sol)) / (2π), digits=2)
# Δωmax3 = round(maximum(map(x -> max(abs.(extrema(x[state_idx]))...), sol3.sol)) / (2π), digits=2)
# Δωmax4 = round(maximum(map(x -> max(abs.(extrema(x[state_idx]))...), sol4.sol)) / (2π), digits=2)
# Δωmax = round(maximum([Δωmax1, Δωmax2, Δωmax3]), digits=2)

#= times network snapshots (for all plots in this script)
#NOTE see `gparguments`: node colors are calculated via simple derivation of phase angles. 
This does not work for `sol.sol.t[1]` and `sol.sol.t[end]`.=#
h = 0.01


exp_nrs = "$prefix,$exp_nr_full,tid_n=$task_id_full_n,tid_w=$task_id_full_w,$exp_nr_BH,tid=$task_id_inertia_n,$exp_nr_Pmech,tid=$task_id_power_n"
string_plot_params = "exp=$exp_nrs,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I,init_fail=$initial_fail"

## fig parameters
CairoMakie.activate!()
set_theme!(theme_minimal())
fontsize = 35
fig = Figure(resolution=(2480,3508), fontsize=fontsize) # DIN A4 at 300 DPI

inset_lims = [-3.0, 3.8, -3.1, 1.2] # [2.495, 3.8, -1.2, 1.2] # x0, x1, y0, y1  # [-3.0, 3.8, -3.1, 1.2] 
selected_lines = [9, 13, 14, 15, 16] # selected_lines = sol4.failures.saveval[1:21]

gb = fig[1, 1] = GridLayout()
gc = fig[2, 1] = GridLayout()

# Tweak relative row size
rowsize!(fig.layout, 1, Relative(0.6))  
# rowsize!(fig.layout, 2, Relative(1))

## network snapshots (4x4 grid)
times = [0.09000000000000001, 6.508362795315281, 8.3, sol4.failures.t[3], sol4.failures.t[4], sol4.failures.t[5]]
xticks_ax_gc21 = round.([0.1, 6.508362795315281, 8.3, sol4.failures.t[3], sol4.failures.t[4], sol4.failures.t[5]], digits=1)

sol4.failures.t[3]

# Create a 2×1 sub-grid in the region gb[1, 3:4] for colorbars
subgrid = gb[1,3:4] = GridLayout(rows = 3, cols = 1)
subgrid[1,1] = Colorbar(fig, width=Relative(0.85), vertical=false,
            colormap=Makie.ColorScheme([colorant"yellow", colorant"red"]), label="line load relative to rating", tellheight=false)

Δωmax = round(f_b_narrow, digits=2)
subgrid[2,1] = Colorbar(fig, width=Relative(0.85), limits = (-Δωmax,Δωmax), ticks=[-Δωmax, 0.0, Δωmax], vertical=false,
            colormap=ColorSchemes.diverging_bkr_55_10_c35_n256, label="node frequency deviation narrow [Hz]", tellheight=false)

# network initial steady state 
gb[1,1] = ax_gb11 = plotnetwork(fig, sol1, times[1]; node_size=10, edge_width=3, show_labels=false)
gb[1,1] = Label(fig, "(a) t = $(round(0.0,digits=2)) s", tellwidth=false, tellheight=false, halign=:right, valign=:top)

# rectangle indicating inset
x0 = inset_lims[1]; x1 = inset_lims[2]; y0 = inset_lims[3]; y1 = inset_lims[4]   
xs = [x0, x1, x1, x0, x0]
ys = [y0, y0, y1, y1, y0]
lines!(ax_gb11, xs, ys; color=:black, linewidth=2)

gb[1,2] = plotnetwork(fig, sol1, times[2]; inset=inset_lims, show_labels=false)
gb[1,2] = Label(fig, "(b) t = $(round(times[2],digits=2)) s \n (first node failure narrow)", tellwidth=false, tellheight=false, halign=:right, valign=:bottom, textsize=fontsize-5)


# network final state
gb[2,1] = plotnetwork(fig, sol1, times[3]; inset=inset_lims, show_labels=true)
gb[2,1] = Label(fig, "(cFn) t = $(round(times[3],digits=2)) s", tellwidth=false, tellheight=false, halign=:right, valign=:bottom)

gb[2,2] = plotnetwork(fig, sol1, times[4]; inset=inset_lims, show_labels=false)
gb[2,2] = Label(fig, "(dFn) t = $(round(times[4],digits=2)) s", tellwidth=false, tellheight=false, halign=:right, valign=:bottom)

gb[2,3] = plotnetwork(fig, sol1, times[5]; inset=inset_lims, show_labels=false)
gb[2,3] = Label(fig, "(eFn) t = $(round(times[5],digits=2)) s", tellwidth=false, tellheight=false, halign=:right, valign=:bottom)

gb[2,4] = plotnetwork(fig, sol1, sol1.sol.t[end] - h; inset=inset_lims, show_labels=false)
gb[2,4] = Label(fig, "(fFn) t → ∞ \n (steady state)", tellwidth=false, tellheight=false, halign=:right, valign=:bottom)




Δωmax = round(f_b_wide, digits=2)
subgrid[3,1] = Colorbar(fig, width=Relative(0.85), limits = (-Δωmax,Δωmax), ticks=[-Δωmax, 0.0, Δωmax], vertical=false,
            colormap=ColorSchemes.diverging_bkr_55_10_c35_n256, label="node frequency deviation wide [Hz]", tellheight=false)

gb[3,1] = plotnetwork(fig, sol4, times[3]; inset=inset_lims, show_labels=false)
gb[3,1] = Label(fig, "(cFw) t = $(round(times[3],digits=2)) s", tellwidth=false, tellheight=false, halign=:right, valign=:bottom)

gb[3,2] = plotnetwork(fig, sol4, times[4]; inset=inset_lims, show_labels=false)
gb[3,2] = Label(fig, "(dFw) t = $(round(times[4],digits=2)) s", tellwidth=false, tellheight=false, halign=:right, valign=:bottom)

gb[3,3] = plotnetwork(fig, sol4, times[5]; inset=inset_lims, show_labels=false)
gb[3,3] = Label(fig, "(eFw) t = $(round(times[5],digits=2)) s", tellwidth=false, tellheight=false, halign=:right, valign=:bottom)

gb[3,4] = plotnetwork(fig, sol4, sol4.sol.t[end] - h; inset=inset_lims, show_labels=false)
gb[3,4] = Label(fig, "(fFw) t → ∞ \n (steady state)", tellwidth=false, tellheight=false, halign=:right, valign=:bottom)


line_colors = distinguishable_colors(length(selected_lines)+1)
deleteat!(line_colors, 2) # delete yellow
# color = line_colors[findfirst(x -> x == i, all_failing_lines_idxs)]

## trajectories
# rating = get_prop(sol1.network, Graphs.edges(sol1.network), :rating)[1]
# gc[1,1] = ax_gc11 = Axis(fig; xlabel="time t in s", title="full failures (solid lines), inertia failures (dashed lines) and power failures (dotted lines) narrow", titlealign=:left)





# # inertia failures
# # vlines!(ax_gc11, sol2.failures_nodes.t; color=:black, linestyle=:dash, linewidth=1) # BUG add Legend/explanation in plot
# for (i, l) in pairs(selected_lines)
#     t, s = seriesforidx(sol2.load_S, l)
#     lines!(ax_gc11, t, s; linestyle=:dash,color=line_colors[i],  linewidth=5)
# end
# for (i, l) in pairs(selected_lines)
#     t, s = seriesforidx(sol2.load_S, l)
#     scatter!(ax_gc11, (t[end], s[end]); marker=:star5,color=line_colors[i],  markersize=35)
# end

# # power failures
# # vlines!(ax_gc11, sol3.failures_nodes.t; color=:black, linestyle=:dot, linewidth=1, label="node failure \n (vertical lines)")
# for (i, l) in pairs(selected_lines)
#     t, s = seriesforidx(sol3.load_S, l)
#     lines!(ax_gc11, t, s; linestyle=:dot, color=line_colors[i],  linewidth=5)
# end
# for (i, l) in pairs(selected_lines)
#     t, s = seriesforidx(sol3.load_S, l)
#     scatter!(ax_gc11, (t[end], s[end]); marker=:star5,color=line_colors[i],  markersize=35)
# end

# explanations in plot
# text!(ax_gc11, 25, 0.45; text="↓ full failure model", textsize=fontsize)
# text!(ax_gc11, 25, 0.10; text="↑ inertia failure model", textsize=fontsize)
# text!(ax_gc11, 14.5, 0.95; text="← power failure model", textsize=fontsize)




gc[1,1] = ax_gc11 = Axis(fig; xlabel="time t in s", ylabel="apparent power flow in p.u.", title="full failure model: narrow frequency bounds (solid) wide frequency bounds (dashed)", titlealign=:left)
# hlines!(ax_gc11, get_prop(sol1.network, Graphs.edges(sol1.network), :rating)[1]; color=:black, linestyle=:dash, linewidth=5, label="line rating") # sol1 as rating is the same
# vlines!(ax_gc11, sol4.failures_nodes.t; color=:black, linewidth=1, label="node failure \n (vertical lines)")

# full failures narrow
# vlines!(ax_gc11, sol1.failures_nodes.t; color=:black, linewidth=1, label="node failure \n (vertical lines)")
for (i, l) in pairs(selected_lines)
    t, s = seriesforidx(sol1.load_S, l)
    lines!(ax_gc11, t, s; label="line $l",color=line_colors[i], linewidth=5)
end
for (i, l) in pairs(selected_lines)
    t, s = seriesforidx(sol1.load_S, l)
    scatter!(ax_gc11, (t[end], s[end]); marker=:star5,color=line_colors[i],  markersize=35)
end

# full failure wide
for (i, l) in pairs(selected_lines)
    t, s = seriesforidx(sol4.load_S, l)
    lines!(ax_gc11, t, s; linestyle=:dash, color=line_colors[i],  linewidth=5)
end
for (i, l) in pairs(selected_lines)
    t, s = seriesforidx(sol4.load_S, l)
    scatter!(ax_gc11, (t[end], s[end]); marker=:star5,color=line_colors[i],  markersize=35)
end

hlines!(ax_gc11, rating; color=:black, linestyle=:dot, linewidth=5, label="line rating")


# text!(ax_gc11, sol4.failures_nodes.t[2]+0.3, 0.2; text="← vertical lines indicate \n     node failures", textsize=fontsize)
# text!(ax_gc11, sol4.failures.t[7]-0.2, 1.9; text="↑ further line failures", textsize=fontsize)
# additional_stars_wide_bounds = []
# for (i, line) in enumerate(sol4.failures.saveval[1:21])
#     if line in selected_lines
#         nothing
#     else
#         push!(additional_stars_wide_bounds, line)
#     end
# end
# for (i, l) in pairs(additional_stars_wide_bounds)
#     t, s = seriesforidx(sol4.load_S, l)
#     scatter!(ax_gc11, (t[end], s[end]); marker=:star5,color=:grey,  markersize=35)
# end

ax_gc11.xticks = xticks_ax_gc21
ax_gc11.xticksvisible = true

xlims!(ax_gc11, -0.1, 18.0)
# axs = [ax_gc11, ax_gc21]
# linkxaxes!(axs...)
# hidexdecorations!.(axs[1:end-1], grid=false)
# for i in axs
#     xlims!(i, -0.1, 25.0)
# end
Legend(gc[1,2], ax_gc11)# plots show the same lines in all plots

vlines!(ax_gc11, xticks_ax_gc21; color=:black, linestyle=:dot, linewidth=1, label="node failure \n (vertical lines)")

datetime = Dates.format(now(), "_yyyymmdd_HHMMSS.s")
# save different versions of script
script_dest = joinpath(MA_DIR, string("WS/trajectories/WS_snapshot_network_different_node_models_script", datetime, ".jl"))
cp(@__FILE__, script_dest; force=true)

save(joinpath(MA_DIR, string("WS/trajectories/WS_snapshots_network+trajectories_different_models_narrow,full_failure_wide_tweak_power_injections=$tweak_power_injections,", string_plot_params, datetime, ".pdf")), fig)
save(joinpath(MA_DIR, string("WS/trajectories/WS_snapshots_network+trajectories_different_models_narrow,full_failure_wide_tweak_power_injections=$tweak_power_injections,", string_plot_params, datetime, ".png")), fig)
fig