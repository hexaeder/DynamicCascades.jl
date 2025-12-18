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
                     ecolorscaling=Observable(1.0), offlinecolor=colorant"lightgray", node_size=55, edge_width=10, 
                     show_angles=false, show_labels=false, inset=nothing) # node size 45
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
                         offlinewidth=10,
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

# small phase-angle disks on top of a network plot
function plot_angle_disks!(ax, c::SolutionContainer, t;
                           disk_size = 40, # 46   # in px
                           line_len  = 0.08)   # in data coords
    
    pos = spring(c.network) # pos = read_pos_or_spring(c.network)
    (nd,) = nd_model(c.network)
    θ_idx = idx_containing(nd, "θ")
    node_idx = map(s -> parse(Int, String(s)[4:end]), nd.syms[θ_idx])

    u = c.sol(t)
    # u0 = c.sol[1] # for phase angle relative to the initial steady state
    θ = fill(0.0, nv(c.network))
    for (si, ni) in zip(θ_idx, node_idx)
        θ[ni] = u[si]
        # θ[ni] = u[si] - u0[si] # for phase angle relative to the initial steady state
    end

    xs = [p[1] for p in pos]
    ys = [p[2] for p in pos]

    # draw disks
    scatter!(ax, xs, ys;
             marker = :circle,
             markersize = disk_size,
             color = :white,
             strokecolor = :black,
             strokewidth = 2)

    # draw centred line + small arrow head at the forward end
    for (x, y, φ) in zip(xs, ys, θ)
        dx = line_len * cos(φ)
        dy = line_len * sin(φ)

        # 1) full line through the disk centre (symmetric)
        x1, y1 = x - dx, y - dy
        x2, y2 = x + dx, y + dy
        lines!(ax, [x1, x2], [y1, y2];
            color = :black, linewidth = 1.5)

        # 2) short arrow segment near the forward end of the line
        frac_start = 0.6      # where the arrow segment starts (60% along the line)
        frac_len   = 0.4      # arrow segment length (40% of the line)

        sx = x + frac_start * dx
        sy = y + frac_start * dy
        vx = frac_len * dx
        vy = frac_len * dy

        arrows!(ax,
            [sx], [sy],        # start position of arrow segment
            [vx], [vy],        # direction vector of arrow segment
            arrowsize = 12,    # head size in px
            linewidth = 1.5,
            color = :black,
        )
    end

    return nothing
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
fontsize = 30
fig = Figure(resolution=(2480,3508), fontsize=fontsize) # DIN A4 at 300 DPI

inset_lims = [2.44, 3.8, -1.2, 0.37] # [2.0, 3.8, -1.2, 0.37] #  # [-3.0, 3.8, -3.1, 1.2] #  # x0, x1, y0, y1  
selected_lines = [9, 12, 13, 14, 15, 16] # selected_lines = sol4.failures.saveval[1:21]

gb = fig[1, 1] = GridLayout()
gc = fig[2, 1] = GridLayout()
gd = fig[3, 1] = GridLayout()

# Tweak relative row size
rowsize!(fig.layout, 1, Relative(0.6))  
# rowsize!(fig.layout, 2, Relative(1))

## network snapshots (4x4 grid)
times = [0.09000000000000001, sol1.failures.t[2], sol1.failures_nodes.t[1], 8.5, 9.3, sol4.failures.t[3], 10.65, 11.0, 11.2, sol4.failures.t[4], 11.7,11.9, sol4.failures.t[5]]
# xticks_ax_gc21 = round.([0.1, sol1.failures_nodes.t[1], 8.3, sol4.failures.t[3], sol4.failures.t[4], sol4.failures.t[5]], digits=1)
xticks_ax_gc21 = round.(times, digits=1)


# Create a 2×1 sub-grid in the region gb[1, 3:4] for colorbars
subgrid = gb[1,7:10] = GridLayout(rows = 3, cols = 1)
subgrid[1,1] = Colorbar(fig, width=Relative(0.85), vertical=false,
            colormap=Makie.ColorScheme([colorant"yellow", colorant"red"]), label="line load relative to rating", tellheight=false)

Δωmax = round(f_b_narrow, digits=2)
subgrid[2,1] = Colorbar(fig, width=Relative(0.85), limits = (-Δωmax,Δωmax), ticks=[-Δωmax, 0.0, Δωmax], vertical=false,
            colormap=ColorSchemes.diverging_bkr_55_10_c35_n256, label="node frequency deviation narrow [Hz]", tellheight=false)

# network initial steady state 
gb[1,1:3] = ax_gb11 = plotnetwork(fig, sol1, times[1]; node_size=10, edge_width=3, show_labels=true)
gb[1,1:3] = Label(fig, "(a) t = $(round(0.0,digits=2)) s", tellwidth=false, tellheight=false, halign=:left, valign=:top)

# rectangle indicating inset
x0 = inset_lims[1]; x1 = inset_lims[2]; y0 = inset_lims[3]; y1 = inset_lims[4]   
xs = [x0, x1, x1, x0, x0]
ys = [y0, y0, y1, y1, y0]
lines!(ax_gb11, xs, ys; color=:black, linewidth=2)

gb[1,4] = plotnetwork(fig, sol1, times[1]; inset=inset_lims, show_labels=true)
gb[1,4] = Label(fig, "(b) t = $(round(0.0,digits=2)) s", tellwidth=false, tellheight=false, halign=:left, valign=:top)

gb[1,5] = plotnetwork(fig, sol1, times[2]; inset=inset_lims, show_labels=false)
gb[1,5] = Label(fig, "(c) t = $(round(times[2],digits=2)) s", tellwidth=false, tellheight=false, halign=:left, valign=:top)

gb[1,6] = plotnetwork(fig, sol1, times[3]; inset=inset_lims, show_labels=false)
gb[1,6] = Label(fig, "(d) t = $(round(times[3],digits=2)) s", tellwidth=false, tellheight=false, halign=:left, valign=:top)

gb[1,1] = Label(fig, "Narrow bounds", font=:bold, textsize=35, tellwidth=false, tellheight=false, halign=:left, valign=:bottom)
gb[2,1] = plotnetwork(fig, sol1, times[4]; inset=inset_lims, show_labels=false)
gb[2,1] = Label(fig, "(e) t = $(round(times[4],digits=2)) s", tellwidth=false, tellheight=false, halign=:left, valign=:top)

gb[2,2] = plotnetwork(fig, sol1, times[5]; inset=inset_lims, show_labels=false)
gb[2,2] = Label(fig, "(f) t = $(round(times[5],digits=2)) s", tellwidth=false, tellheight=false, halign=:left, valign=:top)

gb[2,3] = plotnetwork(fig, sol1, times[6]; inset=inset_lims, show_labels=false)
gb[2,3] = Label(fig, "(g) t = $(round(times[6],digits=2)) s", tellwidth=false, tellheight=false, halign=:left, valign=:top)

gb[2,4] = plotnetwork(fig, sol1, times[7]; inset=inset_lims, show_labels=false)
gb[2,4] = Label(fig, "(h) t = $(round(times[7],digits=2)) s", tellwidth=false, tellheight=false, halign=:left, valign=:top)

gb[2,5] = plotnetwork(fig, sol1, times[8]; inset=inset_lims, show_labels=false)
gb[2,5] = Label(fig, "(i) t = $(round(times[8],digits=2)) s", tellwidth=false, tellheight=false, halign=:left, valign=:top)

gb[2,6] = plotnetwork(fig, sol1, times[9]; inset=inset_lims, show_labels=false)
gb[2,6] = Label(fig, "(j) t = $(round(times[9],digits=2)) s", tellwidth=false, tellheight=false, halign=:left, valign=:top)

gb[2,7] = plotnetwork(fig, sol1, times[10]; inset=inset_lims, show_labels=false)
gb[2,7] = Label(fig, "(k) t = $(round(times[10],digits=2)) s", tellwidth=false, tellheight=false, halign=:left, valign=:top)

gb[2,8] = plotnetwork(fig, sol1, times[11]; inset=inset_lims, show_labels=false)
gb[2,8] = Label(fig, "(l) t = $(round(times[11],digits=2)) s", tellwidth=false, tellheight=false, halign=:left, valign=:top)

gb[2,9] = plotnetwork(fig, sol1, times[12]; inset=inset_lims, show_labels=false)
gb[2,9] = Label(fig, "(m) t = $(round(times[12],digits=2)) s", tellwidth=false, tellheight=false, halign=:left, valign=:top)

gb[2,10] = plotnetwork(fig, sol1, sol1.sol.t[end] - h; inset=inset_lims, show_labels=false)
gb[2,10] = Label(fig, "(n) t → ∞ \n (steady state)", tellwidth=false, tellheight=false, halign=:left, valign=:top)


Δωmax = round(f_b_wide, digits=2)
subgrid[3,1] = Colorbar(fig, width=Relative(0.85), limits = (-Δωmax,Δωmax), ticks=[-Δωmax, 0.0, Δωmax], vertical=false,
            colormap=ColorSchemes.diverging_bkr_55_10_c35_n256, label="node frequency deviation wide [Hz]", tellheight=false)

gb[2,1] = Label(fig, "Wide bounds", font=:bold, textsize=35, tellwidth=false, tellheight=false, halign=:left, valign=:bottom)
gb[3,1] = plotnetwork(fig, sol4, times[4]; inset=inset_lims, show_labels=false)
gb[3,1] = Label(fig, "(e) t = $(round(times[4],digits=2)) s", tellwidth=false, tellheight=false, halign=:left, valign=:top)

gb[3,2] = plotnetwork(fig, sol4, times[5]; inset=inset_lims, show_labels=false)
gb[3,2] = Label(fig, "(f) t = $(round(times[5],digits=2)) s", tellwidth=false, tellheight=false, halign=:left, valign=:top)

gb[3,3] = plotnetwork(fig, sol4, times[6]; inset=inset_lims, show_labels=false)
gb[3,3] = Label(fig, "(g) t = $(round(times[6],digits=2)) s", tellwidth=false, tellheight=false, halign=:left, valign=:top)

gb[3,4] = plotnetwork(fig, sol4, times[7]; inset=inset_lims, show_labels=false)
gb[3,4] = Label(fig, "(h) t = $(round(times[7],digits=2)) s", tellwidth=false, tellheight=false, halign=:left, valign=:top)

gb[3,5] = plotnetwork(fig, sol4, times[8]; inset=inset_lims, show_labels=false)
gb[3,5] = Label(fig, "(i) t = $(round(times[8],digits=2)) s", tellwidth=false, tellheight=false, halign=:left, valign=:top)

gb[3,6] = plotnetwork(fig, sol4, times[9]; inset=inset_lims, show_labels=false)
gb[3,6] = Label(fig, "(j) t = $(round(times[9],digits=2)) s", tellwidth=false, tellheight=false, halign=:left, valign=:top)

gb[3,7] = plotnetwork(fig, sol4, times[10]; inset=inset_lims, show_labels=false)
gb[3,7] = Label(fig, "(k) t = $(round(times[10],digits=2)) s", tellwidth=false, tellheight=false, halign=:left, valign=:top)

gb[3,8] = plotnetwork(fig, sol4, times[11]; inset=inset_lims, show_labels=false)
gb[3,8] = Label(fig, "(l) t = $(round(times[11],digits=2)) s", tellwidth=false, tellheight=false, halign=:left, valign=:top)

gb[3,9] = plotnetwork(fig, sol4, times[12]; inset=inset_lims, show_labels=false)
gb[3,9] = Label(fig, "(m) t = $(round(times[12],digits=2)) s", tellwidth=false, tellheight=false, halign=:left, valign=:top)

gb[3,10] = plotnetwork(fig, sol4, sol4.sol.t[end] - h; inset=inset_lims, show_labels=false)
gb[3,10] = Label(fig, "(n) t → ∞ \n (steady state)", tellwidth=false, tellheight=false, halign=:left, valign=:top)

# Keep row 1 height unchanged, squeeze rows 2 and 3
rowsize!(gb, 1, Auto())           # leave row 1 as it is
rowsize!(gb, 2, Relative(0.25))   # shrink row 2
rowsize!(gb, 3, Relative(0.25))   # shrink row 3

###
### Trajectories power flow
###
line_colors = distinguishable_colors(length(selected_lines)+1)
deleteat!(line_colors, 2) # delete yellow

gc[1,1] = ax_gc11 = Axis(fig;
    xlabel = "time t in s",
    ylabel = "apparent power flow in p.u.",
    title = "Full failure model: narrow frequency bounds (solid) wide frequency bounds (dashed)",
    titlealign = :left,
    xticklabelsize = 10,
    yticklabelsize = 45,
    xlabelsize = 45,
    ylabelsize = 45,
    titlesize = 45,
)

# full failures narrow
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
    lines!(ax_gc11, t, s; linestyle=:dash, color=line_colors[i], linewidth=5)
end
for (i, l) in pairs(selected_lines)
    t, s = seriesforidx(sol4.load_S, l)
    scatter!(ax_gc11, (t[end], s[end]); marker=:star5,color=line_colors[i], markersize=35)
end

rating = get_prop(sol1.network, Graphs.edges(sol1.network), :rating)[1]
hlines!(ax_gc11, rating; color=:black, linestyle=:dot, linewidth=5, label="line rating")
vlines!(ax_gc11, times; color=:black, linestyle=:dot, linewidth=1)

ax_gc11.xticks = xticks_ax_gc21
ax_gc11.xticksvisible = true

xlims!(ax_gc11, -0.1, 18.0)
Legend(gc[1,2], ax_gc11, labelsize=35) # plots show the same lines in all plots



###
### Trajectories phase angles
###
selected_nodes = [4, 5, 6, 82] 
node_colors = distinguishable_colors(length(selected_nodes)+1)
deleteat!(node_colors, 2) # delete yellow

gd[1,1] = ax_gd11 = Axis(fig;
    xlabel = "time t in s",
    ylabel = "phase angle",
    title = "Full failure model: narrow frequency bounds (solid) wide frequency bounds (dashed)",
    titlealign = :left,
    xticklabelsize = 45,
    yticklabelsize = 45,
    xlabelsize = 45,
    ylabelsize = 45,
    titlesize = 45,
)


(nd, p, overload_cb) = nd_model(sol1.network)
state_idx = idx_containing(nd, "θ") # array: indices of θ-states
node_idx = map(s -> parse(Int, String(s)[4:end]), nd.syms[state_idx]) # array: indices of vertices that are generators

# full failures narrow
for (i, l) in pairs([state_idx[findfirst(x -> x == i, node_idx)] for i in selected_nodes])
    t, s = seriesforidx(sol1.sol, l)
    node_idx = selected_nodes[i]
    lines!(ax_gd11, t, s; label="Node $node_idx", color=node_colors[i], linewidth=5)
    scatter!(ax_gd11, (t[end], s[end]); color=node_colors[i], marker=:star5, markersize=25)
end

node_idx = map(s -> parse(Int, String(s)[4:end]), nd.syms[state_idx]) # array: indices of vertices that are generators
# full failure wide
for (i, l) in pairs([state_idx[findfirst(x -> x == i, node_idx)] for i in selected_nodes])
    t, s = seriesforidx(sol4.sol, l)
    node_idx = selected_nodes[i]
    lines!(ax_gd11, t, s; linestyle=:dash, color=node_colors[i], linewidth=5)
    scatter!(ax_gd11, (t[end], s[end]); color=node_colors[i], marker=:star5, markersize=25)
end

vlines!(ax_gd11, sol1.failures_nodes.t; color=:black, linewidth=2)
vlines!(ax_gd11, sol4.failures_nodes.t; color=:black, linestyle=:dash, linewidth=2)

xlims!(ax_gd11, -0.1, 20)
ylims!(ax_gd11, -1, 2.0)
Legend(gd[1,2], ax_gd11, labelsize=35) # plots show the same lines in all plots

datetime = Dates.format(now(), "_yyyymmdd_HHMMSS.s")
# save different versions of script
script_dest = joinpath(MA_DIR, string("WS/trajectories/WS_snapshot_network_different_node_models_script", datetime, ".jl"))
cp(@__FILE__, script_dest; force=true)

save(joinpath(MA_DIR, string("WS/trajectories/WS_snapshots_network+trajectories_different_models_narrow,full_failure_wide_tweak_power_injections=$tweak_power_injections,", string_plot_params, datetime, ".pdf")), fig)
save(joinpath(MA_DIR, string("WS/trajectories/WS_snapshots_network+trajectories_different_models_narrow,full_failure_wide_tweak_power_injections=$tweak_power_injections,", string_plot_params, datetime, ".png")), fig)
fig