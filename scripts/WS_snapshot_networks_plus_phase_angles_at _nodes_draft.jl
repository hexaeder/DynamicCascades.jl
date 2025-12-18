"""
Sketchy GPT plot. Switched to Inkscape.
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


using CairoMakie
using Colors
using Graphs
using NetworkDynamics
using NetworkLayout: spring

###############  SCHEMATIC: 4-node inset plots + connectors ending at box borders ###############
# - static (CairoMakie)
# - node insets: phase angle θ(t) for narrow (solid black) + wide (dashed black)
# - colored connectors connect BOXES and STOP at box borders
# - "make lines longer" is done by MOVING THE BOXES (spread factor), not extending lines
# - grey stubs: solid near box, dashed outward (subnetwork indication)
#
# REQUIREMENTS:
#   sol1 = narrow SolutionContainer
#   sol4 = wide   SolutionContainer
#
# USAGE:
#   fig = plot_schematic_with_node_insets(sol1, sol4; outpath="schematic_boxes.png")
#   display(fig)
###############################################################################################


# --- extract θ(t) for one node id from a SolutionContainer
function theta_series(c::SolutionContainer, nodeid::Int;
                      tmin=0.0, tmax=20.0, npts=800)
    (nd,) = nd_model(c.network)
    θ_state_idx = idx_containing(nd, "θ")
    θ_node_idx  = map(s -> parse(Int, String(s)[4:end]), nd.syms[θ_state_idx])

    k = findfirst(==(nodeid), θ_node_idx)
    k === nothing && error("Node $nodeid not found among θ states.")
    si = θ_state_idx[k]

    trange = range(tmin, tmax; length=npts)
    θ = Float32[c.sol(t)[si] for t in trange]
    return trange, θ
end

# --- same palette logic as in your script for powerflow curves
function powerflow_color_map(selected_lines::Vector{Int})
    cols = distinguishable_colors(length(selected_lines) + 1)
    deleteat!(cols, 2) # delete yellow
    return Dict(l => cols[i] for (i,l) in pairs(selected_lines))
end

# --- robust finite min/max (avoid NaNs/Infs)
function finite_minmax(v1, v2)
    vals = Float64[]
    append!(vals, Float64.(v1[isfinite.(v1)]))
    append!(vals, Float64.(v2[isfinite.(v2)]))
    isempty(vals) && return nothing
    return minimum(vals), maximum(vals)
end

# --------------------- BBox helpers (no GeometryBasics import) ---------------------
function bb_l_r_b_t(bb)
    l = bb.origin[1]
    b = bb.origin[2]
    r = bb.origin[1] + bb.widths[1]
    t = bb.origin[2] + bb.widths[2]
    return l, r, b, t
end

bb_center(bb) = begin
    l,r,b,t = bb_l_r_b_t(bb)
    Point2f((l+r)/2, (b+t)/2)
end

bb_size(bb) = begin
    l,r,b,t = bb_l_r_b_t(bb)
    (r-l, t-b)
end

function bb_with_center(bb, newc_any)
    newc = Point2f(newc_any[1], newc_any[2])
    w, h = bb_size(bb)
    return BBox(newc[1]-w/2, newc[1]+w/2, newc[2]-h/2, newc[2]+h/2)
end

function bb_side_point(bb, dir::Point2f; pad=6f0)
    l,r,b,t = bb_l_r_b_t(bb)
    c = bb_center(bb)
    d = dir / (norm(dir) + 1f-12)

    ts = Float32[]
    if abs(d[1]) > 1e-9
        push!(ts, (l - c[1]) / d[1]); push!(ts, (r - c[1]) / d[1])
    end
    if abs(d[2]) > 1e-9
        push!(ts, (b - c[2]) / d[2]); push!(ts, (t - c[2]) / d[2])
    end
    ts = filter(x -> x > 0, ts)
    isempty(ts) && return c

    tmin = minimum(ts)
    p = c + tmin*d
    return p - pad*d
end

# --- connector between boxes, ends at the box borders
function connect_bboxes!(scene, bbA, bbB; color=:black, linewidth=6, pad=8f0)
    cA = bb_center(bbA)
    cB = bb_center(bbB)
    dirAB = cB - cA
    pA = bb_side_point(bbA,  dirAB; pad=pad)
    pB = bb_side_point(bbB, -dirAB; pad=pad)
    lines!(scene, [pA[1], pB[1]], [pA[2], pB[2]];
           color=color, linewidth=linewidth, space=:pixel)
    return nothing
end

# --- grey "external network" stub: solid near box, dashed outward
function grey_stub!(scene, bb, dir::Point2f;
                    solid_len=80f0, dashed_len=180f0,
                    linewidth=4, color=RGBAf(0.6,0.6,0.6,1.0),
                    pad=8f0)
    d = dir / (norm(dir) + 1f-12)
    p0 = bb_side_point(bb, d; pad=pad)
    p1 = p0 + solid_len*d
    p2 = p1 + dashed_len*d

    lines!(scene, [p0[1], p1[1]], [p0[2], p1[2]];
           color=color, linewidth=linewidth, linestyle=:solid, space=:pixel)
    lines!(scene, [p1[1], p2[1]], [p1[2], p2[2]];
           color=color, linewidth=linewidth, linestyle=:dash, space=:pixel)
    return nothing
end

# --------------------- layout helpers ---------------------
function spread_boxes(node_bbox_px::Dict{Int,V}; spread=1.0, shifts=Dict{Int,Point2f}()) where {V}
    centers = Dict(k => bb_center(bb) for (k,bb) in node_bbox_px)
    c0 = Point2f(mean(p[1] for p in values(centers)),
                 mean(p[2] for p in values(centers)))

    out = Dict{Int,Any}()
    for (k,bb) in node_bbox_px
        ck = centers[k]
        ck_spread = c0 + Float32(spread) * (ck - c0)
        ck_shift  = get(shifts, k, Point2f(0,0))
        out[k] = bb_with_center(bb, ck_spread + ck_shift)
    end
    return out
end

# HARD guarantee: nothing is truncated on left/bottom.
# This uses the empty space on the right/top by shifting everything inward.
function force_boxes_inside!(node_bbox_px::Dict{Int,Any}; margin=60f0)
    minx = +Inf
    miny = +Inf
    for bb in values(node_bbox_px)
        l,r,b,t = bb_l_r_b_t(bb)
        minx = min(minx, l)
        miny = min(miny, b)
    end

    dx = minx < margin ? (margin - minx) : 0f0
    dy = miny < margin ? (margin - miny) : 0f0

    if dx != 0f0 || dy != 0f0
        for (k,bb) in node_bbox_px
            c = bb_center(bb)
            node_bbox_px[k] = bb_with_center(bb, c + Point2f(dx, dy))
        end
    end
    return node_bbox_px
end

# --------------------- main plotting function ---------------------
function plot_schematic_with_node_insets(sol_narrow::SolutionContainer,
                                         sol_wide::SolutionContainer;
                                         selected_nodes = [4,5,6,82],
                                         selected_lines = [9,12,13,14,15,16],
                                         tmin=0.0, tmax=20.0, npts=800,
                                         resolution = (2200, 1400),
                                         outpath::Union{Nothing,String}=nothing)

    CairoMakie.activate!()
    set_theme!(theme_minimal())

    lcol = powerflow_color_map(selected_lines)

    fig = Figure(resolution=resolution, fontsize=22)

    # Use an LScene as a pure canvas (no axis clipping)
    canvas = LScene(fig[1,1], show_axis=false)

    # ---------------- YOU EDIT THESE BY HAND ----------------
    node_bbox_px_base = Dict(
        82 => BBox(  110,  560,  760, 1120),
        5  => BBox(  110,  560,  380,  740),
        6  => BBox( 540,  990,  120,  480),
        4  => BBox(1030, 1480,  260,  620),
    )

    # shorter connections than before (your "factor 4 instead of 5")
    spread = 1.4

    # node 4 further down: negative y = down (pixel coords)
    shifts = Dict(
        4 => Point2f(0, -60),
    )

    node_bbox_px = spread_boxes(node_bbox_px_base; spread=spread, shifts=shifts)

    # >>> THIS is the truncation fix <<<
    force_boxes_inside!(node_bbox_px; margin=80f0)

    colored_box_connections = [
        (82, 5,  9),
        (5,  4, 12),
        (5,  6, 14),
        (6,  4, 13),
    ]

    # TWO grey lines going left from boxes 82, 5, 6
    external_grey_stubs = [
        (82, Point2f(-1,  0.20)), (82, Point2f(-1, -0.20)),
        (5,  Point2f(-1,  0.20)), (5,  Point2f(-1, -0.20)),
        (6,  Point2f(-1,  0.20)), (6,  Point2f(-1, -0.20)),

    ]
    # --------------------------------------------------------

    # draw grey stubs first
    for (n, d) in external_grey_stubs
        grey_stub!(canvas.scene, node_bbox_px[n], d;
                   solid_len=80f0, dashed_len=180f0,
                   linewidth=4, color=RGBAf(0.6,0.6,0.6,1.0))
    end

    # draw colored connectors (end at box borders)
    for (a, b, lid) in colored_box_connections
        connect_bboxes!(canvas.scene, node_bbox_px[a], node_bbox_px[b];
                        color=get(lcol, lid, :black), linewidth=6)
    end

    # inset axes with θ(t)
    for n in selected_nodes
        bb = node_bbox_px[n]
        ax_in = Axis(fig; bbox=bb)
        hidedecorations!(ax_in); hidespines!(ax_in)

        # black frame
        poly!(ax_in, Rect(0,0,1,1);
              color = RGBAf(0,0,0,0),
              strokecolor = :black,
              strokewidth = 2,
              space = :relative)

        t1, θ1 = theta_series(sol_narrow, n; tmin=tmin, tmax=tmax, npts=npts)
        t2, θ2 = theta_series(sol_wide,   n; tmin=tmin, tmax=tmax, npts=npts)

        m1 = isfinite.(θ1); m2 = isfinite.(θ2)
        lines!(ax_in, t1[m1], θ1[m1]; color=:black, linewidth=2)
        lines!(ax_in, t2[m2], θ2[m2]; color=:black, linewidth=2, linestyle=:dash)

        xlims!(ax_in, tmin, tmax)

        mm = finite_minmax(θ1, θ2)
        if mm === nothing
            ylims!(ax_in, -1, 1)
        else
            ymin, ymax = mm
            if ymin == ymax
                ymin -= 1e-3; ymax += 1e-3
            end
            pad = 0.06 * (ymax - ymin)
            ylims!(ax_in, ymin - pad, ymax + pad)
        end

        text!(ax_in, "Node $n";
              position=(0.02, 0.95), space=:relative,
              align=(:left,:top), color=:black, fontsize=18)
    end

    outpath !== nothing && save(outpath, fig)
    return fig
end


# ---- USAGE:
fig = plot_schematic_with_node_insets(sol1, sol4; outpath="schematic_insets.png")
fig
