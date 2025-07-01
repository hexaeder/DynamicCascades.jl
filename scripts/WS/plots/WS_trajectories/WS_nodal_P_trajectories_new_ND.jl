"""
Watts-Strogatz-Network: Related to the curve in WS_lines+nodes_uebergang.jl with
frequency bound f_b=0.03. Plotting all frequency and power flow trajectories of all
nodes and lines that fail for I=[0.2, 3.0, 30.0].
"""

include(abspath(@__DIR__, "..", "helpers_jarray.jl"))


using DynamicCascades
using Colors
using CairoMakie # for normal plots
CairoMakie.activate!()

###
### save solution objects
###
initial_fail = 78
task_id = 1720
task_id_array = [1720, 1723, 1728]

exp_name_date = "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250321_171511.976"
exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)

for task_id in task_id_array
    sol, nw = simulate(exp_data_dir, task_id, initial_fail;
        tspan=(0., 35.),
        solverargs = (;dtmax=0.01),
        verbose = true);

    Serialization.serialize(joinpath(exp_data_dir, "trajectories", "task_id=$task_id,initial_fail=$initial_fail.sol"), sol)
end


###
### plot trajectories
###

# load solution objects
sol1720 = Serialization.deserialize(joinpath(exp_data_dir, "trajectories", "task_id=1720,initial_fail=78.sol"));
sol1723 = Serialization.deserialize(joinpath(exp_data_dir, "trajectories", "task_id=1723,initial_fail=78.sol"));
sol1728 = Serialization.deserialize(joinpath(exp_data_dir, "trajectories", "task_id=1728,initial_fail=78.sol"));

# 78: initial trigger line
lines_task_id_1720 = [78] # sol1720.failures.saveval
lines_task_id_1723 = [78] # sol1723.failures.saveval
lines_task_id_1728 = [78, 199, 194, 131, 147, 146, 14, 106, 135, 149, 148] # sol1728.failures.saveval
nodes_task_id_1720 = [29, 98, 92, 58] # sol1720.failures_nodes.saveval
nodes_task_id_1723 = [29, 98] # sol1723.failures_nodes.saveval
nodes_task_id_1728 = [98, 59, 62, 61, 60] # sol1728.failures_nodes.saveval

all_failing_lines_idxs = [14, 78, 106, 131, 135, 146, 147, 148, 149, 194, 199]
all_failing_nodes_idxs = [29, 58, 59, 60, 61, 62, 92, 98]
# all_failing_lines_idxs = Int64[1:200...]
# all_failing_nodes_idxs = Int64[1:100...]

node_colors = distinguishable_colors(length(all_failing_nodes_idxs)+1)
deleteat!(node_colors, 2) # delete yellow
line_colors = distinguishable_colors(length(all_failing_lines_idxs)+1)
deleteat!(line_colors, 2) # delete yellow

################################################################################
############################ Line and nodes ####################################
################################################################################
fontsize = 35
titlesize = (fontsize+5)
linewidth = 3.5
# fig = Figure(resolution=(2100,1500), fontsize= fontsize)
# xlim_30 = 36.0
fig = Figure(size=(3100,1500), fontsize= fontsize)
xlim_30 = 35.0

# FREQUENCIES ########################################################################
# frequencies of failed gen nodes I=0.2 s^2
sol = sol1720
fig[1,1] = ax = Axis(fig; ylabel="Frequency [Hz]", title=L"Inertia $I=0.2$ $s^2$", titlealign = :left, titlesize = titlesize)
for i in all_failing_nodes_idxs
    x, y = remove_zero_tail!(deepcopy(sol.t), deepcopy(map(first, sol(sol.t, idxs=vidxs(i, :f)))))
    color = node_colors[findfirst(x -> x == i, all_failing_nodes_idxs)]
    lines!(ax, x, y; label="Node $i", color=color, linewidth=linewidth) # consider plotdensity=Int(1e5)
    scatter!(ax, (x[end], y[end]); color=color, marker=:star5, markersize=25)
end
xlims!(ax, 0, 1.1)
ax.xticks = [0.0, 0.25, 0.5, 0.75, 1.0]
ax.yticks = [-0.03, -0.02, -0.01, 0.00, 0.01, 0.02, 0.03]
# axislegend(ax, position = :rt, nbanks = 3)

# frequencies of failed gen nodes I=3.0 s^2
sol = sol1723
fig[1,2] = ax = Axis(fig; ylabel="Frequency [Hz]", title=L"Inertia $I=3.0$ $s^2$", titlealign = :left, titlesize = titlesize)
for i in all_failing_nodes_idxs
    x, y = remove_zero_tail!(deepcopy(sol.t), deepcopy(map(first, sol(sol.t, idxs=vidxs(i, :f)))))
    color = node_colors[findfirst(x -> x == i, all_failing_nodes_idxs)]
    lines!(ax, x, y; label="Node $i", color=color, linewidth=linewidth) # consider plotdensity=Int(1e5)
    scatter!(ax, (x[end], y[end]); color=color, marker=:star5, markersize=25)
end
xlims!(ax, 0, 4.1)
ax.xticks = [0.1, 1, 2, 3, 4]
ax.yticks = [-0.03, -0.02, -0.01, 0.00, 0.01, 0.02, 0.03]
# axislegend(ax, position = :rt, nbanks = 2)

# frequencies of failed gen nodes I=30.0 s^2
sol = sol1728
fig[1,3] = ax = Axis(fig; ylabel="Frequency [Hz]", title=L"Inertia $I=30.0$ $s^2$", titlealign = :left, titlesize = titlesize)
for i in all_failing_nodes_idxs
    x, y = remove_zero_tail!(deepcopy(sol.t), deepcopy(map(first, sol(sol.t, idxs=vidxs(i, :f)))))
    color = node_colors[findfirst(x -> x == i, all_failing_nodes_idxs)]
    lines!(ax, x, y; label="Node $i", color=color, linewidth=linewidth) # consider plotdensity=Int(1e5)
    scatter!(ax, (x[end], y[end]); color=color, marker=:star5, markersize=25)
end
xlims!(ax, -1., xlim_30)
# ax.xticks = [0.0, 5., 10., 15., 20., 25., 30., 35.]
ax.yticks = [-0.03, -0.02, -0.01, 0.00, 0.01, 0.02, 0.03]
# axislegend(ax, position = :rt, nbanks = 2)

# FLOWS ########################################################################
# failed power flows I=0.2 s^2
sol = sol1720
fig[3,1] = ax = Axis(fig; xlabel="Time [s]", ylabel="App. power flow [p.u.]")
for i in all_failing_lines_idxs
    x, y = remove_zero_tail!(deepcopy(sol.t), deepcopy(map(first, sol(sol.t, idxs=eidxs(i, :S)))))
    color = line_colors[findfirst(x -> x == i, all_failing_lines_idxs)]
    lines!(ax, x, y; label="Node $i", color=color, linewidth=linewidth) # consider plotdensity=Int(1e5)
    scatter!(ax, (x[end], y[end]); color=color, marker=:star5, markersize=25)
end
xlims!(ax, 0, 1.1)
ax.xticks = [0.0, 0.25, 0.5, 0.75, 1.0]
# hlines!(ax, 2.1; color=:red, linestyle=:dash, linewidth=linewidth, label="Rating")
# axislegend(ax, position = (0.9,0.92), nbanks = 6, labelsize=22)

# failed power flows I=3.0 s^2
sol = sol1723
fig[3,2] = ax = Axis(fig; xlabel="Time [s]", ylabel="App. power flow [p.u.]")
for i in all_failing_lines_idxs
    x, y = remove_zero_tail!(deepcopy(sol.t), deepcopy(map(first, sol(sol.t, idxs=eidxs(i, :S)))))
    color = line_colors[findfirst(x -> x == i, all_failing_lines_idxs)]
    lines!(ax, x, y; label="Node $i", color=color, linewidth=linewidth) # consider plotdensity=Int(1e5)
    scatter!(ax, (x[end], y[end]); color=color, marker=:star5, markersize=25)
end
xlims!(ax, 0, 4.1)
# ax.xticks = [0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5]

# axislegend(ax, position = :rt, nbanks = 3)

# failed power flows I=30.0 s^2
sol = sol1728
fig[3,3] = ax = Axis(fig; xlabel="Time [s]", ylabel="App. power flow [p.u.]")
for i in all_failing_lines_idxs
    x, y = remove_zero_tail!(deepcopy(sol.t), deepcopy(map(first, sol(sol.t, idxs=eidxs(i, :S)))))
    color = line_colors[findfirst(x -> x == i, all_failing_lines_idxs)]
    lines!(ax, x, y; label="Node $i", color=color, linewidth=linewidth) # consider plotdensity=Int(1e5)
    scatter!(ax, (x[end], y[end]); color=color, marker=:star5, markersize=25)
end
xlims!(ax, -1., xlim_30)
# ax.xticks = [0.0, 5., 10., 15., 20., 25., 30., 35.]
# axislegend(ax, position = :rt, nbanks = 3)


# ΔP TRAJECTORIES ##############################################################

# ΔP trajectories I=0.2 s^2
sol = sol1720
task_id = 1720
fig[2,1] = ax = Axis(fig; xlabel="Time [s]", ylabel="ΔP [p.u.]")
for i in all_failing_nodes_idxs
    x, y = remove_zero_tail!(deepcopy(sol.t), deepcopy(map(first, sol(sol.t, idxs=vidxs(i, :ΔP)))))
    color = node_colors[findfirst(x -> x == i, all_failing_nodes_idxs)]
    lines!(ax, x, y; label="Node $i", color=color, linewidth=linewidth) # consider plotdensity=Int(1e5)
    scatter!(ax, (x[end], y[end]); color=color, marker=:star5, markersize=25)
end
xlims!(ax, 0, 1.1)
ax.xticks = [0.0, 0.25, 0.5, 0.75, 1.0]
# hlines!(ax, 2.1; color=:red, linestyle=:dash, linewidth=linewidth, label="Rating")
# axislegend(ax, position = (0.9,0.92), nbanks = 6, labelsize=22)

# ΔP trajectories I=3.0 s^2
sol = sol1723
task_id = 1723
fig[2,2] = ax = Axis(fig; xlabel="Time [s]", ylabel="ΔP [p.u.]")
for i in all_failing_nodes_idxs
    x, y = remove_zero_tail!(deepcopy(sol.t), deepcopy(map(first, sol(sol.t, idxs=vidxs(i, :ΔP)))))
    color = node_colors[findfirst(x -> x == i, all_failing_nodes_idxs)]
    lines!(ax, x, y; label="Node $i", color=color, linewidth=linewidth) # consider plotdensity=Int(1e5)
    scatter!(ax, (x[end], y[end]); color=color, marker=:star5, markersize=25)
end
xlims!(ax, 0, 4.1)
# ax.xticks = [0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5]

# ΔP trajectories I=30.0 s^2
sol = sol1728
task_id = 1728
fig[2,3] = ax = Axis(fig; xlabel="Time [s]", ylabel="ΔP [p.u.]")
for i in all_failing_nodes_idxs
    x, y = remove_zero_tail!(deepcopy(sol.t), deepcopy(map(first, sol(sol.t, idxs=vidxs(i, :ΔP)))))
    color = node_colors[findfirst(x -> x == i, all_failing_nodes_idxs)]
    lines!(ax, x, y; label="Node $i", color=color, linewidth=linewidth) # consider plotdensity=Int(1e5)
    scatter!(ax, (x[end], y[end]); color=color, marker=:star5, markersize=25)
end
xlims!(ax, -1., xlim_30)
# ax.xticks = [0.0, 5., 10., 15., 20., 25., 30., 35.]get_ΔP_ΔS_at_node(99, sol, network)

CairoMakie.save(joinpath(exp_data_dir, "trajectories", "WS_traj,lines+nodes,task_ids=$task_id_array,init_fail=$initial_fail,_DeltaP_failing_lines+nodes_new_ND.pdf"),fig)
CairoMakie.save(joinpath(exp_data_dir, "trajectories", "WS_traj,lines+nodes,task_ids=$task_id_array,init_fail=$initial_fail,_DeltaP_failing_lines+nodes_new_ND.png"),fig)
fig
