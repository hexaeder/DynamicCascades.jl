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

# using NetworkDynamicsInspector

using Graphs
using Unitful
using Statistics
using GraphMakie
 
include(abspath(@__DIR__, "paper_plots_helpers_and_parameters.jl"))

###
### save solution objects
###
initial_fail = 78
task_id = 1720
task_id_array = [1720, 1723, 1728]

exp_name_date = "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250321_171511.976"
exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)

# # NOTE For resimulating set `res_tol=1e-6,` in ND_model.jl
# for task_id in task_id_array
#     sol, nw = simulate(exp_name_date, task_id, initial_fail;
#         tspan=(0., 35.),
#         solverargs = (;dtmax=0.01),
#         verbose = true);

#     Serialization.serialize(joinpath(exp_data_dir, "trajectories_paper", "task_id=$task_id,initial_fail=$initial_fail.sol"), sol)
# end

###
### plot trajectories
###

# load solution objects
sol1720 = Serialization.deserialize(joinpath(exp_data_dir, "trajectories_paper", "task_id=1720,initial_fail=78.sol"));
sol1723 = Serialization.deserialize(joinpath(exp_data_dir, "trajectories_paper", "task_id=1723,initial_fail=78.sol"));
sol1728 = Serialization.deserialize(joinpath(exp_data_dir, "trajectories_paper", "task_id=1728,initial_fail=78.sol"));

# # 78: initial trigger line
# lines_task_id_1720 = [78] # sol1720.failures.saveval
# lines_task_id_1723 = [78] # sol1723.failures.saveval
# lines_task_id_1728 = [78, 199, 194, 131, 147, 146, 14, 106, 135, 149, 148] # sol1728.failures.saveval
# nodes_task_id_1720 = [29, 98, 92, 58] # sol1720.failures_nodes.saveval
# nodes_task_id_1723 = [29, 98] # sol1723.failures_nodes.saveval
# nodes_task_id_1728 = [98, 59, 62, 61, 60] # sol1728.failures_nodes.saveval

all_failing_lines_idxs = [14, 78, 106, 131, 135, 146, 147, 148, 149, 194, 199]
all_failing_nodes_idxs = [29, 58, 59, 60, 61, 62, 92, 98]

# node_colors = distinguishable_colors(9); deleteat!(node_colors, 2) # delete yellow
# line_colors = distinguishable_colors(12); deleteat!(line_colors, 2) # delete yellow

colors = distinguishable_colors(20)
deleteat!(colors, 2) # delete yellow
line_colors = colors[1:11]
node_colors = colors[12:19]

line_colors[1:3]
################################################################################
############################ Line and nodes ####################################
################################################################################
# NOTE different approaches
# # This does not work
# lines!(ax, sol.t, sol(sol.t, idxs=vidxs(all_failing_nodes_idxs, :ω)))
# lines!(ax, sol.t, sol(sol.t, idxs=vidxs(1, :ω)))
# lines!(ax, sol.t, sol(sol.t, idxs=VIndex(1,:ω)))
# # This also works
# for i in all_failing_nodes_idxs
#     lines!(sol.t, (sol(sol.t, idxs=VIndex(i,2)).u)./(2*π);
#         label="Node $i",
#         color=node_colors[findfirst(x -> x == i, all_failing_nodes_idxs)], 
#         linewidth=linewidth)    
# end
# And this works: `lines!(ax, sol, idxs=vidxs(i, :f), plotdensity=Int(1e5))`

fontsize = titlesize = 36
linewidth = 4
patchsize = (35,9)
fig = Figure(size=(2100,1000), fontsize= fontsize)

# FREQUENCIES ########################################################################
# frequencies of failed gen nodes I=0.2 s^2
sol = sol1720;
fig[1,1] = ax11 = Axis(fig; ylabel="Frequency [Hz]", title=L"Inertia $I=0.2$ $s^2$", titlealign = :left, titlesize = titlesize)
hlines!(ax11, 0.03; color=:black, linestyle=:dash, linewidth=linewidth, label="frequency bound")
hlines!(ax11, -0.03; color=:black, linestyle=:dash, linewidth=linewidth)
for i in all_failing_nodes_idxs
    x, y = remove_zero_tail!(deepcopy(sol.t), deepcopy(map(first, sol(sol.t, idxs=vidxs(i, :f)))))
    color = node_colors[findfirst(x -> x == i, all_failing_nodes_idxs)]
    lines!(ax11, x, y; color=color, linewidth=linewidth) # consider plotdensity=Int(1e5)
    scatter!(ax11, (x[end], y[end]); color=color, marker=:star5, markersize=25)
end
xlims!(ax11, 0, 1.1)
ax11.yticks = [-0.03, -0.02, -0.01, 0.00, 0.01, 0.02, 0.03];
ax11.xticks = [0.0, 0.25, 0.5, 0.75, 1.0]


# frequencies of failed gen nodes I=3.0 s^2
sol = sol1723
fig[1,2] = ax12 = Axis(fig; title=L"Inertia $I=3.0$ $s^2$", titlealign = :left, titlesize = titlesize)
hlines!(ax12, 0.03; color=:black, linestyle=:dash, linewidth=linewidth)
hlines!(ax12, -0.03; color=:black, linestyle=:dash, linewidth=linewidth)
for i in all_failing_nodes_idxs
    x, y = remove_zero_tail!(deepcopy(sol.t), deepcopy(map(first, sol(sol.t, idxs=vidxs(i, :f)))))
    color = node_colors[findfirst(x -> x == i, all_failing_nodes_idxs)]
    lines!(ax12, x, y; label="Node $i", color=color, linewidth=linewidth) # consider plotdensity=Int(1e5)
    scatter!(ax12, (x[end], y[end]); color=color, marker=:star5, markersize=25)
end
xlims!(ax12, 0, 4.1)
ax12.xticks = [0.1, 1, 2, 3, 4];
ax12.yticks = [-0.03, -0.02, -0.01, 0.00, 0.01, 0.02, 0.03];


# frequencies of failed gen nodes I=30.0 s^2
sol = sol1728
fig[1,3] = ax13 = Axis(fig; title=L"Inertia $I=30.0$ $s^2$", titlealign = :left, titlesize = titlesize)
hlines!(ax13, 0.03; color=:black, linestyle=:dash, linewidth=linewidth)
hlines!(ax13, -0.03; color=:black, linestyle=:dash, linewidth=linewidth)
for i in all_failing_nodes_idxs
    x, y = remove_zero_tail!(deepcopy(sol.t), deepcopy(map(first, sol(sol.t, idxs=vidxs(i, :f)))))
    color = node_colors[findfirst(x -> x == i, all_failing_nodes_idxs)]
    lines!(ax13, x, y; label="Node $i", color=color, linewidth=linewidth) # consider plotdensity=Int(1e5)
    scatter!(ax13, (x[end], y[end]); color=color, marker=:star5, markersize=25)
end
xlims!(ax13, -1., 34.5)
ax13.yticks = [-0.03, -0.02, -0.01, 0.00, 0.01, 0.02, 0.03];
ax13.xticks = [0.0, 5., 10., 15., 20., 25., 30., 35.];
axislegend(ax11, position = (0.99,0.85), labelsize=fontsize-4, patchsize = patchsize)


# FLOWS ########################################################################

# rating
df_config = DataFrame(CSV.File(joinpath(RESULTS_DIR, exp_name_date, "config.csv")))
network = import_system_wrapper(df_config, 1)
rating = get_prop(network, Graphs.edges(network), :rating)[1]

# failed power flows I=0.2 s^2
sol = sol1720
fig[2,1] = ax21 = Axis(fig; ylabel="Apparent power flow [p.u.]")
hlines!(ax21, rating; color=:black, linestyle=:dot, linewidth=linewidth, label="line rating")
scatter!(ax21, (NaN, NaN); label="line initially removed", color=line_colors[2], marker=:star5, markersize=25)
for i in all_failing_lines_idxs
    x, y = remove_zero_tail!(deepcopy(sol.t), deepcopy(map(first, sol(sol.t, idxs=eidxs(i, :S)))))
    color = line_colors[findfirst(x -> x == i, all_failing_lines_idxs)]
    lines!(ax21, x, y; color=color, linewidth=linewidth) # consider plotdensity=Int(1e5)
    scatter!(ax21, (x[end], y[end]); color=color, marker=:star5, markersize=25)

end
xlims!(ax21, 0, 1.1)
ax21.xticks = [0.0, 0.25, 0.5, 0.75, 1.0]


# failed power flows I=3.0 s^2
sol = sol1723
fig[2,2] = ax22 = Axis(fig; xlabel="Time [s]")
hlines!(ax22, rating; color=:black, linestyle=:dot, linewidth=linewidth)
for i in all_failing_lines_idxs
    x, y = remove_zero_tail!(deepcopy(sol.t), deepcopy(map(first, sol(sol.t, idxs=eidxs(i, :S)))))
    color = line_colors[findfirst(x -> x == i, all_failing_lines_idxs)]
    lines!(ax22, x, y; label="Node $i", color=color, linewidth=linewidth) # consider plotdensity=Int(1e5)
    scatter!(ax22, (x[end], y[end]); color=color, marker=:star5, markersize=25)
end
xlims!(ax22, 0, 4.1)


# failed power flows I=30.0 s^2
sol = sol1728
fig[2,3] = ax23 = Axis(fig)
hlines!(ax23, rating; color=:black, linestyle=:dot, linewidth=linewidth)
for i in all_failing_lines_idxs
    x, y = remove_zero_tail!(deepcopy(sol.t), deepcopy(map(first, sol(sol.t, idxs=eidxs(i, :S)))))
    color = line_colors[findfirst(x -> x == i, all_failing_lines_idxs)]
    lines!(ax23, x, y; label="Node $i", color=color, linewidth=linewidth) # consider plotdensity=Int(1e5)
    scatter!(ax23, (x[end], y[end]); color=color, marker=:star5, markersize=25)
end
xlims!(ax23, -1., 34.5)
ax23.xticks = [0.0, 5., 10., 15., 20., 25., 30., 35.];
axislegend(ax21, position = (0.99,0.88), labelsize=fontsize-4, patchsize = patchsize)

linkyaxes!(ax11,ax12,ax13)
linkyaxes!(ax21,ax22,ax23)
linkxaxes!(ax11, ax21); linkxaxes!(ax12, ax22); linkxaxes!(ax13, ax23);
# hideydecorations!(ax12); hideydecorations!(ax13); hideydecorations!(ax22); hideydecorations!(ax23)
ax12.yticklabelsvisible = false; ax13.yticklabelsvisible = false; ax22.yticklabelsvisible = false; ax23.yticklabelsvisible = false;

ax11.xticklabelsvisible = false; ax12.xticklabelsvisible = false; ax13.xticklabelsvisible = false


CairoMakie.save(joinpath(MA_DIR, "WS_traj,lines+nodes,task_ids=$task_id_array,init_fail=$initial_fail,_new_ND.pdf"),fig)
# CairoMakie.save(joinpath(MA_DIR, "WS_traj,lines+nodes,task_ids=$task_id_array,init_fail=$initial_fail,_new_ND.png"),fig)
fig
