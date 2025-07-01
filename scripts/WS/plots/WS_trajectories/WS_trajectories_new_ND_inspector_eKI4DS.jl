"""
"""



using DynamicCascades
using NetworkDynamics
using NetworkDynamicsInspector
using WGLMakie # for inspector
# # using Bonito # for using plot pane and memorizing plots
# # Bonito.set_cleanup_time!(720)
using Serialization

initial_fail = 78
task_id = 1728
task_id_array = [1720, 1723, 1728]

exp_name_date = "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250321_171511.976"
exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)

# load solution objects
sol1720 = Serialization.deserialize(joinpath(exp_data_dir, "trajectories", "task_id=1720,initial_fail=78_eKI4DS_example.sol"));
sol1723 = Serialization.deserialize(joinpath(exp_data_dir, "trajectories", "task_id=1723,initial_fail=78_eKI4DS_example.sol"));
sol1728 = Serialization.deserialize(joinpath(exp_data_dir, "trajectories", "task_id=1728,initial_fail=78_eKI4DS_example.sol"));

# 78: initial trigger line
lines_task_id_1720 = [78] # sol1720.failures.saveval
lines_task_id_1723 = [78] # sol1723.failures.saveval
lines_task_id_1728 = [78, 199, 194, 131, 147, 146, 14, 106, 135, 149, 148] # sol1728.failures.saveval
nodes_task_id_1720 = [29, 98, 92, 58] # sol1720.failures_nodes.saveval
nodes_task_id_1723 = [29, 98] # sol1723.failures_nodes.saveval
nodes_task_id_1728 = [98, 59, 62, 61, 60] # sol1728.failures_nodes.saveval

all_failing_lines_idxs = [14, 78, 106, 131, 135, 146, 147, 148, 149, 194, 199]
all_failing_nodes_idxs = [29, 58, 59, 60, 61, 62, 92, 98]

# call for interactive inspection
inspect(sol1728; restart=true, reset=true)
# dump_app_state()

# initial state
set_sol!(sol1728) # optional if after inspect(sol)
set_state!(; t=0.0, tmin=0.0, tmax=35.0)
set_graphplot!(; nstate=[:ω], estate=[:S], nstate_rel=false, estate_rel=false, ncolorrange=(-0.18849556f0, 0.18849556f0), ecolorrange=(0.0f0, 2.099972f0))
define_timeseries!([
    (; selcomp=[], states=Symbol[], rel=false),
])

# cascade
set_sol!(sol1728) # optional if after inspect(sol)
set_state!(; t=0.02875924404272802, tmin=0.0, tmax=35.0)
set_graphplot!(; nstate=[:ω], estate=[:S], nstate_rel=false, estate_rel=false, ncolorrange=(-0.18849556f0, 0.18849556f0), ecolorrange=(0.0f0, 2.099972f0))
define_timeseries!([
    (; selcomp=[VIndex(98), VIndex(59), VIndex(60), VIndex(61), VIndex(62)], states=[:ω, :ωmax], rel=false),
    (; selcomp=[EIndex(78), EIndex(199), EIndex(194), EIndex(131), EIndex(147), EIndex(14), EIndex(106), EIndex(135), EIndex(149), EIndex(148)], states=[:S, :rating], rel=false),
])