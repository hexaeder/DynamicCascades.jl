includet(abspath(@__DIR__, "..", "..", "..","helpers_jarray.jl"))

using DynamicCascades
using NetworkDynamics
using Serialization

###
### inspect solutions
###
using NetworkDynamicsInspector
using WGLMakie # for inspector
# using Bonito # for using plot pane and memorizing plots
# Bonito.set_cleanup_time!(720)


function inspect_failing(dir, task_id, initial_fail)
    sol = Serialization.deserialize(joinpath(dir, "task_id=$task_id,initial_fail=$initial_fail.sol"));
    nw = NetworkDynamics.extract_nw(sol)

    vindices = [i for i in 1:nv(nw) if sol(sol.t[end], idxs=vidxs(i, :node_swing_stat))[1] != 1]
    eindices = [i for i in 1:ne(nw) if sol(sol.t[end], idxs=eidxs(i, :line_stat))[1] == 0]
    
    # set marker for failed vertices
    for i in vindices
        set_marker!(nw[VIndex(i)], :xcross)
    end

    inspect(sol; restart=true, reset=true)
    set_sol!(sol) # optional if after inspect(sol)
    set_graphplot!(; nstate=[:ω], estate=[:S], nstate_rel=false, estate_rel=false)
    define_timeseries!([
        (; selcomp=[VIndex(i) for i in vindices], states=[:ω, :ωmax], rel=false),
        (; selcomp=[EIndex(i) for i in eindices], states=[:S, :rating], rel=false),
    ])

    return sol
end

# dump_app_state()

## plot only failing vs all
# using CairoMakie
# CairoMakie.activate!()
# df_config = DataFrame(CSV.File(joinpath(exp_data_dir, "config.csv")));
# plot_simulation(dir, dir, df_config, task_id, initial_fail)

# # all lines and nodes 
# define_timeseries!([
#     (; selcomp=[VIndex(i) for i in 1:nv(nw)], states=[:ω, :ωmax], rel=false),
#     (; selcomp=[EIndex(i) for i in 1:ne(nw)], states=[:S, :rating], rel=false),
# ])

# Load the data
exp_name_date = "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250321_171511.976"
exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)


###
### "line graph"
###
sub_dir = "general_investigations"
initial_fail = 35
task_id = 68 # 140
# sol = inspect_failing(exp_data_dir, sub_dir, task_id, initial_fail)

# ωmax = 0.251
# Find steady state...
# Shutdown line 35 at t=0.1
# Vertex 58 tripped at t=1.6799790081136325
sol = simulate(exp_name_date, 68, initial_fail;
    gen_model=SwingDynLoadModel,
    tspan=(0., 40.),
    solverargs = (;dtmax=0.01),
    verbose = true);

# ωmay = 0.503
# Find steady state...
# Shutdown line 35 at t=0.1
# Line 142 tripped at t=2.3283051398978425
# Vertex 58 tripped at t=2.941264178376854
sol2 = simulate(exp_name_date, 140, initial_fail;
    gen_model=SwingDynLoadModel,
    tspan=(0., 40.),
    solverargs = (;dtmax=0.01),
    verbose = true);   

inspect(sol)
set_state!(; t=0.06009656910048066, tmin=0.0, tmax=4.810359414780986)
set_graphplot!(; nstate=[:ω], estate=[:S], nstate_rel=false, estate_rel=false, ncolorrange=(-0.25132743f0, 0.25132743f0), ecolorrange=(0.0f0, 1.8353579f0))
define_timeseries!([
    (; selcomp=[VIndex(58), VIndex(60)], states=[:ω, :ωmax], rel=false),
    (; selcomp=[EIndex(35), EIndex(142)], states=[:S, :rating], rel=false),
    ])
    
set_sol!(sol) # optional if after inspect(sol)



exp_name_date = "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250321_171511.976"
exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
sub_dir = "trajectories_MA_after_up2"
dir = joinpath(exp_data_dir, sub_dir)

sol = simulate(exp_name_date, task_id, initial_fail;
    gen_model=SwingDynLoadModel,
    tspan=(0., 40.),
    solverargs = (;dtmax=0.01),
    verbose = true);
Serialization.serialize(joinpath(dir, "task_id=$task_id,initial_fail=$initial_fail.sol"), sol)

"""
See org/Trajectories WS-Networks
"""


###
### inertia variation (task_id_array = [1720, 1723, 1728])
###
sub_dir = "trajectories_MA"
initial_fail = 78
# low inertia
"""
 - Adjacent nodes of initially triggered line fail soon.
 - There are 3 other nodes directly connected to one of the failing nodes, two of them fail, a third one doesn't.
"""
task_id = 1720 
inspect_failing(dir, task_id, initial_fail)
# high inertia
"""
 - Due to inertia, the two adjacent nodes (see I=0.2) don't fail.
 - Slowly line loadings increase more and more leading eventually to multiple line failures that then cause node failures again.
"""
task_id = 1728 
inspect_failing(dir, task_id, initial_fail)


###
### f_b variation 1 
###
exp_name_date = "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250321_171511.976"
exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
sub_dir = "different_node_failure_models"
node_model = "BH+Pmech(standard_model)_after_up"
dir = joinpath(exp_data_dir, sub_dir, node_model)
initial_fail = 95

# narrow bound
task_id = 2265
sol = simulate(exp_name_date, task_id, initial_fail;
    gen_model=SwingDynLoadModel,
    tspan=(0., 40.),
    solverargs = (;dtmax=0.01),
    verbose = true);
Serialization.serialize(joinpath(dir, "task_id=$task_id,initial_fail=$initial_fail.sol"), sol)
sol = Serialization.deserialize(joinpath(dir, "task_id=$task_id,initial_fail=$initial_fail.sol"));
nw = NetworkDynamics.extract_nw(sol)
sol = inspect_failing(dir, task_id, initial_fail)
set_sol!(sol) # optional if after inspect(sol)
set_state!(; t=1.6726953360654804, tmin=0.0, tmax=11.306688830104303)
define_timeseries!([
    (; selcomp=[VIndex(33), VIndex(35), VIndex(37), VIndex(72)], states=[:ω, :ωmax], rel=false),
    (; selcomp=[EIndex(3), EIndex(88), EIndex(89), EIndex(90), EIndex(91), EIndex(95), EIndex(97)], states=[:S, :rating], rel=false),
])
set_graphplot!(; nstate=[:ω], estate=[:S], nstate_rel=false, estate_rel=false, ncolorrange=(-0.50265485f0, 0.50265485f0), ecolorrange=(0.0f0, 2.095634f0))


# wider bound
task_id = 2310
sol = inspect_failing(dir, task_id, initial_fail)
set_sol!(sol) # optional if after inspect(sol)
set_state!(; t=4.926033627210094, tmin=0.0, tmax=10.487046632124352)
set_graphplot!(; nstate=[:ω], estate=[:P], nstate_rel=false, estate_rel=false, ncolorrange=(-0.28274333f0, 0.28274333f0), ecolorrange=(-1.9670131f0, 1.9670131f0))
define_timeseries!([
    (; selcomp=[VIndex(33), VIndex(35), VIndex(37), VIndex(72)], states=[:ω, :ωmax], rel=false),
    (; selcomp=[EIndex(3), EIndex(88), EIndex(89), EIndex(90), EIndex(91), EIndex(95), EIndex(97), EIndex(94)], states=[:S, :rating], rel=false),
])


###
### different node (and line) failure models
###

## no failures after initial perturbation
task_id = 2265
initial_fail = 95
exp_name_date = "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250321_171511.976"
exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
sub_dir = "different_node_failure_models"
node_model = "no_failures"
dir = joinpath(exp_data_dir, sub_dir, node_model)
sol = simulate(exp_name_date, task_id, initial_fail;
    gen_model=SwingDynLoadModel,
    trip_lines=:none,
    trip_nodes=:none,
    tspan=(0., 40.),
    solverargs = (;dtmax=0.01),
    verbose = true);
Serialization.serialize(joinpath(dir, "task_id=$task_id,initial_fail=$initial_fail.sol"), sol)

sol = Serialization.deserialize(joinpath(dir, "task_id=$task_id,initial_fail=$initial_fail.sol"));
nw = NetworkDynamics.extract_nw(sol)
inspect(sol; restart=true, reset=true)
set_sol!(sol) # optional if after inspect(sol)
define_timeseries!([
    (; selcomp=[VIndex(33), VIndex(35), VIndex(37), VIndex(72)], states=[:ω, :ωmax], rel=false),
    (; selcomp=[EIndex(3), EIndex(88), EIndex(89), EIndex(90), EIndex(91), EIndex(95), EIndex(97), EIndex(94)], states=[:S, :rating], rel=false),
])


## change_to_BH_only
initial_fail = 95
exp_name_date = "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250321_171511.976"
exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
sub_dir = "different_node_failure_models"
node_model = "change_to_BH_only"
dir = joinpath(exp_data_dir, sub_dir, node_model)
# sol = simulate(exp_name_date, task_id, initial_fail;
#     gen_model=SwingDynLoadModel_change_to_BH_only,
#     tspan=(0., 40.),
#     solverargs = (;dtmax=0.01),
#     verbose = true);
# Serialization.serialize(joinpath(dir, "task_id=$task_id,initial_fail=$initial_fail.sol"), sol)

# narrow bounds
task_id = 2265
sol = Serialization.deserialize(joinpath(dir, "task_id=$task_id,initial_fail=$initial_fail.sol"));
nw = NetworkDynamics.extract_nw(sol)
inspect(sol; restart=true, reset=true)
set_sol!(sol) # optional if after inspect(sol)
set_state!(; t=0.0, tmin=0.0, tmax=14.847645429362881)
define_timeseries!([
    (; selcomp=[VIndex(4), VIndex(33), VIndex(35), VIndex(37), VIndex(72)], states=[:ω, :ωmax], rel=false),
    (; selcomp=[EIndex(3), EIndex(6), EIndex(11), EIndex(12), EIndex(88), EIndex(89), EIndex(90), EIndex(91), EIndex(94), EIndex(95), EIndex(97)], states=[:S, :rating], rel=false),
])

# wide bounds
task_id = 2310
inspect_failing(dir, task_id, initial_fail);
dump_app_state()



## change_Pmech_only # BUG Macht keinen Sinn, weil ω nicht auf null gesetzt wird im CB?
task_id = 2265
initial_fail = 95
exp_name_date = "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250321_171511.976"
exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
sub_dir = "different_node_failure_models"
node_model = "change_Pmech_only"
dir = joinpath(exp_data_dir, sub_dir, node_model)
# sol = simulate(exp_name_date, task_id, initial_fail;
#     gen_model=SwingDynLoadModel_change_Pmech_only,
#     tspan=(0., 40.),
#     solverargs = (;dtmax=0.01),
#     verbose = true);
# Serialization.serialize(joinpath(dir, "task_id=$task_id,initial_fail=$initial_fail.sol"), sol)

# wide bounds
task_id = 2310
sol = Serialization.deserialize(joinpath(dir, "task_id=$task_id,initial_fail=$initial_fail.sol"));
nw = NetworkDynamics.extract_nw(sol)
inspect(sol; restart=true, reset=true)
set_sol!(sol) # optional if after inspect(sol)
define_timeseries!([
    (; selcomp=[VIndex(33), VIndex(35), VIndex(37), VIndex(72)], states=[:ω, :ωmax], rel=false),
    (; selcomp=[EIndex(3), EIndex(88), EIndex(89), EIndex(90), EIndex(91), EIndex(95), EIndex(97), EIndex(94)], states=[:S, :rating], rel=false),
])

# narrow bounds
task_id = 2265
sol = Serialization.deserialize(joinpath(dir, "task_id=$task_id,initial_fail=$initial_fail.sol"));
nw = NetworkDynamics.extract_nw(sol)
inspect(sol; restart=true, reset=true)
set_sol!(sol) # optional if after inspect(sol)
define_timeseries!([
    (; selcomp=[VIndex(33), VIndex(35), VIndex(37), VIndex(72)], states=[:ω, :ωmax], rel=false),
    (; selcomp=[EIndex(3), EIndex(88), EIndex(89), EIndex(90), EIndex(91), EIndex(95), EIndex(97), EIndex(94)], states=[:S, :rating], rel=false),
])

###
### Node failure and inertia reduction leads to another line failure
###
# wider bound 
task_id = 2310
exp_name_date = "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250321_171511.976"
exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
sub_dir = "different_node_failure_models"
node_model = "BH+Pmech(standard_model)"
dir = joinpath(exp_data_dir, sub_dir, node_model)
initial_fail = 95
sol = Serialization.deserialize(joinpath(dir, "task_id=$task_id,initial_fail=$initial_fail.sol"));
nw = NetworkDynamics.extract_nw(sol)
set_sol!(sol) # optional if after inspect(sol)
set_state!(; t=4.458474846100491, tmin=1.5438983304386071, tmax=6.035238928078193)
set_graphplot!(; nstate=[:ω], estate=[:P], nstate_rel=false, estate_rel=false, ncolorrange=(-0.28274333f0, 0.28274333f0), ecolorrange=(-1.9670131f0, 1.9670131f0))
define_timeseries!([
    (; selcomp=[VIndex(33), VIndex(35), VIndex(37)], states=[:ω, :ωmax], rel=false),
    (; selcomp=[EIndex(3), EIndex(88), EIndex(89), EIndex(90), EIndex(91), EIndex(95), EIndex(97), EIndex(94)], states=[:S, :rating], rel=false),
])


###
### f_b variation 2 (for wider f_b the line adjacent to the failing node fails as well, for narrower bounds this line does not fail.)
###
exp_name_date = "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250321_171511.976"
exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
sub_dir = "general_investigations"
initial_fail = 35
dir = joinpath(exp_data_dir, sub_dir)

# narrow bound
task_id = 68
inspect_failing(dir, task_id, initial_fail)

# wider bound
task_id = 140
inspect_failing(dir, task_id, initial_fail)



###
### f_b variation 3 (for wider f_b the line adjacent to the failing node fails as well, this line triggers large cascade, 
### Pmech of the failing node (v33) does not have effect (v33 is isolated after line failure))
exp_name_date = "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250321_171511.976"
exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
sub_dir = "trajectories_braess"
initial_fail = 93
dir = joinpath(exp_data_dir, sub_dir)

# narrow bound
task_id = 1446
inspect_failing(dir, task_id, initial_fail)

# wider bound
task_id = 1617
inspect_failing(dir, task_id, initial_fail)



task_id = 1446
initial_fail = 2
sol = Serialization.deserialize(joinpath(exp_data_dir, "trajectories_braess", "task_id=$task_id,initial_fail=$initial_fail.sol"));
vindices = [i for i in 1:100 if sol(sol.t[end], idxs=vidxs(i, :node_swing_stat))[1] != 1]
eindices = [i for i in 1:200 if sol(sol.t[end], idxs=eidxs(i, :line_stat))[1] == 0]

inspect(sol; restart=true, reset=true)
set_sol!(sol) # optional if after inspect(sol)
set_graphplot!(; nstate=[:ω], estate=[:S], nstate_rel=false, estate_rel=false)




