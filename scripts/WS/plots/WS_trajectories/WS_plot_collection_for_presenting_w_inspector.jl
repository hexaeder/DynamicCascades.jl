include(abspath(@__DIR__, "..", "..", "..","helpers_jarray.jl"))

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
dir = joinpath(exp_data_dir, sub_dir)
initial_fail = 35
task_id = 140

# sol = simulate(exp_name_date, task_id, initial_fail;
#     gen_model=SwingDynLoadModel,
#     verbose = true);
# Serialization.serialize(joinpath(dir, "task_id=$task_id,initial_fail=$initial_fail.sol"), sol)
# task_id = 68
# ωmax = 0.251
# Shutdown line 35 at t=0.1
# Vertex 58 tripped at t=1.6799804655331712
# Terminated on steady state at 110.88068502409718
# ---
# task_id = 140
# ωmay = 0.503
# Shutdown line 35 at t=0.1
# Line 142 tripped at t=2.328306445619492
# Vertex 58 tripped at t=2.9412654518564687
# Terminated on steady state at 120.92007289341954

task_id = 68
sol68 = Serialization.deserialize(joinpath(dir, "task_id=$task_id,initial_fail=$initial_fail.sol"));
task_id = 140
sol140 = Serialization.deserialize(joinpath(dir, "task_id=$task_id,initial_fail=$initial_fail.sol"));
 
inspect_wrapper(sol68; which_trajectories=:failing)
inspect_wrapper(sol140; which_trajectories=:failing)

set_state!(; t=0.06009656910048066, tmin=0.0, tmax=4.0)
define_timeseries!([
    (; selcomp=[VIndex(58), VIndex(60)], states=[:ω, :ωmax], rel=false),
    (; selcomp=[EIndex(35), EIndex(142)], states=[:S, :rating], rel=false),
    ])


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
inspect_wrapper(dir, task_id, initial_fail)
# high inertia
"""
 - Due to inertia, the two adjacent nodes (see I=0.2) don't fail.
 - Slowly line loadings increase more and more leading eventually to multiple line failures that then cause node failures again.
"""
task_id = 1728 
inspect_wrapper(dir, task_id, initial_fail)


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
sol = inspect_wrapper(dir, task_id, initial_fail)
set_sol!(sol) # optional if after inspect(sol)
set_state!(; t=1.6726953360654804, tmin=0.0, tmax=11.306688830104303)
define_timeseries!([
    (; selcomp=[VIndex(33), VIndex(35), VIndex(37), VIndex(72)], states=[:ω, :ωmax], rel=false),
    (; selcomp=[EIndex(3), EIndex(88), EIndex(89), EIndex(90), EIndex(91), EIndex(95), EIndex(97)], states=[:S, :rating], rel=false),
])
set_graphplot!(; nstate=[:ω], estate=[:S], nstate_rel=false, estate_rel=false, ncolorrange=(-0.50265485f0, 0.50265485f0), ecolorrange=(0.0f0, 2.095634f0))


# wider bound
task_id = 2310
sol = inspect_wrapper(dir, task_id, initial_fail)
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
inspect_wrapper(dir, task_id, initial_fail);
dump_app_state()



## change_Pmech_only
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
node_model = "full_failure"
dir = joinpath(exp_data_dir, sub_dir, node_model)
initial_fail = 95

# sol = simulate(exp_name_date, task_id, initial_fail;
#     gen_model=SwingDynLoadModel,
#     verbose = true);
# Serialization.serialize(joinpath(dir, "task_id=$task_id,initial_fail=$initial_fail.sol"), sol)

t_stop_CBs = 4.85
# sol_t_stop_CBs = simulate(exp_name_date, task_id, initial_fail;
#     gen_model=SwingDynLoadModel,
#     tspan=(0., 100.),
#     t_stop_CBs = t_stop_CBs,
#     verbose = true);
# Serialization.serialize(joinpath(dir, "task_id=$task_id,initial_fail=$initial_fail,t_stop_CBs=$t_stop_CBs.sol"), sol_t_stop_CBs)

sol = Serialization.deserialize(joinpath(dir, "task_id=$task_id,initial_fail=$initial_fail.sol"));
sol_t_stop_CBs = Serialization.deserialize(joinpath(dir, "task_id=$task_id,initial_fail=$initial_fail,t_stop_CBs=$t_stop_CBs.sol"));

CairoMakie.save(joinpath(MA_DIR, "braessness", "task_id=$task_id,initial_fail=$initial_fail,t_stop_CBs=$t_stop_CBs.pdf"), plot_simulation(sol_t_stop_CBs))
CairoMakie.save(joinpath(MA_DIR, "braessness", "task_id=$task_id,initial_fail=$initial_fail.pdf"), plot_simulation(sol))

nw = NetworkDynamics.extract_nw(sol)
inspect(sol; restart=true, reset=true)
set_state!(; t=4.458474846100491, tmin=1.5438983304386071, tmax=6.035238928078193)
set_graphplot!(; nstate=[:ω], estate=[:P], nstate_rel=false, estate_rel=false, ncolorrange=(-0.28274333f0, 0.28274333f0), ecolorrange=(-1.9670131f0, 1.9670131f0))
define_timeseries!([
    (; selcomp=[VIndex(33), VIndex(35), VIndex(37)], states=[:ω, :ωmax], rel=false),
    (; selcomp=[EIndex(3), EIndex(88), EIndex(89), EIndex(90), EIndex(91), EIndex(95), EIndex(97), EIndex(94)], states=[:S, :rating], rel=false),
])


inspect_wrapper(sol_t_stop_CBs; which_trajectories = :all)
# set_sol!(sol_t_stop_CBs)

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
inspect_wrapper(dir, task_id, initial_fail)

# wider bound
task_id = 140
inspect_wrapper(dir, task_id, initial_fail)



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
inspect_wrapper(dir, task_id, initial_fail)

# wider bound
task_id = 1617
inspect_wrapper(dir, task_id, initial_fail)



task_id = 1446
initial_fail = 2
sol = Serialization.deserialize(joinpath(exp_data_dir, "trajectories_braess", "task_id=$task_id,initial_fail=$initial_fail.sol"));
vindices = [i for i in 1:100 if sol(sol.t[end], idxs=vidxs(i, :node_swing_stat))[1] != 1]
eindices = [i for i in 1:200 if sol(sol.t[end], idxs=eidxs(i, :line_stat))[1] == 0]

inspect(sol; restart=true, reset=true)
set_sol!(sol) # optional if after inspect(sol)
set_graphplot!(; nstate=[:ω], estate=[:S], nstate_rel=false, estate_rel=false)


##################################################################################
######## Braessness vs Braessness: Example trajectories (lines + nodes) ##########
##################################################################################

###
### example_trajectories RTS #####################################################
###

# choose inertia
I = 3 # scaling factor

# choose `f_b_narrow` and `f_b_wide`
f_b_narrow = 0.3
f_b_wide = 1.0
N_nodes = 73 # nv(import_system(:rtsgmlc))
failures = "line+node"

### plot coordinates in order to find ArrayTaskID and initially_failed_line
###


## full failure
exp_name_date = "RTS_exp04_variation_frequency+inertia_PIK_HPC_20250616_213442.721"
# write_failures_f_b_to_df(exp_name_date, freq_bounds, I)
# write_network_measures_to_df(exp_name_date)
braessness_lines, braessness_nodes = get_braessness(exp_name_date, f_b_narrow, f_b_wide, I);
x_ρ, x_dist = get_network_measures(exp_name_date);
exp_nr_full = exp_name_date[8:9]

## inertia failure
exp_name_date = "RTS_exp05_variation_frequency+inertia_inertia_failure_PIK_HPC_20250701_174359.744"
# write_failures_f_b_to_df(exp_name_date, freq_bounds, I)
# write_network_measures_to_df(exp_name_date)
braessness_lines_BH, braessness_nodes_BH = get_braessness(exp_name_date, f_b_narrow, f_b_wide, I);
x_ρ_BH, x_dist_BH = get_network_measures(exp_name_date);
exp_nr_BH = exp_name_date[8:9]

## power failure
exp_name_date = "RTS_exp06_variation_frequency+inertia_power_failure_PIK_HPC_20250701_175708.899"
# write_failures_f_b_to_df(exp_name_date, freq_bounds, I)
# write_network_measures_to_df(exp_name_date)
braessness_lines_Pmech, braessness_nodes_Pmech = get_braessness(exp_name_date, f_b_narrow, f_b_wide, I);
x_ρ_Pmech, x_dist_Pmech = get_network_measures(exp_name_date);
exp_nr_Pmech = exp_name_date[8:9]

prefix = exp_name_date[1:7]
exp_nrs = "$prefix,$exp_nr_full,$exp_nr_BH,$exp_nr_Pmech"


xs = braessness_lines_Pmech .+ braessness_nodes_Pmech;
ys = braessness_lines_BH .+ braessness_nodes_BH;
zs = braessness_lines .+ braessness_nodes;

fig = plot_braessness_vs_braessness_2D_count(xs, ys, zs, failures, f_b_narrow, f_b_wide, I; markersize=5, plot_coordinates = true)
# CairoMakie.save(joinpath(MA_DIR, "braessness", "exp=$exp_nrs,braessness_vs_braessness_lines_and_nodes_logcount_colorscale,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I,coordinates_plotted.pdf"),fig)
# CairoMakie.save(joinpath(MA_DIR, "braessness", "exp=$exp_nrs,braessness_vs_braessness_lines_and_nodes_logcount_colorscale,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I,coordinates_plotted.png"),fig)


exp_name_date_full = "RTS_exp04_variation_frequency+inertia_PIK_HPC_20250616_213442.721"
exp_name_date_inertia = "RTS_exp05_variation_frequency+inertia_inertia_failure_PIK_HPC_20250701_174359.744"
exp_name_date_power = "RTS_exp06_variation_frequency+inertia_power_failure_PIK_HPC_20250701_175708.899"


exp_name_date = exp_name_date_inertia
task_id = 790
initial_fail = 27

gen_model = Serialization.deserialize(joinpath(RESULTS_DIR, exp_name_date, "exp.params"))[:gen_model]
dir = joinpath(RESULTS_DIR, exp_name_date, "example_trajectories", "sims", "gen_model=$gen_model,task_id=$task_id,initial_fail=$initial_fail.sol")
sol = Serialization.deserialize(dir);
inspect_wrapper(sol; which_trajectories = :failing);
describe_failures(sol)
# set_sol!(sol)

"""
full
 - wide: some islands with ∑P_i < 0
 - narrow: islands with ∑P_i > 0
inertia
 - wide: 
 - narrow:
power
 - wide: 
 - narrow: strong ∑P_i < 0
"""
b_power_failure = 6
b_inertia_failure = 15
# x-Axis: wide:    ArrayTaskID=790, initially_failed_line=27 
#         narrow:  ArrayTaskID=212, initially_failed_line=27
# y-Axis: wide:    ArrayTaskID=790, initially_failed_line=27 
#         narrow:  ArrayTaskID=212, initially_failed_line=27
# z-Axis: wide:    ArrayTaskID=790, initially_failed_line=27, Braessness=17 
#         narrow:  ArrayTaskID=212, initially_failed_line=27




###
### example_trajectories WS #################################################
###


#= 
## API
# sol-Object
describe_failures(sol)

# Inspector
set_state!(; t=0.0, tmin=0.0, tmax=tmax)
dump_app_state()
=#


I = 7.5
f_b_narrow = 0.03
f_b_wide = 0.15

### plot coordinates in order to find ArrayTaskID and initially_failed_line
###

## full failure
exp_name_date = "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250321_171511.976"
# write_failures_f_b_to_df(exp_name_date, freq_bounds, I)
# write_network_measures_to_df(exp_name_date)
braessness_lines, braessness_nodes = get_braessness(exp_name_date, f_b_narrow, f_b_wide, I);
x_ρ, x_dist = get_network_measures(exp_name_date);
exp_nr_full = exp_name_date[11:12]

## inertia failure
exp_name_date = "WS_k=4_exp11_vary_I_only_lines_and_nodes_change_to_BH_complement_f_bPIK_HPC_K_=3,N_G=32_20250420_145304.04"
# write_failures_f_b_to_df(exp_name_date, freq_bounds, I)
# write_network_measures_to_df(exp_name_date)
braessness_lines_BH, braessness_nodes_BH = get_braessness(exp_name_date, f_b_narrow, f_b_wide, I);
x_ρ_BH, x_dist_BH = get_network_measures(exp_name_date);
exp_nr_BH = exp_name_date[11:12]

## power failure
exp_name_date = "WS_k=4_exp12_vary_I_only_lines_and_nodes_change_Pmech_complement_f_b_PIK_HPC_K_=3,N_G=32_20250718_160303.054"
write_failures_f_b_to_df(exp_name_date, freq_bounds, I)
# write_network_measures_to_df(exp_name_date)
braessness_lines_Pmech, braessness_nodes_Pmech = get_braessness(exp_name_date, f_b_narrow, f_b_wide, I);
# x_ρ_Pmech, x_dist_Pmech = get_network_measures(exp_name_date);
exp_nr_Pmech = exp_name_date[11:12]

prefix = exp_name_date[1:10]
exp_nrs = "$prefix,$exp_nr_full,$exp_nr_BH,$exp_nr_Pmech"


xs = braessness_lines_Pmech .+ braessness_nodes_Pmech
ys = braessness_lines_BH .+ braessness_nodes_BH
zs = braessness_lines .+ braessness_nodes

fig = plot_braessness_vs_braessness_2D_count(xs, ys, zs, failures, f_b_narrow, f_b_wide, I; markersize=5, plot_coordinates = true)
# CairoMakie.save(joinpath(MA_DIR, "braessness", "exp=$exp_nrs,braessness_vs_braessness_lines_and_nodes_logcount_colorscale,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I,coordinates_plotted.pdf"),fig)
# CairoMakie.save(joinpath(MA_DIR, "braessness", "exp=$exp_nrs,braessness_vs_braessness_lines_and_nodes_logcount_colorscale,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I,coordinates_plotted.png"),fig)



exp_name_date_full = "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250321_171511.976"
exp_name_date_inertia = "WS_k=4_exp11_vary_I_only_lines_and_nodes_change_to_BH_complement_f_bPIK_HPC_K_=3,N_G=32_20250420_145304.04"
exp_name_date_power = "WS_k=4_exp12_vary_I_only_lines_and_nodes_change_Pmech_complement_f_b_PIK_HPC_K_=3,N_G=32_20250718_160303.054"


function WS_inspect_wrapper(exp_name_date, task_id, initial_fail;
        start_inspector= true, kwargs...)

    if exp_name_date[11:12] == "04"
        gen_model = SwingDynLoadModel
    elseif exp_name_date[11:12] == "11"
        gen_model = SwingDynLoadModel_change_to_BH_only
    elseif exp_name_date[11:12] == "12"
        gen_model = SwingDynLoadModel_change_Pmech_only
    end

    sol = Serialization.deserialize(joinpath(RESULTS_DIR, exp_name_date, "example_trajectories", "sims", "gen_model=$gen_model,task_id=$task_id,initial_fail=$initial_fail.sol"));
    start_inspector ? inspect_wrapper(sol; kwargs...) : nothing
    return sol
end



""" --- 01 --- high power, inertia, full (53,54,48) """
# b_power_failure = 53
# b_inertia_failure = 54
# get_task_ids_and_failed_line_from_coordinates(exp_name_date_power, exp_name_date_inertia, exp_name_date_full, b_power_failure, b_inertia_failure, f_b_narrow, f_b_wide, I;
#     run_save_simulation = true,
#     generate_plot = true)

function adapt_inspector01()
    set_state!(; t=0.0, tmin=0.0, tmax=40)
    define_timeseries!([
        (; selcomp=[VIndex(3), VIndex(4), VIndex(5), VIndex(6), VIndex(7), VIndex(8), VIndex(72), VIndex(73), VIndex(77), VIndex(82)], states=[:ω, :ωmax], rel=false),
        (; selcomp=[EIndex(9), EIndex(12), EIndex(13), EIndex(14), EIndex(15), EIndex(16)], states=[:S, :rating], rel=false)])
end

""" power: 53 """
# x-Axis: narrow:  ArrayTaskID=2229, initially_failed_line=15
#          wide:   ArrayTaskID=2256, initially_failed_line=15
sol = WS_inspect_wrapper(exp_name_date_power, 2229, 15; which_trajectories=:failing); adapt_inspector01(); describe_failures(sol) # narrow
sol = WS_inspect_wrapper(exp_name_date_power, 2256, 15; which_trajectories=:failing, tmax=40); describe_failures(sol) # wide

""" inertia 54 """
# y-Axis: narrow:  ArrayTaskID=2229, initially_failed_line=15
#          wide:   ArrayTaskID=2256, initially_failed_line=15
sol = WS_inspect_wrapper(exp_name_date_inertia, 2229, 15; which_trajectories=:failing); adapt_inspector01(); describe_failures(sol) # narrow
sol = WS_inspect_wrapper(exp_name_date_inertia, 2256, 15; which_trajectories=:failing, tmax=40); describe_failures(sol)  # wide

""" full 48 """
# z-Axis: narrow:  ArrayTaskID=7584, initially_failed_line=15
#          wide:   ArrayTaskID=7755, initially_failed_line=15, Braessness=48
sol = WS_inspect_wrapper(exp_name_date_full, 7584, 15; which_trajectories=:failing); adapt_inspector01(); describe_failures(sol)  # narrow
sol = WS_inspect_wrapper(exp_name_date_full, 7755, 15; which_trajectories=:failing, tmax=40); describe_failures(sol) # wide


"""
Tweaked power injections
"""
initial_fail = 15
exp_name_date = exp_name_date_full
# task_id = 7584 # narrow
task_id = 7755 # wide 

# exp_name_date = exp_name_date_inertia
# task_id = 2229 # narrow
# task_id = 2256 # wide

# exp_name_date = exp_name_date_power
# task_id = 2229 # narrow
# task_id = 2256 # wide

df_config = DataFrame(CSV.File(joinpath(RESULTS_DIR, exp_name_date, "config.csv")))
_,_,_,_,_,_,_,_,_,_,_,_,freq_bound,trip_lines,trip_nodes,_,_ = get_network_args_stripped(df_config, task_id)
network = import_system_wrapper(df_config, task_id)
df_config[!, :filepath_graph] = replace.(df_config[!, :filepath_graph],"/home/brandner" => "/home/brandner/nb_data/HU_Master/2122WS/MA")
graph = loadgraph(df_config[task_id,:filepath_graph])

# Adapt v5
Pmech_v5 = get_prop(network, 5, :Pmech)
Pload_v5 = get_prop(network, 5, :Pload)
Pmech_v5_new = Pmech_v5 - 2*Pmech_v5
Pload_v5_new = Pload_v5 - Pmech_v5
set_prop!(network, 5, :Pmech, Pmech_v5_new)
set_prop!(network, 5, :Pload, Pload_v5_new)

# Restore global power balance
v_restore_id = 35
Pmech_v_restore = get_prop(network, v_restore_id, :Pmech)
Pload_v_restore = get_prop(network, v_restore_id, :Pload)
Pmech_v_restore_new = Pmech_v_restore + 2*Pmech_v5
Pload_v_restore_new = Pload_v_restore + Pmech_v5
set_prop!(network, v_restore_id, :Pmech, Pmech_v_restore_new)
set_prop!(network, v_restore_id, :Pload, Pload_v_restore_new)

if exp_name_date[11:12] == "04"
    gen_model = SwingDynLoadModel
elseif exp_name_date[11:12] == "11"
    gen_model = SwingDynLoadModel_change_to_BH_only
elseif exp_name_date[11:12] == "12"
    gen_model = SwingDynLoadModel_change_Pmech_only
end

# Check global power balance
# sum(get_prop(network, 1:nv(network), :Pmech))
# sum(get_prop(network, 1:nv(network), :Pload))
@assert isapprox(
    sum(get_prop(network, 1:nv(network), :Pmech)),
    sum(get_prop(network, 1:nv(network), :Pload)), atol=1e-8)

# This checks power balance again
sol = simulate(network;
    graph=graph,
    gen_model = gen_model,
    x_static= steadystate(network; graph=graph, zeroidx=1),
    initial_fail=initial_fail,
    trip_lines=trip_lines,
    trip_nodes=trip_nodes,
    freq_bound=freq_bound
    );

###
### set node positions in inspector
### 
nw = extract_nw(sol)
exp_name_date_full = "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250321_171511.976"
df = CSV.read(joinpath(RESULTS_DIR, exp_name_date_full, "graph_positions_spring.csv"), DataFrame)
positions = collect(zip(df.x, df.y))

for i in 1:nv(nw)
    set_position!(nw, VIndex(i), positions[i])
end
inspect_wrapper(sol; which_trajectories = :failing); # adapt_inspector01()


""" --- 02 --- full > power, inertia (31,18,41) """
# get_task_ids_and_failed_line_from_coordinates(exp_name_date_power, exp_name_date_inertia, exp_name_date_full, 31, 18, f_b_narrow, f_b_wide, I; run_save_simulation = false, generate_plot = false)
function adapt_inspector02()
    set_state!(; t=0.0, tmin=0.0, tmax=30)
    define_timeseries!([
        (; selcomp=[VIndex(1), VIndex(21), VIndex(31), VIndex(33), VIndex(36), VIndex(38), VIndex(39), VIndex(47), VIndex(51), VIndex(52), VIndex(57), VIndex(59), VIndex(61), VIndex(74), VIndex(77), VIndex(78), VIndex(79), VIndex(80), VIndex(99)], states=[:ω, :ωmax], rel=false),
        (; selcomp=[EIndex(1), EIndex(2), EIndex(3), EIndex(34), EIndex(35), EIndex(60), EIndex(69), EIndex(86), EIndex(92), EIndex(96), EIndex(98), EIndex(103), EIndex(126), EIndex(127), EIndex(136), EIndex(174), EIndex(175), EIndex(177), EIndex(178)], states=[:S, :rating], rel=false)])
end

""" power: 31
 - narrow: local power disbalance mitigated => trajectories return to fixpoint
 - wide: mainly further nodes fail
"""
# x-Axis: wide:    ArrayTaskID=150, initially_failed_line=2 
#         narrow:  ArrayTaskID=123, initially_failed_line=2
sol = WS_inspect_wrapper(exp_name_date_power, 123, 2; which_trajectories=:failing); adapt_inspector02(); describe_failures(sol) # narrow
sol = WS_inspect_wrapper(exp_name_date_power, 150, 2; which_trajectories=:failing, tmax=50); describe_failures(sol) # wide

""" inertia 18
 - narrow:
 - wide:
"""
# y-Axis: wide:    ArrayTaskID=150, initially_failed_line=2 
#         narrow:  ArrayTaskID=123, initially_failed_line=2
sol = WS_inspect_wrapper(exp_name_date_inertia, 123, 2; which_trajectories=:failing); adapt_inspector02(); describe_failures(sol) # narrow
sol = WS_inspect_wrapper(exp_name_date_inertia, 150, 2; which_trajectories=:failing, tmax=50); describe_failures(sol) # wide

""" full 41
 - narrow: 
 - wide: 
"""
# z-Axis: wide:    ArrayTaskID=501, initially_failed_line=2, Braessness=41 
#         narrow:  ArrayTaskID=330, initially_failed_line=2
sol = WS_inspect_wrapper(exp_name_date_full, 330, 2; which_trajectories=:failing); adapt_inspector02(); describe_failures(sol) # narrow
sol = WS_inspect_wrapper(exp_name_date_full, 501, 2; which_trajectories=:failing, tmax=50); describe_failures(sol) # wide



"""
"""
b_inertia_failure = -39
b_full_failure = 17
""" inertia -39
 - narrow:
 - wide:
"""
# x-Axis: wide:    ArrayTaskID=474, initially_failed_line=79 
#         narrow:  ArrayTaskID=447, initially_failed_line=79
WS_inspect_wrapper(exp_name_date_inertia, 447, 79; which_trajectories=:failing); # narrow
WS_inspect_wrapper(exp_name_date_inertia, 474, 79; which_trajectories=:all); # wide

""" full 17
 - narrow:
 - wide: 
"""
# y-Axis: wide:    ArrayTaskID=1617, initially_failed_line=79 
#         narrow:  ArrayTaskID=1446, initially_failed_line=79
WS_inspect_wrapper(exp_name_date_full, 1446, 79; which_trajectories=:failing); # narrow
WS_inspect_wrapper(exp_name_date_full, 1617, 79; which_trajectories=:all);  # wide

""" power: 19
 - narrow: 
 - wide: 
"""
# z-Axis: wide:    ArrayTaskID=474, initially_failed_line=79, Braessness=19 
#         narrow:  ArrayTaskID=447, initially_failed_line=79
WS_inspect_wrapper(exp_name_date_power, 447, 79; which_trajectories=:failing); # narrow
WS_inspect_wrapper(exp_name_date_power, 474, 79; which_trajectories=:all); # wide



###
### For comparing with `scripts/WS_snapshot_networks_plot.jl`
###

initial_fail = 15
task_id = 7584
exp_name_date = "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250321_171511.976"
sol = simulate(exp_name_date, task_id, initial_fail;
    gen_model=SwingDynLoadModel,
    verbose = true);


## power: 53
# x-Axis: narrow:  ArrayTaskID=2229, initially_failed_line=15
#          wide:   ArrayTaskID=2256, initially_failed_line=15
inspect_wrapper(sol; which_trajectories = :failing); adapt_inspector01(); describe_failures(sol) # narrow
sol = WS_inspect_wrapper(exp_name_date_power, 2256, 15; which_trajectories=:failing, tmax=40); describe_failures(sol) # wide

## inertia 54
# y-Axis: narrow:  ArrayTaskID=2229, initially_failed_line=15
#          wide:   ArrayTaskID=2256, initially_failed_line=15
sol = WS_inspect_wrapper(exp_name_date_inertia, 2229, 15; which_trajectories=:failing); adapt_inspector01(); describe_failures(sol) # narrow
sol = WS_inspect_wrapper(exp_name_date_inertia, 2256, 15; which_trajectories=:failing, tmax=40); describe_failures(sol)  # wide

## full 48
# z-Axis: narrow:  ArrayTaskID=7584, initially_failed_line=15
#          wide:   ArrayTaskID=7755, initially_failed_line=15, Braessness=48
sol = WS_inspect_wrapper(exp_name_date_full, 7584, 15; which_trajectories=:failing); adapt_inspector01(); describe_failures(sol)  # narrow
sol = WS_inspect_wrapper(exp_name_date_full, 7755, 15; which_trajectories=:failing, tmax=40); describe_failures(sol) # wide



