include(abspath(@__DIR__, "..", "helpers_jarray.jl"))


using DynamicCascades
using NetworkDynamics
using Serialization


# Load the data
exp_name_date = "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250321_171511.976"
exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)


###
### get large diff in number of failures for different f_b
###
df_all_failures = DataFrame(CSV.File(joinpath(exp_data_dir, "all_failures.csv")))
select!(df_all_failures, Not([:τ, :K, :α, :failure_modes, :filepath_steady_state, :init_pert, :k, :N_nodes, :σ, :μ]))

# Filter only for inertia_values = 7.5
df = filter(:inertia_values => ==(7.5), df_all_failures)
df = filter(:freq_bounds => in([0.03, 0.15]), df)


# Compute the sum of failures for each row
df[!, :sum_avg_failures] .= df.avg_line_failures .+ df.avg_node_failures;
df[!, :diff] .= NaN;

# Calculate the difference of failures for different f_b
for i in 1:size(df, 1)
    if iseven(i)
        df[i, :diff] = df[i, :sum_avg_failures] - df[i-1, :sum_avg_failures]
    end
end


# calculate number of line and node failures at end of cascade
N_nodes = nv(nw)
N_lines = ne(nw)
vindices = [i for i in 1:N_nodes if sol(sol.t[end], idxs=vidxs(i, :node_swing_stat))[1] == 0]
eindices = [i for i in 1:N_lines if sol(sol.t[end], idxs=eidxs(i, :line_stat))[1] == 0]
number_of_node_failures = length(vindices)
number_of_line_failures = length(eindices) - 1

###
### rerun multiple simulations
###
task_id_array = [1446, 1617]
initial_fail_array = [1, 2, 15, 93]
for task_id in task_id_array
    for initial_fail in initial_fail_array
        println("task_id=$task_id, initial_fail=$initial_fail")
        sol, nw = simulate(exp_data_dir, task_id, initial_fail;
            solverargs = (;dtmax=0.01),
            verbose = true);
        
        Serialization.serialize(joinpath(exp_data_dir, "trajectories_braess", "task_id=$task_id,initial_fail=$initial_fail.sol"), sol)
        Serialization.serialize(joinpath(exp_data_dir, "trajectories_braess", "task_id=$task_id,initial_fail=$initial_fail.nw"), nw)
    end
end

###
### inspect solutions
###
using NetworkDynamicsInspector
using WGLMakie # for inspector
# # using Bonito # for using plot pane and memorizing plots
# # Bonito.set_cleanup_time!(720)

task_id = 1446
initial_fail = 15
sol = Serialization.deserialize(joinpath(exp_data_dir, "trajectories_braess", "task_id=$task_id,initial_fail=$initial_fail.sol"));
inspect(sol; restart=true, reset=true)
set_sol!(sol) # optional if after inspect(sol)
set_graphplot!(; nstate=[:ω], estate=[:S], nstate_rel=false, estate_rel=false)
vindices = [i for i in 1:N_nodes if sol(sol.t[end], idxs=vidxs(i, :node_swing_stat))[1] == 0]
eindices = [i for i in 1:N_lines if sol(sol.t[end], idxs=eidxs(i, :line_stat))[1] == 0]
define_timeseries!([
    (; selcomp=[VIndex(i) for i in vindices], states=[:ω, :ωmax], rel=false),
    (; selcomp=[EIndex(i) for i in eindices], states=[:S, :rating], rel=false),
])